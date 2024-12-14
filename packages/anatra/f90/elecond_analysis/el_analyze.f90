!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_ctrl
  use mod_xtcio
  use mod_dcdio
  use mod_traj
  use mod_com
  use xdr, only: xtcfile

  ! subroutines
  !
  public  :: analyze

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, trajopt, traj)
!-----------------------------------------------------------------------
      implicit none

      ! formal arguments
      type(s_input),  intent(in)    :: input 
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_trajopt),intent(in)    :: trajopt 
      type(s_traj),   intent(inout) :: traj

      ! local variables
      type(s_dcd)            :: dcd
      type(xtcfile)          :: xtc

      integer                :: iunit, io, itraj, ixyz
      integer                :: istep, jstep, ij, tab, istep_tot, iatm, jatm
      integer                :: nmol, natm, natm_full, nstep, nstep_tot
      integer                :: nend, nt_range
      integer                :: trajtype
      integer                :: ithr, nthreads
      real(8)                :: prd, rit(3), ri0(3), rjt(3), rj0(3)
      real(8)                :: dri(3), drj(3)
      logical                :: is_end
      character(len=MaxChar) :: fname

      type(s_traj)           :: qt 
      real(8), allocatable   :: tcf_thr(:, :), tcf(:), tcf_red(:)
      real(8), allocatable   :: count(:), count_red(:)


      nt_range = option%nt_range
      !
      !
      nthreads = 1
      !$omp parallel
      !$ nthreads = omp_get_num_threads()
      !$omp end parallel

      ! Get trajectory type
      !
      call get_trajtype(input%ftraj(1), trajtype)

      ! Get Total number of steps
      !
      if (trajtype == TrajTypeDCD) then
        call get_total_step_from_dcd(input%ftraj, nstep_tot)
      else if (trajtype == TrajTypeXTC) then
        call get_total_step_from_xtc(input%ftraj, nstep_tot)
      end if

      call setup_traj_from_args(trajopt, nstep_tot, qt, trajid = 1)

      natm  = qt%natm
      nstep = nstep_tot

      ! Get Trajectory weighted with charge 
      !
      write(iw,*)
      write(iw,'("Analyze> Load trajectories")')
      istep_tot = 0
      do itraj = 1, input%ntraj

        call open_trajfile(input%ftraj(itraj), trajtype, iunit, dcd, xtc)
        call init_trajfile(trajtype, iunit, dcd, xtc, natm_full)

        is_end = .false.
        istep  = 0
        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          call read_trajfile_oneframe(trajtype, iunit, istep, dcd, xtc, is_end)

          if (is_end) exit

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ", i0)') istep_tot
          end if

          call send_coord_to_traj(1, trajtype, dcd, xtc, traj)
          qt%coord(:, :, istep_tot) = traj%coord(:, :, 1) 

        end do

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, iunit, dcd, xtc)

      end do

      do istep = 1, nstep_tot
        do iatm = 1, natm 
          qt%coord(1:3, iatm, istep) = qt%charge(iatm) * qt%coord(1:3, iatm, istep)
        end do
      end do

      ! Calculate TCF or Displacement 
      !
      write(iw,*)
      write(iw,'("Analyze> Calculate TCF")')

      allocate(tcf_thr(0:option%nt_range, 0:nthreads - 1), tcf(0:option%nt_range))
      allocate(tcf_red(0:option%nt_range))
      allocate(count(0:option%nt_range))
      allocate(count_red(0:option%nt_range))

      tcf_thr = 0.0d0
      tcf     = 0.0d0
      tcf_red = 0.0d0
      count   = 0.0d0

      if (option%calctype == CalcTypeCOORD) then
        ithr = 0
        tcf_thr = 0.0d0

        do istep = 1, nstep
          write(iw,'("Step : ", i0)') istep

          nend = min(istep + nt_range * option%nt_sparse, nstep_tot) 

          tcf_red = 0.0d0
          !$omp parallel private(iatm, jatm, jstep, tab, prd,          &
          !$omp        &         ri0, rit, rj0, rjt, dri, drj),        &
          !$omp        & shared (qt),                                  & 
          !$omp        & default(shared),                              & 
          !$omp        & reduction(+:tcf_red), reduction(+:count_red)
          !$omp do
          do iatm = 1, natm
            ri0(1:3) = qt%coord(1:3, iatm, istep)

            ! Self-part
            !
            !do jstep = istep, min(istep + option%tcfrange + 1, nstep_tot)

            tab =  - 1
            do jstep = istep, nend, option%nt_sparse
              !tab = jstep - istep
              tab = tab + 1

              rit(1:3) = qt%coord(1:3, iatm, jstep)
              dri(1:3) = rit(1:3) - ri0(1:3)
              prd      = dot_product(dri(:), dri(:))
           
              !tcf_thr(tab, ithr) = tcf_thr(tab, ithr) + prd
              !tcf(tab) = tcf(tab) + prd
              tcf_red(tab)    = tcf_red(tab)   + prd
              count_red(tab)  = count_red(tab) + 1.0d0
            end do

            ! Distinct-part
            !
            if (iatm < natm) then
              do jatm = iatm + 1, natm
                tab = - 1
                rj0(1:3) = qt%coord(1:3, jatm, istep)
             
                !do jstep = istep, min(istep + option%tcfrange + 1, nstep_tot)
                do jstep = istep, nend, option%nt_sparse
             
                  !tab      = jstep - istep
                  tab      = tab + 1
             
                  rit(1:3) = qt%coord(1:3, iatm, jstep)
                  rjt(1:3) = qt%coord(1:3, jatm, jstep)
             
                  dri(1:3) = rit(1:3) - ri0(1:3)
                  drj(1:3) = rjt(1:3) - rj0(1:3)
             
                  prd = dot_product(dri(:), drj(:))
             
                  tcf_red(tab) = tcf_red(tab) + 2.0d0 * prd
                end do
             
              end do
            end if
          end do
          !$omp end do
          !$omp end parallel

          tcf = tcf + tcf_red

        end do

        !do istep = 0, option%tcfrange
        !  !tcf(istep) = sum(tcf_thr(istep, 0:nthreads - 1))
        !  write(iw,'("Step = ", i0)') istep
        !  do ithr = 0, nthreads - 1
        !    write(iw,'("Threads tch_thr = ", i5, e15.7)') ithr, tcf_thr(istep, ithr)
        !    tcf(istep) = tcf(istep) + tcf_thr(istep, ithr)
        !  end do
        !end do

        !do ithr = 0, nthreads - 1
        !  write(iw,'("THR ",i5,2x,e15.7)') ithr, tcf_thr(0, ithr)
        !end do 

      else if (option%calctype == CalcTypeVEL) then
        write(iw,'("calctype = VEL is not supported yet. Sorry.")')
        stop
      end if

      ! Average
      !
      do istep = 0, option%nt_range
        jstep = istep * option%nt_sparse
        tcf(istep) = tcf(istep) / (nstep_tot - jstep)
      end do

      ! Output
      !
      write(fname,'(a,".eltcf")') trim(output%fhead)
      call open_file(fname, io)
      do istep = 0, option%nt_range
        write(io,'(2e15.7)') option%dt_out * istep, tcf(istep)
      end do
      close(io)

      ! Memory deallocation
      !
      deallocate(tcf)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
