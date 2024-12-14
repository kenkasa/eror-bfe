!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj
  use mod_com

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj(2)


      integer                :: trajtype, iunit, itraj
      integer                :: i, ig, igc, imol
      integer                :: istep, istep_tot, nstep_tot
      integer                :: natm, ngrid, nmol(2)
      integer                :: io_ts
      real(8)                :: sumd, zave, dlen
      real(8)                :: dg, xsta, cossq, sorder_ave
      character(len=MaxChar) :: fout
      logical                :: is_end

      type(s_com)   :: com(2)
      type(s_dcd)   :: dcd
      type(xtcfile) :: xtc

      real(8), allocatable :: dp(:, :), angle(:), cost(:), theta(:)
      real(8), allocatable :: sorder(:)
      real(8), allocatable :: distr(:)


      if (option%xcoord == XcoordModeCOST) then
        dg    = option%dcost
        ngrid = 2.0d0 / dble(option%dcost)
        xsta  = -1.0d0
      else if (option%xcoord == XcoordModeTHETA) then
        dg    = option%dtheta
        ngrid = 180.0d0 / dble(option%dtheta)
        xsta  = 0.0d0
      end if

      ! Get trajectory types
      !
      call get_trajtype(input%ftraj(1), trajtype)


      ! Setup molecule
      !
      do i = 1, 2
        call get_com(option%mode(i),       &
                     traj(i),              &
                     com(i),               &
                     setup = .true.,       &
                     calc_coord = .false., &
                     myrank = 0)

        nmol(i) = com(i)%nmol
      end do

      if (nmol(1) /= nmol(2)) then
        write(iw,'("Analyze> Error.")')
        write(iw,'("Number of molecules in two dcd should be the same.")')
        stop
      end if

      ! allocate memory
      !

      allocate(dp(3, nmol(1)), angle(nmol(1)), cost(nmol(1)))
      allocate(theta(nmol(1)), sorder(nmol(1)))
      allocate(distr(0:ngrid))

      ! calculate histogram
      !
      write(iw,*)
      write(iw,'("Analyze> Start orientation analysis")')
      !
      distr      = 0.0d0
      sorder_ave = 0.0d0
      istep_tot  = 0

      write(fout,'(a,".angle")') trim(output%fhead)
      open(UnitOUT, file=trim(fout))
      do itraj = 1, input%ntraj
        call open_trajfile(input%ftraj(itraj), trajtype, iunit, dcd, xtc)
        call init_trajfile(trajtype, iunit, dcd, xtc, natm)

        is_end = .false.
        istep  = 0

        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          call read_trajfile_oneframe(trajtype, iunit, istep, dcd, xtc, is_end)

          if (is_end) exit

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ",i0)') istep_tot
          end if

          do i = 1, 2
            call send_coord_to_traj(1, trajtype, dcd, xtc, traj(i))
            call get_com(option%mode(i),       &
                         traj(i),              &
                         com(i),               &
                         setup = .false.,      &
                         calc_coord = .true.,  &
                         myrank = 1)
          end do

!$omp parallel private(imol, dlen, zave, cossq), &
!$omp        & shared(dp, cost, theta, angle),   &
!$omp        & default(shared),                  &
!$omp        & reduction(+:sorder_ave)
!$omp do
!
          do imol = 1, nmol(1)
            dp(:, imol) = com(2)%coord(:, imol, 1) &
                        - com(1)%coord(:, imol, 1)

            dlen = sqrt(dot_product(dp(1:3, imol), dp(1:3, imol)))

            if (option%judgeup == JudgeUpModeNMOLUP) then
              if (imol <= option%nmolup) then
                cost(imol) =  dp(3, imol) / dlen 
              else
                cost(imol) = -dp(3, imol) / dlen
              end if
            else if (option%judgeup == JudgeUpModeCOORD) then
              zave = com(1)%coord(3, imol, 1) + com(2)%coord(3, imol, 1)
              zave = 0.5d0 * zave

              if (zave >= 0.0d0) then
                cost(imol) =  dp(3, imol) / dlen
              else
                cost(imol) = -dp(3, imol) / dlen
              end if
            else
              cost(imol) =  dp(3, imol) / dlen 
            endif

            theta(imol) = acos(cost(imol)) / PI * 180.0d0

            cossq = cost(imol) ** 2
            sorder(imol) = 1.5d0 * cossq - 0.5d0

            sorder_ave = sorder_ave + sorder(imol)

            if (option%xcoord == XcoordModeCOST) then
              angle(imol) = cost(imol)
            else if (option%xcoord == XcoordModeTHETA) then
              angle(imol) = theta(imol)
            end if
          end do
!$omp     end do
!$omp     end parallel

          do imol = 1, nmol(1)
            ig  = nint((angle(imol) - xsta)/ dg)
            igc = (angle(imol) - xsta)/ dg

            if (ig < 0 .or. ig > ngrid) then
              write(iw,'("istep = ", i0)') istep_tot
              write(iw,'("imol  = ", i0)') imol
              write(iw,'("ig = ",i0)') ig
              write(iw,'("cost  = ",f15.7)') cost(imol)
            else
              if (igc == ngrid - 1) then
                distr(igc)     = distr(igc)     + 0.5d0
                distr(igc + 1) = distr(igc + 1) + 0.5d0
              else
                distr(ig) = distr(ig) + 1.0d0
              end if
            end if
            
          end do

          ! Print out Time-series
          !
          write(UnitOUT,'(f20.10)', advance='no') traj(1)%dt * istep_tot
          do imol = 1, nmol(1) 
            write(UnitOUT,'(f20.10)',advance='no') angle(imol)
          end do
          write(UnitOUT,*)

        end do

        istep_tot = istep_tot - 1
        call close_trajfile(trajtype, iunit, dcd, xtc)

      end do
      close(UnitOUT)

      nstep_tot = istep_tot

      do ig = 0, ngrid
        distr(ig) = distr(ig) / (nstep_tot * dg * nmol(1))
      end do

      !sorder_ave = sorder_ave / (nstep_tot * nmol(1))

      ! normalize
      !
      sumd = 0.0d0
      do ig = 0, ngrid
        sumd = sumd + distr(ig) * dg
      end do
      distr = distr / sumd

      !write(iw,'("S = (3/2)*<cos(T)**2> - 1/2 = ",f20.10)') sorder_ave
      !Note:  bugfix of sorder is not done yet
      

      ! generate files
      !
      ! - dipole
      !write(fout,'(a,".angle")') trim(output%fhead)
      !open(UnitOUT, file=trim(fout))
      !if (option%xcoord == XcoordModeCOST) then 
      !  do istep = 1, traj(1)%nstep 
      !    write(UnitOUT,'(f20.10)', advance='no') traj(1)%dt * istep
      !    do id = 1, nmol(1) 
      !      write(UnitOUT,'(f20.10)',advance='no') cost(id, istep)
      !    end do
      !    write(UnitOUT,*)
      !  end do
      !else if (option%xcoord == XcoordModeTHETA) then
      !  do istep = 1, traj(1)%nstep 
      !    write(UnitOUT,'(f20.10)', advance='no') traj(1)%dt * istep
      !    do id = 1, nmol(1) 
      !      write(UnitOUT,'(f20.10)',advance='no') theta(id, istep)
      !    end do
      !    write(UnitOUT,*)
      !  end do
      !end if
      !close(UnitOUT)
      !
      ! - cost distribution 
      !
      write(fout,'(a,".hist")') trim(output%fhead) 
      open(UnitOUT, file=trim(fout))
      do ig = 0, ngrid
        write(UnitOUT,'(2f20.10)') &
          xsta + dble(ig) * dg, distr(ig)
      end do
      close(UnitOUT)

      ! Deallocate memory
      !
      deallocate(dp, angle, cost, theta, distr)
      deallocate(sorder)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
