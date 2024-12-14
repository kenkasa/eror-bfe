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
      type(s_traj),   intent(inout) :: traj(1)

      integer, parameter     :: UP = 1, LOW = 2

      integer                :: trajtype, iunit, itraj
      integer                :: i, ig, ig1, ig2, imol
      integer                :: istep, istep_tot, nstep_tot
      integer                :: natm, nmol(1)
      integer                :: ngrid
      integer                :: istate, nup, nlow, nup_tot, nlow_tot
      real(8)                :: zsta, dz, zave, sumd
      real(8)                :: zcom(2), dev(2), delta, delta_sq
      character(len=MaxChar) :: fout
      logical                :: is_end

      type(s_com)   :: com(1)
      type(s_dcd)   :: dcd
      type(xtcfile) :: xtc

      integer, allocatable :: state_uplow(:)
      real(8), allocatable :: distr_dev(:), distr_rmsd(:)


      zsta  = 0.0d0
      ngrid = option%ngrid
      dz    = option%dz

      ! Get trajectory types
      !
      call get_trajtype(input%ftraj(1), trajtype)


      ! Setup molecule
      !
      !do i = 1, 2
      i = 1
      call get_com(option%mode(i),       &
                   traj(i),              &
                   com(i),               &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)

      nmol(i) = com(i)%nmol
      !end do

      ! allocate memory
      !

      allocate(state_uplow(nmol(1)))
      allocate(distr_dev(0:ngrid), distr_rmsd(0:ngrid))

      ! calculate histogram
      !
      write(iw,*)
      write(iw,'("Analyze> Start Z-deviation analysis")')
      !
      distr_dev  = 0.0d0
      distr_rmsd = 0.0d0
      istep_tot  = 0

      nup_tot    = 0
      nlow_tot   = 0

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

          i = 1
          call send_coord_to_traj(1, trajtype, dcd, xtc, traj(i))
          call get_com(option%mode(i),       &
                       traj(i),              &
                       com(i),               &
                       setup = .false.,      &
                       calc_coord = .true.,  &
                       myrank = 1)

          ! Calculate Zcom
          !
          zcom        = 0.0d0
          nup         = 0
          nlow        = 0
          state_uplow = 0
          do imol = 1, nmol(1)
            
            if (option%judgeup == JudgeUpModeNMOLUP) then

              if (imol <= option%nmolup) then
                nup               = nup  + 1
                zcom(UP)          = zcom(UP)   + com(1)%coord(3, imol, 1)
                state_uplow(imol) = UP 
              else
                nlow              = nlow + 1
                zcom(LOW)  = zcom(LOW) + com(1)%coord(3, imol, 1)
                state_uplow(imol) = LOW 
              end if

            else if (option%judgeup == JudgeUpModeCOORD) then

              zave = com(1)%coord(3, imol, 1) 

              if (zave >= 0.0d0) then
                nup               = nup  + 1
                zcom(UP)          = zcom(UP)  + com(1)%coord(3, imol, 1)
                state_uplow(imol) = UP 
              else
                nlow              = nlow + 1
                zcom(LOW)         = zcom(LOW) + com(1)%coord(3, imol, 1) 
                state_uplow(imol) = LOW 
              end if
            else
              nup               = nup  + 1
              zcom(UP)          = zcom(UP)  + com(1)%coord(3, imol, 1)
              state_uplow(imol) = UP 
            endif

          end do

          nup_tot  = nup_tot  + nup
          nlow_tot = nlow_tot + nlow

          zcom(UP)  = zcom(UP)  / dble(nup)
          if (nlow /= 0) then
            zcom(LOW) = zcom(LOW) / dble(nlow)
          end if

          ! Calculate Histogram
          !
          dev = 0.0d0

          do imol = 1, nmol(1)
            istate   = state_uplow(imol)

            delta_sq = (com(1)%coord(3, imol, 1) - zcom(istate))**2
            delta    = sqrt(delta_sq)

            ! distr_dev
            !
            ig       = nint((delta - zsta) / dz)
            if (ig >= 0 .and. ig <= ngrid) then
              distr_dev(ig) = distr_dev(ig) + 1.0d0
            end if

            ! rmsd (histogram generated after)
            !
            dev(istate) = dev(istate) + delta_sq
          end do

          ! distr_rmsd
          !
          dev(UP)  = sqrt(dev(UP)  / dble(nup))
          dev(LOW) = sqrt(dev(LOW) / dble(nlow))
          
          ig1 = nint((dev(UP)  - zsta) / dz)
          ig2 = nint((dev(LOW) - zsta) / dz)
          if (ig1 >= 0 .and. ig1 <= ngrid) then
            distr_rmsd(ig1) = distr_rmsd(ig1) + 1.0d0
          end if

          if (ig2 >= 0 .and. ig2 <= ngrid) then 
            distr_rmsd(ig2) = distr_rmsd(ig2) + 1.0d0
          end if

        end do

        istep_tot = istep_tot - 1
        call close_trajfile(trajtype, iunit, dcd, xtc)

      end do

      nstep_tot = istep_tot

      do ig = 0, ngrid
        distr_dev(ig)  = distr_dev(ig)  / (nstep_tot * dz * nmol(1))
        distr_rmsd(ig) = distr_rmsd(ig) / (nstep_tot * dz)
      end do

      distr_dev(0)  = distr_dev(0)  * 2.0d0
      distr_rmsd(0) = distr_rmsd(0) * 2.0d0 

      !sorder_ave = sorder_ave / (nstep_tot * nmol(1))

      ! normalize
      !
      sumd = distr_dev(0) * 0.5d0 * dz
      do ig = 1, ngrid
        sumd = sumd + distr_dev(ig) * dz
      end do
      distr_dev = distr_dev / sumd

      sumd = distr_rmsd(0) * 0.5d0 * dz
      do ig = 1, ngrid
        sumd = sumd + distr_rmsd(ig) * dz
      end do
      distr_rmsd = distr_rmsd / sumd

      ! Output
      !
      write(fout,'(a,".distr_dev")') trim(output%fhead) 
      open(UnitOUT, file=trim(fout))
      do ig = 0, ngrid
        write(UnitOUT,'(2f20.10)') &
          zsta + dble(ig) * dz, distr_dev(ig)
      end do
      close(UnitOUT)

      write(fout,'(a,".distr_rmsd")') trim(output%fhead) 
      open(UnitOUT, file=trim(fout))
      do ig = 0, ngrid
        write(UnitOUT,'(2f20.10)') &
          zsta + dble(ig) * dz, distr_rmsd(ig)
      end do
      close(UnitOUT)

      ! Deallocate memory
      !
      deallocate(state_uplow)
      deallocate(distr_dev, distr_rmsd)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
