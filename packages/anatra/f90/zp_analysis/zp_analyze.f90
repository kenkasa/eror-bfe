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
      type(s_traj),   intent(inout) :: traj

      type(s_com)   :: com
      type(s_dcd)   :: dcd
      type(xtcfile) :: xtc

      integer :: iunit, trajtype, itraj
      integer :: istep, istep_tot, nstep_tot
      integer :: iatm, imol, ixyz, iz, ierr
      integer :: natm, nmol, nz
      real(8) :: boxave(3), vol, dv, zsta, zmin, zmax
      real(8) :: m, d(3), d2, z
      logical :: is_end

      real(8), allocatable :: weight(:) 
      real(8), allocatable :: gz(:) 


      allocate(weight(traj%natm))


      ! Get trajectory type
      !
      call get_trajtype(input%ftraj(1), trajtype)

      ! Setup molecule
      !
      call get_com(option%mode,          &
                   traj,                 &
                   com,                  &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)

      nmol = com%nmol

      ! Get weight
      !
      if (option%denstype == DensityTypeNUMBER) then
        weight = 1.0d0
      else if (option%denstype == DensityTypeELECTRON) then
        do iatm = 1, traj%natm
          weight(iatm) = - traj%charge(iatm) &
                         + get_atomicnum(traj%mass(iatm), ierr)
          if (ierr /= 0) then
            write(iw,'("Analyze> Error.")')
            write(iw,'("Failed to assign atomic number of ",a)') &
                    trim(traj%atmname(iatm))
          end if
        end do
      end if

      ! Get boxsize at first step for determining # of grids
      !
      call open_trajfile(input%ftraj(1), trajtype, iunit, dcd, xtc)
      call init_trajfile(trajtype, iunit, dcd, xtc, natm)
      call read_trajfile_oneframe(trajtype, iunit, 1, dcd, xtc, is_end)
      call send_coord_to_traj(1, trajtype, dcd, xtc, traj)

      boxave(:) = traj%box(:, 1)
      zmax      = boxave(3)
      zsta      = - 0.5d0 * zmax
      nz        = zmax / option%dz + 1

      call close_trajfile(trajtype, iunit, dcd, xtc)
      boxave(:) = 0.0d0

      ! Allocate memory
      !
      allocate(gz(0:nz))

      gz        = 0.0d0
      is_end    = .false.
      istep_tot = 0
      do itraj = 1, input%ntraj
        call open_trajfile(input%ftraj(itraj), trajtype, iunit, dcd, xtc)
        call init_trajfile(trajtype, iunit, dcd, xtc, natm)

        istep = 0
        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          call read_trajfile_oneframe(trajtype, iunit, istep, dcd, xtc, is_end)

          if (is_end) exit

          call send_coord_to_traj(1, trajtype, dcd, xtc, traj)

          call get_com(option%mode,         &
                       traj,                &
                       com,                 &
                       setup = .false.,     &
                       calc_coord = .true., & 
                       myrank = 1)

          boxave(:) = boxave(:) + traj%box(:, 1)

          do imol = 1, nmol
            if (option%centertype == CenterTypeZERO) then
              iz  = (com%coord(3, imol, 1) - zsta) / option%dz
            else if (option%centertype == CenterTypeHALF) then
              iz  = (com%coord(3, imol, 1)     &
                    - traj%box(3, 1) * 0.5d0   &
                    - zsta) / option%dz

            end if

            if (iz >= 0 .and. iz <= nz) &
              gz(iz) = gz(iz) + weight(imol)
          end do
        end do ! step

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, iunit, dcd, xtc)

      end do   ! traj

      nstep_tot = istep_tot

      boxave(:) = boxave(:) / dble(nstep_tot)
      dv        = option%dz * boxave(1) * boxave(2)
      vol       = boxave(1) * boxave(2) * boxave(3)

      gz = gz / (dv * nstep_tot)

      ! generate files
      !
      ! - Z-profile 
      open(10, file=trim(output%fzp))
      do iz = 0, nz
        z = zsta + option%dz * iz 
        write(10,'(2f20.10)') z, gz(iz) 
      end do
      close(10)

      deallocate(weight, gz)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
