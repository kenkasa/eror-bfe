!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_xtcio
  use mod_dcdio
  use mod_traj
  use mod_com
  use xdr, only: xtcfile

  ! subroutines
  !
  public :: analyze
  public :: update_distr

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj

      type(s_dcd)   :: dcd
      type(xtcfile) :: xtc
      type(s_com)   :: com

      integer :: itraj, istep, imol
      integer :: istep_tot, ig
      integer :: nstep_tot, ngrid
      integer :: natm, nstep, nmol
      integer :: iunit
      integer :: trajtype
      logical :: is_end

      real(8), allocatable :: dipole(:, :, :)
      real(8), allocatable :: distr(:) 


      ! Setup
      !

      ! Get trajectory type
      !
      call get_trajtype(input%ftraj(1), trajtype)

      call get_com(option%mode,          &
                   traj,                 &
                   com,                  &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)

      nmol  = com%nmol
      ngrid = 2.0d0 / dble(option%dcost)

      allocate(distr(0:ngrid))
      distr = 0.0d0

      istep_tot = 0
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

          call send_coord_to_traj(1, trajtype, dcd, xtc, traj)

          call get_com(option%mode,         &
                       traj,                &
                       com,                 &
                       setup = .false.,     &
                       calc_coord = .true., & 
                       myrank = 1)

          call update_distr(option, traj, com, ngrid, distr)

        end do

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, iunit, dcd, xtc)

      end do

      nstep_tot = istep_tot
      distr(:)  = distr(:) / (nstep_tot * option%dcost * nmol)

      ! generate files
      !
      call open_file(trim(output%fdpori), iunit)
      do ig = 0, ngrid
        write(iunit,'(2f20.10)') &
          dble(ig) * option%dcost - 1.0d0, distr(ig)
      end do
      close(iunit)

      ! Deallocate memory
      !
      deallocate(distr)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine update_distr(option, traj, com, ngrid, distr)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_traj),   intent(in)    :: traj
      type(s_com),    intent(in)    :: com
      integer,        intent(in)    :: ngrid
      real(8),        intent(inout) :: distr(0:ngrid)

      integer              :: iatm, id, imol, ig, nmol
      real(8)              :: c, dlen, ct
      real(8), allocatable :: dp(:, :)


      nmol = com%nmol

      if (.not. allocated(dp)) &
        allocate(dp(3, nmol))

      dp = 0.0d0
      do iatm = 1, traj%natm
        c  = traj%charge(iatm)
        id = com%molid(iatm)

        dp(:, id) = dp(:, id) &
          + c * (traj%coord(:, iatm, 1) - com%coord(:, id, 1))
      end do

      do imol = 1, nmol
        dlen = sqrt(dot_product(dp(1:3, imol), dp(1:3, imol)))
        if (abs(dlen) < 1.0d-6) then
          write(iw,'("dlen = ", e20.10)') dlen
        end if

        ct = dp(3, imol) / dlen
        if (imol <= option%nmolup) then
          ct = ct
        else
          ct = -ct
        end if

        ig = nint((ct + 1.0d0) / option%dcost)

        if (ig >= 0 .and. ig <= ngrid) then
          distr(ig) = distr(ig) + 1.0d0
        end if
      end do

    end subroutine update_distr
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
