!=======================================================================
module mod_tcf
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_cv
  use mod_bootstrap
  use mod_random
  use mod_analyze_str

  ! constants
  !

  ! structures
  !

  ! subroutines
  !
  public :: get_transtcf

  contains

!-----------------------------------------------------------------------
    subroutine get_transtcf(option, cvinfo, state) 
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_cvinfo), intent(in)    :: cvinfo
      type(s_state),  intent(inout) :: state(cvinfo%nfile)

      integer :: ifile, istep, jstep, it
      integer :: is, js
      integer :: nstep, nt_range, nt, nstate

      real(8), allocatable :: hist(:, :, :), norm(:, :)


      nstep    = state(1)%nstep
      nt_range = option%nt_range
      nstate   = option%nstate

      ! Memory allocation (if needed) 
      !
      do ifile = 1, cvinfo%nfile
        if (.not. allocated(state(ifile)%hist)) then
          allocate(state(ifile)%hist(0:nt_range, 0:nstate, 0:nstate))
          allocate(state(ifile)%norm(0:nt_range, 0:nstate))
          state(ifile)%hist = 0.0d0
        end if
      end do

      allocate(hist(0:nt_range, 0:nstate, 0:nstate))
      allocate(norm(0:nt_range, 0:nstate))
      hist = 0.0d0
      norm = 0.0d0

      ! Calculate tcf
      !
      do ifile = 1, cvinfo%nfile
        hist = 0.0d0
        norm = 0.0d0
        do istep = 1, nstep - nt_range - 1

          is = state(ifile)%data(istep)
          nt = min(istep + nt_range, nstep)

          it = - 1
          do jstep = istep, nt 
            it = it + 1
            js = state(ifile)%data(jstep)
            hist(js, is, it) = hist(js, is, it) + 1.0d0
            norm(is, it)     = norm(is, it)     + 1.0d0 
          end do

        end do

        state(ifile)%hist(:, :, :) = hist(:, :, :)
        state(ifile)%norm(:, :)    = norm(:, :)

      end do

      ! Memory deallocation
      !
      deallocate(hist, norm)

    end subroutine get_transtcf 
!-----------------------------------------------------------------------


end module mod_tcf
!=======================================================================
