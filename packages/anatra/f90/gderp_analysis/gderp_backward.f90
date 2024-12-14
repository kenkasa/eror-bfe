!=======================================================================
module mod_gderp_backward
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_gderp_ctrl
  use mod_gderp_str

  ! structures
  !

  ! subroutines
  !
  public :: eval_kernel_backward
  public :: propagate_backward

  contains
!-----------------------------------------------------------------------
    subroutine eval_kernel_backward(gde)
!-----------------------------------------------------------------------
      implicit none

      type(s_gde), intent(inout) :: gde

      integer              :: i, j
      integer              :: nstep, ncut
      real(8)              :: dt
      real(8)              :: diff, cum, ri


      ! Setup variables
      !
      nstep = gde%nstep
      dt    = gde%dt
      ncut  = gde%ncut 

      ! Calculate Memory kernel
      !
      gde%kernel    = 0.0d0
      gde%kernel(1) = - gde%p(1) / (dt * gde%pdot(0))

      do i = 2, ncut
        cum = - gde%p(i) / dt

        ri  = 0.0d0
        do j = 1, i - 1
          ri = ri - gde%kernel(j) * gde%pdot(i - j)
        end do

        cum           = cum + ri
        gde%kernel(i) = cum / gde%pdot(0)
      end do

    end subroutine eval_kernel_backward
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine propagate_backward(gde)
!-----------------------------------------------------------------------
      implicit none

      type(s_gde), intent(inout) :: gde

      integer              :: i, j
      integer              :: nstep, ncut, nextend
      real(8)              :: dt
      real(8)              :: diff, cum, ri
      real(8)              :: k1, fact1, fact2

      ! Setup variables
      !
      nstep   = gde%nstep
      dt      = gde%dt
      ncut    = gde%ncut
      nextend = gde%nextend

      k1      = gde%kernel(1)
      fact1   = 1.0d0 / (1.0d0 + k1)
      fact2   = k1 / (1.0d0 + k1)

      ! Extend dP(t)/dt & P(t)
      !

      !   Initialize
      !
      gde%p_ext       = 0.0d0
      gde%pdot_ext    = 0.0d0

      gde%p_ext(0)    = gde%p(0)
      gde%pdot_ext(0) = gde%pdot(0) 

      !   Propagate
      !
      gde%p_ext(1) = fact2 * gde%p_ext(0)

      do i = 2, nextend

        cum = gde%p_ext(i - 1) * fact2

        ri  = 0.0d0
        do j = 2, min(i, ncut)
          ri = ri - dt * gde%kernel(j) * gde%pdot_ext(i - j)
        end do
        ri = ri * fact1 

        cum                 = cum + ri
        gde%p_ext(i)        = cum
        gde%pdot_ext(i - 1) = (gde%p_ext(i) - gde%p_ext(i - 1)) / dt 
      end do

    end subroutine propagate_backward
!-----------------------------------------------------------------------

end module mod_gderp_backward
!=======================================================================
