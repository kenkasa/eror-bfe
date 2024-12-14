!=======================================================================
module mod_gderp_forward
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
  public :: eval_kernel_forward
  public :: propagate_forward

  contains
!-----------------------------------------------------------------------
    subroutine eval_kernel_forward(gde)
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
      gde%kernel(0) = - gde%p(1) / (dt * gde%pdot(1))

      do i = 2, ncut + 1
        cum = - gde%p(i) / dt

        ri  = 0.0d0
        do j = 0, i - 2 
          ri = ri - gde%kernel(j) * gde%pdot(i - j)
        end do

        cum               = cum + ri
        gde%kernel(i - 1) = cum / gde%pdot(1)
      end do

    end subroutine eval_kernel_forward
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine propagate_forward(gde)
!-----------------------------------------------------------------------
      implicit none

      type(s_gde), intent(inout) :: gde

      integer              :: i, j
      integer              :: nstep, ncut, nextend
      real(8)              :: dt
      real(8)              :: diff, cum, ri

      ! Setup variables
      !
      nstep   = gde%nstep
      dt      = gde%dt
      ncut    = gde%ncut
      nextend = gde%nextend 

      ! Extend dP(t)/dt & P(t)
      !

      !   Initialize
      !
      gde%p_ext       = 0.0d0
      gde%pdot_ext    = 0.0d0

      gde%p_ext(0)    = gde%p(0)
      gde%p_ext(1)    = gde%p(1)
      gde%p_ext(2)    = gde%p(2)

      gde%pdot_ext(0) = gde%pdot(0) 
      gde%pdot_ext(1) = gde%pdot(1)

      !   Propagate
      !
      do i = 2, nextend - 1 

        cum = - gde%p_ext(i)

        ri  = 0.0d0
        do j = 1, min(i - 1, ncut)
          ri = ri - dt * gde%kernel(j) * gde%pdot_ext(i - j)
        end do

        cum              = cum + ri
        gde%pdot_ext(i)  = cum / (gde%kernel(0) * dt)
        gde%p_ext(i + 1) = gde%p_ext(i) + gde%pdot_ext(i) * dt

      end do

    end subroutine propagate_forward
!-----------------------------------------------------------------------

end module mod_gderp_forward
!=======================================================================
