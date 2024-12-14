!=======================================================================
module mod_gderp_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_movave_ctrl
  use mod_movave_analyze
  use mod_gderp_ctrl
  use mod_gderp_str
  use mod_gderp_forward
  use mod_gderp_backward

  ! structures
  !

  ! subroutines
  !
  public  :: gderp_analyze
  private :: eval_timeconst
  private :: gderp_write

  contains
!-----------------------------------------------------------------------
    subroutine gderp_analyze(input, output, option, movave_option)
!-----------------------------------------------------------------------
       implicit none

       type(s_input),                   intent(in)  :: input
       type(s_output),                  intent(in)  :: output
       type(s_option),                  intent(in)  :: option
       type(s_movave_option),           intent(in)  :: movave_option

       type(s_movave) :: movave
       type(s_gde)    :: gde

       integer        :: nstep, ncut, nextend
       real(8)        :: dt, t_cut


       ! Read input data & Perform moving average 
       !
       call movave_analyze(input, output, movave_option, movave)

       ! Write averaged data for user check
       !
       call movave_write(output, movave)

       ! Setup
       !
       nstep   = movave%ngrid
       dt      = option%dt
       ncut    = nint(option%t_cut    / dt)
       nextend = nint(option%t_extend / dt)

       if (ncut >= nstep - 3) &
         ncut = nstep - 3 

       allocate(gde%p       (0:nstep-1),   &
                gde%pdot    (0:nstep-1),   &
                gde%kernel  (0:nextend),   &
                gde%p_ext   (0:nextend),   &
                gde%pdot_ext(0:nextend),   &
                gde%tau_p   (0:nextend),   &
                gde%tau_k   (0:nextend))

       gde%dt       = dt
       gde%nstep    = nstep
       gde%ncut     = ncut
       gde%nextend  = nextend
       gde%p        = movave%data
       gde%pdot     = movave%deriv
       gde%kernel   = 0.0d0
       gde%p_ext    = 0.0d0
       gde%pdot_ext = 0.0d0
       gde%tau_p    = 0.0d0
       gde%tau_k    = 0.0d0

       if (option%integrator == IntegratorTypeFORWARD) then
         call eval_kernel_forward(gde)
         call propagate_forward  (gde)
       else
         call eval_kernel_backward(gde)
         call propagate_backward  (gde)
         !write(iw,'("Gderp_Analyze> Error.")')
         !write(iw,'("Sorry, only integrator = forward is supported currently.")')
         !stop
       end if

       call eval_timeconst(option%integrator, gde)
       call gderp_write(output, gde)

    end subroutine gderp_analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine eval_timeconst(integrator, gde)
!-----------------------------------------------------------------------
      implicit none

      integer,     intent(in)    :: integrator
      type(s_gde), intent(inout) :: gde

      integer              :: i, j, ishift
      integer              :: nstep, ncut, nextend
      real(8)              :: dt
      real(8)              :: diff, cum, ri


      ! Setup variables
      !
      nstep   = gde%nstep
      dt      = gde%dt
      ncut    = gde%ncut
      nextend = gde%nextend 

      ! Calculate time constant for kernel
      !
      if (integrator == IntegratorTypeFORWARD) then
        ishift = 0
      else if (integrator == IntegratorTypeBACKWARD) then
        ishift = 1
      end if

      gde%tau_k = 0.0d0
      ri        = 0.0d0
      do i = 0, ncut
        ri           = ri + gde%kernel(i + ishift) * dt 
        gde%tau_k(i) = ri 
      end do
      gde%tau_k(ncut+1:nextend) = ri

      ! Calculate time constant for P(t)
      !
      gde%tau_p = 0.0d0
      ri        = 0.0d0
      do i = 0, nextend
        ri           = ri + gde%p_ext(i) * dt
        gde%tau_p(i) = ri
      end do

    end subroutine eval_timeconst
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine gderp_write(output, gde)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in) :: output
      type(s_gde),    intent(in) :: gde 

      integer                :: i, ncmp
      real(8)                :: dt, time
      character(len=MaxChar) :: fwrite


      dt   = gde%dt
      ncmp = 5 
      write(fwrite,'(a,".gderp")') trim(output%fhead)
      open(10, file = trim(fwrite))
        write(10,'("# 1: time")')
        write(10,'("# 2: P(t)")')
        write(10,'("# 3: tau_P(t)")')
        write(10,'("# 4: Kernel(t)")')
        write(10,'("# 5: tau_K(t)")')
        do i = 0, gde%nextend
          time = dt * i
          write(10, '(e15.7,2x)', advance = 'no') time
          write(10, '(e15.7,2x)', advance = 'no') gde%p_ext (i)
          write(10, '(e15.7,2x)', advance = 'no') gde%tau_p (i)
          write(10, '(e15.7,2x)', advance = 'no') gde%kernel(i) / gde%kernel(0)
          write(10, '(e15.7,2x)', advance = 'no') gde%tau_k (i)
          write(10,*)
        end do
      close(10)

    end subroutine gderp_write
!-----------------------------------------------------------------------

end module mod_gderp_analyze
!=======================================================================
