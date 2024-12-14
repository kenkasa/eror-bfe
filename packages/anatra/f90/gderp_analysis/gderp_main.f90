!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_gderp_ctrl
  use mod_gderp_analyze
  use mod_movave_ctrl
  use mod_movave_analyze

  type(s_input)         :: input
  type(s_output)        :: output
  type(s_option)        :: option
  type(s_movave_option) :: movave_option
  !type(s_trajopt) :: trajopt
  !type(s_traj)    :: traj

  call show_title
  call show_usage
  call read_ctrl(input, output, option, movave_option)
  call gderp_analyze(input, output, option, movave_option)

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("    Generalized Diffusion Equation Analysis")')
  write(6,'("         for Returning Probability")')
  write(6,*)
  write(6,'("==================================================")') 

end subroutine show_title
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine show_usage
!-----------------------------------------------------------------------
  use mod_const

  implicit none

  character(len=MaxChar) :: f_ctrl


  call getarg(1, f_ctrl)

  if (trim(f_ctrl) == "-h") then
    write(iw,'("&input_param")')
    write(iw,'(" fts = ""data""   ! Input file name")')
    write(iw,'("/")')

    write(iw,'("&output_param")')
    write(iw,'(" fhead = ""out""    ! Output file header")')
    write(iw,'("/")')

    write(iw,'("&option_param")')
    write(iw,'("  dt           = 0.001d0  ! time grid")')
    write(iw,'("  t_cut        = 10.0d0   ! cutoff time for kernel")')
    write(iw,'("  t_extend     = 100.0d0  ! extended time region for P(t)")')
    write(iw,'("/")')

    write(iw,'("&movave_option_param")')
    write(iw,'("  t_sep        = 0.0d0    ! boundary of two different region")')
    write(iw,'("  nregion      = 2        ! # of diffrent regions")')
    write(iw,'("  npoint       = 3 5      ! # of data points used for each region (odd number)")')
    write(iw,'("  include_zero = .false.  ! whether origin is changed or not by the average")')
    write(iw,'("/")')
    write(iw,*)

    stop
  end if

end subroutine show_usage
!-----------------------------------------------------------------------
