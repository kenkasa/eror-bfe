!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_setup
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option 


  call show_title
  call show_usage
  call read_ctrl(input, output, option)

  if (option%shuffle) then
    call analyze_shuffle(input, output, option)
  else
    call analyze(input, output, option)
  end if
  call termination("ss_analysis")


end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("           Structural Sampling Analysis")')
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
    write(iw,'(" fdcd = ""inp.dcd"" ! input dcd file")')
    write(iw,'("/")')
    write(iw,*)
    write(iw,'("&output_param")')
    write(iw,'(" fdcd = ""inp.dcd"" ! output dcd file")')
    write(iw,'("/")')
    write(iw,*)
    write(iw,'("&option_param")')
    write(iw,'(" nsample = 100")')
    write(iw,'(" iseed   = 3141592")')
    write(iw,'("/")')
    stop
  end if


end subroutine show_usage
!-----------------------------------------------------------------------

