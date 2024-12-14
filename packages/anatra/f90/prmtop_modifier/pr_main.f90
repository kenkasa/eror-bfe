!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option


  call show_title
  call read_ctrl(input, output, option)
  call analyze(input, output, option)
  call termination("Prmtop modifier")

end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("                Prmtop Modifier")')
  write(6,*)
  write(6,'("==================================================")') 

end subroutine show_title
!-----------------------------------------------------------------------
