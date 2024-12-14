!=======================================================================
program main 
!=======================================================================
  use mod_const
  use mod_ctrl
  !use mod_setup
  use mod_dcdio
  use mod_xtcio
  use xdr, only: xtcfile
  use mod_analyze

  type(s_input)   :: input
  type(s_output)  :: output
  type(s_option)  :: option 
  type(s_dcd)     :: dcdin
  type(xtcfile)   :: xtcin

  call show_title
  call show_usage
  call read_ctrl(input, output, option)
  !call setup(input, dcdin, xtcin)
  call analyze(input, output, option)
  !call dealloc_dcd(dcdin)
  call termination('ed_analysis')


end program main 
!=======================================================================

!-----------------------------------------------------------------------
subroutine show_title
!-----------------------------------------------------------------------
  implicit none

  write(6,'("==================================================")')
  write(6,*)
  write(6,'("            Extract DCD Analysis")')
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
    write(iw,'(" ftraj = ""inp.dcd"" ! input trajectory file (dcd or xtc)")')
    write(iw,'(" fts   = ""ts.dat""   ! time-series data")')
    write(iw,'("/")')
    write(iw,*)
    write(iw,'("&output_param")')
    write(iw,'(" ftraj = ""inp.dcd"" ! output dcd file")')
    write(iw,'(" fts   = ""ts.dat""  ! output time-series data")')
    write(iw,'("/")')
    write(iw,*)
    write(iw,'("&option_param")')
    write(iw,'(" dt   = 0.1d0 ")')
    write(iw,'(" ndim = 2  ! dimension of time-series data")')
    write(iw,'(" react_range = -5.0d0 5.0d0 1.0d0 0.0d0")')
    write(iw,'("/")')
    stop
  end if


end subroutine show_usage
!-----------------------------------------------------------------------

