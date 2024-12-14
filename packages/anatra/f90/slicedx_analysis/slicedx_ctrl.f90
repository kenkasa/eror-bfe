!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_grid3d
  use mod_input
  use mod_output
  use mod_cv
  use mod_bootstrap
  implicit none

  ! constants
  !
  integer, parameter, public :: MaxDim   = 3
  integer, parameter, public :: MaxState = 100 
  integer, parameter, public :: MaxNcell = 100 

  ! structures
  !

  type :: s_option
    real(8) :: zval = -1.0d+100
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, cvinfo, mpl2d, bootopt)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_cvinfo),  intent(out) :: cvinfo
      type(s_mpl2d),   intent(out) :: mpl2d
      type(s_bootopt), intent(out) :: bootopt

      character(len=MaxChar)       :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input  (iunit, input)
        call read_ctrl_output (iunit, output)
        call read_ctrl_option (iunit, option)
      close(iunit)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      real(8) :: zval = 1.0d+100 

      integer :: i, j, k 

      namelist /option_param/ &
        zval


      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("zval              = ", f15.7)')  zval 

      ! setup option variables
      !
      option%zval = zval 

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
