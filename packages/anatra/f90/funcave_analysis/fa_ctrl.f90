!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_cv
  use mod_bootstrap
  implicit none

  ! constants
  !

  ! structures
  !
  type :: s_option
    logical :: use_bootstrap = .false.
    real(8) :: xsta          = 0.0d0
    real(8) :: dx            = 0.1d0
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, cvinfo, bootopt)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option 
      type(s_cvinfo),  intent(out) :: cvinfo
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

        if (option%use_bootstrap) then
          call read_ctrl_bootstrap(iunit, bootopt)
        end if

      close(iunit)

      open(iunit, file=trim(input%flist))
        call read_ctrl_cvinfo  (iunit, cvinfo)
      close(iunit) 

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer, parameter :: ndim_max = 3 
!
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      logical :: use_bootstrap = .false.
      real(8) :: xsta          = 0.0d0
      real(8) :: dx            = 0.1d0

      integer :: i, j
      integer :: iopt, ierr

      namelist /option_param/ &
              use_bootstrap,  &
              xsta,           &
              dx

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("use_bootstrap = ", a)') get_tof(use_bootstrap)
      write(iw,'("xsta          = ", f20.10)') xsta 
      write(iw,'("dx            = ", f20.10)') dx

      !iopt = get_opt(mode, CoMMode, ierr)
      !if (ierr /= 0) then
      !  write(iw,'("Read_Ctrl_Option> Error.")')
      !  write(iw,'("mode = ",a," is not available.")') trim(mode)
      !  stop
      !end if
      !option%mode          = iopt

      option%use_bootstrap = use_bootstrap
      option%xsta          = xsta
      option%dx            = dx

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

end module
!=======================================================================
