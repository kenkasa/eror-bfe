!=======================================================================
module mod_gderp_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_movave_ctrl
  implicit none

  ! constants
  !
  integer,      parameter, public :: IntegratorTypeFORWARD  = 1
  integer,      parameter, public :: IntegratorTypeBACKWARD = 2 
  character(*), parameter, public :: IntegratorTypes(2) &
                                       = (/'FORWARD   ',&
                                           'BACKWARD  '/)

  ! structures
  !
  type :: s_option
    real(8) :: dt         = 0.001d0
    real(8) :: t_cut      = 10.0d0
    real(8) :: t_extend   = 150.0d0
    integer :: integrator = IntegratorTypeFORWARD 
  end type s_option

  ! subroutines
  !
  public :: read_ctrl
  public :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, movave_option)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter                   :: iunit = 10

      type(s_input),           intent(out) :: input
      type(s_output),          intent(out) :: output
      type(s_option),          intent(out) :: option
      type(s_movave_option),   intent(out) :: movave_option

      character(len=MaxChar)       :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input          (iunit, input)
        call read_ctrl_output         (iunit, output)
        call read_ctrl_option         (iunit, option)
        call read_movave_ctrl_option  (iunit, movave_option, .true.)
      close(iunit)

      ! Several variables are transferred to movave structure
      movave_option%dt    = option%dt
      movave_option%t_sta = 0.0d0

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,         intent(in)  :: iunit
      type(s_option),  intent(out) :: option 

      real(8)                :: dt           = 0.1d0
      real(8)                :: t_cut        = 0.0d0 
      real(8)                :: t_extend     = 0.0d0
      integer                :: npoint(2)    = 5
      character(len=MaxChar) :: integrator   = "forward"

      integer                :: iopt, ierr

      namelist /option_param/ dt,          &
                              t_cut,       &
                              t_extend,    &
                              integrator

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("dt           = ", f15.7)')     dt
      write(iw,'("t_cut        = ", f15.7)')     t_cut 
      write(iw,'("t_extend     = ", f15.7)')     t_extend
      write(iw,'("integrator   = ", a)')         trim(integrator)

      option%dt           = dt
      option%t_cut        = t_cut
      option%t_extend     = t_extend

      iopt = get_opt(integrator, IntegratorTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("integrator = ",a," is not available.")') &
          trim(integrator)
        stop
      end if
      option%integrator   = iopt

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
