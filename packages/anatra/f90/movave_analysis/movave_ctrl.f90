!=======================================================================
module mod_movave_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  implicit none

  ! constants
  !
  integer, parameter :: MaxSep = 10

  ! structures
  !
  type :: s_movave_option
    real(8) :: dt                    = 0.1d0
    real(8) :: t_sta                 = 0.0d0
    real(8) :: t_sep(MaxSep - 1)     = 0.0d0
    integer :: nregion               = 1
    integer :: npoint(MaxSep)        = 5
    logical :: include_zero = .false.
  end type s_movave_option

  ! subroutines
  !
  public :: read_movave_ctrl
  public :: read_movave_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_movave_ctrl(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter                   :: iunit = 10

      type(s_input),           intent(out) :: input
      type(s_output),          intent(out) :: output
      type(s_movave_option),   intent(out) :: option

      character(len=MaxChar)       :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Movave_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input          (iunit, input)
        call read_ctrl_output         (iunit, output)
        call read_movave_ctrl_option  (iunit, option)
      close(iunit)

    end subroutine read_movave_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_movave_ctrl_option(iunit, option, external_use)
!-----------------------------------------------------------------------
      implicit none
!
      integer,                intent(in)  :: iunit
      type(s_movave_option),  intent(out) :: option 
      logical, optional,      intent(in)  :: external_use

      real(8)                :: dt               = 0.1d0
      real(8)                :: t_sta            = 0.0d0 
      real(8)                :: t_sep(MaxSep -1) = 0.0d0
      integer                :: nregion          = 1
      integer                :: npoint(MaxSep)   = 5
      logical                :: include_zero     = .true.
      logical                :: extr


      namelist /movave_option_param/ dt,          &
                                     t_sta,       &
                                     t_sep,       &
                                     nregion,     &
                                     npoint,      &
                                     include_zero

      rewind iunit
      read(iunit, movave_option_param)

      extr = .false.
      if (present(external_use)) & 
        extr = external_use

      write(iw,*)
      write(iw,'(">> Option section parameters")')

      if (.not. extr) then 
        write(iw,'("dt           = ", f15.7)')     dt
        write(iw,'("t_sta        = ", f15.7)')     t_sta
      end if

      write(iw,'("t_sep        = ", 10(f15.7))') t_sep(1:2)
      write(iw,'("nregion      = ", i0)')        nregion
      write(iw,'("npoint       = ", 10(i0,2x))') npoint(1:3)
      write(iw,'("include_zero = ", a)')         get_tof(include_zero)

      option%dt           = dt
      option%t_sta        = t_sta
      option%t_sep        = t_sep
      option%nregion      = nregion
      option%npoint       = npoint
      option%include_zero = include_zero 

    end subroutine read_movave_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
