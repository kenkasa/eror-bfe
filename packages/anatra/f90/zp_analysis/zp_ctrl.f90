!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_traj
  use mod_com

  implicit none

  ! constants
  !
  integer,      parameter, public :: DensityTypeNUMBER   = 1
  integer,      parameter, public :: DensityTypeELECTRON = 2 
  character(*), parameter, public :: DensityType(2)      = (/'NUMBER   ',&
                                                             'ELECTRON '/)

  integer,      parameter, public :: CenterTypeZERO      = 1
  integer,      parameter, public :: CenterTypeHALF      = 2 
  character(*), parameter, public :: CenterTypes(2)      = (/'ZERO',&
                                                             'HALF'/)

  ! structures
  !
  type :: s_option
    real(8) :: dz         = 0.1d0
    integer :: mode       = CoMModeRESIDUE
    integer :: denstype   = DensityTypeNUMBER
    integer :: centertype = CenterTypeZERO 
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, trajopt)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_trajopt), intent(out) :: trajopt

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
        call read_ctrl_trajopt(iunit, trajopt)
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

      real(8)                :: dz         = 0.1d0
      character(len=MaxChar) :: mode       = "RESIDUE"
      character(len=MaxChar) :: denstype   = "NUMBER"
      character(len=MaxChar) :: centertype = "ZERO" 

      integer                :: iopt, ierr

      namelist /option_param/ dz, mode, denstype, centertype

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("dz         = ", f15.7)') dz 
      write(iw,'("mode       = ", a)')     trim(mode)
      write(iw,'("denstype   = ", a)')     trim(denstype)
      write(iw,'("centertype = ", a)')     trim(centertype)


      iopt = get_opt(mode, CoMMode, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("mode = ",a," is not available.")') trim(mode)
        stop
      end if
      option%mode     = iopt

      iopt = get_opt(denstype, DensityType, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("mode = ",a," is not available.")') trim(mode)
        stop
      end if
      option%denstype = iopt

      iopt = get_opt(centertype, CenterTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("mode = ",a," is not available.")') trim(mode)
        stop
      end if
      option%centertype = iopt

      option%dz         = dz

      ! Combination check
      !
      if (option%mode /= ComModeATOM .and. &
          option%denstype == DensityTypeELECTRON ) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("Combination mode = ATOM and denstype = ELECTIONS")')
        write(iw,'("is not available.")')
      end if 

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
