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

  ! structures
  !
  type :: s_option
    real(8) :: dr            = 0.1d0
    integer :: mode(2)       = (/CoMModeRESIDUE, CoMModeRESIDUE/)
    logical :: identical     = .false.
    logical :: normalize     = .false.
    logical :: separate_self = .false. 
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

      real(8)                :: dr            = 0.1d0
      character(len=MaxChar) :: mode(2)       = (/"RESIDUE", "RESIDUE"/)
      logical                :: identical     = .false.
      logical                :: normalize     = .false.
      logical                :: separate_self = .false.

      integer                :: itraj
      integer                :: iopt, ierr

      namelist /option_param/ dr, mode, identical, normalize, &
                              separate_self

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("dr            = ", f15.7)')   dr 
      write(iw,'("mode          = ", a,2x,a)')  trim(mode(1)), trim(mode(2))
      write(iw,'("identical     = ", a)')       get_tof(identical)
      write(iw,'("normalize     = ", a)')       get_tof(normalize)
      write(iw,'("separate_self = ", a)')       get_tof(separate_self)

      do itraj = 1, 2
        iopt = get_opt(mode(itraj), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(itraj))
          stop
        end if
        option%mode(itraj) = iopt
      end do
                           
      option%dr            = dr
      option%identical     = identical
      option%normalize     = normalize
      option%separate_self = separate_self

      ! Combination check
      !

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
