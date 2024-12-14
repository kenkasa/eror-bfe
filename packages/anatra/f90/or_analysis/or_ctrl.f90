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
  integer,      parameter, public :: XcoordModeCOST  = 1
  integer,      parameter, public :: XcoordModeTHETA = 2 
  character(*), parameter, public :: XcoordMode(2) = (/'COST      ',&
                                                       'THETA     '/)

  integer,      parameter, public :: JudgeUpModeNMOLUP = 1
  integer,      parameter, public :: JudgeUpModeCOORD  = 2 
  integer,      parameter, public :: JudgeUpModeNONE   = 3 
  character(*), parameter, public :: JudgeUpMode(3) = (/'NMOLUP   ',&
                                                        'COORD    ',&
                                                        'NONE     '/)

  ! structures
  !
  type :: s_option
    integer :: judgeup  = JudgeUpModeNMOLUP 
    integer :: mode(2)  = (/CoMModeRESIDUE, CoMModeRESIDUE/)
    integer :: xcoord   = XcoordModeCOST
    real(8) :: dcost    = 0.1d0
    real(8) :: dtheta   = 1.0d0
    integer :: nmolup   = 0
    integer :: nmollow  = 0
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

      character(len=MaxChar) :: judgeup  = "NMOLUP"
      character(len=MaxChar) :: mode(2)  = (/"RESIDUE", "RESIDUE"/) 
      character(len=MaxChar) :: xcoord   = "COST"
      real(8)                :: dcost    = 0.1d0
      real(8)                :: dtheta   = 1.0d0
      integer                :: nmolup   = 0
      integer                :: nmollow  = 0

      integer                :: itraj
      integer                :: iopt, ierr

      namelist /option_param/ judgeup, mode, xcoord, dcost, dtheta, &
                              nmolup, nmollow


      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("judgeup  = ", a)')      trim(judgeup)
      write(iw,'("nmolup   = ", i0)')     nmolup 
      write(iw,'("mode     = ", a,2x,a)') trim(mode(1)), trim(mode(2))
      write(iw,'("xcoord   = ", a)')      trim(xcoord)
      write(iw,'("dcost    = ", f15.7)')  dcost
      write(iw,'("dtheta   = ", f15.7)')  dtheta

      iopt = get_opt(judgeup, JudgeUpMode, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("judgeup = ",a," is not available.")') trim(judgeup)
        stop
      end if
      option%judgeup = iopt

      do itraj = 1, 2
        iopt = get_opt(mode(itraj), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(itraj))
          stop
        end if
        option%mode(itraj) = iopt
      end do

      iopt = get_opt(xcoord, XcoordMode, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("xcoord = ",a," is not available.")') trim(xcoord)
        stop
      end if
      option%xcoord   = iopt

      option%dcost    = dcost
      option%dtheta   = dtheta
      option%nmolup   = nmolup
      option%nmollow  = nmollow


    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
