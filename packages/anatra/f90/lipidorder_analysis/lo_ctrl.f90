!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_traj
  implicit none

  ! constants
  !
  !integer,      parameter, public :: ParameterAAA = 1
  !integer,      parameter, public :: ParameterBBB = 2 
  !character(*), parameter, public :: Parameter(2) = (/'AAAAA     ',&
  !                                                    'AAAAA     '/)

  ! structures
  !
  ! not used in this script
  type :: s_option
    integer :: carbon_id(MaxTraj)  
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

      integer :: carbon_id(MaxTraj)

      integer :: i
      !integer                :: iopt, ierr

      namelist /option_param/ carbon_id 


      ! Initialize
      !
      do i = 1, MaxTraj 
        carbon_id(i) = i
      end do

      ! Read
      !
      rewind iunit
      read(iunit, option_param)

      ! Output 
      !
      !write(iw,*)
      !write(iw,'(">> Option section parameters")')
      !write(iw,'("carbon_id   = ", i0)') carbon_id 

      option%carbon_id = carbon_id 

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
