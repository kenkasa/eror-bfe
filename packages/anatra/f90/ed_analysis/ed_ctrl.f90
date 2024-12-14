!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  implicit none

  ! constants
  !
  integer, parameter, public :: Nstate   = 1 
  integer, parameter, public :: REACTIVE = 1 
  integer, parameter, public :: OTHERS   = 2 
  integer, parameter, public :: StateInfo(Nstate) = (/1/)

  ! structures
  !
  type :: s_option
    integer              :: ndim         = 1
    real(8)              :: dt           = 1.0d0
    real(8), allocatable :: state_def(:, :, :)
    real(8), allocatable :: react_range(:, :)
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option

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
      integer, parameter :: ndim_max = 3 
!
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      integer :: ndim                       = 1
      real(8) :: dt                         = 1.0d0
      real(8) :: react_range(2, ndim_max)   = 0.0d0

      integer :: i, j
      integer :: iopt, ierr

      namelist /option_param/ ndim, dt, react_range

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("ndim          = ", i0)') ndim
      write(iw,*)
      write(iw,'("react_range   =")')
      do i = 1, ndim
        write(iw,'(es15.7, " <= component ",i0," < ",es15.7)') &
          react_range(1, i), i, react_range(2, i) 
      end do

      !iopt = get_opt(mode, CoMMode, ierr)
      !if (ierr /= 0) then
      !  write(iw,'("Read_Ctrl_Option> Error.")')
      !  write(iw,'("mode = ",a," is not available.")') trim(mode)
      !  stop
      !end if
      !option%mode          = iopt
                           
      option%ndim          = ndim
      option%dt            = dt

      allocate(option%react_range(1:2, ndim))
      allocate(option%state_def(2, ndim, Nstate))

      do i = 1, ndim
        option%react_range(1:2, i)   = react_range(1:2, i)
      end do

      do i = 1, ndim
        option%state_def(1:2, i, 1) = react_range(1:2, i)
      end do

      ! Combination check
      !

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

end module
!=======================================================================
