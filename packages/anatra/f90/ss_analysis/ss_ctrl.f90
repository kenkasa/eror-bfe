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

  ! structures
  !
  type :: s_option
    integer              :: nsample     =  1
    integer              :: iseed       = -1
    logical              :: duplicate   = .false.
    logical              :: out_rst7    = .false.
    logical              :: use_allsnap = .false.
    logical              :: shuffle     = .false.
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

      integer :: nsample     = 100
      integer :: iseed       = -1
      logical :: duplicate   = .false.
      logical :: out_rst7    = .false.
      logical :: use_allsnap = .false.
      logical :: shuffle     = .false.

      integer :: i, j
      integer :: iopt, ierr

      namelist /option_param/ nsample,      &
                              iseed,        &
                              duplicate,    &
                              out_rst7,     &
                              use_allsnap,  &
                              shuffle

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("nsample       = ", i0)') nsample 
      write(iw,'("iseed         = ", i0)') iseed 
      write(iw,'("duplicate     = ", a)')  get_tof(duplicate) 
      write(iw,'("out_rst7      = ", a)')  get_tof(out_rst7) 
      write(iw,'("use_allsnap   = ", a)')  get_tof(use_allsnap) 
      write(iw,'("shuffle       = ", a)')  get_tof(shuffle) 
      write(iw,*)

      !iopt = get_opt(mode, CoMMode, ierr)
      !if (ierr /= 0) then
      !  write(iw,'("Read_Ctrl_Option> Error.")')
      !  write(iw,'("mode = ",a," is not available.")') trim(mode)
      !  stop
      !end if
      !option%mode          = iopt
                           
      option%nsample       = nsample 
      option%iseed         = iseed
      option%duplicate     = duplicate
      option%use_allsnap   = use_allsnap
      option%out_rst7      = out_rst7
      option%shuffle       = shuffle

      ! Combination check
      !

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

end module
!=======================================================================
