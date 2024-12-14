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
    real(8) :: dx       = 0.1d0
    integer :: xsta     = 0.0d0 
    integer :: nx       = 100
    integer :: ncol     = 1
    integer :: data_sta = 0
    integer :: data_end = 0
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
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      real(8)                :: dx       = 0.1d0
      real(8)                :: xsta     = 0.0d0 
      integer                :: nx       = 100
      integer                :: ncol     = 1
      integer                :: data_sta = 0
      integer                :: data_end = 0


      namelist /option_param/ dx, xsta, nx, ncol, &
                              data_sta, data_end

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("dx       = ", f15.7)')   dx
      write(iw,'("xsta     = ", f15.7)')   xsta 
      write(iw,'("nx       = ", i0)')      nx
      write(iw,'("ncol     = ", i0)')      ncol
      write(iw,'("data_sta = ", i0)')      data_sta 
      write(iw,'("data_end = ", i0)')      data_end

      option%dx       = dx
      option%xsta     = xsta
      option%nx       = nx
      option%ncol     = ncol
      option%data_sta = data_sta
      option%data_end = data_end

    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
