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

  integer,      parameter, public :: CalcTypeCOORD     = 1
  integer,      parameter, public :: CalcTypeVEL       = 2 
  character(*), parameter, public :: CalcTypes(2)      = (/'COORD     ',&
                                                           'VEL       '/)


  ! structures
  !
  type :: s_option
    integer :: mode      = CoMModeRESIDUE
    real(8) :: dt        = 1.0d0
    integer :: calctype  = CalcTypeCOORD
    integer :: tcfrange  = 0 
    logical :: bigtraj   = .false.
    real(8) :: t_sparse  = -1.0d0
    real(8) :: t_range   = -1.0d0

    ! prepared after reading namelist
    !
    real(8) :: dt_out    = 0.0d0
    integer :: nstep     = 0
    integer :: nt_sparse = 0
    integer :: nt_range  = 0

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
      type(s_input_check)          :: input_check

      ! Get control file name
      !
      call getarg(1, f_ctrl)

      ! Setup input checker
      !
      call setup_input_check(input_check)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input  (iunit, input)
        call read_ctrl_output (iunit, output)
        call read_ctrl_option (iunit, option)
        call read_ctrl_trajopt(iunit, trajopt)

        call check_input(input, option, input_check) 
      close(iunit)

    end subroutine read_ctrl
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine setup_input_check(ic)
!-----------------------------------------------------------------------
      implicit none

      type(s_input_check), intent(out) :: ic


      call initialize_input_check(ic)

      ic%require(InputFTRAJ) = .true.

    end subroutine setup_input_check
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      logical                :: bigtraj   = .false.
      integer                :: tcfrange  = 1000
      character(len=MaxChar) :: calctype  = 'COORD'
      real(8)                :: dt        = 1.0d0
      real(8)                :: t_sparse  = -1.0d0
      real(8)                :: t_range   = -1.0d0

      integer                :: iopt, ierr

      namelist /option_param/ &
        dt,       &
        calctype, &
        bigtraj,  &
        t_sparse, &
        t_range,  &
        tcfrange

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("calctype  = ", a)')     trim(calctype)
      write(iw,'("bigtraj   = ", a)')     get_tof(bigtraj)
      write(iw,'("tcfrange  = ", i0)')    tcfrange
      write(iw,'("dt        = ", f15.7)') dt
      write(iw,'("t_sparse  = ", f15.7)') t_sparse 
      write(iw,'("t_range   = ", f15.7)') t_range

      ! Get calctype 
      !
      iopt = get_opt(calctype, CalcTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("calctype = ",a," is not available.")') trim(calctype)
        stop
      end if
      option%calctype  = iopt

      option%mode      = CoMModeATOM 

      ! Namelist => Option
      !
      option%tcfrange  = tcfrange 
      option%dt        = dt
      option%bigtraj   = bigtraj
      option%t_sparse  = t_sparse 
      option%t_range   = t_range

      ! Check combinations
      !

      ! Convert from real to integer
      !
      if (option%t_sparse < 0.0d0) then
        option%t_sparse = option%dt
      end if

      option%nt_sparse = nint(option%t_sparse / option%dt)
      option%dt_out    = option%dt * option%nt_sparse

      if (option%t_range > 0.0d0) then
        option%nt_range = nint(option%t_range / option%dt)
      else
        option%nt_range = option%tcfrange
      end if

      if (option%nt_sparse > 1) then
        option%nt_range = nint(dble(option%nt_range) / dble(option%nt_sparse))  
      end if

      write(iw,'("dt_out    = ", f15.7)') option%dt_out
      write(iw,'("nt_sparse = ", i0)')    option%nt_sparse



    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine check_input(input, option, ic)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),       intent(in) :: input
      type(s_option),      intent(in) :: option
      type(s_input_check), intent(in) :: ic

      integer :: i


      do i = 1, NInputParm

        if (ic%require(InputFTRAJ)) then
          if (.not. input%specify_ftraj) then
            write(iw,*)
            write(iw,'("Check_Input> Error.")')
            write(iw,'("ftraj or flist_traj : not specified in input_param")')
            stop
          end if
        end if

      end do

    end subroutine check_input
!-----------------------------------------------------------------------

!!-----------------------------------------------------------------------
!    subroutine check_combination(input, option)
!!-----------------------------------------------------------------------
!      implicit none
!
!      type(s_input)  :: input
!      type(s_option) :: option
!
!      
!      if (option%use_cond) then
!        if (trim(input%fts) == "") then
!          write(iw,'("Check_Combination> Error.")')
!          write(iw,'("fcv should be specified in input_param")')
!          write(iw,'("if use_cond = .true. in option_param.")')
!          stop
!        end if
!      end if
!
!    end subroutine check_combination
!!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
