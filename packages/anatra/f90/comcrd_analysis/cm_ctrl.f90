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

  integer,      parameter, public :: CoMFormatTypeXYZ        = 1
  integer,      parameter, public :: CoMFormatTypeTIMESERIES = 2 
  character(*), parameter, public :: CoMFormatType(2)        = (/'XYZ       ',&
                                                             'TIMESERIES'/)


  ! structures
  !
  type :: s_option
    integer :: mode      = CoMModeRESIDUE
    integer :: comformat = CoMFormatTypeXYZ 
    integer :: msddim    = 3
    integer :: msdrange  = 1000
    real(8) :: dt        = 1.0d0
    real(8) :: dr        = 0.1d0     ! for vhfcalc
    real(8) :: rmax      = 50.0d0    ! for vhfcalc
    real(8) :: t_sparse  = -1.0d0    ! for vhfcalc
    real(8) :: t_range   = -1.0d0    ! for vhfcalc
    real(8) :: rcrange(2)= 0.0d0
    logical :: msdcalc   = .false.
    logical :: ngpcalc   = .false.
    logical :: vhfcalc   = .false.
    logical :: comcalc   = .false. 
    logical :: onlyz     = .false.
    logical :: unwrap    = .false.
    logical :: use_cond  = .false.

    ! prepared after reading namelists
    !
    real(8) :: dt_out    = 0.0d0
    integer :: nstep     = 0
    integer :: nt_sparse = 0
    integer :: nt_range  = 0
    integer :: nr        = 0
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

      character(len=MaxChar) :: mode      = "RESIDUE"
      character(len=MaxChar) :: comformat = "XYZ"
      logical                :: unwrap    = .false. 
      logical                :: msdcalc   = .false.
      logical                :: ngpcalc   = .false.
      logical                :: vhfcalc   = .false.
      logical                :: comcalc   = .false.
      logical                :: onlyz     = .false.
      logical                :: use_cond  = .false.
      integer                :: msddim    = 3 
      integer                :: msdrange  = 1000
      real(8)                :: dt        = 1.0d0
      real(8)                :: dr        = 0.1d0
      real(8)                :: rmax      = 50.0d0
      real(8)                :: t_sparse  = -1.0d0
      real(8)                :: t_range   = -1.0d0
      real(8)                :: rcrange(2)= 0.0d0

      integer                :: iopt, ierr

      namelist /option_param/ mode,        &
                              comformat,   &
                              unwrap,      &
                              msdcalc,     &
                              ngpcalc,     &
                              vhfcalc,     &
                              comcalc,     &
                              onlyz,       &
                              use_cond,    &
                              msddim,      &
                              msdrange,    &
                              dt,          &
                              dr,          &
                              rmax,        &
                              t_sparse,    &
                              t_range,     &
                              rcrange 

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("mode      = ", a)')     trim(mode)
      write(iw,'("comformat = ", a)')     trim(comformat)
      write(iw,'("unwrap    = ", a)')     get_tof(unwrap)
      write(iw,'("msdcalc   = ", a)')     get_tof(msdcalc)
      write(iw,'("ngpcalc   = ", a)')     get_tof(ngpcalc)
      write(iw,'("vhfcalc   = ", a)')     get_tof(vhfcalc)
      write(iw,'("comcalc   = ", a)')     get_tof(comcalc)
      write(iw,'("onlyz     = ", a)')     get_tof(onlyz)
      write(iw,'("use_cond  = ", a)')     get_tof(use_cond)
      write(iw,'("msddim    = ", i0)')    msddim 
      write(iw,'("msdrange  = ", i0)')    msdrange
      write(iw,'("dt        = ", f15.7)') dt
      write(iw,'("dr        = ", f15.7)') dr
      write(iw,'("rmax      = ", f15.7)') rmax
      write(iw,'("t_sparse  = ", f15.7, " : Used for VHF calculation")') t_sparse 
      write(iw,'("t_range   = ", f15.7, " : Used for VHF calculation")') t_range
      write(iw,'("rcrange   = ", e15.7, 2x, e15.7)') rcrange(1), rcrange(2)


      ! Get mode
      !
      iopt = get_opt(mode, CoMMode, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("mode = ",a," is not available.")') trim(mode)
        stop
      end if
      option%mode      = iopt

      ! Get comformat
      !
      iopt = get_opt(comformat, CoMFormatType, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("comformat = ",a," is not available.")') trim(comformat)
        stop
      end if
      option%comformat = iopt

      ! Namelist => Option
      !
      option%msddim    = msddim
      option%msdrange  = msdrange
      option%comcalc   = comcalc
      option%msdcalc   = msdcalc
      option%ngpcalc   = ngpcalc
      option%vhfcalc   = vhfcalc
      option%comcalc   = comcalc
      option%onlyz     = onlyz
      option%unwrap    = unwrap
      option%use_cond  = use_cond
      option%rcrange   = rcrange
      option%dt        = dt
      option%dr        = dr
      option%rmax      = rmax
      option%t_sparse  = t_sparse
      option%t_range   = t_range

      ! Check combinations
      !
      if (option%ngpcalc .and. .not. option%msdcalc) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("msdcalc should be true if ngpcalc = .true.")')
        stop
      end if

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
        option%nt_range = option%msdrange
      end if

      if (option%nt_sparse > 1) then
        option%nt_range = nint(dble(option%nt_range) / dble(option%nt_sparse))  
      end if

      option%nr = nint(option%rmax / option%dr)

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

!-----------------------------------------------------------------------
    subroutine check_combination(input, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input)  :: input
      type(s_option) :: option

      
      if (option%use_cond) then
        if (trim(input%fts) == "") then
          write(iw,'("Check_Combination> Error.")')
          write(iw,'("fcv should be specified in input_param")')
          write(iw,'("if use_cond = .true. in option_param.")')
          stop
        end if
      end if

    end subroutine check_combination
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
