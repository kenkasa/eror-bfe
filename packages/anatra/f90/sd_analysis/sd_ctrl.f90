!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_cv
  use mod_traj
  use mod_com
  implicit none

  ! constants
  !
  integer,      parameter, public :: NdimMax = 3

  integer,      parameter, public :: CenterTypeZERO = 1
  integer,      parameter, public :: CenterTypeHALF = 2
  character(*), parameter, public :: CenterTypes(2) = (/'ZERO', &
                                                        'HALF'/)

  !integer,      parameter, public :: CoMModeRESIDUE = 1
  !integer,      parameter, public :: CoMModeWHOLE   = 2 
  !character(*), parameter, public :: CoMMode(2)     = (/'RESIDUE   ',&
  !                                                      'WHOLE     '/)

  integer,      parameter, public :: Nstate   = 1
  integer,      parameter, public :: REACTIVE = 1
  integer,      parameter, public :: OTHERS   = 2
  integer,      parameter, public :: StateInfo(Nstate) = (/1/)

  ! structures
  !

  type :: s_option
    integer :: mode                   = CoMModeRESIDUE
    integer :: ng3(3)                 = (/100, 100, 100/)
    real(8) :: del(3)                 = (/0.1d0, 0.1d0, 0.1d0/) 
    real(8) :: origin(3)              = 0.0d0

    logical :: use_pbcwrap            = .false.
    integer :: centertype             = CenterTypeZERO 

    logical :: use_spline             = .false.
    integer :: spline_resolution      = 4

    logical :: is_whole_snap          = .true.
    integer :: nstep_whole            = 0

    logical :: use_restriction        = .false.
    logical :: normalize_reacsnap     = .false.
    integer :: ndim                   = 1

    logical :: use_weight             = .false.

    real(8) :: count_threshold        = 1.0d-10

    real(8), allocatable :: react_range(:, :)  
    real(8), allocatable :: state_def(:, :)
  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, trajopt, cvinfo)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_trajopt), intent(out) :: trajopt 
      type(s_cvinfo),  intent(out) :: cvinfo 

      integer                      :: i, iunit_cv, iunit_fwl, nfile
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

      if (option%use_restriction) then
        if (trim(input%flist) == "") then
          allocate(cvinfo%fcv(1))
          cvinfo%nfile  = 1
          cvinfo%fcv(1) = trim(input%fcv)
        else
          nfile = 0
          call open_file(input%flist, iunit_cv)
          do while (.true.)
            read(iunit_cv,*,end=100)
            nfile = nfile + 1
          end do

100       rewind iunit_cv
          cvinfo%nfile = nfile
          allocate(cvinfo%fcv(nfile))

          do i = 1, nfile
            read(iunit_cv,'(a)') cvinfo%fcv(i)
          end do

          close (iunit_cv)

        end if
      end if

      if (option%use_weight) then
        if (trim(input%flist_weight) /= "") then
          call open_file(input%flist_weight, iunit_fwl)
          call read_ctrl_cvinfo_weight(iunit_fwl, cvinfo)
          close(iunit_fwl)
        else
          write(iw,'("Read_Ctrl> Error.")')
          write(iw,'("flist_weight should be specified if use_weight = .true.")')
          stop
        end if
      end if


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


      character(len=MaxChar) :: mode                    = 'RESIDUE' 
      integer                :: ng3(3)                  = (/50, 50, 50/)
      real(8)                :: del(3)                  = (/1.0d0, 1.0d0, 1.0d0/) 
      real(8)                :: origin(3)               = (/0.0d0, 0.0d0, 0.0d0/)
      ! used for pbc wrap
      logical                :: use_pbcwrap             = .false.
      character(len=MaxChar) :: centertype              = 'ZERO'
      ! used for restricted sdf
      logical                :: use_restriction         = .false.
      integer                :: ndim                    = 1
      real(8)                :: react_range(2, NdimMax) = 0.0d0  
      ! used for spline
      logical                :: use_spline              = .false.
      integer                :: spline_resolution       = 4
      ! used for conditional                            
      logical                :: is_whole_snap           = .true.
      logical                :: normalize_reacsnap      = .false.
      integer                :: nstep_whole             = 0

      logical                :: use_weight              = .false.

      real(8)                :: count_threshold         = 1.0d-10


      integer                :: iopt, ierr
      integer                :: i, j, k 


      namelist /option_param/ mode,                   &
                              ng3,                    &
                              del,                    &
                              origin,                 &
                              use_pbcwrap,            &
                              centertype,             &
                              use_restriction,        &
                              ndim,                   &
                              react_range,            &
                              use_spline,             &
                              spline_resolution,      &
                              normalize_reacsnap,     &
                              is_whole_snap,          & 
                              nstep_whole,            &
                              use_weight,             &
                              count_threshold

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("mode              = ", a)')               trim(mode)
      write(iw,'("ng3               = ", 3(i0,2x))')        (ng3(i),    i = 1, 3)
      write(iw,'("del               = ", 3(f15.7,2x))')     (del(i),    i = 1, 3)
      write(iw,'("origin            = ", 3(f15.7,2x))')     (origin(i), i = 1, 3)

      if (use_pbcwrap) then
        write(iw,'("use_pbcwrap       = ", a)')             get_tof(use_pbcwrap)
        write(iw,'("centertype        = ", a)')             trim(centertype)
      end if

      write(iw,'("use_restriction   = ", a)')               get_tof(use_restriction)
      write(iw,'("normalize_reacsnap= ", a)')               get_tof(normalize_reacsnap)
      if (use_restriction) then
        write(iw,'("ndim             = ", i0)')             ndim
        do i = 1, ndim
          write(iw,'(es15.7, " <= component ",i0," < ",es15.7)') &
                  react_range(1, i), i, react_range(2, i)
        end do 
      end if

      write(iw,'("use_spline        = ", a)')               get_tof(use_spline)
      write(iw,'("spline_resolution = ", i0)')              spline_resolution
      write(iw,'("is_whole_snap     = ", a)')               get_tof(is_whole_snap)
      if (.not. is_whole_snap) then
        write(iw,'("nstep_whole       = ", i0)')            nstep_whole 
      end if

      write(iw,'("use_weight          = ", a)')             get_tof(use_weight)

      write(iw,'("count_threshold     = ", e15.7)')         count_threshold

      iopt = get_opt(mode, CoMMode, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("mode = ",a," is not available.")') trim(mode)
        stop
      end if
      option%mode              = iopt

      iopt = get_opt(centertype, CenterTypes, ierr)
      if (ierr /= 0) then
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("center = ",a," is not available.")') trim(centertype)
        stop
      end if
      option%centertype        = iopt


      allocate(option%react_range(1:2, ndim))
      allocate(option%state_def(2, ndim))

      option%ng3                   = ng3
      option%del                   = del
      option%origin                = origin
      option%use_pbcwrap           = use_pbcwrap
      option%use_restriction       = use_restriction
      option%normalize_reacsnap    = normalize_reacsnap
      option%ndim                  = ndim
      do i = 1, ndim
        option%react_range(1:2, i) = react_range(1:2, i)
        option%state_def(1:2, i)   = react_range(1:2, i) 
      end do 
      option%use_spline            = use_spline
      option%spline_resolution     = spline_resolution
      option%is_whole_snap         = is_whole_snap 
      option%nstep_whole           = nstep_whole
      option%use_weight            = use_weight
      option%count_threshold       = count_threshold

      ! Combination check
      !
      if (use_pbcwrap) then
        write(iw,*)
        write(iw,'("Read_Ctrl_Option> Remark")')
        write(iw,'("use_pbcwrap = .true. detected.")')
        write(iw,'("Please do not use this option if your trajectory is rotated by fitting.")')
      end if

      if (.not. is_whole_snap .and. use_weight) then
        write(iw,*)
        write(iw,'("Read_Ctrl_Option> Error.")')
        write(iw,'("is_whole_snap = .false. and use_weight = .true. can not be used at the same time.")')
        stop
      end if


    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

end module
!=======================================================================
