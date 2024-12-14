!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_grid3d
  use mod_input
  use mod_output
  use mod_cv
  use mod_bootstrap
  implicit none

  ! constants
  !
  integer, parameter, public :: MaxDim   = 3
  integer, parameter, public :: MaxState = 100 
  integer, parameter, public :: MaxNcell = 100 

  ! structures
  !

  type :: s_option
    logical :: use_bootstrap          = .false.
    integer :: ndim                   = MaxDim
    integer :: nstart                 = 1 
    integer :: xyzcol(MaxDim)         = (/2, 3, 4/) 
    integer :: ng3(MaxDim)            = 100
    real(8) :: del(MaxDim)            = (/0.1d0, 0.1d0, 0.1d0/) 
    real(8) :: origin(MaxDim)         = 0.0d0
    real(8) :: temperature            = 298.0d0
    logical :: use_spline             = .false.
    integer :: spline_resolution      = 4
    logical :: space3D                = .false.
    integer :: nstep_system           = 0
    real(8) :: box_system(3)          = (/0.0d0, 0.0d0, 0.0d0/)

    ! define region of reference point
    logical :: use_refrange           = .false.
    real(8) :: refrange(2, 3)         = 0.0d0

    ! define cutoff range of calculating PMF
    logical :: use_cutrange           = .false.
    real(8) :: cutrange(2, 3)         = 0.0d0

    ! define range of each state
    logical :: calc_statepop              = .false.
    integer :: nstate                     = 2 
    real(8) :: staterange(2, 3, MaxState) = 0.0d0 

    ! generate matplotlib script for 2d-case
    logical :: gen_script_mpl2d       = .false.

    ! use radial coordinate
    logical :: use_radial             = .false.
    integer :: radcol                 = 0

    ! calculate g_integral = int_{glower}^{gupper} g(x) (available only for ndim = 1)
    logical :: calc_gintegral         = .false.
    real(8) :: gint_range(2)          = 0.0d0
    real(8) :: gint_refrange(2)       = 0.0d0

    ! use voronoi tessellation (ugly implementation) 
    logical :: use_voronoi            = .false.
    logical :: use_cellsearch         = .true.
    integer :: vr_ncell               = 5
    real(8) :: vr_separation(3)       = (/1.0d0, 1.0d0, 1.0d0/)
    real(8) :: vr_normvec(3)          = (/1.0d0, 1.0d0, 1.0d0/)
    real(8), allocatable :: vr_cellpos(:, :)

    ! skip calculation of PMF
    logical :: skip_calc              = .false.


  end type s_option

  ! subroutines
  !
  public  :: read_ctrl
  private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, cvinfo, mpl2d, bootopt)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_cvinfo),  intent(out) :: cvinfo
      type(s_mpl2d),   intent(out) :: mpl2d
      type(s_bootopt), intent(out) :: bootopt

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

        if (option%use_bootstrap) then
          call read_ctrl_bootstrap(iunit, bootopt)
        end if

        if (option%gen_script_mpl2d .and. option%ndim == 2) then
          call read_ctrl_matplotlib2d(iunit, mpl2d)
        end if
      close(iunit)

      open(iunit, file=trim(input%flist))
        call read_ctrl_cvinfo (iunit, cvinfo)
      close(iunit)

      if (trim(input%flist_weight) /= "") then
        open(iunit, file=trim(input%flist_weight))
          call read_ctrl_cvinfo_weight (iunit, cvinfo)
        close(iunit)
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

      logical :: use_bootstrap              = .false.
      integer :: ndim                       = 2
      integer :: nstart                     = 1 
      integer :: xyzcol(MaxDim)             = (/1, 2, 3/) 
      integer :: ng3(MaxDim)                = (/1, 1, 1/)
      real(8) :: del(MaxDim)                = (/1.0d0, 1.0d0, 1.0d0/) 
      real(8) :: origin(MaxDim)             = (/0.0d0, 0.0d0, 0.0d0/)
      real(8) :: temperature                = 298.0d0
      logical :: use_spline                 = .false.
      integer :: spline_resolution          = 4
      logical :: space3D                    = .false.
      integer :: nstep_system               = 0
      real(8) :: box_system(3)              = (/0.0d0, 0.0d0, 0.0d0/)
      logical :: use_refrange               = .false.
      real(8) :: refrange(2, 3)             = 0.0d0
      logical :: use_cutrange               = .false.
      real(8) :: cutrange(2, 3)             = 0.0d0
      logical :: calc_statepop              = .false.
      integer :: nstate                     = 2
      real(8) :: staterange(6*MaxState)     = 0.0d0
      logical :: gen_script_mpl2d           = .false.
      logical :: use_radial                 = .false.
      integer :: radcol                     = 0
      logical :: calc_gintegral             = .false.
      real(8) :: gint_range(2)              = 0.0d0
      real(8) :: gint_refrange(2)           = 0.0d0
      logical :: use_voronoi                = .false.
      logical :: use_cellsearch             = .false.
      integer :: vr_ncell                   = 5 
      real(8) :: vr_separation(3)           = (/1.0d0, 1.0d0, 1.0d0/)
      real(8) :: vr_normvec(3)              = (/1.0d0, 1.0d0, 1.0d0/)
      real(8) :: vr_cellpos(3, MaxNcell)    = 0.0d0 
      logical :: skip_calc                  = .false.


      integer :: i, j, k 


      namelist /option_param/ use_bootstrap, ndim, nstart, xyzcol, ng3, del, origin,  &
                              temperature, use_spline, spline_resolution,             &
                              space3D, nstep_system, box_system, use_refrange,        &
                              refrange, use_cutrange, cutrange,                       &
                              calc_statepop, nstate, staterange,                      &
                              gen_script_mpl2d, skip_calc, use_radial, radcol,        &
                              calc_gintegral, gint_range, gint_refrange,              &
                              use_voronoi, use_cellsearch, vr_ncell,                  &
                              vr_separation, vr_normvec, vr_cellpos

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("use_bootstrap     = ", a)')               get_tof(use_bootstrap) 
      write(iw,'("ndim              = ", i0)')              ndim 
      write(iw,'("nstart            = ", i0)')              nstart 
      write(iw,'("xyzcol            = ", 3(i0,2x))')        (xyzcol(i), i = 1, ndim) 
      write(iw,'("ng3               = ", 3(i0,2x))')        (ng3(i),    i = 1, ndim)
      write(iw,'("del               = ", 3(f15.7,2x))')     (del(i),    i = 1, ndim)
      write(iw,'("origin            = ", 3(f15.7,2x))')     (origin(i), i = 1, ndim)
      write(iw,'("temperature       = ", f15.7)')           temperature
      write(iw,'("use_spline        = ",a)')                get_tof(use_spline)
      write(iw,'("spline_resolution = ",i0)')               spline_resolution
      if (ndim == 2) then
        write(iw,'("gen_script_mpl2d  = ",a)')              get_tof(gen_script_mpl2d)
      else if (ndim == 3) then                                   
        write(iw,'("space3D           = ", a)')             get_tof(space3D)
        if (space3D) then
          write(iw,'("nstep_system      = ", i0)')          nstep_system 
          write(iw,'("box_system        = ", 3(f15.7,2x))') (box_system(i), i = 1, ndim)
        end if
      end if

      write(iw,'("use_radial        = ",a)')              get_tof(use_radial)
      if (use_radial) then
        write(iw,'("radcol            = ",i0)')           radcol
      end if

      write(iw,'("calc_gintegral    = ",a)')              get_tof(calc_gintegral)
      if (calc_gintegral) then
        write(iw,'("gint_range        = ",2f15.7)')       gint_range(1), gint_range(2)
        write(iw,'("gint_refrange     = ",2f15.7)')       gint_refrange(1), gint_refrange(2)
      end if

      write(iw,'("use_voronoi       = ", a)')             get_tof(use_voronoi)
      if (use_voronoi) then
        write(iw,'("use_cellsearch    = ", a)')             get_tof(use_cellsearch)
        write(iw,'("vr_ncell          = ", i0)')            vr_ncell
        write(iw,'("vr_separation     = ", 3(f15.7,2x))')   (vr_separation(i), i = 1, 3)
        write(iw,'("vr_normvec        = ", 3(f15.7,2x))')   (vr_normvec(i),    i = 1, 3)

        if (.not. use_cellsearch) then
          write(iw,'("vr_cellpos : ")')
          do i = 1, vr_ncell
            write(iw,'(i3,3f20.10)') i, (vr_cellpos(j, i), j = 1, 3) 
          end do
          write(iw,*)
        end if

      end if

      write(iw,'("skip_calc         = ",a)')              get_tof(skip_calc)

      ! memory allocation
      !
      allocate(option%vr_cellpos(3, vr_ncell))

      ! setup option variables
      !
      option%xyzcol = 0
      option%ng3    = 1
      option%del    = 1.0d0

      option%use_bootstrap               = use_bootstrap
      option%ndim                        = ndim
      option%nstart                      = nstart
      option%xyzcol(1:ndim)              = xyzcol(1:ndim)
      option%ng3(1:ndim)                 = ng3(1:ndim)
      option%del(1:ndim)                 = del(1:ndim)
      option%origin                      = origin
      option%temperature                 = temperature
      option%box_system                  = box_system
      option%use_refrange                = use_refrange
      option%refrange                    = refrange
      option%use_cutrange                = use_cutrange
      option%cutrange                    = cutrange

      option%calc_statepop               = calc_statepop
      option%nstate                      = nstate

      option%use_spline                  = use_spline
      option%spline_resolution           = spline_resolution
      option%space3D                     = space3D
      option%nstep_system                = nstep_system
      option%box_system                  = box_system
      option%gen_script_mpl2d            = gen_script_mpl2d
      option%use_radial                  = use_radial
      option%radcol                      = radcol
      option%calc_gintegral              = calc_gintegral
      option%gint_range                  = gint_range
      option%gint_refrange               = gint_refrange
      option%use_voronoi                 = use_voronoi
      option%use_cellsearch              = use_cellsearch
      option%vr_ncell                    = vr_ncell
      option%vr_separation               = vr_separation
      option%vr_normvec                  = vr_normvec
      option%vr_cellpos(1:3, 1:vr_ncell) = vr_cellpos(1:3, 1:vr_ncell)
      option%skip_calc                   = skip_calc
    
      k = 0
      do i = 1, nstate
        do j = 1, ndim
          option%staterange(1, j, i) = staterange(k + 1)
          option%staterange(2, j, i) = staterange(k + 2)
          k = k + 2 
        end do
      end do
      
      ! combination check
      !
      !if (option%ndim > 2) then
      !  write(iw,'("Read_Ctrl_Option> Error.")')
      !  write(iw,'("ndim > 2 is currently not supported.")')
      !  write(iw,'("sorry for the inconvenience.")')
      !  stop
      !end if

      if (option%use_radial) then
        if (option%radcol == 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("radcol (column id of radial coordinate) &
                     &should be specified")')
          stop
        end if
      end if

      if (option%use_voronoi) then
        !if (option%ndim /= 2) then
        !  write(iw,'("Read_Ctrl_Option> Error.")')
        !  write(iw,'("voronoi tessellation is available &
        !             &only if ndim = 2")')
        !  stop
        !end if

        if (.not. option%use_spline) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("voronoi tesselation is available &
                     &only if use_spline = .true.")')
        end if
      end if



    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
