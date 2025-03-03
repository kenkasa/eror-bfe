!=======================================================================
module mod_input
!=======================================================================
  use mod_util
  use mod_const

  implicit none

  ! parameters
  !
  integer, parameter :: NInputParm        = 17

  integer, parameter :: InputFXYZ         = 1
  integer, parameter :: InputFTRAJ        = 2
  integer, parameter :: InputFPDB         = 3
  integer, parameter :: InputFDXS         = 4
  integer, parameter :: InputFTS          = 5
  integer, parameter :: InputFCV          = 6
  integer, parameter :: InputFWEIGHT      = 7
  integer, parameter :: InputFDX          = 8
  integer, parameter :: InputFCHARGE      = 9
  integer, parameter :: InputFLIST        = 10
  integer, parameter :: InputFLIST_CV     = 11
  integer, parameter :: InputFLIST_TRAJ   = 12
  integer, parameter :: InputFLIST_WEIGHT = 13
  integer, parameter :: InputFPRMTOP      = 14
  integer, parameter :: InputFPRMTOP2     = 15
  integer, parameter :: InputFANAPARM     = 16
  integer, parameter :: InputFANAPARM2    = 17

  ! structures
  !
  type :: s_input
    character(len=MaxChar) :: fxyz(MaxTraj)
    logical                :: specify_fxyz

    character(len=MaxChar) :: fdcd(MaxTraj)
    !logical                :: specify_fdcd

    character(len=MaxChar) :: ftraj(MaxTraj)  ! can be used both for dcd and xtc
    logical                :: specify_ftraj 

    character(len=MaxChar) :: fpdb(MaxTraj)
    logical                :: specify_fpdb

    character(len=MaxChar) :: fcvs(MaxTraj)
    !logical                :: specify_fcvs

    character(len=MaxChar) :: fdxs(MaxTraj)
    logical                :: specify_fdxs

    character(len=MaxChar) :: fts
    logical                :: specify_fts

    character(len=MaxChar) :: fcv
    logical                :: specify_fcv

    character(len=MaxChar) :: fweight
    logical                :: specify_fweight

    character(len=MaxChar) :: fdx
    logical                :: specify_fdx

    character(len=MaxChar) :: fcharge
    logical                :: specify_fcharge

    character(len=MaxChar) :: flist         ! TODO: should be changed to flist_cv for public release
    logical                :: specify_flist

    character(len=MaxChar) :: flist_cv
    logical                :: specify_flist_cv

    character(len=MaxChar) :: flist_traj
    logical                :: specify_flist_traj

    character(len=MaxChar) :: flist_weight
    logical                :: specify_flist_weight

    character(len=MaxChar) :: fprmtop
    logical                :: specify_fprmtop

    character(len=MaxChar) :: fprmtop2
    logical                :: specify_fprmtop2

    character(len=MaxChar) :: fanaparm
    logical                :: specify_fanaparm

    character(len=MaxChar) :: fanaparm2
    logical                :: specify_fanaparm2

    logical                :: present_ftraj
    integer                :: ntraj
  end type s_input

  type :: s_input_check
    logical :: require(NInputParm)
  end type s_input_check

  ! subroutines
  !
  public :: read_ctrl_input
  public :: initialize_input_check

  contains
!-----------------------------------------------------------------------
    subroutine read_ctrl_input(iunit, input, myrank)
!-----------------------------------------------------------------------
      implicit none
!
      integer,           intent(in)  :: iunit
      type(s_input),     intent(out) :: input
      integer, optional, intent(in)  :: myrank

      character(len=MaxChar) :: fdcd0 = "", fdcd1 = ""
      character(len=MaxChar) :: fdcd2 = ""
      character(len=MaxChar) :: fxyz(1:MaxTraj) = ""
      character(len=MaxChar) :: fdcd(1:MaxTraj) = ""
      character(len=MaxChar) :: ftraj(1:MaxTraj) = ""
      character(len=MaxChar) :: fpdb(1:MaxTraj) = ""
      character(len=MaxChar) :: fcvs(1:MaxTraj) = ""
      character(len=MaxChar) :: fdxs(1:MaxTraj) = ""
      character(len=MaxChar) :: fts = "", fcv = "", fweight = "", fdx = "", fcharge = ""

      character(len=MaxChar) :: flist         = ""  ! TODO: should be changed to flist_cv for public release 
      character(len=MaxChar) :: flist_cv      = ""
      character(len=MaxChar) :: flist_traj    = ""
      character(len=MaxChar) :: flist_weight  = ""
      character(len=MaxChar) :: fprmtop       = ""
      character(len=MaxChar) :: fprmtop2      = ""
      character(len=MaxChar) :: fanaparm      = ""
      character(len=MaxChar) :: fanaparm2     = ""

      integer :: i, ndcd, ntraj, nxyz, npdb, ndx
      integer :: irank
      integer :: iunit_tmp

      namelist /input_param/ fdcd0,           &
                             fdcd1,           &
                             fdcd2,           &
                             fxyz,            &
                             fdcd,            &
                             ftraj,           &
                             fpdb,            &
                             fdxs,            &
                             fts,             &
                             fcv,             &
                             fweight,         &
                             fdx,             &
                             fcharge,         &
                             flist,           &
                             flist_cv,        &
                             flist_traj,      &
                             flist_weight,    &
                             fprmtop,         &
                             fprmtop2,        &
                             fanaparm,        &
                             fanaparm2


      if (present(myrank)) then
        irank = myrank
      else
        irank = 0
      end if

      ! Initialize
      !
      fdcd0        = ""
      fdcd1        = ""
      fxyz         = ""
      fdcd         = ""
      ftraj        = ""
      fpdb         = ""
      fdxs         = ""
      fts          = ""
      fcv          = ""
      fweight      = ""
      fdx          = ""
      fcharge      = ""
      flist        = ""
      flist_cv     = ""
      flist_traj   = ""
      flist_weight = ""
      fprmtop      = ""
      fprmtop2     = ""
      fanaparm     = ""
      fanaparm2    = ""

      input%specify_fxyz         = .false.
      input%specify_ftraj        = .false.
      input%specify_fpdb         = .false.
      !specify_fcvs         = .false.
      input%specify_fdxs         = .false.
      input%specify_fts          = .false.
      input%specify_fcv          = .false.
      input%specify_fweight      = .false.
      input%specify_fdx          = .false.
      input%specify_fcharge      = .false.
      input%specify_flist        = .false.
      input%specify_flist_cv     = .false.
      input%specify_flist_traj   = .false.
      input%specify_flist_weight = .false.
      input%specify_fprmtop      = .false.
      input%specify_fprmtop2     = .false.
      input%specify_fanaparm     = .false.
      input%specify_fanaparm2    = .false.


      ! Read INPUT_PARAM namelist
      !
      rewind iunit
      read(iunit,input_param)

      ! Print out input variables
      !
      if (irank == 0) then
        write(iw,*)
        write(iw,'(">> Input section parameters")')
      end if

      ! Setup fdcd
      ! TODO: should be removed for public release
      !
      ndcd = 0
      do i = 1, MaxTraj
        if (trim(fdcd(i)) /= "") then
          !write(iw,'("fdcd = ", a)') trim(fdcd(i))
          ndcd = ndcd + 1
        end if
      end do

      ! Setup ftraj 
      !
      ! - From ftraj directly
      !
      input%specify_ftraj = is_specified(ftraj(1))
      input%present_ftraj = .false.
      ntraj = 0
      do i = 1, MaxTraj
        if (trim(ftraj(i)) /= "") then
          !write(iw,'("ftraj = ",a)') trim(ftraj(i))
          ntraj = ntraj + 1
        end if
      end do

      if (ntraj /= 0) then
        input%present_ftraj = .true.
        input%ntraj         = ntraj
      end if


      ! - From flist_traj
      !
      input%specify_flist_traj = is_specified(flist_traj)

      if (trim(flist_traj) /= "") then

        !if (ntraj /= 0) then
        !  write(iw,*)
        !  write(iw,'("Read_Ctrl_Input> Remark")')
        !  write(iw,'("Trajectory list specified in ftraj is replaced by that in flist_traj.")')
        !end if

        ntraj = 0
        call open_file(flist_traj, iunit_tmp)

        do while (.true.)
          read(iunit_tmp,*,end=100)
          ntraj = ntraj + 1
        end do

100     rewind iunit_tmp

        do i = 1, ntraj
          read(iunit_tmp,'(a)')      ftraj(i)
          !write(iw,'("ftraj = ",a)') trim(ftraj(i))
        end do

        input%present_ftraj = .true.
        input%ntraj         = ntraj

        close(iunit_tmp)
      end if

      ! Setup fxyz
      !
      input%specify_fxyz = is_specified(fxyz(1))
      nxyz = 0
      do i = 1, MaxTraj
        if (trim(fxyz(i)) /= "") then
          !write(iw,'("fxyz = ", a)') trim(fxyz(i))
          nxyz = nxyz + 1
        end if
      end do

      ! Setup fpdb
      !
      input%specify_fpdb = is_specified(fpdb(1))
      npdb = 0
      do i = 1, MaxTraj
        if (trim(fpdb(i)) /= "") then
          !write(iw,'("fpdb = ", a)') trim(fpdb(i))
          npdb = npdb + 1
        end if
      end do

      ! Setup fdx
      !
      input%specify_fdxs = is_specified(fdxs(1))
      ndx = 0
      do i = 1, MaxTraj
        if (trim(fdxs(i)) /= "") then
          !write(iw,'("fdxs = ", a)') trim(fdxs(i))
          ndx = ndx + 1
        end if
      end do

      input%specify_fts       = is_specified(fts)
      input%specify_fcv       = is_specified(fcv)
      input%specify_fweight   = is_specified(fweight)
      input%specify_fdx       = is_specified(fdx)
      input%specify_fcharge   = is_specified(fcharge)
      input%specify_flist     = is_specified(flist)
      input%specify_fprmtop   = is_specified(fprmtop)
      input%specify_fprmtop2  = is_specified(fprmtop2)
      input%specify_fanaparm  = is_specified(fanaparm)
      input%specify_fanaparm2 = is_specified(fanaparm2)


      ! Print out
      !
      if (irank == 0) then

        if (input%specify_ftraj) then
          if (input%specify_flist_traj) then
            write(iw,'("Trajectory list specified in ftraj is &
                       &replaced by that in flist_traj.")')
          else
            do i = 1, MaxTraj
              if (trim(ftraj(i)) /= "") & 
                write(iw,'("ftraj = ", a)') trim(ftraj(i))
            end do
          end if
        end if

        if (input%specify_flist_traj) then
          input%specify_ftraj = .true.
          write(iw,'("flist_traj   = ", a)') trim(flist_traj) 
          do i = 1, MaxTraj
            if (trim(ftraj(i)) /= "") &
              write(iw,'("  ftraj      = ", a)') trim(ftraj(i))
          end do
        end if

        if (input%specify_fxyz) then
          do i = 1, MaxTraj
            if (trim(fxyz(i)) /= "") &
              write(iw,'("fxyz     = ", a)') trim(fxyz(i))
          end do
        end if

        if (input%specify_fpdb) then
          do i = 1, MaxTraj
            if (trim(fpdb(i)) /= "") &
              write(iw,'("fpdb     = ", a)') trim(fpdb(i))
          end do
        end if

        if (input%specify_fdxs) then
          do i = 1, MaxTraj
            if (trim(fdxs(i)) /= "") &
              write(iw,'("fdxs     = ", a)') trim(fdxs(i))
          end do
        end if

        if (input%specify_fts) & 
          write(iw,'("fts          = ", a)') trim(fts)
        if (input%specify_fcv) &
          write(iw,'("fcv          = ", a)') trim(fcv)
        if (input%specify_fweight) &
          write(iw,'("fweight      = ", a)') trim(fweight)
        if (input%specify_fdx) &
          write(iw,'("fdx          = ", a)') trim(fdx)
        if (input%specify_fcharge) &
          write(iw,'("fcharge      = ", a)') trim(fcharge)
        if (input%specify_flist) &
          write(iw,'("flist        = ", a)') trim(flist)
        !if (trim(flist_traj) /= "") &
        !  write(iw,'("flist_traj   = ", a)') trim(flist_traj)
        if (input%specify_flist_weight) &
          write(iw,'("flist_weight = ", a)') trim(flist_weight)
        if (input%specify_fprmtop) &
          write(iw,'("fprmtop      = ", a)') trim(fprmtop)
        if (input%specify_fprmtop2) &
          write(iw,'("fprmtop2     = ", a)') trim(fprmtop2)
        if (input%specify_fanaparm) &
          write(iw,'("fanaparm     = ", a)') trim(fanaparm)
        if (input%specify_fanaparm2) &
          write(iw,'("fanaparm2    = ", a)') trim(fanaparm2)
      end if
      ! 
      !

      ! Send values of input parameters to input structure
      !
      input%fxyz          = fxyz
      input%fdcd          = fdcd
      input%ftraj         = ftraj
      input%fpdb          = fpdb
      input%fts           = fts
      input%fcv           = fcv
      input%fweight       = fweight
      input%fdx           = fdx
      input%fdxs          = fdxs
      input%fcharge       = fcharge
      input%flist         = flist
      input%flist_traj    = flist_traj
      input%flist_weight  = flist_weight
      input%fprmtop       = fprmtop
      input%fprmtop2      = fprmtop2
      input%fanaparm      = fanaparm 
      input%fanaparm2     = fanaparm2

      ! For backward compatibility
      ! fdcd0 => input%fdcd(1), fdcd1 => input%fdcd(2)
      ndcd = 0
      do i = 1, MaxTraj
        if (trim(fdcd(i)) /= "") then
          ndcd = ndcd + 1
        end if
      end do

      if (ndcd == 0) then
        if (trim(fdcd0) /= "") then
          if (irank == 0) then
            write(iw,'("fdcd = ", a)') trim(fdcd0)
          end if
          input%fdcd(1) = fdcd0
        end if

        if (trim(fdcd1) /= "") then
          if (irank == 0) then
            write(iw,'("fdcd = ", a)') trim(fdcd1)
          end if
          input%fdcd(2) = fdcd1
        end if

        if (trim(fdcd2) /= "") then
          if (irank == 0) then
            write(iw,'("fdcd = ", a)') trim(fdcd2)
          end if
          input%fdcd(3) = fdcd2
        end if
      end if

      if (trim(fdx) /= "") then
        if (trim(fdxs(1)) == "") then
          fdxs(1) = trim(fdx)
        end if
      end if

    end subroutine read_ctrl_input
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine initialize_input_check(ic)
!-----------------------------------------------------------------------
      implicit none
      
      type(s_input_check), intent(out) :: ic

      integer :: i

      do i = 1, NInputParm
        ic%require(i) = .false.
      end do

    end subroutine initialize_input_check
!-----------------------------------------------------------------------

end module mod_input
!=======================================================================
