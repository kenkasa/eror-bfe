!=======================================================================
module mod_analyze
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_dcdio
  use mod_xtcio
  use mod_traj
  use xdr, only: xtcfile


  ! constants
  !

  ! structures
  !
  type :: s_cv
    integer :: nstep
    integer :: ndim
    real(8), allocatable :: data(:,:)
  end type s_cv

  type :: s_state
    integer :: nstep
    integer :: ndim
    logical :: is_reacted = .false.
    integer, allocatable :: data(:) 
    real(8), allocatable :: hist(:)
  end type s_state

  ! subroutines
  !
  public :: analyze
  public :: read_ts
  public :: get_state
  public :: extract_dcd2dcd
  public :: extract_dcd2xtc
  public :: extract_xtc2xtc

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)   :: input
      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option

      type(s_cv)    :: cv, weight
      type(s_state) :: state
      integer       :: istep, next
      integer       :: iunit_weight
      integer       :: trajtype, trajtype_out
      logical       :: present_traj
      logical       :: present_weight


      ! Get input/output trajectory type
      !
      if (trim(input%ftraj(1)) /= "") then
        present_traj = .true.
      else
        present_traj = .false.
      end if

      if (present_traj) then
        call get_trajtype(input%ftraj(1), trajtype)
        call get_trajtype(output%ftraj, trajtype_out)
      end if

      ! Read cv files
      !
      write(iw,*)
      write(iw,'("Analyze> Read Time-Series Data file")')
      call read_ts(input%fts, option%ndim, cv)

      ! Read weight file if present
      !
      present_weight = .false.
      if (trim(input%fweight) /= "") then
        write(iw,*)
        write(iw,'("Analyze> Read Weight file")')
        call read_ts(input%fweight, 1, weight)
        present_weight = .true.
      end if

      ! Get state
      !
      write(iw,*)
      write(iw,'("Analyze> Get State")')
      call get_state(option%ndim,      &
                     option%state_def, &
                     cv,               &
                     state)

      next = 0
      do istep = 1, cv%nstep
        if (state%data(istep) == REACTIVE) then
          write(iw,'(i0,2x,f15.7)') istep, cv%data(1, istep)
          next = next + 1
        end if
      end do

      if (next == 0) then
        write(iw,'("No frames for reactive state.")')
        return
      end if

      ! Extract weight file (if ftraj is not specified)
      !
      if (.not. present_traj) then 

        if (present_weight) then
          call open_file(output%fweight, iunit_weight)
          do istep = 1, cv%nstep
            if (state%data(istep) == REACTIVE) then
              write(iunit_weight,'(i0,2x,e20.10)') &
                istep, weight%data(1, istep)
            end if
          end do
          close(iunit_weight)
        end if

        return

      end if 

      !  Extract snapshots
      !
      if (trajtype == TrajTypeDCD) then

        if (trajtype_out == TrajTypeDCD) then

          call extract_dcd2dcd(input,          &
                               output,         &
                               option,         &
                               cv,             &
                               weight,         &
                               state,          &
                               present_weight, &
                               next)

        else if (trajtype_out == TrajTypeXTC) then

          call extract_dcd2xtc(input,          &
                               output,         &
                               option,         &
                               cv,             &
                               weight,         &
                               state,          &
                               present_weight, &
                               next)

        end if

      else if (trajtype == TrajTypeXTC) then

        if (trajtype_out == TrajTypeDCD) then

          write(iw,'("Analyze> Error.")')
          write(iw,'("Sorry, convert xtc => dcd is not supported currently.")')
          stop

        else if (trajtype_out == TrajTypeXTC) then

          call extract_xtc2xtc(input,          &
                               output,         &
                               option,         &
                               cv,             &
                               weight,         &
                               state,          &
                               present_weight, &
                               next)

        end if

      end if


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine read_ts(fcv, ndim, cv)
!-----------------------------------------------------------------------
      implicit none

      character(len=MaxChar), intent(in)  :: fcv
      integer,                intent(in)  :: ndim
      type(s_cv),             intent(out) :: cv

      integer                :: istep, icv
      integer                :: nstep
      character(len=MaxChar) :: cdum


      open(10,file=trim(fcv))
      nstep = 0
      do while(.true.)
        read(10,'(a)',end=100) cdum
        nstep = nstep + 1
      end do

100   rewind 10

      cv%nstep = nstep
      allocate(cv%data(ndim, nstep)) 
      do istep = 1, nstep
        read(10,*) cdum, (cv%data(icv, istep), icv = 1, ndim) 
      end do
      close(10)

      cv%ndim = ndim
      

    end subroutine read_ts
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state(ndim, state_def, cv, state)
!-----------------------------------------------------------------------
      implicit none

      integer,                intent(in)  :: ndim
      real(8),                intent(in)  :: state_def(2, ndim, Nstate)
      type(s_cv),             intent(in)  :: cv
      type(s_state),          intent(out) :: state

      integer :: istep, istate, icv, ia
      integer :: nstep
      logical :: is_assigned
      real(8) :: wrk(ndim) 


      nstep       = cv%nstep

      ! for global
      allocate(state%data(nstep))

      do istep = 1, cv%nstep
        wrk(:) = cv%data(:, istep) 
        is_assigned = .false.
        do istate = 1, nstate
          if (is_assigned) then
            exit
          else
            ia = 0
            do icv = 1, ndim 
              if (wrk(icv) >= state_def(1, icv, istate) &
                .and. wrk(icv) < state_def(2, icv, istate)) then
                ia = ia + 1
              end if
            end do

            if (ia == ndim) then 
              is_assigned       = .true.
              state%data(istep) = StateInfo(istate)
            end if

          end if
        end do

        if (.not. is_assigned) then
          state%data(istep) = OTHERS 
        end if
      end do
      
      state%nstep = nstep
      state%ndim  = ndim

    end subroutine get_state 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine extract_dcd2dcd(input,          &
                               output,         &
                               option,         &
                               cv,             &
                               weight,         &
                               state,          &
                               present_weight, &
                               next)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in) :: input
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      type(s_cv),     intent(in) :: cv
      type(s_cv),     intent(in) :: weight
      type(s_state),  intent(in) :: state
      logical,        intent(in) :: present_weight
      integer,        intent(in) :: next

      type(s_dcd) :: dcdin, dcdout
      integer     :: iunit_in, iunit_out, iunit_weight
      integer     :: itraj, istep, istep_tot, jstep_tot, iwrite
      integer     :: natm, nstep
      logical     :: traj_reactive = .false.


      call dcd_open(input%ftraj(1), iunit_in)
      call dcd_open(output%ftraj,   iunit_out)

      call dcd_read_header(iunit_in, dcdin)
      call dcd_close(iunit_in, input%ftraj(1))

      if (present_weight) then
        call open_file(output%fweight, iunit_weight)
      end if

      natm = dcdin%natm

      call alloc_dcd(natm, 1, dcdin)
      call alloc_dcd(natm, 1, dcdout)

      dcdout%natm       = natm
      dcdout%nstep      = next 
      dcdout%dcdinfo    = dcdin%dcdinfo
      dcdout%dcdinfo(1) = next
      dcdout%dcdinfo(4) = next
       
      call dcd_write_header(iunit_out, dcdout)

      istep_tot = 0
      iwrite    = 0
      do itraj = 1, input%ntraj
        call dcd_open(input%ftraj(itraj), iunit_in)
        call dcd_read_header(iunit_in, dcdin)

        nstep = dcdin%nstep

        ! Check the existence of reactive state in trajectory
        !
        jstep_tot     = istep_tot
        traj_reactive = .false.
        do istep = 1, nstep
          jstep_tot = jstep_tot + 1
          if (state%data(jstep_tot) == REACTIVE) then
            traj_reactive = .true.
            exit
          end if
        end do

        if (.not. traj_reactive) then
          write(iw,'("Skip reading ",a)') trim(input%ftraj(itraj))
          istep_tot = istep_tot + nstep
        else
          do istep = 1, nstep
            istep_tot = istep_tot + 1
         
            call read_dcd_oneframe(iunit_in, dcdin)
         
            if (state%data(istep_tot) == REACTIVE) then
              write(iw,'("Step ", i8, " : Extracted")') istep_tot
         
              dcdout%box(1:3, 1)           = dcdin%box(1:3, 1)
              dcdout%coord(1:3, 1:natm, 1) = dcdin%coord(1:3, 1:natm, 1)
              call write_dcd_oneframe(iunit_out, 1, dcdout)
         
              if (present_weight) then
                write(iunit_weight,'(i0,2x,e20.10)') &
                  istep_tot, weight%data(1, istep_tot)
              end if
            end if
         
          end do
        end if

        call dcd_close(iunit_in, input%ftraj(itraj)) 

      end do

      call dcd_close(iunit_out)
      call dealloc_dcd(dcdin)
      call dealloc_dcd(dcdout)

      if (present_weight) then
        close(iunit_weight)
      end if 

    end subroutine extract_dcd2dcd
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine extract_dcd2xtc(input,          &
                               output,         &
                               option,         &
                               cv,             &
                               weight,         &
                               state,          &
                               present_weight, &
                               next)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in) :: input
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      type(s_cv),     intent(in) :: cv
      type(s_cv),     intent(in) :: weight
      type(s_state),  intent(in) :: state
      logical,        intent(in) :: present_weight
      integer,        intent(in) :: next

      real(8), parameter :: ang2nm = 0.1d0

      type(s_dcd)   :: dcdin
      type(xtcfile) :: xtcout 
      integer       :: iunit_in, iunit_out, iunit_weight
      integer       :: itraj, istep, istep_tot, jstep_tot, iwrite, nstep
      integer       :: natm
      real(8)       :: box(3, 3)
      logical       :: traj_reactive = .false.

      real(4), allocatable :: coord(:, :)


      ! Setup Input DCD / Output XTC
      !
      call dcd_open(input%ftraj(1), iunit_in)
      call dcd_read_header(iunit_in, dcdin)
      natm = dcdin%natm
      call alloc_dcd(natm, 1, dcdin)
      call xtcout%init(trim(output%ftraj), 'w')

      call dcd_close(iunit_in)

      if (present_weight) then
        call open_file(output%fweight, iunit_weight)
      end if

      ! Allocate memory
      !
      allocate(coord(1:3, 1:natm))

      istep_tot = 0
      iwrite    = 0
      do itraj = 1, input%ntraj
        call dcd_open(input%ftraj(itraj), iunit_in)
        call dcd_read_header(iunit_in, dcdin)
        nstep = dcdin%nstep

        ! Check the existence of reactive state in trajectory
        !
        jstep_tot     = istep_tot
        traj_reactive = .false. 
        do istep = 1, nstep
          jstep_tot = jstep_tot + 1
          if (state%data(jstep_tot) == REACTIVE) then
            traj_reactive = .true.
            exit
          end if
        end do

        if (.not. traj_reactive) then
          write(iw,'("Skip reading ",a)') trim(input%ftraj(itraj))
          istep_tot = istep_tot + nstep
        else

          do istep = 1, nstep
            istep_tot = istep_tot + 1
         
            call read_dcd_oneframe(iunit_in, dcdin)
         
            if (state%data(istep_tot) == REACTIVE) then
              write(iw,'(" Step ",i8, ": Extracted")') istep 
             
              box = 0.0
              box(1, 1) = dcdin%box(1, 1) * ang2nm
              box(2, 2) = dcdin%box(2, 1) * ang2nm
              box(3, 3) = dcdin%box(3, 1) * ang2nm
              coord(1:3, 1:natm) = dcdin%coord(1:3, 1:natm, 1) * ang2nm
         
              call xtcout%write(natm,               &
                                0,                  &
                                0.0,                &
                                real(box),          &
                                coord(1:3, 1:natm), &
                                real(1000.0d0))
         
              if (present_weight) then
                write(iunit_weight,'(i0,2x,e20.10)') &
                  istep_tot, weight%data(1, istep_tot)
              end if
            end if
          end do

        end if

        call dcd_close(iunit_in)

      end do

      if (present_weight) then
        close(iunit_weight)
      end if 

      ! Deallocate memory
      !
      call xtcout%close
      call dealloc_dcd(dcdin)

      deallocate(coord)

    end subroutine extract_dcd2xtc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine extract_xtc2xtc(input,          &
                               output,         &
                               option,         &
                               cv,             & 
                               weight,         &
                               state,          &
                               present_weight, &
                               next)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in) :: input
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      type(s_cv),     intent(in) :: cv
      type(s_cv),     intent(in) :: weight
      type(s_state),  intent(in) :: state
      logical,        intent(in) :: present_weight
      integer,        intent(in) :: next

      type(xtcfile) :: xtcin, xtcout
      integer       :: iunit_weight
      integer       :: itraj, istep, istep_tot, jstep
      integer       :: natm
      logical       :: is_end


      ! Setup Input XTC / Output XTC
      !
      call xtcin%init(trim(input%ftraj(1)))
      call xtcout%init(trim(output%ftraj), 'w')
      natm = xtcin%natoms

      if (present_weight) then
        call open_file(output%fweight, iunit_weight)
      end if

      ! Start extract
      !
      istep_tot = 0
      do itraj = 1, input%ntraj
       
        is_end = .false.
        do while (.not. is_end)
          call xtcin%read

          if (xtcin%STAT /= 0) then
            is_end = .true.
            exit
          else
            if (state%data(istep_tot + 1) == REACTIVE) then
              write(iw,'(" Step ",i8, ": Extracted")') istep_tot + 1

              call xtcout%write(xtcin%natoms,   &
                                istep_tot + 1,  &
                                xtcin%time,     &
                                xtcin%box,      &
                                xtcin%pos,      &
                                xtcin%prec)
              
              if (present_weight) then
                write(iunit_weight,'(i0,2x,e20.10)') &
                  istep_tot, weight%data(1, istep_tot + 1)
              end if

            end if

            istep_tot = istep_tot + 1 

          end if

        end do

        call xtcin%close

      end do

      ! Start extract
      !
!      jstep = 0 
!      do istep = 1, cv%nstep
!
!        call xtcin%read
!
!        if (state%data(istep) == REACTIVE) then
!          write(iw,'(" Step ",i8, ": Extracted")') istep
!
!          jstep = jstep + 1
!          call xtcout%write(xtcin%natoms,   &
!                            jstep,          &
!                            xtcin%time,     &
!                            xtcin%box,      &
!                            xtcin%pos,      &
!                            xtcin%prec)
!                            
!        end if
!
!      end do

      ! Deallocate memory
      !
      !call xtcin%close

      if (present_weight) then
        close(iunit_weight)
      end if

      call xtcout%close

    end subroutine extract_xtc2xtc
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================