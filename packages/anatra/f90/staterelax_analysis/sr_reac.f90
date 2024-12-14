!=======================================================================
module mod_reac
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_cv
  use mod_bootstrap
  use mod_random
  use mod_analyze_str

  ! constants
  !

  ! structures
  !

  ! subroutines
  !
  public :: get_reaction_info
  public :: check_transient
  public :: state_average
  public :: search_reaction

  contains

!-----------------------------------------------------------------------
    subroutine get_reaction_info(option, cvinfo, state)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_cvinfo), intent(in)    :: cvinfo
      type(s_state),  intent(inout) :: state(cvinfo%nfile)

      integer :: ifile
      integer :: nstep_sparse
      integer :: nBevents, nUBevents


      ! Average state
      !
!$omp parallel private(ifile, nstep_sparse), default(shared)
!$omp do
      do ifile = 1, cvinfo%nfile
        call state_average(option%nt_transient, state(ifile), nstep_sparse)
      end do
!$omp end do
!$omp end parallel

      ! Search reaction
      !
!$omp parallel private(ifile, nBevents, nUBevents),  &
!$omp          shared(state, nstep_sparse),          &
!$omp          default(shared)
!$omp do
      do ifile = 1, cvinfo%nfile
        call search_reaction(state(ifile), nstep_sparse, nBevents, nUBevents)
        state(ifile)%nBevents   = nBevents
        state(ifile)%nUBevents = nUBevents 
      end do
!$omp end do
!$omp end parallel

    end subroutine get_reaction_info
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine check_transient(output, option, cvinfo, state)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_cvinfo), intent(in)    :: cvinfo
      type(s_state),  intent(inout) :: state(cvinfo%nfile)

      integer                :: ntr_sta, ntr_interval, nntr
      integer                :: nfile

      integer                :: ifile, itr, ntr
      integer                :: istep, jstep, imin, imax
      integer                :: iunit
      real(8)                :: r, rcurr, dt_sparse
      integer                :: nstep_sparse
      character(len=MaxChar) :: fname

      integer, allocatable :: nBevents(:, :), nUBevents(:, :)
      integer, allocatable :: nstay(:, :)
      real(8), allocatable :: krate(:)


      ! Setup parameters
      !
      ntr_sta      = option%ntr_sta
      ntr_interval = option%ntr_interval
      nntr         = option%nntr
      nfile        = cvinfo%nfile

      ! Memory allocation
      !
      allocate(nBevents(0:nntr, nfile), nUBevents(0:nntr, nfile))
      allocate(nstay(0:nntr, nfile))
      allocate(krate(0:nntr))
      
      nBevents  = 0
      nUBevents = 0

!$omp parallel private(ifile, itr, ntr, nstep_sparse),    &
!$omp          shared(state, nBevents, nUBevents, nstay), &
!$omp          default(shared)
!$omp do
      do ifile = 1, nfile
        do itr = 0, nntr
          ntr = ntr_sta + itr * ntr_interval

          ! Average states using moving average scheme
          !
          call state_average(ntr, state(ifile), nstep_sparse)

          ! Search Binding/Unbinding
          !
          call search_reaction(state(ifile),          &
                               nstep_sparse,          &
                               nBevents(itr, ifile),  &
                               nUBevents(itr, ifile))

          ! Calculate time in which complex stays in Reaction zone
          !
          call calculate_staytime(state(ifile),        &
                                  option%reaczone_id,  &
                                  nstay(itr, ifile))

        end do
      end do
!$omp end do
!$omp end parallel

      ! Calculate 1st-order rate constant (krate) 
      !
      do itr = 0, nntr
        ntr = ntr_sta + itr * ntr_interval
        dt_sparse  = ntr * option%dt
        krate(itr) = dble(sum(nBevents(itr, :))) / (sum(nstay(itr, :)) * dt_sparse)
      end do

      ! Output Dt-dependency of krate
      !
      write(fname, '(a,".kcheck")') trim(output%fhead)
      call open_file(fname, iunit)
      do itr = 0, nntr
        ntr       = ntr_sta + itr * ntr_interval
        dt_sparse = ntr * option%dt

        write(iunit,'(f15.7,2x)', advance='no') dt_sparse
        write(iunit,'(f15.7)')                  krate(itr)
      end do
      close(iunit)

      ! Memory Deallocation
      !
      deallocate(nBevents, nUBevents)
      deallocate(nstay, krate)


    end subroutine check_transient
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine state_average(ntr, state, nstep_sparse)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)    :: ntr
      type(s_state),  intent(inout) :: state
      integer,        intent(out)   :: nstep_sparse

      integer :: istep, jstep, nstep
      integer :: imin, imax
      real(8) :: r


      state%avedata = 0.0d0
      nstep         = state%nstep

      jstep = 0
      do istep = ntr, nstep - ntr, ntr
        jstep = jstep + 1

        imin  = istep - ntr / 2
        imax  = istep + ntr / 2
        r     = sum(state%data_ts(imin:imax)) / dble(ntr + 1)

        state%avedata(jstep) = r
      end do

      imin = nstep - ntr / 2
      imax = nstep
      r    = sum(state%data_ts(imin:imax)) / dble(ntr /2 + 1)

      state%avedata(jstep + 1) = r
      nstep_sparse             = jstep + 1


    end subroutine state_average
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine search_reaction(state, nstep_sparse, nBevents, nUBevents)
!-----------------------------------------------------------------------
      implicit none

      type(s_state), intent(in)  :: state
      integer,       intent(in)  :: nstep_sparse
      integer,       intent(out) :: nBevents
      integer,       intent(out) :: nUBevents

      integer :: istep, jstep
      real(8) :: rcurr, rnext
      logical :: is_changed
     

      nBevents   = 0
      nUBevents  = 0

      istep      = 0
      is_changed = .true.
      rcurr      = state%data_ts(1)
      do while (is_changed)
        is_changed = .false.
        do jstep = istep + 1, nstep_sparse
          rnext = state%avedata(jstep)

          if (rcurr < 0.0d0 .and. rnext >= 0.0d0) then
            nBevents   = nBevents + 1
            rcurr      = rnext
            istep      = jstep
            is_changed = .true.
            exit
          else if (rcurr >= 0.0d0 .and. rnext < 0.0d0) then
            nUBevents  = nUBevents + 1
            rcurr      = rnext
            istep      = jstep
            is_changed = .true.
            exit
          end if
        end do

      end do

    end subroutine search_reaction
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calculate_staytime(state, target_id, nstay)
!-----------------------------------------------------------------------
      implicit none

      type(s_state), intent(in)  :: state
      integer,       intent(in)  :: target_id
      integer,       intent(out) :: nstay

      integer :: istep, nstep

      
      ! Setup parameters
      !
      nstep = state%nstep

      ! Calculate
      !
      nstay = 0
      do istep = 1, nstep
        if (state%data(istep) == target_id) then
          nstay = nstay + 1
        end if
      end do

    end subroutine calculate_staytime
!-----------------------------------------------------------------------

end module mod_reac
!=======================================================================
