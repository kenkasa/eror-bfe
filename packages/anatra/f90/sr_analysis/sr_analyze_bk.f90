!=======================================================================
module mod_analyze
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
  use mod_reac
  use mod_tcf

  ! constants
  !

  ! structures
  !

  ! subroutines
  !
  public  :: analyze 
  private :: analyze_bootstrap
  private :: calc_kine_bootstrap
  private :: calc_kineave_bootstrap

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, cvinfo, bootopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option
      type(s_cvinfo),  intent(in)    :: cvinfo
      type(s_bootopt), intent(inout) :: bootopt

      type(s_cv),    allocatable :: cv(:)
      type(s_state), allocatable :: state(:)


      integer :: ifile, istep, itr, ntr, ncut
      integer :: nBevents, nUBevents
      integer :: nreact, nreact_final, nstay, nstay_bound
      integer :: nstep_sparse
      integer :: nbound, nubound
      integer :: nthreads, thread_id, nwork, ifinish, iper
      real(8) :: rpinit, kr, t, dth
      real(8) :: tau_r, tau_d
      real(8) :: progress

      integer, allocatable :: nbcheck(:, :), nubcheck(:, :)
      real(8), allocatable :: rp(:), norm_tot(:) 
      character(len=MaxChar) :: frp


      ! Memory allocation
      !
      allocate(cv(cvinfo%nfile), state(cvinfo%nfile))

      ! Read cv files
      !
      write(iw,*)
      write(iw,'("Analyze> Read CV file")')

      iper = 1
      do ifile = 1, cvinfo%nfile

        ! Check progress
        !
        progress = ifile / dble(cvinfo%nfile)
        if (progress >= iper * 0.1d0) then
          write(iw,'("  Progress : ", f6.2, "%")') progress * 100.0d0
          iper = iper + 1
        end if

        ! Read
        !
        call read_cv(cvinfo%fcv(ifile),         &
                     option%ndim,               &
                     cv(ifile),                 &
                     nstep_read = option%nstep)
      end do

      ! Get state info
      !
      write(iw,*)
      write(iw,'("Analyze> Get State")')

!$omp parallel private(ifile, nthreads, thread_id, nwork, progress, ifinish), & 
!$omp          default(shared)
      nthreads  = omp_get_num_threads()
      thread_id = omp_get_thread_num()
      nwork     = cvinfo%nfile / nthreads + 1
      ifinish   = 0
!$omp do
      do ifile = 1, cvinfo%nfile

        call get_state(option, cv(ifile), state(ifile))

        if (thread_id == 0) then
          ifinish = ifinish + 1
          progress = ifinish / dble(nwork) * 100.0d0
          write(iw,'("  Progress : ", f6.2, "%")') progress
        end if

      end do
!$omp end do
!$omp end parallel

      ! Check transient time scale
      !
      if (option%check_timescale) then
        write(iw,*)
        write(iw,'("Analyze> Check transient time-scale")')
        call check_transient(output, option, cvinfo, state)
        call termination("sr_analysis")
      end if

      ! Count # of binding/unbinding events
      ! (stored in state%nBevents and state%nUBevents)
      !
      write(iw,*)
      write(iw,'("Analyze> Analyze reaction events")')
      call get_reaction_info(option, cvinfo, state)

      ! Calculate transition time-correlation functions 
      !
      write(iw,*)
      write(iw,'("Analyze> Get transition TCF")')
      call get_transtcf(option, cvinfo, state)

      write(iw,*)
      write(iw,'("Analyze> Get kinetic properties")')
      call get_kineprop(option, output, cvinfo, bootopt, state)

      !deallocate(rp)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_state(option, cv, state)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),         intent(in)  :: option
      type(s_cv),             intent(in)  :: cv
      type(s_state),          intent(out) :: state 

      integer                :: istep, istate, icv, ia
      integer                :: nstate, ndim, nstep
      logical                :: is_assigned, is_quenched
      real(8)                :: wrk(option%ndim)
      character(len=MaxChar) :: fname


      nstate = option%nstate
      ndim   = option%ndim
      nstep  = cv%nstep

      state%quench_step = nstep

      ! for global
      allocate(state%data(nstep), state%data_ts(nstep), &
               state%avedata(nstep))

      is_quenched = .false.
      do istep = 1, cv%nstep
        wrk(:) = cv%data(:, istep) 
        is_assigned = .false.
        do istate = 1, nstate
          if (is_assigned) then
            exit
          else
            ia = 0
            do icv = 1, ndim 
              if (wrk(icv) >= option%state_def(1, icv, istate) &
                .and. wrk(icv) < option%state_def(2, icv, istate)) then
                ia = ia + 1
              end if
            end do

            if (ia == ndim) then 
              is_assigned       = .true.
              state%data(istep) = istate 
            end if

          end if
        end do

        if (.not. is_assigned) then
          state%data(istep) = 0 
        end if

        ! Check quenching
        !
        if (option%use_quench) then
          if (.not. is_quenched) then
            ia = 0
            do icv = 1, ndim
              if (wrk(icv) >= option%state_def(1, icv, nstate + 1) &
                .and. wrk(icv) < option%state_def(2, icv, nstate + 1)) then
                ia = ia + 1
              end if
            end do

            if (ia == ndim) then
              is_quenched       = .true.
              state%quench_step = istep
            end if

          end if
        end if


      end do

      state%data_ts = OFF 
      do istep = 1, nstep
        if (state%data(istep) == option%bound_id) then 
          state%data_ts(istep) = ON 
        end if
      end do
      
      state%nstep = nstep
      state%ndim  = ndim

    end subroutine get_state 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_kineprop(option, output, cvinfo, bootopt, state)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),  intent(in)    :: option
      type(s_output),  intent(in)    :: output
      type(s_cvinfo),  intent(in)    :: cvinfo
      type(s_bootopt), intent(inout) :: bootopt
      type(s_state),   intent(inout) :: state(cvinfo%nfile)

      ! Local variables
      !
      type(s_booteach) :: beach
      type(s_bootave)  :: bave

      integer :: nstate, nstep, ntrial, nsample, nt_range, nfile
      integer :: ifile, istep, itrial, is, iunit
      character(len=MaxChar) :: fname

      integer, allocatable :: nstay(:, :)


      ! Prepare parameters 
      !
      nstate   = option%nstate
      nstep    = state(1)%nstep
      ntrial   = bootopt%ntrial
      nsample  = bootopt%nsample
      nt_range = option%nt_range
      nfile    = cvinfo%nfile

      if (.not. option%use_bootstrap) then
        ntrial  = 1
        nsample = cvinfo%nfile 
      end if

      ! Allocate memory
      !
      allocate(nstay(0:nstate, 1:nfile))

      allocate(beach%rand(nsample, ntrial))
      allocate(beach%taud(ntrial), beach%kd(ntrial))
      allocate(beach%tauins(ntrial), beach%kins(ntrial))
      allocate(beach%hist(0:nt_range, 0:nstate, 0:nstate, ntrial))

      allocate(bave%hist      (0:nt_range, 0:nstate, 0:nstate))
      allocate(bave%hist_err  (0:nt_range, 0:nstate, 0:nstate))
      allocate(bave%hist_stdev(0:nt_range, 0:nstate, 0:nstate))
      allocate(bave%cumm      (0:nt_range))

      if (option%calc_2nd) then
        allocate(beach%hist2(0:nt_range, 0:nstate, ntrial))
        allocate(beach%pa3(0:nt_range, ntrial))
        allocate(beach%taud2(0:nstate, ntrial))
        allocate(beach%taupa3(ntrial))

        allocate(bave%hist2(0:nt_range, 0:nstate))
        allocate(bave%pa3(0:nt_range))
        allocate(bave%taud2(0:nstate))
        allocate(bave%taud2_err(0:nstate))
        allocate(bave%cumm2(0:nt_range, 0:nstate))
        allocate(bave%cummpa3(0:nt_range))
      end if

      ! Setup Bootstrap (if use_bootstrap = .true.) 
      !
      if (option%use_bootstrap) then

        write(iw,*)
        write(iw,'("Get_Kineprop> Setup Bootstrap")')

        ! Generate random seed
        !
        call get_seed(bootopt%iseed)
        call initialize_random(bootopt%iseed)

        ! Generate random numbers
        !
        do itrial = 1, ntrial
          call get_random_integer(nsample,              &
                                  1,                    &
                                  nfile,                &
                                  bootopt%duplicate,    &
                                  beach%rand(1, itrial))
        end do

      end if

      ! Calculate time spent in each region for each traj  
      !
      nstay = 0
!$omp parallel private(ifile, istep, is), shared(nstay), &
!$omp          default(shared)
!$omp do
      do ifile = 1, nfile
        do istep = 1, nstep
          is = state(ifile)%data(istep)
          nstay(is, ifile) = nstay(is, ifile) + 1
        end do
      end do
!$omp end do
!$omp end parallel

      ! Calculate kinetic properties
      !
      write(iw,*)
      write(iw,'("Get_Kineprop> Calculate average values")')
      call get_kinetics_average(option, cvinfo, bootopt, state, nstay, beach, bave)

      call output_kinetics(option, output, bave)

    end subroutine get_kineprop
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_kinetics_average(option, cvinfo, bootopt, state, nstay, beach, bave)
!-----------------------------------------------------------------------
      implicit none

      type(s_option),   intent(in)    :: option
      type(s_cvinfo),   intent(in)    :: cvinfo
      type(s_bootopt),  intent(in)    :: bootopt
      type(s_state),    intent(in)    :: state(cvinfo%nfile)
      integer,          intent(in)    :: nstay(0:option%nstate, cvinfo%nfile)
      type(s_booteach), intent(inout) :: beach
      type(s_bootave),  intent(inout) :: bave

      integer :: nstate, ntrial, nsample, nt_range
      integer :: itrial, isample, itraj, istep, is, js, rid
      integer :: nthreads, thread_id, nwork, ifinish
      integer :: nev, nuev
      real(8) :: t, p1, p2, f, rval
      real(8) :: dev, val, sterr
      real(8) :: progress

      real(8), allocatable :: norm_tot(:, :)


      ! Setup parameters
      !
      nstate   = option%nstate
      nt_range = option%nt_range
      ntrial   = bootopt%ntrial
      nsample  = bootopt%nsample
      rid      = option%reaczone_id

      if (.not. option%use_bootstrap) then
        ntrial  = 1
        nsample = cvinfo%nfile 
      end if

      ! Memory allocation
      !
      allocate(norm_tot(0:nt_range, 0:nstate))
      norm_tot = 0.0d0


      beach%hist = 0.0d0
      beach%taud = 0.0d0
      beach%kd   = 0.0d0
      beach%kins = 0.0d0

      beach%hist2 = 0.0d0

!$omp parallel private(itrial, isample, itraj, istep, is, js, nev, nuev, norm_tot, t, &
!$omp                  nthreads, thread_id, nwork, progress, ifinish),   &
!$omp          shared(nstay), default(shared)
      nthreads  = omp_get_num_threads()
      thread_id = omp_get_thread_num()
      nwork     = ntrial / nthreads + 1
      ifinish   = 0
!$omp do 
      do itrial = 1, ntrial
        
        if (option%use_bootstrap .and. thread_id == 0) then
          ifinish  = ifinish + 1
          progress = ifinish / dble(nwork) * 100.0d0
          write(iw,'("  Progress : ", f6.2, "%")') progress
        end if

        nev      = 0
        nuev     = 0
        norm_tot = 0.0d0

        do isample = 1, nsample

          if (option%use_bootstrap) then
            itraj = beach%rand(isample, itrial)
          else
            itraj = isample
          end if

          do is = 0, nstate

            do istep = 0, nt_range
              norm_tot(istep, is) &
                = norm_tot(istep, is) + state(itraj)%norm(istep, is)
            end do

            do js = 0, nstate
              do istep = 0, nt_range
                beach%hist(istep, js, is, itrial)       &
                  = beach%hist(istep, js, is, itrial)   &
                  + state(itraj)%hist(istep, js, is)
              end do
            end do

          end do

          ! Integrate TCFs
          !
          call trape_integral(0,                               &
                              nt_range,                        &
                              option%dt,                       &
                              beach%hist(:, rid, rid, itrial), &
                              beach%taud(itrial))

          nev  = nev  + state(itraj)%nBevents
          nuev = nuev + nstay(option%reaczone_id, itraj)

          ! 2nd-order transition
          !

          if (option%calc_2nd) then
            do is = 0, nstate
              do istep = 0, nt_range
                beach%hist2(istep, is, itrial)     &
                  = beach%hist2(istep, is, itrial) &
                  + state(itraj)%hist2(istep, is)
              end do
            end do
          end if

        end do

        do is = 0, nstate
          if (norm_tot(0, is) /= 0.0d0) then
            do js = 0, nstate
              do istep = 0, nt_range
                beach%hist(istep, js, is, itrial)      &
                  = beach%hist(istep, js, is, itrial)  &
                    / norm_tot(istep, is)
              end do
            end do
          end if
        end do

        beach%kins(itrial)   = dble(nev) / (dble(nuev) * option%dt)
        beach%tauins(itrial) = 1.0d0 / beach%kins(itrial)


        if (option%calc_2nd) then

          ! P_A<=i<=A (t)
          !
          do is = 0, nstate
            do istep = 0, nt_range
              beach%hist2(istep, is, itrial) &
                = beach%hist2(istep, is, itrial) / norm_tot(istep, option%reaczone_id)
            end do
          end do

          ! P_A<=A<=A (t)
          !
          beach%pa3(:, itrial) &
            = beach%hist2(:, option%reaczone_id, itrial)

          ! P_A<=i<=A (t) / t
          !
          beach%hist2(0, option%reaczone_id, itrial) &
            = beach%hist2(0, option%reaczone_id, itrial) / option%dt

          do is = 0, nstate
            do istep = 1, nt_range
              t = option%dt * istep
              beach%hist2(istep, is, itrial) = beach%hist2(istep, is, itrial) / t
            end do
          end do

        end if

      end do
!$omp end do
!$omp end parallel

      ! Integrate TCFs for each trial
      !
      do itrial = 1, ntrial
        call trape_integral(0, nt_range, option%dt, beach%hist(:, rid, rid, itrial), beach%taud(itrial))
        beach%kd(itrial) = 1.0d0 / beach%taud(itrial)
      end do

      if (option%calc_2nd) then

        do itrial = 1, ntrial

          ! P_A<i<A (t) / t
          do is = 1, nstate
            call trape_integral(0, nt_range, option%dt, beach%hist2(:, is, itrial), beach%taud2(is, itrial))
          end do

          ! P_A<A<A (t)
          call trape_integral(0, nt_range, option%dt, beach%pa3(:, itrial), beach%taupa3(itrial))

        end do

      end if

      ! MC boostrap average 
      !
      if (option%use_bootstrap) then

        bave%hist   = 0.0d0
        bave%taud   = 0.0d0
        bave%kd     = 0.0d0
        bave%tauins = 0.0d0
        bave%kins   = 0.0d0

        if (option%calc_2nd) then
          bave%hist2     = 0.0d0
          bave%taud2     = 0.0d0
          bave%taud2_err = 0.0d0
          bave%pa3       = 0.0d0
          bave%taupa3    = 0.0d0
        end if

        do itrial = 1, ntrial
          bave%hist(:, :, :) = bave%hist(:, :, :) + beach%hist(:, :, :, itrial)
          bave%taud          = bave%taud          + beach%taud(itrial)
          !bave%tauins        = bave%tauins        + beach%tauins(itrial)
          bave%kins          = bave%kins          + beach%kins(itrial)

          if (option%calc_2nd) then
            bave%hist2(:, :) = bave%hist2(:, :)   + beach%hist2(:, :, itrial)
            bave%taud2(:)    = bave%taud2(:)      + beach%taud2(:, itrial)
            bave%pa3(:)      = bave%pa3(:)        + beach%pa3(:, itrial)
            bave%taupa3      = bave%taupa3        + beach%taupa3(itrial)
          end if

        end do

        bave%hist(:, :, :) = bave%hist(:, :, :) / ntrial
        bave%taud          = bave%taud          / ntrial
        bave%kins          = bave%kins          / ntrial

        bave%kd            = 1.0d0 / bave%taud
        bave%tauins        = 1.0d0 / bave%kins

        if (option%calc_2nd) then
          bave%hist2       = bave%hist2  / ntrial
          bave%taud2       = bave%taud2  / ntrial
          bave%pa3         = bave%pa3    / ntrial
          bave%taupa3      = bave%taupa3 / ntrial
        end if

        ! Running integral of averaged TCFs

        ! P_A<A (t)
        call trape_integral(0, nt_range, option%dt, bave%hist(:, rid, rid), rval, cumm = bave%cumm)

        if (option%calc_2nd) then

          ! P_A<i<A (t) / t
          do is = 0, nstate 
            call trape_integral(0, nt_range, option%dt, bave%hist2(:, is), rval, cumm = bave%cumm2(:, is))
          end do

          ! P_A<A<A (t)
          call trape_integral(0, nt_range, option%dt, bave%pa3, rval, cumm = bave%cummpa3)

        end if

        ! Error analysis
        !   Note: StdDev of Bootstrap is StdErr
        ! - taud
        !
        call standard_deviation(ntrial, bave%taud, beach%taud, bave%taud_stdev)
        bave%taud_err = bave%taud_stdev 

        ! - kd
        !
        bave%kd_stdev = bave%taud_stdev / bave%taud**2 
        bave%kd_err   = bave%kd_stdev

        ! - kins
        !
        call standard_deviation(ntrial, bave%kins, beach%kins, bave%kins_stdev)
        bave%kins_err = bave%kins_stdev

        ! - tauins
        !
        bave%tauins_stdev = bave%kins_stdev / bave%kins**2 
        bave%tauins_err   = bave%tauins_stdev

        ! - hist
        !
        do is = 0, nstate
          do js = 0, nstate
            do istep = 0, nt_range
              call standard_deviation(ntrial,                         &
                                      bave%hist(istep, js, is),       &
                                      beach%hist(istep, js, is, :),   &
                                      bave%hist_stdev(istep, js, is))
              bave%hist_err = bave%hist_stdev
            end do
          end do
        end do

        if (option%calc_2nd) then

          ! - taud2
          !
          do is = 1, nstate
            call standard_deviation(ntrial,             &
                                    bave%taud2(is),     &
                                    beach%taud2(is, :), &
                                    bave%taud2_err(is))
          end do

          ! - taupa3
          !
          call standard_deviation(ntrial,         &
                                  bave%taupa3,    &
                                  beach%taupa3,   &
                                  bave%taupa3_err)

        end if

      else
        bave%taud         = beach%taud(1)
        bave%taud_stdev   = 0.0d0

        bave%kd           = beach%kd(1)
        bave%taud_stdev   = 0.0d0

        bave%tauins       = beach%tauins(1)
        bave%tauins_stdev = 0.0d0

        bave%kins         = beach%kins(1)
        bave%kins_stdev   = 0.0d0

        bave%hist         = beach%hist(:, :, :, 1)
        bave%hist_stdev   = 0.0d0

        if (option%calc_2nd) then
          bave%hist2   = beach%hist2(:, :, 1)
          bave%pa3     = beach%pa3(:, 1)
          bave%taud2   = beach%taud2(:, 1)
          bave%taupa3  = beach%taupa3(1)

          bave%taud_err   = 0.0d0
          bave%taupa3_err = 0.0d0

          bave%cumm2   = 0.0d0
          bave%cummpa3 = 0.0d0
        end if

        bave%cumm    = 0.0d0
        bave%cumm(0) = 0.5d0 * option%dt &
                     * (bave%hist(0, rid, rid) + bave%hist(1, rid, rid))
        do istep = 1, nt_range - 1
          bave%cumm(istep) = bave%cumm(istep - 1) &
            + 0.5d0 * option%dt &
              * (bave%hist(istep, rid, rid) + bave%hist(istep + 1, rid, rid))
        end do

        bave%cumm(nt_range) = bave%cumm(nt_range - 1) &
          + option%dt * bave%hist(nt_range, rid, rid)

      end if


    end subroutine get_kinetics_average
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine output_kinetics(option, output, bave)
!-----------------------------------------------------------------------
      implicit none 

      type(s_option),  intent(in) :: option
      type(s_output),  intent(in) :: output
      type(s_bootave), intent(in) :: bave

      integer                :: iunit, istep, is, js, k
      real(8)                :: t
      character(len=MaxChar) :: fname
    

      ! PROP
      !
      write(fname,'(a,".prop")') trim(output%fhead)
      call open_file(fname, iunit)
      
      write(iunit,'("property      value          err")')
      write(iunit,'("taud    ", 2f15.7)') bave%taud,   bave%taud_stdev
      write(iunit,'("tauins  ", 2f15.7)') bave%tauins, bave%tauins_stdev
      write(iunit,'("kd      ", 2f15.7)') bave%kd,     bave%kd_stdev
      write(iunit,'("kins    ", 2f15.7)') bave%kins,   bave%kins_stdev

      if (option%calc_2nd) then
        write(iunit,*)
        do is = 1, option%nstate
          if (is /= option%bound_id) then
            write(iunit,'("taud2 ",i0,1x,2f15.7)') is, bave%taud2(is), bave%taud2_err(is)
          end if
        end do

        write(iunit,'("taupa3  ", 2f15.7)') bave%taupa3, bave%taupa3_err
      end if

      close(iunit)

      ! Pret(t)
      !
      write(fname,'(a,".pret")') trim(output%fhead)
      call open_file(fname, iunit)
      is = option%reaczone_id
      do istep = 0, option%nt_range
        t = istep * option%dt 
        write(iunit,'(f15.7)', advance='no') t
        write(iunit,'(f15.7)', advance='no') bave%hist(istep, is, is)
        write(iunit,'(f15.7)', advance='no') bave%hist_err(istep, is, is)
        write(iunit,'(f15.7)', advance='no') bave%cumm(istep)

        if (option%calc_2nd) then
          do js = 1, option%nstate
            if (js /= option%bound_id) then
              write(iunit,'(f15.7)', advance='no') bave%hist2(istep, js)
              write(iunit,'(f15.7)', advance='no') bave%cumm2(istep, js)
            end if
          end do
        end if
        write(iunit,*)
      end do
      close(iunit)

      ! Transition probability matrix
      !
      write(fname,'(a,".tpm_to_reac")') trim(output%fhead)
      call open_file(fname, iunit)
      is = option%reaczone_id
      write(iunit,'("# col. ",i5," : time")') 1
      k = 1 
      do js = 0, option%nstate
        if (js /= is) then
          k = k + 1
          write(iunit,'("# col. ",i5," : state ",i5, " to ", i5)') k, js, is
          if (js == 0) then
            write(iunit,'("#   Note: state 0 is a region which is not defined in fstate file")')
          end if
        end if
      end do

      do istep = 0, option%nt_range
        write(iunit,'(f15.7)', advance='no') option%dt * istep
        do js = 0, option%nstate
          if (js /= is) then 
            write(iunit,'(f15.7)', advance='no') bave%hist(istep, is, js)
          end if
        end do
        write(iunit,*)
      end do

      close(iunit)

      write(fname,'(a,".tpm_from_reac")') trim(output%fhead)
      call open_file(fname, iunit)
      is = option%reaczone_id
      write(iunit,'("# col. ",i5," : time")') 1
      k = 1 
      do js = 0, option%nstate
        if (js /= is) then
          k = k + 1
          write(iunit,'("# col. ",i5," : state ",i5, " to ", i5)') k, is, js
          if (js == 0) then
            write(iunit,'("#   Note: state 0 is a region which is not defined in fstate file")')
          end if
        end if
      end do

      do istep = 0, option%nt_range
        write(iunit,'(f15.7)', advance='no') option%dt * istep
        do js = 0, option%nstate
          if (js /= is) then 
            write(iunit,'(f15.7)', advance='no') bave%hist(istep, js, is)
          end if
        end do
        write(iunit,*)
      end do

      close(iunit)

      ! 2nd-order transition
      !
      if (option%calc_2nd) then
        write(fname,'(a,".pret2nd")') trim(output%fhead)
        call open_file(fname, iunit)
        do istep = 0, option%nt_range
          write(iunit,'(f15.7)', advance='no') option%dt * istep
          write(iunit,'(f15.7)', advance='no') bave%pa3(istep)
          write(iunit,'(f15.7)', advance='no') bave%cummpa3(istep)
          write(iunit,*)
        end do
        close(iunit)
      end if


      ! Generate input for RP rate 
      !
      write(fname,'(a,".rprate_input")') trim(output%fhead)
      call open_file(fname, iunit)
    
      write(iunit,'("&option_param")')
      k = 1 
      if (option%calc_2nd) then
        k = 2 
      end if

      write(iunit,'("  rporder     = ", i0)') k
      write(iunit,'("  timeunit    = 1.0d-9 ! please modify it by yourself")')
      write(iunit,'("  kins        = ", f15.7)') bave%kins
      write(iunit,'("  kins_err    = ", f15.7)') bave%kins_err
      write(iunit,'("  taud        = ", f15.7)') bave%taud
      write(iunit,'("  taud_err    = ", f15.7)') bave%taud_err
      write(iunit,'("  Kstar       = VALUE1 ! please modify it by yourself")')
      write(iunit,'("  Kstar_err   = VALUE2 ! please modify it by yourself")')

      if (option%calc_2nd) then
        write(iunit,'("  taupa3      = ", f15.7)') bave%taupa3
        write(iunit,'("  taupa3_err  = ", f15.7)') bave%taupa3_err
      end if

      write(iunit,'("/")')

      close(iunit)


    end subroutine output_kinetics
!-----------------------------------------------------------------------


end module mod_analyze
!=======================================================================
