!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_ctrl
  use mod_xtcio
  use mod_dcdio
  use mod_traj
  use mod_com
  use xdr, only: xtcfile

  ! subroutines
  !
  public  :: analyze
  private :: calc_msd 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      ! formal arguments
      type(s_input),  intent(in)    :: input 
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj

      ! local variables
      type(s_dcd)            :: dcd
      type(xtcfile)          :: xtc
      type(s_com)            :: com, com_all

      integer                :: iunit, itraj
      integer                :: istep, istep_tot
      integer                :: nmol, natm, nstep_tot
      integer                :: trajtype
      logical                :: is_end

      real(8), allocatable   :: msd(:, :), msdave(:)
      real(8), allocatable   :: msd_stdev(:)


      ! Get trajectory type
      !
      call get_trajtype(input%ftraj(1), trajtype)

      ! Get Total number of steps
      !
      if (trajtype == TrajTypeDCD) then
        call get_total_step_from_dcd(input%ftraj, nstep_tot)
      else if (trajtype == TrajTypeXTC) then
        call get_total_step_from_xtc(input%ftraj, nstep_tot)
      end if

      call get_com(option%mode,          &
                   traj,                 &
                   com,                  &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 0)

      call get_com(option%mode,          &
                   traj,                 &
                   com_all,              &
                   setup = .true.,       &
                   calc_coord = .false., &
                   myrank = 1)

      nmol          = com_all%nmol
      com_all%nstep = nstep_tot

      allocate(com_all%coord(1:3, nmol, nstep_tot))

      ! Get CoM
      !
      write(iw,*)
      write(iw,'("Analyze> Get CoM coordinates")')
      istep_tot = 0
      do itraj = 1, input%ntraj

        call open_trajfile(input%ftraj(itraj), trajtype, iunit, dcd, xtc)
        call init_trajfile(trajtype, iunit, dcd, xtc, natm)

        is_end = .false.
        istep  = 0
        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          call read_trajfile_oneframe(trajtype, iunit, istep, dcd, xtc, is_end)

          if (is_end) exit

          if (mod(istep_tot, 100) == 0) then
            write(iw,'("Step ", i0)') istep_tot
          end if

          call send_coord_to_traj(1, trajtype, dcd, xtc, traj)


          call get_com(option%mode,         &
                       traj,                &
                       com,                 &
                       setup = .false.,     &
                       calc_coord = .true., & 
                       myrank = 1)

          com_all%coord(1:3, 1:nmol, istep_tot) &
            = com%coord(1:3, 1:nmol, 1)

        end do

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, iunit, dcd, xtc)

      end do

      ! Calculate MSD
      !
      if (option%msdcalc) then
        allocate(msd(0:option%msdrange, nmol))
        allocate(msdave(0:option%msdrange))
        allocate(msd_stdev(0:option%msdrange))

       
        write(iw,*)
        write(iw,'("Analyze> Calculate MSD")')

        if (option%use_cond) then
          write(iw,'("use_cond = .true. >> Calculate conditional MSD")')
          call calc_msd_cond(com_all, input, option, msd, msdave, msd_stdev)
        else
          call calc_msd(com_all, option, msd, msdave, msd_stdev)
        end if

        write(iw,'("Finished")')
        write(iw,*)
      end if

      call generate_comfile(output, option, traj, com_all)

      call generate_msdfile(output,   &
                            option,   &
                            traj,     &
                            com_all,  &
                            msd,      &
                            msdave,   &
                            msd_stdev) 

      ! Deallocate memory
      !
      if (allocated(msd)) &
        deallocate(msd, msdave)


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_msd(com, option, msd, msdave, msd_stdev)
!-----------------------------------------------------------------------
      implicit none

      type(s_com),    intent(in)    :: com
      type(s_option), intent(in)    :: option
      real(8),        intent(inout) :: msd(0:option%msdrange, com%nmol)
      real(8),        intent(inout) :: msdave(0:option%msdrange)
      real(8),        intent(inout) :: msd_stdev(0:option%msdrange)

      integer :: imol, istep, jstep, ij, nend
      integer :: nmol, nstep, msdrange, msddim 
      real(8) :: d(3), d2, dev

      real(8), allocatable :: crd(:, :, :)


      ! setup variables
      !
      nmol     = com%nmol
      nstep    = com%nstep
      msdrange = option%msdrange
      msddim   = option%msddim

      ! allocate memory 
      !
      allocate(crd(3, nstep, nmol))

      ! prepare modified array of CoM 
      !
      do istep = 1, nstep
        do imol = 1, nmol
          crd(1:3, istep, imol) = com%coord(1:3, imol, istep)
        end do
      end do

      ! calculate MSD
      !
      msd = 0.0d0

!$omp parallel private(imol, istep, nend, ij, d, d2), &
!$omp        & default(shared)
!$omp do
      do imol = 1, nmol
        do istep = 1, nstep - 1

          nend = istep + msdrange

          if (nend > nstep) &
            nend = nstep

          do jstep = istep, nend
            ij            = jstep - istep
            d(1:3)        = crd(1:3, jstep, imol) - crd(1:3, istep, imol)
            d2            = dot_product(d(1:msddim), d(1:msddim))

            msd(ij, imol) = msd(ij, imol) + d2 

          end do
        end do

        ! average over time steps for each molecule
        !
        do istep = 0, msdrange - 1
          msd(istep, imol) = msd(istep, imol) / (nstep - istep)
        end do

      end do
!$omp end do
!$omp end parallel

      ! average over molecules
      !
      do istep = 0, msdrange - 1
        msdave(istep) = sum(msd(istep, 1:nmol)) / dble(nmol)
      end do

      do istep = 0, msdrange - 1
        dev = 0.0d0
        do imol = 1, nmol
          dev = dev + (msd(istep, imol) - msdave(istep)) ** 2
        end do
        dev              = sqrt(dev / dble(nmol - 1))
        msd_stdev(istep) = dev
      end do

      ! deallocate memory
      !
      deallocate(crd)

    end subroutine calc_msd
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_msd_cond(com, input, option, msd, msdave, msd_stdev)
!-----------------------------------------------------------------------
      implicit none

      type(s_com),    intent(in)    :: com
      type(s_input),  intent(in)    :: input
      type(s_option), intent(in)    :: option
      real(8),        intent(inout) :: msd(0:option%msdrange, com%nmol)
      real(8),        intent(inout) :: msdave(0:option%msdrange)
      real(8),        intent(inout) :: msd_stdev(0:option%msdrange)

      integer :: imol, istep, jstep, ij, nend
      integer :: nmol, nstep, msdrange, msddim
      integer :: iunit
      character(len=MaxChar) :: dum
      real(8) :: val, d(3), d2, dev

      integer, allocatable :: nij(:, :)
      real(8), allocatable :: crd(:, :, :)
      real(8), allocatable :: ts(:, :) 


      ! setup variables
      !
      nmol     = com%nmol
      nstep    = com%nstep
      msdrange = option%msdrange
      msddim   = option%msddim

      ! allocate memory 
      !
      allocate(crd(3, nstep, nmol))
      allocate(ts(nstep, nmol))
      allocate(nij(0:nstep, nmol))

      ! prepare modified array of CoM 
      !
      do istep = 1, nstep
        do imol = 1, nmol
          crd(1:3, istep, imol) = com%coord(1:3, imol, istep)
        end do
      end do

      ! load time-series data
      !
      call open_file(trim(input%fts), iunit)

      do istep = 1, nstep
        read(iunit,*) dum, (ts(istep, imol), imol = 1, nmol) 
      end do

      close(iunit)

      ! calculate MSD
      !
      msd = 0.0d0
      nij = 0

!$omp parallel private(imol, istep, val, nend, ij, d, d2), &
!$omp        & default(shared)
!$omp do
      do imol = 1, nmol
        do istep = 1, nstep - 1

          val  = ts(istep, imol) 
          nend = istep + msdrange

          if (val < option%rcrange(1) .or. val > option%rcrange(2)) &
            cycle

          if (nend > nstep) &
            nend = nstep

          do jstep = istep, nend
            ij            = jstep - istep
            d(1:3)        = crd(1:3, jstep, imol) - crd(1:3, istep, imol)
            d2            = dot_product(d(1:msddim), d(1:msddim))

            msd(ij, imol) = msd(ij, imol) + d2
            nij(ij, imol) = nij(ij, imol) + 1

          end do
        end do

        ! average over time steps for each molecule
        !
        if (nij(0, imol) == 0) then
          msd(istep, imol) = 0.0d0
        else
          do istep = 0, msdrange - 1
!            msd(istep, imol) = msd(istep, imol) / (nstep - istep)
            msd(istep, imol) = msd(istep, imol) / dble(nij(istep, imol))
          end do
        end if

      end do
!$omp end do
!$omp end parallel

      ! average over molecules
      !
      do istep = 0, msdrange - 1
!        msdave(istep) = sum(msd(istep, 1:nmol)) / dble(nmol)
        msdave(istep) = sum(msd(istep, 1:nmol) * nij(istep, 1:nmol)) &
                      / dble(sum(nij(istep, 1:nmol)))
      end do

      do istep = 0, msdrange - 1
        dev = 0.0d0
        do imol = 1, nmol
          dev = dev + (msd(istep, imol) - msdave(istep)) ** 2
        end do
        dev              = sqrt(dev / dble(nmol - 1))
        msd_stdev(istep) = dev
      end do

      ! deallocate memory
      !
      deallocate(crd)
      deallocate(ts)

    end subroutine calc_msd_cond
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine generate_comfile(output, option, traj, com)
!-----------------------------------------------------------------------
      implicit none

      ! formal arguments
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      type(s_traj),   intent(in) :: traj
      type(s_com),    intent(in) :: com

      ! local variables
      integer                :: istep, imol, ixyz, nmol 
      character(len=MaxChar) :: fname


      if (.not. option%outcom) &
        return

      nmol = com%nmol

      if (trim(output%fcom) /= "") then
        write(fname,'(a)') trim(output%fcom)
      else
        write(fname,'(a,".com")') trim(output%fhead)
      end if

      open(UnitOUT, file=trim(fname))

      if (option%onlyz) then

        do istep = 1, com%nstep
          write(UnitOUT,'(f20.10)',advance='no') traj%dt * istep
          do imol = 1, nmol
            write(UnitOUT,'(f20.10)',advance='no') com%coord(3, imol, istep)
          end do
          write(UnitOUT,*)
        end do

      else
        if (option%comformat == CoMFormatTypeXYZ) then

          do istep = 1, com%nstep
            write(UnitOUT,'(i0)') nmol 
            write(UnitOUT,*)
            do imol = 1, nmol 
              write(UnitOUT,'("Ar ",3f20.10)') &
                (com%coord(ixyz, imol, istep), ixyz = 1, 3)
            end do
          end do

        else if (option%comformat == CoMFormatTypeTIMESERIES) then

          do istep = 1, com%nstep
            write(UnitOUT,'(f20.10)',advance='no') traj%dt * istep
            do imol = 1, nmol
              write(UnitOUT,'(3f20.10)',advance='no') & 
                (com%coord(ixyz, imol, istep), ixyz = 1, 3)
            end do
            write(UnitOUT,*)
          end do

        end if 
      end if

      close(UnitOUT)

    end subroutine generate_comfile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine generate_msdfile(output,    &
                                option,    &
                                traj,      &
                                com,       &
                                msd,       &
                                msdave,    &
                                msd_stdev)
!-----------------------------------------------------------------------
      implicit none

      ! formal arguments
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      type(s_traj),   intent(in) :: traj
      type(s_com),    intent(in) :: com 
      real(8),        intent(in) :: msd(0:option%msdrange, com%nmol)
      real(8),        intent(in) :: msdave(0:option%msdrange)
      real(8),        intent(in) :: msd_stdev(0:option%msdrange)

      ! local variables
      integer                :: istep, imol, ixyz, nmol 
      character(len=MaxChar) :: fname


      nmol = com%nmol

      if (.not. option%msdcalc) &
        return

      ! MSD (each)
      !
      if (trim(output%fmsd) /= "") then
        write(fname,'(a)') trim(output%fmsd)
      else
        write(fname,'(a,".msd")') trim(output%fhead)
      end if

      open(UnitOUT, file=trim(fname))
      do istep = 0, option%msdrange - 1
        write(UnitOUT,'(f20.10)', advance='no') traj%dt * istep
        do imol = 1, nmol
          write(UnitOUT,'(f20.10)', advance='no') msd(istep, imol)
        end do
        write(UnitOUT,*) 
      end do
      close(UnitOUT)

      ! MSD (average)
      !
      if (trim(output%fmsdave) /= "") then
        write(fname,'(a)') trim(output%fmsdave)
      else
        write(fname,'(a,".msdave")') trim(output%fhead)
      end if

      open(UnitOUT, file=trim(fname))
      do istep = 0, option%msdrange - 1
        write(UnitOUT,'(3f20.10)') traj%dt * istep,  &
                                   msdave(istep),    &
                                   msd_stdev(istep)
      end do
      close(UnitOUT)


    end subroutine generate_msdfile
!-----------------------------------------------------------------------


end module mod_analyze
!=======================================================================
