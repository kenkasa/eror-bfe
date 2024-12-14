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
      real(8), allocatable   :: ngp(:, :), ngpave(:)
      real(8), allocatable   :: msd_stdev(:)
      real(8), allocatable   :: vhf(:, :)


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
        allocate(msd(0:option%nt_range, nmol))
        allocate(msdave(0:option%nt_range))
        allocate(msd_stdev(0:option%nt_range))
        allocate(ngp(0:option%nt_range, nmol))
        allocate(ngpave(0:option%nt_range))

        !allocate(msd(0:option%msdrange, nmol))
        !allocate(msdave(0:option%msdrange))
        !allocate(msd_stdev(0:option%msdrange))
        !allocate(ngp(0:option%msdrange, nmol))
        !allocate(ngpave(0:option%msdrange))

        write(iw,*)
        write(iw,'("Analyze> Calculate MSD")')

        if (option%use_cond) then
          write(iw,'("use_cond = .true. >> Calculate conditional MSD")')
          call calc_msd_cond(com_all, input, option, msd, msdave, msd_stdev)
        else
          call calc_msd(com_all, option, msd, msdave, msd_stdev, ngp, ngpave)
        end if

        write(iw,'("Finished")')
        write(iw,*)
      end if

      if (option%vhfcalc) then
        allocate(vhf(0:option%nr, 0:option%nt_range))
        call calc_vhf(com_all, option, vhf) 
      end if

      call generate_comfile(output,   &
                            option,   &
                            traj,     &
                            com_all)

      call generate_msdfile(output,   &
                            option,   &
                            traj,     &
                            com_all,  &
                            msd,      &
                            msdave,   &
                            msd_stdev)

      call generate_ngpfile(output,   &
                            option,   &
                            traj,     &
                            ngpave)

      call generate_vhffile(output,   &
                            option,   &
                            vhf)

      ! Deallocate memory
      !
      if (allocated(msd)) &
        deallocate(msd, msdave, ngp, ngpave)

      if (allocated(vhf)) &
        deallocate(vhf)


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine calc_msd(com, option, msd, msdave, msd_stdev, ngp, ngpave)
!-----------------------------------------------------------------------
      implicit none

      type(s_com),    intent(in)    :: com
      type(s_option), intent(in)    :: option
      real(8),        intent(inout) :: msd(0:option%nt_range, com%nmol)
      real(8),        intent(inout) :: msdave(0:option%nt_range)
      real(8),        intent(inout) :: msd_stdev(0:option%nt_range)
      real(8),        intent(inout) :: ngp(0:option%nt_range, com%nmol)
      real(8),        intent(inout) :: ngpave(0:option%nt_range)

      integer :: imol, istep, jstep, ij, nend
      integer :: nmol, nstep, nt_range, msddim 
      real(8) :: d(3), d2, d4, dev

      real(8), allocatable :: crd(:, :, :)


      ! setup variables
      !
      nmol     = com%nmol
      nstep    = com%nstep
      nt_range = option%nt_range
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
      ngp = 0.0d0

!$omp parallel private(imol, istep, jstep, nend, ij, d, d2, d4), &
!$omp        & default(shared)
!$omp do
      do imol = 1, nmol
        do istep = 1, nstep - 1

          !nend = istep + msdrange

          !if (nend > nstep) &
          !  nend = nstep

          nend = min(istep + nt_range * option%nt_sparse, &
                     nstep)
          

          ij = -1
          do jstep = istep, nend, option%nt_sparse
            !ij            = jstep - istep
            ij            = ij + 1
            d(1:3)        = crd(1:3, jstep, imol) - crd(1:3, istep, imol)
            d2            = dot_product(d(1:msddim), d(1:msddim))
            d4            = d2 * d2

            msd(ij, imol) = msd(ij, imol) + d2
            ngp(ij, imol) = ngp(ij, imol) + d4 

          end do
        end do

        ! average over time steps for each molecule
        !
        do istep = 0, option%nt_range !msdrange - 1
          jstep = istep * option%nt_sparse

          !msd(istep, imol) = msd(istep, imol) / (nstep - istep)
          !ngp(istep, imol) = ngp(istep, imol) / (nstep - istep)
          msd(istep, imol) = msd(istep, imol) / (nstep - jstep)
          ngp(istep, imol) = ngp(istep, imol) / (nstep - jstep)
        end do

        do istep = 1, option%nt_range
          ngp(istep, imol) = 3.0d0*ngp(istep, imol) / (5.0d0*msd(istep, imol)**2.0) - 1.0
        end do

      end do
!$omp end do
!$omp end parallel

      ! average over molecules
      !
      do istep = 0, option%nt_range !msdrange - 1
        msdave(istep) = sum(msd(istep, 1:nmol)) / dble(nmol)
        ngpave(istep) = sum(ngp(istep, 1:nmol)) / dble(nmol)
      end do

      do istep = 0, option%nt_range !msdrange - 1
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
      real(8),        intent(inout) :: msd(0:option%nt_range, com%nmol)
      real(8),        intent(inout) :: msdave(0:option%nt_range)
      real(8),        intent(inout) :: msd_stdev(0:option%nt_range)

      integer :: imol, istep, jstep, ij, nend
      integer :: nmol, nstep, nt_range, msddim
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
      nt_range = option%nt_range
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

!$omp parallel private(imol, istep, jstep, val, nend, ij, d, d2), &
!$omp        & default(shared)
!$omp do
      do imol = 1, nmol
        do istep = 1, nstep - 1

          val  = ts(istep, imol) 
          !nend = istep + msdrange

          if (val < option%rcrange(1) .or. val > option%rcrange(2)) &
            cycle

          !if (nend > nstep) &
          !  nend = nstep
          nend = min(istep + nt_range * option%nt_sparse, &
                     nstep)


          ij = - 1
          do jstep = istep, nend, option%nt_sparse
            ij            = ij + 1
            !ij            = jstep - istep
            d(1:3)        = crd(1:3, jstep, imol) - crd(1:3, istep, imol)
            d2            = dot_product(d(1:msddim), d(1:msddim))

            msd(ij, imol) = msd(ij, imol) + d2
            nij(ij, imol) = nij(ij, imol) + 1

          end do
        end do

        ! average over time steps for each molecule
        !
        if (nij(0, imol) == 0) then
          msd(:, imol) = 0.0d0
        else
          do istep = 0, option%nt_range !msdrange - 1
!            msd(istep, imol) = msd(istep, imol) / (nstep - istep)
            msd(istep, imol) = msd(istep, imol) / dble(nij(istep, imol))
          end do
        end if

      end do
!$omp end do
!$omp end parallel

      ! average over molecules
      !
      do istep = 0, option%nt_range !msdrange - 1
!        msdave(istep) = sum(msd(istep, 1:nmol)) / dble(nmol)
        msdave(istep) = sum(msd(istep, 1:nmol) * nij(istep, 1:nmol)) &
                      / dble(sum(nij(istep, 1:nmol)))
      end do

      do istep = 0, option%nt_range !msdrange - 1
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
    subroutine calc_vhf(com, option, vhf)
!-----------------------------------------------------------------------
      implicit none

      type(s_com),    intent(in)    :: com
      type(s_option), intent(in)    :: option
      real(8),        intent(inout) :: vhf(0:option%nr, 0:option%nt_range)

      integer :: imol, istep, jstep, it, ij, nend
      integer :: nmol, nstep
      integer :: rind
      real(8) :: d(3), d2, dev
      real(8) :: dr, vsum
      integer :: nt_range, nt_sparse, nr

      real(8), allocatable :: norm(:)
      real(8), allocatable :: crd(:, :, :)


      ! Setup variables
      !
      nmol      = com%nmol
      nstep     = com%nstep
      nt_range  = option%nt_range
      nt_sparse = option%nt_sparse
      nr        = option%nr
      dr        = option%dr

      ! Allocate memory 
      !
      allocate(crd(3, nstep, nmol))
      allocate(norm(0:nt_range))

      ! Prepare modified array of CoM 
      !
      do istep = 1, nstep
        do imol = 1, nmol
          crd(1:3, istep, imol) = com%coord(1:3, imol, istep)
        end do
      end do

      ! Calculate MSD
      !
      vhf  = 0.0d0
      norm = 0.0d0

!$omp parallel private(imol, istep, jstep, nend, it, rind, d, d2), &
!$omp        & default(shared), &
!$omp        & reduction(+:vhf), reduction(+:norm) 
!$omp do
      do imol = 1, nmol
        do istep = 1, nstep - 1

          nend = min(istep + nt_range * nt_sparse, nstep)

          it = -1
          do jstep = istep, nend, option%nt_sparse
            it = it + 1
            d(1:3)        = crd(1:3, jstep, imol) - crd(1:3, istep, imol)
            d2            = dot_product(d(1:3), d(1:3))

            rind          = nint(sqrt(d2) / dr) 
            vhf(rind, it) = vhf(rind, it) + 1.0d0 
            norm(it)      = norm(it)      + 1.0d0
          end do
        end do

      end do
!$omp end do
!$omp end parallel

      ! Normalize 
      !
      do istep = 0, nt_range 
        vhf(:, istep) = vhf(:, istep) / (norm(istep) * dr)
      end do
      vhf(0,:) = vhf(0, :) * 2.0d0 

      ! Test 
      !
      !do istep = 0, nt_range
      !  vsum = vhf(0, istep) * 0.5d0 * dr
      !  do rind = 1, nr
      !    vsum = vsum + vhf(rind, istep) * dr 
      !  end do
      !  write(iw,'("Integral of vhf at step ", i0, " = ", f20.10)') istep, vsum
      !end do

      ! deallocate memory
      !
      deallocate(crd, norm)

    end subroutine calc_vhf
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


      if (.not. option%comcalc) &
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
      real(8),        intent(in) :: msd(0:option%nt_range, com%nmol)
      real(8),        intent(in) :: msdave(0:option%nt_range)
      real(8),        intent(in) :: msd_stdev(0:option%nt_range)

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
      do istep = 0, option%nt_range !option%msdrange - 1
        write(UnitOUT,'(f20.10)', advance='no') option%dt_out * istep !traj%dt * istep
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
      do istep = 0, option%nt_range !option%msdrange - 1
        write(UnitOUT,'(3f20.10)') option%dt_out * istep, & !traj%dt * istep,  &
                                   msdave(istep),    &
                                   msd_stdev(istep)
      end do
      close(UnitOUT)


    end subroutine generate_msdfile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine generate_ngpfile(output,    &
                                option,    &
                                traj,      &
                                ngpave)
!-----------------------------------------------------------------------
      implicit none

      ! formal arguments
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      type(s_traj),   intent(in) :: traj
      real(8),        intent(in) :: ngpave(0:option%nt_range)

      ! local variables
      integer                :: istep
      character(len=MaxChar) :: fname


      if (.not. option%ngpcalc) &
        return

      write(fname,'(a,".ngpave")') trim(output%fhead)

      open(UnitOUT, file=trim(fname))
      do istep = 1, option%nt_range !option%msdrange - 1
        write(UnitOUT,'(2f20.10)') option%dt_out * istep, & !traj%dt * istep,  &
                                   ngpave(istep)
      end do
      close(UnitOUT)


    end subroutine generate_ngpfile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine generate_vhffile(output,    &
                                option,    &
                                vhf)
!-----------------------------------------------------------------------
      implicit none

      ! formal arguments
      type(s_output), intent(in) :: output
      type(s_option), intent(in) :: option
      real(8),        intent(in) :: vhf(0:option%nr, 0:option%nt_range)

      ! local variables
      integer                :: istep, rind
      real(8)                :: r
      character(len=MaxChar) :: fname


      if (.not. option%vhfcalc) &
        return

      ! VHF 
      !
      write(fname,'(a, ".vhf")') trim(output%fhead)

      open(UnitOUT, file=trim(fname))
      do rind = 0, option%nr
        r = option%dr * dble(rind)
        write(UnitOUT,'(f20.10)', advance='no') r
        do istep = 0, option%nt_range
          write(UnitOUT,'(f20.10)', advance='no') vhf(rind, istep)
        end do
        write(UnitOUT,*)
      end do
      close(UnitOUT)

    end subroutine generate_vhffile
!-----------------------------------------------------------------------


end module mod_analyze
!=======================================================================
