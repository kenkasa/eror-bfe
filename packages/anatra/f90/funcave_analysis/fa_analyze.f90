!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl
  use mod_bootstrap
  use mod_random

  ! constants
  !

  ! structures
  !
  type :: s_booteach
    integer, allocatable :: rand(:, :)
    real(8), allocatable :: func(:, :)

  end type s_booteach

  type :: s_bootave
    real(8), allocatable :: ave(:)
    real(8), allocatable :: err(:) 
  end type s_bootave

  ! subroutines
  !
  public  :: analyze

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, cvinfo)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option
      type(s_cvinfo),  intent(in)    :: cvinfo

      type(s_cv),    allocatable :: cv(:)

      integer                :: nstep, nfile
      integer                :: ifile, istep
      real(8)                :: val, ave, dev, sterr
      real(8)                :: x
      character(len=MaxChar) :: fout

      real(8), allocatable :: norm(:)
      real(8), allocatable :: func_ave(:), func_stdev(:), func_sterr(:) 


      allocate(cv(cvinfo%nfile))

      ! read cv files
      !
      write(iw,*)
      write(iw,'("Analyze> Read CV file")')
      do ifile = 1, cvinfo%nfile
        call read_cv(cvinfo%fcv(ifile), 1, cv(ifile))
      end do

      nstep = cv(1)%nstep
      nfile = cvinfo%nfile

      allocate(func_ave(nstep), func_stdev(nstep), func_sterr(nstep))
      allocate(norm(nstep))

      func_ave   = 0.0d0
      func_stdev = 0.0d0
      func_sterr = 0.0d0
      norm       = 0.0d0
      
      do istep = 1, nstep
        ave   = 0.0d0
        dev   = 0.0d0
        sterr = 0.0d0

        do ifile = 1, nfile
          if (.not. isnan(cv(ifile)%data(1, istep))) then
            norm(istep) = norm(istep) + 1.0d0
            ave = ave + cv(ifile)%data(1, istep)  
          end if 
        end do
        !ave = ave / dble(nfile)
        ave = ave / norm(istep) 

        do ifile = 1, nfile
          if (.not. isnan(cv(ifile)%data(1, istep))) then
            dev = dev + (cv(ifile)%data(1, istep) - ave)**2
          end if 
        end do
        !dev   = sqrt(dev / (nfile - 1))
        !sterr = dev / sqrt(dble(nfile))
        dev   = sqrt(dev / (norm(istep) - 1.0d0))
        sterr = dev / sqrt(norm(istep))

        func_ave(istep)   = ave
        func_stdev(istep) = dev
        func_sterr(istep) = sterr
      end do

      write(fout, '(a,".ave")') trim(output%fhead)
      open(10,file=trim(fout))
      do istep = 1, nstep 
        x = option%xsta + option%dx * (istep - 1)
        write(10,'(e15.7,2x)',advance='no') x
        write(10,'(e15.7,2x)',advance='no') func_ave(istep)
        write(10,'(e15.7,2x)',advance='no') func_stdev(istep)
        write(10,'(e15.7,2x)')              func_sterr(istep)
      end do
      close(10)

      deallocate(func_ave, func_stdev, func_sterr)
      deallocate(norm)

    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine analyze_bootstrap(input, output, option, cvinfo, bootopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option
      type(s_cvinfo),  intent(in)    :: cvinfo
      type(s_bootopt), intent(inout) :: bootopt 

      type(s_booteach) :: beach
      type(s_bootave)  :: bave

      type(s_cv),    allocatable :: cv(:)

      integer                :: ntrial, nsample
      integer                :: nstep, nfile

      integer                :: ifile, isample, istep, ig, itrial
      real(8)                :: val, ave, dev, sterr
      real(8)                :: x
      character(len=MaxChar) :: fout

      real(8), allocatable :: func_ave(:), func_stdev(:), func_sterr(:) 


      ! Read CV files
      !
      allocate(cv(cvinfo%nfile))

      write(iw,*)
      write(iw,'("Analyze> Read CV file")')
      do ifile = 1, cvinfo%nfile
        call read_cv(cvinfo%fcv(ifile), 1, cv(ifile))
      end do

      nstep = cv(1)%nstep
      nfile = cvinfo%nfile

      ! Setup Bootstrap
      !
      ntrial  = bootopt%ntrial
      nsample = bootopt%nsample

      allocate(beach%func(nstep, ntrial))
      allocate(bave%ave(nstep), bave%err(nstep))

      ! Generate random numbers 
      !
      call get_seed(bootopt%iseed)
      call initialize_random(bootopt%iseed)

      allocate(beach%rand(nsample, ntrial))

      do itrial = 1, ntrial
        call get_random_integer(nsample, 1, nfile, bootopt%duplicate, beach%rand(1, itrial)) 
      end do

      ! Run Bootstrap
      !
      allocate(func_ave(nstep))

      beach%func = 0.0d0

      !$omp parallel private(itrial, isample, istep, val, func_ave) &
      !$omp          default(shared)
      !$omp do
      do itrial = 1, ntrial

        if (mod(itrial, 10) == 0) then
          write(iw,'("trial : ", i0)') itrial
        end if

        func_ave = 0.0d0

        do isample = 1, nsample
          ifile = beach%rand(isample, itrial)

          do istep = 1, nstep
            val             = cv(ifile)%data(1, istep) 
            func_ave(istep) = func_ave(istep) + val 
          end do

        end do

        func_ave                    = func_ave / dble(nsample)
        beach%func(1:nstep, itrial) = func_ave(1:nstep) 

      end do

      !$omp end do
      !$omp end parallel

      ! Calculate average & error
      !
      bave%ave = 0.0d0
      bave%err = 0.0d0

      do istep = 1, nstep
        bave%ave(istep) = sum(beach%func(istep, 1:ntrial)) / dble(ntrial)
      end do

      do istep = 1, nstep
        dev = 0.0d0
        ave = bave%ave(istep)
        do itrial = 1, ntrial
          val = beach%func(istep, itrial)
          dev = dev + (val - ave)**2 
        end do
        bave%err(istep) = sqrt(dev / dble(ntrial - 1))
      end do 

      ! Output
      !
      write(fout, '(a,".ave")') trim(output%fhead)
      open(10,file=trim(fout))
      do istep = 1, nstep 
        x = option%xsta + option%dx * (istep - 1)
        write(10,'(e15.7,2x)',advance='no') x
        write(10,'(e15.7,2x)',advance='no') bave%ave(istep)
        write(10,'(e15.7,2x)')              bave%err(istep)
      end do
      close(10)

      ! Deallocate
      !
      deallocate(func_ave)

    end subroutine analyze_bootstrap
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
