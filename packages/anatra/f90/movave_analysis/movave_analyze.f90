!=======================================================================
module mod_movave_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_movave_ctrl

  ! structures
  !
  type :: s_movave
    integer              :: ngrid
    real(8), allocatable :: grid(:)
    real(8), allocatable :: data(:)
    real(8), allocatable :: deriv(:)
  end type s_movave

  ! subroutines
  !
  public :: movave_analyze
  public :: movave_write

  contains
!-----------------------------------------------------------------------
    subroutine movave_analyze(input, output, option, movave, finp)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),                    intent(in)  :: input
      type(s_output),                   intent(in)  :: output
      type(s_movave_option),            intent(in)  :: option
      type(s_movave),                   intent(out) :: movave
      character(len=MaxChar), optional, intent(in)  :: finp

      integer                :: i, isep, istep, ngrid, nsplit(0:MaxSep)
      integer                :: reg_sta, reg_end
      integer                :: ista, iend, iorg
      real(8)                :: rdum, diff
      character(len=MaxChar) :: fread

      integer, allocatable   :: range(:, :)
      real(8), allocatable   :: input_data(:)


      ! Setup input file name
      !
      fread = trim(input%fts)
      if (present(finp)) &
        fread = finp
     
      ! Read input file 
      !

      !   Get # of grids
      !
      ngrid = 0
      open(10,file=trim(fread))
      do while (.true.) 
        read(10,*,end=100)
        ngrid = ngrid + 1 
      end do

100   continue

      !   Memory Alloc. for storing input data
      !
      movave%ngrid = ngrid
      allocate(input_data  (0:ngrid - 1), &
               movave%grid (0:ngrid - 1), &
               movave%data (0:ngrid - 1), &
               movave%deriv(0:ngrid - 1))

      rewind 10

      !   Read input data
      !
      input_data  = 0.0d0
      movave%grid = 0.0d0
      do istep = 0, ngrid - 1 
        read(10,*) movave%grid(istep), input_data(istep)
      end do

      close(10)


      ! Setup part 
      !

      !   Setup starting point for averaging
      !
      iorg = 1
      if (option%include_zero) then
        iorg = 0
      end if

      !   Setup boundary point
      !

      if (option%nregion == 1) then
        nsplit(1) = ngrid - 1
        if (option%include_zero) then
          nsplit(0) = iorg
        end if
      else
        do i = 1, option%nregion - 1
          nsplit(i) = nint((option%t_sep(i) - option%t_sta) / option%dt) 
        end do
        nsplit(0)              = iorg
        nsplit(option%nregion) = ngrid - 1 
      end if

      !   Setup average range for each point
      !
      allocate(range(1:3, 0:ngrid - 1))
      range = 0

      !   ... Initialize
      do istep = 0, ngrid - 1
        range(1:2, istep) = istep  ! 1: start , 2: end
        range(3,   istep) = 1      ! # of points used for averaging
      end do

      !   ... Determine the average range for each point
      do isep = 0, option%nregion - 1 
        reg_sta = nsplit(isep) 
        reg_end = nsplit(isep + 1) 

        do istep = reg_sta, reg_end 
          ista = istep - option%npoint(isep + 1) / 2
          iend = istep + option%npoint(isep + 1) / 2
       
          if (ista < 0) then
            ista = 0
          end if
       
          if (iend > ngrid - 1) then
            iend = ngrid - 1
          end if
       
          range(1, istep) = ista
          range(2, istep) = iend
          range(3, istep) = iend - ista + 1
       
        end do
      end do

      ! Perform Averaging
      !
      movave%data  = 0.0d0
      movave%deriv = 0.0d0

      if (.not.option%include_zero) then
        movave%data(0) = input_data(0)
      end if

      do istep = iorg, ngrid - 1
        ista = range(1, istep)
        iend = range(2, istep)
        movave%data(istep) &
          = sum(input_data(ista:iend)) / range(3, istep)
      end do

      do istep = 0, ngrid - 2
        diff = movave%data(istep + 1) - movave%data(istep)
        movave%deriv(istep) = diff / option%dt 
      end do
      diff = movave%data(ngrid) - movave%data(ngrid - 1)
      movave%deriv(ngrid - 1) = diff / option%dt

      ! Memory Dealloc.
      !
      deallocate(input_data, range)                 

    end subroutine movave_analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine movave_write(output, movave)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in) :: output
      type(s_movave), intent(in) :: movave

      integer                :: istep
      character(len=MaxChar) :: fwrite


      write(fwrite,'(a,".movave")') trim(output%fhead)
      open(10, file = trim(fwrite))
        do istep = 0, movave%ngrid - 1
          write(10,'(3(e15.7,2x))') &
            movave%grid (istep),   &
            movave%data (istep),   &
            movave%deriv(istep)
        end do
      close(10)

    end subroutine movave_write
!-----------------------------------------------------------------------

end module mod_movave_analyze
!=======================================================================
