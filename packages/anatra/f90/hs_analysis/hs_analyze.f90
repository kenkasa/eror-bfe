!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_ctrl

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)   :: input
      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option

      integer                :: istep, ix, icol, np
      integer                :: ista, iend, npt
      real(8)                :: x, ave 
      character(len=MaxChar) :: line

      real(8), allocatable   :: dat(:, :), distr(:) 

      np = 0
      open(10,file=trim(input%fts))
      do while (.true.) 
        read(10,*,end=100) line
        np = np + 1 
      end do

100   continue
      allocate(dat(option%ncol, np), distr(0:option%nx))
      rewind 10
      do istep = 1, np
        read(10,*) x, (dat(icol, istep), icol = 1, option%ncol) 
      end do
      close(10)

      if (option%data_sta == 0 .and. option%data_end == 0) then
        ista = 1
        iend = np
        npt  = np
      else
        ista = option%data_sta
        iend = option%data_end
        npt  = iend - ista + 1
      end if

      distr = 0.0d0
      do istep = ista, iend 
        do icol = 1, option%ncol
          !ix        = (dat(icol, istep) - option%xsta) / option%dx
          !if (ix >= 0 .and. ix <= option%nx) then
          !  distr(ix) = distr(ix) + 1.0d0
          !end if
          
          ix        = nint((dat(icol, istep) - option%xsta) / option%dx)
          if (ix > 0 .and. ix < option%nx) then
            distr(ix) = distr(ix) + 1.0d0
          else if (ix == 0 .or. ix == option%nx) then
            distr(ix) = distr(ix) + 2.0d0
          end if
        end do
      end do

      distr = distr / (option%ncol * npt * option%dx)
      ave   = sum(dat(:, :)) / dble(npt * option%ncol)

      open(10,file=trim(output%fdist))
      write(10,'("# Average value : ",f20.10)') ave
      do ix = 0, option%nx
        x = option%xsta + ix * option%dx
        write(10,'(2f20.10)') x, distr(ix)
      end do 
      close(10)


    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
