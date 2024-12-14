!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_prmtop
  !use mod_potential

  ! constants
  !

  ! subroutines
  !
  public  :: analyze

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option)
!-----------------------------------------------------------------------
      implicit none

      ! parameters
      real(8), parameter :: charge_const = 332.05221729d0 

      ! formal arguments
      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option

      type(s_prmtop)         :: prmtop
      integer                :: io_i, io_o, io_c, iatm, natm, ierr
      real(8)                :: cval
      character(len=MaxChar) :: fout, line

      real(8), allocatable   :: charge_original(:), charge_modified(:)


      ! Read parameter file
      !
      write(iw,'("")')
      write(iw,'("Analyze> Read potential parameter files")')
      call read_prmtop(input%fprmtop, prmtop)
      
      natm = prmtop%natom

      ! Allocate memory
      !
      allocate(charge_original(1:natm), charge_modified(1:natm))

      ! Load charge info.
      !
      do iatm = 1, natm
        charge_original(1:natm) = prmtop%charge(1:natm)
        charge_modified(1:natm) = prmtop%charge(1:natm)
      end do

      ! Read charge file
      !
      call open_file(option%fcharge, io_c)
      do while (.true.)
        read(io_c,*, end=100) iatm, cval
        charge_modified(iatm) = sqrt(charge_const) * cval
      end do
100   close(io_c)

      ! Generate modified prmtop
      !
      write(fout,'(a,".prmtop")') trim(output%fhead) 

      call open_file(input%fprmtop, io_i)
      call open_file(fout,          io_o)

      do while (.true.)
        read(io_i,'(a)', end=101) line

        if (line(1:12) == '%FLAG CHARGE') then
          write(io_o,'(a)') line(1:80)
          read(io_i, '(a)') line
          write(io_o,'(a)') line(1:80)
          do iatm = 1, natm
            write(io_o,'(e16.8)',advance='no') charge_modified(iatm)
            if (mod(iatm, 5) == 0) &
              write(io_o,*)
          end do
          if (mod(natm, 5) /= 0) &
            write(io_o,*)

          call seek_prmtop_next_flag(io_i, line, ierr)
          write(io_o,'(a)') line(1:80)
          write(iw,'(a)') line(1:80)
        else
          write(io_o,'(a)') line(1:80)
        end if
      end do

101   continue

      close(io_i)
      close(io_o)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
