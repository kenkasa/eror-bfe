!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj(2)

      type(s_dcd)   :: dcd 
      type(xtcfile) :: xtc 
      type(s_com)   :: com(2)

      integer :: iunit, trajtype, natm
      integer :: i, istep, istep_tot, nstep_tot
      integer :: iatm, ires, ires_prev, id, ixyz, iz, ierr
      integer :: jstep, ij, nend, itraj, nmol_large
      integer :: imol, jmol, ig
      integer :: imsta, imend, jmsta, jmend, nself
      integer :: nmol(2), nr, npair, npair_distinct
      real(8) :: boxave(3), vol, dv, zsta, zmin, zmax
      real(8) :: m, d(3), d2, r, fourpi, fact
      real(8) :: grsum
      logical :: is_end

      character(len=MaxChar) :: fname

      real(8), allocatable :: hist(:), hist_self(:), hist_distinct(:)
      real(8), allocatable :: gr(:), gr_self(:), gr_distinct(:)
      real(8), allocatable :: cumgr(:)


!$    write(iw,'("Analyze> OpenMP parallelization is activated")')

      ! Get trajectory types
      !
      call get_trajtype(input%ftraj(1), trajtype)

      ! Setup molecule
      !
      do i = 1, 2
        call get_com(option%mode(i),       &
                     traj(i),              &
                     com(i),               &
                     setup = .true.,       &
                     calc_coord = .false., &
                     myrank = 0)

        nmol(i) = com(i)%nmol
      end do

      ! get boxsize at first step for determining # of grids
      !
      call open_trajfile(input%ftraj(1), trajtype, iunit, dcd, xtc)
      call init_trajfile(trajtype, iunit, dcd, xtc, natm)
      call read_trajfile_oneframe(trajtype, iunit, 1, dcd, xtc, is_end)
      call send_coord_to_traj(1, trajtype, dcd, xtc, traj(1))

      boxave(:) = traj(1)%box(:, 1)
      zmax      = boxave(3)
      nr        = zmax / option%dr + 1

      call close_trajfile(trajtype, iunit, dcd, xtc)
      boxave(:) = 0.0d0

      ! get constants for normalization
      !
      fourpi = 4.0d0 * PI
      if (option%identical) then
        npair = nmol(1) * (nmol(1) - 1) / 2
        imsta = 1
        imend = nmol(1) - 1
        jmend = nmol(1)
      else
        npair = nmol(1) * nmol(2)
        imsta = 1
        imend = nmol(1)
        jmend = nmol(2)
      end if

      npair_distinct = nmol(1) * (nmol(1) - 1) / 2

      ! allocate memory
      !
      allocate(gr(0:nr), gr_self(0:nr), gr_distinct(0:nr))
      allocate(hist(0:nr), hist_self(0:nr), hist_distinct(0:nr))
      allocate(cumgr(0:nr))

      ! calculate RDF 
      !
      gr          = 0.0d0
      gr_self     = 0.0d0
      gr_distinct = 0.0d0

      hist          = 0.0d0
      hist_self     = 0.0d0
      hist_distinct = 0.0d0

      nself       = 0

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
            write(iw,'("Step ",i0)') istep_tot
          end if

          boxave(:) = boxave(:) + traj(1)%box(:, 1)

          do i = 1, 2
            call send_coord_to_traj(1, trajtype, dcd, xtc, traj(i))
            call get_com(option%mode(i),       &
                         traj(i),              &
                         com(i),               &
                         setup = .false.,      &
                         calc_coord = .true.,  &
                         myrank = 1)

          end do

          call update_hist(option,        &
                           traj,          &
                           com,           &
                           nr,            &
                           imsta,         &
                           imend,         &
                           jmend,         &
                           hist,          &
                           hist_self,     &
                           hist_distinct)

        end do ! step

        istep_tot = istep_tot - 1

        call close_trajfile(trajtype, iunit, dcd, xtc)

      end do   ! traj

      nstep_tot = istep_tot
      boxave(:) = boxave(:) / dble(nstep_tot)
      vol       = boxave(1) * boxave(2) * boxave(3)

      call normalize_hist(option,          &
                          vol,             &
                          nstep_tot,       &
                          nmol,            &
                          nr,              &
                          npair,           &
                          npair_distinct,  &
                          hist,            &
                          hist_self,       &
                          hist_distinct,   &
                          gr,              &
                          gr_self,         &
                          gr_distinct)

      ! for checking normalization
      if (option%separate_self) then
        grsum = 0.0d0
        do ig = 1, nr
          r = dble(ig) * option%dr
          fact = fourpi * r * r * option%dr * gr_self(ig) 
          grsum = grsum + fact 
        end do
        write(iw,'("grsum = ", f20.10)') grsum
      end if

      ! calculate cumlative  
      !
      cumgr = 0.0d0
      grsum = 0.0d0
      do ig = 1, nr
        r = option%dr * dble(ig)
        grsum     = grsum + fourpi * r * r * option%dr * gr(ig)
        cumgr(ig) = grsum 
      end do

      ! Generate files
      !
      ! - Output RDF

      call gen_file_name(output%fhead, "rd", fname, output%frd)
      call open_file(fname, iunit)

      if (option%separate_self) then
        do ig = 0, nr 
          r = option%dr * dble(ig) 
          write(iunit,'(4f20.10)') r, gr(ig), gr_distinct(ig), gr_self(ig) 
        end do
      else
        do ig = 0, nr
          r = option%dr * dble(ig)
          write(iunit,'(3f20.10)') r, gr(ig), cumgr(ig)
        end do
      end if

      close(iunit)

      ! - Output binary format RDF
      !
      call gen_file_name(output%fhead, "rdbin", fname, output%frdbin)
      call open_file(fname, iunit, "unformatted")

      write(iunit) option%identical
      write(iunit) option%separate_self
      !write(iunit) traj(1)%nstep
      write(iunit) nstep_tot
      write(iunit) nmol
      write(iunit) vol
      write(iunit) option%dr
      write(iunit) nr
      write(iunit) hist, hist_distinct, hist_self

      close(iunit)
     
      ! Deallocate memory
      !
      deallocate(gr, gr_self, gr_distinct)
      deallocate(hist, hist_distinct, hist_self)
      deallocate(cumgr)


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine update_hist(option,        &
                           traj,          &
                           com,           &
                           nr,            &
                           imsta,         &
                           imend,         &
                           jmend,         &
                           hist,          &
                           hist_self,     &
                           hist_distinct)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)    :: option
      type(s_traj),   intent(in)    :: traj(2)
      type(s_com),    intent(in)    :: com(2)
      integer,        intent(in)    :: nr
      integer,        intent(in)    :: imsta
      integer,        intent(in)    :: imend
      integer,        intent(in)    :: jmend
      real(8),        intent(inout) :: hist(0:nr)
      real(8),        intent(inout) :: hist_self(0:nr)
      real(8),        intent(inout) :: hist_distinct(0:nr)

      integer :: imol, jmol, jmsta, ig
      real(8) :: d(3), r


!$omp parallel private(imol, jmol, jmsta, d, r, ig), shared(traj, com) default(shared), &
!$omp        & reduction(+:hist), reduction(+:hist_self), reduction(+:hist_distinct)
!$omp do
      do imol = imsta, imend
      
        if (option%identical) then
          jmsta = imol + 1      
        else
          jmsta = 1
        end if
      
        do jmol = jmsta, jmend
          d(1:3) = com(2)%coord(1:3, jmol, 1) &
                 - com(1)%coord(1:3, imol, 1)
          d(1:3) = d(1:3) - traj(1)%box(1:3, 1) &
                          * nint(d(1:3) / traj(1)%box(1:3, 1))
          r      = sqrt(dot_product(d, d))
          ig     = nint(r/option%dr)
      
          if (ig >= 0 .and. ig <= nr) then
            hist(ig) = hist(ig) + 1.0d0
            !gr(ig) = gr(ig) + 1.0d0
            if (imol == jmol) then
              !nself       = nself + 1
              hist_self(ig) = hist_self(ig) + 1.0d0
              !gr_self(ig)   = gr_self(ig) + 1.0d0
            else
              hist_distinct(ig) = hist_distinct(ig) + 1.0d0 
              !gr_distinct(ig)   = gr_distinct(ig) + 1.0d0
            end if
          end if
      
        end do
      
      end do
!$omp end do
!$omp end parallel


    end subroutine update_hist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine normalize_hist(option,          &
                              vol,             &
                              nstep,           &
                              nmol,            &
                              nr,              &
                              npair,           &
                              npair_distinct,  &
                              hist,            &
                              hist_self,       &
                              hist_distinct,   &
                              gr,              &
                              gr_self,         &
                              gr_distinct)
!-----------------------------------------------------------------------
      implicit none

      type(s_option), intent(in)  :: option
      real(8),        intent(in)  :: vol
      integer,        intent(in)  :: nstep
      integer,        intent(in)  :: nmol(2)
      integer,        intent(in)  :: nr
      integer,        intent(in)  :: npair
      integer,        intent(in)  :: npair_distinct
      real(8),        intent(in)  :: hist(0:nr)
      real(8),        intent(in)  :: hist_self(0:nr)
      real(8),        intent(in)  :: hist_distinct(0:nr)
      real(8),        intent(out) :: gr(0:nr)
      real(8),        intent(out) :: gr_self(0:nr)
      real(8),        intent(out) :: gr_distinct(0:nr)

      integer :: ig
      real(8) :: r, fact, fourpi


      fourpi = 4.0d0 * PI

      if (option%normalize) then
        do ig = 1, nr
          r      = dble(ig) * option%dr
          fact   = vol / (fourpi * r * r * option%dr * npair * nstep)
          gr(ig) = hist(ig) * fact

          if (option%separate_self) then
            fact   = 1.0d0 / (fourpi * r * r * option%dr * nmol(1) * nstep)
            gr_self(ig) = hist_self(ig) * fact
            if (npair_distinct == 0) then
              gr_distinct(ig) = 0.0d0
            else
              fact   = vol / (fourpi * r * r * option%dr * npair_distinct * nstep)
              gr_distinct(ig) = hist_distinct(ig) * fact 
            end if
          end if
        end do
      else
        do ig = 1, nr
          r      = dble(ig) * option%dr
          fact   = nmol(2) / (fourpi * r * r * option%dr * npair * nstep)
          gr(ig) = hist(ig) * fact 
          if (option%separate_self) then
            fact   = 1.0d0 / (fourpi * r * r * option%dr * nmol(1) * nstep)
            gr_self(ig) = hist_self(ig) * fact 
            if (npair_distinct == 0) then
              gr_distinct(ig) = 0.0d0
            else
              fact   = (nmol(1) - 1) / (fourpi * r * r * option%dr * npair_distinct * nstep)
              gr_distinct(ig) = hist_distinct(ig) * fact
            end if
          end if
        end do

      end if


    end subroutine normalize_hist
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
