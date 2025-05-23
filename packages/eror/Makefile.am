ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = erdst slvfe
nobase_pkgdata_SCRIPTS = \
	tools/AMBER/gen_structure tools/AMBER/gen_input \
	tools/GROMACS/gen_structure tools/GROMACS/gen_input \
	tools/NAMD/gen_structure tools/NAMD/gen_input \
	tools/LAMMPS/gen_structure tools/LAMMPS/gen_input \
	tools/CHARMM/gen_structure tools/CHARMM/gen_input \
	tools/GENESIS/gen_structure tools/GENESIS/gen_input \
	tools/GENESIS/gen_structure_AMBER tools/GENESIS/gen_structure_CHARMM \
	tools/bootstrap \
	tools/refs_converge

#	tools/GENESIS/gen_structure_GROMACS # <- In GENESIS, GROMACS is unsupported officially so not included in SCRIPTS


install-exec-hook:
	cd $(DESTDIR)$(bindir) && rm -f ermod$(EXEEXT) && $(LN_S) erdst$(EXEEXT) ermod$(EXEEXT)

uninstall-hook:
	cd $(DESTDIR)$(bindir) && rm -f ermod$(EXEEXT)

dist_pkglibexec_SCRIPTS = $(top_srcdir)/vmdplugins/libexec/*.so
subdir=../anatra/f90/lib

anatra_sources =     \
	anatra_util.f90    \
	anatra_const.f90   \
	anatra_output.f90  \
	anatra_anaparm.f90 \
	anatra_prmtop.f90  \
	anatra_grid3d.f90  \
	anatra_xyz.f90     \
	anatra_dcd.f90     \
	anatra_traj.f90    \
	anatra_fitting.f90 \
	anatra_pot.f90     \
        anatra_com.f90     \
	anatra_input.f90   \
        anatra_lib.f90     \
	anatra_energy.f90 



ersources = enganal.F90 engmain.F90 realcal.F90 spline.F90 fft_iface.F90 insertion.F90 mpiproc.F90 engproc.F90 recpcal.F90 setconf.F90 trajectory.F90 utility.F90 bestfit.F90 randgen.F90

EXTRA_DIST = parameters_er.example parameters_fe.example
EXTRA_DIST += extensions/*
EXTRA_DIST += trajectory-formats/*
EXTRA_DIST += tools/AMBER/*
EXTRA_DIST += tools/GROMACS/*
EXTRA_DIST += tools/NAMD/*
EXTRA_DIST += tools/LAMMPS/*
EXTRA_DIST += tools/CHARMM/*
EXTRA_DIST += tools/GENESIS/*
EXTRA_DIST += tools/mol_dissection
EXTRA_DIST += tools/refs_converge
EXTRA_DIST += tools/bootstrap
EXTRA_DIST += vmdplugins/*.C vmdplugins/*.c vmdplugins/*.h vmdplugins/*.cxx vmdplugins/*.hxx vmdplugins/Makefile vmdplugins/LICENSE
CLEANFILES = *.mod

if USE_FFTW
FFTFLAGS = -DFFTW
FFTLIBS = $(FFTWLIBS)
endif

if USE_MKL
FFTFLAGS = -DMKL
endif

erdst_SOURCES = $(anatra_sources) $(ersources)
erdst_CFLAGS = -DINSTALL_PLUGIN_PATH=\"$(pkglibexecdir)\"
erdst_FCFLAGS = $(VMDPLUGINFLAGS) $(FFTFLAGS)
erdst_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS) -ldl -lm 

if USE_VMDPLUGIN
erdst_SOURCES += vmdfio.c xoshiro256ss.c
endif

if USE_TRANSFER
erdst_FCFLAGS += -DHAVE_TRANSFER
else
erdst_SOURCES += hash_real.c
endif

if USE_FLUSH
erdst_FCFLAGS += -DHAVE_FLUSH
endif

extra_flags = 

# module dependencies
$(top_builddir)/anatra-energy.$(OBJEXT): $(top_srcdir)/anatra_energy.f90 \
    $(top_builddir)/anatra-lib.$(OBJEXT)     \
    $(top_builddir)/anatra-pot.$(OBJEXT)     \
    $(top_builddir)/anatra-com.$(OBJEXT)     \
    $(top_builddir)/anatra-traj.$(OBJEXT)    \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-lib.$(OBJEXT): $(top_srcdir)/anatra_lib.f90 \
    $(top_builddir)/anatra-com.$(OBJEXT)     \
    $(top_builddir)/anatra-fitting.$(OBJEXT) \
    $(top_builddir)/anatra-traj.$(OBJEXT)    \
    $(top_builddir)/anatra-grid3d.$(OBJEXT)  \
    $(top_builddir)/anatra-prmtop.$(OBJEXT)  \
    $(top_builddir)/anatra-input.$(OBJEXT)   \
    $(top_builddir)/anatra-output.$(OBJEXT)  \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-fitting.$(OBJEXT): $(top_srcdir)/anatra_fitting.f90 \
    $(top_builddir)/anatra-traj.$(OBJEXT)    \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-anaparm.$(OBJEXT): $(top_srcdir)/anatra_anaparm.f90 \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-output.$(OBJEXT): $(top_srcdir)/anatra_output.f90 \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-grid3d.$(OBJEXT): $(top_srcdir)/anatra_grid3d.f90 \
    $(top_builddir)/anatra-const.$(OBJEXT) 

$(top_builddir)/anatra-xyz.$(OBJEXT): $(top_srcdir)/anatra_xyz.f90 \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-input.$(OBJEXT): $(top_srcdir)/anatra_input.f90 \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-traj.$(OBJEXT): $(top_srcdir)/anatra_com.f90 \
    $(top_builddir)/anatra-traj.$(OBJEXT)    \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-pot.$(OBJEXT): $(top_srcdir)/anatra_pot.f90 \
    $(top_builddir)/anatra-anaparm.$(OBJEXT)  \
    $(top_builddir)/anatra-prmtop.$(OBJEXT)   \
    $(top_builddir)/anatra-traj.$(OBJEXT)     \
    $(top_builddir)/anatra-const.$(OBJEXT)    \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-prmtop.$(OBJEXT): $(top_srcdir)/anatra_prmtop.f90 \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-com.$(OBJEXT): $(top_srcdir)/anatra_traj.f90 \
    $(top_builddir)/anatra-const.$(OBJEXT)   \
    $(top_builddir)/anatra-util.$(OBJEXT) 

$(top_builddir)/anatra-const.$(OBJEXT): $(top_srcdir)/anatra_const.f90 \
    $(top_builddir)/anatra-util.$(OBJEXT) 



$(top_builddir)/erdst-enganal.$(OBJEXT): $(top_srcdir)/enganal.F90 \
    $(top_builddir)/erdst-engproc.$(OBJEXT) \
    $(top_builddir)/erdst-setconf.$(OBJEXT) \
    $(top_builddir)/erdst-insertion.$(OBJEXT) \
    $(top_builddir)/erdst-trajectory.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-engproc.$(OBJEXT): $(top_srcdir)/engproc.F90 \
    $(top_builddir)/erdst-insertion.$(OBJEXT) \
    $(top_builddir)/erdst-realcal.$(OBJEXT) \
    $(top_builddir)/erdst-recpcal.$(OBJEXT) \
    $(top_builddir)/erdst-bestfit.$(OBJEXT) \
    $(top_builddir)/erdst-utility.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-randgen.$(OBJEXT): $(top_srcdir)/randgen.F90
$(top_builddir)/erdst-realcal.$(OBJEXT): $(top_srcdir)/realcal.F90 \
    $(top_builddir)/erdst-engmain.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT)
$(top_builddir)/erdst-recpcal.$(OBJEXT): $(top_srcdir)/recpcal.F90 \
    $(top_builddir)/erdst-spline.$(OBJEXT) \
    $(top_builddir)/erdst-fft_iface.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-insertion.$(OBJEXT): $(top_srcdir)/insertion.F90 \
    $(top_builddir)/erdst-setconf.$(OBJEXT) \
    $(top_builddir)/erdst-bestfit.$(OBJEXT) \
    $(top_builddir)/erdst-trajectory.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-randgen.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-setconf.$(OBJEXT): $(top_srcdir)/setconf.F90 \
    $(top_builddir)/erdst-utility.$(OBJEXT) \
    $(top_builddir)/erdst-trajectory.$(OBJEXT) \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT) \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-bestfit.$(OBJEXT): $(top_srcdir)/bestfit.F90 \
    $(top_builddir)/erdst-mpiproc.$(OBJEXT)
$(top_builddir)/erdst-trajectory.$(OBJEXT): $(top_srcdir)/trajectory.F90 \
    $(top_builddir)/erdst-utility.$(OBJEXT)
$(top_builddir)/erdst-utility.$(OBJEXT): $(top_srcdir)/utility.F90 \
    $(top_builddir)/erdst-engmain.$(OBJEXT)
$(top_builddir)/erdst-mpiproc.$(OBJEXT): $(top_srcdir)/mpiproc.F90 \
    $(top_builddir)/erdst-engmain.$(OBJEXT)

# slvfe part
slvfe_SOURCES = slvfe.F90 sfemain.F90 sfecorrect.F90 engmain.F90 $(anatra_sources)
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
slvfe_LINK = $(FC_CANONICAL) $(FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
slvfe.$(OBJEXT): slvfe.F90 sfemain.$(OBJEXT) sfecorrect.$(OBJEXT) engmain.$(OBJEXT)
sfemain.$(OBJEXT): sfemain.F90 engmain.$(OBJEXT)
sfecorrect.$(OBJEXT): sfecorrect.F90 sfemain.$(OBJEXT)
engmain.$(OBJEXT): engmain.F90

# these are necessary to compile slvfe's sources without mpi
$(top_builddir)/slvfe.$(OBJEXT): $(top_srcdir)/slvfe.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/sfemain.$(OBJEXT): $(top_srcdir)/sfemain.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/sfecorrect.$(OBJEXT): $(top_srcdir)/sfecorrect.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
$(top_builddir)/engmain.$(OBJEXT): $(top_srcdir)/engmain.F90
	$(FC_CANONICAL) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	    $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS) $(slvfe_FCFLAGS) -c -o $@ $<
