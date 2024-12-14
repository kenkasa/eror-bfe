ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = ermod slvfe
nobase_pkgdata_SCRIPTS = \
	tools/AMBER/gen_structure tools/AMBER/gen_input \
	tools/gromacs/gen_structure tools/gromacs/gen_input \
	tools/NAMD/gen_structure tools/NAMD/gen_input \
	tools/MARBLE/gen_structure tools/MARBLE/gen_input \
	tools/MODYLAS/gen_structure tools/MODYLAS/gen_input \
	tools/LAMMPS/gen_structure tools/LAMMPS/gen_input \
	tools/CHARMM/gen_structure tools/CHARMM/gen_input \
	tools/GENESIS/gen_structure_CHARMM  tools/GENESIS/gen_structure_AMBER \
	tools/GENESIS/gen_structure_GROMACS tools/GENESIS/gen_input \
	tools/refs_converge tools/refs_converge_py2 tools/refs_converge_py3

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

ersources = enganal.F90 engmain.F90 realcal.F90 spline.F90 fft_iface.F90 mpiproc.F90 engproc.F90 recpcal.F90 insertion.F90 setconf.F90 trajectory.F90 utility.F90 bestfit.F90

EXTRA_DIST = parameters_er.example parameters_fe.example
EXTRA_DIST += extensions/*
EXTRA_DIST += trajectory-formats/*
EXTRA_DIST += tools/AMBER/*
EXTRA_DIST += tools/gromacs/*
EXTRA_DIST += tools/NAMD/*
EXTRA_DIST += tools/MARBLE/*
EXTRA_DIST += tools/MODYLAS/*
EXTRA_DIST += tools/LAMMPS/*
EXTRA_DIST += tools/CHARMM/*
EXTRA_DIST += tools/GENESIS/*
EXTRA_DIST += tools/mol_dissection
EXTRA_DIST += tools/refs_converge tools/refs_converge_py2 tools/refs_converge_py3
EXTRA_DIST += vmdplugins/*.C vmdplugins/*.c vmdplugins/*.h vmdplugins/*.cxx vmdplugins/*.hxx vmdplugins/Makefile vmdplugins/LICENSE
CLEANFILES = *.mod

if USE_FFTW
FFTFLAGS = -DFFTW
FFTLIBS = $(FFTWLIBS)
endif

if USE_MKL
FFTFLAGS = -DMKL
endif

ermod_SOURCES = $(anatra_sources) $(ersources)
ermod_CFLAGS = -DINSTALL_PLUGIN_PATH=\"$(pkglibexecdir)\"
ermod_FCFLAGS = $(VMDPLUGINFLAGS) $(FFTFLAGS)
ermod_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS) $(FFTLIBS) $(MPILIBS) -ldl -lm 

if USE_VMDPLUGIN
ermod_SOURCES += vmdfio.c
endif

if USE_TRANSFER
ermod_FCFLAGS += -DHAVE_TRANSFER
else
ermod_SOURCES += hash_real.c
endif

if USE_FLUSH
ermod_FCFLAGS += -DHAVE_FLUSH
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

$(top_builddir)/ermod-enganal.$(OBJEXT): $(top_srcdir)/enganal.F90 \
    $(top_builddir)/ermod-engproc.$(OBJEXT) \
    $(top_builddir)/ermod-setconf.$(OBJEXT) \
    $(top_builddir)/ermod-insertion.$(OBJEXT) \
    $(top_builddir)/ermod-trajectory.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-engproc.$(OBJEXT): $(top_srcdir)/engproc.F90 \
    $(top_builddir)/ermod-insertion.$(OBJEXT) \
    $(top_builddir)/ermod-realcal.$(OBJEXT) \
    $(top_builddir)/ermod-recpcal.$(OBJEXT) \
    $(top_builddir)/ermod-bestfit.$(OBJEXT) \
    $(top_builddir)/ermod-utility.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-realcal.$(OBJEXT): $(top_srcdir)/realcal.F90 \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-recpcal.$(OBJEXT): $(top_srcdir)/recpcal.F90 \
    $(top_builddir)/ermod-spline.$(OBJEXT) \
    $(top_builddir)/ermod-fft_iface.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-insertion.$(OBJEXT): $(top_srcdir)/insertion.F90 \
    $(top_builddir)/ermod-setconf.$(OBJEXT) \
    $(top_builddir)/ermod-bestfit.$(OBJEXT) \
    $(top_builddir)/ermod-trajectory.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-setconf.$(OBJEXT): $(top_srcdir)/setconf.F90 \
    $(top_builddir)/ermod-utility.$(OBJEXT) \
    $(top_builddir)/ermod-trajectory.$(OBJEXT) \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT) \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-bestfit.$(OBJEXT): $(top_srcdir)/bestfit.F90 \
    $(top_builddir)/ermod-mpiproc.$(OBJEXT)
$(top_builddir)/ermod-trajectory.$(OBJEXT): $(top_srcdir)/trajectory.F90 \
    $(top_builddir)/ermod-utility.$(OBJEXT)
$(top_builddir)/ermod-utility.$(OBJEXT): $(top_srcdir)/utility.F90 \
    $(top_builddir)/ermod-engmain.$(OBJEXT)
$(top_builddir)/ermod-mpiproc.$(OBJEXT): $(top_srcdir)/mpiproc.F90 \
    $(top_builddir)/ermod-engmain.$(OBJEXT)

# slvfe part
slvfe_SOURCES = slvfe.F90 sfemain.F90 sfecorrect.F90
slvfe_LDADD = $(LAPACK_LIBS) $(BLAS_LIBS)
slvfe_LINK = $(FC_CANONICAL) $(FCFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
slvfe.$(OBJEXT): slvfe.F90 sfemain.$(OBJEXT) sfecorrect.$(OBJEXT)
sfemain.$(OBJEXT): sfemain.F90
sfecorrect.$(OBJEXT): sfecorrect.F90 sfemain.$(OBJEXT)

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
