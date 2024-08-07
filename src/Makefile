#------------------------------------------------------------------------------#
# Makefile for Monte Carlo radiative transfer code MCFOST
#------------------------------------------------------------------------------#

IS_DARWIN := $(shell uname | tr '[a-z]' '[A-Z]' 2>&1 | grep -c DARWIN)
cpu =  $(shell uname -m)

KNOWN_SYSTEM=no

LIB= $(MCFOST_INSTALL)/lib

ifndef SYSTEM
    SYSTEM=ifort
endif

ifeq ($(gfortran), yes)
    SYSTEM=gfortran
endif

ifeq ($(ifort), yes)
    SYSTEM=ifort
endif

ifeq ($(ifx), yes)
    SYSTEM=ifx
endif

ifeq ($(release), yes)
    RELEASE=yes
endif

ifndef OPENMP
    OPENMP=yes
endif

# faster compilation
ifeq ($(dev), yes)
    DEV=yes
endif

ifeq ($(debug), yes)
    DEBUG=yes
endif

ifeq ($(mono), yes)
    OPENMP=no
endif

ifeq ($(DEBUG), yes)
    OPENMP=no
endif

ifeq ($(debugp), yes)
    DEBUGP=yes
endif

ifeq ($(DEBUGP), yes)
    DEBUG=yes
    OPENMP=yes
endif

ifeq ($(profiling), yes)
    PROFILING=yes
endif

ifeq ($(static), yes)
    STATIC=yes
endif

ifdef xgboost
    XGBOOST=$(xgboost)
endif

ifndef XGBOOST
    ifeq ($(MCFOST_XGBOOST), yes)
        XGBOOST = yes
    else
        XGBOOST = no
    endif
endif

#------------------------------------------------------------------------------#
# System dependent compiler flags

ifeq ($(SYSTEM), ifort)
    STATIC_FLAGS = -static-intel -qopenmp-link static
    ifeq ($(IS_DARWIN),1)
        ARCH = -axSSSE3,SSE4.1,SSE4.2,AVX,CORE-AVX2 -mmacosx-version-min=10.13 -mdynamic-no-pic
    else
        ARCH = -axSSE2,SSSE3,SSE4.1,SSE4.2,AVX,CORE-AVX2,CORE-AVX512
        STATIC_FLAGS += -static-libstdc++
    endif
    FC= ifort
    FFLAGS= -fpp -O3 -fp-model fast=2 -traceback
    DBLFLAG= -r8
    DEBUGFLAG= -check all -g -fpe0 -traceback -no-ftz \
               -warn uninitialized -warn unused -warn truncated_source
    KNOWN_SYSTEM= yes
    FOMPFLAGS= -fopenmp
    IPOFLAGS= -ipo
    CC= icc
    CXX=icpc -std=c++11
    CCFLAGS= -O3  -fp-model fast=2 -traceback
    COMPFLAGS= -fopenmp
    LIBCXX= -cxxlib
    LIBS= $(LIB)/ifort
endif

ifeq ($(SYSTEM), ifx)
    STATIC_FLAGS = -static-intel -qopenmp-link static
    ifeq ($(IS_DARWIN),1)
        ARCH = -axSSSE3,SSE4.1,SSE4.2,AVX,CORE-AVX2 -mmacosx-version-min=10.13 -mdynamic-no-pic
    else
        ARCH = -axSSSE3,SSE4.1,SSE4.2,AVX,CORE-AVX2,CORE-AVX512
        STATIC_FLAGS += -static-libstdc++
    endif
    FC= ifx
    FFLAGS= -fpp -O3 -fp-model fast=2 -traceback
    DBLFLAG= -r8
    DEBUGFLAG= -check all -C -g -fpe0 -traceback -no-ftz \
               -warn uninitialized -warn unused -warn truncated_source
    KNOWN_SYSTEM= yes
    FOMPFLAGS= -qopenmp
    IPOFLAGS= -ipo
    CC= icx
    CXX=icpx -std=c++11
    CCFLAGS= -O3 -fp-model fast=2 -traceback
    COMPFLAGS= -qopenmp
    LIBCXX= -cxxlib
    LIBS= $(LIB)/ifx
endif

ifeq ($(SYSTEM), gfortran)
    FC= gfortran
    FFLAGS+= -m64 -cpp -O3 -funroll-loops -ffinite-math-only -fno-trapping-math \
	     -ffree-line-length-132 -fbacktrace
    STATIC_FLAGS= -static-libgcc -static-libstdc++ -static-libgfortran -static-libquadmath
    DBLFLAG= -fdefault-real-8 -fdefault-double-8
    DEBUGFLAG= -g -fbounds-check -Wunused-variable -Wunused-function -Wuninitialized -fbacktrace
    KNOWN_SYSTEM= yes
    FOMPFLAGS= -fopenmp
    CC= gcc
    CXX= g++
    CCFLAGS= -Wall -pedantic -m64 -O3 -funroll-loops -ffinite-math-only -fno-trapping-math
    COMPFLAGS= -fopenmp
    ifeq ($(cpu),arm64)
        LIBCXX= -lstdc++ -Wl,-ld_classic # for github action
    else
        LIBCXX= -lstdc++
    endif
    ifneq ($(shell $(CC) -v 2>&1 | grep -c "clang"), 0)
        LIBCXX= -lc++
    endif
    LIBS= $(LIB)/gfortran
endif

#------------------------------------------------------------------------------#
# Set other optional flags depending on system settings

# Force static flags when compiling the releases, and march=native if not a release
ifeq ($(RELEASE), yes)
    STATIC = yes
    FFLAGS += $(ARCH)
    CCFLAGS += $(ARCH)
else
    ifeq ($(SYSTEM), gfortran)
	FFLAGS += -march=native
        CCFLAGS += -march=native
    else
	FFLAGS += -xHOST #-march=native
        CCFLAGS += -xHOST #-march=native
    endif
endif

ifeq ($(DEBUG), yes)
    FFLAGS+= ${DEBUGFLAG}
    FFLAGS:= $(FFLAGS:-O3=-O0)
endif

ifeq ($(PROFILING), yes)
    FFLAGS+= -g3
endif

ifeq ($(DEV), yes)
    FFLAGS:= $(FFLAGS:-O3=-O1)
endif

ifeq ($(OPENMP), yes)
    FFLAGS+= ${FOMPFLAGS}
    CCLAGS+= ${FOMPFLAGS}
endif

ifeq ($(IPO), yes)
    FFLAGS+= ${IPOFLAGS}
    CCFLAGS+= ${IPOFLAGS}
endif


#-- optional compilation with xgboost, as there are some issues on some systems
ifeq ($(XGBOOST), yes)
    XGBOOST_WRAPPER = xgboost_wrapper
    LXGBOOST= -lxgboost -ldmlc -lrabit
else
    XGBOOST_WRAPPER = no_xgboost_wrapper
endif

#-- HDF5 library
ifndef HDF5_DIR
    HDF5_LIB= $(LIBS)
    HDF5_INC= $(MCFOST_INSTALL)/include/hdf5
endif
ifdef HDF5_DIR
    HDF5_LIB= $(HDF5_DIR)/lib
    HDF5_INC= $(HDF5_DIR)/include
endif

#------------------------------------------------------------------------------#
# Library flags

INCLUDE+= -I$(MCFOST_INSTALL)/include \
          -I$(MCFOST_INSTALL)/include/voro++ \
          -I$(MCFOST_INSTALL)/include/$(SYSTEM) \
          -I$(HDF5_INC)
LDFLAGS+= -L$(LIBS) -L$(HDF5_LIB) -lcfitsio -lvoro++ $(LIBCXX) \
          -lsprng $(LXGBOOST) -lhdf5_fortran -lhdf5 -lz -ldl
ifeq ($(STATIC), yes)
    LDFLAGS+= $(STATIC_FLAGS)
endif


#------------------------------------------------------------------------------#
# Sources and objects
#gas transfer. (Only atoms at the moment)
gdir=gas
VPATH = $(gdir)/

SOURCES = mcfost_env.f90 parameters.f90 constants.f90 healpix_mod.f90 sha.f90 messages.f90 \
          operating_system.f90 random_numbers.f90 utils.f90 sort.f90 fits_utils.f90 \
          grains.f90 cylindrical_grid.f90 spherical_grid.f90 kdtree2.f90 \
          elements_type.f90 Voronoi.f90 grid.f90 wavelengths.f90 stars.f90 \
          read_DustEM.f90 Temperature.f90 density.f90 read_opacity.f90 scattering.f90 \
          read_fargo3d.f90 hdf5_utils.f90 utils_hdf5.f90 read_athena++.f90	\
          readVTK.f90 read_idefix.f90 read_pluto.f90 voigts.f90 	\
          coated_sphere.f90 dust_prop.f90 molecular_emission.f90 PAH.f90            \
          read1d_models.f90 read_spherical_grid.f90	\
          input.f90 benchmarks.f90 atom_type.f90  \
          wavelengths_gas.f90 broad.f90                             \
          read_param.f90 dust_ray_tracing.f90 uplow.f90 abo.f90     \
          occupation_probability.f90 lte.f90                        \
          radiation_field.f90 thermal_emission.f90 diffusion.f90                    \
          io_prodimo.f90 disk_physics.f90 gas_contopac.f90 		     \
          opacity_atom.f90 MRW.f90 optical_depth.f90 ML_prodimo.f90                                          \
          output.f90 mem.f90  init_mcfost.f90 io_phantom_infiles.f90                \
          io_phantom_utils.f90 mess_up_SPH.f90 read_phantom.f90     \
          read_gadget2.f90 fsolve.f90 integrate.f90 reconstruct_from_moments.f90 \
	  SPH2mcfost.f90 mhd2mcfost.f90 dust_transfer.f90          \
          mol_transfer.f90 collision_atom.f90 io_atom.f90           \
          electron_density.f90 see.f90 escape.f90                             \
          atom_transfer.f90

OBJECTS= $(SOURCES:.f90=.o)

#---------	---------------------------------------------------------------------#
# Compile
main : mcfost

#include .depend

all : mcfost libmcfost.a

mcfost :  $(OBJECTS) voro++_wrapper.o  $(XGBOOST_WRAPPER).o mcfost.o
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) voro++_wrapper.o $(XGBOOST_WRAPPER).o mcfost.o $(LDFLAGS)

libmcfost.a : $(OBJECTS) voro++_wrapper.o no_xgboost_wrapper.o mcfost2phantom.o
	rm -f mcfost.o xgboost_wrapper.o libmcfost.a
	ar vr libmcfost.a *.o
	ranlib libmcfost.a

mcfost_env.o : .make_lastsystem mcfost_env.f90
	mkdir -p ../utils ; rm -rf ../utils/Version ;
	grep required_utils_version mcfost_env.f90 | awk '{print $$6}' > ../utils/Version
	${MAKE} checkmcfost
	$(FC) $(FFLAGS) -c -DPOINTER_SIZE=8 mcfost_env.f90 $(INCLUDE)

io_phantom_utils.o :  io_phantom_utils.f90
	$(FC) $(FFLAGS) -c $(DBLFLAG) io_phantom_utils.f90

voro++_wrapper.o :  voro++_wrapper.cpp
	$(CXX) $(CCFLAGS) -c voro++_wrapper.cpp -I$(MCFOST_INSTALL)/include

$(XGBOOST_WRAPPER).o :  $(XGBOOST_WRAPPER).cpp
	$(CXX) $(CCFLAGS) -c $(XGBOOST_WRAPPER).cpp -I$(MCFOST_INSTALL)/include

%.o : %.f90
	$(FC) $(FFLAGS) -c -DPOINTER_SIZE=8 $< -o $@ $(INCLUDE)

io_prodimo.o : sha.f90

utils.o : sha.f90

Voronoi.o : voro++_wrapper.o

ifeq ($(MCFOST_GIT),1)
sha.f90 : ../.git/objects mcfost_env.o
	printf "module sha\n  implicit none\n" > sha.f90
	git rev-parse HEAD | awk \
		'{print "  character(len=40), parameter :: sha_id = \""$$1"\""}' \
		>> sha.f90
	printf "end module sha\n" >> sha.f90
endif

operating_system.f90 : sha.o
	printf "module os\n  implicit none\n" > operating_system.f90
	uname -s | awk \
		'{print "  character(len=7), parameter :: operating_system = \""$$1"\""}' \
		>> operating_system.f90
	uname -m | awk \
		'{print "  character(len=7), parameter :: architecture = \""$$1"\""}' \
		>> operating_system.f90
	printf "end module os\n" >> operating_system.f90

mcfost2prodimo :  mcfost mcfost2prodimo.o
	$(FF) $(FFLAGS_EXE) -o mcfost2prodimo mcfost2prodimo.f90 $(LIBS_FITS)

#------------------------------------------------------------------------------#
# Release

.PHONY : release_utils fosti dmz wardlaw

release_doc:
	pushd . ; cd ../docs ; make html ; cd _build ; tar czf html.tgz html ; scp -r html.tgz $(IPAG_SERVER):/user/publicdir/pintec/mcfost/docs ; ssh $(IPAG_SERVER) "cd /user/publicdir/pintec/mcfost/docs ; tar xzf html.tgz ; rm -f html.tgz" ; popd

release_utils :
	tar czvf mcfost_utils.tgz -C../utils .
	openssl sha1 mcfost_utils.tgz | awk '{print $2}' > mcfost_utils.sha1
	scp -r mcfost_utils.tgz mcfost_utils.sha1 ../utils/Version $(IPAG_SERVER):/user/publicdir/pintec/mcfost_utils/
	rm -f mcfost_utils.tgz mcfost_utils.sha1

#------------------------------------------------------------------------------#
# Cleaning

.PHONY : clean cleandata cleandata_dust cleanall

clean :
	rm -rf *.o $(gdir)/*.o *.mod _dust_prop*.tmp _voronoi.tmp *.i90 operating_system.f90 system.o libmcfost.a
	if test $(MCFOST_GIT) = 1 ; then rm -rf sha.f90; fi

cleandata :
	rm -rf data_* seed_*

cleandata_dust :
	rm -rf _dust_prop*.tmp

cleanall : clean cleandata
	rm -rf mcfost

cleanGasmod :
	rm -rf $(gdir)/*.o $(gdir)/*.mod

#------------------------------------------------------------------------------#
# Testing

LASTSYSTEM = ${shell if [ -e .make_lastsystem ]; then cat .make_lastsystem; fi}

.PHONY : checkmcfost checksys err

checkmcfost:
   ifdef MCFOST_INSTALL
	@echo "MCFOST directory is "$(MCFOST_INSTALL);
   else
	@echo; echo "ERROR: cannot find MCFOST installation directory"
	@echo "- set this using MCFOST_INSTALL"; echo; ${MAKE} err
   endif

err:
	$(error aborting);

.make_lastsystem : checksys

checksys:
   ifeq ($(KNOWN_SYSTEM), yes)
	@echo ""
	@echo "Compiling MCFOST for $(SYSTEM) system..........."
	@echo ""
       ifneq ($(SYSTEM),$(LASTSYSTEM))
	    @echo system changed from ${LASTSYSTEM} to ${SYSTEM}
	    @${MAKE} clean
	    @echo $(SYSTEM) > .make_lastsystem
	    @${MAKE} mcfost_env.o
       endif
   else
	@echo ""
	@echo "make: ERROR: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo "=> set the environment variable SYSTEM to one listed "
	@echo "   in src/Makefile and try again"
	@echo ""
	@${MAKE} err;
endif
