#------------------------------------------------------------------------------#
# Makefile for Monte Carlo radiative transfer code MCFOST
#------------------------------------------------------------------------------#
# Linux compiling server
LINUX_SERVER= ozstar

KNOWN_SYSTEM= no

LIB= $(MCFOST_INSTALL)/lib

ifndef SYSTEM
    SYSTEM= ifort
endif

ifeq ($(gfortran), yes)
    SYSTEM=gfortran
endif

ifndef OPENMP
    OPENMP= yes
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

#------------------------------------------------------------------------------#
# System dependent compiler flags

ifeq ($(SYSTEM), ifort)
    ifeq ($(shell uname | tr '[a-z]' '[A-Z]' 2>&1 | grep -c DARWIN),1)
        ARCH= -axSSSE3,SSE4.1,SSE4.2,AVX,CORE-AVX2 -mmacosx-version-min=10.12 -mdynamic-no-pic -static-intel -qopenmp-link static
    else
        ARCH= -axSSE2,SSSE3,SSE4.1,SSE4.2,AVX,CORE-AVX2,CORE-AVX512 -static
    endif
    FC= ifort
    FFLAGS= -fpp -O3 -no-prec-div -fp-model fast=2 -traceback $(ARCH)
    DBLFLAG= -r8
    DEBUGFLAG= -check all -C -g -fpe0 -traceback -no-ftz \
               -warn uninitialized -warn unused -warn truncated_source
    KNOWN_SYSTEM= yes
    FOMPFLAGS= -fopenmp
    IPOFLAGS= -ipo
    CC= icc
    CCFLAGS= -O3 -no-prec-div -fp-model fast=2 -traceback $(ARCH)
    COMPFLAGS= -fopenmp
    LIBCXX= -cxxlib
    LIBS= $(LIB)/ifort
endif

ifeq ($(SYSTEM), xeonphi)
    FC= ifort
    FFLAGS= -fpp -O3 -no-prec-div -traceback -mmic -static-intel
    DBLFLAG= -r8
    DEBUGFLAG= -check all -C -g -fpe0 -traceback -no-ftz \
               -warn uninitialized -warn unused -warn truncated_source
    KNOWN_SYSTEM= yes
    FOMPFLAGS= -qopenmp
    IPOFLAGS= -ipo
    CC= icc
    CCFLAGS= -O3 -no-prev-div -traceback -mmic
    COMPFLAGS= -qopenmp
    LIBCXX= -cxxlib
    LIBS= $(LIB)/xeon-phi
endif

ifeq ($(SYSTEM), gfortran)
    FC= gfortran
    FFLAGS+= -m64 -cpp -O3 -funroll-loops -ffinite-math-only \
             -fno-trapping-math -ffree-line-length-none -fbacktrace  -static-libgcc -static-libstdc++
    DBLFLAG= -fdefault-real-8 -fdefault-double-8
    DEBUGFLAG= -g -fbounds-check -Wunused-variable -Wunused-function \
               -Wuninitialized -fbacktrace
    KNOWN_SYSTEM= yes
    FOMPFLAGS= -fopenmp
    CC= gcc
    CXX= g++
    CCFLAGS= -Wall -pedantic -m64 -O3 -funroll-loops -ffinite-math-only \
             -fno-trapping-math
    COMPFLAGS= -fopenmp
    LIBCXX= -lstdc++
    ifneq ($(shell $(CC) -v 2>&1 | grep -c "clang"), 0)
	LIBCXX= -lc++
    endif
    LIBS= $(LIB)/gfortran
endif

#------------------------------------------------------------------------------#
# Set other optional flags depending on system settings

ifeq ($(DEBUG), yes)
    FFLAGS+= ${DEBUGFLAG}
    FFLAGS:= $(FFLAGS:-O3=-O0)
endif

ifeq ($(DEV), yes)
   FFLAGS:= $(FFLAGS:-O3=-O1)
endif

ifeq ($(OPENMP), yes)
    FFLAGS+= ${FOMPFLAGS}
    CCFLAGS+= ${COMPFLAGS}
endif

ifeq ($(IPO), yes)
    FFLAGS+= ${IPOFLAGS}
    CCFLAGS+= ${IPOFLAGS}
endif

#------------------------------------------------------------------------------#
# Library flags

INCLUDE+= -I$(MCFOST_INSTALL)/include -I$(MCFOST_INSTALL)/include/voro++

LDFLAGS+= -L$(LIBS) $(LIBCXX) -lcfitsio -lvoro++ -lsprng

#------------------------------------------------------------------------------#
# Sources and objects

SOURCES = mcfost_env.f90 parameters.f90 constants.f90 sha.f90 messages.f90 operating_system.f90 \
	random_numbers.f90 utils.f90 fits_utils.f90 grains.f90 \
	cylindrical_grid.f90 spherical_grid.f90 kdtree2.f90 Voronoi.f90 grid.f90 \
	wavelengths.f90 stars.f90 read_DustEM.f90 Temperature.f90 density.f90 read_opacity.f90 \
	scattering.f90  coated_sphere.f90 dust_prop.f90 molecular_emission.f90 PAH.f90 \
	input.f90 benchmarks.f90 read_param.f90 dust_ray_tracing.f90 radiation_field.f90 \
	thermal_emission.f90 diffusion.f90 io_prodimo.f90 disk_physics.f90 output.f90 mem.f90 \
	optical_depth.f90 init_mcfost.f90 \
	io_phantom_infiles.f90 io_phantom_utils.f90 read_phantom.f90 read_gadget2.f90 SPH2mcfost.f90 \
	dust_transfer.f90 mol_transfer.f90

#------------------------------------------------------------------------------#
#Atomic line transfer
Atom_dir = Atoms
SOURCES_ATOM = $(Atom_dir)/math.f90 $(Atom_dir)/atom_type.f90 $(Atom_dir)/zeeman.f90 \
	$(Atom_dir)/getline.f90 $(Atom_dir)/constant.f90 $(Atom_dir)/voigtfunctions.f90 \
	$(Atom_dir)/uplow.f90 $(Atom_dir)/grid_type.f90 $(Atom_dir)/getlambda.f90 \
	$(Atom_dir)/spectrum_type.f90 $(Atom_dir)/barklem.f90 $(Atom_dir)/readatom.f90 \
	$(Atom_dir)/collision.f90 $(Atom_dir)/lte.f90 $(Atom_dir)/solvene.f90 \
	$(Atom_dir)/broad.f90  $(Atom_dir)/profilefunctions.f90 \
	$(Atom_dir)/rayleigh_scattering.f90 $(Atom_dir)/thomson_scattering.f90 \
	$(Atom_dir)/hydrogen.f90  $(Atom_dir)/opacity.f90  $(Atom_dir)/planck.f90 \
	$(Atom_dir)/metal.f90 $(Atom_dir)/project.f90 $(Atom_dir)/accelerate.f90

#$(Atom_dir)/atomictransfer.f90 $(Atom_dir)/readatmos.f90

SOURCES+= $(SOURCES_ATOM)
OBJECTS= $(SOURCES:.f90=.o)

#------------------------------------------------------------------------------#
# Compile
main : mcfost

#include .depend

all : mcfost libmcfost.a

mcfost :  $(OBJECTS) voro++_wrapper.o mcfost.o
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) voro++_wrapper.o mcfost.o $(LDFLAGS)

libmcfost.a : $(OBJECTS) voro++_wrapper.o mcfost2phantom.o
	rm -rf libtmp
	mkdir -p libtmp
	cp $(LIBS)/libcfitsio.a libtmp
	cp $(LIBS)/libsprng.a libtmp
	cp $(LIBS)/libvoro++.a libtmp
	cd libtmp; \
	ar x libcfitsio.a; \
	ar x libsprng.a; \
	ar x libvoro++.a
	rm -f mcfost.o io_phantom_utils.o io_phantom_infiles.o libtmp/pmlcg.o \
		libmcfost.a
	ar vr libmcfost.a *.o libtmp/*.o
	ranlib libmcfost.a

mcfost_env.o : .make_lastsystem mcfost_env.f90
	rm -rf ../utils/Version ; grep required_utils_version mcfost_env.f90 | awk \
		'{print $$6}' > ../utils/Version
	${MAKE} checkmcfost
	$(FC) $(FFLAGS) -c -DPOINTER_SIZE=8 mcfost_env.f90 $(INCLUDE)

io_phantom_utils.o :  io_phantom_utils.f90
	$(FC) $(FFLAGS) -c $(DBLFLAG) io_phantom_utils.f90

voro++_wrapper.o :  voro++_wrapper.cpp
	$(CC) $(CCFLAGS) -c voro++_wrapper.cpp $(INCLUDE)

%.o : %.f90
	$(FC) $(FFLAGS) -c -DPOINTER_SIZE=8 $< -o $@ $(INCLUDE)

io_prodimo.o : sha.f90

utils.o : sha.f90

Voronoi.o : voro++_wrapper.o

ifeq ($(MCFOST_GIT),1)
sha.f90 : ../.git/objects mcfost_env.o
	printf "module sha\n  implicit none\n" > sha.f90
	git log --pretty=format:"%H %d" | grep \(HEAD | awk \
		'{print "  character(len=40), parameter :: sha_id = \""$$1"\""}' \
		>> sha.f90
	printf "end module sha\n" >> sha.f90
endif

operating_system.f90 : sha.o
	printf "module os\n  implicit none\n" > operating_system.f90
	uname | awk \
		'{print "  character(len=7), parameter :: operating_system = \""$$1"\""}' \
		>> operating_system.f90
	printf "end module os\n" >> operating_system.f90

mcfost2prodimo :  mcfost mcfost2prodimo.o
	$(FF) $(FFLAGS_EXE) -o mcfost2prodimo mcfost2prodimo.f90 $(LIBS_FITS)

#------------------------------------------------------------------------------#
# Release

.PHONY : release release_mac release_utils fosti dmz wardlaw

release : main
	- git tag -d release
	git tag release
	git push github
	\cp -f mcfost $(MCFOST_INSTALL)/bin/
	rm -rf web
	mkdir web web/linux web/macos
	cp mcfost mcfost_update
	tar czvf mcfost_bin.tgz mcfost_update ; rm -rf mcfost_update
	openssl sha1 mcfost_bin.tgz | awk '{print $$2}' > web/macos/mcfost.sha1
	mv mcfost_bin.tgz web/macos/
	ssh $(LINUX_SERVER) "cd mcfost/src ; git pull ; git checkout release ; export MCFOST_GIT=1 ; export MAKEFLAGS= ; make clean mcfost ; \cp -f mcfost ~/mcfost_cigri ; \cp -f mcfost mcfost_update ; rm -f mcfost.sha1 ; tar czf mcfost_bin.tgz mcfost_update --remove-files ; sha1sum mcfost_bin.tgz  | awk '{print $$1}' > mcfost.sha1"
	scp $(LINUX_SERVER):mcfost/src/mcfost_bin.tgz web/linux/
	scp $(LINUX_SERVER):mcfost/src/mcfost.sha1 web/linux/
	cp `grep mcfost_version mcfost_env.f90 | awk '{print "ref"$$6".para ref"$$6"_multi.para ref"$$6"_3D.para"}' ` web
	grep mcfost_release mcfost_env.f90 | awk '{print $$6}' | sed s/\"//g > web/version.txt
	git log --pretty=format:"%cd  %h %d %s" --date=short  --graph | sed s/"dmz\/master, phare\/master, github\/master, bitbucket\/master"// | sed s/"()"// | sed s/", ,"/","/ > web/history.txt
	mkdir -p yorick ; cp ~/yorick/init/mcfost_struct.i ~/yorick/init/mcfost_utils.i yorick/
	cp -r yorick web/
	rsync -Pur  web/* horus:/user/publicdir/pintec/mcfost
	ssh horus "cd  /user/publicdir/pintec/mcfost/macos ; tar xzf mcfost_bin.tgz ; \cp -f mcfost_update mcfost ; tar czf mcfost.tgz mcfost --remove-files ; rm -f mcfost_update ; cd  /user/publicdir/pintec/mcfost/linux ; tar xzf mcfost_bin.tgz ; \cp -f mcfost_update mcfost ; tar czf mcfost.tgz mcfost --remove-files ; rm -f mcfost_update"
	rm -rf yorick web
	pushd . ; cd ../docs ; make html ; cd _build ; tar czf html.tgz html ; scp -r html.tgz horus:/user/publicdir/pintec/mcfost/docs ; ssh horus "cd /user/publicdir/pintec/mcfost/docs ; tar xzf html.tgz ; rm -f html.tgz" ; popd

release_mac : main
	\cp -f mcfost $(MCFOST_INSTALL)/bin/
	rm -rf web
	mkdir web web/macos
	cp mcfost mcfost_update
	tar czvf mcfost_bin.tgz mcfost_update ; rm -rf mcfost_update
	openssl sha1 mcfost_bin.tgz | awk '{print $$2}' > web/macos/mcfost.sha1
	mv mcfost_bin.tgz web/macos/
	cp `grep mcfost_version mcfost_env.f90 | awk '{print "ref"$$6".para ref"$$6"_multi.para ref"$$6"_3D.para"}' ` web
	grep mcfost_release mcfost_env.f90 | awk '{print $$6}' | sed s/\"//g > web/version.txt
	git log --pretty=format:"%cd  %h  %d %s" --date=short --graph --branches > web/history.txt
	mkdir -p yorick ; cp ~/yorick/init/mcfost_struct.i ~/yorick/init/mcfost_utils.i yorick/
	cp -r yorick web/
	rsync -Pur -r web/* horus:/user/publicdir/pintec/mcfost
	rm -rf yorick web

release_utils :
	tar czvf mcfost_utils.tgz -C../utils .
	openssl sha1 mcfost_utils.tgz | awk '{print $2}' > mcfost_utils.sha1
	scp -r mcfost_utils.tgz mcfost_utils.sha1 ../utils/Version horus:/user/publicdir/pintec/mcfost_utils/
	rm -f mcfost_utils.tgz mcfost_utils.sha1

fosti :
	rsync -Pur *.f90 Makefile fosti:mcfost/src
	ssh fosti "cd mcfost/src ; make clean mcfost fostino=yes; \cp -f mcfost ~/mcfost_cigri/"

dmz :
	rsync -Pur *.f90 Makefile dmz:mcfost/src
	ssh dmz "cd mcfost/src ; make clean mcfost"

wardlaw :
	rsync -Pur *.f90 Makefile wardlaw:mcfost/src
	ssh wardlaw "source .login ; cd mcfost/src ; make clean mcfost ; \cp -f mcfost ~/bin/"

#------------------------------------------------------------------------------#
# Cleaning

.PHONY : clean cleandata cleandata_dust cleanall

clean :
	rm -rf *.o $(Atom_dir)/*.o *.mod _dust_prop*.tmp _voronoi.tmp *.i90 operating_system.f90 system.o libmcfost.a
	if test $(MCFOST_GIT) = 1 ; then rm -rf sha.f90; fi

cleandata :
	rm -rf data* seed*

cleandata_dust :
	rm -rf _dust_prop*.tmp

cleanall : clean cleandata
	rm -rf mcfost

#------------------------------------------------------------------------------#
# Backup

.PHONY : save

save :
	#git push phare
	git push bitbucket --mirror
	git push github --mirror
	#git push dmz

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
       endif
   else
	@echo ""
	@echo "make: ERROR: value of SYSTEM = $(SYSTEM) not recognised..."
	@echo "=> set the environment variable SYSTEM to one listed "
	@echo "   in src/Makefile and try again"
	@echo ""
	@${MAKE} err;
   endif
