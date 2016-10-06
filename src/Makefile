#*************************************************************
# Makefile code transfert radiatif Monte Carlo MCFOST
# Multi-systemes
# Versions sequentielle et parallele (OpenMP)
#*************************************************************
# 1er mac est en SSSE3
# nouveau mac en SSE4.1 (SSE4.2 pour core i7)
# -mmacosx-version-min=10.4 requis sur Snow Leopard

# Choix du type de compilateur
ARCH = $(shell uname | tr '[a-z]' '[A-Z]')
PLATFORM =  $(shell uname -p)

# Default is ifort
COMPILO = ifort
ifeq ($(g95),yes)
COMPILO = g95
endif
ifeq ($(gfortran),yes)
COMPILO = gfortran
endif
ifeq ($(gfortran32),yes)
COMPILO = gfortran32
endif

BINARY = mcfost

FFLAGS_R8 = -r8

# Linux ou Mac x86, ifort 64 bits
ifeq ($(COMPILO),ifort)
FC_SER = ifort -fpp -traceback
FC_OMP = ifort -fpp -openmp -traceback

CC_SER = icc -traceback
CC_OMP = icc -openmp -traceback

OPTS = -mdynamic-no-pic -O3 -no-prec-div #-fast (because -fast implies -xHOST and -ipo)
OPTS = -O3 -no-prec-div
ifeq ($(ipo),yes)
OPTS2 = $(OPTS) -ipo
OPTS = $(OPTS2)
endif
DEBUG_OPTS = -check all  -C -g  -fpe0 -traceback -no-ftz  -warn uninitialized -warn unused -warn truncated_source # -no-ftz permet eviter bug dans conversion db -> sl

PREC = #
FFLAGS =-c $(PREC)

ifeq ($(ARCH),DARWIN)
MACH =  -axSSSE3,SSE4.1,SSE4.2  #-tpp7 #-xW -tune pn4 # - arch i386 ou x86_64
FFLAGS_EXE = $(PREC) -static-intel -mmacosx-version-min=10.6 #-static-intel #-static-libgcc #problem -static-intel avec Xcode 3.2.2
endif

ifeq ($(ARCH),LINUX)
MACH =  -axSSE2,SSSE3,SSE4.1,SSE4.2,AVX # for froggy, dmz and ifort >=13
ifeq ($(fostino),yes)
MACH = -xW -axP # for fostino and ifort 10
endif
FFLAGS_EXE = $(PREC)  -static #-static-libcxa -i-static
endif

LIBS = $(MCFOST_INSTALL)/lib/ifort/

# Modification de quelques variables pour les xeon phi
ifeq ($(phi),yes)
LIBS = $(MCFOST_INSTALL)/lib/xeon-phi/
OPTS = -O3 -no-prec-div
MACH = -mmic
FFLAGS_EXE = $(PREC)  -static-intel
BINARY = mcfost.phi
endif

LIBS_FITS = -L $(LIBS) -lcfitsio
LIBS_VORO =  -L$(LIBS) -lvoro++ -cxxlib

# SPRNG
FLAGS_SPRNG = $(FLAGS) -DPOINTER_SIZE=8
endif

# Linux ou Mac x86, g95 64 bits
ifeq ($(COMPILO),g95)
FC_SER = g95 -cpp
FC_OMP = g95 -cpp
OPTS = -O3 -ffast-math #-mtune=athlon -mcpu=athlon -mfpmath=sse -msse
DEBUG_OPTS =  -g -fbounds-check -ftrace=full -fimplicit-none #-g -C -ggdb -ftrace=full
# -Wall -Wextra -Wimplicit-none -Werror=113,115,137,146,147,159,163 -ftrace=full -fbounds-check -freal=nan -finteger=-2109876543 -fmodule-private
MACH = #-march=athlon
PREC =
FFLAGS = -c $(PREC)
FFLAGS_EXE = $(PREC)

LIBS = $(MCFOST_INSTALL)/lib/g95/

# SPRNG
FLAGS_SPRNG = $(FLAGS) -DPOINTER_SIZE=8
endif

# Linux ou Mac x86, gfortran 32 bits
ifeq ($(COMPILO),gfortran32)
FC_SER = gfortran -m32 -cpp
FC_OMP = gfortran -m32 -cpp -fopenmp
OPTS = -O3 -ffast-math  -funroll-loops
FFLAGS_R8 = -fdefault-real-8 -fdefault-double-8
#-mtune=athlon -mcpu=athlon -mfpmath=sse -msse
DEBUG_OPTS =  -g -fbounds-check -Wunused-variable -Wunused-function -Wuninitialized -Wconversion -Warray-temporaries
#-g -C -ggdb -ftrace=full  -fimplicit-none
# -Wall -Wextra -Wimplicit-none -Werror=113,115,137,146,147,159,163 -ftrace=full -fbounds-check -freal=nan -finteger=-2109876543 -fmodule-private
MACH = #-march=athlon
PREC =
FFLAGS = -c $(PREC)
FFLAGS_EXE = $(PREC)

LIBS = $(MCFOST_INSTALL)/lib/gfortran32/

# SPRNG
FLAGS_SPRNG = $(FLAGS)
endif


# Linux ou Mac x86, gfortran64
ifeq ($(COMPILO),gfortran)
FC_SER = gfortran -m64 -cpp
FC_OMP = gfortran -m64 -cpp -fopenmp
CC_SER = g++ -Wall -ansi -pedantic -m64
CC_OMP = g++ -Wall -ansi -pedantic -m64
OPTS = -O3 -funroll-loops -ffinite-math-only -fno-trapping-math #-ffast-math : -funsafe-math-optimizations make mcfost stall, rest not much faster than 03
#-mtune=athlon -mcpu=athlon -mfpmath=sse -msse
DEBUG_OPTS =  -g -fbounds-check -Wunused-variable -Wunused-function -Wuninitialized #-Wconversion -Warray-temporaries
#-Wunused-variable -Wunused-function -Wuninitialized -Wconversion -Warray-temporaries
#-g -C -ggdb -ftrace=full  -fimplicit-none
# -Wall -Wextra -Wimplicit-none -Werror=113,115,137,146,147,159,163 -ftrace=full -fbounds-check -freal=nan -finteger=-2109876543 -fmodule-private
FFLAGS_R8 = -fdefault-real-8 -fdefault-double-8

MACH = #-march=athlon
PREC =
FFLAGS = -c $(PREC)
FFLAGS_EXE = $(PREC)

LIBS = $(MCFOST_INSTALL)/lib/gfortran64/
LIBS_VORO =  -L$(LIBS) -lvoro++ -lstdc++

FLAGS_SPRNG = $(FLAGS) -DPOINTER_SIZE=8
endif

LIBS_FITS = -L$(LIBS) -lcfitsio

# Numerical Recipes
MODS_NR = -I$(LIBS)/nr
LIBS_NR =  -L$(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L$(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE =  -L$(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L$(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
INCLUDE_SPRNG = -I$(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L$(LIBS) -lsprng

#Voro++
INCLUDE_VORO = -I$(MCFOST_INSTALL)/include/voro++

FF = $(FC_OMP)  $(MACH) $(OPTS)
CC = $(CC_OMP)  $(MACH) $(OPTS)
ifeq ($(debugp), yes)
FF = $(FC_OMP)  $(DEBUG_OPTS)
CC = $(CC_OMP)
endif
ifeq ($(debug), yes)
FF = $(FC_SER) $(DEBUG_OPTS)
CC = $(CC_SER)
endif
ifeq ($(mono), yes)
FF = $(FC_SER)  $(MACH) $(OPTS)
CC = $(CC_SER)  $(MACH) $(OPTS)
endif
ifeq ($(prof), yes)
FF = $(FC_SER) -O0 -g
CC = $(CC_SER) -O0 -g
endif
ifeq ($(dev), yes)
FF = $(FC_SER) -O1
CC = $(CC_SER) -O1
endif
ifeq ($(devp), yes)
FF = $(FC_OMP) -O1
CC = $(CC_OMP) -O1
endif


SOURCES = modules.f90 sha.f90 system.f90 utils.f90 mem.f90 \
	cylindrical_grid.f90 spherical_grid.f90 Voronoi.f90 grid.f90 \
	input.f90 output.f90 benchmarks.f90 read_param.f90 io_prodimo.f90 \
	init_mcfost.f90 scattering.f90 coated_sphere.f90 dust_ray_tracing.f90 \
	dust.f90 molecules.f90  optical_depth.f90 disk_physics.f90 \
	stars.f90 density.f90 PAH.f90 diffusion.f90 thermal_emission.f90 \
	io_phantom_infiles.f90 io_phantom_utils.f90 read_phantom.f90 \
	read_gadget2.f90 SPH2mcfost.f90 \
	dust_transfer.f90 mol_transfer.f90 \

CC_SOURCES = voro++_wrapper.cpp

OBJECTS = $(SOURCES:.f90=.o)

main : mcfost

#*************************************************************
# Code séquentiel et parallèle
#*************************************************************
mcfost : $(OBJECTS) voro++_wrapper.o mcfost.o
	$(FF) $(FFLAGS_EXE) -o $(BINARY)  $(OBJECTS) voro++_wrapper.o mcfost.o   $(LIBS_NR_EQDIFF) $(LIBS_SORT) $(LIBS_NR) $(LIBS_SPRNG) $(LIBS_FITS) $(LIBS_VORO) $(INCLUDE_SPRNG) $(INCLUDE_VORO)

libmcfost.a : checkmcfost $(OBJECTS) voro++_wrapper.o mcfost2phantom.o
	mkdir ../lib/artmp
	cp ../lib/lib/libcfitsio.a ../lib/lib/libsprng.a ../lib/lib/libvoro++.a ../lib/lib/nr/libnr.a ../lib/lib/nr/eq_diff/libnr_eq_diff.a ../lib/lib/nr/sort/libnr_sort.a ../lib/artmp
	ar x ../lib/artmp/libcfitsio.a
	ar x ../lib/artmp/libnr.a
	ar x ../lib/artmp/libnr_eq_diff.a
	ar x ../lib/artmp/libnr_sort.a
	ar x ../lib/artmp/libsprng.a
	ar x ../lib/artmp/libvoro++.a
	ar r ../lib/artmp/libmcfost.a *.o
	mv ../lib/artmp/libmcfost.a .
	rm -r ../lib/artmp
	ranlib libmcfost.a

# Fichiers necessaires aux suivants
modules.o : modules.f90
	rm -rf ../utils/Version ; grep required_utils_version modules.f90 | awk '{print $$6}' > ../utils/Version
	${MAKE} checkmcfost
	$(FF) $(FFLAGS) $(FLAGS_SPRNG) modules.f90 $(INCLUDE_SPRNG)

# Fichier compile sans optimisation
mem.o : mem.f90  modules.o
	$(FF) $(FFLAGS) $(FLAGS_SPRNG) mem.f90 $(INCLUDE_SPRNG)

io_phantom_utils.o : io_phantom_utils.f90
	$(FF) $(FFLAGS) $(FFLAGS_R8) io_phantom_utils.f90

voro++_wrapper.o : voro++_wrapper.cpp
	$(CC) $(FFLAGS) voro++_wrapper.cpp $(INCLUDE_VORO)

# defaut
%.o : %.f90  modules.o
	$(FF) $(FFLAGS) $(FLAGS_SPRNG) $< -o $@ $(MODS_NR) $(INCLUDE_SPRNG)



io_prodimo.o : sha.f90
utils.o :      sha.f90

# Le fichier sha.f90 n'est cree que sur mon mac
ifeq ($(MCFOST_GIT),1)
sha.f90 : ../.git/objects modules.o
	printf "module sha\n  implicit none\n" > sha.f90
	git log --pretty=format:"%H %d" | grep \(HEAD | awk '{print "  character(len=40), parameter :: sha_id = \""$$1"\""}' >> sha.f90
	printf "end module sha\n" >> sha.f90
endif

system.f90 : sha.o
	printf "module system\n  implicit none\n" > system.f90
	uname | awk '{print "  character(len=7), parameter :: operating_system = \""$$1"\""}' >> system.f90
	printf "end module system\n" >> system.f90

mcfost2prodimo : mcfost mcfost2prodimo.o
	$(FF) $(FFLAGS_EXE) -o mcfost2prodimo mcfost2prodimo.f90 $(LIBS_FITS)

#voronoi : Voronoi.o
#	$(FF) $(FFLAGS_EXE) -o voronoi Voronoi.o utils.o modules.o $(LIBS_SPRNG) $(INCLUDE_SPRNG)

#*************************************************************
# Release
#*************************************************************

release : main
	- git tag -d release
	\cp -f mcfost $(MCFOST_INSTALL)/bin/
	rm -rf web
	mkdir web web/linux web/macos
	cp mcfost mcfost_update
	tar czvf mcfost_bin.tgz mcfost_update ; rm -rf mcfost_update
	openssl sha1 mcfost_bin.tgz | awk '{print $$2}' > web/macos/mcfost.sha1
	mv mcfost_bin.tgz web/macos/
	rsync -Pu *.f90 Makefile dmz:mcfost/src
	ssh dmz "cd mcfost/src ; make clean mcfost ; \cp -f mcfost ~/mcfost_cigri ; \cp -f mcfost mcfost_update ; rm -f mcfost.sha1 ; tar czf mcfost_bin.tgz mcfost_update --remove-files ; sha1sum mcfost_bin.tgz  | awk '{print $$1}' > mcfost.sha1"
	scp dmz:mcfost/src/mcfost_bin.tgz web/linux/
	scp dmz:mcfost/src/mcfost.sha1 web/linux/
	cp `grep mcfost_version modules.f90 | awk '{print "ref"$$6".para ref"$$6"_multi.para ref"$$6"_3D.para"}' ` web
	grep mcfost_release modules.f90 | awk '{print $$6}' | sed s/\"//g > web/version.txt
	git log --pretty=format:"%cd  %h %d %s" --date=short  --graph | sed s/"riri\/master, phare\/master, github\/master, dropbox\/master, bitbucket\/master"// | sed s/"()"// | sed s/", ,"/","/ > web/history.txt
	cp -r Doc web/
	openssl sha1 web/Doc/MCFOSTManual.pdf | awk '{print $$2}' > web/Doc/MCFOSTManual.sha1
	mkdir -p yorick ; cp ~/yorick/init/mcfost_struct.i ~/yorick/init/mcfost_utils.i yorick/
	cp -r yorick web/
	rsync -Pur  web/* horus:/user/publicdir/pintec/mcfost
	rm -rf yorick web
	git tag release


release_mac : main
	\cp -f mcfost $(MCFOST_INSTALL)/bin/
	rm -rf web
	mkdir web web/macos
	cp mcfost mcfost_update
	tar czvf mcfost_bin.tgz mcfost_update ; rm -rf mcfost_update
	openssl sha1 mcfost_bin.tgz | awk '{print $$2}' > web/macos/mcfost.sha1
	mv mcfost_bin.tgz web/macos/
	cp `grep mcfost_version modules.f90 | awk '{print "ref"$$6".para ref"$$6"_multi.para ref"$$6"_3D.para"}' ` web
	grep mcfost_release modules.f90 | awk '{print $$6}' | sed s/\"//g > web/version.txt
	git log --pretty=format:"%cd  %h  %d %s" --date=short --graph --branches > web/history.txt
	cp -r Doc web/
	openssl sha1 web/Doc/MCFOSTManual.pdf | awk '{print $$2}' > web/Doc/MCFOSTManual.sha1
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

#*************************************************************
# Nettoyage
#*************************************************************
clean :
	rm -rf *.o  *.mod _dust_prop*.tmp *.i90 system.f90 system.o
	if test $(MCFOST_GIT) ; then rm -rf sha.f90; fi

cleandata :
	rm -rf data* seed*

cleandata_dust :
	rm -rf _dust_prop*.tmp

cleanall : clean cleandata
	rm -rf mcfost

#*************************************************************
# Sauvergarde
#*************************************************************
save :
	git push phare
	git push bitbucket
	git push github
	git push riri

#*************************************************************
# test
#*************************************************************
checkmcfost:
   ifdef MCFOST_INSTALL
	@echo "MCFOST directory is "$(MCFOST_INSTALL);
   else
	@echo; echo "ERROR: cannot find MCFOST installation directory - set this using MCFOST_INSTALL"; echo; ${MAKE} err
   endif

err:
	$(error aborting);
