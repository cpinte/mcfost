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
ifeq ($(g95),yes)
COMPILO = g95
endif
ifeq ($(gfortran),yes)
COMPILO = gfortran
endif
ifeq ($(sun),yes)
COMPILO = sunf95
endif
ifeq ($(ifort32),yes)
COMPILO = ifort32
endif

# Linux, ifort
ifeq ($(ARCH),LINUX)
PLATFORM =  $(shell uname -m)
ifeq ($(PLATFORM),i686)
FC_SER = ifort -fpp  -traceback
FC_OMP = ifort -fpp  -openmp -traceback
OPTS = -fast  #-ipo -no-prec-div -fast
# fpe0 bug avec openmp
DEBUG_OPTS =  -check all -traceback  -fpstkchk  -C -g -d2  -fpe0
#PARALLEL = -openmp
MACH = #-tpp7 #-xW -tune pn4
PREC = #-pc80 -prec_div -mp
# Option de conservation de la precision
# -mp : maintain precision fait prendre 1 facteur 2 dans le tps d'execution
# mais sinon donne des resultats faux
# -pc80 -prec_div ne fait pas perdre plus de tps
FFLAGS = -c  $(PREC) 
FFLAGS_EXE = $(PREC)  -static #-static-libcxa -i-static

LIBS = $(MCFOST_INSTALL)/lib/ifort/
LIBS_PGPLOT = -L $(LIBS) -lpgplot -L /usr/X11R6/lib/ -lX11 -static
LIBS_FITS = -L $(LIBS) -lcfitsio
LIBS_SYST = -L $(LIBS) -lappel_syst -Vaxlib

# Numerical recipes
MODS_NR = -I $(LIBS)/nr 
LIBS_NR = -L $(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L $(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE = -L $(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG = $(FLAGS) 
INCLUDE_SPRNG = -I $(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L $(LIBS) -lsprng

# Linux, sunf95
ifeq ($(COMPILO),sunf95)
FC_SER = sunf95 -w -fpp
FC_OMP = sunf95 -w -fpp  -xopenmp -O3
OPTS =  -O3  #-mtune=athlon -mcpu=athlon -mfpmath=sse -msse 
DEBUG_OPTS = -g -C -ggdb
#PARALLEL = -openmp
MACH = #-march=athlon
PREC = 
FFLAGS = -c $(PREC) 
FFLAGS_EXE = $(PREC) 

LIBS = $(MCFOST_INSTALL)/lib/sunf95/
LIBS_PGPLOT = -L $(LIBS) -lpgplot
LIBS_FITS = -L $(LIBS) -lcfitsio
LIBS_SYST = -L $(LIBS) -lappel_syst

# Numerical Recipes
MODS_NR = -M$(LIBS)/nr
LIBS_NR =  -L $(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L $(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE =  -L $(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG = $(FLAGS) 
INCLUDE_SPRNG = -I$(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L $(LIBS) -lsprng
endif

endif

ifeq ($(PLATFORM),x86_64)
FC_SER = ifort -fpp -traceback
FC_OMP = ifort -fpp  -openmp -traceback
OPTS = -O3 -ipo -no-prec-div #-fast (because -fast implies -xP, cannot be overwritten)
# fpe0 bug avec openmp
DEBUG_OPTS =  -check all -traceback  -fpstkchk  -C -g -d2  -fpe0
#PARALLEL = -openmp
MACH =  -xW -axP #-tpp7 #-xW -tune pn4
PREC = #-pc80 -prec_div -mp
# Option de conservation de la precision
# -mp : maintain precision fait prendre 1 facteur 2 dans le tps d'execution
# mais sinon donne des resultats faux
# -pc80 -prec_div ne fait pas perdre plus de tps
FFLAGS = -c  $(PREC) 
FFLAGS_EXE = $(PREC)  -static #-static-libcxa -i-static

LIBS = $(MCFOST_INSTALL)/lib/ifort/
LIBS_PGPLOT = -L $(LIBS) -lpgplot -L /usr/X11R6/lib/ -lX11 -static
LIBS_FITS = -L $(LIBS) -lcfitsio
LIBS_SYST = -L $(LIBS) -lappel_syst -Vaxlib

# Numerical recipes
MODS_NR = -I $(LIBS)/nr 
LIBS_NR = -L $(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L $(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE = -L $(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG = $(FLAGS) -DPOINTER_SIZE=8
INCLUDE_SPRNG = -I $(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L $(LIBS) -lsprng
endif


ifeq ($(PLATFORM),ia64)
FC_SER = ifort -fpp 
FC_OMP = ifort -fpp  -openmp
OPTS = -O3 #-ipo -no-prec-div 
DEBUG_OPTS =  -check all -traceback  -fpstkchk  -C -g -d2  -fpe0
MACH = -tpp2 #-xW -tune pn4
PREC = #-pc80 -prec_div -mp
FFLAGS = -c  $(PREC) 
FFLAGS_EXE = $(PREC)  -static-libcxa -i-static

LIBS = $(MCFOST_INSTALL)/lib/ifort/
LIBS_PGPLOT = -L $(LIBS) -lpgplot -L /usr/X11R6/lib/ -lX11 -static
LIBS_FITS = -L $(LIBS) -lcfitsio
LIBS_SYST = -L $(LIBS) -lappel_syst -Vaxlib

# Numerical recipes
MODS_NR = -I $(LIBS)/nr 
LIBS_NR = -L $(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L $(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE = -L $(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG = $(FLAGS) -DPOINTER_SIZE=8 # 64 bits Itanium
INCLUDE_SPRNG = -I $(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L $(LIBS) -lsprng
endif

# Linux ppc64
ifeq ($(PLATFORM),ppc64)
FC_SER = xlf90 -qsuffix=f=f90
FC_OMP = xlf90_r -qsuffix=f=f90 -qsmp=omp -qnosave #-qreport=smplist
OPTS =   -O3 -qstrict -qtune=ppc970 -qarch=ppc970 -qcache=auto
#-qhot fait des optimisations sur les boucles -> decalage des indices dans la sousroutine opacite
DEBUG_OPTS = -g -qsigtrap -qcheck  -qnooptimize #-qflttrap=ov:und:zero:inv:imp:en
MACH = -qarch=auto -qtune=auto
FFLAGS = -c 
FFLAGS_EXE =  

LIBS = $(MCFOST_INSTALL)/lib/
LIBS_PGPLOT = -L $(LIBS) -lpgplot -lX11
LIBS_FITS =  -L $(LIBS) -lcfitsio 
COMP_FITS = -qextname=ftclos:ftcrhd:ftfiou:ftgerr:ftgiou:ftgmsg:ftinit:ftphpr:ftpkye:ftppre:ftgknj:ftgpve:ftopen:ftpprd:ftpprj:ftgpvj:ftmahd:ftpkys:ftpkyj
LIBS_SYST =  -L $(LIBS) -lappel_syst

# Numerical Recipes
MODS_NR = -I $(LIBS)/nr
LIBS_NR = -L $(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L $(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE = -L $(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG = -qsuffix=f=f90:cpp=f90
INCLUDE_SPRNG = -I $(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L $(LIBS) -lsprng
#profiling -pg  puis utiliser gprof
endif

endif


# AIX , xlf90 , Power3 (hal,...)
# Marche avec GNU-make
ifeq ($(ARCH),AIX)
FC_SER = xlf90 -qsuffix=f=f90
FC_OMP = xlf90_r -qsuffix=f=f90 -qsmp=omp -qnosave #-qreport=smplist
OPTS =   -O3 -qstrict -lmass #-qsclk=micro
#-qhot fait des optimisations sur les boucles -> decalage des indices dans la sousroutine opacite
DEBUG_OPTS = -g -qsigtrap -qcheck  -qnooptimize #-qflttrap=ov:und:zero:inv:imp:en
MACH = -qarch=pwr3 -qtune=pwr3
FFLAGS = -c 
FFLAGS_EXE =  -bmaxdata:0x80000000 -bmaxstack:0x20000000

LIBS = $(MCFOST_INSTALL)/lib/
LIBS_PGPLOT = -L $(LIBS) -lpgplot -lX11
LIBS_FITS =  -L $(LIBS) -lcfitsio
LIBS_SYST =  -L $(LIBS) -lappel_syst

# Numerical Recipes
MODS_NR = -I $(LIBS)/nr
LIBS_NR = -L $(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L $(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE = -L $(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG = -qsuffix=f=f90:cpp=f90
INCLUDE_SPRNG = -I $(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L $(LIBS) -lsprng
#profiling -pg  puis utiliser gprof
endif

# SUN, AMD64
# -DPOINTER_SIZE=8 obligatoire pour sprng 
# a cause conversion db -> sl peut renvoyer rn = 1.0
# A priori -f pour cfitsio sur sun (mais pas besoin)
# il faut tout compiler avec  -xarch=amd64 (y compris code C)
ifeq ($(ARCH),SUNOS)
FC_SER = f95 -fpp  
FC_OMP = f95 -fpp  -xopenmp
OPTS = -fast 
# fpe0 bug avec openmp
DEBUG_OPTS = -C -g -ftrap=%all
#PARALLEL = -openmp
MACH =  -xarch=amd64
PREC = 
# Option de conservation de la precision
# -mp : maintain precision fait prendre 1 facteur 2 dans le tps d'execution
# mais sinon donne des resultats faux
# -pc80 -prec_div ne fait pas perdre plus de tps
FFLAGS = -c  $(PREC) 
FFLAGS_EXE = $(PREC) 

LIBS = $(MCFOST_INSTALL)/lib/
LIBS_PGPLOT = -L $(LIBS) -lpgplot -L /usr/X11R6/lib/ -lX11
LIBS_FITS = -L $(LIBS) -lcfitsio -lnsl -lsocket
LIBS_SYST = -L $(LIBS) -lappel_syst

# Numerical recipes
MODS_NR = -M$(LIBS)/nr 
LIBS_NR = -L$(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L$(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE = -L$(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
#FLAGS_SPRNG = $(FLAGS) 
FLAGS_SPRNG = $(FLAGS) -DPOINTER_SIZE=8 # 64 bits 
INCLUDE_SPRNG = -I$(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L$(LIBS) -lsprng
endif

# Les macs
ifeq ($(ARCH),DARWIN)
PLATFORM =  $(shell uname -m)

ifeq ($(PLATFORM),i386)
FC_SER = ifort -fpp -traceback
FC_OMP = ifort -fpp -openmp -traceback #-stand f03 #-diag-enable sv3
# fast = -O3 -ipo -no-prec-div -mdynamic-no-pic
OPTS = -O3 #-fp-model fast=2 -no-prec-div -no-prec-sqrt -mdynamic-no-pic -ftz -auto -ipo-jobs2
# traceback does not slow down the calculations
# fpe0 bug avec openmp
DEBUG_OPTS = -check all  -C -g  -fpe0  -traceback -no-ftz # -no-ftz permet eviter bug dans conversion db -> sl 
#PARALLEL = -openmp
MACH =  -axSSSE3,SSE4.1,SSE4.2  #-tpp7 #-xW -tune pn4 # - arch i386 ou x86_64
PREC = #-pc80 -prec_div -mp
# Option de conservation de la precision
# -mp : maintain precision fait prendre 1 facteur 2 dans le tps d'execution
# mais sinon donne des resultats faux
# -pc80 -prec_div ne fait pas perdre plus de tps
FFLAGS =-c $(PREC)
FFLAGS_EXE = $(PREC) -static-intel -mmacosx-version-min=10.4 #-static-intel #-static-libgcc #problem -static-intel avec Xcode 3.2.2

LIBS = $(MCFOST_INSTALL)/lib/ifort/
LIBS_PGPLOT = -L $(LIBS) -lpgplot -L /usr/X11R6/lib/ -lX11 -static
LIBS_FITS = -L $(LIBS) -lcfitsio
LIBS_SYST = -L $(LIBS) -lappel_syst 

# Numerical recipes
MODS_NR = -I $(LIBS)/nr 
LIBS_NR = -L $(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L $(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE = -L $(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG =$(FLAGS) -DPOINTER_SIZE=8
INCLUDE_SPRNG = -I $(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L $(LIBS) -lsprng
endif


ifeq ($(PLATFORM),Power Macintosh)
NAME = $(shell uname -n)

# g95 par defaut
FC_SER = g95 -cpp
FC_OMP = g95 -cpp
OPTS = -O3 -ffast-math #-mtune=athlon -mcpu=athlon -mfpmath=sse -msse 
DEBUG_OPTS = -g -C -ggdb -ftrace=full
#PARALLEL = -openmp
MACH = #-march=athlon
PREC = 
FFLAGS = -c $(PREC) 
FFLAGS_EXE = $(PREC) 
LIBS = /Users/gaspardduchene/Work/Science/Disks/MC_simus/mcfost/lib
LIBS_PGPLOT = -lpgplot
LIBS_FITS = -lcfitsio
LIBS_SYST = -lappel_syst
# Numerical Recipes

MODS_NR = -I $(LIBS)/lib/nr
LIBS_NR =  -lnr
LIBS_NR_EQDIFF = -lnr_eq_diff
LIBS_SPLINE =  -lnr_splin
LIBS_SORT = -lnr_sort
LIBS_LU = -lnr_LU

# SPRNG
FLAGS_SPRNG = $(FLAGS) 
INCLUDE_SPRNG = -I $(LIBS)/include
LIBS_SPRNG =   -lsprng

# xlf sur la machine gagmac2 = xlfmac
ifeq ($(NAME),gagmac2.obs.ujf-grenoble.fr)
FC_SER = xlf90 -F/Users/fost/xlflib/xlf.cfg:xlf90_s
FC_OMP = xlf90_r -F/Users/fost/xlflib/xlf.cfg:xlf90_r_s -qsmp=omp -qnosave #-qreport=smpl\
ist
OPTS =  -qsuffix=f=f90 -O3 -qstrict #-lmass #-qsclk=micro
#-qhot fait des optimisations sur les boucles -> decalage des indices dans la sousroutine\
 opacite
DEBUG_OPTS = -g -qsigtrap -qcheck  -qnooptimize #-qflttrap=ov:und:zero:inv:imp:en
MACH = -qarch=auto -qtune=auto
FFLAGS = -c $(OPTS) $(MACH)
FFLAGS_EXE = $(OPTS) $(MACH) #-bmaxdata:0x80000000 -bmaxstack:0x20000000

LIBS = $(MCFOST_INSTALL)/lib/
LIBS_PGPLOT = -L $(LIBS) -lpgplot -lX11
LIBS_FITS =  -L $(LIBS) -lcfitsio 
COMP_FITS = -qextname=ftclos:ftcrhd:ftfiou:ftgerr:ftgiou:ftgmsg:ftinit:ftphpr:ftpkye:ftppre:ftgknj:ftgpve:ftopen:ftpprd:ftpprj:ftgpvj:ftmahd:ftpkys:ftpkyj
LIBS_SYST =  -L $(LIBS) -lappel_syst

# Numerical Recipes
MODS_NR = -I $(LIBS)/nr
LIBS_NR = -L $(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L $(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE = -L $(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG = -qsuffix=f=f90:cpp=f90
INCLUDE_SPRNG = -I $(MCFOST_INSTALL)/include/
LIBS_SPRNG =  -L $(LIBS) -lsprng
endif

endif

endif	

# Linux ou Mac x86, g95
ifeq ($(COMPILO),g95)
FC_SER = g95 -cpp 
FC_OMP = g95 -cpp
OPTS = -O3 -ffast-math #-mtune=athlon -mcpu=athlon -mfpmath=sse -msse 
DEBUG_OPTS =  -g -fbounds-check -ftrace=full -fimplicit-none #-g -C -ggdb -ftrace=full
# -Wall -Wextra -Wimplicit-none -Werror=113,115,137,146,147,159,163 -ftrace=full -fbounds-check -freal=nan -finteger=-2109876543 -fmodule-private 

#PARALLEL = -openmp
MACH = #-march=athlon
PREC = 
FFLAGS = -c $(PREC) 
FFLAGS_EXE = $(PREC) 

LIBS = $(MCFOST_INSTALL)/lib/g95/
LIBS_PGPLOT = -L$(LIBS) -lpgplot
LIBS_FITS = -L$(LIBS) -lcfitsio
LIBS_SYST = -L$(LIBS) -lappel_syst

# Numerical Recipes
MODS_NR = -I$(LIBS)/nr
LIBS_NR =  -L$(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L$(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE =  -L$(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L$(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
ifeq ($(ARCH),DARWIN) # g95 seulement en 32bits sur mac
FLAGS_SPRNG = $(FLAGS) 
endif
ifeq ($(ARCH),LINUX) # version 64 bits de g95 sur linux
FLAGS_SPRNG = $(FLAGS) -DPOINTER_SIZE=8
endif
INCLUDE_SPRNG = -I$(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L$(LIBS) -lsprng
endif

# Linux ou Mac x86, gfortran
ifeq ($(COMPILO),gfortran)
FC_SER = gfortran -cpp 
FC_OMP = gfortran -cpp -fopenmp
OPTS = -O3 -ffast-math #-mtune=athlon -mcpu=athlon -mfpmath=sse -msse 
DEBUG_OPTS =  -g -fbounds-check -ftrace=full -fimplicit-none #-g -C -ggdb -ftrace=full
# -Wall -Wextra -Wimplicit-none -Werror=113,115,137,146,147,159,163 -ftrace=full -fbounds-check -freal=nan -finteger=-2109876543 -fmodule-private 

#PARALLEL = -openmp
MACH = #-march=athlon
PREC = 
FFLAGS = -c $(PREC) 
FFLAGS_EXE = $(PREC) 

LIBS = $(MCFOST_INSTALL)/lib/gfortran/
LIBS_PGPLOT = -L$(LIBS) -lpgplot
LIBS_FITS = -L$(LIBS) -lcfitsio
LIBS_SYST = -L$(LIBS) -lappel_syst

# Numerical Recipes
MODS_NR = -I$(LIBS)/nr
LIBS_NR =  -L$(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L$(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE =  -L$(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L$(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG = $(FLAGS) 
INCLUDE_SPRNG = -I$(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L$(LIBS) -lsprng
endif


# force ifort en 32bits
ifeq ($(COMPILO),ifort32)
FC_SER = /opt/intel/Compiler/11.1/084/bin/ia32/ifort -fpp -traceback
FC_OMP = /opt/intel/Compiler/11.1/084/bin/ia32/ifort -fpp  -openmp -traceback
# fast = -O3 -ipo -no-prec-div -mdynamic-no-pi
OPTS = -fast  # -fast #-fp-model fast=2 -no-prec-div -no-prec-sqrt -mdynamic-no-pic -ftz -auto -ipo-jobs2
#OPTS = -ipo -mdynamic-no-pic -O3 -no-prec-div # -fast #-fp-model fast=2 -no-prec-div -no-prec-sqrt -mdynamic-no-pic -ftz -auto -ipo-jobs2
# fpe0 bug avec openmp
DEBUG_OPTS =  -check all  -C -g  -fpe0  -traceback -no-ftz # -no-ftz permet eviter bug dans conversion db -> sl 
#PARALLEL = -openmp
MACH = #-tpp7 #-xW -tune pn4
PREC = #-pc80 -prec_div -mp
# Option de conservation de la precision
# -mp : maintain precision fait prendre 1 facteur 2 dans le tps d'execution
# mais sinon donne des resultats faux
# -pc80 -prec_div ne fait pas perdre plus de tps
FFLAGS = -c  $(PREC) 
FFLAGS_EXE = $(PREC)  -static-intel

LIBS = $(MCFOST_INSTALL)/lib/ifort32/
LIBS_PGPLOT = -L $(LIBS) -lpgplot -L /usr/X11R6/lib/ -lX11 -static
LIBS_FITS = -L $(LIBS) -lcfitsio
LIBS_SYST = -L $(LIBS) -lappel_syst -Vaxlib

# Numerical recipes
MODS_NR = -I $(LIBS)/nr 
LIBS_NR = -L $(LIBS)/nr -lnr
LIBS_NR_EQDIFF = -L $(LIBS)/nr/eq_diff -lnr_eq_diff
LIBS_SPLINE = -L $(LIBS)/nr/spline -lnr_splin
LIBS_SORT = -L $(LIBS)/nr/sort -lnr_sort
LIBS_LU = -L $(LIBS)/nr/LU -lnr_LU

# SPRNG
FLAGS_SPRNG = $(FLAGS)  
INCLUDE_SPRNG = -I $(MCFOST_INSTALL)/include
LIBS_SPRNG =  -L $(LIBS) -lsprng
endif


FF = $(FC_OMP)  $(MACH) $(OPTS)
ifeq ($(debugp), yes)
FF = $(FC_OMP)  $(DEBUG_OPTS) 
endif
ifeq ($(debug), yes)
FF = $(FC_SER) $(DEBUG_OPTS) 
endif
ifeq ($(mono), yes)
FF = $(FC_SER)  $(MACH) $(OPTS)
endif
ifeq ($(prof), yes)
FF = $(FC_SER) -O0 -pg
endif
ifeq ($(dev), yes)
FF = $(FC_SER) -O1
endif
ifeq ($(devp), yes)
FF = $(FC_OMP) -O1
endif

SOURCES = modules.f90 sha.f90 mem.f90 utils.f90 input.f90 output.f90 \
	benchmarks.f90 	read_param.f90 io_prodimo.f90 init_mcfost.f90 \
	scattering.f90 coated_sphere.f90 dust_ray_tracing.f90 dust.f90 \
	molecules.f90 grid.f90 optical_depth.f90 disk_physics.f90 \
	stars.f90 density.f90 PAH.f90 diffusion.f90 \
	thermal_emission.f90 dust_transfer.f90\
	mol_transfer.f90 \
	mcfost.f90

OBJECTS = $(SOURCES:.f90=.o)

main : mcfost

#*************************************************************
# Code séquentiel et parallèle
#*************************************************************
mcfost : $(OBJECTS)  
	$(FF) $(FFLAGS_EXE) -o mcfost $(OBJECTS) $(COMP_FITS) $(LIBS_SYST) $(LIBS_NR_EQDIFF) $(LIBS_SPLINE) $(LIBS_SORT) $(LIBS_NR) $(LIBS_SPRNG) $(LIBS_FITS) $(INCLUDE_SPRNG)

# Fichiers necessaires aux suivants
modules.o : modules.f90
	$(FF) $(FFLAGS) $(FLAGS_SPRNG) modules.f90 $(INCLUDE_SPRNG) $(COMP_FITS)

# Fichier compile sans optimisation
mem.o : mem.f90  modules.o
	$(FF) $(FFLAGS) $(FLAGS_SPRNG) mem.f90 $(INCLUDE_SPRNG) $(COMP_FITS)

# defaut
%.o : %.f90  modules.o 
	$(FF) $(FFLAGS) $(FLAGS_SPRNG) $< -o $@ $(MODS_NR) $(INCLUDE_SPRNG) $(COMP_FITS)

io_prodimo.o : sha.f90
utils.o :      sha.f90

# Le fichier sha.f90 n'est cree que sur mon mac
ifeq ($(ARCH),DARWIN)
sha.f90 : ../.git/COMMIT_EDITMSG
	printf "module sha\n  implicit none\n" > sha.f90
	git log --pretty=format:"%H %d" | grep HEAD | awk '{print "  character(len=40), parameter :: sha_id = \""$$1"\""}' >> sha.f90
	printf "end module sha\n" >> sha.f90
endif

mcfost2prodimo : mcfost mcfost2prodimo.o
	$(FF) $(FFLAGS_EXE) -o mcfost2prodimo mcfost2prodimo.f90 $(COMP_FITS) $(LIBS_FITS)

#*************************************************************
# Release
#*************************************************************

release : main
	\cp -f mcfost $(MCFOST_INSTALL)/bin/
	rm -rf web
	mkdir web web/linux_64bits web/macos_intel_64bits
	cp mcfost web/macos_intel_64bits
	rsync -Pu *.f90 Makefile fosti:mcfost/src
	ssh fosti "cd mcfost/src ; make ; \cp -f mcfost mcfost_cigri"
	scp fosti:mcfost/src/mcfost web/linux_64bits/
	rm -rf ../utils/Version
	grep mcfost_version modules.f90 | awk '{print $$6}' | sed s/\"//g > ../utils/Version
	grep mcfost_release modules.f90 | awk '{print $$6}' | sed s/\"//g > web/version.txt
	git log --pretty=format:"%cd  %h %d %s" --date=short  --graph --branches > web/history.txt
	scp -r web/* horus:/user/publicdir/pintec/mcfost	
	rm -rf web

release_mac : main
	\cp -f mcfost $(MCFOST_INSTALL)/bin/
	rm -rf web
	mkdir web web/macos_intel_64bits
	cp mcfost web/macos_intel_64bits
	rm -rf ../utils/Version
	grep mcfost_version modules.f90 | awk '{print $$6}' | sed s/\"//g > ../utils/Version
	grep mcfost_release modules.f90 | awk '{print $$6}' | sed s/\"//g > web/version.txt
	git log --pretty=format:"%cd  %h  %d %s" --date=short --graph --branches > web/history.txt
	scp -r web/* horus:/user/publicdir/pintec/mcfost
	rm -rf web

release_utils : 
	tar czvf mcfost_utils.tgz -C../utils .
	scp -r mcfost_utils.tgz ../utils/Version horus:/user/publicdir/pintec/mcfost_utils/
	rm -f mcfost_utils.tgz 

#*************************************************************
# Nettoyage
#*************************************************************
clean : 
	rm -rf sha.f90 *.o *.exe *.mod outils/*.o outils/*.exe core*

cleandata : 
	rm -rf data* seed* density.fits.gz volume.fits.gz grid.fits.gz

cleanall : clean cleandata
	rm -rf mcfost