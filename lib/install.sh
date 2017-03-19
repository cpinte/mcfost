#!/bin/sh

if [ "$SYSTEM" = "ifort" ] ; then
    echo "Building MCFOST's libraries with ifort"
    export CC=icc

elif [ "$SYSTEM" = "gfortran" ] ; then
    echo "Building MCFOST's libraries with gfortran"
    export CC=gcc

else
    echo "Unknown system to build mcfost:"
    echo $SYSTEM
    echo ""
    echo "Please choose ifort or gfortran"
    echo "Exiting"
    exit
fi

rm -rf lib include
mkdir lib
mkdir include

# SPRNG
# wget http://sprng.cs.fsu.edu/Version2.0/sprng2.0b.tar.gz
tar xzvf sprng2.0b.tar.gz
\cp -f $SYSTEM/make.CHOICES sprng2.0
\cp -f $SYSTEM/make.INTEL sprng2.0/SRC
if [ $(uname | tr '[a-z]' '[A-Z]' 2>&1 | grep -c DARWIN) -eq 1 ]; then
  \cp macos/insertmenu.mac sprng2.0/SRC/insertmenu
fi
cd sprng2.0
make -B
\cp lib/libsprng.a ../lib
\cp include/*.h ../include
cd ..
rm -rf sprng2.0

# cfitsio
# g77 ou f77 needed by configure to set up the fortran wrapper in Makefile
# wget ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3030.tar.gz
tar xzvf cfitsio3030.tar.gz
cd cfitsio
if [ "$SYSTEM" = "gfortran" ] ; then
    export CFLAGS="-m64"
    export FC="gfortran"
fi
./configure
make
\cp libcfitsio.a ../lib
cd ..
rm -rf cfitsio

# voro++
# Original voro++ can be obtained from
# svn checkout --username anonsvn https://code.lbl.gov/svn/voro/trunk voro
# Downloading last version tested with mcfost :
git clone https://cpinte@bitbucket.org/cpinte/voro.git
if [ "$SYSTEM" = "ifort" ] ; then
    \cp -f  ifort/config.mk voro
fi
cd voro
make
\cp src/libvoro++.a ../lib
mkdir -p ../include/voro++
\cp src/*.hh ../include/voro++/
cd ..
rm -rf voro

# Numerical recipes
mkdir lib/nr lib/nr/eq_diff lib/nr/spline lib/nr/sort
cd nr
if [ "$SYSTEM" = "ifort" ] ; then
    ./compile_ifort.sh
elif [ "$SYSTEM" = "gfortran" ] ; then
    ./compile_gfortran.sh
fi
\cp libnr.a *.mod ../lib/nr
\cp eq_diff/libnr_eq_diff.a eq_diff/*.mod ../lib/nr/eq_diff
\cp spline/libnr_splin.a ../lib/nr/spline
\cp sort/libnr_sort.a ../lib/nr/sort
./clean.sh
cd ..

# Put in directory
mkdir -p $MCFOST_INSTALL/include
\cp -r include/* $MCFOST_INSTALL/include
mkdir -p $MCFOST_INSTALL/lib/$SYSTEM
\cp -r lib/* $MCFOST_INSTALL/lib/$SYSTEM
