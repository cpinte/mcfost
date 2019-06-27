#!/bin/bash
set -eu

#-- Preliminary checks
for comm in svn make tar git-lfs
do
    command -v $comm
    if [ $? != 0 ] ; then echo "error: $comm command not found"; exit 1; fi
done


du *tar.gz | while read line
do
   typeset -i size=$(echo $line | awk '{print $1}')
   if (($size<1000)); then
       git lfs pull
   fi
done

if [ ! $# = 0 ]; then SYSTEM=$1 ; fi

set +u # personalized error messages
if [ -z $SYSTEM ]; then echo "error: SYSTEM need to be set (choose: ifort, gfortran or xeon-phi)." ; exit 1 ; fi
if [ -z $MCFOST_INSTALL ]; then echo "error: MCFOST_INSTALL needs to point to a directory." ; exit 1 ; fi
set -u

if [ $SYSTEM = "ifort" ] ; then
    echo "Building MCFOST's libraries with ifort"
    export CC=icc
    export FC=ifort
    export CXX=icpc
elif [ $SYSTEM = "gfortran" ] ; then
    echo "Building MCFOST's libraries with gfortran"
    export CC=gcc
    export FC=gfortran
    export CXX=g++
    export CFLAGS="-m64"
elif [ $SYSTEM = "xeon-phi" ] ; then
    echo "Building MCFOST's libraries with ifort for Xeon-Phi"
    export CC=icc
    export FC=ifort
    export CXX=icpc
    export CFLAGS=-mmic
else
    echo "Unknown system to build mcfost: "$SYSTEM"\nPlease choose ifort or gfortran\ninstall.sh <system>\nExiting" ; exit 1
fi


#-- Clean previous files if any
rm -rf lib include sprng2.0 cfitsio voro
mkdir lib include
pushd .

#-------------------------------------------
# SPRNG
#-------------------------------------------
echo "Compiling sprng ..."
#Original version from http://sprng.cs.fsu.edu/Version2.0/sprng2.0b.tar.gz
tar xzvf sprng2.0b.tar.gz
\cp -f $SYSTEM/make.CHOICES sprng2.0
\cp -f $SYSTEM/make.INTEL sprng2.0/SRC
if [ $(uname | tr '[a-z]' '[A-Z]' 2>&1 | grep -c DARWIN) -eq 1 ]; then \cp -f macos/insertmenu.mac sprng2.0/SRC/insertmenu ; fi

cd sprng2.0
make -B
\cp lib/libsprng.a ../lib
\cp include/*.h ../include
cd ~1
echo "Done"

#-------------------------------------------
# cfitsio
#-------------------------------------------
echo "Compiling cfitsio ..."
# g77 ou f77 needed by configure to set up the fortran wrapper in Makefile
# Original version from ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3420.tar.gz
tar xzvf cfitsio3420.tar.gz
cd cfitsio
./configure --enable-ssse3

#--- Tweaking Makefile
if [ "$SYSTEM" = "xeon-phi" ] ; then \cp -f ../Makefile_cfitsio.xeon_phi Makefile ; fi
# We do not want to have to link with libcurl
cat Makefile | sed s/-DCFITSIO_HAVE_CURL=1// | sed s/-lcurl// >> Makefile.tmp && \mv -f Makefile.tmp Makefile

make
\cp libcfitsio.a ../lib
cd ~1
echo "Done"

#-------------------------------------------
# voro++
#-------------------------------------------
echo "Compiling voro++ ..."
# Downloading last version tested with mcfost : git clone git@bitbucket.org:cpinte/voro.git
# Original voro++ can be obtained from svn checkout --username anonsvn https://code.lbl.gov/svn/voro/trunk voro
svn checkout --username anonsvn --password anonsvn https://code.lbl.gov/svn/voro/trunk voro
if [ "$SYSTEM" = "ifort" ] ; then
    \cp -f  ifort/config.mk voro
elif [ "$SYSTEM" = "xeon-phi" ] ; then
    \cp -f  ifort/config.mk voro # To be tested
elif [ "$SYSTEM" = "gfortran" ] ; then
    \cp -f  gfortran/config.mk voro
fi

cd voro
make
\cp src/libvoro++.a ../lib
mkdir -p ../include/voro++ ; \cp src/*.hh ../include/voro++/
cd ~1
echo "Done"

#---------------------------------------------
# hdf5 : you can skip hdf5 (slow to compile)
# and use system library if prefered
# In that case, you need to define HDF5ROOT
# for the mcfost Makefile
#---------------------------------------------
if [ skip_hdf5 != "yes" ] ; then
    echo "Compiling hdf5 ..."
    tar xjvf hdf5-1.10.5.tar.bz2 ; mv hdf5-1.10.5 hdf5
    cd hdf5
    ./configure --prefix=$MCFOST_INSTALL/hdf5/$SYSTEM --enable-fortran
    make install
    cd ~1
    echo "Done"
fi

#-- Put in final directory
echo "Installing MCFOST libraries in "$MCFOST_INSTALL/lib/$SYSTEM
mkdir -p $MCFOST_INSTALL/include
\cp -r include/*.h include/voro++ $MCFOST_INSTALL/include/
mkdir -p $MCFOST_INSTALL/lib/$SYSTEM
\cp -r lib/*.a $MCFOST_INSTALL/lib/$SYSTEM/

#-- Final cleaning
rm -rf lib include sprng2.0 cfitsio voro hdf5
popd
