#!/bin/bash
for comm in svn make tar
do
    command -v $comm
    if [ $? != 0 ] ; then echo "error: $comm command not found"; exit 1; fi
done
set -eu

#preliminary checks
if [ ! $# = 0 ]; then SYSTEM=$1 ; fi

set +u # personalized error messages
if [ -z $SYSTEM ]; then echo "error: SYSTEM need to be set (choose: ifort, gfortran or xeon-phi)." ; exit 1 ; fi
if [ -z $MCFOST_INSTALL ]; then echo "error: MCFOST_INSTALL needs to point to a directory." ; exit 1 ; fi
set -u

if [ $SYSTEM = "ifort" ] ; then
    echo "Building MCFOST's libraries with ifort" ; export CC=icc
elif [ $SYSTEM = "gfortran" ] ; then
    echo "Building MCFOST's libraries with gfortran" ; export CC=gcc
elif [ $SYSTEM = "xeon-phi" ] ; then
    echo "Building MCFOST's libraries with ifort for Xeon-Phi" ; export CC=icc
else
    echo "Unknown system to build mcfost: "$SYSTEM"\nPlease choose ifort or gfortran\ninstall.sh <system>\nExiting" ; exit 1
fi

rm -rf lib include sprng2.0 cfitsio voro # Clean previous files if any


echo "Installing MCFOST libraries in "$MCFOST_INSTALL/lib/$SYSTEM
mkdir lib include
pushd .

#-------------------------------------------
# SPRNG
#-------------------------------------------
# wget http://sprng.cs.fsu.edu/Version2.0/sprng2.0b.tar.gz
tar xzvf sprng2.0b.tar.gz
\cp -f $SYSTEM/make.CHOICES sprng2.0
\cp -f $SYSTEM/make.INTEL sprng2.0/SRC
if [ $(uname | tr '[a-z]' '[A-Z]' 2>&1 | grep -c DARWIN) -eq 1 ]; then \cp -f macos/insertmenu.mac sprng2.0/SRC/insertmenu ; fi

cd sprng2.0
make -B
\cp lib/libsprng.a ../lib
\cp include/*.h ../include
cd ~1

#-------------------------------------------
# cfitsio
#-------------------------------------------
# g77 ou f77 needed by configure to set up the fortran wrapper in Makefile
# wget ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3420.tar.gz
tar xzvf cfitsio3420.tar.gz
if [ "$SYSTEM" = "ifort" ] ; then
    export CC="icc" ; export FC="ifort"
elif [ "$SYSTEM" = "gfortran" ] ; then
    export CFLAGS="-m64" ; export FC="gfortran"
elif [ "$SYSTEM" = "xeon-phi" ] ; then
    export CFLAGS=-mmic
fi

cd cfitsio
./configure --enable-ssse3
# Tweaking Makefile
if [ "$SYSTEM" = "xeon-phi" ] ; then \cp -f ../Makefile_cfitsio.xeon_phi Makefile ; fi
cat Makefile | sed s/-DCFITSIO_HAVE_CURL=1// | sed s/-lcurl// >> Makefile.tmp && \mv -f Makefile.tmp Makefile # We do not want to have to link with libcurl
make
\cp libcfitsio.a ../lib
cd ~1

#-------------------------------------------
# voro++
#-------------------------------------------
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

# Put in final directory
mkdir -p $MCFOST_INSTALL/include
\cp -r include/*.h include/voro++ $MCFOST_INSTALL/include/
mkdir -p $MCFOST_INSTALL/lib/$SYSTEM
\cp -r lib/*.a $MCFOST_INSTALL/lib/$SYSTEM/

# Final cleaning
rm -rf lib include sprng2.0 cfitsio voro
popd
