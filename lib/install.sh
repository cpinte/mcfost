#--------------------------------------------------------
# This script downloads and installs the libraries
# required by mcfost:
#  - sprng2
#  - cfitsio
#  - voro++
#  - xgboost and dependencies
#  - hdf5
#
# It is likely that you some of the libaries already
# available on your system. You are free to use them
# if you wish a "cleaner" installation. The versions
# provided here have been tested with mcfost.
#
# The libraries will be installed in $MCFOST_INSTALL
#
# Please refer to the licenses for each library.
#--------------------------------------------------------
#!/bin/bash
set -eu

#-- Preliminary checks
for comm in svn make tar
do
    command -v $comm
    if [ $? != 0 ] ; then echo "error: $comm command not found"; exit 1; fi
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
rm -rf lib include sprng2.0 cfitsio voro xgboost
mkdir lib include
mkdir include/$SYSTEM
mkdir include/hdf5
pushd .

#-- Downloading libraries
wget -N http://sprng.org/Version2.0/sprng2.0b.tar.gz
wget -N http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.47.tar.gz
if [ skip_hdf5 != "yes" ] ; then
    wget -N https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.bz2
fi
svn checkout --username anonsvn --password anonsvn https://code.lbl.gov/svn/voro/trunk voro
git clone --recursive https://github.com/dmlc/xgboost


#-------------------------------------------
# SPRNG
#-------------------------------------------
echo "Compiling sprng ..."
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
tar xzvf cfitsio-3.47.tar.gz
mv cfitsio-3.47 cfitsio
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
if [ "$SYSTEM" = "ifort" ] ; then
    \cp -f  ifort/config.mk voro
elif [ "$SYSTEM" = "xeon-phi" ] ; then
    \cp -f  ifort/config.mk voro # To be tested
elif [ "$SYSTEM" = "gfortran" ] ; then
    \cp -f  gfortran/config.mk voro
fi

cd voro
svn up -r604
make
\cp src/libvoro++.a ../lib
mkdir -p ../include/voro++ ; \cp src/*.hh ../include/voro++/
cd ~1
echo "Done"


#-------------------------------------------
# xgboost
#-------------------------------------------
echo "Compiling xgboost ..."
cd xgboost
git checkout v0.90
#-- we remove the test for the moment even if this works for gfortran
#if [ "$SYSTEM" = "ifort" ] ; then
\cp ../ifort/xgboost/base.h include/xgboost/base.h
#fi
\cp ../ifort/xgboost/rabit/Makefile rabit/ # g++ was hard-coded in the Mekefile
make -j
\cp dmlc-core/libdmlc.a rabit/lib/librabit.a lib/libxgboost.a ../lib
\cp -r dmlc-core/include/dmlc rabit/include/rabit include/xgboost ../include
#We will need the .h file when linking mcfost, the path is hard-coded in the the xgboost files
mkdir -p  $MCFOST_INSTALL/src/common
\cp -r src/common/*.h $MCFOST_INSTALL/src/common
cd ~1
echo "Done"


#---------------------------------------------
# hdf5 : you can skip hdf5 (slow to compile)
# and use system library if prefered
# In that case, you need to define HDF5ROOT
# for the mcfost Makefile
#---------------------------------------------
if [ skip_hdf5 != "yes" ] ; then
    wget -N https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.bz2
    echo "Compiling hdf5 ..."
    tar xjvf hdf5-1.10.5.tar.bz2 ; mv hdf5-1.10.5 hdf5
    mkdir -p $HOME/hdf5_install_tmp
    cd hdf5
    ./configure --prefix=$HOME/hdf5_install_tmp --enable-fortran --disable-shared
    make -j -l6 install
    cd ~1
    \cp $HOME/hdf5_install_tmp/lib/libhdf5.a $HOME/hdf5_install_tmp/lib/libhdf5_fortran.a lib/
    \cp $HOME/hdf5_install_tmp/include/*.h include/hdf5/
    \cp $HOME/hdf5_install_tmp/include/*.mod include/$SYSTEM
    echo "Done"
fi

#-- Put in final directory
echo "Installing MCFOST libraries in "$MCFOST_INSTALL/lib/$SYSTEM
mkdir -p $MCFOST_INSTALL/include
\cp -r include/* $MCFOST_INSTALL/include/
mkdir -p $MCFOST_INSTALL/lib/$SYSTEM
\cp -r lib/*.a $MCFOST_INSTALL/lib/$SYSTEM/

#-- Final cleaning
./clean.sh
popd
