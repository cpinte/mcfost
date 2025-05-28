#!/usr/bin/env bash

#--------------------------------------------------------
#
# This script downloads and installs the libraries
# required by MCFOST:
#
#  - SPRNG
#  - CFITSIO
#  - Voro++
#  - XGBoost and dependencies (optional)
#  - astrochem and dependencies (optional)
#  - HDF5
#
# It is likely that you some of the libaries already
# available on your system. You are free to use them
# if you wish a "cleaner" installation. The versions
# provided here have been tested with MCFOST.
#
# The libraries will be installed in $MCFOST_INSTALL
#
# Please refer to the licenses for each library.
#
# Variables
# ---------
#
# This build script reads the following variables from
# the environment in which it is run.
#
# SYSTEM         : the compiler system. Options are
#                  "ifort" or "gfortran".
#
# MCFOST_INSTALL : the git repository directory for
#                  MCFOST. The libraries will be
#                  installed under this directory.
#
# MCFOST_XGBOOST   : a boolean flag ("yes" or "no") to
#                  compile the XGBoost libraries.
#                  If "yes", you also need to compile MCFOST
#                  with MCFOST_XGBOOST=yes to use XGBoost.
#                  Optional, default value is "no".
#
# MCFOST_ASTROCHEM : a boolean flag ("yes" or "no") to
#                  compile the astrochem libraries.
#                  If "yes", you also need to compile MCFOST
#                  with MCFOST_ASTROCHEM=yes to use astrochem.
#                  Optional, default value is "no".
#
# SKIP_HDF5      : a boolean flag ("yes" or "no") to
#                  skip compiling the HDF5 libraries.
#                  If "yes", you must provide HDF5_DIR
#                  when compiling MCFOST. Optional,
#                  default value is "no".
#
#--------------------------------------------------------

set -eu

#-- Preliminary checks
for comm in make tar wget; do
    if ! command -v $comm >/dev/null 2>&1; then
        echo "error: $comm command not found"
        exit 1
    fi
done

if [ ! $# = 0 ]; then SYSTEM=$1; fi

set +u # personalized error messages
if [ -z "$SYSTEM" ]; then
    echo "error: SYSTEM need to be set (choose: ifort, ifx, or gfortran)."
    exit 1
fi
if [ -z "$MCFOST_INSTALL" ]; then
    echo "error: MCFOST_INSTALL needs to point to a directory."
    exit 1
fi
set -u

if [ "$SYSTEM" = "ifort" ]; then
    echo "Building MCFOST's libraries with ifort"
    export CC=icc
    export FC=ifort
    export CXX=icpc
    export CFLAGS=""
elif [ "$SYSTEM" = "ifx" ]; then
    echo "Building MCFOST's libraries with ifx"
    export CC=icx
    export FC=ifx
    export CXX=icpx
    export CFLAGS="-std=gnu90"
elif [ "$SYSTEM" = "gfortran" ]; then
    echo "Building MCFOST's libraries with gfortran"
    export CC=gcc
    export FC=gfortran
    export CXX=g++
    export CFLAGS="-m64"
else
    echo "Unknown system to build MCFOST: $SYSTEM"
    echo "Please choose ifort, ifx, or gfortran"
    echo "install.sh <system>"
    echo "Exiting"
    exit 1
fi


# Detecting operating system
if [ "$(uname | tr '[:lower:]' '[:upper:]' 2>&1 | grep -c DARWIN)" -eq 1 ]; then
    os="macos"
else
    os="linux"
fi


#-- Check if SKIP_HDF5 is set, if not, set to 'no'
if [ -z ${SKIP_HDF5+x} ]; then SKIP_HDF5=no; fi

#-- Check if MCFOST_ASTROCHEM is set, if not, set to 'no'
if [ -z ${MCFOST_ASTROCHEM+x} ]; then MCFOST_ASTROCHEM=no; fi

#-- Check if MCFOST_XGBOOST is set, if not, set to 'no'
if [ -z ${MCFOST_XGBOOST+x} ]; then MCFOST_XGBOOST=no; fi

#-- Clean previous files if any
rm -rf lib include sprng2.0 cfitsio voro xgboost
mkdir lib include
mkdir "include/$SYSTEM"
mkdir include/hdf5
pushd .

#-- Downloading libraries
wget -N http://sprng.org/Version2.0/sprng2.0b.tar.gz
wget -N http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.3.0.tar.gz
git clone https://github.com/cpinte/voro
if [ "$SKIP_HDF5" != "yes" ]; then
    wget -N https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.bz2
fi
if [ "$MCFOST_XGBOOST" = "yes" ]; then
    git clone --recursive https://github.com/dmlc/xgboost
    cd xgboost
    git checkout v0.90
    # We need to install rabit manually as the submodule was removed
    \rm -rf rabit
    git clone https://github.com/dmlc/rabit
    cd ..
fi
if [ "$MCFOST_ASTROCHEM" = "yes" ]; then
    git clone https://github.com/cpinte/astrochem
fi

#-------------------------------------------
# SPRNG
#-------------------------------------------
echo "Compiling SPRNG ..."
tar xzvf sprng2.0b.tar.gz
\cp -f "$SYSTEM/make.CHOICES" sprng2.0
\cp -f "$SYSTEM/make.INTEL" sprng2.0/SRC
if [ "$SYSTEM" = "gfortran" ] ; then
    \cp -f "$SYSTEM/primes_64.c" sprng2.0/SRC # recent gfortran has an issue with having duplicate symbols between primes_32 and primes_64
fi
if [ "$SYSTEM" = "ifx" ]; then
    patch -p1 -d sprng2.0 < ifx/sprng2.0b.patch
fi
if [ "$os" = "macos" ]; then
    \cp -f macos/insertmenu.mac sprng2.0/SRC/insertmenu
fi

cd sprng2.0
make -B
\cp lib/libsprng.a ../lib
\cp include/*.h ../include
cd ~1
echo "Done"

#-------------------------------------------
# CFITSIO
#-------------------------------------------
echo "Compiling CFITSIO ..."
tar xzvf cfitsio-4.3.0.tar.gz
mv cfitsio-4.3.0 cfitsio
cd cfitsio
./configure --enable-ssse3 --disable-curl

make
\cp libcfitsio.a ../lib
cd ~1
echo "Done"

#-------------------------------------------
# Voro++
#-------------------------------------------
echo "Compiling Voro++ library ..."
if [ "$SYSTEM" = "ifort" ]; then
    \cp -f ifort/config.mk voro
elif [ "$SYSTEM" = "ifx" ]; then
    \cp -f ifx/config.mk voro
elif [ "$SYSTEM" = "gfortran" ]; then
    \cp -f gfortran/config.mk voro
fi

cd voro
# Allowing for up to 1e8 particles (1e7 by default)
\cp -f ../voro++/config.hh src/
cd src
pwd
make libvoro++.a
cd -
\cp src/libvoro++.a ../lib
mkdir -p ../include/voro++
\cp src/*.hh ../include/voro++/
cd ~1
echo "Done"

#-------------------------------------------
# XGBoost
#-------------------------------------------
if [ "$MCFOST_XGBOOST" = "yes" ]; then
    echo "Compiling XGBoost ..."
    cd xgboost
    \cp ../ifort/xgboost/rabit/Makefile rabit/ # g++ is hard-coded in the Makefile

    if [ "$os" = "linux" ]; then
        # Forcing g++ for now as there are some issues with ifort 2020 on linux
        CC_old=$CC
        FC_old=$FC
        CXX_old=$CXX
        CFLAGS_old=$CFLAGS

        export CC=gcc
        export FC=gfortran
        export CXX=g++
        export CFLAGS="-m64"
    fi

    #-- we remove the ifort test for the moment even if the default works for gfortran
    #if [ "$SYSTEM" = "ifort" ] ; then
    \cp ../ifort/xgboost/base.h include/xgboost/base.h
    #fi

    make -j
    \cp dmlc-core/libdmlc.a rabit/lib/librabit.a lib/libxgboost.a ../lib
    \cp -r dmlc-core/include/dmlc rabit/include/rabit include/xgboost ../include
    # We will need the .h file when linking MCFOST, the path is hard-coded in the the xgboost files
    mkdir -p "$MCFOST_INSTALL/src/common"
    \cp -r src/common/*.h "$MCFOST_INSTALL/src/common"
    cd ~1
    echo "Done"

    if [ "$os" = "linux" ]; then
        # Restoring flags
        export CC=$CC_old
        export FC=$FC_old
        export CXX=$CXX_old
        export CFLAGS=$CFLAGS_old
    fi
else
    echo "Skipping XGBoost ..."
    echo "Make sure that MCFOST_XGBOOST is not set yes when compiling MCFOST"
fi

#---------------------------------------------
# Astrochem (does not compile with sundials 6)
#---------------------------------------------
if [ "$MCFOST_ASTROCHEM" = "yes" ]; then
    conda create --name astrochem -c conda-forge sundials=5.7.0 python cython numpy matplotlib h5py autoconf automake libtool
    conda activate astrochem
    cd astrochem
    autoupdate
    ./bootstrap
    ./configure CPPFLAGS="-I$CONDA_PREFIX/include" LDFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib -L/$CONDA_PREFIX/lib"  --prefix=$CONDA_PREFIX
    make
    \cp ./src/.libs/libastrochem.a ../lib
    \cp ./src/libastrochem.h ../include
    cd ~
    conda deactivate
else
    echo "Skipping Astrochem"
fi


#---------------------------------------------
# HDF5
#---------------------------------------------
if [ "$SKIP_HDF5" != "yes" ]; then
    wget -N https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.bz2
    echo "Compiling HDF5 ..."
    tar xjvf hdf5-1.10.5.tar.bz2
    cd hdf5-1.10.5
    mkdir -p "$HOME/hdf5_install_tmp"
    ./configure --prefix="$HOME/hdf5_install_tmp" --enable-fortran --disable-shared
    export CFLAGS="-m64 -fpermisive"
    make install
    cd ~1
    \cp "$HOME/hdf5_install_tmp/lib/libhdf5.a" "$HOME/hdf5_install_tmp/lib/libhdf5_fortran.a" lib/
    \cp "$HOME"/hdf5_install_tmp/include/*.h include/hdf5/
    \cp "$HOME"/hdf5_install_tmp/include/*.mod "include/$SYSTEM"
    echo "Done"
else
    echo "Skipping HDF5 ..."
    echo "Make sure to set HDF5_DIR when compiling MCFOST"
fi

#-- Put in final directory
echo "Installing MCFOST libraries in $MCFOST_INSTALL/lib/$SYSTEM"
mkdir -p "$MCFOST_INSTALL/include"
\cp -r include/* "$MCFOST_INSTALL/include/"
mkdir -p "$MCFOST_INSTALL/lib/$SYSTEM"
\cp -r lib/*.a "$MCFOST_INSTALL/lib/$SYSTEM"

#-- Final cleaning
./clean.sh
popd
