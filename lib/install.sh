#!/bin/sh

#preliminary checks

if [ -z $MCFOST_INSTALL ]; then
    echo "error: MCFOST_INSTALL needs to be set."
    exit 1
fi
if [ ! -d $MCFOST_INSTALL ]; then
    echo "error: MCFOST_INSTALL needs to point to an existing directory."
    exit 1
fi

if [ ! $# = 0 ]; then
    SYSTEM=$1
fi

if [ -z $SYSTEM ]; then
    echo "error: SYSTEM need to be set (choose: ifort, gfortran or xeon-phi)"
    exit 1
fi

if [ $SYSTEM = "ifort" ] ; then
    echo "Building MCFOST's libraries with ifort"
    export CC=icc

elif [ $SYSTEM = "gfortran" ] ; then
    echo "Building MCFOST's libraries with gfortran"
    export CC=gcc

elif [ $SYSTEM = "xeon-phi" ] ; then
    echo "Building MCFOST's libraries with ifort for Xeon-Phi"
    export CC=icc

else
    echo "Unknown system to build mcfost."
    echo $SYSTEM
    echo ""
    echo "Please choose ifort or gfortran"
    echo "install.sh <system>"
    echo "Exiting"
    exit
fi

#clear previous files
rm -rf lib include
mkdir lib
mkdir include

ANCHOR=$(pwd)



# SPRNG
# wget http://sprng.cs.fsu.edu/Version2.0/sprng2.0b.tar.gz
tar xzvf sprng2.0b.tar.gz \
&& \cp -f $SYSTEM/make.CHOICES sprng2.0 \
&& \cp -f $SYSTEM/make.INTEL sprng2.0/SRC \
&& if [ $(uname | tr '[a-z]' '[A-Z]' 2>&1 | grep -c DARWIN) -eq 1 ]; then
  \cp -f macos/insertmenu.mac sprng2.0/SRC/insertmenu
fi
cd sprng2.0 \
&& make -B \
&& \cp lib/libsprng.a ../lib \
&& \cp include/*.h ../include

cd $ANCHOR
rm -rf sprng2.0



# cfitsio
# g77 ou f77 needed by configure to set up the fortran wrapper in Makefile
# wget ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3420.tar.gz
tar xzvf cfitsio3420.tar.gz \
&& cd cfitsio \
&& if [ "$SYSTEM" = "ifort" ] ; then
    export CC="icc"
    export FC="ifort"
elif [ "$SYSTEM" = "gfortran" ] ; then
    export CFLAGS="-m64"
    export FC="gfortran"
elif [ "$SYSTEM" = "xeon-phi" ] ; then
    export CFLAGS=-mmic
fi \
&& ./configure --enable-ssse3 \
&& if [ "$SYSTEM" = "xeon-phi" ] ; then
    \cp -f ../Makefile_cfitsio.xeon_phi Makefile
fi
# We do not want to have to link with libcurl
cat Makefile | sed s/-DCFITSIO_HAVE_CURL=1// | sed s/-lcurl// >> Makefile.tmp && \mv -f Makefile.tmp Makefile \
&& make \
&& \cp libcfitsio.a ../lib

cd $ANCHOR
rm -rf cfitsio



# voro++

# Downloading last version tested with mcfost :
# git clone git@bitbucket.org:cpinte/voro.git
# Original voro++ can be obtained from
echo
echo "password is 'anonsvn'"
echo
svn checkout --username anonsvn https://code.lbl.gov/svn/voro/trunk voro \
&& if [ "$SYSTEM" = "ifort" ] ; then
    \cp -f  ifort/config.mk voro
elif [ "$SYSTEM" = "xeon-phi" ] ; then
    \cp -f  ifort/config.mk voro # To be tested
elif [ "$SYSTEM" = "gfortran" ] ; then
    \cp -f  gfortran/config.mk voro
fi \
&& cd voro \
&& make \
&& \cp src/libvoro++.a ../lib \
&& mkdir -p ../include/voro++ \
&& \cp src/*.hh ../include/voro++/ \

cd $ANCHOR
rm -rf voro



# Put in final directory
mkdir -p $MCFOST_INSTALL/include
\cp -r include/*.h include/voro++ $MCFOST_INSTALL/include/
mkdir -p $MCFOST_INSTALL/lib/$SYSTEM
\cp -r lib/*.a $MCFOST_INSTALL/lib/$SYSTEM/
