export CC=gcc

rm -rf lib include
mkdir lib
mkdir include

# SPRNG
# gcc-2.95 needed
# wget http://sprng.cs.fsu.edu/Version2.0/sprng2.0b.tar.gz
tar xzvf sprng2.0b.tar.gz
\cp -f gfortran/make.CHOICES sprng2.0
\cp -f gfortran/make.INTEL sprng2.0/SRC
if [ $(uname | tr '[a-z]' '[A-Z]' 2>&1 | grep -c DARWIN) -eq 1 ]; then
  \cp insertmenu.mac sprng2.0/SRC/insertmenu
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
export CFLAGS="-m64"
export FC="gfortran"
./configure
make
\cp libcfitsio.a ../lib
cd ..
rm -rf cfitsio

# voro++
tar xzvf voro++_594.tgz
cd trunk
make
\cp src/libvoro++.a ../lib
mkdir -p ../include/voro++
\cp src/*.hh ../include/voro++/
cd ..
rm -rf trunk

# Numerical recipes
mkdir lib/nr lib/nr/eq_diff lib/nr/spline lib/nr/sort
cd nr
./compile_gfortran.sh
\cp libnr.a *.mod ../lib/nr
\cp eq_diff/libnr_eq_diff.a eq_diff/*.mod ../lib/nr/eq_diff
\cp spline/libnr_splin.a ../lib/nr/spline
\cp sort/libnr_sort.a ../lib/nr/sort
./clean.sh
cd ..

# Put in directory
mkdir -p $MCFOST_INSTALL/include
\cp -r include/* $MCFOST_INSTALL/include
mkdir -p $MCFOST_INSTALL/lib/gfortran
\cp -r lib/* $MCFOST_INSTALL/lib/gfortran
