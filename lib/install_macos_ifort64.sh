export CC=cc # cfitsio does not compile with icc
export CFLAGS=-mmacosx-version-min=10.6

mkdir lib
mkdir include

# SPRNG
# gcc-2.95 needed
#wget http://sprng.cs.fsu.edu/Version2.0/sprng2.0b.tar.gz
tar xzvf sprng2.0b.tar.gz
\cp -f  macos/ifort64/make.CHOICES sprng2.0
\cp -f macos/ifort64/make.IFORT sprng2.0/SRC
cd sprng2.0
make -B
mv lib/libsprng.a ../lib
mv include/*.h ../include
cd ..
rm -rf sprng2.0

# cfitsio
# g77 ou f77 needed by configure to set up the fotran wrapper in Makefile
#wget ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3030.tar.gz
tar xzvf cfitsio3030.tar.gz
cd cfitsio
./configure
make
\cp libcfitsio.a ../lib
cd ..
rm -rf cfitsio

# Numerical recipes
mkdir lib/nr lib/nr/eq_diff lib/nr/spline lib/nr/sort
cd nr
./compile_linux.com
\cp libnr.a *.mod ../lib/nr
\cp eq_diff/libnr_eq_diff.a eq_diff/*.mod ../lib/nr/eq_diff
\cp spline/libnr_splin.a ../lib/nr/spline
\cp sort/libnr_sort.a ../lib/nr/sort
./clean.com
cd ..

cp -r include $MCFOST_INSTALL
mkdir $MCFOST_INSTALL/lib
cp -r lib $MCFOST_INSTALL/lib/ifort
