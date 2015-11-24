rm -rf lib include
mkdir lib
mkdir include

# SPRNG
# gcc-2.95 needed
#wget http://sprng.cs.fsu.edu/Version2.0/sprng2.0b.tar.gz
tar xzvf sprng2.0b.tar.gz
cp -f  linux/ifort64/make.CHOICES sprng2.0
cp  linux/ifort64/make.IFORT64 sprng2.0/SRC
cd sprng2.0
make
mv lib/libsprng.a ../lib
mv include/*.h ../include
cd ..
rm -rf sprng2.0

# cfitsio
# g77 ou f77 needed by configure to set up the fotran wrapper in Makefile
# wget ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3030.tar.gz
tar xzvf cfitsio3030.tar.gz
cd cfitsio
./configure
make
cp libcfitsio.a ../lib
cd ..
rm -rf cfitsio

# voro++
tar xzvf voro++-0.4.6.tar.gz
\cp -f  macos/ifort64/config.mk voro++-0.4.6
cd voro++-0.4.6
make
\cp src/libvoro++.a ../lib
\cp src/voro++.hh ../include
cd ..
rm -rf voro++-0.4.6

# Numerical recipes
mkdir lib/nr lib/nr/eq_diff lib/nr/spline lib/nr/sort
cd nr
./compile_linux.com
cp libnr.a *.mod ../lib/nr
cp eq_diff/libnr_eq_diff.a eq_diff/*.mod ../lib/nr/eq_diff
cp spline/libnr_splin.a ../lib/nr/spline
cp sort/libnr_sort.a ../lib/nr/sort
./clean.com
cd ..

mkdir -p $MCFOST_INSTALL/include
\cp -r include/* $MCFOST_INSTALL/include
mkdir -p $MCFOST_INSTALL/lib/ifort
\cp -r lib/* $MCFOST_INSTALL/lib/ifort
