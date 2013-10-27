mkdir lib
mkdir include

# SPRNG
# gcc-2.95 needed
#wget http://sprng.cs.fsu.edu/Version2.0/sprng2.0b.tar.gz
tar xzvf sprng2.0b.tar.gz
cp -f  linux/ifort_xeon_phi/make.CHOICES sprng2.0
cp  linux/ifort_xeon_phi/make.IFORT64 sprng2.0/SRC
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
export CC=icc
export CFLAGS=-mmic
./configure
\cp -f ../Makefile_cfitsio.xeon_phi Makefile
make
cp libcfitsio.a ../lib
cd ..
rm -rf cfitsio

# Numerical recipes
mkdir lib/nr lib/nr/eq_diff lib/nr/spline lib/nr/sort
cd nr
./compile_xeon_phi.com
cp libnr.a *.mod ../lib/nr
cp eq_diff/libnr_eq_diff.a eq_diff/*.mod ../lib/nr/eq_diff
cp spline/libnr_splin.a ../lib/nr/spline
cp sort/libnr_sort.a ../lib/nr/sort
./clean.com
cd ..

cp -r include $MCFOST_INSTALL
mkdir $MCFOST_INSTALL/lib
cp -r lib $MCFOST_INSTALL/lib/xeon_phi
