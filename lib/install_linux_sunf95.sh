mkdir lib
mkdir include

# SPRNG
# gcc-2.95 needed
#wget http://sprng.cs.fsu.edu/Version2.0/sprng2.0a.tgz
tar xzvf sprng2.0a.tgz
cp -f  linux/sunf95/make.CHOICES sprng2.0
cp  linux/sunf95/make.SUNF95 sprng2.0/SRC
cd sprng2.0
make
mv lib/libsprng.a ../lib
mv include/*.h ../include
cd ..
rm -rf sprng2.0

# cfitsio
# g77 ou f77 needed by configure to set up the fotran wrapper in Makefile 
#wget ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3006.tar.gz
tar xzvf cfitsio3006.tar.gz
cd cfitsio
./configure
make
cp libcfitsio.a ../lib
cd ..
rm -rf cfitsio

# appels systemes
cd appels_systeme/linux/sunf95/
make
cp libappel_syst.a ../../../lib
make clean
cd ../../..

# Numerical recipes
mkdir lib/nr lib/nr/eq_diff lib/nr/spline
cd nr
./compile_sunf95.com
cp libnr.a *.mod ../lib/nr
cp eq_diff/libnr_eq_diff.a eq_diff/*.mod ../lib/nr/eq_diff
cp spline/libnr_splin.a ../lib/nr/spline
./clean.com
cd ..

cp -r include $MCFOST_INSTALL
mkdir $MCFOST_INSTALL/lib
cp -r lib $MCFOST_INSTALL/lib/sunf95