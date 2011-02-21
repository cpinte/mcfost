rm -rf *.o *.mod *.a
rm -rf eq_diff/*.o eq_diff/*.mod eq_diff/*.a
rm -rf generateur/*.o generateur/*.mod  generateur/*.a
rm -rf spline/*.o spline/*.mod  spline/*.a
rm -rf sort/*.o sort/*.mod  sort/*.a
