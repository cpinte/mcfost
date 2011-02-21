ifort -c -w -O3 -tpp6 -mp -pc80 -prec_div nrtype.f90 nr.f90 nrutil.f90
ar rv libnr.a *.o
ranlib libnr.a	
cd eq_diff
ifort -c -w -O3 -tpp6 -mp -pc80 -prec_div odeint.f90  rk4.f90  rkck.f90 rkdumb.f90  rkqs.f90 -I ..
ar rv libnr_eq_diff.a *.o
ranlib libnr_eq_diff.a
cd ../generateur
ifort -c -w -O3 -tpp6 -mp -pc80 -prec_div  ran_state.f90 ran0.f90 ran1.f90 ran2.f90 -I ..
ar rv libnr_random.a *.o
ranlib libnr_random.a
cd ../spline
ifort -c -w -O3 -tpp6 -mp -pc80 -prec_div spline.f90 splint.f90 locate.f90 tridag.f90 -I ..
ar rv  libnr_splin.a *.o
ranlib libnr_splin.a
cd ../sort
ifort -c -w -O3 -tpp6 -mp -pc80 -prec_div indexx.f90 -I ..
ar rv  libnr_sort.a *.o
ranlib libnr_sort.a
cd ../LU
ifort -c -w -O3 -tpp6 -mp -pc80 -prec_div lubksb.f90  ludcmp.f90 -I ..
ar rv libnr_LU.a *.o
ranlib libnr_LU.a