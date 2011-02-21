f95 -c  -O3 -C  nrtype.f90 nr.f90 nrutil.f90 -xarch=amd64
ar rv libnr.a *.o
ranlib libnr.a
cd eq_diff
f95 -c -O3  odeint.f90  rk4.f90  rkck.f90 rkdumb.f90  rkqs.f90 -M.. -xarch=amd64
ar rv libnr_eq_diff.a *.o
ranlib libnr_eq_diff.a
cd ../generateur
f95 -c -O3  ran_state.f90 ran0.f90 ran1.f90 ran2.f90 -M.. -xarch=amd64
ar rv libnr_random.a *.o
ranlib libnr_random.a
cd ../spline
f95 -c -O3  spline.f90 splint.f90 locate.f90 tridag.f90 -M.. -xarch=amd64
ar rv  libnr_splin.a *.o
ranlib libnr_splin.a
cd ../sort
f95 -c -O3  indexx.f90 -M.. -xarch=amd64
ar rv  libnr_sort.a *.o
ranlib libnr_sort.a
