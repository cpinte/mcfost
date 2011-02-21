sunf95 -c  nrtype.f90 nr.f90 nrutil.f90
ar rv libnr.a *.o
ranlib libnr.a	
cd eq_diff
sunf95 -c  odeint.f90  rk4.f90  rkck.f90 rkdumb.f90  rkqs.f90 -M..
ar rv libnr_eq_diff.a *.o
ranlib libnr_eq_diff.a
cd ../generateur
sunf95 -c   ran_state.f90 ran0.f90 ran1.f90 ran2.f90 -M..
ar rv libnr_random.a *.o
ranlib libnr_random.a
cd ../spline
sunf95 -c  spline.f90 splint.f90 locate.f90 tridag.f90 -M..
ar rv  libnr_splin.a *.o
ranlib libnr_splin.a
cd ../sort
sunf95 -c indexx.f90  -M..
ar rv  libnr_sort.a *.o
ranlib libnr_sort.a
