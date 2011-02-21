xlf90 -qsuffix=f=f90 -O3 -qstrict -qarch=auto -qtune=auto -c  nrtype.f90 nr.f90 nrutil.f90
ar rv libnr.a nrtype.o nr.o nrutil.o
ranlib libnr.a	
cd eq_diff
xlf90  -qsuffix=f=f90 -O3 -qstrict -qarch=auto -qtune=auto  -c  odeint.f90  rk4.f90  rkck.f90 rkdumb.f90  rkqs.f90 -I ..
ar rv libnr_eq.a *.o
mv libnr_eq.a libnr_eq_diff.a
ranlib libnr_eq_diff.a
cd ../generateur
xlf90  -qsuffix=f=f90 -O3 -qstrict -qarch=auto -qtune=auto  -c   ran_state.f90 ran0.f90 ran1.f90 ran2.f90 -I ..
ar rv libnr_rand.a *.o
mv libnr_rand.a libnr_random.a
ranlib libnr_random.a
cd ../spline
xlf90 -qsuffix=f=f90 -O3 -qstrict -qarch=auto -qtune=auto -c spline.f90 splint.f90 tridag.f90 locate.f90  -I .. 
ar rv libnr_splin.a spline.o splint.o tridag.o locate.o  
ranlib libnr_splin.a
cd ../sort
xlf90 -qsuffix=f=f90 -O3 -qstrict -qarch=auto -qtune=auto -c indexx.f90 -I ..
ar rv  libnr_sort.a *.o
ranlib libnr_sort.a
