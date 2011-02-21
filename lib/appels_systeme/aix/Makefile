main : 
	xlf90 -qsuffix=f=f90 -c appel_syst_aix.f90
	ar rv libappel_syst.a appel_syst_aix.o
	ranlib libappel_syst.a

clean : 
	rm -f appel_syst_aix.o libappel_syst.a