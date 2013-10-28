program BIGCRUNCH
  !***********************************************************
  ! Code de transfert radiatif Monte-Carlo : MCFOST
  !***********************************************************
  ! Transfert dans la poussiere : lumiere diffusee et emission
  ! thermique et dans le gaz : raies moleculaires
  ! Code parallèle OpenMP
  !
  ! F. Menard, G. Duchene et C. Pinte
  ! Grenoble, Exeter, Grenoble, Santiago
  !
  !
  ! Quelques remarques :
  ! - Vive le 64bits !!
  ! - il faut augmenter la stack size sur les architectures x86
  ! et AIX pour les lib numerical recipes. Passage en
  ! allocation dynamique dans MCFOST en cas de pb
  ! - pas de pointeur en OpenMP sur Sun x86_64
  ! - ifort n'aime pas le retypage a la volee des tableaux
  !
  !***********************************************************
  use parametres
  use init_mcfost
  use dust_transfer
  use mol_transfer

  implicit none

  integer :: time
  real :: cpu_time_begin, cpu_time_end

  ! debut de l'execution
  call system_clock(time_begin,count_rate=time_tick,count_max=time_max)
  call cpu_time(cpu_time_begin)

  call initialisation_mcfost()

  !*************************************************************************
  ! Here it comes ...
  !*************************************************************************


  ! Transfert radiatif dans le continu
  ldust_transfer = .true.
  if (ldust_transfer) then
     call transfert_poussiere()
  endif

  ! Emission moleculaire ...
  if (lemission_mol) then
     call mol_line_transfer()
  endif

  ! Temps d'execution
  call system_clock(time_end)
  if (time_end < time_begin) then
     time=int((time_end + (1.0 * time_max)- time_begin)/time_tick)
  else
     time=int((time_end - time_begin)/time_tick)
  endif
  write (*,'(" Processing complete in ", I3, "h", I3, "m", I3, "s")')  time/3600, mod(time/60,60), mod(time,60)
  call cpu_time(cpu_time_end)
  time = int(cpu_time_end - cpu_time_begin)
  write (*,'(" CPU time used          ", I3, "h", I3, "m", I3, "s")')  time/3600, mod(time/60,60), mod(time,60)


end program BIGCRUNCH
