program BIGCRUNCH
  !***********************************************************
  ! Code de transfert radiatif Monte-Carlo : MCFOST
  !***********************************************************
  ! Transfert dans la poussiere : lumiere diffusee et emission
  ! thermique et dans le gaz : raies moleculaires
  ! Code parall�le OpenMP
  !
  ! F. Menard, G. Duchene et C. Pinte
  ! Grenoble, Exeter, Grenoble, Santiago
  !
  !***********************************************************
  use parametres
  use init_mcfost
  use dust_transfer
  use mol_transfer

  implicit none

  integer :: itime
  real :: time, cpu_time_begin, cpu_time_end

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
     time=(time_end + (1.0 * time_max)- time_begin)/real(time_tick)
  else
     time=(time_end - time_begin)/real(time_tick)
  endif
  if (time > 60) then
     itime = real(time)
     write (*,'(" Processing complete in ", I3, "h", I3, "m", I3, "s")')  itime/3600, mod(itime/60,60), mod(itime,60)
  else
     itime = int(time)
     write (*,'(" Processing complete in ", F5.2, "s")')  time
  endif
  call cpu_time(cpu_time_end)
  time = cpu_time_end - cpu_time_begin
  if (time > 60) then
     itime = int(time)
     write (*,'(" CPU time used          ", I3, "h", I3, "m", I3, "s")')  itime/3600, mod(itime/60,60), mod(itime,60)
  else
     write (*,'(" CPU time used          ", F5.2, "s")')  time
  endif

end program BIGCRUNCH
