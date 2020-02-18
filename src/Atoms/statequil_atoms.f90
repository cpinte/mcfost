MODULE statequil_atoms

	use atmos_type, only					: atmos, nTotal_atom
	use atom_type
	use spectrum_type, only					: NLTEspec
	use constant
	use opacity
	use math, only							: locate, any_nan_infinity_matrix, any_nan_infinity_vector, is_nan_infinity, integrate_x
	use parametres
	use accelerate
	use collision, only						: CollisionRate !future deprecation
	use impact, only						: Collision_Hydrogen
	use getlambda, only						: hv
 
	use mcfost_env, only					: dp
	use constantes, only					: tiny_dp

	IMPLICIT NONE
 
	real(kind=dp), parameter :: prec_pops = 1d-100 !1d-8
	character(len=15), parameter :: invpop_file = "inversion_populations.txt"

	CONTAINS

 
	SUBROUTINE fill_collision_matrix(id, icell)
		integer, intent(in) :: icell, id
		type (AtomType), pointer :: atom
		integer :: nact
  
		do nact=1, atmos%NactiveAtoms
			atom => atmos%ActiveAtoms(nact)%ptr_atom
    
			CALL collision_matrix_atom(id, icell, atom)
  	 
			atom => NULL()
		enddo
  
	RETURN 
	END SUBROUTINE fill_collision_matrix


 !Occupa prob ?
 !The matrix itself could kept in memory, and we would just have to multiply it by ne
	SUBROUTINE collision_matrix_atom(id, icell, atom)
		integer, intent(in) :: icell, id
		type (AtomType), pointer, intent(inout) :: atom
  

    
		atom%C(:,:,id) = 0d0
		if (atom%ID=="H") then
			atom%C(:,:,id) = Collision_Hydrogen(icell)
		else
! 			write(*,*) "Collision RH not set for large atoms"
! 			stop
! 			if (id==1) write(*,*) " Using RH routine for collision atom ", atom%ID
			atom%C(:,:,id) = CollisionRate(icell, atom)
		end if
     
		if (any_nan_infinity_matrix(atom%C(:,:,id))>0) then
			write(*,*) atom%C(:,:,id)
			write(*,*) "Bug at compute_collision_matrix", id, icell, atom%ID, atom%n(:,icell)
			stop
		endif

  
	RETURN 
	END SUBROUTINE collision_matrix_atom
 
	SUBROUTINE initGamma_atom(id, atom)
		integer, intent(in) :: id
		type (AtomType), pointer, intent(inout) :: atom
		integer :: kc
  
		atom%Gamma(:,:,id) = atom%C(:,:,id)
  
	RETURN 
	END SUBROUTINE initGamma_atom
 
	SUBROUTINE initGamma(id)
	! ------------------------------------------------------ !
	!for each active atom, allocate and init Gamma matrix
	!deallocated in freeAtom()
	! Collision matrix has to be allocated
	! if already allocated, simply set Gamma to its new icell
	! value.
	!
	! n(l)*Cl->lp = n(lp)*Clp->l
	! e.g.,
	! Cij = Cji * (nj/ni)^star, with Cji = C(ij) = 
	! colision rate from j to i.
	!
	! (ij = j->i = (i-1)*N + j)
	! (ji = i->j = (j-1)*N + i)
	! ------------------------------------------------------ !
		integer :: nact, Nlevel, lp, l, ij, ji, nati, natf
		integer, intent(in) :: id
		type (AtomType), pointer :: atom
  
		nati = 1; natf = atmos%Nactiveatoms
		do nact = nati, natf
			atom => atmos%ActiveAtoms(nact)%ptr_atom

			CALL initGamma_atom(id, atom)

			NULLIFY(atom)
		enddo
	RETURN
	END SUBROUTINE initGamma

 
	SUBROUTINE remove_delta(id, atom)
	!Remove therm in delta(l,l')
	!need to be done once Gamma is know, so not in FillGamma_mu which is the 
	! wavelength-ray building of Gamma.
	integer, intent(in) :: id
	type (AtomType), intent(inout) :: atom
	integer :: l
  
	do l = 1, atom%Nlevel
		atom%Gamma(l,l,id) = 0.0_dp
		atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
	end do
  
	RETURN
	END SUBROUTINE remove_delta

 
	SUBROUTINE calc_rates(id, icell, iray, n_rayons)
		integer, intent(in) :: id, icell, iray, n_rayons
		integer :: nact
		
		do nact=1, atmos%NactiveAtoms
		
			call calc_rates_atom(id, icell, iray, atmos%ActiveAtoms(nact)%ptr_atom, n_rayons)
		
		enddo
	
	RETURN
	END SUBROUTINE 
	
	SUBROUTINE init_rates_atom(id, atom)
		integer, intent(in) :: id
		type (AtomType), intent(inout) :: atom
		integer :: k
		
		do k=1, atom%Nline
		
			!atom%lines(k)%Jbar(id) = 0.0
			atom%lines(k)%Rij(id) = 0.0
			atom%lines(k)%Rji(id) = atom%lines(k)%Aji
		
		enddo
		
		do k=1, atom%Ncont
		
			atom%continua(k)%Rij(id) = 0.0
			atom%continua(k)%Rji(id) = 0.0
		
		enddo
		
	
	RETURN
	END SUBROUTINE init_rates_atom
 
 	!occupa prob
	SUBROUTINE calc_rates_atom(id, icell, iray, atom, n_rayons)

		integer, intent(in) :: icell, id, n_rayons, iray
		type (AtomType), intent(inout) :: atom
		integer :: kr, kc, i, j, l, Nblue, Nred, Nl
		real(kind=dp) :: Ieff, Jnu_ray, aJnu_ray, etau, wphi
		integer :: dk
		integer, parameter :: dl = 5 !Do continuum integral dl by dl bins 
								!In this case, integrate from 1, Nl+dl,dl to take the last point !

  
		tr_loop : do kr=1, atom%Ntr
	 		kc = atom%at(kr)%ik
   
			SELECT CASE (atom%at(kr)%trtype)
   
			CASE ("ATOMIC_LINE") !Rji initialized to Aji !!
				i = atom%lines(kc)%i; j = atom%lines(kc)%j
				Nblue = atom%lines(kc)%Nblue; Nred = atom%lines(kc)%Nred
				Nl = atom%lines(kc)%Nlambda
				dk = atom%lines(kc)%dk(iray,id)
				
				!Integration over frequencies for this direction
				!Trapezoidal rule ? Here included in the weight, except if we use hv
				!and in this case we integrate I*phi*dv (rectangle)
				
				Jnu_ray = 0.0
				wphi = 0.0
				!Do we need line%w_lam, since the line is linearly sampled in dv ?
				do l=1, Nl !Fast relatively
				
					wphi = wphi + atom%lines(kc)%phi(l,icell)*1e3*hv!*atom%lines(kc)%w_lam(l)
					
					!-> to check
					!write(*,*) "wl = ", atom%lines(kc)%w_lam(l)*1e-3

!-> Only for MALI scheme
! 					Ieff = NLTEspec%I(Nblue+l-1-dk,iray,id)*NLTEspec%etau(Nblue+l-1-dk,iray,id) + &
!        		 NLTEspec%Psi(Nblue+l-1-dk,iray,id) * atom%eta(Nblue+l-1-dk,iray,id)
       		 
!        		 		Ieff = NLTEspec%I(Nblue+l-1-dk,iray,id) * NLTEspec%etau(Nblue+l-1-dk,iray,id) + &
!        		 		NLTEspec%Psi(Nblue+l-1-dk,iray,id) * NLTEspec%S(Nblue+l-1-dk,iray,id)


					etau = dexp(-ds(iray,id) * NLTEspec%chi(Nblue+l-1-dk,iray,id))
       		 		Ieff = NLTEspec%I(Nblue+l-1-dk,iray,id) * etau  + (1.0_dp - etau) * NLTEspec%S(Nblue+l-1-dk,iray,id)
      		 
					Jnu_ray = Jnu_ray + Ieff * atom%lines(kc)%phi(l,icell)*(1e3 * hv)!*atom%lines(kc)%w_lam(l)
					
				enddo
				!Renormalise profile ?
! 				write(*,*) wphi
				if (wphi < 0.95) then
					write(*,*) "icell = ", icell, " id = ", id
					write(*,*) " --> Beware, profile not well normalized for line ", i, j, " area = ", wphi
				endif
     			Jnu_ray = Jnu_ray / n_rayons / wphi

				!Integration over directions
				
				!atom%lines(kc)%Jbar(id) = atom%lines(kc)%Jbar(id) + Jnu_ray
				atom%lines(kc)%Rji(id) = atom%lines(kc)%Rji(id) + Jnu_ray * atom%lines(kc)%Bji
				atom%lines(kc)%Rij(id) = atom%lines(kc)%Rij(id) + Jnu_ray * atom%lines(kc)%Bij
				
! write(*,*) "id = ", id, " icell = ", icell
! write(*,*) " Ray = ", iray
! if (i==1 .and. j==4) &				
! 				write(*,*) "line ", "Aji = ", atom%lines(kc)%Aji, 'Rij = ', atom%lines(kc)%Rij(id), 'Rji = ', atom%lines(kc)%Rji(id)

        
			CASE ("ATOMIC_CONTINUUM")
				i = atom%continua(kc)%i; j = atom%continua(kc)%j
				Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
				Nl = atom%continua(kc)%Nlambda
! if (i==1 .and. j==4) &
! 				write(*,*) "Nlambda = ", Nl, " iray = ", iray
				
				!Integration over frequencies for this direction
				!Trapezoidal rule, in the weight definition.
				
				
				Jnu_ray = 0.0
				aJnu_ray = 0.0
				do l=1, Nl !Can takes time
    
! 					Ieff = NLTEspec%I(Nblue+l-1,iray,id)*NLTEspec%etau(Nblue+l-1,iray,id) + &
!        				 NLTEspec%Psi(Nblue+l-1,iray,id) * atom%eta(Nblue+l-1,iray,id)
       				 
! 					Ieff = NLTEspec%I(Nblue+l-1,iray,id) * NLTEspec%etau(Nblue+l-1,iray,id) + &
! 					NLTEspec%Psi(Nblue+l-1,iray,id) * NLTEspec%S(Nblue+l-1,iray,id)

					etau = dexp(-ds(iray,id) * NLTEspec%chi(Nblue+l-1,iray,id))
					Ieff = NLTEspec%I(Nblue+l-1,iray,id) * etau + (1.0_dp - etau) * NLTEspec%S(Nblue+l-1,iray,id)
       		 
					aJnu_ray = aJnu_ray + Ieff * atom%continua(kc)%alpha(l)*atom%continua(kc)%w_lam(l)

					Jnu_ray = Jnu_ray + (Ieff  + &
					atom%continua(kc)%twohnu3_c2(l)) * atom%continua(kc)%gij(l,icell) * atom%continua(kc)%alpha(l)*atom%continua(kc)%w_lam(l)
! if (i==1 .and. j==4) then					
! 					write(*,*) "i=",i, "j=",j, l, "lam=", NLTEspec%lambda(nblue+l-1)
! 					write(*,*) '2hnu3_c2 = ', atom%continua(kc)%twohnu3_c2(l), " gij = ", atom%continua(kc)%gij(l,icell)
! 					write(*,*) "alpha = ", atom%continua(kc)%alpha(l), " Ieff = ", Ieff
! 					write(*,*) " Wlam = ", atom%continua(kc)%w_lam(l), " eta = ", atom%eta(Nblue+l-1,iray,id)
! 					write(*,*) " exp(-dtau) = ", NLTEspec%etau(Nblue+l-1,iray,id), " Psi = ", NLTEspec%psi(Nblue+l-1,iray,id)
! 					write(*,*) "aJ = ", aJnu_ray * fourPI_h / n_rayons, Jnu_ray * fourPI_h / n_rayons
! 
! endif

				enddo

				!Integration over directions

				atom%continua(kc)%Rij(id) = atom%continua(kc)%Rij(id) + fourpi_h * aJnu_ray / n_rayons
				atom%continua(kc)%Rji(id) = atom%continua(kc)%Rji(id) + fourpi_h * Jnu_ray / n_rayons

! if (i==1 .and. j==4) &
! 				write(*,*) "Cont ", 'Rij = ', atom%continua(kc)%Rij(id), 'Rji = ', atom%continua(kc)%Rji(id)

! if (iray==n_rayons) then
! 	write(*,*) " Integration over frequencies and directions done", id, icell
! 	stop
! endif
			CASE DEFAULT
    
				CALL Error("Unkown transition type", atom%at(kr)%trtype)
     
			END SELECT
  
		end do tr_loop

 
	RETURN
	END SUBROUTINE calc_rates_atom

	!occupa prob
	SUBROUTINE rate_matrix_atom(id, icell, atom, switch_to_lte)
	!Gamma init to collision rates
		integer, intent(in) :: icell, id
		type (AtomType), intent(inout) :: atom
		logical, optional, intent(in) :: switch_to_lte
		integer :: kr, kc, i, j, l, Nblue, Nred, Nl

		if (present(switch_to_lte)) then
			if (switch_to_lte) return
		end if 
  
		tr_loop : do kr=1, atom%Ntr
			kc = atom%at(kr)%ik
   
			SELECT CASE (atom%at(kr)%trtype)
   
			CASE ("ATOMIC_LINE")
				i = atom%lines(kc)%i; j = atom%lines(kc)%j

! write(*,*) atom%Gamma(:,:,id) - atom%C(:,:,id)
! stop	
				atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + atom%lines(kc)%Rji(id)
				atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + atom%lines(kc)%Rij(id)
				     
     !remove diagonal
        
			CASE ("ATOMIC_CONTINUUM")


				i = atom%continua(kc)%i; j = atom%continua(kc)%j


				atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + atom%continua(kc)%Rij(id)
				atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + atom%continua(kc)%Rji(id)
				
       
			CASE DEFAULT
    
				CALL Error("Unkown transition type", atom%at(kr)%trtype)
     
			END SELECT
  
		end do tr_loop

 
	RETURN
	END SUBROUTINE rate_matrix_atom

 
 !also remove delta term. Has to be done after Gamma is known 
 !after wavelength-ray integration
	SUBROUTINE SEE_atom(id, icell,atom, dM)
	! --------------------------------------------------------------------!
	! For atom atom solves for the Statistical Equilibrium Equations (SEE) 
	!
	! We solve for :
	!  Sum_l' Gamma_l'l n_l' = 0 (m^-3 s^-1)
	!
	! which is equivalent in matrix notation to:
	!
	! GG(l,l') dot n_l' = 0 with l' is on the columns and l on the rows
	!
	! In particular for a 2-level atom the two equations are:
	!
	! n1*G_11 + n2*G_21 = 0
	! n1*G12 + n2*G22 = 0,
	! with one of these equations has to be replaced by
	! n1 + n2 = N
	!
	! For numerical stability, the row with the largest populations is 
	! removed (i.e., it is replaced by the mass conservation law).
	!
	! This matrix Gamma for atoms (or later molecules) is compared to 
	! molecules A matrix in mcfost:
	! diag(Gamma) = -diag(A); offdiag(Gamma) = -offdiag(A)
	! A = transpose(Gamma)
	! --------------------------------------------------------------------!
  
		integer, intent(in) :: icell, id
		type(AtomType), intent(inout) :: atom
		integer :: lp, imaxpop, l
		real(kind=dp), intent(out) :: dM
		real(kind=dp), dimension(atom%Nlevel) :: ndag
		real(kind=dp) :: n0
		real(kind=dp), dimension(atom%Nlevel, atom%Nlevel) :: Aij
		
		n0 = 0.0_dp
		if (atom%ID=="H") n0 = atmos%nHmin(icell)
		if (n0 >= ntotal_atom(icell,atom)) call error("Error n0")
  

  !Gamma is known = wavelength and ray integrated
		CALL remove_delta(id, atom)
  !store old populations
		ndag(:) = atom%n(:,icell) /( ntotal_atom(icell,atom) - n0 )
  
		imaxpop = locate(atom%n(:,icell), maxval(atom%n(:,icell)))
		!imaxpop = atom%Nlevel
  
		atom%n(:,icell) = 0d0
  !use normalised now
		atom%n(imaxpop,icell) = 1d0 !(nTotal_atom(icell, atom)-n0)
  
  !Sum_l'_imaxpop * n_l' = N
		atom%Gamma(:,imaxpop,id) = 1d0 !all columns of the last row for instance
  !(G11 G21)  (n1)  (0)
  !(       ) .(  ) =( )
  !(1    1 )  (n2)  (N)
  
  !Y a peut être un transpose ici  par rapport à MCFOST, atom%Gamma.T ?
		Aij = transpose(atom%Gamma(:,:,id))

		if ((any_nan_infinity_matrix(atom%Gamma(:,:,id))>0)) then
			write(*,*) "BUG Gamma", id, icell 
			write(*,*) atom%Gamma(:,:,id)
			write(*,*) "n = ", atom%n(:,icell)
			write(*,*) "ndag=", ndag
			stop
		end if
  !write(*,*) atom%ID, " Rate matrix (G.T)", Aij
		CALL GaussSlv(Aij, atom%n(:,icell),atom%Nlevel)
		if ((any_nan_infinity_vector(atom%n(:,icell))>0)) then
			write(*,*) "BUG pops", id, icell, atom%n(:,icell)
			stop
		end if
  
		dM = 0.0_dp
		do l=1,atom%Nlevel
   !if (atom%n(l,icell) < prec_pops) then
			if (atom%n(l,icell)/maxval(atom%n(:,icell)) < prec_pops) then
    !n(l)/ntot < prec_pops
				write(*,*) icell, atom%n(l,icell), " minor populations level ", l, " of atom ", atom%ID
				atom%n(l,icell) = 0.0_dp
			else !denormalise
				dM = max(dM, dabs(atom%n(l,icell) - ndag(l))/atom%n(l,icell))
				atom%n(l,icell) = atom%n(l,icell) * ( ntotal_atom(icell,atom) - n0 )
			endif
		enddo
  
! 		write(*,*) " ***** Inverted pops *****"
! 		write(*,*) id, icell, atom%ID
! 		write(*,*) (atom%n(l,icell)/ntotal_atom(icell,atom) ,l=1,atom%Nlevel)
! 		write(*,*) (atom%n(l,icell),l=1,atom%Nlevel)
! 		write(*,*) " ***** ************* *****"

	RETURN
	END SUBROUTINE SEE_atom
 
	!occupa prob
 	!Define a negative Tex for inversion of Temperature ? because the ratio of ni/nj exists
 	!even with inversion of pops.
 	!For negative Tex (or Tion), the old value is used, and the transition does not contribute to
 	!the global convergence of the atom.
 	!Take absolute value instead ??
	SUBROUTINE calc_delta_Tex_atom(icell, atom, dT, Tex, Tion)
		integer, intent(in) :: icell
		type(AtomType), intent(inout) :: atom
		real(kind=dp), intent(out) :: dT, Tex, Tion
		integer :: nact, kr, kc, i, j
		real(kind=dp) :: deltaE_k, ratio, Tdag, gij, dT_loc, Tion_loc, Tex_loc
		real(kind=dp) :: wi, wj
		real :: sign
		
  
		dT = 0.0_dp !for all transitions of one atom
		Tex = 0.0
		Tion = 0.0
		
		tr_loop : do kr=1, atom%Ntr
			kc = atom%at(kr)%ik
   
			SELECT CASE (atom%at(kr)%trtype)
   
			CASE ("ATOMIC_LINE")
				i = atom%lines(kc)%i; j = atom%lines(kc)%j
				if (atom%ID=="H") then! .or. atom%ID=="He") then
					wi = 1.0
					wj = 1.0
				else
					wi = 1.0
					wj = 1.0
				endif


				Tdag = atom%lines(kc)%Tex(icell)
				deltaE_k = (atom%E(j)-atom%E(i)) / KBOLTZMANN
                
				ratio = dlog(wi * atom%n(j,icell) * atom%lines(kc)%gij / (atom%n(i,icell)*wj))
				

				if (ratio /= 0.0) then
					sign = real(ratio / dabs(ratio))
      
					atom%lines(kc)%Tex(icell) = -deltaE_k / ratio
					Tex_loc = atom%lines(kc)%Tex(icell)
       
					if (atom%lines(kc)%Tex(icell) < 0) then
						write(*,*) "Tex negative ( njgij > ni) ", wi * atom%n(j,icell) * atom%lines(kc)%gij, atom%n(i,icell)*wj
						write(*,*) "ratio = ", ratio, " diff = ", (atom%n(i,icell) - atom%n(j,icell) * atom%lines(kc)%gij)/(1d-50 + atom%n(i,icell))
						write(*,*) "icell = ", icell, " :: Te = ", atmos%T(icell)
						write(*,*) " --> Setting to old value", Tdag, atom%ID, " line (i,j) = ", i, j
						atom%lines(kc)%Tex(icell) = Tdag
						cycle tr_loop !do not count this line for max relative change
        				!stop
       				endif
       
       				dT_loc = dabs(Tdag-Tex_loc)/(1d-50 + Tex_loc)
       				dT = max(dT, dT_loc)
       				Tex = max(Tex, Tex_loc)

     			 endif


        
    		CASE ("ATOMIC_CONTINUUM")

     			i = atom%continua(kc)%i; j = atom%continua(kc)%j
				wi = 1.0
           
      			Tdag = atom%continua(kc)%Tex(icell)

				deltaE_k = (atom%E(j)-atom%E(i)) / KBOLTZMANN
				!ratio = dlog( atom%n(j,icell)*atom%g(i) / ( atom%n(i,icell)*atom%g(j) ) )
				!Doesnt make sens for continua Tex = -deltaE_k / ratio
      
				!at threshold
				!i.e., deltaE is hnu0
				gij = atom%nstar(i,icell)/(1d-50 + atom%nstar(j,icell) ) * dexp(-hc_k/atom%continua(kc)%lambda0/atmos%T(icell))
				ratio = log( atom%n(i,icell)  / ( atom%n(j,icell) * gij ) )
					
				if (ratio /= 0.0) then
					sign = real(ratio / dabs(ratio))
					!ionisation temperature
					atom%continua(kc)%Tex(icell) = deltaE_k / ratio
					Tion_loc = atom%continua(kc)%Tex(icell)
					
					if (atom%continua(kc)%Tex(icell) < 0) then
						write(*,*) "Tion negative (njgij > ni) ", wi * atom%n(j,icell) * gij, atom%n(i,icell)*wj
						write(*,*) "ratio = ", ratio, " diff = ", (atom%n(i,icell) - atom%n(j,icell) * gij)/(1d-50 + atom%n(i,icell))
						write(*,*) "icell = ", icell, " :: Te = ", atmos%T(icell)
						write(*,*) " --> Setting to old value", Tdag, atom%ID, " cont (i,j) = ", i, j
						atom%continua(kc)%Tex(icell) = Tdag
						cycle tr_loop !do not count this cont for max relative change
        				!stop
       				endif
       				
       				dT_loc = dabs(Tdag-Tion_loc)/(Tion_loc + 1d-50)
       				dT = max(dT, dT_loc)		
       				Tion = max(Tion_loc, Tion)			

				endif
    
 
    		CASE DEFAULT
    		
     			CALL Error("Unkown transition type", atom%at(kr)%trtype)
     
  			 END SELECT
  
  		end do tr_loop
 
	RETURN
	END SUBROUTINE calc_delta_Tex_atom
 
 	!occupa prob
 	!Define nefative Tex ? 
	SUBROUTINE calc_Tex_atom(icell, atom)
	!
	! For lines:
	!
	! n(j)gij / n(i) = exp(-dE / kTex)
	!
	! -> Tex = -dE/k * (log(nj*gij) - log(ni))**-1
	!
	!
	! For continua I define %Tex = Tion = hnu/k * 1 / log(2hnu3/c2/Snu_cont + 1), 
	! the ionisation temperature. If Tion=Tle, Snu_cont(T=Tlte) = Bnu(Tlte) (gij = f(Tlte))
	!
		integer, intent(in) :: icell
		type(AtomType), intent(inout) :: atom
		integer :: nact, kr, kc, i, j
		real(kind=dp) :: deltaE_k, ratio, Tex, gij, wi, wj
		real :: sign
		
		wi = 1.0
		wj = 1.0
  
		tr_loop : do kr=1, atom%Ntr
			kc = atom%at(kr)%ik
   
			SELECT CASE (atom%at(kr)%trtype)
   
			CASE ("ATOMIC_LINE")
				i = atom%lines(kc)%i; j = atom%lines(kc)%j
				
				deltaE_k = (atom%E(j)-atom%E(i)) / KBOLTZMANN
				ratio = dlog(atom%n(j,icell) * atom%lines(kc)%gij) - dlog(atom%n(i,icell))
				if (ratio /= 0.0) then
					!sign negative means positive Tex, not included yet
					sign = real(ratio/abs(ratio))

					!should be de = hnu0
					atom%lines(kc)%Tex(icell) = -deltaE_k / ratio
					if (atom%lines(kc)%Tex(icell) < 0) then
						write(*,*) ratio, "Tex negative (njgij > ni) ", wi * atom%n(j,icell) * atom%lines(kc)%gij, atom%n(i,icell)*wj
						write(*,*) "icell = ", icell, " :: Te = ", atmos%T(icell)
						write(*,*) " --> Setting to abs()", atom%ID, " line (i,j) = ", i, j
						atom%lines(kc)%Tex(icell) = abs(atom%lines(kc)%Tex(icell))
       				endif
				endif
        
			CASE ("ATOMIC_CONTINUUM")
				i = atom%continua(kc)%i; j = atom%continua(kc)%j
     
     
				deltaE_k = (atom%E(j)-atom%E(i)) / KBOLTZMANN
				!ratio = dlog( atom%n(j,icell)*atom%g(i) / ( atom%n(i,icell)*atom%g(j) ) )
				!Doesn't make sens for continua Tex = -deltaE_k / ratio
      
				!at threshold
				!i.e., deltaE is hnu0
				gij = atom%nstar(i,icell)/(1d-50 + atom%nstar(j,icell) ) * dexp(-hc_k/atom%continua(kc)%lambda0/atmos%T(icell))
				ratio = log( atom%n(i,icell)  / ( atom%n(j,icell) * gij ) )
					
				if (ratio /= 0.0) then
					sign = real(ratio / dabs(ratio))
					!ionisation temperature
					atom%continua(kc)%Tex(icell) = deltaE_k / ratio
					if (atom%lines(kc)%Tex(icell) < 0) then
						write(*,*) ratio, "Tion negative (njgij > ni) ", wi * atom%n(j,icell) * gij, atom%n(i,icell)*wj
						write(*,*) "icell = ", icell, " :: Te = ", atmos%T(icell)
						write(*,*) " --> Setting to abs()", atom%ID, " cont (i,j) = ", i, j
						atom%continua(kc)%Tex(icell) = abs(atom%continua(kc)%Tex(icell))
       				endif
				endif
    
 
			CASE DEFAULT
    
				CALL Error("Unkown transition type", atom%at(kr)%trtype)
     
			END SELECT
  
		end do tr_loop

 
	RETURN
	END SUBROUTINE calc_Tex_atom
 
	SUBROUTINE calc_Tex(icell) !for all atoms
		integer, intent(in) :: icell
		type(AtomType), pointer :: atom
		integer :: nact
  
		do nact=1, atmos%NactiveAtoms
			atom => atmos%ActiveAtoms(nact)%ptr_atom
   
			CALL calc_Tex_atom(icell, atom)
  
			atom => NULL()
		enddo
 
	RETURN
	END SUBROUTINE calc_Tex
 
	SUBROUTINE update_populations(id, icell, delta,verbose, nit)
	! --------------------------------------------------------------------!
	! Performs a solution of SEE for each atom.
	!
	! to do: implements Ng's acceleration iteration here 
	! --------------------------------------------------------------------!

		integer, intent(in) :: id, icell, nit
		logical, intent(in) :: verbose
		type(AtomType), pointer :: atom
		integer :: nat, nati, natf
		logical :: accelerate = .false.
		real(kind=dp) :: dM, dT, dTex, dpop, Tion, Tex
		real(kind=dp), intent(out) :: delta
  
		dpop = 0.0_dp
		dTex = 0.0_dp
		
		if (verbose) write(*,*) " --> niter = ", nit, " id = ", id, " icell = ", icell

		nati = 1; natf = atmos%Nactiveatoms
		do nat=nati,natf !loop over each active atoms
  
			atom => atmos%ActiveAtoms(nat)%ptr_atom
			CALL SEE_atom(id, icell, atom, dM)
   !compute dTex and new values of Tex
			CALL calc_delta_Tex_atom(icell, atom, dT, Tex, Tion)
   

			!For one atoms = for all transitions
			if (verbose) then
				write(*,*) atom%ID, " -> dM = ", real(dM), " -> dT = ", real(dT)
				write(*,*) "   <::> Tex (K) = ", Tex, " Tion (K) = ", Tion, ' Te (K) = ', atmos%T(icell)
			endif
   
   !max for all atoms
			dpop = max(dpop, dM)
   ! and all transitions...
			dTex = max(dTex, dT)
			   
			atom => NULL()
		enddo
  
		!flag the one with a "*"
		delta = dTex
		!delta = dM

		!Compare all atoms
		if (verbose) then
			write(*,*) icell, "    >> *max(dT) = ", real(dTex)
			write(*,*) icell, "    >> max(dpops) = ", real(dpop)
		endif
 
	RETURN
	END SUBROUTINE update_populations
 
	SUBROUTINE Gamma_LTE(id,icell)
	! ------------------------------------------------------------------------- !
	! Fill the rate matrix Gamma, whose elements are Gamma(lp,l) is the rate
	! of transition from level lp to l.
	! At initialisation, Gamma(lp,l) = C(J,I), the collisional rates from upper
	! level j to lower level i.
	!
	! This is the LTE case where Gamma(l',l) = C(l'l)
	! Gamma(l',l) = Cl'l - delta(l,l')Sum_l" (Cll").
	!
	! This Gamma is frequency and angle independent. 
	! FOR ALL ATOMS
	! ------------------------------------------------------------------------- !
		integer, intent(in) :: id, icell
		integer :: nact, kr, l, lp, nati, natf
		type (AtomType), pointer :: atom

	!nati = (1. * (id-1)) / NLTEspec%NPROC * atmos%Nactiveatoms + 1
	!natf = (1. * id) / NLTEspec%NPROC * atmos%Nactiveatoms
		nati = 1; natf = atmos%Nactiveatoms
		do nact=nati,natf !loop over each active atoms
   			atom => atmos%ActiveAtoms(nact)%ptr_atom
!    do l=1,atom%Nlevel
!     do lp=1,atom%Nlevel   
!       atom%Gamma(lp, l) = atom%Ckij(icell, (l-1)*atom%Nlevel+lp) !lp->l; C_kij = C(j->i)
!     end do
!   end do

   !fill the diagonal here, delta(l',l)
   !because for each G(i,i); Cii, Rii is 0.
   !and Gii = -Sum_j Cij + Rij = -sum_j Gamma(i,j)
   !diagonal of Gamma, Gamma(col,col) is sum_row Gamma(col,row)
   !G(1,1) = - (G12 + G13 + G14 ...)
   !G(2,2) = - (G21 + G23 + G24 ..) first index is for column and second row
   
   !wavelength and ray indep can remove it now
			do l = 1, atom%Nlevel
				atom%Gamma(l,l,id) = 0d0
				atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
			enddo

   
			NULLIFY(atom)
  		end do !loop over atoms

	RETURN
	END SUBROUTINE Gamma_LTE


END MODULE statequil_atoms

