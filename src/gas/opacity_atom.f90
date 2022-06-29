!No dissolution and no mag yet!
module Opacity_atom

   use atom_type
   use grid
   use parametres
   use broad, Only               : line_damping
   use voigts, only              : voigt
   use gas_contopac, only        : H_bf_Xsection, alloc_gas_contopac, background_continua_lambda, &
                                     dealloc_gas_contopac, hnu_k
   use wavelengths, only         :  n_lambda
   use wavelengths_gas, only     : Nlambda_max_line, Nlambda_max_cont, n_lambda_cont, tab_lambda_cont, tab_lambda_nm, Nlambda_max_line_vel
   use constantes, only          : c_light
   use molecular_emission, only  : v_proj, ds
   use utils, only               : linear_1D_sorted
   !$ use omp_lib

   implicit none

   !local profile for cell id in direction iray for all atoms and b-b trans
   real(kind=dp), allocatable :: Itot(:,:,:), psi(:,:,:), phi_loc(:,:,:,:,:)
   real(kind=dp), allocatable :: eta_atoms(:,:,:), Uji_down(:,:,:,:), chi_up(:,:,:,:), chi_down(:,:,:,:), chi_tot(:)
   integer, parameter 		   :: NvspaceMax = 151
   logical 		               :: lnon_lte_loop                                      

   contains

   !could be parralel
   subroutine alloc_atom_opac(N,x)
      integer, intent(in) :: N
      real(kind=dp), dimension(N) :: x
      integer :: nat, kr, icell, nb, nr
      type(AtomType), pointer :: atm
      real(kind=dp) :: vth
      integer(kind=8) :: mem_loc

      mem_loc = 0

      ! call alloc_gas_contopac(N,x)
      !-> on a small grid and interpolated later
      call alloc_gas_contopac(n_lambda_cont,tab_lambda_cont)

      do nat=1, n_atoms
         atm => atoms(nat)%p
         !if commong gauss prof add in mem_loc
         do kr=1,atm%nline
               if (.not.atm%lines(kr)%lcontrib) cycle
               if (atm%lines(kr)%Voigt) then
                  allocate(atm%lines(kr)%v(atm%lines(kr)%Nlambda),atm%lines(kr)%phi(atm%lines(kr)%Nlambda,n_cells))
                  allocate(atm%lines(kr)%a(n_cells))
                  mem_loc = mem_loc + sizeof(atm%lines(kr)%a)+sizeof(atm%lines(kr)%phi)+sizeof(atm%lines(kr)%v)
               endif
         enddo

         do icell=1, n_cells
            if (icompute_atomRT(icell) <= 0) cycle
            do kr=1,atm%nline
               if (.not.atm%lines(kr)%lcontrib) cycle
               if (atm%lines(kr)%Voigt) then
                  nb = atm%lines(kr)%nb; nr = atm%lines(kr)%nr
                  atm%lines(kr)%a(icell) = line_damping(icell,atm%lines(kr))
                  !tmp because of vbroad, recomputed after in m/s
                  vth = vbroad(T(icell),atm%weight, vturb(icell))
                  atm%lines(kr)%v(:) = c_light * (x(nb:nr)-atm%lines(kr)%lambda0)/atm%lines(kr)%lambda0 / vth !unitless
                  atm%lines(kr)%phi(:,icell) = Voigt(atm%lines(kr)%Nlambda, atm%lines(kr)%a(icell), atm%lines(kr)%v(:)) / (vth * sqrtpi)
               endif
            enddo
         enddo

         do kr=1,atm%nline
            if (.not.atm%lines(kr)%lcontrib) cycle
            if (atm%lines(kr)%Voigt) then
               nb = atm%lines(kr)%nb; nr = atm%lines(kr)%nr
               !-> will not change during the non-LTE loop.
               atm%lines(kr)%v(:) = c_light * (x(nb:nr)-atm%lines(kr)%lambda0)/atm%lines(kr)%lambda0 !m/s
            endif
         enddo

         do kr = 1, atm%Ncont
            if (.not.atm%continua(kr)%lcontrib) cycle
            ! allocate(atm%continua(kr)%twohnu3_c2(atm%continua(kr)%Nlambda))
            ! atm%continua(kr)%twohnu3_c2(:) = twohc/x(atm%continua(kr)%Nb:atm%continua(kr)%Nr)**3
            ! allocate(atm%continua(kr)%alpha(atm%continua(kr)%Nlambda))
            ! if (atm%continua(kr)%hydrogenic) then
            !    atm%continua(kr)%alpha(:) = H_bf_Xsection(atm%continua(kr), x(atm%continua(kr)%Nb:atm%continua(kr)%Nr))
            ! else
            !    atm%continua(kr)%alpha(:) = linear_1D_sorted(size(atm%continua(kr)%alpha_file),&
            !    atm%continua(kr)%lambda_file,atm%continua(kr)%alpha_file,atm%continua(kr)%Nlambda,x(atm%continua(kr)%Nb:atm%continua(kr)%Nr))
            ! endif
            !-> on a small grid and interpolated later
            allocate(atm%continua(kr)%twohnu3_c2(atm%continua(kr)%Nlambdac))
            atm%continua(kr)%twohnu3_c2(:) = twohc/tab_lambda_cont(atm%continua(kr)%Nbc:atm%continua(kr)%Nrc)**3
            allocate(atm%continua(kr)%alpha(atm%continua(kr)%Nlambdac))
            if (atm%continua(kr)%hydrogenic) then
               atm%continua(kr)%alpha(:) = H_bf_Xsection(atm%continua(kr), tab_lambda_cont(atm%continua(kr)%Nbc:atm%continua(kr)%Nrc))
            else
               atm%continua(kr)%alpha(:) = linear_1D_sorted(size(atm%continua(kr)%alpha_file),&
               atm%continua(kr)%lambda_file,atm%continua(kr)%alpha_file,atm%continua(kr)%Nlambdac,tab_lambda_cont(atm%continua(kr)%Nbc:atm%continua(kr)%Nrc))
            endif

         enddo

         atm => null()
      enddo

      write(*,'("allocate "(1F6.4)" GB for line profiles")') real(mem_loc) / 1024.0/1024./1024.0


      return
   end subroutine alloc_atom_opac

   subroutine dealloc_atom_opac()
      integer :: nat, kr
      type(AtomType), pointer :: atm

      call dealloc_gas_contopac()
  
      do nat=1, n_atoms
         atm => atoms(nat)%p

         !first allocation in io_atom
         if (allocated(atm%vg)) deallocate(atm%vg,atm%phig)

         do kr=1,atm%nline
            if (allocated( atm%lines(kr)%a)) deallocate(atm%lines(kr)%a )
            if (allocated( atm%lines(kr)%v)) deallocate(atm%lines(kr)%v,atm%lines(kr)%phi)
         enddo

         do kr = 1, atm%Ncont
            if (allocated(atm%continua(kr)%twohnu3_c2)) then 
               deallocate(atm%continua(kr)%twohnu3_c2)
               deallocate(atm%continua(kr)%alpha)
            endif
         enddo

         atm => null()
      enddo


      return
   end subroutine dealloc_atom_opac

   subroutine contopac_atom_loc(icell,N,lambda,chi,snu)
      integer, intent(in) :: icell, N
      real(kind=dp), intent(in), dimension(N) :: lambda
      real(kind=dp), intent(inout), dimension(N) :: chi, Snu
      real(kind=dp), dimension(N_lambda_cont) :: chic, snuc
      integer :: la, lac, i0
      real(kind=dp) :: w

      !init continuous opacity with background gas continuum.
      call background_continua_lambda(icell, n_lambda_cont, tab_lambda_cont, chic, Snuc)
      !Snu = Snu + scat(lambda, icell) * Jnu(:,icell)
      !accumulate b-f
      call opacity_atom_bf_loc(icell, n_lambda_cont, tab_lambda_cont, chic, Snuc)

      chi(1) = chic(1)
      snu(1) = snuc(1)

      !linear interpolation
      i0 = 2
      do la=1, N
         loop_i : do lac=i0, n_lambda_cont
            if (tab_lambda_cont(lac) > lambda(la)) then
               w = (lambda(la) - tab_lambda_cont(lac-1)) / (tab_lambda_cont(lac) - tab_lambda_cont(lac-1))
               chi(la) = (1.0_dp - w) * chic(lac-1)  + w * chic(lac)
               snu(la) = (1.0_dp - w) * snuc(lac-1)  + w * snuc(lac)
               i0 = lac
               exit loop_i
            endif
         enddo loop_i 
      enddo
      chi(N) = chic(n_lambda_cont)
      snu(N) = snuc(n_lambda_cont)

      return
   end subroutine contopac_atom_loc

   subroutine opacity_atom_bf_loc(icell,N,lambda,chi,Snu)
   !to do: remove lambda dep since it must be consistent with Nr, Nb
      integer, intent(in) :: icell, N
      real(kind=dp), intent(in), dimension(N) :: lambda
      real(kind=dp), intent(inout), dimension(N) :: chi, Snu
      integer :: nat, Nred, Nblue, kr, i, j, idelem
      type(AtomType), pointer :: atom
      real(kind=dp) :: wi, wj, chi_ion, Diss, nn, gij, ni_on_nj_star
      real(kind=dp), dimension(N) :: ehnukt

      ehnukt(:) = exp(hnu_k/T(icell))!exphckT(icell)**(lambda_base/lambda(:))

      atom_loop : do nat = 1, N_atoms
         atom => Atoms(nat)%p
         idelem = atom%periodic_table

         tr_loop : do kr = 1,atom%Ncont

            if (.not.atom%continua(kr)%lcontrib) cycle

            !beware Nc here, assumed tab_lambda_cont
            Nred = atom%continua(kr)%Nrc; Nblue = atom%continua(kr)%Nbc
            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            !ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
            ni_on_nj_star = atom%nstar(i,icell)/(atom%nstar(j,icell) + 1d-100)

            gij = ni_on_nj_star * exp(-hc_k/T(icell)/atom%continua(kr)%lambda0)
            if ((atom%n(i,icell) - atom%n(j,icell) * gij) <= 0.0_dp) then
               cycle tr_loop
            endif

            !should be the same as directly exp(-hc_kT/lambda)
            chi(Nblue:Nred) = chi(Nblue:Nred) + atom%continua(kr)%alpha(:) * (atom%n(i,icell) - &
               ni_on_nj_star * ehnukt(Nblue:Nred) * atom%n(j,icell))

            Snu(Nblue:Nred) = Snu(Nblue:Nred) + atom%n(j,icell) * atom%continua(kr)%alpha(:) * atom%continua(kr)%twohnu3_c2(:) *& 
               ni_on_nj_star * ehnukt(Nblue:Nred)

         end do tr_loop

         atom => NULL()
      end do atom_loop


      return
   end subroutine opacity_atom_bf_loc

   function calc_vloc(icell,u,v,w,x,y,z,x1,y1,z1)
      !computes the local mean velocity in direction (u,v,w) of the cell icell
      !by averaging between points (x,y,z) and (x1,y1,z1).
      real(kind=dp) :: calc_vloc
      integer, intent(in) :: icell
      real(kind=dp), intent(in) :: x, y, z, x1, y1, z1
      real(kind=dp), intent(in) :: u, v, w

      calc_vloc = 0.5 * (v_proj(icell,x1,y1,z1,u,v,w) + v_proj(icell,x,y,z,u,v,w))
      return
   end function calc_vloc

   subroutine opacity_atom_bb_loc(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l_void_before,l_contrib, &
         iterate,N,lambda,chi,Snu)
   !to do: remove lambda dep since it must be consistent with Nr, Nb
      integer, intent(in) :: id, icell, iray,N
      logical, intent(in) :: iterate
      real(kind=dp), intent(in), dimension(N) :: lambda
      real(kind=dp), intent(inout), dimension(N) :: chi, Snu
      real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l_void_before,l_contrib
      integer :: nat, Nred, Nblue, kr, i, j, Nlam
      real(kind=dp) :: dv
      type(AtomType), pointer :: atom
      real(kind=dp), dimension(Nlambda_max_line) :: phi0
      ! real(kind=dp), dimension(Nlambda_max_line_vel) :: phi0

      dv = 0.0_dp
      if (lnon_lte_loop.and..not.iterate) then !not iterate but non-LTE
         dv = calc_vloc(icell,u,v,w,x,y,z,x1,y1,z1)
      endif

      atom_loop : do nat = 1, N_Atoms
         atom => Atoms(nat)%p

         tr_loop : do kr = 1,atom%Nline

            if (.not.atom%lines(kr)%lcontrib) cycle

            !if .not.labs (image or LTE), Nred and Nblue includes the maxium extension of line
            !due to velocity fields (and not only the natural width + damping)
            Nred = atom%lines(kr)%Nr; Nblue = atom%lines(kr)%Nb
            Nlam = atom%lines(kr)%Nlambda
            i = atom%lines(kr)%i; j = atom%lines(kr)%j

            if ((atom%n(i,icell) - atom%n(j,icell)*atom%lines(kr)%gij) <= 0.0_dp) cycle tr_loop
            
            if (abs(dv)>atom%lines(kr)%vmax) then
               !move the profile to the red edge up to Nover_sup
               !change Nlam ??
               if (dv > 0) then
                  ! Nblue = Nred
                  Nred = atom%lines(kr)%Nover_sup
                  Nblue =  Nred - Nlam + 1
               !move to the blue edge down to Nover_inf
               else
                  ! Nred = Nblue
                  Nblue =  atom%lines(kr)%Nover_inf
                  Nred = Nlam + Nblue - 1
               endif
               Nlam = Nred - Nblue + 1
            endif

            phi0(1:Nlam) = profile_art_i(atom%lines(kr),id,icell,iterate,Nlam,lambda(Nblue:Nred),&
                                 x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
            !to interpolate the profile we need to find the index of the first lambda on the grid and then increment


            chi(Nblue:Nred) = chi(Nblue:Nred) + &
               hc_fourPI * atom%lines(kr)%Bij * phi0(1:Nlam) * (atom%n(i,icell) - atom%lines(kr)%gij*atom%n(j,icell))

            Snu(Nblue:Nred) = Snu(Nblue:Nred) + &
               hc_fourPI * atom%lines(kr)%Aji * phi0(1:Nlam) * atom%n(j,icell)


            if ((iterate.and.atom%active)) then
               phi_loc(1:Nlam,atom%ij_to_trans(i,j),atom%activeindex,iray,id) = phi0(1:Nlam)
            endif


            !if (lmagnetized) then
            !endif


         end do tr_loop

         atom => null()

      end do atom_loop


      return
   end subroutine opacity_atom_bb_loc


   subroutine xcoupling(id, icell, iray, lupdate_psi)
   !test: update_psi in case of local subiterations.
   !background not iterated (included in I)
   !to do add bound-free
   !to do add background
      integer, intent(in) :: id, icell, iray
      logical, intent(in) :: lupdate_psi
      integer :: nact, j, i, kr, Nb, Nr, la, Nl
      integer :: i0, j0, la0, N1, N2
      type (AtomType), pointer :: aatom
      real(kind=dp) :: gij, chicc, wl, ni_on_nj_star, wphi, twohnu3_c2, wt
      real(kind=dp) :: term1(Nlambda_max_cont), term2(Nlambda_max_cont), term3(Nlambda_max_cont)
      real(kind=dp), dimension(Nlambda_max_line) :: phi0, wei_line

      Uji_down(:,:,:,id) = 0.0_dp
      chi_down(:,:,:,id) = 0.0_dp
      chi_up(:,:,:,id)   = 0.0_dp

      eta_atoms(:,:,id) = 0.0_dp
      !-> missing LTE here
      if (lupdate_psi) chi_tot(:) = 1d-50

      aatom_loop : do nact=1, Nactiveatoms

         aatom => ActiveAtoms(nact)%p

         cont_loop : do kr = 1, aatom%Ncont

            if (.not.aatom%continua(kr)%lcontrib) cycle cont_loop

            j = aatom%continua(kr)%j
            i = aatom%continua(kr)%i
            Nb = aatom%continua(kr)%Nbc; Nr = aatom%continua(kr)%Nrc
            Nl = Nr-Nb+1
            N1 = aatom%continua(kr)%Nb; N2 = aatom%continua(kr)%Nr

            !ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
            ni_on_nj_star = aatom%nstar(i,icell)/(aatom%nstar(j,icell) + 1d-100)


            gij = ni_on_nj_star * exp(-hc_k/T(icell)/aatom%continua(kr)%lambda0)

      
            if (aatom%n(i,icell) - gij*aatom%n(j,icell) <= 0.0_dp) cycle cont_loop

            !The wavelength integration weight cannot be used in that loop because we interpolate after.
            freq_loop : do la=1, Nl
               gij = ni_on_nj_star * exp(hnu_k(Nb-1+la)/T(icell))
               term1(la) = fourpi_h * aatom%continua(kr)%alpha(la) * (aatom%n(i,icell) - gij*aatom%n(j,icell))
               term2(la) = aatom%continua(kr)%alpha(la) * aatom%continua(kr)%twohnu3_c2(la) * gij
               term3(la) = aatom%continua(kr)%alpha(la) * aatom%continua(kr)%twohnu3_c2(la) * gij * aatom%n(j,icell)
            enddo freq_loop
            !linear interpolation + adding wavelength integration weight
            !wl = 0.5*(tab_lambda_nm(N1+1)-tab_lambda_nm(N1)) / tab_lambda_nm(N1)
            chi_down(N1,j,nact,id) = chi_down(N1,j,nact,id) + term1(1)! * wl
            chi_up(n1,i,nact,id) = chi_up(n1,i,nact,id) + term1(1)! * wl
            Uji_down(n1,j,nact,id) = Uji_down(N1,j,nact,id) + term2(1)
            eta_atoms(N1,nact,id) = eta_atoms(N1,nact,id) + term3(1)
            i0 = 2
            do la=N1,N2
               ! if (la>1) then
               !    wl = 0.5*(tab_lambda_nm(la+1)-tab_lambda_nm(la-1)) / tab_lambda_nm(la)
               ! !otherwise, wl is the previous one, first point of the grid
               ! endif
               loop_i : do la0=i0, Nl
                  if (tab_lambda_cont(Nb+la0-1) > tab_lambda_nm(la)) then
                     wt = (tab_lambda_nm(la) - tab_lambda_cont(Nb+la0-2)) / (tab_lambda_cont(Nb+la0-1) - tab_lambda_cont(Nb+la0-2))

                     chi_down(la,j,nact,id) = chi_down(la,j,nact,id) + ( (1.0_dp - wt) *term1(la0-1)  + wt * term1(la0) )!*wl
                     chi_up(la,i,nact,id) = chi_up(la,i,nact,id) + ( (1.0_dp - wt) *term1(la0-1)  + wt * term1(la0) )!*wl
                     Uji_down(la,j,nact,id) = Uji_down(la,j,nact,id) + (1.0_dp - wt) *term2(la0-1)  + wt * term2(la0)
                     eta_atoms(la,nact,id) = eta_atoms(la,nact,id) + (1.0_dp - wt) *term3(la0-1)  + wt * term3(la0)

                     i0 = la0
                     exit loop_i
                  endif
               enddo loop_i
            enddo
            ! wl = 0.5*(tab_lambda_nm(N2)-tab_lambda_nm(N2-1)) / tab_lambda_nm(N2-1)
            chi_down(n2,j,nact,id) = chi_down(n2,j,nact,id) + term1(Nl)! * wl
            chi_up(n2,i,nact,id) = chi_up(n2,i,nact,id) + term1(Nl)! * wl
            Uji_down(n2,j,nact,id) = Uji_down(n2,j,nact,id) + term2(Nl)
            eta_atoms(n2,nact,id) = eta_atoms(n2,nact,id) + term3(Nl)

         enddo cont_loop

         do la=1, n_lambda
             if (la==1) then
                wl = 0.5*(tab_lambda_nm(la+1)-tab_lambda_nm(la)) / tab_lambda_nm(la)
             elseif (la==n_lambda) then
                wl = 0.5*(tab_lambda_nm(la)-tab_lambda_nm(la-1)) / tab_lambda_nm(la-1)
             else
                wl = 0.5*(tab_lambda_nm(la+1)-tab_lambda_nm(la-1)) / tab_lambda_nm(la)
             endif    
             chi_down(la,:,:,id) =  chi_down(la,:,:,id)  * wl
             chi_up(la,:,:,id) = chi_up(la,:,:,id) * wl     
         enddo


         line_loop : do kr = 1, aatom%Nline

            if (.not.aatom%lines(kr)%lcontrib) cycle line_loop

            j = aatom%lines(kr)%j
            i = aatom%lines(kr)%i

            if (aatom%n(i,icell) - aatom%lines(kr)%gij*aatom%n(j,icell) <= 0.0_dp) cycle line_loop

            Nb = aatom%lines(kr)%Nb; Nr = aatom%lines(kr)%Nr
            Nl = Nr - Nb + 1

            wphi = 0.0
            do la=1,Nl
               if (la==1) then
                  wl = 0.5*(tab_lambda_nm(Nb+1)-tab_lambda_nm(Nb)) * c_light / aatom%lines(kr)%lambda0
               elseif (la==Nl) then
                  wl = 0.5*(tab_lambda_nm(Nr)-tab_lambda_nm(Nr-1)) * c_light / aatom%lines(kr)%lambda0
               else
                  wl = 0.5*(tab_lambda_nm(Nb+la)-tab_lambda_nm(Nb+la-2)) * c_light / aatom%lines(kr)%lambda0
               endif
               wei_line(la) = wl
               phi0(la) = phi_loc(la,kr,nact,iray,id)
               wphi = wphi + wl * phi0(la)
            enddo

            freq2_loop : do la=1, Nl
               wl = wei_line(la)

               Uji_down(Nb-1+la,j,nact,id) = Uji_down(Nb-1+la,j,nact,id) + hc_fourPI * aatom%lines(kr)%Aji * phi0(la)!/wphi

               !there is a hc_fourpi factor that simplifies here, because integral is over dnu/hnu dOmega = dv/hc dOmega * hc/4pi
               !dv dOmega/4pi which is whtat is contained in wl (for dv) and what the angular integration provides (dOmega/4pi)
               chicc = wl * aatom%lines(kr)%Bij * (aatom%n(i,icell) - aatom%lines(kr)%gij*aatom%n(j,icell)) * phi0(la)/wphi

               chi_down(Nb-1+la,j,nact,id) = chi_down(Nb-1+la,j,nact,id) + chicc
               chi_up(Nb-1+la,i,nact,id) = chi_up(Nb-1+la,i,nact,id) + chicc

               eta_atoms(Nb-1+la,nact,id) = eta_atoms(Nb-1+la,nact,id) + &
                  hc_fourPI * aatom%lines(kr)%Aji * aatom%n(j,icell) * phi0(la)!/wphi


            enddo freq2_loop
            if (lupdate_psi) then
               chi_tot(Nb:Nr) = chi_tot(Nb:Nr) + hc_fourPI * aatom%lines(kr)%Bij * (aatom%n(i,icell) - aatom%lines(kr)%gij*aatom%n(j,icell)) * phi0(1:Nl)
            endif

         enddo line_loop

         aatom => null()
      enddo aatom_loop

      if (lupdate_psi) then
         psi(:,1,id) = (1.0_dp - exp(-ds(iray,id)*chi_tot))/chi_tot
      endif

    return
   end subroutine xcoupling

   function profile_art(line,id,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
      ! phi = Voigt / sqrt(pi) / vbroad(icell)
      integer, intent(in)                    :: id, icell, N
      type (AtomicLine), intent(in)          :: line
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer 											:: Nvspace
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: norm, vth
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                                dv, omegav_mean
      integer                                :: Nred, Nblue, i, j, nv
      real(kind=dp), dimension(N)            :: u0, profile_art, u1, u0sq

      Nvspace = NvspaceMax
      i = line%i; j = line%j
      Nred = line%Nr; Nblue = line%Nb
      vth = vbroad(T(icell),line%Atom%weight, vturb(icell))

      u0(:) = (lambda - line%lambda0)/line%lambda0  * ( c_light/vth )

      v0 = v_proj(icell,x,y,z,u,v,w)
      if (lvoronoi) then
         omegav(1) = v0
         Nvspace = 1
         omegav_mean = v0
      else

         Omegav = 0.0
         omegav(1) = v0
         v1 = v_proj(icell,x1,y1,z1,u,v,w)

         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/vth*20.)),NvspaceMax)

         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
         enddo
         omegav(Nvspace) = v1
         omegav_mean = sum(omegav(1:Nvspace))/real(Nvspace,kind=dp)
      endif
      !in non-LTE:
      !the actual cell icell_nlte must be centered on 0 (moving at vmean).
      !the other cells icell crossed must be centered in v(icell) - vmean(icell_nlte) 
      if (lsubstract_avg) then!labs == .true.
         omegav(1:Nvspace) = omegav(1:Nvspace) - omegav_mean
      endif


      if (line%voigt) then
         u1(:) = u0(:) - omegav(1)/vth
         profile_art(:) = Voigt(N, line%a(icell), u1(:))
         do nv=2, Nvspace
            u1(:) = u0(:) - omegav(nv)/vth
            profile_art(:) = profile_art(:) + Voigt(N, line%a(icell), u1(:))
         enddo

      else
         !u1 = (u0 - omegav(nv)/vth)**2
         u0sq(:) = u0(:)*u0(:)
         u1(:) = u0sq(:) + (omegav(1)/vth)*(omegav(1)/vth) - 2*u0(:) * omegav(1)/vth
         profile_art(:) = exp(-u1(:))
         do nv=2, Nvspace
            u1(:) = u0sq(:) + (omegav(nv)/vth)*(omegav(nv)/vth) - 2*u0(:) * omegav(nv)/vth
            profile_art(:) = profile_art(:) + exp(-u1(:))
         enddo
      endif

      norm = Nvspace * vth * sqrtpi
      profile_art(:) = profile_art(:) / norm

      return
   end function profile_art

   function profile_art_i(line,id,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
      integer, intent(in)                    :: icell, N, id
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer 											:: Nvspace
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: norm, vth
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                                dv, omegav_mean
      type (AtomicLine), intent(in)          :: line
      integer                                :: Nred, Nblue, i, j, nv
      real(kind=dp), dimension(N)            :: uloc, u0, profile_art_i, u1, u0sq

      Nvspace = NvspaceMax
      i = line%i; j = line%j
      Nred = line%Nr; Nblue = line%Nb
      vth = vbroad(T(icell),line%Atom%weight, vturb(icell))

      u0(:) = (lambda - line%lambda0)/line%lambda0  * ( c_light/vth )

      v0 = v_proj(icell,x,y,z,u,v,w)
      if (lvoronoi) then
         omegav(1) = v0
         Nvspace = 1
         omegav_mean = v0
      else

         Omegav = 0.0
         omegav(1) = v0
         v1 = v_proj(icell,x1,y1,z1,u,v,w)

         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/vth*20.)),NvspaceMax)

         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
         enddo
         omegav(Nvspace) = v1
         omegav_mean = sum(omegav(1:Nvspace))/real(Nvspace,kind=dp)
      endif
      if (lsubstract_avg) then!labs == .true.
         omegav(1:Nvspace) = omegav(1:Nvspace) - omegav_mean
      endif


      if (line%voigt) then
         u1(:) = u0(:) - omegav(1)/vth
         uloc(:) = line%v(:) / vth
         profile_art_i(:) = linear_1D_sorted(N,uloc(:),line%phi(:,icell),N,u1(:))
         do nv=2, Nvspace
            u1(:) = u0(:) - omegav(nv)/vth
            profile_art_i(:) = profile_art_i(:) + linear_1D_sorted(N,uloc(:),line%phi(:,icell),N,u1(:))
         enddo

      else
         !u1 = (u0 - omegav(nv)/vth)**2
         u0sq(:) = u0(:)*u0(:)
         u1(:) = u0sq(:) + (omegav(1)/vth)*(omegav(1)/vth) - 2*u0(:) * omegav(1)/vth
         profile_art_i(:) = exp(-u1(:)) / sqrtpi / vth
         do nv=2, Nvspace
            u1(:) = u0sq(:) + (omegav(nv)/vth)*(omegav(nv)/vth) - 2*u0(:) * omegav(nv)/vth
            profile_art_i(:) = profile_art_i(:) + exp(-u1(:)) / sqrtpi / vth
         enddo
      endif

      profile_art_i(:) = profile_art_i(:) / Nvspace

      return
   end function profile_art_i


   subroutine local_intensity(id,icell,iray,iloc,ii)
   !computes the contribution of a cell to the emergent intensity at that cell (ii)
   ! and the total intensity to that cell (iloc)
      integer, intent(in) :: id, icell, iray
      real(kind=dp) :: x, y, z, x1, y1, z1, u, v, w, l_void_before,l_contrib!intent(in) ? 
      real(kind=8), intent(out) :: ii(:), iloc(:)
      real(kind=dp) :: chi0(n_lambda)
      integer :: nat, i, j, kr, Nr, Nb, Nlam
      real(kind=dp), dimension(Nlambda_max_line) :: phi0
      type(AtomType), pointer :: atom

      iloc(:) = Itot(:,1,id)
      call contopac_atom_loc(icell,n_lambda,tab_lambda_nm,chi0,ii(:))
      !add lines in the ref of the cell

      x = 0.0_dp; y = x; z = x; x1 = 0.0_dp; y1 = x1; z1 = x1
      u = 0.0_dp; v = 0.0_dp; w = 0.0_dp
      l_void_before = 0.0_dp; l_contrib = 0.0_dp

      atom_loop : do nat = 1, N_Atoms
         atom => Atoms(nat)%p

         tr_loop : do kr = 1,atom%Nline

            if (.not.atom%lines(kr)%lcontrib) cycle

            Nr = atom%lines(kr)%Nr; Nb = atom%lines(kr)%Nb
            Nlam = atom%lines(kr)%Nlambda
            i = atom%lines(kr)%i; j = atom%lines(kr)%j

            if ((atom%n(i,icell) - atom%n(j,icell)*atom%lines(kr)%gij) <= 0.0_dp) cycle tr_loop
            

            phi0(1:Nlam) = profile_art_i(atom%lines(kr),id,icell,.true.,Nlam,tab_lambda_nm(Nb:Nr),&
                                 x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)


            chi0(Nb:Nr) = chi0(Nb:Nr) + &
               hc_fourPI * atom%lines(kr)%Bij * phi0(1:Nlam) * (atom%n(i,icell) - atom%lines(kr)%gij*atom%n(j,icell))

            ii(Nb:Nr) = ii(Nb:Nr) + hc_fourPI * atom%lines(kr)%Aji * phi0(1:Nlam) * atom%n(j,icell)

         end do tr_loop

         atom => null()

      end do atom_loop

      !or ii(:) = ii(:) * psi(:,1,id) !(1.0 - exp(-dtau))/chi
      ii(:) = ii(:)/chi0(:) * (1.0_dp - exp(-chi0*ds(iray,id)))

      return
   end subroutine local_intensity
 
   subroutine write_opacity_emissivity_bin(Nlambda,lambda)
      !not para to be able to write while computing!
      integer, intent(in) :: Nlambda
      real(kind=dp), intent(in) :: lambda(Nlambda)
      integer :: unit, unit2, status = 0
      integer :: alloc_status, id, icell, m, Nrec
      real(kind=dp), allocatable, dimension(:,:,:) :: chi_tmp, eta_tmp, rho_tmp
      character(len=11) :: filename_chi="chi_map.bin"
      character(len=50) :: filename_eta="eta_map.bin"  
      character(len=18) :: filename_rho="magnetoopt_map.bin"  
      
      call warning("opacity emissivity map bug here because of molecular emission and r=0!")

      write(*,*) " Writing emissivity and opacity map (rest frame) ..."
      ! if (lmagnetized) then
      !    Nrec = 4
      !    allocate(rho_tmp(Nlambda,n_cells,Nrec-1),stat=alloc_status)
      !    if (alloc_status /= 0) call error("Cannot allocate rho_tmp !")
      !    call warning("  Polarized case with an observer // to z!")
      ! else
         Nrec = 1
      ! endif
      allocate(chi_tmp(Nlambda, n_cells, Nrec),stat=alloc_status)
      if (alloc_status /= 0) call error("Cannot allocate chi_tmp !")
      allocate(eta_tmp(Nlambda, n_cells, Nrec),stat=alloc_status)
      if (alloc_status /= 0) call error("Cannot allocate eta_tmp !")
  
  
      call ftgiou(unit,status)
      call ftgiou(unit2,status)
      open(unit, file=trim(filename_chi),form="unformatted",status='unknown',access="sequential",iostat=status)
      open(unit2, file=trim(filename_eta),form="unformatted",status='unknown',access="sequential",iostat=status) 
      id = 1
      do icell=1, n_cells
         !$ id = omp_get_thread_num() + 1
         if (icompute_atomRT(icell) > 0) then
            call contopac_atom_loc(icell,Nlambda,lambda,chi_tmp(:,icell,1),eta_tmp(:,icell,1))
            call opacity_atom_bb_loc(id,icell,1,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,&
               0d0,0d0,.false.,Nlambda,lambda,chi_tmp(:,icell,1), eta_tmp(:,icell,1))
            ! do m=2,Nrec
            !    !m=1, unpolarized already filled.
            !    !etaQUV, chiQUV only for Q(=1), U, V
            !    eta_tmp(:,m,icell) = etaQUV_p(:,m-1,id)
            !    chi_tmp(:,m,icell) = chiQUV_p(:,m-1,id)
            !    !only QUV for rho_p
            !    rho_tmp(:,m-1,icell) = rho_p(:,m-1,id)
            ! enddo
         else
            chi_tmp(:,icell,:) = 0.0
            eta_tmp(:,icell,:) = 0.0
         endif
         !for each cell, write all wavelengths and all records (only one if not pol, or 4 Stokes)
         write(unit,iostat=status) chi_tmp(:,icell,:)
         write(unit2,iostat=status) eta_tmp(:,icell,:)
         !if lmagnetized write rho_tmp(:,icell,Nrec)
      enddo
  
   
      ! if (lmagnetized) then
      !   write(unit, iostat=status) rho_tmp
      !   deallocate(rho_tmp)
      ! endif
      deallocate(chi_tmp, eta_tmp)
      close(unit)
      close(unit2)                 
  
      !free unit
      call ftfiou(unit, status)
      call ftfiou(unit2, status)
    
      return
   end subroutine write_opacity_emissivity_bin

end module Opacity_atom
