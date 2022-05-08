!No dissolution and no mag yet!
module Opacity_atom

   use atom_type
   use grid
   use parametres
   use broad, Only               : line_damping
   use voigts, only              : voigt
   use gas_contopac, only        : H_bf_Xsection, alloc_gas_contopac, background_continua_lambda, &
                                     dealloc_gas_contopac, hnu_k!exphckT, lambda_base,
   use wavelengths_gas, only     : Nlambda_max_line, n_lambda_cont, tab_lambda_cont
   use constantes, only          : c_light
   use molecular_emission, only  : v_proj
   use utils, only               : linear_1D_sorted
   !$ use omp_lib

   implicit none

   !local profile for cell id in direction iray for all atoms and b-b trans
   real(kind=dp), allocatable :: Itot(:,:,:), psi(:,:,:), phi_loc(:,:,:,:,:)!move in see
   real(kind=dp) :: vlabs

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
  
      allocate(psi(N,1,nb_proc))
      allocate(Itot(N,1,nb_proc))

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
                  !tmp because of vbroad!
                  vth = vbroad(T(icell),atm%weight, vturb(icell))
                  atm%lines(kr)%v(:) = c_light * (x(nb:nr)-atm%lines(kr)%lambda0)/atm%lines(kr)%lambda0 / vth
                  atm%lines(kr)%phi(:,icell) = Voigt(atm%lines(kr)%Nlambda, atm%lines(kr)%a(icell), atm%lines(kr)%v(:)) / (vth * sqrtpi)
               endif
            enddo
         enddo

         do kr=1,atm%nline
            if (.not.atm%lines(kr)%lcontrib) cycle
            if (atm%lines(kr)%Voigt) then
               nb = atm%lines(kr)%nb; nr = atm%lines(kr)%nr
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
  
      if (allocated(psi)) then
         deallocate(psi)
      endif
      deallocate(Itot)

      do nat=1, n_atoms
         atm => atoms(nat)%p

         !first allocation in io_atom
         if (allocated(atm%vg)) deallocate(atm%vg,atm%phig)

         do kr=1,atm%nline
            if (allocated( atm%lines(kr)%a)) deallocate(atm%lines(kr)%a )
            if (allocated( atm%lines(kr)%v)) deallocate(atm%lines(kr)%v,atm%lines(kr)%phi)
         enddo

         do kr = 1, atm%Ncont
            deallocate(atm%continua(kr)%twohnu3_c2)
            deallocate(atm%continua(kr)%alpha)
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
      if (lambda(N)==tab_lambda_cont(n_lambda_cont)) then
         chi(N) = chic(n_lambda_cont)
         snu(N) = snuc(n_lambda_cont)
      endif


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
      type(AtomType), pointer :: atom
      real(kind=dp), dimension(Nlambda_max_line) :: phi0


      !average speed of the cell icell with iterate = (lsubstract_avg.and.nbr_cell==1) in direction iray.
      !always 0 if labs is .false.
      vlabs = 0.0_dp

      atom_loop : do nat = 1, N_Atoms
         atom => Atoms(nat)%p

         tr_loop : do kr = 1,atom%Nline

            if (.not.atom%lines(kr)%lcontrib) cycle

            !if .not.labs (image or LTE), Nred and Nblue includes the maxium extension of line
            !due to velocity fields (and not only the natural width + damping)
            Nred = atom%lines(kr)%Nr; Nblue = atom%lines(kr)%Nb
            Nlam = atom%lines(kr)%Nlambda
            i = atom%lines(kr)%i; j = atom%lines(kr)%j

            if ((atom%n(i,icell) - atom%n(j,icell)*atom%lines(kr)%gij) <= 0.0_dp) then
               cycle tr_loop
            endif
            
            if (iterate) then
               vlabs = calc_vloc(icell,u,v,w,x,y,z,x1,y1,z1)
               write(*,*) id, icell, " vlabs = ", vlabs * 1d-3
            else
            !-> no. in non-LTE (vlabs /= 0) if vlabs=0.0001 m/s the profile is completely shifted wich is wrong
            !compute the maximum relative shift (max(abs(v(icell)-v'(icell))))
               if (vlabs < 0) then
                  Nblue = atom%lines(kr)%Nover_inf
                  Nred = Nlam - 1 + Nblue
               elseif (vlabs > 0) then
                  Nred = atom%lines(kr)%Nover_sup
                  Nblue = Nred - Nlam + 1
               !vlabs == 0.0 (or labs) does not change Nblue and Nred
               endif
            endif

            phi0(1:Nlam) = profile_art(atom%lines(kr),icell,iterate,Nlam,lambda(Nblue:Nred),&
                                 x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
            !to interpolate the profile we need to find the index of the first lambda on the grid and then increment




            chi(Nblue:Nred) = chi(Nblue:Nred) + &
               hc_fourPI * atom%lines(kr)%Bij * phi0(1:Nlam) * (atom%n(i,icell) - atom%lines(kr)%gij*atom%n(j,icell))

            Snu(Nblue:Nred) = Snu(Nblue:Nred) + &
               hc_fourPI * atom%lines(kr)%Aji * phi0(1:Nlam) * atom%n(j,icell)


            !Store profile for integration or accumulate radiative rates here ?
            if (iterate) then
               phi_loc(1:Nlam,atom%ij_to_trans(i,j),nat,iray,id) = phi0(1:Nlam)
            endif


            !if (lmagnetized) then
            !endif


         end do tr_loop

         atom => null()

      end do atom_loop


      return
   end subroutine opacity_atom_bb_loc

   function profile_art(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
      ! phi = Voigt / sqrt(pi) / vbroad(icell)
      integer, intent(in)                    :: icell, N
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer, parameter 							:: NvspaceMax = 151                                      
      integer 											:: Nvspace
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: norm, vth
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                                dv, omegav_mean
      type (AtomicLine), intent(in)          :: line
      integer                                :: Nred, Nblue, i, j, nv
      real(kind=dp), dimension(N)            :: u0, profile_art, u1, u0sq
      ! real(kind=dp), dimension(N,NvspaceMax) :: u1


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
         !omegav_mean should be close to vlabs
         write(*,*) "vlabs=", vlabs*1d-3, omegav_mean * 1d-3
         stop
      else
         !recentre the cell on the speed of the non-lte cell,
         !so that the kinematics is computed with respect to a non-moving non-lte cell.
         !it is always 0 in LTE or image (when labs = .false.)
         omegav(1:Nvspace) = omegav(1:Nvspace) - vlabs
         ! write(*,*) "vlabs = ", vlabs * 1d-3
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

   function profile_art_i(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
      integer, intent(in)                    :: icell, N
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer, parameter 							:: NvspaceMax = 151                                      
      integer 											:: Nvspace
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: norm, vth
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                                dv, omegav_mean
      type (AtomicLine), intent(in)          :: line
      integer                                :: Nred, Nblue, i, j, nv
      real(kind=dp), dimension(N)            :: uloc, u0, profile_art_i, u1, u0sq
      ! real(kind=dp), dimension(N,NvspaceMax) :: u1


      Nvspace = NvspaceMax
      i = line%i; j = line%j
      Nred = line%Nr; Nblue = line%Nb
      vth = vbroad(T(icell),line%Atom%weight, vturb(icell))

      u0(:) = (lambda - line%lambda0)/line%lambda0  * ( c_light/vth )
      uloc(:) = line%v(:) / vth

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
         !omegav_mean should be close to vlabs
         write(*,*) "vlabs=", vlabs*1d-3, omegav_mean * 1d-3
         stop
      else
         !recentre the cell on the speed of the non-lte cell,
         !so that the kinematics is computed with respect to a non-moving non-lte cell.
         !it is always 0 in LTE or image (when labs = .false.)
         omegav(1:Nvspace) = omegav(1:Nvspace) - vlabs
         ! write(*,*) "vlabs = ", vlabs * 1d-3
      endif


      if (line%voigt) then
         u1(:) = u0(:) - omegav(1)/vth
         profile_art_i(:) = linear_1D_sorted(N,uloc(:),line%phi(:,icell),N,u1(:))
         do nv=2, Nvspace
            u1(:) = u0(:) - omegav(nv)/vth
            profile_art_i(:) = profile_art_i(:) + linear_1D_sorted(N,uloc(:),line%phi(:,icell),N,u1(:))
         enddo

      else
         !u1 = (u0 - omegav(nv)/vth)**2
         u0sq(:) = u0(:)*u0(:)
         u1(:) = u0sq(:) + (omegav(1)/vth)*(omegav(1)/vth) - 2*u0(:) * omegav(1)/vth
         profile_art_i(:) = exp(-u1(:))
         do nv=2, Nvspace
            u1(:) = u0sq(:) + (omegav(nv)/vth)*(omegav(nv)/vth) - 2*u0(:) * omegav(nv)/vth
            profile_art_i(:) = profile_art_i(:) + exp(-u1(:))
         enddo
      endif

      profile_art_i(:) = profile_art_i(:) / Nvspace

      return
   end function profile_art_i

  !no level dissolution yet
  !fills also eta_atoms
  !check conditions on negative opac
!   subroutine cross_coupling_terms(id, icell, iray)
!     integer, intent(in) :: id, icell, iray
!     integer :: nact, j, i, kr, kc, Nb, Nr, la, Nl, icell_d
!     type (AtomType), pointer :: aatom
!     real(kind=dp) :: gij, wi, wj, chicc, wl,  wphi, ni_on_nj_star

!     !for one ray
!     Uji_down(:,:,:,id) = 0.0_dp
!     chi_down(:,:,:,id) = 0.0_dp
!     chi_up(:,:,:,id)   = 0.0_dp

!     eta_atoms(:,:,id) = 0.0_dp

!     aatom_loop : do nact=1, Nactiveatoms
!        aatom => ActiveAtoms(nact)%ptr_atom

!        cont_loop : do kr = aatom%Ntr_line+1, aatom%Ntr

!           kc = aatom%at(kr)%ik

!           j = aatom%continua(kc)%j
!           i = aatom%continua(kc)%i
!           Nb = aatom%continua(kc)%Nb; Nr = aatom%continua(kc)%Nr
!           Nl = Nr - Nb + 1

!           ! 					ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
!           ni_on_nj_star = aatom%nstar(i,icell)/(aatom%nstar(j,icell) + 1d-100)


!           icell_d = 1
!           if (ldissolve) then
!              if (aatom%ID=="H") icell_d = icell
!           endif

!           gij = ni_on_nj_star * exp(-hc_k/T(icell)/aatom%continua(kc)%lambda0)

!           if (aatom%n(i,icell) - gij*aatom%n(j,icell) <= 0.0_dp) then
!              cycle cont_loop
!           endif


!           freq_loop : do la=1, Nl
!              if (la==1) then
!                 wl = 0.5*(lambda(Nb+1)-lambda(Nb)) / lambda(Nb)
!              elseif (la==Nl) then
!                 wl = 0.5*(lambda(Nb+la-1)-lambda(Nb+la-2)) / lambda(Nb+la-1)
!              else
!                 wl = 0.5*(lambda(Nb+la)-lambda(Nb+la-2)) / lambda(Nb+la-1)
!              endif

!              gij = ni_on_nj_star * exp(-hc_k/T(icell)/lambda(Nb+la-1))

!              !small inversions
!              !chicc = wl * fourpi_h * aatom%continua(kc)%alpha_nu(la,icell_d) * abs(aatom%n(i,icell) - gij*aatom%n(j,icell))
!              chicc = wl * fourpi_h * aatom%continua(kc)%alpha_nu(la,icell_d) * (aatom%n(i,icell) - gij*aatom%n(j,icell))
!              ! 						if (chicc < 0.0) chicc = 0.0_dp !should not happend

!              Uji_down(Nb+la-1,j,nact,id) = Uji_down(Nb+la-1,j,nact,id) + &
!                   aatom%continua(kc)%alpha_nu(la,icell_d) * (twohc/lambda(Nb+la-1)**3) * gij

!              chi_down(Nb+la-1,j,nact,id) = chi_down(Nb+la-1,j,nact,id) + chicc

!              chi_up(Nb+la-1,i,nact,id) = chi_up(Nb+la-1,i,nact,id) + chicc
!              !check consistency with nlte b-f and how negative opac is handled
!              !if (chicc > 0.0_dp) &
!              eta_atoms(Nb+la-1,nact,id) = eta_atoms(Nb+la-1,nact,id) + &
!                   aatom%continua(kc)%alpha_nu(la,icell_d) * (twohc/lambda(Nb+la-1)**3)  * gij * aatom%n(j,icell)
!           enddo freq_loop

!        enddo cont_loop

!        !for each line eventually
!        wi = 1.0; wj = 1.0
!        if (ldissolve) then
!           if (aatom%ID=="H") then
!              !nn
!              wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1),hydrogen%n(1,icell))!1 for H
!              wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1),hydrogen%n(1,icell))
!           endif
!        endif

!        line_loop : do kr = 1, aatom%Ntr_line

!           kc = aatom%at(kr)%ik


!           j = aatom%lines(kc)%j
!           i = aatom%lines(kc)%i

!           if (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell) <= 0.0_dp) then
!              cycle line_loop
!           endif

!           Nb = aatom%lines(kc)%Nblue; Nr = aatom%lines(kc)%Nred

!           Nl = Nr-dk_min+dk_max-Nb+1
!           wphi = 0.0
!           freq2_loop : do la=1, Nl
!              if (la==1) then
!                 wl = 0.5*1d3*hv
!                 !wl = 0.5*(lambda(Nb+1)-lambda(Nb)) * clight / aatom%lines(kc)%lambda0
!              elseif (la==Nl) then
!                 wl = 0.5*1d3*hv
!                 !wl = 0.5*(lambda(Nr)-lambda(Nr-1)) * clight / aatom%lines(kc)%lambda0
!              else
!                 wl = 1d3*hv
!                 !wl = 0.5*(lambda(Nb+la+1)-lambda(Nb+la-1)) * clight / aatom%lines(kc)%lambda0
!              endif


!              Uji_down(Nb+dk_min-1+la,j,nact,id) = Uji_down(Nb+dk_min-1+la,j,nact,id) + &
!                   hc_fourPI * aatom%lines(kc)%Aji * aatom%lines(kc)%phi_loc(la,iray,id)

!              !small inversions
!              ! 					if (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell) >= 0.0_dp) then

!              chi_down(Nb+dk_min-1+la,j,nact,id) = chi_down(Nb+dk_min-1+la,j,nact,id) + &
!                   wl * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(la,iray,id) * &
!                   (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))

!              chi_up(Nb+dk_min-1+la,i,nact,id) = chi_up(Nb+dk_min-1+la,i,nact,id) + &
!                   wl * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(la,iray,id) * &
!                   (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))


!              ! 					endif

!              eta_atoms(Nb+dk_min-1+la,nact,id) = eta_atoms(Nb+dk_min-1+la,nact,id) + &
!                   hc_fourPI * aatom%lines(kc)%Aji * aatom%lines(kc)%phi_loc(la,iray,id) * aatom%n(j,icell)


!              wphi = wphi + wl * aatom%lines(kc)%phi_loc(la,iray,id)
!           enddo freq2_loop

!           chi_down(Nb+dk_min:Nr+dk_max,j,nact,id) = chi_down(Nb+dk_min:Nr+dk_max,j,nact,id) / wphi
!           chi_up(Nb+dk_min:Nr+dk_max,i,nact,id) = chi_up(Nb+dk_min:Nr+dk_max,i,nact,id) / wphi

!        enddo line_loop


!     enddo aatom_loop

!     return
!   end subroutine cross_coupling_terms

  !TO DO merge with opacity_atom_loc()
!   subroutine opacity_atom_zeeman_loc(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l_void_before,l_contrib,  iterate)

!     integer, intent(in) :: id, icell, iray
!     logical, intent(in) :: iterate
!     real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l_void_before,l_contrib
!     integer :: nact, Nred, Nblue, kc, kr, i, j, m
!     type(AtomType), pointer :: aatom
!     integer :: dk0, dk1, Nlam
!     real(kind=dp) :: wi, wj, adamp, dv_dl, chil, etal
!     real(kind=dp), dimension(Nlambda_max_line+2*dk_max) :: phi0
!     real(kind=dp), dimension(Nlambda_max_line+2*dk_max,3) :: phiZ, psiZ

!     chiQUV_p(:,:,id) = 0.0_dp
!     etaQUV_p(:,:,id) = 0.0_dp
!     rho_p(:,:,id) = 0.0_dp

!     atom_loop : do nact = 1, Natom
!        aatom => Atoms(nact)%ptr_atom

!        tr_loop : do kr = 1,aatom%Ntr_line

!           kc = aatom%at(kr)%ik

!           Nred = aatom%lines(kc)%Nred; Nblue = aatom%lines(kc)%Nblue
!           i = aatom%lines(kc)%i;j = aatom%lines(kc)%j


!           Nblue = Nblue + dk_min
!           Nred = Nred + dk_max
!           Nlam = Nred - Nblue + 1

!           wj = 1.0; wi = 1.0
!           if (ldissolve) then
!              if (aatom%ID=="H") then
!                 !nn
!                 wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1),hydrogen%n(1,icell)) !1 for H
!                 wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1),hydrogen%n(1,icell))
!              endif
!           endif

!           if ((aatom%n(i,icell)*wj/wi - aatom%n(j,icell)*aatom%lines(kc)%gij) <= 0.0_dp) then
!              cycle tr_loop
!           endif

!           !if not init, bug when line is not polarizable why
!           phiz = 0.0
!           psiz = 0.0
!           !!fixed at the moment
!           call local_profile_zv(aatom%lines(kc),icell,iterate,Nlam,&
!           lambda(Nblue:Nred),phi0(1:Nlam),phiZ(1:Nlam,:), psiZ(1:Nlam,:), x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)

!           etal = hc_fourPI * aatom%lines(kc)%Aji * aatom%n(j,icell)
!           chil = hc_fourPI * aatom%lines(kc)%Bij * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))

!           ! 				if ((aatom%n(i,icell)*wj/wi - aatom%n(j,icell)*aatom%lines(kc)%gij) > 0.0_dp) then


!           chi(Nblue:Nred,id) = chi(Nblue:Nred,id) + chil * phi0(1:Nlam)
!           eta(Nblue:Nred,id)= eta(Nblue:Nred,id) + etal * phi0(1:Nlam)
          
!           if (minval(chi(Nblue:Nred,id)) < 0) then
!           	write(*,*) id, icell, chil
!           	write(*,*) "phi=", phi0(1:Nlam)
!           	call error("chi neg!")
!           endif

!           do m=1,3
!              chiQUV_p(Nblue:Nred,m,id) = chiQUV_p(Nblue:Nred,m,id) + chil * phiz(1:Nlam,m)
!              etaQUV_p(Nblue:Nred,m,id) = etaQUV_p(Nblue:Nred,m,id) + etal * phiz(1:Nlam,m)
!              rho_p(Nblue:Nred,m,id) = rho_p(Nblue:Nred,m,id) + chil * psiz(1:Nlam,m)
!           	 if (any_nan_infinity_vector(phiz(:,m)) > 0) then
!           	 	write(*,*) "(phiz)", icell, m, Nlam, chil, etal
!           	 	write(*,*) phiz(:,m)
!           	 	stop
!           	 endif
!           	 if (any_nan_infinity_vector(psiz(:,m)) > 0) then
!           	 	write(*,*) "(psiz)", icell, m, Nlam, chil, etal
!           	 	write(*,*) psiz(:,m)
!           	 	stop
!           	 endif
! !           	 if (aatom%lines(kc)%polarizable) then
! !           	 write(*,*) "m=", m, size(phiz)
! !           	 write(*,*) "phi0=",maxval(phi0)
! !           	 write(*,*) "phiz(m)", maxval(phiz(:,m))
! !           	 write(*,*) "psiz(m)=", maxval(psiz(:,m))
! !           	 endif
!           enddo

!           ! 				else !neg or null
!           ! 					eta(Nblue:Nred,id)= eta(Nblue:Nred,id) + etal * phi0(1:Nlam)
!           ! 					do m=1,3
!           ! ! 						chiQUV_p(Nblue:Nred,m,id) = chiQUV_p(Nblue:Nred,m,id) + chil * phiz(1:Nlam,m)
!           ! 						etaQUV_p(Nblue:Nred,m,id) = etaQUV_p(Nblue:Nred,m,id) + etal * phiz(1:Nlam,m)
!           ! ! 						rho_p(Nblue:Nred,m,id) = rho_p(Nblue:Nred,m,id) + chil * psiz(1:Nlam,m)
!           ! 					enddo
!           ! 				endif


!        end do tr_loop

!        aatom => NULL()

!     end do atom_loop


!     return
!   end subroutine opacity_atom_zeeman_loc
  
!   subroutine solve_stokes_mat(id,icell,la,Q,P)
! 	!Fill the modified Stokes matrix K'
! 	!and computes the matrix R  and the vector P
! 	!such that dot(I,R) = P.
! 	!Then, solve for IR = P, and get P the new
! 	!Stokes vector = (I,Q,U,V)^t
! 	!
! 	!Q is the atenutation factor for the local point
! 	!For wavelength index la.
!   	integer, intent(in) :: id, la,icell
!   	real(kind=dp), intent(in) :: Q
!   	real(kind=dp), intent(inout) :: P(4)
!   	real(kind=dp) :: R(4,4)
 
!   	!    / eta/chi  \
!   	!    | etaQ/chi |
!   	!S = | etaU/chi |
!   	!    \ etaV/chi / 	
!   	!P = Svect * Q
  	
!   	!Absorption-dispersion matrix
!   	!    /  chi  chiQ  chiU   chiV \
!   	!    |                         |
!   	!    |  chiQ  chi  rhoV  -rhoU |
!   	!KK =|                         |
!   	!    |  chiU -rhoV  chi   rhoQ |
!   	!    |                         |
!   	!    \  chiV  rhoU  -rhoQ  chi /
  	
!   	!KK' = K/chi - eye(4); diag(KK') = 0
!   	!R = eye(4) + Q*KK'; diag(R) = 1.0
!   	!

! 	!P(1) already filled with Snu * Q
! 	P(2) = etaQUV_p(la,1,id)/chi(la,id) * Q!SQ exp(-tau) (1 - exp(-dtau))
! 	P(3) = etaQUV_p(la,2,id)/chi(la,id) * Q!SU exp(-tau) (1 - exp(-dtau))
! 	P(4) = etaQUV_p(la,3,id)/chi(la,id) * Q!SV exp(-tau) (1 - exp(-dtau))
	
!     if (any_nan_infinity_vector(P) > 0) then
!         write(*,*) "(nan/inf P)", " icell=",icell, " la=",la, " lam=",lambda(la)
!         write(*,*) 'Q=',Q
!         write(*,*) 'P=',P
!         write(*,*) 'chi=',chi(la,id)
!     	stop
!     endif

! ! 	R = 0.0	 
! 	!eye(4)
! 	R(1,1) = 1.0
! 	R(2,2) = 1.0
! 	R(3,3) = 1.0
! 	R(4,4) = 1.0
	
! 	R(1,2) = chiQUV_p(la,1,id)/chi(la,id) * Q
! 	R(1,3) = chiQUV_p(la,2,id)/chi(la,id) * Q
! 	R(1,4) = chiQUV_p(la,3,id)/chi(la,id) * Q
	
! 	R(2,1) = R(1,2)
! 	R(2,3) = rho_p(la,3,id)/chi(la,id) * Q
! 	R(2,4) = -rho_p(la,2,id)/chi(la,id) * Q
	
! 	R(3,1) = R(1,3)
! 	R(3,2) = -R(2,3)
! 	R(3,4) = rho_p(la,1,id)/chi(la,id) * Q
	
! 	R(4,1) = R(1,4)
! 	R(4,2) = -R(2,4)
! 	R(4,3) = -R(3,4)

!     if (any_nan_infinity_matrix(R) > 0) then
!         write(*,*) "(2)", icell, la, lambda(la),R
!         stop
!     endif

! ! 	call solve_lin(R,P,4,.true.)
! 	call invert_4x4(R)
!     if (any_nan_infinity_matrix(R) > 0) then
!         write(*,*) "(3)",icell, la, lambda(la), R
!         stop
!     endif

! 	P = matmul(R,P)
! 	if (any_nan_infinity_vector(P) > 0) then
!         write(*,*) "(4)", icell, la, lambda(la), P
!         stop
! 	endif


! !     Stokes_Q(la,id) += P(2)
! !     Stokes_U(la,id) += P(3)
! !     Stokes_V(la,id) += P(4)

!   	return
!   end subroutine

   subroutine write_opacity_emissivity_bin(Nlambda,lambda)
      !not para to be able to write while computing!
      integer, intent(in) :: Nlambda
      real(kind=dp), intent(in) :: lambda(Nlambda)
      integer :: unit, unit2, status = 0
      integer :: alloc_status, id, icell, m, Nrec
      integer(kind=8) :: Nsize
      real(kind=dp), allocatable, dimension(:,:,:) :: chi_tmp, eta_tmp, rho_tmp
      character(len=50) :: filename_chi="chi_map.bin"
      character(len=50) :: filename_eta="eta_map.bin"  
      character(len=50) :: filename_rho="magnetoopt_map.bin"  
      
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
         if (icompute_atomRT(icell)) then
            call background_continua_lambda(icell, Nlambda, lambda, chi_tmp(:,icell,1), eta_tmp(:,icell,1))
            !Snu = Snu + scat(lambda, icell) * Jnu(:,icell)
            !accumulate b-f and b-b un chi and Sny
            call opacity_atom_bf_loc(icell,Nlambda,lambda,chi_tmp(:,icell,1), eta_tmp(:,icell,1))
            call opacity_atom_bb_loc(id,icell,1,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,-1d0,&
               0d0,0d0,.true.,Nlambda,lambda,chi_tmp(:,icell,1), eta_tmp(:,icell,1))
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
