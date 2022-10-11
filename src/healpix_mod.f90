! Few routines taken from: (free)
! Copyright (C) 1997-2013 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday,
!                          Matthias Bartelmann, Hans K. Eriksen,
!                          Frode K. Hansen, Martin Reinecke
!
! Routines adapted for mcfost (mostly from C distribution)
!
! This routine does not do the mesh refinement of the pixels.
! This is just a module with a couple of routines needed for pixelisation.
! It defines a structure healpix_map which contains the informations needed for the
! pixelisation of the sphere with non-uniform pixel sizes (adaptive pixelisation).
!
!
! It also provides routines to easy find the heirs of a pixel or to find the
! ancestor and the brothers and the parents of a pixel.
!
!
module healpix_mod

	use mcfost_env, only	: dp
	use constantes, only	: pi
	use parametres, only	: n_cells

	implicit none

	integer, private, dimension(0:1023) :: pix2y, pix2x
	integer, private, dimension(0:127) :: x2pix, y2pix

  ! coordinate of the lowest corner of each face
	integer, private, parameter, dimension(0:11) :: jrll = (/ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 /) ! in unit of nside
	integer, private, parameter, dimension(0:11) :: jpll = (/ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 /) ! in unit of nside/2

	type healpix_map
		integer, dimension(:), allocatable :: p
		integer, dimension(:), allocatable :: l
		integer :: npix !for this cell
		!npix is sum_{j=1,12*4**lbase} [ 4**(l(j) - lbase) ]
		!if (l(j) == lbase) sum is 12*4**lbase
	end type healpix_map

	type (healpix_map), dimension(:), allocatable :: healpix_cell

	contains

	!find common ancestor in resolution lend (< lstart) and siblings to a pixel
	!in resolution lstart.
	!subroutine
	!
	!return
	!end subroutine

	!C indexes
	!lend must be > lstart and Nheirs = 4**(dl)
	subroutine healpix_pix_heirs(Nheirs, pix, lstart, lend, p_children)
	!For parent pixel pix in the resolution parameter lstart, gives p_children
	!heirs (Nheirs = 4**(lstart-lend)) in resolution lend (lend > lstart).
	!(if dl=1 the heirs are the children)
		integer, intent(in) :: pix, lstart, lend, Nheirs
		integer, dimension(4) :: children
		integer :: dl, Npix, l, k, kpix, i
		integer, dimension(Nheirs), intent(out) :: p_children
		integer, dimension(Nheirs) :: pix_tmp

		dl = lend - lstart
		!if (dl <= 0) then
		!	write(*,*) " Error!"
		!endif
		Npix = 4**dl !number of heirs
		!Nheirs = Npix !!

		p_children(1:Npix) = -1
		p_children(1) = pix
		pix_tmp(1:Npix) = -1

		do l=lstart+1,lend
			k = 1
			infinity : do
				kpix = p_children(k)
				if (kpix == -1) exit infinity

				children(:) = healpix_children_nested_index(kpix)

				do i=1,4

					pix_tmp(i + (k-1) * 4) = children(i)

				enddo

				k = k + 1

			enddo infinity

			p_children(1:4*k) = pix_tmp(1:4*k)
			pix_tmp(1:4*k) = -1

		enddo

	return
	end subroutine healpix_pix_heirs

	subroutine healpix_pix_heirs_new(pix, lstart, lend, p_children)
	!For parent pixel pix in the resolution parameter lstart, gives p_children
	!heirs (Nheirs = 4**(lstart-lend)) in resolution lend (lend > lstart).
		integer, intent(in) :: pix, lstart, lend
		integer :: Nheirs !number of children in the highest resolution for pix in the lowest resolution
		integer, dimension(4) :: children
		integer :: l, k, kpix, i
		integer, dimension(4**(lstart - lend)), intent(out) :: p_children
		integer, dimension(4**(lstart - lend)) :: pix_tmp

		Nheirs = 4**(lstart - lend)

		p_children(1:Nheirs) = -1
		p_children(1) = pix
		pix_tmp(1:Nheirs) = -1

		do l=lstart+1,lend
			k = 1
			infinity : do
				kpix = p_children(k)
				if (kpix == -1) exit infinity

				children(:) = healpix_children_nested_index(kpix)

				do i=1,4

					pix_tmp(i + (k-1) * 4) = children(i)

				enddo

				k = k + 1

			enddo infinity

			p_children(1:4*k) = pix_tmp(1:4*k)
			pix_tmp(1:4*k) = -1

		enddo

	return
	end subroutine healpix_pix_heirs_new

	function mean_number_of_pixels()
		integer :: npix_mean
		integer :: icell, nn
		real :: mean_number_of_pixels

		npix_mean = 0
		nn = 0
		do icell=1,n_cells
			if (allocated(healpix_cell(icell)%p)) then
				npix_mean = npix_mean + healpix_cell(icell)%npix
				nn = nn + 1
			endif
		enddo

		mean_number_of_pixels = real(npix_mean) / real(nn)

	return
	end function mean_number_of_pixels

	subroutine allocate_healpix_map_cell(l,hpm_cell)
		integer, intent(in) :: l
		type (healpix_map), intent(inout) :: hpm_cell
		integer :: Npix, pp

		Npix = 12*4**l

		allocate(hpm_cell%p(Npix)); hpm_cell%p(:) = -1
		do pp=1,Npix
			hpm_cell%p(pp) = pp-1 !C indexes
		enddo
		allocate(hpm_cell%l(Npix)); hpm_cell%l(:) = l !same resolution
		 hpm_cell%npix = 0 !computed after refinement

	return
	end subroutine allocate_healpix_map_cell

	subroutine deallocate_healpix_map_cell(hpm_cell)
		type (healpix_map), intent(inout) :: hpm_cell

		!not allocated if cell is empty
		if (allocated(hpm_cell%p)) then
			deallocate(hpm_cell%p, hpm_cell%l)
		endif
		hpm_cell%npix = 0

	return
	end subroutine deallocate_healpix_map_cell

	subroutine init_healpix_map(l,mask)
		integer, intent(in) :: l
		integer :: icell, Npix
		logical, intent(in), optional :: mask(n_cells)


		Npix = 12*4**l

		allocate(healpix_cell(n_cells))

		if (present(mask)) then

			do icell=1,n_cells

				if (mask(icell)) then
					call allocate_healpix_map_cell(l,healpix_cell(icell))
				endif

			enddo

		else

			do icell=1,n_cells

				call allocate_healpix_map_cell(l,healpix_cell(icell))

			enddo

		endif

	return
	end subroutine init_healpix_map

	subroutine deallocate_healpix_map()
		integer :: icell

		do icell=1,n_cells

			call deallocate_healpix_map_cell(healpix_cell(icell))

		enddo

		deallocate(healpix_cell)

	return
	end subroutine deallocate_healpix_map

	subroutine init_pix2xy_and_xy2pix()
		integer :: kpix, jpix, ix, iy, ip, id
		integer :: i, j, k

		!first pix2x and 2y arrays

		pix2y(:) = 0
		pix2x(:) = 0

		do kpix=0,1023
			jpix = kpix
			ix = 0
			iy = 0
			ip = 1                     ! bit position (in x and y)
			infinity : do
				if ( jpix == 0) exit infinity   ! go through all the bits
				id = modulo(jpix,2)       ! bit value (in kpix), goes in ix
				jpix = jpix/2
				ix = id*ip+ix

				id = modulo(jpix,2)       !bit value (in kpix), goes in iy
				jpix = jpix/2
				iy = id*ip+iy

				ip = 2*ip             ! next bit (in x and y)
			enddo infinity
			pix2x(kpix) = ix          ! in 0,31
			pix2y(kpix) = iy          ! in 0,31
		enddo

		x2pix(:) = 0
		y2pix(:) = 0

		!now x and y 2pix
		do i=0,127					  !for converting x,y into
			j  = i                    !pixel numbers
			k  = 0
			ip = 1
			infinity_b : do
				if (j == 0) then
					x2pix(i) = k
					exit infinity_b
				else
					id = modulo(j, 2)
					j  = j/2
					k  = ip*id+k
					ip = ip*4
				endif
			enddo infinity_b
		enddo

		y2pix(:) = 2 * x2pix(:)


	return
	end subroutine init_pix2xy_and_xy2pix

	subroutine pix2xyf(l, pix, ix, iy, iface)
	!C index pix
		integer, intent(in) :: l, pix
		integer :: nside, npface
		integer, intent(out) :: ix, iy, iface
		integer :: ipf_tmp, scale, ip_low, smax, i
		integer :: ip_trunc, ip_med, ip_hi

		nside = 2**l
		npface = nside*nside
		iface = pix / npface
		scale = 1

		ix = 0
		iy = 0
		ipf_tmp = pix

		if (nside <= 8192) then !l<=13
			smax = 2
		else
			smax = 5
		endif
!
!
! 		do i=0,smax-1
! 			ip_low = iand(ipf_tmp,1023)
! 			ix = ix + scale * pix2x(ip_low)
! 			iy = iy + scale * pix2y(ip_low)
! 			scale = scale * 32
! 			ipf_tmp = ipf_tmp / 1024
! 		enddo
!
! 		ix = ix + scale * pix2x(ipf_tmp)
! 		iy = iy + scale * pix2y(ipf_tmp)

		if (nside <= 8192) then
			ip_low = iand(ipf_tmp,1023)   ! content of the last 10 bits
			ip_trunc =    ipf_tmp/1024    ! truncation of the last 10 bits
			ip_med = iand(ip_trunc,1023)  ! content of the next 10 bits
       		ip_hi  =      ip_trunc/1024   ! content of the high weight 10 bits

			ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low)
			iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low)
! 			write(*,*) "ix", ix, ip_hi, ip_med, ip_low
! 			write(*,*) "iy", iy

		else
			ix = 0
			iy = 0
			scale = 1
			do i=0, smax-1
				ip_low = iand(ipf_tmp,1023)
				ix = ix + scale * pix2x(ip_low)
				iy = iy + scale * pix2y(ip_low)
				scale = scale * 32
				ipf_tmp   = ipf_tmp/1024
			enddo
			ix = ix + scale * pix2x(ipf_tmp)
			iy = iy + scale * pix2y(ipf_tmp)
		endif

	return
	end subroutine pix2xyf

	subroutine xyf2pix(l,ix,iy,face,pix)
		integer, intent(in) :: l, ix, iy, face
		integer, intent(out) :: pix
		integer, parameter :: scale_factor = 128*128
		integer :: nside, ix_tmp, iy_tmp, smax, scale, npface
		integer :: ix_low, iy_low, i


		nside = 2**l

		ix_tmp = ix
		iy_tmp = iy
		if (nside > 8192) then
			smax = 4
		else
			smax = 1
		endif

		scale = 1
		pix = 0
		npface = nside*nside

		do i=0, smax-1
			ix_low = iand(ix_tmp,127)
			iy_low = iand(iy_tmp,127)
			pix = pix + (x2pix(ix_low) + y2pix(iy_low)) * scale
			scale = scale * scale_factor
			ix_tmp = ix_tmp / 128
			iy_tmp = iy_tmp / 128

		enddo
		pix = pix + (x2pix(ix_tmp) + y2pix(iy_tmp)) * scale
		pix = pix + face * npface ! in {0, 12*nside**2 - 1}

	return
	end subroutine xyf2pix

	subroutine healpix_nested_pix(l, mu, phi, pix)
		integer, intent(in) :: l
		real(kind=dp), intent(in) :: mu, phi
		integer, intent(out) :: pix
		real(kind=dp) :: za, tt,tp, tmp, theta
		integer :: nside, iface, ix, iy
		integer :: jm, jp, ifp, ifm, ntt

		nside = 2**l

		za = abs(mu)
		theta =  acos(mu) !in 0 and pi
		tt = modulo(phi, 2.0*pi) * 2.0/pi

		if (za <= 2.0/3.0) then
			jp = int(nside*(0.5 + tt - mu*0.75))
			jm = int(nside*(0.5 + tt + mu*0.75))

			!fin the face
			ifp = jp / nside
			ifm = jm / nside

			if (ifp == ifm) then !faces 4 to 7
				iface = iand(ifp,3)+4
! 				iface = ior(ifp,4)
			elseif (ifp < ifm) then !faces 0 to 3
				iface = iand(ifp,3)
! 				iface = ifp
			else !faces 8 to 11
				iface = iand(ifm,3) + 8
! 				iface = ifm + 8
			endif

			ix = iand(jm, nside-1)
			iy = nside - iand(jp,nside-1) - 1

		else
			ntt = int(tt)
			if (ntt >=4) ntt = 3
			tp = tt - real(ntt)
			tmp = nside * sqrt(3 * (1.0-za))

			jp = int(tp*tmp)
			jm = int((1.0 - tp)*tmp)
! 			if (jp>=nside) jp = nside-1
! 			if (jm>=nside) jm = nside-1
			jp = min(nside-1,jp) !edge points
			jm = min(nside-1,jm)

			if (mu>=0) then
				iface = ntt
				ix = nside - jm - 1
				iy = nside - jp - 1
			else
				iface = ntt + 8
				ix = jp
				iy = jm
			endif

		endif

		!!write(*,*) "mod=", tt," face_num=", iface

		call xyf2pix(l,ix,iy,iface,pix)

! 		if (za <= 2.0/3.0) then !equator
! 			jp = int(nside*(0.5 + tt - mu*0.75))
! 			jm = int(nside*(0.5 + tt + mu*0.75))
!
! 			!fin the face
! 			ifp = jp / nside
! 			ifm = jm / nside
!
! 			if (ifp == ifm) then !faces 4 to 7
! 				iface = iand(ifp,3)+4
! 			elseif (ifp < ifm) then !faces 0 to 3
! 				iface = iand(ifp,3)
! 			else !faces 8 to 11
! 				iface = iand(ifm,3) + 8
! 			endif
!
! 			ix = iand(jm, nside-1)
! 			iy = nside - iand(jp,nside-1) - 1
!
! 		else !polar region, abs(mu) > 2/3
! 			ntt = int(tt)
! 			if (ntt >= 4) ntt = 3
! 			tp = tt - real(ntt)
! 			if (mu > 0.0) then
! 				tmp = sqrt(6.0) * sin(theta * 0.5)
! 			else
! 				tmp = sqrt(6.0) * cos(theta * 0.5)
! 			endif
!
! 			jp = int(nside * tp * tmp)
! 			jm = int(nside * (1.0 - tp)*tmp)
! 			jp = min(nside-1,jp) !edge points
! 			jm = min(nside-1,jm)
!
! 			if (mu >= 0) then
! 				iface = ntt
! 				ix = nside-jm-1
! 				iy = nside-jp-1
! 			else
! 				iface = ntt + 8
! 				ix = jp
! 				iy = jm
! 			endif
!
! 		endif
!
! 		if (nside <= 8192) then
! 			ix_low = iand(ix, 127)
! 			iy_low = iand(iy, 127)
! 			ipf =     x2pix1(ix_low) + y2pix1(iy_low) &
!             & + (x2pix1(ix/128) + y2pix1(iy/128)) * 16384
! 		else
! 			scale = 1
! 			scale_factor = 16384 ! 128*128
!        		ipf = 0
!        		ismax = 1 ! for nside in [2^14, 2^20]
!        		if (nside >  1048576 ) ismax = 3
!       		do i=0, ismax
!           		ix_low = iand(ix, 127) ! last 7 bits
!           		iy_low = iand(iy, 127) ! last 7 bits
!           		ipf = ipf + (x2pix1(ix_low)+y2pix1(iy_low)) * scale
!           		scale = scale * scale_factor
!           		ix  =     ix / 128 ! truncate out last 7 bits
!           		iy  =     iy / 128
!       		 enddo
!        		ipf =  ipf + (x2pix1(ix)+y2pix1(iy)) * scale
!     	endif
!     	pix = ipf + face_num* int(nside,MKD) * nside    ! in {0, 12*nside**2 - 1}

	return
	end subroutine healpix_nested_pix

	!pix from 0 to Npix - 1 !C indexes
	subroutine healpix_nested_mu_and_phi (l, pix, mu, phi)
		integer, intent(in) :: l, pix
		real(kind=dp), intent(out) :: mu, phi
		integer :: nside, nl4, npix, iface, ix, iy, nr
		integer :: jr, kshift, jp, npface, ipf, jrt, jpt
		real :: fact2, fact1

		nside = 2**l
		nl4 = nside * 4
		npix = 12*nside*nside
		npface = nside * nside
		!!iface = pix / npface
		ipf = modulo(pix,npface)

		fact2 = 4.0/real(npix)

		call pix2xyf(l, ipf, ix, iy, iface)
		iface = pix / npface !otherwise does not work with ipf
		jrt = ix + iy
		jpt = ix - iy
		jr = jrll(iface)*nside - jrt - 1

		if (jr < nside) then
			nr = jr
			mu = 1.0 - nr*nr*fact2
			kshift = 0
		elseif (jr > 3*nside) then
			nr = nl4 - jr
			mu = nr*nr*fact2 - 1.0
			kshift = 0
		else
			fact1 = lshift(nside,1) * fact2
			nr = nside
			mu = (2*nside-jr)*fact1
			kshift = iand((jr - nside),1)
		endif

		jp = (jpll(iface)*nr + jpt + 1 + kshift) / 2
		if (jp > nl4) jp = jp - nl4
		if (jp < 1) jp = jp + nl4
		phi = (jp - (kshift+1)*0.5) * 0.5 * pi / real(nr)

! 		jp  = jpll(iface)*nr + jpt
! 		if (jp < 0) jp = jp + 2 * nl4
! 		phi = jp * 0.25 * pi / nr


	return
	end subroutine healpix_nested_mu_and_phi


	function int_pow(a,n)
		integer, intent(in) :: a, n
		integer :: int_pow

		int_pow = a**n

	return
	end function int_pow

	!largest_integer_smaller_than_x
	function healpix_listx(x)!largest_integer_smaller_than_x(x)
  		real(kind=dp), intent(in) :: x
  		integer(kind=8) healpix_listx

  		healpix_listx = int(ceiling(x)) - 1

  	return
	end function healpix_listx

	function healpix_children_nested_index(pix)
  		!given the index of the parent pixel pix
  		!computes the 4 children pixel indexes
  		integer, parameter :: n_children = 4
  		integer, intent(in) :: pix
  		integer, dimension(n_children) :: healpix_children_nested_index
  		integer :: n

  		do n=0,n_children-1
  			healpix_children_nested_index(n+1) = 4*pix + n
  		enddo

  		return
	end function healpix_children_nested_index

	function healpix_parent_nested_index(m)
		integer, intent(in) :: m
		integer :: healpix_parent_nested_index

  	!m=0,1,2,3 have the same parent
  	!m=4,5,6,7 have the same parent...
		healpix_parent_nested_index = int(real(m)/4.0)

		return
	end function healpix_parent_nested_index

	function healpix_npix(l)
		integer, intent(in) :: l
		integer(kind=8) :: healpix_npix

  	!Npix = 12 * Nside**2
  	!with Nside = 2**l, the resolution of the grid
		healpix_npix = 12 * 4**l

		return
	end function healpix_npix

	function healpix_weight(l)
		integer, intent(in) :: l
		real(kind=dp) :: healpix_weight

  	!weight = area/4pi
  	!area is 4pi / Npix
  	!Npix = 12 * Nside**2
  	!Nside = 2**l
  	!area = pi/3/Nside/Nside = pi/3/4**l
  	!weight = (pi/3)/(4pi)/4**l
  	!weight = (1./12.)/4**l
		healpix_weight = (1.0_dp/12.0_dp) / real(4**l,kind=dp)

		return
	end function healpix_weight

	function healpix_angular_resolution(l)
  	!in degrees
		integer, intent(in) :: l
		real(kind=dp) :: healpix_angular_resolution

  	!4pi * healpix_weight is actually the solid angle it's the "square" resolution
  	!healpix_angular_resolution = 180.0 * sqrt( 4.0 * pi * healpix_weight(l) )/ pi
		healpix_angular_resolution = 180.0 * sqrt(pi/real(3*4**l)) / pi

		return
	end function healpix_angular_resolution

    !pix from 1 to Npix (fortran indexes)
	subroutine healpix_ring_mu_and_phi(l, pix, mu, phi)
		integer, intent(in) :: l, pix
		integer :: p, i, j, pp
		real(kind=dp), intent(out) :: mu, phi !cos(theta) and azimuth from the north pole
		integer(kind=8) :: Npix, Nside, North_cap
		real(kind=dp) :: four_on_Npix, ph, fodd, s

		Npix = healpix_npix(l)
		Nside = 2**l
		North_cap = 2 * Nside * (Nside-1)

  	!C index
		p = pix - 1

  	!better to do that out of the subroutine if called many times
  	!if pix in a do loop (do pix=1,Npix) should never happen
		if (pix > Npix) then
			write(*,*) pix, Npix, l
			write(*,*) "(Healpix_ring_mu_and_phi) pixel index larger than Npix"
			stop
		endif

		four_on_Npix = 4.0_dp / real(Npix,kind=dp)

		if (p < North_cap) then !north
			ph = 5d-1 * real(1 + p,kind=dp)

  		! i = largest_integer_small_than_x(ph)
  			i = rshift(int(1 + floor(sqrt(real(1+2*p,kind=dp)))),1)

			j = (p+1) - 2 * i * (i-1)

			mu = 1.0_dp - real(i*i,kind=dp) * four_on_Npix
			phi = (real(j,kind=dp) - 0.5_dp) * 0.5_dp * pi / real(i,kind=dp)

		elseif (p < (Npix - North_cap)) then !equator
			pp = p - North_cap
			i = pp/(4*Nside) + Nside !integer
			j = mod(pp,4*Nside) + 1

		!1 if i+Nside odd, 0.5 otherwise
! 			if (mod(i+Nside,2)==0) then
			if (iand(i+Nside,1_8)==0) then
				fodd = 0.5
			else !odd
				fodd = 1.0
			endif
			s = 0.5 * mod(i-Nside+1,2)
! 			write(*,*) "fodd", fodd, s
! 			fodd = s!, should give the same result

		!!see, it is simply 4/3
		!!(2*Nside * 4.0 / Npix ) * (Nside <<1) = 4/3 whatever Nside
		!!(4/Npix) * (Nside << 1) = 2/3/Nside whatever Nside
! 		mu = four_on_Npix * real( ( 2 * Nside - i) * lshift(Nside,1) )
			mu = 4.0_dp / 3.0_dp - 2.0_dp * real(i,kind=dp) / real(3*Nside,kind=dp)
			phi = (real(j,kind=dp) - fodd) * pi / real(2 * Nside,kind=dp)

		else !south
			pp = Npix - p
			i = rshift(int(1 + floor(sqrt(real(2*pp-1,kind=dp)))),1)
			j = 4 * i + 1 - (pp - 2*i*(i-1))

			mu = -1.0_dp + real(i*i,kind=dp) * four_on_Npix
			phi = (real(j,kind=dp) - 0.5_dp) * 0.5_dp * pi / real(i,kind=dp)

		endif

	return
	end subroutine healpix_ring_mu_and_phi

	subroutine healpix_ring_northern_pixels(l, pixs)
		integer, intent(in) :: l
		integer, dimension(:), intent(out) :: pixs
		integer :: p
		integer(kind=8) :: Npix, Nside, North_cap, Nupper

		Npix = healpix_npix(l)
		Nside = 2**l
		North_cap = 2 * Nside * (Nside-1)

		Nupper = Npix - North_cap
		if (size(pixs) /= Nupper) then
			write(*,*) "error Nupper size (pixs)"
		endif

		do p=0,Npix-1

			if (p < North_cap) then !north
				pixs(p+1) = p+1
			elseif (p < (Npix - North_cap)) then !equator
				pixs(p+1) = p+1
			endif

		enddo

	return
	end subroutine healpix_ring_northern_pixels

	subroutine healpix_sphere(l,mu,xmu,ymu)
		integer, intent(in) :: l
  		real(kind=dp), dimension(12*4**l), intent(out) :: xmu, ymu, mu
  		real(kind=dp) :: z, phi
  		integer :: pix


		call init_pix2xy_and_xy2pix
  		do pix=1,healpix_npix(l)

!   			call healpix_ring_mu_and_phi(l,pix,z,phi)
  			call healpix_nested_mu_and_phi (l, pix-1, z, phi)
  			mu(pix) = z
  			xmu(pix) = sqrt(1.0_dp - z*z) * cos(phi)
  			ymu(pix) = sqrt(1.0_dp - z*z) * sin(phi)

  		enddo

  	return
	end subroutine healpix_sphere


end module healpix_mod
