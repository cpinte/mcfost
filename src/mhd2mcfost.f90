!Interface between MCFOST and grid-base codes.
!Includes interface between stellar atmosphere codes.

module mhd2mcfost

    use constantes
    use grid, only : cell_map, vfield_x, vfield_y, vfield_z, alloc_atomrt_grid, nHtot, ne, v_char, lmagnetized, vturb, T, icompute_atomRT, lcalc_ne
    use parametres
    use messages
    use utils
    use sph2mcfost, only : SPH_to_Voronoi, Hydro_to_Voronoi_atomic
    use sort, only : find_kth_smallest_inplace
    use elements_type

    implicit none

    real(kind=dp) :: Taccretion
    logical :: laccretion_shock

    contains


    subroutine read_pluto()
        !read pluto format in hdf5
    
        call error("Pluto interface not available yet!")
    
        return
    end subroutine read_pluto

    subroutine setup_mhd_to_mcfost()
        !here lmhd_voronoi ist true even with model ascii !
        integer                                  :: n_points,k,icell
        integer                                  :: Nread
        real(kind=dp), allocatable, dimension(:) :: x,y,z,h!allocated in reading
        real(kind=dp), allocatable, dimension(:) :: vx,vy,vz,mass_gas, mass_ne_on_massgas, T_tmp, vt_tmp, dz
        real(kind=dp), allocatable, dimension(:) :: rho, hydro_grainsizes
        real(kind=dp), allocatable, dimension(:,:) :: rhodust, massdust
        integer,       allocatable, dimension(:) :: is_ghost
    
        integer, parameter                       :: Nheader = 3 !Add more, for ascii file
        integer                                  :: syst_status, acspot, alloc_status
        character(len=512)                       :: inputline, FormatLine, cmd
        logical                                  :: check_previous_tesselation
        integer,  allocatable, dimension(:)      :: particle_id
        real(kind=dp), dimension(6)              :: hydro_limits
        integer                                  :: ndusttypes !, voroindex, N_fixed_ne = 0
        real, parameter                          :: limit_factor = 1.005!, Lextent = 1.01

    
        write(FormatLine,'("(1"A,I3")")') "A", 512
    
        !There will be an error if lphantom_file is true. Because density_files and density_file
        !stores pluto's model name. But also the filename from phantom.. So to date, the two
        !codes cannot be merged.
        !This is to be able to use read previous tesselation
    
        lmagnetized = .false.
        lfix_star = .true.
        lphantom_file = .false.
        lignore_dust = .true.
        lrandomize_voronoi = .false.
        check_previous_tesselation = (.not.lrandomize_voronoi)
        n_points = 0 ! to avoid compiler warning
    
        if (lignore_dust) then
           ndusttypes = 0
           if (allocated(rhodust)) deallocate(rhodust,massdust)
        else
           call error("Dust not handled yet for pluto models!")
        endif
    
        if (lpluto_file) then
           write(*,*) "Voronoi tesselation on Pluto model..."
           !read and convert to mcfost units
           call read_pluto() ! to date empty
        elseif (lmodel_ascii) then

           cmd = "wc -l "//trim(density_file)//" > ntest.txt"
           call appel_syst(cmd,syst_status)
           open(unit=1,file="ntest.txt",status="old")
           read(1,*) N_points
           close(unit=1)
           N_points = N_points - Nheader
           write(*,*) " Input model has ", N_points," grid points !"
    
           N_points = N_points + n_etoiles
    
           open(unit=1,file=density_file, status="old")
           call read_line(1, FormatLine, inputline, Nread)
           !unused
        !    read(inputline(1:Nread),*) rotation_law
    
           !!lmagnetoaccr = .false.
           !!lspherical_velocity = .false.
           lvelocity_file = .false.
           lvfield_sphere_coord = .false.
           lvfield_cyl_coord = .true. !in mcfost, need lvelocity_file .false.
           !-> .false. with Voronoi
           !read T shock and if accretion spots
           call read_line(1, FormatLine, inputline, Nread)
           read(inputline(1:Nread),*) Taccretion, acspot
           laccretion_shock =  (acspot == 1)
           if (Taccretion==0.0_dp) Taccretion = -1.0_dp
    
           call read_line(1, FormatLine, inputline, Nread)
        !    read(inputline(1:Nread),*) thetai, thetao, tilt
        !    !for compatibility
        !    tilt = tilt * pi / 180.0
        !    thetai = thetai * pi / 180.0
        !    thetao = thetao * pi / 180.
    
           allocate(h(n_points), stat=alloc_status)
           if (alloc_status > 0) then
              call error("Allocation error smoothing length h")
           endif
           !cut cells larger than 3*h
    
           allocate(particle_id(n_points), stat=alloc_status)
           if (alloc_status > 0) then
              call error("Allocation error particle_id (mhd_to_mcfost)")
           endif
           particle_id(:) = 0
    
    
           allocate(x(n_points), y(n_points), z(n_points), stat=alloc_status)
           if (alloc_status > 0) then
              call error("Allocation error x, y, z")
           endif
    
           allocate(vx(n_points), vy(n_points), vz(n_points), stat=alloc_status)
           if (alloc_status > 0) then
              call error("Allocation error vx, vy, vz")
           endif
    
           allocate(T_tmp(n_points), vt_tmp(n_points), dz(n_points), mass_gas(n_points), mass_ne_on_massgas(n_points), &
                stat=alloc_status)
           if (alloc_status > 0) then
              call error("Allocation error T_tmp")
           endif
    
           icell = 0
           Nread = 0
           do icell=1, n_points
              particle_id(icell) = icell
              if (lmagnetized) then
                 call error("Magnetic field not available with Voronoi!")
              else
                 call read_line(1, FormatLine, inputline, Nread)
                 read(inputline(1:Nread),*) x(icell), y(icell), z(icell), T_tmp(icell), mass_gas(icell),&
                      mass_ne_on_massgas(icell), vx(icell), vy(icell), vz(icell), vt_tmp(icell), dz(icell), h(icell)
              endif
           enddo
           !density_file    
        else
           call error("lpluto_file or lmodel_ascii required for lmhd_voronoi!")
        end if

        hydro_limits(:) = 0
    
        k = 1
        hydro_limits(1) = find_kth_smallest_inplace(k,real(x))*limit_factor
        hydro_limits(3) = find_kth_smallest_inplace(k,real(y))*limit_factor
        hydro_limits(5) = find_kth_smallest_inplace(k,real(z))*limit_factor
    
        k = n_points
        hydro_limits(2) = find_kth_smallest_inplace(k,real(x))*limit_factor
        hydro_limits(4) = find_kth_smallest_inplace(k,real(y))*limit_factor
        hydro_limits(6) = find_kth_smallest_inplace(k,real(z))*limit_factor
        if (n_etoiles > 0) then
            write(*,*) "# Model limits (Rstar) :"
            write(*,*) "x =", hydro_limits(1)/etoile(1)%r, hydro_limits(2)/etoile(1)%r
            write(*,*) "y =", hydro_limits(3)/etoile(1)%r, hydro_limits(4)/etoile(1)%r
            write(*,*) "z =", hydro_limits(5)/etoile(1)%r, hydro_limits(6)/etoile(1)%r
        endif
    
        !also work with grid-based code
        !massdust, rhodust, hydro_grainsizes not allocated if ndusttypes = 0 !
        call sph_to_voronoi(n_points-n_etoiles, ndusttypes, particle_id, x, y, z, h, vx, vy, vz, &
             T_tmp, mass_gas, massdust, rho, rhodust, hydro_grainsizes, hydro_limits, check_previous_tesselation, is_ghost)
        ! -> correction for small density applied on mass_gas directly inside
    
        call hydro_to_Voronoi_atomic(n_points,T_tmp,vt_tmp,mass_gas,mass_ne_on_massgas,dz)
    ! 	call empty_cells
    
        !deallocating temporary variables from input file.
        deallocate(h,vx,vy,vz,mass_gas, mass_ne_on_massgas, x,y,z,T_tmp, vt_tmp, dz)
    
        return
      end subroutine setup_mhd_to_mcfost
    


    subroutine read_spheregrid_ascii(filename)
    ! ------------------------------------------- !
    ! Read from ascii file a model to be used.
    ! ------------------------------------------- !
        character(len=*), intent(in)	:: filename
        integer, parameter :: Nhead = 3 !Add more
        character(len=50) :: rotation_law
        integer :: icell, Nread, syst_status, N_points, k, i, j, acspot, N_fixed_ne = 0
        character(len=512) :: inputline, FormatLine, cmd
        real(kind=dp) :: rr, zz, pp, Vmod

        call alloc_atomrt_grid
        call read_abundance
    
        lvfield_sphere_coord = .false.
        lVoronoi = .false.
        !deactivated at the moment. I'll put back Zeeman pol later
        lmagnetized = .false.  
        lvelocity_file = .false.!needed with the new flag
        lvfield_cyl_coord = .false. !in mcfost, need lvelocity_file .false.
    
        write(FormatLine,'("(1"A,I3")")') "A", 512
    
        !could add header with magnetic field and so on
        !location of spots + lmagnetoaccretion flags if other kind of models with the use of spots
        ! + Tschok
    
        cmd = "wc -l "//trim(filename)//" > ntest.txt"
        call appel_syst(cmd,syst_status)
        open(unit=1,file="ntest.txt",status="old")
        read(1,*) N_points
        close(unit=1)
        !-N headers lines
        write(*,*) "Found ", N_points - Nhead, " points and grid has", n_cells, " points"
        if (N_points - Nhead/= n_cells) then
           call error( "Should read a model for the exact same grid as mcfost !" )
        end if
    
        open(unit=1,file=filename, status="old")
        call read_line(1, FormatLine, inputline, Nread)
        read(inputline(1:Nread),*) rotation_law
    
    
        select case (rotation_law)
        case ("magneto-accretion")
            lvfield_cyl_coord = .true.
           write(*,*) " Velocity law is ", trim(rotation_law), lvfield_cyl_coord
           if (lmagnetized) then
               call warning(" --> CHeck projection with lmagnetoaccr and Bfield!")
           endif
        case ("spherical_vector")
            lvfield_sphere_coord = .true.
           write(*,*) " Velocity law is ", trim(rotation_law), lvfield_sphere_coord
        case default
           write(*,*) " Velocity law ", rotation_law," not handled yet"
           stop
        end select
    
    
        !read T shock and if accretion spots
        call read_line(1, FormatLine, inputline, Nread)
        read(inputline(1:Nread),*) Taccretion, acspot
        laccretion_shock = .false.
        if (acspot==1) laccretion_shock = .true.
        if (Taccretion==0.0_dp) Taccretion = -1.0_dp
 
        call read_line(1, FormatLine, inputline, Nread)
        ! read(inputline(1:Nread),*) thetai, thetao, tilt
        ! thetai = thetai * pi/180.
        ! thetao = thetao * pi/180.
        ! tilt = tilt * pi/180.
    
        do i=1, n_rad
           do j=j_start,nz !j_start = -nz in 3D
              do k=1, n_az
                 if (j==0) then !midplane
                    !icell = cell_map(i,1,k)
                    cycle
                 else
                    icell = cell_map(i,j,k)
                 end if
                 Nread = 0
                 if (lmagnetized) then
                    stop
                    ! call getnextline(1, "#", FormatLine, inputline, Nread)
    
                    ! !In case of no polarisation, but magnetic field is present in the file, it is better to
                    ! !read the mandatory variables and put the magnetic field at the end of the file
                    ! read(inputline(1:Nread),*) rr, zz, pp, T(icell), nHtot(icell), ne(icell), &
                    !      vR(icell), V2(icell), Vphi(icell), vturb(icell), icompute_atomRT(icell), &
                    !      Bmag(icell), gammab(icell), chib(icell)
                    ! !BR(icell), B2(icell), Bphi(icell)
                    ! !Bmag(icell), gammab(icell), chib(icell)
                 else
                    call read_line(1, FormatLine, inputline, Nread)
                    read(inputline(1:Nread),*) rr, zz, pp, T(icell), nHtot(icell), ne(icell), &
                         vfield_x(icell), vfield_y(icell), vfield_z(icell), vturb(icell), icompute_atomRT(icell)

                 end if !magnetized
              end do
           end do
        end do
        close(unit=1)
        
        !rho -> nH
        nHtot = nHtot * 1d3 / masseH / wght_per_H
    
        write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT>0)), " density zones"
        write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT==0)), " transparent zones"
        write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT<0)), " dark zones"
    
        Vmod = maxval(sqrt(vfield_x**2+vfield_y(:)**2+vfield_z(:)**2))
    
        v_char = Vmod
    
        !to change according to what B quantities are read!!
        ! if (lmagnetized) then
        !    gammaB = gammaB * pi / 180.0
        !    chiB = chiB * pi / 180.0
        !    B_char = maxval(abs(Bmag))
        !    write(*,*)  "Typical Magnetic field modulus (G)", B_char * 1d4
        !    if (B_char <= 0.0_dp) then
        !       !    		deallocate(BR,Bphi)
        !       !    		if (allocated(B_z)) deallocate(B_z)
        !       !    		if (allocated(Btheta)) deallocate(Btheta)
        !       deallocate(Bmag, gammab,chib)
        !       lmagnetized = .false.
        !    endif
        ! endif
    
    
        lcalc_ne = .false.
        icell_loop : do icell=1,n_cells
           !check that in filled cells there is electron density otherwise we need to compute it
           !from scratch.
           if (icompute_atomRT(icell) > 0) then
    
              if (ne(icell) <= 0.0_dp) then
                 write(*,*) "  ** No electron density found in the model! ** "
                 lcalc_ne = .true.
                 exit icell_loop
              endif
    
           endif
        enddo icell_loop
        N_fixed_ne = size(pack(icompute_atomRT,mask=(icompute_atomRT==2)))
        if (N_fixed_ne > 0) then
           write(*,'("Found "(1I5)" cells with fixed electron density values! ("(1I3)" %)")') &
                N_fixed_ne, nint(real(N_fixed_ne) / real(n_cells) * 100)
        endif
    
        write(*,*) "Maximum/minimum velocities in the model (km/s):"
        write(*,*) " V1 = ", 1e-3 * maxval(abs(vfield_x)), 1d-3*minval(abs(vfield_x),mask=icompute_atomRT>0)
        write(*,*) " V2 = ",  1d-3 * maxval(abs(vfield_y)), 1d-3*minval(abs(vfield_y),mask=icompute_atomRT>0)
        write(*,*) "Vphi = ",  1d-3 * maxval(abs(vfield_z)), 1d-3*minval(abs(vfield_z),mask=icompute_atomRT>0)
    
    
        write(*,*) "Typical line extent due to V fields (km/s):"
        write(*,*) v_char/1d3
    
        write(*,*) "Maximum/minimum turbulent velocity (km/s):"
        write(*,*) maxval(vturb)/1d3, minval(vturb, mask=icompute_atomRT>0)/1d3
    
        write(*,*) "Maximum/minimum Temperature in the model (K):"
        write(*,*) real(maxval(T)), real(maxval(T,mask=icompute_atomRT>0))
        write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
        write(*,*) real(maxval(nHtot)), real(minval(nHtot,mask=icompute_atomRT>0))
        if (.not.lcalc_ne) then
           write(*,*) "Maximum/minimum ne density in the model (m^-3):"
           write(*,*) real(maxval(ne)), real(minval(ne,mask=icompute_atomRT>0))
        endif
    
        return
      end subroutine read_spheregrid_ascii

      subroutine read_marcs()

        return
      end subroutine read_marcs
    

end module mhd2mcfost