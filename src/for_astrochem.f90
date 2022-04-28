

  n_cells ! size(data_grid)

  call compute_column(1, column) ! particle/cm^-2

  G(:) = compute_UV_field() ! Habing

  ! unit conversions for astrochem :
  gas_density(:) = densite_gaz(1:n_cells) * masse_mol_gaz / m3_to_cm3 ! nH2/m**3 --> g/cm**3

  do icell=1,n_cells
     dust_mass_density(icell) = sum(densite_pouss(:,icell) * M_grain(:)) ! M_grain en g
     dust_particle_dens(icell,:) = densite_pouss(:,icell) * m3_to_cm3    ! particle/cm^-3
  enddo ! icell

  ! No conversion
  Tdust(1:n_cells)  ! K

  r_grain ! in micron

  r_grid ! in au
