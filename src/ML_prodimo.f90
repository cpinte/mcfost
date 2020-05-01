module ML_ProDiMo

  use mcfost_env
  use wavelengths
  use prodimo, only : n_phot_envoyes_ISM ! todo : move out of prodimo module
  use constantes

  implicit none

  public :: init_ML, save_J_ML, xgb_compute_features, xgb_predict_Tgas, xgb_predict_abundance

  private

  integer, parameter :: n_lambda_ML = 39
  integer :: n_features

  real(kind=sp), dimension(:,:), allocatable, save :: J_ML
  real(kind=sp), dimension(:,:), allocatable, save :: feature_Tgas
  real(kind=sp), dimension(:,:), allocatable, save :: feature_abundance

  interface
     subroutine xgb_predict(model_name, feature, nrow, nfea, output) bind(C, name='predict')

       use, intrinsic :: iso_c_binding

       character(len=1,kind=c_char), dimension(*) :: model_name
       real(c_float), dimension(nrow, nfea), intent(in) :: feature
       integer(c_int), intent(in), value :: nrow, nfea

       real(c_float), dimension(nrow), intent(out) :: output

     end subroutine xgb_predict
  end interface

contains

  pure function str_f2c (f_string) result (c_string)
    ! Converts a fortran string to C string
    ! From http://fortranwiki.org/fortran/show/Generating+C+Interfaces
    use, intrinsic :: iso_c_binding, only: c_char, c_null_char
    implicit none

    character(len=*), intent(in) :: f_string
    character(len=1,kind=c_char), dimension(len_trim(f_string)+1) :: c_string
    integer :: n, i

    n = len_trim(f_string)
    do i = 1, n
       c_string(i) = f_string(i:i)
    enddo
    c_string(n + 1) = c_null_char

  end function str_f2c

!-----------------------------------------------------------------------------------------

  subroutine init_ML()

    use messages

    integer :: alloc_status

    n_features = 45 ! When no spatial information is used
    write(*,*) "Initializing ML with", n_features, "features"

    ! This should never happen because we fix the wavelenght grid, but we leave it just in case
    if (n_lambda2 /= n_lambda_ML) then
       write(*,*) n_lambda2, "wavelenght bins !!!"
       call error("Incorrect number of wavelength bins for xgboost")
    endif

    alloc_status = 0
    allocate(J_ML(n_lambda2,n_cells), stat=alloc_status)
    if (alloc_status /= 0) call error("Allocation J_ML")

    allocate(feature_Tgas(n_features, n_cells), stat=alloc_status)
    if (alloc_status /= 0) call error("Allocation feature_Tgas")

    allocate(feature_abundance(n_features+1, n_cells), stat=alloc_status)
    if (alloc_status /= 0) call error("Allocation feature_abundance")

    ! Todo : this is ugly, I re-use a variable from io_prodimo.f90
    allocate(n_phot_envoyes_ISM(n_lambda2,nb_proc),  stat=alloc_status)
    if (alloc_status /= 0) call error('Allocation error n_phot_envoyes_ISM')
    n_phot_envoyes_ISM = 0.0

    return

  end subroutine init_ML

!-----------------------------------------------------------------------------------------

  subroutine save_J_ML(lambda, lISM)
    ! sauvegarde le champ de radiation pour ProDiMo
    ! avant et apres le calcul du champ ISM

    use dust_ray_tracing, only : n_phot_envoyes
    use radiation_field, only : xJ_abs, xN_abs
    use prodimo, only : chi_ISM, R_ISM, T_ISM_stars, Wdil, Tcmb
    use parametres, only : Rmax
    use temperature, only : E_disk
    use cylindrical_grid, only : volume
    use wavelengths, only : tab_lambda
    use stars, only : E_stars
    use utils, only : Blambda

    integer, intent(in) :: lambda
    logical, intent(in) :: lISM

    integer :: icell
    real(kind=dp) :: n_photons_envoyes, energie_photon, facteur
    real :: wl

    ! Note: this is a slow loop as we are swapping dimensions
    if (.not.lISM) then
       ! Step2
       n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
       energie_photon = hp * c_light**2 / 2. * (E_stars(lambda) + E_disk(lambda)) / n_photons_envoyes &
            * tab_lambda(lambda) * 1.0e-6  !lambda.F_lambda

       do icell=1, n_cells
          facteur = energie_photon / volume(icell)
          J_ML(lambda,icell) = facteur * sum(xJ_abs(icell,lambda,:))
       enddo

       ! reset for ISM radiation
       xJ_abs(:,lambda,:) = 0.0
    else ! Champs ISM
       n_photons_envoyes = sum(n_phot_envoyes_ISM(lambda,:))
       if (n_photons_envoyes > 0.) then ! We test if there were ISM packets
          wl = tab_lambda(lambda) * 1e-6
          energie_photon = (chi_ISM * 1.71 * Wdil * Blambda(wl,T_ISM_stars) + Blambda(wl,TCmb)) * wl & !lambda.F_lambda
               * (4.*pi*R_ISM**2) / n_photons_envoyes / pi

          do icell=1, n_cells
             facteur = energie_photon / volume(icell)
             J_ML(lambda,icell) =  J_ML(lambda,icell) +  facteur * sum(xJ_abs(icell,lambda,:))
          enddo
       endif
    endif

    return

  end subroutine save_J_ML

!-----------------------------------------------------------------------------------------

  subroutine xgb_compute_features()

    use grains, only : grain, r_grain, n_grains_tot
    use optical_depth, only : compute_column
    use density, only : densite_pouss
    use cylindrical_grid, only : r_grid, z_grid
    use molecular_emission, only : densite_gaz, Tcin
    use temperature, only : Tdust
    use Voronoi_grid, only : Voronoi

    real, dimension(0:3,n_cells) :: N_grains
    logical, dimension(n_grains_tot) :: mask_not_PAH

    integer, parameter :: n_directions = 4
    real, dimension(n_cells,n_directions) :: CD

    real(kind=dp) :: N
    integer :: icell, i

    !--- Moments de la distribution de grain
    mask_not_PAH(:) = .not.grain(:)%is_PAH
    do icell=1, n_cells
       N = sum(densite_pouss(:,icell),mask=mask_not_PAH)
       N_grains(0,icell) = N
       if (N > 0) then
          N_grains(1,icell) = sum(densite_pouss(:,icell) * r_grain(:),mask=mask_not_PAH) / N
          N_grains(2,icell) = sum(densite_pouss(:,icell) * r_grain(:)**2,mask=mask_not_PAH) / N
          N_grains(3,icell) = sum(densite_pouss(:,icell) * r_grain(:)**3,mask=mask_not_PAH) / N
       else
          N_grains(1,icell) = 0.0
          N_grains(2,icell) = 0.0
          N_grains(3,icell) = 0.0
       endif
    enddo
    ! part.cm^-3 --> part.m^-3
    N_grains(0,:) = N_grains(0,:) /  (cm_to_m**3)

    if (n_features == 51) then
       !--- Column density
       write(*,'(a31, $)') "Computing column densities ..."
       call compute_column(1, CD)
       write(*,*) " Done"

       if (lVoronoi) then
          feature_Tgas(1,:) = sqrt(Voronoi(:)%xyz(1)**2 + Voronoi(:)%xyz(2)**2)
          feature_Tgas(2,:) = abs(Voronoi(:)%xyz(3))
       else
          feature_Tgas(1,:) = r_grid(:)
          feature_Tgas(2,:) = z_grid(:)
       endif
       feature_Tgas(3,:) = Tdust(:)
       feature_Tgas(4,:) = densite_gaz(:) * masse_mol_gaz / m3_to_cm3 ! g.cm^3
       feature_Tgas(5:43,:) = J_ML(:,:)
       feature_Tgas(44:47,:) = N_grains(:,:)
       do i=1,n_directions
          feature_Tgas(48+i-1,:) = CD(:,i)  ! CD(n_cells, n_directions)
       enddo
    else if (n_features == 45) then
       feature_Tgas(1,:) = Tdust(:)
       feature_Tgas(2,:) = densite_gaz(:) * masse_mol_gaz / m3_to_cm3 ! g.cm^3
       feature_Tgas(3:41,:) = J_ML(:,:)
       feature_Tgas(42:45,:) = N_grains(:,:)
    endif

    feature_Tgas = log10(max(feature_Tgas,tiny_real))

    return

  end subroutine xgb_compute_features

!-----------------------------------------------------------------------------------------

  subroutine xgb_predict_Tgas()

    use fits_utils, only : cfitsWrite
    use molecular_emission, only : Tcin
    use utils, only : in_dir
    use messages

    character(len=512) :: filename, dir
    integer :: ios

    call xgb_compute_features()

    ! Reading xgboosr file
    filename = "Tgas.xgb"
    dir = in_dir(filename, ML_dir,  status=ios)
    if (ios /=0) then
       call error("xgboost file cannot be found: "//trim(filename))
    else
       filename = trim(dir)//trim(filename)
       write(*,*) "Reading "//trim(filename)
    endif

    ! Predict Tgas
    call xgb_predict(str_f2c(filename), feature_Tgas, n_cells, n_features, Tcin) ! A terme remplacer par un Path

    ! Prepare the features for the abundance prediction
    feature_abundance(1:n_features,:) = feature_Tgas
    feature_abundance(n_features+1,:) = Tcin ! this is still log10(Tcin) here

    ! We predicted log10(Tcin)
    Tcin = 10**Tcin

    filename = trim(data_dir)//"/Tgas_ML.fits"
    if (.not.lVoronoi) then
       call cfitsWrite(trim(filename),Tcin,[n_rad,nz])
    else
       call cfitsWrite(trim(filename),Tcin,[n_cells])
    endif

    return

  end subroutine xgb_predict_Tgas

!-----------------------------------------------------------------------------------------

  subroutine xgb_predict_abundance(imol)

    use fits_utils, only : cfitsWrite
    use molecular_emission, only : tab_abundance, mol
    use utils, only : in_dir
    use messages

    integer, intent(in) :: imol

    character(len=512) :: filename, dir
    integer :: ios
    real :: factor
    logical :: lfactor

    ! TODO : Some molecules are given with in different units, we need to adapt the code (?????)
    factor = 1.0 ; lfactor = .false.
    select case(trim(mol(imol)%name))
    case("CO")
       filename = "12co.xgb"
    case("13C16O")
       filename = "12co.xgb"
       factor = 1.428571428e-2 ; lfactor = .true.
    case("C18O")
       filename = "12co.xgb"
       factor = 2e-3 ; lfactor = .true.
    case("HCO+")
       filename = "hco+.xgb"
    case default
       call error("Selected molecule does not have xgboost training data")
    end select

    dir = in_dir(filename, ML_dir,  status=ios)
    if (ios /=0) then
       call error("xgboost file cannot be found: "//trim(filename))
    else
       filename = trim(dir)//trim(filename)
       write(*,*) "Reading "//trim(filename)
       if (lfactor) write(*,*) "Correcting abundance by factor ", factor
    endif

    call xgb_predict(str_f2c(filename) , feature_abundance, n_cells, n_features+1, tab_abundance)
    if (lfactor) tab_abundance = tab_abundance * factor

    tab_abundance = max(0.0, tab_abundance) ! preventing slightly negative abundances

    filename = trim(data_dir2(imol))//"/abundance_ML.fits"
    if (.not.lVoronoi) then
       call cfitsWrite(trim(filename),tab_abundance,[n_rad,nz])
    else
       call cfitsWrite(trim(filename),tab_abundance,[n_cells])
    endif

    return

  end subroutine xgb_predict_abundance

end module ML_ProDiMo
