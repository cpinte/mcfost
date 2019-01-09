MODULE simple_models

 use atmos_type, only : atmos, init_atomic_atmos
 
 ! MCFOST's modules
 use input
 use parametres
 use grid
 use density
 use constantes
 use fits_utils, only : print_error
 
 IMPLICIT NONE

 CONTAINS
  
  SUBROUTINE uniform_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! constant on the grid
   ! Values width idk are from a Solar model atmosphere by 
   ! Fontenla et al. namely the FAL C model.
  ! ----------------------------------------------------------- !
   double precision, dimension(n_cells) :: nHtot, T, ne
   double precision :: Vkepmax
   !idk = 10
!    nHtot =  2.27414200581936d16
!    T = 45420d0
!    ne = 2.523785d16

   !idk = 75
   T=7590d0
   ne = 4.446967d20
   nHtot = 1.259814d23

    !idk = 81   
!     T=9400d0
!     ne = 3.831726d21
!     nHtot = 1.326625d23

    !idk = 0 
!     T = 100000d0
!     ne = 1.251891d16
!     nHtot = 1.045714d16
    
   !!more or less the same role as init_molecular_disk
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   atmos%moving = .false. !force to be static for this case
   !atmos%vturb = 0d0
   !tmos%vturb = 9.506225d3 !m/s !idk=10
   atmos%vturb = 1.696164d3 !idk = 75
   !atmos%vturb = 1.806787d3 !idk=81
   !atmos%vturb = 10.680960d3 !idk=0
   atmos%velocity_law = 0 !keplerian = 0
   atmos%v_char = maxval(atmos%vturb)
   write(*,*) " >-< Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3

  RETURN
  END SUBROUTINE uniform_law_model
  
  SUBROUTINE prop_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! proportional to dust density and temperature.
  ! ----------------------------------------------------------- !
   double precision, dimension(n_cells) :: nHtot, T, ne
   double precision :: Vkepmax
   
   nHtot = 1d14 * densite_gaz/MAXVAL(densite_gaz) + 1d12
   T = 3000d0+Tdust
   ne = 0d0
   
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   if (lstatic) atmos%moving = .false.
   atmos%velocity_law = 0 !keplerian
   atmos%v_char = maxval(atmos%vturb)
   if (atmos%moving) then
    Vkepmax = dsqrt(Ggrav * sum(etoile%M) * Msun_to_kg * Rmin**2 / &
                ((Rmin**2 + minval(z_grid)**2)**1.5 * AU_to_m) )

    atmos%v_char = atmos%v_char + Vkepmax
    write(*,*) " >-< Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3
    if (atmos%velocity_law == 0) then
      if (.not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
      atmos%Vmap = 0d0
    end if
   else
    write(*,*) " >-< Model is static"
   end if

  RETURN
  END SUBROUTINE prop_law_model
  
  SUBROUTINE spherically_symmetric_shells_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple spherically symmetric model of shells
   ! in expansion or in contraction.
   ! building
  ! ----------------------------------------------------------- !
  integer :: izone, i, j, k, icell, n_zones=1 !Only one component firstly
  double precision, dimension(n_cells) :: Vr, Tr,vturbr, nHtotr, ner
  type(disk_zone_type) :: dz ! to define the properties of the model
  double precision :: r, nH0, T0, vt, rcyl, z, r0, Mdot, yr_to_sec = 3.154d7
  double precision :: vinf, v0, beta, b=1d0
  
  ! -- Parameters -- !
  nH0 = 1d12 !m^-3
  T0 = 2362d0
  Mdot = 1d-4 * Msun_to_kg / yr_to_sec ! Msun/yr -> kg/s
  beta = 5d-1
  vt = 5d0
  v0 = 14.86d0 !km/s
  vinf = 40d0 !km/s
  r0 = Rmin !/Rmin
  ! ---------------- !
  
  ! ---- Model------ !
  Tr = 0d0
  nHtotr = 0d0
  Vr = 0d0
  vturbr = 0d0
  ner = 0d0
  ! ---------------- !
  
  main : do izone=1, n_zones
   !dz = disk_zone(izone) ! for now the parameters are hardcoded
    do i=1, n_rad
     do j=j_start,nz
      do k=1, n_az
       if (j==0) then !midplane
        icell = cell_map(i,1,k)
        rcyl = r_grid(icell) !AU
        z = 0.0_dp
       else
        icell = cell_map(i,j,k)
        rcyl = r_grid(icell)
        z = z_grid(icell)/z_scaling_env
       end if
       r = sqrt(rcyl**2 + z**2)
       if (r >= Rmin .and. r<= 10) then
        Tr(icell) = T0 * dsqrt(r0/r)
        Vr(icell) = 1d3 * (v0 + (vinf-v0) * (1d0 - b*r0/r)**beta) !m/s
        nHtotr(icell) = Mdot / (4*PI * (r*AU_to_m)**2 * Vr(icell)) * 1d3/masseH!kg*1d3=g
        vturbr(icell) = vt*1d3 !m/s
        write(*,*) icell, r, nHtotr(icell), Tr(icell), Vr(icell)/1d3 
       end if
      end do !k
     end do !j
    end do !i
  end do main !over n_zones
    
  CALL init_atomic_atmos(n_cells, Tr, ner, nHtotr)
  if (lstatic) atmos%moving = .false.
  atmos%vturb = vturbr
  atmos%v_char = maxval(atmos%vturb)
  if (atmos%moving) then
   if (.not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
   atmos%Vmap = Vr !and set keyword lkeplerian=.false. and linfall=.true.
  				  ! otherwise vinfall/vkep = cte = expected.
   atmos%velocity_law = -1
   atmos%v_char = atmos%v_char + maxval(abs(atmos%Vmap))
   write(*,*) " >-< Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3
  else
   write(*,*) " >-< Model is static"
  end if
  RETURN
  END SUBROUTINE spherically_symmetric_shells_model

  SUBROUTINE spherically_symmetric_star_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple spherically symmetric model of a star
   ! building
  ! ----------------------------------------------------------- !
  use math, only : interp1D, locate
  integer :: izone, i, j, k, icell, n_zones=1, n
  double precision, dimension(n_cells) :: Vr, Tr,vturbr, nHtotr, ner
  double precision, dimension(n_cells) :: Scale_star_disk
  type(disk_zone_type) :: dz ! to define the properties of the model
  integer, parameter :: NDEP = 56 !82!56
  double precision, parameter :: Rstarmarcs=165d0, f = 1d0 !1d0 to no change density
  double precision, dimension(NDEP) :: Tmarcs, nHmarcs, Rmarcs, vtmarcs, nemarcs
  double precision :: r, vt, rcyl, z, v0, alpha, b = 1d0, rin, rmin, rmax
  double precision :: Mdot, yr_to_sec = 3.154d7, Runit, rs, T0, Rlimit, r0 !Vr = v0 at r0

  ! Dwarf 1Rsun
!   vtmarcs = 2d0
!   data Rmarcs / 0.99985071, 0.9998732 , 0.99989189, 0.9999077 , 0.99992117,&
!        0.99993283, 0.99993806, 0.99994307, 0.99994799, 0.99995291,&
!        0.99995809, 0.99996381, 0.99997061, 0.9999789 , 0.99998872,&
!        1.        , 1.00001276, 1.00002696, 1.00004244, 1.00005907,&
!        1.00007659, 1.00009474, 1.00011331, 1.00013214, 1.00015114,&
!        1.00017   , 1.000189  , 1.00020786, 1.00022657, 1.00024529,&
!        1.000264  , 1.00028243, 1.00030086, 1.00031914, 1.00033743,&
!        1.00035557, 1.00037371, 1.00039171, 1.00040957, 1.00042743,&
!        1.00044529, 1.00046286, 1.00048057, 1.00049814, 1.00051571,&
!        1.00053314, 1.00056871, 1.00060414, 1.00063914, 1.000674  ,&
!        1.00070871, 1.000743  , 1.000777  , 1.00081071, 1.00084386,&
!        1.00087671/
!   data Tmarcs / 10170.6,  9908.4,  9653.9,  9398.7,  9137.1,  8859.4,  8706.6,&
!         8540.9,  8356.7,  8146.3,  7900. ,  7606.6,  7269.3,  6953.7,&
!         6700. ,  6473.7,  6275.1,  6099.8,  5945.9,  5811.1,  5693.2,&
!         5590.3,  5499.8,  5420.9,  5351.6,  5289.7,  5235.4,  5186.7,&
!         5141.7,  5102. ,  5064.9,  5030.7,  4998.8,  4968.8,  4940.3,&
!         4912.9,  4886.3,  4860.5,  4835.3,  4810.4,  4785.8,  4761.5,&
!         4737.6,  4713.7,  4690. ,  4666.3,  4620. ,  4573.9,  4528.2,&
!         4483. ,  4438.2,  4393.9,  4350.2,  4307.3,  4266.1,  4232.3 / !K
!   data nHmarcs /1.22249032e+23, 1.18052880e+23, 1.14965225e+23, 1.12769465e+23,&
!        1.11384905e+23, 1.10773320e+23, 1.10803050e+23, 1.11074865e+23,&
!        1.11622744e+23, 1.12510391e+23, 1.13814256e+23, 1.15576810e+23,&
!        1.17526237e+23, 1.18320448e+23, 1.17131255e+23, 1.14561749e+23,&
!        1.10599188e+23, 1.05430446e+23, 9.92763728e+22, 9.24427607e+22,&
!        8.52354022e+22, 7.79515956e+22, 7.08334266e+22, 6.40380385e+22,&
!        5.76673621e+22, 5.17808572e+22, 4.63827708e+22, 4.14756512e+22,&
!        3.70420852e+22, 3.30370533e+22, 2.94427177e+22, 2.62195802e+22,&
!        2.33353627e+22, 2.07582118e+22, 1.84583976e+22, 1.64078892e+22,&
!        1.45812040e+22, 1.29545580e+22, 1.15071403e+22, 1.02194143e+22,&
!        9.07396665e+21, 8.05550786e+21, 7.15002240e+21, 6.34561834e+21,&
!        5.63040374e+21, 4.99503496e+21, 3.92883856e+21, 3.08812164e+21,&
!        2.42557131e+21, 1.90377044e+21, 1.49324406e+21, 1.17063301e+21,&
!        9.17419863e+20, 7.18909588e+20, 5.63592499e+20, 4.40978216e+20 / !m^-3
!   data nemarcs / 1.44230356e+20, 1.37966709e+20, 1.33260966e+20, 1.29789561e+20,&
!        1.27401761e+20, 1.26007604e+20, 1.25731678e+20, 1.25728648e+20,&
!        1.26064498e+20, 1.26795391e+20, 1.28007209e+20, 1.29754839e+20,&
!        1.31760166e+20, 1.32542674e+20, 1.31150963e+20, 1.28228284e+20,&
!        1.23768333e+20, 1.17970903e+20, 1.11077079e+20, 1.03419587e+20,&
!        9.53531133e+19, 8.72006350e+19, 7.92351946e+19, 7.16289026e+19,&
!        6.45022949e+19, 5.79192876e+19, 5.18822076e+19, 4.63925806e+19,&
!        4.14330091e+19, 3.69542408e+19, 3.29334288e+19, 2.93290378e+19,&
!        2.61039662e+19, 2.32208518e+19, 2.06484084e+19, 1.83546164e+19,&
!        1.63126232e+19, 1.44928504e+19, 1.28733289e+19, 1.14329297e+19,&
!        1.01518670e+19, 9.01262129e+18, 7.99968770e+18, 7.09956145e+18,&
!        6.29949007e+18, 5.58876833e+18, 4.39592142e+18, 3.45526982e+18,&
!        2.71405686e+18, 2.13022404e+18, 1.67095425e+18, 1.30988793e+18,&
!        1.02653312e+18, 8.04451630e+17, 6.30659452e+17, 4.93431958e+17 /
   ! AGB 165Rsun
   vtmarcs = 1d0
   data rmarcs / 0.97802597, 0.97863203, 0.97914286, 0.97961905, 0.98014719,&
       0.98107359, 0.98176623, 0.98274459, 0.98408658, 0.98588745,&
       0.98813853, 0.99070996, 0.99334805, 0.99581385, 0.99802597,&
       1.        , 1.00178961, 1.00345455, 1.00504935, 1.00661645,&
       1.00819394, 1.00980952, 1.01148052, 1.01322944, 1.01506494,&
       1.01697835, 1.01898701, 1.02107359, 1.0232381 , 1.02547186,&
       1.02775758, 1.03009524, 1.03248485, 1.03490909, 1.03736797,&
       1.03987013, 1.04240693, 1.04497835, 1.04758442, 1.05022511,&
       1.05290043, 1.05561905, 1.05837229, 1.06116883, 1.06400866,&
       1.06688312, 1.07294372, 1.07916883, 1.0855671 , 1.09212121,&
       1.09878788, 1.10554113, 1.11238095, 1.11904762, 1.12554113,&
       1.13177489 /
   data Tmarcs / 8043.4, 7862.7, 7661.1, 7406.5, 7005. , 6332. , 6010.9, 5718.6,&
       5452.5, 5211.2, 4992.4, 4794.7, 4616.1, 4455.2, 4310. , 4178.5,&
       4058.8, 3949.4, 3849.1, 3757.2, 3673.1, 3596.2, 3526.1, 3462.4,&
       3404.9, 3352.7, 3306. , 3263.6, 3225.9, 3191.7, 3161.4, 3132.9,&
       3106.6, 3082. , 3058.9, 3036.9, 3015.7, 2995.1, 2974.9, 2955. ,&
       2935.4, 2915.8, 2896.4, 2877. , 2857.6, 2839.3, 2799.2, 2759.8,&
       2719.9, 2679.8, 2639.7, 2599.7, 2560.1, 2521.5, 2485.5, 2452. /
   data nHmarcs /1.64344109d21, 1.68138055d21, 1.72654658d21, 1.78737016d21,&
       1.89095091d21, 2.06739951d21, 2.14026737d21, 2.18724003d21,&
       2.19988652d21, 2.16977584d21, 2.09449913d21, 1.98549845d21,&
       1.86565793d21, 1.75545283d21, 1.65970086d21, 1.57840202d21,&
       1.50673859d21, 1.44049509d21, 1.37726265d21, 1.31463243d21,&
       1.25079778d21, 1.18395207d21, 1.11469749d21, 1.04303407d21,&
       9.69564001d20, 8.94889506d20, 8.20215012d20, 7.46744945d20,&
       6.75683732d20, 6.08838015d20, 5.45244252d20, 4.86709084d20,&
       4.32931404d20, 3.83790769d20, 3.39347400d20, 2.99239970d20,&
       2.63227593d20, 2.31009162d20, 2.02283570d20, 1.76749711d20,&
       1.54046256d20, 1.33992540d20, 1.16287459d20, 1.00690125d20,&
       8.69596533d19, 7.49153799d19, 5.52169710d19, 4.02941163d19,&
       2.91170307d19, 2.08546592d19, 1.48084340d19, 1.04604514d19,&
       7.36507312d18, 5.19529729d18, 3.69578526d18, 2.65636448d18 / 
   data nemarcs / 7.70270855e+16, 6.16081894e+16, 4.73464123e+16, 3.32099744e+16,&
       1.79599683e+16, 5.36241700e+15, 2.75093478e+15, 1.41980852e+15,&
       7.50660423e+14, 4.18769725e+14, 2.57225042e+14, 1.80971116e+14,&
       1.44839637e+14, 1.24400362e+14, 1.08895923e+14, 9.47811526e+13,&
       8.14445204e+13, 6.89924481e+13, 5.76558598e+13, 4.75960227e+13,&
       3.88855320e+13, 3.15601491e+13, 2.55323265e+13, 2.06447580e+13,&
       1.67560007e+13, 1.36618985e+13, 1.12193013e+13, 9.26781705e+12,&
       7.72138418e+12, 6.47204225e+12, 5.46873471e+12, 4.63302956e+12,&
       3.93784120e+12, 3.35825354e+12, 2.86979733e+12, 2.45652178e+12,&
       2.10248203e+12, 1.79918365e+12, 1.53871644e+12, 1.31475603e+12,&
       1.12169821e+12, 9.55604841e+11, 8.12465744e+11, 6.89550467e+11,&
       5.83975899e+11, 4.94885049e+11, 3.49829676e+11, 2.45489907e+11,&
       1.70614396e+11, 1.17544121e+11, 8.03946103e+10, 5.47182420e+10,&
       3.72034067e+10, 2.53839830e+10, 1.75485352e+10, 1.23147419e+10 /  
  !FAL C 1Rsun
!   data Rmarcs / 0.99851 , 0.99853 , 0.998551, 0.998571, 0.998592, 0.998612,&
!        0.998633, 0.998653, 0.998674, 0.998694, 0.998715, 0.998735,&
!        0.998756, 0.998786, 0.998817, 0.998868, 0.998919, 0.99897 ,&
!        0.999021, 0.999072, 0.999123, 0.999224, 0.999324, 0.999424,&
!        0.999525, 0.999624, 0.999705, 0.999776, 0.999847, 0.999927,&
!        1.000027, 1.000135, 1.000235, 1.000335, 1.000435, 1.000555,&
!        1.000694, 1.000878, 1.001063, 1.001258, 1.001446, 1.001657,&
!        1.001865, 1.002037, 1.002226, 1.002388, 1.002514, 1.002614,&
!        1.002728, 1.002838, 1.002936, 1.003001, 1.00303 , 1.003045,&
!        1.003049, 1.003051, 1.003053, 1.003054, 1.003055, 1.003056,&
!        1.003056, 1.003057, 1.003058, 1.003058, 1.003059, 1.003059,&
!        1.00306 , 1.00306 , 1.003061, 1.003062, 1.003062, 1.003063,&
!        1.003064, 1.003065, 1.003066, 1.003068, 1.00307 , 1.003071,&
!        1.003075, 1.003078, 1.003081, 1.003083 /
!   data Tmarcs / 9400.,   9140.,   8860.,   8540.,   8220.,   7900.,   7590.,&
!          7280.,   6980.,   6720.,   6520.,   6340.,   6180.,   5980.,&
!          5790.,   5580.,   5410.,   5270.,   5150.,   5060.,   4990.,&
!          4880.,   4780.,   4690.,   4610.,   4540.,   4510.,   4500.,&
!          4520.,   4560.,   4680.,   4900.,   5160.,   5380.,   5570.,&
!          5760.,   5950.,   6180.,   6370.,   6570.,   6740.,   6910.,&
!          7080.,   7220.,   7410.,   7600.,   7780.,   7970.,   8273.,&
!          8635.,   8988.,   9228.,   9358.,   9458.,   9587.,   9735.,&
!          9983.,  10340.,  10850.,  11440.,  12190.,  13080.,  14520.,&
!         16280.,  17930.,  20420.,  24060.,  27970.,  32150.,  36590.,&
!         41180.,  45420.,  49390.,  53280.,  60170.,  66150.,  71340.,&
!         75930.,  83890.,  90820.,  95600., 100000. /
!   data nHmarcs /1.32662462e+23, 1.31275567e+23, 1.30020496e+23, 1.29180598e+23,&
!        1.28402892e+23, 1.27401250e+23, 1.25981369e+23, 1.24607714e+23,&
!        1.22812486e+23, 1.20282003e+23, 1.16619660e+23, 1.12638909e+23,&
!        1.08385618e+23, 1.01434202e+23, 9.45184976e+22, 8.21612134e+22,&
!        7.04214461e+22, 5.97136909e+22, 5.01483483e+22, 4.17101990e+22,&
!        3.43691248e+22, 2.30240577e+22, 1.52160292e+22, 9.93341605e+21,&
!        6.41610951e+21, 4.10450602e+21, 2.84430487e+21, 2.05290447e+21,&
!        1.47210490e+21, 1.00590598e+21, 6.20960707e+20, 3.65754559e+20,&
!        2.29178909e+20, 1.48043185e+20, 9.80419662e+19, 6.15684103e+19,&
!        3.66282003e+19, 1.85672803e+19, 9.77814031e+18, 5.24351028e+18,&
!        3.04612025e+18, 1.73738021e+18, 1.04788518e+18, 7.11571167e+17,&
!        4.76695155e+17, 3.45149141e+17, 2.71927128e+17, 2.26075118e+17,&
!        1.83065106e+17, 1.50485091e+17, 1.28173075e+17, 1.16667062e+17,&
!        1.12308055e+17, 1.09119049e+17, 1.04571040e+17, 1.00765034e+17,&
!        9.65480272e+16, 9.21940221e+16, 8.72100175e+16, 8.24780142e+16,&
!        7.75050116e+16, 7.24630094e+16, 6.57624071e+16, 5.91966052e+16,&
!        5.41514040e+16, 4.80215027e+16, 4.12490015e+16, 3.58761009e+16,&
!        3.14878205e+16, 2.78867302e+16, 2.49529201e+16, 2.27414201e+16,&
!        2.09965070e+16, 1.95244380e+16, 1.73522158e+16, 1.58073849e+16,&
!        1.46666937e+16, 1.37824624e+16, 1.24723003e+16, 1.15182218e+16,&
!        1.09401812e+16, 1.04571358e+16/
!   data nemarcs /3.831726e+21, 2.952398e+21, 2.194603e+21, 1.531806e+21,&
!        1.041290e+21, 6.869667e+20, 4.446967e+20, 2.789487e+20,&
!        1.719226e+20, 1.099800e+20, 7.641340e+19, 5.456858e+19,&
!        4.028244e+19, 2.749041e+19, 1.950519e+19, 1.329338e+19,&
!        9.816085e+18, 7.584832e+18, 5.984860e+18, 4.803401e+18,&
!        3.880555e+18, 2.557346e+18, 1.680047e+18, 1.099721e+18,&
!        7.167248e+17, 4.648930e+17, 3.294358e+17, 2.441017e+17,&
!        1.818187e+17, 1.307890e+17, 9.114873e+16, 7.544027e+16,&
!        8.363294e+16, 9.742137e+16, 1.111659e+17, 1.229867e+17,&
!        1.308237e+17, 1.340907e+17, 1.280117e+17, 1.220873e+17,&
!        1.137706e+17, 1.031897e+17, 9.612962e+16, 9.185673e+16,&
!        8.838377e+16, 8.445100e+16, 8.095357e+16, 7.846423e+16,&
!        7.561758e+16, 7.168647e+16, 6.651616e+16, 6.141942e+16,&
!        5.841143e+16, 5.801540e+16, 6.265075e+16, 6.599262e+16,&
!        6.848673e+16, 6.967801e+16, 6.991762e+16, 6.918466e+16,&
!        6.747161e+16, 6.520094e+16, 6.137226e+16, 5.700496e+16,&
!        5.337548e+16, 4.860836e+16, 4.291177e+16, 3.804151e+16,&
!        3.394113e+16, 3.043595e+16, 2.748732e+16, 2.523785e+16,&
!        2.344633e+16, 2.192829e+16, 1.969550e+16, 1.811689e+16,&
!        1.694766e+16, 1.603707e+16, 1.467464e+16, 1.366348e+16,&
!        1.304293e+16, 1.251891e+16/
!   data vtmarcs / 1.806787 ,  1.785877 ,  1.760494 ,  1.746632 ,  1.733641 ,&
!         1.716951 ,  1.696164 ,  1.677151 ,  1.651916 ,  1.621151 ,&
!         1.579327 ,  1.541506 ,  1.506637 ,  1.448816 ,  1.395714 ,&
!         1.29188  ,  1.195529 ,  1.093189 ,  0.995005 ,  0.8929698,&
!         0.8043706,  0.6274118,  0.5728015,  0.5217747,  0.5555816,&
!         0.6536559,  0.7225354,  0.8043894,  0.8897392,  0.9947371,&
!         1.168133 ,  1.368657 ,  1.533542 ,  1.77223  ,  2.000838 ,&
!         2.254545 ,  2.595622 ,  3.053286 ,  3.552633 ,  4.0603   ,&
!         4.504405 ,  4.926045 ,  5.296849 ,  5.590575 ,  5.932102 ,&
!         6.207403 ,  6.431674 ,  6.614389 ,  6.82316  ,  7.044414 ,&
!         7.268274 ,  7.399466 ,  7.452542 ,  7.492713 ,  7.547351 ,&
!         7.591289 ,  7.642009 ,  7.696761 ,  7.762707 ,  7.82944  ,&
!         7.905294 ,  7.98737  ,  8.10573  ,  8.234054 ,  8.342731 ,&
!         8.494678 ,  8.700418 ,  8.889256 ,  9.065828 ,  9.230186 ,&
!         9.380583 ,  9.506225 ,  9.614243 ,  9.717166 ,  9.899261 ,&
!        10.04315  , 10.15883  , 10.25485  , 10.40904  , 10.53191  ,&
!        10.61132  , 10.68096  /

  ! -- Parameters -- !
  alpha = 3
  v0 = 100d0 !km/s
  T0 = 0d3 ! for envelope
  r0 = Rmarcs(1)!1d0
  ! ---------------- !

  ! ---- Model------ !
  Tr = 0d0
  nHtotr = 0d0
  Vr = 0d0
  vturbr = 0d0
  ner = 0d0
  ! ---------------- !

  main : do izone=1, n_zones
   dz = disk_zone(izone) ! for now the parameters are hardcoded
   rin = dz%rin
   rmin = dz%rmin
   rmax = dz%rmax
   Rlimit = maxval(sqrt(r_grid**2 + z_grid**2))!rmax
   Runit = (rmarcs(NDEP) - rmarcs(1)) / (Rlimit - rmin) !Rstar/AU

    do i=1, n_rad
     do j=j_start,nz
      do k=1, n_az
       if (j==0) then !midplane
        icell = cell_map(i,1,k)
        rcyl = r_grid(icell) !AU
        z = 0.0_dp
       else
        icell = cell_map(i,j,k)
        rcyl = r_grid(icell)
        z = z_grid(icell)/z_scaling_env
       end if
       r = sqrt(rcyl**2 + z**2) !ne va pas jusqu'à Rmax why ?
       rs = (r - Rmin) * Runit + rmarcs(1)
       !Scale_star_disk(icell) = Rsun * Rstarmarcs * Runit
       Scale_star_disk(icell) = Rsun*Rstarmarcs*rs/r !Rstar/AU * Rstar(in m) results in m/AU
       if (r >= Rmin .and. r <= Rlimit) then
        !interp cmass on the new R grid
        Tr(icell) = dexp(interp1D(log(rmarcs),LOG(Tmarcs),log(rs)))
        nHtotr(icell) = dexp(interp1D(log(rmarcs),LOG(f*nHmarcs),log(rs)))
        ner(icell) = dexp(interp1D(log(rmarcs),LOG(f*nemarcs),log(rs)))
        vturbr(icell) = 1d3 * dexp(interp1D(log(rmarcs),LOG(vtmarcs),log(rs))) !m/s 
        Vr(icell) = 1d3 * v0 * (rs/r0)**alpha !m/s, with rs in Rstar
 
        !Tr(icell) = T0 * sqrt(Rmin/r)
        !Vr(icell) = 1d3 * v0 * (1d0-b*Rmin/r)**alpha !m/s
        !nHtotr(icell) = 1d3/masseH * Mdot / (4*pi * (r*AU_to_m)**2 * Vr(icell)) !m^-3
        
        write(*,*) icell, rs, ner(icell), nHtotr(icell), Tr(icell), Vr(icell), vturbr(icell)
       end if   
      end do !k
     end do !j
    end do !i
  end do main !over n_zones

  CALL init_atomic_atmos(n_cells, Tr, ner, nHtotr)
  
  atmos%Scale2disk = Scale_star_disk
  write(*,*) "min/max distance in Rstar per AU in km/AU", 1d-3*minval(atmos%Scale2disk), &
  														1d-3*maxval(atmos%Scale2disk)
  if (lstatic) atmos%moving = .false.
  atmos%vturb = vturbr
  atmos%v_char = maxval(atmos%vturb)
  if (atmos%moving) then
   if (.not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
   atmos%Vmap = Vr !and set keyword lkeplerian=.false. and linfall=.true.
  				  ! otherwise vinfall/vkep = cte = expected.
   atmos%velocity_law = -1
   atmos%v_char = atmos%v_char + maxval(abs(atmos%Vmap))
   write(*,*) " >-< Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3
  else
   write(*,*) " >-< Model is static"
  end if

  RETURN
  END SUBROUTINE spherically_symmetric_star_model

 END MODULE simple_models