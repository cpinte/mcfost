module mcfost_env

  use iso_fortran_env

  implicit none

  real, parameter :: mcfost_version = 4.0
  character(8), parameter :: mcfost_release = "4.0.00"
  real, parameter :: required_utils_version = 4.0

  character(len=128) :: web_server    = "http://ipag.osug.fr/public/pintec/"
  character(len=128) :: webpage       = "/mcfost/"
  character(len=128) :: utils_webpage = "/mcfost_utils/"
  character(len=128) :: doc_webpage   = "https://mcfost.readthedocs.io/"

  ! Système
  integer :: nb_proc
  integer, parameter :: cache_line_size = 64 ! 64 bytes = 16 floats = 8 double, from Core 2 Duo to i7 + Xeon Phi

  integer, parameter :: sp = selected_real_kind(p=6,r=37)
  integer, parameter :: dp = selected_real_kind(p=13,r=200)
  integer, parameter :: lp = logical_kinds(1) ! smallest logical available
  integer, parameter :: limite_stack = 5000000

  real :: max_mem = 4. ! GBytes maximum size for 1 array (mcfost can have 2 arrays of this size)
  logical :: low_mem_scattering, low_mem_th_emission, low_mem_th_emission_nLTE, lMueller_pos_multi

  character(len=512) :: mcfost_utils, my_mcfost_utils, data_dir, root_dir, tmp_dir, basename_data_dir, seed_dir
  character(len=512) :: lambda_filename, band, model_pah, pah_grain, cmd_opt
  character(len=512), dimension(100) :: data_dir2, basename_data_dir2
  character(len=512), dimension(:), allocatable :: search_dir, dust_dir, mol_dir, star_dir, lambda_dir, ML_dir

  integer :: time_begin, time_end, time_tick, time_max

end module mcfost_env
