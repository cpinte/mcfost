subroutine appel_syst(commande, status)

  implicit none

  character(len=128), intent(in) :: commande
  integer, intent(out) :: status

  call system(commande, status)
      
  return

end subroutine appel_syst
