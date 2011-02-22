subroutine appel_syst(commande, status)

  implicit none

  character(len=512), intent(in) :: commande
  integer, intent(out) :: status
  integer :: system

  status = system(commande)
      
  return

end subroutine appel_syst
