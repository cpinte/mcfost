

module escape

    use parametres
    use grid
    use spectrum_type, only : Istar_cont, lambda_cont
    use utils, only : progress_bar, Bpnu
    use molecular_emission, only : v_proj
    use naleat, only  : seed, stream, gtype
    use mcfost_env, only  : time_begin, time_end, time_tick, time_max

   !$ use omp_lib
#include "sprng_f.h"

    subroutine mean_velocity_gradient()


        return
    end subroutine mean_velocity_gradient

    subroutine escape_prob
    !Probabilistic solution of the non-LTE transfer.

    end subroutine escape_prob

    subroutine averaged_sobolev()
    ! A Sobolev method with averaged quantities
    ! dv/ds -> <dv/ds>

        return
    end subroutine averaged_sobolev

    subroutine sobolev()
    !Large velocity gradient / supersonic approximation

        return
    end subroutine sobolev


end module escape