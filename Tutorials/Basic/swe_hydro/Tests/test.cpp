
#include <AMReX.H>
#include <AMReX_Utility.H>
//#include <stdlib.h>

#include "test_f.H"

using namespace amrex;

int main(int argc, char *argv[]) {
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    bool passed = true;

    test_cons_to_prim(&passed);
    test_rhoh_from_p(&passed);
    test_p_from_rhoh(&passed);
    test_p_from_rho_eps(&passed);
    test_calc_gamma_up(&passed);
    test_calc_gamma_down(&passed);
    test_gr_sources(&passed);
    test_W_swe(&passed);
    test_swe_from_comp(&passed);
    test_calc_gamma_swe(&passed);
    test_comp_from_swe(&passed);
    test_gr_swe_flux(&passed);
    test_f_of_p(&passed);
    test_zbrent(&passed);
    test_gr_comp_flux(&passed);

    BL_PROFILE_VAR_STOP(pmain);
    amrex::Finalize();

    if (!passed) {
        exit(EXIT_FAILURE);
    }
    return 0;
}
