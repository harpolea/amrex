
#include <AMReX.H>
#include <AMReX_Utility.H>
//#include <stdlib.h>

#include "test_f.H"

using namespace amrex;

int main(int argc, char *argv[]) {
    amrex::Initialize(argc,argv);

    BL_PROFILE_VAR("main()", pmain);

    bool passed;

    test_cons_to_prim(&passed);

    BL_PROFILE_VAR_STOP(pmain);
    amrex::Finalize();

    if (!passed) {
        exit(EXIT_FAILURE);
    }
    return 0;
}
