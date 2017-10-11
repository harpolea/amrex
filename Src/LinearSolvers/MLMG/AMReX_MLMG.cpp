
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>

namespace amrex {

MLMG::MLMG (MLLinOp& a_lp)
    : linop(a_lp),
      namrlevs(a_lp.NAMRLevels()),
      finest_amr_lev(a_lp.NAMRLevels()-1)
{}

MLMG::~MLMG ()
{}

void
MLMG::solve (const Vector<MultiFab*>& a_sol, const Vector<MultiFab const*>& a_rhs,
             Real a_tol_real, Real a_tol_abs)
{
    AMREX_ASSERT(a_sol[0]->nGrow() > 0);
    AMREX_ASSERT(namrlevs <= a_sol.size());
    AMREX_ASSERT(namrlevs <= a_rhs.size());

    linop.prepareForSolve();

    sol = a_sol;

    rhs.resize(namrlevs);
    for (int alev = 0; alev < namrlevs; ++alev)
    {
        rhs[alev].define(a_rhs[alev]->boxArray(), a_rhs[alev]->DistributionMap(), 1, 0);
        MultiFab::Copy(rhs[alev], *a_rhs[alev], 0, 0, 1, 0);        
    }

    const auto& amrrr = linop.AMRRefRatio();

    for (int falev = finest_amr_lev; falev > 0; --falev)
    {
        amrex::average_down(*sol[falev], *sol[falev-1], 0, 1, amrrr[falev-1]);
        amrex::average_down( rhs[falev],  rhs[falev-1], 0, 1, amrrr[falev-1]);
    }

    const int nc = 1;
    int ng = 0;
    linop.make(res, nc, ng);
    linop.make(rescor, nc, ng);
    ng = 1;
    linop.make(cor, nc, ng);

    cor_hold.resize(namrlevs-1);
    for (int alev = 1; alev < finest_amr_lev; ++alev)
    {
        cor_hold[alev].define(cor[alev][0].boxArray(),
                              cor[alev][0].DistributionMap(),
                              cor[alev][0].nComp(),
                              cor[alev][0].nGrow());
    }

    //
    // TODO: We need to fill the fine amr level ghost cells by interploating from the coarse
    // Example: ml_prolongation.f90
    // We can use the fillpatch meta data.
    //

    //
    // TODO: enforce solvability if appropriate
    //
    
    for (int alev = finest_amr_lev; alev >= 0; --alev)
    {
        //
        // TODO: compute residues res given sol and rhs
        // Example: compute_defect.f90
        // 
        if (alev < finest_amr_lev)
        {
            //
            // TODO: (1) Fix crse level residue at crse/fine boundary
            //       (2) retrict the fine res down to this crse level
            // Example: cc_ml_resid.f90
            //
        }
    }

    //
    // TODO: compute the intial inf-norm of res and rhs
    // Example: ml_norm.f90
    //

    //
    // TODO: need a multi-levle covergence test function
    // Example: ml_cc.f90
    // 
    if (false) // replace with the covergence test
    {
        if (verbose >= 1) {
            amrex::Print() << "MLMG: No iterations needed\n";
        }
    }
    else
    {
//        for (int iter = 0; iter < max_iters; ++iter)
        for (int iter = 0; iter < 10; ++iter)
        {
            oneIter();

            // test convergence
        }
    }
}

void
MLMG::oneIter ()
{
    // if converged?
    //    return

    computeResidual(finest_amr_lev, 0);

    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        miniCycle(alev);

        MultiFab::Add(*sol[alev], cor[alev][0], 0, 0, 1, 0);

        computeResWithCrseSolFineCor(alev-1,alev);

        if (alev != finest_amr_lev) {
            MultiFab::Copy(cor_hold[alev], cor[alev][0], 0, 0, 1, 1); // save it for the up cycle
        }
    }

    // coarest amr level
    {    
        // enforce solvability if appropriate

        mgCycle ();

        MultiFab::Add(*sol[0], cor[0][0], 0, 0, 1, 0);
    }

    for (int alev = 1; alev <= finest_amr_lev; ++alev)
    {
        interpCorrection(alev);

        MultiFab::Add(*sol[alev], cor[alev][0], 0, 0, 1, 0);

        if (alev != finest_amr_lev) {
            MultiFab::Add(cor_hold[alev], cor[alev][0], 0, 0, 1, 1);
        }

        computeResWithCrseCorFineCor(alev-1, alev);

        miniCycle (alev);

        MultiFab::Add(*sol[alev], cor[alev][0], 0, 0, 1, 0);

        if (alev != finest_amr_lev) {
            MultiFab::Add(cor[alev][0], cor_hold[alev], 0, 0, 1, 1);            
        }
    }

    for (int alev = finest_amr_lev; alev > 0; --alev)
    {
        // restrict sol[alev] to sol[alev-1]
    }

    // ...
}

void
MLMG::computeResidual (int alev, int mlev)
{
    MultiFab& x = *sol[alev];
    const MultiFab& b = rhs[alev];
    MultiFab& r = res[alev][mlev];

    if (alev > 0) {
        linop.updateBC(alev, *sol[alev-1]);
    }
    linop.residual(alev, mlev, r, x, b, MLLinOp::BCMode::Inhomogeneous);
    VisMF::Write(r, "res");
    amrex::Abort("xxxxx");
}

void
MLMG::computeResWithCrseSolFineCor (int calev, int falev)
{
    MultiFab& crse_sol = *sol[calev];
    const MultiFab& crse_rhs = rhs[calev];
    MultiFab& crse_res = res[calev][0];
    MultiFab& fine_cor = cor[falev][0];
    MultiFab& fine_res = res[falev][0];
    MultiFab& fine_rescor = rescor[falev][0];

    // compute crse and fine, and then construct crse/fine residual
}

void
MLMG::computeResWithCrseCorFineCor (int calev, int falev)
{
    const MultiFab& crse_cor = cor[calev][0];
    MultiFab& fine_cor = cor[falev][0];
    MultiFab& fine_res = res[falev][0];
    MultiFab& fine_rescor = rescor[falev][0];

    // interp correction to supply boundary conditions
    // compute residual
}

void
MLMG::miniCycle (int alev)
{
    Vector<MultiFab>& xs = cor[alev];
    Vector<MultiFab>& bs = res[alev];
    Vector<MultiFab>& rs = rescor[alev];

    for (auto& x : xs) x.setVal(0.0);

    for (int i = 0; i < nu1; ++i) {
        int mglev = 0;
        // linop.smooth(alev, mglev, xs[mglev], bs[mglev]);
    }

    // for ref ratio of 4 ...
    
}

void
MLMG::mgCycle ()
{
    Vector<MultiFab>& xs = cor[0];
    Vector<MultiFab>& bs = res[0];
    Vector<MultiFab>& rs = rescor[0];

    for (auto& x : xs) x.setVal(0.0);

    for (int i = 0; i < nu1; ++i) {
        int mglev = 0;
        // linop.smooth(alev, mglev, xs[mglev], bs[mglev]);
    }

    // compute defect 
    
}

void
MLMG::interpCorrection (int alev)
{
    const MultiFab& crse_cor = cor[alev-1][0];
    MultiFab& fine_cor = cor[alev][0];
    fine_cor.setVal(0.0);
    // prolongation 
}

}

