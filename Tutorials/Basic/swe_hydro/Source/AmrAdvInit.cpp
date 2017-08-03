
#include "AmrAdv.H"
#include "AmrAdv_F.H"

using namespace amrex;

void
AmrAdv::InitData ()
{
    if (restart_chkfile.empty())
    {
    	const Real time = 0.0;
    	InitFromScratch(time);
    	AverageDown();

    	if (plot_int > 0) {
    	    WritePlotFile();
    	}
    }
    else
    {
	       InitFromCheckpoint();
    }

    for (int lev = 0; lev < maxLevel(); lev++) {

        //phi_old[lev]->FillBoundary(geom[lev].periodicity());
        //fill_physbc(*phi_old[lev], geom[lev]);

        //phi_new[lev]->FillBoundary(geom[lev].periodicity());
        //fill_physbc(*phi_new[lev], geom[lev]);
    }
}

void AmrAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
				      const DistributionMapping& dm)
{
    int ncomp;
    const int nghost = 6;//phi_new[lev]->nGrow();

    if (lev > max_swe_level) {
        ncomp = 5;
    } else {
        ncomp = 3;
    }

    phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
    phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	       flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }

    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();
    Real cur_time = t_new[lev];

    MultiFab& state = *phi_new[lev];

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        const int* lo  = box.loVect();
        const int* hi  = box.hiVect();

    	initdata(lev, cur_time, ARLIM_3D(lo), ARLIM_3D(hi),
             BL_TO_FORTRAN_3D(state[mfi]),
             ZFILL(dx),
    		 ZFILL(prob_lo), &ncomp, &alpha0, &M, &R,
             p, &nlayers, &gamma, rho);
    }
    //state.FillBoundary(geom[lev].periodicity());
    //fill_physbc(state, geom[lev]);

    //phi_old[lev]->FillBoundary(geom[lev].periodicity());
    //fill_physbc(*phi_old[lev], geom[lev]);
}
