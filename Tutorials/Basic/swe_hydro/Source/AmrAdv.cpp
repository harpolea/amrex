
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <AmrAdv.H>
#include <AmrAdvBC.H>
#include <AmrAdvMesh.H>
#include <AmrAdv_F.H>

using namespace amrex;

AmrAdv::AmrAdv ()
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    int nlevs_max = maxLevel() + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= maxLevel(); ++lev) {
	       nsubsteps[lev] = MaxRefRatio(lev-1);
    }

    t_new.resize(nlevs_max, 0.0);
    t_old.resize(nlevs_max, -1.e100);
    dt.resize(nlevs_max, 1.e100);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);

    flux_reg.resize(nlevs_max+1);
}

AmrAdv::~AmrAdv ()
{
    delete[] rho;
    delete[] p;
}

void
AmrAdv::ReadParameters ()
{
    {
    	ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
    	pp.query("max_step", max_step);
    	pp.query("stop_time", stop_time);
    }

    {
    	ParmParse pp("amr"); // Traditionally, these have prefix, amr.

    	pp.query("regrid_int", regrid_int);

    	pp.query("check_file", check_file);
    	pp.query("check_int", check_int);

    	pp.query("plot_file", plot_file);
    	pp.query("plot_int", plot_int);

    	pp.query("restart", restart_chkfile);
    }

    {
    	ParmParse pp("adv");

    	pp.query("cfl", cfl);

    	pp.query("do_reflux", do_reflux);

        pp.query("max_swe_level", max_swe_level);
    }

    {
    	ParmParse pp("swe");

        pp.query("max_swe_level", max_swe_level);
        pp.query("nlayers", nlayers);

        Array<Real> rho_arr, p_arr;

        int n = pp.countval("rho");
        if (n > 0) {
            pp.getarr("rho", rho_arr, 0, n);
        }
        n = pp.countval("p");
        if (n > 0) {
            pp.getarr("p", p_arr, 0, n);
        }

        rho = new Real[nlayers];
        p = new Real[nlayers];

        for (int i = 0; i < nlayers; i++) {
            rho[i] = rho_arr[i];
            p[i] = p_arr[i];
        }
    }

}

void
AmrAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				const DistributionMapping& dm)
{
    int ncomp = phi_new[lev-1]->nComp();
    const int nghost = phi_new[lev-1]->nGrow();

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

    FillCoarsePatch(lev, time, *phi_new[lev], 0, ncomp);

    //if (lev == max_swe_level+1) {
        // do conversion from swe to comp
        //comp_from_swe(phi_new[lev], phi_new[lev-1], p, rho, lo, hi, 4, 3, gamma, gamma_up, dz);
    //}
}


void
AmrAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
		     const DistributionMapping& dm)
{
    const int ncomp = phi_new[lev]->nComp();
    const int nghost = phi_new[lev]->nGrow();

#if __cplusplus >= 201402L
    auto new_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
    auto old_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
#else
    std::unique_ptr<MultiFab> new_state(new MultiFab(ba, dm, ncomp, nghost));
    std::unique_ptr<MultiFab> old_state(new MultiFab(ba, dm, ncomp, nghost));
#endif

    FillPatch(lev, time, *new_state, 0, ncomp);

    std::swap(new_state, phi_new[lev]);
    std::swap(old_state, phi_old[lev]);

    t_new[lev] = time;
    t_old[lev] = time - 1.e200;

    if (lev > 0 && do_reflux) {
	       flux_reg[lev].reset(new FluxRegister(ba, dm, refRatio(lev-1), lev, ncomp));
    }
}

void
AmrAdv::ClearLevel (int lev)
{
    phi_new[lev].reset(nullptr);
    phi_old[lev].reset(nullptr);
    flux_reg[lev].reset(nullptr);
}

void
AmrAdv::InitFromScratch (Real time)
{
    MakeNewGrids(time);
}

void
AmrAdv::AverageDown ()
{

    for (int lev = finest_level-1; lev >= 0; --lev)
    {
        // first check to see if we're in the danger zone
        if (lev == max_swe_level) {
            // convert first
            int n_swe_comp = phi_new[lev]->nComp();
            int nghost = phi_new[lev]->nGrow();
            BoxArray ba = grids[lev+1];
            std::unique_ptr<MultiFab> phi_swe(new MultiFab(ba, dmap[lev+1], n_swe_comp, nghost));

            // TODO: this doesn't work as haven't defined gamma_up_mf

            swe_from_comp_wrapper(lev+1, *phi_swe, *phi_new[lev+1]);

            amrex::average_down(*phi_swe, *phi_new[lev],
 			     geom[lev+1], geom[lev],
 			     0, n_swe_comp, refRatio(lev));

        } else {
            amrex::average_down(*phi_new[lev+1], *phi_new[lev],
 			     geom[lev+1], geom[lev],
 			     0, phi_new[lev]->nComp(), refRatio(lev));
        }
    }
}

void
AmrAdv::AverageDownTo (int crse_lev)
{

    if (crse_lev == max_swe_level) {
        // convert first
        int n_swe_comp = phi_new[crse_lev]->nComp();
        int nghost = phi_new[crse_lev]->nGrow();
        BoxArray ba = grids[crse_lev+1];
        std::unique_ptr<MultiFab> phi_swe(new MultiFab(ba, dmap[crse_lev+1], n_swe_comp, nghost));

        swe_from_comp_wrapper(crse_lev+1, *phi_swe, *phi_new[crse_lev+1]);

        amrex::average_down(*phi_swe, *phi_new[crse_lev],
             geom[crse_lev+1], geom[crse_lev],
             0, n_swe_comp, refRatio(crse_lev));
    } else {
        amrex::average_down(*phi_new[crse_lev+1], *phi_new[crse_lev],
    			 geom[crse_lev+1], geom[crse_lev],
    			 0, phi_new[crse_lev]->nComp(), refRatio(crse_lev));
    }

}

long
AmrAdv::CountCells (int lev)
{
    const int N = grids[lev].size();

    long cnt = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:cnt)
#endif
    for (int i = 0; i < N; ++i)
    {
        cnt += grids[lev][i].numPts();
    }

    return cnt;
}

void
AmrAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
    	Array<MultiFab*> smf;
    	Array<Real> stime;
    	GetData(0, time, smf, stime);

    	AmrAdvPhysBC physbc;
    	amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
    				     geom[lev], physbc);
    }
    else
    {
    	Array<MultiFab*> cmf, fmf;
    	Array<Real> ctime, ftime;
    	GetData(lev-1, time, cmf, ctime);
    	GetData(lev  , time, fmf, ftime);

    	AmrAdvPhysBC cphysbc, fphysbc;
    	Interpolater* mapper = &cell_cons_interp;

    	int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaries
    	int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    	Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

        if (lev-1 == max_swe_level) {
            // do comp from swe conversion
            BoxArray ba = grids[lev-1];
            Array<MultiFab*> comp_mf;

            for (int i = 0; i< cmf.size(); i++) {
                comp_mf.push_back(new MultiFab(ba, dmap[lev-1], 5, phi_new[lev]->nGrow()));
                comp_from_swe_wrapper(lev-1, *cmf[i], *comp_mf[i]);
            }

        	amrex::FillPatchTwoLevels(mf, time, comp_mf, ctime,
                           fmf, ftime,
        				   0, icomp, ncomp, geom[lev-1], geom[lev],
        				   cphysbc, fphysbc, refRatio(lev-1),
        				   mapper, bcs);
        } else {
            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
        				   0, icomp, ncomp, geom[lev-1], geom[lev],
        				   cphysbc, fphysbc, refRatio(lev-1),
        				   mapper, bcs);
        }
    }
}

void
AmrAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Array<MultiFab*> cmf;
    Array<Real> ctime;
    GetData(lev-1, time, cmf, ctime);

    if (cmf.size() != 1) {
	       amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    AmrAdvPhysBC cphysbc, fphysbc;
    Interpolater* mapper = &cell_cons_interp;

    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaryies
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

    // NOTE: 0 index here as getdata using the vector push_back operation, so it puts the data from lev-1 into index 0 element of multifab array.
    amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				 cphysbc, fphysbc, refRatio(lev-1),
				 mapper, bcs);

    if (lev-1 == max_swe_level) {
        Array<MultiFab*> swe_mf;
        GetData(lev, time, swe_mf, ctime);
        // do comp from swe conversion
        comp_from_swe_wrapper(lev, *swe_mf[0], mf);
    }
}

void
AmrAdv::GetData (int lev, Real time, Array<MultiFab*>& data, Array<Real>& datatime)
{
    data.clear();
    datatime.clear();

    const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

    if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
    	data.push_back(phi_new[lev].get());
    	datatime.push_back(t_new[lev]);
    }
    else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
    	data.push_back(phi_old[lev].get());
    	datatime.push_back(t_old[lev]);
    }
    else
    {
    	data.push_back(phi_old[lev].get());
    	data.push_back(phi_new[lev].get());
    	datatime.push_back(t_old[lev]);
    	datatime.push_back(t_new[lev]);
    }
}

void AmrAdv::comp_from_swe_wrapper(int lev, MultiFab& swe_mf, MultiFab& comp_mf) {

    const int n_cons_comp = 5;
    const int n_swe_comp = 3;
    const Real* dx = geom[lev].CellSize();


    for (MFIter mfi(swe_mf, true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const FArrayBox& U_swe      =   swe_mf[mfi];
        FArrayBox& U_comp      =   comp_mf[mfi];

        comp_from_swe(BL_TO_FORTRAN_3D(U_comp),
            BL_TO_FORTRAN_3D(U_swe),
            p, rho,
            bx.loVect(), bx.hiVect(),
            &n_cons_comp, &n_swe_comp,
            &gamma, dx, &alpha0, &M, &R);
    }
}

void AmrAdv::swe_from_comp_wrapper(int lev, MultiFab& swe_mf, MultiFab& comp_mf) {

    const int n_cons_comp = 5;
    const int n_swe_comp = 3;
    const Real* dx = geom[lev].CellSize();

    for (MFIter mfi(swe_mf, true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const FArrayBox& U_swe      =   swe_mf[mfi];
        FArrayBox& U_comp      =   comp_mf[mfi];
        FArrayBox& U_prim      =   comp_mf[mfi];
        FArrayBox& p_comp      =   comp_mf[mfi];

        // do prim conversion here first
        cons_to_prim(BL_TO_FORTRAN_3D(U_comp),
            BL_TO_FORTRAN_3D(U_prim),
            BL_TO_FORTRAN_3D(p_comp),
            bx.loVect(), bx.hiVect(),
            &n_cons_comp, &gamma, &alpha0, &M, &R, dx);

        // then calculate swe
        swe_from_comp(BL_TO_FORTRAN_3D(U_prim),
            BL_TO_FORTRAN_3D(U_swe),
            BL_TO_FORTRAN_3D(p_comp),
            p,
            bx.loVect(), bx.hiVect(),
            &n_cons_comp, &n_swe_comp,
            &alpha0, &M, &R,
            dx);
    }
}
