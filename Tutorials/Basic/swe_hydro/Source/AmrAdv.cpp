
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_BCUtil.H>

#include "AmrAdv.H"
#include "AmrAdvBC.H"
#include "AmrAdv_F.H"

using namespace amrex;

AmrAdv::AmrAdv ()
{
    ReadParameters();

    //Initialize();
    InitAmrCore();

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
    }

    {
    	ParmParse pp("swe");

        pp.query("max_swe_level", max_swe_level);
        pp.query("nlayers", nlayers);

        std::cout << "nlayers are: " << nlayers << '\n';

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

        pp.query("gr", gr);
        pp.query("M", M);
        pp.query("R", R);
        alpha0 = sqrt(1.0 - 2.0 * M / R);
        pp.query("gamma", gamma);
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

    phi_new[lev]->FillBoundary(geom[lev].periodicity());
    fill_physbc(*phi_new[lev], geom[lev]);

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

    phi_new[lev]->FillBoundary(geom[lev].periodicity());
    fill_physbc(*phi_new[lev], geom[lev]);
    phi_old[lev]->FillBoundary(geom[lev].periodicity());
    fill_physbc(*phi_old[lev], geom[lev]);
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
            int nghost = phi_new[lev+1]->nGrow();
            BoxArray ba = grids[lev+1];
            std::unique_ptr<MultiFab> phi_swe(new MultiFab(ba,
                    dmap[lev+1], n_swe_comp, nghost));

            std::cout << "averagedowning here \n\n\n";

            swe_from_comp_wrapper(lev+1, *phi_swe, *phi_new[lev+1]);

            amrex::average_down(*phi_swe, *phi_new[lev],
 			     geom[lev+1], geom[lev],
 			     0, n_swe_comp, refRatio(lev));

                 phi_new[lev]->FillBoundary(geom[lev].periodicity());
                 fill_physbc(*phi_new[lev], geom[lev]);

        } else {
            amrex::average_down(*phi_new[lev+1], *phi_new[lev],
 			     geom[lev+1], geom[lev],
 			     0, phi_new[lev]->nComp(), refRatio(lev));

                 phi_new[lev]->FillBoundary(geom[lev].periodicity());
                 fill_physbc(*phi_new[lev], geom[lev]);
        }
    }
}

void
AmrAdv::AverageDownTo (int crse_lev)
{

    if (crse_lev == max_swe_level) {
        // convert first
        int n_swe_comp = phi_new[crse_lev]->nComp();
        int nghost = phi_new[crse_lev+1]->nGrow();
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

        mf.FillBoundary(geom[lev].periodicity());
        fill_physbc(mf, geom[lev]);
    }
    else
    {
    	Array<MultiFab*> cmf, fmf;
    	Array<Real> ctime, ftime;
    	GetData(lev-1, time, cmf, ctime);
    	GetData(lev  , time, fmf, ftime);

    	AmrAdvPhysBC cphysbc, fphysbc;
    	Interpolater* mapper = &cell_cons_interp;

        int lo_bc[] = {INT_DIR, INT_DIR, HOEXTRAP}; // outflow boundaries
        int hi_bc[] = {INT_DIR, INT_DIR, HOEXTRAP};
    	Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

        if (lev-1 == max_swe_level) {
            // do comp from swe conversion
            BoxArray ba = grids[lev-1];
            Array<MultiFab*> comp_mf;

            for (int i = 0; i< cmf.size(); i++) {
               //comp_mf.push_back(new MultiFab(grids[lev],
                //    dmap[lev], phi_new[lev]->nComp(),
                //    phi_new[lev]->nGrow()));
                comp_mf.push_back(new MultiFab(ba,
                     dmap[lev-1], phi_new[lev]->nComp(),
                     phi_new[lev]->nGrow()));
                comp_from_swe_wrapper(lev-1, *cmf[i], *comp_mf[i]);
                comp_mf[i]->FillBoundary(geom[lev].periodicity());

                fill_physbc(*comp_mf[i], geom[lev]);
            }

        	amrex::FillPatchTwoLevels(mf, time, comp_mf, ctime,
                           fmf, ftime,
        				   0, icomp, phi_new[lev]->nComp(),
                           geom[lev-1], geom[lev],
        				   cphysbc, fphysbc, refRatio(lev-1),
        				   mapper, bcs);
        } else {
            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
        				   0, icomp, ncomp, geom[lev-1], geom[lev],
        				   cphysbc, fphysbc, refRatio(lev-1),
        				   mapper, bcs);
        }

        mf.FillBoundary(geom[lev].periodicity());
        fill_physbc(mf, geom[lev]);
    }
}

void
AmrAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    std::cout << "Filling coarse patch level " << lev << '\n';
    BL_ASSERT(lev > 0);

    Array<MultiFab*> cmf;
    Array<Real> ctime;
    GetData(lev-1, time, cmf, ctime);

    if (cmf.size() != 1) {
	       amrex::Abort("FillCoarsePatch: how did this happen?");
    }

    AmrAdvPhysBC cphysbc, fphysbc;
    Interpolater* mapper = &cell_cons_interp;

    //int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaries
    //int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    int lo_bc[] = {INT_DIR, INT_DIR, HOEXTRAP}; // outflow boundaries
    int hi_bc[] = {INT_DIR, INT_DIR, HOEXTRAP};
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

    //amrex::FillDomainBoundary(mf, geom[lev-1], bcs);
    mf.FillBoundary(geom[lev].periodicity());
    fill_physbc(mf, geom[lev]);
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

    std::cout << "*********comp from swe*********\n\n";

    const int n_cons_comp = 5;
    const int n_swe_comp = 3;
    int nghost = swe_mf.nGrow();
    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    for (MFIter mfi(swe_mf, true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const FArrayBox& U_swe      =   swe_mf[mfi];
        FArrayBox& U_comp      =   comp_mf[mfi];

        comp_from_swe(BL_TO_FORTRAN_3D(U_comp),
            BL_TO_FORTRAN_3D(U_swe),
            p, rho,
            bx.loVect(), bx.hiVect(),
            &n_cons_comp, &n_swe_comp,
            &gamma, ZFILL(dx), &alpha0, &M, &R, &nghost, ZFILL(prob_lo));
    }
    comp_mf.FillBoundary(geom[lev].periodicity());
    fill_physbc(comp_mf, geom[lev]);
}

void AmrAdv::swe_from_comp_wrapper(int lev, MultiFab& swe_mf, MultiFab& comp_mf) const {

    std::cout << "*********swe from comp*********\n\n";

    const int n_cons_comp = 5;
    const int n_swe_comp = 3;
    const Real* dx = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    int n_com_comp = comp_mf.nComp();
    int nghost = comp_mf.nGrow();
    BoxArray ba = grids[lev];
    MultiFab U_prim_mf(ba, dmap[lev], n_com_comp, nghost);
    MultiFab p_mf(ba, dmap[lev], 1, nghost);

    for (MFIter mfi(comp_mf, true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        FArrayBox& U_swe      =   swe_mf[mfi];
        FArrayBox& U_comp      =   comp_mf[mfi];
        FArrayBox& U_prim      =   U_prim_mf[mfi];
        FArrayBox& p_comp      =   p_mf[mfi];

        // do prim conversion here first
        cons_to_prim(BL_TO_FORTRAN_3D(U_comp),
            BL_TO_FORTRAN_3D(U_prim),
            BL_TO_FORTRAN_3D(p_comp),
            bx.loVect(), bx.hiVect(),
            &n_cons_comp, &gamma, &alpha0, &M, &R, ZFILL(dx), ZFILL(prob_lo));

        // then calculate swe
        swe_from_comp(BL_TO_FORTRAN_3D(U_prim),
            BL_TO_FORTRAN_3D(U_swe),
            BL_TO_FORTRAN_3D(p_comp),
            p,
            bx.loVect(), bx.hiVect(),
            &n_cons_comp, &n_swe_comp,
            &alpha0, &M, &R,
            ZFILL(dx), ZFILL(prob_lo));
    }
    swe_mf.FillBoundary(geom[lev].periodicity());
    fill_physbc(swe_mf, geom[lev]);

    //swe_mf.FillBoundary(0, n_swe_comp, geom[lev].periodicity());
}

void AmrAdv::InitAmrMesh (int max_level_in,
        const Array<int>& n_cell_in,
        std::vector<int> refrat) {

    AmrCore::InitAmrMesh(max_level_in, n_cell_in, refrat);

    std::cout << ref_ratio.size() << '\n';

    //int m = maxLevel();
    //int s = max_swe_level+1;
    //int n_swe_layers = std::min(m, s);
    // rehack ref ratios of swe layers
    for (int i = 0; i < std::min(max_level, max_swe_level); i++) {
        ref_ratio[i][2] = 1;
    }

    for (int i = 0; i < max_level; i++) {
        std::cout << "ref ratio 2: " << ref_ratio[i] << '\n';
    }

    ParmParse pp("amr");

    Array<int> n_cell(BL_SPACEDIM);
    Array<int> n_cell_swe(BL_SPACEDIM);
	if (n_cell_in[0] == -1)
	{
	    pp.getarr("n_cell",n_cell,0,BL_SPACEDIM);
	}
	else
	{
	    for (int i = 0; i < BL_SPACEDIM; i++) {
            n_cell[i] = n_cell_in[i];
        }
	}

    for (int i = 0; i < BL_SPACEDIM; i++) {
        n_cell_swe[i] = n_cell[i];
    }

    //std::cout << "nlayers: " << nlayers << '\n';
    if (0 <= max_swe_level) {
        n_cell_swe[2] = nlayers;
    }
    IntVect lo(IntVect::TheZeroVector()), hi(n_cell_swe);
    hi -= IntVect::TheUnitVector();
    Box index_domain(lo,hi);

    for (int i = 0; i <= max_level; i++)
    {
        geom[i].define(index_domain);
        //std::cout << "index domain " << index_domain << '\n';
        if (i < max_level) index_domain.refine(ref_ratio[i]);


    }

    Real offset[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        const Real delta = Geometry::ProbLength(i)/(Real)n_cell[i];
        offset[i]        = Geometry::ProbLo(i) + delta*lo[i];
    }
    CoordSys::SetOffset(offset);

    {
	// chop up grids to have more grids than the number of procs
	pp.query("refine_grid_layout", refine_grid_layout);
    }

    for (int i = 0; i < BL_SPACEDIM; i++) {
        std::cout << "n_cell " << i << ' ' << n_cell[i] << '\n';
        std::cout << "n_cell_swe " << i << ' ' << n_cell_swe[i] << '\n';
    }

    std::cout << "printing geometry " << geom[0] << '\n';
    std::cout << "printing geometry domain " << geom[0].Domain() << '\n';

}
