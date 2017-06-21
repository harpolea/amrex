
#include <AMReX_ParmParse.H>

#include <AmrAdv.H>
#include <AmrAdv_F.H>

using namespace amrex;

void
AmrAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    static bool first = true;
    static Array<Real> phierr;

    const int ncomp = phi_new[lev]->nComp();

    if (first)
    {
    	first = false;
    	ParmParse pp("adv");
        if (ncomp < 5) {
            int n = pp.countval("swe_phierr");
        	if (n > 0) {
        	    pp.getarr("swe_phierr", phierr, 0, n);
        	}
        } else {
            int n = pp.countval("phierr");
            if (n > 0) {
                pp.getarr("phierr", phierr, 0, n);
            }
        }
    }

    if (lev >= phierr.size()) return;

    const int clearval = TagBox::CLEAR;
    const int   tagval = TagBox::SET;

    const Real* dx      = geom[lev].CellSize();
    const Real* prob_lo = geom[lev].ProbLo();

    const MultiFab& state = *phi_new[lev];


#ifdef _OPENMP
#pragma omp parallel
#endif
    {
            Array<int>  itags;

    	for (MFIter mfi(state,true); mfi.isValid(); ++mfi)
    	{
    	    const Box&  tilebx  = mfi.tilebox();

                TagBox&     tagfab  = tags[mfi];

    	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
    	    // So we are going to get a temporary integer array.
    	    tagfab.get_itags(itags, tilebx);

                // data pointer and index space
    	    int*        tptr    = itags.dataPtr();
    	    const int*  tlo     = tilebx.loVect();
    	    const int*  thi     = tilebx.hiVect();
    	    state_error(tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
    			BL_TO_FORTRAN_3D(state[mfi]),
    			&tagval, &clearval,
    			ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()),
    			ZFILL(dx), ZFILL(prob_lo), &time, &phierr[lev], &ncomp);
    	    //
    	    // Now update the tags in the TagBox.
    	    //
    	    tagfab.tags_and_untags(itags, tilebx);
    	}
    }
}
