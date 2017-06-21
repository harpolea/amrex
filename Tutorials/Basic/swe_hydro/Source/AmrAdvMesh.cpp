
#include <AMReX.H>
#include <AMReX_Cluster.H>
#include <AMReX_ParmParse.H>

#include <AmrAdv.H>
#include <AmrAdvBC.H>
#include <AmrAdvMesh.H>

using namespace amrex;

AmrAdvMesh::AmrAdvMesh() : AmrMesh() {
	ParmParse pp("swe");

    pp.query("max_swe_level", max_swe_level);
    pp.query("nlayers", nlayers);
}


void
AmrAdvMesh::MakeNewGrids (int lbase, Real time, int& new_finest, Array<BoxArray>& new_grids)
{
    BL_ASSERT(lbase < max_level);

    // Add at most one new level
    int max_crse = std::min(finest_level, max_level-1);

    if (new_grids.size() < max_crse+2) new_grids.resize(max_crse+2);

    //
    // Construct problem domain at each level.
    //
    Array<IntVect> bf_lev(max_level); // Blocking factor at each level.
    Array<IntVect> rr_lev(max_level);
    Array<Box>     pc_domain(max_level);  // Coarsened problem domain.

    for (int i = 0; i <= max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            bf_lev[i][n] = std::max(1,blocking_factor[i+1]/ref_ratio[i][n]);
    }
    for (int i = lbase; i < max_crse; i++)
    {
        for (int n=0; n<BL_SPACEDIM; n++)
            rr_lev[i][n] = (ref_ratio[i][n]*bf_lev[i][n])/bf_lev[i+1][n];
    }
    for (int i = lbase; i <= max_crse; i++) {
	pc_domain[i] = amrex::coarsen(Geom(i).Domain(),bf_lev[i]);
    }
    //
    // Construct proper nesting domains.
    //
    Array<BoxList> p_n(max_level);      // Proper nesting domain.
    Array<BoxList> p_n_comp(max_level); // Complement proper nesting domain.

    BoxList bl(grids[lbase]);
    bl.simplify();
    bl.coarsen(bf_lev[lbase]);
    p_n_comp[lbase].complementIn(pc_domain[lbase],bl);
    p_n_comp[lbase].simplify();
    p_n_comp[lbase].accrete(n_proper);
    if (Geometry::isAnyPeriodic()) {
	       ProjPeriodic(p_n_comp[lbase], Geometry(pc_domain[lbase]));
    }
    p_n[lbase].complementIn(pc_domain[lbase],p_n_comp[lbase]);
    p_n[lbase].simplify();
    bl.clear();

    for (int i = lbase+1; i <= max_crse; i++)
    {
        p_n_comp[i] = p_n_comp[i-1];

        // Need to simplify p_n_comp or the number of grids can too large for many levels.
        p_n_comp[i].simplify();

        p_n_comp[i].refine(rr_lev[i-1]);
        p_n_comp[i].accrete(n_proper);

	if (Geometry::isAnyPeriodic()) {
	    ProjPeriodic(p_n_comp[i], Geometry(pc_domain[i]));
	}

        p_n[i].complementIn(pc_domain[i],p_n_comp[i]);
        p_n[i].simplify();
    }
    //
    // Now generate grids from finest level down.
    //
    new_finest = lbase;

    for (int levc = max_crse; levc >= lbase; levc--)
    {
        int levf = levc+1;
        //
        // Construct TagBoxArray with sufficient grow factor to contain
        // new levels projected down to this level.
        //
        int ngrow = 0;

        if (levf < new_finest)
        {
            BoxArray ba_proj(new_grids[levf+1]);

            ba_proj.coarsen(ref_ratio[levf]);
            ba_proj.growcoarsen(n_proper, ref_ratio[levc]);

            BoxArray levcBA = grids[levc];

            while (!levcBA.contains(ba_proj))
            {
                levcBA.grow(1);
                ++ngrow;
            }
        }
        TagBoxArray tags(grids[levc],dmap[levc],n_error_buf[levc]+ngrow);

        //
        // Only use error estimation to tag cells for the creation of new grids
        //      if the grids at that level aren't already fixed.
        //

        if ( ! (useFixedCoarseGrids() && levc < useFixedUpToLevel()) ) {
	    ErrorEst(levc, tags, time, ngrow);
	}

        //
        // If new grids have been constructed above this level, project
        // those grids down and tag cells on intersections to ensure
        // proper nesting.
        //
        // NOTE: this loop replaces the previous code:
        //      if (levf < new_finest)
        //          tags.setVal(ba_proj,TagBox::SET);
        // The problem with this code is that it effectively
        // "buffered the buffer cells",  i.e., the grids at level
        // levf+1 which were created by buffering with n_error_buf[levf]
        // are then coarsened down twice to define tagging at
        // level levc, which will then also be buffered.  This can
        // create grids which are larger than necessary.
        //
        if (levf < new_finest)
        {
            int nerr = n_error_buf[levf];

            BoxList bl_tagged(new_grids[levf+1]);
            bl_tagged.simplify();
            bl_tagged.coarsen(ref_ratio[levf]);
            //
            // This grows the boxes by nerr if they touch the edge of the
            // domain in preparation for them being shrunk by nerr later.
            // We want the net effect to be that grids are NOT shrunk away
            // from the edges of the domain.
            //
            for (BoxList::iterator blt = bl_tagged.begin(), End = bl_tagged.end();
                 blt != End;
                 ++blt)
            {
                for (int idir = 0; idir < BL_SPACEDIM; idir++)
                {
                    if (blt->smallEnd(idir) == Geom(levf).Domain().smallEnd(idir))
                        blt->growLo(idir,nerr);
                    if (blt->bigEnd(idir) == Geom(levf).Domain().bigEnd(idir))
                        blt->growHi(idir,nerr);
                }
            }
            Box mboxF = amrex::grow(bl_tagged.minimalBox(),1);
            BoxList blFcomp;
            blFcomp.complementIn(mboxF,bl_tagged);
            blFcomp.simplify();
            bl_tagged.clear();

            const IntVect& iv = IntVect(AMREX_D_DECL(nerr/ref_ratio[levf][0],
                                               nerr/ref_ratio[levf][1],
                                               nerr/ref_ratio[levf][2]));
            blFcomp.accrete(iv);
            BoxList blF;
            blF.complementIn(mboxF,blFcomp);
            BoxArray baF(blF);
            blF.clear();
            baF.grow(n_proper);
            //
            // We need to do this in case the error buffering at
            // levc will not be enough to cover the error buffering
            // at levf which was just subtracted off.
            //
            for (int idir = 0; idir < BL_SPACEDIM; idir++)
            {
                if (nerr > n_error_buf[levc]*ref_ratio[levc][idir])
                    baF.grow(idir,nerr-n_error_buf[levc]*ref_ratio[levc][idir]);
            }

            baF.coarsen(ref_ratio[levc]);

            tags.setVal(baF,TagBox::SET);
        }
        //
        // Buffer error cells.
        //
        tags.buffer(n_error_buf[levc]+ngrow);

        if (useFixedCoarseGrids())
        {
    	    if (levc>=useFixedUpToLevel())
    	    {
    		          tags.setVal(GetAreaNotToTag(levc), TagBox::CLEAR);
    	    }
    	    else
    	    {
    		          new_finest = std::max(new_finest,levf);
    	    }
        }

        //
        // Coarsen the taglist by blocking_factor/ref_ratio.
        //
        int bl_max = 0;
        for (int n=0; n<BL_SPACEDIM; n++) {
            bl_max = std::max(bl_max,bf_lev[levc][n]);
        }
        if (bl_max >= 1) {
            tags.coarsen(bf_lev[levc]);
        } else {
            amrex::Abort("blocking factor is too small relative to ref_ratio");
        }
        //
        // Remove or add tagged points which violate/satisfy additional
        // user-specified criteria.
        //
	ManualTagsPlacement(levc, tags, bf_lev);
        //
        // Map tagged points through periodic boundaries, if any.
        //
        tags.mapPeriodic(Geometry(pc_domain[levc]));
        //
        // Remove cells outside proper nesting domain for this level.
        //
        tags.setVal(p_n_comp[levc],TagBox::CLEAR);
        //
        // Create initial cluster containing all tagged points.
        //
	std::vector<IntVect> tagvec;
	tags.collate(tagvec);
        tags.clear();

        if (tagvec.size() > 0)
        {
            //
            // Created new level, now generate efficient grids.
            //
            if ( !(useFixedCoarseGrids() && levc<useFixedUpToLevel()) ) {
                new_finest = std::max(new_finest,levf);
	    }
            //
            // Construct initial cluster.
            //
            ClusterList clist(&tagvec[0], tagvec.size());
            clist.chop(grid_eff);
            BoxDomain bd;
            bd.add(p_n[levc]);
            clist.intersect(bd);
            bd.clear();
            //
            // Efficient properly nested Clusters have been constructed
            // now generate list of grids at level levf.
            //
            BoxList new_bx;
            clist.boxList(new_bx);
            new_bx.refine(bf_lev[levc]);
            new_bx.simplify();
            BL_ASSERT(new_bx.isDisjoint());

	    if (new_bx.size()>0) {
    		if ( !(Geom(levc).Domain().contains(BoxArray(new_bx).minimalBox())) ) {
    		// Chop new grids outside domain, note that this is likely to result in
    		//  new grids that violate blocking_factor....see warning checking below
    		    new_bx = amrex::intersect(new_bx,Geom(levc).Domain());
    		}
	    }

            IntVect largest_grid_size;
            for (int n = 0; n < BL_SPACEDIM; n++) {
                largest_grid_size[n] = max_grid_size[levf] / ref_ratio[levc][n];
            }
            //
            // Ensure new grid boxes are at most max_grid_size in index dirs.
            //
            new_bx.maxSize(largest_grid_size);

#ifdef BL_FIX_GATHERV_ERROR
	      int wcount = 0, iLGS = largest_grid_size[0];

              while (new_bx.size() < 64 && wcount++ < 4)
              {
                  iLGS /= 2;
		  amrex::Print() << "BL_FIX_GATHERV_ERROR:  using iLGS = " << iLGS
				 << "   largest_grid_size was:  " << largest_grid_size[0] << '\n'
				 << "BL_FIX_GATHERV_ERROR:  new_bx.size() was:   "
				 << new_bx.size() << '\n';

                  new_bx.maxSize(iLGS);

		  amrex::Print() << "BL_FIX_GATHERV_ERROR:  new_bx.size() now:   "
				 << new_bx.size() << '\n';
	      }
#endif
            //
            // Refine up to levf.
            //
            new_bx.refine(ref_ratio[levc]);
            BL_ASSERT(new_bx.isDisjoint());

	    if (new_bx.size()>0) {
    		if ( !(Geom(levf).Domain().contains(BoxArray(new_bx).minimalBox())) ) {
    		    new_bx = amrex::intersect(new_bx,Geom(levf).Domain());
    		}
#if 0
// Let's not check this, because of the hack that uses blocking factor larger than max grid size.
		if (ParallelDescriptor::IOProcessor()) {
		    for (int d=0; d<BL_SPACEDIM; ++d) {
			bool ok = true;
			for (BoxList::const_iterator bli = new_bx.begin(); bli != new_bx.end(); ++bli) {
			    int len = bli->length(d);
			    int bf = blocking_factor[levf];
			    ok &= (len/bf) * bf == len;
			}
			if (!ok) {
			    amrex::Warning("WARNING: New grids violate blocking factor near upper boundary");
			}
		    }
		}
#endif
	    }

            if(levf > useFixedUpToLevel()) {
              new_grids[levf].define(new_bx);
	         }
        }
    }

    for (int lev = lbase+1; lev <= new_finest; ++lev) {
        if (new_grids[lev].empty())
        {
            if (!(useFixedCoarseGrids() && lev<useFixedUpToLevel()) ) {
                amrex::Abort("AmrMesh::MakeNewGrids: how did this happen?");
            }
        }
        else if (refine_grid_layout)
        {
            ChopGrids(lev,new_grids[lev],ParallelDescriptor::NProcs());
            if (new_grids[lev] == grids[lev]) {
                new_grids[lev] = grids[lev]; // to avoid dupliates
            }
        }
    }
}

void
AmrAdvMesh::MakeNewGrids (Real time)
{
    // define coarse level BoxArray and DistributionMap
    {
	finest_level = 0;

	const BoxArray& ba = MakeBaseGrids();
	DistributionMapping dm(ba);

	MakeNewLevelFromScratch(0, time, ba, dm);

	SetBoxArray(0, ba);
	SetDistributionMap(0, dm);
    }

    if (max_level > 0) // build fine levels
    {
    	Array<BoxArray> new_grids(max_level+1);
    	new_grids[0] = grids[0];
    	do
    	{
    	    int new_finest;

    	    // Add (at most) one level at a time.
    	    MakeNewGrids(finest_level,time,new_finest,new_grids);

    	    if (new_finest <= finest_level) break;
    	    finest_level = new_finest;

    	    DistributionMapping dm(new_grids[new_finest]);

                MakeNewLevelFromScratch(new_finest, time, new_grids[finest_level], dm);

    	    SetBoxArray(new_finest, new_grids[new_finest]);
    	    SetDistributionMap(new_finest, dm);
    	}
    	while (finest_level < max_level);

    	// Iterate grids to ensure fine grids encompass all interesting junk.
    	for (int it=0; it<4; ++it)  // try at most 4 times
    	{
    	    for (int i = 1; i <= finest_level; ++i) {
    		          new_grids[i] = grids[i];
    	    }

    	    int new_finest;
    	    MakeNewGrids(0, time, new_finest, new_grids);

    	    if (new_finest < finest_level) break;
    	    finest_level = new_finest;

    	    bool grids_the_same = true;
    	    for (int lev = 1; lev <= new_finest; ++lev) {
        		if (new_grids[lev] != grids[lev]) {
        		    grids_the_same = false;
        		    DistributionMapping dm(new_grids[lev]);

                    MakeNewLevelFromScratch(lev, time, new_grids[lev], dm);

        		    SetBoxArray(lev, new_grids[lev]);
        		    SetDistributionMap(lev, dm);
        		}
    	    }
    	    if (grids_the_same) break;
    	}
    }
}


void AmrAdvMesh::InitAmrMesh (int max_level_in,
        const Array<int>& n_cell_in,
        std::vector<int> refrat) {

    AmrMesh::InitAmrMesh(max_level_in, n_cell_in, refrat);

    // rehack ref ratios of swe layers
    for (int i = 0; i < max_swe_level; i++) {

        ref_ratio[i][2] = 1;
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
	    for (int i = 0; i < BL_SPACEDIM; i++) n_cell[i] = n_cell_in[i];
	}

    std::cout << "nlayers: " << nlayers << '\n';
    n_cell_swe[2] = nlayers;
    IntVect lo(IntVect::TheZeroVector()), hi(n_cell);

    for (int i = 0; i <= max_level; i++)
    {
        if (i <= max_swe_level) {
            hi = IntVect(n_cell_swe);
        }

        hi -= IntVect::TheUnitVector();
        Box index_domain(lo,hi);
        geom[i].define(index_domain);
        if (i < max_level) index_domain.refine(ref_ratio[i]);
    }

    Real offset[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        const Real delta = Geometry::ProbLength(i)/(Real)n_cell[i];
        offset[i]        = Geometry::ProbLo(i) + delta*lo[i];
    }
    CoordSys::SetOffset(offset);


}
