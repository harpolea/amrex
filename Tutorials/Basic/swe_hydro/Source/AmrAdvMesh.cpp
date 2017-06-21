
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
