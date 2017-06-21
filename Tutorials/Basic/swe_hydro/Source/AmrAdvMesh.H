#ifndef AMR_Adv_Mesh_H_
#define AMR_Adv_Mesh_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrMesh.H>

class AmrAdvMesh
    : public amrex::AmrMesh
{
public:
    AmrAdvMesh();

    void
    MakeNewGrids (amrex::Real time);

    void
    MakeNewGrids (int lbase, amrex::Real time, int& new_finest, amrex::Array<amrex::BoxArray>& new_grids);

    void InitAmrMesh (int max_level_in,
            const amrex::Array<int>& n_cell_in,
            std::vector<int> refrat = std::vector<int>());

private:

    int max_swe_level = -1;
    int nlayers = 1;

};

#endif
