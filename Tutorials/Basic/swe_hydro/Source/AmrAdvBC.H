#ifndef AmrAdvBC_H_
#define AmrAdvBC_H_

#include <AMReX_PhysBCFunct.H>
#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>

// In this tutorial, we assume the boundaries are all periodic.
// Thus there are no physical boundaries. And FillBoundary does nothing.

class AmrAdvPhysBC
    : public amrex::PhysBCFunctBase
{
public:
    AmrAdvPhysBC () {}
    virtual ~AmrAdvPhysBC () {}
    virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override {}
    using amrex::PhysBCFunctBase::FillBoundary;

    virtual void FillBoundary (amrex::MultiFab& mf, amrex::Real time) override;

    virtual void FillBoundary (amrex::MultiFab& mf);

private:
    amrex::Geometry m_geom;
    amrex::BCRec    m_bcr;
    std::unique_ptr<amrex::BndryFunctBase> m_bc_func;
};

#endif
