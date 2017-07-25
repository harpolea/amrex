#ifndef AMR_Adv_H_
#define AMR_Adv_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_BoxArray.H>

class AmrAdv
    : public amrex::AmrCore
{
public:
    AmrAdv ();
    virtual ~AmrAdv();

    void InitData ();
    void Evolve ();

    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					  const amrex::DistributionMapping& dm) override;
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
					 const amrex::DistributionMapping& dm) override;
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
			      const amrex::DistributionMapping& dm) override;
    virtual void ClearLevel (int lev) override;
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) override;

    virtual void InitFromScratch (amrex::Real time) override;
    //virtual amrex::BoxArray MakeBaseGrids () const override;

    void fill_physbc (amrex::MultiFab& phi, const amrex::Geometry& geom) const;


protected:

    virtual void InitAmrMesh (int max_level_in,
                const amrex::Array<int>& n_cell_in,
                std::vector<int> refrat = std::vector<int>()) override;

private:

    void ReadParameters();

    void InitFromCheckpoint ();

    void AverageDown ();
    void AverageDownTo (int crse_lev);
    long CountCells (int lev);

    void ComputeDt ();
    amrex::Real EstTimeStep (int lev, bool local=false) const;
    void timeStep (int lev, amrex::Real time, int iteration);
    void Advance (int lev, amrex::Real time, amrex::Real dt, int iteration, int ncycle);

    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);
    void GetData (int lev, amrex::Real time, amrex::Array<amrex::MultiFab*>& data, amrex::Array<amrex::Real>& datatime);
    void comp_from_swe_wrapper(int lev, amrex::MultiFab& swe_mf, amrex::MultiFab& comp_mf);
    void swe_from_comp_wrapper(int lev, amrex::MultiFab& swe_mf, amrex::MultiFab& comp_mf) const;

    void average_down(const amrex::MultiFab& S_fine, amrex::MultiFab& S_crse,
                       int scomp, int ncomp, const amrex::IntVect& ratio);

    std::string PlotFileName (int lev) const;
    amrex::Array<const amrex::MultiFab*> PlotFileMF () const;
    amrex::Array<std::string> PlotFileVarNames () const;
    void WritePlotFile () const;
    void
    WriteMultiLevelPlotfile (const std::string& plotfilename,
                             int nlevels,
                             const amrex::Array<const amrex::MultiFab*>& mf,
                             const amrex::Array<std::string>& varnames,
                             const amrex::Array<amrex::Geometry>& geom, Real time, const amrex::Array<int>& level_steps,
                             const amrex::Array<amrex::IntVect>& ref_ratio,
                             const std::string &versionName = "amrex_swe_v1",
                             const std::string &levelPrefix = "Level_",
                             const std::string &mfPrefix = "Cell") const;

    amrex::Array<int> istep;      // which step?
    amrex::Array<int> nsubsteps;  // how many substeps on each level?

    amrex::Array<amrex::Real> t_new;
    amrex::Array<amrex::Real> t_old;
    amrex::Array<amrex::Real> dt;

    amrex::Array<std::unique_ptr<amrex::MultiFab> > phi_new;
    amrex::Array<std::unique_ptr<amrex::MultiFab> > phi_old;

    amrex::Array<std::unique_ptr<amrex::FluxRegister> > flux_reg;

    // runtime parameters

    int max_step   = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    amrex::Real cfl = 0.7;

    int regrid_int = 2;

    int do_reflux = 1;

    std::string restart_chkfile;

    std::string check_file {"chk"};
    std::string plot_file {"plt"};
    int check_int = -1;
    int plot_int = -1;

    int max_swe_level = -1;
    int nlayers = 1;

    Real *rho;
    Real *p;

    amrex::Real alpha0, M, R, gamma;

    bool gr = false;

};

#endif
