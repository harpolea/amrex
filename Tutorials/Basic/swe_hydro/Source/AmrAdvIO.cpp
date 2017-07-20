
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include "AmrAdv.H"

using namespace amrex;

std::string
AmrAdv::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

Array<const MultiFab*>
AmrAdv::PlotFileMF () const
{
    Array<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	       r.push_back(phi_new[i].get());
    }
    return r;
}

Array<std::string>
AmrAdv::PlotFileVarNames () const
{
    if (phi_new[0]->nComp() < 5) {

        return {"h", "hu", "hv"};
    } else {
        return {"rho", "rhou", "rhov", "rhow", "E"};
    }
}

void
AmrAdv::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
				    Geom(), t_new[0], istep, refRatio());
}

void
AmrAdv::InitFromCheckpoint ()
{
    amrex::Abort("AmrAdv::InitFromCheckpoint: todo");
}

void
AmrAdv::WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
                         const Array<const MultiFab*>& mf,
                         const Array<std::string>& varnames,
                         const Array<Geometry>& geom, Real time, const Array<int>& level_steps,
                         const Array<IntVect>& ref_ratio,
                         const std::string &versionName,
                         const std::string &levelPrefix,
                         const std::string &mfPrefix)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int finest_level = nlevels-1;

    //
    // Only let 64 CPUs be writing at any one time.
    //
    int saveNFiles(VisMF::GetNOutFiles());
    VisMF::SetNOutFiles(64);

    bool callBarrier(true);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);

    if (ParallelDescriptor::IOProcessor()) {
      std::string HeaderFileName(plotfilename + "/Header");
      std::ofstream HeaderFile(HeaderFileName.c_str(), std::ofstream::out   |
	                 std::ofstream::trunc | std::ofstream::binary);
      if( ! HeaderFile.good()) {
          FileOpenFailed(HeaderFileName);
      }

      Array<BoxArray> boxArrays(nlevels);
      for(int level(0); level < boxArrays.size(); ++level) {
	         boxArrays[level] = mf[level]->boxArray();
      }

      WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                 geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix);
    }

    for (int level = 0; level <= finest_level; ++level)
    {
        MultiFab* data;
        std::unique_ptr<MultiFab> mf_tmp;
        if (mf[level]->nGrow() > 0) {
            mf_tmp.reset(new MultiFab(mf[level]->boxArray(),
                                      mf[level]->DistributionMap(),
                                      mf[level]->nComp(), 0));
            MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
            data = mf_tmp.get();
        } else {
            MultiFab::Copy(*data, *mf[level], 0, 0, mf[level]->nComp(), 0);
            //data = mf[level];
        }

        if (level > max_swe_level) {
            int n_swe_comp = phi_new[max_swe_level]->nComp();
            int nghost = phi_new[max_swe_level]->nGrow();
            BoxArray ba = grids[level];
            std::unique_ptr<MultiFab> phi_swe(new MultiFab(ba,
                    dmap[level], n_swe_comp, nghost));

            swe_from_comp_wrapper(level, *phi_swe, *data);
            MultiFab::Copy(*data, *phi_swe, 0, 0, n_swe_comp, 0);
        }
	    VisMF::Write(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

    VisMF::SetNOutFiles(saveNFiles);
}
