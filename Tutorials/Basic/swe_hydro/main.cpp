
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFabUtil.H>


#include <array>

#include <fstream>
#include <iomanip>

#include "myfunc_F.H"
#include <MyAmr.H>

using namespace amrex;

void main_main ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void advance (MultiFab& old_data, MultiFab& new_data,
	      Real dt, const Geometry& geom, Real gamma,
          MultiFab& gamma_up, MultiFab& beta, MultiFab& alpha,
          Real alpha0, Real M, Real R)
{
    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    old_data.FillBoundary(geom.periodicity());

    int Ncomp = old_data.nComp();
    int ng_p = old_data.nGrow();

    const Real* dx = geom.CellSize();

    // Advance the solution one grid at a time
    for ( MFIter mfi(old_data); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        update_data(bx.loVect(), bx.hiVect(), &Ncomp,
                   BL_TO_FORTRAN_ANYD(old_data[mfi]),
                   BL_TO_FORTRAN_ANYD(new_data[mfi]),
                   dx, dt, &gamma,
                   BL_TO_FORTRAN_ANYD(gamma_up[mfi]),
                   BL_TO_FORTRAN_ANYD(beta[mfi]),
                   BL_TO_FORTRAN_ANYD(alpha[mfi]),
                   alpha0, M, R);
    }
}

void main_main ()
{
    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    // BL_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, nsteps, plot_int, is_periodic[BL_SPACEDIM], nlayers;

    Real alpha0, M, R, gamma, min_compressible_resolution;

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be writtenq
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Default nsteps to 0, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        // stellar mass
        M = 1.0;
        pp.query("M", M);

        // stellar radius
        R = 100.0;
        pp.query("R", R);

        // lapse function
        alpha0 = sqrt(1.0 - 2 * M / R);
        pp.query("alpha0", alpha0);

        // adiabatic index
        gamma = 5.0/3.0;
        pp.query("gamma", gamma);

        // nlayers
        pp.get("nlayers", nlayers);

        // get models
        min_compressible_resolution = 5.0e-3;
        pp.query("min_compressible_resolution", min_compressible_resolution);
    }

    //Print() << R << ' ' << alpha0 << ' ' << M << '\n';
    //return;

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(IntVect(AMREX_D_DECL(0,0,0)));
        IntVect dom_hi(IntVect(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

        // This defines the physical size of the box.  Right now the box is [-1,1] in each direction.
        RealBox real_box;
        for (int n = 0; n < BL_SPACEDIM; n++) {
            real_box.setLo(n,-1.0);
            real_box.setHi(n, 1.0);
        }

        // This says we are using Cartesian coordinates
        int coord = 0;

        // This sets the boundary conditions to be doubly or triply periodic
        int is_periodic[BL_SPACEDIM];
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            is_periodic[i] = 1;
        }

        // This defines a Geometry object
        geom.define(domain,&real_box,coord,is_periodic);
    }

    // Nghost = number of ghost cells for each array
    // Need 6 ghosts for slope limited rk3
    int Nghost = 6;

    // Ncomp = number of components for each array
    int Ncomp  = 2 + BL_SPACEDIM;

    // time = starting time in the simulation
    Real time = 0.0;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two data multifabs; one will store the old state, the other the new.
    MultiFab data_old(ba, dm, Ncomp, Nghost);
    MultiFab data_new(ba, dm, Ncomp, Nghost);

    data_old.setVal(0.0);
    data_new.setVal(0.0);

    // make metric multifabs
    MultiFab gamma_up(ba, dm, 9, Nghost);
    MultiFab beta(ba, dm, 3, Nghost);
    MultiFab alpha(ba, dm, 1, Nghost);

    gamma_up.setVal(0.0);
    beta.setVal(0.0);
    alpha.setVal(0.0);

    // Initialize data_new by calling a Fortran routine.
    // MFIter = MultiFab Iterator
    for ( MFIter mfi(data_new); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        const Real* dx = geom.CellSize();

        if (std::min(std::min(dx[0], dx[1]), dx[2]) > min_compressible_resolution) {
            init_swe_data(BL_TO_FORTRAN_ANYD(data_new[mfi]),
                      bx.loVect(), bx.hiVect(), &Ncomp,
                     geom.CellSize(), geom.ProbLo(), geom.ProbHi(),
                     BL_TO_FORTRAN_ANYD(gamma_up[mfi]),
                     BL_TO_FORTRAN_ANYD(beta[mfi]),
                     BL_TO_FORTRAN_ANYD(alpha[mfi]),
                     alpha0, M, R);
        } else {
            init_comp_data(BL_TO_FORTRAN_ANYD(data_new[mfi]),
                      bx.loVect(), bx.hiVect(), &Ncomp,
                     geom.CellSize(), geom.ProbLo(), geom.ProbHi(),
                     BL_TO_FORTRAN_ANYD(gamma_up[mfi]),
                     BL_TO_FORTRAN_ANYD(beta[mfi]),
                     BL_TO_FORTRAN_ANYD(alpha[mfi]),
                     alpha0, M, R);
        }
    }

    // compute the time step
    const Real* dx = geom.CellSize();
    Real dt = 0.9*dx[0] / (2.0*BL_SPACEDIM);

    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
#if (BL_SPACEDIM == 2)
        WriteSingleLevelPlotfile(pltfile, data_new, {"D", "Sx", "Sy", "tau"}, geom, time, 0);
#else
        WriteSingleLevelPlotfile(pltfile, data_new, {"D", "Sx", "Sy", "Sz", "tau"}, geom, time, 0);
#endif
    }

    for (int n = 1; n <= nsteps; ++n)
    {
        MultiFab::Copy(data_old, data_new, 0, 0, Ncomp, 0);

        // new_data = old_data + dt * (something)
        advance(data_old, data_new, dt, geom, gamma, gamma_up, beta, alpha, alpha0, M, R);
        time = time + dt;

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && n%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,5);
#if (BL_SPACEDIM == 2)
            WriteSingleLevelPlotfile(pltfile, data_new, {"D", "Sx", "Sy", "tau"}, geom, time, n);
#else
            WriteSingleLevelPlotfile(pltfile, data_new, {"D", "Sx", "Sy", "Sz", "tau"}, geom, time, n);
#endif
        }
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
