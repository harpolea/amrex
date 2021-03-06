#ifndef AMREX_EBAMRUTIL_H_
#define AMREX_EBAMRUTIL_H_

#include <AMReX_TagBox.H>
#include <AMReX_MultiFab.H>

namespace amrex {

    void TagCutCells (TagBoxArray& tags, const MultiFab& state);

}

#endif
