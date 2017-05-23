#ifndef MYFUNC_F_H__
#define MUFUNC_F_H__

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void init_data(amrex_real* data, const int* dlo, const int* dhi,
                  const int* lo, const int* hi,
                  const int* Ncomp,
                  const amrex_real* dx, const amrex_real* prob_lo, const amrex_real* prob_hi);

    void compute_flux(amrex_real* data, const int* philo, const int* phihi,
                      const int* lo, const int* hi,
                      const int* Ncomp,
                      amrex_real* fluxx, const int* fxlo, const int* fxhi,
                      amrex_real* fluxy, const int* fylo, const int* fyhi,
#if (BL_SPACEDIM == 3)
                      amrex_real* fluxz, const int* fzlo, const int* fzhi,
#endif
                      const amrex_real* dx, amrex_real* g, amrex_real* dt);

    void update_data(const int* lo, const int* hi, const int* Ncomp,
                    amrex_real* pold, const int* polo, const int* pohi,
                    amrex_real* pnew, const int* pnlo, const int* pnhi,
                    amrex_real* fluxx, const int* fxlo, const int* fxhi,
                    amrex_real* fluxy, const int* fylo, const int* fyhi,
#if (BL_SPACEDIM == 3)
                    amrex_real* fluxz, const int* fzlo, const int* fzhi,
#endif
                    const amrex_real* dx, amrex_real dt); // note dt is passed by value

#ifdef __cplusplus
}
#endif

#endif
