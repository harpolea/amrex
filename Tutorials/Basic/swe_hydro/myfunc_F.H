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
                  const amrex_real* dx, const amrex_real* prob_lo,
                  const amrex_real* prob_hi,
                  amrex_real* gamma_up, const int* glo, const int* ghi,
                  amrex_real* beta, const int* blo, const int* bhi,
                  amrex_real* alpha, const int* alo, const int* ahi,
                  amrex_real alpha0, amrex_real M, amrex_real R);

    void update_data(const int* lo, const int* hi, const int* Ncomp,
                    amrex_real* pold, const int* polo, const int* pohi,
                    amrex_real* pnew, const int* pnlo, const int* pnhi,
                    const amrex_real* dx, amrex_real dt,
                    const amrex_real* gamma,
                    amrex_real* gamma_up, const int* glo, const int* ghi,
                    amrex_real* beta, const int* blo, const int* bhi,
                    amrex_real* alpha, const int* alo, const int* ahi,
                    amrex_real alpha0, amrex_real M, amrex_real R); // note dt is passed by value

#ifdef __cplusplus
}
#endif

#endif
