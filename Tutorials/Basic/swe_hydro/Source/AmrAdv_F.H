#ifndef _Adv_F_H_
#define _Adv_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void initdata(const int& level, const amrex_real& time,
		  const int* lo, const int* hi,
		  BL_FORT_FAB_ARG_3D(state),
		  const amrex_real* dx, const amrex_real* problo,
          const int* Ncomp, const amrex_real* alpha0,
          const amrex_real* M, const amrex_real* R,
          const amrex_real* p, const int* nlayers,
          const amrex_real* gamma, const amrex_real* rho_swe);

    void get_face_velocity(const int& level, const amrex_real& time,
               const BL_FORT_FAB_ARG_3D(state),
			   AMREX_D_DECL(BL_FORT_FAB_ARG(xvel),
				  BL_FORT_FAB_ARG(yvel),
				  BL_FORT_FAB_ARG(zvel)),
			   const int* Ncomp);

    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
		     const BL_FORT_FAB_ARG_3D(state),
		     const int* tagval, const int* clearval,
		     const int* lo, const int* hi,
		     const amrex_real* dx, const amrex_real* problo,
		     const amrex_real* time, const amrex_real* phierr, const int* Ncomp);

    void advect(const amrex_real& time, const int* lo, const int*hi,
		const BL_FORT_FAB_ARG_3D(statein),
		BL_FORT_FAB_ARG_3D(stateout),
		AMREX_D_DECL(BL_FORT_FAB_ARG_3D(fx),
		       BL_FORT_FAB_ARG_3D(fy),
		       BL_FORT_FAB_ARG_3D(fz)),
		const amrex_real* dx, const amrex_real& dt, const int* Ncomp,
        const bool* gr, const amrex_real* alpha0, const amrex_real* M, const amrex_real* R,const amrex_real* prob_lo);

    void comp_from_swe(BL_FORT_FAB_ARG_3D(U_comp),
        const BL_FORT_FAB_ARG_3D(U_swe),
        amrex_real* p, amrex_real* rho,
        const int* lo, const int* hi,
        const int* n_cons_comp, const int* n_swe_comp,
        amrex_real* gamma,
        const amrex_real* dx, const amrex_real* alpha0,
        const amrex_real* M, const amrex_real* R, const int * nghost,
        const amrex_real* prob_lo);

    void swe_from_comp(BL_FORT_FAB_ARG_3D(U_prim),
        const BL_FORT_FAB_ARG_3D(U_swe),
        BL_FORT_FAB_ARG_3D(p_comp),
        amrex_real* p_swe,
        const int* lo, const int* hi,
        const int* n_cons_comp, const int* n_swe_comp,
        const amrex_real* alpha0, const amrex_real* M,
        const amrex_real* R, const amrex_real* dx,
        const amrex_real* prob_lo);

    void cons_to_prim(const BL_FORT_FAB_ARG_3D(U),
        BL_FORT_FAB_ARG_3D(U_prim),
        BL_FORT_FAB_ARG_3D(p),
        const int* lo, const int* hi, const int* Ncomp,
        const amrex_real* gamma, const amrex_real* alpha0,
        const amrex_real* M,
        const amrex_real* R, const amrex_real* dx,const amrex_real* prob_lo);

#ifdef __cplusplus
}
#endif

#endif
