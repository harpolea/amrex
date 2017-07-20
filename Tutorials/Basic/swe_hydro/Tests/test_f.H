#ifndef _Test_F_H
#define _Test_F_H

#ifdef __cplusplus
extern "C"
{
#endif

    void test_cons_to_prim(bool * passed);
    void test_rhoh_from_p(bool * passed);
    void test_p_from_rhoh(bool * passed);
    void test_p_from_rho_eps(bool * passed);
    void test_calc_gamma_up(bool * passed);
    void test_calc_gamma_down(bool * passed);
    void test_gr_sources(bool * passed);
    void test_W_swe(bool * passed);
    void test_swe_from_comp(bool * passed);
    void test_calc_gamma_swe(bool * passed);
    void test_comp_from_swe(bool * passed);
    void test_gr_swe_flux(bool * passed);
    void test_f_of_p(bool * passed);
    void test_zbrent(bool * passed);
    void test_gr_comp_flux(bool * passed);

#ifdef __cplusplus
}
#endif

#endif
