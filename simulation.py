import numpy as np
from scipy.integrate import ode
import lib

y0, t0 = [1.0], 0

# Notations
# glc is glucose
# atp ia ATP
# glc6p is glucose-6-P
# fru6p is fructose-6-P


arg = {

    "glc_ext" : 0,
    "glc_cyt" : 0,
    "atp_cyt" : 0,
    "adp_cyt" : 0,
    "amp_cyt" : 0,

    "atp_mit": 0,
    "adp_mit": 0,
    "amp_mit": 0,
    "pi_cyt" : 0,
    "pi_mit" : 0,


    "glc6p"   : 0,
    "fru6p"   : 0,
    "fru26p"  : 0,
    "fru16bp" : 0,
    "fru16p"  : 0,

    "grap"    : 0,
    "dhap"    : 0,
    "bpg13"   : 0,
    "pg3"     : 0,
    "pg2"     : 0,
    "pep"     : 0,
    "pyr_cyt" : 0,
    "lac"     : 0,
    "lac_ext" : 0,
    "cr"      : 0,
    "crp"     : 0,

    "mal_cyt" : 0,
    "oa_cyt"  : 0,
    "mal_mit":  0,
    "oa_mit":   0,
    "asp_cyt" : 0,
    "asp_mit" : 0,
    "akg_cyt" : 0,
    "akg_mit" : 0,
    "glu_cyt" : 0,
    "glu_mit" : 0,

    "nad_cyt" : 0,
    "nadh_cyt": 0,
    "nad_mit" : 0,
    "nadh_mit": 0,
    "fad"     : 0,
    "fadh2"   : 0,

    "dhap_cyt": 0,
    "g3p_cyt" : 0,

    "Q"       : 0,
    "QH2"     : 0,

    "atp_consumption" : {
        "Vmax" : 1,
        "Km_atp" : 1,
        "activation" : 0,
    },

    "atp/atp_axchanger" : {
        "Vmax" : 5.4 * 10**-5,
        "S_Vmm" : 0.3,

    },

    "atp_syntase" : {

        "dG0" : 30500,
        "n" : 3,
        "k" : 3,
    },

    "mitgly3pdehyd": {
        "Vmax_g3pdh": 6.4 * 10**4,
        "Vmax_Q" : 3.2 * 10**6,
        "Em_dhap_g3p" : 190, # mV
        "Em_FAD_g3p"  : 210, # mV
        "Em_Q"        : 200, # mV !!!!!!!! нет данных, значние взято от балды

        "Km_dhap": 0.17, # сзято от цитозольного фермента
        "Km_g3p": 0.3,   # сзято от цитозольного фермента
    },

    "cytgly3pdehyd": {
        "Vmax": 3.2 * 10**4,
        "Keq": 3257.3,
        "Km_dhap": 0.17,
        "Km_nadh": 0.01,
        "Km_g3p": 0.3,
        "Km_nad": 0.03,
    },

    "mal_akg_carrier": {
        "Vmax" : 32,
        "Km_mal_cyt" : 1.36,
        "Km_akg_cyt" : 0.1,
        "Km_mal_mit" : 0.71,
        "Km_akg_mit" : 0.2,
    },

    "asp_glu_carrier": {
        "Vmax" : 3200,
        "Km_asp_mit" : 0.05,
        "Km_glu_cyt" : 2.8,
        "Km_asp_cyt" : 0.05,
        "Km_glu_mit" : 2.8,
    },

    "mito_membrane" : {
        "Vmm" : -200,
        "h_cyt" : 1, # !!!!!
        "h_mit" : 1, # !!!!!

    },

    "asp_aminotrans" : {
        "Vmax" : 32,
        "Keq"  : 0.147,

    },

    "malatdehyd" : {
        "Vmax" : 10**4,
        "Keq"  : 10**-4,
        "Km_nad" : 0.05,
        "Km_mal" : 0.77,
        "Km_oa"  : 0.04,
        "Km_nadh" : 0.05,
    },

    "creatinekinase" : {
        "Vmax" : 0.0135,
        "Keq"  : 7.0
    },

    "MCT" : {
        "Vmax" : 5.0,
        "Keq_mct" : 1.737,
        "Km_lac" : 1.1,
        "Km_lac_ex" : 1.1,
    },

    "LDG" : {
        "Vmax" : 10**5,
        "Keq_ldg" : 8400,
        "Km_pyr"  : 0.36,
        "Km_nadh" : 0.043,
        "Km_lac"  : 4.2,
        "Km_nad"  : 0.088,
    },

    "pyruvatekinase": {
        "Vmax" : 23.76,
        "Km_pep" : 0.074,
        "Km_adp" : 0.42,
        "Ki_atp" : 4.4,
    },

    "enolase": {
        "Vmax" : 216000,
        "Keq_eno" : 0.5,
        "Km_pg2"  : 0.05,
        "Km_pep"  : 0.15,
    },

    "p-gricerate_mutase": {
        "Vmax" : 14400,
        "Keq_pgm" : 0.1814,
        "Km_pg3"  : 0.22,
        "Km_pg2"  : 0.28,

    },

    "p-glyceratekinase": {
        "Vmax" : 396,
        "Keq_pgk" : 1310,
        "Km_bpg13" : 0.063,
        "Km_adp"   : 0.42,
        "Km_pg3"   : 0.67,
        "Km_atp"   : 0.25,
    },

    "grap_dehydr" : {
        "Vmax" : 72000,
        "Keq_gapdh" : 0.0868,
        "Km_nad"    : 0.01, # 0.027
        "Km_grap"   : 0.101,
        "Km_pi"     : 3.9,
        "Km_nadh"   : 0.008,
        "Km_bpg13"  : 0.0035,


    },

    "triosep-isomerase": {
        "Vmax" : 10**6,
        "Keq_tri": 0.0545,
        "Km_dhap" : 0.84,
        "Km_grap" : 1.65,
    },

    "aldolase" : {
        "Vmax" : 46.8,
        "Keq_aldo" : 0.0976,
        "Km_fru16p" : 0.003,
        "Km_grap"   : 0.08,
        "Km_dhap"   : 0.03,

    },

    "fru-2,6-bisphosphatase" : {
        "Vmax" : 0.052,
        "Km_fru26p" : 0.07,
        "Ki_fru6p" : 0.02,

    },

    "phosphofructokinase2": {
        "Vmax" : 0.0026,
        "Km_fru6p" : 0.027,
        "Km_atp"   : 0.055,
        "Ka_amp"   : 0.073,
        "Ka_adp"   : 0.056,

    },

    "fru-1,6-bisphosphatase" : {
        "Vmax" : 0.455,
        "Km_fru6p" : 0.132,
    },

    "phosphofructokinase1": {
        "Vmax" : 49.6,
        "Km_fru6p" : 0.111,
        "Km_atp"   : 0.04,
        "n" : 1.8,
        "Ki_atp" : 1.2,
        "K0" : 0.55,
        "Ka_fru26p" : 0.0042,
        "nfru26p" : 5.5,
        "Ka_fru26p" : 0.005,

    },

    "glc6p_isomerase" : {
        "Vmax" : 24.4,
        "Keq_g6piso": 0.5157,
        "Km_glc6p" : 0.593,
        "Km_fru6p": 0.095,
    },

    "hexokinase" : {
        "Vmax" : 9.36,
        "Km_glc" : 0.043,
        "Km_atp" : 0.37,
        "Kiatp_glc6p" : 0.074,
        "Ki_glc6p" : 0.1,


    },

    "glc_trs" : {
        "Vmax" : 0.72,
        "Km_glc_cyt" : 2.87,
        "Km_glc_ext" : 2.87,
    }


}


def model_equations(t, y, arg):

    # Glycolysis

    # Glucose transporter
    vglc_transp = lib.getVglucosetransporter(arg)
    # Hexokinase
    vhexokinase = lib.getVhexokinase(arg)
    glucose6p_isomerase = lib.getVglucose6p_isomerase(arg)
    phosphofructokinase1 = lib.getVphosphofructokinase1(arg)
    fru16bisphosphatase = lib.getVfru16bisphosphatase(arg)
    aldolase = lib.getValdolase(arg)
    triosep_isomerase = lib.getVtriosep_isomerase(arg)
    grap_dehydr = lib.getVgrap_dehydrogenase(arg)
    ph_glyceratekinase = lib.getVphosphoglyceratekinase(arg)
    ph_glyceratemutase = lib.getVphosphoglyceratemutase(arg)
    enolase = lib.getVenolase(arg)
    pyruvatekinase = lib.getVpyruvatekinase(arg)
    ldg = lib.getVlactatedehydrogenase(arg)
    # Monocarboxilate transporter
    mct = lib.getVmonocarboxilatetransporter(arg)
    creatinekinase = lib.getVcreatinekinase(arg)
    cyt_malatdehydrogenase = lib.getVmalatdehydrogenase(arg, mode="cyt")
    mito_malatdehydrogenase = lib.getVmalatdehydrogenase(arg, mode="mit")
    cyt_asp_aminotrans = lib.getVaspartateaminotransferase(arg, mode="mit")
    mito_asp_aminotrans = lib.getVaspartateaminotransferase(arg, mode="mit")

    asp_glu_carrier = lib.getVasp_glu_carrier(arg)
    mal_akg_carrier = lib.getVmal_akg_carrier(arg)
    cytg3pdehyd = lib.getVcytg3pdehyd(arg)

    atp_syntase = lib.getVatp_syntase(arg)

    mitg3pdehydFADH2, mitg3pdehydQH2 = lib.getVmitg3pdehyd(arg)
    atp_atp_axchanger = lib.getVatp_atp_axchanger(arg)

    atp_consumption = lib.getVatp_consumption(arg)



    return [vglc_transp, vhexokinase, glucose6p_isomerase]

model_equations(0, 0, arg)

# r = ode(model_equations).set_integrator('dopri5', method='bdf')
#
# r.set_initial_value(y0, t0).set_f_params([0.2])  #.set_jac_params([-0.2])
# t1 = 10
# dt = 0.1
# while r.successful() and r.t < t1:
#     # r.integrate(r.t + dt)
#     print(r.t+dt, r.integrate(r.t+dt))