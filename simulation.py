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

    "glc_ext" : 1.0,
    "glc_cyt" : 1.0,
    "atp_cyt" : 1.0,
    "adp_cyt" : 1.0,
    "amp_cyt" : 1.0,

    "atp_mit": 1.0,
    "adp_mit": 1.0,
    "amp_mit": 1.0,
    "pi_cyt" : 1.0,
    "pi_mit" : 1.0,
    "gtp_mit" : 1.0,
    "gdp_mit" : 1.0,

    "glc6p"   : 1.0,
    "fru6p"   : 1.0,
    "fru26p"  : 1.0,
    "fru16bp" : 1.0,
    "fru16p"  : 1.0,

    "grap"    : 1.0,
    "dhap"    : 1.0,
    "bpg13"   : 1.0,
    "pg3"     : 1.0,
    "pg2"     : 1.0,
    "pep"     : 1.0,
    "pyr_cyt" : 1.0,
    "pyr_mit" : 1.0,
    "lac"     : 1.0,
    "lac_ext" : 1.0,
    "cr"      : 1.0,
    "crp"     : 1.0,

    "mal_cyt" : 1.0,
    "oa_cyt"  : 1.0,
    "mal_mit":  1.0,
    "oa_mit":   1.0,
    "asp_cyt" : 1.0,
    "asp_mit" : 1.0,
    "akg_cyt" : 1.0,
    "akg_mit" : 1.0,
    "glu_cyt" : 1.0,
    "glu_mit" : 1.0,
    "suc"     : 1.0,
    "fum"     : 1.0,

    "nad_cyt" : 1.0,
    "nadh_cyt": 1.0,
    "nad_mit" : 1.0,
    "nadh_mit": 1.0,
    "fad"     : 1.0,
    "fadh2"   : 1.0,

    "dhap_cyt": 1.0,
    "g3p_cyt" : 1.0,

    "Q"       : 1.0,
    "QH2"     : 1.0,
    "cytc_ox" : 1.0,
    "cytc_red" : 1.0,
    "O2_mit"  : 1.0,

    "K_mit"   : 1.0,
    "K_cyt"   : 1.0,
    "Na_mit"  : 1.0,
    "Na_cyt"  : 1.0,
    "H+_cyt"  : 1.0,
    "H+_mit"  : 1.0,
    "Ca_cyt"  : 1.0,
    "Ca_mit"  : 1.0,

    "CoA"     : 1.0,
    "ACoA"    : 1.0,
    "sucCoA" : 1.0,
    "fad_pdhg" : 1.0,
    "fadh2_pdhg" : 1.0,

    "citr" : 1.0,
    "isocitr" : 1.0,

    "mal_dehydr" : {
        "Vmax" : 3.2 * 10**4,
        "Keq"  : 0.0001,
        "Km_nad" : 0.06,
        "Km_mal" : 0.145,
        "Km_oa"  : 0.017,
        "Km_nadh" : 0.044,
    },

    "fumarase" : {
        "Vmax" : 6.4 * 10**7,
        "Keq"  : 4.4,
        "Km_fum" : 0.14,
        "Km_mal" : 0.3,
    },

    "suc_dehydr": {
        "Vmax_succdh" : 1.6*10**5,
        "Vmax_nadh" : 10**12,
        "Km_suc"    : 1.6,
        "Ki_mal"    : 2.2,
        "Km_nad"    : 1.0, # !!!!!!! значение свято от балды
        "Em_FAD-succdh" : 100,
    },

    "sucCoAsyntase" : {
        "Vmax" : 1.92 * 10**4,
        "Keq"  : 3.8,
        "Amax_P" : 1.2,
        "Km_P"   : 2.5, # 0.72 другое значение указвнное в статье !!!!!!!!!!!!
        "n_P"    : 3,
        "Km_sucCoA" : 0.041,
        "Km_adp"   : 0.25,
        "Km_suc"   : 1.6,
        "Km_CoA"   : 0.056,
        "Km_atp"   : 0.017,
        "Km_sucCoA_G" : 0.086,
        "Km_gdp"  : 0.007,
        "Km_gtp"  : 0.036,
        "Km_suc_G" : 0.49,
        "Km_CoA_G" : 0.036,

    },

    "akg_dehydr" : {
        "Vmax_nad" : 134.4,
        "Vmax_fad" : 1e4,

        "Km1" : 2.5,
        "Km2" : 2.5, #  !!!!! величина взята от балды !!!!!
        "Ki_ca" : 1, # !!!!! величина взята от балды !!!!!


        "Km_nad" : 0.021,
        "Ki_nadh" : 0.0045,
        "Km_CoA"  : 0.013,
        "Km_SucCoA" : 0.0045,
        "Km_fad" : 0.00001,

        "Em_fad" : 297,
        "Em_nad"  : 297, # !!!!! величина взята от балды !!!!!

    },

    "isocit_dehydr" : {
        "Vmax" : 64,
        "n_isocit" : 1.9,
        "Km1_isocit" : 0.11,
        "Km2_isocit" : 0.06,
        "Ka_Ca" : 0.0074,
        "n_Ca"  : 2,
        "Km_nad" : 0.091,
        "Ki_nadh" : 0.041,
    },

    "aconitase" : {
        "Vmax" : 1.6 * 10**6,
        "Keq"  : 0.067,
        "Km_cit" : 0.48,
        "Km_isocit" : 0.12,
    },

    "citrate_syntase" : {
        "Vmax" : 1.28 * 10**3,
        "Km_oxa" : 0.0045,
        "Ki_cit" : 3.7,
        "Km_accoa" : 0.005,
        "Ki_CoA"  : 0.025,
    },

    "pyr_dehyd_comp" : {
        "Vmax_pdhc_fad" : 13.1,
        "Vmax_pdhc_nad" : 1e4, # !!!!!!
        "Amax_Ca"  : 1.7,
        "Ka_Ca"    : 10**-3,
        "Km_pyr"   : 0.068,
        "Km_nad"   : 0.041,
        "Km_CoA"   : 0.0047,
        "Ki_AcoA"  : 0.0004,
        "Km_fad"   : 0.00001,
        "Em_fad"   : 297,
        "Em_nad"   : 300, # !!!!!!! нет данных
    },

    "pyr_exchanger" : {
        "Vmax" : 128,
        "Km_pyr_cit" : 0.15,
        "Km_pyr_mit" : 0.15,

    },

    "complex4" : {
        "Vmax" : 32.5,
        "Km_O2" : 0.001,
        "Km_cytc" : 0.001,
        "n"       : 2,
        "dGh"     : 1,

    },

    "complex3" : {
        "Vmax" : 2.25*10**4,
        "n" : 2,
    },

    "complex1" : {
        "Vmax" : 2.25,
    },

    "ca_h_pump" : {
        "Km_ca" : 0.01,
        "n_H"   : 3,
    },

    "ca_na_pump" : {
        "Km_na" : 8,
        "Km_ca" : 0.0096,
        "n_Na"  : 3,
        "n"     : 2.8,
    },

    "calcium_ed" : {
        "P_RMC"    : 2*10**-6, # m/s
        "Ki_cacyt" : 0.0001,
        "P_Mcu"    : 2 * 10**-2, # m/s
        "Km_Mcu"   : 19.2,
        "n"        : 0.6,
        "Ka"       : 0.0003,
        "n_a"      : 5,
    },

    "phos_pump" : {
        "Vmax"  : 43.3485,
    },

    "Na_pump": {
        "Vmax": 5*10**-3,
    },

    "K_pump" : {
        "Vmax" : 7.5 * 10**-4,
    },

    "protons_ed" : {
        "P_H_mit" : 2 * 10**-4, # m/s
    },

    "sodium_ed" : {
        "P_Namit" : 10**-10 # m/s
    },

    "potassium_ed" : {
        "P_Kmit" : 2 * 10**-10, # m/s !!!!!
    },

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
        "Am"    : 3.7 * 10**-5, # cm**2
        "Cmit"  : 0.9 * 10**-6, #F/cm**2

        "Em_N"  : 1,     # !!!!! нет данных
        "Em_Q"  : 1,     # !!!!! нет данных
        "Em_cytc" : 1,  # !!!!! нет данных

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

class Simulator():

    def __init__(self, arg):

        self.arg = arg

    def model_equations(self, t, y):
        arg = self.arg

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

        potassium_current_ed = lib.getIed(arg, ion="K")
        sodium_current_ed = lib.getIed(arg, ion="Na")
        protons_current_ed = lib.getIed(arg, ion="H+")

        phos_pump = lib.getVpumps(arg, ion="pi")

        calcium_ed = lib.getIca_ed(arg)
        ca_na_pump = lib.getIca_na_pump(arg)
        ca_h_pump = lib.getIca_h_pump(arg)

        complex1 = lib.getVcomplex1(arg)
        complex3 = lib.getVcomplex3(arg)
        complex4 = lib.getVcomplex4(arg)

        pyr_exchanger = lib.getVpyr_exchanger(arg)
        pyr_dehyd_compACoA, pyr_dehyd_compFad = lib.getVpyr_dehydrogenase_complex(arg)
        citrate_syntase = lib.getVcitratesyntase(arg)
        aconitase = lib.getVaconitase(arg)
        isocit_dehydr = lib.getVisocit_dehydrogenase(arg)
        akg_dehydr_fad, akg_dehydr_nad  = lib.getVakg_dehydrogenase(arg)

        sucCoAsyntase_atp = lib.getVsucCoAsyntase(arg, mode="atp")
        sucCoAsyntase_gtp = lib.getVsucCoAsyntase(arg, mode="gtp")

        v_succdh_fad, v_succdh = lib.getVsuc_dehydrydrogenase(arg)
        fumarase = lib.getVfumarase(arg)

        mal_dehydr = lib.getVmal_dehydr(arg)


        # calculate balans of currents
        # update mitochondrial potential


        return []

simulalor = Simulator(arg)



simulalor.model_equations(0, 0)

