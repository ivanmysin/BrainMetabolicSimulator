import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
import lib

# Notations
metabolites = []
metabolites.append({"idx" : 0, "full" : "Extracellular glucose", "short" : "glc_ext", "rest" : 2.48 })
metabolites.append({"idx" : 1, "full" : "Cytosolic glucose", "short" : "glc", "rest" : 1.2})
metabolites.append({"idx" : 2, "full" : "Cytosolic ATP", "short" : "atp_cyt", "rest" : 2.2})
metabolites.append({"idx" : 3, "full" : "Cytosolic ADP", "short" : "adp_cyt", "rest" : 0.01  })
metabolites.append({"idx" : 4, "full" : "Cytosolic AMP", "short" : "amp_cyt", "rest" : 0.002})
metabolites.append({"idx" : 5, "full" : "Mitochondrial ATP", "short" : "atp_mit", "rest" : 2.2 }) # такая же как в цит
metabolites.append({"idx" : 6, "full" : "Mitochondrial ADP", "short" : "adp_mit", "rest" : 0.01}) # такая же как в цит
metabolites.append({"idx" : 7, "full" : "Mitochondrial AMP", "short" : "amp_mit", "rest" : 0.002}) # такая же как в цит
metabolites.append({"idx" : 8, "full" : "Cytosolic inorganic phosphate", "short" : "pi_cyt", "rest" : 0.5}) # !!!!!!!
metabolites.append({"idx" : 9, "full" : "Mitochondrial inorganic phosphate", "short" : "pi_mit", "rest" : 0.5}) # !!!!!!!
metabolites.append({"idx" : 10, "full" : "Mitochondrial GTP", "short" : "gtp_mit", "rest" : 1.0}) # !!!!!!!
metabolites.append({"idx" : 11, "full" : "Mitochondrial GDP", "short" : "gdp_mit", "rest" : 0.01}) # !!!!!!!
metabolites.append({"idx" : 12, "full" : "Glucose-6-phosphate", "short" : "glc6p", "rest" : 0.1})
metabolites.append({"idx" : 13, "full" : "Fructose-6-phosphate", "short" : "fru6p", "rest" : 0.03})
metabolites.append({"idx" : 14, "full" : "Fructose-2,6-bisphosphate", "short" : "fru26p", "rest" : 0.04}) # !!!!!!!
metabolites.append({"idx" : 15, "full" : "Fructose-1,6-bisphosphate", "short" : "fru16p", "rest" : 0.04})
metabolites.append({"idx" : 16, "full" : "Glycerol 3-phosphate", "short" : "grap", "rest" : 0.01})
metabolites.append({"idx" : 17, "full" : "Dihydroxyacetone phosphate", "short" : "dhap", "rest" : 0.05})
metabolites.append({"idx" : 18, "full" : "1,3-Bisphosphoglycerate", "short" : "bpg13", "rest" : 0.05}) # !!!!!!!!!!
metabolites.append({"idx" : 19, "full" : "3-Phosphoglycerate", "short" : "pg3", "rest" : 0.1})
metabolites.append({"idx" : 20, "full" : "Phosphoenolpyruvate", "short" : "pep", "rest" : 0.01})
metabolites.append({"idx" : 21, "full" : "Cytosolic pyruvate", "short" : "pyr_cyt", "rest" : 0.15})
metabolites.append({"idx" : 22, "full" : "Mitochondrial pyruvate", "short" : "pyr_mit", "rest" : 0.04})
metabolites.append({"idx" : 23, "full" : "Cytosolic lactate", "short" : "lac", "rest" : 0.6})
metabolites.append({"idx" : 24, "full" : "Extracellular lactate", "short" : "lac_ext", "rest" : 0.6})
metabolites.append({"idx" : 25, "full" : "Creatine", "short" : "cr", "rest" : 0.01 }) # !!!!!!!!
metabolites.append({"idx" : 26, "full" : "Creatine phosphate", "short" : "crp", "rest" : 4.9})
metabolites.append({"idx" : 27, "full" : "Cytosolic malate", "short" : "mal_cyt", "rest" : 2.0})  # !!!!!!
metabolites.append({"idx" : 28, "full" : "Cytosolic oxaloacetate", "short" : "oa_cyt", "rest" : 0.01}) # !!!!!!
metabolites.append({"idx" : 29, "full" : "Mitochondrial malate", "short" : "mal_mit", "rest" : 2.0})
metabolites.append({"idx" : 30, "full" : "Mitochondrial oxaloacetate", "short" : "oa_mit", "rest" : 0.01})
metabolites.append({"idx" : 31, "full" : "Cytosolic aspartate", "short" : "asp_cyt", "rest" : 1.0}) # !!!!!!!!!!
metabolites.append({"idx" : 32, "full" : "Mitochondrial aspartate", "short" : "asp_mit", "rest" : 1.0}) # !!!!!!!!!!
metabolites.append({"idx" : 33, "full" : "Cytosolic alpha-ketoglutarate", "short" : "akg_cyt", "rest" : 1.0}) # !!!!!!!!!!
metabolites.append({"idx" : 34, "full" : "Mitochondrial alpha-ketoglutarate", "short" : "akg_mit", "rest" : 0.05})
metabolites.append({"idx" : 35, "full" : "Cytosolic glutamate", "short" : "glu_cyt", "rest" : 10.0})
metabolites.append({"idx" : 36, "full" : "Mitochondrial glutamate", "short" : "glu_mit", "rest" : 10.0})
metabolites.append({"idx" : 37, "full" : "Succinate", "short" : "suc", "rest" : 1.2})
metabolites.append({"idx" : 38, "full" : "Fumarate", "short" : "fum", "rest" : 0.05})
metabolites.append({"idx" : 39, "full" : "Cytosolic NAD+", "short" : "nad_cyt", "rest" : 0.001}) # !!!!!!!!!!
metabolites.append({"idx" : 40, "full" : "Cytosolic NADH", "short" : "nadh_cyt", "rest" : 0.006})
metabolites.append({"idx" : 41, "full" : "Mitochondrial NAD+", "short" : "nad_mit", "rest" : 0.01}) # !!!!!!!!!!
metabolites.append({"idx" : 42, "full" : "Mitochondrial NADH", "short" : "nadh_mit", "rest" : 0.12})
metabolites.append({"idx" : 43, "full" : "FAD of succinate dehydrogenase", "short" : "fad_sucdh", "rest" : 0.001}) # !!!!!!!!!!
metabolites.append({"idx" : 44, "full" : "FADH2 of succinate dehydrogenase", "short" : "fadh2_sucdh", "rest" : 0.001}) # !!!!!!!!!!
metabolites.append({"idx" : 45, "full" : "FAD of glycerol-3-phosphate dehydrogenase", "short" : "fad_g3dh", "rest" : 0.001})# !!!!!!!!!!
metabolites.append({"idx" : 46, "full" : "FADH2 of glycerol-3-phosphate dehydrogenase", "short" : "fadh2_g3dh", "rest" : 0.001})# !!!!!!!!!!
metabolites.append({"idx" : 47, "full" : "Coenzyme Q oxidized", "short" : "Q", "rest" : 0.01}) # !!!!!!!!!!
metabolites.append({"idx" : 48, "full" : "Coenzyme QH2 reduced", "short" : "QH2", "rest" : 0.01})# !!!!!!!!!!
metabolites.append({"idx" : 49, "full" : "Cytochrome c oxidized", "short" : "cytc_ox", "rest" : 0.01}) # !!!!!!!!!!
metabolites.append({"idx" : 50, "full" : "Cytochrome c reduced", "short" : "cytc_red", "rest" : 0.01}) # !!!!!!!!!!
metabolites.append({"idx" : 51, "full" : "Mitochondrial oxigen (O2)", "short" : "o2_mit", "rest" : 0.028})
metabolites.append({"idx" : 52, "full" : "Mitochondrial potassium (K)", "short" : "k_mit", "rest" : 10.0}) # !!!!!!!!!!
metabolites.append({"idx" : 53, "full" : "Cytosolic potassium (K)", "short" : "k_cyt", "rest" : 140.0})
metabolites.append({"idx" : 54, "full" : "Mitochondrial sodium (Na)", "short" : "na_mit", "rest" : 7.0}) # !!!!!!!!!!
metabolites.append({"idx" : 55, "full" : "Cytosolic sodium (Na)", "short" : "na_cyt", "rest" : 7.0})
metabolites.append({"idx" : 56, "full" : "Cytosolic proton (H+)", "short" : "h+_cyt", "rest" : 10**-7}) # !!!!!!!!!!
metabolites.append({"idx" : 57, "full" : "Mitochondrial proton (H+)", "short" : "h+_mit", "rest" : 10**-7}) # !!!!!!!!!!
metabolites.append({"idx" : 58, "full" : "Cytosolic calcium (Ca)", "short" : "ca_cyt", "rest" : 10**-7})
metabolites.append({"idx" : 59, "full" : "Mitochondrial calcium (Ca)", "short" : "ca_mit", "rest" : 10**-7})
metabolites.append({"idx" : 60, "full" : "CoA", "short" : "coa", "rest" : 0.02})
metabolites.append({"idx" : 61, "full" : "Acetyl-CoA", "short" : "acoa", "rest" : 1.0})
metabolites.append({"idx" : 62, "full" : "Succinyl-CoA", "short" : "succoa", "rest" : 0.001})
metabolites.append({"idx" : 63, "full" : "Citrate", "short" : "citr", "rest" : 1.3})
metabolites.append({"idx" : 64, "full" : "Isocitrate", "short" : "isocitr", "rest" : 0.02})
metabolites.append({"idx" : 65, "full" : "Voltage on mitochondrial membrane", "short" : "Vmm", "rest" : -135.0})
metabolites.append({"idx" : 66, "full" : "FAD of pyruvate dehydrogenase complex", "short" : "fad_pdhc", "rest" : 0.001}) # !!!!!!!!!!
metabolites.append({"idx" : 67, "full" : "FADH2 of pyruvate dehydrogenase complex", "short" : "fadh2_pdhc", "rest" : 0.001}) # !!!!!!!!!!
metabolites.append({"idx" : 68, "full" : "FAD of alpha-ketoglutarate dehydrogenase complex", "short" : "fad_akgdhc", "rest" : 0.001}) # !!!!!!!!!!
metabolites.append({"idx" : 69, "full" : "FADH2 of alpha-ketoglutarate dehydrogenase complex", "short" : "fadh2_akgdhc", "rest" : 0.001}) # !!!!!!!!!!
metabolites.append({"idx" : 70, "full" : "CO2", "short" : "co2", "rest" : 0.0}) # !!!!!!!!!!
metabolites.append({"idx" : 71, "full" : "2-Phosphoglycerate", "short" : "pg2", "rest" : 0.02})
metabolites.append({"idx" : 72, "full" : "Extracellular pyruvate", "short" : "pyr_ext", "rest" : 0.5})

global_params = {
    "Cmm" : 0.9 * 10**-6 * 3.7 * 10**-5, #farad,  capacity of mitochondrial membrane, 0.9 * 10**-6 F/cm^2, square 3.7 * 10**-5 cm^2
    "Volume_cyt_mit" : 0.07,             # part of mitochondrial volume in total volume of neuron, (Jolivet, 2015)
    "Volume_extracellular2cell" : 0.444, # ratio of extracellular volume to total volume of neurons, (Jolivet, 2015)
}



enzyme_params = {
    # "mal_dehydr" : {
    #     "Vmax" : 3.2 * 10**4,
    #     "Keq"  : 0.0001,
    #     "Km_nad" : 0.06,
    #     "Km_mal" : 0.145,
    #     "Km_oa"  : 0.017,
    #     "Km_nadh" : 0.044,
    # },

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
        "Em_FAD" : 100,
    },

    "sucCoAsyntase_4atp" : {
        "Vmax" : 1.92 * 10**4,
        "Keq"  : 3.8,
        "Amax_P" : 1.2,
        "Km_P"   : 2.5, # 0.72 другое значение указвнное в статье !!!!!!!!!!!!
        "n_P"    : 3,
        "Km_sucCoA" : 0.041,
        "Km_ndp"   : 0.25,
        "Km_suc"   : 1.6,
        "Km_CoA"   : 0.056,
        "Km_ntp"   : 0.017,
        # "Km_sucCoA_G" : 0.086,
        # "Km_gdp"  : 0.007,
        # "Km_gtp"  : 0.036,
        # "Km_suc_G" : 0.49,
        # "Km_CoA_G" : 0.036,

    },

    "sucCoAsyntase_4gtp": {
        "Vmax": 1.92 * 10 ** 4,
        "Keq": 3.8,
        "Amax_P": 1.2,
        "Km_P": 2.5,  # 0.72 другое значение указвнное в статье !!!!!!!!!!!!
        "n_P": 3,

        # "Km_sucCoA": 0.041,
        # "Km_adp": 0.25,
        # "Km_suc": 1.6,
        # "Km_CoA": 0.056,
        # "Km_atp": 0.017,
        "Km_sucCoA": 0.086,
        "Km_ndp": 0.007,
        "Km_ntp": 0.036,
        "Km_suc": 0.49,
        "Km_CoA": 0.036,
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
        "Km_pyr_cyt" : 0.15,
        "Km_pyr_mit" : 0.15,

        "Volume_cyt_mit" : global_params["Volume_cyt_mit"]

    },

    "complex4" : {
        "Vmax" : 32.5,
        "Km_O2" : 0.001,
        "Km_cytc" : 0.001,
        "n"       : 2,
        "dGh"     : 1,
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],
    },

    "complex3" : {
        "Vmax" : 2.25*10**4,
        "n" : 2,
        "Em_Q": 1,  # !!!!!!!!
        "Em_cytc" : 1, # !!!!!!!
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],
    },

    "complex1" : {
        "Vmax" : 2.25,
        "Em_N": 1,  # !!!!! нет данных
        "Em_Q": 1,  # !!!!! нет данных
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],
    },

    "ca_h_pump" : {
        "Km_ca" : 0.01,
        "n_H"   : 3,
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],
    },

    "ca_na_pump" : {
        "Km_na" : 8,
        "Km_ca" : 0.0096,
        "n_Na"  : 3,
        "n"     : 2.8,
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],
    },

    "calcium_ed" : {
        "P_RMC"    : 2*10**-6, # m/s
        "Ki_cacyt" : 0.0001,
        "P_Mcu"    : 2 * 10**-2, # m/s
        "Km_Mcu"   : 19.2,
        "n"        : 0.6,
        "Ka"       : 0.0003,
        "n_a"      : 5,
        "Am": 3.7 * 10 ** -5,  # cm**2
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],
    },

    "phos_pump" : {
        "Vmax"  : 43.3485,
        "is_simport" : True,
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],

    },

    "Na_pump": {
        "Vmax": 5*10**-3,
        "is_simport" : False,
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],

    },

    "K_pump" : {
        "Vmax" : 7.5 * 10**-4,
        "is_simport" : False,
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],

    },

    "protons_ed" : {
        "P" : 2 * 10**-4, # m/s
        "Am": 3.7 * 10 ** -5,  # cm**2
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],
    },

    "sodium_ed" : {
        "P" : 10**-10, # m/s
        "Am": 3.7 * 10 ** -5,  # cm**2
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],
    },

    "potassium_ed" : {
        "P" : 2 * 10**-10, # m/s !!!!!
        "Am": 3.7 * 10 ** -5,  # cm**2
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],
    },

    "atp_consumption" : {
        "Vmax" : 1,
        "Km_atp" : 1,
        "activation" : 0,
    },

    "atp/adp_axchanger" : {
        "Vmax" : 5.4 * 10**-5,
        "S_Vmm" : 0.3,
        "Volume_cyt_mit": global_params["Volume_cyt_mit"],
        "Cmm": global_params["Cmm"],

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
        "Keq" : 1.737,
        "Km_lac_cyt" : 1.1,
        "Km_lac_ext" : 1.1,
        "Volume_extracellular2cell" : global_params["Volume_extracellular2cell"],
    },

    "LDG" : {
        "Vmax" : 10**5,
        "Keq" : 8400,
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
        "Keq" : 0.5,
        "Km_pg2"  : 0.05,
        "Km_pep"  : 0.15,
    },

    "p-gricerate_mutase": {
        "Vmax" : 14400,
        "Keq" : 0.1814,
        "Km_pg3"  : 0.22,
        "Km_pg2"  : 0.28,

    },

    "p-glyceratekinase": {
        "Vmax" : 396,
        "Keq" : 1310,
        "Km_bpg13" : 0.063,
        "Km_adp"   : 0.42,
        "Km_pg3"   : 0.67,
        "Km_atp"   : 0.25,
    },

    "grap_dehydr" : {
        "Vmax" : 72000,
        "Keq" : 0.0868,
        "Km_nad"    : 0.01, # 0.027
        "Km_grap"   : 0.101,
        "Km_pi"     : 3.9,
        "Km_nadh"   : 0.008,
        "Km_bpg13"  : 0.0035,


    },

    "triosep-isomerase": {
        "Vmax" : 10**6,
        "Keq": 0.0545,
        "Km_dhap" : 0.84,
        "Km_grap" : 1.65,
    },

    "aldolase" : {
        "Vmax" : 46.8,
        "Keq" : 0.0976,
        "Km_fru16p" : 0.003,
        "Km_grap"   : 0.08,
        "Km_dhap"   : 0.03,

    },

    "fru-2,6-bisphosphatase" : {
        "Vmax" : 0.052,
        "Km" : 0.07,
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
        "Km" : 0.132,
    },

    "phosphofructokinase1": {
        "Vmax" : 49.6,
        "Km_fru6p" : 0.111,
        "Km_atp"   : 0.04,
        "n" : 1.8,
        "Ki_atp" : 1.2,
        "K0" : 0.55,
        "Ka_fru26p" : 0.0042,
        "n_fru26p" : 5.5,
        "Ka_fru26p" : 0.005,

    },

    "glc6p_isomerase" : {
        "Vmax" : 24.4,
        "Keq": 0.5157,
        "Km_glc6p" : 0.593,
        "Km_fru6p": 0.095,
    },

    "hexokinase" : {
        "Vmax" : 9.36,
        "Km_glc" : 0.043,
        "Km_atp" : 0.37,
        "Ki_atp" : 0.074,
        "Ki_glc6p" : 0.1,


    },

    "glc_trs" : {
        "Vmax" : 0.72,
        "Km_glc_cyt" : 2.87,
        "Km_glc_ext" : 2.87,
        "Volume_extracellular2cell": global_params["Volume_extracellular2cell"],
    },

    "glc_diffussion" : {
        "env_glc_level" : 4.5,
        "D" : 0.01, ## !!!!!!!! значение от балды
    },

    "lac_diffussion": {
        "env_lac_level": 0.6,
        "D": 0.01,  ## !!!!!!!! значение от балды
    },

    "pyr_diffussion": {
        "env_pyr_level": 5.0,
        "D": 0.01,  ## !!!!!!!! значение от балды
    },

    "oxigen_diffussion": {
        "env_o2_level": 8.35,
        "D": 0.01,  ## !!!!!!!! значение от балды
    },

}

class Simulator():

    def __init__(self, enzyme_params, metabolites):

        # enzyme_params
        self.enzymes = []

        self.enzymes.append(lib.Oxigen_diffusion(51, enzyme_params["oxigen_diffussion"]))
        self.enzymes.append(lib.Pyruvate_diffusion(72, enzyme_params["pyr_diffussion"]))
        self.enzymes.append(lib.Lactate_diffusion(24, enzyme_params["lac_diffussion"]))
        self.enzymes.append(lib.Glucose_diffusion(0, enzyme_params["glc_diffussion"]))

        self.enzymes.append( lib.GlucoseTransporter(0, 1, enzyme_params["glc_trs"]) )
        self.enzymes.append( lib.Hexokinase(1, 2, 12, 6, enzyme_params["hexokinase"]) )

        # self.enzymes.append( lib.Glucose6phosphate_isomerase(12, 13, enzyme_params["glc6p_isomerase"]) )
        # self.enzymes.append( lib.Phosphofructokinase_type1(13, 2, 15, 3, 8, 14, enzyme_params["phosphofructokinase1"] ) )
        # self.enzymes.append( lib.Fructose16_bisphosphatase(15,13, 8, enzyme_params["fru-1,6-bisphosphatase"] ) )
        # self.enzymes.append( lib.Phosphofructokinase_type2(13, 2, 14, 3, 4, enzyme_params["phosphofructokinase2"] ) )
        # self.enzymes.append( lib.Fructose26_bisphosphatase(14, 13, 8, enzyme_params["fru-2,6-bisphosphatase"] ) )
        # self.enzymes.append( lib.Aldolase(15, 16, 17, enzyme_params["aldolase"] ) )
        # self.enzymes.append( lib.Triosophosphate_isomerase(16, 17, enzyme_params["triosep-isomerase"] ) )
        # self.enzymes.append( lib.Glyceraldehyde_3_phosphate_dehydrogenase(16, 8, 39, 18, 40, enzyme_params["grap_dehydr"] ) )
        # self.enzymes.append( lib.Phosphoglycerate_kinase(18, 3, 19, 2, enzyme_params["p-glyceratekinase"] ) )
        # self.enzymes.append( lib.Phosphoglycerate_mutase(19, 71, enzyme_params["p-gricerate_mutase"] ) )
        # self.enzymes.append( lib.Enolase(71, 20, enzyme_params["enolase"] ) )
        # self.enzymes.append( lib.Pyruvate_kinase(20, 3, 21, 2, enzyme_params["pyruvatekinase"] ) )
        # self.enzymes.append( lib.Lactate_dehydrogenase(21, 40, 23, 39, enzyme_params["LDG"] ) )
        # self.enzymes.append( lib.Monocarboxilate_transporter(24, 23, enzyme_params["MCT"] ) )
        # self.enzymes.append( lib.Creatine_kinase(2, 25, 3, 26,  enzyme_params["creatinekinase"] ) )

        # self.enzymes.append( lib.Malate_dehydrogenase(27, 28, 39, 40,  enzyme_params["malatdehyd"] ) ) # Cytosolic enzyme
        # self.enzymes.append( lib.Malate_dehydrogenase(29, 30, 41, 42,  enzyme_params["malatdehyd"] ) ) # Mitochondrial enzyme
        # self.enzymes.append( lib.Aspartate_aminotransferase(31, 33, 28, 35, enzyme_params["asp_aminotrans"]))  # Cytosolic enzyme
        # self.enzymes.append( lib.Aspartate_aminotransferase(32, 34, 30, 36,  enzyme_params["asp_aminotrans"] ) ) # Mitochondrial enzyme
        # self.enzymes.append( lib.Aspartate_glutamate_carrier(32, 35, 56, 31, 36, 57, 65, enzyme_params["asp_glu_carrier"] ) )
        # self.enzymes.append( lib.Malate_alphaketoglutarate_carrier(27, 34, 29, 33, enzyme_params["mal_akg_carrier"] ) )
        # self.enzymes.append( lib.Glycerol_3phosphate_dehydrogenase_cytosolic(17, 40, 19, 39, enzyme_params["cytgly3pdehyd"] ) )
        # self.enzymes.append( lib.Glycerol_3phosphate_dehydrogenase_mitochondrial(19, 17, 45, 46, 47, 48, enzyme_params["mitgly3pdehyd"] ) )
        # self.enzymes.append( lib.ATP_synthetase(5, 6, 9, 56, 57, 65, enzyme_params["atp_syntase"] ) )
        # self.enzymes.append( lib.ATP_ADP_axchanger(5, 3, 6, 2, 65, enzyme_params["atp/adp_axchanger"] ) )
        # self.enzymes.append( lib.ATP_consumption(2, 3, 8, enzyme_params["atp_consumption"] ) )
        # self.enzymes.append( lib.Passive_efflux_ion(53, 52, 65, enzyme_params["potassium_ed"]))  # Efflux for potassium (K)
        # self.enzymes.append( lib.Passive_efflux_ion(55, 54, 65, enzyme_params["sodium_ed"]))  # Efflux for sodium (Na)
        # self.enzymes.append( lib.Passive_efflux_ion(56, 57, 65, enzyme_params["protons_ed"]))  # Efflux for protons (H)
        # self.enzymes.append( lib.Pump(53, 52, 56, 57, enzyme_params["K_pump"]))  # Pump for potassium (K)
        # self.enzymes.append( lib.Pump(55, 54, 56, 57, enzyme_params["Na_pump"]))  # Pump for sodium (Na)
        # self.enzymes.append( lib.Pump(8, 9, 56, 57, enzyme_params["phos_pump"]))  # Pump for inorganic phosphate
        # self.enzymes.append( lib.Calcium_effux(58, 59, 65, enzyme_params["calcium_ed"]) )
        #
        # # ca_mit, na_cyt, ca_cyt, na_mit, Vmm
        # self.enzymes.append( lib.Ca_Na_pump(59, 55, 58, 54, 65, enzyme_params["ca_na_pump"]) )
        #
        # #  ca_mit, h_cyt, ca_cyt, h_mit, Vmm
        # self.enzymes.append( lib.Ca_H_pump(59, 56, 58, 57, 65, enzyme_params["ca_h_pump"]) )
        #
        # # h_cyt, h_mit, q, qh2, nad, nadh, Vmm
        # self.enzymes.append( lib.Complex1(56, 57, 47, 48, 41, 42, 65, enzyme_params["complex1"]) )
        #
        # # h_cyt, h_mit, q, qh2, cytc_ox, cytc_red, Vmm,
        # self.enzymes.append( lib.Complex3(56, 57, 47, 48, 49, 50, 65, enzyme_params["complex3"]) )
        #
        # # h_cyt, h_mit, cytc_ox, cytc_red, o2, Vmm
        # self.enzymes.append( lib.Complex4(56, 57, 49, 50, 51, 65, enzyme_params["complex4"]) )
        #
        # # pyr_cyt, pyr_mit, h_cyt, h_mit
        # self.enzymes.append( lib.Pyruvate_exchanger(21, 22, 56, 57, enzyme_params["pyr_exchanger"]) )
        #
        # # pyr, CoA, acCoA, fad_pdhc, fadh2_pdhc, nad, nadh, ca
        # self.enzymes.append( lib.Pyruvate_dehydrogenase_complex(22, 60, 61, 66, 67, 41, 42, 59, enzyme_params["pyr_dehyd_comp"]) )
        #
        # # oa, acCoA, CoA, cit
        # self.enzymes.append( lib.Citrate_synthetase(30, 61, 60, 63, enzyme_params["citrate_syntase"]) )
        #
        # # citr, isocitr
        # self.enzymes.append( lib.Aconitase(63, 64, enzyme_params["aconitase"]) )
        #
        # # isocitr, nad, akg, nadh, ca,
        # self.enzymes.append( lib.Isocitrate_dehydrogenase(64, 41, 34, 42, 59, enzyme_params["isocit_dehydr"]) )
        #
        # # ca, akg, nadh, nad, CoA, sucCoA, fad, fadh2
        # self.enzymes.append( lib.Alpha_ketoglutarate_dehydrogenase(59, 34, 42, 41, 60, 62, 68, 69, enzyme_params["akg_dehydr"]) )
        #
        # #  sucCoA, pi, suc, CoA, adp, atp
        # self.enzymes.append( lib.Succinil_CoA_synthetase(62, 9, 37, 60, 6, 5, enzyme_params["sucCoAsyntase_4atp"]) )
        #
        # #  sucCoA, pi, suc, CoA, gdp, gtp
        # self.enzymes.append( lib.Succinil_CoA_synthetase(62, 9, 37, 60, 11, 10, enzyme_params["sucCoAsyntase_4gtp"]) )
        #
        # #  suc, fad, fadh2, fum, q, qh2, mal,
        # self.enzymes.append( lib.Succinate_dehydrydrogenase(37, 43, 44, 38, 47, 48, 29, enzyme_params["suc_dehydr"]) )
        #
        # # mal, fum,
        # self.enzymes.append( lib.Fumarase(29, 38, enzyme_params["fumarase"]) )


    def run_model(self, t, y):
        dydt = np.zeros_like(y) # [0.0 for _ in range(len(y))]

        y[np.isnan(y)] = 0
        y[(y < 0.0000001) & (y > -0.0000001)] = 0


        for enzyme in self.enzymes:
            dydt = enzyme.update(y, dydt)

        # dydt[ ] = -0.01 * y [] # to remove of somebody metabolite

        return dydt


simulalor = Simulator(enzyme_params, metabolites)

y0 = [metabolites[idx]["rest"] for idx in range(len(metabolites))]



# simulalor.run_model(0, y0)
sol = solve_ivp(simulalor.run_model, [0, 0.1], y0, method="LSODA")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

for idx in range(5):  # len(metabolites)
    ax.plot(sol.t,sol.y[idx, :], linewidth = 1, label=metabolites[idx]["full"])

plt.legend()
plt.show()

