import numpy as np
from scipy.integrate import odeint, solve_ivp
import lib

# Notations

metabolites = []

metabolites.append({"idx" : 0, "full" : "External glucose", "short" : "glc_ext"})
metabolites.append({"idx" : 1, "full" : "Cytosolic glucose", "short" : "glc"})
metabolites.append({"idx" : 2, "full" : "Cytosolic ATP", "short" : "atp_cyt"})
metabolites.append({"idx" : 3, "full" : "Cytosolic ADP", "short" : "adp_cyt"})
metabolites.append({"idx" : 4, "full" : "Cytosolic AMP", "short" : "amp_cyt"})
metabolites.append({"idx" : 5, "full" : "Mitochondrial ATP", "short" : "atp_mit"})
metabolites.append({"idx" : 6, "full" : "Mitochondrial ADP", "short" : "adp_mit"})
metabolites.append({"idx" : 7, "full" : "Mitochondrial AMP", "short" : "amp_mit"})
metabolites.append({"idx" : 8, "full" : "Cytosolic inorganic phosphate", "short" : "pi_cyt"})
metabolites.append({"idx" : 9, "full" : "Mitochondrial inorganic phosphate", "short" : "pi_mit"})
metabolites.append({"idx" : 10, "full" : "Mitochondrial GTP", "short" : "gtp_mit"})
metabolites.append({"idx" : 11, "full" : "Mitochondrial GDP", "short" : "gdp_mit"})
metabolites.append({"idx" : 12, "full" : "Glucose-6-phosphate", "short" : "glc6p"})
metabolites.append({"idx" : 13, "full" : "Fructose-6-phosphate", "short" : "fru6p"})
metabolites.append({"idx" : 14, "full" : "Fructose-2,6-bisphosphate", "short" : "fru26p"})
metabolites.append({"idx" : 15, "full" : "Fructose-1,6-bisphosphate", "short" : "fru16p"})
metabolites.append({"idx" : 16, "full" : "Glycerol 3-phosphate", "short" : "grap"})
metabolites.append({"idx" : 17, "full" : "Dihydroxyacetone phosphate", "short" : "dhap"})
metabolites.append({"idx" : 18, "full" : "1,3-Bisphosphoglycerate", "short" : "bpg13"})
metabolites.append({"idx" : 19, "full" : "3-Phosphoglycerate", "short" : "pg3"})
metabolites.append({"idx" : 20, "full" : "Phosphoenolpyruvate", "short" : "pep"})
metabolites.append({"idx" : 21, "full" : "Cytosolic pyruvate", "short" : "pyr_cyt"})
metabolites.append({"idx" : 22, "full" : "Mitochondrial pyruvate", "short" : "pyr_mit"})
metabolites.append({"idx" : 23, "full" : "Cytosolic lactate", "short" : "lac"})
metabolites.append({"idx" : 24, "full" : "External lactate", "short" : "lac_ext"})
metabolites.append({"idx" : 25, "full" : "Creatine", "short" : "cr"})
metabolites.append({"idx" : 26, "full" : "Creatine phosphate", "short" : "crp"})
metabolites.append({"idx" : 27, "full" : "Cytosolic malate", "short" : "mal_cyt"})
metabolites.append({"idx" : 28, "full" : "Cytosolic oxaloacetate", "short" : "oa_cyt"})
metabolites.append({"idx" : 29, "full" : "Mitochondrial malate", "short" : "mal_mit"})
metabolites.append({"idx" : 30, "full" : "Mitochondrial oxaloacetate", "short" : "oa_mit"})
metabolites.append({"idx" : 31, "full" : "Cytosolic aspartate", "short" : "asp_cyt"})
metabolites.append({"idx" : 32, "full" : "Mitochondrial aspartate", "short" : "asp_mit"})
metabolites.append({"idx" : 33, "full" : "Cytosolic alpha-ketoglutarate", "short" : "akg_cyt"})
metabolites.append({"idx" : 34, "full" : "Mitochondrial alpha-ketoglutarate", "short" : "akg_mit"})
metabolites.append({"idx" : 35, "full" : "Cytosolic glutamate", "short" : "glu_cyt"})
metabolites.append({"idx" : 36, "full" : "Mitochondrial glutamate", "short" : "glu_mit"})
metabolites.append({"idx" : 37, "full" : "Succinate", "short" : "suc"})
metabolites.append({"idx" : 38, "full" : "Fumarate", "short" : "fum"})
metabolites.append({"idx" : 39, "full" : "Cytosolic NAD+", "short" : "nad_cyt"})
metabolites.append({"idx" : 40, "full" : "Cytosolic NADH", "short" : "nadh_cyt"})
metabolites.append({"idx" : 41, "full" : "Mitochondrial NAD+", "short" : "nad_mit"})
metabolites.append({"idx" : 42, "full" : "Mitochondrial NADH", "short" : "nadh_mit"})
metabolites.append({"idx" : 43, "full" : "FAD of succinate dehydrogenase", "short" : "fad_sucdh"})
metabolites.append({"idx" : 44, "full" : "FADH2 of succinate dehydrogenase", "short" : "fadh2_sucdh"})
metabolites.append({"idx" : 45, "full" : "FAD of glycerol-3-phosphate dehydrogenase", "short" : "fad_g3dh"})
metabolites.append({"idx" : 46, "full" : "FADH2 of glycerol-3-phosphate dehydrogenase", "short" : "fadh2_g3dh"})
metabolites.append({"idx" : 47, "full" : "Coenzyme Q oxidized", "short" : "Q"})
metabolites.append({"idx" : 48, "full" : "Coenzyme QH2 reduced", "short" : "QH2"})
metabolites.append({"idx" : 49, "full" : "Cytochrome c oxidized", "short" : "cytc_ox"})
metabolites.append({"idx" : 50, "full" : "Cytochrome c reduced", "short" : "cytc_red"})
metabolites.append({"idx" : 51, "full" : "Mitochondrial oxigen (O2)", "short" : "o2_mit"})
metabolites.append({"idx" : 52, "full" : "Mitochondrial potassium (K)", "short" : "k_mit"})
metabolites.append({"idx" : 53, "full" : "Cytosolic potassium (K)", "short" : "k_cyt"})
metabolites.append({"idx" : 54, "full" : "Mitochondrial sodium (Na)", "short" : "na_mit"})
metabolites.append({"idx" : 55, "full" : "Cytosolic sodium (Na)", "short" : "na_cyt"})
metabolites.append({"idx" : 56, "full" : "Cytosolic proton (H+)", "short" : "h+_cyt"})
metabolites.append({"idx" : 57, "full" : "Mitochondrial proton (H+)", "short" : "h+_mit"})
metabolites.append({"idx" : 58, "full" : "Cytosolic calcium (Ca)", "short" : "ca_cyt"})
metabolites.append({"idx" : 59, "full" : "Mitochondrial calcium (Ca)", "short" : "ca_mit"})
metabolites.append({"idx" : 60, "full" : "CoA", "short" : "coa"})
metabolites.append({"idx" : 61, "full" : "Acetyl-CoA", "short" : "acoa"})
metabolites.append({"idx" : 62, "full" : "Succinyl-CoA", "short" : "succoa"})
metabolites.append({"idx" : 63, "full" : "Citrate", "short" : "citr"})
metabolites.append({"idx" : 64, "full" : "Isocitrate", "short" : "isocitr"})
metabolites.append({"idx" : 65, "full" : "Voltage on mitochondrial membrane", "short" : "Vmm"})
metabolites.append({"idx" : 66, "full" : "FAD of pyruvate dehydrogenase complex", "short" : "fad_pdhc"})
metabolites.append({"idx" : 67, "full" : "FADH2 of pyruvate dehydrogenase complex", "short" : "fadh2_pdhc"})
metabolites.append({"idx" : 68, "full" : "FAD of alpha-ketoglutarate dehydrogenase complex", "short" : "fad_akgdhc"})
metabolites.append({"idx" : 69, "full" : "FADH2 of alpha-ketoglutarate dehydrogenase complex", "short" : "fadh2_akgdhc"})
metabolites.append({"idx" : 70, "full" : "CO2", "short" : "co2"})
metabolites.append({"idx" : 71, "full" : "2-Phosphoglycerate", "short" : "pg2"})

"""
enzymes = []
enzymes.append({"reagends":[27, 39], "pruducts":[28, 40], "short":"mal_dehydr", "full":"Cytosolic malate dehydrogenase", "func": lib.getVmal_dehydr })
enzymes.append({"reagends":[29, 41], "pruducts":[30, 42], "short":"mal_dehydr", "full":"Mitochondrial malate dehydrogenase",  "func": lib.getVmal_dehydr })
enzymes.append({"reagends":[38], "pruducts":[29], "short":"fumarase", "full":"Fumarase",  "func": lib.getVfumarase })
enzymes.append({"reagends":[37, 43], "pruducts":[38, 44], "short":"suc_dehydr_s1", "full":"Succinate dehydrogenase stage 1",  "func":lib.getVsuc_dehydrydrogenase_stage1  })
enzymes.append({"reagends":[37, 47], "pruducts":[38, 48], "short":"suc_dehydr_s2", "full":"Succinate dehydrogenase stage 2",  "func":lib.getVsuc_dehydrydrogenase_stage2  })
enzymes.append({"reagends":[62, 6, 9], "pruducts":[37, 60, 5], "short":"sucCoAsyntase_atp", "full":"Succinil-CoA synthetase with ATP",  "func": lib.getVsucCoAsyntase })
enzymes.append({"reagends":[62, 11, 9], "pruducts":[37, 60, 10], "short":"sucCoAsyntase_gtp", "full":"Succinil-CoA synthetase with GTP",  "func": lib.getVsucCoAsyntase })
enzymes.append({"reagends":[34, 68, 60], "pruducts":[62, 69], "short":"akg_dehydr_s1", "full":"Alpha-ketoglutarate dehydrogenase stage 1",  "func": lib.getVakg_dehydrogenase_stage1 })
enzymes.append({"reagends":[69, 41], "pruducts":[68, 42], "short":"akg_dehydr_s2", "full":"Alpha-ketoglutarate dehydrogenase stage 2",  "func": lib.getVakg_dehydrogenase_stage2 })
enzymes.append({"reagends":[64, 41], "pruducts":[34, 42], "short":"isocit_dehydr", "full":"Isocitrate dehydrogenase",  "func": lib.getVisocit_dehydrogenase })
enzymes.append({"reagends":[63], "pruducts":[64], "short":"aconitase", "full":"Aconitase",  "func": lib.getVaconitase })
enzymes.append({"reagends":[30, 61], "pruducts":[63], "short":"citrate_syntase", "full":"Citrate syntase",  "func": lib.getVcitratesyntase })
enzymes.append({"reagends":[22, 60, 66], "pruducts":[61, 69, 70], "short":"pyr_dehyd_comp_s1", "full":"Pyruvate dehydrogenase stage 1",  "func": lib.getVpyr_dehydrogenase_complex_stage1 })
enzymes.append({"reagends":[69, 41], "pruducts":[42, 68], "short":"pyr_dehyd_comp_s2", "full":"Pyruvate dehydrogenase stage 2",  "func": lib.getVpyr_dehydrogenase_complex_stage2 })
enzymes.append({"reagends":[21, 56], "pruducts":[22, 57], "short":"pyr_exchanger", "full":"Pyruvate exchanger",  "func": lib.getVpyr_exchanger })
enzymes.append({"reagends":[50, 51, 57], "pruducts":[49, 56], "short":"complex4", "full":"Complex IV",  "func": lib.getVcomplex4})
enzymes.append({"reagends":[48, 49, 57], "pruducts":[47, 50, 56], "short":"complex3", "full":"Complex III",  "func": lib.getVcomplex3})
enzymes.append({"reagends":[42, 47, 57], "pruducts":[41, 48, 56], "short":"complex1", "full":"Complex I",  "func": lib.getVcomplex1})
enzymes.append({"reagends":[59, 56], "pruducts":[58, 57], "short":"ca_h_pump", "full":"Ca - H pump",  "func": lib.getIca_h_pump})
enzymes.append({"reagends":[59], "pruducts":[58], "short":"calcium_ed", "full":"calcium_ed",  "func": lib.getIca_ed})
enzymes.append({"reagends":[2], "pruducts":[3, 8], "short":"atp_consumption", "full":"ATP consumption",  "func": lib.getVatp_consumption})
enzymes.append({"reagends":[5, 3], "pruducts":[2, 6], "short":"atp/adp_axchanger", "full":"ATP/ADP axchanger",  "func": lib.getVatp_consumption})
enzymes.append({"reagends":[5], "pruducts":[6, 9], "short":"atp_syntase", "full":"ATP synthetase",  "func": lib.getVatp_syntase})
enzymes.append({"reagends":[16, 45], "pruducts":[17, 46], "short":"mitgly3pdehyd_s1", "full":"Mitochondrial glycerol-3-phosphate dehydrogenase stage 1",  "func": lib.getVmitg3pdehyd_stage1})
enzymes.append({"reagends":[46, 47], "pruducts":[45, 48], "short":"mitgly3pdehyd_s2", "full":"Mitochondrial glycerol-3-phosphate dehydrogenase stage 2",  "func": lib.getVmitg3pdehyd_stage2})
enzymes.append({"reagends":[17, 40], "pruducts":[16, 39], "short":"cytgly3pdehyd", "full":"Cytosolic glycerol-3-phosphate dehydrogenase",  "func": lib.getVcytg3pdehyd})
enzymes.append({"reagends":[27, 34], "pruducts":[33, 29], "short":"mal_akg_carrier", "full":"Malate/alpha-ketoglutarate carrier",  "func": lib.getVmal_akg_carrier})
enzymes.append({"reagends":[32, 36, 56], "pruducts":[31,35, 57], "short":"asp_glu_carrier", "full":"Aspartate/glutamate carrier",  "func": lib.getVasp_glu_carrier})
enzymes.append({"reagends":[32, 34], "pruducts":[30, 36], "short":"asp_aminotrans", "full":"Mitochondrial aminotranferase",  "func": lib.getVaspartateaminotransferase})
enzymes.append({"reagends":[31, 33], "pruducts":[28, 35], "short":"asp_aminotrans", "full":"Cytosolic aminotranferase",  "func": lib.getVaspartateaminotransferase})
enzymes.append({"reagends":[2, 25], "pruducts":[26, 3], "short":"creatinekinase", "full":"Creatine kinase",  "func": lib.getVcreatinekinase})
enzymes.append({"reagends":[23], "pruducts":[24], "short":"MCT", "full":"Monocarboxylate transporter",  "func": lib.getVmonocarboxilatetransporter})
enzymes.append({"reagends":[21, 40], "pruducts":[24, 39], "short":"LDG", "full":"Lactate dehydrogenase",  "func": lib.getVlactatedehydrogenase})

"""




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
        "Em_Q": 1,  # !!!!!!!!
        "Em_cytc" : 1, # !!!!!!!
    },

    "complex1" : {
        "Vmax" : 2.25,
        "Em_N": 1,  # !!!!! нет данных
        "Em_Q": 1,  # !!!!! нет данных
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
        "Am": 3.7 * 10 ** -5,  # cm**2
        "Cmm": 0.9 * 10 ** -6,  # F/cm**2
    },

    "phos_pump" : {
        "Vmax"  : 43.3485,
        "is_simport" : True,
    },

    "Na_pump": {
        "Vmax": 5*10**-3,
        "is_simport" : False,
    },

    "K_pump" : {
        "Vmax" : 7.5 * 10**-4,
        "is_simport" : False,
    },

    "protons_ed" : {
        "P" : 2 * 10**-4, # m/s
        "Am": 3.7 * 10 ** -5,  # cm**2
        "Cmm": 0.9 * 10 ** -6,  # F/cm**2
    },

    "sodium_ed" : {
        "P" : 10**-10, # m/s
        "Am": 3.7 * 10 ** -5,  # cm**2
        "Cmm": 0.9 * 10 ** -6,  # F/cm**2
    },

    "potassium_ed" : {
        "P" : 2 * 10**-10, # m/s !!!!!
        "Am": 3.7 * 10 ** -5,  # cm**2
        "Cmm": 0.9 * 10 ** -6,  # F/cm**2
    },

    "atp_consumption" : {
        "Vmax" : 1,
        "Km_atp" : 1,
        "activation" : 0,
    },

    "atp/adp_axchanger" : {
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
    }
}

class Simulator():

    def __init__(self, enzyme_params, metabolites):

        # enzyme_params
        self.enzymes = []

        self.enzymes.append( lib.GlucoseTransporter(0, 1, enzyme_params["glc_trs"]) )
        self.enzymes.append( lib.Hexokinase(1, 2, 12, 6, enzyme_params["hexokinase"]) )
        self.enzymes.append( lib.Glucose6phosphate_isomerase(12, 13, enzyme_params["glc6p_isomerase"]) )
        self.enzymes.append( lib.Phosphofructokinase_type1(13, 2, 15, 3, 8, 14, enzyme_params["phosphofructokinase1"] ) )
        self.enzymes.append( lib.Fructose16_bisphosphatase(15,13, 8, enzyme_params["fru-1,6-bisphosphatase"] ) )
        self.enzymes.append( lib.Phosphofructokinase_type2(13, 2, 14, 3, 4, enzyme_params["phosphofructokinase2"] ) )
        self.enzymes.append( lib.Fructose26_bisphosphatase(14, 13, 8, enzyme_params["fru-2,6-bisphosphatase"] ) )
        self.enzymes.append( lib.Aldolase(15, 16, 17, enzyme_params["aldolase"] ) )
        self.enzymes.append( lib.Triosophosphate_isomerase(16, 17, enzyme_params["triosep-isomerase"] ) )
        self.enzymes.append( lib.Glyceraldehyde_3_phosphate_dehydrogenase(16, 8, 39, 18, 40, enzyme_params["grap_dehydr"] ) )
        self.enzymes.append( lib.Phosphoglycerate_kinase(18, 3, 19, 2, enzyme_params["p-glyceratekinase"] ) )
        self.enzymes.append( lib.Phosphoglycerate_mutase(19, 71, enzyme_params["p-gricerate_mutase"] ) )
        self.enzymes.append( lib.Enolase(71, 20, enzyme_params["enolase"] ) )
        self.enzymes.append( lib.Pyruvate_kinase(20, 3, 21, 2, enzyme_params["pyruvatekinase"] ) )
        self.enzymes.append( lib.Lactate_dehydrogenase(21, 40, 23, 39, enzyme_params["LDG"] ) )
        self.enzymes.append( lib.Monocarboxilate_transporter(24, 23, enzyme_params["MCT"] ) )
        self.enzymes.append( lib.Creatine_kinase(2, 25, 3, 26,  enzyme_params["creatinekinase"] ) )
        self.enzymes.append( lib.Malate_dehydrogenase(27, 28, 39, 40,  enzyme_params["malatdehyd"] ) ) # Cytosolic enzyme
        self.enzymes.append( lib.Malate_dehydrogenase(29, 30, 41, 42,  enzyme_params["malatdehyd"] ) ) # Mitochondrial enzyme
        self.enzymes.append( lib.Aspartate_aminotransferase(31, 33, 28, 35, enzyme_params["asp_aminotrans"]))  # Cytosolic enzyme
        self.enzymes.append( lib.Aspartate_aminotransferase(32, 34, 30, 36,  enzyme_params["asp_aminotrans"] ) ) # Mitochondrial enzyme
        self.enzymes.append( lib.Aspartate_glutamate_carrier(32, 35, 56, 31, 36, 57, 65, enzyme_params["asp_glu_carrier"] ) )
        self.enzymes.append( lib.Malate_alphaketoglutarate_carrier(27, 34, 29, 33, enzyme_params["mal_akg_carrier"] ) )
        self.enzymes.append( lib.Glycerol_3phosphate_dehydrogenase_cytosolic(17, 40, 19, 39, enzyme_params["cytgly3pdehyd"] ) )
        self.enzymes.append( lib.Glycerol_3phosphate_dehydrogenase_mitochondrial(19, 17, 45, 46, 47, 48, enzyme_params["mitgly3pdehyd"] ) )
        self.enzymes.append( lib.ATP_synthetase(5, 6, 9, 56, 57, 65, enzyme_params["atp_syntase"] ) )
        self.enzymes.append( lib.ATP_ADP_axchanger(5, 3, 6, 2, 65, enzyme_params["atp/adp_axchanger"] ) )
        self.enzymes.append( lib.ATP_consumption(2, 3, 8, enzyme_params["atp_consumption"] ) )
        self.enzymes.append( lib.Passive_efflux_ion(53, 52, 65, enzyme_params["potassium_ed"]))  # Efflux for potassium (K)
        self.enzymes.append( lib.Passive_efflux_ion(55, 54, 65, enzyme_params["sodium_ed"]))  # Efflux for sodium (Na)
        self.enzymes.append( lib.Passive_efflux_ion(56, 57, 65, enzyme_params["protons_ed"]))  # Efflux for protons (H)
        self.enzymes.append( lib.Pump(53, 52, 56, 57, enzyme_params["K_pump"]))  # Pump for potassium (K)
        self.enzymes.append( lib.Pump(55, 54, 56, 57, enzyme_params["Na_pump"]))  # Pump for sodium (Na)
        self.enzymes.append( lib.Pump(8, 9, 56, 57, enzyme_params["phos_pump"]))  # Pump for inorganic phosphate
        self.enzymes.append( lib.Calcium_effux(58, 59, 65, enzyme_params["calcium_ed"]) )

        # ca_mit, na_cyt, ca_cyt, na_mit, Vmm
        self.enzymes.append( lib.Ca_Na_pump(59, 55, 58, 54, 65, enzyme_params["ca_na_pump"]) )

        #  ca_mit, h_cyt, ca_cyt, h_mit, Vmm
        self.enzymes.append( lib.Ca_H_pump(59, 56, 58, 57, 65, enzyme_params["ca_h_pump"]) )

        # h_cyt, h_mit, q, qh2, nad, nadh, Vmm
        self.enzymes.append( lib.Complex1(56, 57, 47, 48, 41, 42, 65, enzyme_params["complex1"]) )

        # h_cyt, h_mit, q, qh2, cytc_ox, cytc_red, Vmm,
        self.enzymes.append( lib.Complex3(56, 57, 47, 48, 49, 50, 65, enzyme_params["complex3"]) )

        # h_cyt, h_mit, cytc_ox, cytc_red, o2, Vmm
        self.enzymes.append( lib.Complex4(56, 57, 49, 50, 51, 65, enzyme_params["complex4"]) )

        # pyr_cyt, pyr_mit, h_cyt, h_mit
        self.enzymes.append( lib.Pyruvate_exchanger(21, 22, 56, 57, enzyme_params["pyr_exchanger"]) )

        # pyr, CoA, acCoA, fad_pdhc, fadh2_pdhc, nad, nadh, ca
        self.enzymes.append( lib.Pyruvate_dehydrogenase_complex(22, 60, 61, 66, 67, 41, 42, 59, enzyme_params["pyr_dehyd_comp"]) )

        # oa, acCoA, CoA, cit
        self.enzymes.append( lib.Citrate_synthetase(30, 61, 60, 63, enzyme_params["citrate_syntase"]) )

        # citr, isocitr
        self.enzymes.append( lib.Aconitase(63, 64, enzyme_params["aconitase"]) )

        # isocitr, nad, akg, nadh, ca,
        self.enzymes.append( lib.Isocitrate_dehydrogenase(64, 41, 34, 42, 59, enzyme_params["isocit_dehydr"]) )

        # ca, akg, nadh, nad, CoA, sucCoA, fad, fadh2
        self.enzymes.append( lib.Alpha_ketoglutarate_dehydrogenase(59, 34, 42, 41, 60, 62, 68, 69, enzyme_params["akg_dehydr"]) )

        #  sucCoA, pi, suc, CoA, adp, atp
        self.enzymes.append( lib.Succinil_CoA_synthetase(62, 9, 37, 60, 6, 5, enzyme_params["sucCoAsyntase_4atp"]) )

        #  sucCoA, pi, suc, CoA, gdp, gtp
        self.enzymes.append( lib.Succinil_CoA_synthetase(62, 9, 37, 60, 11, 10, enzyme_params["sucCoAsyntase_4gtp"]) )

        #  suc, fad, fadh2, fum, q, qh2, mal,
        self.enzymes.append( lib.Succinate_dehydrydrogenase(37, 43, 44, 38, 47, 48, 29, enzyme_params["suc_dehydr"]) )

        # mal, fum,
        self.enzymes.append( lib.Fumarase(29, 38, enzyme_params["fumarase"]) )

    def run_model(self, t, y):
        dydt = [0.0 for _ in range(len(y))]

        for enzyme in self.enzymes:
            dydt = enzyme.update(y, dydt)
        return dydt


simulalor = Simulator(enzyme_params, metabolites)

y0 = [1.0 for _ in range(len(metabolites))]
y0[65] = -200


simulalor.run_model(0, y0)
# sol = solve_ivp(simulalor.run_model, [0, 1000], y0)


