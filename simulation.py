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
        self.arg["glc_ext"] = y[0]
        self.arg["glc_cyt"] = y[1]
        self.arg["atp_cyt"] = y[2]
        self.arg["adp_cyt"] = y[3]
        self.arg["amp_cyt"] = y[4]

        self.arg["atp_mit"] = y[5]
        self.arg["adp_mit"] = y[6]
        self.arg["amp_mit"] = y[7]
        self.arg["pi_cyt"] = y[8]
        self.arg["pi_mit"] = y[9]
        self.arg["gtp_mit"] = y[10]
        self.arg["gdp_mit"] = y[11]

        self.arg["glc6p"] = y[12]
        self.arg["fru6p"] = y[13]
        self.arg["fru26p"] = y[14]
        self.arg["fru16bp"] = y[15]
        self.arg["fru16p"] = y[16]

        self.arg["grap"] = y[17]
        self.arg["dhap"] = y[18]
        self.arg["bpg13"] = y[19]
        self.arg["pg3"] = y[20]
        self.arg["pg2"] = y[21]
        self.arg["pep"] = y[22]
        self.arg["pyr_cyt"] = y[23]
        self.arg["pyr_mit"] = y[24]
        self.arg["lac"] = y[25]
        self.arg["lac_ext"] = y[26]
        self.arg["cr"] = y[27]
        self.arg["crp"] = y[28]

        self.arg["mal_cyt"] = y[29]
        self.arg["oa_cyt"] = y[30]
        self.arg["mal_mit"] = y[31]
        self.arg["oa_mit"] = y[32]
        self.arg["asp_cyt"] = y[33]
        self.arg["asp_mit"] = y[34]
        self.arg["akg_cyt"] = y[35]
        self.arg["akg_mit"] = y[36]
        self.arg["glu_cyt"] = y[37]
        self.arg["glu_mit"] = y[38]
        self.arg["suc"] = y[39]
        self.arg["fum"] = y[40]

        self.arg["nad_cyt"] = y[41]
        self.arg["nadh_cyt"] = y[42]
        self.arg["nad_mit"] = y[43]
        self.arg["nadh_mit"] = y[44]
        self.arg["fad"] = y[45]
        self.arg["fadh2"] = y[46]

        self.arg["dhap_cyt"] = y[47]
        self.arg["g3p_cyt"] = y[48]

        self.arg["Q"] = y[49]
        self.arg["QH2"] = y[50]
        self.arg["cytc_ox"] = y[51]
        self.arg["cytc_red"] = y[52]
        self.arg["O2_mit"] = y[53]

        self.arg["K_mit"] = y[54]
        self.arg["K_cyt"] = y[55]
        self.arg["Na_mit"] = y[56]
        self.arg["Na_cyt"] = y[57]
        self.arg["H+_cyt"] = y[58]
        self.arg["H+_mit"] = y[59]
        self.arg["Ca_cyt"] = y[60]
        self.arg["Ca_mit"] = y[61]

        self.arg["CoA"] = y[62]
        self.arg["ACoA"] = y[63]
        self.arg["sucCoA"] = y[64]
        self.arg["fad_pdhg"] = y[65]
        self.arg["fadh2_pdhg"] = y[66]

        self.arg["citr"] = y[67]
        self.arg["isocitr"] = y[68]

        self.arg["mito_membrane"]["Vmm"] = y[69]



        arg = self.arg


        # Glycolysis


        vglc_transp = lib.getVglucosetransporter(arg)                     # Glucose transporter
        vhexokinase = lib.getVhexokinase(arg)                             # Hexokinase
        glucose6p_isomerase = lib.getVglucose6p_isomerase(arg)            # Glucose-6-phosphate isomerase
        phosphofructokinase1 = lib.getVphosphofructokinase1(arg)          # Phosphofructokinase 1
        phosphofructokinase2 = lib.getVphosphofructokinase2(arg)          # Phosphofructokinase 2

        fru16bisphosphatase = lib.getVfru16bisphosphatase(arg)            # Fructose-1,6-bisphosphotase
        fru26bisphosphatase = lib.getVfru26bisphosphatase(arg)            # Fructose-2,6-bisphosphatase
        aldolase = lib.getValdolase(arg)                                  # Aldolase
        triosep_isomerase = lib.getVtriosep_isomerase(arg)                # Triosophosphateisomerase
        grap_dehydr = lib.getVgrap_dehydrogenase(arg)                     # Glyverolphosphatedehydrogenase
        ph_glyceratekinase = lib.getVphosphoglyceratekinase(arg)          # Phosphoglyceratekinase
        ph_glyceratemutase = lib.getVphosphoglyceratemutase(arg)          # Phosphoglyceratemutase
        enolase = lib.getVenolase(arg)                                    # Enolase
        pyruvatekinase = lib.getVpyruvatekinase(arg)                      # Pyruvatekinase
        ldg = lib.getVlactatedehydrogenase(arg)                           # Lactate dedydrogenase
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

        dydt = [0.0 for _ in range(70)]

        dydt[0] = 0                             # external glucose
        dydt[1] = vglc_transp - vhexokinase     # cytosole glucose

        # cytosole ATP
        dydt[2] =  -vhexokinase-phosphofructokinase1-phosphofructokinase2+ph_glyceratekinase+pyruvatekinase-creatinekinase

        dydt[3] = vhexokinase + phosphofructokinase1+phosphofructokinase2 - ph_glyceratekinase-pyruvatekinase+creatinekinase                   # cytosole ADP

        dydt[4] = 0           # cytosole AMP


        dydt[5] = 0        # Mitochondrial ATP
        dydt[6] = 0        # Mitochondrial ADP
        dydt[7] = 0         # Mitochondrial AMP
        dydt[8] = vhexokinase + phosphofructokinase1 + fru16bisphosphatase + fru26bisphosphatase - grap_dehydr        # cytosole nonorganic phosphate
        dydt[9] = 0          # mitochondrial nonorganic phosphate
        dydt[10] = 0         # Mitochondrial GTF
        dydt[11] = 0            # "Mitochondrial GDF

        dydt[12] = vhexokinase - glucose6p_isomerase     # glukose-6-phosphate
        dydt[13] = glucose6p_isomerase - phosphofructokinase1 + fru16bisphosphatase - phosphofructokinase2 + fru26bisphosphatase    # fructose-6-phosphate
        dydt[14] = phosphofructokinase2 - fru26bisphosphatase     # fructose-2,6-phosphate
        dydt[15] = phosphofructokinase1-fru16bisphosphatase-aldolase      # fructose-1,6-bisphosphate

        # dydt[16] = phosphofructokinase1-aldolase      # fructose-1,6-phosphate fru16p

        dydt[17] = aldolase + triosep_isomerase - grap_dehydr     # Glycerol phosphate   "grap"
        dydt[18] = aldolase - triosep_isomerase       #   hap
        dydt[19] = grap_dehydr - ph_glyceratekinase       #   bpg13
        dydt[20] = ph_glyceratekinase - ph_glyceratemutase      #   pg
        dydt[21] = ph_glyceratemutase - enolase      #   pg2
        dydt[22] = enolase - pyruvatekinase      # Phosphoenol pyruvate "pep"
        dydt[23] = pyruvatekinase-ldg       # Cytosole Pyruvate pyr_cyt
        dydt[24] = 0       # Mitochondrial pyruvate  "pyr_mit"
        dydt[25] = ldg+mct       # Cytosole lactate lac
        dydt[26] = 0       # Extracellular lactate    lac_ext
        dydt[27] = -creatinekinase       # Creatine cr
        dydt[28] = creatinekinase       # Creatine phosphate "crp"

        dydt[29] = 0       # Cytosole malate   mal_cyt
        dydt[30] = 0       # Cytosole oxaloacetate oa_cyt
        dydt[31] = 0       # Mitochondrial malate mal_mit
        dydt[32] = 0       # Mitochondrial oxaloacetate   oa_mit
        dydt[33] = 0       # Cytosole aspartate asp_cyt
        dydt[34] = 0       # Mitochondrial aspartate asp_mit
        dydt[35] = 0       # Cytosole alpha-ketoglutorate akg_cyt
        dydt[36] = 0       # Mitochondrial alpha-ketoglutorate akg_mit
        dydt[37] = 0       # Cytosole glutamate glu_cyt
        dydt[38] = 0       # Mitochondrial glutamate glu_mit
        dydt[39] = 0       # Succinate suc
        dydt[40] = 0       # Fumarate fum

        dydt[41] = -grap_dehydr       # Cytosole NAD+ nad_cyt
        dydt[42] = grap_dehydr       # Cytosole NADH nadh_cyt
        dydt[43] = ldg       # Mitochondrial NAD+ nad_mit
        dydt[44] = -ldg       # Mitochondrial NADH nadh_mit
        dydt[45] = 0       # Mitochondrial FAD   fad
        dydt[46] = 0       # Mitochondrial FADH2 fadh2

        dydt[47] = 0       #  dhap_cyt
        dydt[48] = 0       #  g3p_cyt


        dydt[49] = 0        # Coensime Q || Q
        dydt[50] = 0        # Coensime Q reduced "QH2"
        dydt[51] = 0        # Cytochrome c oxidized cytc_ox
        dydt[52] = 0        # cytochrome c reduced cytc_red
        dydt[53] = 0        # Mitochondrial Oxigen O2_mit

        dydt[54] = 0         # Mitochondrial potassium K_mit
        dydt[55] = 0         # Cytosole potassium  K_cyt
        dydt[56] = 0         # Mitochondrial sodium Na_mit
        dydt[57] = 0         # Cytosole sodium Na_cyt
        dydt[58] = 0         # Cytosole proton H+_cyt
        dydt[59] = 0         # Mitochondrial proton H+_mit
        dydt[60] = 0         # Cytosole calcium Ca_cyt
        dydt[61] = 0         # Mitochondrial calcium Ca_mit

        dydt[62] = 0          # Coensime CoA CoA
        dydt[63] = 0          # Acetile coensime CoA ACoA
        dydt[64] = 0          # Succinile CoA sucCoA
        dydt[65] = 0          # fad_pdhg
        dydt[66] = 0          # fadh2_pdhg

        dydt[67] = 0           # Citrate citr
        dydt[68] = 0           # Isocitrate isocitr

        dydt[69] = 0            # Voltage on mitochondrial membrane  mito_membrane-Vmm



        return dydt


simulalor = Simulator(arg)

y0 = [1.0 for _ in range(70)]
y0[-1] = -200


simulalor.model_equations(0, y0)

sol = solve_ivp(simulalor.model_equations, [0, 1000], y0)


"""