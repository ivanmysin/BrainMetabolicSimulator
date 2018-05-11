
import numpy as np
from scipy import constants as const
from scipy.constants import physical_constants

F = physical_constants["Faraday constant"][0]
R = const.R
T = 310 # 37 grad Celcium

def getVglucosetransporter(arg):

    glc_diff = arg["glc_ext"] - arg["glc_cyt"]
    glc_cyt_Km_ratio = arg["glc_cyt"]/arg["glc_trs"]["Km_glc_cyt"]
    glc_ext_Km_ratio = arg["glc_ext"]/arg["glc_trs"]["Km_glc_ext"]
    V = arg["glc_trs"]["Vmax"] * glc_diff / (1 + glc_cyt_Km_ratio + glc_ext_Km_ratio)

    return V

def getVhexokinase(arg):

    Vglc = arg["glc_cyt"] / (arg["glc_cyt"] + arg["hexokinase"]["Km_glc"])
    Vatp = arg["atp_cyt"] / (arg["atp_cyt"] + arg["hexokinase"]["Km_atp"]*(1 + arg["glc6p"]/arg["hexokinase"]["Kiatp_glc6p"]) )
    Inh = 1 - arg["glc6p"]/ arg["hexokinase"]["Ki_glc6p"]
    V = arg["hexokinase"]["Vmax"] * Vglc * Vatp * Inh

    return V

def getVglucose6p_isomerase(arg):

    tmp1 = arg["glc6p"] - arg["glc6p"] / arg["glc6p_isomerase"]["Keq_g6piso"]
    tmp2 = 1 + arg["glc6p"]/arg["glc6p_isomerase"]["Km_glc6p"] + arg["fru6p"] / arg["glc6p_isomerase"]["Km_fru6p"]
    V = arg["glc6p_isomerase"]["Vmax"] * tmp1 / tmp2
    return V

def getVphosphofructokinase1(arg):

    fru__n = arg["fru6p"]**arg["phosphofructokinase1"]["nfru26p"]
    tmp1 = fru__n / (fru__n + arg["phosphofructokinase1"]["Ka_fru26p"]**arg["phosphofructokinase1"]["nfru26p"] )
    Vfru = arg["fru6p"] / (arg["fru6p"] + arg["phosphofructokinase1"]["Km_fru6p"]*(1 - arg["phosphofructokinase1"]["K0"]*tmp1) )
    Vatp = arg["atp_cyt"] / (arg["atp_cyt"] + arg["phosphofructokinase1"]["Km_atp"])
    Vatp__n = 1 - arg["atp_cyt"]**arg["phosphofructokinase1"]["n"]/(arg["atp_cyt"]**arg["phosphofructokinase1"]["n"] + arg["phosphofructokinase1"]["Ki_atp"]**arg["phosphofructokinase1"]["n"])
    Vfru26p = arg["fru26p"] / (arg["fru26p"] + arg["phosphofructokinase1"]["Ka_fru26p"])
    V = arg["phosphofructokinase1"]["Vmax"] * Vfru * Vatp * Vatp__n * Vfru26p

    return V

def getVfru16bisphosphatase(arg):

    V = arg["fru-1,6-bisphosphatase"]["Vmax"]*arg["fru16bp"] / (arg["fru16bp"] + arg["fru-1,6-bisphosphatase"]["Km_fru6p"])

    return V


def getVphosphofructokinase2(arg):

    Vfru = arg["fru6p"] / (arg["fru6p"] + arg["phosphofructokinase2"]["Km_fru6p"])
    Vatp = arg["atp_cyt"] / (arg["atp_cyt"] + arg["phosphofructokinase2"]["Km_atp"])
    Vamp =  arg["amp_cyt"] / (arg["amp_cyt"] + arg["phosphofructokinase2"]["Ka_amp"])
    Vadp = arg["adp_cyt"] / (arg["adp_cyt"] + arg["phosphofructokinase2"]["Ka_adp"])

    V = Vfru * Vatp * Vamp * Vadp
    return V


def getVfru26bisphosphatase(arg):

    tmp = 1 + arg["fru6p"] / arg["fru-2,6-bisphosphatase"]["Ki_fru6p"]
    V = arg["fru-2,6-bisphosphatase"]["Vmax"] * arg["fru26p"] / (arg["fru26p"] + arg["fru-2,6-bisphosphatase"]["Km_fru26p"]*tmp)
    return V

def getValdolase(arg):

    tmp1 = arg["fru16p"] - arg["grap"]*arg["dhap"]/arg["aldolase"]["Keq_aldo"]
    tmp2 = 1 + arg["fru16p"]/arg["aldolase"]["Km_fru16p"]
    tmp3 = 1 + arg["grap"]/arg["aldolase"]["Km_grap"]
    tmp4 = 1 + arg["dhap"] / arg["aldolase"]["Km_dhap"]

    V = arg["aldolase"]["Vmax"] * tmp1 / (tmp2 + tmp3*tmp4 - 1)

    return V

def getVtriosep_isomerase(arg):

    tmp1 = arg["dhap"] - arg["dhap"]/arg["triosep-isomerase"]["Keq_tri"]
    tmp2 = 1 + arg["dhap"] / arg["triosep-isomerase"]["Km_dhap"] + arg["grap"]/arg["triosep-isomerase"]["Km_grap"]

    V = arg["triosep-isomerase"]["Vmax"] * tmp1 * tmp2

    return V

def getVgrap_dehydrogenase(arg):

    tmp1 = arg["nad_cyt"]*arg["grap"]*arg["pi"] - arg["bpg13"]*arg["nadh_cyt"]/arg["grap_dehydr"]["Keq_gapdh"]

    tmp2 = 1 + arg["nad_cyt"]/arg["grap_dehydr"]["Km_nad"]
    tmp3 = 1 + arg["grap"] / arg["grap_dehydr"]["Km_grap"]
    tmp4 = 1 + arg["pi"] / arg["grap_dehydr"]["Km_pi"]

    tmp5 = 1 + arg["nadh_cyt"] / arg["grap_dehydr"]["Km_nadh"]
    tmp6 = 1 + arg["bpg13"] / arg["grap_dehydr"]["Km_bpg13"]

    V = arg["grap_dehydr"]["Vmax"]*tmp1/(tmp2*tmp3*tmp4 + tmp5*tmp6 - 1)

    return V

def getVphosphoglyceratekinase(arg):

    tmp1 = arg["bpg13"]*arg["adp_cyt"] - arg["pg3"]*arg["atp_cyt"]/arg["p-glyceratekinase"]["Keq_pgk"]
    tmp2 = 1 + arg["bpg13"] / arg["p-glyceratekinase"]["Km_bpg13"]
    tmp3 = 1 + arg["adp_cyt"] / arg["p-glyceratekinase"]["Km_adp"]
    tmp4 = 1 +  arg["pg3"] / arg["p-glyceratekinase"]["Km_pg3"]
    tmp5 = 1 + arg["atp_cyt"] / arg["p-glyceratekinase"]["Km_atp"]

    V = arg["p-glyceratekinase"]["Vmax"] * tmp1 / (tmp2*tmp3 + tmp4*tmp5 - 1)

    return V


def getVphosphoglyceratemutase(arg):

    tmp1 = arg["pg3"] - arg["pg2"]/arg["p-gricerate_mutase"]["Keq_pgm"]
    tmp2 = 1 + arg["pg3"]/arg["p-gricerate_mutase"]["Km_pg3"]
    tmp3 = 1 + arg["pg2"] / arg["p-gricerate_mutase"]["Km_pg2"]

    V = arg["p-gricerate_mutase"]["Vmax"] * tmp1 / (tmp2 + tmp3 - 1)

    return V

def getVenolase(arg):

    tmp1 = arg["pg2"] - arg["pep"]/ arg["enolase"]["Keq_eno"]
    tmp2 = 1 + arg["pg2"] / arg["enolase"]["Km_pg2"]
    tmp3 = 1 + arg["pep"] / arg["enolase"]["Km_pep"]

    V = arg["enolase"]["Vmax"] * tmp1 / (tmp2 + tmp3 - 1)

    return V

def getVpyruvatekinase(arg):
    Vpep = arg["pep"] / (arg["pep"] + arg["pyruvatekinase"]["Km_pep"] )
    tmp = 1 + arg["atp_cyt"]/arg["pyruvatekinase"]["Ki_atp"]
    Vadp = arg["adp_cyt"] / (arg["adp_cyt"] + arg["pyruvatekinase"]["Km_adp"] * tmp)

    V = arg["pyruvatekinase"]["Vmax"] * Vpep * Vadp

    return V

def getVlactatedehydrogenase(arg):

    tmp1 = arg["pyr_cyt"] * arg["nadh_cyt"] - arg["lac"]*arg["nad_cyt"]/arg["LDG"]["Keq_ldg"]
    tmp2 = 1 + arg["pyr_cyt"] / arg["LDG"]["Km_pyr"]
    tmp3 = 1 + arg["nadh_cyt"] / arg["LDG"]["Km_nadh"]
    tmp4 = 1 + arg["lac"] / arg["LDG"]["Km_lac"]
    tmp5 = 1 + arg["nad_cyt"] / arg["LDG"]["Km_nad"]

    V = arg["LDG"]["Vmax"] * tmp1 / (tmp2*tmp3 + tmp4*tmp5 - 1)

    return V

def getVmonocarboxilatetransporter(arg):

    tmp1 = arg["lac"] - arg["lac_ext"] / arg["MCT"]["Keq_mct"]
    tmp2 = 1 + arg["lac"] / arg["MCT"]["Km_lac"]
    tmp3 = 1 + arg["lac_ext"] / arg["MCT"]["Km_lac_ex"]

    V = arg["MCT"]["Vmax"] * tmp1 / (tmp2 + tmp3 - 1)

    return V

def getVcreatinekinase(arg):

    tmp = arg["atp_cyt"]*arg["cr"] / arg["creatinekinase"]["Keq"]
    V = arg["creatinekinase"]["Vmax"]*(arg["adp_cyt"]*arg["crp"] - tmp)

    return V

def getVmalatdehydrogenase(arg, mode="mit"):

    if mode == "cyt":
        mal = arg["mal_cyt"]
        oa = arg["oa_cyt"]
        nad = arg["nad_cyt"]
        nadh = arg["nadh_cyt"]

    elif mode == "mit":
        mal = arg["mal_mit"]
        oa = arg["oa_mit"]
        nad = arg["nad_mit"]
        nadh = arg["nadh_mit"]

    tmp1 = mal*nad - oa*nadh/arg["malatdehyd"]["Keq"]
    tmp2 = 1 + mal / arg["malatdehyd"]["Km_mal"]
    tmp3 = 1 + nad / arg["malatdehyd"]["Km_nad"]
    tmp4 = 1 + oa / arg["malatdehyd"]["Km_oa"]
    tmp5 = 1 + nadh / arg["malatdehyd"]["Km_nadh"]

    V = arg["malatdehyd"]["Vmax"] * tmp1 / (tmp2*tmp3 + tmp4*tmp5 - 1)

    return V

def getVaspartateaminotransferase(arg, mode="mit"):

    if mode == "cyt":
        asp = arg["asp_cyt"]
        oa = arg["oa_cyt"]
        akg = arg["nad_cyt"]
        glu = arg["nadh_cyt"]
    elif mode == "mit":
        asp = arg["asp_mit"]
        oa = arg["oa_mit"]
        akg = arg["nad_mit"]
        glu = arg["nadh_mit"]

    tmp = oa * glu / arg["asp_aminotrans"]["Keq"]
    V = arg["asp_aminotrans"]["Vmax"] * (asp * akg - tmp)

    return V

def getVasp_glu_carrier(arg):

    dG = -arg["mito_membrane"]["Vmm"] + 1000 * R * T / F * np.log(arg["mito_membrane"]["h_cyt"] / arg["mito_membrane"]["h_mit"])

    Keq = np.exp(F * dG/(1000 * R * T) )

    tmp1 = arg["asp_mit"] * arg["glu_cyt"] - arg["asp_cyt"] * arg["glu_mit"] / Keq
    tmp2 = (arg["asp_mit"] + arg["asp_glu_carrier"]["Km_asp_mit"]) * (arg["glu_cyt"] +  arg["asp_glu_carrier"]["Km_glu_cyt"])
    tmp3 = (arg["asp_cyt"] + arg["asp_glu_carrier"]["Km_asp_cyt"]) * (arg["glu_mit"] +  arg["asp_glu_carrier"]["Km_glu_mit"])

    V = arg["asp_glu_carrier"]["Vmax"] * tmp1 / (tmp2 + tmp3)

    return V

def getVmal_akg_carrier(arg):

    tmp1 = arg["mal_cyt"]*arg["akg_mit"] - arg["mal_mit"]*arg["akg_cyt"]
    tmp2 = (arg["mal_cyt"] + arg["mal_akg_carrier"]["Km_mal_cyt"]) * (arg["akg_mit"] + arg["mal_akg_carrier"]["Km_akg_mit"])
    tmp3 = (arg["mal_mit"] + arg["mal_akg_carrier"]["Km_mal_mit"]) * (arg["akg_cyt"] + arg["mal_akg_carrier"]["Km_akg_cyt"])

    V = arg["mal_akg_carrier"]["Vmax"] * tmp1 / (tmp2 + tmp3)

    return V
