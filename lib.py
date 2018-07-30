
import numpy as np
from scipy import constants as const
from scipy.constants import physical_constants

F = physical_constants["Faraday constant"][0]
R = const.R
T = 310 # 37 grad Celcium


class Eznyme:

    def __init__(self):
        pass

    def update(self):
        pass



class GlucoseTransporter(Eznyme):

    def __init__(self, glc_ext, glc_cyt, Km_glc_cyt, Km_glc_ext, Vmax):
        self.glc_ext_idx = glc_ext
        self.glc_cyt_idx = glc_cyt
        self.Km_glc_cyt = Km_glc_cyt
        self.Km_glc_ext = Km_glc_ext
        self.Vmax = Vmax

    def update(self, metabolits, dydt):

        glc_ext = metabolits[self.glc_ext_idx]
        glc_cyt = metabolits[self.glc_cyt_idx]

        glc_diff = glc_ext - glc_cyt

        glc_cyt_Km_ratio = glc_cyt / self.Km_glc_cyt
        glc_ext_Km_ratio = glc_ext / self.Km_glc_cyt

        V = self.Vmax * glc_diff / (1 + glc_cyt_Km_ratio + glc_ext_Km_ratio)

        dydt[self.glc_ext_idx] -= V
        dydt[self.glc_cyt_idx] += V

        return dydt

class Hexokinase(Eznyme):

    def __init__(self, glc_cyt, atp_cyt, glc6p, adp_cyt, Vmax, Km_glc, Km_atp, Ki_atp, Ki_glc6p):

        self.glc_cyt_idx = glc_cyt
        self.atp_cyt_idx = atp_cyt
        self.glc6p_idx = glc6p
        self.adp_cyt_idx = adp_cyt
        self.Vmax = Vmax
        self.Km_glc = Km_glc
        self.Km_atp = Km_atp
        self.Ki_atp = Ki_atp
        self.Ki_glc6p = Ki_glc6p


    def update(self, metabolits, dydt):
        glc_cyt = metabolits[self.glc_cyt_idx]
        atp_cyt = metabolits[self.adp_cyt_idx]
        glc6p = metabolits[self.glc6p_idx]

        Vglc = glc_cyt / (glc_cyt + self.Km_glc)
        Vatp = atp_cyt / (atp_cyt + self.Km_atp * (1 + glc6p / self.Ki_atp) )
        Inh = 1 - glc6p/ self.Ki_glc6p
        V = self.Vmax * Vglc * Vatp * Inh

        dydt[self.glc_cyt_idx] -= V
        dydt[self.atp_cyt_idx] -= V

        dydt[self.glc6p_idx] += V
        dydt[self.adp_cyt_idx] += V

        return dydt


class Glucose6phosphate_Isomerase(Eznyme):

    def __init__(self, glc6p, fru6p, Vmax, Keq, Km_glc6p, Km_fru6p):

        self.glc6p_idx = glc6p
        self.fru6p_idx = fru6p
        self.Vmax = Vmax
        self.Keq = Keq
        self.Km_glc6p = Km_glc6p
        self.Km_fru6p = Km_fru6p


    def update(self, metabolits, dydt):
        glc6p = metabolits[self.glc6p_idx]
        fru6p = metabolits[self.fru6p_idx]

        tmp1 = glc6p - fru6p / self.Keq
        tmp2 = 1 + glc6p / self.Km_glc6p + fru6p / self.Km_fru6p
        V = self.Vmax * tmp1 / tmp2

        dydt[self.glc6p_idx] += V
        dydt[self.fru6p_idx] -= V

        return dydt

def Phosphofructokinase_type1(Eznyme):

    def __init__(self, fru6p, atp_cyt, fru16p, adp_cyt, pi_cyt, fru26p, Vmax, Km_fru6p, Km_atp, n, Ki_atp, K0, Ka_fru26p, n_fru26p, Kafru26p):

        self.fru6p_idx = fru6p
        self.atp_cyt_idx = atp_cyt
        self.fru16p_idx = fru16p
        self.adp_cyt_idx = adp_cyt
        self.pi_cyt_idx = pi_cyt
        self.fru26p_idx = fru26p
        self.Vmax = Vmax
        self.Km_fru6p = Km_fru6p
        self.Km_atp = Km_atp
        self.n = n
        self.Ki_atp = Ki_atp
        self.K0 = K0
        self.Ka_fru26p = Ka_fru26p
        self.n_fru26p = n_fru26p
        self.Kafru26p = Kafru26p


    def update(self, metabolits, dydt):

        fru6p = metabolits[self.fru6p_idx]
        atp_cyt = metabolits[self.atp_cyt_idx]
        fru26p = metabolits[self.fru16p_idx]

        fru__n = fru6p**self.n_fru26p
        tmp1 = fru__n / (fru__n + self.Ka_fru26p**self.nfru26p )

        Vfru = fru6p / (fru6p + self.Km_fru6p * (1 - self.K0 * tmp1) )
        Vatp = atp_cyt / (atp_cyt + self.Km_atp)
        Vatp__n = 1 - atp_cyt**self.n /(atp_cyt**self.n + self.Ki_atp**self.n)

        Vfru26p = fru26p / (fru26p + self.Ka_fru26p)
        V = self.Vmax * Vfru * Vatp * Vatp__n * Vfru26p

        dydt[self.fru6p_idx] -= V
        dydt[self.atp_cyt_idx] -= V


        dydt[self.adp_cyt_idx] += V
        dydt[self.fru16p_idx] += V

        return dydt
#####################################################################################################
class Fructose16_bisphosphatase(Eznyme):
    def __init__(self, fru16p, fru6p, pi_cyt, Vmax, Km):
        self.fru16p_idx = fru16p
        self.fru6p_idx = fru6p
        self.pi_cyt_idx = pi_cyt
        self.Vmax = Vmax
        self.Km = Km



    def update(self, metabolits, dydt):

        fru16bp = metabolits[self.fru16p_idx]

        V = self.Vmax * fru16bp / (fru16bp + self.Km)

        dydt[self.fru16p_idx] -= V
        dydt[self.fru6p_idx] += V
        dydt[self.pi_cyt_idx] += V

        return dydt
#####################################################################################################

class Phosphofructokinase_type2(Eznyme):

    def __init__(self, fru6p, atp_cyt, fru26p, adp_cyt, amp, Vmax, Km_fru6p, Km_atp, Ka_amp, Ka_adp):
        self.fru6p_idx = fru6p
        self.atp_cyt_idx = atp_cyt
        self.fru26p_idx = fru26p
        self.adp_cyt_idx = adp_cyt
        self.amp_idx = amp
        self.Vmax = Vmax
        self.Km_fru6p = Km_fru6p
        self.Km_atp = Km_atp
        self.Ka_amp = Ka_amp
        self.Ka_adp = Ka_adp

    def update(self, metabolits, dydt):
        fru6p = metabolits[self.fru6p_idx]
        atp_cyt = metabolits[self.atp_cyt_idx]
        amp_cyt = metabolits[self.amp_idx]
        adp_cyt = metabolits[self.adp_cyt_idx]

        Vfru = fru6p / (fru6p + self.Km_fru6p)
        Vatp = atp_cyt / (atp_cyt + self.Km_atp)
        Vamp = amp_cyt / (amp_cyt + self.Ka_amp)
        Vadp = adp_cyt / (adp_cyt + self.Ka_adp)

        V = Vfru * Vatp * Vamp * Vadp

        dydt[self.fru6p_idx] -= V
        dydt[self.atp_cyt_idx] -= V

        dydt[self.adp_cyt_idx] += V
        dydt[self.fru26p_idx] += V

        return dydt
######################################################################################################

class Fructose26_bisphosphatase(Eznyme):

    def __init__(self, fru26p, fru6p, pi_cyt, Vmax, Km, Ki_fru6p):
        self.fru26p_idx = fru26p
        self.fru6p_idx = fru6p
        self.pi_cyt_idx = pi_cyt
        self.Vmax = Vmax
        self.Km = Km
        self.Ki_fru6p = Ki_fru6p



    def update(self, metabolits, dydt):
        fru26p = metabolits[self.fru26p_idx]
        fru6p = metabolits[self.fru6p_idx]

        tmp = 1 + fru6p / self.Ki_fru6p
        V = self.Vmax * fru26p / (fru26p + self.Km_fru26p * tmp)

        dydt[self.fru26p_idx] -= V

        dydt[self.fru6p_idx] += V
        dydt[self.pi_cyt_idx] += V

        return dydt


#######################################################################################################

class Aldolase(Eznyme):

    def __init__(self, fru16p, grap, dhap, Vmax, Keq, Km_fru16p, Km_grap, Km_dhap):
        self.fru16p_idx = fru16p
        self.grap_idx = grap
        self.dhap_idx = dhap
        self.Vmax = Vmax
        self.Keq = Keq
        self.Km_fru16p = Km_fru16p
        self.Km_grap = Km_grap
        self.Km_dhap = Km_dhap



    def update(self, metabolits, dydt):
        fru16p = metabolits[self.fru16p_idx]
        grap = metabolits[self.grap_idx]
        dhap = metabolits[self.dhap_idx]


        tmp1 = fru16p - grap * dhap / self.Keq
        tmp2 = 1 + fru16p / self.Km_fru16p
        tmp3 = 1 + grap / self.Km_grap
        tmp4 = 1 + dhap / self.Km_dhap

        V = self.Vmax * tmp1 / (tmp2 + tmp3*tmp4 - 1)

        dydt[self.fru16p_idx] -= V

        dydt[self.grap_idx] += V
        dydt[self.dhap_idx] += V

        return dydt
#######################################################################################################

class Triosofosphate_Isomerase(Eznyme):

    def __init__(self, grap, dhap, Vmax, Keq, Km_grap, Km_dhap):
        self.grap_idx = grap
        self.dhap_idx = dhap
        self.Vmax = Vmax
        self.Keq = Keq
        self.Km_grap = Km_grap
        self.Km_dhap = Km_dhap

    def update(self, metabolits, dydt):
        grap = metabolits[self.grap_idx]
        dhap = metabolits[self.dhap_idx]

        tmp1 = dhap - grap / self.Keq
        tmp2 = 1 + dhap / self.Km_dhap + grap / self.Km_grap

        V = self.Vmax * tmp1 * tmp2

        dydt[self.dhap_idx] -= V
        dydt[self.grap_idx] += V

        return grap

######################################################################################################

class Glyceraldehyde_3_phosphate_dehydrogenase(Eznyme):

    def __init__(self, grap, pi, nad, bpg13, nadh, Vmax, Km_nad, Km_grap, Km_pi, Km_nadh, Km_bpg13):
        self.grap_idx = grap
        self.pi_idx = pi
        self.nad_idx = nad
        self.bpg13_idx = bpg13
        self.nadh_idx = nadh
        self.Vmax = Vmax
        self.Km_nad = Km_nad
        self.Km_grap = Km_grap
        self.Km_pi = Km_pi
        self.Km_nadh = Km_nadh
        self.Km_bpg13 = Km_bpg13

    def update(self, metabolits, dydt):

        grap = metabolits[self.grap_idx]
        pi = metabolits[self.pi_idx]
        nad = metabolits[self.nad_idx]
        bpg13 = metabolits[self.bpg13_idx]
        nadh = metabolits[self.nadh_idx]


        tmp1 = nad * grap * pi - bpg13 * nadh / self.Keq

        tmp2 = 1 + nad / self.Km_nad
        tmp3 = 1 + grap / self.Km_grap
        tmp4 = 1 + pi / self.Km_pi

        tmp5 = 1 + nadh / self.Km_nadh
        tmp6 = 1 + bpg13 / self.Km_bpg13

        V = self.Vmax * tmp1 / (tmp2*tmp3*tmp4 + tmp5*tmp6 - 1)

        dydt[self.grap_idx] -= V
        dydt[self.pi_idx] -= V
        dydt[self.nad_idx] -= V

        dydt[self.bpg13_idx] += V
        dydt[self.nadh_idx] += V

        return dydt

#####################################################################################################

class Phosphoglycerate_kinase(Eznyme):

    def __init__(self, bpg13, adp, pg3, atp, Vmax, Keq, Km_bpg13, Km_adp, Km_pg3, Km_atp):

        self.bpg13_idx = bpg13
        self.adp_idx = adp
        self.atp_idx = atp
        self.pg3_idx = pg3
        self.Vmax = Vmax
        self.Keq = Keq
        self.Km_bpg13 = Km_bpg13
        self.Km_adp = Km_adp
        self.Km_pg3 = Km_pg3
        self.Km_atp = Km_atp

    def update(self, metabolits, dydt):

        bpg13 = metabolits[self.bpg13_idx]
        adp_cyt = metabolits[self.adp_idx ]
        pg3 = metabolits[self.pg3_idx]
        atp_cyt = metabolits[self.atp_idx ]

        tmp1 = bpg13 * adp_cyt - pg3 * atp_cyt / self.Keq_pgk
        tmp2 = 1 + bpg13 / self.Km_bpg13
        tmp3 = 1 + adp_cyt / self.Km_adp
        tmp4 = 1 +  pg3 / self.Km_pg3
        tmp5 = 1 + atp_cyt / self.Km_atp
        V = self.Vmax * tmp1 / (tmp2*tmp3 + tmp4*tmp5 - 1)

        dydt[self.bpg13_idx] -= V
        dydt[self.adp_idx] -= V

        dydt[self.atp_idx] += V
        dydt[self.pg3_idx] += V

        return dydt

#######################################################################################################

class Phosphoglycerate_mutase(Eznyme):

    def __init__(self, pg3, pg2, params):

        self.pg3_idx = pg3
        self.pg2_idx = pg2

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_pg2 = params["Km_pg2"]


    def update(self, metabolits, dydt):
        pg2 = metabolits[self.pg2_idx]
        pg3 = metabolits[self.pg3_idx]

        tmp1 = pg3 - pg2 / self.Keq
        tmp2 = 1 + pg3 / self.Km_pg3
        tmp3 = 1 + pg2 / self.Km_pg2

        V = self.Vmax * tmp1 / (tmp2 + tmp3 - 1)

        dydt[self.pg3_idx] -= V

        dydt[self.pg2_idx] += V

        return dydt

#######################################################################################################

class Enolase(Eznyme):

    def __init__(self, pg2, pep, params):
        self.pg2_idx = pg2
        self.pep_idx = pep

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_pg2 = params["Km_pg2"]
        self.Km_pep = params["Km_pep"]

    def update(self, metabolits, dydt):
        pg2 = metabolits[self.pg2_idx]
        pep = metabolits[self.pep_idx]

        tmp1 = pg2 - pep / self.Keq
        tmp2 = 1 + pg2 / self.Km_pg2
        tmp3 = 1 + pep / self.Km_pep

        V = self.Vmax * tmp1 / (tmp2 + tmp3 - 1)

        dydt[self.pg2_idx] -= V
        dydt[self.pep_idx] += V

        return dydt
######################################################################################################

class Pyruvate_kinase(Eznyme):

    def __init__(self, pep, adp, pyr, atp, params):

        self.pep_idx = pep
        self.adp_cyt_idx = adp
        self.atp_cyt_idx = atp
        self.pyr_cyt_idx = pyr

        self.Vmax = params["Vmax"]
        self.Km_pep = params["Km_pep"]
        self.Km_adp = params["Km_adp"]
        self.Ki_atp = params["Ki_atp"]

    def update(self, metabolits, dydt):

        pep = metabolits[self.pep_idx]
        adp = metabolits[self.adp_cyt_idx]
        atp = metabolits[self.atp_cyt_idx]

        Vpep = pep / (pep + self.Km_pep )
        tmp = 1 + atp / self.Ki_atp
        Vadp = adp / ( adp + self.Km_adp * tmp)

        V = self.Vmax * Vpep * Vadp

        dydt[self.pep_idx] -= V
        dydt[self.adp_cyt_idx] -= V

        dydt[self.atp_cyt_idx] += V
        dydt[self.pyr_cyt_idx] += V

        return dydt

########################################################################################################
class Lactate_dehydrogenase(Eznyme):

    def __init__(self, pyr, nadh, lac, nad, params):

        self.pyr_idx = pyr
        self.nadh_idx = nadh
        self.lac_idx = lac
        self.nad_idx = nad

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_pyr = params["Km_pyr"]
        self.Km_nadh = params["Km_nadh"]
        self.Km_nad = params["Km_nad"]
        self.Km_lac = params["Km_lac"]


    def update(self, metabolits, dydt):

        pyr = metabolits[self.pyr_idx]
        nadh = metabolits[self.nadh_idx]
        nad = metabolits[self.nad_idx]
        lac = metabolits[self.lac_idx]

        tmp1 = pyr * nadh - lac * nad / self.Keq
        tmp2 = 1 + pyr / self.Km_pyr
        tmp3 = 1 + self.nadh / self.Km_nadh
        tmp4 = 1 + lac / self.Km_lac
        tmp5 = 1 + nad / self.Km_nad

        V = self.Vmax * tmp1 / (tmp2*tmp3 + tmp4*tmp5 - 1)

        dydt[self.pyr_idx] -= V
        dydt[self.nadh_idx] -= V

        dydt[self.lac_idx] += V
        dydt[self.nad_idx] += V

        return dydt

########################################################################################################

class Monocarboxilate_transporter(Eznyme):

    def __init__(self, lac_ext, lac_cyt, params):
        self.lac_ext_idx = lac_ext
        self.lac_cyt_idx = lac_cyt

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_lac_cyt = params["Km_lac_cyt"]
        self.Km_lac_ext = params["Km_lac_ext"]


    def update(self, metabolits, dydt):

        lac_ext = metabolits[self.lac_ext_idx]
        lac_cyt = metabolits[self.lac_cyt_idx]


        tmp1 = lac_cyt - lac_ext / self.Keq
        tmp2 = 1 + lac_cyt / self.Km_lac_cyt
        tmp3 = 1 + lac_ext / self.Km_lac_ext

        V = self.Vmax * tmp1 / (tmp2 + tmp3 - 1)

        dydt[self.lac_cyt_idx] -= V
        dydt[self.lac_ext_idx] += V

        return dydt

#########################################################################################################

class Creatine_kinase(Eznyme):

    def __init__(self, atp, cr, adp, crp, params):

        self.atp_idx = atp
        self.cr_idx = cr
        self.adp_idx = adp
        self.crp_idx = crp

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]


    def update(self, metabolits, dydt):

        atp = metabolits[self.atp_idx]
        adp = metabolits[self.adp_idx]
        cr = metabolits[self.cr_idx]
        crp = metabolits[self.crp_idx]


        tmp = atp * cr / self.Keq
        V = self.Vmax * (adp * crp - tmp)

        dydt[self.atp_idx] -= V
        dydt[self.cr_idx] -= V

        dydt[self.crp_idx] += V
        dydt[self.adp_idx] += V

        return dydt
#########################################################################################################


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

def getVcytg3pdehyd(arg):
     tmp1 = arg["dhap_cyt"] * arg["nadh_cyt"] - arg["g3p_cyt"] * arg["nad_cyt"] / arg["cytgly3pdehyd"]["Keq"]

     tmp2 = 1 + arg["dhap_cyt"]/arg["cytgly3pdehyd"]["Km_dhap"]
     tmp3 = 1 + arg["nadh_cyt"] / arg["cytgly3pdehyd"]["Km_nadh"]

     tmp4 = 1 + arg["g3p_cyt"] / arg["cytgly3pdehyd"]["Km_g3p"]
     tmp5 = 1 + arg["nad_cyt"] / arg["cytgly3pdehyd"]["Km_nad"]

     V = arg["cytgly3pdehyd"]["Vmax"] * tmp1 / (tmp2*tmp3 + tmp4*tmp5 - 1)

     return V

def getVmitg3pdehyd_stage1(arg):

    Keq_g3pdh = np.exp( (2*arg["mitgly3pdehyd"]["Em_dhap_g3p"] - 2*arg["mitgly3pdehyd"]["Em_FAD_g3p"])*F / 1000 / R / T )
    tmp1 = arg["dhap_cyt"] * arg["fad"] - arg["g3p_cyt"]*arg["fadh2"] / Keq_g3pdh

    tmp2 = 1 + arg["dhap_cyt"] / arg["mitgly3pdehyd"]["Km_dhap"]
    tmp3 = 1 + arg["g3p_cyt"] / arg["mitgly3pdehyd"]["Km_g3p"]

    Vg3pdh = arg["mitgly3pdehyd"]["Vmax_g3pdh"] * tmp1 / (tmp2 + tmp3)

    return Vg3pdh

def getVmitg3pdehyd_stage2(arg):

    Keq_fad_Q = np.exp( (2*arg["mitgly3pdehyd"]["Em_FAD_g3p"] + 2*arg["mitgly3pdehyd"]["Em_Q"] )*F / 1000 / R / T )
    Vqh2 = arg["mitgly3pdehyd"]["Vmax_Q"] * (arg["fadh2"] * arg["Q"] - arg["fad"]*arg["QH2"] / Keq_fad_Q  )

    return Vqh2


def getVatp_syntase(arg):

    dG = -arg["mito_membrane"]["Vmm"] + R * T * np.log( arg["mito_membrane"]["h_cyt"] / arg["mito_membrane"]["h_mit"]  ) / (1000 * F)
    Vmax = 1.8 * 10**-16 * dG**arg["atp_syntase"]["n"]

    U = arg["mito_membrane"]["Vmm"] * F / (1000 * R * T)

    Keq = np.exp( arg["atp_syntase"]["dG0"]/(R * T) - arg["atp_syntase"]["k"] * U ) * (( arg["mito_membrane"]["h_cyt"] / arg["mito_membrane"]["h_mit"]  )**arg["atp_syntase"]["k"])


    V = Vmax * (arg["adp_mit"] * arg["pi_mit"] - arg["atp_mit"] / Keq  )


    return V


def getVatp_atp_axchanger(arg):

    U = arg["mito_membrane"]["Vmm"] * F / (1000 * R * T)

    tmp1 = 1 - np.exp(U) * arg["atp_cyt"] * arg["adp_mit"] / ( arg["adp_cyt"] * arg["atp_mit"])
    tmp2 = 1 + arg["atp_cyt"] / arg["adp_cyt"]  * np.exp(arg["atp/atp_axchanger"]["S_Vmm"] * U)
    tmp3 = 1 + arg["atp_mit"] / arg["adp_mit"]


    V = arg["atp/atp_axchanger"]["Vmax"] * tmp1 / (tmp2 * tmp3)

    return V

def getVatp_consumption(arg):
    V = arg["atp_consumption"]["Vmax"] * arg["atp_cyt"] / (arg["atp_cyt"] + arg["atp_consumption"]["Km_atp"] ) * (1 + arg["atp_consumption"]["activation"])

    return V

def getIed(arg, ion="K"):
    if ion == "K":
        Ion_in = arg["K_cyt"]
        Ion_mit = arg["K_mit"]
        P = arg["potassium_ed"]["P_Kmit"]
    elif ion == "Na":
        Ion_in = arg["Na_cyt"]
        Ion_mit = arg["Na_mit"]
        P = arg["sodium_ed"]["P_Namit"]
    elif ion == "H+":
        Ion_in = arg["H+_cyt"]
        Ion_mit = arg["H+_mit"]
        P = arg["protons_ed"]["P_H_mit"]



    U = arg["mito_membrane"]["Vmm"] * F / (1000 * R * T)
    tmp = (Ion_in - Ion_mit*np.exp(U)) / (1 - np.exp(U))
    I = arg["mito_membrane"]["Am"] * P * U * F * tmp

    return I


def getVpumps(arg, ion="K"):

    if ion == "K":
        ion_in = arg["K_cyt"]
        ion_mit = arg["K_mit"]
        Vmax = arg["K_pump"]["Vmax"]

    elif ion == "Na":
        ion_in = arg["Na_cyt"]
        ion_mit = arg["Na_mit"]
        Vmax = arg["Na_pump"]["Vmax"]

    elif ion == "pi":
        ion_in = arg["pi_cyt"]
        ion_mit = arg["pi_mit"]
        Vmax = arg["phos_pump"]["Vmax"]

    if ion=="K" or ion=="Na":
        V = Vmax * (ion_in*arg["H+_mit"] - ion_mit*arg["H+_cyt"])
    else:
        V = Vmax * (ion_in * arg["H+_cyt"] - ion_mit * arg["H+_mit"])

    return V


def getIca_ed(arg):

    U = arg["mito_membrane"]["Vmm"] * F / (1000 * R * T)
    tmp1 = (arg["Ca_cyt"] - arg["Ca_mit"] * np.exp(2 * U)) / (1 - np.exp(2 * U))
    tmp2 = arg["calcium_ed"]["P_RMC"]*(1 - arg["Ca_cyt"]/(arg["Ca_cyt"] + arg["calcium_ed"]["Ki_cacyt" ]))
    tmp3 = arg["calcium_ed"]["P_Mcu"] * arg["Ca_cyt"]**arg["calcium_ed"]["n"] / (arg["Ca_cyt"]**arg["calcium_ed"]["n"] + arg["calcium_ed"]["Km_Mcu"]**arg["calcium_ed"]["n"])
    tmp4 = arg["Ca_cyt"]**arg["calcium_ed"]["n_a"] / (arg["Ca_cyt"]**arg["calcium_ed"]["n_a"] + arg["calcium_ed"]["Ka"]**arg["calcium_ed"]["n_a"])

    I = arg["mito_membrane"]["Am"] * 2 * U * F * tmp1 * (tmp2 + tmp3 * tmp4)
    return I


def getIca_na_pump(arg):

    Keq = np.exp(-0.001 * arg["mito_membrane"]["Vmm"] * F / R / T)
    tmp1 = arg["Ca_mit"] / (arg["Ca_mit"] + arg["ca_na_pump"]["Km_ca"])
    tmp2 = arg["Na_cyt"]**arg["ca_na_pump"]["n"] / (arg["Na_cyt"]**arg["ca_na_pump"]["n"] + arg["ca_na_pump"]["Km_na"]**arg["ca_na_pump"]["n"])
    tmp3 = arg["Ca_mit"] * arg["Na_cyt"]**arg["ca_na_pump"]["n_Na"] - arg["Ca_cyt"] * arg["Na_mit"]**arg["ca_na_pump"]["n_Na"] / Keq

    I = tmp1 * tmp2 * tmp3

    return I

def getIca_h_pump(arg):
    Keq = np.exp(-0.001 * arg["mito_membrane"]["Vmm"] * F / R / T)

    tmp1 =  arg["Ca_mit"] / (arg["Ca_mit"] + arg["ca_h_pump"]["Km_ca"])
    tmp2 =  arg["Ca_mit"] * arg["H+_cyt"]**arg["ca_h_pump"]["n_H"] - arg["Ca_cyt"] * arg["H+_mit"]**arg["ca_h_pump"]["n_H"] / Keq

    I = tmp1 * tmp2

    return I

def getVcomplex1(arg):
    U = arg["mito_membrane"]["Vmm"] * F / (1000 * R * T)
    Keq = np.exp(2*arg["mito_membrane"]["Em_N"] + 2*arg["mito_membrane"]["Em_Q"]  + 4*U ) * (arg["mito_membrane"]["h_mit"]/arg["mito_membrane"]["h_cyt"])**4
    V = arg["complex1"]["Vmax"] * (arg["nadh_mit"] * arg["Q"] - arg["nad_mit"] * arg["QH2"] / Keq)

    return V

def getVcomplex3(arg):
    U = arg["mito_membrane"]["Vmm"] * F / (1000 * R * T)

    tmp = (arg["mito_membrane"]["h_mit"]/arg["mito_membrane"]["h_cyt"])**4
    Keq = np.exp(-2 * arg["mito_membrane"]["Em_Q"] + 2 * arg["mito_membrane"]["Em_cytc"] + 2*U) * tmp

    n =  arg["complex3"]["n"]
    V  = arg["complex3"]["Vmax"]*(arg["QH2"] * arg["cytc_ox"]**n - arg["Q"] * arg["cytc_red"]**n  / Keq)
    return V

def getVcomplex4(arg):

    tmp1 = arg["cytc_red"]**arg["complex4"]["n"] / (arg["cytc_red"]**arg["complex4"]["n"] + arg["complex4"]["Km_cytc"]**arg["complex4"]["n"])
    tmp2 = arg["O2_mit"] / (arg["O2_mit"] + arg["complex4"]["Km_O2"])
    tmp3 = (np.exp(-0.001 * arg["complex4"]["dGh"] * F / R / T  ))**2

    V = arg["complex4"]["Vmax"] * tmp1*tmp2*tmp3

    return V

def getVpyr_exchanger(arg):

    tmp1 = arg["pyr_cyt"] * arg["H+_cyt"] - arg["pyr_mit"] * arg["H+_mit"]
    tmp2 = 1 + arg["pyr_cyt"]/arg["pyr_exchanger"]["Km_pyr_cit"]
    tmp3 = 1 + arg["pyr_mit"]/arg["pyr_exchanger"]["Km_pyr_mit"]

    V = arg["pyr_exchanger"]["Vmax"] * tmp1 / (tmp2 * tmp3)

    return V

def getVpyr_dehydrogenase_complex_stage1(arg):

    tmp1 = 1 + arg["pyr_dehyd_comp"]["Amax_Ca"] * arg["Ca_mit"] / (arg["Ca_mit"] + arg["pyr_dehyd_comp"]["Ka_Ca"] )
    tmp2 = arg["pyr_mit"] / (arg["pyr_mit"] + arg["pyr_dehyd_comp"]["Km_pyr"] )
    tmp3 = arg["fad_pdhg"] / (arg["fad_pdhg"] + arg["pyr_dehyd_comp"]["Km_fad"])
    tmp4 = arg["CoA"] / (arg["CoA"] + arg["pyr_dehyd_comp"]["Km_CoA"]*(1 + arg["ACoA"] / arg["pyr_dehyd_comp"]["Ki_AcoA"]))
    pyr_dehyd_compACoA = arg["pyr_dehyd_comp"]["Vmax_pdhc_fad"] * tmp1 * tmp2 * tmp3 * tmp4

    return pyr_dehyd_compACoA

def getVpyr_dehydrogenase_complex_stage2(arg):
    Keq = np.exp( 2*(arg["pyr_dehyd_comp"]["Em_fad"] + arg["pyr_dehyd_comp"]["Em_nad"])*F*0.001 / R / T  )
    tmp5 = arg["fadh2_pdhg"] * arg["nad_mit"] - arg["fad_pdhg"]*arg["nadh_mit"] / Keq
    pyr_dehyd_compFad = arg["pyr_dehyd_comp"]["Vmax_pdhc_nad"] * tmp5 / (arg["nad_mit"] + arg["pyr_dehyd_comp"]["Km_nad"])

    return pyr_dehyd_compFad

def getVcitratesyntase(arg):

    tmp1 = arg["oa_mit"] / ( arg["oa_mit"] + arg["citrate_syntase"]["Km_oxa"]*(1 + arg["citr"]/arg["citrate_syntase"]["Ki_cit"]))
    tmp2 = arg["ACoA"] / (arg["ACoA"] + arg["citrate_syntase"]["Km_accoa"]*(1 + arg["CoA"]/arg["citrate_syntase"]["Ki_CoA"]))

    V = arg["citrate_syntase"]["Vmax"] * tmp1 * tmp2
    return V

def getVaconitase(arg):

    tmp1 = arg["citr"] - arg["isocitr"] /  arg["aconitase"]["Keq"]
    tmp2 = 1 + arg["citr"]/arg["aconitase"]["Km_cit"] + arg["isocitr"]/arg["aconitase"]["Km_isocit"]

    V = arg["aconitase"]["Vmax"] * tmp1 / tmp2
    return V

def getVisocit_dehydrogenase(arg):

    Km_isocitr = arg["isocit_dehydr"]["Km1_isocit"] / (1 + (arg["Ca_mit"]/arg["isocit_dehydr"]["Ka_Ca"])**arg["isocit_dehydr"]["n_Ca"]) + arg["isocit_dehydr"]["Km2_isocit"]
    tmp1 = arg["isocitr"]**arg["isocit_dehydr"]["n_isocit"] / (arg["isocitr"]**arg["isocit_dehydr"]["n_isocit"] + Km_isocitr**arg["isocit_dehydr"]["n_isocit"])
    tmp2 = arg["nad_mit"] / (arg["nad_mit"] + arg["isocit_dehydr"]["Km_nad"] * (1 + arg["nadh_mit"]/arg["isocit_dehydr"]["Ki_nadh"]) )
    V = arg["isocit_dehydr"]["Vmax"] * tmp1 * tmp2
    return V


def getVakg_dehydrogenase_stage1(arg):

    Km = (arg["akg_dehydr"]["Km1"]/(1 + arg["Ca_mit"]/arg["akg_dehydr"]["Ki_ca"]) + arg["akg_dehydr"]["Km2"] ) * (1 + arg["nadh_mit"] / arg["akg_dehydr"]["Ki_nadh"])
    tmp1 = arg["akg_mit"] / (arg["akg_mit"] + Km)
    tmp2 = arg["fad"] / (arg["fad"] + arg["akg_dehydr"]["Km_fad"])
    tmp3 = arg["CoA"] / (arg["CoA"] + arg["akg_dehydr"]["Km_CoA"]*(1 + arg["sucCoA"]/arg["akg_dehydr"]["Km_SucCoA"]) )
    Vfad = arg["akg_dehydr"]["Vmax_fad"] * tmp1 * tmp2 * tmp3

    return Vfad

def getVakg_dehydrogenase_stage2(arg):

    Keq = np.exp(0.002 * F * (arg["akg_dehydr"]["Em_fad"] + arg["akg_dehydr"]["Em_nad"]) / R / T )
    tmp4 = arg["fadh2"] * arg["nad_mit"] - arg["fad"]*arg["nadh_mit"] / Keq
    tmp5 = arg["nad_mit"] + arg["akg_dehydr"]["Km_nad"] * (1 + arg["nadh_mit"]/arg["akg_dehydr"]["Ki_nadh"])
    Vnad = arg["akg_dehydr"]["Vmax_nad"] * tmp4 / tmp5

    return Vnad

def getVsucCoAsyntase(arg, mode="gtp"):

    if mode == "gtp":
        Km_sucCoA = arg["sucCoAsyntase"]["Km_sucCoA_G"]
        Km_suc = arg["sucCoAsyntase"]["Km_suc_G"]
        Km_CoA = arg["sucCoAsyntase"]["Km_CoA_G"]
        nuc3P = arg["gtp_mit"]
        nuc2P = arg["gtp_mit"]
        Km_nuc3P = arg["sucCoAsyntase"]["Km_gdp"]
        Km_nuc2P = arg["sucCoAsyntase"]["Km_gtp"]

    elif mode == "atp":
        Km_sucCoA = arg["sucCoAsyntase"]["Km_sucCoA"]

        Km_suc = arg["sucCoAsyntase"]["Km_suc"]
        Km_CoA = arg["sucCoAsyntase"]["Km_CoA"]

        nuc3P = arg["atp_mit"]
        nuc2P = arg["adp_mit"]
        Km_nuc3P = arg["sucCoAsyntase"]["Km_atp"]
        Km_nuc2P = arg["sucCoAsyntase"]["Km_atp"]

    tmp1 = arg["sucCoAsyntase"]["Amax_P"]**arg["sucCoAsyntase"]["n_P"] * arg["pi_mit"]**arg["sucCoAsyntase"]["n_P"]
    tmp1 /= (arg["pi_mit"]**arg["sucCoAsyntase"]["n_P"] + arg["sucCoAsyntase"]["Km_P"]**arg["sucCoAsyntase"]["n_P"])
    tmp1 += 1

    tmp2 = arg["sucCoA"] * arg["gdp_mit"] * arg["pi_mit"] - arg["suc"]*arg["CoA"]*arg["gtp_mit"]/arg["sucCoAsyntase"]["Keq"]

    vsuccoa = 1 + arg["sucCoA"] / Km_sucCoA
    vnuc3P = 1 + nuc3P/Km_nuc3P
    vp = 1 + arg["pi_mit"] / arg["sucCoAsyntase"]["Km_P"]

    Vsuc = 1 + arg["suc"] / Km_suc
    Vcoa = 1 + arg["CoA"] / Km_CoA
    vnuc2P = 1 + nuc2P / Km_nuc2P

    V = arg["sucCoAsyntase"]["Vmax"] * tmp1 * tmp2 / (vsuccoa*vnuc3P*vp + Vsuc*Vcoa*vnuc2P - 1)

    return V


def getVsuc_dehydrydrogenase_stage1(arg):

    Keq_succdh = 1.0 # !!!!!!! np.exp( (25-arg["suc_dehydr"]["Em_FAD-succdh"]) * F / R / T )
    tmp1 = arg["suc"] * arg["Q"] - arg["fum"] * arg["QH2"] / Keq_succdh
    tmp2 = arg["suc"] + arg["suc_dehydr"]["Km_suc"]*(1 + arg["mal_mit"]/arg["suc_dehydr"]["Ki_mal"])
    v_succdh_fad = arg["suc_dehydr"]["Vmax_succdh"] * tmp1 / tmp2
    return v_succdh_fad

def getVsuc_dehydrydrogenase_stage2(arg):

    Keq_pdhc_fad_nad = 1.0 # !!! np.exp(-arg["suc_dehydr"]["Em_FAD-succdh"] * F / R / T) ## !!!!!! add - before Em
    tmp3 = arg["fadh2"] * arg["nad_mit"] - arg["fad"] * arg["nadh_mit"] / Keq_pdhc_fad_nad
    v_succdh = arg["suc_dehydr"]["Vmax_nadh"] * tmp3 / (arg["suc_dehydr"]["Km_nad"])
    return v_succdh


def getVfumarase(arg):

    tmp1 = arg["fum"] - arg["mal_mit"] / arg["fumarase"]["Keq"]
    tmp2 = 1 + arg["fum"]/arg["fumarase"]["Km_fum"] - arg["mal_mit"] / arg["fumarase"]["Km_mal"]

    V = arg["fumarase"]["Vmax"] * tmp1 / tmp2

    return V

def getVmal_dehydr(arg):

    tmp1 = arg["mal_mit"] * arg["nad_mit"] - arg["oa_mit"]*arg["nadh_mit"] / arg["mal_dehydr"]["Keq"]

    tmp2 = (1 + arg["mal_mit"] / arg["mal_dehydr"]["Km_mal"]) * (1 + arg["nad_mit"] / arg["mal_dehydr"]["Km_nad"])
    tmp3 = (1 + arg["oa_mit"] / arg["mal_dehydr"]["Km_oa"] ) * (1 + arg["nadh_mit"] / arg["mal_dehydr"]["Km_nadh"])

    V = arg["mal_dehydr"]["Vmax"] * tmp1 / (tmp2 + tmp3 - 1)

    return V



