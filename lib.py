
import numpy as np
from scipy import constants as const
from scipy.constants import physical_constants

F = physical_constants["Faraday constant"][0]
R = const.R
T = 310 # 37 grad Celcium


class Enzyme:

    def __init__(self):
        pass

    def update(self):
        pass

    def print_reag(self, metabilites):
        attrs = dir(self)

        print("#################################################")

        enzyme_name = str(self)
        # enzyme_name = enzyme_name.split("'")[1]
        print(enzyme_name)
        for a in attrs:
            if "idx" in a:
                v = self.__getattribute__(a)
                try:
                    idx = int(v)
                    print(metabilites[idx]["full"])
                except TypeError:
                    pass
                except ValueError:
                    pass


class Glucose_diffusion(Enzyme):
    def __init__(self, glc_ext, params):
        self.glc_ext_idx = glc_ext

        self.env_glc_level = params["env_glc_level"]
        self.D = params["D"]

    def update(self, metabolites, dydt):
        glc_ext = metabolites[self.glc_ext_idx]

        V = self.D * (self.env_glc_level - glc_ext)
        dydt[self.glc_ext_idx] += V
        return dydt

class Lactate_diffusion(Enzyme):
    def __init__(self, lac_ext, params):
        self.lac_ext_idx = lac_ext

        self.env_lac_level = params["env_lac_level"]
        self.D = params["D"]

    def update(self, metabolites, dydt):
        lac_ext = metabolites[self.lac_ext_idx]

        V = self.D * (self.env_lac_level - lac_ext)
        dydt[self.lac_ext_idx] += V
        return dydt

class Pyruvate_diffusion(Enzyme):
    def __init__(self, pyr_ext, params):
        self.pyr_ext_idx = pyr_ext

        self.env_pyr_level = params["env_pyr_level"]
        self.D = params["D"]

    def update(self, metabolites, dydt):
        lac_ext = metabolites[self.pyr_ext_idx]

        V = self.D * (self.env_pyr_level - lac_ext)
        dydt[self.pyr_ext_idx] += V
        return dydt

class Oxigen_diffusion(Enzyme):
    def __init__(self, o2_mit, params):
        self.o2_mit_idx = o2_mit

        self.env_o2_level = params["env_o2_level"]
        self.D = params["D"]

    def update(self, metabolites, dydt):
        o2_mit = metabolites[self.o2_mit_idx]

        V = self.D * (self.env_o2_level - o2_mit)
        dydt[self.o2_mit_idx] += V
        return dydt

########################################################################################################################
class GlucoseTransporter(Enzyme):

    def __init__(self, glc_ext, glc_cyt, params):
        self.glc_ext_idx = glc_ext
        self.glc_cyt_idx = glc_cyt
        self.Km_glc_cyt = params["Km_glc_cyt"]
        self.Km_glc_ext = params["Km_glc_ext"]
        self.Vmax = params["Vmax"]
        self.Volume_extracellular2cell = params["Volume_extracellular2cell"]

    def update(self, metabolites, dydt):

        glc_ext = metabolites[self.glc_ext_idx]
        glc_cyt = metabolites[self.glc_cyt_idx]

        glc_diff = glc_ext - glc_cyt

        glc_cyt_Km_ratio = glc_cyt / self.Km_glc_cyt
        glc_ext_Km_ratio = glc_ext / self.Km_glc_ext

        V = self.Vmax * glc_diff / (1 + glc_cyt_Km_ratio + glc_ext_Km_ratio)

        dydt[self.glc_ext_idx] -= V / self.Volume_extracellular2cell
        dydt[self.glc_cyt_idx] += V

        return dydt
########################################################################################################

class Hexokinase(Enzyme):

    def __init__(self, glc_cyt, atp_cyt, glc6p, adp_cyt, params):

        self.glc_cyt_idx = glc_cyt
        self.atp_cyt_idx = atp_cyt
        self.glc6p_idx = glc6p
        self.adp_cyt_idx = adp_cyt

        self.Vmax = params["Vmax"]
        self.Km_glc = params["Km_glc"]
        self.Km_atp = params["Km_atp"]
        self.Ki_atp = params["Ki_atp"]
        self.Ki_glc6p = params["Ki_glc6p"]

    def update(self, metabolites, dydt):
        glc_cyt = metabolites[self.glc_cyt_idx]
        atp_cyt = metabolites[self.atp_cyt_idx]
        glc6p = metabolites[self.glc6p_idx]

        # print(atp_cyt)

        Vglc = glc_cyt / (glc_cyt + self.Km_glc)
        Vatp = atp_cyt / (atp_cyt + self.Km_atp * (1 + glc6p / self.Ki_atp) )
        # atp_cyt / (atp_cyt + self.Km_atp * (1 + atp_cyt / self.Ki_atp) ) #

        # Inh = 1 - glc6p / self.Ki_glc6p
        V = self.Vmax * Vglc * Vatp #* Inh

        dydt[self.glc_cyt_idx] -= V
        dydt[self.atp_cyt_idx] -= V

        dydt[self.glc6p_idx] += V
        dydt[self.adp_cyt_idx] += V

        return dydt
#######################################################################################################

class Glucose6phosphate_isomerase(Enzyme):

    def __init__(self, glc6p, fru6p, params):

        self.glc6p_idx = glc6p
        self.fru6p_idx = fru6p

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_glc6p = params["Km_glc6p"]
        self.Km_fru6p = params["Km_fru6p"]

    def update(self, metabolites, dydt):
        glc6p = metabolites[self.glc6p_idx]
        fru6p = metabolites[self.fru6p_idx]

        tmp1 = glc6p - fru6p / self.Keq
        tmp2 = 1 + glc6p / self.Km_glc6p + fru6p / self.Km_fru6p
        V = self.Vmax * tmp1 / tmp2

        dydt[self.glc6p_idx] -= V
        dydt[self.fru6p_idx] += V

        return dydt
#################################################################################################

class Phosphofructokinase_type1(Enzyme):

    def __init__(self, fru6p, atp_cyt, fru16p, adp_cyt, pi_cyt, fru26p, params):

        self.fru6p_idx = fru6p
        self.atp_cyt_idx = atp_cyt
        self.fru16p_idx = fru16p
        self.adp_cyt_idx = adp_cyt
        self.pi_cyt_idx = pi_cyt
        self.fru26p_idx = fru26p

        self.Vmax = params["Vmax"]
        self.Km_fru6p = params["Km_fru6p"]
        self.Km_atp = params["Km_atp"]
        self.n = params["n"]
        self.Ki_atp = params["Ki_atp"]
        self.K0 = params["K0"]
        self.Ka_fru26p = params["Ka_fru26p"]
        self.n_fru26p = params["n_fru26p"]
        self.Ka_fru26p = params["Ka_fru26p"]

    def update(self, metabolites, dydt):

        fru6p = metabolites[self.fru6p_idx]
        atp_cyt = metabolites[self.atp_cyt_idx]
        fru26p = metabolites[self.fru16p_idx]

        fru__n = fru26p**self.n_fru26p
        # print (fru__n)
        tmp1 = fru__n / (fru__n + self.Ka_fru26p**self.n_fru26p )

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
class Fructose16_bisphosphatase(Enzyme):
    def __init__(self, fru16p, fru6p, pi_cyt, params):
        self.fru16p_idx = fru16p
        self.fru6p_idx = fru6p
        self.pi_cyt_idx = pi_cyt

        self.Vmax = params["Vmax"]
        self.Km = params["Km"]

    def update(self, metabolites, dydt):

        fru16bp = metabolites[self.fru16p_idx]

        V = self.Vmax * fru16bp / (fru16bp + self.Km)

        dydt[self.fru16p_idx] -= V
        dydt[self.fru6p_idx] += V
        dydt[self.pi_cyt_idx] += V

        return dydt
#####################################################################################################

class Phosphofructokinase_type2(Enzyme):

    def __init__(self, fru6p, atp_cyt, fru26p, adp_cyt, amp, params):
        self.fru6p_idx = fru6p
        self.atp_cyt_idx = atp_cyt
        self.fru26p_idx = fru26p
        self.adp_cyt_idx = adp_cyt
        self.amp_idx = amp

        self.Vmax = params["Vmax"]
        self.Km_fru6p = params["Km_fru6p"]
        self.Km_atp = params["Km_atp"]
        self.Ka_amp = params["Ka_amp"]
        self.Ka_adp = params["Ka_adp"]

    def update(self, metabolites, dydt):
        fru6p = metabolites[self.fru6p_idx]
        atp_cyt = metabolites[self.atp_cyt_idx]
        amp_cyt = metabolites[self.amp_idx]
        adp_cyt = metabolites[self.adp_cyt_idx]

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

class Fructose26_bisphosphatase(Enzyme):

    def __init__(self, fru26p, fru6p, pi_cyt, params):
        self.fru26p_idx = fru26p
        self.fru6p_idx = fru6p
        self.pi_cyt_idx = pi_cyt

        self.Vmax = params["Vmax"]
        self.Km = params["Km"]
        self.Ki_fru6p = params["Ki_fru6p"]




    def update(self, metabolites, dydt):
        fru26p = metabolites[self.fru26p_idx]
        fru6p = metabolites[self.fru6p_idx]

        tmp = 1 + fru6p / self.Ki_fru6p
        V = self.Vmax * fru26p / (fru26p + self.Km * tmp)

        dydt[self.fru26p_idx] -= V

        dydt[self.fru6p_idx] += V
        dydt[self.pi_cyt_idx] += V

        return dydt

#######################################################################################################

class Aldolase(Enzyme):

    def __init__(self, fru16p, grap, dhap, params):
        self.fru16p_idx = fru16p
        self.grap_idx = grap
        self.dhap_idx = dhap

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_fru16p = params["Km_fru16p"]
        self.Km_grap = params["Km_grap"]
        self.Km_dhap = params["Km_dhap"]



    def update(self, metabolites, dydt):
        fru16p = metabolites[self.fru16p_idx]
        grap = metabolites[self.grap_idx]
        dhap = metabolites[self.dhap_idx]


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

class Triosophosphate_isomerase(Enzyme):

    def __init__(self, grap, dhap, params):
        self.grap_idx = grap
        self.dhap_idx = dhap
        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_grap = params["Km_grap"]
        self.Km_dhap = params["Km_dhap"]

    def update(self, metabolites, dydt):
        grap = metabolites[self.grap_idx]
        dhap = metabolites[self.dhap_idx]

        tmp1 = dhap - grap / self.Keq # !!!! поменяны местами метаболиты
        tmp2 = 1 + dhap / self.Km_dhap + grap / self.Km_grap

        V = self.Vmax * tmp1 / tmp2

        dydt[self.dhap_idx] -= V
        dydt[self.grap_idx] += V

        return dydt

######################################################################################################

class Glyceraldehyde_3_phosphate_dehydrogenase(Enzyme):

    def __init__(self, grap, pi, nad, bpg13, nadh, params):
        self.grap_idx = grap
        self.pi_idx = pi
        self.nad_idx = nad
        self.bpg13_idx = bpg13
        self.nadh_idx = nadh

        self.Vmax = params["Vmax"]
        self.Km_nad = params["Km_nad"]
        self.Km_grap = params["Km_grap"]
        self.Km_pi = params["Km_pi"]
        self.Km_nadh = params["Km_nadh"]
        self.Km_bpg13 = params["Km_bpg13"]
        self.Keq = params["Keq"]

    def update(self, metabolites, dydt):

        grap = metabolites[self.grap_idx]
        pi = metabolites[self.pi_idx]
        nad = metabolites[self.nad_idx]
        bpg13 = metabolites[self.bpg13_idx]
        nadh = metabolites[self.nadh_idx]


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

class Phosphoglycerate_kinase(Enzyme):

    def __init__(self, bpg13, adp, pg3, atp, params):

        self.bpg13_idx = bpg13
        self.adp_idx = adp
        self.atp_idx = atp
        self.pg3_idx = pg3


        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_bpg13 = params["Km_bpg13"]
        self.Km_adp = params["Km_adp"]
        self.Km_pg3 = params["Km_pg3"]
        self.Km_atp = params["Km_atp"]

    def update(self, metabolites, dydt):

        bpg13 = metabolites[self.bpg13_idx]
        adp_cyt = metabolites[self.adp_idx ]
        pg3 = metabolites[self.pg3_idx]
        atp_cyt = metabolites[self.atp_idx ]

        tmp1 = bpg13 * adp_cyt - pg3 * atp_cyt / self.Keq
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

class Phosphoglycerate_mutase(Enzyme):

    def __init__(self, pg3, pg2, params):

        self.pg3_idx = pg3
        self.pg2_idx = pg2

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_pg2 = params["Km_pg2"]
        self.Km_pg3 = params["Km_pg3"]


    def update(self, metabolites, dydt):
        pg2 = metabolites[self.pg2_idx]
        pg3 = metabolites[self.pg3_idx]

        tmp1 = pg3 - pg2 / self.Keq  #  поменяны местами метаболиты !!!!!
        tmp2 = 1 + pg3 / self.Km_pg3
        tmp3 = 1 + pg2 / self.Km_pg2

        V = self.Vmax * tmp1 / (tmp2 + tmp3 - 1)

        dydt[self.pg3_idx] -= V

        dydt[self.pg2_idx] += V

        return dydt

#######################################################################################################

class Enolase(Enzyme):

    def __init__(self, pg2, pep, params):
        self.pg2_idx = pg2
        self.pep_idx = pep

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_pg2 = params["Km_pg2"]
        self.Km_pep = params["Km_pep"]

    def update(self, metabolites, dydt):
        pg2 = metabolites[self.pg2_idx]
        pep = metabolites[self.pep_idx]

        tmp1 = pg2 - pep / self.Keq
        tmp2 = 1 + pg2 / self.Km_pg2
        tmp3 = 1 + pep / self.Km_pep

        V = self.Vmax * tmp1 / (tmp2 + tmp3 - 1)

        dydt[self.pg2_idx] -= V
        dydt[self.pep_idx] += V

        return dydt
######################################################################################################

class Pyruvate_kinase(Enzyme):

    def __init__(self, pep, adp, pyr, atp, params):

        self.pep_idx = pep
        self.adp_cyt_idx = adp
        self.atp_cyt_idx = atp
        self.pyr_cyt_idx = pyr

        self.Vmax = params["Vmax"]
        self.Km_pep = params["Km_pep"]
        self.Km_adp = params["Km_adp"]
        self.Ki_atp = params["Ki_atp"]

    def update(self, metabolites, dydt):

        pep = metabolites[self.pep_idx]
        adp = metabolites[self.adp_cyt_idx]
        atp = metabolites[self.atp_cyt_idx]

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
class Lactate_dehydrogenase(Enzyme):

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


    def update(self, metabolites, dydt):

        pyr = metabolites[self.pyr_idx]
        nadh = metabolites[self.nadh_idx]
        nad = metabolites[self.nad_idx]
        lac = metabolites[self.lac_idx]

        tmp1 = pyr * nadh - lac * nad / self.Keq
        tmp2 = 1 + pyr / self.Km_pyr
        tmp3 = 1 + nadh / self.Km_nadh
        tmp4 = 1 + lac / self.Km_lac
        tmp5 = 1 + nad / self.Km_nad

        V = self.Vmax * tmp1 / (tmp2*tmp3 + tmp4*tmp5 - 1)

        dydt[self.pyr_idx] -= V
        dydt[self.nadh_idx] -= V

        dydt[self.lac_idx] += V
        dydt[self.nad_idx] += V

        return dydt

########################################################################################################

class Monocarboxilate_transporter(Enzyme):

    def __init__(self, lac_ext, lac_cyt, params):
        self.lac_ext_idx = lac_ext
        self.lac_cyt_idx = lac_cyt

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_lac_cyt = params["Km_lac_cyt"]
        self.Km_lac_ext = params["Km_lac_ext"]
        self.Volume_extracellular2cell = params["Volume_extracellular2cell"]



    def update(self, metabolites, dydt):
        lac_ext = metabolites[self.lac_ext_idx]
        lac_cyt = metabolites[self.lac_cyt_idx]

        tmp1 = lac_cyt - lac_ext / self.Keq
        tmp2 = 1 + lac_cyt / self.Km_lac_cyt
        tmp3 = 1 + lac_ext / self.Km_lac_ext

        V = self.Vmax * tmp1 / (tmp2 + tmp3 - 1)

        dydt[self.lac_cyt_idx] -= V
        dydt[self.lac_ext_idx] += V / self.Volume_extracellular2cell

        return dydt

#########################################################################################################

class Creatine_kinase(Enzyme):

    def __init__(self, atp, cr, adp, crp, params):

        self.atp_idx = atp
        self.cr_idx = cr
        self.adp_idx = adp
        self.crp_idx = crp

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]


    def update(self, metabolites, dydt):

        atp = metabolites[self.atp_idx]
        adp = metabolites[self.adp_idx]
        cr = metabolites[self.cr_idx]
        crp = metabolites[self.crp_idx]


        tmp = atp * cr / self.Keq
        V = self.Vmax * (adp * crp - tmp)

        dydt[self.atp_idx] -= V
        dydt[self.cr_idx] -= V

        dydt[self.crp_idx] += V
        dydt[self.adp_idx] += V

        return dydt
#########################################################################################################


class Malate_dehydrogenase(Enzyme):

    def __init__(self, mal, oa, nad, nadh, params):

        self.mal_idx = mal
        self.oa_idx = oa
        self.nad_idx = nad
        self.nadh_idx = nadh

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_nad = params["Km_nad"]
        self.Km_nadh = params["Km_nadh"]
        self.Km_mal = params["Km_mal"]
        self.Km_oa = params["Km_oa"]

    def update(self, metabolites, dydt):

        mal = metabolites[self.mal_idx]
        oa = metabolites[self.oa_idx]
        nad = metabolites[self.nad_idx]
        nadh = metabolites[self.nadh_idx]

        tmp1 = mal * nad - oa * nadh / self.Keq
        tmp2 = 1 + mal / self.Km_mal
        tmp3 = 1 + nad / self.Km_nad
        tmp4 = 1 + oa / self.Km_oa
        tmp5 = 1 + nadh / self.Km_nadh

        V = self.Vmax * tmp1 / (tmp2*tmp3 + tmp4*tmp5 - 1)

        dydt[self.mal_idx] -= V
        dydt[self.nad_idx] -= V

        dydt[self.nadh_idx] += V
        dydt[self.oa_idx] += V

        return dydt
########################################################################################################

class Aspartate_aminotransferase(Enzyme):

    def __init__(self, asp, akg, oa, glu, params):
        self.oa_idx = oa
        self.asp_idx = asp
        self.akg_idx = akg
        self.glu_idx = glu

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]

    def update(self, metabolites, dydt):

        asp = metabolites[self.asp_idx]
        oa = metabolites[self.oa_idx]
        akg = metabolites[self.akg_idx]
        glu = metabolites[self.glu_idx]

        tmp = oa * glu / self.Keq
        V = self.Vmax * (asp * akg - tmp)

        dydt[self.asp_idx] -= V
        dydt[self.akg_idx] -= V

        dydt[self.akg_idx] += V
        dydt[self.glu_idx] += V

        return dydt
#############################################################################################################

class Aspartate_glutamate_carrier(Enzyme):

    def __init__(self, asp_mit, glu_cyt, h_cyt, asp_cyt, glu_mit, h_mit, Vmm, params):
        self.asp_mit_idx = asp_mit
        self.glu_cyt_idx = glu_cyt
        self.h_cyt_idx = h_cyt

        self.asp_cyt_idx = asp_cyt
        self.glu_mit_idx = glu_mit
        self.h_mit_idx = h_mit
        self.Vmm_idx = Vmm

        self.Vmax = params["Vmax"]
        self.Km_asp_mit = params["Km_asp_mit"]
        self.Km_glu_cyt = params["Km_glu_cyt"]

        self.Km_asp_cyt = params["Km_asp_cyt"]
        self.Km_glu_mit = params["Km_glu_mit"]

    def update(self, metabolites, dydt):
        asp_mit = metabolites[self.asp_mit_idx]
        glu_cyt = metabolites[self.glu_cyt_idx]
        h_cyt = metabolites[self.h_cyt_idx]

        asp_cyt = metabolites[self.asp_cyt_idx]
        glu_mit = metabolites[self.glu_mit_idx]
        h_mit = metabolites[self.h_mit_idx]
        Vmm = metabolites[self.Vmm_idx]

        dG = -Vmm + 1000 * R * T / F * np.log( h_cyt / h_mit)

        Keq = np.exp(F * dG/(1000 * R * T) )

        tmp1 = asp_mit * glu_cyt - asp_cyt * glu_mit / Keq
        tmp2 = (asp_mit + self.Km_asp_mit) * (glu_cyt +  self.Km_glu_cyt)
        tmp3 = (asp_cyt + self.Km_asp_cyt) * (glu_mit +  self.Km_glu_mit)

        V = self.Vmax * tmp1 / (tmp2 + tmp3)

        dydt[self.asp_mit_idx] -= V
        dydt[self.glu_cyt_idx] -= V
        dydt[self.h_cyt_idx] -= V

        dydt[self.h_mit_idx] += V
        dydt[self.asp_cyt_idx] += V
        dydt[self.glu_mit_idx] += V

        return dydt

#################################################################################################################

class Malate_alphaketoglutarate_carrier(Enzyme):

    def __init__(self, mal_cyt, akg_mit, mal_mit, akg_cyt, params):

        self.mal_cyt_idx = mal_cyt
        self.akg_mit_idx = akg_mit
        self.mal_mit_idx = mal_mit
        self.akg_cyt_idx = akg_cyt

        self.Vmax = params["Vmax"]
        self.Km_mal_cyt = params["Km_mal_cyt"]
        self.Km_akg_mit = params["Km_akg_mit"]
        self.Km_mal_mit = params["Km_mal_mit"]
        self.Km_akg_cyt = params["Km_akg_cyt"]

    def update(self, metabolites, dydt):

        mal_cyt = metabolites[self.mal_cyt_idx]
        akg_mit = metabolites[self.akg_mit_idx]
        mal_mit = metabolites[self.mal_mit_idx]
        akg_cyt = metabolites[self.akg_cyt_idx]

        tmp1 = mal_cyt * akg_mit - mal_mit * akg_cyt
        tmp2 = ( mal_cyt + self.Km_mal_cyt) * (akg_mit + self.Km_akg_mit)
        tmp3 = ( mal_mit + self.Km_mal_mit) * (akg_cyt + self.Km_akg_cyt)

        V = self.Vmax * tmp1 / (tmp2 + tmp3)

        dydt[self.mal_cyt_idx] -= V
        dydt[self.akg_mit_idx] -= V

        dydt[self.mal_mit_idx] += V
        dydt[self.akg_cyt_idx] += V

        return dydt

###############################################################################################################

class Glycerol_3phosphate_dehydrogenase_cytosolic(Enzyme):

    def __init__(self, dhap, nadh, g3p, nad, params):
        self.dhap_idx = dhap
        self.nadh_idx = nadh
        self.nad_idx = nad
        self.g3p_idx = g3p

        self.Vmax = params["Vmax"]
        self.Keq = params["Keq"]
        self.Km_dhap = params["Km_dhap"]
        self.Km_nadh = params["Km_nadh"]
        self.Km_nad = params["Km_nad"]
        self.Km_g3p = params["Km_g3p"]

    def update(self, metabolites, dydt):
        dhap = metabolites[self.dhap_idx]
        nadh = metabolites[self.nadh_idx]
        nad = metabolites[self.nad_idx]
        g3p = metabolites[self.g3p_idx]

        tmp1 = dhap * nadh - g3p * nad / self.Keq
        tmp2 = 1 + dhap / self.Km_dhap
        tmp3 = 1 + nadh / self.Km_nadh
        tmp4 = 1 + g3p / self.Km_g3p
        tmp5 = 1 + nad / self.Km_nad

        V = self.Vmax * tmp1 / (tmp2*tmp3 + tmp4*tmp5 - 1)

        dydt[self.dhap_idx] -= V
        dydt[self.nadh_idx] -= V

        dydt[self.nad_idx] += V
        dydt[self.g3p_idx] += V

        return dydt


###########################################################################################################

class Glycerol_3phosphate_dehydrogenase_mitochondrial(Enzyme):

    def __init__(self, g3p, dhap, fad_g3dh, fadh2_g3dh, Q, QH2, params):
        self.dhap_idx = dhap
        self.fad_g3dh_idx = fad_g3dh
        self.fadh2_g3dh_idx = fadh2_g3dh
        self.g3p_idx = g3p
        self.q_idx = Q
        self.qh2_idx = QH2

        self.Vmax_g3pdh = params["Vmax_g3pdh"]
        self.Vmax_Q = params["Vmax_Q"]

        self.Em_FAD_g3p = params["Em_FAD_g3p"]
        self.Em_Q = params["Em_Q"]
        self.Em_dhap_g3p = params["Em_dhap_g3p"]

        self.Km_dhap = params["Km_dhap"]
        self.Km_g3p = params["Km_g3p"]

    def update(self, metabolites, dydt):

        dhap = metabolites[self.dhap_idx]
        g3p = metabolites[self.g3p_idx]

        fad = metabolites[self.fad_g3dh_idx]
        fadh2 = metabolites[self.fadh2_g3dh_idx]

        Q = metabolites[self.q_idx]
        QH2 = metabolites[self.qh2_idx]

        Keq_g3pdh = np.exp( (2*self.Em_dhap_g3p - 2*self.Em_FAD_g3p)* F / 1000 / R / T )

        tmp1 = g3p * fad - dhap * fadh2 / Keq_g3pdh
        tmp2 = 1 + dhap / self.Km_dhap
        tmp3 = 1 + g3p / self.Km_g3p

        Vg3pdh = self.Vmax_g3pdh * tmp1 / (tmp2 + tmp3)

        dydt[self.g3p_idx] -= Vg3pdh
        dydt[self.fad_g3dh_idx] -= Vg3pdh

        dydt[self.dhap_idx] += Vg3pdh
        dydt[self.fadh2_g3dh_idx] += Vg3pdh

        Keq_fad_Q = np.exp( (2 * self.Em_FAD_g3p + 2 * self.Em_Q ) * F / 1000 / R / T )
        Vqh2 = self.Vmax_Q * ( fadh2 * Q - fad * QH2 / Keq_fad_Q  )

        dydt[self.fadh2_g3dh_idx] -= Vqh2
        dydt[self.q_idx] -= Vqh2

        dydt[self.fad_g3dh_idx] += Vqh2
        dydt[self.qh2_idx] += Vqh2


        return dydt

#########################################################################################################

class ATP_synthetase(Enzyme):

    def __init__(self, atp, adp, pi, h_cyt, h_mit, Vmm, params):
        self.atp_idx = atp
        self.adp_idx = adp
        self.pi_idx = pi
        self.h_cyt_idx = h_cyt
        self.h_mit_idx = h_mit
        self.Vmm_idx = Vmm

        self.n = params["n"]
        self.k = params["k"]
        self.dG0 = params["dG0"]


    def update(self, metabolites, dydt):
        atp = metabolites[self.atp_idx]
        adp = metabolites[self.adp_idx]
        pi = metabolites[self.pi_idx]
        h_cyt = metabolites[self.h_cyt_idx]
        h_mit = metabolites[self.h_mit_idx]
        Vmm = metabolites[self.Vmm_idx]

        dG = -Vmm + R * T * np.log( h_cyt / h_mit ) / (1000 * F)
        Vmax = 1.8 * 10**-16 * dG**self.n

        U = Vmm * F / (1000 * R * T)

        Keq = np.exp( self.dG0 /(R * T) - self.k * U ) * (( h_cyt / h_mit )**self.k)

        V = Vmax * (adp * pi - atp / Keq  )

        dydt[self.atp_idx] -= V
        dydt[self.adp_idx] += V
        dydt[self.pi_idx] += V

        return dydt

##########################################################################################################

class ATP_ADP_axchanger(Enzyme):

    def __init__(self, atp_mit, adp_cyt, adp_mit, atp_cyt, Vmm, params):
        self.atp_mit_idx = atp_mit
        self.adp_cyt_idx = adp_cyt
        self.adp_mit_idx = adp_mit
        self.atp_cyt_idx = atp_cyt
        self.Vmm_idx = Vmm

        self.Vmax = params["Vmax"]
        self.S_Vmm = params["S_Vmm"]
        self.Volume_cyt_mit = params["Volume_cyt_mit"]
        self.Cmm = params["Cmm"]


    def update(self, metabolites, dydt):

        atp_mit = metabolites[self.atp_mit_idx]
        adp_cyt = metabolites[self.adp_cyt_idx]
        adp_mit = metabolites[self.adp_mit_idx]
        atp_cyt = metabolites[self.atp_cyt_idx]
        Vmm = metabolites[self.Vmm_idx]

        U = Vmm * F / (1000 * R * T)

        tmp1 = 1 - np.exp(U) * atp_cyt * adp_mit / (adp_cyt * atp_mit)
        tmp2 = 1 + atp_cyt / adp_cyt  * np.exp( self.S_Vmm * U)
        tmp3 = 1 + atp_mit / adp_mit

        V = self.Vmax * tmp1 / (tmp2 * tmp3)

        dydt[self.atp_mit_idx] -= V / self.Volume_cyt_mit
        dydt[self.adp_cyt_idx] -= V / (1 - self.Volume_cyt_mit)

        dydt[self.adp_mit_idx] += V / self.Volume_cyt_mit
        dydt[self.atp_cyt_idx] += V / (1 - self.Volume_cyt_mit)

        dydt[self.Vmm_idx] += V * F / self.Cmm
        return dydt

############################################################################################################

class ATP_consumption(Enzyme):

    def __init__(self, atp, adp, pi, params):

        self.atp_idx = atp
        self.adp_idx = adp
        self.pi_idx = pi

        self.Vmax = params["Vmax"]
        self.Km_atp = params["Km_atp"]
        self.activation = params["activation"]


    def update(self, metabolites, dydt):
        atp = metabolites[self.atp_idx]

        V = self.Vmax * atp / ( atp + self.Km_atp) * (1 + self.activation)

        dydt[self.atp_idx] -= V
        dydt[self.adp_idx] += V
        dydt[self.pi_idx] += V

        return dydt
############################################################################################################

class Passive_efflux_ion(Enzyme):

    def __init__(self, ion_in, ion_out, Vmm, params):
        self.ion_in_idx = ion_in
        self.ion_out_idx = ion_out

        self.Vmm_idx = Vmm

        self.Am = params["Am"]
        self.P = params["P"]

        self.Cmm = params["Cmm"]
        self.Volume_cyt_mit = params["Volume_cyt_mit"]

    def update(self, metabolites, dydt):
        ion_in = metabolites[self.ion_in_idx]
        ion_out = metabolites[self.ion_out_idx]
        Vmm = metabolites[self.Vmm_idx]

        U = Vmm * F / (1000 * R * T)
        tmp = (ion_in - ion_out * np.exp(U)) / (1 - np.exp(U))
        I = self.P * U * F * tmp * self.Am

        dydt[self.ion_in_idx] -= I / F / (1 - self.Volume_cyt_mit)
        dydt[self.ion_out_idx] += I / F / self.Volume_cyt_mit


        dydt[self.Vmm_idx] += I / self.Cmm

        return dydt
############################################################################################################
class Pump(Enzyme):

    def __init__(self, ion_cyt, ion_mit, h_cyt, h_mit, params):

        self.ion_cyt_idx = ion_cyt
        self.ion_mit_idx = ion_mit
        self.h_cyt_idx = h_cyt
        self.h_mit_idx = h_mit

        self.Vmax = params["Vmax"]
        self.is_simport = params["is_simport"]
        self.Volume_cyt_mit = params["Volume_cyt_mit"]


    def update(self, metabolites, dydt):

        ion_cyt = metabolites[self.ion_cyt_idx]
        ion_mit = metabolites[self.ion_mit_idx]
        h_cyt = metabolites[self.h_cyt_idx]
        h_mit = metabolites[self.h_mit_idx]

        if self.is_simport:
            V = self.Vmax * (ion_cyt * h_cyt - ion_mit * h_mit)

            dydt[self.h_cyt_idx] -= V / (1 - self.Volume_cyt_mit)
            dydt[self.h_mit_idx] += V / self.Volume_cyt_mit
        else:
            V = self.Vmax * (ion_cyt * h_mit - ion_mit * h_cyt)

            dydt[self.h_mit_idx] -= V / self.Volume_cyt_mit
            dydt[self.h_cyt_idx] += V  / (1 - self.Volume_cyt_mit)

        dydt[self.ion_cyt_idx] -= V  / (1 - self.Volume_cyt_mit)
        dydt[self.ion_mit_idx] += V / self.Volume_cyt_mit

        return dydt

###########################################################################################################

class Calcium_effux(Enzyme):
    def __init__(self, ca_cyt, ca_mit, Vmm, params):

        self.ca_cyt_idx = ca_cyt
        self.ca_mit_idx = ca_mit
        self.Vmm_idx = Vmm

        self.P_RMC = params["P_RMC"]
        self.Ki_cacyt = params["Ki_cacyt"]
        self.P_Mcu = params["P_Mcu"]
        self.n = params["n"]
        self.n_a = params["n_a"]
        self.Ka = params["Ka"]
        self.Am = params["Am"]
        self.Km_Mcu = params["Km_Mcu"]
        self.Volume_cyt_mit = params["Volume_cyt_mit"]
        self.Cmm = params["Cmm"]

    def update(self, metabolites, dydt):

        ca_cyt = metabolites[self.ca_cyt_idx]
        ca_mit = metabolites[self.ca_mit_idx]
        Vmm = metabolites[self.Vmm_idx]

        U = Vmm * F / (1000 * R * T)
        tmp1 = ( ca_cyt - ca_mit * np.exp(2 * U) ) / (1 - np.exp(2 * U))

        tmp2 = self.P_RMC * (1 - ca_cyt/(ca_cyt + self.Ki_cacyt))
        tmp3 = self.P_Mcu * ca_cyt**self.n / ( ca_cyt**self.n + self.Km_Mcu**self.n)
        tmp4 = ca_cyt**self.n_a / (ca_cyt**self.n_a + self.Ka**self.n_a)

        I = self.Am * 2 * U * F * tmp1 * (tmp2 + tmp3 * tmp4)

        dydt[self.ca_cyt_idx] -= I * F / (1 - self.Volume_cyt_mit)
        dydt[self.ca_mit_idx] += I * F / self.Volume_cyt_mit

        dydt[self.Vmm_idx] += I / self.Cmm
        return dydt

###############################################################################################

class Ca_Na_pump(Enzyme):

    def __init__(self, ca_mit, na_cyt, ca_cyt, na_mit, Vmm, params):
        self.ca_mit_idx = ca_mit
        self.na_cyt_idx = na_cyt
        self.ca_cyt_idx = ca_cyt
        self.na_mit_idx = na_mit
        self.Vmm_idx = Vmm

        self.n_Na = params["n_Na"]
        self.n = params["n"]
        self.Km_ca = params["Km_ca"]
        self.Km_na = params["Km_na"]
        self.Volume_cyt_mit = params["Volume_cyt_mit"]
        self.Cmm = params["Cmm"]

    def update(self, metabolites, dydt):
        Vmm = metabolites[self.Vmm_idx]

        ca_mit = metabolites[self.ca_mit_idx]
        na_cyt = metabolites[self.na_cyt_idx]


        ca_cyt = metabolites[self.ca_cyt_idx]
        na_mit = metabolites[self.na_mit_idx]

        Keq = np.exp(-0.001 * Vmm * F / R / T)
        tmp1 = ca_mit / (ca_mit + self.Km_ca)
        tmp2 = na_cyt**self.n / (na_cyt**self.n + self.Km_na**self.n)
        tmp3 = ca_mit * na_cyt**self.n_Na - ca_cyt * na_mit**self.n_Na / Keq

        V = tmp1 * tmp2 * tmp3

        dydt[self.ca_mit_idx] -= V  / self.Volume_cyt_mit
        dydt[self.na_cyt_idx] -= 3 * V / (1 - self.Volume_cyt_mit)

        dydt[self.na_mit_idx] += 3 * V  / self.Volume_cyt_mit
        dydt[self.ca_cyt_idx] += 3 * V / (1 - self.Volume_cyt_mit)

        dydt[self.Vmm_idx] += V * F / self.Cmm

        return dydt
###############################################################################################

class Ca_H_pump(Enzyme):

    def __init__(self,  ca_mit, h_cyt, ca_cyt, h_mit, Vmm, params):
        
        self.ca_mit_idx = ca_mit
        self.h_cyt_idx = h_cyt
        self.ca_cyt_idx = ca_cyt
        self.h_mit_idx = h_mit
        self.Vmm_idx = Vmm

        self.Km_ca = params["Km_ca"]
        self.n_H = params["n_H"]

        self.Volume_cyt_mit = params["Volume_cyt_mit"]
        self.Cmm = params["Cmm"]


    def update(self, metabolites, dydt):

        Vmm = metabolites[self.Vmm_idx]

        ca_mit = metabolites[self.ca_mit_idx]
        h_cyt = metabolites[self.h_cyt_idx]

        ca_cyt = metabolites[self.ca_cyt_idx]
        h_mit = metabolites[self.h_mit_idx]

        Keq = np.exp(-0.001 * Vmm * F / R / T)
        tmp1 = ca_mit / ( ca_mit + self.Km_ca)
        tmp2 =  ca_mit * h_cyt**self.n_H - ca_cyt * h_mit**self.n_H / Keq

        V = tmp1 * tmp2

        dydt[self.ca_mit_idx] -= V  / self.Volume_cyt_mit
        dydt[self.h_cyt_idx] -= 3 * V / (1 - self.Volume_cyt_mit)

        dydt[self.h_mit_idx] += 3 * V  / self.Volume_cyt_mit
        dydt[self.ca_cyt_idx] += 3 * V / (1 - self.Volume_cyt_mit)

        dydt[self.Vmm_idx] += V * F / self.Cmm

        return dydt
##############################################################################################

class Complex1(Enzyme):

    def __init__(self, h_cyt, h_mit, q, qh2, nad, nadh, Vmm, params):
        self.h_cyt_idx = h_cyt
        self.h_mit_idx = h_mit
        self.q_idx = q
        self.qh2_idx = qh2
        self.nad_idx = nad
        self.nadh_idx = nadh
        self.Vmm_idx = Vmm

        self.Em_N = params["Em_N"]
        self.Em_Q = params["Em_Q"]
        self.Vmax = params["Vmax"]

        self.Volume_cyt_mit = params["Volume_cyt_mit"]
        self.Cmm = params["Cmm"]

    def update(self, metabolites, dydt):

        Vmm = metabolites[self.Vmm_idx]
        h_cyt = metabolites[self.h_cyt_idx]
        h_mit = metabolites[self.h_mit_idx]
        Q = metabolites[self.q_idx]
        QH2 = metabolites[self.qh2_idx]
        nadh = metabolites[self.nadh_idx]
        nad = metabolites[self.nad_idx]

        U = Vmm * F / (1000 * R * T)
        Keq = np.exp(2 * self.Em_N + 2*self.Em_Q  + 4*U ) * (h_mit/h_cyt)**4
        V = self.Vmax * (nadh * Q - nad * QH2 / Keq)

        dydt[self.nadh_idx] -= V
        dydt[self.q_idx] -= V
        dydt[self.h_mit_idx] -= 4 * V / self.Volume_cyt_mit

        dydt[self.h_cyt_idx] += 4 * V  / (1 - self.Volume_cyt_mit)
        dydt[self.qh2_idx] += V
        dydt[self.nad_idx] += V

        dydt[self.Vmm_idx] -= 4 * V * F / self.Cmm  # !!!!!!

        return dydt
###################################################################################

class Complex3(Enzyme):

    def __init__(self, h_cyt, h_mit, q, qh2, cytc_ox, cytc_red, Vmm, params):
        self.h_cyt_idx = h_cyt
        self.h_mit_idx = h_mit
        self.q_idx = q
        self.qh2_idx = qh2
        self.cytc_ox_idx = cytc_ox
        self.cytc_red_idx = cytc_red
        self.Vmm_idx = Vmm

        self.Em_Q = params["Em_Q"]
        self.Em_cytc = params["Em_cytc"]
        self.n = params["n"]
        self.Vmax = params["Vmax"]
        self.Volume_cyt_mit = params["Volume_cyt_mit"]
        self.Cmm = params["Cmm"]

    def update(self, metabolites, dydt):
        Vmm = metabolites[self.Vmm_idx]
        h_cyt = metabolites[self.h_cyt_idx]
        h_mit = metabolites[self.h_mit_idx]
        Q = metabolites[self.q_idx]
        QH2 = metabolites[self.qh2_idx]
        cytc_ox = metabolites[self.cytc_ox_idx]
        cytc_red = metabolites[self.cytc_red_idx]

        U = Vmm * F / (1000 * R * T)

        tmp = (h_mit / h_cyt)**4
        Keq = np.exp(-2 * self.Em_Q + 2 * self.Em_cytc + 2*U) * tmp

        V  = self.Vmax * (QH2 * cytc_ox**self.n - Q * cytc_red**self.n  / Keq)

        dydt[self.qh2_idx] -= V
        dydt[self.cytc_ox_idx] -= 2 * V
        dydt[self.h_mit_idx] -= 2 * V / self.Volume_cyt_mit

        dydt[self.h_cyt_idx] += 2 * V  / (1 - self.Volume_cyt_mit)
        dydt[self.cytc_red_idx] += 2* V
        dydt[self.q_idx] += V

        dydt[self.Vmm_idx] -= 2 * V * F / self.Cmm # !!!!!!

        return dydt

class Complex4(Enzyme):

    def __init__(self, h_cyt, h_mit, cytc_ox, cytc_red, o2, Vmm, params):

        self.h_cyt_idx = h_cyt
        self.h_mit_idx = h_mit
        self.o2_idx = o2

        self.cytc_ox_idx = cytc_ox
        self.cytc_red_idx = cytc_red
        self.Vmm_idx = Vmm

        self.Km_O2 = params["Km_O2"]
        self.dGh = params["dGh"]
        self.n = params["n"]
        self.Vmax = params["Vmax"]
        self.Km_cytc = params["Km_cytc"]
        self.Volume_cyt_mit = params["Volume_cyt_mit"]
        self.Cmm = params["Cmm"]

    def update(self, metabolites, dydt):

        # Vmm = metabolites[self.Vmm_idx]
        # h_cyt = metabolites[self.h_cyt_idx]
        # h_mit = metabolites[self.h_mit_idx]
        o2 = metabolites[self.o2_idx]


        cytc_red = metabolites[self.cytc_red_idx]

        tmp1 = cytc_red**self.n / (cytc_red**self.n + self.Km_cytc**self.n)

        tmp2 = o2 / (o2 + self.Km_O2)
        tmp3 = (np.exp(-0.001 * self.dGh * F / R / T  ))**2

        V = self.Vmax * tmp1*tmp2*tmp3

        dydt[self.cytc_red_idx] -= 2*V
        dydt[self.o2_idx] -= V
        dydt[self.h_mit_idx] -= 2*V  / self.Volume_cyt_mit

        dydt[self.h_cyt_idx] += 2*V / (1 - self.Volume_cyt_mit)
        dydt[self.cytc_ox_idx] += 2*V

        dydt[self.Vmm_idx] -= 2 * V * F / self.Cmm # !!!!!!

        return dydt
########################################################################################################################

class Pyruvate_exchanger(Enzyme):
    def __init__(self, pyr_cyt, pyr_mit, h_cyt, h_mit, params):
        self.pyr_cyt_idx = pyr_cyt
        self.pyr_mit_idx = pyr_mit
        self.h_cyt_idx = h_cyt
        self.h_mit_idx = h_mit

        self.Km_pyr_mit = params["Km_pyr_mit"]
        self.Km_pyr_cyt = params["Km_pyr_cyt"]
        self.Vmax = params["Vmax"]
        self.Volume_cyt_mit = params["Volume_cyt_mit"]

    def update(self, metabolites, dydt):

        h_cyt = metabolites[self.h_cyt_idx]
        h_mit = metabolites[self.h_mit_idx]
        pyr_cyt = metabolites[self.pyr_cyt_idx]
        pyr_mit = metabolites[self.pyr_mit_idx]

        tmp1 = pyr_cyt * h_cyt - pyr_mit * h_mit
        tmp2 = 1 + pyr_cyt / self.Km_pyr_cyt
        tmp3 = 1 + pyr_mit / self.Km_pyr_mit

        V = self.Vmax * tmp1 / (tmp2 * tmp3)

        dydt[self.pyr_cyt_idx] -= V / (1 - self.Volume_cyt_mit)
        dydt[self.pyr_mit_idx] += V / self.Volume_cyt_mit

        dydt[self.h_mit_idx] -= V / self.Volume_cyt_mit
        dydt[self.h_cyt_idx ] += V / (1 - self.Volume_cyt_mit)

        return dydt
########################################################################################################################

class Pyruvate_dehydrogenase_complex(Enzyme):

    def __init__(self, pyr, CoA, acCoA, fad_pdhc, fadh2_pdhc, nad, nadh, ca, params):
        self.pyr_idx = pyr
        self.CoA_idx = CoA
        self.acCoA_idx = acCoA
        self.fad_pdhc_idx = fad_pdhc
        self.fadh2_pdhc_idx = fadh2_pdhc
        self.nad_idx = nad
        self.nadh_idx = nadh
        self.ca_idx = ca

        self.Amax_Ca = params["Amax_Ca"]
        self.Ka_Ca = params["Ka_Ca"]
        self.Ki_AcoA = params["Ki_AcoA"]
        self.Vmax_pdhc_fad = params["Vmax_pdhc_fad"]
        self.Em_fad = params["Em_fad"]
        self.Em_nad = params["Em_nad"]
        self.Km_nad = params["Km_nad"]
        self.Vmax_pdhc_nad = params["Vmax_pdhc_nad"]
        self.Vmax_pdhc_nad = params["Vmax_pdhc_nad"]
        self.Km_pyr = params["Km_pyr"]
        self.Km_fad = params["Km_fad"]
        self.Km_CoA = params["Km_CoA"]


    def update(self, metabolites, dydt):

        pyr = metabolites[self.pyr_idx]
        CoA = metabolites[self.CoA_idx]
        acCoA = metabolites[self.acCoA_idx]
        fad_pdhc = metabolites[self.fad_pdhc_idx]
        fadh2_pdhc = metabolites[self.fadh2_pdhc_idx]
        nad = metabolites[self.nad_idx]
        nadh = metabolites[self.nadh_idx]
        ca = metabolites[self.ca_idx]

        tmp1 = 1 + self.Amax_Ca * ca / ( ca + self.Ka_Ca )
        tmp2 = pyr / (pyr + self.Km_pyr )
        tmp3 = fad_pdhc / (fad_pdhc + self.Km_fad)

        tmp4 = CoA / (CoA + self.Km_CoA*(1 + acCoA / self.Ki_AcoA))
        pyr_dehyd_compACoA = self.Vmax_pdhc_fad * tmp1 * tmp2 * tmp3 * tmp4

        dydt[self.pyr_idx] -= pyr_dehyd_compACoA
        dydt[self.CoA_idx] -= pyr_dehyd_compACoA
        dydt[self.fad_pdhc_idx] -= pyr_dehyd_compACoA

        dydt[self.acCoA_idx] += pyr_dehyd_compACoA
        dydt[self.fadh2_pdhc_idx] += pyr_dehyd_compACoA

        Keq = np.exp( 2*( self.Em_fad + self.Em_nad)*F*0.001 / R / T  )
        tmp5 = fadh2_pdhc * nad - fad_pdhc * nadh / Keq
        pyr_dehyd_compFad = self.Vmax_pdhc_nad * tmp5 / (nad + self.Km_nad)

        dydt[self.fadh2_pdhc_idx] -= pyr_dehyd_compFad
        dydt[self.nad_idx] -= pyr_dehyd_compFad
        dydt[self.nadh_idx] += pyr_dehyd_compFad
        dydt[self.fad_pdhc_idx] += pyr_dehyd_compFad

        return dydt
########################################################################################################################

class Citrate_synthetase(Enzyme):
    def __init__(self, oa, acCoA, CoA, cit, params):
        self.oa_idx = oa
        self.acCoA_idx = acCoA
        self.CoA_idx = CoA
        self.cit_idx = cit

        self.Km_oxa = params["Km_oxa"]
        self.Ki_cit = params["Ki_cit"]
        self.Ki_CoA = params["Ki_CoA"]
        self.Km_accoa = params["Km_accoa"]
        self.Vmax = params["Vmax"]

    def update(self, metabolites, dydt):
        oa = metabolites[self.oa_idx]
        acCoA = metabolites[self.acCoA_idx]
        CoA = metabolites[self.CoA_idx]
        cit = metabolites[self.cit_idx]

        tmp1 = oa / ( oa + self.Km_oxa * (1 + cit / self.Ki_cit))
        tmp2 = acCoA / (acCoA + self.Km_accoa*(1 + CoA / self.Ki_CoA))

        V = self.Vmax * tmp1 * tmp2

        dydt[self.oa_idx] -= V
        dydt[self.acCoA_idx] -= V
        dydt[self.cit_idx] += V

        return dydt
########################################################################################################################

class Aconitase(Enzyme):

    def __init__(self, citr, isocitr, params):

        self.citr_idx = citr
        self.isocitr_idx = isocitr

        self.Keq = params["Keq"]
        self.Km_cit = params["Km_cit"]
        self.Km_isocit = params["Km_isocit"]
        self.Vmax = params["Vmax"]

    def update(self, metabolites, dydt):
        citr = metabolites[self.citr_idx]
        isocitr = metabolites[self.isocitr_idx]

        tmp1 = citr - isocitr /  self.Keq
        tmp2 = 1 + citr / self.Km_cit + isocitr / self.Km_isocit

        V = self.Vmax * tmp1 / tmp2

        dydt[self.citr_idx] -= V
        dydt[self.isocitr_idx] += V

        return dydt
########################################################################################################################

class Isocitrate_dehydrogenase(Enzyme):

    def __init__(self, isocitr, nad, akg, nadh, ca, params):
        self.isocitr_idx = isocitr
        self.nad_idx = nad
        self.akg_idx = akg
        self.nadh_idx = nadh
        self.ca_idx = ca

        self.Km1_isocit = params["Km1_isocit"]
        self.Km2_isocit = params["Km2_isocit"]
        self.Ka_Ca = params["Ka_Ca"]
        self.n_Ca = params["n_Ca"]
        self.n_isocit = params["n_isocit"]
        self.Ki_nadh = params["Ki_nadh"]
        self.Km_nad = params["Km_nad"]
        self.Vmax = params["Vmax"]

    def update(self, metabolites, dydt):
        isocitr = metabolites[self.isocitr_idx]
        nad = metabolites[self.nad_idx]
        nadh = metabolites[self.nadh_idx]
        ca = metabolites[self.ca_idx]

        Km_isocitr = self.Km1_isocit / (1 + (ca / self.Ka_Ca)**self.n_Ca) + self.Km2_isocit
        tmp1 = isocitr**self.n_isocit / (isocitr**self.n_isocit + Km_isocitr**self.n_isocit)
        tmp2 = nad / (nad + self.Km_nad * (1 + nadh / self.Ki_nadh) )
        V = self.Vmax * tmp1 * tmp2

        dydt[self.isocitr_idx] -= V
        dydt[self.nad_idx] -= V

        dydt[self.nadh_idx] += V
        dydt[self.akg_idx] += V

        return dydt
########################################################################################################################

class Alpha_ketoglutarate_dehydrogenase(Enzyme):

    def __init__(self, ca, akg, nadh, nad, CoA, sucCoA, fad, fadh2, params):
        self.ca_idx = ca
        self.nadh_idx = nadh
        self.nad_idx = nad
        self.CoA_idx = CoA
        self.sucCoA_idx = sucCoA
        self.akg_idx = akg
        self.fad_idx = fad
        self.fadh2_idx = fadh2


        self.Ki_ca = params["Ki_ca"]
        self.Ki_nadh = params["Ki_nadh"]
        self.Km2 = params["Km2"]
        self.Km1 = params["Km1"]
        self.Km_fad = params["Km_fad"]
        self.Km_CoA = params["Km_CoA"]
        self.Vmax_fad = params["Vmax_fad"]
        self.Km_SucCoA = params["Km_SucCoA"]
        self.Em_nad = params["Em_nad"]
        self.Em_fad = params["Em_fad"]
        self.Km_nad = params["Km_nad"]
        self.Vmax_nad = params["Vmax_nad"]

    def update(self, metabolites, dydt):

        ca = metabolites[self.ca_idx]
        nadh = metabolites[self.nadh_idx]
        nad = metabolites[self.nad_idx]
        CoA = metabolites[self.CoA_idx]
        sucCoA = metabolites[self.sucCoA_idx]
        akg = metabolites[self.akg_idx]
        fad = metabolites[self.fad_idx]
        fadh2 = metabolites[self.fadh2_idx]

        Km = ( self.Km1 / (1 + ca / self.Ki_ca) + self.Km2 ) * (1 + nadh / self.Ki_nadh)
        tmp1 = akg / (akg + Km)

        tmp2 = fad / (fad + self.Km_fad)
        tmp3 = CoA / (CoA + self.Km_CoA *(1 + sucCoA / self.Km_SucCoA) )
        Vfad = self.Vmax_fad * tmp1 * tmp2 * tmp3

        dydt[self.akg_idx] -= Vfad
        dydt[self.CoA_idx] -= Vfad
        dydt[self.fad_idx] -= Vfad

        dydt[self.fadh2_idx] += Vfad
        dydt[self.sucCoA_idx] += Vfad

        Keq = np.exp(0.002 * F * ( self.Em_fad + self.Em_nad) / R / T )
        tmp4 = fadh2 * nad - fad * nadh / Keq
        tmp5 = nad + self.Km_nad * (1 + nadh / self.Ki_nadh)
        Vnad = self.Vmax_nad * tmp4 / tmp5

        dydt[self.fadh2_idx] -= Vnad
        dydt[self.nad_idx] -= Vnad

        dydt[self.fad_idx] += Vnad
        dydt[self.nadh_idx] += Vnad
        return dydt
########################################################################################################################

class Succinil_CoA_synthetase(Enzyme):

    def __init__(self, sucCoA, pi, suc, CoA, ndp, ntp, params):

        self.sucCoA_idx = sucCoA
        self.pi_idx = pi
        self.suc_idx = suc
        self.CoA_idx = CoA
        self.ndp_idx = ndp # nucleotide diphosphate ADP or GDP
        self.ntp_idx = ntp # nucleotide triphosphate ATP or GDP

        self.Km_sucCoA = params["Km_sucCoA"]
        self.Km_suc = params["Km_suc"]
        self.Km_CoA = params["Km_CoA"]

        self.Km_nuc3P = params["Km_ndp"]
        self.Km_nuc2P = params["Km_ntp"]

        self.Amax_P = params["Amax_P"]
        self.n_P = params["n_P"]
        self.Km_P = params["Km_P"]
        self.Keq = params["Keq"]
        self.Vmax = params["Vmax"]

    def update(self, metabolites, dydt):

        sucCoA = metabolites[self.sucCoA_idx]
        pi = metabolites[self.pi_idx]
        suc = metabolites[self.suc_idx]
        CoA = metabolites[self.CoA_idx]
        ndp = metabolites[self.ndp_idx]            # nucleotide diphosphate ADP or GDP
        ntp = metabolites[self.ntp_idx]            # nucleotide triphosphate ATP or GDP


        tmp1 = self.Amax_P**self.n_P * pi**self.n_P
        tmp1 /= ( pi**self.n_P + self.Km_P**self.n_P)
        tmp1 += 1

        tmp2 = sucCoA * ndp * pi - suc * CoA * ntp / self.Keq

        vsuccoa = 1 + sucCoA / self.Km_sucCoA
        vnuc2P = 1 + ndp / self.Km_nuc2P
        vp = 1 + pi / self.Km_P

        Vsuc = 1 + suc / self.Km_suc
        Vcoa = 1 + CoA / self.Km_CoA
        vnuc3P = 1 + ntp / self.Km_nuc3P

        V = self.Vmax * tmp1 * tmp2 / (vsuccoa*vnuc2P*vp + Vsuc*Vcoa*vnuc3P - 1)

        dydt[self.sucCoA_idx] -= V
        dydt[self.ndp_idx] -= V
        dydt[self.pi_idx] -= V

        dydt[self.suc_idx] += V
        dydt[self.CoA_idx] += V
        dydt[self.ntp_idx] += V

        return dydt

########################################################################################################################

class Succinate_dehydrydrogenase(Enzyme):
    def __init__(self, suc, fad, fadh2, fum, q, qh2, mal, params):
        self.suc_idx = suc
        self.fad_idx = fad
        self.fadh2_idx = fadh2
        self.fum_idx = fum
        self.q_idx = q
        self.qh2_idx = qh2
        self.mal_idx = mal

        self.Em_FAD = params["Em_FAD"]
        self.Ki_mal = params["Ki_mal"]
        self.Km_suc = params["Km_suc"]
        self.Vmax_succdh = params["Vmax_succdh"]


    def update(self, metabolites, dydt):
        suc = metabolites[self.suc_idx]
        fad = metabolites[self.fad_idx]
        fadh2 = metabolites[self.fadh2_idx]
        fum = metabolites[self.fum_idx]
        Q = metabolites[self.q_idx]
        QH2 = metabolites[self.qh2_idx]
        mal = metabolites[self.mal_idx]

        Keq_succdh = 1.0 # !!!!!!! np.exp( (25 - self.Em_FAD) * F / R / T )
        tmp1 = suc * Q - fum * QH2 / Keq_succdh
        tmp2 = suc + self.Km_suc * (1 + mal / self.Ki_mal )
        v_succdh_fad = self.Vmax_succdh * tmp1 / tmp2
        # dydt[] -= v_succdh_fad


        Keq_pdhc_fad_nad = 1.0 # !!! np.exp(- self.Em_FAD * F / R / T) ## !!!!!! added minus before Em
        # tmp3 = fadh2 * nad - fad * nadh / Keq_pdhc_fad_nad
        # v_succdh = self.Vmax_nadh * tmp3 / self.Km_nad


        return dydt

########################################################################################################################

class Fumarase(Enzyme):
    def __init__(self, mal, fum, params):
        self.mal_idx = mal
        self.fum_idx = fum

        self.Keq = params["Keq"]
        self.Km_fum = params["Km_fum"]
        self.Km_mal = params["Km_mal"]
        self.Vmax = params["Vmax"]

    def update(self, metabolites, dydt):

        mal = metabolites[self.mal_idx]
        fum = metabolites[self.fum_idx]

        tmp1 = fum - mal / self.Keq
        tmp2 = 1 + fum / self.Km_fum - mal / self.Km_mal

        V = self.Vmax * tmp1 / tmp2

        dydt[self.fum_idx] -= V
        dydt[self.mal_idx] += V

        return dydt
########################################################################################################################



