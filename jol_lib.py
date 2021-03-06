import numpy as np
import sympy as sym
from sympy.functions.special.delta_functions import Heaviside

from scipy import constants as const
from scipy.constants import physical_constants


log = sym.functions.elementary.exponential.log
exp = sym.functions.elementary.exponential.exp
sqrt = sym.sqrt
cos = sym.cos
pi = sym.pi




# F = physical_constants["Faraday constant"][0]
# R = const.R

R   = 8.314510       # J mol-1 K-1
F   = 9.64853e04     # C mol-1
T   = 310  # 37 grad Celcium
RTF = 26.7300



def get_Fin(t):
    return 1.0

def get_adp(atp, qAK, A):

    # thr = 0.00001
    # atp_new = atp # sym.Piecewise((thr, atp < thr), (atp, atp >= thr))
    # if type(atp) is np.ndarray:
    #     sqrt = np.sqrt
    #     # adp = 0.5 * atp * ( -qAK + sqrt(qAK * qAK + 4 * qAK * (A / atp - 1)) )
    # else:
    #     sqrt = sym.sqrt
    adp = 0.5 * atp * ( -qAK + sqrt(qAK * qAK + 4 * qAK * (A / atp - 1)) )

    return adp

def get_damp_datp(atp, qAK, A):

    # if type(atp) is np.ndarray:
    #     sqrt = np.sqrt
    # else:
    #     sqrt = sym.sqrt

    # thr = 0.00001
    #atp_new = atp # sym.Piecewise((thr, atp < thr), (atp, atp >= thr))

    sqrt_ux = sqrt( qAK**2 + 4*qAK * (A / atp - 1) )
    damp_datp = -1.0 + 0.5*qAK - 0.5*sqrt_ux + qAK * A/(atp * sqrt_ux)

    return damp_datp
    
def get_blood_flow():

    # lastF0 = F0 * builtin1(tstim - lag)
    # taurise = 0.1 * 1 / builtin2(0)
    #
    # if CBF == 'on':
    #     if t <= (tstim + t1) & & t >= (t1 + lag):
    #         y = F0 * builtin1(t - t1 - lag);
    #     elif t > (tstim + t1):
    #         y = F0 + (lastF0 - F0) * exp(-(t - tstim - t1) / 5);
    # else:
    #     y = F0 * (1 + 0.1 * builtin3(t - t1, lag, taurise))

    return 0.012


def builtin1():
    pass

def builtin2():
    pass

def builtin3():
    pass


class BaseRate:
    def __init__(self, params):
        pass
    def update(self):
        pass

class SodiumLeak(BaseRate):
    def __init__(self, na, Vpl, params):

        self.na_idx = na
        self.Vpl_idx = Vpl
        self.SmVx = params["SmVx"]
        self.gna_leak = params["gna_leak"]
        self.Cmpl = params["Cmpl"]
        self.na_ext = params["na_ext"]
        self.Vcpl_g = params["Vcpl_g"]


    def update(self, y, dydt, t):

        na = y[self.na_idx]

        if self.Vcpl_g is None:
            Vpl = y[self.Vpl_idx]
        else:
            Vpl = self.Vcpl_g

        Ena = RTF * log(self.na_ext / na)
        Ina = self.gna_leak * (Ena  - Vpl)
        Jl = self.SmVx / F * Ina
        dydt[self.na_idx] += Jl

        if self.Vcpl_g is None:
            dydt[self.Vpl_idx] += Ina / self.Cmpl

        return dydt


class PotassiumLeak(BaseRate):

    def __init__(self, kext, Vpl, params):
        self.kext_idx = kext
        self.Vpl_idx = Vpl
        self.SmVx = params["SmVx"]
        self.gk_leak = params["gk_leak"]
        self.Cmpl = params["Cmpl"]
        self.Kcyt = params["Kcyt"]
        self.Vcpl_g = params["Vcpl_g"]
        self.rex = params["rex"]


    def update(self, y, dydt, t):

        kext = y[self.kext_idx]

        if self.Vcpl_g is None:
            Vpl = y[self.Vpl_idx]
        else:
            Vpl = self.Vcpl_g

        EK = RTF * log(kext / self.Kcyt)
        IK = self.gk_leak * (Vpl - EK)
        Jl = self.SmVx / F * IK
        dydt[self.kext_idx] += Jl / self.rex

        if self.Vcpl_g is None:
            dydt[self.Vpl_idx] -= IK / self.Cmpl

        return dydt


########################################################################################################################
class NaKATPase(BaseRate):
    def __init__(self, nax, Vpl, atp, kext, params):
        self.nax_idx = nax
        self.Vpl_idx = Vpl
        self.atp_idx = atp
        self.kext_idx = kext

        self.kx = params["kx"]
        self.SmVx = params["SmVx"]
        self.Cmpl = params["Cmpl"]
        self.Km = params["Km"]
        self.Na0 = params["Na0"]
        self.rex = params["rex"]

    def update(self, y, dydt, t):
        nax = y[self.nax_idx]
        atp = y[self.atp_idx]
        kext = y[self.kext_idx]

        Jpump = self.SmVx * self.kx * atp * nax / (1 + atp / self.Km) * kext

        dydt[self.nax_idx] -= 3 * Jpump

        dydt[self.kext_idx] -= 2 * Jpump / self.rex

        if not (self.Cmpl is None):
            dydt[self.atp_idx] -= Jpump
            dIPump = F * self.kx * atp * (nax - self.Na0) / (1 + atp / self.Km)
            dydt[self.Vpl_idx] -= dIPump / self.Cmpl
        else:
            dydt[self.atp_idx] -= 7/4 * Jpump

        return dydt

########################################################################################################################
class NaKATPase2(BaseRate):
    def __init__(self, nax, Vpl, atp, kext, params):
        self.nax_idx = nax
        self.Vpl_idx = Vpl
        self.atp_idx = atp
        self.kext_idx = kext

        self.kx = params["kx"]
        self.SmVx = params["SmVx"]
        self.Cmpl = params["Cmpl"]
        self.KmNa = params["KmNa"]
        self.KmK = params["KmK"]
        self.KmATP = params["KmATP"]

        self.rex = params["rex"]

        #self.Na0 = params["Na0"]

    def update(self, y, dydt, t):
        nax = y[self.nax_idx]
        atp = y[self.atp_idx]
        kext = y[self.kext_idx]

        atp = 1.5

        Vatp = 1 / (1 + self.KmATP / atp)
        Vna = (1 / (1 + self.KmNa / nax))**3
        Vk = (1 / (1 + self.KmK/kext))**2
        Jpump = self.SmVx * self.kx * Vk * Vna * Vatp

        dydt[self.nax_idx] -= 3 * Jpump
        dydt[self.kext_idx] -= 2 * Jpump / self.rex
        dydt[self.atp_idx] -= Jpump

        if not (self.Cmpl is None):
            dIPump = F * Jpump / self.SmVx  # is Vext         # F * self.kx * atp * (nax - self.Na0) / (1 + atp / self.Km)
            dydt[self.Vpl_idx] -= dIPump / self.Cmpl
        # else:
        #     dydt[self.atp_idx] -= 7 / 4 * Jpump

        return dydt
########################################################################################################################
class GlucoseTransport(BaseRate):
    def __init__(self, glcx, glcy, params):
        self.glcx_idx = glcx
        self.glcy_idx = glcy

        self.Tmax = params["Tmax"]
        self.Km = params["Km"]
        self.rxy = params["rxy"]

    def update(self, y, dydt, t):
        glcx = y[self.glcx_idx]
        glcy = y[self.glcy_idx]

        Jglc = self.Tmax * (glcx / (glcx + self.Km) - glcy / (glcy + self.Km))

        dydt[self.glcx_idx] -= Jglc / self.rxy
        dydt[self.glcy_idx] += Jglc

        return dydt
########################################################################################################################

class HexokinasePhosphofructoKinase(BaseRate):
    def __init__(self, glcx, atp, gap, params):
        self.glcx_idx = glcx
        self.atp_idx = atp
        self.gap_idx = gap

        self.kx = params["kx"]
        self.Km = params["Km"]
        self.Ki_atp = params["Ki_atp"]
        self.nH = params["nH"]

    def update(self, y, dydt, t):
        glcx = y[self.glcx_idx]
        atp = y[self.atp_idx]

        Jhkpfk = self.kx * atp * glcx / (glcx + self.Km) / ( 1 + (atp/self.Ki_atp)**self.nH )

        # print(2 * Jhkpfk)
        dydt[self.glcx_idx] -= Jhkpfk
        dydt[self.atp_idx] -= 2 * Jhkpfk
        dydt[self.gap_idx] += 2 * Jhkpfk

        return dydt
########################################################################################################################

class PhosphoglycerateKinase(BaseRate):
    def __init__(self, gap, atp, nadh_cyt, pep, params):
        self.gap_idx = gap
        self.pep_idx = pep
        self.atp_idx = atp
        self.nadh_cyt_idx = nadh_cyt

        self.N = params["N"]
        self.kx = params["kx"]

        self.qAK = params["qAK"]
        self.A = params["A"]

    def update(self, y, dydt, t):
        gap = y[self.gap_idx]
        atp = y[self.atp_idx]
        nadh_cyt = y[self.nadh_cyt_idx]

        adp = get_adp(atp, self.qAK, self.A)

        Jpgk = self.kx * gap * adp * (self.N - nadh_cyt) / nadh_cyt
        # print(Jpgk)
        dydt[self.gap_idx] -= Jpgk
        dydt[self.nadh_cyt_idx] += Jpgk
        dydt[self.pep_idx] += Jpgk
        dydt[self.atp_idx] += Jpgk

        return dydt
########################################################################################################################
class PyruvateKinase(BaseRate):
    def __init__(self, pep, atp, pyr, params):
        self.pep_idx = pep
        self.atp_idx = atp
        self.pyr_idx = pyr
        self.kx = params["kx"]
        self.qAK = params["qAK"]
        self.A = params["A"]

    def update(self, y, dydt, t):
        pep = y[self.pep_idx]
        atp = y[self.atp_idx]

        adp = get_adp(atp, self.qAK, self.A)
        Jpk = self.kx * pep * adp

        dydt[self.pep_idx] -= Jpk
        dydt[self.atp_idx] += Jpk
        dydt[self.pyr_idx] += Jpk

        return dydt
########################################################################################################################
class LactateDehydrogenase(BaseRate):

    def __init__(self, lac, pyr, nadh_cyt, params):

        self.lac_idx = lac
        self.pyr_idx = pyr
        self.nadh_cyt_idx = nadh_cyt

        self.kxplus = params["kxplus"]
        self.kxminus = params["kxminus"]
        self.N = params["N"]

    def update(self, y, dydt, t):
        lac = y[self.lac_idx]
        pyr = y[self.pyr_idx]
        nadh_cyt = y[self.nadh_cyt_idx]

        Jldg = self.kxplus*pyr*nadh_cyt - self.kxminus*lac*(self.N - nadh_cyt)

        dydt[self.pyr_idx] -= Jldg
        dydt[self.nadh_cyt_idx] -= Jldg
        dydt[self.lac_idx] += Jldg

        return dydt
########################################################################################################################
class LactateTransport(BaseRate):
    def __init__(self, lacx, lacy, params):
        self.lacx_idx = lacx
        self.lacy_idx = lacy

        self.Tmax = params["Tmax"]
        self.Km = params["Km"]
        self.rxy = params["rxy"]

    def update(self, y, dydt, t):
        lacx = y[self.lacx_idx]
        lacy = y[self.lacy_idx]

        Jlac = self.Tmax * (lacx / (lacx + self.Km) - lacy / (lacy + self.Km))

        dydt[self.lacx_idx] -= Jlac / self.rxy
        dydt[self.lacy_idx] += Jlac

        return dydt
########################################################################################################################
class TCA(BaseRate):
    def __init__(self, pyr, nadh_mit, atp, params):
        self.pyr_idx = pyr
        self.nadh_mit_idx = nadh_mit
        self.atp_idx = atp

        self.kx = params["kx"]
        self.Km_pyr = params["Km_pyr"]
        self.Km_nad = params["Km_nad"]
        self.N = params["N"]

    def update(self, y, dydt, t):

        pyr = y[self.pyr_idx]
        nadh_mit = y[self.nadh_mit_idx]

        Jtca = self.kx * pyr / (pyr + self.Km_pyr) * (self.N - nadh_mit) / (self.N - nadh_mit + self.Km_nad)

        dydt[self.pyr_idx] -= Jtca
        dydt[self.nadh_mit_idx] += 4 * Jtca
        # dydt[self.atp_idx] += Jtca

        return dydt
########################################################################################################################

class ETC(BaseRate):
    def __init__(self, o2, atp, nadh, params):
        self.o2_idx = o2
        self.atp_idx = atp
        self.nadh_idx = nadh

        self.kx = params["kx"]
        self.Km_o2 = params["Km_o2"]
        self.Km_adp = params["Km_adp"]
        self.Km_nadh = params["Km_nadh"]
        self.qAK = params["qAK"]
        self.A = params["A"]

    def update(self, y, dydt, t):
        o2 = y[self.o2_idx]
        atp = y[self.atp_idx]
        nadh = y[self.nadh_idx]
        adp = get_adp(atp, self.qAK, self.A)

        Jetc = self.kx * o2 / (o2 + self.Km_o2) * adp / (adp + self.Km_adp) * nadh / (nadh + self.Km_nadh)

        dydt[self.o2_idx] -= 0.6 * Jetc
        dydt[self.nadh_idx] -= Jetc
        dydt[self.atp_idx] += 3.6 * Jetc

        return dydt

########################################################################################################################

class NADHShuttle(BaseRate):
    def __init__(self, Rcyt, Rmit, params):
        self.Rcyt_idx = Rcyt
        self.Rmit_idx = Rmit

        self.Tmax = params["Tmax"]
        self.Km_cyt = params["Km_cyt"]
        self.Km_mit = params["Km_mit"]
        self.N = params["N"]

    def update(self, y, dydt, t):
        Rcyt = y[self.Rcyt_idx]
        Rmit = y[self.Rmit_idx]

        Rx_minus = Rcyt / (self.N - Rcyt)
        Rx_plus = (self.N - Rmit) / Rmit

        Jshuttle = self.Tmax * Rx_minus / (Rx_minus + self.Km_cyt) *  Rx_plus / (Rx_plus + self.Km_mit)

        dydt[self.Rcyt_idx] -= Jshuttle
        dydt[self.Rmit_idx] += Jshuttle

        return dydt
########################################################################################################################
class CreatineKinase(BaseRate):
    def __init__(self, pcr, atp, params):
        self.pcr_idx = pcr
        self.atp_idx = atp

        self.kxplus = params["kxplus"]
        self.kxminus = params["kxminus"]
        self.Creatinefull = params["Creatinefull"]
        self.qAK = params["qAK"]
        self.A = params["A"]

    def update(self, y, dydt, t):
        pcr = y[self.pcr_idx]
        atp = y[self.atp_idx]
        adp = get_adp(atp, self.qAK, self.A)

        Jck = self.kxplus * adp * pcr  - self.kxminus * atp * (self.Creatinefull - pcr)

        dydt[self.pcr_idx] -= Jck
        dydt[self.atp_idx] += Jck

        return dydt
########################################################################################################################

class OxygenExchange(BaseRate):
     def __init__(self, o2c, o2x, params):
         self.o2c_idx = o2c
         self.o2x_idx = o2x

         self.PScapVx = params["PScapVx"]
         self.Ko2 = params["Ko2"]
         self.HbOP = params["HbOP"]
         self.nh = params["nh"]
         self.rxy = params["rxy"]

     def update(self, y, dydt, t):

         o2c = y[self.o2c_idx]
         o2x = y[self.o2x_idx]

         Jcxo2m = self.PScapVx * ( self.Ko2 / (self.HbOP/o2c - 1)**self.nh - o2x )

         dydt[self.o2c_idx] -= Jcxo2m / self.rxy
         dydt[self.o2x_idx] += Jcxo2m

         return dydt
########################################################################################################################
class CappilaryFlow(BaseRate):
    def __init__(self, varc, params):

        self.varc_idx = varc

        self.Volcap = params["Volcap"]
        self.vara = params["vara"]

    def update(self, y, dydt, t):

        varc = y[self.varc_idx]
        BF = get_blood_flow()
        Jcap = 2 * BF / self.Volcap * (self.vara - varc)

        dydt[self.varc_idx] += Jcap

        return dydt
########################################################################################################################
class LeakCurrent(BaseRate):
    def __init__(self, Vpl, nan, kext, params):
        self.Vpl_idx = Vpl
        self.nan_idx = nan
        self.kext_idx = kext

        self.gl = params["gl"]

        self.Cmpl = params["Cmpl"]
        self.gKpas = params["gKpas"]
        # self.EK = params["EK"]
        self.gNan = params["gNan"]
        self.nae = params["nae"]
        self.Kcyt = params["Kcyt"]

    def update(self, y, dydt, t):
         Vpl = y[self.Vpl_idx]
         nan = y[self.nan_idx]
         kext = y[self.kext_idx]

         EK = -80.0 # RTF * log(kext / self.Kcyt)
         ENa =  RTF * log(self.nae / nan)

         El = self.gKpas * EK / (self.gKpas + self.gNan) + self.gNan * ENa / (self.gKpas + self.gNan)

         Il = self.gl * (El - Vpl)

         dydt[self.Vpl_idx] += Il / self.Cmpl

         return dydt
########################################################################################################################
class SodiumCurrent(BaseRate):
    def __init__(self, Vpl, hNa, nan, params):
        self.Vpl_idx = Vpl
        self.hNa_idx = hNa

        self.nan_idx = nan

        self.gmax = params["gmax"]
        self.Cmpl = params["Cmpl"]
        self.nae = params["nae"]
        self.SmVn = params["SmVn"]

    def update(self, y, dydt, t):

        Vpl = y[self.Vpl_idx]
        hNa = y[self.hNa_idx]
        nan = y[self.nan_idx]

        mNa = self._get_minf(Vpl)
        Ina = self.gmax * mNa**3 * hNa * (RTF * log(self.nae / nan) - Vpl)

        dydt[self.Vpl_idx] += Ina / self.Cmpl
        dydt[self.hNa_idx] = 4 * self._get_dhdt(Vpl, hNa)
        dydt[self.nan_idx] += Ina * self.SmVn / F

        return dydt

    def _get_minf(self, V):

        alpha_m = -0.1*(V + 33) / ( exp( -0.1 * (V + 33) ) - 1 )
        beta_m = 4 * exp( (V + 58) / -12  )
        m_inf = alpha_m / (alpha_m + beta_m )
        return m_inf

    def _get_dhdt(self, V, h):
        alpha_h = 0.07 * exp(-0.1 * (V + 50))
        beta_h = 1 / ( exp(-0.1*(V + 20)) + 1)

        tau_h = 0.001 / (alpha_h + beta_h )
        h_inf = alpha_h / (alpha_h + beta_h )

        dhdt = (h_inf - h) / tau_h
        return dhdt
########################################################################################################################

class PotassiumCurrent(BaseRate):
    def __init__(self, Vpl, nK, kext, params):
        self.Vpl_idx = Vpl
        self.nK_idx = nK
        self.kext_idx = kext

        self.gmax = params["gmax"]
        self.Cmpl = params["Cmpl"]
        self.Kcyt = params["Kcyt"]
        self.SmVn = params["SmVn"]
        self.ren = params["ren"]

    def update(self, y, dydt, t):

        Vpl = y[self.Vpl_idx]
        nK = y[self.nK_idx]
        kext = y[self.kext_idx]

        EK = RTF * log(kext / self.Kcyt ) # -80.0

        IK = self.gmax * nK**4 * (EK - Vpl)

        dydt[self.Vpl_idx] += IK / self.Cmpl
        dydt[self.nK_idx] = 4 * self._get_dndt(Vpl, nK)
        dydt[self.kext_idx] -= IK * self.SmVn / F / self.ren

        return dydt

    def _get_dndt(self, V, n):
        alpha_n = -0.01 * (V + 34) / ( exp(-0.1*(V + 34)) -1 )
        beta_n = 0.125 * exp(-0.04 * (V + 44) )
        tau_n = 0.001 / (alpha_n + beta_n)
        n_inf = alpha_n / (alpha_n + beta_n)

        dndt = (n_inf - n) / tau_n
        return dndt
########################################################################################################################
class CalciumCurrent(BaseRate):
    def __init__(self,  Vpl,  ca, params):
        self.Vpl_idx = Vpl
        self.ca_idx = ca

        self.gmax = params["gmax"]
        self.Cmpl = params["Cmpl"]
        self.ECa = params["ECa"]
        self.SmVn = params["SmVn"]

    def update(self, y, dydt, t):
        Vpl = y[self.Vpl_idx]

        mCa = 1 / (1 + exp( (Vpl + 20) /-9))
        ICa = self.gmax * mCa**2 * (self.ECa - Vpl)

        dydt[self.Vpl_idx] += ICa / self.Cmpl
        dydt[self.ca_idx] += ICa * self.SmVn / F

        return dydt


class AHPCurrent(BaseRate):

    def __init__(self, Vpl, ca, kext, params):
        self.Vpl_idx = Vpl
        self.ca_idx = ca
        self.kext_idx = kext

        self.gmax = params["gmax"]
        self.Kcyt = params["Kcyt"]
        self.KD = params["KD"]
        self.Cmpl = params["Cmpl"]
        self.SmVn = params["SmVn"]
        self.ren = params["ren"]

    def update(self, y, dydt, t):
        Vpl = y[self.Vpl_idx]
        ca = y[self.ca_idx]

        kext = y[self.kext_idx]

        EK = RTF * log(kext / self.Kcyt ) #  -80.0 #

        Iahp = self.gmax * ca / (ca + self.KD) * (EK - Vpl)

        dydt[self.Vpl_idx] += Iahp / self.Cmpl
        dydt[self.kext_idx] -= Iahp * self.SmVn / F / self.ren

        return dydt
########################################################################################################################

class CalciumDecay(BaseRate):

    def __init__(self, ca, params):
        self.ca_idx = ca

        self.tau = params["tau"]
        self.Ca0 = params["Ca0"]


    def update(self, y, dydt, t):

        ca = y[self.ca_idx]
        dydt[self.ca_idx] -= (ca - self.Ca0) / self.tau

        return dydt

########################################################################################################################

class VenousVolume(BaseRate):

    def __init__(self, Vv, params):
        self.Vv_idx = Vv
        self.F0 = params["F0"]
        self.Vv0 = params["Vv0"]
        self.alphav = params["alphav"]
        self.tauv = params["tauv"]

    def update(self, y, dydt, t):
        Vv = y[self.Vv_idx]

        BF = get_blood_flow()

        dydt[self.Vv_idx] = (BF-self.F0*(Vv / self.Vv0)**(1/self.alphav)) / (1 + self.F0 * self.tauv/ self.Vv0*(Vv / self.Vv0)**(-1/2))

        return dydt

class DeoxyhemoglobinRate(BaseRate):
    def __init__(self, dHb, o2c, Vv, params):
        self.dHb_idx = dHb
        self.o2c_idx = o2c
        self.Vv_idx = Vv

        self.F0 = params["F0"]
        self.Vv0 = params["Vv0"]
        self.alphav = params["alphav"]
        self.tauv = params["tauv"]
        self.o2a = params["o2a"]

    def update(self, y, dydt, t):
        dHb = y[self.dHb_idx]
        o2c = y[self.o2c_idx]
        Vv = y[self.Vv_idx]
        dVvdt = dydt[self.Vv_idx]

        BF = get_blood_flow()

        Fout = self.F0 * ( (Vv / self.Vv0)**(1 / self.alphav) + dVvdt * self.tauv / self.Vv0 * (Vv / self.Vv0) ** (-1 / 2))
        #Fout =      F0 * ( (Y(24) / Vv0) ^  (1 / alphav     ) + dY(24) *     tauv /      Vv0 * (Y(24) / Vv0) ^ (-1 / 2));
        dydt[self.dHb_idx] = 2*BF*( self.o2a - o2c) - Fout * dHb / Vv

        return dydt

class ATP_Consumption(BaseRate):
    def __init__(self, atpx, params):
        self.atp_idx = atpx
        self.Vmax = params["Vmax"]
        self.Km = params["Km"]

    def update(self, y, dydt, t):
        atpx = y[self.atp_idx]
        Jatpc = self.Vmax * atpx / (atpx + self.Km)
        dydt[self.atp_idx] -= Jatpc

        return dydt

class Stimulation(BaseRate):
    def __init__(self, Vpl, nan, nag, params):
        self.Vpl_idx = Vpl
        self.nan_idx = nan
        self.nag_idx = nag

        self.ge = params["ge"]
        self.tau = params["tau"]
        self.frequency = params["frequency"]
        self.Cmpl = params["Cmpl"]
        self.Eext = params["Eext"]
        self.SmVn = params["SmVn"]
        self.SmVg = params["SmVg"]
        self.glia_na = params["glia_na"]
        self.gmax_na = params["gmax_na"]

    def update(self, y, dydt, t):
        Vpl = y[self.Vpl_idx]

        # s = cos(2 * pi * t * self.frequency)
        #
        # piece_wise_conditions = tuple([(0, s < 0.9), (t, s >= 0.9)])
        #
        # tstim = sym.Piecewise(*piece_wise_conditions)
        # g = self.ge * exp(-tstim / self.tau)
        # dVdt = g * (self.Eext - Vpl) / self.Cmpl

        tstart = 10
        tend = 30

        tstim = t - tstart
        g = self.ge * exp( -tstim / self.tau )
        Isyne =  g * (self.Eext - Vpl)
        dVdtstim = Isyne / self.Cmpl

        dVdt = sym.Piecewise( (0, t < tstart), (dVdtstim, (t >= tstart) & (t <= tend) ), (0, t > tend)  )
        dNandt =  sym.Piecewise( (0, t < tstart), (self.SmVn / F * 2/3 * Isyne,  (t >= tstart) & (t <= 30)), (0, t > tend)  )
        dNagdt =  sym.Piecewise( (0, t < tstart), (self.SmVg / F * 2/3 * self.glia_na * g, (t >= tstart) & (t <= tend)), (0, t > tend)  )

        dydt[self.Vpl_idx] += dVdt
        dydt[self.nan_idx] += dNandt
        dydt[self.nag_idx] += dNagdt

        return dydt





def get_model(short_names, agents, params, glob_params):

    metabls = sym.symbols(short_names, real=True) # , nonnegative=True
    t = sym.symbols("t")

    # thr = 0.00001
    # metabls[14] = sym.Piecewise((sym.Float(0), metabls[14] < thr), (metabls[14], metabls[14] >= thr))
    # metabls[15] = sym.Piecewise((sym.Float(0), metabls[15] < thr), (metabls[15], metabls[15] >= thr))

    fsx_expr = [0 for _ in range(len(metabls))]

    for en in agents:
        fsx_expr = en.update(metabls, fsx_expr, t)

    # fsx_expr[14] -= params["vATPases n"]["v"]
    # fsx_expr[15] -= params["vATPases g"]["v"]
    # fsx_expr[15] += params["Jpump0"]["vPumpg0"]


    # fsx_expr[14] = sym.Piecewise((0.0, metabls[14] < 0.01), (dvatpn, metabls[14] >= 0.01) )
    # fsx_expr[15] = sym.Piecewise((0.0, metabls[15] < 0.01), (dvatpg, metabls[15] >= 0.01) )

    fsx_expr[14] = fsx_expr[14] / (1 - get_damp_datp(metabls[14], glob_params["qAK"], glob_params["A"]))
    fsx_expr[15] = fsx_expr[15] / (1 - get_damp_datp(metabls[15], glob_params["qAK"], glob_params["A"]))

    fsx_expr[12] = fsx_expr[12] / (1 - glob_params["MVF"])
    fsx_expr[13] = fsx_expr[13] / (1 - glob_params["MVF"])
    fsx_expr[31] = fsx_expr[31] / glob_params["MVF"]
    fsx_expr[32] = fsx_expr[32] / glob_params["MVF"]


    model = sym.lambdify([t, metabls], fsx_expr, modules=['numpy', 'sympy'])
    jacobian_expr = sym.Matrix(fsx_expr).jacobian(metabls)
    jacobian = sym.lambdify([t, metabls], jacobian_expr)

    simulator = Simulator(model, jacobian, short_names)

    return simulator.run_model, simulator.jacobian


class Simulator:
    def __init__(self, model, jacobian, short_names):

        self.model = model
        self.jacobian = jacobian
        self.short_names = short_names

    def _remove_negative(self, y):
        y = np.asarray(y)
        y[ np.isnan(y) ] = 0
        belowzero_idx = y < 0 # 1e-7
        belowzero_idx[27] = False

        # print("################")
        # for i in np.arange(0, y.size, dtype=int)[belowzero_idx]:
        #     print(self.short_names[i], " ", y[i])


        y[belowzero_idx] = 1e-10
        return y

    def run_model(self, t, y):
        y = self._remove_negative(y)
        dydt = self.model(t, y)
        dydt = np.asarray(dydt)
        return dydt

    def jacobian(self, t, y):
        y = self._remove_negative(y)
        jac = self.jacobian(t, y)
        jac = np.asarray(jac)
        return jac
