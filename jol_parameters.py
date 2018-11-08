


# Physical constants and parameters
R       =   8.314510      # J mol-1 K-1
F       =   9.64853e04     # C mol-1
RTF     =   26.73          # mV (corresponds to body temperature 37°C)

FourCIN = 1
# Volumes
Ve      =   0.2
Vcap    =   0.0055
Vg      =   0.25
Vn      =   0.45
MVF     =   0.07           # mitochondrial volume fraction

# Surface to Volume ratios and sodium conductances
SmVn    =   2.5e04         # cm-1
SmVg    =   2.5e04         # cm-1
gNan    =   0.0155         # mS cm-2
gNag    =   0.00623        # mS cm-2
gKpas   =   0.232          # mS cm-2

Nae     =   150            # mmol/L
Vm      =   -70            # mV

qAK = 0.92
totalAdenine = 2.212  # mmol / L

na_ext = 150.0     # Extracellular sodium concentration
Vcpl_g = -70.0     # Membrane potetial of astrocyte
Cmpl = 1e-03
N =  0.212
Volcap = 0.0055
EK = -80.0

# "Ve": Ve,
# "Vcap": Vcap,
# "Vg": Vg,
# "Vn": Vn,
# "MVF": MVF,  # mitochondrial volume fraction

ren = Ve / Vn
reg = Ve / Vg
rce = Vcap / Ve
rcn = Vcap / Vn
rcg = Vcap / Vg

glob_params = {
    "qAK": qAK,
    "A": totalAdenine,  # mmol / L
}



# Varibles y in equstions dydt = f(y, t)
vars = []
vars.append({"idx": 0, "short" : "nan", "rest" : 8.0, "full" : "Neuronal Na"})
vars.append({"idx": 1, "short" : "nag", "rest" : 15.0, "full" : "Glyal Na"})
vars.append({"idx": 2, "short" : "glcn", "rest" : 1.2, "full" : "Neuronal glucose"})
vars.append({"idx": 3, "short" : "glcg", "rest" : 1.9, "full" : "Glyal glucose"})
vars.append({"idx": 4, "short" : "g3pn", "rest" : 0.0046, "full" : "Neuronal glyceraldegyde-3-phosphate"})
vars.append({"idx": 5, "short" : "g3pg", "rest" : 0.0046, "full" : "Glyal glyceraldegyde-3-phosphate"})
vars.append({"idx": 6, "short" : "pepn", "rest" : 0.015, "full" : "Neuronal phosphoenolpuruvate"})
vars.append({"idx": 7, "short" : "pepg", "rest" : 0.015, "full" : "Glyal phosphoenolpuruvate"})
vars.append({"idx": 8, "short" : "pyrn", "rest" : 0.17, "full" : "Neuronal pyruvate"})
vars.append({"idx": 9, "short" : "pyrg", "rest" : 0.17, "full" : "Glyal pyruvate"})
vars.append({"idx": 10, "short" : "lacn", "rest" : 0.6, "full" : "Neuronal lactate"})
vars.append({"idx": 11, "short" : "lacg", "rest" : 0.6, "full" : "Glyal lactate"})
vars.append({"idx": 12, "short" : "nadh_cyt_n", "rest" : 0.006, "full" : "Neuronal cytosolic NADH"})
vars.append({"idx": 13, "short" : "nadh_cyt_g", "rest" : 0.1, "full" : "Glyal cytosolic NADH"})
vars.append({"idx": 14, "short" : "atpn", "rest" : 2.2, "full" : "Neuronal ATP"})
vars.append({"idx": 15, "short" : "atpg", "rest" : 2.2, "full" : "Glyal ATP"})
vars.append({"idx": 16, "short" : "crn", "rest" : 4.9, "full" : "Neuronal phosphocreatine"})
vars.append({"idx": 17, "short" : "crg", "rest" : 4.9, "full" : "Glyal phosphocreatine"})
vars.append({"idx": 18, "short" : "o2n", "rest" : 0.028, "full" : "Neuronal oxygen"})
vars.append({"idx": 19, "short" : "o2g", "rest" : 0.028, "full" : "Glyal oxygen"})
vars.append({"idx": 20, "short" : "o2c", "rest" : 7.0, "full" : "Capillary oxygen"})
vars.append({"idx": 21, "short" : "glcc", "rest" : 4.5, "full" : "Capillary glucose"})
vars.append({"idx": 22, "short" : "lacc", "rest" : 0.55, "full" : "Capillary lactate"})
vars.append({"idx": 23, "short" : "Volven", "rest" : 0.02, "full" : "Venous volume"})
vars.append({"idx": 24, "short" : "dHb", "rest" : 0.058, "full" : "Deoxyhemoglobin"})
vars.append({"idx": 25, "short" : "glce", "rest" : 2.48, "full" : "Extracellular glucose"})
vars.append({"idx": 26, "short" : "lace", "rest" : 0.6, "full" : "Extracellular lactate"})
vars.append({"idx": 27, "short" : "Vpl", "rest" : -73.0, "full" : "Neuronal membrane potential"})
vars.append({"idx": 28, "short" : "hgate", "rest" : 0.99, "full" : "h gate variable"})
vars.append({"idx": 29, "short" : "ngate", "rest" : 0.02, "full" : "n gate variable"})
vars.append({"idx": 30, "short" : "can", "rest" : 0.00005, "full" : "Neuronal calcium"})
vars.append({"idx": 31, "short" : "nadh_mit_n", "rest" : 0.0, "full" : "Neuronal mitochondrial NADH"})
vars.append({"idx": 32, "short" : "nadh_mit_g", "rest" : 0.0, "full" : "Glyal mitochondrial NADH"})


params = {
    "sodium leak n" : {
        "gna_leak": 0.0155,
        "na_ext" : na_ext,
        "SmVx" : SmVn,
        "Vcpl_g" : None,
        "Cmpl"  : Cmpl,
    },

    "sodium leak g": {
        "gna_leak": 0.00623,
        "na_ext": na_ext,
        "Vcpl_g" : Vcpl_g,
        "SmVx" : SmVg,
        "Cmpl"  : Cmpl,
    },

    "Na/K-ATPase n" : {
        "kx"  : 2.49e-06,       # cm (mmol/L)-1 sec-1
        "Km"  : 0.5,            # mmol/L
        "vPumpg0" :  0.0708,         # mmol/L sec-1
        "SmVx" : SmVn,
        "Cmpl"  : Cmpl,
    },
    "Na/K-ATPase g": {
        "kx": 4.64e-07,  # cm (mmol/L)-1 sec-1
        "Km": 0.5,  # mmol/L
        "vPumpg0": 0.0708,  # mmol/L sec-1
        "SmVx": SmVg,
        "Cmpl": Cmpl,
    },

    # GLC exchange constants
    "GLC_exchange en" : {
        "Tmax"   : 0.041,        # mmol/L sec-1
        "Km" :  8,              # mmol/L
    },

    "GLC_exchange ce": {
        "Tmax": 0.239,  # mmol/L sec-1
        "Km": 8,  # mmol/L
    },

    "GLC_exchange eg": {
        "Tmax": 0.147,  # mmol/L sec-1
        "Km": 8,  # mmol/L
    },

    "GLC_exchange cg": {
        "Tmax": 0.00164,  # mmol/L sec-1
        "Km": 8,  # mmol/L
    },

    # Hexokinase-phosphofructokinase system
    "Hexokinase-phosphofructokinase n" : {
        "kx" : 0.050435,       # sec-1
        "Ki_atp"   :  1,              # mmol/L
        "nH"      :   4,              #
        "Km"      :   0.05,           # mmol/L
    },

    "Hexokinase-phosphofructokinase g": {
        "kx": 0.185,  # sec-1
        "Ki_atp": 1,  # mmol/L
        "nH": 4,  #
        "Km": 0.05,  # mmol/L
    },

    # Phosphoglycerate kinase
    "Phosphoglycerate kinase n" : {
        "kx" : 3.97,           # L/mmol sec-1
        "N"  : N,          # mmol/L
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },

    "Phosphoglycerate kinase g": {
        "kx": 401.7,  # L/mmol sec-1
        "N": N,  # mmol/L
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },

    # Pyruvate kinase
    "Pyruvate kinase n" : {
        "kx" : 36.7,         # L/mmol/sec
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },

    "Pyruvate kinase g": {
        "kx": 135.2,  # L/mmol/sec
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },

    # Lactate exchange constants
    "Lactate exchange ne" : {
        "Tmax" : FourCIN*23.5, # mmol/L sec-1
        "Km" : 0.74,       # mmol/L
    },

    "Lactate exchange ge": {
        "Tmax": 107,  # mmol/L sec-1
        "Km": 3.5,  # mmol/L
    },

    "Lactate exchange ec": {
        "Tmax": 0.30,  # mmol/L sec-1
        "Km": 1.0,  # mmol/L # and leegsma-vogt01.pdf
    },

    "Lactate exchange gc": {
        "Tmax": 2.43e-03,  # mmol/L sec-1
        "Km": 1.0,  # mmol/L # see cremer79.pdf
    },


    # Lactate dehydrogenase
    "Lactate dehydrogenase n" : {
        "kxplus" : 78.1,           # L/mmol/sec
        "kxminus" : 0.768,          # L/mmol/sec
        "N": N,  # mmol/L
    },

    "Lactate dehydrogenase g": {
        "kxplus": 1.71,  # L/mmol/sec
        "kxminus": 0.099,  # L/mmol/sec
        "N": N,  # mmol/L
    },

    # NADH Shuttles
    "NADH Shuttles n" : {
        "Tmax" : 9446,           # mmol/L/sec
        "Km_cyt" : 4.653e-08,      # Total sum NAD+ and NADH
        "Km_mit" : 3.666e+05,      #
    },

    "NADH Shuttles g": {
        "Tmax": 134.2,  # mmol/L/sec
        "Km_cyt": 2.614e-04,  #
        "Km_mit": 9.620e+03,  #
    },

    # Constant ATP rates
    "ATP rates n" : {
        "vATPasesn" : 0.1388, #0.0832;     # mmol/L/sec
    },

    "ATP rates g": {
        "vATPasesg": 0.1351,  # 0.1061;     # mmol/L/sec
    },

    # Mitochondrial respiration
    "TCA n" : {
        "kx": 0.147,  # mmol/L/sec
        "N": N,  # mmol/L
        "Km_pyr" : 0.04,
        "Km_nad" : 0.409,
    },

    "TCA g": {
        "kx": 5.31,  # mmol/L/sec
        "N": N,  # mmol/L
        "Km_nad" : 40.3,
        "Km_pyr" : 0.04,
    },

    "Mitochondrial respiration n" : {
        "Km_nadh"  : 0.04,           # mmol/L
        "Km_o2" : 0.001,          # mmol/L
        "Km_adp" : 3.328e-03,       # mmol/L
        "kx" : 0.1610,    # mmol/L/sec
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },

    "Mitochondrial respiration g": {
        "Km_nadh": 0.04,  # mmol/L
        "Km_o2": 0.001,  # mmol/L
        "Km_adp": 4.989e-04,  # mmol/L
        "kx": 0.0627,  # mmol/L/sec
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },

    # Creatine kinase
    "Creatine kinase n" : {
        "kxplus" : 8.85e-02,       # L/mmol/sec
        "kxminus" : 0.00057,        # L/mmol/sec
        "Creatinefull"      : 10.0,               # mmol/L
        "qAK" : qAK,
        "A" : totalAdenine,  # mmol / L
    },

    "Creatine kinase g": {
        "kxplus": 1.0e-03,  # L/mmol/sec
        "kxminus": 0.000007,  # L/mmol/sec
        "Creatinefull": 10.0,  # mmol/L
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },


    # Oxygen exchange constants
    "Oxygen exchange n" : {
        "PScapVx" : 1.6608,    #1.79;       # sec-1
        "Ko2" : 0.0361,      # mmol/L
        "HbOP" : 8.6,        # mmol/L
        "nh" : 2.73,       #
        "Volx" : Vn,
    },

    "Oxygen exchange g": {
        "PScapVx": 0.8736,  # 0.94;       # sec-1
        "Ko2": 0.0361,  # mmol/L
        "HbOP": 8.6,  # mmol/L
        "nh": 2.73,  #
        "Volx" : Vg,
    },

    # Blood flow contributions
    "Blood flow oxygen" : {
        "vara" : 8.35,       # O2a mmol/L
        "Volcap" : Volcap,
    },

    "Blood flow glucose": {
        "vara": 4.75,  # GLCa mmol/L
        "Volcap" : Volcap,

    },
    "Blood flow lactate": {
        "vara": 0.498,  # LACa mmol/L
        "Volcap" : Volcap,
    },



    # Venous flow
    "Venous flow" : {
        "tauv" : 35,         # secs
        "alphav" : 0.5,        #
    },

    # Neuronal membrane parameters
    "leak current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gl" : 0.02,           # mS/cm2
        "gKpas" : 0.2035,
        "EK"    : EK,
        "gNan" : 0.0136,
        "nae"  : na_ext,
    },

    "sodium current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gmax" : 40,   # mS/cm2
        "nae"  : na_ext,

    },

    "potassium current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gmax" : 18,             # mS/cm2
        "EK": EK,  # mV
    },

    "calcium current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gmax" : 0.02,           # mS/cm2
        "ECa" : 120,            # mV
    },

    "AHP current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gmax" : 6.5,            # mS/cm2
        "KD" : 30e-03,         # mmol/L
        "EK" :  EK,            # mV
    },

    "calcium decay" : {
        "tau" : 150e-03,        # sec
        "Ca0"   : 0.5e-04,        # mmol/L
    },



}





"""


# Ratio of excitatory conductance
glia        =   50;         # mV

# Neuronal parameters






Ee      =   0;              # mV
Ei      =   -80;            # mV
phih    =   4;              #
phin    =   4;              #

# # This is auto
man!

kCKnps = 4.33e-002;
kCKgps = 1.35e-003;
KmNADn = 4.09e-001;
KmNADg = 4.03e+001;
KmNADHn = 4.44e-002;
KmNADHg = 2.69e-002;
kLDHnps = 7.23e+001;
kLDHgps = 1.59;
MnCyto = 4.9e-008;
MgCyto = 2.5e-004;
MnMito = 3.93e+005;
MgMito = 1.06e+004;
KmADPn = 3.41e-003;
KmADPg = 4.83e-004;
TMaxLACgc = 2.59e-003;

########################################################################################################################
# # Constrained parameters

TMaxGLCen = 0.041;
TMaxGLCce = 0.239;
TMaxGLCeg = 0.147;
TMaxGLCcg = 0.0016;
kHKPFKn = 0.0504;
kHKPFKg = 0.185;
kPGKn = 3.97;
kPGKg = 401.7;
kPKn = 36.7;
kPKg = 135.2;
kCKnms = 0.00028;
kCKgms = 0.00001;
PScapVn = 1.66;
PScapVg = 0.87;
VMaxMitooutn = 0.164;
VMaxMitooutg = 0.064;
TnNADH = 10330;
TgNADH = 150;
VMaxMitoinn = 0.1303;
VMaxMitoing = 5.7;
vATPasesn = 0.1695;
vATPasesg = 0.1404;
kLDHnms = 0.72;
kLDHgms = 0.071;
TMaxLACne = FourCIN * 24.3;
TMaxLACge = 106.1;
TMaxLACec = 0.25;
LACa = 0.506;

# #

kPumpn = 2.2e-06;
gNan = 0.0136;
gNag = 0.0061;
gKpas = 0.2035;
kPumpg = 4.5e-07;
vPumpg0 = 0.0687;


"""