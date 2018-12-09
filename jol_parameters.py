

FourCIN = 1
# Volumes
Ve   =   0.2
Vcap =   0.0055
Vg   =   0.25
Vn   =   0.45
MVF  =   0.07           # mitochondrial volume fraction

# Surface to Volume ratios and sodium conductances
SmVn = 2.5e04         # cm-1    0.734259   #
SmVg = 2.5e04         # cm-1
gNan = 0.014      # 0.0136   # mS cm-2
gNag = 0.0061     # mS cm-2
gKpasn = 0.21    # 0.232    #          # mS cm-2
gKpasg = 0.02

qAK = 0.92
totalAdenine = 2.212  # mmol / L

Nae = 120.0     # Extracellular sodium concentration,   mmol/L
Kcyt = 80.0    #

Vcpl_g = -70.0     # Membrane potetial of astrocyte
Cmpl = 1e-03
N = 0.212 # 0.212
Volcap = 0.0055
# EK = -112.258 # -80.0

# Parameters for venouse volume
F0 = 0.012 # !!!!
Vv0 = 0.021
alphav = 0.5
tauv = 35
O2a = 8.35

ren = Ve / Vn
reg = Ve / Vg
rce = Vcap / Ve
rcn = Vcap / Vn
rcg = Vcap / Vg

glob_params = {
    "qAK": qAK,
    "A": totalAdenine,  # mmol / L
    "MVF" : MVF,
}

################################################################
params = {
    "sodium leak n" : {
        "gna_leak": gNan,
        "na_ext" : Nae,
        "SmVx" : SmVn,
        "Vcpl_g" : None,
        "Cmpl"  : Cmpl,
    },

    "sodium leak g": {
        "gna_leak": gNag,
        "na_ext": Nae,
        "Vcpl_g" : Vcpl_g,
        "SmVx" : SmVg,
        "Cmpl"  : Cmpl,
    },

    "potassiun leak n" : {
        "gk_leak"  : gKpasn,
        "Kcyt"   : Kcyt,
        "SmVx"   : SmVn,
        "Vcpl_g": None,
        "Cmpl": Cmpl,
        "rex": ren,
    },

    "potassiun leak g": {
        "gk_leak": gKpasg,
        "Kcyt": Kcyt,
        "Vcpl_g" : Vcpl_g,
        "SmVx" : SmVg,
        "Cmpl": Cmpl,
        "rex": ren,
    },

    "Na/K-ATPase n" : {
        "kx"  : 2.2000e-06,       # cm (mmol/L)-1 sec-1
        "Km"  : 0.5,            # mmol/L
        "SmVx" : SmVn,
        "Cmpl"  : Cmpl,
        "Na0": 7.9717,
        "rex": ren,
    },
    "Na/K-ATPase g": {
        "kx": 4.5e-07,  # cm (mmol/L)-1 sec-1
        "Km": 0.5,  # mmol/L
        "SmVx": SmVg,
        "Cmpl": None,
        "Na0": None,
        "rex": reg,
    },

    "Na/K-ATPase2 n": {
        "kx"  : 0.8, # 10000 * 2.2e-06, # 1000 * 2.2000e-06,    # cm (mmol/L)-1 sec-1
        "KmK": 3.5,  # mmol/L
        "KmNa": 10,  # mmol/L
        "KmATP": 170,  # mmol/L
        "SmVx": SmVn,
        "Cmpl": Cmpl,
        "Imax" : 0.013*100, # mA/(cm^2 )
        "rex" : ren,
    },
    "Na/K-ATPase2 g": {
        "kx": 0.8, # 10000 * 4.5e-07,  # cm (mmol/L)-1 sec-1
        "KmK": 3.5,  # mmol/L
        "KmNa": 10,  # mmol/L
        "KmATP": 170,  # mmol/L
        "SmVx": SmVg,
        "Cmpl": None,
        "Imax" : None,
        "rex": reg,
    },

    # GLC exchange constants
    "GLC_exchange en" : {
        "Tmax"   : 0.041,        # mmol/L sec-1
        "Km" : 8,              # mmol/L
        "rxy" : ren,
    },

    "GLC_exchange ce": {
        "Tmax": 0.239,  # mmol/L sec-1
        "Km": 8,  # mmol/L
        "rxy" : rce,
    },

    "GLC_exchange eg": {
        "Tmax": 0.147,  # mmol/L sec-1
        "Km": 8,  # mmol/L
        "rxy" : reg,
    },

    "GLC_exchange cg": {
        "Tmax": 0.0016,  # mmol/L sec-1
        "Km": 8,  # mmol/L
        "rxy" : rcg,
    },

    # Hexokinase-phosphofructokinase system
    "Hexokinase-phosphofructokinase n" : {
        "kx" : 0.0504,       # sec-1
        "Ki_atp"   :  1,              # mmol/L
        "nH"      :   4,              #
        "Km"      :   0.05,           # mmol/L
    },

    "Hexokinase-phosphofructokinase g": {
        "kx": 0.185,  # sec-1
        "Ki_atp": 1,  # 0.5 mmol/L # !!!!!!!!
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
    "Lactate exchange en" : {
        "Tmax" : 0.243, # FourCIN*24.3, # mmol/L sec-1
        "Km" : 0.74,       # mmol/L
        "rxy" : ren,
    },

    "Lactate exchange eg": {
        "Tmax": 1.061, # 106.1,  # mmol/L sec-1
        "Km": 3.5,  # mmol/L
        "rxy": reg,
    },

    # "Lactate exchange ce": {
    #     "Tmax": 0.25,  # mmol/L sec-1
    #     "Km": 1.0,  # mmol/L # and leegsma-vogt01.pdf
    #     "rxy": rce,
    # },
    #
    # "Lactate exchange cg": {
    #     "Tmax": 0.0026,  # mmol/L sec-1
    #     "Km": 1.0,  # mmol/L # see cremer79.pdf
    #     "rxy": rcg,
    # },


    # Lactate dehydrogenase
    "Lactate dehydrogenase n" : {
        "kxplus" : 72.3,           # L/mmol/sec
        "kxminus" : 0.144,          # L/mmol/sec
        "N": N,  # mmol/L
    },

    "Lactate dehydrogenase g": {
        "kxplus": 1.59,  # L/mmol/sec
        "kxminus": 0.071,  # L/mmol/sec
        "N": N,  # mmol/L
    },

    # NADH Shuttles
    "NADH Shuttles n" : {
        "Tmax" : 10330,           # mmol/L/sec
        "Km_cyt" : 4.9e-08,      # Total sum NAD+ and NADH
        "Km_mit" : 393000,      #
        "N": N,  # mmol/L
    },

    "NADH Shuttles g": {
        "Tmax": 150,  # mmol/L/sec
        "Km_cyt": 22.5e-04,  #
        "Km_mit": 10600,  #
        "N": N,  # mmol/L
    },

    # Constant ATP rates
    # "vATPases n" : {
    #     "v" : 0.1695, #0.0832;     # mmol/L/sec
    # },
    #
    # "vATPases g": {
    #     "v": 0.1404,  # 0.1061;     # mmol/L/sec
    # },

    "ATP consumption n" : {
        "Vmax" : 0.1695,
        "Km"   : 1.0,
    },

    "ATP consumption g": {
        "Vmax": 0.1404,
        "Km": 1.0,
    },

    # Mitochondrial respiration
    "TCA n" : {
        "kx": 0.1303,  # mmol/L/sec
        "N": N,  # mmol/L
        "Km_pyr" : 0.04,
        "Km_nad" : 0.409,
    },

    "TCA g": {
        "kx": 5.7,  # mmol/L/sec
        "N": N,  # mmol/L
        "Km_nad" : 40.3,
        "Km_pyr" : 0.04,
    },

    "Mitochondrial respiration n" : {
        "Km_nadh"  : 0.0444,           # mmol/L
        "Km_o2" : 0.001,          # mmol/L
        "Km_adp" : 0.0034,       # mmol/L
        "kx" : 0.164,    # mmol/L/sec
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },

    "Mitochondrial respiration g": {
        "Km_nadh": 0.0269,  # mmol/L
        "Km_o2": 0.001,  # mmol/L
        "Km_adp": 4.8300e-04,  # mmol/L
        "kx": 0.064,  # mmol/L/sec
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },

    # Creatine kinase
    "Creatine kinase n" : {
        "kxplus" : 0.0433,       # L/mmol/sec
        "kxminus" : 2.8e-04,        # L/mmol/sec
        "Creatinefull"      : 10.0,               # mmol/L
        "qAK" : qAK,
        "A" : totalAdenine,  # mmol / L
    },

    "Creatine kinase g": {
        "kxplus": 0.0014,  # L/mmol/sec
        "kxminus": 0.00001,  # L/mmol/sec
        "Creatinefull": 10.0,  # mmol/L
        "qAK": qAK,
        "A": totalAdenine,  # mmol / L
    },

    # Oxygen exchange constants
    "Oxygen exchange n" : {
        "PScapVx" : 1.66,    #1.79;       # sec-1
        "Ko2" : 0.0361,      # mmol/L
        "HbOP" : 8.6,        # mmol/L
        "nh" : 1 / 2.73,       #
        "rxy": rcn,
    },

    "Oxygen exchange g": {
        "PScapVx": 0.87,  # 0.94;       # sec-1
        "Ko2": 0.0361,  # mmol/L
        "HbOP": 8.6,  # mmol/L
        "nh": 1 / 2.73,  #
        "rxy": rcg,
    },

    # Blood flow contributions
    "Blood flow oxygen" : {
        "vara" : O2a,  # O2a mmol/L
        "Volcap" : Volcap,
    },

    "Blood flow glucose": {
        "vara": 4.75,  # GLCa mmol/L
        "Volcap" : Volcap,

    },
    "Blood flow lactate": {
        "vara": 0,# !!! 0.506,  # LACa mmol/L
        "Volcap" : Volcap,
    },

    # Venous flow
    "Venous flow" : {
        "tauv" : tauv,         # secs
        "alphav" : alphav,     #
        "F0" : F0,
        "Vv0" : Vv0,
    },

    "Dexyhemoglobin rate" : {
        "tauv": tauv,  # secs
        "alphav": alphav,  #
        "F0": F0,
        "Vv0": Vv0,
        "o2a" : O2a,

    },

    # Neuronal membrane parameters
    "leak current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gl" : 0.02,           # mS/cm2
        "gKpas" : gKpasn,
        # "EK"    : EK,
        "gNan" : gNan,
        "nae"  : Nae,
        "Kcyt": Kcyt,  # mV
    },

    "sodium current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gmax" : 40,   # mS/cm2
        "nae"  : Nae,
        "SmVn": SmVn,
    },

    "potassium current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gmax" : 18,             # mS/cm2
        "Kcyt": Kcyt,  # mV
        "SmVn": SmVn,
        "ren" : ren,
    },

    "calcium current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gmax" : 0.02,           # mS/cm2
        "ECa" : 120,            # mV
        "SmVn" : SmVn,
    },

    "AHP current" : {
        "Cmpl" : Cmpl,          # mF/cm2
        "gmax" : 6.5,            # mS/cm2
        "KD" : 30e-03,         # mmol/L
        "Kcyt": Kcyt,  # mV
        "SmVn": SmVn,
        "ren" : ren,
    },

    "calcium decay" : {
        "tau" : 150e-03,        # sec
        "Ca0"   : 0.5e-04,        # mmol/L
    },

    # "Jpump0" : {
    #     "vPumpg0" : 0.75 * 0.0687,
    # },


    "stimulation" : {
        "ge" : 0.08,  # 0.0599,
        "tau" : 10.0, # 2.5
        "frequency" : 10.0,
        "Cmpl" : Cmpl,          # mF/cm2
        "Eext" : 0.0,
        "SmVn" : SmVn,
        "SmVg" : SmVg,
        "glia_na" : 50.0,
        "gmax_na" : 40,
    },

}


# Varibles y in equstions dydt = f(y, t)
vars = []
vars.append({"idx": 0, "short" : "nan", "rest" : 8.0, "full" : "Neuronal Na+", "fullrus" : r'Нейрональный $Na^+$', "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 1, "short" : "nag", "rest" : 15.0, "full" : "Glyal Na+", "fullrus" : r"Астроцитарный $Na^+$", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 2, "short" : "glcn", "rest" : 1.2, "full" : "Neuronal glucose", "fullrus" : "Нейрональная глюкоза", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 3, "short" : "glcg", "rest" : 1.9, "full" : "Glyal glucose", "fullrus" : "Астроцитарная глюкоза", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 4, "short" : "g3pn", "rest" : 0.0046, "full" : "Neuronal glyceraldegyde-3-phosphate", "fullrus" : "Нейрональный глицеральдегид 3-фосфат", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 5, "short" : "g3pg", "rest" : 0.0046, "full" : "Glyal glyceraldegyde-3-phosphate", "fullrus" : "Астроцитарный глицеральдегид 3-фосфат", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 6, "short" : "pepn", "rest" : 0.015, "full" : "Neuronal phosphoenolpuruvate", "fullrus" : "Нейрональный фосфоенолпируват", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 7, "short" : "pepg", "rest" : 0.015, "full" : "Glyal phosphoenolpuruvate", "fullrus" : "Астроцитарный фосфоенолпируват", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 8, "short" : "pyrn", "rest" : 0.17, "full" : "Neuronal pyruvate", "fullrus" : "Нейрональный пируват", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 9, "short" : "pyrg", "rest" : 0.17, "full" : "Glyal pyruvate", "fullrus" : "Астроцитарный пируват", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 10, "short" : "lacn", "rest" : 0.6, "full" : "Neuronal lactate", "fullrus" : "Нейрональный лактат", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 11, "short" : "lacg", "rest" : 0.6, "full" : "Glyal lactate", "fullrus" : "Астроцитарный лактат", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 12, "short" : "nadh_cyt_n", "rest" : 0.006, "full" : "Neuronal cytosolic NADH", "fullrus" : "Нейрональный цитозольный НАДН", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 13, "short" : "nadh_cyt_g", "rest" : 0.1, "full" : "Glyal cytosolic NADH", "fullrus" : "Астроцитарный цитозольный НАДН", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 14, "short" : "atpn", "rest" : 2.2, "full" : "Neuronal ATP", "fullrus" : "Нейрональный АТФ", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 15, "short" : "atpg", "rest" : 2.2, "full" : "Glyal ATP", "fullrus" : "Астроцитарный АТФ", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 16, "short" : "crn", "rest" : 4.9, "full" : "Neuronal phosphocreatine", "fullrus" : "Нейроранальный фосфокреатин", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 17, "short" : "crg", "rest" : 4.9, "full" : "Glyal phosphocreatine", "fullrus" : "Астроцитарный фосфокреатин", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 18, "short" : "o2n", "rest" : 0.028, "full" : "Neuronal oxygen", "fullrus" : "Нейрональный кислород", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 19, "short" : "o2g", "rest" : 0.028, "full" : "Glyal oxygen", "fullrus" : "Астроцитарный кислород", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 20, "short" : "o2c", "rest" : 7.0, "full" : "Capillary oxygen", "fullrus" : "Капилярный кислород", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 21, "short" : "glcc", "rest" : 4.5, "full" : "Capillary glucose", "fullrus" : "Капилярная глюкоза", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 22, "short" : "lacc", "rest" : 0.55, "full" : "Capillary lactate", "fullrus" : "Капилярный лактат", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 23, "short" : "Volven", "rest" : 0.02, "full" : "Venous volume", "fullrus" : "Венозный объем", "units":" ", "unitsrus":" "})
vars.append({"idx": 24, "short" : "dHb", "rest" : 0.058, "full" : "Deoxyhemoglobin", "fullrus" : "Дезоксигемоглобин", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 25, "short" : "glce", "rest" : 2.48, "full" : "Extracellular glucose", "fullrus" : "Экстраклеточная глюкоза", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 26, "short" : "lace", "rest" : 0.6, "full" : "Extracellular lactate", "fullrus" : "Экстраклеточный лактат", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 27, "short" : "Vpl", "rest" : -73.0, "full" : "Neuronal membrane potential", "fullrus" : "Потенциал мембраны нейрона", "units":"mV", "unitsrus":"мВ"})
vars.append({"idx": 28, "short" : "hgate", "rest" : 0.99, "full" : "h gate variable", "fullrus" : "h", "units":" ", "unitsrus":" "})
vars.append({"idx": 29, "short" : "ngate", "rest" : 0.02, "full" : "n gate variable", "fullrus" : "n", "units":" ", "unitsrus":" "})
vars.append({"idx": 30, "short" : "can", "rest" : 0.00005, "full" : "Neuronal calcium", "fullrus" : r"Нейрональный $Ca^{2+}$", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 31, "short" : "nadh_mit_n", "rest" : 0.0, "full" : "Neuronal mitochondrial NADH", "fullrus" : "Нейрональный митохондриальный НАДН", "units":"mM", "unitsrus":"мМ"})
vars.append({"idx": 32, "short" : "nadh_mit_g", "rest" : 0.0, "full" : "Glyal mitochondrial NADH", "fullrus" : "Астроцитарный митохондриальный НАДН", "units":"mM", "unitsrus":"мМ"})

vars.append({"idx": 33, "short" : "kex", "rest" : 3.0, "full" : "Extracellular K+", "fullrus" : r"Экстраклеточный $K^+$", "units":"mM", "unitsrus":"мМ"})







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