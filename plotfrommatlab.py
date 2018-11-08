import numpy as np
import matplotlib.pyplot as plt

simulated_data = []

simulated_data.append({"name" : "Neuronal Na"})
simulated_data.append({"name" : "Glyal Na"})
simulated_data.append({"name" : "Neuronal glucose"})
simulated_data.append({"name" : "Glyal glucose"})
simulated_data.append({"name" : "Neuronal glyceraldegyde-3-phosphate"})
simulated_data.append({"name" : "Glyal glyceraldegyde-3-phosphate"})
simulated_data.append({"name" : "Neuronal phosphoenolpuruvate"})
simulated_data.append({"name" : "Glyal phosphoenolpuruvate"})
simulated_data.append({"name" : "Neuronal pyruvate"})
simulated_data.append({"name" : "Glyal pyruvate"})
simulated_data.append({"name" : "Neuronal lactate"})
simulated_data.append({"name" : "Glyal lactate"})
simulated_data.append({"name" : "Neuronal cytosolic NADH"})
simulated_data.append({"name" : "Glyal cytosolic NADH"})
simulated_data.append({"name" : "Neuronal ATP"})
simulated_data.append({"name" : "Glyal ATP"})
simulated_data.append({"name" : "Neuronal phosphocreatine"})
simulated_data.append({"name" : "Glyal phosphocreatine"})
simulated_data.append({"name" : "Neuronal oxygen"})
simulated_data.append({"name" : "Glyal oxygen"})
simulated_data.append({"name" : "Capillary oxygen"})
simulated_data.append({"name" : "Capillary glucose"})
simulated_data.append({"name" : "Capillary lactate"})
simulated_data.append({"name" : "Venous volume"})
simulated_data.append({"name" : "Deoxyhemoglobin"})
simulated_data.append({"name" : "Extracellular glucose"})
simulated_data.append({"name" : "Extracellular lactate"})
simulated_data.append({"name" : "Neuronal membrane potential"})
simulated_data.append({"name" : "h gate variable"})
simulated_data.append({"name" : "n gate variable"})
simulated_data.append({"name" : "Neuronal calcium"})
simulated_data.append({"name" : "Neuronal mitochondrial NADH"})
simulated_data.append({"name" : "Glyal mitochondrial NADH"})


path = "/home/ivan/coding/Jol_code4matlab/Jolivet et al PLoS Comput Biol 2015/Data files/control/"
data_control = np.loadtxt(path + "raw_0020.data")
path = "/home/ivan/coding/Jol_code4matlab/Jolivet et al PLoS Comput Biol 2015/Data files/expr2/"

data_exper1 = np.loadtxt(path + "raw_0020.data")

# data = np.append(data[0, :].reshape(1, -1), data, axis=0 )   #  cat(1,DATA(1,:),DATA);
# data = np.append(data[0, :].reshape(1, -1), data, axis=0 )   #  cat(1,DATA(1,:),DATA);
# data[0, 0] = -100
# data[1, 0] = 0
indexes4polot = [27, 14, 10, 16]
tc = data_control[:, 0]
texp = data_exper1[:, 0]

for idx, var in enumerate(simulated_data):
    if not (idx in indexes4polot):
        continue

    fig, ax = plt.subplots(ncols=1, nrows=2, sharex=True, sharey=True)
    ax[0].set_title( var["name"] )
    ax[0].plot(tc, data_control[:, idx+8], linewidth=2, label="Контроль", color="green"  ) # var["name"]
    ax[0].legend()
    ax[1].plot(texp, data_exper1[:, idx+8], linewidth=2, label="Эксперимент", color="blue") # var["name"]
    ax[1].set_xlim(0, 25)
    ax[1].legend()

plt.show()



















"""
%% Balance equations

dY(1)   =   vLeakNan-3*vPumpn+vnstim; # Neuronal Na"
dY(2)   =   vLeakNag-3*vPumpg+vgstim; # Glyal Na
dY(3)   =   vGLCen-vHKPFKn;           # Neuronal glucose
dY(4)   =   vGLCcg+vGLCeg-vHKPFKg;    # Glyal glucose
dY(5)   =   2*vHKPFKn-vPGKn;          # Neuronal glyceraldegyde-3-phosphate
dY(6)   =   2*vHKPFKg-vPGKg;          # Glyal glyceraldegyde-3-phosphate
dY(7)   =   vPGKn-vPKn;               # Neuronal phosphoenolpuruvate
dY(8)   =   vPGKg-vPKg;               # Glyal phosphoenolpuruvate
dY(9)   =   vPKn-vLDHn-vMitoinn;      # Neuronal pyruvate
dY(10)  =   vPKg-vLDHg-vMitoing;      # Glyal pyruvate
dY(11)  =   vLDHn-vLACne;             # Neuronal lactate
dY(12)  =   vLDHg-vLACge-vLACgc;      # Glyal lactate
dY(13)  =   1/(1-MVF)*(vPGKn-vLDHn-vShuttlen);  # Neuronal cytosolic NADH
dY(14)  =   1/(1-MVF)*(vPGKg-vLDHg-vShuttleg);  # Glyal cytosolic NADH
AMPATP
dY(15)  =   (-2*vHKPFKn+vPGKn+vPKn-vATPasesn-vPumpn+3.6*vMitooutn+vCKn)/(1-dAMPdATPn);   # Neuronal ATP
dY(16)  =   (-2*vHKPFKg+vPGKg+vPKg-vATPasesg-7/4*vPumpg+3/4*vPumpg0+3.6*vMitooutg+vCKg)/(1-dAMPdATPg); # Glyal ATP
dY(17)  =   -vCKn;  # Neuronal phosphocreatine
dY(18)  =   -vCKg;  # Glyal phosphocreatine
dY(19)  =   vO2mcn-0.6*vMitooutn; # Neuronal oxygen
dY(20)  =   vO2mcg-0.6*vMitooutg; # Glyal oxygen

%% Capillaries

dY(21)  =   vO2c-1/rcn*vO2mcn-1/rcg*vO2mcg;   # Capillary oxygen
dY(22)  =   vGLCc-1/rce*vGLCce-1/rcg*vGLCcg;  # Capillary glucose   
dY(23)  =   vLACc+1/rce*vLACec+1/rcg*vLACgc;  # Capillary lactate   
dY(24)  =   (BF-F0*(Y(24)/Vv0)^(1/alphav))/(1+F0*tauv/Vv0*(Y(24)/Vv0)^(-1/2));  # Venous volume
Fout    =   F0*((Y(24)/Vv0)^(1/alphav)+dY(24)*tauv/Vv0*(Y(24)/Vv0)^(-1/2));
dY(25)  =   BF*(O2a-O2cbar)-Fout*Y(25)/Y(24);  # Deoxyhemoglobin

%% Extracellular compartment

dY(26)  =   vGLCce-1/reg*vGLCeg-1/ren*vGLCen;  # Extracellular glucose
dY(27)  =   1/ren*vLACne+1/reg*vLACge-vLACec;  # Extracellular lactate

%% Hodgkin-Huxley equations

dY(28)  =   1/Cm*(-IL-INa-IK-ICa-ImAHP-dIPump+Isyne+Isyni);  # Neuronal membrane potential
dY(29)  =   phih*(hinf-Y(29))/tauh;                          # h gate variable
dY(30)  =   phin*(ninf-Y(30))/taun;                          # n gate variable
dY(31)  =   -SmVn/F*ICa-(Y(31)-Ca0)/tauCa;                   # Neuronal calcium


%% Mitochondrial NADH

dY(32)  =   1/MVF*(4*vMitoinn-vMitooutn+vShuttlen);  # Neuronal mitochondrial NADH
dY(33)  =   1/MVF*(4*vMitoing-vMitooutg+vShuttleg);  # Glyal mitochondrial NADH
"""
