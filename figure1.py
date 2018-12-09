import numpy as np
import sympy as sym
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import jol_lib as lib
from jol_parameters import params, vars, glob_params

plt.rc('axes', linewidth=2)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
plt.rc('font', size=18)
plt.rc('lines', linewidth=1)
plt.rc('lines', markersize=4)
plt.rc('lines', color="black")
# RTF = lib.RTF
# params["Hexokinase-phosphofructokinase n"]["kx"] *= 0.7
# params["Hexokinase-phosphofructokinase g"]["kx"] *= 0.7

agents = []
agents.append( lib.SodiumLeak(0, 27, params["sodium leak n"]) )
agents.append( lib.SodiumLeak(1, 27, params["sodium leak g"]) )

agents.append( lib.PotassiumLeak(33, 27, params["potassiun leak n"]) )
agents.append( lib.PotassiumLeak(33, 27, params["potassiun leak g"]) )

agents.append( lib.NaKATPase2(0, 27, 14, 33, params["Na/K-ATPase2 n"]) )
agents.append( lib.NaKATPase2(1, 27, 15, 33, params["Na/K-ATPase2 g"]) )

## agents.append( lib.NaKATPase(0, 27, 14, 33, params["Na/K-ATPase n"]) )
## agents.append( lib.NaKATPase(1, 27, 15, 33, params["Na/K-ATPase g"]) )

agents.append( lib.ATP_Consumption(14, params["ATP consumption n"]) )
agents.append( lib.ATP_Consumption(15, params["ATP consumption g"]) )

agents.append( lib.GlucoseTransport(25, 2, params["GLC_exchange en"]) )
agents.append( lib.GlucoseTransport(21, 25, params["GLC_exchange ce"]) )
agents.append( lib.GlucoseTransport(25, 3, params["GLC_exchange eg"]) )
agents.append( lib.GlucoseTransport(21, 3, params["GLC_exchange cg"]) )
agents.append( lib.HexokinasePhosphofructoKinase(2, 14, 4, params["Hexokinase-phosphofructokinase n"])  )
agents.append( lib.HexokinasePhosphofructoKinase(3, 15, 5, params["Hexokinase-phosphofructokinase g"])  )
agents.append( lib.PhosphoglycerateKinase(4, 14, 12, 6, params["Phosphoglycerate kinase n"]) )
agents.append( lib.PhosphoglycerateKinase(5, 15, 13, 7, params["Phosphoglycerate kinase g"]) )
agents.append( lib.PyruvateKinase(6, 14, 8, params["Pyruvate kinase n"]) )
agents.append( lib.PyruvateKinase(7, 15, 9, params["Pyruvate kinase g"]) )
agents.append( lib.LactateDehydrogenase(10, 8, 12, params["Lactate dehydrogenase n"]) )
agents.append( lib.LactateDehydrogenase(11, 9, 13, params["Lactate dehydrogenase g"]) )
agents.append( lib.LactateTransport(26, 10, params["Lactate exchange en"]) )
agents.append( lib.LactateTransport(26, 11, params["Lactate exchange eg"]) )
agents.append( lib.LactateTransport(22, 26, params["Lactate exchange ce"]) )
agents.append( lib.LactateTransport(22, 11, params["Lactate exchange cg"]) )
agents.append( lib.TCA(8, 31, 14, params["TCA n"]) )
agents.append( lib.TCA(9, 32, 15, params["TCA g"]) )
agents.append( lib.ETC(18, 14, 31, params["Mitochondrial respiration n"]) )
agents.append( lib.ETC(19, 15, 32, params["Mitochondrial respiration g"]) )
agents.append( lib.NADHShuttle(12, 31, params["NADH Shuttles n"]) )
agents.append( lib.NADHShuttle(13, 32, params["NADH Shuttles g"]) )

agents.append( lib.CreatineKinase(16, 14, params["Creatine kinase n"]) )
agents.append( lib.CreatineKinase(17, 15, params["Creatine kinase g"]) )


agents.append( lib.OxygenExchange(20, 18, params["Oxygen exchange n"]) )
agents.append( lib.OxygenExchange(20, 19, params["Oxygen exchange g"]) )
agents.append( lib.CappilaryFlow(20, params["Blood flow oxygen"]) )
agents.append( lib.CappilaryFlow(21, params["Blood flow glucose"]) )
agents.append( lib.CappilaryFlow(22, params["Blood flow lactate"]) )

# # agents.append( lib.LeakCurrent(27, 0, 33, params["leak current"]) )

agents.append( lib.SodiumCurrent(27, 28, 0, params["sodium current"]) )
agents.append( lib.PotassiumCurrent(27, 29, 33, params["potassium current"]) )
agents.append( lib.CalciumCurrent(27, 30, params["calcium current"]) )
agents.append( lib.AHPCurrent(27, 30, 33, params["AHP current"]) )
agents.append( lib.CalciumDecay(30, params["calcium decay"]) )
agents.append( lib.VenousVolume(23, params["Venous flow"]) )
agents.append( lib.DeoxyhemoglobinRate(24, 20, 23, params["Dexyhemoglobin rate"]) )

agents.append( lib.Stimulation(27, 0, 1, params["stimulation"]) )


# y0 = np.loadtxt("initstate.data", dtype=np.float64)  #  np.asarray( [vars[idx]["rest"] for idx in range(len(vars))] )  #
# y0 = np.append(y0, [3.0])
y0 = np.load("init_vars.npy")  #

short_names = [met["short"] for met in vars ]
run_model, jacobian = lib.get_model(short_names, agents, params, glob_params)



y = np.empty( (len(vars), 0), dtype=float)
t = np.array([], dtype=float)

one_simulation_time = 50.0
for i in range(1):
    sol = solve_ivp(run_model, [0, one_simulation_time], y0, method="LSODA", jac=jacobian, rtol=1e-3, atol=1e-6)  # min_step=0.0001
    y0 = sol.y[:, -1]

    y = np.append(y, sol.y, axis=1)
    t = np.append(t, sol.t + one_simulation_time * i)



adpn = lib.get_adp(y[14, :], 0.92, 2.212)
adpg = lib.get_adp(y[15, :], 0.92, 2.212)

ext_g_input = np.zeros_like(t)
ext_g_input[(t>=10)&(t<=30)] = params["stimulation"]["ge"] * np.exp(-(t[(t>=10)&(t<=30)]-10) / params["stimulation"]["tau"])



indxes = [27, 0, 1, 30, 33] # , 33, 14, 15, 12, 13, 31, 32, 18, 20] # list( range(len(vars)) ) #
long_names = [vars[idx]["fullrus"] for idx in range(len(vars))]
units = [vars[idx]["unitsrus"] for idx in range(len(vars))]

fig1, axes = plt.subplots(nrows=len(indxes)+1, ncols=1, figsize=(10, 15), constrained_layout=True )

axes[0].set_title("Возбуждающий сигнал")
axes[0].plot( t, ext_g_input, linewidth=2, color="black")
axes[0].set_ylim(0, 1.2*params["stimulation"]["ge"])
axes[0].broken_barh([(10, 20)], (- 2, 1000), color="gray", alpha=0.2)
axes[0].set_xlim(0, t[-1])
axes[0].set_ylabel(r"$мкСм/см^2$")

for i, idx in enumerate(indxes):
    i += 1
    axes[i].plot( t, y[idx, :], linewidth=2, color="black") # , label=long_names[idx]
    axes[i].set_title(long_names[idx])
    # axes[i].legend()
    axes[i].set_xlim(0, t[-1])
    ymin = np.min(y[idx, :])
    ymin -= 0.2*np.abs(ymin)
    ymax = 1.2 * np.max(y[idx, :])
    axes[i].set_ylim(ymin, ymax)
    axes[i].broken_barh( [(10, 20)], (ymin-2, 1000), color="gray", alpha=0.2 )
    axes[i].set_ylabel(units[idx])

axes[i].set_xlabel("Время, сек")

fig1.savefig("./figures/fig2.png", dpi=300)



indxes = [2, 3 , 12, 13, 31, 32, 14, 15, -1, -2]
fig2, axes = plt.subplots(nrows=len(indxes)//2, ncols=2, constrained_layout=True, figsize=(15, 15) )


for i, idx in enumerate(indxes):

    if i%2 == 0:
        j = 0
    else:
        j = 1

    i = i // 2


    if idx > 0:
        ymet = y[idx, :]
    elif idx == -1:
        ymet = adpn
        long_names[idx] = "Нейрональный АДФ"
        units[idx] = "мМ"
    elif idx == -2:
        ymet == adpg
        long_names[idx] = "Астроцитарный АДФ"
        units[idx] = "мМ"

    axes[i, j].plot( t, ymet, linewidth=2, color="black") # , label=long_names[idx]
    axes[i, j].set_title(long_names[idx])
    axes[i, j].set_xlim(0, t[-1])

    ymin = 0.9*np.min(ymet)
    # ymin -= 0.2*np.abs(ymin)
    ymax = 1.1 * np.max(ymet)

    axes[i, j].set_ylim(ymin, ymax)
    axes[i, j].broken_barh( [(10, 20)], (ymin-2, 1000), color="gray", alpha=0.2 )
    axes[i, j].set_ylabel(units[idx])

axes[i, 0].set_xlabel("Время, сек")
axes[i, 1].set_xlabel("Время, сек")




fig2.savefig("./figures/fig3.png", dpi=300)

plt.show()








