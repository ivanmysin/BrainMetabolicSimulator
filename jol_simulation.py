import numpy as np
import sympy as sym
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import jol_lib as lib
from jol_parameters import params, vars, glob_params


agents = []
agents.append( lib.SodiumLeak(0, 27, params["sodium leak n"])  )
agents.append( lib.SodiumLeak(1, 27, params["sodium leak g"])  )
agents.append( lib.NaKATPase(0, 27, 14, params["Na/K-ATPase n"])  )
agents.append( lib.NaKATPase(1, 27, 15, params["Na/K-ATPase g"])  )
agents.append( lib.GlucoseTransport(25, 2, params["GLC_exchange en"])  )
agents.append( lib.GlucoseTransport(21, 25, params["GLC_exchange ce"])  )
agents.append( lib.GlucoseTransport(25, 3, params["GLC_exchange eg"])  )
agents.append( lib.GlucoseTransport(21, 3, params["GLC_exchange cg"])  )
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

agents.append( lib.LeakCurrent(27, 0, params["leak current"]) )
agents.append( lib.SodiumCurrent(27, 28, 0, params["sodium current"]) )
agents.append( lib.PotassiumCurrent(27, 29, params["potassium current"]) )
agents.append( lib.CalciumCurrent(27, 30, params["calcium current"]) )
agents.append( lib.AHPCurrent(27, 30, params["AHP current"]) )
agents.append( lib.CalciumDecay(30, params["calcium decay"]) )

agents.append( lib.VenousVolume(23, params["Venous flow"]) )
agents.append( lib.DeoxyhemoglobinRate(24, 20, 23, params["Dexyhemoglobin rate"]) )


short_names = [met["short"] for met in vars ]
run_model, jacobian = lib.get_model(short_names, agents, params, glob_params)
y0 = np.asarray( [vars[idx]["rest"] for idx in range(len(vars))] )
sol = solve_ivp(run_model, [0, 100], y0, method="LSODA", jac=jacobian, rtol=1e-3, atol=1e-6 ) #min_step=0.0001

y = np.copy(sol.y)
t = np.copy(sol.t)

indxes = [18, 19, 10, 11, 31, 32]
long_names = [vars[idx]["full"] for idx in range(len(vars))]
for idx in indxes:
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(t, y[idx, :], linewidth=2, label=vars[idx]["full"])
    plt.legend()
plt.show()


# for en in agents:
#
#     run_model, jacobian = lib.get_model(short_names, [en], params, glob_params)
#
#     sol = solve_ivp(run_model, [0, 20], y0, method="LSODA", min_step=0.001, rtol=0.001,   atol=1e-6) # min_step=0.001,  LSODA
#
#
#     print("#################################################")
#     attrs = dir(en)
#     indxes = []
#     enzyme_name = str(en)
#     enzyme_name = enzyme_name.split(" ")[0]
#     enzyme_name = enzyme_name[5:]
#
#     print(enzyme_name)
#     for a in attrs:
#         if "idx" in a:
#             v = en.__getattribute__(a)
#             try:
#                 idx = int(v)
#                 print(vars[idx]["full"])
#                 indxes.append(idx)
#             except TypeError:
#                 pass
#             except ValueError:
#                 pass
#
#     for idx in indxes:
#         fig = plt.figure()
#         ax = fig.add_subplot(1, 1, 1)
#         ax.plot(sol.t,sol.y[idx, :], linewidth = 2, label=vars[idx]["full"])
#         plt.legend()
#
#     for idx, met in enumerate(vars):
#         if idx in indxes:
#             continue
#         res = np.sum( sol.y[idx, :] != y0[idx] )
#         if res != 0:
#             print("Error in metabolite %s with idx %d" % ( vars[idx]["full"], vars[idx]["idx"]) )
#
#     plt.show()





