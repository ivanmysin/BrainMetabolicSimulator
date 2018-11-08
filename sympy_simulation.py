import numpy as np
import sympy as sym
import sympy_lib as lib
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from parameters import *

np.set_printoptions(precision=5, suppress=True)

enzymes = []

enzymes.append( lib.Oxigen_diffusion(51, enzyme_params["oxigen_diffussion"]) )
# enzymes.append( lib.Pyruvate_diffusion(71, enzyme_params["pyr_diffussion"]) )
enzymes.append( lib.Lactate_diffusion(24, enzyme_params["lac_diffussion"]) )
enzymes.append( lib.Glucose_diffusion(0, enzyme_params["glc_diffussion"]) )

enzymes.append( lib.GlucoseTransporter(0, 1, enzyme_params["glc_trs"]) )
enzymes.append( lib.Hexokinase(1, 2, 12, 3, enzyme_params["hexokinase"]) )
enzymes.append( lib.Glucose6phosphate_isomerase(12, 13, enzyme_params["glc6p_isomerase"]) )
enzymes.append( lib.Phosphofructokinase_type1(13, 2, 15, 3, 8, 14, enzyme_params["phosphofructokinase1"] ) )
enzymes.append( lib.Fructose16_bisphosphatase(15,13, 8, enzyme_params["fru-1,6-bisphosphatase"] ) )
enzymes.append( lib.Phosphofructokinase_type2(13, 2, 14, 3, 4, enzyme_params["phosphofructokinase2"] ) )
enzymes.append( lib.Fructose26_bisphosphatase(14, 13, 8, enzyme_params["fru-2,6-bisphosphatase"] ) )
enzymes.append( lib.Aldolase(15, 16, 17, enzyme_params["aldolase"] ) )
enzymes.append( lib.Triosophosphate_isomerase(16, 17, enzyme_params["triosep-isomerase"] ) )
enzymes.append( lib.Glyceraldehyde_3_phosphate_dehydrogenase(16, 8, 39, 18, 40, enzyme_params["grap_dehydr"] ) )
enzymes.append( lib.Phosphoglycerate_kinase(18, 3, 19, 2, enzyme_params["p-glyceratekinase"] ) )
enzymes.append( lib.Phosphoglycerate_mutase(19, 70, enzyme_params["p-gricerate_mutase"] ) )
enzymes.append( lib.Enolase(70, 20, enzyme_params["enolase"] ) )
enzymes.append( lib.Pyruvate_kinase(20, 3, 21, 2, enzyme_params["pyruvatekinase"] ) )
enzymes.append( lib.Lactate_dehydrogenase(21, 40, 23, 39, enzyme_params["LDG"] ) )


enzymes.append( lib.Monocarboxilate_transporter(24, 23, enzyme_params["MCT"] ) )
# # enzymes.append( lib.Creatine_kinase(2, 25, 3, 26,  enzyme_params["creatinekinase"] ) ) # Делает АТФ отрицательным !!!!

enzymes.append( lib.Malate_dehydrogenase(27, 28, 39, 40,  enzyme_params["malatdehyd"] ) ) # Cytosolic enzyme
enzymes.append( lib.Malate_dehydrogenase(29, 30, 41, 42,  enzyme_params["malatdehyd"] ) ) # Mitochondrial enzyme
enzymes.append( lib.Aspartate_aminotransferase(31, 33, 28, 35, enzyme_params["asp_aminotrans"]))  # Cytosolic enzyme
enzymes.append( lib.Aspartate_aminotransferase(32, 34, 30, 36,  enzyme_params["asp_aminotrans"] ) ) # Mitochondrial enzyme
enzymes.append( lib.Aspartate_glutamate_carrier(32, 35, 56, 31, 36, 57, 65, enzyme_params["asp_glu_carrier"] ) ) # !!!!!!!!
enzymes.append( lib.Malate_alphaketoglutarate_carrier(27, 34, 29, 33, enzyme_params["mal_akg_carrier"] ) )
enzymes.append( lib.Glycerol_3phosphate_dehydrogenase_cytosolic(17, 40, 19, 39, enzyme_params["cytgly3pdehyd"] ) )
enzymes.append( lib.Glycerol_3phosphate_dehydrogenase_mitochondrial(19, 17, 45, 46, 47, 48, enzyme_params["mitgly3pdehyd"] ) )
enzymes.append( lib.ATP_synthetase(5, 6, 9, 56, 57, 65, enzyme_params["atp_syntase"] ) )
enzymes.append( lib.ATP_ADP_axchanger(5, 3, 6, 2, 65, enzyme_params["atp/adp_axchanger"] ) )
enzymes.append( lib.ATP_consumption(2, 3, 8, enzyme_params["atp_consumption"] ) )

enzymes.append( lib.Passive_efflux_ion(53, 52, 65, enzyme_params["potassium_ed"]))  # Efflux for potassium (K)
enzymes.append( lib.Passive_efflux_ion(55, 54, 65, enzyme_params["sodium_ed"]))  # Efflux for sodium (Na)
enzymes.append( lib.Passive_efflux_ion(56, 57, 65, enzyme_params["protons_ed"]))  # Efflux for protons (H)
enzymes.append( lib.Pump(53, 52, 56, 57, enzyme_params["K_pump"]))  # Pump for potassium (K)
enzymes.append( lib.Pump(55, 54, 56, 57, enzyme_params["Na_pump"]))  # Pump for sodium (Na)
enzymes.append( lib.Pump(8, 9, 56, 57, enzyme_params["phos_pump"]))  # Pump for inorganic phosphate
enzymes.append( lib.Calcium_effux(58, 59, 65, enzyme_params["calcium_ed"]) )
#
#         # # ca_mit, na_cyt, ca_cyt, na_mit, Vmm
enzymes.append( lib.Ca_Na_pump(59, 55, 58, 54, 65, enzyme_params["ca_na_pump"]) )
# # #         #
# # #         # #  ca_mit, h_cyt, ca_cyt, h_mit, Vmm
enzymes.append( lib.Ca_H_pump(59, 56, 58, 57, 65, enzyme_params["ca_h_pump"]) )
# #         #
# #         # # h_cyt, h_mit, q, qh2, nad, nadh, Vmm
enzymes.append( lib.Complex1(56, 57, 47, 48, 41, 42, 65, enzyme_params["complex1"]) )
# # #         #
# # #         # # h_cyt, h_mit, q, qh2, cytc_ox, cytc_red, Vmm,
enzymes.append( lib.Complex3(56, 57, 47, 48, 49, 50, 65, enzyme_params["complex3"]) )
# # #         #
# # #         # # h_cyt, h_mit, cytc_ox, cytc_red, o2, Vmm
enzymes.append( lib.Complex4(56, 57, 49, 50, 51, 65, enzyme_params["complex4"]) )
# # #         #
# #         # # pyr_cyt, pyr_mit, h_cyt, h_mit
enzymes.append( lib.Pyruvate_exchanger(21, 22, 56, 57, enzyme_params["pyr_exchanger"]) )
# #
# #         # # pyr, CoA, acCoA, fad_pdhc, fadh2_pdhc, nad, nadh, ca
enzymes.append( lib.Pyruvate_dehydrogenase_complex(22, 60, 61, 66, 67, 41, 42, 59, enzyme_params["pyr_dehyd_comp"]) )
# #         #
# #         # # oa, acCoA, CoA, cit
enzymes.append( lib.Citrate_synthetase(30, 61, 60, 63, enzyme_params["citrate_syntase"]) )
# #         #
# #         # # citr, isocitr
enzymes.append( lib.Aconitase(63, 64, enzyme_params["aconitase"]) )
# #         #
# #         # # isocitr, nad, akg, nadh, ca,
enzymes.append( lib.Isocitrate_dehydrogenase(64, 41, 34, 42, 59, enzyme_params["isocit_dehydr"]) )
# #         #
# #         # # ca, akg, nadh, nad, CoA, sucCoA, fad, fadh2
enzymes.append( lib.Alpha_ketoglutarate_dehydrogenase(59, 34, 42, 41, 60, 62, 68, 69, enzyme_params["akg_dehydr"]) )
# #         #
# #         # #  sucCoA, pi, suc, CoA, adp, atp
enzymes.append( lib.Succinil_CoA_synthetase(62, 9, 37, 60, 6, 5, enzyme_params["sucCoAsyntase_4atp"]) )
# #         #
# #         # #  sucCoA, pi, suc, CoA, gdp, gtp
enzymes.append( lib.Succinil_CoA_synthetase(62, 9, 37, 60, 11, 10, enzyme_params["sucCoAsyntase_4gtp"]) )
# #         #
# #         # #  suc, fad, fadh2, fum, q, qh2, mal,
enzymes.append( lib.Succinate_dehydrydrogenase(37, 43, 44, 38, 47, 48, 29, enzyme_params["suc_dehydr"]) )
# #         #
# #         # # mal, fum,
enzymes.append( lib.Fumarase(29, 38, enzyme_params["fumarase"]) )

short_names = [met["short"] for met in metabolites ]
# print(len(names), len(set(names)))
metabls = sym.symbols(short_names, real=True) # , nonnegative=True

fsx_expr = [ 0 for _ in range(len(short_names)) ]

for en in enzymes:
    fsx_expr = en.update(metabls, fsx_expr)

y0 = np.asarray( [metabolites[idx]["rest"] for idx in range(len(metabolites))] )

long_names = [metabolites[idx]["full"] for idx in range(len(metabolites))]

t = sym.symbols("t")
model = sym.lambdify([t, metabls], fsx_expr, "numpy")
jacobian_expr = sym.Matrix(fsx_expr).jacobian(metabls)
jacobian = sym.lambdify([t, metabls], jacobian_expr, "numpy")

def run_model(t, y):
    y = np.asarray(y)
    dydt = np.asarray( model(t, y) )
    return dydt

def get_jacobian(t, y):
    J = np.asarray(jacobian(t, y))
    return J



sol = solve_ivp(run_model, [0, 100], y0, method="LSODA", jac=get_jacobian) #  , min_step=0.001, rtol=0.0001,   atol=1e-6

y = np.copy(sol.y)
t = np.copy(sol.t)

# print( type(y) )
# diff = (y[:, -1] - y0) / y0
#
# diff_idx = np.argsort(diff)
#
# for idx in diff_idx:
#     print(diff[idx], metabolites[idx]["full"])

# y0 = y[:, -1]

# for idx in range(y0.size):
#     print(metabolites[idx]["full"] + " %f" % y0[idx])

# print(y0[70])
#
# y0[70] = 0 #  0.0001
# y0[56] = 10**-3
# y0[57] = 10**-5
# #y0[55] = 0.07
# y0[58] = 0.0001
#
# sol = solve_ivp(run_model, [0, 100], y0, method="LSODA", jac=get_jacobian)
#
# y = np.append(y, sol.y, axis=1)
# t = np.append(t, sol.t + t[-1])


indxes = [55, 58, 40, 42]

for idx in indxes:
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(t, y[idx, :], linewidth=2, label=metabolites[idx]["full"])
    plt.legend()
plt.show()

