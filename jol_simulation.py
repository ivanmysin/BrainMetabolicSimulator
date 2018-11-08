import numpy as np
import sympy as sym
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
agents.append( lib.LactateTransport(10, 26, params["Lactate exchange ne"]) )
agents.append( lib.LactateTransport(11, 26, params["Lactate exchange ge"]) )
agents.append( lib.LactateTransport(26, 22, params["Lactate exchange ec"]) )
agents.append( lib.LactateTransport(11, 22, params["Lactate exchange gc"]) )
agents.append( lib.TCA(8, 31, 14, params["TCA n"]) )
agents.append( lib.TCA(9, 32, 15, params["TCA g"]) )
agents.append( lib.ETC(18, 14, 31, params["Mitochondrial respiration n"]) )
agents.append( lib.ETC(19, 15, 32, params["Mitochondrial respiration g"]) )
agents.append( lib.NADHShuttle(12, 31, params["NADH Shuttles n"]) )
agents.append( lib.NADHShuttle(13, 32, params["NADH Shuttles n"]) )
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
agents.append( lib.CalciumCurrent(27, params["calcium current"]) )
agents.append( lib.AHPCurrent(27, 30, params["AHP current"]) )
agents.append( lib.CalciumDecay(30, params["calcium decay"]) )


short_names = [met["short"] for met in vars ]
metabls = sym.symbols(short_names, real=True) # , nonnegative=True
t = sym.symbols("t")


fsx_expr = [ 0 for _ in range(len(short_names)) ]

for en in agents:
    fsx_expr = en.update(metabls, fsx_expr, t)


fsx_expr[14] = fsx_expr[14] / (1 - lib.get_damp_datp(metabls[14], glob_params["qAK"], glob_params["A"]))
fsx_expr[14] = fsx_expr[15] / (1 - lib.get_damp_datp(metabls[15], glob_params["qAK"], glob_params["A"]))




y0 = np.asarray( [vars[idx]["rest"] for idx in range(len(vars))] )
