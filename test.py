import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
sqrt = np.sqrt
import jol_parameters as p
import jol_lib as l
import sympy as sym
qAK = 0.01 # 0.92
A =2.212

atp = sym.symbols("atp") # np.linspace(0.01, 2.2, 1000)

# damp_datp = l.get_damp_datp(atp, qAK, A)
adp = l.get_adp(atp, qAK, A)

fun = sym.lambdify([atp], adp, modules=['numpy', 'sympy'])



a = fun(0.0000000001)

print(a)

# plt.plot(atp, damp_datp)
# plt.show()



# RTF = 26.73
# R       =   8.314510       # J mol-1 K-1
# F       =   9.64853e04     # C mol-1
# T = 310  # 37 grad Celcium
#
#
# Kext = 3.0
# Kcyt = 200.0
# Nacyt = 8.0
# Naext = 140.0
#
# gKpas = 0.232
# gNapas = 0.0136
#
#
# tmp = (gKpas * Kext + gNapas * Naext) / (gKpas * Kcyt + gNapas * Nacyt)
# E = RTF * np.log( tmp )
# print(E)
#
# V = -40
#
# U = V / RTF
# IK = U * gKpas * F * ( Kext - Kcyt * np.exp(U) ) / (1 - np.exp(U))
#
# print(    )




# EK = RTF * np.log(Kex / p.Kcyt)
# ENa = RTF * np.log(p.na_ext / Nacyt)
#
# El = p.gKpas * EK / (p.gKpas + p.gNan) + p.gNan * ENa / (p.gKpas + p.gNan)
#
# print(El)


"""
alpha = 0.1
beta = 10 # 0.015
gamma = 0.0225
delta = 0.02

def get_jac(t, y):

    J = np.zeros( [y.size, y.size], dtype=np.float)

    J[0, 0] = alpha - beta * y[1]
    J[0, 1] = -y[0] * beta

    J[1, 0] = delta * y[1]
    J[1, 1] = delta*y[0] - gamma

    return J

def model(t, z):
    x, y = z[0], z[1]
    dxdt =  x * (alpha - beta * y)
    dydt = -y * (gamma - delta * x)
    return [dxdt, dydt]

# res = get_jac(0, np.array( [0.1, 0.2] ) )
# print(res)


y0 = np.array([1.0, 1.0])
# solve ODE
sol = solve_ivp(model, [0, 1000], y0, method="LSODA", jac=get_jac)
print (sol.t.size)

X, Y = sol.y[0, :],sol.y[1, :]
plt.plot(sol.t, X,'k',linewidth = 5)
plt.plot(sol.t, Y,'r',linewidth = 5)
myleg = plt.legend(['X','Y'],loc='upper right',prop = {'size':28,'weight':'bold'}, bbox_to_anchor=(1,0.9))
plt.show()

"""


"""
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt


def get_jac(t, y):

    J = np.zeros( [y.size, y.size], dtype=np.float)

    J[0, 0] = 0.08*y[0] + 5
    J[0, 1] = -1

    J[1, 0] = 0.02
    J[1, 1] = -0.1

    print (J)
    return J



# function that returns dy/dt
def model(t, y):
    a = 0.1
    b = 0.2
    I = 5.5


    v = y[0]
    u = y[1]

    if v >= 30:
        y[0] = -65
        y[1] += 2

        return [0, 0]

    dvdt = 0.04 * v**2 + 5*v + 140 - u + I
    dudt = a*(b*v - u)

    return [dvdt, dudt]

# initial condition
y0 = [-65, -6.5]

# time points
t = np.linspace(0, 1000, 1000)

# solve ODE
sol = solve_ivp(model, [0, 1000], y0, method="LSODA", jac=get_jac)
print (sol.t.size)
# plot results
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.plot(sol.t,sol.y[0, :],'k',linewidth = 5)
plt.plot(sol.t,sol.y[1, :],'r',linewidth = 5)
myleg = plt.legend(['v','u'],loc='upper right',prop = {'size':28,'weight':'bold'}, bbox_to_anchor=(1,0.9))

plt.show()
"""
