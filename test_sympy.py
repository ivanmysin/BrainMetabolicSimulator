import sympy as sym
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

t = sym.symbols("t")
tau = 10.5

tstim = t - 10
g = 0.0599 * 7000 * sym.exp(-tstim / tau)

dVdt = sym.Piecewise((0, t < 10), (g, (t >= 10) & (t <= 30) ), (0, t > 30)  )

f = sym.lambdify(t, dVdt, "numpy")

t_vec = np.linspace(0, 100, 1000)
g_vec = f(t_vec)
# print(g_vec)
plt.plot(t_vec, g_vec)
plt.show()


"""
alpha = 0.1
beta = 0.015 
gamma = 0.0225
delta = 0.02

dxdt, dydt, x, y, t = symbols("dxdt, dydt, x, y, t")

dxdt = x * (alpha - beta * y)
dydt = -y * (gamma - delta * x)

func = lambdify([t, [x, y] ], [dxdt, dydt], "numpy")

J = Matrix([dxdt, dydt]).jacobian([x, y])
Jfunc = lambdify([t, [x, y]], J)

y0 = np.array([1.0, 1.0])

# # solve ODE
sol = solve_ivp(func, [0, 1000], y0, method="LSODA", jac=Jfunc)
print (sol.t.size)

X, Y = sol.y[0, :],sol.y[1, :]
plt.plot(sol.t, X,'k',linewidth = 2)
plt.plot(sol.t, Y,'r',linewidth = 2)
myleg = plt.legend(['X','Y'],loc='upper right',prop = {'size':28,'weight':'bold'}, bbox_to_anchor=(1,0.9))
plt.show()

"""



