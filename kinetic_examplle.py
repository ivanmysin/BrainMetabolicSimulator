from operator import mul
from functools import reduce
import sympy as sym
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def prod(seq):
    return reduce(mul, seq) if seq else 1


def mk_exprs_symbs(rxns, names):
    # create symbols for reactants
    symbs = sym.symbols(names, real=True, nonnegative=True)
    # map between reactant symbols and keys in r_stoich, net_stoich
    c = dict(zip(names, symbs))
    f = {n: 0 for n in names}
    k = []

    for coeff, r_stoich, net_stoich in rxns:
        k.append(sym.S(coeff))
        r = k[-1]*prod([c[rk]**p for rk, p in r_stoich.items()])  # EXERCISE: c[rk]**p
        for net_key, net_mult in net_stoich.items():
            f[net_key] += net_mult*r  # EXERCISE: net_mult*r

    return [f[n] for n in names], symbs, tuple(k)

reactions = [
    # (coeff, r_stoich, net_stoich)
    ('k1', {'A': 1}, {'B': 1, 'A': -1}),
    ('k2', {'B': 1, 'C': 1}, {'A': 1, 'B': -1}),
    ('k3', {'B': 2}, {'B': -1, 'C': 1})
]
names = 'A B C'.split()

sym.init_printing()
ydot, y, k = mk_exprs_symbs(reactions, names)
print(ydot)

t = sym.symbols('t')  # not used in this case.
f = sym.lambdify((y, t) + k, ydot)

# %exercise exercise_lambdify_jac.py
J = sym.Matrix(ydot).jacobian(y)  # EXERCISE: jacobian
J_cb = sym.lambdify((y, t) + k, J)  # EXERCISE: (y, t) + k

tout = np.logspace(-6, 6)
k_vals = (0.04, 1e4, 3e7)  # from the literature
y0 = [1, 0, 0]

yout, info = odeint(f, y0, tout, k_vals, full_output=True, Dfun=J_cb)
plt.loglog(tout, yout)
plt.legend(names)
print("The Jacobian was evaluated %d times." % info['nje'][-1])
plt.show()