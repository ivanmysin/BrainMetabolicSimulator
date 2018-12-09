import numpy as np

import matplotlib.pyplot as plt
sqrt = np.sqrt
import jol_parameters as p
import jol_lib as l
plt.rc('axes', linewidth=2)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)
plt.rc('font', size=18)
plt.rc('lines', linewidth=1)
plt.rc('lines', markersize=4)
plt.rc('lines', color="black")


qAK = p.qAK
A = p.totalAdenine

atp = np.linspace(0.1, A, 1000)

damp_datp = 1 - l.get_damp_datp(atp, qAK, A)
adp = l.get_adp(atp, qAK, A)

fig, ax = plt.subplots(ncols=2, nrows=1, constrained_layout=True, figsize=(10, 5))
ax[0].plot(atp, damp_datp, linewidth=2, color="black")
ax[0].set_xlabel("АТФ, мМ")
ax[0].set_ylabel(r"$1 - \frac{d[АТФ]}{d[АМФ]}$", fontsize=24)

ax[1].plot(atp, adp, linewidth=2, color="black")
ax[1].set_xlabel("АТФ, мМ")
ax[1].set_ylabel("АДФ, мМ")

fig.savefig("./figures/fig6.png", dpi=100)
plt.show()