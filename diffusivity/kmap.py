#!/usr/bin/env python
# Time-stamp: "2015-12-03 16:58:29 marine"

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

import eos
import kmodel

# Light elements:
# name, z: valence charge, x: atomic fraction
# f1,f2,f3 are fitting parameters
IRON = {'name': 'Fe', 'f1': 5.26e-9, 'f2': 1.24, 'f3': -3.21}
SILICIUM = {'name': 'Si', 'z': 1, 'x': 22.5,
            'f1': 3.77e-8, 'f2': 1.48, 'f3': -3.10}
CARBON = {'name': 'C', 'z': 1, 'x': 30}
OXYGENn = {'name': 'O', 'z': 0.5, 'x': 23.2}
SULFUR = {'name': 'S', 'z': 0.5, 'x': 19.4}

MURNAGHAN = {'T0': 1812., 'rho0': 7010., 'K': 130.*1e9, 'KPrim': 4.,
                        'alphaPrim': 130.*np.log(2.)/(360.-135.),
                        'alpha': 1.e-5 / np.exp(-0.400485 * 135./130.),
                        'Cp': 800., 'cond': 140.}


N = 100
p = np.linspace(130e9, 450e9, N)
t = np.linspace(2000., 5000., N)

P, T = np.meshgrid(p, t)

conductivity = np.zeros(np.shape(P))


for i in range(0, N):
    for j in range(0, N):
        f = eos.fval(P[i, j], T[i, j], MURNAGHAN)
        res = kmodel.resistivity(f, T[i, j], El=SILICIUM, theta0=417., gamma=1.52, Ir=IRON)
        conductivity[i, j] = kmodel.k(res, T[i, j], f, L=2.44e-8)

# adiabatic profile
fig = plt.figure(1)
ax = fig.add_subplot(111)
im = ax.contourf(P/1e9, T, conductivity, 20)

P0, T0 = 135e9, 3750.

for i in range(0, 10):
    Tad = []
    T0 = 2000. + i*300.
    for p_i in p:
        # print Tadiabatic(p, T0, P0, MURNAGHAN)
        Tad = np.append(Tad, eos.Tadiabatic(p_i, T0, P0, MURNAGHAN))
    ax.plot(p/1e9, Tad, color="black")
    ## ax.annotate("isentropic profile", xy=(p[int(N/2)]/1e9, Tad[int(N/2)]), xycoords='data',
    ##             xytext=(-30, 30), textcoords='offset points',
    ##             arrowprops=dict(arrowstyle='->'))

cbar = plt.colorbar(im, ax=ax)

ax.axis([130., 450., 2000., 5000.])
ax.set_xlabel("Pressure (GPa)")
ax.set_xlabel("Temperature (K)")

plt.savefig("map_conductivity.pdf")
plt.show()
