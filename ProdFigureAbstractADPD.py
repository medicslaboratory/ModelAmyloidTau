# Date : 21 September 2022
# Autor : Éléonore Chamberland

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import ticker
import os

from scipy.integrate import solve_ivp

import equations as eqns
from InitialConditions import InitialConditions
import parameters as param

p = param.Parameters()

AgeStart = 30

y0 = InitialConditions(AgeStart)

AgeEnd = 80
decades = int((AgeEnd - AgeStart) / 10)

max_step = 0.1
maxstepstr = str(max_step).replace('.', '')
rtol = 1e-10  # Default value : 1e-3
rtolstr = "{:.0e}".format(rtol)
method = "BDF"
# method = "Radau"
# method = "LSODA"
sol = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0, method=method, max_step=max_step, rtol=rtol)

"""Generate the figure"""
fig, axs = plt.subplots(nrows=3, ncols=3, sharex="all")  # If no text under the figure, add: layout="constrained"
fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)

# Making a list for Label names in the plot
labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$', r'$\tau$',
             '$F_i$', '$F_o$', '$N$', '$A$', '$M_{NA}$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$',
             r'$\hat{M}_{anti}$',
             r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

graphs = [0, 1, 2, 3, 4, 6, 7, 8, 17]

k = np.argmax(sol.t >= ((AgeStart+0.5)*365))  # First indice to skip the first half year.

ax = axs.flat
j = 0
for i in graphs:
    ax[j].plot(sol.t[k:] / 365, sol.y[i, k:])
    ax[j].grid()
    if i >= 6:
        ax[j].set_xlabel('Age (years)')
    ax[j].set_ylabel(labelname[i])
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    ax[j].yaxis.set_major_formatter(formatter)
    j = j+1

plt.subplots_adjust(hspace=.8, wspace=.8)
plt.tight_layout()

"""Write the full name of the variables"""
plt.subplots_adjust(bottom=0.16)
VariablesDef = [r'$A\beta^{i}$: Intracellular beta-amyloid monomer',
                r'$A\beta_{m}^{o}$: Extracellular beta-amyloid monomer',
                r'$A\beta_{o}^{o}$: Extracellular beta-amyloid oligomer',
                r'$A\beta_{p}^{o}$: Extracellular beta-amyloid plaques',
                r'$GSK 3$: Glycogen synthase kinase-3', r'$\tau$: Hyperphosphorylated tau protein',
                r'$F_i$: Intracellular NFTs', r'$F_o$: Extracellular NFTs',
                r'$N$: Neurons', r'$A$: Astrocytes', r'$M_{NA}$: Resting microglia',
                r'$M_{pro}$ Proinflammatory microglia', r'$M_{anti}$: Anti-inflammatory microglia',
                r'$\hat{M}_{pro}$: Proinflammatory macrophages',
                r'$\hat{M}_{anti}$: Anti-inflammatory macrophages', r'$T_{\beta}$: TGF-$\beta$',
                r'$I_{10}$: Interleukin-10', r'$T_{\alpha}$: TNF-$\alpha$', '$P$: MCP-1']
VariablesDefText = "; ".join([VariablesDef[i] for i in graphs]) + "."
plt.text(0.03, 0.08, VariablesDefText, fontsize=9, ha='left', va='top', transform=plt.gcf().transFigure, wrap=True)

if p.S == 0:
    sex = "F"
else:  # p.S == 1:
    sex = "M"

if p.AP == 1:
    APOE = "+"
else:  # p.AP == 0:
    APOE = "-"

number = 1
my_path = os.path.abspath('FiguresAbstractADPD')
FigName = "FigureAbstract_" + f"{number:02}" + "_" + method + "_APOE" + APOE + "_" + sex + "_" + \
          str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + ".png"

plt.savefig(os.path.join(my_path, FigName), dpi=180)

plt.show()
