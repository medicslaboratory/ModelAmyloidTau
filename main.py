# Date : 17 January 2022
# Author : Éléonore Chamberland

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import os
from scipy.integrate import solve_ivp
from PIL import Image
from PIL import PngImagePlugin
import equations as eqns
from InitialConditions import InitialConditions
import parameters as param

p = param.Parameters()

AgeStart = 30

y0 = InitialConditions(AgeStart)

AgeEnd = 40
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
fig, axs = plt.subplots(nrows=4, ncols=5, sharex="all")  # If no text under the figure, add: layout="constrained"
fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)

# Making a list for Label names in the plot
labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$', r'$\tau$',
             '$F_i$', '$F_o$', '$N$', '$A$', '$M_{NA}$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$', r'$\hat{M}_{anti}$',
             r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

k = np.argmax(sol.t >= ((AgeStart+0.5)*365))  # First indice to skip the first half year.

i = 0
for ax in axs.flat:
    if i < 19:
        # # """Plot tout."""
        ax.plot(sol.t / 365, sol.y[i, :])  # , '.-', ms=2
        # # """Plot sans la première demie année."""
        # ax.plot(sol.t[k:] / 365, sol.y[i, k:])
        ax.grid()
        if i >= 14:
            ax.set_xlabel('Age (years)')
        ax.set_ylabel(labelname[i])
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
        # notation will be used if exp <= -1 or exp >= 1.
        ax.yaxis.set_major_formatter(formatter)
    i = i+1
axs.flat[14].tick_params('x', labelbottom=True)

axs[3, 4].remove()

"""Pour afficher certains les graphs de concentration intraneuronale."""
# """Add the plot of ABi/N in the graph for ABi"""
# axs[0, 0].set_ylabel(labelname[0], color='b')
# ax1 = axs[0, 0].twinx()
# ax1.plot(sol.t / 365, sol.y[0, :]/sol.y[8, :], "g-")  # ABi/N
# ax1.set_ylabel(r"$A \beta^{i} / N$", color='g')
# formatter1 = ticker.ScalarFormatter(useMathText=True)
# formatter1.set_scientific(True)
# formatter1.set_powerlimits((-1, 1))
# ax1.yaxis.set_major_formatter(formatter1)
#
# """Add the plot of GSK3/N in the graph for GSK3"""
# axs[0, 4].set_ylabel(labelname[4], color='b')
# ax2 = axs[0, 4].twinx()
# ax2.plot(sol.t / 365, sol.y[4, :] / sol.y[8, :], "g-")  # G/N
# ax2.set_ylabel(r"$GSK3/N$", color='g')
# formatter2 = ticker.ScalarFormatter(useMathText=True)
# formatter2.set_scientific(True)
# formatter2.set_powerlimits((-1, 1))
# ax2.yaxis.set_major_formatter(formatter2)
#
# """Add the plot of tau/N in the graph for tau"""
# axs[1, 0].set_ylabel(labelname[5], color='b')
# ax3 = axs[1, 0].twinx()
# ax3.plot(sol.t / 365, sol.y[5, :] / sol.y[8, :], "g-")  # tau/N
# ax3.set_ylabel(r"$\tau/N$", color='g')
# formatter3 = ticker.ScalarFormatter(useMathText=True)
# formatter3.set_scientific(True)
# formatter3.set_powerlimits((-1, 1))
# ax3.yaxis.set_major_formatter(formatter3)

"""Pour afficher un graph à la dernière position 
(doit retirer le 'axs[3, 4].remove()' et 'axs.flat[14].tick_params('x', labelbottom=True)' plus haut)."""
# ax = axs[3, 4]
# dABpodt = p.kappa_ABooABpo * (sol.y[2, :] ** 2) - ((p.d_MantiABpo * sol.y[12, :] + p.d_hatMantiABpo * sol.y[14, :])
#                                                    * (1 + p.AP * p.delta_APdp) * (sol.y[3, :] / (sol.y[3, :]
#                                                                                                  + p.K_ABpo)))
# ax.plot(sol.t / 365, dABpodt)  # , '.-', ms=2
# dNdtFi = - p.d_FiN * (1 / (1 + np.exp(- p.n * (sol.y[6, :] - p.K_Fi)))) * sol.y[8, :]
# ax.plot(sol.t / 365, dNdtFi)  # , '.-', ms=2
# ax.grid()
# ax.set_xlabel('Age (years)')
# # ax.set_ylabel("dABpo/dt")
# ax.set_ylabel("dN/dt du à F_i")
# ax2 = axs[3, 4].twinx()
# axs[3, 4].plot(sol.t / 365, sol.y[12, :], "b-")  # M_anti
# axs[3, 4].set_ylabel("Manti", color='b')
# ax2.plot(sol.t / 365, sol.y[16, :], "g-")  # I_10
# ax2.set_ylabel("I10", color='g')

plt.subplots_adjust(hspace=.2, wspace=.2)
plt.tight_layout()

# plt.tight_layout(rect=(0, 0.1, 1, 1))

"""Write the initial values used"""
plt.subplots_adjust(bottom=0.16)
icNameValue = [str(labelname[i]) + " = " + "{:.2e}".format(y0[i]) for i in np.arange(19)]
initcond = "Initial conditions used (in g/mL) : \n" + ", ".join(icNameValue[:10]) + ", \n" + ", ".join(icNameValue[10:])
plt.text(0.03, 0.08, initcond, fontsize=9, ha='left', va='top', transform=plt.gcf().transFigure)  # , wrap=True

"""Write the full name of the variables"""
# plt.subplots_adjust(bottom=0.16)
# VariablesDef = r'$A \beta^{i}$: Intracellular beta-amyloid monomer; ' \
#                r'$A \beta_{m}^{o}$: Extracellular beta-amyloid monomer; ' \
#                r'$A \beta_{o}^{o}$: Extracellular beta-amyloid oligomer; ' \
#                r'$A \beta_{p}^{o}$: Extracellular beta-amyloid plaques; ' \
#                r'$GSK 3$: Glycogen synthase kinase-3; ' + '\n' + r'$\tau$: Hyperphosphorylated tau protein; ' \
#                r'$F_i$: Intracellular NFTs; $F_o$: Extracellular NFTs; ' \
#                r'$N$: Neurons; $A$: Astrocytes; $M_{NA}$: Resting microglia; ' \
#                r'$M_{pro}$ Proinflammatory microglia; $M_{anti}$: Anti-inflammatory microglia; ' + '\n' + \
#                r'$\hat{M}_{pro}$: Proinflammatory macrophages; ' \
#                r'$\hat{M}_{anti}$: Anti-inflammatory macrophages; $T_{\beta}$: TGF-$\beta$; ' \
#                r'$I_{10}$: Interleukin-10; $T_{\alpha}$: TNF-$\alpha$; $P$: MCP-1.'
# plt.text(0.03, 0.08, VariablesDef, fontsize=9, ha='left', va='top', transform=plt.gcf().transFigure)

"""Save the plot as a .png file"""
if p.S == 0:
    sex = "F"
else:  # p.S == 1:
    sex = "M"

if p.AP == 1:
    APOE = "+"
else:  # p.AP == 0:
    APOE = "-"

number = 5
date = "22-09-22"
my_path = os.path.abspath('Figures')
FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + method + "_APOE" + APOE + "_" + sex + "_" + \
          str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + "_.png"
# d_TaN=001on365ETprodGmultNsurN0ET2*transfosETDegGavecN
while os.path.exists(os.path.join(my_path, FigName)):
    number = number+1
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + method + "_APOE" + APOE + "_" + sex + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + ".png"

plt.savefig(os.path.join(my_path, FigName), dpi=180)

"""Add information to the figure."""
FigInfos = {"max_step": str(max_step),
            "Début": "Début intégration",  # "Ignore la première demi-année."
            "Modification(s)": "d_TaN = 0.01/365. p.lambda_InsG * (p.Ins_0 / p.Ins(t, p.S)) * (y[8]/p.N_0). transfo aggrédation *2. Ajout '- (y[4] / y[8]) * abs(dydt[8])' à eqn GSK3."}

im = Image.open("Figures/" + FigName)
Infos = PngImagePlugin.PngInfo()
for x in im.info:
    Infos.add_text(x, str(im.info[x]))
for x in FigInfos:
    Infos.add_text(x, FigInfos[x])
im.save("Figures/" + FigName, "png", pnginfo=Infos)

# im = Image.open("Figures/" + FigName)
# print(im.info)

# # To access the infos by typing the following in the console :
# from PIL import Image
# im = Image.open("Path of the figure")  # "Path of the figure" is, for example, "Figures/Figure_22-09-08_..._01.png".
# im.info

# fig, axs = plt.subplots(1, 2)
# dNdtFi = - p.d_FiN * (1 / (1 + np.exp(- p.n * (sol.y[6, :] - p.K_Fi)))) * sol.y[8, :]
# dNdtTa = - p.d_TaN * (sol.y[17, :] / (sol.y[17, :] + p.K_Ta)) * (1 / (1 + (sol.y[16, :] / p.K_I10))) * sol.y[8, :]
#
# axs[0].plot(sol.t / 365, dNdtFi)
# axs[1].plot(sol.t / 365, dNdtTa)

# fig2, axs = plt.subplots(1, 3, sharex="all", layout="constrained")
# ax1, ax2, ax3 = axs.flat
# ABiperN = sol.y[0, :] / sol.y[8, :]
# GperN = sol.y[4, :] / sol.y[8, :]
# tauperN = sol.y[5, :] / sol.y[8, :]
# ax1.plot(sol.t / 365, ABiperN)
# ax1.set_ylabel(r"$A\beta^i/N$")
# ax1.set_xlabel("Age (years)")
#
# ax2.plot(sol.t / 365, GperN)
# ax2.set_ylabel(r"$GSK3/N$")
# ax2.set_xlabel("Age (years)")
#
# ax3.plot(sol.t / 365, tauperN)
# ax3.set_ylabel(r"$\tau/N$")
# ax3.set_xlabel("Age (years)")
#
# # formatter = ticker.ScalarFormatter(useMathText=True)
# # formatter.set_scientific(True)
# # formatter.set_powerlimits((-1, 1))
# # # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
# # # notation will be used if exp <= -1 or exp >= 1.
# # ax1.yaxis.set_major_formatter(formatter)
# # ax2.yaxis.set_major_formatter(formatter)
# # ax3.yaxis.set_major_formatter(formatter)

plt.show()

