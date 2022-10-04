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

AgeEnd = 80
decades = int((AgeEnd - AgeStart) / 10)

max_step = 0.01
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
#
# """Add the plot of F_i/N in the graph for F_i"""
# axs[1, 1].set_ylabel(labelname[6], color='b')
# ax4 = axs[1, 1].twinx()
# ax4.plot(sol.t / 365, sol.y[6, :] / sol.y[8, :], "g-")  # tau/N
# ax4.set_ylabel(r"$F_i/N$", color='g')
# formatter4 = ticker.ScalarFormatter(useMathText=True)
# formatter4.set_scientific(True)
# formatter4.set_powerlimits((-1, 1))
# ax4.yaxis.set_major_formatter(formatter4)

"""Pour afficher un graph à la dernière position 
(doit retirer le 'axs[3, 4].remove()' et 'axs.flat[14].tick_params('x', labelbottom=True)' ci-dessous)."""
axs.flat[14].tick_params('x', labelbottom=True)
axs[3, 4].remove()
# ax = axs[3, 4]
# # # dABpodt = p.kappa_ABooABpo * (sol.y[2, :] ** 2) - ((p.d_MantiABpo * sol.y[12, :] + p.d_hatMantiABpo * sol.y[14, :])
# # #                                                    * (1 + p.AP * p.delta_APdp) * (sol.y[3, :] / (sol.y[3, :]
# # #                                                                                                  + p.K_ABpo)))
# # # ax.plot(sol.t / 365, dABpodt)  # , '.-', ms=2
# # # dNdtFi = - p.d_FiN * (1 / (1 + np.exp(- p.n * (sol.y[6, :] - p.K_Fi)))) * sol.y[8, :]
# # # ax.plot(sol.t / 365, dNdtFi)  # , '.-', ms=2
# ax.grid()
# # # ax.set_xlabel('Age (years)')
# # # ax.set_ylabel("dABpo/dt")
# # # ax.set_ylabel("dN/dt du à F_i")
# ax2 = axs[3, 4].twinx()
# # # axs[3, 4].plot(sol.t / 365, sol.y[12, :], "b-")  # M_anti
# # # axs[3, 4].set_ylabel("Manti", color='b')
# # # ax2.plot(sol.t / 365, sol.y[16, :], "g-")  # I_10
# # # ax2.set_ylabel("I10", color='g')
# axs[3, 4].plot(sol.t / 365, sol.y[11, :], "b-")  # M_pro
# # axs[3, 4].set_ylabel("M_pro", color='b')
# axs[3, 4].plot(sol.t / 365, sol.y[13, :], "r-")  # hat{M}_pro
# axs[3, 4].set_ylabel("M_pro (blue); hat{M}_pro (red)")
# ax2.plot(sol.t / 365, sol.y[17, :], "g-")  # T_alpha
# ax2.set_ylabel("T_a", color='g')

plt.subplots_adjust(hspace=.2, wspace=.2)
plt.tight_layout()

# """Peut-être utiliser cela plutôt que les deux lignes précédentes si overlap. Ne l'utilise pas de base, car produit
# des graphs plus petits (plus distancés)."""
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

number = 3
date = "22-10-04"
CommentModif = "_n=15&K_TaM*2e2&K_Fi=1.708e-10&kappa_PMhat*1e-1&K_P*1e2&Fo0=5e-17&CondInitModif"
# &kappa_TaMhatanti*1e-2
# "_n=10ETd_FiNe-5ETFo0=0ETK_TaM*2e2ETK_Fi*1e-6ETkappa_MproTaFadok98minETkappa_PMhat*1e-1_3"
# ETkappaTalphaHallswoth94
# ETAjustCondInitMantiethatMantipourTaetI10egaux
# conserve : d_TaN=001on365ETprodGmultNsurN0ET2*transfosETDegGavecN
# À partir de 22-09-28_5..., conserve : d_TaN=7e-5on365 et F_i0ettau0equilibre
# À partir de 22-09-29_8..., conserve : kappa_MproTaFadok98min # (impact aussi kappa_MhatproTa)
my_path = os.path.abspath('Figures')
FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + method + "_APOE" + APOE + "_" + sex + "_" + \
          str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + CommentModif + ".png"
while os.path.exists(os.path.join(my_path, FigName)):
    number = number+1
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + method + "_APOE" + APOE + "_" + sex + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + \
              CommentModif + ".png"

plt.savefig(os.path.join(my_path, FigName), dpi=180)

"""Add information to the figure."""
FigInfos = {"max_step": str(max_step),
            "Début": "Début intégration",  # "Ignore la première demi-année."
            "Modification(s)": "d_TaN = 7e-5 / 365. p.lambda_InsG * (p.Ins_0 / p.Ins(t, p.S)) * (y[8]/p.N_0). "
                               "transfo aggrédation *2. Ajout '- (y[4] / y[8]) * abs(dydt[8])' à eqn GSK3. "
                               "n = 15. "
                               "[F_i]_0 = y0[6] = équilibre, et tau aussi. "
                               # "d_FiN = 1 / (20 * 365) (Kril et al. 2002). "
                               "d_FiN = 1 / (2.51 * 365) (original). "
                               "[F_o]_0 = 5e-17. " 
                               # "K_Fi = 0.1 * (0.6 * (6e-3 * self.rho_cerveau)) * 1e-6 ~= 3.708e-10. "
                               "K_Fi = 1.708e-10. "
                               "K_TaM = 2.24e-12 * 2e2, impact aussi K_TaMhat. "
                               "kappa_MhatproTa = (1.5e-9 / 18 * 24) / (4e6 * m_Mhat), min value from Fadok98 (impact aussi kappa_MproTa). "
                               # "Cond init Manti et hat{M}_anti pour que cond init Ta et I_10 soit égales."
                               "kappa_PMhat = 0.33 * 1e-1."
                               # "kappa_TaA = 0.92 / 100e-9 * 1e-2. "
                               # "kappa_TaManti = 4.8 * 1e-1. "
                               # "kappa_TaMhatanti = 1 / (10 / 24) * 1e-2. "
                               "K_P = 6.23e-10 * 1e2. "
                               # "K_ABpo = (1.11 + 0.53) / 527.4 / 1000 * 1e-4. "
                               "[M_anti]_0 = 1e-5. "
                               # "[AB_p^o]_0 = 1e-29. "
                               "[AB_m^o]_0, [AB_o^o]_0 et [AB_p^o]_0 à l'équilibre. "
                               # "[AB_m^o]_0 = 4e-11, [AB_o^o]_0 = 6e-17 et [AB_p^o]_0 = 5e-28. "
                               # "self.K_Fo = 11 * ((1000 * 72500) / Avogadro) * 1000 * 5. "
                               # "K_ABooM = 0.060 / 527.4 / 1000 * 1e-5. "
                               # "[M_NA]_0 = (1 - 0.67/100) * M_0 - (y0[11] + y0[12]), [hat{M}_pro]_0 = [hat{M}_anti]_0 = 0.67/100 * M_0 * 0.5. Impact aussi plaques et cytokines. "
            }

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

"""Trace les graphs des pertes neuronales par chaque cause (F_i et TNFa) dans une figure indépendante."""
# dNdtFi = - p.d_FiN * (1 / (1 + np.exp(- p.n * (sol.y[6, :] - p.K_Fi)))) * sol.y[8, :]
# dNdtTa = - p.d_TaN * (sol.y[17, :] / (sol.y[17, :] + p.K_Ta)) * (1 / (1 + (sol.y[16, :] / p.K_I10))) * sol.y[8, :]
#
# # # Pour un graph.
# fig, ax1 = plt.subplots(1, 1)
# # # Un graph : Deux axes x.
# ax1.plot(sol.t / 365, dNdtFi, "b-")  # loss by F_i
# ax1.set_ylabel(r"$dN/dt$ par $F_i$", color="b")
# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True)
# formatter.set_powerlimits((-1, 1))
# # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
# # notation will be used if exp <= -1 or exp >= 1.
# ax1.yaxis.set_major_formatter(formatter)
# ax2 = ax1.twinx()
# ax2.plot(sol.t / 365, dNdtTa, "g-")  # loss by TNFa
# ax2.set_ylabel(r"$dN/dt$ par $T_\alpha$", color='g')
# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True)
# formatter.set_powerlimits((-1, 1))
# # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
# # notation will be used if exp <= -1 or exp >= 1.
# ax2.yaxis.set_major_formatter(formatter)
# ax1.set_xlabel("Age (years)")
# ax1.grid()
#
# # # Un graph : Une axe x. Recommande moins.
# # ax1.plot(sol.t / 365, dNdtFi, "b-", label=r"Par $F_i$")  # loss by F_i
# # ax1.plot(sol.t / 365, dNdtTa, "g-", label=r"Par $T_\alpha$")  # loss by TNFa
# # ax1.legend()
# # ax1.set_ylabel(r"$dN/dt$ par facteur")
# # ax1.set_xlabel("Age (years)")
# # ax1.grid()
#
# # # Pour deux graphs distincts
# # fig, axs = plt.subplots(1, 2)
# # axs[0].plot(sol.t / 365, dNdtFi)
# # axs[0].set_ylabel(r"$dN/dt$ par $F_i$")
# # axs[1].plot(sol.t / 365, dNdtTa)
# # axs[1].set_ylabel(r"$dN/dt$ par $T_i$")
#
#  plt.tight_layout()

""" Autre"""
fig, ax = plt.subplots(1, 1)
M_activ = p.kappa_FoM * (sol.y[7, :] / (sol.y[7, :] + p.K_Fo)) * sol.y[10, :] + \
          p.kappa_ABooM * (sol.y[2, :] / (sol.y[2, :] + p.K_ABooM)) * sol.y[10, :]
epsilon_Ta = sol.y[17, :] / (sol.y[17, :] + p.K_TaAct)
epsilon_I10 = sol.y[16, :] / (sol.y[16, :] + p.K_I10Act)
ax.plot(sol.t / 365, ((p.beta * epsilon_Ta) / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ, "b-", label="Activ pro")
ax.plot(sol.t / 365, (epsilon_I10 / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ, "r-", label="Activ anti")
ax.set_ylabel("Rates")
ax.plot(sol.t / 365, M_activ, "k", label="M_activ")
ax.plot(sol.t / 365, p.kappa_FoM * (sol.y[7, :] / (sol.y[7, :] + p.K_Fo)) * sol.y[10, :], label="Mactiv F_o")
ax.plot(sol.t / 365, p.kappa_ABooM * (sol.y[2, :] / (sol.y[2, :] + p.K_ABooM)) * sol.y[10, :], label="Mactiv ABoo")
ax.grid()
ax.set_xlabel('Age (years)')
# ax2 = ax.twinx()
AntiToPro = p.kappa_TaManti * (sol.y[17, :] / (sol.y[17, :] + p.K_TaM)) * sol.y[12, :]
ProToAnti = p.kappa_TbMpro * (sol.y[15, :] / (sol.y[15, :] + p.K_TbM)) * sol.y[11, :]
ax.plot(sol.t / 365, AntiToPro, "g-", label="rate anti -> pro")
ax.plot(sol.t / 365, ProToAnti, "m-", label="rate pro -> anti")
# ax2.plot(sol.t / 365, AntiToPro, "g-", label="rate anti -> pro")
# ax2.set_ylabel("Rate anti to pro", color='g')
ax.legend()
plt.tight_layout()

plt.show()

