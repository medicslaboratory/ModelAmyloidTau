# Date : 17 January 2022
# Autor : Éléonore Chamberland

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

AgeEnd = 31
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
# fig = plt.figure()
fig, axs = plt.subplots(nrows=4, ncols=5, sharex="all")
fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)

# Making a list for Label names in the plot
labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$', r'$\tau$',
             '$F_i$', '$F_o$', '$N$', 'A', '$M$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$', r'$\hat{M}_{anti}$',
             r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

plt.subplots_adjust(hspace=.8, wspace=.8)

# for i in range(0, 19):
#     ax = fig.add_subplot(4, 5, i + 1)
#
#     ax.plot(sol.t / 365, sol.y[i, :])
#     if decades > 5:
#         indentxaxis = int(decades / 2)
#     elif decades < 1:
#         indentxaxis = AgeEnd - AgeStart
#     else:
#         indentxaxis = decades
#     ax.set_xticks(np.linspace(AgeStart, AgeEnd, indentxaxis + 1))
#     ax.set(xlabel='Age (years)', ylabel=labelname[i])
#     formatter = ticker.ScalarFormatter(useMathText=True)
#     formatter.set_scientific(True)
#     formatter.set_powerlimits((-1, 1))
#     # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
#     # notation will be used if exp <= -1 or exp >= 1.
#     ax.yaxis.set_major_formatter(formatter)

# axs.tolist()
# axs = [item for sublist in axs for item in sublist]

i = 0
for ax in axs.flat:
    if i < 19:
        ax.plot(sol.t / 365, sol.y[i, :])  # , '.-', ms=2
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

plt.tight_layout()

# Write the initial values used
plt.subplots_adjust(bottom=0.16)
icNameValue = [str(labelname[i]) + "= " + "{:.2e}".format(y0[i]) for i in np.arange(19)]
initcond = "Initial conditions used (in g/mL) : \n" + ", ".join(icNameValue[:10]) + ", \n" + ", ".join(icNameValue[10:])
plt.text(0.03, 0.08, initcond, fontsize=9, ha='left', va='top', transform=plt.gcf().transFigure)  # , wrap=True

# Save the plot as a .png file
if p.S == 0:
    sex = "F"
else:  # p.S == 1:
    sex = "M"

if p.AP == 1:
    APOE = "+"
else:  # p.AP == 0:
    APOE = "-"

number = 1
date = "22-09-21"
my_path = os.path.abspath('Figures')
FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + method + "_APOE" + APOE + "_" + sex + "_" + \
          str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + "test.png"
while os.path.exists(os.path.join(my_path, FigName)):
    number = number+1
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + method + "_APOE" + APOE + "_" + sex + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + ".png"

plt.savefig(os.path.join(my_path, FigName), dpi=180)

"""Add information to the figure."""
FigInfos = {"max_step": str(max_step),
            "Début": "Début intégration",  # "Ignore la première demi-année."
            "Modification(s)": "Normal."}

im = Image.open("Figures/" + FigName)
Infos = PngImagePlugin.PngInfo()
for x in im.info:
    Infos.add_text(x, str(im.info[x]))
for x in FigInfos:
    Infos.add_text(x, FigInfos[x])
im.save("Figures/" + FigName, "png", pnginfo=Infos)

# im = Image.open("Figures/" + FigName)
# print(im.info)

# ## To access the infos by typing the following in the console :
# from PIL import Image
# im = Image.open("Path of the figure")  # "Path of the figure" is, for example, "Figures/Figure_22-09-08_..._01.png".
# im.info

plt.show()

