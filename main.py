# Date : 7 Octobre 2022
# Author : Éléonore Chamberland

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import os
from scipy.integrate import solve_ivp
from PIL import Image
from PIL import PngImagePlugin
import datetime
import equations as eqns
from InitialConditions import InitialConditions
import parameters as param
import figuresfunctions as ff

p = param.Parameters()

AgeStart = 30

y0 = InitialConditions(AgeStart)

AgeEnd = 80
decades = int((AgeEnd - AgeStart) / 10)

max_step = 0.5
maxstepstr = str(max_step).replace('.', '')
# rtol = np.ones(19) * 1e-10
# rtol[3] = 1e-29
# rtolstr = "Array"
rtol = 1e-10  # Default value : 1e-3
rtolstr = "{:.0e}".format(rtol)

atol = np.ones(19) * 1e-15
atol[3] = 1e-35
atolstr = "array"

method = "BDF"
# method = "Radau"
# method = "LSODA"
sol = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0, method=method, max_step=max_step, rtol=rtol, atol=atol)

number = 3

# 22-10-05_20 pas mal, perte environ -11,5%
# 22-10-06_15 : ~-6.178%
# 22-10-06_21 : perte environ ~-10.77%
# 22-10-06_33 : perte environ -10.628%  n=15&K_TaM*2e2&K_Fi=1.25e-10&kappa_PMhat*1e-2&K_P*1e2&Fo0=5e-17&lambda_ABiLindstrom21&d_TaNmodif*10
# 22-10-06_34 : perte environ -11.3947%  " " K_Fi=1.2e-10 " "

# 22-10-07_01 (même que 22-10-06_33) (K_Fi=1.2e-10): -11.394726350603838 %
# 22-10-07_02 (K_Fi=1.25e-10) : -10.628227731913643 %
# 22-10-07_03 (comme _02, mais APOE-) : -10.627498374162005 %

CommentModif = "comme02"
# CommentModif = "n=15&K_TaM*2e2&K_Fi=1.2e-10&kappa_PMhat*1e-2&K_P*1e2&Fo0=5e-17&lambda_ABiLindstrom21&d_TaNmodif*10"
# CommentModif = "_n=15&K_TaM*2e2&K_Fi=1.25e-10&kappa_PMhat*1e-2&K_P*1e2&Fo0=5e-17&CondInitModif&lambda_ABiLindstrom21" \
#                "&d_TaNmodif*10_atol=array"
# &kappa_TaMhatanti*1e-2   &d_TaNmodif*10 &d_TaNmodif*5 &K_Fo=16*-
# "_n=10ETd_FiNe-5ETFo0=0ETK_TaM*2e2ETK_Fi*1e-6ETkappa_MproTaFadok98minETkappa_PMhat*1e-1_3"
# ETkappaTalphaHallswoth94
# ETAjustCondInitMantiethatMantipourTaetI10egaux
# conserve : d_TaN=001on365ETprodGmultNsurN0ET2*transfosETDegGavecN
# À partir de 22-09-28_5..., conserve : d_TaN=7e-5on365 et F_i0ettau0equilibre
# À partir de 22-09-29_8..., conserve : kappa_MproTaFadok98min # (impact aussi kappa_MhatproTa)

solt = sol.t / 365

FigName = ff.main_figure(solt, sol.y, y0, AgeStart, AgeEnd, method, maxstepstr, rtolstr, atolstr, number, CommentModif)


"""Add information to the figure."""
FigInfos = {"max_step": str(max_step),
            "Début": "Début intégration",  # "Ignore la première demi-année."
            "Modification(s)": "d_TaN = 7.26e-3 / 365 * 10. "
                               "p.lambda_InsG * (p.Ins_0 / p.Ins(t, p.S)) * (y[8]/p.N_0). "
                               # "transfo aggrédation *2. "
                               "Ajout '- (y[4] / y[8]) * abs(dydt[8])' à eqn GSK3. "
                               "n = 15. "
                               "[F_i]_0 = y0[6] = équilibre, et tau aussi. "
                               # "d_FiN = 1 / (20 * 365) (Kril et al. 2002). "
                               "d_FiN = 1 / (2.51 * 365). "
                               "[F_o]_0 = 5e-17. " 
                               # "K_Fi = 0.1 * (0.6 * (6e-3 * self.rho_cerveau)) * 1e-6 ~= 3.708e-10. "
                               "K_Fi = 1.25e-10. "
                               # "K_Ta = 2.24e-12 (pour perte neurones)."
                               "K_TaM = 2.24e-12 * 2e2, impact aussi K_TaMhat. "
                               "kappa_MhatproTa = (1.5e-9 / 18 * 24) / (4e6 * m_Mhat), min value from Fadok98 (impact aussi kappa_MproTa). "
                               # "Cond init Manti et hat{M}_anti pour que cond init Ta et I_10 soit égales."
                               "kappa_PMhat = 0.33 * 1e-3."
                               # "kappa_TaA = 0.92 / 100e-9 * 1e-2. "
                               # "kappa_TaManti = 4.8 * 1e-1. "
                               # "kappa_TaMhatanti = 1 / (10 / 24) * 1e-2. "
                               "K_P = 6.23e-10 * 1e2. "
                               "K_ABpo = (1.11 + 0.53) / 527.4 / 1000 * 1e-18. "
                               # "beta = 1e-1. "
                               "[M_anti]_0 = 1e-4, [M_pro]_0 = 1e-12. "
                               "[hat{M}_pro]_0 = [hat{M}_anti]_0 = 0. "
                               # "[AB_p^o]_0 = 2.5e-25. "
                               # "[AB_m^o]_0 et [AB_o^o]_0 à l'équilibre. "
                               "[AB_m^o]_0, [AB_o^o]_0 et [AB_p^o]_0 à l'équilibre. "
                               # "[AB_m^o]_0 = 4e-11, [AB_o^o]_0 = 6e-17 et [AB_p^o]_0 = 5e-28. "
                               # "self.K_Fo = 11 * ((1000 * 72500) / Avogadro) * 1000 * 5. "
                               # "K_ABooM = 0.060 / 527.4 / 1000 * 1e-5. "
                               # "[M_NA]_0 = (1 - 0.67/100) * M_0 - (y0[11] + y0[12]), [hat{M}_pro]_0 = [hat{M}_anti]_0 = 0.67/100 * M_0 * 0.5. Impact aussi plaques et cytokines. "
                               "lambda_ABi = 1.4157348479999998e-09 (Lindstrom21); impact aussi lambda_ABmo. "
                               # "kappa_FoM = kappa_ABooM = TotalMaxActivRateM * 1/2. "
                               # "kappa_MproTa = kappa_MhatproTa * 0.5. "
                               # "TotalMaxActivRateM = 0.2141 * 0.5. "
                               # "K_Fo = 16 * ((1000 * 72500) / Avogadro) * 1000."
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

"""Rate figures"""

fig_neuloss = ff.fig_neuronal_loss_rates(sol.t, sol.y)

fig_microactiv = ff.fig_microglia_activation_rates(sol.t, sol.y)

fig_astroactiv = ff.fig_astrocyte_activation_rates(sol.t, sol.y)

"""Save the figures of rates"""
# today = datetime.date.today()
# date = today.strftime("%y-%m-%d")
# my_path = os.path.abspath('Figures')
# # Save the neuronal loss figure
# fig_neuloss.savefig(os.path.join(my_path, FigName[:19] + "NeuronalLoss"), dpi=180)
# # Save the microglia activation figure
# fig_microactiv.savefig(os.path.join(my_path, FigName[:19] + "MicrogliaActiv"), dpi=180)
# # Save the astroglia activation figure
# fig_astroactiv.savefig(os.path.join(my_path, FigName[:19] + "AstrogliaActiv"), dpi=180)

NeuronalLossInPercent = (sol.y[8, -1] - sol.y[8, 0]) / sol.y[8, 0] * 100
print("Concentration de neurones initiale: ", sol.y[8, 0], "\nConcentration de neurones finale: ", sol.y[8, -1],
      "\nPourcentage de perte: ", NeuronalLossInPercent)

plt.show()
