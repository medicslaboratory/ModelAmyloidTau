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
import ratefiguresfunctions as rff


AgeStart = 30
AgeEnd = 80
decades = int((AgeEnd - AgeStart) / 10)

max_step = np.inf  # Default
maxstepstr = "Default"
# max_step = 1
# maxstepstr = str(max_step).replace('.', '')
# rtol = np.ones(19) * 1e-10
# rtol[3] = 1e-29
# rtolstr = "Array"
# rtol = 1e-10  # Default value : 1e-3
rtol = 1e-5
# rtol = 1e-10
rtolstr = "{:.0e}".format(rtol)

# # atol default value : 1e-6
# atol = np.ones(19) * 1e-15
# atol[3] = 1e-20
# atolstr = "array"
atol = 1e-22
# atol = 1e-26
atolstr = "{:.0e}".format(atol)

method = "BDF"
# method = "Radau"
# method = "LSODA"


# number = 1
# CommentModif = "xi=5e-1"
# CommentModif = "Inst=Ins0"
# CommentModif = "maxstep" + maxstepstr + "_rtol" + rtolstr + \
#                   "_atol" + atolstr + "_layouttight"
# CommentModif = "maxstep" + "Default" + "_rtol" + rtolstr + \
#                   "_atol" + atolstr  # + "_layouttight_modifplacelegend"

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

# p = param.Parameters(Sex=0, APOE_status=0)
# y0 = InitialConditions(p, AgeStart)
# sol = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0, method=method, max_step=max_step, args=[0, 0],
#                 rtol=rtol, atol=atol)
#
# solt = sol.t / 365
#
# FigName = ff.main_figure(solt, sol.y, p, y0, AgeStart, AgeEnd, method, maxstepstr, rtolstr, atolstr, number,
#                          CommentModif="", IntraneuronalConcentration=True)

"""main figure, 4 sub-models"""
number = 1
CommentModif = "vf"
# CommentModif = "xi=5e-1"
# CommentModif = "Inst=Ins0"

sols = ff.run_4_models(AgeStart, AgeEnd, method, max_step, rtol, atol, InsVar=True, xi=0.5)

FigName = ff.main_figure_4_models(sols, AgeStart, AgeEnd, method, max_step, rtol, atol, number,
                                  CommentModif=CommentModif, SkipFirstHalfYear=False, orientation="paysage")

# ff.figure_intracellular_concentrations(sols, AgeStart, AgeEnd, method, maxstepstr, rtolstr, atolstr, number,
#                                        CommentModif=Commentmodif)

"""Add information to the figure."""
# FigInfos = {"max_step": str(max_step),
#             "Début": "Début intégration",  # "Ignore la première demi-année."
#             "Modification(s)": "d_TaN = 7.26e-3 / 365 * 10. "
#                                "p.lambda_InsG * (p.Ins_0 / p.Ins(t, p.S)) * (y[8]/p.N_0). "
#                                # "transfo aggrédation *2. "
#                                "Ajout '- (y[4] / y[8]) * abs(dydt[8])' à eqn GSK3. "
#                                "n = 15. "
#                                "[F_i]_0 = y0[6] = équilibre, et tau aussi. "
#                                # "d_FiN = 1 / (20 * 365) (Kril et al. 2002). "
#                                "d_FiN = 1 / (2.51 * 365). "
#                                "[F_o]_0 = 5e-17. "
#                                # "K_Fi = 0.1 * (0.6 * (6e-3 * self.rho_cerveau)). "  # * 1e-6 ~= 3.708e-10
#                                "K_Fi = 1.25e-10. "
#                                # "K_Ta = 2.24e-12 (pour perte neurones)."
#                                "K_TaM = 2.24e-12 * 2e2, impact aussi K_TaMhat. "
#                                "kappa_MhatproTa = (1.5e-9 / 18 * 24) / (4e6 * m_Mhat), min value from Fadok98 (impact aussi kappa_MproTa). "
#                                # "Cond init Manti et hat{M}_anti pour que cond init Ta et I_10 soit égales."
#                                "kappa_PMhat = 1 / 3 * 1e-2."
#                                # "kappa_TaA = 0.92 / 100e-9 * 1e-2. "
#                                # "kappa_TaManti = 4.8 * 1e-1. "
#                                # "kappa_TaMhatanti = 1 / (10 / 24) * 1e-2. "
#                                "K_P = 6.23e-10 * 1e2. "
#                                # "K_ABpo = (1.11 + 0.53) / 527.4 / 1000 * 1e-18. "
#                                # "beta = 1e-1. "
#                                "[M_anti]_0 = 1e-4, [M_pro]_0 = 1e-12. "
#                                "[hat{M}_pro]_0 = [hat{M}_anti]_0 = 0. "
#                                # "[AB_p^o]_0 = 2.5e-25. "
#                                # "[AB_m^o]_0 et [AB_o^o]_0 à l'équilibre. "
#                                "[AB_m^o]_0, [AB_o^o]_0 et [AB_p^o]_0 à l'équilibre. "
#                                # "[AB_m^o]_0 = 4e-11, [AB_o^o]_0 = 6e-17 et [AB_p^o]_0 = 5e-28. "
#                                # "self.K_Fo = 11 * ((1000 * 72500) / Avogadro) * 1000 * 5. "
#                                # "K_ABooM = 0.428 / 496.7 / 1000 * 1e2. "
#                                "K_ABooM = 0.060 / 527.4 / 1000 * 1.5e2. "  # * 2e2
#                                # "[M_NA]_0 = (1 - 0.67/100) * M_0 - (y0[11] + y0[12]), [hat{M}_pro]_0 = [hat{M}_anti]_0 = 0.67/100 * M_0 * 0.5. Impact aussi plaques et cytokines. "
#                                "lambda_ABi = 1.415734848e-06 / 2 (Lindstrom21); impact aussi lambda_ABmo et lambda_AABmo. "
#                                # "kappa_FoM = kappa_ABooM = TotalMaxActivRateM * 1/2. "
#                                # "kappa_MproTa = kappa_MhatproTa * 0.5. "
#                                # "TotalMaxActivRateM = 0.2141 * 0.5. "
#                                # "K_Fo = 16 * ((1000 * 72500) / Avogadro) * 1000. "
#                                # "kappa_ABmoABoo = kappa_ABmoABoo_max, pas _min."
#                                # "kappa_MhatantiTb = (kappa_MhatantiTb_max + kappa_MhatantiTb_min) / 2; impact aussi kappa_MantiTb. "
#             }
#
# im = Image.open("Figures/" + FigName)
# Infos = PngImagePlugin.PngInfo()
# for x in im.info:
#     Infos.add_text(x, str(im.info[x]))
# for x in FigInfos:
#     Infos.add_text(x, FigInfos[x])
# im.save("Figures/" + FigName, "png", pnginfo=Infos)

# im = Image.open("Figures/" + FigName)
# print(im.info)

# # To access the infos by typing the following in the console :
# from PIL import Image
# im = Image.open("Path of the figure")  # "Path of the figure" is, for example, "Figures/Figure_22-09-08_..._01.png".
# im.info


"""Figure of comparison of models (base vs Without insulin variation)"""
# sol_Fpos_InsVarT = ff.run_1_model(0, 1, AgeStart, AgeEnd, method, max_step, rtol, atol, InsVar=True)
# sol_Fpos_InsVarF = ff.run_1_model(0, 1, AgeStart, AgeEnd, method, max_step, rtol, atol, InsVar=False)
#
# number = 1
# CommentModif = "F_APOE+_InsVar_portrait"
# # # CommentModif = "F_APOE+_InsVar_paysage"
# #
# # FigName = ff.main_figure_comparison(sol_Fpos_InsVarT.t, sol_Fpos_InsVarT.y, sol_Fpos_InsVarF.t, sol_Fpos_InsVarF.y,
# #                                     ["$Ins(t)$ original", "$Ins(t) = Ins_0$"], AgeStart, AgeEnd, method, number,
# #                                     CommentModif=CommentModif, color1="k", color2="darkorange")  #, orientation="paysage")
# FigName = ff.main_figure_comparison_many([sol_Fpos_InsVarT.t, sol_Fpos_InsVarF.t], [sol_Fpos_InsVarT.y, sol_Fpos_InsVarF.y],
#                                          ["$Ins(t)$ original", "$Ins(t) = Ins_0$"], ["k", "darkorange"], AgeStart,
#                                          AgeEnd, method, number, CommentModif=CommentModif)  #, orientation="paysage")

"""Figure of comparison of models (base vs 0<xi<1)"""
# xis = [1, 0.8, 0.5, 0.3]
# solts = []
# solys = []
# for i in range(len(xis)):
#     sol = ff.run_1_model(0, 1, AgeStart, AgeEnd, method, max_step, rtol, atol, InsVar=True, xi=xis[i])  # F,+
#     # sol = ff.run_1_model(1, 0, AgeStart, AgeEnd, method, max_step, rtol, atol, InsVar=True, xi=xis[i])  # H,-
#     solts.append(sol.t)
#     solys.append(sol.y)
#
# number = 1
# xi_text = "{:.0e}".format(xis[0])
# labels = [r"$\xi = 1$ (base)"]
# for xi in xis[1:]:
#     xi_text = xi_text + "_{:.0e}".format(xi)
#     labels.append(r"$\xi = {:.1f}$".format(xi))
#
# # CommentModif = "H_APOE-_xi=" + xi_text + "_portrait"
# CommentModif = "F_APOE+_xi=" + xi_text + "_portrait"
# # CommentModif = "F_APOE+_xi=" + xi_text + "_paysage"
#
# colors_list = ["k", "b", "g", "orange", "purple"]
# alphas = [1, 0.8, 0.7, 0.6, 0.5]
#
# # FigName = ff.main_figure_comparison(sol_Fpos_InsVarT.t, sol_Fpos_InsVarT.y, sol_Fpos_InsVarF.t, sol_Fpos_InsVarF.y,
# #                                     [r"$\xi = 1$ (base)", r"$\xi = 0.5$ (modifié)"], AgeStart, AgeEnd, method, number,
# #                                     CommentModif=CommentModif, color1="k", color2="darkorange")  #, orientation="paysage")
# FigName = ff.main_figure_comparison_many(solts, solys, labels, colors_list[:len(xis)], AgeStart, AgeEnd, method, number,
#                                          CommentModif=CommentModif, alphas=alphas)  # , orientation="paysage")


"""Figure of a single variable."""
# FigNameABpo = ff.fig_one_variable_all_models(3, sols, AgeStart, AgeEnd, method, number,
#                                              CommentModif=CommentModif, SkipFirstHalfYear=False)

# FigNameNeu = ff.fig_one_variable_all_models(8, sols, AgeStart, AgeEnd, method, number,
#                                             CommentModif=CommentModif, SkipFirstHalfYear=False)


"""Rate figures"""
# today = datetime.date.today()
# date = today.strftime("%y-%m-%d")
# my_path = os.path.abspath('Figures')
#
# for combi in [[0, 0], [0, 1], [1, 0], [1, 1]]:
#     sex, APOE4 = combi
#     p = param.Parameters(Sex=sex, APOE_status=APOE4)
#     y0 = InitialConditions(p, AgeStart)
#     sol = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0, method=method, max_step=max_step,
#                     args=[sex, APOE4], rtol=rtol, atol=atol)
#
#     # fig_neuloss_V2 = rff.fig_neuronal_loss_rates_V2(sol.t, sol.y, p)
#     fig_astroactiv = rff.fig_astrocyte_activation_rates(sol.t, sol.y, p)
#
#     if p.S == 0:
#         sex = "F"
#     else:  # p.S == 1:
#         sex = "M"
#     if p.AP == 1:
#         APOE = "+"
#     else:  # p.AP == 0:
#         APOE = "-"
#
#     # fig_neuloss_V2.savefig(os.path.join(my_path, "Figure_" + date + "_" + f"{number:02}" + "_NeuronalLossV2_" + sex +
#     #                                     "_APOE" + APOE + "_" + method + "_" + str(AgeEnd - AgeStart).replace(".", "") +
#     #                                     "y_" + CommentModif + ".png"), dpi=180)
#     fig_astroactiv.savefig(os.path.join(my_path, "Figure_" + date + "_" + f"{number:02}" + "_AstrogliaActiv_" + sex +
#                                         "_APOE" + APOE + "_" + method + "_" + str(AgeEnd - AgeStart).replace(".", "") +
#                                         "y_" + CommentModif + ".png"), dpi=180)


# sex = 1
# APOE4 = 0
# p = param.Parameters(Sex=sex, APOE_status=APOE4)
# y0 = InitialConditions(p, AgeStart)
# sol = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0, method=method, max_step=max_step, args=[sex, APOE4],
#                 rtol=rtol, atol=atol)

# fig_neuloss_V2 = rff.fig_neuronal_loss_rates_V2(sol.t, sol.y, p)

# fig_neuloss = rff.fig_neuronal_loss_rates(sol.t, sol.y, p)
#
# fig_microactiv = rff.fig_microglia_activation_rates(sol.t, sol.y, p)
# fig_microactiv_V2 = rff.fig_microglia_activation_rates_V2(sol.t, sol.y, p)
#
# fig_astroactiv = rff.fig_astrocyte_activation_rates(sol.t, sol.y, p)


"""Save the figures of rates"""
# today = datetime.date.today()
# date = today.strftime("%y-%m-%d")
# my_path = os.path.abspath('Figures')
#
# if p.S == 0:
#     sex = "F"
# else:  # p.S == 1:
#     sex = "M"
# if p.AP == 1:
#     APOE = "+"
# else:  # p.AP == 0:
#     APOE = "-"
#
# FigId = "Figure_" + date + "_" + f"{number:02}" + "_" + method + "_APOE" + APOE + "_" + sex + "_" + \
#               str(AgeEnd - AgeStart).replace(".", "") + "y"

# # Save the neuronal loss figure
# fig_neuloss_V2.savefig(os.path.join(my_path, FigId + "_NeuronalLossV2_" + CommentModif + ".png"), dpi=180)
# fig_neuloss.savefig(os.path.join(my_path, FigId + "_NeuronalLoss_" + CommentModif + ".png"), dpi=180)
# fig_neulossV2.savefig(os.path.join(my_path, "Figure_" + date + "_" + f"{number:02}" + "_NeuronalLossV2_" + sex +
#                                  "_APOE" + APOE + "_" + method + "_" + str(AgeEnd - AgeStart).replace(".", "") +
#                                  "y_" + CommentModif + ".png"), dpi=180)
# # Save the microglia activation figure
# fig_microactiv.savefig(os.path.join(my_path, FigId + "_MicrogliaActiv_" + CommentModif + ".png"), dpi=180)
# fig_microactiv_V2.savefig(os.path.join(my_path, "Figure_" + date + "_" + f"{number:02}" + "_MicrogliaActiv_" + sex +
#                                        "_APOE" + APOE + "_" + method + "_" + str(AgeEnd - AgeStart).replace(".", "") +
#                                        "y_" + CommentModif + ".png"), dpi=180)
# # Save the astroglia activation figure
# fig_astroactiv.savefig(os.path.join(my_path, FigId + "_AstrogliaActiv_" + CommentModif + ".png"), dpi=180)


"""Print neuronal loss in percent"""
# NeuronalLossInPercent = (sol.y[8, -1] - sol.y[8, 0]) / sol.y[8, 0] * 100
# print("Concentration de neurones initiale: ", sol.y[8, 0], "\nConcentration de neurones finale: ", sol.y[8, -1],
#       "\nPourcentage de perte: ", NeuronalLossInPercent)

"""Microglia activation"""
# s = 0
# labels = ["F, APOE4-", "F, APOE4+", "M, APOE4-", "M, APOE4+"]
# for sol in sols:  # [sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos]:
#     MicrogliaActivInPercent = (sol.y[10, -1] - sol.y[10, 0]) / sol.y[10, 0] * 100
#     # print(labels[s])
#     # print("Concentration de neurones initiale: ", sol.y[8, 0], "\nConcentration de neurones finale: ", sol.y[8, -1],
#     #       "\nPourcentage de perte: ", NeuronalLossInPercent)
#     print(labels[s] + ": " + str(MicrogliaActivInPercent) + " %")
#     s = s + 1
#
plt.show()
