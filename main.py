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

"""main figure, 4 sub-models"""
number = 1

"""Base"""
# CommentModif = "vf"
# sols = ff.run_4_models(AgeStart, AgeEnd, method, max_step, rtol, atol, InsVar=True, xi=1)

# FigName = ff.main_figure_4_models(sols, AgeStart, AgeEnd, method, max_step, rtol, atol, number,
#                                   CommentModif=CommentModif, SkipFirstHalfYear=False)  #, orientation="paysage")
#
# ff.figure_intracellular_concentrations(sols, AgeStart, AgeEnd, method, maxstepstr, rtolstr, atolstr, number,
#                                        CommentModif=CommentModif)
# ff.figure_intracellular_concentrations_SexDiff(sols, AgeStart, AgeEnd, method, maxstepstr, rtolstr, atolstr, number,
#                                                CommentModif=CommentModif)

"""xi = 0.5"""
# CommentModif = "xi=5e-1"
# sols_xi = ff.run_4_models(AgeStart, AgeEnd, method, max_step, rtol, atol, InsVar=True, xi=0.5)
# _ = ff.main_figure_4_models(sols_xi, AgeStart, AgeEnd, method, max_step, rtol, atol, number,
#                             CommentModif=CommentModif, SkipFirstHalfYear=False)
# ff.figure_intracellular_concentrations(sols_xi, AgeStart, AgeEnd, method, maxstepstr, rtolstr, atolstr, number,
#                                        CommentModif=CommentModif)
# ff.figure_intracellular_concentrations_SexDiff(sols_xi, AgeStart, AgeEnd, method, maxstepstr, rtolstr, atolstr, number,
#                                                CommentModif=CommentModif)

"""Sans variation de l'insuline"""
# CommentModif = "Inst=Ins0"
# sols_ins = ff.run_4_models(AgeStart, AgeEnd, method, max_step, rtol, atol, InsVar=False, xi=1)
# _ = ff.main_figure_4_models(sols_ins, AgeStart, AgeEnd, method, max_step, rtol, atol, number,
#                             CommentModif=CommentModif, SkipFirstHalfYear=False)
# ff.figure_intracellular_concentrations(sols_ins, AgeStart, AgeEnd, method, maxstepstr, rtolstr, atolstr, number,
#                                        CommentModif=CommentModif)
# ff.figure_intracellular_concentrations_SexDiff(sols_ins, AgeStart, AgeEnd, method, maxstepstr, rtolstr, atolstr, number,
#                                                CommentModif=CommentModif)

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
# _ = ff.main_figure_comparison_many([sol_Fpos_InsVarT.t, sol_Fpos_InsVarF.t], [sol_Fpos_InsVarT.y, sol_Fpos_InsVarF.y],
#                                    ["$Ins(t)$ original", "$Ins(t) = Ins_0$"], ["k", "darkorange"], AgeStart,
#                                    AgeEnd, method, number, CommentModif=CommentModif)  #, orientation="paysage")

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
# CommentModif = "F_APOE+_xi=" + xi_text + "_portrait_vf"
# # CommentModif = "F_APOE+_xi=" + xi_text + "_paysage"
#
# colors_list = ["k", "b", "g", "orange", "purple"]
# alphas = [1, 0.8, 0.7, 0.6, 0.5]
#
# _ = ff.main_figure_comparison_many(solts, solys, labels, colors_list[:len(xis)], AgeStart, AgeEnd, method, number,
#                                    CommentModif=CommentModif, alphas=alphas)  # , orientation="paysage")


"""Figure of a single variable."""
# CommentModif = "vf"
#
# FigNameABpo = ff.fig_one_variable_all_models(3, sols, AgeStart, AgeEnd, method, number,
#                                              CommentModif=CommentModif, SkipFirstHalfYear=False)
#
# FigNameNeu = ff.fig_one_variable_all_models(8, sols, AgeStart, AgeEnd, method, number,
#                                             CommentModif=CommentModif, SkipFirstHalfYear=False)


"""Rate figures"""
# CommentModif = "vf"
#
# today = datetime.date.today()
# date = today.strftime("%y-%m-%d")
# my_path = os.path.abspath('Figures')
#
# InsVar = True
# # xi = 0.5
# xi = 1
# for combi in [[0, 0], [0, 1], [1, 0], [1, 1]]:
#     sex, APOE4 = combi
#     sol, p = ff.run_1_model(sex, APOE4, AgeStart, AgeEnd, method, max_step, rtol, atol,
#                             InsVar=InsVar, xi=xi, return_p=True)
#     # p = param.Parameters(Sex=sex, APOE_status=APOE4)
#     # y0 = InitialConditions(p, AgeStart)
#     # sol = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0, method=method, max_step=max_step,
#     #                 args=[sex, APOE4, InsVar], rtol=rtol, atol=atol)
#
#     fig_neuloss_V2 = rff.fig_neuronal_loss_rates_V2(sol.t, sol.y, p)
#     fig_astroactiv = rff.fig_astrocyte_activation_rates(sol.t, sol.y, p)
#     fig_microactiv_V2 = rff.fig_microglia_activation_rates_V2(sol.t, sol.y, p)
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
#     fig_neuloss_V2.savefig(os.path.join(my_path, "Figure_" + date + "_" + f"{number:02}" + "_NeuronalLossV2_" + sex +
#                                         "_APOE" + APOE + "_" + method + "_" + str(AgeEnd - AgeStart).replace(".", "") +
#                                         "y_" + CommentModif + ".png"), dpi=180)
#     fig_astroactiv.savefig(os.path.join(my_path, "Figure_" + date + "_" + f"{number:02}" + "_AstrogliaActiv_" + sex +
#                                         "_APOE" + APOE + "_" + method + "_" + str(AgeEnd - AgeStart).replace(".", "") +
#                                         "y_" + CommentModif + ".png"), dpi=180)
#     fig_microactiv_V2.savefig(os.path.join(my_path, "Figure_" + date + "_" + f"{number:02}" + "_MicrogliaActiv_" + sex +
#                                            "_APOE" + APOE + "_" + method + "_" + str(AgeEnd - AgeStart).replace(".", "")
#                                            + "y_" + CommentModif + ".png"), dpi=180)


"""Print neuronal loss in percent"""
# NeuronalLossInPercent = (sol.y[8, -1] - sol.y[8, 0]) / sol.y[8, 0] * 100
# print("Concentration de neurones initiale: ", sol.y[8, 0], "\nConcentration de neurones finale: ", sol.y[8, -1],
#       "\nPourcentage de perte: ", NeuronalLossInPercent)


"""Microglia activation"""


def print_microglia_activ(sols, model=""):
    if model != "":
        titre = "Pourcentage d'activation des microglies pour le modèle " + model
        print(titre)
    else:
        print("Pourcentage d'activation des microglies")
    s = 0
    labels = ["F, APOE4-", "F, APOE4+", "M, APOE4-", "M, APOE4+"]
    for sol in sols:  # [sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos]:
        # MicrogliaActivInPercent = abs(sol.y[10, -1] - sol.y[10, 0]) / sol.y[10, 0] * 100
        #   # abs([M_NA]_f + [M_NA]_0) / [M_NA]_0 * 100
        # # Donne pas exactement le même résultat que :
        MicrogliaActivInPercent = (sol.y[11, -1] + sol.y[12, -1]) / sol.y[10, 0] * 100
        #   # ([M_pro]_f + [M_anti]_f) / [M_NA]_0 * 100
        print(labels[s] + ": " + str(MicrogliaActivInPercent) + " %" + " ({:.2f}%)".format(MicrogliaActivInPercent))
        s = s + 1


# print_microglia_activ(sols, "base")
# print_microglia_activ(sols_ins, "Ins(t)=Ins_0")
# print_microglia_activ(sols_xi, "xi=0.5")

"""show figures"""
plt.show()
