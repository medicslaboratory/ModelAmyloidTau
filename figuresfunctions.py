# Date : 7 Octobre 2022
# Author : Éléonore Chamberland

import os
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from scipy.integrate import solve_ivp
import parameters as param
from InitialConditions import InitialConditions
import equations as eqns
import datetime

# p = param.Parameters()


def run_1_model(Sex, APOE_status, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=True, xi=1):
    p = param.Parameters(Sex=Sex, APOE_status=APOE_status, xi=xi)
    y0 = InitialConditions(p, AgeStart)
    sol = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0, method=methodstr,
                    args=[Sex, APOE_status, InsVar, xi], max_step=max_step, rtol=rtol, atol=atol)
    return sol


def run_4_models(AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=True, xi=1):
    sol_Fneg = run_1_model(0, 0, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=InsVar, xi=xi)

    sol_Fpos = run_1_model(0, 1, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=InsVar, xi=xi)

    sol_Mneg = run_1_model(1, 0, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=InsVar, xi=xi)

    sol_Mpos = run_1_model(1, 1, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=InsVar, xi=xi)

    # p_Fneg = param.Parameters(Sex=0, APOE_status=0)  # F ; APOE4-
    # y0_Fneg = InitialConditions(p_Fneg, AgeStart)
    # sol_Fneg = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Fneg, method=methodstr, args=[0, 0, InsVar],
    #                      max_step=max_step, rtol=rtol, atol=atol)
    #
    # p_Fpos = param.Parameters(Sex=0, APOE_status=1)  # F ; APOE4+
    # y0_Fpos = InitialConditions(p_Fpos, AgeStart)
    # sol_Fpos = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Fpos, method=methodstr, args=[0, 1, InsVar],
    #                      max_step=max_step, rtol=rtol, atol=atol)
    #
    # p_Mneg = param.Parameters(Sex=1, APOE_status=0)  # M ; APOE4-
    # y0_Mneg = InitialConditions(p_Mneg, AgeStart)
    # sol_Mneg = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Mneg, method=methodstr, args=[1, 0, InsVar],
    #                      max_step=max_step, rtol=rtol, atol=atol)
    #
    # p_Mpos = param.Parameters(Sex=1, APOE_status=1)  # M ; APOE+
    # y0_Mpos = InitialConditions(p_Mpos, AgeStart)
    # sol_Mpos = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Mpos, method=methodstr, args=[1, 1, InsVar],
    #                      max_step=max_step, rtol=rtol, atol=atol)

    return sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos


def main_figure(solt, soly, p, y0, AgeStart, AgeEnd, methodstr, maxstepstr, rtolstr, atolstr, number, CommentModif="",
                IntraneuronalConcentration=False, SkipFirstHalfYear=False, color="b"):
    """
    Figure principale du modèle présentant les graphiques de chaque variable en fonction du temps.
    La figure est générée puis enregistrée dans le dossier "Figure".

    :param color:
    :param solt: Array des temps (en années). (ndarray, shape (n_points,) ; Time points).
    :param soly: Array des arrays des valeurs de chaques paramètres.
                (ndarray, shape (n, n_points) ; Values of the solution at t).
    :param y0: Array ou liste des conditions initiales.
    :param AgeStart: int. Âge du début de la simulation (en années).
    :param AgeEnd: int. Âge de fin de la simulation (en années).
    :param methodstr: str. Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param maxstepstr: str. Pas de temps maximal permis pour l'intégration du modèle. Sera ajouté au nom
                d'enregistrement de la figure.
    :param rtolstr: str. "rtol" utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param atolstr: str. "atol" utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param number: int. Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1.
    :param CommentModif: str, optionnel. Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure.
    :param IntraneuronalConcentration: bool, optionnel. Si l'on affiche les graphiques des concentrations
                intraneuronales pour ABi, GSK3, tau et F_i sur leur graphique respectif. Valeur par défaut est "False".
    :param SkipFirstHalfYear: bool, optionnel. Si l'on ignore la première demie année pour l'affichage les valeurs.
                Valeur par défaut est "False".
    :return: FigName : Nom de la figure utilisé pour l'enregistrement.
    """
    fig, axs = plt.subplots(nrows=4, ncols=5, sharex="all")  # If no text under the figure, add: layout="constrained"
    fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)

    # Making a list for Label names in the plot
    labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$',
                 r'$\tau$', '$F_i$', '$F_o$', '$N$', '$A$', '$M_{NA}$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$',
                 r'$\hat{M}_{anti}$', r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

    k = np.argmax(solt >= (AgeStart+0.5))  # First indice of "solt" to skip the first half year.

    i = 0
    for ax in axs.flat:
        if i < 19:
            ax.plot(solt, soly[i, :], color=color)  # , '.-', ms=2

            if SkipFirstHalfYear:  # """Plot sans la première demie année."""
                ax.plot(solt[k:], soly[i, k:], color=color)

            ax.grid()
            if i >= 14:
                # ax.set_xlabel('Age (years)')
                ax.set_xlabel('Âge (années)')
            ax.set_ylabel(labelname[i])
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
            # notation will be used if exp <= -1 or exp >= 1.
            ax.yaxis.set_major_formatter(formatter)
        i = i+1

    """Pour afficher certains les graphs de concentration intraneuronale."""
    if IntraneuronalConcentration:
        """Add the plot of ABi/N in the graph for ABi"""
        axs[0, 0].set_ylabel(labelname[0], color=color)
        ax1 = axs[0, 0].twinx()
        ax1.plot(solt, soly[0, :]/soly[8, :], "g-")  # ABi/N
        ax1.set_ylabel(r"$A \beta^{i} / N$", color='g')
        formatter1 = ticker.ScalarFormatter(useMathText=True)
        formatter1.set_scientific(True)
        formatter1.set_powerlimits((-1, 1))
        ax1.yaxis.set_major_formatter(formatter1)

        """Add the plot of GSK3/N in the graph for GSK3"""
        axs[0, 4].set_ylabel(labelname[4], color=color)
        ax2 = axs[0, 4].twinx()
        ax2.plot(solt, soly[4, :] / soly[8, :], "g-")  # G/N
        ax2.set_ylabel(r"$GSK3/N$", color='g')
        formatter2 = ticker.ScalarFormatter(useMathText=True)
        formatter2.set_scientific(True)
        formatter2.set_powerlimits((-1, 1))
        ax2.yaxis.set_major_formatter(formatter2)

        """Add the plot of tau/N in the graph for tau"""
        axs[1, 0].set_ylabel(labelname[5], color=color)
        ax3 = axs[1, 0].twinx()
        ax3.plot(solt, soly[5, :] / soly[8, :], "g-")  # tau/N
        ax3.set_ylabel(r"$\tau/N$", color='g')
        formatter3 = ticker.ScalarFormatter(useMathText=True)
        formatter3.set_scientific(True)
        formatter3.set_powerlimits((-1, 1))
        ax3.yaxis.set_major_formatter(formatter3)

        """Add the plot of F_i/N in the graph for F_i"""
        axs[1, 1].set_ylabel(labelname[6], color=color)
        ax4 = axs[1, 1].twinx()
        ax4.plot(solt, soly[6, :] / soly[8, :], "g-")  # tau/N
        ax4.set_ylabel(r"$F_i/N$", color='g')
        formatter4 = ticker.ScalarFormatter(useMathText=True)
        formatter4.set_scientific(True)
        formatter4.set_powerlimits((-1, 1))
        ax4.yaxis.set_major_formatter(formatter4)

    """Pour afficher un graph à la dernière position 
    (doit retirer le 'axs[3, 4].remove()' et 'axs.flat[14].tick_params('x', labelbottom=True)' ci-dessous)."""
    axs.flat[14].tick_params('x', labelbottom=True)
    axs[3, 4].remove()
    # ax = axs[3, 4]
    # # # dABpodt = p.kappa_ABooABpo * (soly[2, :] ** 2) - ((p.d_MantiABpo * soly[12, :] + p.d_hatMantiABpo
    # # #                                                    * soly[14, :]) * (1 + p.AP * p.delta_APdp) *
    # # #                                                   (soly[3, :] / (soly[3, :] + p.K_ABpo)))
    # # # ax.plot(solt, dABpodt)  # , '.-', ms=2
    # # # dNdtFi = - p.d_FiN * (1 / (1 + np.exp(- p.n * (soly[6, :] - p.K_Fi)))) * soly[8, :]
    # # # ax.plot(solt, dNdtFi)  # , '.-', ms=2
    # ax.grid()
    # # # # ax.set_xlabel('Age (years)')
    # # # ax.set_xlabel('Âge (années)')
    # # # ax.set_ylabel("dABpo/dt")
    # # # ax.set_ylabel("dN/dt du à F_i")
    # ax2 = axs[3, 4].twinx()
    # # # axs[3, 4].plot(solt, soly[12, :], "b-")  # M_anti
    # # # axs[3, 4].set_ylabel("Manti", color='b')
    # # # ax2.plot(solt, soly[16, :], "g-")  # I_10
    # # # ax2.set_ylabel("I10", color='g')
    # axs[3, 4].plot(solt, soly[11, :], "b-")  # M_pro
    # # axs[3, 4].set_ylabel("M_pro", color='b')
    # axs[3, 4].plot(solt, soly[13, :], "r-")  # hat{M}_pro
    # axs[3, 4].set_ylabel("M_pro (blue); hat{M}_pro (red)")
    # ax2.plot(solt, soly[17, :], "g-")  # T_alpha
    # ax2.set_ylabel("T_a", color='g')

    plt.subplots_adjust(hspace=.2, wspace=.2)
    plt.tight_layout()

    """Peut-être utiliser cela plutôt que les deux lignes précédentes si overlap. Ne l'utilise pas de base, car
    produit des graphs plus petits (plus distancés)."""
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

    today = datetime.date.today()
    date = today.strftime("%y-%m-%d")

    if CommentModif != "":
        CommentModif = "_" + CommentModif

    my_path = os.path.abspath('Figures')
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_APOE" + APOE + "_" + sex + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + \
              "_atol" + atolstr + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number+1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_APOE" + APOE + "_" + sex + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + \
                  "_atol" + atolstr + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=180)

    return FigName


def main_figure_4_models(sols, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, number, CommentModif="",
                         SkipFirstHalfYear=False, orientation="portrait"):
    """
    Figure principale du modèle présentant les graphiques de chaque variable en fonction du temps, pour les quatre
    possibilités: APOE4 +/- et sexe (F/M)
    La figure est générée puis enregistrée dans le dossier "Figure".

    :param sols: tuple. Solution de l'intégration du modèle pour les quatres possibilités, dans l'ordre suivant :
                sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos
    :param AgeStart: int. Âge du début de la simulation (en années).
    :param AgeEnd: int. Âge de fin de la simulation (en années).
    :param methodstr: str. Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param max_step: str. Pas de temps maximal permis pour l'intégration du modèle. Sera ajouté au nom
                d'enregistrement de la figure.
    :param rtol: str. "rtol" utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param atol: str. "atol" utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param number: int. Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1.
    :param CommentModif: str, optionnel. Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure.
    :param SkipFirstHalfYear: bool, optionnel. Si l'on ignore la première demie année pour l'affichage les valeurs.
                Valeur par défaut est "False".
    :return: FigName : Nom de la figure utilisé pour l'enregistrement.
    """
    sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos = sols

    p_Fneg = param.Parameters(Sex=0, APOE_status=0)  # F ; APOE4-
    y0_Fneg = InitialConditions(p_Fneg, AgeStart)
    # sol_Fneg = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Fneg, method=methodstr, args=[0, 0],
    #                      max_step=max_step, rtol=rtol, atol=atol)
    #
    p_Fpos = param.Parameters(Sex=0, APOE_status=1)  # F ; APOE4+
    y0_Fpos = InitialConditions(p_Fpos, AgeStart)
    # sol_Fpos = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Fpos, method=methodstr, args=[0, 1],
    #                      max_step=max_step, rtol=rtol, atol=atol)
    #
    p_Mneg = param.Parameters(Sex=1, APOE_status=0)  # M ; APOE4-
    y0_Mneg = InitialConditions(p_Mneg, AgeStart)
    # sol_Mneg = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Mneg, method=methodstr, args=[1, 0],
    #                      max_step=max_step, rtol=rtol, atol=atol)
    #
    p_Mpos = param.Parameters(Sex=1, APOE_status=1)  # M ; APOE4+
    y0_Mpos = InitialConditions(p_Mpos, AgeStart)
    # sol_Mpos = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Mpos, method=methodstr, args=[1, 1],
    #                      max_step=max_step, rtol=rtol, atol=atol)

    if orientation == "portrait":
        fig, axs = plt.subplots(nrows=5, ncols=4, sharex="all", figsize=(9, 10), layout="tight")  # "constrained"
    elif orientation == "paysage":
        fig, axs = plt.subplots(nrows=4, ncols=5, sharex="all", layout="tight")
        fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)
    else:
        print("Paramètre <orientation> non valide.")
        return

    # Making a list for Label names in the plot
    labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$',
                 r'$\tau$', '$F_i$', '$F_o$', '$N$', '$A$', '$M_{NA}$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$',
                 r'$\hat{M}_{anti}$', r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

    colors = ["hotpink", "purple", "royalblue", "darkblue"]
    labels = ["F, APOE4-", "F, APOE4+", "M, APOE4-", "M, APOE4+"]
    s = 0
    for sol in [sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos]:
        solt = sol.t / 365
        soly = sol.y
        k = np.argmax(solt >= (AgeStart+0.5))  # First indice of "solt" to skip the first half year.
        i = 0
        for ax in axs.flat:
            if i < 19:
                if i < 18:
                    ax.plot(solt, soly[i, :], color=colors[s])  # , '.-', ms=2

                    if SkipFirstHalfYear:  # """Plot sans la première demie année."""
                        ax.plot(solt[k:], soly[i, k:], color=colors[s])

                else:  # i == 18
                    ax.plot(solt, soly[i, :], color=colors[s], label=labels[s])  # , '.-', ms=2

                    if SkipFirstHalfYear:  # """Plot sans la première demie année."""
                        ax.plot(solt[k:], soly[i, k:], color=colors[s], label=labels[s])

                if s == 3:  # Si l'on est à la dernière option.
                    ax.grid()
                    if orientation == "portrait":
                        indice = 15
                    else:  # orientation == "paysage"
                        indice = 14
                    if i >= indice:
                        # ax.set_xlabel('Age (years)')
                        ax.set_xlabel('Âge (années)')
                    ax.set_ylabel(labelname[i])
                    formatter = ticker.ScalarFormatter(useMathText=True)
                    formatter.set_scientific(True)
                    formatter.set_powerlimits((-1, 1))
                    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
                    # notation will be used if exp <= -1 or exp >= 1.
                    ax.yaxis.set_major_formatter(formatter)
            i = i + 1
        s = s + 1
    handles, labels = axs.flat[18].get_legend_handles_labels()

    if orientation == "portrait":
        fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.81, 0.06))  # (0.82, 0.06)
        axs.flat[15].tick_params('x', labelbottom=True)
    else:  # orientation == "paysage"
        fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.85, 0.08))
        axs.flat[14].tick_params('x', labelbottom=True)  # Ajout des graduations au 14e graphique (hat{M}_anti)

    axs.flat[19].remove()  # Retire le dernier graphique

    # plt.subplots_adjust(hspace=.2, wspace=.2)
    # plt.tight_layout()

    """Save the plot as a .png file"""

    today = datetime.date.today()
    date = today.strftime("%y-%m-%d")

    if CommentModif != "":
        CommentModif = "_" + CommentModif

    my_path = os.path.abspath('Figures')
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_AllInOne_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number+1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_AllInOne_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=180)

    s = 0
    for sol in [sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos]:
        NeuronalLossInPercent = (sol.y[8, -1] - sol.y[8, 0]) / sol.y[8, 0] * 100
        # print(labels[s])
        # print("Concentration de neurones initiale: ", sol.y[8, 0], "\nConcentration de neurones finale: ", sol.y[8, -1],
        #       "\nPourcentage de perte: ", NeuronalLossInPercent)
        print(labels[s] + ": " + str(NeuronalLossInPercent) + " %")
        s = s + 1

    s = 0

    for y0 in [y0_Fneg, y0_Fpos, y0_Mneg, y0_Mpos]:
        icNameValue = [str(labelname[i]) + " = " + "{:.5e}".format(y0[i]) for i in np.arange(19)]
        print(labels[s])
        initcond = "Initial conditions used (in g/mL): " + ", ".join(icNameValue[:])
        print(initcond)
        s = s + 1

    return FigName


def fig_one_variable_all_models(variable, sols, AgeStart, AgeEnd, methodstr, number, CommentModif="",
                                SkipFirstHalfYear=False):

    sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos = sols

    # Making a list for Label names in the plot
    labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$',
                 r'$\tau$', '$F_i$', '$F_o$', '$N$', '$A$', '$M_{NA}$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$',
                 r'$\hat{M}_{anti}$', r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

    colors = ["hotpink", "purple", "royalblue", "darkblue"]
    labels = ["F, APOE4-", "F, APOE4+", "M, APOE4-", "M, APOE4+"]

    if variable == 3:  # Si ABpo
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex="all", figsize=(6.4, 7))

        if SkipFirstHalfYear:  # """Plot sans la première demie année."""
            k = np.argmax(sol_Fpos.t / 365 >= (AgeStart + 0.5))  # First indice of "solt" to skip the first half year.
            ax1.plot(sol_Fpos.t[k:] / 365, sol_Fpos.y[variable, k:], color=colors[1], label=labels[1])
            k = np.argmax(sol_Mpos.t / 365 >= (AgeStart + 0.5))  # First indice of "solt" to skip the first half year.
            ax1.plot(sol_Mpos.t[k:] / 365, sol_Mpos.y[variable, k:], color=colors[3], label=labels[3])
        else:
            ax1.plot(sol_Fpos.t / 365, sol_Fpos.y[variable, :], color=colors[1], label=labels[1])
            ax1.plot(sol_Mpos.t / 365, sol_Mpos.y[variable, :], color=colors[3], label=labels[3])

        # ax.set_xlabel('Age (years)')
        # ax1.set_xlabel('Âge (années)')
        ax1.legend()
        # ax1.set_ylabel(labelname[variable] + ", APOE4+")
        ax1.set_ylabel(labelname[variable])
        ax1.grid()
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
        # notation will be used if exp <= -1 or exp >= 1.
        ax1.yaxis.set_major_formatter(formatter)

        if SkipFirstHalfYear:  # """Plot sans la première demie année."""
            k = np.argmax(sol_Fneg.t / 365 >= (AgeStart + 0.5))  # First indice of "solt" to skip the first half year.
            ax2.plot(sol_Fneg.t[k:] / 365, sol_Fneg.y[variable, k:], color=colors[0], label=labels[0])
            k = np.argmax(sol_Mneg.t / 365 >= (AgeStart + 0.5))  # First indice of "solt" to skip the first half year.
            ax2.plot(sol_Mneg.t[k:] / 365, sol_Mneg.y[variable, k:], color=colors[2], label=labels[2])
        else:
            ax2.plot(sol_Fneg.t / 365, sol_Fneg.y[variable, :], color=colors[0], label=labels[0])
            ax2.plot(sol_Mneg.t / 365, sol_Mneg.y[variable, :], color=colors[2], label=labels[2])

        # ax2.set_ylabel(labelname[variable] + ", APOE4-")
        ax2.set_ylabel(labelname[variable])
        ax2.set_xlabel("Âge (années)")
        # ax2.set_xlabel("Age (years)")
        ax2.legend()
        ax2.grid()
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
        # notation will be used if exp <= -1 or exp >= 1.
        ax2.yaxis.set_major_formatter(formatter)

    else:
        fig, ax = plt.subplots(nrows=1, ncols=1, sharex="all")  # If no text under the figure, add: layout="constrained"
        # fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)
        s = 0
        for sol in [sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos]:
            solt = sol.t / 365
            soly = sol.y
            k = np.argmax(solt >= (AgeStart+0.5))  # First indice of "solt" to skip the first half year.

            ax.plot(solt, soly[variable, :], color=colors[s], label=labels[s])  # , '.-', ms=2

            if SkipFirstHalfYear:  # """Plot sans la première demie année."""
                ax.plot(solt[k:], soly[variable, k:], color=colors[s], label=labels[s])

            ax.set_ylabel(labelname[variable])
            # ax.set_xlabel('Age (years)')
            ax.set_xlabel('Âge (années)')
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
            # notation will be used if exp <= -1 or exp >= 1.
            ax.yaxis.set_major_formatter(formatter)
            s = s + 1

        ax.legend()
        ax.grid()
        # handles, labels = ax.get_legend_handles_labels()
        # fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.85, 0.08))

    plt.tight_layout()

    """Save the plot as a .png file"""

    today = datetime.date.today()
    date = today.strftime("%y-%m-%d")

    VariablesName = ["ABi", "ABmo", "ABoo", "ABpo", "GSK3", "tau", "NFTi", "NFTo", "Neurones", "Astrocytes",
                     "MicrogliesNA", "MicrogliesPro", "MicrogliesAnti", "MacrophagesPro", "MacrophagesAnti", "TGFb",
                     "IL10", "TNFa", "MCP1"]

    if CommentModif != "":
        CommentModif = "_" + CommentModif

    my_path = os.path.abspath('Figures')
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_" + VariablesName[variable] + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number + 1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_" + VariablesName[variable] + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=180)

    return FigName


def main_figure_comparison(solt1, soly1, solt2, soly2, labels, AgeStart, AgeEnd, methodstr, number, CommentModif="", SkipFirstHalfYear=False, color1="k", color2="b", orientation="portrait"):
    """
    Figure principale du modèle présentant les graphiques de chaque variable en fonction du temps.
    La figure est générée puis enregistrée dans le dossier "Figure".

    :param solt1: Array des temps (en jours) de la première série. (ndarray, shape (n_points,) ; Time points).
    :param soly1: Array des arrays des valeurs de chaque variable pour la première série.
                (ndarray, shape (n, n_points) ; Values of the solution at t).
    :param solt2: Array des temps (en années) de la deuxième série. (ndarray, shape (n_points,) ; Time points).
    :param soly2: Array des arrays des valeurs de chaque variable pour la deuxième série.
                (ndarray, shape (n, n_points) ; Values of the solution at t).
    :param labels: Liste des labels à associer [label de la première série, label de la deuxième série].
    :param AgeStart: int. Âge du début de la simulation (en années).
    :param AgeEnd: int. Âge de fin de la simulation (en années).
    :param methodstr: str. Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param maxstepstr: str. Pas de temps maximal permis pour l'intégration du modèle. Sera ajouté au nom
                d'enregistrement de la figure.
    :param rtolstr: str. "rtol" utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param atolstr: str. "atol" utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param number: int. Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1.
    :param CommentModif: str, optionnel. Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure.
    :param SkipFirstHalfYear: bool, optionnel. Si l'on ignore la première demie année pour l'affichage les valeurs.
                Valeur par défaut est "False".
    :return: FigName : Nom de la figure utilisé pour l'enregistrement.
    """
    solt1 = solt1 / 365
    solt2 = solt2 / 365

    if orientation == "portrait":
        fig, axs = plt.subplots(nrows=5, ncols=4, sharex="all", figsize=(9, 10), layout="tight")  # "constrained"
    elif orientation == "paysage":
        fig, axs = plt.subplots(nrows=4, ncols=5, sharex="all", layout="tight")
        fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)
    else:
        print("Paramètre <orientation> non valide.")
        return

    # Making a list for Label names in the plot
    labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$',
                 r'$\tau$', '$F_i$', '$F_o$', '$N$', '$A$', '$M_{NA}$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$',
                 r'$\hat{M}_{anti}$', r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

    k1 = np.argmax(solt1 >= (AgeStart+0.5))  # First indice of "solt1" to skip the first half year.
    k2 = np.argmax(solt2 >= (AgeStart + 0.5))  # First indice of "solt2" to skip the first half year.

    i = 0
    for ax in axs.flat:
        if i < 19:
            if SkipFirstHalfYear:  # """Plot sans la première demie année."""
                ax.plot(solt1[k1:], soly1[i, k1:], color=color1, label=labels[0])
                ax.plot(solt2[k2:], soly2[i, k2:], color=color2, label=labels[1])
            else:
                ax.plot(solt1, soly1[i, :], color=color1, label=labels[0])  # , '.-', ms=2
                ax.plot(solt2, soly2[i, :], color=color2, alpha=0.8, label=labels[1])


            ax.grid()
            if orientation == "portrait":
                indice = 15
            else:  # orientation == "paysage"
                indice = 14
            if i >= indice:
                # ax.set_xlabel('Age (years)')
                ax.set_xlabel('Âge (années)')
            ax.set_ylabel(labelname[i])
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
            # notation will be used if exp <= -1 or exp >= 1.
            ax.yaxis.set_major_formatter(formatter)
        i = i+1
    axs.flat[19].remove()  # Retire le dernier graphique

    handles, labels = axs.flat[18].get_legend_handles_labels()

    if orientation == "portrait":
        # fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.80, 0.07))
        fig.legend(handles, labels, loc='best', bbox_to_anchor=(0.75, 0., 0.5, 0.25))
        axs.flat[15].tick_params('x', labelbottom=True)
    else:  # orientation == "paysage"
        fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.85, 0.09))
        axs.flat[14].tick_params('x', labelbottom=True)  # Ajout des graduations au 14e graphique (hat{M}_anti)

    """Save the plot as a .png file"""
    today = datetime.date.today()
    date = today.strftime("%y-%m-%d")

    if CommentModif != "":
        CommentModif = "_" + CommentModif

    my_path = os.path.abspath('Figures')
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_Comparaison_" + methodstr + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number+1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_Comparaison_" + methodstr + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=180)

    return FigName


def main_figure_comparison_many(solst, solsy, labels, colors, AgeStart, AgeEnd, methodstr, number, CommentModif="", SkipFirstHalfYear=False, orientation="portrait", alphas=[1, 0.8]):
    """
    Figure principale du modèle présentant les graphiques de chaque variable en fonction du temps.
    La figure est générée puis enregistrée dans le dossier "Figure".

    :param solt1: Array des temps (en jours) de la première série. (ndarray, shape (n_points,) ; Time points).
    :param soly1: Array des arrays des valeurs de chaque variable pour la première série.
                (ndarray, shape (n, n_points) ; Values of the solution at t).
    :param solt2: Array des temps (en années) de la deuxième série. (ndarray, shape (n_points,) ; Time points).
    :param soly2: Array des arrays des valeurs de chaque variable pour la deuxième série.
                (ndarray, shape (n, n_points) ; Values of the solution at t).
    :param labels: Liste des labels à associer [label de la première série, label de la deuxième série].
    :param AgeStart: int. Âge du début de la simulation (en années).
    :param AgeEnd: int. Âge de fin de la simulation (en années).
    :param methodstr: str. Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param maxstepstr: str. Pas de temps maximal permis pour l'intégration du modèle. Sera ajouté au nom
                d'enregistrement de la figure.
    :param rtolstr: str. "rtol" utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param atolstr: str. "atol" utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param number: int. Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1.
    :param CommentModif: str, optionnel. Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure.
    :param SkipFirstHalfYear: bool, optionnel. Si l'on ignore la première demie année pour l'affichage les valeurs.
                Valeur par défaut est "False".
    :return: FigName : Nom de la figure utilisé pour l'enregistrement.
    """
    if len(solst) != len(solsy):
        print("Les tailles de solst et solsy ne correspondent pas.")
        return
    if len(solst) != len(labels):
        print("La taille des solutions (solst et solsy) ne correspond pas à celle de labels.")
        return
    if len(solst) != len(colors):
        print("La taille des solutions (solst et solsy) ne correspond pas à celle de colors.")
        return

    nbmodel = len(solst)

    for i in range(nbmodel):
        solst[i] = solst[i] / 365

    if orientation == "portrait":
        fig, axs = plt.subplots(nrows=5, ncols=4, sharex="all", figsize=(9, 10), layout="tight")  # "constrained"
    elif orientation == "paysage":
        fig, axs = plt.subplots(nrows=4, ncols=5, sharex="all", layout="tight")
        fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)
    else:
        print("Paramètre <orientation> non valide.")
        return

    # Making a list for Label names in the plot
    labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$',
                 r'$\tau$', '$F_i$', '$F_o$', '$N$', '$A$', '$M_{NA}$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$',
                 r'$\hat{M}_{anti}$', r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

    ks = np.zeros(nbmodel)
    for i in range(nbmodel):
        ks[i] = np.argmax(solst[i] >= (AgeStart+0.5))  # First indice of "solst[i]" to skip the first half year.

    i = 0
    for ax in axs.flat:
        if i < 19:
            if SkipFirstHalfYear:  # """Plot sans la première demie année."""
                for j in range(nbmodel):
                    ax.plot(solst[j], solsy[j][i, ks[j]:], color=colors[j], alpha=alphas[j], label=labels[j], linewidth=1.5)
            else:
                for j in range(nbmodel):
                    ax.plot(solst[j], solsy[j][i, :], color=colors[j], alpha=alphas[j], label=labels[j], linewidth=1.5)  # , '.-', ms=2

            ax.grid()
            if orientation == "portrait":
                indice = 15
            else:  # orientation == "paysage"
                indice = 14
            if i >= indice:
                # ax.set_xlabel('Age (years)')
                ax.set_xlabel('Âge (années)')
            ax.set_ylabel(labelname[i])
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
            # notation will be used if exp <= -1 or exp >= 1.
            ax.yaxis.set_major_formatter(formatter)
        i = i+1
    axs.flat[19].remove()  # Retire le dernier graphique

    handles, labels = axs.flat[18].get_legend_handles_labels()

    if orientation == "portrait":
        # fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.80, 0.07))
        fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.875, 0.11))
        axs.flat[15].tick_params('x', labelbottom=True)
    else:  # orientation == "paysage"
        # fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.85, 0.09))
        fig.legend(handles, labels, loc='center', bbox_to_anchor=(0.9, 0.125))
        axs.flat[14].tick_params('x', labelbottom=True)  # Ajout des graduations au 14e graphique (hat{M}_anti)

    """Save the plot as a .png file"""
    today = datetime.date.today()
    date = today.strftime("%y-%m-%d")

    if CommentModif != "":
        CommentModif = "_" + CommentModif

    my_path = os.path.abspath('Figures')
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_Comparaison_" + methodstr + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number+1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_Comparaison_" + methodstr + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=180)

    return FigName
