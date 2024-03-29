"""
Author : Éléonore Chamberland
Date : 7 Octobre 2022
"""

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

# fig_dir = 'Figures'
fig_dir = 'FiguresPresentationADPD'
# fig_dir = "."
dpi = 250

xlabel = "Âge (années)"
# xlabel = "Age (years)"


def run_1_model(Sex, APOE_status, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=True, xi=1, return_p=False):
    """Runs 1 variant of the model.

    :param int Sex: Sex value: 0 for women, 1 for men.
    :param int APOE_status: APOE4 status: 0 for non bearing APOE4 allele, 1 otherwise.
    :param float AgeStart: Age to start the model at, in years.
    :param float AgeEnd: Age at which we end running the model, in years.
    :param str methodstr: Method of integration of the `scipy.integrate.solve_ivp` method. Usually ``"BDF"``.
    :param float max_step: Maximal time step of integration.
    :param float rtol: `rtol` of the `scipy.integrate.solve_ivp` method.
    :param float atol: `atol` the `scipy.integrate.solve_ivp` method
    :param InsVar: Default: ``True``. Whether we want the insulin concentration to vary or not.
    :type InsVar: bool, optional
    :param float xi: xi value of the model: 0 < xi <= 1. Will be multiplied to `TotalMaxActivRateM`, so to
        `kappa_FoM` and `kappa_ABooM`.
    :param return_p: Default: ``False``. Whether to return the `parameters.Parameters` object or not.
    :type return_p: bool, optional
    :return: sol, (p).
        sol: The resolved model. The output of the `scipy.integrate.solve_ivp` method;
        p: `parameters.Parameters` object, returned if `return_p` is set to ``True``.
    """
    p = param.Parameters(Sex=Sex, APOE_status=APOE_status, xi=xi)
    y0 = InitialConditions(p, AgeStart)
    sol = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0, method=methodstr,
                    args=[Sex, APOE_status, InsVar, xi], max_step=max_step, rtol=rtol, atol=atol)
    if return_p:
        return sol, p
    else:
        return sol


def run_4_models(AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=True, xi=1):
    """Runs the solutions for the four (4) options of the model, for the sex, and APOE4 status.

    :param float AgeStart: Age to start the model at, in years.
    :param float AgeEnd: Age at which we end running the model, in years.
    :param str methodstr: Method of integration of the `scipy.integrate.solve_ivp` method. Usually ``"BDF"``.
    :param float max_step: Maximal time step of integration.
    :param float rtol: `rtol` of the `scipy.integrate.solve_ivp` method.
    :param float atol: `atol` the `scipy.integrate.solve_ivp` method
    :param InsVar: Default: ``True``. Whether we want the insulin concentration to vary or not.
    :type InsVar: bool, optional
    :param float xi: xi value of the model: 0 < xi <= 1. Will be multiplied to `TotalMaxActivRateM`, so to
        `kappa_FoM` and `kappa_ABooM`.
    :return: sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos.
        sol_Fneg: Solution for female APOE4 negative;
        sol_Fpos: Solution for female APOE4 positive;
        sol_Mneg: Solution for male APOE4 negative;
        sol_Mpos: Solution for male APOE4 positive.
    """
    sol_Fneg = run_1_model(0, 0, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=InsVar, xi=xi)

    sol_Fpos = run_1_model(0, 1, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=InsVar, xi=xi)

    sol_Mneg = run_1_model(1, 0, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=InsVar, xi=xi)

    sol_Mpos = run_1_model(1, 1, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, InsVar=InsVar, xi=xi)

    return sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos


def main_figure(solt, soly, p, y0, AgeStart, AgeEnd, methodstr, maxstepstr, rtolstr, atolstr, number, CommentModif="",
                IntraneuronalConcentration=False, SkipFirstHalfYear=False, color="b"):
    """
    Figure principale du modèle présentant les graphiques de chaque variable en fonction du temps.
    La figure est générée, puis enregistrée dans le dossier `fig_dir`.

    :param solt: Array des temps (en années). (ndarray, shape (n_points,) ; Time points).
    :type solt: :py:`numpy.ndarray`
    :param soly: Array des arrays des valeurs de chaques paramètres.
                (ndarray, shape (n, n_points) ; Values of the solution at t).
    :type soly: :py:`numpy.ndarray`
    :param y0: Array ou liste des conditions initiales.
    :type y0: :py:`numpy.ndarray` or list
    :param int AgeStart: Âge du début de la simulation (en années).
    :param int AgeEnd: Âge de fin de la simulation (en années).
    :param str methodstr: Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param str maxstepstr: Pas de temps maximal permis pour l'intégration du modèle. Sera ajouté au nom
                d'enregistrement de la figure.
    :param str rtolstr: ``rtol`` utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param str atolstr: ``atol`` utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param int number: Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1 jusqu'à l'obtention d'un nom qui n'est pas
                déjà existant.
    :param CommentModif: Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure. Default: ``""``.
    :param CommentModif: str, optional
    :param IntraneuronalConcentration: Si l'on affiche les graphiques des concentrations intraneuronales pour ABi,
        GSK3, tau et F_i sur leur graphique respectif. Default: ``False``.
    :type IntraneuronalConcentration: bool, optional
    :param SkipFirstHalfYear: Si l'on ignore la première demie année pour l'affichage les valeurs. Default: ``False``.
    :type SkipFirstHalfYear: bool, optional
    :param color: Color to plot the graphs lines.
    :type color: str, optional
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
                ax.set_xlabel(xlabel)
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
        ax1.plot(solt, soly[0, :]/soly[8, :], "g-", alpha=0.7)  # ABi/N
        ax1.set_ylabel(r"$A \beta^{i} / N$", color='g')
        formatter1 = ticker.ScalarFormatter(useMathText=True)
        formatter1.set_scientific(True)
        formatter1.set_powerlimits((-1, 1))
        ax1.yaxis.set_major_formatter(formatter1)

        """Add the plot of GSK3/N in the graph for GSK3"""
        axs[0, 4].set_ylabel(labelname[4], color=color)
        ax2 = axs[0, 4].twinx()
        ax2.plot(solt, soly[4, :] / soly[8, :], "g-", alpha=0.7)  # G/N
        ax2.set_ylabel(r"$GSK3/N$", color='g')
        formatter2 = ticker.ScalarFormatter(useMathText=True)
        formatter2.set_scientific(True)
        formatter2.set_powerlimits((-1, 1))
        ax2.yaxis.set_major_formatter(formatter2)

        """Add the plot of tau/N in the graph for tau"""
        axs[1, 0].set_ylabel(labelname[5], color=color)
        ax3 = axs[1, 0].twinx()
        ax3.plot(solt, soly[5, :] / soly[8, :], "g-", alpha=0.7)  # tau/N
        ax3.set_ylabel(r"$\tau/N$", color='g')
        formatter3 = ticker.ScalarFormatter(useMathText=True)
        formatter3.set_scientific(True)
        formatter3.set_powerlimits((-1, 1))
        ax3.yaxis.set_major_formatter(formatter3)

        """Add the plot of F_i/N in the graph for F_i"""
        axs[1, 1].set_ylabel(labelname[6], color=color)
        ax4 = axs[1, 1].twinx()
        ax4.plot(solt, soly[6, :] / soly[8, :], "g-", alpha=0.7)  # tau/N
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
    # # # ax.set_xlabel(xlabel)
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

    my_path = os.path.abspath(fig_dir)
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_APOE" + APOE + "_" + sex + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + \
              "_atol" + atolstr + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number+1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_APOE" + APOE + "_" + sex + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y_maxstep" + maxstepstr + "_rtol" + rtolstr + \
                  "_atol" + atolstr + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=dpi)

    return FigName


def figure_intracellular_concentrations(sols, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, number, CommentModif=""):
    """Trace et sauvegarde le tracé des concentrations intracellulaires, i.e. de {A\beta^{i}}/N, GSK3/N, \tau/N, {F_i}/N.

    Voir aussi `figure_intracellular_concentrations_SexDiff`.

    :param tuple sols: Solutions des 4 sous-modèles. Résultats de la fonction `run_4_models`: (sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos).
    :param int AgeStart: Âge du début de la simulation (en années).
    :param int AgeEnd: Âge de fin de la simulation (en années).
    :param str methodstr: Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param str max_step: Pas de temps maximal permis pour l'intégration du modèle. Sera ajouté au nom
                d'enregistrement de la figure.
    :param str rtol: ``rtol`` utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param str atol: ``atol`` utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param int number: Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1 jusqu'à l'obtention d'un nom qui n'est pas
                déjà existant.
    :param CommentModif: Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure. Default: ``""``.
    :param CommentModif: str, optional
    :return: None
    """

    sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos = sols

    fig, axs = plt.subplots(nrows=2, ncols=2, sharex="all", figsize=(8, 7), layout="tight")
    axs = axs.flat

    variablenames = [r"$A \beta^{i} / N$", r"$GSK3/N$", r"$\tau/N$", r"$F_i/N$"]
    colors = ["purple", "hotpink", "darkblue", "royalblue"]
    labels = ["F, APOE4+", "F, APOE4-", "M, APOE4+", "M, APOE4-"]
    alphas = [1, 0.8, 1, 0.8]
    s = 0
    for sol in [sol_Fpos, sol_Fneg, sol_Mpos,  sol_Mneg]:
        solt = sol.t / 365
        soly = sol.y

        axs[0].plot(solt, (soly[0, :] / soly[8, :]), color=colors[s], label=labels[s])  # ABi/N
        axs[1].plot(solt, (soly[4, :] / soly[8, :]), color=colors[s], label=labels[s], alpha=alphas[s])  # G/N
        axs[2].plot(solt, (soly[5, :] / soly[8, :]), color=colors[s], label=labels[s], alpha=alphas[s])  # tau/N
        axs[3].plot(solt, (soly[6, :] / soly[8, :]), color=colors[s], label=labels[s], alpha=alphas[s])  # Fi/N
        s = s + 1

    i = 0
    for ax in axs:
        ax.grid()
        ax.set_ylabel(variablenames[i])
        if i >= 2:
            ax.set_xlabel(xlabel)
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
        # notation will be used if exp <= -1 or exp >= 1.
        ax.yaxis.set_major_formatter(formatter)
        i += 1
    # handles, labels = axs[3].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.81, 0.06))  # (0.82, 0.06)
    axs[3].legend()

    """Save the plot as a .png file"""
    today = datetime.date.today()
    date = today.strftime("%y-%m-%d")
    if CommentModif != "":
        CommentModif = "_" + CommentModif
    my_path = os.path.abspath(fig_dir)
    FigName = "Figure_" + date + "_" + "{:0>2d}".format(number) + "_Intracellular_" + methodstr + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + "_mstep=" + max_step + "_rtol" + \
              rtol + "_atol" + atol + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number + 1
        FigName = "Figure_" + date + "_" + "{:0>2d}".format(number) + "_Intracellular_" + methodstr + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + "_mstep=" + max_step + "_rtol" + \
                  rtol + "_atol" + atol + ".png"
    plt.savefig(os.path.join(my_path, FigName), dpi=dpi)


def figure_intracellular_concentrations_SexDiff(sols, AgeStart, AgeEnd, methodstr, max_step, rtol, atol, number, CommentModif=""):
    """Trace et sauvegarde le tracé des concentrations intracellulaires, i.e. de {A\beta^{i}}/N, GSK3/N, \tau/N, {F_i}/N.
    Traces les courbes pour les différents sexes dans des figures distictes.

    Voir aussi `figure_intracellular_concentrations`.

    :param tuple sols: Solutions des 4 sous-modèles. Résultats de la fonction `run_4_models`: (sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos).
    :param int AgeStart: Âge du début de la simulation (en années).
    :param int AgeEnd: Âge de fin de la simulation (en années).
    :param str methodstr: Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param str max_step: Pas de temps maximal permis pour l'intégration du modèle. Sera ajouté au nom
                d'enregistrement de la figure.
    :param str rtol: ``rtol`` utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param str atol: ``atol`` utilisé pour l'intégration du modèle. Sera ajouté au nom d'enregistrement de la figure.
    :param int number: Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1 jusqu'à l'obtention d'un nom qui n'est pas
                déjà existant.
    :param CommentModif: Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure. Default: ``""``.
    :param CommentModif: str, optional
    :return: None
    """

    sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos = sols

    fig, axs = plt.subplots(nrows=4, ncols=2, sharex="all", figsize=(8, 9), layout="tight")

    variablenames = [r"$A \beta^{i} / N$", r"$GSK3/N$", r"$\tau/N$", r"$F_i/N$"]
    colors = ["purple", "hotpink", "darkblue", "royalblue"]
    labels = ["F, APOE4+", "F, APOE4-", "M, APOE4+", "M, APOE4-"]
    alphas = [1, 0.8, 1, 0.8]
    s = 0
    for sol in [sol_Fpos, sol_Fneg]:
        solt = sol.t / 365
        soly = sol.y

        axs[0, 0].plot(solt, soly[0, :] / soly[8, :], color=colors[s], label=labels[s])  # ABi/N
        # y = soly[0, :] / soly[8, :]
        axs[1, 0].plot(solt, soly[4, :] / soly[8, :], color=colors[s], label=labels[s], alpha=alphas[s])  # G/N
        axs[2, 0].plot(solt, soly[5, :] / soly[8, :], color=colors[s], label=labels[s], alpha=alphas[s])  # tau/N
        axs[3, 0].plot(solt, (soly[6, :] / soly[8, :]), color=colors[s], label=labels[s], alpha=alphas[s])  # Fi/N
        axs[0, 0].set_ylabel(variablenames[0])
        axs[1, 0].set_ylabel(variablenames[1])
        axs[2, 0].set_ylabel(variablenames[2])
        axs[3, 0].set_ylabel(variablenames[3])
        s = s + 1

    for sol in [sol_Mpos,  sol_Mneg]:
        solt = sol.t / 365
        soly = sol.y

        axs[0, 1].plot(solt, soly[0, :] / soly[8, :], color=colors[s], label=labels[s])  # ABi/N
        # y = soly[0, :] / soly[8, :]
        axs[1, 1].plot(solt, soly[4, :] / soly[8, :], color=colors[s], label=labels[s], alpha=alphas[s])  # G/N
        axs[2, 1].plot(solt, soly[5, :] / soly[8, :], color=colors[s], label=labels[s], alpha=alphas[s])  # tau/N
        axs[3, 1].plot(solt, (soly[6, :] / soly[8, :]), color=colors[s], label=labels[s], alpha=alphas[s])  # Fi/N
        axs[0, 1].set_ylabel(variablenames[0])
        axs[1, 1].set_ylabel(variablenames[1])
        axs[2, 1].set_ylabel(variablenames[2])
        axs[3, 1].set_ylabel(variablenames[3])
        s = s + 1

    i = 0
    for ax in axs.flat:
        ax.grid()
        if i >= 6:
            ax.set_xlabel(xlabel)
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-1, 1))
        # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
        # notation will be used if exp <= -1 or exp >= 1.
        ax.yaxis.set_major_formatter(formatter)
        if i == 2:
            ax.set_ylim(1.186e-4, 1.189e-4)
        if i == 3:
            ax.set_ylim(3.571e-5, 3.575e-5)
        if i in [2,3]:
            ax.ticklabel_format(useOffset=False)
        i += 1
    # handles, labels = axs[3].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.81, 0.06))  # (0.82, 0.06)
    axs[3, 0].legend()
    axs[3, 1].legend()

    """Save the plot as a .png file"""
    today = datetime.date.today()
    date = today.strftime("%y-%m-%d")
    if CommentModif != "":
        CommentModif = "_" + CommentModif
    my_path = os.path.abspath(fig_dir)
    FigName = "Figure_" + date + "_" + "{:0>2d}".format(number) + "_IntracellularSex_" + methodstr + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + "_mstep=" + max_step + "_rtol" + \
              rtol + "_atol" + atol + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number + 1
        FigName = "Figure_" + date + "_" + "{:0>2d}".format(number) + "_IntracellularSex_" + methodstr + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + "_mstep=" + max_step + "_rtol" + \
                  rtol + "_atol" + atol + ".png"
    plt.savefig(os.path.join(my_path, FigName), dpi=dpi)


def main_figure_4_models(sols, AgeStart, AgeEnd, methodstr, number, CommentModif="",
                         SkipFirstHalfYear=False, orientation="portrait"):
    """
    Figure principale du modèle présentant les graphiques de chaque variable en fonction du temps, pour les quatre
    possibilités : APOE4 +/- et sexe (F/M)
    La figure est générée puis enregistrée dans le dossier `fig_dir` (défini dans le haut de ce fichier).

    Imprime aussi les pourcentages perte neuronale et les conditions initiales de chaque sous-modèle.

    :param sols: tuple. Solution de l'intégration du modèle pour les quatres possibilités, dans l'ordre suivant :
                sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos
    :param AgeStart: int. Âge du début de la simulation (en années).
    :param AgeEnd: int. Âge de fin de la simulation (en années).
    :param methodstr: str. Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param int number: Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1 jusqu'à l'obtention d'un nom qui n'est pas
                déjà existant.
    :param CommentModif: Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure. Default: ``""``.
    :param CommentModif: str, optional
    :param SkipFirstHalfYear: Si l'on ignore la première demie année pour l'affichage les valeurs. Default: ``False``.
    :type SkipFirstHalfYear: bool, optional
    :param orientation: Orientation de la figure: "portrait" (default) ou "paysage".
    :type orientation: str, optional
    :return: FigName : Nom de la figure utilisé pour l'enregistrement.
    """
    sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos = sols

    if orientation == "portrait":
        fig, axs = plt.subplots(nrows=5, ncols=4, sharex="all", figsize=(9, 10), layout="tight")  # "constrained"
    elif orientation == "paysage":
        fig, axs = plt.subplots(nrows=4, ncols=5, sharex="all", layout="tight")
        fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)
    else:
        print("Paramètre `orientation` non valide.")
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
                    ax.plot(solt, soly[i, :], color=colors[s])  # , '.-', ms=2, , alpha=0.8

                    if SkipFirstHalfYear:  # """Plot sans la première demie année."""
                        ax.plot(solt[k:], soly[i, k:], color=colors[s])

                else:  # i == 18
                    ax.plot(solt, soly[i, :], color=colors[s], label=labels[s])  # , '.-', ms=2, , alpha=0.8

                    if SkipFirstHalfYear:  # """Plot sans la première demie année."""
                        ax.plot(solt[k:], soly[i, k:], color=colors[s], label=labels[s])

                if s == 3:  # Si l'on est à la dernière option.
                    ax.grid()
                    if orientation == "portrait":
                        indice = 15
                    else:  # orientation == "paysage"
                        indice = 14
                    if i >= indice:
                        ax.set_xlabel(xlabel)
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

    my_path = os.path.abspath(fig_dir)
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_AllInOne_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number+1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_AllInOne_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=dpi)

    s = 0
    for sol in [sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos]:
        NeuronalLossInPercent = (sol.y[8, -1] - sol.y[8, 0]) / sol.y[8, 0] * 100
        # print(labels[s])
        # print("Concentration de neurones initiale: ", sol.y[8, 0], "\nConcentration de neurones finale: ", sol.y[8, -1],
        #       "\nPourcentage de perte: ", NeuronalLossInPercent)
        print(labels[s] + ": " + str(NeuronalLossInPercent) + " %")
        s = s + 1

    s = 0
    for sol in [sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos]:
        y0 = sol.y[:, 0]
        icNameValue = [str(labelname[i]) + " = " + "{:.5e}".format(y0[i]) for i in np.arange(19)]
        print(labels[s])
        initcond = "Initial conditions used (in g/mL): " + ", ".join(icNameValue[:])
        print(initcond)
        s = s + 1

    return FigName


def fig_one_variable_all_models(variable, sols, AgeStart, AgeEnd, methodstr, number, CommentModif="",
                                SkipFirstHalfYear=False):
    """Trace et sauvegarde le graphique d'une variable en fonction du temps pour les 4 sous-modèles (APOE+/-, F/M).

    :param int variable: Numéro correspondant de la varible entre 0 et 18 (voir le fichier equations.py au besoin).
    :param tuple sols: Solution de l'intégration du modèle pour les quatres possibilités, dans l'ordre suivant :
                sol_Fneg, sol_Fpos, sol_Mneg, sol_Mpos. Résultat de la fonction `run_4_models`.
    :param AgeStart: int. Âge du début de la simulation (en années).
    :param AgeEnd: int. Âge de fin de la simulation (en années).
    :param methodstr: str. Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param int number: Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1 jusqu'à l'obtention d'un nom qui n'est pas
                déjà existant.
    :param CommentModif: Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure. Default: ``""``.
    :param CommentModif: str, optional
    :param SkipFirstHalfYear: Si l'on ignore la première demie année pour l'affichage les valeurs. Default: ``False``.
    :type SkipFirstHalfYear: bool, optional
    :return: FigName : Nom de la figure utilisé pour l'enregistrement.
    """
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

        # ax.set_xlabel(xlabel)
        # ax1.set_xlabel(xlabel)
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
        ax2.set_xlabel(xlabel)
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
            ax.set_xlabel(xlabel)
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

    my_path = os.path.abspath(fig_dir)
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_" + VariablesName[variable] + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number + 1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_" + methodstr + "_" + VariablesName[variable] + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=dpi)

    return FigName


def main_figure_comparison(solt1, soly1, solt2, soly2, labels, AgeStart, AgeEnd, methodstr, number, CommentModif="",
                           SkipFirstHalfYear=False, color1="k", color2="b", orientation="portrait"):
    """
    Figure principale du modèle présentant les graphiques de chaque variable en fonction du temps, pour deux variantes
    du modèle.
    La figure est générée puis enregistrée dans le dossier `fig_dir`.

    Voir aussi `main_figure_comparison_many`.

    :param solt1: Array des temps (en jours) de la première série. (ndarray, shape (n_points,) ; Time points).
    :type solt1: :py:`numpy.ndarray`
    :param soly1: Array des arrays des valeurs de chaque variable pour la première série.
                (ndarray, shape (n, n_points) ; Values of the solution at t).
    :type soly1: :py:`numpy.ndarray`
    :param solt2: Array des temps (en années) de la deuxième série. (ndarray, shape (n_points,) ; Time points).
    :type solt2: :py:`numpy.ndarray`
    :param soly2: Array des arrays des valeurs de chaque variable pour la deuxième série.
                (ndarray, shape (n, n_points) ; Values of the solution at t).
    :type soly2: :py:`numpy.ndarray`
    :param list[str] labels: Liste des labels à associer [label de la première série, label de la deuxième série].
    :param int AgeStart: int. Âge du début de la simulation (en années).
    :param int AgeEnd: int. Âge de fin de la simulation (en années).
    :param str methodstr: Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param int number: int. Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1 jusqu'à l'obtention d'un nom qui n'est pas
                déjà existant.
    :param CommentModif: Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure. Default: ``""``.
    :param CommentModif: str, optional
    :param SkipFirstHalfYear: Si l'on ignore la première demie année pour l'affichage les valeurs. Default: ``False``.
    :type SkipFirstHalfYear: bool, optional
    :param str, optional color1: Couleur de tracé de la première variante du modèle. Défault: ``"k"`` (noir).
    :param str, optional color2: Couleur de tracé de la deuxième variante du modèle. Défault: ``"b"`` (bleu).
    :param orientation: Orientation de la figure: "portrait" (default) ou "paysage".
    :type orientation: str, optional
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
                ax.set_xlabel(xlabel)
            ax.set_ylabel(labelname[i])
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((-1, 1))
            # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10,
            #   scientific notation will be used if exp <= -1 or exp >= 1.
            ax.yaxis.set_major_formatter(formatter)
        i = i+1
    axs.flat[19].remove()  # Retire le dernier graphique

    handles, labels = axs.flat[18].get_legend_handles_labels()

    if orientation == "portrait":
        # fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.80, 0.07))
        # fig.legend(handles, labels, loc='best', bbox_to_anchor=(0.75, 0., 0.5, 0.25))
        fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.81, 0.06))
        axs.flat[15].tick_params('x', labelbottom=True)
    else:  # orientation == "paysage"
        fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.85, 0.09))
        axs.flat[14].tick_params('x', labelbottom=True)  # Ajout des graduations au 14e graphique (hat{M}_anti)

    """Save the plot as a .png file"""
    today = datetime.date.today()
    date = today.strftime("%y-%m-%d")

    if CommentModif != "":
        CommentModif = "_" + CommentModif

    my_path = os.path.abspath(fig_dir)
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_Comparaison_" + methodstr + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number+1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_Comparaison_" + methodstr + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=dpi)

    return FigName


def main_figure_comparison_many(solst, solsy, labels, colors, AgeStart, AgeEnd, methodstr, number, CommentModif="",
                                SkipFirstHalfYear=False, orientation="portrait", alphas=(1, 0.8)):
    """
    Figure principale du modèle présentant les graphiques de chaque variable en fonction du temps, pour plusieurs
    variantes du modèle.
    La figure est générée puis enregistrée dans le dossier `fig_dir`.

    :param solst: Liste des array des temps (en jours). [(ndarray, shape (n_points,) ; Time points)].
    :type solst: list[:py:`numpy.ndarray`]
    :param solsy: Liste des array des arrays des valeurs de chaque variable.
                (ndarray, shape (n, n_points) ; Values of the solution at t).
    :type solsy: list[:py:`numpy.ndarray`]
    :param list[str] labels: Liste des labels à associer à chaque variante du modèle.
    :param list[str] colors: Liste des couleurs à associer à chaque variante du modèle.
    :param int AgeStart: Âge du début de la simulation (en années).
    :param int AgeEnd: Âge de fin de la simulation (en années).
    :param str methodstr: Nom de la méthode utilisée pour l'intégration du modèle. Sera ajouté au nom d'enregistrement
                de la figure.
    :param int number: Nombre pour identifier la figue. Sera mis dans son nom d'enregistrement. Si ce nombre est déjà
                utilisé pour la date du jour, il sera incrémenté de 1 jusqu'à l'obtention d'un nom qui n'est pas
                déjà existant.
    :param CommentModif: Commentaire pour le titre de la figure résumant les modifications, s'il y a
                lieu. Sera ajouté à la fin du nom d'enregistrement de la figure. Default: ``""``.
    :param CommentModif: str, optional
    :param SkipFirstHalfYear: Si l'on ignore la première demie année pour l'affichage les valeurs. Default: ``False``.
    :type SkipFirstHalfYear: bool, optional
    :param orientation: Orientation de la figure: "portrait" (default) ou "paysage".
    :type orientation: str, optional
    :param alphas: Tuple ou liste des alphas à utiliser pour les tracés, valeurs qui définissent la transparence du
        tracé. Chaque valeur doit être entre 0 et 1: 0 -> transparent, 1 -> opaque. Doit avoir une valeur par variante
        du modèle.
        Default: ``(1, 0.8)``.
    :type alphas: tuple(float) or list(float), optional
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
                    ax.plot(solst[j], solsy[j][i, ks[j]:], color=colors[j], alpha=alphas[j], label=labels[j],
                            linewidth=1.5)
            else:
                for j in range(nbmodel):
                    ax.plot(solst[j], solsy[j][i, :], color=colors[j], alpha=alphas[j], label=labels[j],
                            linewidth=1.5)  # , '.-', ms=2

            ax.grid()
            if orientation == "portrait":
                indice = 15
            else:  # orientation == "paysage"
                indice = 14
            if i >= indice:
                ax.set_xlabel(xlabel)
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

    my_path = os.path.abspath(fig_dir)
    FigName = "Figure_" + date + "_" + f"{number:02}" + "_Comparaison_" + methodstr + "_" + \
              str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"
    while os.path.exists(os.path.join(my_path, FigName)):
        number = number+1
        FigName = "Figure_" + date + "_" + f"{number:02}" + "_Comparaison_" + methodstr + "_" + \
                  str(AgeEnd - AgeStart).replace(".", "") + "y" + CommentModif + ".png"

    plt.savefig(os.path.join(my_path, FigName), dpi=dpi)

    return FigName
