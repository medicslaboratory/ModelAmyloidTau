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


def main_figure(solt, soly, p, y0, AgeStart, AgeEnd, methodstr, maxstepstr, rtolstr, atolstr, number, CommentModif="",
                IntraneuronalConcentration=False, SkipFirstHalfYear=False):
    """
    Figure principale du modèle présentant les graphiques de chaque variable en fonction du temps.
    La figure est générée puis enregistrée dans le dossier "Figure".

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
            ax.plot(solt, soly[i, :])  # , '.-', ms=2

            if SkipFirstHalfYear:  # """Plot sans la première demie année."""
                ax.plot(solt[k:], soly[i, k:])

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
    if IntraneuronalConcentration:
        """Add the plot of ABi/N in the graph for ABi"""
        axs[0, 0].set_ylabel(labelname[0], color='b')
        ax1 = axs[0, 0].twinx()
        ax1.plot(solt, soly[0, :]/soly[8, :], "g-")  # ABi/N
        ax1.set_ylabel(r"$A \beta^{i} / N$", color='g')
        formatter1 = ticker.ScalarFormatter(useMathText=True)
        formatter1.set_scientific(True)
        formatter1.set_powerlimits((-1, 1))
        ax1.yaxis.set_major_formatter(formatter1)

        """Add the plot of GSK3/N in the graph for GSK3"""
        axs[0, 4].set_ylabel(labelname[4], color='b')
        ax2 = axs[0, 4].twinx()
        ax2.plot(solt, soly[4, :] / soly[8, :], "g-")  # G/N
        ax2.set_ylabel(r"$GSK3/N$", color='g')
        formatter2 = ticker.ScalarFormatter(useMathText=True)
        formatter2.set_scientific(True)
        formatter2.set_powerlimits((-1, 1))
        ax2.yaxis.set_major_formatter(formatter2)

        """Add the plot of tau/N in the graph for tau"""
        axs[1, 0].set_ylabel(labelname[5], color='b')
        ax3 = axs[1, 0].twinx()
        ax3.plot(solt, soly[5, :] / soly[8, :], "g-")  # tau/N
        ax3.set_ylabel(r"$\tau/N$", color='g')
        formatter3 = ticker.ScalarFormatter(useMathText=True)
        formatter3.set_scientific(True)
        formatter3.set_powerlimits((-1, 1))
        ax3.yaxis.set_major_formatter(formatter3)

        """Add the plot of F_i/N in the graph for F_i"""
        axs[1, 1].set_ylabel(labelname[6], color='b')
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
    # # # ax.set_xlabel('Age (years)')
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


def fig_neuronal_loss_rates(solt, soly, p):
    dNdtFi = - p.d_FiN * (1 / (1 + np.exp(- p.n * (soly[6, :] - p.K_Fi)))) * soly[8, :]
    dNdtTa = - p.d_TaN * (soly[17, :] / (soly[17, :] + p.K_Ta)) * (1 / (1 + (soly[16, :] / p.K_I10))) * soly[8, :]

    """Pour un graph."""
    # fig, ax1 = plt.subplots(1, 1)
    # plt.suptitle("Taux de perte neuronale")

    """Un graph : Deux axes des x."""
    # ax1.plot(solt, dNdtFi + dNdtTa, "r-", label=r"Total")
    # ax1.plot(solt, dNdtFi, "b-", label=r"Par $F_i$")  # loss by F_i
    # ax1.set_ylabel(r"$dN/dt$ total et par $F_i$", color="k")
    # # ax1.set_ylabel(r"$dN/dt$ par $F_i$", color="b")
    # formatter = ticker.ScalarFormatter(useMathText=True)
    # formatter.set_scientific(True)
    # formatter.set_powerlimits((-1, 1))
    # # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # # notation will be used if exp <= -1 or exp >= 1.
    # ax1.yaxis.set_major_formatter(formatter)
    # ax2 = ax1.twinx()
    # ax2.plot(solt, dNdtTa, "g-", label=r"Par $T_\alpha$")  # loss by TNFa
    # ax2.set_ylabel(r"$dN/dt$ par $T_\alpha$", color='g')
    # formatter = ticker.ScalarFormatter(useMathText=True)
    # formatter.set_scientific(True)
    # formatter.set_powerlimits((-1, 1))
    # # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # # notation will be used if exp <= -1 or exp >= 1.
    # ax2.yaxis.set_major_formatter(formatter)
    # ax1.set_xlabel("Age (years)")
    # ax1.legend(loc='center left')
    # ax1.grid()

    """Un graph : Une axe des x. Recommande moins."""
    # ax1.plot(solt, dNdtFi + dNdtTa, "r-", label=r"Total")
    # ax1.plot(solt, dNdtFi, "b-", label=r"Par $F_i$")  # loss by F_i
    # ax1.plot(solt, dNdtTa, "g-", label=r"Par $T_\alpha$")  # loss by TNFa
    # ax1.legend()
    # ax1.set_ylabel(r"$dN/dt$ par facteur")
    # ax1.set_xlabel("Age (years)")
    # formatter = ticker.ScalarFormatter(useMathText=True)
    # formatter.set_scientific(True)
    # formatter.set_powerlimits((-1, 1))
    # # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # # notation will be used if exp <= -1 or exp >= 1.
    # ax1.yaxis.set_major_formatter(formatter)
    # ax1.grid()

    """Pour deux graphs superposés"""
    fig, axs = plt.subplots(2, 1, sharex="all")
    plt.suptitle("Taux de perte neuronale")
    axs[0].plot(solt, dNdtFi + dNdtTa, "r-", label=r"Total")  # total loss
    axs[0].plot(solt, dNdtFi, "b-", label=r"Par $F_i$")  # loss by F_i
    axs[0].set_ylabel(r"$dN/dt$ total et par $F_i$")
    axs[0].legend()
    # axs[0].plot(solt, dNdtFi)
    # axs[0].set_ylabel(r"$dN/dt$ par $F_i$")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    axs[0].yaxis.set_major_formatter(formatter)
    axs[0].grid()
    axs[1].plot(solt, dNdtTa)  # loss by T_a
    axs[1].set_ylabel(r"$dN/dt$ par $T_\alpha$")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    axs[1].yaxis.set_major_formatter(formatter)
    axs[1].grid()

    plt.tight_layout()

    return fig


def fig_microglia_activation_rates(solt, soly, p):
    """
    Graphique des taux d'activation des microglies.
    :param solt:
    :param soly:
    :return:
    """
    fig, ax = plt.subplots(1, 1)
    plt.suptitle("Taux d'activation pour les microglies")
    M_activ = p.kappa_FoM * (soly[7, :] / (soly[7, :] + p.K_Fo)) * soly[10, :] + \
              p.kappa_ABooM * (soly[2, :] / (soly[2, :] + p.K_ABooM)) * soly[10, :]
    epsilon_Ta = soly[17, :] / (soly[17, :] + p.K_TaAct)
    epsilon_I10 = soly[16, :] / (soly[16, :] + p.K_I10Act)
    ax.plot(solt, ((p.beta * epsilon_Ta) / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ, "b-",
            label="Taux activ. pro-infl.")
    ax.plot(solt, (epsilon_I10 / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ, "r-",
            label="Taux activ. anti-infl.")
    ax.set_ylabel("Taux")
    ax.plot(solt, M_activ, "k", label=r"$M_{activ}$")
    ax.plot(solt, p.kappa_FoM * (soly[7, :] / (soly[7, :] + p.K_Fo)) * soly[10, :],
            label=r"$M_{activ}$ par $F_o$")
    ax.plot(solt, p.kappa_ABooM * (soly[2, :] / (soly[2, :] + p.K_ABooM)) * soly[10, :],
            label=r"$M_{activ}$ par ${A\beta}_{o}^{o}$")
    ax.grid()
    ax.set_xlabel('Age (years)')
    # ax2 = ax.twinx()
    AntiToPro = p.kappa_TaManti * (soly[17, :] / (soly[17, :] + p.K_TaM)) * soly[12, :]
    ProToAnti = p.kappa_TbMpro * (soly[15, :] / (soly[15, :] + p.K_TbM)) * soly[11, :]
    ax.plot(solt, AntiToPro, "g-", label=r"Taux conversion anti $\rightarrow$ pro")
    ax.plot(solt, ProToAnti, "m-", label=r"Taux conversion pro $\rightarrow$ anti")
    # ax2.plot(solt, AntiToPro, "g-", label="rate anti -> pro")
    # ax2.set_ylabel("Rate anti to pro", color='g')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    ax.yaxis.set_major_formatter(formatter)
    ax.legend()
    plt.tight_layout()

    return fig


def fig_astrocyte_activation_rates(solt, soly, p):
    """
    Graphiques des taux d'activation des astrocytes.
    :param solt:
    :param soly:
    :return:
    """
    dAdtABpo = p.kappa_ABpoA * soly[3, :] * (p.A_max - soly[9, :])
    dAdtTa = p.kappa_TaA * soly[17, :] * (p.A_max - soly[9, :])

    """Un graphique"""
    fig, ax1 = plt.subplots(1, 1)
    plt.suptitle("Taux d'activation des astrocytes")

    """Un graph : Deux axes des x."""
    ax1.plot(solt, dAdtABpo + dAdtTa, "r-", label=r"Total")  # total activation
    ax1.plot(solt, dAdtTa, "b-", label=r"Par $T_\alpha$")  # activ by T_a
    ax1.set_ylabel(r"$dA/dt$ total et par $T_\alpha$", color="k")
    ax1.legend(loc='center left')
    # ax1.set_ylabel(r"$dA/dt$ par $T_\alpha$", color="b")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    ax1.yaxis.set_major_formatter(formatter)
    ax2 = ax1.twinx()
    ax2.plot(solt, dAdtABpo, "g-", label=r"Par $A\beta_p^o$")  # activ by ABpo
    ax2.set_ylabel(r"$dN/dt$ par $A\beta_p^o$", color='g')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    ax2.yaxis.set_major_formatter(formatter)
    ax1.set_xlabel("Age (years)")
    ax1.grid()

    plt.tight_layout()

    return fig


def main_figure_4_models(AgeStart, AgeEnd, methodstr, max_step, rtol, atol, number, CommentModif="",
                         SkipFirstHalfYear=False):

    p_Fneg = param.Parameters(Sex=0, APOE_status=0)  # F ; APOE-
    y0_Fneg = InitialConditions(p_Fneg, AgeStart)
    sol_Fneg = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Fneg, method=methodstr, args=[0, 0],
                         max_step=max_step, rtol=rtol, atol=atol)

    p_Fpos = param.Parameters(Sex=0, APOE_status=1)  # F ; APOE+
    y0_Fpos = InitialConditions(p_Fpos, AgeStart)
    sol_Fpos = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Fpos, method=methodstr, args=[0, 1],
                         max_step=max_step, rtol=rtol, atol=atol)

    p_Mneg = param.Parameters(Sex=1, APOE_status=0)  # M ; APOE-
    y0_Mneg = InitialConditions(p_Mneg, AgeStart)
    sol_Mneg = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Mneg, method=methodstr, args=[1, 0],
                         max_step=max_step, rtol=rtol, atol=atol)

    p_Mpos = param.Parameters(Sex=1, APOE_status=1)  # M ; APOE+
    y0_Mpos = InitialConditions(p_Mpos, AgeStart)
    sol_Mpos = solve_ivp(eqns.ODEsystem, [365 * AgeStart, 365 * AgeEnd], y0_Mpos, method=methodstr, args=[1, 1],
                         max_step=max_step, rtol=rtol, atol=atol)

    fig, axs = plt.subplots(nrows=4, ncols=5, sharex="all")  # If no text under the figure, add: layout="constrained"
    fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)

    # Making a list for Label names in the plot
    labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$',
                 r'$\tau$', '$F_i$', '$F_o$', '$N$', '$A$', '$M_{NA}$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$',
                 r'$\hat{M}_{anti}$', r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

    colors = ["hotpink", "purple", "royalblue", "darkblue"]
    labels = ["F, APOE-", "F, APOE+", "M, APOE-", "M, APOE+"]
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
                    if i >= 14:
                        ax.set_xlabel('Age (years)')
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
    fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.85, 0.08))

    axs.flat[14].tick_params('x', labelbottom=True)  # Ajout des graduations au 14e graphique (hat{M}_anti)
    axs[3, 4].remove()  # Retire le dernier graphique

    plt.subplots_adjust(hspace=.2, wspace=.2)
    plt.tight_layout()

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
