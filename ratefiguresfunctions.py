# Date : 13 Octobre 2022
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


def fig_neuronal_loss_rates(solt, soly, p):
    solt = solt / 365

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
    # # ax1.set_xlabel("Age (years)")
    # ax1.set_xlabel('Âge (années)')
    # ax1.legend(loc='center left')
    # ax1.grid()

    """Un graph : Une axe des x. Recommande moins."""
    # ax1.plot(solt, dNdtFi + dNdtTa, "r-", label=r"Total")
    # ax1.plot(solt, dNdtFi, "b-", label=r"Par $F_i$")  # loss by F_i
    # ax1.plot(solt, dNdtTa, "g-", label=r"Par $T_\alpha$")  # loss by TNFa
    # ax1.legend()
    # ax1.set_ylabel(r"$dN/dt$ par facteur")
    # # ax1.set_xlabel("Age (years)")
    # ax1.set_xlabel('Âge (années)')
    # formatter = ticker.ScalarFormatter(useMathText=True)
    # formatter.set_scientific(True)
    # formatter.set_powerlimits((-1, 1))
    # # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # # notation will be used if exp <= -1 or exp >= 1.
    # ax1.yaxis.set_major_formatter(formatter)
    # ax1.grid()

    """Pour deux graphs superposés"""
    fig, axs = plt.subplots(2, 1, sharex="all", layout="tight", figsize=(6.4, 7))
    # plt.suptitle("Taux de perte neuronale")
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
    axs[1].plot(solt, dNdtTa, "g-")  # loss by T_a
    axs[1].set_ylabel(r"$dN/dt$ par $T_\alpha$")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    axs[1].yaxis.set_major_formatter(formatter)
    axs[1].grid()
    axs[1].set_xlabel("Âge (années)")

    # plt.tight_layout()

    return fig


def fig_neuronal_loss_rates_V2(solt, soly, p):
    solt = solt / 365

    RatedNdtFi = - p.d_FiN * (1 / (1 + np.exp(- p.n * (soly[6, :] - p.K_Fi))))
    RatedNdtTa = - p.d_TaN * (soly[17, :] / (soly[17, :] + p.K_Ta)) * (1 / (1 + (soly[16, :] / p.K_I10)))
    dNdtFi = RatedNdtFi * soly[8, :]
    dNdtTa = RatedNdtTa * soly[8, :]

    """Pour graphs superposés"""
    fig, axs = plt.subplots(3, 1, sharex="all", layout="tight", figsize=(6.4, 7.3))
    # plt.suptitle("Taux de perte neuronale")
    axs[0].plot(solt, dNdtFi + dNdtTa, "r-", label=r"$\frac{dN}{dt} (t)$ (g/mL/j)")  # total loss
    axs[0].plot(solt, dNdtFi, "-", color="mediumblue", label=r"$\frac{dN}{dt}$ par $F_i$ (g/mL/j)")  # loss by F_i
    # axs[0].plot(solt, RatedNdtFi, "-", color="dodgerblue", label=r"Taux de $dN/dt$ par $F_i$ (/j)")
    axs[0].set_ylabel(r"Taux")
    axs[0].legend()
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    axs[0].yaxis.set_major_formatter(formatter)
    axs[0].grid()

    # axs[1].plot(solt, dNdtFi + dNdtTa, "r-", label=r"$\frac{dN}{dt} (t)$ (g/mL/j)")  # total loss
    # axs[1].plot(solt, dNdtFi, "-", color="mediumblue", label=r"$\frac{dN}{dt}$ par $F_i$ (g/mL/j)")  # loss by F_i
    axs[1].plot(solt, RatedNdtFi, "-", color="dodgerblue", label=r"Taux de $\frac{dN}{dt}$ par $F_i$ (/j)")
    axs[1].set_ylabel(r"Taux")
    axs[1].legend()
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    axs[1].yaxis.set_major_formatter(formatter)
    axs[1].grid()

    # axs[0].set_ylim(-2.35e-4, -2.0e-4)
    # axs[1].set_ylim(-5.48e-4, -5.44e-4)
    #
    # axs[0].spines.bottom.set_visible(False)
    # axs[1].spines.top.set_visible(False)
    # # axs[0].xaxis.tick_top()
    # # axs[0].tick_params(labeltop=False)  # don't put tick labels at the top
    # axs[1].xaxis.tick_bottom()

    # d = .5  # proportion of vertical to horizontal extent of the slanted line
    # kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
    #               linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    # axs[0].plot([0, 1], [0, 0], transform=axs[0].transAxes, **kwargs)
    # axs[1].plot([0, 1], [1, 1], transform=axs[1].transAxes, **kwargs)

    axs[2].plot(solt, dNdtTa, "g-", label=r"$\frac{dN}{dt}$ par $T_\alpha$ (g/mL/j)")  # loss by T_a
    axs[2].plot(solt, RatedNdtTa, "-", color="limegreen", label=r"Taux de $\frac{dN}{dt}$ par $T_\alpha$ (/j)")
    axs[2].set_ylabel(r"Taux")
    axs[2].legend()
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    axs[2].yaxis.set_major_formatter(formatter)
    axs[2].grid()
    axs[2].set_xlabel("Âge (années)")

    axs[0].set_ylim(-2.5e-4, -2.0e-4)
    axs[1].set_ylim(-5.457621565e-4 + -8.05e-14, -5.457621565e-4 + -3.0e-14)
    axs[2].set_ylim(-1.25e-5, 0.05e-5)

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
    # ax.set_xlabel('Age (years)')
    ax.set_xlabel('Âge (années)')
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


def fig_microglia_activation_rates_V2(solt, soly, p):
    """
    Graphique des taux d'activation des microglies.
    :param solt:
    :param soly:
    :return:
    """
    solt = solt / 365

    fig, axs = plt.subplots(3, 1, sharex="all", layout="tight", figsize=(6.4, 7.5))
    # plt.suptitle("Taux d'activation pour les microglies")

    M_activ = p.kappa_FoM * (soly[7, :] / (soly[7, :] + p.K_Fo)) * soly[10, :] + \
              p.kappa_ABooM * (soly[2, :] / (soly[2, :] + p.K_ABooM)) * soly[10, :]
    epsilon_Ta = soly[17, :] / (soly[17, :] + p.K_TaAct)
    epsilon_I10 = soly[16, :] / (soly[16, :] + p.K_I10Act)

    # axs[0].plot(solt, ((p.beta * epsilon_Ta) / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ, "b-",
    #         label="Taux activ. pro-infl.")
    # axs[0].plot(solt, (epsilon_I10 / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ, "r-",
    #         label="Taux activ. anti-infl.")
    axs[0].plot(solt, M_activ, "k", label=r"$M_{activ}$")
    axs[0].plot(solt, p.kappa_FoM * (soly[7, :] / (soly[7, :] + p.K_Fo)) * soly[10, :],
            label=r"$M_{activ}$ par $F_o$")
    axs[0].plot(solt, p.kappa_ABooM * (soly[2, :] / (soly[2, :] + p.K_ABooM)) * soly[10, :],
            label=r"$M_{activ}$ par ${A\beta}_{o}^{o}$")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    axs[0].yaxis.set_major_formatter(formatter)
    axs[0].set_ylabel("Taux")
    axs[0].grid()
    axs[0].legend()
    # ax.set_xlabel('Age (years)')
    # axs[0].set_xlabel('Âge (années)')

    axs[1].plot(solt, ((p.beta * epsilon_Ta) / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ, "b-",
                label="Taux activ. pro-infl.")
    axs[1].plot(solt, (epsilon_I10 / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ, "r-",
                label="Taux activ. anti-infl.")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    axs[1].yaxis.set_major_formatter(formatter)
    axs[1].set_ylabel("Taux")
    axs[1].grid()
    axs[1].legend()

    AntiToPro = p.kappa_TaManti * (soly[17, :] / (soly[17, :] + p.K_TaM)) * soly[12, :]
    ProToAnti = p.kappa_TbMpro * (soly[15, :] / (soly[15, :] + p.K_TbM)) * soly[11, :]

    axs[2].plot(solt, AntiToPro, "g-", label=r"Taux conversion anti $\rightarrow$ pro")
    axs[2].plot(solt, ProToAnti, "m-", label=r"Taux conversion pro $\rightarrow$ anti")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    # formatter.set_powerlimits((-1, 1)): For a number representable as a * 10^{exp} with 1<abs(a)<=10, scientific
    # notation will be used if exp <= -1 or exp >= 1.
    axs[2].yaxis.set_major_formatter(formatter)
    axs[2].set_ylabel("Taux")
    axs[2].grid()
    axs[2].legend()
    axs[2].set_xlabel('Âge (années)')

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
    # ax1.set_xlabel("Age (years)")
    ax1.set_xlabel('Âge (années)')
    ax1.grid()

    plt.tight_layout()

    return fig