import math
import numpy as np
import parameters as param
p = param.Parameters()


def InitialConditions(AgeStart=30):
    """
    Function that defines de vector of initial conditions of the model, in g/mL.
    :param AgeStart: Age for which the integration of the model will start. # TODO: Si utilise toujours pas, à retirer.
    :return: y0 : a vector with the initial conditions.
    """
    y0 = np.zeros(19)

    """AB^i (Amyloid-beta monomer intracell.)"""
    # 1e-6 changé pour éviter saut départ...
    y0[0] = p.lambda_ABi * (1 + p.AP * p.delta_APi) / p.d_ABi

    """AB_m^o (Amyloid-beta monomer extracell.)"""
    # y0[1] = 1e-11   # Todo ou  ??
    # y0[1] = (p.lambda_ABmo * (1 + p.AP * p.delta_APm) + p.lambda_AABmo) / (p.d_ABmo(AgeStart*365) + p.kappa_ABmoABoo *
    #                                                                        (1 + p.AP * p.delta_APmo))
    # A = p.kappa_ABmoABoo * (1 + p.AP * p.delta_APmo)
    # B = p.d_ABmo(AgeStart*365)
    # C = -(p.lambda_ABmo * (1 + p.AP * p.delta_APm) + p.lambda_AABmo)
    # y0[1] = (-B + math.sqrt((B**2) - (4 * A * C))) / (2 * A)
    # Prendre + donne nég. Avec AgeStart=30, = 4.233652303020852e-11.
    # Par expérience, pour pas avoir de différence au départ, prendre (22-09-07_..._12 vs _13):
    y0[1] = 0

    """AB_o^o (Amyloid-beta oligomers extracell.)"""
    # y0[2] = 5e-13  # 0
    # y0[2] = p.kappa_ABmoABoo * (1 + p.AP * p.delta_APmo) * (y0[1] ** 2) / (p.d_ABoo + p.kappa_ABooABpo)
    # A = p.kappa_ABooABpo
    # B = p.d_ABoo
    # C = - p.kappa_ABmoABoo * (1 + p.AP * p.delta_APmo) * (y0[1] ** 2)
    # y0[2] = (-B + math.sqrt((B**2) - (4 * A * C))) / (2 * A)
    # Prendre + donne nég. Avec AgeStart=30, = 6.791648132775663e-17.
    # Par expérience, pour pas avoir de différence au départ, prendre (22-09-07_..._12 vs _13):
    y0[2] = 0

    """AB_p^o (Amyloid-beta plaque extracell.)"""
    y0[3] = 1e-18

    """G (GSK3)"""
    # y0[4] = p.lambda_InsG / p.d_G  # ~= 0.16168
    y0[4] = p.G_0   # ~ 5.344e-05 (for women)

    """tau (tau proteins)"""
    # y0[5] = (p.lambda_tau + p.lambda_Gtau)/p.d_tau
    # y0[5] = 6e-7
    # (p.lambda_tau + p.lambda_Gtau * (y0[4] / p.G_0))/p.d_tau => avant fig _39
    # = 2.57e-5    # Hao: Concentration of tau proteins is, in health, 137 pg/ml and, in AD, 490 pg/ml
    # A = p.kappa_tauFi
    # B = p.d_tau
    # C = -(p.lambda_tau + p.lambda_Gtau)
    # y0[5] = (-B + math.sqrt((B**2) - (4 * A * C))) / (2 * A)
    # Prendre "+" donne nég. Donc, prend "-" et donne ~ 6.4947e-07.
    # Si prend lambda_Gtau * 1e-2 => ~6.6886e-09
    y0[5] = 1e-10

    """F_i (NFT inside the neurons)"""
    # y0[6] = p.kappa_tauFi * (y0[5] ** 2) / p.d_Fi
    y0[6] = 0

    """F_o (NFT outside the neurons)"""
    y0[7] = 0

    """N (Living neurons)"""
    y0[8] = p.N_0

    """A (Activated astrocytes)"""
    # y0[9] = 0.14 /2
    # T_alpha_0 = 2.79e-5
    # Q = (p.kappa_ABpoA * y0[3] + p.kappa_TaA * T_alpha_0)
    # y0[9] = (Q * p.A_max) / (Q + p.d_A)
    y0[9] = 0  # p.A_0/1000

    """M_NA (Resting microglia)"""
    if p.S == 0:  # women
        y0[10] = 3.811e-2
    else:  # men
        y0[10] = 3.193e-2

    """M_pro (Proinflammatory microglia)"""
    y0[11] = 1e-12  # y0[10] * (p.beta / (p.beta + 1))

    """M_anti (Anti-inflammatory microglia)"""
    # y0[12] = 1e-12  # y0[10] * (1 / (p.beta + 1))
    # Changement 1e-12 pour 1e-5 (fig 22-09-09_..._03 vs _01) pas vrm amélioration

    """hat{M}_pro (Proinflammatory macrophages)"""
    # y0[13] = p.Mprohateq/1.5  # 0
    y0[13] = 1e-12
    # y0[13] = p.Mhatmax/2

    """hat{M}_anti (Anti-inflammatory macrophages)"""
    # y0[14] = 1e-9  # ou (p.kappa_TB * y0[15])/p.d_Mantihat, si y0[15] defini avant # Hao: 0
    y0[14] = 1e-12
    # y0[14] = p.Mhatmax/3

    # Changement 1e-12 pour 1e-5 (fig 22-09-09_..._03 vs _01) pas vrm amélioration

    """T_{beta} (TGF-beta)"""
    # y0[15] = p.K_TbM * 1e-7  # = 5.9e-18 (i.e. même diff ordre de grandeur entre [Tb]_0 et K_TbM que [Ta]_0 et K_TaM).
    #                          # Figure_22-09-13_..._05...: Modif pas utile, atteint équilibre avant de réaugmenter.
    y0[15] = (p.kappa_MantiTb * y0[12] + p.kappa_MhatantiTb * y0[14])/p.d_Tb  # = 3.774586223657386e-22
    # 1.0e-6

    """I_10 (IL-10 = Interleukin 10)"""
    y0[16] = (p.kappa_MantiI10 * y0[12] + p.kappa_MhatantiI10 * y0[14]) / p.d_I10  # = 2.6173735792085045e-17
    # Hao : 1.0e-5

    """T_{alpha} (TNF-alpha)"""
    # y0[17] = 1e-12
    y0[17] = (p.kappa_MproTa * y0[11] + p.kappa_MhatproTa * y0[13]) / p.d_Ta  # nouv kappas (Fadok98): 7.31e-21
    # 1.1620103844701855e-19 (av)
    # Hao 2e-5

    """P (MCP-1)"""
    y0[18] = (p.kappa_MproP * y0[11] + p.kappa_MhatproP * y0[13] + p.kappa_AP * y0[9]) / p.d_P  # 3.975362086617884e-19
    # Hao: 5e-9

    return y0
