# Date : 17 January 2022
# Autor : Éléonore Chamberland

# This file contains the equations of the model.
# It defines a function that generates the ODE system of the model.

# import math
import numpy as np
import parameters as param


def ODEsystem(t, y, Sex, APOE_status, InsVar=True, xi=1):
    """
    Function that defines the equations of the model.
    :param t:
    :param y:
        y[0]   AB^i (Amyloid-beta monomer inside the neurons)
        y[1]   Amyloid-beta monomer outside the neurons
        y[2]   Amyloid-beta oligomer outside
        y[3]   Amyloid-beta plaque outside the neurons
        y[4]   GSK-3
        y[5]   tau proteins
        y[6]   F_i (NFT inside the neurons)
        y[7]   F_o (NFT outside the neurons)
        y[8]   Living neurons
        y[9]   Astrocytes
        y[10]  M (Microglia)
        y[11]  M_pro (Proinflammatory microglia)
        y[12]  M_anti (anti-inflammatory microglias)
        y[13]  M_pro^hat (M_pro macrophages)
        y[14]  M_anti^hat (M_anti macrophages)
        y[15]  T_{beta} (TGF-beta)
        y[16]  I_10 (IL-10 = Interleukin 10)
        y[17]  T_{alpha} (TNF-alpha)
        y[18]  P (MCP-1)
    """
    p = param.Parameters(Sex, APOE_status, xi)

    print(t/365, "ans")

    dydt = np.zeros(19)

    # Living neurons (N)
    # dydt[8] = -p.d_FiN * (y[6] / (y[6] + p.K_Fi)) * y[8] \
    #           - p.d_TaN * (y[17] / (y[17] + p.K_Ta)) * (1 / (1 + (y[16] / p.K_I10))) * y[8]
    dydt[8] = -p.d_FiN * (1 / (1 + np.exp(- p.n * (y[6] - p.K_Fi) / p.K_Fi))) * y[8] \
              - p.d_TaN * (y[17] / (y[17] + p.K_Ta)) * (1 / (1 + (y[16] / p.K_I10))) * y[8]
    # Division ("/ p.K_Fi") dans sigmoïde utile? Retiré, ok? NON!
    # dydt[8] = -p.d_FiN * (1 / (1 + np.exp(- p.n * (y[6] - p.K_Fi)))) * y[8] \
    #           - p.d_TaN * (y[17] / (y[17] + p.K_Ta)) * (1 / (1 + (y[16] / p.K_I10))) * y[8]
    # dydt[8] = - p.d_TaN * (y[17] / (y[17] + p.K_Ta)) * (1 / (1 + (y[16] / p.K_I10))) * y[8]

    # Amyloid-beta monomer inside the neurons (AB^i)
    dydt[0] = p.lambda_ABi * (1 + p.AP * p.delta_APi) * (y[8] / p.N_0) - p.d_ABi * y[0] - (y[0] / y[8]) * abs(dydt[8])

    # Amyloid-beta monomer outside the neurons (AB_m^o)
    dydt[1] = (y[0] / y[8]) * abs(dydt[8]) + p.lambda_ABmo * (1 + p.AP * p.delta_APm) * (y[8] / p.N_0) \
              + p.lambda_AABmo * (y[9] / p.A_0) - p.kappa_ABmoABoo * (1 + p.AP * p.delta_APmo) * (y[1] ** 2) \
              - p.d_ABmo(t) * y[1]

    # Amyloid-beta oligomers outside (AB_o^o)
    dydt[2] = p.kappa_ABmoABoo * (1 + p.AP * p.delta_APmo) * (y[1] ** 2) - p.kappa_ABooABpo * (y[2] ** 2) \
              - p.d_ABoo * y[2]

    # Amyloid-beta plaque outside the neurons (AB_p^o)
    # Deuxième terme problématique pour les plaques (22-09-08_..._14 vs _15 ; 15 retire le terme), mais ne règle pas les
    # autres. Donc, vrai problème ailleurs.
    dydt[3] = p.kappa_ABooABpo * (y[2] ** 2) - ((p.d_MantiABpo * y[12] + p.d_hatMantiABpo * y[14])
                                                * (1 + p.AP * p.delta_APdp) * (y[3] / (y[3] + p.K_ABpo)))

    # Glycogen synthase kinase-type 3 (GSK-3) (G)
    # dydt[4] = p.lambda_InsG * (p.Ins(t, p.S) / p.Ins_0) - p.d_G * y[4]
    # TODO: À approuver: Augmentation de l'activité de GSK3 quand la concentration d'insuline diminue
    #  (Jolivat08; DOI: 10.1002/jnr.21787), donc inverse division.
    #  Fait que augmentation de l'activité de la GSK3 (22-09-08_..._08 vs _09). Bien!
    #  À confirmer, déjà ajusté dans Latex.
    # Ajout du "* (y[8]/p.N_0)" et "- (y[4] / y[8]) * abs(dydt[8])" - OK.
    if InsVar:
        dydt[4] = p.lambda_InsG * (p.Ins_0 / p.Ins(t, p.S)) * (y[8] / p.N_0) - p.d_G * y[4] - (y[4] / y[8]) * abs(dydt[8])  # ****
    # # Ajout - (y[4] / y[8]) * abs(dydt[8]), pas impact.
    # dydt[4] = p.lambda_InsG * (p.Ins_0 / p.Ins(t, p.S)) - p.d_G * y[4] - (y[4] / y[8]) * abs(dydt[8])
    else:  # Insuline cste
        dydt[4] = p.lambda_InsG * (y[8] / p.N_0) - p.d_G * y[4] - (y[4] / y[8]) * abs(dydt[8])

    # tau proteins (tau)
    dydt[5] = p.lambda_tau * (y[8] / p.N_0) + p.lambda_Gtau * (y[4] / p.G_0) \
              - p.kappa_tauFi * (y[5] ** 2) * (y[8] / p.N_0) - (y[5] / y[8]) * abs(dydt[8]) - p.d_tau * y[5]
    # # Ajout "* (y[8] / p.N_0)" au 2e terme.
    # dydt[5] = p.lambda_tau * (y[8] / p.N_0) + p.lambda_Gtau * (y[4] / p.G_0) * (y[8] / p.N_0)\
    #           - p.kappa_tauFi * (y[5] ** 2) * (y[8] / p.N_0) - (y[5] / y[8]) * abs(dydt[8]) - p.d_tau * y[5]
    # # Retire "* (y[8] / p.N_0)" au 3e terme.
    # dydt[5] = p.lambda_tau * (y[8] / p.N_0) + p.lambda_Gtau * (y[4] / p.G_0) \
    #           - p.kappa_tauFi * (y[5] ** 2) - (y[5] / y[8]) * abs(dydt[8]) - p.d_tau * y[5]

    # NFT inside the neurons (F_i)
    dydt[6] = p.kappa_tauFi * (y[5] ** 2) * (y[8] / p.N_0) - (y[6] / y[8]) * abs(dydt[8]) - p.d_Fi * y[6]
    # dydt[6] = p.kappa_tauFi * (y[5] ** 2) - (y[6] / y[8]) * abs(dydt[8]) - p.d_Fi * y[6]

    # NFT outside the neurons (F_o)
    dydt[7] = (y[6] / y[8]) * abs(dydt[8]) - p.kappa_MFo * (y[12] / (y[12] + p.K_Manti)) * y[7] - p.d_Fo * y[7]

    # Astrocytes (A)
    dydt[9] = (p.kappa_ABpoA * y[3] + p.kappa_TaA * y[17]) * (p.A_max - y[9]) - p.d_A * y[9]

    # Microglia (M)
    # dydt[10] = p.kappa_FoM * (y[7] / (y[7] + p.K_Fo)) * (p.M_max - y[10]) \
    #            + p.lambda_ABpoM * (y[3] / (y[3] + p.K_ABpoM)) * (p.M_max - y[10]) - p.d_M * y[10]
    # dydt[10] = p.kappa_FoM * (y[7] / (y[7] + p.K_Fo)) * (p.M_max - y[10]) \
    #            + p.kappa_ABooM * (y[2] / (y[2] + p.K_ABooM)) * (p.M_max - y[10]) - p.d_M * y[10]

    M_activ = p.kappa_FoM * (y[7] / (y[7] + p.K_Fo)) * y[10] + p.kappa_ABooM * (y[2] / (y[2] + p.K_ABooM)) * y[10]
    # 22-09-14_..._05 : Sans terme pour AB_o^o ; 22-09-14_..._06 : Sans terme pour F_o.

    # Resting microglia (not activated) (M_NA)
    dydt[10] = p.d_Mpro * y[11] + p.d_Manti * y[12] - M_activ

    epsilon_Ta = y[17] / (y[17] + p.K_TaAct)
    epsilon_I10 = y[16] / (y[16] + p.K_I10Act)

    # Proinflammatory microglia (M_pro)
    dydt[11] = ((p.beta * epsilon_Ta) / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ \
               - p.kappa_TbMpro * (y[15] / (y[15] + p.K_TbM)) * y[11] \
               + p.kappa_TaManti * (y[17] / (y[17] + p.K_TaM)) * y[12] \
               - p.d_Mpro * y[11]
    # 22-09-14_..._07 : Sans transfert anti -> pro
    # dydt[11] = ((p.beta * epsilon_Ta) / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ \
    #            - p.kappa_TbMpro * (y[15] / (y[15] + p.K_TbM)) * y[11] \
    #            - p.d_Mpro * y[11]

    # Anti-inflammatory microglia (M_anti)
    # Sans terme de création : 22-09-09_..._04. Tout de même saut en M_anti, hat{M}_anti, T_b et I_10, donc
    # peu concluant.
    dydt[12] = (epsilon_I10 / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ \
               + p.kappa_TbMpro * (y[15] / (y[15] + p.K_TbM)) * y[11] \
               - p.kappa_TaManti * (y[17] / (y[17] + p.K_TaM)) * y[12] - p.d_Manti * y[12]
    # 22-09-14_..._07 : Sans transfert anti -> pro
    # dydt[12] = (epsilon_I10 / (p.beta * epsilon_Ta + epsilon_I10)) * M_activ \
    #             + p.kappa_TbMpro * (y[15] / (y[15] + p.K_TbM)) * y[11] \
    #             - p.d_Manti * y[12]

    # Proinflammatory macrophages (hat{M}_pro)
    dydt[13] = p.kappa_PMhat * (y[18] / (y[18] + p.K_P)) * (p.Mhatmax - (y[13] + y[14])) * \
               ((p.beta * epsilon_Ta) / (p.beta * epsilon_Ta + epsilon_I10)) \
               - p.kappa_TbMhatpro * (y[15] / (y[15] + p.K_TbMhat)) * y[13] \
               + p.kappa_TaMhatanti * (y[17] / (y[17] + p.K_TaMhat)) * y[14] - p.d_Mhatpro * y[13]

    # Anti-inflammatory macrophages (hat{M}_anti)
    # Sans terme de création : 22-09-09_..._04. Tout de même saut en M_anti, hat{M}_anti, T_b et I_10, donc
    # peu concluant.
    dydt[14] = p.kappa_PMhat * (y[18] / (y[18] + p.K_P)) * (p.Mhatmax - (y[13] + y[14])) * \
               (epsilon_I10 / (p.beta * epsilon_Ta + epsilon_I10)) \
               + p.kappa_TbMhatpro * (y[15] / (y[15] + p.K_TbMhat)) * y[13] \
               - p.kappa_TaMhatanti * (y[17] / (y[17] + p.K_TaMhat)) * y[14] - p.d_Mhatanti * y[14]
    # dydt[14] = p.kappa_TbMhatpro * (y[15] / (y[15] + p.K_TbMhat)) * y[13] \
    #            - p.kappa_TaMhatanti * (y[17] / (y[17] + p.K_TaMhat)) * y[14] - p.d_Mhatanti * y[14]

    # TGF-Beta = Transforming growth factor beta (T_beta)
    dydt[15] = p.kappa_MantiTb * y[12] + p.kappa_MhatantiTb * y[14] - p.d_Tb * y[15]

    # IL-10 = Interleukin 10 (I_10)
    dydt[16] = p.kappa_MantiI10 * y[12] + p.kappa_MhatantiI10 * y[14] - p.d_I10 * y[16]

    # TNF-alpha = Tumor necrosis factor alpha (T_alpha)
    dydt[17] = p.kappa_MproTa * y[11] + p.kappa_MhatproTa * y[13] - p.d_Ta * y[17]

    # MCP-1 (P)
    dydt[18] = p.kappa_MproP * y[11] + p.kappa_MhatproP * y[13] + p.kappa_AP * y[9] - p.d_P * y[18]

    return dydt
