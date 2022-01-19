# Date : 17 January 2022
# Autor : Éléonore Chamberland

# This file contains the equations of the model.
# It defines a function that generates the ODE system of the model.

# import math
import numpy as np
import parameters as param


def ODEsystem(t, y):
    """
    Function that defines the equations of the model.
    :param t:
    :param y:
        y[0]   AB^i (Amyloid-beta monomer inside the neurons)
        y[1]   Amyloid-beta monomer outside the neurons
        y[2]   Amyloid-beta plaque outside the neurons
        y[3]   GSK3
        y[4]   tau proteins
        y[5]   F_i (NFT inside the neurons)
        y[6]   F_o (NFT outside the neurons)
        y[7]   Living neurons
        y[8]   Astrocytes
        y[9]   Amyloid-beta oligomers outside
        y[10]  M (Microglia)
        y[11]  M_1 (Proinflammatory microglia)
        y[12]  M_2 (anti-inflammatory microglias)
        y[13]  M_1^hat (M_1 macrophages)
        y[14]  M_2^hat (M_2 macrophages)
        y[15]  T_{beta} (TGF-beta)
        y[16]  I_10 (IL-10 = Interleukin 10)
        y[17]  T_{alpha} (TNF-alpha)
        y[18]  P (MCP-1)
    """
    p = param.Parameters()

    dydt = np.zeros(19)

    # Amyloid-beta monomer inside the neurons (AB^i)
    # TODO: Pas certaine du dernier terme.
    #  Si on se fie à Hao, au deuxième terme devrait aussi avoir "*(y[7] / p.N_0)"?
    dydt[0] = p.lambda_ABi * (y[7] / p.N_0) - p.d_ABi * y[0] - (y[0] / y[7]) * abs(dydt[7])

    # Amyloid-beta monomer outside the neurons (AB_m^o)
    # TODO: Pas certaine du premier terme (modif en fonction du dernier terme précécent).
    dydt[1] = (y[0] / y[7]) * abs(dydt[7]) - p.d_ABmo * y[1] + p.lambda_ABmo * (y[7] / p.N_0)

    # Amyloid-beta plaque outside the neurons (AB_p^o)
    dydt[2] = (p.lambda_AABpo * (y[8] / p.A_0) + p.lambda_ABooABpo * y[9] * (1 + p.AP * p.delta_AP)
               - (p.d_M2hatABpo * y[14] + p.d_MABpo * (y[11] + p.theta * y[12])) * (y[2] / (y[2] + p.K_ABpo)))

    # Glycogen synthase kinase-type 3 (GSK-3) (G)
    dydt[3] = p.lambda_ABiG * y[0] - p.d_G * y[3]

    # tau proteins (tau)
    dydt[4] = (y[7] / p.N_0) * (p.lambda_tau + p.lambda_Gtau * y[3] - p.d_tau * y[4])

    # NFT inside the neurons (F_i)
    dydt[5] = p.lambda_tauFi * y[4] * (y[7] / p.N_0) - p.d_Fi * y[5] - (y[5] / y[7]) * abs(dydt[7])

    # NFT outside the neurons (F_o)
    # Avant: = (y[5] / p.N_0) * abs(dydt[7]) - p.d_Fo * y[6] # TODO: Revoir pas clair dans thèse.
    dydt[6] = (y[5] / y[7]) * abs(dydt[7]) - p.d_Fo * y[6]

    # Living neurons (N)
    dydt[7] = (-p.d_FiN * (y[5] / (y[5] + p.K_Fi)) * y[7]
               - p.d_TaN * (y[17] / (y[17] + p.K_Ta)) * (1 / (1 + (y[16] / p.K_I10))) * y[7])

    # Astrocytes (A)
    dydt[8] = p.lambda_ABpoA * y[2] / p.W_A + p.lambda_TaA * y[17] / p.W_A - p.d_A * y[8]

    # Amyloid Beta Oligomers outside (AB_o^o)
    dydt[9] = -p.d_ABoo * y[9] + (p.lambda_ABmoABoo * y[1] - p.lambda_ABooABpo * y[9]) * (1 + p.AP * p.delta_AP)

    # Microglia (M)
    dydt[10] = (p.lambda_FoM * (y[6] / (y[6] + p.K_Fo)) + p.lambda_ABpoM * (y[2] / (y[2] + p.K_ABpo))
                - p.d_M * y[10] - p.lambda_MM1 * y[10])

    # Proinflammatory microglia (M_1)
    dydt[11] = y[10] * (p.beta / (p.beta + 1)) * p.lambda_MM1 - p.d_M1 * y[11]

    # Anti-inflammatory microglia (M_2)
    dydt[12] = y[10] * (1 / (p.beta + 1)) * p.lambda_MM1 + p.lambda_TBM2 * y[15] - p.d_M2 * y[12]

    # Proinflammatory macrophages (M_1^hat)
    dydt[13] = (y[18] / (y[18] + p.K_P)) * (p.M1hateq - y[13]) * p.lambda_PM1hat - p.d_M1hat * y[13]

    # Anti-inflammatory macrophages (M_2^hat)
    dydt[14] = p.lambda_TB * y[15] - p.d_M2hat * y[14]

    # TGF-Beta = Transforming growth factor beta (T_beta)
    dydt[15] = p.lambda_M1TB * y[11] + p.lambda_M1hatTB * y[13] - p.d_TB * y[15]

    # IL-10 = Interleukin 10 (I_10)
    # Erreur : dydt[16] = p.lambda_M2I10 * y[12] / p.K_M2 - p.lambda_M2I10 * y[16]
    # TODO: Pourquoi 2 cstes plutôt qu'une?
    #   p.lambda_M2I10/p.K_M2 = 0.3924 ; Hao (1 cste) = 6.67e-3
    dydt[16] = p.lambda_M2I10 * y[12] / p.K_M2 - p.d_I10 * y[16]

    # TNF-alpha = Tumor necrosis factor alpha (T_alpha)
    # #TODO: Pourquoi 2 cstes plutôt qu'une?
    #   p.lambda_M1hatTa/p.K_M1hat = 2.675 ; Hao (1 cste) = 3e-2
    #   p.lambda_M1Ta/p.K_M1 = 3,56667 ; Hao (1 cste) = 3e-2
    dydt[17] = p.lambda_M1Ta * y[11] / p.K_M1 + p.lambda_M1hatTa * y[13] / p.K_M1hat - p.d_Ta * y[17]

    # MCP-1 (P)
    dydt[18] = p.lambda_AP * y[8] - p.d_P * y[18]

    return dydt
