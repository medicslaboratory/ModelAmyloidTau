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
        y[2]   Amyloid-beta oligomers outside
        y[3]   Amyloid-beta plaque outside the neurons
        y[4]   GSK3
        y[5]   tau proteins
        y[6]   F_i (NFT inside the neurons)
        y[7]   F_o (NFT outside the neurons)
        y[8]   Living neurons
        y[9]   Astrocytes
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

    # Living neurons (N)
    dydt[8] = -p.d_FiN * (y[6] / (y[6] + p.K_Fi)) * y[8] \
              - p.d_TaN * (y[17] / (y[17] + p.K_Ta)) * (1 / (1 + (y[16] / p.K_I10))) * y[8]

    # Amyloid-beta monomer inside the neurons (AB^i)
    dydt[0] = p.lambda_ABi * (y[8] / p.N_0) - p.d_ABi * y[0] - (y[0] / y[8]) * abs(dydt[8])

    # Amyloid-beta monomer outside the neurons (AB_m^o)
    # TODO: Retiré "+ p.lambda_ABmo * (y[8] / p.N_0)" => _15 (vs _14)
    dydt[1] = (y[0] / y[8]) * abs(dydt[8]) - p.d_ABmo * y[1] - p.lambda_ABmoABoo * y[1] * (1 + p.AP * p.delta_AP)

    # Amyloid Beta oligomers outside (AB_o^o)
    dydt[2] = - p.d_ABoo * y[2] + (p.lambda_ABmoABoo * y[1] - p.lambda_ABooABpo * y[2]) * (1 + p.AP * p.delta_AP)

    # Amyloid-beta plaque outside the neurons (AB_p^o)
    dydt[3] = p.lambda_AABpo * (y[9] / p.A_0) + p.lambda_ABooABpo * y[2] * (1 + p.AP * p.delta_AP) \
              - (p.d_M2hatABpo * y[14] + p.d_MABpo * (y[11] + p.theta * y[12])) * (y[3] / (y[3] + p.K_ABpo))

    # print(str(p.lambda_AABpo * (y[9] / p.A_0)) + ",   " + str(p.lambda_ABooABpo * y[2] * (1 + p.AP * p.delta_AP)) + ",   " + str(- p.d_M2hatABpo * y[14] * (y[3] / (y[3] + p.K_ABpo))) + ",   " + str(- p.d_MABpo * (y[11] + p.theta * y[12]) * (y[3] / (y[3] + p.K_ABpo))))

    # Glycogen synthase kinase-type 3 (GSK-3) (G)
    dydt[4] = p.lambda_ABiG * y[0] - p.d_G * y[4]

    # tau proteins (tau)
    dydt[5] = (y[8] / p.N_0) * p.lambda_tau + p.lambda_Gtau * y[4] - p.d_tau * y[5]

    # NFT inside the neurons (F_i)
    dydt[6] = p.lambda_tauFi * y[5] * (y[8] / p.N_0) - p.d_Fi * y[6] - (y[6] / y[8]) * abs(dydt[8])

    # NFT outside the neurons (F_o)
    # Avant: = (y[6] / p.N_0) * abs(dydt[8]) - p.d_Fo * y[7]
    dydt[7] = (y[6] / y[8]) * abs(dydt[8]) - p.d_Fo * y[7]

    # Astrocytes (A)
    # Avant : p.lambda_ABpoA * y[3] / p.W_A + p.lambda_TaA * y[17] / p.W_A - p.d_A * y[9]
    #         où p.W_A = 1
    dydt[9] = p.lambda_ABpoA * y[3] + p.lambda_TaA * y[17] - p.d_A * y[9]

    # Microglia (M)
    # TODO:
    #  Dans Hao, ajout en fonction de ABO et non AB_out (ici plaque). Peut-être changer y[3] -> y[2]?
    #  Mais revoir comment arranger les unitées/termes pour que être certain que les deux premiers termes fonctionnent.
    dydt[10] = p.lambda_FoM * (y[7] / (y[7] + p.K_Fo)) + p.lambda_ABpoM * (y[3] / (y[3] + p.K_ABpoM)) \
               - p.d_M * y[10]
    # print(str(p.lambda_FoM * (y[7] / (y[7] + p.K_Fo))) + ",   " + str(p.lambda_ABpoM * (y[3] / (y[3] + p.K_ABpoM))) + ",   " + str(-p.d_M * y[10]))

    # Proinflammatory microglia (M_1)
    # TODO: Revoir modif. Retiré "- p.d_M1 * y[11]", car déjà pris en compte dans l'eq pour M
    #  Ajout terme conversion M1 -> M2
    #  Modif "+ y[10] * (p.beta / (p.beta + 1)) * p.lambda_MM1" par "+ dydt[10] * (p.beta / (p.beta + 1))" probleme
    #  potentiel si Neurons (donc F_o, donc M) diminu trop rapidement (changement fait à partir de fig _17)
    dydt[11] = dydt[10] * (p.beta / (p.beta + 1)) - p.lambda_TBM2 * (y[15] / (y[15] + p.K_TB)) * y[11]

    # Anti-inflammatory microglia (M_2)
    # TODO: Revoir modif. Retiré "- p.d_M2 * y[12]", car déjà pris en compte dans l'eq pour M
    #  Modif terme conversion M1 -> M2 : "+ p.lambda_TBM2 * y[15]" -> celui de Hao (fait + de sens)
    #  Modif "+ y[10] * (1 / (p.beta + 1)) * p.lambda_MM1" par "+ dydt[10] * (1 / (p.beta + 1))" problème potentiel si
    #   Neurons (donc F_o, donc M) diminu trop rapidement (changement fait à partir de fig _17)
    dydt[12] = dydt[10] * (1 / (p.beta + 1)) + p.lambda_TBM2 * (y[15] / (y[15] + p.K_TB)) * y[11]

    # Avec modif faites sur eqns ci-avant, rendue à graph ...ModifEqns_8

    # Proinflammatory macrophages (M_1^hat)
    # TODO: Thèse : Dit que la forme serait la même que dans Hao, mais pas vraiment. Le premier terme devrait être en
    #  relation avec M^hat que nous n'avons pas ici.
    #  dydt[13] = (y[18] / (y[18] + p.K_P)) * (p.M1hateq - y[13]) * p.lambda_PM1hat - p.d_M1hat * y[13]
    dydt[13] = (y[18] / (y[18] + p.K_P)) * (p.M1hateq - y[13]) - p.d_M1hat * y[13]

    # Anti-inflammatory macrophages (M_2^hat)
    # TODO: À discuter, voir note (OneNote)
    dydt[14] = p.lambda_TB * y[15] - p.d_M2hat * y[14]

    # TGF-Beta = Transforming growth factor beta (T_beta)
    # TODO: Pourquoi M_1 et M_1^hat? Dans Hao, il est produit par M_2 et M_2^hat, et les constantes sont les mêmes.
    #         dydt[15] = p.lambda_M1TB * y[12] + p.lambda_M1hatTB * y[14] - p.d_TB * y[15] (graph ..._10)
    # dydt[15] = p.lambda_M1TB * y[11] + p.lambda_M1hatTB * y[13] - p.d_TB * y[15]
    dydt[15] = p.lambda_M2TB * y[12] + p.lambda_M2hatTB * y[14] - p.d_TB * y[15]

    # IL-10 = Interleukin 10 (I_10)
    # Erreur : dydt[16] = p.lambda_M2I10 * y[12] / p.K_M2 - p.lambda_M2I10 * y[16]
    # TODO: Pourquoi 2 cstes plutôt qu'une?
    #  p.lambda_M2I10 * y[12] / p.K_M2 -> p.lambda_M2I10 * y[12]
    #  Grosse différence d'ordre de grandeur: (Cependant, ne change pas vrm les résultats).
    #  p.lambda_M2I10/p.K_M2 = 6.67e-3/0.017 = 0.3924 vs. Hao (1 cste) = 6.67e-3.
    dydt[16] = p.lambda_M2I10 * y[12] - p.d_I10 * y[16]

    # TNF-alpha = Tumor necrosis factor alpha (T_alpha)
    # #TODO: Pourquoi 2 cstes plutôt qu'une?
    #   p.lambda_M1hatTa/p.K_M1hat = 2.675 ; Hao (1 cste) = 3e-2
    #   p.lambda_M1Ta/p.K_M1 = 3.56667 ; Hao (1 cste) = 3e-2
    # Avant : dydt[17] = p.lambda_M1Ta * y[11] / p.K_M1 + p.lambda_M1hatTa * y[13] / p.K_M1hat - p.d_Ta * y[17]
    # (modif: graph ..._11 vs ..._9)
    dydt[17] = p.lambda_M1Ta * y[11] + p.lambda_M1hatTa * y[13] - p.d_Ta * y[17]

    # MCP-1 (P)
    dydt[18] = p.lambda_AP * y[9] - p.d_P * y[18]

    return dydt
