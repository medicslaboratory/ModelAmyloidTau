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
        y[11]  M_pro (Proinflammatory microglia)
        y[12]  M_anti (anti-inflammatory microglias)
        y[13]  M_pro^hat (M_pro macrophages)
        y[14]  M_anti^hat (M_anti macrophages)
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
    dydt[1] = (y[0] / y[8]) * abs(dydt[8]) - p.d_ABmo(t) * y[1] + p.lambda_ABmo * (1 + p.AP * p.delta_APm) * \
              (y[8] / p.N_0) - p.kappa_ABmoABoo * y[1] * (1 + p.AP * p.delta_APmo) + p.lambda_AABmo * (y[9] / p.A_0)

    # Amyloid Beta oligomers outside (AB_o^o)
    dydt[2] = - p.d_ABoo * y[2] + p.kappa_ABmoABoo * y[1] * (1 + p.AP * p.delta_APmo) \
              - p.kappa_ABooABpo * y[2] * (1 + p.AP * p.delta_APop)

    # Amyloid-beta plaque outside the neurons (AB_p^o)
    # TODO (Simon) Revoir si (1 + p.AP * p.delta_AP) a rapport ici
    # Changé Manti -> Mpro et Mantihat -> Mprohat, change pas vrm figure (_37)
    dydt[3] = p.kappa_ABooABpo * y[2] * (1 + p.AP * p.delta_APop) \
              - (p.d_MproABpo * y[11] + p.d_MprohatABpo * y[13]) * (y[3] / (y[3] + p.K_ABpo))

    # Glycogen synthase kinase-type 3 (GSK-3) (G)
    dydt[4] = p.lambda_InsG * (p.Ins / p.Ins_0) - p.d_G * y[4]

    # tau proteins (tau)
    dydt[5] = p.lambda_tau * (y[8] / p.N_0) + p.lambda_Gtau * (y[4] / p.G_0) - p.d_tau * y[5]

    # NFT inside the neurons (F_i)
    dydt[6] = p.kappa_tauFi * y[5] * (y[8] / p.N_0) - p.d_Fi * y[6] - (y[6] / y[8]) * abs(dydt[8])

    # NFT outside the neurons (F_o)
    # Avant: = (y[6] / p.N_0) * abs(dydt[8]) - p.d_Fo * y[7]
    dydt[7] = (y[6] / y[8]) * abs(dydt[8]) - p.d_Fo * y[7]

    # Astrocytes (A)
    # Avant : p.lambda_ABpoA * y[3] / p.W_A + p.lambda_TaA * y[17] / p.W_A - p.d_A * y[9]
    #         où p.W_A = 1
    dydt[9] = p.kappa_ABpoA * y[3] * (p.A_max - y[9]) + p.kappa_TaA * y[17] * (p.A_max - y[9]) - p.d_A * y[9]

    # Microglia (M)
    # TODO:
    #  Dans Hao, ajout en fonction de ABO et non AB_out (ici plaque). Peut-être changer y[3] -> y[2]?
    #  Mais revoir comment arranger les unitées/termes pour que être certain que les deux premiers termes fonctionnent.
    dydt[10] = p.kappa_FoM * (y[7] / (y[7] + p.K_Fo)) * (p.M_max - y[10]) + p.lambda_ABpoM * (y[3] / (y[3] + p.K_ABpoM)) * (p.M_max - y[10]) - p.d_M * y[10]

    # Proinflammatory microglia (M_pro)
    # TODO: Revoir modif. Retiré "- p.d_Mpro * y[11]", car déjà pris en compte dans l'eq pour M
    #  Ajout terme conversion Mpro -> Manti
    #  Modif "+ y[10] * (p.beta / (p.beta + 1)) * p.lambda_MMpro" par "+ dydt[10] * (p.beta / (p.beta + 1))" probleme
    #  potentiel si Neurons (donc F_o, donc M) diminu trop rapidement (changement fait à partir de fig _17)
    dydt[11] = dydt[10] * (p.beta / (p.beta + 1)) - p.kappa_TBManti * (y[15] / (y[15] + p.K_TB)) * y[11]

    # Anti-inflammatory microglia (M_anti)
    # TODO: Revoir modif. Retiré "- p.d_Manti * y[12]", car déjà pris en compte dans l'eqn pour M
    #  Modif terme conversion Mpro -> Manti : "+ p.lambda_TBManti * y[15]" -> celui de Hao (fait + de sens)
    #  Modif "+ y[10] * (1 / (p.beta + 1)) * p.lambda_MMpro" par "+ dydt[10] * (1 / (p.beta + 1))" problème potentiel si
    #   Neurons (donc F_o, donc M) diminu trop rapidement (changement fait à partir de fig _17)
    dydt[12] = dydt[10] * (1 / (p.beta + 1)) + p.kappa_TBManti * (y[15] / (y[15] + p.K_TB)) * y[11]

    # Avec modif faites sur eqns ci-avant, rendue à graph ...ModifEqns_8

    # Proinflammatory macrophages (M_pro^hat)
    # TODO: Thèse : Dit que la forme serait la même que dans Hao, mais pas vraiment. Le premier terme devrait être en
    #  relation avec M^hat (= M_pro^hat + M_anti^hat).
    dydt[13] = p.kappa_PMprohat * (y[18] / (y[18] + p.K_P)) * (p.Mprohateq - y[13]) - p.d_Mprohat * y[13]

    # Anti-inflammatory macrophages (M_anti^hat)
    # TODO: À discuter, voir note (OneNote)
    dydt[14] = p.kappa_TB * y[15] - p.d_Mantihat * y[14]

    # TGF-Beta = Transforming growth factor beta (T_beta)
    # TODO: Pourquoi M_pro et M_pro^hat? Dans Hao, il est produit par M_anti et M_anti^hat, et les constantes sont les mêmes.
    #         dydt[15] = p.kappa_MproTB * y[12] + p.kappa_MprohatTB * y[14] - p.d_TB * y[15] (graph ..._10)
    dydt[15] = p.kappa_MproTB * y[11] + p.kappa_MprohatTB * y[13] - p.d_TB * y[15]
    # dydt[15] = p.kappa_MantiTB * y[12] + p.kappa_MantihatTB * y[14] - p.d_TB * y[15]

    # IL-10 = Interleukin 10 (I_10)
    # Erreur : dydt[16] = p.lambda_MantiI10 * y[12] / p.K_Manti - p.lambda_MantiI10 * y[16]
    # TODO: Pourquoi 2 cstes plutôt qu'une?
    #  p.kappa_MantiI10 * y[12] / p.K_Manti -> p.kappa_MantiI10 * y[12]
    #  Différence d'ordre de grandeur: (Cependant, ne change pas vrm les résultats).
    #  p.kappa_MantiI10/p.K_Manti = 6.67e-3/0.017 = 0.3924 vs. Hao (1 cste) = 6.67e-3.
    dydt[16] = p.kappa_MantiI10 * y[12] - p.d_I10 * y[16]

    # TNF-alpha = Tumor necrosis factor alpha (T_alpha)
    # #TODO: Pourquoi 2 cstes plutôt qu'une?
    #   p.lambda_MprohatTa/p.K_Mprohat = 2.675 ; Hao (1 cste) = 3e-2
    #   p.lambda_MproTa/p.K_Mpro = 3.56667 ; Hao (1 cste) = 3e-2
    # Avant : dydt[17] = p.kappa_MproTa * y[11] / p.K_Mpro + p.kappa_MprohatTa * y[13] / p.K_Mprohat - p.d_Ta * y[17]
    # (modif: graph ..._11 vs ..._9)
    dydt[17] = p.kappa_MproTa * y[11] + p.kappa_MprohatTa * y[13] - p.d_Ta * y[17]

    # MCP-1 (P)
    dydt[18] = p.kappa_AP * y[9] - p.d_P * y[18]

    return dydt
