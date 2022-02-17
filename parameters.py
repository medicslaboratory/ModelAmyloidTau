# Date : 17 January 2022
# Autor : Éléonore Chamberland

# This file contains the equations of the model
# It generates a class for the parameters of the model.
# Recall : 1 cm^3 = 1 mL (for water)

import math


class Parameters():
    """
    Class defining the parameters of the model.
    """

    def __init__(self):
        self.N_0 = 0.14
        """Reference density of neuron (g/cm^3) (= g/mL) [Value from Hao]"""

        # TODO: Grande différence avec Hao : self.d_ABi = 9.51
        # self.d_ABi = (math.log(2)) / (self.ABihalf / 24) = 1.76
        # self.ABihalf = 9.4 : Half-life of Amyloid Beta42 inside the neurons (hour) (semble ok)
        self.d_ABi = (math.log(2)) / (9.4 / 24)
        """Degradation rate of Amyloid Beta42 inside (/day)"""

        # TODO: Pas sure de la méthode...
        #   Seyed: self.lambda_ABi = self.d_ABi * 1.0e-6  # self.ABi = 10**(-6)
        #   Je prends plutôt la valeur de Hao en attendant (9.51e-6).
        self.lambda_ABi = 9.51e-6
        """Creation rate of Amyloid Beta42 inside (g/mL/day)"""

        self.d_ABmo = self.d_ABi
        """Degradation rate of Amyloid Beta42 monomer outside (/day)"""

        # TODO Pas sure de la méthode...
        #  Selon thèse : = self.d_ABmo * ABmo = self.d_ABmo * 1.5 * ABi = self.d_ABmo * 1.5 * 10**(-6)
        #                           = 9.51e-6 * 1.5 * 10e-6 = 1.4265e-11
        self.lambda_ABmo = self.d_ABmo * 1.5e-5  # self.ABmo = 1.5 * 10**(-6)
        """Creation rate of  Amyloid Beta42 monomer outside (g/mL/day)"""

        self.A_0 = 0.14
        """Reference density of astrocytes (g/cm^3) (= g/mL) [Value from Hao]"""

        self.lambda_AABpo = 8e-11  # (1 / 10) * 8e-10
        """Creation rate of Amyloid Beta42 plaque outside by astrocytes (g/mL/day)"""
        # TODO : Seyed (méthode de Hao): = (1 / 10) * self.lambda_NABpo = (1 / 10) * 8 * 10 ** (-11) = 8e-12
        #   self.lambda_NABpo : Production rate of Amyloid Beta42 plaque outside by neuron (g/mL/day) [Hao: lambda_N = 8e-9]
        #   Hao: 8e-10 g/mL/day (lambda_A)
        #   8e-12 => LSODA_2 vs 8e-10 => LSODA_3 : bonne différence! En attendant: entre-deux : 8e-11 (LSODA_4)
        #   Grande variation plus vraie avec modif finale du modèle (si prend autre valeur, fait faire des variations
        #   bizarres à la fin de Ab^i et GSK3).

        self.lambda_ABooABpo = 25.09
        """Creation rate of Amyloid Beta42 plaque outside by Amyloid Beta oligomers outside (/day)"""
        # TODO: Pas clair d'où vient la valeur...
        #   Devrait être en /day (j'ai fait la modif g/ml/day -> /day)
        #   Testé avec 9e-7 (_5) et change pas grand chose (LSODA_4 vs _5)

        # TODO : Sayed put a value for AP and delta_AP
        self.AP = 1
        """AP equals to 1 if one has the APOE4 gene and 0 otherwise"""

        self.delta_AP = 0.25
        """This constant quantifies the impact of the APOE4 gene."""

        self.theta = 0.9
        """Relative clearance power of amyloid-beta by M_2 compared to M_1 [Value from Hao]"""

        self.d_M2hatABpo = 1e-5  # self.theta * 1.0e-2
        """Degradation rate of Amyloid Beta42 plaque outside by M_2^hat (/day)"""
        # TODO: Revoir
        #   Seyed : 4e-7
        #   Hao : estime 10^-2, mais lui a plutôt le terme : d_M2hatABpo * (M_1^hat + theta * M_2^hat)
        #   Prenons: theta * valeur de Hao = 0.9 * 10^-2 (= 9e-3)
        #  A pas mal d'impact, voir graph _28 (9e-3) vs _29 (1e-5). Je conserve 1e-5, les courbes sont plus douces.

        self.d_MABpo = (1 / 5) * (self.d_M2hatABpo / self.theta)
        """Degradation rate of Amyloid Beta42 plaque outside by microglias (M) (/day) [Computation from Hao : 
        (1 / 5) * self.d_M2hatABpo ; ajout de "/ self.theta" pour compenser pour la modification de self.d_M2hatABpo]"""

        self.K_ABpo = 7e-3
        """Concentration of Amyloid Beta42 plaque outside at which the destruction of AB_p^o by M_1, M_2 and M_2^hat 
        rate is half maximal (Michaelis-Menten constant) (g/mL) \n 
        [Value of Hao] : self.K_ABpo = 10**(3) * self.ABo = 10**(3) * 7*10**(-6) = 7 * 10**(-3) """

        self.d_ABoo = (1 / 10) * self.d_ABmo  # = 0.1 * 1.76
        """Degradation rate of Amyloid Beta42 oligomer outside (/day) [Computation method of Hao]"""

        self.lambda_ABmoABoo = 1/5
        """Creation rate of Amyloid Beta42 oligomer outside by Amyloid Beta42 monomer outside (/day)"""
        # TODO: Pas trop sure du pourquoi du calcul de Seyed. : = 5 * (1 / 25) * self.d_ABoo
        #  Hao : The ratio of soluble AO to total AB_out is approximately 1/25 => lambda_{A_O} = 1/25 * d_{A_O}
        #  Ce calcul d'applique pas ici.
        #  Devrait être selon le nbr de mono pour 1 oligo en moyenne. Posons 1/5 pour le moment.

        self.lambda_ABiG = 0.25
        """Creation rate of GSK-3 by Amyloid Beta42 inside neurons (/day)"""
        # TODO : Revoir valeur.

        self.d_G = 0.408
        """Degradation rate of GSK-3 (/day)"""
        # Seyed : self.d_G = (math.log(2)) / (self.GSK3half / 24) où self.GSK3half : The half-life of GSK3 = (41 ± 4)h
        # De la même source que Seyed : (0.017 ± 0.02)/h = (0.408 ± 0.48)/day  # TODO: Valeur chez souris; pour humain?

        self.lambda_tau = 2.63e-11  # 26.3 * 10 ** (-12)
        """Creation rate of tau in health (g/ml/day)"""
        # Hao : 8.1 × 10^−11 g/ml/day

        self.d_tau = (math.log(2)) / 23  # = 0.03014...
        """Degradation rate of tau proteins (/day) \n 
            self.d_tau = (math.log(2)) / self.Tauhalf \n 
            self.Tauhalf : The half-life of tau proteins in humain = 23 days"""
        # Hao : 0.277/day (half-live = 60h)

        # TODO: Ceci est une valeur posée par Seyed
        self.lambda_Gtau = 0.25
        """Creation rate of tau by GSK3 (/day)"""

        self.d_Fi = 1.0e-2 * self.d_tau
        """Degradation rate of intracellular NFT (/day) [Computation method of Hao]"""

        self.d_Fo = 1.0e-3 * self.d_tau  # 1.0e-1 * self.d_tau
        """Degradation rate of extracellular NFT (/day) [Computation method of Hao]"""
        # TODO: Revoir. J'ai modifié "1.0e-1" pour "1.0e-3". Ralentis la décroissance de F_o, donc de M. (fig _17)

        self.lambda_tauFi = 0.6 * self.d_Fo
        """Production rate of NFT by tau (/day) [Computation method of Hao]"""
        # 60% of the hyperphosphorylated tau become NFT

        self.d_FiN = 3e-5  # ((4 + 4 * 1) / (3 + 2 * 1)) * (math.log(2) / 3650)
        """Degradation rate of neurons by F_i (/day)"""
        # self.d_FiN = ((4 + 4 * self.gamma) / (3 + 2 * self.gamma)) * self.d_N  = 3.038e-4  [Hao: 3.4e-4]
        # self.gamma : I_10 inhibition ratio = 1 (Hao)
        # self.d_N : the death rate of neuron (/day) = ln2/10ans = ln2/3650
        # TODO: Revoir d_N (Seyed semble pas sur)
        #  J'ai modifié pour 3e-5, fait bien ralentir décroissance de neurones et fait plus de sens. (fig _17)

        self.K_Fi = 0.7 * 490e-12  # 0.7 * 490 * 10 ** (-12) = 3.42999e-10
        """Half-saturation of intracellular NFTs (g/mL)"""
        # These: Assuming that in AD, 70% of hyperphosphorylated tau proteins (whose concentration in disease is
        # 490 pg/ml) are in NFT form.
        # self.K_Fi = 0.7 * self.htau = 0.7 * 490 * 10 ** (-12) = 3.43e-10   [Hao : 3.36e-10]
        # self.htau : concentration of hyperphosphorylated tau in disease (g/ml) = 490 * 10**(-12)

        self.d_TaN = (1 / 2) * self.d_FiN
        """Degradation rate of neurons by T_alpha (TNF-alpha) (/day) [Computation method of Hao]"""

        self.K_Ta = 4e-5
        """Half-saturation of T_alpha (TNF-alpha) (g/mL)"""
        # TODO: Seyed a mis 2.5e-5 et dit que c'est la valeur de Hao16, mais est plutôt de 4e-5 g/mL ?

        self.K_I10 = 2.5e-6  # 2.5 * 10 ** (-6)
        """Half-saturation of IL-10 (g/mL) [Value of Hao]"""

        # self.W_A = 1  # 10 ** (-12)
        # """?? (g/astrocyte)"""  #TODO ?? Retiré de l'eq des astrocytes.

        self.lambda_ABpoA = 1.793
        """Creation rate of astrocytes by Amyloid Beta42 plaque outside (/day) [Value from Hao lambda_{A A_beta^o}]"""

        self.lambda_TaA = 1.54
        """Production/activation rate of astrocytes by TNF-alpha (/day) [Value from Hao lambda_{A T_alpha}]"""

        self.d_A = (math.log(2) / 600) * 0.2
        """Death rate of astrocytes (/day) [Computation method and value of Hao]"""
        # self.d_A = (math.log(2) / self.Astrocyteshalf) * (1 / 10)  # TODO: Pk *0.1 ? Retiré pour l'instant.
        # self.Astrocyteshalf : Half-life of astrocytes (day) = 600  [Valeur de Hao]
        # (math.log(2) / 600) = 0.001155 = 1.155e-3 /day
        # Hao : 1.2e-3 /day

        self.lambda_FoM = 1e-3  # 2e-2 ; changement pour 1e-3 à fig ..._26
        """Creation rate of microglias by F_o (NFT) (/day) [Value from Hao (lambda_{MF} en /day)]. 
        Devrait être en g/mL/day pour que le terme ait les bonnes unitées."""  # TODO Unitées

        # TODO: Revoir. Thèse: Valeur prise pour être < K_Fi = 3.36e-10 (car Hao "more NFT reside within neurons than
        #  outside of them"). Hao : K_{F_o} = 2.58e-11
        self.K_Fo = 1.0e-11
        """Average of extracellular NFTs (g/mL)"""

        # TODO: Revoir. He puts : (0.015 * 0.047 - 2e-2 * 1.0e-11) / 1.0e-6 = 705 ; pk? Non...
        #  La valeur de Hao (qui était pour ABO...) = 2.3e-3 /day
        self.lambda_ABpoM = 2.3e-3  # Diminuer la valeur ne semble pas avoir bcp d'impact.
        """Creation rate of microglias by Amyloid Beta42 plaque outside (/day). 
        Devrait être en g/mL/day pour que le terme ait les bonnes unitées."""  # TODO Unitées

        self.K_ABpoM = 1e-7
        """Concentration of Amyloid Beta42 plaque outside at which the creation rate of M by AB_p^o is half maximal 
        (g/mL) [I take the value in Hao K_{A_O}]."""

        self.d_M = 0.015
        """Degradation rate of microglias (/day) [Here he takes it equal to the death rate of M_1 and M_2 from Hao]"""

        # self.lambda_MM1 = 9.3e-3
        # """Creation rate of M_1 by microglias (/day)"""
        # # TODO: Pas certaine de la provenance de cette valeur (voir thèse...)

        # self.d_M1 = 0.015
        # """Death rate of M_1 (proinflammatory microglia) (/day)"""

        self.beta = 10
        """Proinflammatory / anti-inflammatory microglia ratio (M_1/M_2 ratio) [Value from Hao]"""

        self.lambda_TBM2 = 6e-3
        """Production rate of M_2 by TGF-beta (/day) [Value from Hao lambda_{M_1 T_beta}]\n
            [The rate (maximal) by which TGF-beta affects the change of phenotype from M1 to M2]"""

        self.K_TB = 2.5e-7
        """Half-saturation of TGF-beta (g/mL) [Value from Hao K_{T_beta}]\n
            [Concentration of TGF-beta for which the convertion of M1 to M2 is half maximal]"""

        # self.d_M2 = 0.015
        # """Death rate of M_2 (anti-inflammatory microglia) (/day)"""

        self.K_P = 5e-9
        """Half-saturation of MCP-1 (g/mL) \n
            [Thèse : "We consider the value for MCP-1 saturation for influx of macrophages as K_P";
            Value in Hao: K_P = 6e-9]"""

        self.M1hateq = 2.5e-2  # 8.64e-7
        """Concentration of M_1^hat at equilibrium (g/mL)"""
        # TODO: On devrait ici avoir la concentration de Proinflammatory macrophage dans le sang.
        #  J'ai pris la valeur de M_0 de Hao (concentration de monocyte dans le sang) divisé par 2 (supposant que
        #   1/2 M_1^hat et 1/2 M_2^hat dans le sang) = (5e-2)/2 = 2.5e-2.

        self.d_M1hat = 0.015
        """Death rate of M_1^hat macrophages (/day) [Value from Hao]"""

        # self.lambda_PM1hat = (0.04 * self.d_M1hat) / 5e-9
        # """Production rate of M_1^hat by MCP-1 (/day)"""
        # # self.lambda_PM1hat = (self.M1hat * self.d_M1hat) / self.K_P
        # # self.M1hat : M_1hat (g/ml) = 0.04  #TODO: pk?
        # # self.K_P : Half-saturation of MCP-1 (g/ml) = 5*10**(-9)  #TODO: pk? (Hao: 6e-9)

        self.d_M2hat = 0.015
        """Death rate of M_2^hat macrophages (/day) [Value from Hao]"""

        self.d_TB = 3.33e2  # 3.33 * 10 ** (2)
        """Degradation rate of TGF-beta (/day) """

        self.lambda_TB = self.d_TB * 2.5e-7
        """Production rate of T_beta (/day) \n
            self.lambda_TB = self.d_TB * self.K_TB \n
            self.K_TB : half-saturation of TGF-beta (T_beta) = 2.5*10**(-7) [Hao]"""
        # TODO: Calcul fait pas de sens et unitées seraient g/mL/day.

        # self.lambda_M1TB = 1.5e-2
        # """Production rate of TGF-beta by M_1 (/day) [Value from Hao (lambda_{T_{beta} M})]"""
        self.lambda_M2TB = 1.5e-2
        """Production rate of TGF-beta by M_2 (/day) [Value from Hao (lambda_{T_{beta} M})]"""

        # self.lambda_M1hatTB = 1.5e-2
        # """Production rate of TGF-beta by M_1^hat (/day) [Value from Hao (lambda_{T_{beta} M^{hat}})]"""
        self.lambda_M2hatTB = 1.5e-2
        """Production rate of TGF-beta by M_2^hat (/day) [Value from Hao (lambda_{T_{beta} M^{hat}})]"""

        # TODO Vérif.. Thèse: (1.2 ± 0.16) × 10^−12 g/mL/day ??
        self.lambda_M2I10 = 6.67e-3
        """Production rate of IL-10 by M_2 (/day) [Value from Hao (lambda_{I_{10} M_2})]"""

        # self.K_M2 = 0.017
        # """Half-saturation of M2 (g/ml) [Value from Hao]"""

        self.d_I10 = 8.32
        """Degradation rate of IL-10 (/day)"""
        # TODO: Hao donne 16.64 /day
        # half-life est 60 minutes

        # TODO Vérif.. These: (3.09 ± 1.8) × 10^−12 g/ml/day. Unités??
        #   Seyed prend : 1.07e-1 /day
        #   Hao : lambda_{T_{alpha} M_1} = 3e-2 /day
        self.lambda_M1hatTa = 3e-2
        """Production rate of TNF-alpha by M_1^hat (/day)"""

        # self.K_M1 = 0.03
        # """Half-saturation of M1 (g/ml) [Value from Hao]"""

        self.lambda_M1Ta = self.lambda_M1hatTa
        """Production rate of TNF-alpha by M_1 (/day)"""
        # Hao : lambda_{T_{alpha} M_1^hat} = 3e-2 /day (ie même val que lambda_M1hatTa)

        # self.K_M1hat = 0.04
        # """Half-saturation of M_1^hat (g/ml) [Value from Hao]"""

        self.d_Ta = 55.45
        """Degradation rate of TNF-alpha (/day) [Value from Hao]"""

        self.d_P = 1.73
        """Degradation rate of  MCP-1 (/day) [Value from Hao]"""

        self.lambda_AP = 6.6e-8
        """Creation rate of MCP-1 by astrocytes (/day)"""
        # Code : self.lambda_AP = (self.K_P*self.d_P)/self.A_0
        # Problème : thèse :
        #       lambda_AP = (d_P * P) / A_0
        #                 = (1.73/day * 3e-10g/ml) / 0.14g/cm^3 = 3e-9 /day
        # TODO: Revoir. En attendant je prends la valeur de Hao

