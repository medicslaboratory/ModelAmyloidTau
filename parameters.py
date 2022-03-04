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

    def d_ABmo(self, t):
        """
        Fonction for degradation rate of extracellular amyloid-beta 42 monomer. The half-life is linear with age,
        with 3.8h at 30 y.o. to 9.4h at 80 y.o.
        :param t: Age of the person (in days)
        :return: The degradation rate of extracellular amyloid-beta 42 monomer.
        """
        halflife = t/1204500 + (197/1320)
        return math.log(2)/halflife

    def __init__(self):
        self.AP = 1
        """AP equals to 1 if one has the APOE4 gene and 0 otherwise"""

        self.N_0 = 0.14
        """Reference density of neuron (g/cm^3) (= g/mL) [Value from Hao]"""
        # Voir Herculano-Houzel. Varie beaucoup selon la région. Dans néocortex pour nous. Nous prendrons N_0 pour une
        # personne ayant environ 30 ans. Peut-être même variation avec sexe.

        self.A_0 = 0.14
        """Reference density of astrocytes (g/cm^3) (= g/mL) [Value from Hao]"""
        # Même idée que N_0 (en prportion avec N_0)

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

        self.lambda_ABi = 9.51e-6
        """Creation rate of Amyloid Beta42 inside (g/mL/day)"""
        # TODO: Pas sure de la méthode...
        #   Seyed: self.lambda_ABi = self.d_ABi * 1.0e-6  # self.ABi = 10**(-6)
        #   Je prends plutôt la valeur de Hao en attendant (9.51e-6).
        # Différence H/F et APOE.

        self.d_ABi = (math.log(2)) / (9.4 / 24)
        """Degradation rate of Amyloid Beta42 inside (/day)"""
        # Grande différence avec Hao : self.d_ABi = 9.51
        # self.d_ABi = (math.log(2)) / (self.ABihalf / 24) = 1.76
        # self.ABihalf = 9.4 : Half-life of Amyloid Beta42 inside the neurons (hour) (semble ok)
        #   (Simon) Half-life à revoir (souris vs humain).

        # Devenu une fonction du temps!
        # self.d_ABmo = self.d_ABi
        # """Degradation rate of Amyloid Beta42 monomer outside (/day)"""

        self.lambda_ABmo = self.d_ABi * 1.5e-6
        """Creation rate of Amyloid Beta42 monomer outside (without APOE allele) (g/mL/day)"""
        # TODO Pas sure de la méthode...
        #  Selon thèse : = self.d_ABmo * ABmo = self.d_ABmo * 1.5 * ABi = self.d_ABmo * 1.5 * 10**(-6)
        #                           = 9.51e-6 * 1.5 * 10e-6 = 1.4265e-11
        # Différence H/F et APOE.

        self.delta_APm = 0.25
        """This constant quantifies the impact of the APOE4 gene on the cration rate of Amyloid Beta42 
        monomer outside."""
        # = (Rate with APOE / self.lambda_ABmo) - 1

        self.kappa_ABmoABoo = 1 / 5
        """Creation rate of Amyloid Beta42 oligomer outside by Amyloid Beta42 monomer outside (/day)"""
        # TODO: Pas trop sure du pourquoi du calcul de Seyed. : = 5 * (1 / 25) * self.d_ABoo
        #  Hao : The ratio of soluble AO to total AB_out is approximately 1/25 => lambda_{A_O} = 1/25 * d_{A_O}
        #  Ce calcul d'applique pas ici.
        #  Devrait être selon le nbr de mono pour 1 oligo en moyenne. Posons 1/5 pour le moment.

        # TODO : Sayed puts a value for delta_AP
        self.delta_APmo = 0.25
        """This constant quantifies the impact of the APOE4 gene on the conversion rate of Amyloid-Beta monomer outside 
        to Amyloid-Beta oligomer outside."""
        # = (Rate with APOE / self.kappa_ABmoABoo) - 1

        self.lambda_AABmo = 1/15 * 8e-11  # (1 / 10) * 8e-10
        """Creation rate of Amyloid Beta42 monomer outside by astrocytes (g/mL/day)"""
        # Seyed (méthode de Hao): = (1 / 10) * lambda_NABpo = (1 / 10) * 8 * 10 ** (-11) = 8e-12
        # lambda_NABpo : Production rate of Amyloid Beta42 plaque outside by neuron (g/mL/day) [Hao: lambda_N = 8e-9]
        # Hao: 8e-10 g/mL/day (lambda_A)
        # Nom modifié, car changé d'équation (plaque -> monomer)

        self.kappa_ABooABpo = 25.09
        """Creation rate of Amyloid Beta42 plaque outside by Amyloid Beta oligomers outside (/day)"""
        # TODO: Revoir relations avec article de Cohen 2013
        # Testé avec 9e-7 (_5) et change pas grand chose (LSODA_4 vs _5)

        # TODO : Sayed puts a value for delta_AP
        self.delta_APop = 0.25
        """This constant quantifies the impact of the APOE4 gene on the conversion rate of Amyloid-Beta oligomer outside 
        to Amyloid-Beta plaque outside."""
        # = (Rate with APOE / self.kappa_ABooABpo) - 1

        # self.theta =  0.9
        # """Relative clearance power of amyloid-beta by M_anti compared to M_pro [Value from Hao]"""

        self.d_MprohatABpo = 4e-7  # self.theta * 1.0e-2
        """Degradation rate of Amyloid Beta42 plaque outside by proinflammatory macrophages (M_pro^hat) (/day)"""
        # TODO: Revoir
        #   Seyed : 4e-7, Dit que vient de Hao, mais
        #   Hao : estime 1e-2, mais lui a plutôt le terme : d_MantihatABpo * (M_pro^hat + theta * M_anti^hat)
        #  A pas mal d'impact, voir graph _28 (9e-3) vs _29 (1e-5). Je conserve 1e-5, les courbes sont plus douces.

        self.d_MproABpo = (1 / 80) * self.d_MprohatABpo
        """Degradation rate of Amyloid Beta42 plaque outside by proinflammatory microglias (M_pro) (/day).\n
        Microgia are 80x slower than macrophage to degrade amyloid-beta plques."""

        self.K_ABpo = 7e-3
        """Concentration of Amyloid Beta42 plaque outside at which the destruction of AB_p^o by M_pro, M_anti and M_anti^hat 
        rate is half maximal (Michaelis-Menten constant) (g/mL) \n 
        [Value of Hao] : self.K_ABpo = 10**(3) * self.ABo = 10**(3) * 7*10**(-6) = 7 * 10**(-3) """

        self.d_ABoo = (1 / 10) * self.d_ABi  # = 0.1 * 1.76
        """Degradation rate of Amyloid Beta42 oligomer outside (/day) [~Computation method of Hao]"""

        self.Ins = 1e-9
        """Concentration of insulin (supposed constant for now) (g/mL)"""

        self.lambda_InsG = 0.1  # 0.25
        """Creation rate of GSK-3 induced by the insulin (g/mL/day)"""
        # If insuline is to high => inhibition of GSK3.
        # TODO Revoir relation et valeur. Pris 0.1 pour avoir une valeur de GSK3 convenable...

        self.Ins_0 = 1e-9
        """Normal concentration of insulin (supposed constant for now) (g/mL)"""

        self.d_G = 0.408
        """Degradation rate of GSK-3 (/day)"""
        # Seyed : self.d_G = (math.log(2)) / (self.GSK3half / 24) où self.GSK3half : The half-life of GSK3 = (41 ± 4)h
        # De la même source que Seyed : (0.017 ± 0.02)/h = (0.408 ± 0.48)/day  # TODO: Valeur chez souris; pour humain?

        self.lambda_tau = 2.63e-12
        """Creation rate of tau in health (g/ml/day) [Value from Sato18]"""
        # Hao : 8.1e−11 g/ml/day

        self.d_tau = math.log(2) / 23
        """Degradation rate of tau proteins (/day) \n 
            d_tau = ln(2) / tauhalf = ln(2) / 23 = 0.03014...\n 
            tauhalf : The half-life of tau proteins in humain = 23 days [Sato18]"""
        # Hao : 0.277/day (half-live = 60h)

        self.lambda_Gtau = 0.25 / 3.1e-6
        """Creation rate of tau by GSK3 (g/mL/day)\n
        0.25 est une valeur posée par Seyed."""
        # TODO Revoir valeur

        self.G_0 = (self.lambda_InsG * (self.Ins / self.Ins_0))/self.d_G
        """Normal concentration of GSK-3 at a normal insulin concentration (g/mL). \n
        Ce calcul correspond à l'équilibre de G (la valeur de G(0))."""
        # TODO Revoir valeur
        #  2.45e-1 est la valeur approximative obtenue pour G_0 (3.1e-6 est l'ancienne cond init de GSK-3)

        self.d_Fi = 1.0e-2 * self.d_tau
        """Degradation rate of intracellular NFT (/day) [Computation method of Hao]"""

        self.d_Fo = 1/5 * self.d_tau  # 1.0e-1 * self.d_tau
        """Degradation rate of extracellular NFT (/day) [Computation method of Hao].\n
        "large" effect : (1/10), "medium" effect : (1/5), and "small" effect : (1/2)
        """
        # TODO: Revoir. J'ai modifié "1.0e-1" pour "1.0e-3". Ralentis la décroissance de F_o, donc de M. (fig _17)
        #   J'ai modifié "1.0e-3" pour "1/5". Augmente la décroissance de F_o (fig_36)

        self.kappa_tauFi = 0.6 * self.d_Fo
        """Production rate of NFT by tau (/day) [Computation method of Hao]"""
        # 60% of the hyperphosphorylated tau become NFT

        # self.W_A = 1  # 10 ** (-12)
        # """?? (g/astrocyte)"""

        self.kappa_ABpoA = 1.793
        """Creation rate of astrocytes by Amyloid Beta42 plaque outside (/day) [Value from Hao lambda_{A A_beta^o}]"""

        self.A_max = self.A_0
        """Maximal density of astrocytes (g/mL)"""
        # TODO Revoir valeur. Actuellement = A(0) de Hao.

        self.kappa_TaA = 1.54
        """Production/activation rate of astrocytes by TNF-alpha (/day) [Value from Hao lambda_{A T_alpha}]"""

        self.d_A = (math.log(2) / 600) * 0.2
        """Death rate of astrocytes (/day) [Computation method and value of Hao]"""
        # self.d_A = (math.log(2) / self.Astrocyteshalf) * (1 / 10)  # TODO: Pk *0.1 ? Retiré pour l'instant.
        # self.Astrocyteshalf : Half-life of astrocytes (day) = 600  [Valeur de Hao]
        # (math.log(2) / 600) = 0.001155 = 1.155e-3 /day
        # Hao : 1.2e-3 /day

        self.kappa_FoM = 1e-3  # 2e-2 ; changement pour 1e-3 à fig ..._26
        """Creation rate of microglias by F_o (NFT) (/day) [Value from Hao (lambda_{MF} en /day)]. 
        Devrait être en g/mL/day pour que le terme ait les bonnes unitées."""  # TODO Unitées

        # TODO: Revoir. Thèse: Valeur prise pour être < K_Fi = 3.36e-10 (car Hao "more NFT reside within neurons than
        #  outside of them"). Hao : K_{F_o} = 2.58e-11
        self.K_Fo = 1.0e-11
        """Average of extracellular NFTs (g/mL)"""

        self.M_max = 0.047
        """Maximal density of microglias (g/mL)"""
        # TODO Revoir valeur. Actuellement = M(0) de Hao.

        # TODO: Revoir. He puts : (0.015 * 0.047 - 2e-2 * 1.0e-11) / 1.0e-6 = 705 ; pk? Non...
        #  La valeur de Hao (qui était pour ABO...) = 2.3e-3 /day
        self.lambda_ABpoM = 2.3e-3  # Diminuer la valeur ne semble pas avoir bcp d'impact.
        """Creation rate of microglias by Amyloid Beta42 plaque outside (/day). 
        Devrait être en g/mL/day pour que le terme ait les bonnes unitées."""  # TODO Unitées (Si reste en /day => kappa)

        self.K_ABpoM = 1e-7
        """Concentration of Amyloid Beta42 plaque outside at which the creation rate of M by AB_p^o is half maximal 
        (g/mL) [I take the value in Hao K_{A_O}]."""

        self.d_M = 0.015
        """Degradation rate of microglias (/day) [Here he takes it equal to the death rate of M_pro and M_anti from Hao]"""

        # self.lambda_MMpro = 9.3e-3
        # """Creation rate of M_pro by microglias (/day)"""
        # # TODO: Pas certaine de la provenance de cette valeur (voir thèse...)

        # self.d_Mpro = 0.015
        # """Death rate of M_pro (proinflammatory microglia) (/day)"""

        self.beta = 10
        """Proinflammatory / anti-inflammatory microglia ratio (M_pro/M_anti ratio) [Value from Hao]"""

        self.kappa_TBManti = 6e-3
        """Production rate of M_anti by TGF-beta (/day) [Value from Hao lambda_{M_pro T_beta}]\n
            [The rate (maximal) by which TGF-beta affects the change of phenotype from Mpro to Manti]"""

        self.K_TB = 2.5e-7
        """Half-saturation of TGF-beta (g/mL) [Value from Hao K_{T_beta}]\n
            [Concentration of TGF-beta for which the convertion of Mpro to Manti is half maximal]"""

        # self.d_Manti = 0.015
        # """Death rate of M_anti (anti-inflammatory microglia) (/day)"""

        self.K_P = 5e-9
        """Half-saturation of MCP-1 (g/mL) \n
            [Thèse : "We consider the value for MCP-1 saturation for influx of macrophages as K_P";
            Value in Hao: K_P = 6e-9]"""

        self.Mprohateq = 5e-2  # 2.5e-2  # 8.64e-7
        """Concentration of M_pro^hat at equilibrium (g/mL)"""
        # TODO: On devrait ici avoir la concentration de Proinflammatory macrophage dans le sang.
        #  J'ai pris la valeur de M_0 de Hao (concentration de monocyte dans le sang) divisé par 2 (supposant que
        #   1/2 M_pro^hat et 1/2 M_anti^hat dans le sang) = (5e-2)/2 = 2.5e-2. => fig _32 vs retour à 5e-2 => _33,
        #   mieux pour Astrocytes (A) et MCP-1 (P) (augmentation plutôt que diminution).

        self.d_Mprohat = 0.015
        """Death rate of M_pro^hat macrophages (/day) [Value from Hao]"""

        self.kappa_PMprohat = (0.04 * self.d_Mprohat) / 5e-9
        """Production rate of M_pro^hat by MCP-1 (/day)"""
        # self.lambda_PMprohat = (self.Mprohat * self.d_Mprohat) / self.K_P
        # self.Mprohat : M_prohat (g/ml) = 0.04  #TODO: pk?
        # self.K_P : Half-saturation of MCP-1 (g/ml) = 5*10**(-9)  #TODO: pk? (Hao: 6e-9)

        self.d_Mantihat = 0.015
        """Death rate of M_anti^hat macrophages (/day) [Value from Hao]"""

        self.d_TB = 3.33e2  # 3.33 * 10 ** (2)
        """Degradation rate of TGF-beta (/day) """

        self.kappa_TB = self.d_TB * 2.5e-7
        """Production rate of T_beta (/day) \n
            self.lambda_TB = self.d_TB * self.K_TB \n
            self.K_TB : half-saturation of TGF-beta (T_beta) = 2.5*10**(-7) [Hao]"""
        # TODO: Calcul fait pas de sens et unitées seraient g/mL/day.

        self.kappa_MproTB = 1.5e-2
        """Production rate of TGF-beta by M_pro (/day) [Value from Hao (lambda_{T_{beta} M})]"""
        # self.kappa_MantiTB = 1.5e-2
        # """Production rate of TGF-beta by M_anti (/day) [Value from Hao (lambda_{T_{beta} M})]"""

        self.kappa_MprohatTB = 1.5e-2
        """Production rate of TGF-beta by M_pro^hat (/day) [Value from Hao (lambda_{T_{beta} M^{hat}})]"""
        # self.kappa_MantihatTB = 1.5e-2
        # """Production rate of TGF-beta by M_anti^hat (/day) [Value from Hao (lambda_{T_{beta} M^{hat}})]"""

        # TODO Vérif.. Thèse: (1.2 ± 0.16) × 10^−12 g/mL/day ??
        self.kappa_MantiI10 = 6.67e-3
        """Production rate of IL-10 by M_anti (/day) [Value from Hao (lambda_{I_{10} M_anti})]"""

        # self.K_Manti = 0.017
        # """Half-saturation of Manti (g/ml) [Value from Hao]"""

        self.d_I10 = 8.32
        """Degradation rate of IL-10 (/day)"""
        # TODO: Hao donne 16.64 /day
        # half-life est 60 minutes

        # TODO Vérif.. These: (3.09 ± 1.8) × 10^−12 g/ml/day. Unités??
        #   Seyed prend : 1.07e-1 /day
        #   Hao : lambda_{T_{alpha} M_pro} = 3e-2 /day
        self.kappa_MprohatTa = 3e-2
        """Production rate of TNF-alpha by M_pro^hat (/day)"""

        # self.K_Mpro = 0.03
        # """Half-saturation of Mpro (g/ml) [Value from Hao]"""

        self.kappa_MproTa = self.kappa_MprohatTa
        """Production rate of TNF-alpha by M_pro (/day)"""
        # Hao : lambda_{T_{alpha} M_pro^hat} = 3e-2 /day (ie même val que lambda_MprohatTa)

        # self.K_Mprohat = 0.04
        # """Half-saturation of M_pro^hat (g/ml) [Value from Hao]"""

        self.d_Ta = 55.45
        """Degradation rate of TNF-alpha (/day) [Value from Hao]"""

        self.d_P = 1.73
        """Degradation rate of  MCP-1 (/day) [Value from Hao]"""

        self.kappa_AP = 6.6e-8
        """Creation rate of MCP-1 by astrocytes (/day)"""
        # Code : self.lambda_AP = (self.K_P*self.d_P)/self.A_0
        # Problème : thèse :
        #       lambda_AP = (d_P * P) / A_0
        #                 = (1.73/day * 3e-10g/ml) / 0.14g/cm^3 = 3e-9 /day
        # TODO: Revoir. En attendant je prends la valeur de Hao

