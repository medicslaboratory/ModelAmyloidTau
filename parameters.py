# Date : 6 septembre 2022
# Autor : Éléonore Chamberland

# This file defines a class for the parameters of the model.
# Recall : 1 cm^3 = 1 mL

import math
from scipy.constants import Avogadro  # Avogadro number


class Parameters():
    """
    Class defining the parameters of the model.
    """

    def d_ABmo(self, t):
        """
        Fonction for the degradation rate of extracellular amyloid-beta 42 monomer. The half-life is linear with age,
        with 3.8h at 30 y.o. to 9.4h at 80 y.o.
        :param t: Age of the person (in days).
        :return: The degradation rate of extracellular amyloid-beta 42 monomer (/day).
        """
        halflife = (7 / 547500) * t + (11 / 600)
        return math.log(2) / halflife

    def Ins(self, t, S):
        """
        Function for the concentration on insulin in the brain in function of age and the sex. From the data for
        peripherical insulin of Bryhni et al. (2010), with a molar mass of insulin of 5808 g/mol (Litwack, 2022).
        We suppose that brain insulin is 10 times smaller than the peripheral one (Gray et al., 2014).
        :param t: Age of the person (in days).
        :param S: Sex of the person (0 for a woman, and 1 for a man).
        :return: Insulin concentration in brain (g/mL).
        """
        if S == 0:  # woman
            return 0.1 * (-4.151e-15 * t + 3.460e-10)
        elif S == 1:  # men
            return 0.1 * (-4.257e-15 * t + 3.763e-10)

    def __init__(self):
        self.AP = 1
        """AP equals to 1 if one has the APOE4 gene and 0 otherwise."""

        self.S = 0
        """Value for the sex. S equals to 0 if the person is a woman, and 1 for a man."""

        self.rho_cerveau = 1.03
        """Density of the brain (g/mL)."""

        if self.S == 0:  # If woman
            self.N_0 = 0.45
            """Reference density of neuron in woman (g/cm^3) (= g/mL)."""
            self.A_0 = 0.10
            """Reference density of astrocytes in woman (g/cm^3) (= g/mL)."""
        elif self.S == 1:  # If man
            self.N_0 = 0.42
            """Reference density of neuron in man (g/cm^3) (= g/mL)."""
            self.A_0 = 0.12
            """Reference density of astrocytes in man (g/cm^3) (= g/mL)."""

        # self.N_0_F = 0.45
        # """Reference density of neuron in woman (g/cm^3) (= g/mL)."""
        #
        # self.N_0_M = 0.42
        # """Reference density of neuron in man (g/cm^3) (= g/mL)."""
        #
        # self.gamma_N = (self.N_0_M / self.N_0_F) - 1
        # """Constant for the differentiation of the sex on the reference density of neuron (N_0)."""
        #
        # self.N_0 = self.N_0_F * (1 + self.S * self.gamma_N)
        # """Reference density of neuron with differentiation for the sex (g/cm^3) (= g/mL).
        # Voir Herculano-Houzel. Varie beaucoup selon la région. Dans néocortex pour nous. Nous prendrons N_0 pour une
        # personne ayant environ 30 ans."""
        #
        # self.A_0_F = 0.10
        # """Reference density of astrocytes in woman (g/cm^3) (= g/mL)."""
        #
        # self.A_0_M = 0.12
        # """Reference density of astrocytes in man (g/cm^3) (= g/mL)."""
        #
        # self.gamma_A = (self.A_0_M / self.A_0_F) - 1
        # """Constant for the differentiation of the sex on the reference density of astrocytes (N_0)."""
        #
        # self.A_0 = self.A_0_F * (1 + self.S * self.gamma_A)
        # """Reference density of astrocytes with differentiation for the sex (g/cm^3) (= g/mL).
        # Même idée que N_0."""

        M_ABm = 4514
        """Molar mass of a peptide (monomer) of AB42 (g/mol)."""

        m_Mhat = 4.990e-9
        """Mass of a macrophage (or microglia) (g)."""

        ##########################################
        # CONSTANTS FOR THE EQUATION FOR NEURONS #
        ##########################################

        self.d_FiN = 1/(2.51 * 365)  # = 1.0915e-3
        # self.d_FiN = 1 / (2.51 * 365) * 1e-2
        """Maximal death rate of neurons induced by F_i (/day)."""

        self.K_Fi = 0.1 * (0.6 * (6e-3 * self.rho_cerveau))  # approx 3.708e-4
        """Concentration of intracellular NFTs (F_i) for which the death rate of neuron induced by F_i is 
        half-maximal (g/mL)."""

        self.n = 2
        """Sigmoid function coefficient (unitless)."""
        # TODO: À déterminer avec le modèle.

        self.d_TaN = 0.01 / 365  # (1 / 2) * 1/(2.51 * 365)
        # self.d_TaN = (1/2) * self.d_FiN
        # TODO: Peut-être changer
        """Maximal death rate of neurons induced by T_alpha (TNF-alpha) (/day)."""

        self.K_Ta = 4.48e-12
        """Concentration of T_alpha (TNF-alpha) for which the death rate of neuron induced by TNF-alpha is 
        half-maximal (g/mL)."""

        self.K_I10 = 2.12e-12
        """Concentration of IL-10 for which the rate of neuron death induced by TNF-alpha is divided by two (g/mL)."""

        #######################################
        # CONSTANTS FOR THE EQUATION FOR AB^i #
        #######################################

        self.lambda_ABi = (1 / 2) * ((5631e-9 - 783e-9) / (50 * 365)) * self.rho_cerveau  # approx. = 1.3681e-10
        """Creation rate of amyloid-beta42 inside neuron from APP (g/mL/day)."""

        self.delta_APi = (8373e-9 - 2178e-9) / (5631e-9 - 783e-9) - 1  # approx. 0.2778
        """Constant that quantifies the impact of the APOE4 gene on the creation rate amyloid from APP pathway inside 
        neurons (unitless). Equals to the rate of creation with the allele divided by the rate without APOE4, minus 1. 
        Here is the calculation after simplifications."""

        self.d_ABi = (math.log(2)) / (1.75 / 24)  # approx. 9.5060
        """Degradation rate of amyloid-beta42 inside neurons (/day)."""

        #########################################
        # CONSTANTS FOR THE EQUATION FOR AB_m^o #
        #########################################

        self.lambda_ABmo = self.lambda_ABi
        # Voir lambda_ABi pour tests.
        """Creation rate of Amyloid Beta42 monomer outside (without APOE allele) (g/mL/day)."""

        self.delta_APm = self.delta_APi
        """This constant quantifies the impact of the APOE4 gene on the creation rate of amyloid-beta42 monomer outside 
        neurons, i.e. on lambda_ABmo (unitless)."""

        self.lambda_AABmo = (1 / 13) * self.lambda_ABmo  # approx. 1.0524e-11
        """Creation rate of amyloid-beta42 monomer outside by astrocytes (g/mL/day)"""

        kappa_ABmoABoo_min = 38 * 1000 * (1 / (2 * M_ABm)) * 86400  # approx. 3.63669e5
        kappa_ABmoABoo_max = 38 * 1000 * (1 / M_ABm) * 86400  # approx. 7.27337e5
        self.kappa_ABmoABoo = kappa_ABmoABoo_min
        """Conversion rate of extracellular amyloid-beta monomer to extracellular amyloid-beta oligomer (mL/g/day)."""
        # TODO: In the interval kappa_ABmoABoo_min to kappa_ABmoABoo_max, à tester.

        self.delta_APmo = 2.7 - 1
        """This constant quantifies the impact of the APOE4 gene on the conversion rate of extracellular amyloid-beta 
        monomer to extracellular amyloid-beta oligomer (unitless)."""

        # self.d_ABmo = self.d_ABi
        # """Degradation rate of Amyloid Beta42 monomer outside (/day)."""
        # Devenu une fonction du temps (voir ci-haut).

        #########################################
        # CONSTANTS FOR THE EQUATION FOR AB_o^o #
        #########################################

        self.kappa_ABooABpo = (3 / 7) * 1e6 * 1000 / (2 * M_ABm)  # approx 4.7471e4
        """Conversion rate of extracellular amyloid-beta42 oligomer to plaques (mL/g/day)."""

        # self.delta_APop = 1  # (no difference)
        # """This constant quantifies the impact of the APOE4 gene on the conversion rate of extracellular
        # amyloid-beta42 oligomer to plaques."""

        self.d_ABoo = 0.3e-3 * 86400  # approx 25.92
        """Degradation rate of extracellular amyloid-beta42 oligomer (/day)."""

        #########################################
        # CONSTANTS FOR THE EQUATION FOR AB_p^o #
        #########################################

        self.d_hatMantiABpo = math.log(2) / 3  # approx 0.2310
        """Degradation rate of extracellular amyloid-beta42 plaque by anti-inflammatory macrophages 
        (hat{M}_pro) (/day)"""

        self.d_MantiABpo = math.log(2) / 0.85  # approx 0.8155
        """Degradation rate of extracellular amyloid-beta42 plaque by anti-inflammatory microglia (M_pro) (/day)."""

        self.delta_APdp = (5 / 20) - 1  # = -0.75
        """This constant quantifies the impact of the APOE4 gene on the degradation rate of Amyloid Beta42 plaque 
            outside by proinflammatory macrophages."""
        # = (Rate with APOE / self.d_MprohatABpo) - 1

        self.K_ABpo = (1.11 + 0.53) / 527.4 / 1000  # approx 3.11e-6
        """Concentration of extracellular amyloid-beta42 plaques at which the degradation rate of AB_p^o by M_anti and 
        hat{M}_anti is half maximal (Michaelis-Menten constant) (g/mL) """

        ########################################
        # CONSTANTS FOR THE EQUATION FOR GSK-3 #
        ########################################

        # self.Ins = fct
        # """Concentration of insulin (g/mL)"""
        # A function of age. See beginning.

        if self.S == 0:  # woman
            self.Ins_0 = 3.006e-11  # 0.1 * (-4.151e-15 * (365 * 30) + 3.460e-10) = 3.0054655e-11
            """Normal concentration of insulin, sex dependent (g/mL). 
            Correspond to the brain concentration at 30 years old."""
        elif self.S == 1:  # men
            self.Ins_0 = 3.296e-11  # 0.1 * (-4.257e-15 * (365 * 30) + 3.763e-10) = 3.2968585e-11
            """Normal concentration of insulin, sex dependent (g/mL). 
            Correspond to the brain concentration at 30 years old."""

        self.d_G = math.log(2) / (41 / 24)  # approx 0.4057
        # Test pour que G soit en équilibre à G_0. Sinon grande différence entre G_0 et équilibre.
        #    d_G = lambda_InsG / G_0 . Résultat "Figure_22-09-08_solve_ivp_Radau_5y_APOE+_F_10.png" Pas beau :(
        #    Okay avec BDF, donne même chose que "Figure_22-09-08_solve_ivp_BDF_5y_APOE+_F_11.png", où modif faite
        #    sur lambda_InsG.
        # self.d_G = self.lambda_InsG / (1104e-12 * 47000 * self.rho_cerveau)
        """Degradation rate of GSK-3 (/day)."""

        if self.S == 0:  # woman
            self.G_0 = 1104e-12 * 47000 * self.rho_cerveau  # approx 5.3445e-05
            """Normal concentration of GSK-3 at a normal insulin concentration (g/mL)."""
        elif self.S == 1:  # men
            self.G_0 = 310e-12 * 47000 * self.rho_cerveau  # = 1.50071e-05
            """Normal concentration of GSK-3 at a normal insulin concentration (g/mL)."""

        # self.lambda_InsG = 0.18e-6 * 47000 * 1440 * 1.077 * 0.005  # approx 0.065602
        # TODO: Test pour que G soit en équilibre à G_0. Sinon grande différence entre G_0 et équilibre.
        #    lambda_InsG = d_G * G_0 . Résultat "Figure_22-09-08_solve_ivp_Radau_1y_APOE+_F_04.png" Pas beau :(
        #    Okay avec BDF. Voir "Figure_22-09-08_solve_ivp_BDF_5y_APOE+_F_11.png" et suivantes
        #    À confirmer, déjà ajusté dans Latex.
        self.lambda_InsG = self.d_G * self.G_0  # F: ~ 2.1685e-05 ; M: ~ 6.0891e-06
        """Creation rate of GSK-3 induced by the insulin (g/mL/day)"""

        ######################################
        # CONSTANTS FOR THE EQUATION FOR tau #
        ######################################

        self.lambda_tau = 26.3e-12
        """Phosphorylation rate of tau in health by other mechanism than GSK-3 (g/ml/day)."""

        self.lambda_Gtau = ((20 / 21) - (20 / 57)) * 1e-6 / 0.5 / 1000 / 1000 * 72500  # approx 8.72e-8
        # self.lambda_Gtau = (((20 / 21) - (20 / 57)) * 1e-6 / 0.5 / 1000 / 1000 * 72500) * 1e-2  # approx 8.72e-8
        # Pourrait voir le *1e-2 comme différence de GSK-3 entre neurones sains et infecté de 1/1e5 plutôt que 1/1e3.
        # 22-09-22 : Fait varier pas mal, mais pas les neurones.
        """Creation rate of tau by GSK3 (g/mL/day)"""

        # G_0, see GSK-3 section.

        self.kappa_tauFi = (100 / 3) * 1e-6 / 19344 * 86400 * 1000  # approx 0.1489
        """Conversion rate of tau in NFT (/day)."""

        self.d_tau = math.log(2) / 5.16  # approx 0.1343
        """Degradation and un-hyperphosphorylation rate of intracellular tau proteins (/day)."""

        ######################################
        # CONSTANTS FOR THE EQUATION FOR F_i #
        ######################################

        self.d_Fi = 1.0e-2 * self.d_tau  # approx 1.343e-3
        """Degradation rate of intracellular NFT (/day)."""

        ######################################
        # CONSTANTS FOR THE EQUATION FOR F_o #
        ######################################

        self.kappa_MFo = 0.4
        """Maximal rate for the degradation of extracellular NFTs by anti-inflammatory microglia (/day)."""

        if self.S == 0:  # woman
            self.K_Manti = (1 / 4) * 3.811e-2  # = 0.0095275
            """Concentration of anti-inflammatory microglia at which the rate of degradation of extracellular NFTs by 
            these cells is half-maximal (g/mL). Sex dependent."""
        elif self.S == 1:  # man
            self.K_Manti = (1 / 4) * 3.193e-2  # = 0.0079825
            """Concentration of anti-inflammatory microglia at which the rate of degradation of extracellular NFTs by 
            these cells is half-maximal (g/mL). Sex dependent."""

        self.d_Fo = 1 / 10 * self.d_tau  # approx 1.343e-2
        """Degradation rate of extracellular NFT (/day)."""

        #################################################
        # CONSTANTS FOR THE EQUATION FOR ASTROCYTES (A) #
        #################################################

        self.A_max = self.A_0
        """Maximal density of astrocytes (g/mL)."""

        self.kappa_TaA = 0.92 / 100e-9  # = 9.2e6
        """Activation rate of astrocytes by TNF-alpha (mL/g/day)."""

        self.kappa_ABpoA = (self.kappa_TaA * 2.24e-12) / (2 * self.K_ABpo)  # approx 3.3136
        """Activation rate of astrocytes by extracellular amyloid-beta42 plaque (mL/g/day)."""

        self.d_A = 0.4
        """Death rate of astrocytes (/day)."""

        #######################################
        # CONSTANTS FOR THE EQUATION FOR M_NA #
        #######################################

        TotalMaxActivRateM = 0.2141

        self.kappa_FoM = TotalMaxActivRateM * 2 / 3  # approx 0.1427
        # TODO: 28.32 * 2? Trop grand! Latex pour autre, changer si change idée.
        """Activation rate of microglia by F_o (NFT) (/day)."""

        self.K_Fo = 11 * ((1000 * 72500) / Avogadro) * 1000  # approx 1.3243e-12
        """Concentration of extracellular NFTs at which the rate of activation of microglia by F_o 
        is half-maximal (g/mL)."""

        self.kappa_ABooM = TotalMaxActivRateM * 1 / 3  # approx 0.07137
        # TODO: 28.32 ? Trop grand! Latex pour autre, changer si change idée.
        """Activation rate of microglia by extracellular amyloid-beta42 oligomer (/day). """

        self.K_ABooM = 0.060 / 527.4 / 1000  # approx 1.1377e-7
        """Concentration of extracellular amyloid-beta42 oligomer at which the rate of activation of microglia by 
        oligomer is half-maximal (g/mL)."""

        self.d_Mpro = 7.67e-4
        """Degradation rate of proinflammatory microglia (/day)."""

        self.d_Manti = 7.67e-4
        """Degradation rate of anti-inflammatory microglia (/day)."""

        ###################################################
        # CONSTANTS FOR THE EQUATION FOR M_pro AND M_anti #
        ###################################################

        self.beta = 1
        """Proinflammatory / anti-inflammatory environnemental ratio (M_pro/M_anti ratio)"""

        self.K_TaAct = 2.24e-12  # = self.K_TaM (def after)
        """Half-saturation constant of TNF-alpha for the activation of microglia to a proinflammatory 
        polarization (g/mL)."""

        self.K_I10Act = 2.12e-12  # = self.K_I10 (def sect. neurones)
        """Half-saturation constant of TNF-alpha for the activation of microglia to a proinflammatory 
        polarization (g/mL)."""

        self.kappa_TbMpro = 4.8
        """Maximal conversion rate of proinflammatory microglia to anti-inflammatory under TGF-beta signaling (/day)."""

        self.K_TbM = 5.9e-11
        """Concentration of TGF-beta for which the conversion of Mpro to Manti is half maximal (g/mL)."""

        self.kappa_TaManti = 4.8
        """Maximal conversion rate of anti-inflammatory microglia to proinflammatory under TNF-alpha 
        signaling (/day)."""

        self.K_TaM = 2.24e-12
        """Concentration of TNF-alpha for which the conversion of Manti to Mpro is half maximal (g/mL)."""

        #############################################################
        # CONSTANTS FOR THE EQUATION FOR hat{M}_pro AND hat{M}_anti #
        #############################################################

        self.kappa_PMhat = 0.33  # * 1e-3
        """Maximal importation rate of macrophage in the brain under MCP-1 signaling (/day)."""

        self.K_P = 6.23e-10
        """Concentration of MCP-1 for which the rate of importation of macrophage is half-maximal (g/mL)."""

        self.Mhatmax = (830 * m_Mhat) / 2e-4  # = 0.0207085
        """Maximal concentration of macrophage in the brain (g/mL)."""

        self.kappa_TbMhatpro = 1 / (10 / 24)  # = 2.4
        """Maximal conversion rate of proinflammatory macrophage to anti-inflammatory under TGF-beta 
        signaling (/day)."""

        self.K_TbMhat = self.K_TbM
        """Concentration of TGF-beta for which the conversion of hat{M}_pro to hat{M}_anti is half maximal (g/mL)."""

        self.kappa_TaMhatanti = 1 / (10 / 24)  # = 2.4
        """Maximal conversion rate of anti-inflammatory macrophage to proinflammatory under TNF-alpha 
        signaling (/day)."""

        self.K_TaMhat = self.K_TaM
        # Test * 1e2: 22-09-13_..._01 (car test cette augmentation sur K_TaM). Pas amélioration des pics.
        """Concentration of TNF-alpha for which the conversion of hat{M}_anti to hat{M}_pro is half maximal (g/mL)."""

        self.d_Mhatpro = 7.67e-4
        """Death rate of proinflammatory macrophages (/day)."""

        self.d_Mhatanti = 7.67e-4
        """Death rate of anti-inflammatory macrophages (/day)."""

        #########################################
        # CONSTANTS FOR THE EQUATION FOR T_beta #
        #########################################

        kappa_MhatantiTb_min = 10 * (47e-12 / 18 * 24) / (4e6 * m_Mhat)  # approx 3.14e-8 /j
        kappa_MhatantiTb_max = 10 * (47e-12 / 18 * 24) / (2e6 * m_Mhat)  # approx 6.28e-8 /j

        self.kappa_MhatantiTb = kappa_MhatantiTb_max
        """Production rate of TGF-beta by hat{M}_pro (/day)."""
        # TODO: In the interval kappa_MhatantiTb_min to kappa_MhatantiTb_max, à tester.

        self.kappa_MantiTb = self.kappa_MhatantiTb
        """Production rate of TGF-beta by M_pro (/day)."""
        # TODO: In the interval kappa_MhatantiTb_min to kappa_MhatantiTb_max, à tester.

        self.d_Tb = math.log(2) / (3 / 1440)  # approx 332.71
        """Degradation rate of TGF-beta (/day)."""

        #######################################
        # CONSTANTS FOR THE EQUATION FOR I_10 #
        #######################################

        # self.kappa_MhatantiI10 = 47 * (52e-9 / 2) / (4e6 * m_Mhat)  # approx 6.12e-5  # Value from DeWaalMalefyt91
        self.kappa_MhatantiI10 = 660e-12 / (2e5 * m_Mhat)  # approx 6.613e-7 # Value from Mia14 (agreement with Fadok98)
        """Production rate of IL-10 by anti-inflammatory macrophages (hat{M}_anti) (/day)."""

        self.kappa_MantiI10 = self.kappa_MhatantiI10
        # Voir kappa_MhatantiI10.
        """Production rate of IL-10 by anti-inflammatory microglia (M_anti) (/day)."""

        self.d_I10 = math.log(2) / (3.556 / 24)  # approx. 4.6782 /j
        """Degradation rate of IL-10 (/day)."""

        ##########################################
        # CONSTANTS FOR THE EQUATION FOR T_alpha #
        ##########################################

        # self.kappa_MhatproTa = 15.9e-9 / (1e6 * m_Mhat)  # approx 3.186e-6  # Value from Hallsworth94
        self.kappa_MhatproTa = (1.5e-9 / 18 * 24) / (2e6 * m_Mhat)  # max Fadok98 : approx 2.004e-7
        # self.kappa_MhatproTa = (1.5e-9 / 18 * 24) / (4e6 * m_Mhat)  # min Fadok9: approx 1.002e-7
        # TODO: Modif Latex selon ce qu'on conserve (les deux options y sont).
        #   Choix de ce qu'on conserve selon test modèle.
        """Production rate of TNF-alpha by proinflammatory macrophages (hat{M}_pro) (/day)."""

        self.kappa_MproTa = self.kappa_MhatproTa
        # Voir kappa_MhatproTa.
        """Production rate of TNF-alpha by proinflammatory microglia (M_pro) (/day)."""

        self.d_Ta = math.log(2) / (18.2 / 1440)  # approx 54.84
        """Degradation rate of TNF-alpha (/day)."""

        ############################################
        # CONSTANTS FOR THE EQUATION FOR MCP-1 (P) #
        ############################################

        self.kappa_MhatproP = 11e-9 / (2e6 * m_Mhat)  # approx 1.102e-6
        """Production rate of MCP-1 by proinflammatory macrophages (hat{M}_pro) (/day)."""

        self.kappa_MproP = self.kappa_MhatproP
        """Production rate of MCP-1 by proinflammatory microglia (M_pro) (/day)."""

        kappa_AP_min = (1 / 10) * self.kappa_MhatproP  # approx 1.1e-7
        kappa_AP_max = (1 / 2) * self.kappa_MhatproP  # approx 5.5e-7
        self.kappa_AP = kappa_AP_min
        """Production rate of MCP-1 by astrocytes (/day)."""
        # TODO: In the interval kappa_AP_min to kappa_AP_max, à tester.

        self.d_P = math.log(2) / (3 / 24)  # approx 5.5452
        """Degradation rate of  MCP-1 (/day)."""
