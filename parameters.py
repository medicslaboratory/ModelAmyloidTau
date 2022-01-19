# Date : 17 January 2022
# Autor : Éléonore Chamberland

# This file contains the equations of the model
# It generates a class for the parameters of the model.

import math


class Parameters():
    """
    Class defining the parameters of the model.
    """

    def __init__(self):
        # Reference density of neuron (neurons/cm^3)
        self.N_0 = 7e7  # 7 * 10 ** (7)

        # Half-life of Amyloid Beta42 inside the neurons (day)
        self.ABihalf = 9.4

        # Degradation rate of Amyloid Beta42 inside (day^-1)
        # self.d_ABi = 9.51  # [Hao]
        self.d_ABi = (math.log(2)) / (self.ABihalf / 24)  # = 1.76

        # Creation rate of Amyloid Beta42 inside (g/ml/day) #TODO Fait pas de sens...
        # self.lambda_ABi = self.d_ABi * 1.0e-6  # self.ABi = 10**(-6)
        self.lambda_ABi = 9.51e-6  # [Hao]

        # Degradation rate of Amyloid Beta42 monomer outside (day^-1)
        self.d_ABmo = self.d_ABi

        # Creation rate of  Amyloid Beta42 monomer outside (g/ml/day) #TODO Fait pas de sens...
        # Selon thèse : = self.d_ABmo * ABmo = self.d_ABmo * 1.5 * ABi = self.d_ABmo * 1.5 * 10**(-6)
        self.lambda_ABmo = self.d_ABmo * 1.5e-6  # self.ABmo = 1.5 * 10**(-6)

        # Reference density of astrocytes (astrocytes/cm^3)
        self.A_0 = 7e7  # 7 * 10 ** 7

        # Creation rate of Amyloid Beta42 plaque outside by astrocytes
        # = (1 / 10) * self.lambda_NABpo
        # self.lambda_NABpo : the production rate of Amyloid Beta42 plaque outside by neuron (g/ml/day)
        self.lambda_AABpo = (1 / 10) * 8e-11  # (1 / 10) * 8 * 10 ** (-11)

        # Creation rate of Amyloid Beta42 plaque outside by Amyloid Beta oligomers outside (g/ml/day)
        self.lambda_ABooABpo = 25.09

        # Put a value for AP and delta_AP #TODO
        # AP stands for 1 if one has the APOE4 gene and 0 otherwise
        self.AP = 1
        # This constant quantifies the impact of the APOE4 gene.
        self.delta_AP = 0.25

        # Degradation rate of Amyloid Beta42 plaque outside by M_2^hat (/day)
        self.d_M2hatABpo = 4e-7  # 4 * 10 ** (-7)

        # Degradation rate of Amyloid Beta42 plaque outside by microglias (M) (/day)
        self.d_MABpo = (1 / 5) * self.d_M2hatABpo

        # Relative clearance power of amyloid-beta by M_2 compared to M_1
        self.theta = 0.9

        # Concentration of Amyloid Beta42 plaque outside at which the process rate is half maximal (g/ml)
        # #TODO revoir...
        # self.K_ABpo = 10**(3) * self.ABo = 10**(3) * 7*10**(-6) = 7 * 10**(-3)
        self.K_ABpo = 7e-3  # 7 * 10 ** (-3)

        # Creation rate of GSK3 by Amyloid Beta42 inside neurons (/day)
        self.lambda_ABiG = 0.25

        # Degradation rate of GSK3 (/day)
        # self.d_G = (math.log(2)) / (self.GSK3half / 24)
        # self.GSK3half : The half-life of GSK3 = 41 +/- 4 h
        self.d_G = (math.log(2)) / (41 / 24)

        # Creation rate of tau in health (g/ml/day)
        self.lambda_tau = 26.3e-12  # 26.3 * 10 ** (-12)

        # Degradation rate of tau proteins (/day)
        # self.d_tau = (math.log(2)) / self.Tauhalf
        # self.Tauhalf : The half-life of tau proteins in humain = 23 days
        self.d_tau = (math.log(2)) / 23

        # Creation rate of tau by GSK3 (/day) # TODO: Ceci est une valeur posée par Seyed
        self.lambda_Gtau = 0.25

        # Degradation rate of intracellular NFT (/day)
        self.d_Fi = 1.0e-2 * self.d_tau  # 10 ** (-2) * self.d_tau

        # Degradation rate of extracellular NFT (/day)
        self.d_Fo = 1.0e-1 * self.d_tau  # 10 ** (-1) * self.d_tau

        # Production rate of NFT by tau (/day)
        # 60% of the hyperphosphorylated tau become NFT
        self.lambda_tauFi = 0.6 * self.d_Fo

        # Degradation rate of neurons by F_i (/day)
        # self.d_FiN = ((4 + 4 * self.gamma) / (3 + 2 * self.gamma)) * self.d_N
        # self.gamma : I_10 inhibition ratio = 1 (Hao)
        # self.d_N : the death rate of neuron (/day) = ln2/10ans = ln2/3650
        #                           correction (2 années bissextiles): ln2/3648
        # TODO: Revoir d_N (Seyed semble pas sur)
        self.d_FiN = ((4 + 4 * 1) / (3 + 2 * 1)) * (math.log(2) / 3648)

        # Half-saturation of intracellular NFTs (g/ml)
        # These: Assuming that in AD, 70% of hyperphosphorylated tau proteins (whose concentration in disease is
        # 490 pg/ml) are in NFT form.
        # self.K_Fi = 0.7 * self.htau
        # self.htau : concentration of hyperphosphorylated tau in disease (g/ml) = 490 * 10**(-12)
        self.K_Fi = 0.7 * 490e-12  # 0.7 * 490 * 10 ** (-12)

        # Degradation rate of neurons by T_alpha (TNF-alpha) (/day)
        self.d_TaN = (1 / 2) * self.d_FiN

        # Half-saturation of T_alpha (TNF-alpha) (g/ml)
        self.K_Ta = 2.5e-5  # 2.5 * 10 ** (-5)

        # half-saturation of IL-10 (g/ml)
        self.K_I10 = 2.5e-6  # 2.5 * 10 ** (-6)

        # ?? (g/astrocyte) #TODO ??
        self.W_A = 1.0e-12  # 10 ** (-12)

        # Creation rate of astrocytes by Amyloid Beta42 plaque outside (/day)
        # Hao16 lambda_AA_beta^o
        self.lambda_ABpoA = 1.793

        # Production/activation rate of astrocytes by TNF-alpha (/day)
        # Hao16 lambda_AT_alpha
        self.lambda_TaA = 1.54

        # Death rate of astrocytes (/day)
        # self.d_A = (math.log(2) / self.Astrocyteshalf) * (1 / 10)  #TODO: Pk *0.1 ?
        # self.Astrocyteshalf : Half-life of astrocytes (day) = 600 #TODO: Valeur sort d'où?
        self.d_A = (math.log(2) / 600) * (1 / 10)

        # Degradation rate of Amyloid Beta42 oligomer outside (/day)
        self.d_ABoo = (1 / 10) * self.d_ABmo

        # Creation rate of Amyloid Beta42 oligomer outside by Amyloid Beta42 monomer outside (/day)
        self.lambda_ABmoABoo = 5 * (1 / 25) * self.d_ABoo

        # Creation rate of microglias by F_o (NFT) (/day)
        self.lambda_FoM = 2e-2  # 2 * 10 ** (-2)

        # Average of extracellular NFTs (g/ml)
        # #TODO: Valeur prise pour être < K_Fi = 3.36e-10 (car more NFT reside within neurons than outside of them)
        self.K_Fo = 1.0e-11  # 10 ** (-11)

        # Creation rate of microglias by Amyloid Beta42 plaque outside (g/ml/day)
        # TODO: Revoir. He puts : (0.015 * 0.047 - 2 * 10 ** (-2) * 10 ** (-11)) / 10 ** (-6) ; pk?
        self.lambda_ABpoM = (0.015 * 0.047 - 2e-2 * 1.0e-11) / 1.0e-6

        # Degradation rate of microglias (/day) [Here he takes it equal to the death rate of macrophage]
        self.d_M = 0.015

        # Creation rate of M_1 by microglias (/day)
        self.lambda_MM1 = 9.3e-3  # 9.3 * 10 ** (-3)

        # Death rate of M_1 (proinflammatory microglia) (/day)
        self.d_M1 = 0.015

        # Proinflammatory / anti-inflammatory microglia ratio (M_1/M_2 ratio)
        self.beta = 10

        # Production rate of M_2 by TGF-beta (/day)
        # [The rate by which TGF-beta affects the change of phenotype from M1 to M2]
        self.lambda_TBM2 = 6e-3  # 6 * 10 ** (-3)

        # Death rate of M_2 (anti-inflammatory microglia) (/day)
        self.d_M2 = 0.015

        # Half-saturation of MCP-1 (g/ml) [We consider the value for MCP-1 saturation for influx of macrophages as K_P]
        self.K_P = 5e-9  # 5 * 10 ** (-9)

        # Concentration of M_1^hat at equilibrium (g/ml)
        self.M1hateq = 8.64e-7  # 8.64 * 10 ** (-7)

        # Death rate of M_1^hat macrophages (/day) [Hao]
        self.d_M1hat = 0.015

        # Production rate of M_1^hat by MCP-1 (/day)
        # self.lambda_PM1hat = (self.M1hat * self.d_M1hat) / self.K_P
        # self.M1hat : M_1hat (g/ml) = 0.04  #TODO: pk?
        # self.K_P : Half-saturation of MCP-1 (g/ml) = 5*10**(-9)  #TODO: pk? (Hao: 6e-9)
        self.lambda_PM1hat = (0.04 * self.d_M1hat) / 5e-9

        # Death rate of M_2^hat macrophages (/day) [Hao]
        self.d_M2hat = 0.015

        # Degradation rate of TGF-beta (/day) [Hao]
        self.d_TB = 3.33e2  # 3.33 * 10 ** (2)

        # Production rate of T_beta (/day)
        # self.lambda_TB = self.d_TB * self.K_TB
        # self.K_TB : half-saturation of TGF-beta (T_beta) = 2.5*10**(-7) [Hao]
        self.lambda_TB = self.d_TB * 2.5e-7

        # Production rate of T_beta by M_1 (/day)
        self.lambda_M1TB = 1.5e-2  # 1.5*10**(-2)

        # Production rate of TGF-beta by M_1^hat (/day)
        self.lambda_M1hatTB = 1.5e-2   # 1.5*10**(-2)

        # Production rate of IL-10 by M_2 (/day) #TODO Vérif.. These: (1.2 ± 0.16) × 10^−12 g/ml/day
        self.lambda_M2I10 = 6.67e-3  # 6.67*10**(-3)

        # Degradation rate of IL-10 (/day)
        self.d_I10 = 8.32

        # Production rate of TNF-alpha by M_1^hat (/day) #TODO Vérif.. These: (3.09 ± 1.8) × 10^−12 g/ml/day
        self.lambda_M1hatTa = 1.07e-1  # 1.07*10**(-1)

        # Production rate of TNF-alpha by M_1 (/day) #TODO Vérif. selon précédente (g/ml/day)
        self.lambda_M1Ta = self.lambda_M1hatTa

        # Degradation rate of TNF-alpha (/day) [Hao]
        self.d_Ta = 55.45

        # Half-saturation of M1 (g/ml) [Hao]
        self.K_M1 = 0.03

        # Half-saturation of M2 (g/ml) [Hao]
        self.K_M2 = 0.017

        # Half-saturation of M_1^hat (g/ml) [Hao]
        self.K_M1hat = 0.04

        # Degradation rate of  MCP-1 (/day) [Hao]
        self.d_P = 1.73

        # Creation rate of MCP-1 by astrocytes (/day)
        # Code : self.lambda_AP = (self.K_P*self.d_P)/self.A_0
        # Problème : thèse :
        #       lambda_AP = (d_P * P) / A_0
        #                 = (1.73/day * 3e-10g/ml) / 0.14g/cm^3 = NOPE
        # #TODO: Revoir.  En attendant prenons la valeur de Hao
        self.lambda_AP = 6.6e-8

