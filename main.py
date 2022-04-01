# Date : 17 January 2022
# Autor : Éléonore Chamberland

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from matplotlib import ticker

import equations as eqns
import parameters as param
p = param.Parameters()

y0 = np.zeros(19)
# The initial conditions in g/mL

# AB^i (Amyloid-beta monomer inside the neurons) # 1e-6 changé pour éviter saut départ...
y0[0] = p.lambda_ABi * (1 + p.AP * p.delta_APi) / p.d_ABi

# AB_m^o (Amyloid-beta monomer outside the neurons)
# y0[1] = 6e-11   # Todo ou  ??
y0[1] = (p.lambda_ABmo * (1 + p.AP * p.delta_APm) + p.lambda_AABmo) / (p.d_ABmo(0) + p.kappa_ABmoABoo *
                                                                       (1 + p.AP * p.delta_APmo))
# AB_o^o (Amyloid-beta oligomers outside) # y0[2] = 5e-13  # 0
y0[2] = p.kappa_ABmoABoo * y0[1] * (1 + p.AP * p.delta_APmo) / (p.d_ABoo + p.kappa_ABooABpo * (1 + p.AP * p.delta_APop))

# AB_p^o (Amyloid-beta plaque outside the neurons)
y0[3] = 0

# G (GSK3)
y0[4] = (p.lambda_InsG * (p.Ins / p.Ins_0))/p.d_G  # = 3.1e-6
# Seyed : 0, mais fait choc à cause du terme "p.lambda_ABiG * y[0]" où p.lambda_ABiG = 0.25

# tau (tau proteins)
# y0[5] = (p.lambda_tau + p.lambda_Gtau)/p.d_tau
y0[5] = 0
# (p.lambda_tau + p.lambda_Gtau * (y0[4] / p.G_0))/p.d_tau => avant fig _39
# = 2.57e-5    # Hao: Concentration of tau proteins is, in health, 137 pg/ml and, in AD, 490 pg/ml

# F_i (NFT inside the neurons)
# y0[6] = p.kappa_tauFi * y0[5] / p.d_Fi
y0[6] = 0

# F_o (NFT outside the neurons)
y0[7] = 0

# N (Living neurons)
y0[8] = p.N_0
# TODO : Seyed:7e7 (LSODA_80y_1) Fait pas vrm de sens... ; Valeur trouvée dans Hao = 0.14 (si juste ça modifié
#  LSODA_80y_2)

# A (Astrocytes)
# y0[9] = 0.14 /2
T_alpha_0 = 2.79e-5
Q = (p.kappa_ABpoA * y0[3] + p.kappa_TaA * T_alpha_0)
y0[9] = (Q * p.A_max) / (Q + p.d_A)

# M (Microglia)
y0[10] = 0
# TODO : Seyed: 0.02 vs Hao: 0.047.

# M_pro (Proinflammatory microglia)
y0[11] = y0[10] * (p.beta / (p.beta + 1))

# M_anti (Anti-inflammatory microglias)
y0[12] = y0[10] * (1 / (p.beta + 1))

# M_pro^hat (Proinflammatory macrophages)
y0[13] = p.Mprohateq/1.5  # 0

# M_anti^hat (Anti-inflammatory macrophages)
y0[14] = 1e-9  # ou (p.kappa_TB * y0[15])/p.d_Mantihat, si y0[15] defini avant # Hao: 0

# T_{beta} (TGF-beta)
y0[15] = (p.kappa_MproTB*y0[11] + p.kappa_MprohatTB*y0[13])/p.d_TB  # 1.0e-6

# I_10 (IL-10 = Interleukin 10)
y0[16] = p.kappa_MantiI10 * y0[12] / p.d_I10  # Hao : 1.0e-5

# T_{alpha} (TNF-alpha)
y0[17] = (p.kappa_MproTa * y0[11] + p.kappa_MprohatTa * y0[13]) / p.d_Ta
# Hao 2e-5
# (source: https://doi-org.acces.bibl.ulaval.ca/10.1002/1097-0029(20000801)50:3<184::AID-JEMT2>3.0.CO;2-H => 75e-12)

# P (MCP-1)
y0[18] = p.kappa_AP*y0[9]/p.d_P  # Hao: 5e-9


annees = 80
decades = int(annees / 10)

sol = solve_ivp(eqns.ODEsystem, [0, 365 * annees], y0, "BDF")  # "BDF" "LSODA" "RK23
method = "solve_ivp_BDF"


"""Generate the figure"""
fig = plt.figure()
fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)

# Making a list for Label names in the plot
labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{o}^{o}$', r'$A \beta_{p}^{o}$', '$GSK 3$', r'$\tau$',
             '$F_i$', '$F_o$', '$Neurons$', 'A', '$M$', '$M_{pro}$', '$M_{anti}$', r'$\hat{M}_{pro}$', r'$\hat{M}_{anti}$',
             r'$T_{\beta}$', '$I_{10}$', r'$T_{\alpha}$', '$P$']

plt.subplots_adjust(hspace=.8, wspace=.8)
for i in range(0, 19):
    ax = fig.add_subplot(4, 5, i + 1)

    # For odeint
    # ax.plot(t, sol[:, i])

    # For solve_ivp
    ax.plot(sol.t / 365, sol.y[i, :])
    if decades > 5:
        indentxaxis = int(decades / 2)
    elif decades < 1:
        indentxaxis = annees
    else:
        indentxaxis = decades
    ax.set_xticks(np.linspace(0, annees, indentxaxis + 1))
    ax.set(xlabel='Time (years)', ylabel=labelname[i])
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(formatter)

# Write the initial values used
plt.subplots_adjust(bottom=0.2)
icNameValue = [str(labelname[i]) + "= " + "{:.2e}".format(y0[i]) for i in np.arange(19)]
initcond = "Initial conditions used (in g/mL) : \n" + ", ".join(icNameValue)
plt.text(0.1, 0.11, initcond, fontsize=9, ha='left', va='top', transform=plt.gcf().transFigure, wrap=True)

# Save the plot as a .png file
my_path = os.path.abspath('Figures')
plt.savefig(os.path.join(my_path, "Figure_" + method + "_" + str(annees) + "y_ModifEqns_45.png"), dpi=180)
# _20 ... _27 : Figures produitent avec Nicolas.
# 34 : 1ere avec modif simon
plt.show()
