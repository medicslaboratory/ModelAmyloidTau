# Date : 17 January 2022
# Autor : Éléonore Chamberland

# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from matplotlib import ticker

import equations as eqns

y0 = np.zeros(19)
# The initial conditions in g/mL
y0[0] = 1e-6  # AB^i (Amyloid-beta monomer inside the neurons) # 1e-6 change pour éviter saut départ...
y0[1] = 1.0e-8  # AB_m^o (Amyloid-beta monomer outside the neurons) # Todo ou  ??
y0[2] = 0  # AB_p^o (Amyloid-beta plaque outside the neurons)
y0[3] = 0  # AB_o^o (Amyloid-beta oligomers outside)
y0[4] = 2.5e-7  # G (GSK3)
# Seyed : 0, mais fait choc à cause du terme "p.lambda_ABiG * y[0]" où p.lambda_ABiG = 0.25
y0[5] = 1.37e-10  # tau (tau proteins)  # Hao: Concentration of tau proteins is, in health, 137 pg/ml and, in AD, 490 pg/ml
y0[6] = 3.36e-10  # F_i (NFT inside the neurons)
y0[7] = 3.36e-11  # F_o (NFT outside the neurons)
y0[8] = 0.14  # N (Living neurons)
# TODO : Seyed:7e7 (LSODA_80y_1) Fait pas vrm de sens... ; Valeur trouvée dans Hao = 0.14 (si juste ça modifié
#  LSODA_80y_2)
y0[9] = 0.14  # A (Astrocytes)
# TODO : Seyed:7e7 (LSODA_80y_1) Fait pas vrm de sens... ; Valeur trouvée dans Hao = 0.14 (LSODA_80y_3)
y0[10] = 0.047  # M (Microglia) #TODO : lui avait 0.02... changé pour valeur trouvée dans Hao
y0[11] = 0.02  # M_1 (Proinflammatory microglia)
y0[12] = 0.02  # M_2 (Anti-inflammatory microglias)
y0[13] = 1e-3  # M_1^hat (Proinflammatory macrophages) # 0
y0[14] = 0  # M_2^hat (Anti-inflammatory macrophages)
y0[15] = 1.0e-6  # T_{beta} (TGF-beta)
y0[16] = 1.0e-5  # I_10 (IL-10 = Interleukin 10)
y0[17] = 75e-12  # 2e-5  # T_{alpha} (TNF-alpha) (source: https://doi-org.acces.bibl.ulaval.ca/10.1002/1097-0029(20000801)50:3<184::AID-JEMT2>3.0.CO;2-H)
y0[18] = 5e-9  # P (MCP-1)

annees = 80
decades = int(annees / 10)

sol = solve_ivp(eqns.ODEsystem, [0, 365 * annees], y0, "LSODA")  # "BDF" "LSODA" "RK23
method = "solve_ivp_LSODA"

# t = np.linspace(0, annees, 365*annees)  # 365*annees pts equidistant between 0 and annees.
# sol = odeint(eqns.ODEsystem, y0, t, tfirst=True)
# method = "odeint"

## Figure
fig = plt.figure()
fig.set_size_inches(35 / 2.54, 20 / 2.54, forward=True)

# Making a list for Label names in the plot
labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{p}^{o}$', r'$A \beta_{o}^{o}$', '$GSK 3$', r'$\tau$',
             '$F_i$', '$F_o$', '$Neurons$', 'A', '$M$', '$M_1$', '$M_2$', r'$\hat{M}_1$', r'$\hat{M}_2$',
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
icNameValue = [str(labelname[i]) + "= " + "{:.2e}".format(y0[i]) for i in np.arange(18)]
initcond = "Initial conditions used (in g/cm$^3$) : \n" + ", ".join(icNameValue)
plt.text(0.1, 0.11, initcond, fontsize=9, ha='left', va='top', transform=plt.gcf().transFigure, wrap=True)

# Save the plot as a .png file
my_path = os.path.abspath('Figures')
plt.savefig(os.path.join(my_path, "Figure_" + method + "_" + str(annees) + "y_ModifEqns_16.png"), dpi=180)

plt.show()
