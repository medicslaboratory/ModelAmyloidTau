# Date : 17 January 2022
# Autor : Éléonore Chamberland

# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from matplotlib import ticker

import equations as eqns

y0 = np.zeros(19)
# The initial conditions in g/ml
y0[0] = 10**(-6)  # AB^i (Amyloid-beta monomer inside the neurons)
y0[1] = 10**(-8)  # AB_m^o (Amyloid-beta monomer outside the neurons) # Todo ou  ??
y0[2] = 0         # AB_p^o (Amyloid-beta plaque outside the neurons)
y0[3] = 0         # G (GSK3)
y0[4] = 1.37e-10  # tau (tau proteins)
y0[5] = 3.36e-10  # F_i (NFT inside the neurons)
y0[6] = 3.36e-11  # F_o (NFT outside the neurons)
y0[7] = 7e7       # N (Living neurons)
y0[8] = 7e7       # A (Astrocytes)
y0[9] = 0         # AB_o^o (Amyloid-beta oligomers outside)
y0[10] = 0.047    # M (Microglia) #TODO : lui avait 0.02... changé pour valeur trouvée dans Hao
y0[11] = 0.02     # M_1 (Proinflammatory microglia)
y0[12] = 0.02     # M_2 (Anti-inflammatory microglias)
y0[13] = 0        # M_1^hat (Proinflammatory macrophages)
y0[14] = 0        # M_2^hat (Anti-inflammatory macrophages)
y0[15] = 1.0e-6   # T_{beta} (TGF-beta)
y0[16] = 1.0e-5   # I_10 (IL-10 = Interleukin 10)
y0[17] = 2e-5     # T_{alpha} (TNF-alpha)
y0[18] = 5e-9     # P (MCP-1)

annees = 5/365
decades = int(annees/10)

sol = solve_ivp(eqns.ODEsystem, [0, 365*annees], y0, "RK23")  # "BDF" "LSODA"
method = "solve_ivp_RK23"

# t = np.linspace(0, annees, 365*annees)  # 365*annees pts equidistant between 0 and annees.
# sol = odeint(eqns.ODEsystem, y0, t, tfirst=True)
# method = "odeint"

## Figure
fig = plt.figure(figsize=(35/2.54, 20/2.54))

# Making a list for Label names in the plot
labelname = [r'$A \beta^{i}$', r'$A \beta_{m}^{o}$', r'$A \beta_{p}^{o}$',
             '$GSK 3$', r'$\tau$', '$F_i$', '$F_o$', '$Neurons$',
             'A', r'$A \beta_{o}^{o}$', '$M$', '$M_1$', '$M_2$',
             r'$\hat{M}_1$', r'$\hat{M}_2$', r'$T_{\beta}$', '$I_{10}$',
             r'$T_{\alpha}$', '$P$']

fig.subplots_adjust(hspace=.8, wspace=.8)
for i in range(0, 19):
    ax = fig.add_subplot(4, 5, i+1)

    # For odeint
    # ax.plot(t, sol[:, i])

    # For solve_ivp
    ax.plot(sol.t / 365, sol.y[i, :])

    ax.set_xticks(np.linspace(0, annees, decades+1))
    ax.set(xlabel='Time (years)', ylabel=labelname[i])
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(formatter)

# Save the plot as a .png file
plt.savefig("Figure_"+method+"_"+str(annees)+"y.png")

plt.show()