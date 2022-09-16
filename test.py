#
# Ignorer ce fichier, pour tests. Pas utile dans le modèle.
#

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import ticker

t = np.linspace(30*365, 80*365, 51)
y = 0.1 * (-4.151e-15 * t + 3.460e-10)
y = y[0] / y

# x = [2, 4, 8, 10]
# y = np.array([1e-18, 2e-18, 3e-18, 4e-18]) + 2.1111640000e-5

fig, ax = plt.subplots(nrows=1, ncols=1, sharex="all")

# ax.plot(x, y, ".-", ms=2)
ax.plot(t, y)
# formatter = ticker.ScalarFormatter()
# formatter.set_scientific(True)
# formatter.set_powerlimits((-1, 1))
# ax.yaxis.set_major_formatter(formatter)
#
# fig.canvas.draw()
#
# offset = ax.yaxis.get_major_formatter().get_offset()
# print(offset)
#
# pluspos = offset.find('+')
# val = False
# if pluspos != -1:
#     epos = offset[pluspos + 1:].find("e")
#     minuspos = offset[pluspos + 1:].find("−")
#     offset = offset[pluspos+1:].replace(offset[pluspos+1:][minuspos], "-")
#     offset = format(float(offset), '.3g')
#     val = float(offset)
#     print(val, type(val))
#
# # formatter = ticker.ScalarFormatter(useOffset=val, useMathText=True)
# # formatter.set_scientific(True)
# # formatter.set_powerlimits((-1, 1))
# # ax.yaxis.set_major_formatter(formatter)

plt.show()
