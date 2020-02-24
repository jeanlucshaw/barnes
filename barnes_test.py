import numpy as np
from os import getcwd
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
import ctypes
from barnes import barnes

def test_func(inx, iny):
    return np.sin(inx) * np.sin(iny)

# Control data set
xg, yg = np.linspace(-3., 3., 30), np.linspace(-3., 3., 30)
XG, YG = np.meshgrid(xg, yg)
zth = test_func(XG, YG)
m, n = xg.size, yg.size
zg = np.zeros_like(XG)


# Test data set scatter
xs = np.random.uniform(-3., 3., 100)
ys = np.random.uniform(-3., 3., 100)
zs = test_func(xs, ys)

# Call to barnes
zg = barnes(xs, ys, zs, xg, yg)

# Visualize results
CMIN, CMAX = -1, 1
_, AX = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 5))

# Barnes interpolation
AX[0].contourf(xg,
               yg,
               zg,
               cmap='jet',
               levels=np.linspace(CMIN, CMAX, 20),
               vmin=CMIN,
               vmax=CMAX)
AX[0].scatter(xs, ys, c=zs, cmap='jet', s=30, edgecolors='k', vmin=CMIN, vmax=CMAX)

AX[1].contourf(XG,
               YG,
               zth,
               cmap='jet',
               levels=np.linspace(CMIN, CMAX, 20),
               vmin=CMIN,
               vmax=CMAX)
AX[1].scatter(xs, ys, c=zs, cmap='jet', s=30, edgecolors='k', vmin=CMIN, vmax=CMAX)

# Axis aesthetics
AX[0].set_title('Barnes')
AX[1].set_title('Control')
AX[1].set(ylim=(-3, 3), xlim=(-3, 3))

plt.show()
