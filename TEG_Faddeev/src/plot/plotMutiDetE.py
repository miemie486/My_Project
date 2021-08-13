import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import csv
from numpy import genfromtxt
from scipy.interpolate import interp1d

x, y = np.zeros(50), np.zeros(50)
xabs, yabs = np.zeros(50), np.zeros(50)

with open("../../data/detEnergyPlot.csv", 'rb') as csvfile:
    data = csv.reader(csvfile)
    i = 0
    for row in data:
        x[i] = abs(np.double(row[0]))
        y[i] = np.double(row[1])
        xabs[i] = abs(np.double(row[0]))
        yabs[i] = abs(np.double(row[1]))
        i += 1

s0 = 0;
f = interp1d(x[s0:i], y[s0:i], kind='cubic')
horizon = np.zeros(100)

fig, ax = plt.subplots()
ax.scatter(x[s0:i], y[s0:i], marker="o", label = 'Original points')
xmesh = np.linspace(x[s0], x[i - 1], 100)
ax.plot(xmesh, f(xmesh), label = 'Curve_original')
ax.plot(xmesh, horizon, label = 'Horizonal line')
ax.set(xlabel='-Energy/MeV', ylabel='Det',
       title='Find the binding energy of triton \n (mesh number of momemtum = 27)')
ax.grid()
ax.legend()

plt.show()
fig.savefig("NP27.png")
