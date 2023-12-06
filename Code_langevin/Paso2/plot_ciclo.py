import numpy as np
import matplotlib.pylab as plt

# 1) lea los datos de datos.dat y almacenelos.

par = np.genfromtxt("parametros.dat")

num = 20

alpha1 = (par[num-1,1]*par[num-1,1])/par[num-1,0]
alpha2 = (par[3*num-1,1]*par[3*num-1,1])/par[3*num-1,0]

figure, ax = plt.subplots()

# Compresión isotérmica

plt.axhline(y=par[0,1], xmin=0.05, xmax=0.2, color='blue')
ax.scatter(par[0:num-1,0], par[0:num-1,1], color='blue', label='Compresión isotérmica')

# Compresión adiabática

x = np.linspace(0.2, 1.0, 100)
y = np.sqrt(x*alpha1)

ax.plot(x, y, color='purple')
ax.scatter(par[num-1:2*num-1,0], par[num-1:2*num-1,1], color='purple',
           label='Compresión adiabática')

# Expansión isotérmica

plt.axhline(y=par[2*num,1], xmin=0.25, xmax=1.0, color='orange')
ax.scatter(par[2*num-1:3*num-1,0], par[2*num-1:3*num-1,1], color='orange',
           label='Expasión isotérmica')

# Expansión adiabática

x = np.linspace(0.05, 0.25, 100)
y = np.sqrt(x*alpha2)

ax.plot(x, y, color='green')
ax.scatter(par[3*num-1:4*num-1,0], par[3*num-1:4*num-1,1], color='green',
           label='Expansión adiabática')

ax.set(xlabel='k', ylabel='T', title="Ciclo de Carnot")
ax.legend()

plt.savefig("Ciclo_Carnot.png")






