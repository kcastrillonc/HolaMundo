import numpy as np
import matplotlib.pylab as plt


# 1) lea los datos de datos.dat y almacenelos.

dat=np.genfromtxt("datosx.dat")

figure, ax = plt.subplots()

ax.plot(dat[:,0],np.power(dat[:,1],2))
ax.set(xlabel='tiempo', ylabel='varianza')

plt.savefig('varianzaF0.png')

figure, ax = plt.subplots()

for ii in range(dat[0,2:].size):
    ax.plot(dat[:,0],dat[:,ii+2])

ax.set(xlabel='Tiempo', ylabel='Posición')

plt.savefig('posicionesF0.png')
