import numpy as np
import matplotlib.pylab as plt


# 1) lea los datos de datos.dat y almacenelos.

dat=np.genfromtxt("datos.dat")

figure, ax = plt.subplots()

ax.plot(dat[:20:,0],np.power(dat[:20:,2],2))
ax.set(xlabel='tiempo', ylabel='varianza')

plt.savefig('varianza.png')
