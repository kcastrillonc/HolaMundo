import numpy as np
import matplotlib.pylab as plt

# 1) lea los datos de datos.dat y almacenelos.

dat=np.genfromtxt("datos_gamma3.dat")

###Varianza x gamma == 3.0

kboltzmann=1.0
T=1.0
gamma=3.0

D=kboltzmann*T/gamma

figure, ax = plt.subplots()

ya = 2*D*dat[0::2,0]

ax.plot(dat[0::2,0],np.power(dat[0::2,1],2))
#ax.plot(dat[0::2,0],ya)
ax.set(xlabel='Tiempo', ylabel='Varianza de x')

plt.savefig('varianzax_(gamma=3)_F0.png')


###Histograma de posiciones para tiempos de 500 en 500

figure, ax = plt.subplots()

for iii in dat[100::1000,0]:
    ax.hist(dat[int(iii),2:], bins = 400, density=True, stacked=True, label='t = '+str(iii))

ax.legend()
ax.set_title('Distribuciones de posición (gamma = 3.0, np = 10000)')

plt.savefig('posiciones(gamma=3)_hist_F0.png')


###Posiciones gamma ==3

figure, ax = plt.subplots()

for ii in range(dat[0,2:].size):
    ax.plot(dat[0::2,0],dat[0::2,ii+2])


ax.set(xlabel='Tiempo', ylabel='Posición')

plt.savefig('posiciones(gamma=3)F0.png')

#######

dat=np.genfromtxt("datos_gamma0.01.dat")

#######

###Varianza x gamma == 0.01

kboltzmann=1.0
T=1.0
gamma=0.01

D=kboltzmann*T/gamma

figure, ax = plt.subplots()

ya = 2*D*dat[0::2,0]

ax.plot(dat[0::2,0],np.power(dat[0::2,1],2))
#ax.plot(dat[0::2,0],ya)
ax.set(xlabel='Tiempo', ylabel='Varianza de x')

plt.savefig('varianzax_(gamma=0.01)_F0.png')


###Velocidades gamma == 0.01

figure, ax = plt.subplots()

for ii in range(dat[0,2:].size):
    ax.plot(dat[1::2,0],dat[1::2,ii+2])

ax.set(xlabel='Tiempo', ylabel='Velocidad')

plt.savefig('velocidades(gamma=0.01)F0.png')
