import numpy as np
import matplotlib.pylab as plt

# 1) lea los datos de datos.dat y almacenelos.

dat=np.genfromtxt("datos_gamma0.5.dat")
par=np.genfromtxt("parametros.dat")

print(par[:,1])

###Varianza x gamma == 3.0

D=par[2,1]/par[1,1]

figure, ax = plt.subplots()

ya = 2*D*dat[0::2,0]

ax.plot(dat[0::2,0],np.power(dat[0::2,1],2))
ax.plot(dat[0::2,0],ya)
ax.set(xlabel='Tiempo', ylabel='Varianza de x')

plt.savefig('varianzax_(gamma='+str(par[1,1])+'_F0.png')

###Velocidades

figure, ax = plt.subplots()

ax.plot(dat[1::2,0],np.power(dat[1::2,1],2))
ax.set(xlabel='Tiempo', ylabel='Varianza de vx')


plt.axhline(y=par[2,1], color='r')

plt.savefig('varianzavx_(gamma='+str(par[1,1])+'_F0.png')

#######

figure, ax = plt.subplots()

xa = np.linspace(np.min(dat[-3,2:]),np.max(dat[-1,2:]),dat[-1,2:].size)
ya = np.exp(-1*np.power(xa,2)/(2*par[2,1]))/(np.sqrt(par[2,1]*2*np.pi))

ax.hist(dat[-1,2:], bins = 400, density=True)
ax.plot(xa,ya)
ax.set(xlabel='vx', ylabel='probabilidad')

plt.savefig('velocidades_final_(gamma='+str(par[1,1])+'_F0.png')

#dat=np.genfromtxt("datos_gamma0.01.dat")

#######





