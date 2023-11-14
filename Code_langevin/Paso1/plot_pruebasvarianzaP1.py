import numpy as np
import matplotlib.pylab as plt

# 1) lea los datos de datos.dat y almacenelos.

par=np.genfromtxt("parametros.dat")
dat=np.genfromtxt("datosP1_k"+str(par[0,1])+"_gamma"
                  +str(par[1,1])+"_T"+str(par[2,1])+".dat")

k = par[0,1]
gamma = par[1,1]
T = par[2,1]

###Varianza x

D=T/gamma
varf=T/k

figure, ax = plt.subplots()

ax.plot(dat[:,0],np.power(dat[:,1],2))
plt.axhline(y=varf, color='r')

ax.set(xlabel='Tiempo', ylabel='Varianza de x', title="Varianza x con resorte, k="+str(k))

plt.savefig("varianzaxP1_k"+str(par[0,1])+
            "_gamma"+str(par[1,1])+"_T"+str(par[2,1])+".png")

#####

figure, ax = plt.subplots()

xa = np.linspace(np.min(dat[-1,2:]),np.max(dat[-1,2:]),dat[-1,2:].size)
ya = np.exp(-1*np.power(xa,2)/(2*varf))/(np.sqrt(varf*2*np.pi))

ax.hist(dat[-1,2:], bins = 350, density=True)
ax.plot(xa,ya, color = 'r')
ax.set(xlabel='vx', ylabel='probabilidad',
       title="Distribución final de x con resorte, k="+str(k))

plt.savefig("posiciones_finalP1_k"+str(par[0,1])+
            "_gamma"+str(par[1,1])+"_T"+str(par[2,1])+".png")

#######





