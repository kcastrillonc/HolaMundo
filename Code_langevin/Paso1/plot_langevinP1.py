import numpy as np
import matplotlib.pylab as plt

# 1) lea los datos de datos.dat y almacenelos.

par=np.genfromtxt("parametros.dat")
dat=np.genfromtxt("datosP1_k"+str(par[0,1])+"_gamma"
                  +str(par[1,1])+"_T"+str(par[2,1])+".dat")

###Histograma de posiciones para tiempos de 500 en 500

#figure, ax = plt.subplots()

#for iii in dat[-100::50,0]:
#    ax.hist(dat[int(iii),2:], bins = 400, density=True, stacked=True, label='t = '+str(iii))

#ax.legend()
#ax.set_title('Distribuciones de posición (gamma ='+str(par[1,1])+', np = 2000)')

#plt.savefig('posiciones(gamma='+str(par[1,1])+'_hist_P1.png')


###Posiciones gamma ==3

figure, ax = plt.subplots()

for ii in range(dat[0,2:].size):
    ax.plot(dat[:,0],dat[:,ii+2])


ax.set(xlabel='Tiempo', ylabel='Posición',
       title="Posiciones con resorte, k="+str(par[0,1]))

plt.savefig("posicionesP1_k"+str(par[0,1])+"_gamma"
                  +str(par[1,1])+"_T"+str(par[2,1])+".png")

##########################################

#dat=np.genfromtxt("datos_gamma0.01.dat")
#par=np.genfromtxt("parametros.dat")

#######

###Velocidades gamma == 0.01

#figure, ax = plt.subplots()

#for ii in range(dat[0,2:].size):
#    ax.plot(dat[1::2,0],dat[1::2,ii+2])

#ax.set(xlabel='Tiempo', ylabel='Velocidad')

#plt.savefig('velocidades(gamma=0.01_F0.png')
