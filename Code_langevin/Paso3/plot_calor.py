import numpy as np
import matplotlib.pylab as plt

# lea los datos de calor.dat y almacenelos.

calor = np.genfromtxt("calor.dat")
calor_avg = np.average(calor, axis = 1)

#print(calor)
#print(calor_avg)

num = 40
x = np.arange(4*num)

##### Plot de dQ

figure, ax = plt.subplots()

# Compresión isotérmica

ax.scatter(x[0:num-1], calor_avg[0:num-1], color='blue', label='Compresión isotérmica')

# Compresión adiabática

ax.scatter(x[num-1:2*num-1], calor_avg[num-1:2*num-1], color='purple',
           label='Compresión adiabática')

# Expansión isotérmica

ax.scatter(x[2*num-1:3*num-1], calor_avg[2*num-1:3*num-1], color='orange',
           label='Expasión isotérmica')

# Expansión adiabática

ax.scatter(x[3*num-1:4*num-1], calor_avg[3*num-1:4*num-1], color='green',
           label='Expansión adiabática')

ax.set(xlabel='Paso', ylabel='dQ', title="dQ entre pasos del Ciclo de Carnot")
ax.legend()

plt.savefig("dQ_Carnot.png")

##### Plot Q (average)

figure, ax = plt.subplots()

y = np.cumsum(calor_avg)

# Compresión isotérmica

ax.scatter(x[0:num-1], y[0:num-1], color='blue', label='Compresión isotérmica')

# Compresión adiabática

ax.scatter(x[num-1:2*num-1], y[num-1:2*num-1], color='purple',
           label='Compresión adiabática')

# Expansión isotérmica

ax.scatter(x[2*num-1:3*num-1], y[2*num-1:3*num-1], color='orange',
           label='Expasión isotérmica')

# Expansión adiabática

ax.scatter(x[3*num-1:4*num-1], y[3*num-1:4*num-1], color='green',
           label='Expansión adiabática')

#ax.set(xlabel='Paso', ylabel='Q', title=, fontsize=16)
ax.set_title('<Q> del Ciclo de Carnot', fontsize = 18.0) # Title
ax.set_ylabel('Q', fontsize = 16.0) # Y label
ax.set_xlabel('Paso', fontsize = 16) # X label

ax.legend()

plt.savefig("Q_Carnot_avg.png")

##### Plot Q

#figure, ax = plt.subplots()

#for ii in range(Q[0,::5000].size):
#    ax.plot(x,Q[:,ii*5000])

#ax.set_title('Q en el Ciclo de Carnot', fontsize = 18.0) # Title
#ax.set_ylabel('Q', fontsize = 16.0) # Y label
#ax.set_xlabel('Paso', fontsize = 16) # X label

#ax.set_ylim([-0.4,0.3])

#plt.savefig('Q_Carnot')

##### Eficiencia por ciclo

#calor = np.genfromtxt("para_eta.dat")
#calor_avg = np.average(calor, axis = 1)

Q = np.cumsum(calor, axis = 0)
Q_avg = y

k1_fin = 0.2  # Constante k al final de Compresión Isotérmica
k2_fin = 1.0  # Constante k al final de Compresión Adiabática

T1 = 0.5                              # Temperatura del baño térmico frío
T2 = np.sqrt(((T1*T1)/k1_fin)*k2_fin) # Temperatura del baño térmico caliente

Q1 = Q[num-1,:] - Q[0,:]           #Cambio de calor en Compresión Isotérmica (Calor sale)
Q2 = Q[3*num-1,:] - Q[2*num-1,:]   #Cambio de calor en Expansión Isotérmica (Calor entra)

W = Q[-1,:]  # Trabajo realizado igual al calor neto transferido en el ciclo.

eta = W/Q2   # Eficiencia igual a la razón entre trabajo realizado y calor absorbido.

Q1_avg = Q_avg[num-1] - Q_avg[0]         # (Calor sale)
Q2_avg = Q_avg[3*num-1] - Q_avg[2*num-1] # (Calor entra)

print(Q1_avg, Q2_avg)

W_avg = Q_avg[-1]  # Trabajo realizado igual al calores cedido al baño caliente.

eta_avg = W_avg/Q2_avg

print(W_avg, eta_avg)

eta_c = 1 - T1/T2  # Eficiencia de Carnot

##### Histogramas de eficiencia

#print(W)
#print(Q2)
#print(eta/eta_c)

dos_eta = [ ]
tres_eta = [ ]
cuatro_eta = [ ]

for ii in range(int(eta.size/2)):
    itr = np.average(eta[ii*2:(ii+1)*2])
    dos_eta.append(itr)

for iii in range(int(eta.size/5)):
    itr = np.average(eta[iii*5:(iii+1)*5])
    tres_eta.append(itr)

for iv in range(int(eta.size/100)):
    itr = np.average(eta[iv*100:(iv+1)*100])
    cuatro_eta.append(itr)

print(eta_c, eta_avg)

eta = np.sort(eta)[1000:-1000]
dos_eta = np.sort(dos_eta)[1000:-1000]
tres_eta = np.sort(tres_eta)[1000:-1000]
cuatro_eta = np.sort(cuatro_eta)[1000:-1000]

figure, ax = plt.subplots()

ax.hist(eta/eta_c, bins = 500, density=True, stacked=True, label='Un ciclo')
ax.hist(dos_eta/eta_c, bins = 500, density=True, stacked=True, label='Dos ciclos')
ax.hist(tres_eta/eta_c, bins = 100, density=True, stacked=True, label='Cinco ciclos')
#ax.hist(cuatro_eta/eta_c, bins = 10, density=True, stacked=True, label='Diez ciclos')

ax.axvline(x=1, color = 'black', linestyle='--', label = 'ηc')
ax.axvline(x=eta_avg, color = 'red', linestyle='-', label = '<η>')

#ax.set_xlim([-500,1000])

ax.set_title('Distribuciones de η según el número de ciclos', fontsize = 18.0) # Title
ax.set_xlabel('η/ηc', fontsize = 16.0) # Y label
ax.set_ylabel('Probabilidad', fontsize = 16) # X label

ax.legend()

plt.savefig('eficiencias.png')

