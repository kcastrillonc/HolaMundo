import numpy as np
from matplotlib import pyplot

def distancia(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)

def main():
    ### Parámetros de la simulación

    tau = 0.53           # Asociado a la viscosidad cinemática
    plot_cada = 50       # Cada cuanto se toma una imagen de la simulación

    Nx = 400             # Número de sitios en dirección x
    Ny = 100             # Número de sitios en dirección y
    Nt = 7000           # Número de pasos de tiempo
    Nu = 9               # Número de velocidades por nodo

    ### Velocidades
    # Se crean arreglos de Nu elementos
    # El primer término corresponderá a la velocidad nula (permanecer en el nodo).
    # El segundo término correspondera a la velocida que únicamente tiene componente y positiva.
    # Los términos siguientes recorreran las demás velocidades en el orden que tienen pasando de una a otra en sentido horario.

    cxd = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])            # Componente x de la dirección de la velocidad
    cyd = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])            # Componente y de la dirección de la velocidad

    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])         # Pesos asociados a cada dirección


    ### Condiciones iniciales

    F = np.ones((Ny, Nx, Nu)) + 0.01*np.random.randn(Ny, Nx, Nu)       # La definición de la cuadrícula. Cada sitio tiene sus Nu velocidades
                                                                       # y se añade un ruido aleatorio en la inicialización.

    F[:, :, 3] = 2.0           # Velocidad inicial del fluido, tiene dirección 3 (hacia la derecha).

    # Obstáculo

    obstaculo = np.full((Ny,Nx), False)          # Se crea un arreglo que cubra la cuadricula.
                                                 # Si hay un obstáculo, tendrá valor True; False en caso contrario.

    for y in range(0, Ny):
        for x in range(0, Nx):
            if(distancia(Nx//4,Ny//2, x, y) < 13):
                obstaculo[y][x] = True

    ### Ciclo Principal

    figure, (ax1,ax2) = pyplot.subplots(2)      # Para plotear más adelante

    for it in range(Nt):
        print(it)                               # Imprimir el paso de tiempo para facilitar debugging

        ## Condiciones de Frontera

        F[:,-1, [6,7,8]] = F[:,-2, [6, 7, 8]]           # Fronteras abiertas en x=Nx
        F[:, 0, [2,3,4]] = F[:, 1, [2, 3, 4]]           # Fronteras abiertas en x=0

        #F[-1,:, [4,5,6]] = F[-2,:, [4, 5, 6]]           # Fronteras abiertas en y=Ny
        #F[0, :, [8,1,2]] = F[1, :, [8, 1, 2]]           # Fronteras abiertas en y=0

        ## Streaming
        # La función roll de la librería numpy permite hacer streaming fácilmente

        for i, cy, cx in zip(range(Nu), cyd, cxd):
            F[:, :, i] = np.roll(F[:, :, i], cy, axis = 0)         # La velocidad se transmite a los vecinos en y
            F[:, :, i] = np.roll(F[:, :, i], cx, axis = 1)         # La velocidad se transmite a los vecinos en x

        ## Streaming en un obstáculo

        bordeObs = F[obstaculo, :]
        bordeObs = bordeObs[:, [0, 5, 6, 7, 8, 1, 2, 3, 4]]        # Se invierte la dirección de las velocidades en el obtaculo

        # Variables del fluido

        rho = np.sum(F, 2)              # Densidad

        ux = np.sum(F*cxd, 2)/rho           # Momento en dirección x
        uy = np.sum(F*cyd, 2)/rho           # Momento en dirección y

        F[obstaculo, :] = bordeObs          #
        ux[obstaculo] = 0                   # Se hace cero la velocidad en el obstáculo
        uy[obstaculo] = 0                   #

        ## Colision

        Feq = np.zeros(F.shape)

        for i, cy, cx, w in zip(range(Nu), cyd, cxd, weights):
            Feq[:, :, i] = rho*w*(1 + 3*(cx*ux + cy*uy) + (9/2)*((cx*ux + cy*uy)**2) - (3/2)*(ux**2 + uy**2))

        ### Evolución temporal

        F = F + (-1/tau)*(F-Feq)


        ### Plotting (solo se puede tener un plot activo al tiempo)

        if(it%plot_cada == 0):

            dvydx = ux[2:, 1:-1] - ux[0:-2, 1:-1]
            dvxdy = uy[1:-1, 2:] - uy[1:-1, 0:-2]

            rotacional = dvydx - dvxdy                            # Cálculo de la vorticidad


            im1 = ax1.imshow(rotacional, cmap="bwr")
            im2 = ax2.imshow(np.sqrt(ux**2+uy**2))

            ax1.set_title("Vorticidad")
            ax2.set_title("Magnitud del campo de velocidades")

            #pyplot.imshow(np.sqrt(ux**2+uy**2))                  # Mostrar la magnitud del campo de velocidades

            #pyplot.show(block=False)

            pyplot.pause(0.01)

            pyplot.cla()




if __name__=="__main__":
    main()
