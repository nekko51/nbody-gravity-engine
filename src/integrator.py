import numpy as np

def norma(vector):
    return np.linalg.norm(vector)

def vecRelativo(vec1,vec2):
    return vec2-vec1

def distancia(vec1,vec2):
    return norma(vecRelativo(vec1,vec2))

def unitario(vector):
    return vector/norma(vector)

def calcularAceleraciones(posiciones,masas,G):
    """Toma un array de vectores posición, un 
    array de escalares masa y la constante 
    de gravitación univarsal. Devuelve
    los valores de las aceleraciones en un array
    de vectores."""
    N=len(masas)
    dim = posiciones.shape[1]
    matrizDistancias = np.zeros((N, dim, N))
    for i in range(N):
        for j in range(N):
            if i!=j:
                deltaR = vecRelativo(posiciones[i],posiciones[j])
                matrizDistancias[i,:,j]=unitario(deltaR)/norma(deltaR)**2
    aceleraciones = matrizDistancias@masas*G
    return aceleraciones

def stepRungeKutta4(posiciones,velocidades,masas,G,dt):
    """Toma un array de vecores posición, un array de vectores velocidad, un array de escalares masa, la constante de gravitación univarsal y el paso de tiempo. Devuelve una dupla con un array de vectores posición y un array de vecores velocidad."""
    k1_v = calcularAceleraciones(posiciones, masas, G)
    k1_r = velocidades
    k2_v = calcularAceleraciones(posiciones + k1_r * 0.5 * dt, masas, G)
    k2_r = velocidades + k1_v * 0.5 * dt
    k3_v = calcularAceleraciones(posiciones + k2_r * 0.5 * dt, masas, G)
    k3_r = velocidades + k2_v * 0.5 * dt
    k4_v = calcularAceleraciones(posiciones + k3_r * dt, masas, G)
    k4_r = velocidades + k3_v * dt
    nextVels = velocidades + (dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
    nextPoss = posiciones + (dt / 6.0) * (k1_r + 2*k2_r + 2*k3_r + k4_r)
    return nextPoss,nextVels

def distanciaMinimaUltimoCuerpo(posiciones):
    """Calcula la distancia mínima entre el último cuerpo (nave)
    y el resto de cuerpos."""
    nave = posiciones[-1]
    otros = posiciones[:-1]
    distancias = np.linalg.norm(otros - nave, axis=1)
    return np.min(distancias)

def extraerPosicionesVectorEstado(vectE, dimensiones, numeroCuerpos):
    num_coords = numeroCuerpos * dimensiones
    flat_pos = vectE[1 : 1 + num_coords]
    return np.array(flat_pos).reshape((numeroCuerpos, dimensiones))

def extraerVelocidadesVectorEstado(vectE, dimensiones, numeroCuerpos):
    inicio_vel = 1 + (numeroCuerpos * dimensiones)
    flat_vel = vectE[inicio_vel:]
    return np.array(flat_vel).reshape((numeroCuerpos, dimensiones))

def calculoOrbitas(initPos,initVel,masas,tiempoSim,dtFactor,G):
    N=len(masas)
    dim = initPos.shape[1]
    historial = [[0,*initPos.flatten(),*initVel.flatten()]]
    stepNumber = 0
    timePassed = 0
    while timePassed < tiempoSim:
        vectorEstado = historial[-1]
        velActual = extraerVelocidadesVectorEstado(vectorEstado,dim,N)
        posActual = extraerPosicionesVectorEstado(vectorEstado,dim,N)
        velNaveNorm = norma(velActual[-1])
        if velNaveNorm > 1e-9: # Si se mueve
            dt = dtFactor*distanciaMinimaUltimoCuerpo(posActual)/velNaveNorm
        else:
            dt = 100.0
        nextPos, nextVel = stepRungeKutta4(posActual,velActual,masas,G,dt)
        stepNumber += 1
        timePassed += dt
        historial.append(
            [timePassed, *nextPos.flatten(), *nextVel.flatten()]
        )
    return np.array(historial)