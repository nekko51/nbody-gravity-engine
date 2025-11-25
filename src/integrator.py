import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def vecRelativo(vec1,vec2):
    return vec2-vec1

def distancia(vec1,vec2):
    return np.linalg.norm(vecRelativo(vec1,vec2))

def unitario(vector):
    return vector/np.linalg.norm(vector)

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
                matrizDistancias[i,:,j]=unitario(deltaR)/np.linalg.norm(deltaR)**2
    aceleraciones = matrizDistancias @ masas * G
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
    flat_pos = vectE[3 : 3 + num_coords]
    return np.array(flat_pos).reshape((numeroCuerpos, dimensiones))

def extraerVelocidadesVectorEstado(vectE, dimensiones, numeroCuerpos):
    inicio_vel = 3 + (numeroCuerpos * dimensiones)
    flat_vel = vectE[inicio_vel:]
    return np.array(flat_vel).reshape((numeroCuerpos, dimensiones))

def calculoHamiltoniano(posiciones,velocidades,masas,G):
    masas = np.array(masas)
    r_vec = posiciones[:, np.newaxis, :] - posiciones[np.newaxis, :, :]
    r_matrix = np.linalg.norm(r_vec, axis=2)
    with np.errstate(divide='ignore'):
        inv_r = 1.0 / r_matrix
    np.fill_diagonal(inv_r, 0)
    potencialTerm = 0.5 * np.sum(masas[:, None] * inv_r * masas[None, :])
    V = -G * potencialTerm
    v_sq = np.sum(velocidades**2, axis=1) 
    T = 0.5 * np.sum(masas * v_sq)
    return T + V

def calculoHamiltonianoNave(posiciones,velocidades,masas,G):
    nave = posiciones[-1]
    otros = posiciones[:-1]
    distancias = np.linalg.norm(otros - nave, axis=1)
    distancias[distancias < 1e-10] = 1e-10
    V = - G * masas[-1] * (1/distancias) @ masas[:-1]
    velNave = velocidades[-1]
    T = masas[-1] * velNave.T @ velNave / 2
    return T + V

def calculoOrbitas(initPos,initVel,masas,tiempoSim,dtFactor,G,debug=False):
    N=len(masas)
    dim = initPos.shape[1]
    H = calculoHamiltoniano(initPos,initVel,masas,G)
    naveH = calculoHamiltonianoNave(initPos,initVel,masas,G)
    historial = [[0,H,naveH,*initPos.flatten(),*initVel.flatten()]]
    # instante de tiempo : hamiltoniano del sistema : hamiltoniano de la nave : posiciones : velocidades
    stepNumber = 0
    timePassed = 0
    while timePassed < tiempoSim:
        vectorEstado = historial[-1]
        velActual = extraerVelocidadesVectorEstado(vectorEstado,dim,N)
        posActual = extraerPosicionesVectorEstado(vectorEstado,dim,N)
        velNaveNorm = (velActual[-1] @ velActual[-1])**.5
        if velNaveNorm > 1e-9:
            dt = dtFactor*distanciaMinimaUltimoCuerpo(posActual)#/velNaveNorm
        else:
            dt = 100.0
        nextPos, nextVel = stepRungeKutta4(posActual,velActual,masas,G,dt)
        stepNumber += 1
        timePassed += dt
        if debug and stepNumber%1e3==0:
            print(
                f"\rtSim: {int(timePassed/3600):7}h.\ttimeStep:{int(stepNumber/1000):7}e3.\tdt: {int(dt):10}", 
                end="")
        # H = calculoHamiltoniano(nextPos,nextVel,masas,G)
        naveH = calculoHamiltonianoNave(nextPos,nextVel,masas,G)
        historial.append(
            [timePassed,H,naveH, *nextPos.flatten(), *nextVel.flatten()]
        )
    if debug: print("")
    return np.array(historial)

def graficarTrayectorias(historial, nombres=None, assets="assets/"):
    """
    Toma el historial devuelto por calculoOrbitas.
    Genera 3 gráficas: Trayectorias, Hamiltoniano Total y Energía Mecánica Nave.
    
    historial: Array numpy (Pasos, 3 + N*dim*2) -> [tiempo, H_tot, H_nave, pos..., vel...]
    """
    dim = 3
    total_cols = historial.shape[1]
    
    # AJUSTE DE INDICE: 
    # Restamos 3 columnas iniciales: Tiempo, H_total, H_nave
    n_cuerpos = int((total_cols - 3) / 2 / dim)

    if nombres is None:
        nombres = [f"Cuerpo {i}" for i in range(n_cuerpos)]

    # ==========================================
    # GRÁFICA 1: TRAYECTORIAS (X-Y)
    # ==========================================
    plt.figure(figsize=(10, 10))
    colores = ['black', 'orange', 'blue', 'red', 'gray', 'green', 'purple']

    for i in range(n_cuerpos):
        # El índice de datos empieza ahora en 3
        idx_x = 3 + i * dim
        idx_y = idx_x + 1
        
        x = historial[:, idx_x]
        y = historial[:, idx_y]
        
        color = colores[i % len(colores)]
        # Estilo: la nave (último cuerpo) va en punteado
        es_nave = (i == n_cuerpos - 1)
        estilo = '--' if es_nave else '-'
        grosor = 1.5 if es_nave else 2
        
        plt.plot(x, y, color=color, linestyle=estilo, linewidth=grosor, label=nombres[i], alpha=0.8)
        plt.plot(x[0], y[0], marker='.', color=color, markersize=5, alpha=0.5)
        
        marker = '*' if i == 0 else 'o'
        plt.plot(x[-1], y[-1], marker=marker, color=color, markersize=10 if i==0 else 7)

    plt.axis('equal')
    plt.grid(True, alpha=0.3)
    plt.xlabel('Posición X (m)')
    plt.ylabel('Posición Y (m)')
    plt.title(f'Trayectorias (t = {historial[-1,0]:.2e} s)')
    plt.legend()
    plt.tight_layout()

    fichero_orbitas = assets + "orbitas.png"
    plt.savefig(fichero_orbitas, dpi=300)
    print(f"Guardado: {fichero_orbitas}")
    # plt.show() # -----------------------------------------------------------------

    # Extraer datos comunes de tiempo
    tiempo = historial[:, 0]

    # ==========================================
    # GRÁFICA 2: HAMILTONIANO TOTAL (SISTEMA)
    # ==========================================
    hamiltoniano = historial[:, 1]
    h0 = hamiltoniano[0]
    # Evitar división por cero si H0 es 0 (raro, pero posible)
    var_relativa = (np.max(hamiltoniano) - np.min(hamiltoniano)) / (np.abs(h0) if h0!=0 else 1)

    plt.figure(figsize=(10, 5))
    plt.plot(tiempo, hamiltoniano, color='purple', linewidth=1.5)
    
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Energía Total (J)')
    plt.title(f'Estabilidad del Integrador (Error Rel: {var_relativa:.2e})')
    plt.grid(True, alpha=0.5, linestyle='--')
    plt.tight_layout()

    fichero_h = assets + "hamiltoniano.png"
    plt.savefig(fichero_h, dpi=300)
    print(f"Guardado: {fichero_h}")
    # plt.show() # -------------------------------------------------------------------

    # ==========================================
    # GRÁFICA 3: ENERGÍA MECÁNICA NAVE
    # ==========================================
    # Esta es la columna nueva (índice 2)
    nave_H = historial[:, 2] 
    
    plt.figure(figsize=(10, 5))
    plt.plot(tiempo, nave_H, color='teal', linewidth=1.5, label='Energía Nave')
    
    # Línea de referencia en 0 (Límite de escape)
    plt.axhline(0, color='red', linestyle='--', alpha=0.7, label='Límite de Escape (E=0)')
    
    plt.xlabel('Tiempo (s)')
    plt.ylabel('Energía Específica')
    plt.title('Evolución de Energía de la Nave')
    plt.legend()
    plt.grid(True, alpha=0.5, linestyle='--')
    
    # Añadimos anotación del estado final
    estado_final = "ESCAPANDO" if nave_H[-1] >= 0 else "ORBITANDO"
    plt.text(tiempo[-1], nave_H[-1], f" {estado_final}", verticalalignment='bottom')

    plt.tight_layout()

    fichero_nave = assets + "energia_nave.png"
    plt.savefig(fichero_nave, dpi=300)
    print(f"Guardado: {fichero_nave}")
    plt.show()

def optimizationStepSlingshot(faseVenus,deltaVelNave,tiempoSim,dtFactor,debug=False):
    G = 6.67430e-11
    masas=[1.98850e30,4.86750e24,5.97220e24,1e6]
    rVenus = 1.08200e11
    vVenus = 3.50200e4
    posicionesIniciales = [
        [0.0, 0.0,0.0],    
        [rVenus*np.cos(faseVenus),rVenus*np.sin(faseVenus),0.0], 
        [1.49600e11,0.0,0.0],
        [1.49500e11,0.0,0.0],
    ]
    velocidadesIniciales = [
        [0.0, 0.0,0.0],
        [-vVenus*np.sin(faseVenus),vVenus*np.cos(faseVenus),0.0],   
        [0.0, 2.97800e4,0.0],
        [0.0, 2.97800e4+deltaVelNave,0.0],
    ]
    posicionesIniciales = np.array(posicionesIniciales)
    velocidadesIniciales = np.array(velocidadesIniciales)
    resultados = calculoOrbitas(posicionesIniciales,velocidadesIniciales,masas,tiempoSim,dtFactor,G,debug)
    varHNave = resultados[-1,2] - resultados[0,2]
    return varHNave, faseVenus, deltaVelNave, resultados

def graficarSuperficieOptimizacion(min_fase,max_fase,min_dv,max_dv,resFase,resDV,dtFactor,tiempoSim,debug=False, assets="assets/"):
    print("Iniciando cálculo de superficie...")
    fases = np.linspace(min_fase, max_fase, resFase)
    deltas = np.linspace(min_dv, max_dv, resDV)
    
    X, Y = np.meshgrid(fases, deltas, indexing='ij')
    Z = np.zeros_like(X) 
    
    total_puntos = resFase * resDV
    contador = 0
    
    for i in range(resFase):
        for j in range(resDV):
            val_fase = X[i, j]
            val_delta = Y[i, j]
            
            resultado = optimizationStepSlingshot(val_fase, val_delta, tiempoSim, dtFactor,debug)

            if isinstance(resultado, (tuple, list, np.ndarray)):
                Z[i, j] = resultado[0] 
            else:
                Z[i, j] = resultado
            
            contador += 1
            tantoporuno = int(total_puntos/100) if total_puntos>100 else 1
            if contador%tantoporuno==0:
                print(f"\rCalculando: {contador:4}/{total_puntos} ({(contador/total_puntos)*100:.1f}%)")
            
    print("\nCálculo terminado. Generando gráfica...")

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, Z, cmap=cm.viridis, edgecolor='none', alpha=0.9)

    ax.set_xlabel('Fase de Venus (radianes)')
    ax.set_ylabel('Delta V Nave (m/s)')
    ax.set_zlabel('Variación Energía (J/kg)')
    ax.set_title('Optimización de Asistencia Gravitatoria (Slingshot)')

    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label='Ganancia de Energía')
    fichero_nave = assets + "energia_nave.png"
    plt.savefig(fichero_nave, dpi=300)
    print(f"Guardado: {fichero_nave}")
    plt.show()
    return X, Y, Z