import numpy as np
import matplotlib.pyplot as plt
import integrator as ing
# 2.242e7 +2*3.5e8/(2.97800e4-5e3)
tiempoSimulacion = 2e7  # segundos primer corte: 8.29223e6
G = 6.67430e-11
masas=[1.98850e30,4.86750e24,5.97220e24,1e6]

rVenus = 1.08200e11
vVenus = 3.50200e4
fVenus = -.5

posicionesIniciales = [
    [0.0, 0.0, 0.0],    
    [rVenus*np.cos(fVenus),rVenus*np.sin(fVenus),0], 
    [1.49600e11,0,0],
    [1.49500e11,0,0],
]

deltavNave = -3.85e3

velocidadesIniciales = [
    [0, 0, 0],
    [-vVenus*np.sin(fVenus),vVenus*np.cos(fVenus),0],   
    [0, 2.97800e4,0],
    [0, 2.97800e4+deltavNave,0],
]

posicionesIniciales = np.array(posicionesIniciales)
velocidadesIniciales = np.array(velocidadesIniciales)

def graficarTrayectorias(historial, nombres=None):
    """
    Toma el historial devuelto por calculoOrbitas y grafica las trayectorias 2D (X-Y).
    historial: Array numpy (Pasos, 1 + N*dim*2)
    nombres: Lista de strings con los nombres de los cuerpos
    """
    dim = 3 
    total_cols = historial.shape[1]
    n_cuerpos = int((total_cols - 1) / 2 / dim)

    if nombres is None:
        nombres = [f"Cuerpo {i}" for i in range(n_cuerpos)]

    plt.figure(figsize=(10, 10))
    colores = ['black', 'orange', 'blue', 'red', 'white', 'green']
    
    for i in range(n_cuerpos):
        idx_x = 1 + i * dim
        idx_y = idx_x + 1
        x = historial[:, idx_x]
        y = historial[:, idx_y]
        color = colores[i % len(colores)]
        estilo = '--' # if 'Nave' in nombres[i] or i == n_cuerpos-1 else '-'
        grosor = 1 if 'Nave' in nombres[i] else 2
        plt.plot(x, y, color=color, linestyle=estilo, linewidth=grosor, label=nombres[i], alpha=0.8)
        plt.plot(x[0], y[0], marker='.', color=color, markersize=5, alpha=0.5)
        marker = '*' if i == 0 else 'o' 
        plt.plot(x[-1], y[-1], marker=marker, color=color, markersize=8 if i==0 else 6)

    plt.axis('equal')
    plt.grid(True, alpha=0.2)
    plt.xlabel('Posición X (m)')
    plt.ylabel('Posición Y (m)')
    plt.title(f'Simulación N-Cuerpos (t = {historial[-1,0]:.2e} s)')
    plt.legend()
    plt.tight_layout()
    plt.show()

# --- USO ---

# 1. Ejecuta tu simulación y guarda el resultado
resultados = ing.calculoOrbitas(posicionesIniciales, velocidadesIniciales, masas, tiempoSimulacion, 1e-2, G)

# 2. Define los nombres (según el orden de tu lista 'masas')
nombres_cuerpos = ["Sol", "Venus", "Tierra", "Nave"]

# 3. Llama a la función
graficarTrayectorias(resultados, nombres_cuerpos)



