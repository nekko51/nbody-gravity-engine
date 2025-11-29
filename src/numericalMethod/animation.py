import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import numericalMethod.integrator as ing

dtFact = 1e-3
tSim = 3e7
_,_,_,historial = ing.optimizationStepSlingshot(5.333,-4000,tSim,dtFact,True)
n_cuerpos = 4
dim = 3
assets = "assets/"


# Configuración del vídeo final
FPS = 120                # Frames por segundo (fluidez)
DURACION_DESEADA = 10   # ¿Cuántos segundos reales quieres que dure el GIF?
n_frames_video = FPS * DURACION_DESEADA

# Extraemos el tiempo original (que tiene dt variable)
t_original = historial[:, 0]
t_inicio = t_original[0]
t_final = t_original[-1]

# Creamos un NUEVO vector de tiempo que avanza exactamente igual paso a paso
t_interpolado = np.linspace(t_inicio, t_final, n_frames_video)

# Creamos una nueva matriz 'historial_smooth' donde guardaremos los datos interpolados
# Tiene el tamaño: (número de frames x número de columnas originales)
historial_smooth = np.zeros((n_frames_video, historial.shape[1]))
historial_smooth[:, 0] = t_interpolado

print("Re-muestreando datos para velocidad constante...")

# Interpolamos cada columna de posición (X, Y, etc.) al nuevo tiempo
for col in range(1, historial.shape[1]):
    # np.interp(donde_quiero_saber, donde_tengo_datos, valor_de_los_datos)
    historial_smooth[:, col] = np.interp(t_interpolado, t_original, historial[:, col])

# AHORA USAREMOS 'historial_smooth' PARA ANIMAR EN LUGAR DE 'historial'
data_anim = historial_smooth 

# ==========================================
# 2. CONFIGURACIÓN DE EJES Y FIGURA
# ==========================================
fig, ax = plt.subplots(figsize=(10, 10))
colores = ['black', 'orange', 'blue', 'red', 'gray', 'green', 'purple']

# Calcular límites fijos basados en los datos (usamos data_anim)
all_x = []
all_y = []
for i in range(n_cuerpos):
    idx_x = 3 + i * dim
    idx_y = idx_x + 1
    all_x.extend(data_anim[:, idx_x])
    all_y.extend(data_anim[:, idx_y])

# Margen del 5%
rango_max = max(max(all_x) - min(all_x), max(all_y) - min(all_y)) / 2 * 1.05
centro_x = (max(all_x) + min(all_x)) / 2
centro_y = (max(all_y) + min(all_y)) / 2

ax.set_xlim(centro_x - rango_max, centro_x + rango_max)
ax.set_ylim(centro_y - rango_max, centro_y + rango_max)
ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.set_xlabel('Posición X (m)')
ax.set_ylabel('Posición Y (m)')
ax.set_title("Trayectorias (Velocidad de simulación constante)")

# Inicializar objetos gráficos
lineas = []
puntos = []

for i in range(n_cuerpos):
    color = colores[i % len(colores)]
    es_nave = (i == n_cuerpos - 1)
    estilo = '--' if es_nave else '-'
    
    ln, = ax.plot([], [], color=color, linestyle=estilo, linewidth=1, alpha=0.6)
    lineas.append(ln)
    
    marker = '*' if i == 0 else 'o'
    size = 12 if i == 0 else 8
    pt, = ax.plot([], [], marker=marker, color=color, markersize=size)
    puntos.append(pt)

time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

# ==========================================
# 3. FUNCIÓN DE ACTUALIZACIÓN
# ==========================================
def update(frame_idx):
    # frame_idx va de 0 a n_frames_video
    
    t_actual = data_anim[frame_idx, 0]
    time_text.set_text(f'Tiempo = {t_actual:.2e} s')
    
    for i in range(n_cuerpos):
        idx_x = 3 + i * dim
        idx_y = idx_x + 1
        
        # Para la estela, tomamos desde el inicio hasta el frame actual
        # data_anim ya tiene el paso correcto, así que la estela crece suavemente
        x_trail = data_anim[:frame_idx, idx_x]
        y_trail = data_anim[:frame_idx, idx_y]
        
        x_curr = data_anim[frame_idx, idx_x]
        y_curr = data_anim[frame_idx, idx_y]
        
        lineas[i].set_data(x_trail, y_trail)
        puntos[i].set_data([x_curr], [y_curr])
        
    return lineas + puntos + [time_text]

# ==========================================
# 4. GENERAR GIF
# ==========================================
# interval = 1000ms / FPS
intervalo_ms = int(1000 / FPS) 

ani = FuncAnimation(fig, update, frames=n_frames_video, interval=intervalo_ms, blit=True)

fichero_gif = assets + "animacion_constante.gif"
print(f"Guardando animación ({n_frames_video} frames)...")

# Usar 'pillow' es estándar, pero si tienes ffmpeg instalado será más rápido
ani.save(fichero_gif, writer='pillow', fps=FPS)

print(f"Guardado: {fichero_gif}")
plt.show()