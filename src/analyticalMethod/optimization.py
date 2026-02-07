import numpy as np
import matplotlib.pyplot as plt

G = 6.67430e-11
M_SUN = 1.989e30
M_VENUS = 4.867e24
R_VENUS_ORBIT = 108.2e9
R_TIERRA = 1.496e11
R_VENUS_BODY = 6051.8e3 

H_ATMOSPHERE = 350e3  
SAFE_RADIUS = R_VENUS_BODY + H_ATMOSPHERE

MU_SUN = G * M_SUN
MU_VENUS = G * M_VENUS


import numpy as np

def calcular_b_desde_fase(v_lanzamiento, fase_inicial_deg):

    energia = (v_lanzamiento**2) / 2 - MU_SUN / R_TIERRA
    if energia >= 0: 
        return None
    
    a = -MU_SUN / (2 * energia)
    h = R_TIERRA * v_lanzamiento
    p = h**2 / MU_SUN
    
    e_val = np.sqrt(max(0, 1 - p/a))
    cos_nu_llegada = np.clip((p / R_VENUS_ORBIT - 1) / e_val, -1, 1)
    nu_llegada = 2 * np.pi - np.arccos(cos_nu_llegada)
    
    angulo_recorrido_sc = nu_llegada - np.pi
    
    E_inicio = np.pi
    E_llegada = 2 * np.pi - 2 * np.arctan(np.sqrt((1-e_val)/(1+e_val)) * np.tan((2*np.pi - nu_llegada)/2))
    
    M_inicio = E_inicio - e_val * np.sin(E_inicio)
    M_llegada = E_llegada - e_val * np.sin(E_llegada)
    
    n = np.sqrt(MU_SUN / a**3)
    t_vuelo = (M_llegada - M_inicio) / n
    
    omega_v = np.sqrt(MU_SUN / R_VENUS_ORBIT**3)
    
    fase_inicial_rad = np.radians(fase_inicial_deg)
    angulo_venus_final = fase_inicial_rad + (omega_v * t_vuelo)
    
    error_angular = angulo_recorrido_sc - angulo_venus_final
    
    error_angular = (error_angular + np.pi) % (2 * np.pi) - np.pi
    
    b = R_VENUS_ORBIT * np.sin(error_angular)
    
    return b


def simulate_analytic_flyby(v_launch_mag, phase, r_launch = 149.6e9, debug=False):
    e_inicial = (v_launch_mag**2)/2 - MU_SUN/R_TIERRA
    impact_parameter_b = calcular_b_desde_fase(v_launch_mag, phase)
    if np.isnan(impact_parameter_b):
        return np.nan
    
    h = r_launch * v_launch_mag
    energy_initial = (v_launch_mag**2)/2 - MU_SUN/r_launch
    a = -MU_SUN / (2 * energy_initial)
    
    term_e = 1 + (2 * energy_initial * h**2) / (MU_SUN**2)
    if term_e < 0: return np.nan
    e = np.sqrt(term_e)
    p = a * (1 - e**2)
    
    cos_nu = (p / R_VENUS_ORBIT - 1) / e
    if abs(cos_nu) > 1.0: return np.nan
    
    nu_intersect = -np.arccos(cos_nu) 
    
    v_perp = MU_SUN / h * (1 + e * np.cos(nu_intersect))
    v_radial = MU_SUN / h * e * np.sin(nu_intersect)
    phi_sc = np.arctan2(v_radial, v_perp)
    v_sc_mag = np.sqrt(v_perp**2 + v_radial**2)
    
    v_venus_mag = np.sqrt(MU_SUN / R_VENUS_ORBIT)
    v_inf_sq = v_sc_mag**2 + v_venus_mag**2 - 2 * v_sc_mag * v_venus_mag * np.cos(phi_sc)
    
    if v_inf_sq < 0: return np.nan
    v_inf_mag = np.sqrt(v_inf_sq)
    
    if abs(impact_parameter_b) < 1.0: 
        impact_parameter_b = 1.0 if impact_parameter_b >= 0 else -1.0
        
    term_hyp = (impact_parameter_b * v_inf_sq) / MU_VENUS
    e_hyp = np.sqrt(1 + term_hyp**2)
    rp = (MU_VENUS / v_inf_sq) * (e_hyp - 1)
    
    if rp < SAFE_RADIUS:
        if debug:
            print('|', end='')
        return np.nan
        
    tan_half_delta = MU_VENUS / (impact_parameter_b * v_inf_sq)
    delta = 2 * np.arctan(tan_half_delta)
    
    vec_v_venus = np.array([v_venus_mag, 0])
    vec_v_sc = np.array([v_sc_mag * np.cos(phi_sc), v_sc_mag * np.sin(phi_sc)])
    vec_v_inf_in = vec_v_sc - vec_v_venus
    
    c, s = np.cos(delta), np.sin(delta)
    R = np.array(((c, s), (-s, c)))
    vec_v_inf_out = np.dot(R, vec_v_inf_in)
    
    vec_v_sc_out = vec_v_venus + vec_v_inf_out
    v_final_mag = np.linalg.norm(vec_v_sc_out)
    
    e_final = (v_final_mag**2)/2 - MU_SUN/R_VENUS_ORBIT
    
    return e_final - e_inicial

def iterative_optimization(v_center = 26000, v_span = 5000, b_center = 45, 
                           b_span = 90, steps=5, grid_res=60, show=True,
                           debug=False):
    best_v = 0
    best_b = 0
    max_e = -np.inf
    
    final_grid = None
    final_extent = None

    print(f"--- INICIANDO OPTIMIZACIÓN ITERATIVA ({steps} PASOS) ---")
    
    for step in range(steps):
        v_min, v_max = v_center - v_span/2, v_center + v_span/2
        b_min, b_max = b_center - b_span/2, b_center + b_span/2
        
        v_vals = np.linspace(v_min, v_max, grid_res)
        b_vals = np.linspace(b_min, b_max, grid_res)
        
        energy_grid = np.zeros((grid_res, grid_res))
        
        for i, b in enumerate(b_vals):
            for j, v in enumerate(v_vals):
                val = simulate_analytic_flyby(v, b, debug=debug)
                energy_grid[i, j] = val if not np.isnan(val) else -np.inf
        
        if np.all(energy_grid == -np.inf):
            print(f"Paso {step+1}: No se encontraron trayectorias seguras en este rango.")
            break
            
        idx = np.unravel_index(np.argmax(energy_grid), energy_grid.shape)
        
        current_max_e = energy_grid[idx]
        current_best_b = b_vals[idx[0]]
        current_best_v = v_vals[idx[1]]
        
        print(f"\nPaso {step+1}:")
        print(f"\tRango V: [{v_min/1000:.1f}, {v_max/1000:.1f}] km/s")
        print(f"\tRango ϕ: [{b_min:.5f}, {b_max:.5f}] º")
        print(f"\t> Óptimo local:")
        print(f"\t\tV={current_best_v/1000:.4f} km/s")
        print(f"\t\tϕ={current_best_b:.5f} º")
        print(f"\t\tE={current_max_e:.2e}")
        
        v_center = current_best_v
        b_center = current_best_b
        
        v_span *= 0.2
        b_span *= 0.2
        
        if step == steps - 1:
            best_v, best_b, max_e = current_best_v, current_best_b, current_max_e
            final_grid = energy_grid
            final_extent = [v_min/1000, v_max/1000, b_min, b_max]

    print("\n--- OPTIMIZACIÓN COMPLETADA ---")
    print(f"Precisión final V: +/- {v_span/2:.2f} m/s")
    print(f"Precisión final ϕ: +/- {b_span/2:.2f} º")

    target_res = 50  # res TIKZ
    bin_size = grid_res // target_res
    p_weight = 4
    v_sub = np.array([np.mean(v_vals[i:i+bin_size]) for i in range(0, grid_res, bin_size)]) / 1000
    b_sub = np.array([np.mean(b_vals[i:i+bin_size]) for i in range(0, grid_res, bin_size)])
    energy_clean = np.where(energy_grid == -np.inf, np.nan, energy_grid)
    with open("assets/datos_malla_power.txt", "w") as f:
        f.write("fase v_lanzamiento energia\n")
        
        for i in range(target_res):
            for j in range(target_res):
                block = energy_clean[i*bin_size:(i+1)*bin_size, j*bin_size:(j+1)*bin_size]

                valid_data = block[~np.isnan(block)]
                
                if valid_data.size == 0:
                    res_val = "nan"
                else:
                    offset = abs(np.min(valid_data)) if np.min(valid_data) < 0 else 0
                    shifted_data = valid_data + offset
                    
                    pow_mean = (np.mean(shifted_data**p_weight))**(1/p_weight)
                    res_val = f"{pow_mean - offset:.5e}"
                
                f.write(f"{b_sub[i]:.5f} {v_sub[j]:.5f} {res_val}\n")
            f.write("\n")
    print(f"Archivo generado usando media de potencia (p={p_weight}).")
    if show:
        plt.figure(figsize=(10, 8))

        grid_masked = np.ma.masked_where(final_grid == -np.inf, final_grid)

        im = plt.imshow(grid_masked, origin='lower', extent=final_extent, aspect='auto', cmap='plasma')
        cbar = plt.colorbar(im)
        cbar.set_label('Energía (J/kg)')

        # plt.plot(best_v/1000, best_b, 'gx', markersize=20, markeredgewidth=4, label='Óptimo Final')

        plt.title(f'Pasos: {steps}')
        plt.xlabel('Velocidad Inicial (km/s)')
        plt.ylabel('Fase de Venus ϕ (º)')
        # plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()

        # fig = plt.figure(figsize=(12, 8))
        # ax = fig.add_subplot(111, projection='3d')

        # x = np.linspace(final_extent[0], final_extent[1], grid_masked.shape[1])
        # y = np.linspace(final_extent[2], final_extent[3], grid_masked.shape[0])
        # X, Y = np.meshgrid(x, y)

        # surf = ax.plot_surface(X, Y, grid_masked, cmap='plasma', edgecolor='none', alpha=0.8)

        # best_z = np.max(grid_masked)
        # ax.scatter(best_v/1000, best_b/1e6, best_z, color='green', s=100, marker='x', 
        #            linewidths=4, label='Óptimo Final', depthshade=False)

        # ax.set_title(f'Superficie 3D - Zoom Final (Nivel {steps})')
        # ax.set_xlabel('Velocidad Inicial (km/s)')
        # ax.set_ylabel('Parámetro de Impacto b (Millones de km)')
        # ax.set_zlabel('Energía (J/kg)')

        # fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label='Energía (J/kg)')

        # ax.view_init(elev=30, azim=225)

        # plt.legend()
        # plt.show()
    

if __name__ == "__main__":
    iterative_optimization(v_center=14001,
                            v_span=28000,
                            b_center=-56,
                            b_span=30,
                           steps=1, grid_res=1000,
                           debug=True,
                           show=False)
    
# v_center=14001,
# v_span=28000,
# b_center=-56,
# b_span=30,