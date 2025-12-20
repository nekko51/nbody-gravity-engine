import numpy as np
import matplotlib.pyplot as plt

G = 6.67430e-11
M_SUN = 1.989e30
M_VENUS = 4.867e24
R_VENUS_ORBIT = 108.2e9 
R_VENUS_BODY = 6051.8e3 

H_ATMOSPHERE = 350e3  
SAFE_RADIUS = R_VENUS_BODY + H_ATMOSPHERE

MU_SUN = G * M_SUN
MU_VENUS = G * M_VENUS

def simulate_analytic_flyby(v_launch_mag, impact_parameter_b):
    r_launch = 149.6e9 
    
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
        return np.nan 
        
    tan_half_delta = MU_VENUS / (impact_parameter_b * v_inf_sq)
    delta = 2 * np.arctan(tan_half_delta)
    
    vec_v_venus = np.array([v_venus_mag, 0])
    vec_v_sc = np.array([v_sc_mag * np.cos(phi_sc), v_sc_mag * np.sin(phi_sc)])
    vec_v_inf_in = vec_v_sc - vec_v_venus
    
    c, s = np.cos(delta), np.sin(delta)
    R = np.array(((c, -s), (s, c)))
    vec_v_inf_out = np.dot(R, vec_v_inf_in)
    
    vec_v_sc_out = vec_v_venus + vec_v_inf_out
    v_final_mag = np.linalg.norm(vec_v_sc_out)
    
    e_final = (v_final_mag**2)/2 - MU_SUN/R_VENUS_ORBIT
    
    return e_final

def iterative_optimization(steps=5, grid_res=60):
    v_center = 30000 
    v_span = 15000   
    
    b_center = 0
    b_span = 400000e3 
    
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
                val = simulate_analytic_flyby(v, b)
                energy_grid[i, j] = val if not np.isnan(val) else -np.inf
        
        if np.all(energy_grid == -np.inf):
            print(f"Paso {step+1}: No se encontraron trayectorias seguras en este rango.")
            return
            
        idx = np.unravel_index(np.argmax(energy_grid), energy_grid.shape)
        
        current_max_e = energy_grid[idx]
        current_best_b = b_vals[idx[0]]
        current_best_v = v_vals[idx[1]]
        
        print(f"Paso {step+1}:")
        print(f"  Rango V: [{v_min/1000:.1f}, {v_max/1000:.1f}] km/s")
        print(f"  Rango b: [{b_min/1e6:.1f}, {b_max/1e6:.1f}] Mkm")
        print(f"  --> Óptimo local: V={current_best_v/1000:.4f} km/s, b={current_best_b/1e6:.4f} Mkm, E={current_max_e:.2e}")
        
        v_center = current_best_v
        b_center = current_best_b
        
        v_span *= 0.2
        b_span *= 0.2
        
        if step == steps - 1:
            best_v, best_b, max_e = current_best_v, current_best_b, current_max_e
            final_grid = energy_grid
            final_extent = [v_min/1000, v_max/1000, b_min/1e6, b_max/1e6]

    print("\n--- OPTIMIZACIÓN COMPLETADA ---")
    print(f"Precisión final V: +/- {v_span/2:.2f} m/s")
    print(f"Precisión final b: +/- {b_span/2000:.2f} km")
    
    plt.figure(figsize=(10, 8))
    
    grid_masked = np.ma.masked_where(final_grid == -np.inf, final_grid)
    
    im = plt.imshow(grid_masked, origin='lower', extent=final_extent, aspect='auto', cmap='plasma')
    cbar = plt.colorbar(im)
    cbar.set_label('Energía (J/kg)')
    
    plt.plot(best_v/1000, best_b/1e6, 'gx', markersize=20, markeredgewidth=4, label='Óptimo Final')
    
    plt.title(f'Zoom Final (Nivel {steps})\nResolución fina alrededor del máximo')
    plt.xlabel('Velocidad Inicial (km/s)')
    plt.ylabel('Parámetro de Impacto b (Millones de km)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

if __name__ == "__main__":
    iterative_optimization(steps=10, grid_res=50)