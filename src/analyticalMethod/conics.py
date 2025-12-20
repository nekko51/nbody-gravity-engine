import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

def calculate_intersections_and_plot(target_radius, mass_central, mass_satellite, position, velocity):
    G = 6.67430e-11
    mu = G * (mass_central + mass_satellite)
    
    r_vec = np.array(position, dtype=float)
    v_vec = np.array(velocity, dtype=float)
    r_mag = np.linalg.norm(r_vec)
    v_mag = np.linalg.norm(v_vec)
    
    h_vec = np.cross(np.append(r_vec, 0), np.append(v_vec, 0))
    h = np.linalg.norm(h_vec)
    
    term1 = np.cross(np.append(v_vec, 0), h_vec) / mu
    term2 = np.append(r_vec, 0) / r_mag
    e_vec_3d = term1 - term2
    e_vec = e_vec_3d[:2] 
    
    e = np.linalg.norm(e_vec)
    
    if abs(e - 1.0) < 1e-9:
        print("Orbit is parabolic.")
        return
    
    p = (h**2) / mu
    a = p / (1 - e**2) 
    omega = np.arctan2(e_vec[1], e_vec[0])
    
    print(f"--- Orbit Data ---")
    print(f"Eccentricity (e):       {e:.4f}")
    print(f"Semi-major axis (a):    {a/1000:.2f} km")
    print(f"Periapsis (rp):         { (a*(1-e))/1000:.2f} km")
    print(f"Apoapsis (ra):          { (a*(1+e))/1000:.2f} km")
    print(f"Target Radius:          {target_radius/1000:.2f} km")

    intersection_data = []
    
    if e > 0:
        cos_nu = (p / target_radius - 1) / e
        
        if abs(cos_nu) <= 1.0:
            nu_sol = np.arccos(cos_nu)
            anomalies = [nu_sol, -nu_sol]
            
            for nu in anomalies:
                theta_global = nu + omega
                x_int = target_radius * np.cos(theta_global)
                y_int = target_radius * np.sin(theta_global)
                
                phi_rad = np.arctan2(e * np.sin(nu), 1 + e * np.cos(nu))
                phi_deg = np.degrees(phi_rad)
                
                intersection_data.append({
                    'coords': (x_int, y_int),
                    'angle_deg': phi_deg
                })
        else:
            print("\nWARNING: Circle does not intersect the ellipse.")
    else:
        if abs(p - target_radius) < 1e-3:
            print("Orbit is identical to the target circle.")
        else:
            print("Concentric orbits with no intersection.")
    
    return intersection_data
