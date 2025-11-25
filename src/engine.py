import numpy as np
import matplotlib.pyplot as plt
import integrator as ing

dtFact = 1e-6
tSim = 3e7

_,_,_,_,hist = ing.optimizationStepSlingshot(5.3328,-4000,tSim,dtFact,True)
ing.graficarTrayectorias(hist)

# barrido fases

# minF, maxF = 5.33235, 5.33245
# fases = np.linspace(minF,maxF,10)
# maximH = []
# for fase in fases:
#     maximH.append(ing.optimizationStepSlingshot(fase,-4000,tSim,dtFact,True)[2])
# maximH = np.array(maximH)
# plt.plot(fases,maximH)
# plt.show()




# minDv, maxDv = -4000.0,-4000.01

# X, Y, Z = ing.graficarSuperficieOptimizacion(minF,maxF,minDv,maxDv,10,5,dtFact,tSim,False)
# idx_flat = np.nanargmax(Z) 
# i_max, j_max = np.unravel_index(idx_flat, Z.shape)

# print(f"Máximo encontrado en índices: [{i_max}, {j_max}]")
# print(f"Valor Máximo de Energía: {Z[i_max, j_max]:.2e} J/kg")


# radio_idx = 2  

# i_start = max(0, i_max - radio_idx)
# i_end   = min(X.shape[0], i_max + radio_idx + 1)

# j_start = max(0, j_max - radio_idx)
# j_end   = min(X.shape[1], j_max + radio_idx + 1)

# X_zoom = X[i_start:i_end, j_start:j_end]
# Y_zoom = Y[i_start:i_end, j_start:j_end]
# Z_zoom = Z[i_start:i_end, j_start:j_end]

# nuevo_min_fase = X_zoom.min()
# nuevo_max_fase = X_zoom.max()
# nuevo_min_dv   = Y_zoom.min()
# nuevo_max_dv   = Y_zoom.max()

# print("-" * 30)
# print("NUEVOS LÍMITES PARA REFINAR LA BÚSQUEDA:")
# print(f"Fase:    [{nuevo_min_fase:.4f}, {nuevo_max_fase:.4f}]")
# print(f"Delta V: [{nuevo_min_dv:.1f}, {nuevo_max_dv:.1f}]")

