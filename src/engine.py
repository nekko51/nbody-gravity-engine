import numpy as np
import matplotlib.pyplot as plt
import integrator as ing

dtFact = 1e-1
tSim = 3e7
assets = 'assets/'

# _,_,_,hist = ing.optimizationStepSlingshot(5.332368421052632,-4000,tSim,dtFact,True)
# if hist:
#   ing.graficarTrayectorias(hist)

# barrido fases

minF, maxF = 5,6 #5.332160664819945, 5.332576177285318
fases = np.linspace(minF,maxF,500)
maximH = []
for fase in fases:
    variaH = ing.optimizationStepSlingshot(fase,-4000,tSim,dtFact,True)[0]
    maximH.append(variaH)
maximH = np.array(maximH)
with open(assets + "faseVsH.txt", "a") as f:
    f.write(f"fases\tvarH\n")
    for i in zip(fases,maximH):
        f.write(f"{i[0]:.15f}\t{i[1]:.2f}\n")
    f.write('\n')
plt.plot(fases,maximH)
plt.show()




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

