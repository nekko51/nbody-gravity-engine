import numpy as np
import matplotlib.pyplot as plt
import integrator as ing

dtFact = 1e-3
tSim = 2e7
assets = 'assets/'

# _,_,_,hist = ing.optimizationStepSlingshot(-.795,20700,tSim,dtFact,True)
# ing.graficarTrayectorias(hist)


ing.graficarSuperficieOptimizacion(-44.7/180*np.pi,
                                   -43.5/180*np.pi,
                                   19500,
                                   24000,
                                   20,20,
                                   1e-2,1e7)

# barrido fases

# minF, maxF = -44.7/180*np.pi, -43.5/180*np.pi #5.332160664819945, 5.332576177285318
# fases = np.linspace(minF,maxF,5)
# maximH = []
# for fase in fases:
#     variaH = ing.optimizationStepSlingshot(fase,-4000,tSim,dtFact,True)[0]
#     maximH.append(variaH)
# maximH = np.array(maximH)
# with open(assets + "faseVsH.txt", "a") as f:
#     f.write(f"fases\tvarH\n")
#     for i in zip(fases,maximH):
#         f.write(f"{i[0]:.15f}\t{i[1]:.2f}\n")
#     f.write('\n')
# plt.plot(fases,maximH)
# plt.show() 

# 20700, -.795

# 20kms - 21kms
# -.8rad - -.78rad
# -45.8 - -44.7
