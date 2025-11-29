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
fases = np.linspace(minF,maxF,5)
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


