import numpy as np
import matplotlib.pyplot as plt
import integrator as ing

dtFact = 1e-6
tSim = 3e7

_,_,_,hist = ing.optimizationStepSlingshot(5.33275,-4000,tSim,dtFact,True)
ing.graficarTrayectorias(hist)


minF, maxF = 5.3326, 5.333
minDv, maxDv = -4000.0,-4000.01

# fases, deltaV, varH = ing.graficarSuperficieOptimizacion(minF,maxF,minDv,maxDv,10,5,dtFact,tSim,False)
# print([i for i in zip(fases, deltaV, varH)])


# historia = ing.optimizationStepSlingshot(0,3e5,2e7,1e-2,True)