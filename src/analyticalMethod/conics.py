import numpy as np
import matplotlib.pyplot as plt

def getEccentricity(totEnergy: float, 
                    angularMomentum: float, 
                    mass: float, 
                    constantK: float) -> float:
    return np.sqrt(1+2*totEnergy*angularMomentum**2/(mass**3*constantK**2))

def constantsFromState(position, velocity, mass, constantK):
    '''energy, angular momentum'''
    T = velocity @ velocity * .5 * mass
    r = np.linalg.norm(position)
    V = - constantK * mass / r
    energy = T + V
    L_vector = np.cross(position, velocity * mass)
    angularMomentum = np.linalg.norm(L_vector)
    return energy, angularMomentum

def dispersionAngle(impactParameter, totEnergy, constantK):
    return 2*np.arctan(constantK/(2*totEnergy*impactParameter))

def getSOI(semimajor, smallMass, bigMass):
    return semimajor*(smallMass/bigMass)**(2/5)

# print(f'soiV: {getSOI(rVenus,mVenus,mSun)/1e6:.2f} e6\nsoiE: {getSOI(rEarth,mEarth,mSun)/1e6:.2f} e6')
# soiV: 6.1623e8
# soiE: 9.2465e8

def trajectoryRadius(angle, orbitPar, eccentricity):
    return orbitPar/(1+eccentricity*np.cos(angle))

def cartesian(radius, angle):
    return radius*np.cos(angle), radius*np.sin(angle)

def polarTrajectoryFromState(position, velocity, 
                             mass, constantK, 
                             measurementNumber = 1000):
    energy, angularMomentum = constantsFromState(position, velocity, mass, constantK)
    eccentricity = getEccentricity(energy, angularMomentum, mass, constantK)
    angle = np.linspace(0, 2*np.pi, measurementNumber)
    radius = trajectoryRadius(angle, angularMomentum**2/mass**2/constantK, eccentricity)
    x, y = cartesian(radius, angle)
    return x, y 

def bolzano(function):
    sign = np.sign(function)
    out = np.floor(np.abs(np.diff(sign)/2))
    return out

def ellipseIntersec(tra1, tra2):
    relativeTraj = np.linalg.norm(tra2-tra1)
    # posiciones = tra1[bolzano(relativeTraj)==1]
    

if __name__ == '__main__':
    G = 6.67430e-11
    mSun = 1.98850e30
    mEarth = 5.97220e24
    rEarth = 1.49600e11
    vEarth = 2.97800e4

    posVec = np.array([rEarth,0,0])
    velVec = np.array([0,vEarth,0])
    X, Y = cartesian(polarTrajectoryFromState(posVec,velVec,mEarth,G*mSun))
    plt.axis('equal')
    plt.plot(X,Y)
    plt.show()