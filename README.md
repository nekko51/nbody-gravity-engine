# Gravitational Assistance Simulator and Trajectory Optimization

This repository contains an astrodynamic simulation engine designed to analyze and optimize gravity assist maneuvers. It focuses specifically on a reduced model of the Solar System, which includes the Sun, Earth, Venus, and the spacecraft whose heliocentric velocity is to be maximized.

## Flight plan

The analysis will be performed using two different methods:

- Analytical calculation: approximation to patched conics and dispersion.
- Numerical simulation: direct integration, using the 
4th order Runge-Kutta algorithm, of the equations of motion for 4 bodies with
adaptive time stepping.

## Computation
### Numerical model
Lagrangian mechanics is used to obtain the equations of motion:

$$
\ddot {\vec r_i}=-\sum_{j\ne i}^nGM_j\frac{\vec r_i-\vec r_j}{\big|\vec r_i-\vec r_j\big|^3}
$$

These differential equations have no analytical solution for the number of 
bodies involved. Therefore, the fourth-order Runge-Kutta integrator (RK4) is used.

Added to this is the need to overcome singularities when carrying out 
the maneuver under study. To handle these as best as possible, an adaptive time step proportional to the distance between the spacecraft and the 
body closest to it is used:

$$
\delta t=\eta r_\text{min}
$$

### Validation of numerical calculation
To ensure the physical accuracy of the simulation, the code
monitors the value of the system's Hamiltonian. This verifies that 
the numerical error of the integrator remains within limits.

---
---
---
After a few manual iterations (at present, optimization 
is not automated), the following orbits are obtained for 
a Venus phase $\vartheta=5.3328\ \text{rad}$ and a push to the spacecraft
that slows it down $\Delta v=-4000\ \text{ms}^{-1}$ relative to Earth, its point of departure:

---
![Trayectoria de Asistencia Gravitacional](assets/orbitas.png)

---

With this same orbit, the change in the spacecraft's Hamiltonian is reflected 
in the environment the instant it approaches Venus:

---
![Trayectoria de Asistencia Gravitacional](assets/energia_nave.png)

---