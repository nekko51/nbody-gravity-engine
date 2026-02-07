# Gravitational Assistance Simulator and Trajectory Optimization

This engine simulates and optimizes gravity assist maneuvers within a reduced Solar System model (Sun, Earth, Venus, and spacecraft). The primary objective is to maximize the probe's final heliocentric velocity by optimizing the angular momentum exchange during the planetary encounter.
Methodology

The simulation employs a dual-stage approach to balance computational efficiency with physical high-fidelity.

## Analytical Approach: Patched Conics

The trajectory is modeled using the patched conics approximation, dividing the flight into distinct two-body problems:

- Heliocentric Transfer: The Sun-probe system is treated as an isolated central force problem. Velocity is determined via the vis-viva equation:

$$
v^2=\mu_\odot\left(\frac{2}{r}-\frac{1}{a}\right)
$$

- SOI Dispersion: Upon entering Venus's Sphere of Influence ($r<r_\text{SOI}​$), the reference frame shifts to the planet. The encounter is modeled as a hyperbolic scattering problem where the deflection angle $\vartheta$ is a function of the impact parameter $s$ and the hyperbolic excess velocity $v_\infty$​.

## Numerical N-Body Simulation

To account for continuous gravitational perturbations, the system integrates the equations of motion derived from the Lagrangian:

$$
  \ddot{\vec{r}}_{\text{ship}} = - \sum_{j \neq \text{ship}} G m_j \frac{\vec{r}_{\text{ship}} - \vec{r}_j}{|\vec{r}_{\text{ship}} - \vec{r}_j|^3}
$$

The integration uses a 4th-order Runge-Kutta (RK4) algorithm. To manage potential singularities and high gradients during close approach, an adaptive time step ($\delta t$) is implemented based on the distance to the nearest body ($r_\text{min}$​): $\delta t=\eta r_\text{min}$​.

## Optimization and Results

An iterative optimization algorithm was used to sweep the initial condition space (launch velocity $v$ and Venus phase $\phi$). Three primary windows of interest were identified analytically:

Window | Launch Velocity $v$ (km/s) | Initial Venus Phase $\phi$ ($^\circ$)
:-----:|:--------------------------:|:-----------------:
$V_1$​  | [24, 28]                   | [-53, -44]
$V_2$  | [12, 26]                   | [-47, -43]
$V_3$  | [8, 12]                    | [-55, -48]

The window $V_2$​ yielded the highest energy gain and was selected for numerical refinement ($V_f$​) with the following intervals: $v\in[19.5,24]$ km/s and $\phi\in[−44.7,−43.5]$ $^\circ$.

## Technical Conclusions

- *Validation of the Patched Conics Method:* The analytical approximation proved to be a highly effective tool for narrowing the search space. The numerical N-body results confirmed that the optimal windows identified analytically remain valid despite the inclusion of multi-body perturbations.

- *Numerical Precision and Conservation:* The conservation of the system's Hamiltonian served as a quality control metric. Monitoring the total energy verified that the RK4 integrator, combined with the adaptive time step, maintained numerical errors within acceptable bounds throughout the maneuver.

- *Sensitivity Analysis:* The results indicate that the efficiency of the gravity assist is significantly more sensitive to the initial orbital phase than to the absolute launch velocity. Precise timing (phase) is critical to achieving the impact parameter required for maximum deflection and energy transfer.

- *Computational Performance:* The implementation of the adaptive step allowed for high resolution at periapsis without the prohibitive cost of a globally small fixed step, effectively resolving the trajectory even at distances below $6.15\times10^6$ m from the planetary center.
