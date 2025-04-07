## Description of Dimensionalization Process in the Program

### Adopted Physical Parameters
To ensure consistency with the physical parameters reported in Reference [1], the program employs various material properties as illustrated in the accompanying figure. 
![](phy.png)

Using the above values, the shear stress applied to the droplet is calculated as $$\tau = \rho_c (\delta u_d)^2$$where $\rho_c = 1000 \rm{kg/m^3}$ is the mass density of the continumm phase fluid.

The implementation incorporates Taylor's frozen turbulence hypothesis to convert spatial sequences along the streamwise direction of the three-dimensional flow field into temporal sequences. The temporal scale $\delta t$ is determined by the streamwise grid spacing $\Delta x$ and a certain value $0.8U_{inf}$ is chosen as the local convection velocity:
$$\delta t = \Delta x/(0.8U_{inf})$$

### Reynolds Number Conversion
The experimental configuration in Reference [1] was established for Taylor-Couette flow, where the Reynolds number is defined as:
$$Re = \frac{\omega_i r_i d}{\nu}$$
where $\omega_i$ represents the angular velocity of the inner cylinder, $r_i$ the inner cylinder radius, $d$ the gap width, and $\nu$ the kinematic viscosity.

For Taylor-Couette flow, the dimensionless torque parameter $G$ characterizes the externally applied torque:
$$G = \frac{\tau}{2\pi L\rho \nu^2}$$
where $\tau$ denotes the applied torque, which equivalently corresponds to the total wall shear stress ($\tau = \tau_w$). Empirical analysis yields the relationship:
$$G = KRe^{1.58}$$

Through substitution of these expressions into the definition of friction Reynolds number $Re_\tau = \frac{u_\tau d}{\nu}$, we establish the conversion relationship between the experimental Reynolds number $Re$ and the friction Reynolds number $Re_\tau$ as reported in Reference [1].

[1] Yi L, Toschi F, Sun C. Global and local statistics in turbulent emulsions. Journal of Fluid Mechanics. 2021;912:A13. doi:10.1017/jfm.2020.1118
