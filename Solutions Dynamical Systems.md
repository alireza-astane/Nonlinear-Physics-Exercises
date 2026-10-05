

# Contents

- [[#1. — Lava dome]]
- [[#2. — Trinity Test]]
- [[#3. — Turbulent puff]]
- [[#5. — Squeezing of a viscous liquid film]]
- [[#6. — Hole Nucleation in Liquid Sheets]]
- [[#7. — Examples of dynamical systems]]
- [[#8. — Logistic population model, with preys]]
- [[#9. — Lotka and Volterra prey-predator model]]
- [[#10. — An example of bifurcation]]
- [[#11. — A simplistic model of gravitary buckling]]
- [[#12. — Stick-slip]]
- [[#13. — A quadratic iterated map]]
- [[#14. — Logistic map]]
- [[#15. — A simple population dynamics]]
- [[#16. — Unforced non-linear oscillators]]
- [[#17. — Forced non-linear pendulum]]
- [[#18. — Parametric pendulum]]
- [[#19. — Two-dimensional jets]]
- [[#20. — Bullard’s dynamo (1955)]]
- [[#21. — Rikitake’s dynamo (1956)]]
- [[#22. — Epidemiology]]

---

# 1. — Lava dome

### 1. Physical Origin of the Governing Equations

Because the lava dome is thin relative to its width, vertical fluid acceleration is negligible, and the vertical momentum balance reduces to the hydrostatic equation:

$$\frac{\partial P}{\partial z} = -\rho g$$

Integrating from a height $z$ up to the free surface $z = h(x,t)$, where atmospheric pressure is zero, yields the pressure field:

$$P(x,z,t) = \rho g(h(x,t) - z)$$

The horizontal driving force per unit volume is the spatial gradient of this pressure:

$$\frac{\partial P}{\partial x} = \rho g \frac{\partial h}{\partial x}$$

Lava is characterized by an extremely high dynamic viscosity $\eta$, meaning the flow operates at a Reynolds number near zero (Stokes flow). Inertial forces are entirely negligible, and the horizontal pressure gradient is balanced exclusively by vertical viscous shear stresses:

$$\eta \frac{\partial^2 u}{\partial z^2} = \frac{\partial P}{\partial x} = \rho g \frac{\partial h}{\partial x}$$

Integrating this twice with respect to $z$, applying a no-slip boundary condition at the ground ($u(0) = 0$) and a zero-shear stress condition at the free surface ($\frac{\partial u}{\partial z}\vert{}_{h} = 0$), produces a parabolic velocity profile:

$$u(x,z,t) = \frac{\rho g}{\eta} \frac{\partial h}{\partial x} \left( \frac{z^2}{2} - h z \right)$$

The volumetric flux $q(x,t)$ (flow rate per unit width) is the integral of the velocity across the film thickness:

$$q(x,t) = \int_{0}^{h} u(x,z,t) dz = -\frac{\rho g}{3\eta} h^3 \frac{\partial h}{\partial x}$$

Local mass conservation dictates that changes in height correspond to the divergence of flux ($\frac{\partial h}{\partial t} + \frac{\partial q}{\partial x} = 0$). Substituting the flux expression yields the non-linear diffusion equation presented in Equation (31):

$$\frac{\partial h}{\partial t} = \frac{\rho g}{3\eta} \frac{\partial}{\partial x} \left( h^3 \frac{\partial h}{\partial x} \right)$$

**Correction to the Mass Conservation Equation:** Equation (32) contains a typographical error in the provided problem statement, expressing $\int_{-\infty}^{\infty} h(x,t) dx = \int_{-\infty}^{\infty} h(x,0) dx = Qt$. The physical principle states that the total volume at time $t$ equals the initial volume _plus_ the total volume injected over time. Assuming the lava dome starts from an empty plane at $t=0$, the initial volume integral is zero, making the correct global mass conservation statement:

$$\int_{-\infty}^{\infty} h(x,t) dx = Qt$$

### 2. Dimensional Analysis and Scaling (Dominant Balance)

Let $\beta = \frac{\rho g}{3\eta}$ represent the physical coefficient. To find a similarity solution, we determine the characteristic height scale $H(t)$ and front position scale $X_f(t)$ by substituting them into our two governing physical constraints.

1. The Non-Linear Diffusion Equation (Eq 1): Replacing derivatives with characteristic scales ($\frac{\partial h}{\partial t} \sim \frac{H}{t}$ and $\frac{\partial}{\partial x} \sim \frac{1}{X_f}$) establishes the fundamental dominant balance required to sustain the spreading motion:
    
    $$\frac{H}{t} \sim \beta \frac{H^4}{X_f^2}$$
    
1. Global Mass Conservation (Eq 2): The total cross-sectional area scales linearly with time:
    
    $$H X_f \sim Q t \implies H \sim \frac{Q t}{X_f}$$
    

Substitute the $H$ scaling into the dominant balance equation:

$$\frac{1}{t} \left( \frac{Q t}{X_f} \right) \sim \beta \frac{1}{X_f^2} \left( \frac{Q t}{X_f} \right)^4 \implies \frac{Q}{X_f} \sim \beta \frac{Q^4 t^4}{X_f^6}$$

Multiply by $X_f^6 / Q$ to isolate the horizontal scale:

$$X_f^5 \sim \beta Q^3 t^4 \implies X_f(t) \propto (\beta Q^3 t^4)^{1/5}$$

Substitute $X_f(t)$ back to find the vertical scale:

$$H(t) \propto \frac{Q t}{(\beta Q^3 t^4)^{1/5}} = \left( \frac{Q^2 t}{\beta} \right)^{1/5}$$

### 3. Transforming the PDE to an ODE

Using the exact scalings derived above, we define the dimensionless similarity variable $\xi$ and the scaled shape profile $f(\xi)$:

$$\xi = \frac{x}{(\beta Q^3 t^4)^{1/5}}$$

$$h(x,t) = \left( \frac{Q^2 t}{\beta} \right)^{1/5} f(\xi)$$

We compute the necessary partial derivatives to substitute into Equation (1):

- **Time derivative:** Applying the product and chain rules to $h(x,t)$:
    
    $$\frac{\partial h}{\partial t} = \frac{1}{5t} \left( \frac{Q^2 t}{\beta} \right)^{1/5} f(\xi) + \left( \frac{Q^2 t}{\beta} \right)^{1/5} f'(\xi) \left( -\frac{4}{5} \frac{\xi}{t} \right) = \frac{1}{5t} \left( \frac{Q^2}{\beta t} \right)^{1/5} \left[ f(\xi) - 4\xi f'(\xi) \right]$$
    
- **Non-linear flux term:**
    
    $$\beta \frac{\partial}{\partial x} \left( h^3 \frac{\partial h}{\partial x} \right) = \frac{\beta}{(\beta Q^3 t^4)^{2/5}} \left( \frac{Q^2 t}{\beta} \right)^{4/5} \frac{d}{d\xi} \left( f^3 f' \right) = \frac{1}{t} \left( \frac{Q^2 t}{\beta} \right)^{1/5} \frac{d}{d\xi} \left( f^3 f' \right)$$
    

Equating the two sides, the time-dependent prefactors cancel out entirely, confirming the validity of the similarity assumption:

$$\frac{1}{5} \left[ f(\xi) - 4\xi f'(\xi) \right] = \frac{d}{d\xi} \left( f^3 f' \right)$$

### 4. Direct Integration and Asymptotic Balance

Multiply the ODE by $-5$:

$$f(\xi) - 4\xi f'(\xi) = 5 \frac{d}{d\xi} \left( f^3 f' \right)$$

Integrate both sides with respect to the spatial coordinate from an arbitrary point $\xi$ to the moving front boundary $\xi_f$:

$$\int_{\xi}^{\xi_f} [f(y) - 4y f'(y)] dy = 5 \left[ f^3 f' \right]_{\xi}^{\xi_f}$$

Apply integration by parts to the term $4y f'(y)$:

$$\int_{\xi}^{\xi_f} 4y f'(y) dy = \Big[ 4y f(y) \Big]_{\xi}^{\xi_f} - 4 \int_{\xi}^{\xi_f} f(y) dy$$

Substitute this back into the left-hand side integral:

$$\int_{\xi}^{\xi_f} f(y) dy - \Big[ 4y f(y) \Big]_{\xi}^{\xi_f} + 4 \int_{\xi}^{\xi_f} f(y) dy = 5 \left( f^3(\xi_f) f'(\xi_f) - f^3(\xi) f'(\xi) \right)$$

$$5 \int_{\xi}^{\xi_f} f(y) dy - 4\xi_f f(\xi_f) + 4\xi f(\xi) = 5 f^3(\xi_f) f'(\xi_f) - 5 f^3(\xi) f'(\xi)$$

Apply the physical boundary conditions at the moving front edge ($\xi = \xi_f$). The lava height tapers to zero ($f(\xi_f) = 0$), and the physical fluid flux out of the front edge is zero ($f^3 f' \big\vert{}_{\xi_f} = 0$).

$$5 \int_{\xi}^{\xi_f} f(y) dy + 4\xi f(\xi) = -5 f^3(\xi) f'(\xi)$$

To achieve a closed-form solution, we evaluate the dominant balance near the advancing front. For self-similar wedge profiles approaching a boundary ($f \to 0$ as $\xi \to \xi_f$), the integral of the height function $\int_{\xi}^{\xi_f} f(y) dy$ is a higher-order infinitesimal term (representing the vanishing sliver of area at the tip) compared to the local terms $\xi f(\xi)$ and $f^3 f'$.

Neglecting the integral term yields the exact leading-order ODE for the similarity profile:

$$4\xi f = -5 f^3 f' \implies f^2 f' = -\frac{4}{5} \xi$$

### 5. Solving the Height Profile and Front Position

Separate variables and integrate the simplified ODE:

$$f^2 df = -\frac{4}{5} \xi d\xi$$

$$\frac{f^3}{3} = -\frac{2}{5} \xi^2 + C$$

Applying the boundary condition $f(\xi_f) = 0$ resolves the integration constant $C = \frac{2}{5} \xi_f^2$.

$$f(\xi) = \left( \frac{6}{5} \right)^{1/3} (\xi_f^2 - \xi^2)^{1/3}$$

To determine the dimensionless front location $\xi_f$, evaluate the global mass conservation integral (Equation 32). Assuming a symmetric injection source at $x=0$, the total volume is twice the integral of the right half:

$$2 \int_{0}^{x_f} h(x,t) dx = Q t$$

Substitute the dimensional scalings and the solved function $f(\xi)$:

$$2 \left( \frac{Q^2 t}{\beta} \right)^{1/5} (\beta Q^3 t^4)^{1/5} \int_{0}^{\xi_f} \left( \frac{6}{5} \right)^{1/3} (\xi_f^2 - \xi^2)^{1/3} d\xi = Q t$$

Cancel $Q t$ from both sides:

$$2 \left( \frac{6}{5} \right)^{1/3} \int_{0}^{\xi_f} (\xi_f^2 - \xi^2)^{1/3} d\xi = 1$$

Transform the integral using $u = \xi / \xi_f$, giving $d\xi = \xi_f du$:

$$2 \left( \frac{6}{5} \right)^{1/3} \xi_f^{5/3} \int_{0}^{1} (1 - u^2)^{1/3} du = 1$$

The definite integral evaluates to a constant via the Beta function: $\int_{0}^{1} (1 - u^2)^{1/3} du = \frac{\sqrt{\pi} \Gamma(4/3)}{2 \Gamma(11/6)} \approx 0.6441$.

$$\xi_f = \left[ 2 \left( \frac{6}{5} \right)^{1/3} (0.6441) \right]^{-3/5} \approx 0.814$$

The full spatiotemporal height profile for the lava dome is mathematically defined for $\vert{}x\vert{} \le x_f(t)$ as:

$$h(x,t) = \left[ \frac{6}{5} \frac{3\eta}{\rho g t} \left( x_f(t)^2 - x^2 \right) \right]^{1/3}$$

Where the spreading front position is:

$$x_f(t) = \xi_f \left( \frac{\rho g}{3\eta} Q^3 t^4 \right)^{1/5}$$

# 2. — Trinity Test


### 1. Scaling Law Derivation

The expansion of the fireball in the first milliseconds of a nuclear detonation is so energetic that ambient atmospheric pressure is negligible. The radius $R(t)$ of the shock wave depends exclusively on the total energy released $E$, the initial density of the undisturbed air $\rho$, and the time $t$ elapsed since detonation.

We assume a scaling law of the form:

$$R = C E^\alpha \rho^\beta t^\gamma$$

where $C$ is a dimensionless constant. To find the exponents, we balance the fundamental SI dimensions (Mass $[M]$, Length $[L]$, and Time $[T]$):

- $[R] = [L]$
    
- $[E] = [M][L]^2[T]^{-2}$ (Joules)
    
- $[\rho] = [M][L]^{-3}$ (kg/m$^3$)
    
- $[t] = [T]$
    

Substitute these into the assumed relationship:

$$[L] = ([M][L]^2[T]^{-2})^\alpha ([M][L]^{-3})^\beta ([T])^\gamma$$

$$[L] = [M]^{\alpha + \beta} [L]^{2\alpha - 3\beta} [T]^{-2\alpha + \gamma}$$

For dimensional consistency, the exponents of each fundamental unit must match on both sides:

- **Mass:** $0 = \alpha + \beta \implies \beta = -\alpha$
    
- **Time:** $0 = -2\alpha + \gamma \implies \gamma = 2\alpha$
    
- **Length:** $1 = 2\alpha - 3\beta$
    

Substitute $\beta = -\alpha$ into the length equation:

$$1 = 2\alpha - 3(-\alpha) = 5\alpha \implies \alpha = \frac{1}{5}$$

Back-substituting yields $\beta = -\frac{1}{5}$ and $\gamma = \frac{2}{5}$.

Reassembling the variables provides the theoretical scaling law governing the blast wave:

$$R(t) = C \left( \frac{E}{\rho} \right)^{1/5} t^{2/5}$$

### 2. Validating the Theory with Real Data

To test this theory against the empirical data from the Trinity Test photographs, we linearize the scaling law by taking the base-10 logarithm of both sides:

$$\log_{10} R = \frac{2}{5} \log_{10} t + \frac{1}{5} \log_{10}\left(\frac{E}{\rho}\right) + \log_{10} C$$

This represents a linear equation $y = mx + b$ where $y = \log_{10} R$ and $x = \log_{10} t$.

According to the theoretical derivation, the slope $m$ must be exactly $\frac{2}{5} = 0.40$.

Using the linear regression results from the experimental data:

- **Observed Slope:** $0.40582$
    
- **Coefficient of Determination ($R^2$):** $0.99209$
    

The observed slope of $0.40582$ deviates from the theoretical $0.40$ by less than $1.5\%$. Furthermore, the $R^2$ value of $0.992$ indicates that over $99\%$ of the variance in the shock wave's radius is perfectly explained by this logarithmic time dependence. The theory exceptionally explains the observed physical results.

### 3. Calculating the Released Energy

We deduce the highly classified yield of the explosion using the empirical y-intercept ($b = 2.77674$). From the linearized equation, the intercept isolates the energy and density components:

$$b = \frac{1}{5} \log_{10}\left(\frac{E}{\rho}\right) + \log_{10} C = \frac{1}{5} \log_{10}\left(\frac{C^5 E}{\rho}\right)$$

Rearranging to solve for the energy $E$:

$$10^{5b} = \frac{C^5 E}{\rho} \implies E = \frac{\rho \cdot 10^{5b}}{C^5}$$

Sir G.I. Taylor's rigorous integration of the fluid dynamics equations behind the shock front defined a constant $K$ such that $C^5 = 1/K$. For air behaving as an ideal gas with a specific heat ratio $\gamma = 1.40$, Taylor determined $K = 0.856$. We use standard atmospheric density at the test site, $\rho = 1.25 \text{ kg/m}^3$.

Substituting the values:

$$E = (1.25) (0.856) 10^{5(2.77674)}$$

$$E = (1.07) \cdot 10^{13.8837}$$

$$E = (1.07) \cdot (7.65 \times 10^{13}) \approx 8.18 \times 10^{13} \text{ Joules}$$

To express this in standard explosive yield, we convert Joules to kilotons of TNT ($1 \text{ kiloton TNT} = 4.184 \times 10^{12} \text{ Joules}$):

$$\text{Yield} = \frac{8.18 \times 10^{13}}{4.184 \times 10^{12}} \approx 19.5 \text{ kilotons of TNT}$$

This physical derivation yields $\approx 19.5 \text{ kt}$, which is staggeringly close to the officially declassified Trinity Test yield of $21 \text{ kt}$.

### 4. Estimating the Fissile Mass

The problem states that the fission of a single Uranium atom releases approximately $215 \text{ MeV}$ of energy. First, we convert the energy per atom into Joules ($1 \text{ eV} = 1.602 \times 10^{-19} \text{ J}$):

$$E_{\text{atom}} = 215 \times 10^6 \text{ eV} \times (1.602 \times 10^{-19} \text{ J/eV}) = 3.444 \times 10^{-11} \text{ Joules/atom}$$

Next, we calculate the total number of atoms $N$ required to produce the macroscopic explosive energy $E$:

$$N = \frac{E}{E_{\text{atom}}} = \frac{8.18 \times 10^{13} \text{ J}}{3.444 \times 10^{-11} \text{ J/atom}} \approx 2.375 \times 10^{24} \text{ atoms}$$

We convert the number of atoms into moles using Avogadro's constant ($N_A = 6.022 \times 10^{23} \text{ mol}^{-1}$):

$$n = \frac{2.375 \times 10^{24}}{6.022 \times 10^{23}} \approx 3.94 \text{ moles}$$

Finally, we determine the physical mass. The molar mass of fissile Uranium-235 is $235 \text{ g/mol}$:

$$\text{Mass} = 3.94 \text{ moles} \times 235 \text{ g/mol} \approx 926 \text{ grams}$$

Therefore, approximately **$0.93 \text{ kg}$** of fissile material was required to undergo fission to produce the observed Trinity blast wave.




![[Pasted image 20260930171053.png]]

![[Pasted image 20260930171129.png]]

![[Pasted image 20260930171253.png]]







# 3. — Turbulent puff

### 1. Conservation of Momentum

As the turbulent puff propagates and entrains ambient air, its mass $M(t)$ increases, but its total momentum $P$ remains constant. Assuming the ambient air has a constant density $\rho$, the mass of the puff scales with its characteristic volume: $M(t) \sim \rho R(t)^3$. The conserved momentum is the product of mass and velocity:

$$P \sim M(t) u(t) \sim \rho R(t)^3 u(t) = \text{constant}$$

Since $P$ and $\rho$ are constants, this gives us our first scaling relationship between velocity and size:

$$u(t) \propto R(t)^{-3}$$

### 2. Kinetic Energy Dissipation

The total kinetic energy of the puff scales as $E \sim M(t) u(t)^2 \sim \rho R(t)^3 u(t)^2$.

Because $P \sim \rho R(t)^3 u(t)$ is constant, we can rewrite the energy purely in terms of momentum and velocity:

$$E \sim P u(t)$$

The power dissipated (the rate of change of kinetic energy) is $dE/dt \sim P \dot{u}(t)$.

The problem states that in the turbulent regime, this dissipation rate depends _only_ on $u(t)$ and $R(t)$. To obtain the correct physical dimensions for power (Mass $\cdot$ Length$^2$ / Time$^3$), the ambient density $\rho$ (Mass / Length$^3$) must be implicitly included. By dimensional analysis, the only combination of $\rho$, $u$, and $R$ that yields the units of power is:

$$\frac{dE}{dt} \sim -\rho R(t)^2 u(t)^3$$

### 3. Kinematic Relationship (Entrainment)

We equate the two expressions for the rate of change of kinetic energy:

$$P \dot{u}(t) \sim -\rho R(t)^2 u(t)^3$$

Substitute the constant momentum $P \sim \rho R(t)^3 u(t)$ into the left side:

$$\left(\rho R(t)^3 u(t)\right) \dot{u}(t) \sim -\rho R(t)^2 u(t)^3$$

Dividing both sides by $\rho R(t)^2 u(t)$ yields:

$$R(t) \dot{u}(t) \sim -u(t)^2$$

Now, we use our momentum scaling law $u \propto R^{-3}$ to find $\dot{u}$:

$$\dot{u} \propto \frac{d}{dt}(R^{-3}) = -3R^{-4}\dot{R}$$

Substitute $u$ and $\dot{u}$ back into the differential relationship $R\dot{u} \sim -u^2$:

$$R(-3R^{-4}\dot{R}) \sim -(R^{-3})^2$$

$$-3R^{-3}\dot{R} \sim -R^{-6}$$

$$\dot{R} \sim R^{-3}$$

_(Note: Because $\dot{R} \sim R^{-3}$ and $u \sim R^{-3}$, this confirms the classic entrainment hypothesis for turbulent puffs: the rate of growth $\dot{R}$ is strictly proportional to the propagation velocity $u$.)_

### 4. Final Scaling Laws with Time

We solve the differential equation for $R(t)$:

$$\frac{dR}{dt} \propto R^{-3}$$

$$R^3 dR \propto dt$$

Integrating both sides with respect to time gives:

$$R^4 \propto t$$

$$R(t) \propto t^{1/4}$$

Finally, substitute $R(t)$ back into the velocity-size relationship $u \propto R^{-3}$:

$$u(t) \propto (t^{1/4})^{-3}$$

$$u(t) \propto t^{-3/4}$$

# 5. — Squeezing of a viscous liquid film

### Physics Overview: Lubrication Theory and Squeezing Flow

This problem deals with the classic **Stefan lubrication flow** (or Reynolds' squeeze-film theory).

When a disk or cylinder approaches a flat surface separated by a thin fluid film ($h \ll R$), the fluid cannot escape instantaneously. To squeeze out the liquid, a huge pressure builds up inside the film.


- **Mass Conservation:** As the cylinder moves down at speed $V = -\dot{h}$, it displaces fluid at a volumetric rate of $\pi r^2 V$ within a radius $r$. This displaced fluid must exit laterally through the cylindrical sidewall of area $2\pi r h$.
    
- **Viscous Resistance (Lubrication Approximation):** Because $h \ll R$, vertical velocity gradients ($\partial u_r / \partial z$) are much larger than radial gradients. The dominant physical balance in the radial momentum equation is between the radial pressure gradient driving the flow and the vertical viscous shear stress resisting it:
    
    $$\frac{\partial P}{\partial r} \sim \eta \frac{\partial^2 u_r}{\partial z^2}$$
    

### 1. Scaling Law Derivation (Reynolds' Law) & Pressure Profile

#### Scaling Law Derivation

- **Mass Conservation Scale:**
    
    The volume of fluid squeezed out per unit time within radius $r$ is:
    
    $$Q(r) \sim r^2 V$$
    
    This volume escapes through a ring of circumference $2\pi r$ and height $h$, giving a mean radial velocity $\bar{v}_r$:
    
    $$Q(r) \sim r h \bar{v}_r \implies \bar{v}_r \sim \frac{r}{h} V$$
    
- **Momentum Balance Scale:**
    
    In the low-Reynolds-number lubrication limit ($Re \ll 1$), viscous shear forces balance the radial pressure gradient:
    
    $$\nabla_r P \sim \eta \frac{\bar{v}_r}{h^2}$$
    
    Substituting the scaling for $\bar{v}_r$:
    
    $$\nabla_r P \sim \eta \frac{r V / h}{h^2} \sim \frac{\eta r}{h^3} V$$
    
    This proves the required scaling law:
    
    $$\nabla_r P \sim \frac{\eta r}{h^3} V$$
    

#### Qualitative Pressure Profile

- $\nabla_r P \propto r$: The pressure gradient increases linearly with radius $r$.
    
- Integrating $\nabla_r P \sim -r$ from the edge ($r = R$, where $P = P_{atm} \approx 0$) inward to the center ($r = 0$) shows that the pressure profile $P(r)$ is a **parabola pointing upwards**:
    
    $$P(r) \propto (R^2 - r^2)$$
    
- The pressure reaches its absolute maximum at the center ($r = 0$) and drops smoothly to atmospheric pressure at the outer boundary ($r = R$).
    

### 2. Condition for 1D Poiseuille Flow Approximation

The flow through a circular slice at radius $r$ can be approximated as a 1D (planar) Poiseuille flow between two parallel plates if:

1. **Thin Film Geometry:** $h \ll r \le R$ (the curvature of the ring is negligible locally, and vertical shear dominates horizontal variations).
    
2. **Low Reynolds Number:** $Re = \frac{\rho \bar{v}_r h}{\eta} \ll 1$ (inertial terms $\rho (\mathbf{u}\cdot\nabla)\mathbf{u}$ are completely negligible compared to viscous terms $\eta \nabla^2 \mathbf{u}$).
    

### 3. Exact Relational Derivation & Retrieval of Reynolds' Law

#### Relationship between Pressure Gradient and Average Velocity

For a steady 1D Poiseuille flow between two stationary plates separated by distance $h$, the velocity profile across the gap $z \in [0, h]$ is parabolic:

$$u_r(z) = -\frac{1}{2\eta} \frac{dP}{dr} z(h - z)$$

The average velocity $\bar{v}_r$ across the height $h$ is:

$$\bar{v}_r = \frac{1}{h} \int_{0}^{h} u_r(z) dz = -\frac{h^2}{12\eta} \frac{dP}{dr}$$

Rearranging gives the exact relation for the pressure gradient:

$$\frac{dP}{dr} = -\frac{12\eta}{h^2} \bar{v}_r$$

#### Relating $\bar{v}_r$ to Crushing Speed $V$

By exact volume conservation of an incompressible fluid inside a cylinder of radius $r$:

$$\pi r^2 V = 2\pi r h \bar{v}_r \implies \bar{v}_r = \frac{r V}{2h}$$

#### Retrieving Reynolds' Law

Substitute $\bar{v}_r = \frac{rV}{2h}$ into the pressure gradient formula:

$$\frac{dP}{dr} = -\frac{12\eta}{h^2} \left( \frac{r V}{2h} \right) = -\frac{6\eta r}{h^3} V$$

Taking the magnitude yields the exact version of Reynolds' law:

$$\nabla_r P = \frac{6\eta r}{h^3} V$$

### 4. Evolution of Separation $h(t)$ Under Constant Weight

Integrate the exact pressure profile $P(r)$ to find the total upward viscous force $F_p$ supporting the cylinder weight $W = mg$:

$$P(r) - P_{atm} = \int_{r}^{R} -\frac{dP}{dr'} dr' = \int_{r}^{R} \frac{6\eta r'}{h^3} V dr' = \frac{3\eta V}{h^3} (R^2 - r^2)$$

The total upward force exerted by the fluid is:

$$F_p = \int_{0}^{R} P(r) \cdot 2\pi r \, dr = \frac{3\eta V}{h^3} 2\pi \int_{0}^{R} (R^2 r - r^3) dr = \frac{3\eta V}{h^3} 2\pi \left( \frac{R^4}{4} \right) = \frac{3\pi \eta R^4}{2 h^3} V$$

Equating the viscous lifting force to the weight $W = mg$ (since $V = -\frac{dh}{dt}$):

$$\frac{3\pi \eta R^4}{2 h^3} \left( -\frac{dh}{dt} \right) = mg$$

Separate variables and integrate:

$$-h^{-3} dh = \frac{2mg}{3\pi \eta R^4} dt$$

$$\int_{h_0}^{h(t)} h^{-3} dh = -\frac{2mg}{3\pi \eta R^4} \int_{0}^{t} dt'$$

$$\frac{1}{2}\left( \frac{1}{h(t)^2} - \frac{1}{h_0^2} \right) = \frac{2mg}{3\pi \eta R^4} t$$

$$\frac{1}{h(t)^2} = \frac{1}{h_0^2} + \frac{4mg}{3\pi \eta R^4} t$$

Solving for $h(t)$:

$$h(t) = \frac{h_0}{\sqrt{1 + \frac{4 m g h_0^2}{3\pi \eta R^4} t}} \propto t^{-1/2} \quad \text{for large } t$$

### 5. Does the Cylinder Touch the Bottom? Real Situations

#### Theoretical Prediction

According to purely continuum fluid mechanics, $h(t) \to 0$ only as $t \to \infty$. The required force to displace the final fluid layer diverges as $1/h^3$, meaning **the cylinder theoretically never touches the bottom in finite time**.

#### What Happens in a Real Situation?

In reality, contact occurs in finite time due to physical effects neglected by continuum lubrication theory:

1. **Surface Roughness:** No surface is perfectly flat. When $h(t)$ reaches the micro-scale height of surface asperities (roughness peaks), solid-to-solid contact occurs.
    
2. **Nanoscale / Molecular Forces:** At ultra-thin gaps ($h \lesssim 10\text{ nm}$), continuum Navier-Stokes breaks down. Van der Waals attractive forces pull the surfaces together, or hydrophobic attraction causes sudden contact.
    
3. **Deformation:** High pressures can elastically deform the cylinder/substrate (elastohydrodynamic lubrication).
    
4. **Dissolved Gases / Cavitation:** Tiny gas bubbles or local pressure drops during flow can cause film rupture.

### Full Derivation of $u_r(z) = -\frac{1}{2\eta} \frac{dP}{dr} z(h - z)$

This velocity profile is not pulled out of nowhere—it is derived directly from the fundamental **Navier-Stokes equation** under lubrication assumptions. Here is every single step:

#### Step A: The Governing Equation (Momentum Balance)

The 2D axisymmetric Navier-Stokes equation in cylindrical coordinates for the radial component of velocity $u_r(r,z,t)$ is:

$$\rho \left( \frac{\partial u_r}{\partial t} + u_r \frac{\partial u_r}{\partial r} + u_z \frac{\partial u_r}{\partial z} \right) = -\frac{\partial P}{\partial r} + \eta \left( \frac{\partial^2 u_r}{\partial r^2} + \frac{1}{r}\frac{\partial u_r}{\partial r} - \frac{u_r}{r^2} + \frac{\partial^2 u_r}{\partial z^2} \right)$$

Because $h \ll R$ (thin film) and $Re \ll 1$ (low Reynolds number):

- **Inertia is negligible:** The entire left side (density $\rho \times \text{acceleration}$) is zero.
    
- **Vertical viscous shear dominates:** Because the gap $h$ is tiny compared to $R$, derivatives with respect to $z$ ($\sim 1/h^2$) are vastly larger than derivatives with respect to $r$ ($\sim 1/R^2$).
    

The complex Navier-Stokes equation simplifies strictly to:

$$0 = -\frac{dP}{dr} + \eta \frac{\partial^2 u_r}{\partial z^2} \implies \frac{\partial^2 u_r}{\partial z^2} = \frac{1}{\eta} \frac{dP}{dr}$$

Notice that $P$ depends only on $r$ (pressure is uniform across the thin vertical gap $z$). Thus, at any fixed radius $r$, the term $\frac{1}{\eta} \frac{dP}{dr}$ acts as a **constant** with respect to $z$.

#### Step B: First Integration

Integrate $\frac{\partial^2 u_r}{\partial z^2} = \frac{1}{\eta} \frac{dP}{dr}$ once with respect to $z$:

$$\frac{\partial u_r}{\partial z} = \frac{1}{\eta} \frac{dP}{dr} z + C_1$$

#### Step C: Second Integration

Integrate a second time with respect to $z$:

$$u_r(z) = \frac{1}{2\eta} \frac{dP}{dr} z^2 + C_1 z + C_2$$

#### Step D: Applying No-Slip Boundary Conditions

The fluid touches two solid surfaces:

1. **At the bottom surface ($z = 0$):** The ground is stationary.
    
    $$u_r(0) = 0 \implies C_2 = 0$$
    
2. **At the top surface ($z = h$):** The cylinder moves vertically at speed $V$, but has **no horizontal velocity**.
    
    $$u_r(h) = 0$$
    

Substitute $z = h$ and $C_2 = 0$ into the equation:

$$0 = \frac{1}{2\eta} \frac{dP}{dr} h^2 + C_1 h$$

Divide by $h$ (since $h \neq 0$) and solve for $C_1$:

$$C_1 = -\frac{1}{2\eta} \frac{dP}{dr} h$$

#### Step E: Assembling the Final Profile

Substitute $C_1$ and $C_2$ back into the expression for $u_r(z)$:

$$u_r(z) = \frac{1}{2\eta} \frac{dP}{dr} z^2 - \left(\frac{1}{2\eta} \frac{dP}{dr} h\right) z$$

Factor out $-\frac{1}{2\eta} \frac{dP}{dr} z$:

$$u_r(z) = -\frac{1}{2\eta} \frac{dP}{dr} z (h - z)$$

This is the exact parabolic Poiseuille flow profile.

### 3. Essential Hydrodynamics Theory Required for This Problem

To master squeeze-film and thin-film fluid problems, you need to understand four core concepts:

#### Theory 1: Low-Reynolds-Number (Stokes) Flow

The Reynolds number measures the ratio of inertial forces to viscous forces: $Re = \frac{\rho U L}{\eta}$.

When $Re \ll 1$ (as in high-viscosity liquids or ultra-thin gaps), fluid inertia is completely irrelevant. The fluid has no "memory" or momentum; if you stop pushing, it stops instantly. The Navier-Stokes equation reduces to a linear balance between pressure and viscosity:

$$\nabla P = \eta \nabla^2 \mathbf{u}$$

#### Theory 2: The Thin-Film (Lubrication) Approximation

When a fluid layer has a vertical thickness $H$ much smaller than its horizontal scale $L$ ($H/L \ll 1$):

1. Pressure is uniform across the thickness: $\frac{\partial P}{\partial z} \approx 0 \implies P = P(r,t)$.
    
2. Vertical velocity gradients dominate completely over horizontal ones:
    
    $$\nabla^2 \mathbf{u} \approx \frac{\partial^2 u_r}{\partial z^2}$$
    

#### Theory 3: Mass Conservation for Incompressible Fluids

For an incompressible fluid ($\nabla \cdot \mathbf{u} = 0$), the rate at which volume is displaced by a solid boundary must equal the volumetric flow rate escaping through the sides.

For a cylinder of radius $r$ moving down at speed $V = -\dot{h}$:

- **Displaced Volume Rate:** $\frac{d}{dt}(\pi r^2 h) = \pi r^2 \frac{dh}{dt} = -\pi r^2 V$
    
- **Escaping Volume Rate through edge of height $h$:** $2\pi r \int_0^h u_r(z) dz = 2\pi r h \bar{v}_r$
    
- **Equating them:** $\pi r^2 V = 2\pi r h \bar{v}_r \implies \bar{v}_r = \frac{r V}{2 h}$
    

#### Theory 4: Planar Poiseuille Flow

Whenever a viscous fluid is driven through a thin gap of height $h$ by a pressure gradient $\frac{dP}{dr}$, the steady velocity profile across the gap is parabolic:

$$u(z) = -\frac{1}{2\eta}\frac{dP}{dr} z(h-z)$$

Integrating this profile across the gap $h$ gives the fundamental **Flux-Pressure relation**:

$$\bar{v}_r = \frac{1}{h}\int_0^h u(z) dz = -\frac{h^2}{12\eta}\frac{dP}{dr}$$



# 6. — Hole Nucleation in Liquid Sheets

When a liquid film or sheet at rest (with density $\rho$, thickness $h$, and surface tension $\gamma$) is punctured, it is violently unstable. The local rupture creates an opening where the liquid-gas interface is suddenly exposed. Because a liquid film possesses excess surface energy, the system naturally seeks to minimize its total surface area.

The unbalanced capillary forces at the edge of the hole pull the liquid back, causing the hole to expand rapidly. As the hole grows, the displaced liquid from the film is swept outward and accumulates into a thickened cylindrical border at the edge, known as a **liquid rim** (or _bourrelet_).

### 1. What is the physical mechanism driving the expansion of the hole?

- **Mechanism:** The primary driving force is **surface tension (capillary force)**.
    
- **Explanation:** A liquid film has two free surfaces (top and bottom), each possessing a surface energy per unit area equal to the surface tension $\gamma$. When a hole is created, the creation of the hole destroys the film locally, but the capillary forces acting along the perimeter of the rim pull outward on the liquid. This converts stored surface energy into the kinetic energy of the rapidly expanding rim.
    

### 2. What is the mass flux entering the rim?

As the rim expands outward at an instantaneous velocity $V(t)$, it continuously sweeps up the stationary liquid film lying in its path.

- **Mass per unit area of the film:** The film has a homogeneous density $\rho$ and thickness $h$, meaning the mass per unit area (surface density) is $\rho h$.
    
- **Rate of area swept:** In a small time interval $dt$, a rim segment of unit length advances by a distance $V(t)dt$. Therefore, the area of the flat film swept up per unit time per unit length is equal to $V(t)$.
    
- **Mass flux ($\dot{m}(t)$):** The instantaneous rate of change of mass per unit length of the rim ($\dot{m}(t)$) is the product of the surface density and the sweeping speed:
    
    $$\dot{m}(t) = \rho h V(t)$$
    

### 3. The Capillary Force and Momentum Conservation Equation

To write the equation of motion for this growing-mass system, we must account for both the forces acting on the rim and the momentum brought in by the newly collected mass.

#### A. The Force Exerted by the Liquid Film

Because the liquid film has **two** free surfaces (an upper and a lower surface), surface tension $\gamma$ pulls inward on the rim from both sides of the film along the entire perimeter.

- Each surface exerts a capillary force of magnitude $\gamma$ per unit length.
    
- Therefore, the total net capillary force $F$ pulling the rim outward into the film per unit length is:
    
    $$F = 2\gamma$$
    

#### B. The Momentum Conservation Equation

The total momentum of a segment of the rim of length $L$ is $P = (m \cdot L) V(t)$, where $m(t)$ is the mass per unit length.

Using Newton's second law in the form of momentum conservation for a variable-mass system (where mass enters from a stationary state, meaning its initial velocity is $v = 0$):

$$\frac{d}{dt} \big( m(t) V(t) \big) = F_{\text{external}}$$

Expanding the derivative using the product rule:

$$m(t) \frac{dV}{dt} + V(t) \frac{dm}{dt} = 2\gamma$$

Recognizing that $\frac{dm}{dt} = \dot{m}(t) = \rho h V(t)$ is the mass flux entering the rim per unit length:

$$m(t) \dot{V}(t) + \rho h V(t)^2 = 2\gamma$$

### 4. Deduce the Expression for the Terminal Velocity

The problem states that **very quickly, the expansion velocity reaches a steady terminal velocity** $V$.

- **Steady-state condition:** At terminal velocity, the acceleration of the rim drops to zero ($\dot{V} = 0$), and the velocity settles to a constant value $V(t) \to V$.
    
- Substituting $\dot{V} = 0$ into our momentum conservation equation:
    
    $$0 + \rho h V^2 = 2\gamma$$
    
- Solving explicitly for the terminal expansion velocity $V$:
    
    $$V = \sqrt{\frac{2\gamma}{\rho h}}$$
    

This fundamental result is known as the **Culick-Taylor velocity** (or Culick's velocity), which dictates how fast a liquid sheet retracts after being ruptured.


# 7. — Examples of dynamical systems
### System 1: $\dot{x} = \lambda - x - \frac{1}{x}$

**1. Fixed Points and Stability**

To find the fixed points $x^*$, set $\dot{x} = 0$:

$$0 = \lambda - x - \frac{1}{x} \implies x^2 - \lambda x + 1 = 0$$

Using the quadratic formula, the fixed points are:

$$x^*_{1,2} = \frac{\lambda \pm \sqrt{\lambda^2 - 4}}{2}$$

Real roots exist only when the discriminant is non-negative: $\lambda^2 - 4 \ge 0$, which yields two regimes: $\lambda \ge 2$ or $\lambda \le -2$.

To determine stability, compute the derivative of the vector field $f(x) = \lambda - x - \frac{1}{x}$:

$$f'(x) = -1 + \frac{1}{x^2}$$

- **Regime $\lambda > 2$:** $x^*_1 = \frac{\lambda + \sqrt{\lambda^2 - 4}}{2} > 1$, so $f'(x^*_1) < 0$ (**Stable**). $x^*_2 = \frac{\lambda - \sqrt{\lambda^2 - 4}}{2} \in (0, 1)$, so $f'(x^*_2) > 0$ (**Unstable**).
    
- **Regime $\lambda < -2$:** $x^*_3 = \frac{\lambda - \sqrt{\lambda^2 - 4}}{2} < -1$, so $f'(x^*_3) < 0$ (**Stable**). $x^*_4 = \frac{\lambda + \sqrt{\lambda^2 - 4}}{2} \in (-1, 0)$, so $f'(x^*_4) > 0$ (**Unstable**).
    

**2. Interesting Values (Bifurcation Points)**

Bifurcations occur where the discriminant vanishes, meaning roots are created or destroyed. These points are **$\lambda_c = 2$** (with $x_c = 1$) and **$\lambda_c = -2$** (with $x_c = -1$).

**3. Bifurcation Diagram**

- **$\lambda < -2$:** Two branches in the $x < 0$ half-plane. A stable branch lies below $x = -1$, and an unstable branch lies between $x = -1$ and $x = 0$.
    
- **$\lambda = -2$:** The branches collide at $x = -1$.
    
- **$-2 < \lambda < 2$:** No fixed points exist.
    
- **$\lambda = 2$:** A new pair of fixed points is born at $x = 1$.
    
- **$\lambda > 2$:** Two branches in the $x > 0$ half-plane. A stable branch lies above $x = 1$, and an unstable branch lies between $x = 0$ and $x = 1$.
    

**4. Classification**

These are two separate **saddle-node bifurcations** occurring at $\lambda = \pm 2$.

**5. Taylor Expansion to Normal Form**

Expand around the bifurcation point $\lambda_c = 2, x_c = 1$. Introduce small perturbations $u$ and $\mu$ such that $x = 1 + u$ and $\lambda = 2 + \mu$:

$$\dot{u} = (2 + \mu) - (1 + u) - \frac{1}{1 + u}$$

Using the Taylor expansion for $\frac{1}{1+u} \approx 1 - u + u^2 - u^3$:

$$\dot{u} \approx 1 + \mu - u - (1 - u + u^2)$$

$$\dot{u} \approx \mu - u^2$$

This is precisely the analytical normal form of a saddle-node bifurcation.

### System 2: $\dot{x} = \left( \frac{2}{1+e^{-\lambda x}} - 1 \right) - x$

**1. Fixed Points and Stability** Simplify the term in brackets using hyperbolic functions:

$$\frac{2}{1+e^{-\lambda x}} - 1 = \frac{1 - e^{-\lambda x}}{1 + e^{-\lambda x}} = \tanh\left(\frac{\lambda x}{2}\right)$$

The system becomes $\dot{x} = \tanh\left(\frac{\lambda x}{2}\right) - x$.

Set $\dot{x} = 0 \implies x = \tanh\left(\frac{\lambda x}{2}\right)$.

- $x^* = 0$ is always a fixed point.
    
- The slope of $\tanh\left(\frac{\lambda x}{2}\right)$ at $x=0$ is $\frac{\lambda}{2}$. If $\frac{\lambda}{2} > 1$ (i.e., $\lambda > 2$), the curve crosses $y = x$ at two additional non-zero symmetric fixed points, $\pm x^*$.
    

For stability, compute $f'(x) = \frac{\lambda}{2} \text{sech}^2\left(\frac{\lambda x}{2}\right) - 1$:

- At $x^* = 0$: $f'(0) = \frac{\lambda}{2} - 1$.
    
    - If $\lambda < 2$, $f'(0) < 0$ (**Stable**).
        
    - If $\lambda > 2$, $f'(0) > 0$ (**Unstable**).
        
- At $\pm x^*$ (for $\lambda > 2$): Geometrically, the $\tanh$ curve must have a slope less than 1 when intersecting $y=x$ from above. Thus, $\frac{\lambda}{2} \text{sech}^2\left(\frac{\lambda x^*}{2}\right) < 1$, meaning $f'(\pm x^*) < 0$ (**Stable**).
    

**2. Interesting Values (Bifurcation Points)**

The bifurcation occurs where the stability of the origin changes: **$\lambda_c = 2$**.

**3. Bifurcation Diagram**

- **$0 < \lambda < 2$:** A single solid line along the axis $x = 0$ (stable).
    
- **$\lambda = 2$:** The origin loses stability.
    
- **$\lambda > 2$:** The axis $x = 0$ becomes dashed (unstable), and two new solid curves (stable) emerge symmetrically, opening to the right like a pitchfork.
    

**4. Classification**

This is a **supercritical pitchfork bifurcation**.

**5. Taylor Expansion to Normal Form**

Expand around $\lambda_c = 2, x_c = 0$. Let $\lambda = 2 + \mu$. The system is $\dot{x} = \tanh\left( x + \frac{\mu}{2}x \right) - x$.

Using the Taylor expansion $\tanh(z) \approx z - \frac{1}{3}z^3$:

$$\dot{x} \approx \left( x + \frac{\mu}{2}x \right) - \frac{1}{3}\left( x + \frac{\mu}{2}x \right)^3 - x$$

Keeping only the lowest-order terms in $x$ and $\mu$:

$$\dot{x} \approx \frac{\mu}{2}x - \frac{1}{3}x^3$$

By scaling the parameter ($\mu' = \mu/2$) and variable ($X = x/\sqrt{3}$), this matches the analytical normal form of a supercritical pitchfork bifurcation ($\dot{x} = r x - x^3$).

### System 3: $\dot{\phi} = -\lambda \sin\phi - 2 \sin 2\phi$

**1. Fixed Points and Stability** Rewrite using the double-angle identity: $\dot{\phi} = -\lambda \sin\phi - 4 \sin\phi \cos\phi = -\sin\phi (\lambda + 4 \cos\phi)$. Setting $\dot{\phi} = 0$ yields two sets of fixed points on the interval $[-\pi, \pi]$:

1. $\sin\phi = 0 \implies \phi^*_1 = 0$, $\phi^*_2 = \pi$ (and $-\pi$).
    
2. $\lambda + 4 \cos\phi = 0 \implies \cos\phi = -\frac{\lambda}{4}$. This yields solutions $\phi^*_{3,4} = \pm \arccos(-\frac{\lambda}{4})$ only if $\lambda \le 4$.
    

Compute $f'(\phi) = -\lambda \cos\phi - 4 \cos 2\phi = -\lambda \cos\phi - 8\cos^2\phi + 4$:

- **At $\phi^*_1 = 0$:** $f'(0) = -\lambda - 4$. Since $\lambda \ge 0$, this is always negative (**Stable**).
    
- **At $\phi^*_2 = \pi$:** $f'(\pi) = \lambda - 8(1) + 4 = \lambda - 4$.
    
    - If $\lambda < 4$, $f'(\pi) < 0$ (**Stable**).
        
    - If $\lambda > 4$, $f'(\pi) > 0$ (**Unstable**).
        
- **At $\phi^*_{3,4} = \pm \arccos(-\frac{\lambda}{4})$** (for $\lambda < 4$):
    
    $f'(\phi^*_{3,4}) = -\lambda(-\frac{\lambda}{4}) - 8(-\frac{\lambda}{4})^2 + 4 = \frac{\lambda^2}{4} - \frac{\lambda^2}{2} + 4 = 4 - \frac{\lambda^2}{4}$.
    
    Since $\lambda < 4$, $4 - \frac{\lambda^2}{4} > 0$ (**Unstable**).
    

**2. Interesting Values (Bifurcation Points)**

The critical parameter where roots collide and stability changes is **$\lambda_c = 4$**.

**3. Bifurcation Diagram**

Focusing on the boundary $\phi = \pi$ (which is continuous with $-\pi$ on a cylinder):

- **$0 \le \lambda < 4$:** $\phi = \pi$ is a solid line (stable). It is flanked symmetrically by two dashed curves representing the unstable fixed points $\pm \arccos(-\frac{\lambda}{4})$, which converge toward $\pi$ as $\lambda$ increases.
    
- **$\lambda = 4$:** The unstable roots collide with the stable root at $\pi$.
    
- **$\lambda > 4$:** The unstable roots disappear, and $\phi = \pi$ becomes a dashed line (unstable).
    

**4. Classification**

This is a **subcritical pitchfork bifurcation** occurring at $\phi = \pi$.

**5. Taylor Expansion to Normal Form**

Expand around the bifurcation point $\lambda_c = 4, \phi_c = \pi$. Let $\phi = \pi + \theta$ (where $\theta$ is small) and $\lambda = 4 + \mu$:

$$\dot{\theta} = -(4+\mu)\sin(\pi+\theta) - 2\sin(2\pi+2\theta)$$

$$\dot{\theta} = (4+\mu)\sin\theta - 2\sin(2\theta)$$

Use Taylor expansions $\sin\theta \approx \theta - \frac{\theta^3}{6}$ and $\sin(2\theta) \approx 2\theta - \frac{(2\theta)^3}{6} = 2\theta - \frac{4\theta^3}{3}$:

$$\dot{\theta} \approx (4+\mu)\left(\theta - \frac{\theta^3}{6}\right) - 2\left(2\theta - \frac{4\theta^3}{3}\right)$$

$$\dot{\theta} \approx 4\theta - \frac{2}{3}\theta^3 + \mu\theta - \mu\frac{\theta^3}{6} - 4\theta + \frac{8}{3}\theta^3$$

Keeping the dominant terms:

$$\dot{\theta} \approx \mu\theta + \left(\frac{8}{3} - \frac{2}{3}\right)\theta^3 = \mu\theta + 2\theta^3$$

This is the analytical normal form of a subcritical pitchfork bifurcation ($\dot{x} = r x + a x^3$ with $a > 0$).


![[Pasted image 20260925154656.png]]



# 8. — Logistic population model, with preys

### 1. Variables $B$ and $A$ in the Predation Model

The term $p(N) = \frac{B N^2}{A^2 + N^2}$ describes the feedback of the population $N$ on predation. This specific functional form is known in ecology as a Holling Type III functional response.

- **The variable $B$:** This models the **maximum predation rate** (or saturation level). As the prey population becomes infinitely large ($N \to \infty$), the predation rate asymptotically approaches $B$. It represents the maximum rate at which predators can consume prey, limited by handling time or predator satiation.
    
- **The variable $A$:** This models the **half-saturation constant**. When the prey population is exactly $N = A$, the predation rate is $p(A) = \frac{B A^2}{2A^2} = \frac{B}{2}$. It represents the prey density at which predators reach half of their maximum consumption rate, determining how quickly the predators "switch on" and focus on this specific prey as its density increases.
    

### 2. Making the Equations Dimensionless

To simplify the differential equation, we scale the population $N$ and time $t$ by characteristic values to create dimensionless variables $x$ and $\tau$. Because the predation term saturates based on the scale $A$, a natural choice for the population scale is $A$. Let:

$$x = \frac{N}{A} \implies N = A x$$

Substitute this into the governing equation:

$$\frac{d(Ax)}{dt} = \frac{R_0 - 1}{T} (Ax) \left(1 - \frac{Ax}{N_{tot}}\right) - \frac{B (Ax)^2}{A^2 + (Ax)^2}$$

$$A \frac{dx}{dt} = \frac{R_0 - 1}{T} A x \left(1 - \frac{A}{N_{tot}}x\right) - B \frac{x^2}{1 + x^2}$$

Divide the entire equation by $B$:

$$\frac{A}{B} \frac{dx}{dt} = \frac{(R_0 - 1) A}{B T} x \left(1 - \frac{A}{N_{tot}}x\right) - \frac{x^2}{1 + x^2}$$

Define the dimensionless time $\tau = \frac{B}{A} t$, which means $\frac{dx}{d\tau} = \frac{A}{B} \frac{dx}{dt}$. We can also group the remaining parameters into two dimensionless constants:

- $r = \frac{(R_0 - 1) A}{B T}$ (the dimensionless intrinsic growth rate)
    
- $q = \frac{N_{tot}}{A}$ (the dimensionless carrying capacity)
    

The fully dimensionless equation is:

$$\frac{dx}{d\tau} = r x \left(1 - \frac{x}{q}\right) - \frac{x^2}{1 + x^2}$$

### 3. Stability of the Fixed Point $N = 0$

In the dimensionless system, $N=0$ corresponds to the fixed point $x^* = 0$. Let the vector field be $f(x) = r x (1 - x/q) - \frac{x^2}{1 + x^2}$.

To find the linear stability, we evaluate the derivative $f'(x)$ at $x = 0$:

$$f'(x) = r - \frac{2r x}{q} - \frac{2x(1+x^2) - 2x^3}{(1+x^2)^2}$$

$$f'(0) = r$$

Returning to the original parameters, the eigenvalue determining stability is $r = \frac{(R_0 - 1) A}{B T}$.

- If **$R_0 > 1$**, then $r > 0$. The fixed point $N=0$ is **unstable**. A small initial population will grow.
    
- If **$R_0 < 1$**, then $r < 0$. The fixed point $N=0$ is **stable**. The population cannot sustain itself and goes extinct.
    

### 4. Analysis of Non-Zero Fixed Points and Stability

To find the non-zero fixed points, we set $f(x) = 0$ and divide by $x$ (since we already analyzed $x=0$):

$$r \left(1 - \frac{x}{q}\right) - \frac{x}{1 + x^2} = 0 \implies r \left(1 - \frac{x}{q}\right) = \frac{x}{1 + x^2}$$

This is solved graphically by finding the intersections of a line $y_1(x) = r(1 - x/q)$ and a curve $y_2(x) = \frac{x}{1 + x^2}$.

- **The curve $y_2(x)$** starts at $(0,0)$, reaches a peak maximum of $1/2$ at $x=1$, and decays to $0$ as $x \to \infty$.
    
- **The line $y_1(x)$** has a y-intercept of $r$ and an x-intercept of $q$.
    

Depending on the values of $r$ and $q$, there are two main qualitative regimes:

1. **One Intersection:** If the line is very steep (low carrying capacity $q$) or very high (high growth rate $r$), it intersects the curve exactly once at some point $x^*$. Because the line crosses the curve from above to below, $y_1'(x^*) < y_2'(x^*)$, meaning $f'(x^*) = y_1' - y_2' < 0$. This single fixed point is **stable**.
    
2. **Three Intersections (Bistability):** For intermediate values of $r$ and large carrying capacities $q$, the line can intersect the curve three times at $x_1 < x_2 < x_3$.
    
    - At $x_1$ (low population), the line crosses from above to below $\implies$ **Stable**.
        
    - At $x_2$ (intermediate population), the line crosses from below to above $\implies$ **Unstable** (acts as a threshold).
        
    - At $x_3$ (high population), the line crosses from above to below $\implies$ **Stable**.
        

### 5. Physical Content of the Model

This is a classic ecological model (most notably applied to the Spruce Budworm outbreak by Ludwig, Jones, and Holling). It captures the non-linear interaction between a population undergoing resource-limited logistic growth and a predator population with a limited consumption capacity.

The most important physical consequence of this model is **bistability and hysteresis**.

- **Endemic State ($x_1$):** A low-population refuge state maintained by the predators. The predators can eat the prey fast enough to keep the population down.
    
- **Epidemic/Outbreak State ($x_3$):** A high-population state limited only by the environment's carrying capacity ($N_{tot}$). The prey population has grown so large that it overwhelms the predators' capacity to eat them (predator satiation).
    

If the environment slowly changes (e.g., the carrying capacity $q$ increases due to forest maturation), the system can undergo a saddle-node bifurcation. The stable lower state $x_1$ collides with the unstable threshold $x_2$ and disappears. This causes a sudden, catastrophic population explosion (an outbreak) where the system rapidly shifts to the high state $x_3$. Even if the environment later returns to its original condition, the population will remain stuck in the outbreak state until conditions degrade drastically enough to destroy the $x_3$ state, demonstrating hysteresis.



![[Pasted image 20260925160101.png]]



# 9. — Lotka and Volterra prey-predator model
### 1. Making the Equations Dimensionless

The Lotka-Volterra model describes the continuous interaction between a prey population $H$ and a predator population $P$ using the following coupled equations:

$$\dot{H} = (r - aP)H$$

$$\dot{P} = (\gamma aH - \mu)P$$

To make the system dimensionless, we identify the characteristic scales of the system from the parameters. The prey population naturally scales with the predator death-to-conversion ratio ($\mu / \gamma a$), and the predator population scales with the prey growth-to-predation ratio ($r / a$).

We introduce dimensionless variables $x$, $y$, and $\tau$:

- **Scaled prey:** $x = \frac{\gamma a}{\mu} H$
    
- **Scaled predators:** $y = \frac{a}{r} P$
    
- **Scaled time:** $\tau = r t$
    

Substituting these into the original prey equation:

$$\frac{d}{d(\tau/r)} \left( \frac{\mu}{\gamma a} x \right) = \left( r - a \left( \frac{r}{a} y \right) \right) \left( \frac{\mu}{\gamma a} x \right)$$

$$r \frac{\mu}{\gamma a} \frac{dx}{d\tau} = r(1 - y) \frac{\mu}{\gamma a} x \implies \frac{dx}{d\tau} = x(1 - y)$$

Substituting into the predator equation:

$$\frac{d}{d(\tau/r)} \left( \frac{r}{a} y \right) = \left( \gamma a \left( \frac{\mu}{\gamma a} x \right) - \mu \right) \left( \frac{r}{a} y \right)$$

$$r \frac{r}{a} \frac{dy}{d\tau} = \mu(x - 1) \frac{r}{a} y \implies \frac{dy}{d\tau} = \frac{\mu}{r} y(x - 1)$$

By defining the single dimensionless parameter $\alpha = \frac{\mu}{r}$ (the ratio of predator death rate to prey growth rate), the fully dimensionless system is:

$$\frac{dx}{d\tau} = x(1 - y)$$

$$\frac{dy}{d\tau} = \alpha y(x - 1)$$

### 2. Determining the Fixed Points

Fixed points occur where the populations are in a steady state, meaning their time derivatives are zero:

$$\frac{dx}{d\tau} = x(1 - y) = 0$$

$$\frac{dy}{d\tau} = \alpha y(x - 1) = 0$$

Solving this system of equations yields two distinct fixed points:

1. **The Extinction State:** If $x = 0$, the second equation requires $y = 0$. This gives the trivial fixed point **$(x^*, y^*) = (0, 0)$**.
    
2. **The Coexistence State:** If $x \neq 0$, the first equation requires $1 - y = 0 \implies y = 1$. Substituting $y = 1$ into the second equation forces $x - 1 = 0 \implies x = 1$. This gives the non-trivial fixed point **$(x^*, y^*) = (1, 1)$**.
    

In the original physical dimensions, these fixed points correspond to $(H, P) = (0, 0)$ and $(H, P) = \left( \frac{\mu}{\gamma a}, \frac{r}{a} \right)$.

### 3. Showing that Trajectories are Periodic

To prove the trajectories are periodic, we must show that the system possesses a conserved quantity (a first integral) that forms closed loops in the phase plane.

Divide the predator equation by the prey equation to eliminate time $\tau$:

$$\frac{dy}{dx} = \frac{\alpha y(x - 1)}{x(1 - y)}$$

This is a separable differential equation. Group the $y$ terms on the left and the $x$ terms on the right:

$$\frac{1 - y}{y} dy = \alpha \frac{x - 1}{x} dx$$

$$\left( \frac{1}{y} - 1 \right) dy = \alpha \left( 1 - \frac{1}{x} \right) dx$$

Integrate both sides:

$$\int \left( \frac{1}{y} - 1 \right) dy = \int \alpha \left( 1 - \frac{1}{x} \right) dx$$

$$\ln y - y = \alpha (x - \ln x) + C$$

Rearrange to define a conserved "energy" function $V(x,y)$:

$$V(x,y) = \alpha(x - \ln x) + (y - \ln y) = \text{Constant}$$

Consider the mathematical properties of the function $f(z) = z - \ln z$ for $z > 0$:

- As $z \to 0$, $-\ln z \to \infty$, so $f(z) \to \infty$.
    
- As $z \to \infty$, $z$ dominates, so $f(z) \to \infty$.
    
- The function has a strict global minimum at $z = 1$, where $f(1) = 1$.
    

Because both $x - \ln x$ and $y - \ln y$ form "bowls" that go to infinity at the boundaries ($0$ and $\infty$) and have unique minimums at $1$, the total function $V(x,y)$ forms a 3D bowl with a global minimum exactly at the coexistence fixed point $(1,1)$.

Any constant energy level $V(x,y) = C > V(1,1)$ corresponds to a horizontal slice of this bowl, which forms a strictly closed, bounded contour. Since the state of the system is trapped on one of these closed contours, the population dynamics perfectly repeat themselves, proving the trajectories are periodic.

### 4. What Selects the Average Population of Preys and Predators?

While the precise amplitude of the oscillations is determined entirely by the initial conditions, the **time-averaged populations** are selected exclusively by the inherent kinetic parameters of the system (growth rates, death rates, and predation efficiency), independent of the starting state.

To prove this, consider the exact time-average over one full period $T$. Start by rewriting the dimensionless equations logarithmically:

$$\frac{1}{x} \frac{dx}{d\tau} = 1 - y \implies \frac{d}{d\tau}(\ln x) = 1 - y$$

$$\frac{1}{y} \frac{dy}{d\tau} = \alpha(x - 1) \implies \frac{d}{d\tau}(\ln y) = \alpha(x - 1)$$

Integrate the prey equation over one period $T$:

$$\int_0^T \frac{d}{d\tau}(\ln x) d\tau = \int_0^T (1 - y) d\tau$$

$$\ln x(T) - \ln x(0) = T - \int_0^T y d\tau$$

Because the trajectory is periodic, the population returns exactly to its starting value, meaning $x(T) = x(0)$. The left side evaluates to zero:

$$0 = T - \int_0^T y d\tau \implies \frac{1}{T} \int_0^T y d\tau = 1$$

The left side is the mathematical definition of the time average, denoted $\langle y \rangle$. We find that $\langle y \rangle = 1$.

Applying the exact same integration to the predator equation yields:

$$0 = \alpha \int_0^T x d\tau - \alpha T \implies \frac{1}{T} \int_0^T x d\tau = 1 \implies \langle x \rangle = 1$$

The time-averaged populations always perfectly equal the coordinates of the non-trivial fixed point. In physical units, the average populations are:

- **Average Prey:** $\langle H \rangle = \frac{\mu}{\gamma a}$
    
- **Average Predator:** $\langle P \rangle = \frac{r}{a}$
    

Thus, the average predator population is selected solely by the prey's reproduction rate $r$ and the predator's hunting efficiency $a$, while the average prey population is selected by the predator's mortality $\mu$ and conversion efficiency $\gamma$.

![[Pasted image 20260925160644.png]]



# 10. — An example of bifurcation

### 1. The Dimensionless Parameter and Qualitative Behavior

To find the dimensionless parameter, we compare the two competing forces in the system: gravity (pulling the ball down to $\theta = 0$) and the centrifugal force (pushing the ball outwards to $\theta = \pi/2$).

The maximum gravitational force component along the ring is of order $mg$. The maximum centrifugal force in the rotating frame is $m\omega^2R$.

We define the dimensionless parameter $\Omega$ as the ratio of these forces:

$$\Omega = \frac{\omega^2 R}{g}$$

**Qualitative Description:**

- **For small $\Omega$ (slow rotation):** Gravity dominates. The centrifugal force is too weak to push the ball up the ring, so the only stable equilibrium position is at the very bottom ($\theta = 0$).
    
- **For large $\Omega$ (fast rotation):** The centrifugal force overpowers gravity near the bottom. The ball is pushed outwards and upwards until the tangential component of gravity perfectly balances the tangential component of the centrifugal force. The bottom ($\theta = 0$) becomes unstable, and two new stable symmetric equilibria emerge at $\theta = \pm \theta^*$.
    

### 2. Mechanical Energy and Equations of Motion

In the rotating frame of reference, the system is conservative if we introduce the centrifugal potential energy.

- **Kinetic Energy ($T$):** The only motion in the rotating frame is along the ring.
    
    $$T = \frac{1}{2}m (R\dot{\theta})^2 = \frac{1}{2}mR^2\dot{\theta}^2$$
    
- **Effective Potential Energy ($V_{eff}$):** This consists of gravitational potential energy and the fictitious centrifugal potential energy. Taking the center of the ring as the zero-point for gravity, and knowing the distance from the rotation axis is $r = R\sin\theta$:
    
    $$V_{eff}(\theta) = -mgR\cos\theta - \frac{1}{2}m\omega^2(R\sin\theta)^2$$
    
    $$V_{eff}(\theta) = -mgR\cos\theta - \frac{1}{2}m\omega^2R^2\sin^2\theta$$
    

The total mechanical energy in the rotating frame is $E = T + V_{eff}$. Since energy is conserved ($\dot{E} = 0$):

$$\frac{d}{dt} \left[ \frac{1}{2}mR^2\dot{\theta}^2 - mgR\cos\theta - \frac{1}{2}m\omega^2R^2\sin^2\theta \right] = 0$$

$$mR^2\dot{\theta}\ddot{\theta} + mgR\sin\theta\dot{\theta} - m\omega^2R^2\sin\theta\cos\theta\dot{\theta} = 0$$

Dividing by $mR^2\dot{\theta}$ (assuming $\dot{\theta} \neq 0$), we get the equation of motion:

$$\ddot{\theta} + \frac{g}{R}\sin\theta - \omega^2\sin\theta\cos\theta = 0$$

### 3. Equilibria, Stability, and the Bifurcation Diagram

Equilibrium positions occur where the angular acceleration is zero ($\ddot{\theta} = 0$), which corresponds to the extrema of the effective potential ($\partial V_{eff} / \partial \theta = 0$):

$$\sin\theta \left( \frac{g}{R} - \omega^2\cos\theta \right) = 0$$

This gives two sets of solutions:

1. $\sin\theta = 0 \implies \theta_0 = 0$ (the bottom) and $\theta_\pi = \pi$ (the top).
    
2. $\cos\theta = \frac{g}{\omega^2 R} = \frac{1}{\Omega} \implies \theta^* = \pm \arccos\left(\frac{1}{\Omega}\right)$. Since $\cos\theta \le 1$, these solutions only exist physically when $\Omega \ge 1$.
    

**Stability:**

Stability is determined by the second derivative of the effective potential, $V_{eff}''(\theta) = mR^2 \left( \frac{g}{R}\cos\theta - \omega^2(\cos^2\theta - \sin^2\theta) \right)$.

- **At $\theta_0 = 0$:** $V_{eff}''(0) = mgR(1 - \Omega)$.
    
    - Stable for $\Omega < 1$.
        
    - Unstable for $\Omega > 1$.
        
- **At $\theta_\pi = \pi$:** $V_{eff}''(\pi) = mgR(-1 - \Omega)$. This is always negative, so the top is always an unstable equilibrium.
    
- __At $\theta^* =arccos(1/\Omega)$(for $\Omega > 1$):
    
    $V_{eff}''(\theta^*) = mgR \left( \frac{1}{\Omega} - \Omega\left(\frac{1}{\Omega^2} - \left(1 - \frac{1}{\Omega^2}\right)\right) \right) = mgR \left( \Omega - \frac{1}{\Omega} \right)$.
    
    Since this equilibrium only exists when $\Omega > 1$, this second derivative is always strictly positive. Thus, when they exist, these new positions are stable.
    

**Bifurcation Diagram:**

The critical parameter is **$\Omega_c = 1$**.

If you graph the equilibrium angle $\theta_{eq}$ on the y-axis against $\Omega$ on the x-axis:

- For $\Omega < 1$, there is a single solid line along $\theta = 0$ (stable).
    
- At $\Omega = 1$, the stable line at $\theta=0$ becomes dashed (unstable), and two new solid branches diverge symmetrically upwards and downwards, following $\theta = \pm \arccos(1/\Omega)$. This shape perfectly defines a **supercritical pitchfork bifurcation**.
    

### 4. Period of Small Oscillations

For small oscillations around the bottom equilibrium ($\theta = 0$), we can use the small-angle approximations $\sin\theta \approx \theta$ and $\cos\theta \approx 1$.

Substituting these into the equation of motion gives:

$$\ddot{\theta} + \left( \frac{g}{R} - \omega^2 \right)\theta = 0$$

$$\ddot{\theta} + \frac{g}{R}(1 - \Omega)\theta = 0$$

This is the equation of a simple harmonic oscillator with an effective angular frequency $\omega_{eff} = \sqrt{\frac{g}{R}(1 - \Omega)}$. The period $T$ is:

$$T = \frac{2\pi}{\omega_{eff}} = 2\pi \sqrt{\frac{R}{g(1 - \Omega)}}$$

**Behavior in the vicinity of $\Omega_c$:**

As $\Omega \to 1^-$, the denominator approaches zero, causing the period to diverge: $T \to \infty$. In complex systems and statistical physics, this divergence of the relaxation time scale as you approach a critical point is known as **critical slowing down**. The restoring force approaches zero, making the system incredibly sluggish in returning to equilibrium.

### 5. Non-linear Expansion and the Normal Form

To find the normal form of the bifurcation, we expand the equation of motion to the lowest non-linear order (up to $O(\theta^3)$) around $\theta=0$.

Using the Taylor series: $\sin\theta \approx \theta - \frac{\theta^3}{6}$ and $\cos\theta \approx 1 - \frac{\theta^2}{2}$.

The term $\sin\theta\cos\theta = \frac{1}{2}\sin(2\theta) \approx \frac{1}{2}\left(2\theta - \frac{(2\theta)^3}{6}\right) = \theta - \frac{2}{3}\theta^3$.

Substitute these into the equation of motion:

$$\ddot{\theta} + \frac{g}{R}\left(\theta - \frac{\theta^3}{6}\right) - \omega^2\left(\theta - \frac{2}{3}\theta^3\right) = 0$$

$$\ddot{\theta} = \left(\omega^2 - \frac{g}{R}\right)\theta + \left(\frac{g}{6R} - \frac{2\omega^2}{3}\right)\theta^3$$

To write this in a normalized form, we factor out $g/R$ and use our dimensionless parameter $\Omega$:

$$\ddot{\theta} = \frac{g}{R}(\Omega - 1)\theta - \frac{g}{R}\left(\frac{2\Omega}{3} - \frac{1}{6}\right)\theta^3$$

Let's introduce a dimensionless time scale $\tau = t\sqrt{g/R}$ to absorb the $g/R$ prefactors. The equation becomes:

$$\theta'' = (\Omega - 1)\theta - \left(\frac{2\Omega}{3} - \frac{1}{6}\right)\theta^3$$

_(where the primes denote derivatives with respect to $\tau$)_

Because the bifurcation occurs right at $\Omega_c = 1$, we can evaluate the coefficient of the non-linear $\theta^3$ term strictly at the critical point to capture the lowest-order structural dynamics. At $\Omega = 1$, the cubic coefficient is $\left(\frac{2(1)}{3} - \frac{1}{6}\right) = \frac{1}{2}$.

This gives us the classic **normal form for a supercritical pitchfork bifurcation**:

$$\theta'' = (\Omega - 1)\theta - \frac{1}{2}\theta^3$$
# 11. — A simplistic model of gravitary buckling
To analyze this gravitary buckling model, we must first recognize that a bifurcation requires competing forces. If the beam were hanging downwards, both gravity and the spring would work together to pull it to $\theta = 0$, resulting in an unconditionally stable system. For buckling to occur, the unperturbed state must be pointing vertically upwards (an inverted pendulum), where the restoring torque of the spring competes against the destabilizing torque of gravity.

Let $\theta = 0$ represent the vertical upright position.

The potential energy of the system $V(\theta)$ consists of the elastic potential energy of the torsion spring and the gravitational potential energy. Setting the pivot $O$ as the reference height, the mass is at a height $h = L\cos\theta$:

$$V(\theta) = \frac{1}{2}k\theta^2 + mgL\cos\theta$$

The generalized force (torque) $\tau$ driving the system is the negative derivative of the potential energy:

$$\tau = -\frac{dV}{d\theta} = -k\theta + mgL\sin\theta$$

### 1. Equilibrium Positions

Equilibrium occurs when the net torque acting on the system is zero ($\tau = 0$), which corresponds to the extrema of the potential energy ($dV/d\theta = 0$):

$$k\theta = mgL\sin\theta$$

$$\frac{\sin\theta}{\theta} = \frac{k}{mgL}$$

Let's define a critical mass $m_c = \frac{k}{gL}$. The equilibrium condition becomes:

$$\frac{\sin\theta}{\theta} = \frac{m_c}{m}$$

Because the function $f(\theta) = \frac{\sin\theta}{\theta}$ has a maximum value of $1$ at $\theta = 0$ and strictly decreases for $\theta \in (0, \pi]$, the number of solutions depends entirely on the ratio of $m_c$ to the variable mass $m$:

- **Trivial Solution:** The equation is always satisfied at $\theta_0 = 0$, regardless of the mass.
    
- **Buckled Solutions:** If $m > m_c$, the ratio $m_c/m$ is less than $1$. The curve $\sin\theta$ intersects the line $\frac{m_c}{m}\theta$ at three points: the origin, and two non-trivial symmetric angles $\theta^* = \pm \theta(m)$.
    

### 2. Stability of Solutions

Stability is determined by the concavity of the potential energy at the equilibrium points, given by the second derivative:

$$V''(\theta) = k - mgL\cos\theta$$

**Stability of the upright position ($\theta_0 = 0$):**

$$V''(0) = k - mgL = gL(m_c - m)$$

- When **$m < m_c$**, $V''(0) > 0$. The upright position is a potential minimum and is **stable**. The spring is stiff enough to counteract gravity.
    
- When **$m > m_c$**, $V''(0) < 0$. The upright position becomes a potential maximum and is **unstable**. Gravity overcomes the spring stiffness.
    

__Stability of the buckled positions ($\theta^*= \theta(m)$):

These solutions only exist when $m > m_c$. At these points, we know $k\theta^* = mgL\sin\theta^*$.

We want to evaluate $V''(\theta^*) = k - mgL\cos\theta^*$.

Since $\theta^*$ is non-zero (and in the range where $\theta < \pi$), we can use the geometric inequality $\cos\theta^* < \frac{\sin\theta^*}{\theta^*}$.

Substituting the equilibrium condition $\frac{\sin\theta^*}{\theta^*} = \frac{k}{mgL}$ into the inequality gives:

$$\cos\theta^* < \frac{k}{mgL} \implies mgL\cos\theta^* < k$$

Therefore, $V''(\theta^*) = k - mgL\cos\theta^* > 0$.

Whenever the buckled states exist ($m > m_c$), they are strictly **stable**.

### 3. Bifurcation Diagram

To trace the equilibrium branches, we plot the equilibrium angle $\theta_{eq}$ as a function of the control parameter, the variable mass $m$.

- For $m < m_c$, there is a single stable branch lying strictly on the horizontal axis ($\theta = 0$).
    
- At the critical point $m = m_c$, this trivial solution loses stability (marked by a dashed line for $m > m_c$).
    
- For $m > m_c$, two new stable branches emerge symmetrically, bending outwards from the axis.
    

To find the shape of these branches near the onset of buckling, we expand the equilibrium condition using the Taylor series $\sin\theta \approx \theta - \frac{\theta^3}{6}$:

$$k\theta \approx mgL\left(\theta - \frac{\theta^3}{6}\right)$$

$$\theta^2 \approx 6\left(1 - \frac{k}{mgL}\right) = 6\left(\frac{m - m_c}{m}\right)$$

Very close to the transition ($m \approx m_c$), the amplitude of the buckling angle grows as $\theta^* \propto \pm\sqrt{m - m_c}$. This square-root scaling reveals a classic **supercritical pitchfork bifurcation**, functionally identical to the spontaneous symmetry breaking seen in continuous phase 
transitions where a mean-field order parameter emerges below a critical temperature.

![[Pasted image 20261004021838.png]]

# 12. — Stick-slip

**1. Dynamical equations** Let $x(t)$ be the position of the mass $m$, and assume the free end of the spring is pulled at a constant velocity $V$, so its position is $Vt$. The elongation of the spring is defined as $\epsilon(t) = Vt - x(t)$. The restoring force of the spring is $F_T = K\epsilon(t)$. The motion is governed by two distinct states depending on the solid friction:

- **Stick phase ("collé"):** The mass is stationary ($\dot{x} = 0$, $\ddot{x} = 0$). The spring elongates at the rate of the pulling velocity ($\dot{\epsilon} = V$). The static friction force exactly balances the spring force ($R_T = K\epsilon$) until the spring force exceeds the maximum static friction threshold, $\mu_s mg$. The equation of motion is simply:
    
    $$0 = K(Vt - x) - R_T$$
    
- **Slip phase ("glissé"):** Once $K\epsilon > \mu_s mg$, the mass breaks free and slides. The kinetic friction force opposes the motion and is given by $R_T = \mu(\dot{x})mg$, where $\mu(\dot{x})$ is the velocity-dependent friction coefficient. Applying Newton's second law in the direction of motion ($\dot{x} > 0$):
    
    $$m\ddot{x} = K(Vt - x) - \mu(\dot{x})mg$$
    

**2. Steady state solutions**

A steady-state solution implies that the mass slides continuously without acceleration, perfectly matching the pulling velocity. Therefore, $\ddot{x}_{ss} = 0$ and $\dot{x}_{ss} = V$.

Substituting these conditions into the slip phase dynamical equation gives:

$$0 = K(Vt - x_{ss}) - \mu(V)mg$$

Solving for the steady-state position $x_{ss}(t)$:

$$x_{ss}(t) = Vt - \frac{\mu(V)mg}{K}$$

Equivalently, this means the steady-state spring elongation $\epsilon_{ss}$ remains constant at $\epsilon_{ss} = \frac{\mu(V)mg}{K}$.

**3. Stability of these solutions**

To analyze the stability of the steady sliding state, we introduce a small perturbation $\delta x(t)$ to the steady-state position:

$$x(t) = x_{ss}(t) + \delta x(t)$$

This implies the velocity is $\dot{x}(t) = V + \delta\dot{x}(t)$ and the acceleration is $\ddot{x}(t) = \delta\ddot{x}(t)$.

We perform a first-order Taylor expansion of the friction coefficient around the steady velocity $V$:

$$\mu(V + \delta\dot{x}) \approx \mu(V) + \mu'(V)\delta\dot{x}$$

_(where $\mu'(V) = \frac{d\mu}{d\dot{x}}\Big\vert{}_V$)_

Substitute the perturbed variables and the expanded friction term into the dynamical equation:

$$m\delta\ddot{x} = K[Vt - (x_{ss} + \delta x)] - [\mu(V) + \mu'(V)\delta\dot{x}]mg$$

We can use the steady-state balance equation ($K(Vt - x_{ss}) = \mu(V)mg$) to cancel out the zeroth-order terms:

$$m\delta\ddot{x} = - K\delta x - \mu'(V)mg\delta\dot{x}$$

Rearranging this yields the equation of a damped harmonic oscillator for the perturbation:

$$m\delta\ddot{x} + [\mu'(V)mg]\delta\dot{x} + K\delta x = 0$$

The stability of the system is entirely determined by the sign of the damping term, $\mu'(V)mg$:

- **Stable:** If $\mu'(V) > 0$ (velocity-strengthening friction), the damping coefficient is positive. Any small perturbation will exponentially decay over time, and the mass will return to steady sliding at velocity $V$.
    
- **Unstable:** If $\mu'(V) < 0$ (velocity-weakening friction), the damping coefficient is negative. Any small perturbation will be exponentially amplified. The steady sliding state is unstable, forcing the system into the cyclic stick-slip oscillations depicted in the elongation graph.

# 13. — A quadratic iterated map
### 1. Fixed Points

A fixed point $u^*$ remains unchanged under the application of the map, meaning it must satisfy the equation $u^* = (u^*)^2 + \mu$. Rearranging this gives a standard quadratic equation:

$$(u^*)^2 - u^* + \mu = 0$$

Using the quadratic formula, the fixed points are:

$$u^* = \frac{1 \pm \sqrt{1 - 4\mu}}{2}$$

Real fixed points exist only when the discriminant is non-negative ($1 - 4\mu \ge 0$), which occurs for **$\mu \le 1/4$**.

Let us define these two branches as:

- $u_+ = \frac{1 + \sqrt{1 - 4\mu}}{2}$
    
- $u_- = \frac{1 - \sqrt{1 - 4\mu}}{2}$
    

### 2. Stability and Bifurcations

The stability of a fixed point $u^*$ is determined by the absolute value of the map's derivative evaluated at that point. For the map $f(u) = u^2 + \mu$, the derivative is $f'(u) = 2u$. A fixed point is stable if $\vert{}f'(u^*)\vert{} < 1$, which translates to $-1/2 < u^* < 1/2$.

**Evaluating $u_+$:**

The derivative is $f'(u_+) = 1 + \sqrt{1 - 4\mu}$. Since the square root is non-negative, $f'(u_+) \ge 1$ for all valid $\mu$. Therefore, **$u_+$ is always unstable**.

**Evaluating $u_-$:**

The derivative is $f'(u_-) = 1 - \sqrt{1 - 4\mu}$. For stability, we require:

$$-1 < 1 - \sqrt{1 - 4\mu} < 1$$

Solving this compound inequality:

- $1 - \sqrt{1 - 4\mu} < 1 \implies -\sqrt{1 - 4\mu} < 0 \implies 1 - 4\mu > 0 \implies \mu < 1/4$
    
- $-1 < 1 - \sqrt{1 - 4\mu} \implies \sqrt{1 - 4\mu} < 2 \implies 1 - 4\mu < 4 \implies -3 < 4\mu \implies \mu > -3/4$
    

Therefore, **$u_-$ is stable for $-3/4 < \mu < 1/4$**.

**Bifurcations:**

- **At $\mu = 1/4$:** The two fixed points merge at $u^* = 1/2$, where $f'(1/2) = 1$. This is a **saddle-node (or tangent) bifurcation**, where a pair of fixed points (one stable, one unstable) is created as $\mu$ decreases through $1/4$.
    
- **At $\mu = -3/4$:** The fixed point $u_-$ reaches $-1/2$, where $f'(-1/2) = -1$. This indicates a **period-doubling (or flip) bifurcation**. As $\mu$ drops below $-3/4$, $u_-$ loses its stability, and a stable period-2 cycle emerges.
    

### 3. Range for a Stable 2-cycle

A 2-cycle consists of two distinct points, $p$ and $q$, where $f(p) = q$ and $f(q) = p$. These points are fixed points of the second iterate of the map, $f^{(2)}(u) = f(f(u))$:

$$f^{(2)}(u) = (u^2 + \mu)^2 + \mu = u^4 + 2\mu u^2 + \mu^2 + \mu$$

To find these points, we solve $f^{(2)}(u) - u = 0$:

$$u^4 + 2\mu u^2 - u + \mu^2 + \mu = 0$$

Because the period-1 fixed points also satisfy this equation, the polynomial $(u^2 - u + \mu)$ must be a factor. Factoring it out using polynomial division gives:

$$(u^2 - u + \mu)(u^2 + u + \mu + 1) = 0$$

The 2-cycle points $p$ and $q$ are the roots of the second factor:

$$u^2 + u + (\mu + 1) = 0$$

These roots exist when the discriminant is positive: $1 - 4(\mu + 1) > 0 \implies \mu < -3/4$. From Vieta's formulas, we know that for these points $p$ and $q$, the product is $pq = \mu + 1$.

The stability of the 2-cycle depends on the derivative of the second iterate map at $p$ (or $q$). Using the chain rule:

$$(f^{(2)})'(p) = f'(f(p)) \cdot f'(p) = f'(q) \cdot f'(p) = (2q)(2p) = 4pq$$

Substitute $pq = \mu + 1$ into the stability condition $\vert{}(f^{(2)})'(p)\vert{} < 1$:

$$-1 < 4(\mu + 1) < 1$$

$$-1/4 < \mu + 1 < 1/4$$

$$-5/4 < \mu < -3/4$$

Thus, there is a stable 2-cycle for the control parameter range **$\mu \in (-1.25, -0.75)$**. At $\mu = -5/4$, another period-doubling bifurcation occurs, leading to a 4-cycle.

# 14. — Logistic map

**1. Fixed points of $f$** A fixed point $x^*$ satisfies the equation $x^* = f(x^*)$. Using the given function $f(x) = ax(1-x)$:

$$x = ax(1-x)$$

$$x - ax + ax^2 = 0$$

$$x(1 - a + ax) = 0$$

This yields two fixed points:

- The trivial fixed point: __$x_0 = 0$
    
- The non-trivial fixed point: $1 - a + ax = 0 \implies$ $x_1 = 1 - \frac{1}{a}$


**2. Stability of the fixed points**

Stability is determined by evaluating the derivative $f'(x) = a - 2ax$ at the fixed points. A fixed point is stable if $\vert{}f'(x^*)\vert{} < 1$.

- For $x_0 = 0$:
    
    $f'(0) = a$.
    
    This fixed point is stable for **$-1 < a < 1$**. (In standard population dynamics contexts where $a \ge 0$, it is stable for $0 \le a < 1$).
    
- __For $x_1 = 1 - \frac{1}{a}$:
    
    $f'(x_1^*) = a - 2a\left(1 - \frac{1}{a}\right) = a - 2a + 2 = 2 - a$.
    
    For stability, we require $\vert{}2 - a\vert{} < 1$:
    
    $$-1 < 2 - a < 1 \implies -3 < -a < -1 \implies 1 < a < 3$$
    
    This fixed point is stable for **$1 < a < 3$**.
    

**3. Finding $g(x)$** The second iterate map is defined as $x_{n+2} = g(x_n)$ where $g(x) = f(f(x))$.

$$g(x) = a(f(x))(1 - f(x))$$

$$g(x) = a(ax(1-x))(1 - ax(1-x))$$

$$g(x) = a^2 x(1-x)(1 - ax + ax^2)$$

**4. Equation obeyed by the fixed points of $g$** The fixed points of $g$ must satisfy $x = g(x)$. Setting this up gives:

$$x = a^2 x(1-x)(1 - ax + ax^2)$$

$$x [1 - a^2(1-x)(1 - ax + ax^2)] = 0$$

Expanding the bracketed polynomial yields a quartic equation in $x$:

$$-a^3 x^4 + 2a^3 x^3 - a^2(a+1) x^2 + (a^2 - 1) x = 0$$

**5. Fixed points of $g$**

The fixed points of $f$ ($x_0^*$ and $x_1^*$) are inherently fixed points of $g$ (since if $f(x)=x$, then $f(f(x))=x$). We can factor out these known roots from the quartic equation $x - g(x) = 0$.

Factoring out $x$ (which corresponds to $x_0^* = 0$) and then using polynomial division to factor out $(ax - a + 1)$ (which corresponds to the root $x_1^* = 1 - 1/a$), we are left with a quadratic factor:

$$a^2 x^2 - a(a+1)x + (a+1) = 0$$

The roots of this quadratic represent the true period-2 cycle. Using the quadratic formula:

$$x_\pm = \frac{a(a+1) \pm \sqrt{a^2(a+1)^2 - 4a^2(a+1)}}{2a^2}$$

$$x_\pm = \frac{a+1 \pm \sqrt{(a+1)(a-3)}}{2a}$$

These two new fixed points of $g$ (the period-2 cycle) emerge when the term under the square root is positive, which occurs for **$a > 3$**.

**6. Stability of the fixed points of $g$**

To study the stability of the 2-cycle $x_\pm$, we evaluate $g'(x)$. By the chain rule, $g'(x_\pm) = f'(f(x_\pm)) \cdot f'(x_\pm)$. Because $x_+$ and $x_-$ map to each other, this simplifies to $g'(x_\pm) = f'(x_-) \cdot f'(x_+)$.

$$g'(x_\pm) = a(1 - 2x_-) \cdot a(1 - 2x_+) = a^2 [1 - 2(x_+ + x_-) + 4(x_+ x_-)]$$

From Vieta's formulas on the quadratic $a^2 x^2 - a(a+1)x + (a+1) = 0$, we know the sum of the roots is $x_+ + x_- = \frac{a+1}{a}$ and the product is $x_+ x_- = \frac{a+1}{a^2}$. Substituting these in:

$$g'(x_\pm) = a^2 \left[ 1 - 2\left(\frac{a+1}{a}\right) + 4\left(\frac{a+1}{a^2}\right) \right]$$

$$g'(x_\pm) = a^2 - 2a(a+1) + 4(a+1) = -a^2 + 2a + 4$$

For the 2-cycle to be stable, we need $\vert{}-a^2 + 2a + 4\vert{} < 1$:

$$-1 < -a^2 + 2a + 4 < 1$$

Solving the right inequality ($-a^2 + 2a + 3 < 0$) gives $a > 3$.

Solving the left inequality ($-a^2 + 2a + 5 > 0$) gives $a < 1 + \sqrt{6}$.

Therefore, the period-2 cycle is stable for the range **$3 < a < 1 + \sqrt{6}$**.

**7. The next bifurcation point** At $a = 1 + \sqrt{6} \approx 3.449$, the derivative $g'(x_\pm)$ crosses $-1$, meaning the period-2 cycle loses its stability. Exactly as we saw with the transition from period-1 to period-2, this triggers another **period-doubling bifurcation**. The system will spawn a stable orbit of period 4 (a 4-cycle). This is the onset of the famous Feigenbaum cascade, where continuous period-doubling eventually leads to deterministic chaos.


# 15. — A simple population dynamics

**1. Extremal point $x_m$** The map is defined as $f(x) = x \exp(r(1-x))$. To find the extrema, we take the first derivative using the product and chain rules and set it to zero:

$$f'(x) = \exp(r(1-x)) + x(-r)\exp(r(1-x)) = \exp(r(1-x))[1 - rx]$$

Setting $f'(x_m) = 0$ yields $1 - r x_m = 0$, giving the extremal point:

$$x_m = \frac{1}{r}$$

To determine if this is a maximum or minimum, evaluate the second derivative at $x_m$:

$$f''(x) = -r\exp(r(1-x))[1 - rx] - r\exp(r(1-x)) = -r\exp(r(1-x))[2 - rx]$$

Substitute $x_m = 1/r$:

$$f''(1/r) = -r\exp\left(r\left(1-\frac{1}{r}\right)\right)\left[2 - r\left(\frac{1}{r}\right)\right] = -r e^{r-1}$$

Since the exponential term $e^{r-1}$ is strictly positive for all real $r$:

- If **$r > 0$**, $f''(1/r) < 0$, meaning $x_m = 1/r$ is a **maximum**.
    
- If **$r < 0$**, $f''(1/r) > 0$, meaning $x_m = 1/r$ is a **minimum**.
    

**2. Fixed points $x^\star$** A fixed point satisfies the condition $x^\star = f(x^\star)$:

$$x^\star = x^\star \exp(r(1-x^\star))$$

This equation has two solutions:

- The trivial fixed point: **$x_0^\star = 0$**.
    
- For $x^\star \neq 0$, we can divide by $x^\star$ to get $1 = \exp(r(1-x^\star))$. Taking the natural logarithm gives $r(1-x^\star) = 0$. Assuming $r \neq 0$, this yields the non-trivial fixed point: **$x_1^\star = 1$**.
    

**3. Stability intervals** A fixed point is stable if the absolute value of the derivative is less than 1, $\vert{}f'(x^\star)\vert{} < 1$.

- **For $x_0^\star = 0$:** $f'(0) = \exp(r(1-0))[1 - 0] = e^r$. Stability requires $e^r < 1$, which gives the interval **$r < 0$**.
    
- **For $x_1^\star = 1$:** $f'(1) = \exp(r(1-1))[1 - r(1)] = 1 - r$. Stability requires $\vert{}1 - r\vert{} < 1 \implies -1 < 1 - r < 1 \implies -2 < -r < 0$. This gives the interval **$0 < r < 2$**.
    

(Note: Point 4 provides geometric context for the upcoming bifurcation calculation, moving directly to step 5).

**5. Period doubling bifurcation parameter $r_p$** A period-doubling (flip) bifurcation occurs when a stable fixed point loses stability by the derivative passing through $-1$ (i.e., $f'(x^\star) = -1$). Using the non-trivial fixed point $x_1^\star = 1$:

$$f'(1) = 1 - r_p = -1$$

**$r_p = 2$**

**6. Superstability condition for a two-cycle** A cycle is superstable if the derivative of the iterated map is zero at the points of the cycle. For a 2-cycle, the chain rule gives $(f^{(2)})'(x) = f'(f(x))f'(x) = 0$. This condition is met if the cycle contains the extremal point of the map, where the first derivative is zero ($f'(x_m)=0$). Therefore, one point of the superstable 2-cycle must be $x_m = 1/r$. The second point of the cycle is $f(x_m)$:

$$f(1/r) = \frac{1}{r}\exp\left(r\left(1-\frac{1}{r}\right)\right) = \frac{1}{r}e^{r-1}$$

For these two points to form a closed 2-cycle, applying the map to the second point must map exactly back to the first point: $f(f(x_m)) = x_m$.

$$f\left(\frac{1}{r}e^{r-1}\right) = \frac{1}{r}$$

Substitute $\frac{1}{r}e^{r-1}$ into the original function $f(x)$:

$$\left( \frac{1}{r}e^{r-1} \right) \exp\left( r\left(1 - \frac{1}{r}e^{r-1}\right) \right) = \frac{1}{r}$$

Multiply both sides by $r$:

$$e^{r-1} \exp\left( r - e^{r-1} \right) = 1$$

Combine the exponentials:

$$\exp(r - 1 + r - e^{r-1}) = 1$$

$$\exp(2r - 1 - e^{r-1}) = 1$$

Taking the natural logarithm of both sides ($\ln(1) = 0$) results directly in the required condition:

**$2r - 1 - e^{r-1} = 0$**

**7. Trivial solution $r_s$ and its corresponding cycle** We are looking for a simple root $r_s$ for the equation $2r - 1 - e^{r-1} = 0$. By inspection, testing **$r_s = 1$**:

$$2(1) - 1 - e^{1-1} = 2 - 1 - 1 = 0$$

Thus, $r_s = 1$ is a solution.

To determine what "2-cycle" this corresponds to, we evaluate the cycle points $x_m$ and $f(x_m)$ at $r=1$:

- $x_m = 1/1 = 1$
    
- $f(x_m) = \frac{1}{1}e^{1-1} = 1$
    

Both points are exactly the same. Therefore, this trivial solution does not correspond to a true oscillating 2-cycle, but rather to the **fixed point $x_1^\star = 1$**. At exactly $r=1$, the fixed point itself is superstable because $f'(1) = 1 - 1 = 0$.


# 16. — Unforced non-linear oscillators

**1. Fixed points and their stability**

Defining state variables $x_1 = \theta$ and $x_2 = \dot{\theta}$, the general oscillator equation $\ddot{\theta} + \theta + \epsilon h(\theta, \dot{\theta}) = 0$ becomes a system of first-order equations:

$$\dot{x}_1 = x_2$$

$$\dot{x}_2 = -x_1 - \epsilon h(x_1, x_2)$$

Fixed points require $\dot{x}_1 = 0$ and $\dot{x}_2 = 0$, meaning $x_2 = 0$ and $-x_1 - \epsilon h(x_1, 0) = 0$.

- **Van der Pol oscillator:** With $h(\theta, \dot{\theta}) = (\theta^2 - 1)\dot{\theta}$, the condition $h(x_1, 0) = 0$ implies the only fixed point is the origin $(0,0)$. The Jacobian evaluated at $(0,0)$ is:
    
    $$J = \begin{pmatrix} 0 & 1 \\ -1 & \epsilon \end{pmatrix}$$
    
    The eigenvalues are $\lambda = \frac{\epsilon \pm \sqrt{\epsilon^2 - 4}}{2}$. For $\epsilon > 0$, the real part is strictly positive, making the origin an **unstable** fixed point (an unstable spiral for $0 < \epsilon < 2$, and an unstable node for $\epsilon \ge 2$).
    
- **Duffing oscillator:** With $h(\theta, \dot{\theta}) = \theta^3$, the equilibrium condition is $-x_1 - \epsilon x_1^3 = 0 \implies x_1(1 + \epsilon x_1^2) = 0$. For $\epsilon > 0$, the only real root is $x_1 = 0$, yielding a fixed point at $(0,0)$. The Jacobian is:
    
    $$J = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$$
    
    The eigenvalues are $\lambda = \pm i$. The origin is a neutrally stable **center** in the linear approximation (and globally stable due to energy conservation).
    

**2. Phase portrait sketches**

- **Van der Pol:** Trajectories originating near the unstable origin spiral outwards, while trajectories starting far away spiral inwards. They all asymptotically converge onto a single, isolated closed trajectory called an attractive **limit cycle**.
    
- **Duffing:** The system is conservative and non-dissipative. The phase portrait consists of a continuous family of closed, concentric orbits surrounding the center at the origin, with no isolated limit cycle.
    

**3. Energy evolution and the limit cycle ($\epsilon \to 0$)**

The energy of the system is $E = \frac{1}{2}\theta^2 + \frac{1}{2}\dot{\theta}^2$. Differentiating with respect to time gives:

$$\dot{E} = \theta\dot{\theta} + \dot{\theta}\ddot{\theta} = \dot{\theta}(\theta + \ddot{\theta})$$

From the general equation, we substitute $\theta + \ddot{\theta} = -\epsilon h(\theta, \dot{\theta})$:

$$\dot{E} = -\epsilon \dot{\theta} h(\theta, \dot{\theta})$$

For the Van der Pol oscillator, this is $\dot{E} = -\epsilon \dot{\theta}^2 (\theta^2 - 1)$. In the weakly nonlinear limit ($\epsilon \to 0$), the trajectory is approximately harmonic over a single cycle: $\theta(t) \approx A \cos t$ and $\dot{\theta}(t) \approx -A \sin t$. We average the energy dissipation over one period $T = 2\pi$:

$$\langle \dot{E} \rangle = \frac{1}{2\pi} \int_0^{2\pi} -\epsilon (-A \sin t)^2 (A^2 \cos^2 t - 1) dt$$

$$\langle \dot{E} \rangle = -\frac{\epsilon A^2}{2\pi} \int_0^{2\pi} (A^2 \sin^2 t \cos^2 t - \sin^2 t) dt$$

Using the standard integrals $\int_0^{2\pi} \sin^2 t dt = \pi$ and $\int_0^{2\pi} \sin^2 t \cos^2 t dt = \pi/4$:

$$\langle \dot{E} \rangle = -\frac{\epsilon A^2}{2\pi} \left( A^2 \frac{\pi}{4} - \pi \right) = -\frac{\epsilon A^2}{8} (A^2 - 4)$$

A stable limit cycle corresponds to a state of zero net energy change per cycle ($\langle \dot{E} \rangle = 0$). For a non-trivial amplitude ($A > 0$), this requires $A^2 - 4 = 0$, giving a limit cycle radius of **$A = 2$**.

**4. Liénard transformation for Van der Pol**

Starting with $\ddot{\theta} + \epsilon(\theta^2 - 1)\dot{\theta} + \theta = 0$, we integrate the damping term to define $F(x) = \frac{x^3}{3} - x$. This allows us to rewrite the equation as a total derivative:

$$\frac{d}{dt}\left[ \dot{\theta} + \epsilon F(\theta) \right] + \theta = 0$$

Let $x = \theta$. To match the requested form $\dot{x} = \epsilon(y - F(x))$, we define the new variable $y$ such that:

$$y = \frac{\dot{x}}{\epsilon} + F(x)$$

Differentiating $y$ with respect to time gives:

$$\dot{y} = \frac{\ddot{x}}{\epsilon} + F'(x)\dot{x} = \frac{1}{\epsilon}\left[ \ddot{\theta} + \epsilon(\theta^2 - 1)\dot{\theta} \right]$$

Substituting the original differential equation ($\ddot{\theta} + \epsilon(\theta^2 - 1)\dot{\theta} = -\theta = -x$) yields:

$$\dot{y} = -\frac{x}{\epsilon}$$

**5. Limit cycle in the $\epsilon \to \infty$ limit (Relaxation oscillations)**

In the $(x,y)$ plane, as $\epsilon \to \infty$, the velocity $\dot{x} = \epsilon(y - F(x))$ becomes exceptionally fast ($O(\epsilon)$) unless the system is on the cubic nullcline $y = F(x)$, while $\dot{y} = -x/\epsilon$ dictates a very slow vertical drift ($O(1/\epsilon)$). The limit cycle behavior separates into two distinct phases:

- **Slow flow:** The trajectory is pinned to the stable outer branches of the cubic curve $y = \frac{x^3}{3} - x$, slowly drifting towards the local extrema located at $x = 1$ (minimum) and $x = -1$ (maximum).
    
- **Fast jumps:** Upon reaching an extremum ($x=\pm 1$), the branch loses stability. The system instantaneously jumps horizontally ($\dot{y} \approx 0$) to the opposite stable branch of the cubic curve.
    
    The limit cycle is a closed loop connecting $(2, 2/3)$ down to $(1, -2/3)$, jumping left to $(-2, -2/3)$, crawling up to $(-1, 2/3)$, and jumping right back to $(2, 2/3)$.
    

**6. Amplitude-dependent frequency of the Duffing oscillator**

Introducing a normalized time $\tau = \omega t$, the derivatives become $\frac{d}{dt} = \omega \frac{d}{d\tau}$. The equation transforms to:

$$\omega^2 \theta'' + \theta + \epsilon \theta^3 = 0$$

We perform a Poincaré-Lindstedt perturbation expansion for both the solution and frequency:

$$\theta(\tau) = \theta_0(\tau) + \epsilon \theta_1(\tau) + \mathcal{O}(\epsilon^2)$$

$$\omega = 1 + \epsilon \omega_1 + \mathcal{O}(\epsilon^2) \implies \omega^2 = 1 + 2\epsilon\omega_1 + \mathcal{O}(\epsilon^2)$$

Substituting these into the equation and collecting terms by powers of $\epsilon$:

- $\mathcal{O}(1)$: $\theta_0'' + \theta_0 = 0$. We select the fundamental solution $\theta_0(\tau) = A \cos \tau$.
    
- $\mathcal{O}(\epsilon)$: $\theta_1'' + \theta_1 = -2\omega_1 \theta_0'' - \theta_0^3$.
    
    Substituting $\theta_0$ into the $\mathcal{O}(\epsilon)$ right-hand side gives:
    
    $$\theta_1'' + \theta_1 = 2\omega_1 A \cos \tau - A^3 \cos^3 \tau$$
    
    Using the identity $\cos^3 \tau = \frac{3}{4}\cos \tau + \frac{1}{4}\cos 3\tau$:
    
    $$\theta_1'' + \theta_1 = \left( 2\omega_1 A - \frac{3}{4}A^3 \right) \cos \tau - \frac{A^3}{4} \cos 3\tau$$
    
    To prevent secular growth (which would invalidate the periodic assumption), the resonant forcing term proportional to $\cos \tau$ must be eliminated:
    
    $$2\omega_1 A - \frac{3}{4}A^3 = 0 \implies \omega_1 = \frac{3}{8}A^2$$
    
    Thus, the frequency depends quadratically on the amplitude: $\omega(A) \approx 1 + \frac{3}{8}\epsilon A^2$.
    

**7. Two-timing expansion for the Van der Pol oscillator**

We introduce two independent time scales: a fast time $\tau = t$ and a slow time $T = \epsilon t$. The time derivative becomes $\frac{d}{dt} = \partial_\tau + \epsilon \partial_T$, and $\frac{d^2}{dt^2} \approx \partial_{\tau\tau} + 2\epsilon \partial_{\tau T}$. Expanding the solution as $\theta = \theta_0(\tau, T) + \epsilon \theta_1(\tau, T)$:

- $\mathcal{O}(1)$: $\partial_{\tau\tau} \theta_0 + \theta_0 = 0$.
    
    The solution is a harmonic oscillator where the amplitude $A(T)$ and phase $\phi(T)$ evolve slowly: $\theta_0(\tau, T) = A(T) \cos(\tau + \phi(T))$.
    
- $\mathcal{O}(\epsilon)$: $\partial_{\tau\tau} \theta_1 + \theta_1 = -2\partial_{\tau T}\theta_0 - (\theta_0^2 - 1)\partial_\tau \theta_0$.
    
    Calculating the components for the right-hand side:
    
    $$-2\partial_{\tau T}\theta_0 = 2A' \sin(\tau+\phi) + 2A\phi' \cos(\tau+\phi)$$
    
    $$-(\theta_0^2 - 1)\partial_\tau \theta_0 = A\sin(\tau+\phi) \left( A^2 \cos^2(\tau+\phi) - 1 \right)$$
    
    Using $\cos^2 x \sin x = \frac{1}{4}\sin x + \frac{1}{4}\sin 3x$, the nonlinear term yields $\left(\frac{A^3}{4} - A\right)\sin(\tau+\phi)$ plus higher harmonics.
    
    Combining all terms driving $\theta_1$:
    
    $$\partial_{\tau\tau} \theta_1 + \theta_1 = \left[ 2A' + \frac{A^3}{4} - A \right]\sin(\tau+\phi) + \left[ 2A\phi' \right]\cos(\tau+\phi) + \text{H.O.T.}$$
    
    To eliminate secular terms, the coefficients of the fundamental harmonics must vanish:
    

1. $2A\phi' = 0 \implies \phi(T) = \text{constant}$.
    
2. $2A' = A - \frac{A^3}{4} \implies A'(T) = \frac{1}{2}A\left(1 - \frac{A^2}{4}\right)$.
    
    This slow-flow differential equation reveals that any initial amplitude $A > 0$ will asymptotically converge to the stable equilibrium point $A = 2$, proving the existence and exact radius of the limit cycle.



**Secular Growth and the Resonant Forcing Term**

In perturbation theory (like the Poincaré-Lindstedt method used for the Duffing oscillator), you are solving a sequence of linear differential equations. At order $\mathcal{O}(\epsilon)$, the equation takes the form of a driven harmonic oscillator:

$$\theta_1'' + \theta_1 = F(\tau)$$

The left side, $\theta_1'' + \theta_1$, describes a system with a natural angular frequency of $\omega_0 = 1$. The right side acts as an external driving force. If this driving force contains terms that oscillate at the system's exact natural frequency (i.e., terms proportional to $\cos \tau$ or $\sin \tau$), the system is driven at **resonance**.

From the theory of linear differential equations, if you drive a harmonic oscillator at resonance, the amplitude of the particular solution grows linearly with time (e.g., $\tau \sin \tau$). These linearly growing terms are called **secular terms**.

Because the Duffing oscillator is a conservative system that oscillates continuously within a bounded region of phase space, its true solution must be strictly periodic. An amplitude that grows to infinity as $\tau \to \infty$ violates this physical reality. Therefore, to ensure the validity of the periodic ansatz, we must mathematically enforce that resonance does not occur in the $\mathcal{O}(\epsilon)$ equation. We do this by setting the net coefficient of any $\cos \tau$ or $\sin \tau$ terms on the right-hand side strictly to zero.

**
![[Pasted image 20261005233017.png]]


# 17. — Forced non-linear pendulum
### 1. Linear Response and Amplitude

The governing equation is $\ddot{\theta} = -\omega_0^2\sin\theta - \mu\dot{\theta} + f\cos(\omega t)$. At the lowest order of $\theta$, we use the small-angle approximation $\sin\theta \approx \theta$. This yields a linear, driven, damped harmonic oscillator equation:

$$\ddot{\theta} + \mu\dot{\theta} + \omega_0^2\theta = f\cos(\omega t)$$

To find the steady-state response, it is mathematically convenient to use complex exponentials. We represent the forcing as the real part of $f e^{i\omega t}$ and assume a complex solution $\underline{\theta}(t) = \underline{\theta_0} e^{i\omega t}$, where the complex amplitude is $\underline{\theta_0} = \theta_0 e^{i\phi}$.

Taking the derivatives:

- $\dot{\underline{\theta}} = i\omega \underline{\theta_0} e^{i\omega t}$
    
- $\ddot{\underline{\theta}} = -\omega^2 \underline{\theta_0} e^{i\omega t}$
    

Substituting these into the linearized equation and canceling the $e^{i\omega t}$ terms:

$$(-\omega^2 + i\mu\omega + \omega_0^2)\underline{\theta_0} = f$$

$$\underline{\theta_0} = \frac{f}{(\omega_0^2 - \omega^2) + i\mu\omega}$$

The physical amplitude $\theta_0$ is the magnitude of this complex amplitude:

$$\theta_0 = \vert{}\underline{\theta_0}\vert{} = \frac{f}{\sqrt{(\omega_0^2 - \omega^2)^2 + \mu^2\omega^2}}$$

This independently confirms the relation provided in the problem statement.

**Plotting $\theta_0(\omega)$:**

- **For $\mu = 0$ (undamped):** The denominator becomes simply $\vert{}\omega_0^2 - \omega^2\vert{}$. The response curve starts at $f/\omega_0^2$ for $\omega=0$, asymptotically diverges to infinity at the resonance frequency ($\omega \to \omega_0$), and approaches zero as $\omega \to \infty$.
    
- **For $\mu \neq 0$ (small damping):** The curve is bounded. It starts near $f/\omega_0^2$, reaches a finite maximum peak of approximately $f/(\mu\omega_0)$ near $\omega \approx \omega_0$ (specifically at the slightly shifted resonance frequency $\omega_{res} = \sqrt{\omega_0^2 - \mu^2/2}$), and decays to zero at high frequencies.
    

### 2. Next Order Expansion ($\mu = 0$)

Setting $\mu = 0$, we expand $\sin\theta$ to the third order (lowest non-linear term): $\sin\theta \approx \theta - \frac{\theta^3}{6}$. The equation of motion becomes:

$$\ddot{\theta} + \omega_0^2\left(\theta - \frac{\theta^3}{6}\right) = f\cos(\omega t)$$

Assuming a solution of the form $\theta(t) = \theta_0\cos(\omega t)$, we evaluate the terms:

- $\ddot{\theta} = -\omega^2\theta_0\cos(\omega t)$
    
- $\theta^3 = \theta_0^3\cos^3(\omega t)$
    

Using the trigonometric identity $\cos^3(\omega t) = \frac{3}{4}\cos(\omega t) + \frac{1}{4}\cos(3\omega t)$, we substitute the terms back into the differential equation:

$$-\omega^2\theta_0\cos(\omega t) + \omega_0^2\left[ \theta_0\cos(\omega t) - \frac{\theta_0^3}{6}\left( \frac{3}{4}\cos(\omega t) + \frac{1}{4}\cos(3\omega t) \right) \right] = f\cos(\omega t)$$

We neglect the non-resonant higher harmonic term proportional to $\cos(3\omega t)$, and balance the coefficients of the primary resonant term $\cos(\omega t)$:

$$-\omega^2\theta_0 + \omega_0^2\theta_0 - \frac{\omega_0^2}{8}\theta_0^3 = f$$

Rearranging to express the relationship between the driving frequency $\omega$ and the amplitude $\theta_0$:

$$(\omega_0^2 - \omega^2)\theta_0 - \frac{\omega_0^2}{8}\theta_0^3 = f$$

$$\omega^2 = \omega_0^2 - \frac{\omega_0^2}{8}\theta_0^2 - \frac{f}{\theta_0}$$

### 3. Plotting $\theta_0(\omega)$ for $\mu = 0$

To visualize this relationship, it is helpful to look at the unforced system first ($f=0$). The "backbone" curve is given by $\omega = \omega_0 \sqrt{1 - \theta_0^2/8}$. Because the non-linear term $-\theta^3/6$ acts as a "softening" spring (the restoring force is weaker than a linear spring at large amplitudes), the resonance frequency decreases as amplitude increases, causing the backbone curve to bend to the left.

With the forcing term $f > 0$ included, the response curve splits into two disconnected branches:

1. One branch emerges from $\omega = 0$, rises, and closely follows the left side of the softening backbone curve, going to infinite amplitude at lower frequencies.
    
2. The second branch comes from $\omega \to \infty$ at zero amplitude, rises as it approaches the resonance region from the right, and goes to infinite amplitude while bounding the right side of the backbone.
    

### 4. Deducing the Diagram for $\mu \neq 0$

When viscous damping is introduced, the two divergent branches connect to form a continuous, rounded curve with a finite peak.

However, because the entire resonance structure is "bent" to the left by the softening non-linearity, the single connected peak leans significantly to the left (towards lower frequencies).

This leftward lean creates a region of bistability. For a specific range of driving frequencies just below $\omega_0$, a vertical line drawn at a constant $\omega$ will intersect the response curve at three distinct amplitudes:

- A large-amplitude stable solution.
    
- A small-amplitude stable solution.
    
- A mid-amplitude unstable solution that connects the two.
    

Physically, this leads to the "jump phenomenon" or hysteresis. If you slowly sweep the driving frequency downward, the amplitude smoothly increases along the upper branch until it drops discontinuously at the edge of the fold. If you sweep upward from low frequencies, the amplitude stays on the lower branch until it jumps discontinuously up to the larger branch.

![[Pasted image 20260923215405.png]]

# 18. — Parametric pendulum 
### 1. Predicting the shape of the amplitude equation

We are looking for a solution of the form $\theta(t) = A(t) e^{i\omega t/2} + A^*(t) e^{-i\omega t/2}$. This ansatz assumes the pendulum oscillates at half the driving frequency (the primary parametric resonance), while $A(t)$ is a "slowly varying" complex amplitude.

To predict the equation for $\dot{A}$, we look at the symmetries of the unforced system ($\lambda = 0, f = 0$). The basic pendulum equation $\ddot{\theta} + \omega_0^2\sin\theta = 0$ is invariant if we flip the angle ($\theta \to -\theta$). Because the driving term $(1 + f\sin\omega t)$ multiplies $\sin\theta$, the full equation $\ddot{\theta} + 2\lambda\dot{\theta} + \omega_0^2\sin\theta(1 + f\sin\omega t) = 0$ is also invariant under $\theta \to -\theta$.

Since $\theta$ is proportional to $A$, the transformation $\theta \to -\theta$ means $A \to -A$. Therefore, the differential equation governing $A$ must be odd in $A$ (meaning it can only contain terms like $A$, $A^*$, $\vert{}A\vert{}^2 A$, etc., but no even powers like $A^2$ or $\vert{}A\vert{}^2$).

### 2. Deriving the amplitude equation step-by-step

We are given that the detuning $\delta$ (where $\omega = 2\omega_0 + \delta$), the damping $\lambda$, and the forcing amplitude $f$ are all small and of the same order of magnitude. Let's call this small order $\epsilon$. Because these perturbations are small, the amplitude $A(t)$ changes very slowly over time. Mathematically, this means:

- $\dot{A}$ is of order $\epsilon$.
    
- $\ddot{A}$ is of order $\epsilon^2$ (and can therefore be neglected).
    
- Any product of two small parameters (like $\lambda \dot{A}$ or $\delta \dot{A}$) is order $\epsilon^2$ and can be neglected.
    

**Step A: Calculate the derivatives of $\theta$**

Substitute $\theta = A e^{i\omega t/2} + A^* e^{-i\omega t/2}$ into the time derivatives:

$$\dot{\theta} = \left(\dot{A} + i\frac{\omega}{2}A\right) e^{i\omega t/2} + \text{c.c.}$$

$$\ddot{\theta} = \left(\ddot{A} + 2i\frac{\omega}{2}\dot{A} - \frac{\omega^2}{4}A\right) e^{i\omega t/2} + \text{c.c.}$$

Now, apply our small-parameter approximations to $\ddot{\theta}$. We drop $\ddot{A}$. We substitute $\omega = 2\omega_0 + \delta$ into the $\dot{A}$ term (keeping only the lowest order $\omega \approx 2\omega_0$) and into the $A$ term (where we must expand it):

$$\frac{\omega^2}{4} = \frac{(2\omega_0 + \delta)^2}{4} = \omega_0^2 + \omega_0\delta + \frac{\delta^2}{4} \approx \omega_0^2 + \omega_0\delta$$

So, the acceleration simplifies to:

$$\ddot{\theta} \approx \left[ i(2\omega_0)\dot{A} - (\omega_0^2 + \omega_0\delta)A \right] e^{i\omega t/2} + \text{c.c.}$$

**Step B: Expand the non-linear pendulum term**

Expand the sine function to the lowest non-linear order: $\sin\theta \approx \theta - \frac{1}{6}\theta^3$.

We need to calculate $\theta^3 = (A e^{i\omega t/2} + A^* e^{-i\omega t/2})^3$. Expanding this binomial gives four terms:

$$\theta^3 = A^3 e^{i 3\omega t/2} + 3A^2 A^* e^{i\omega t/2} + 3A(A^*)^2 e^{-i\omega t/2} + (A^*)^3 e^{-i 3\omega t/2}$$

We only care about the "secular" or "resonant" terms that oscillate at the fundamental frequency $e^{i\omega t/2}$. The higher harmonics ($e^{i 3\omega t/2}$) oscillate too fast to affect the slow drift of the amplitude and average out to zero.

The resonant part of $\sin\theta$ is therefore:

$$\left( A - \frac{3}{6}\vert{}A\vert{}^2 A \right) e^{i\omega t/2} = \left( A - \frac{1}{2}\vert{}A\vert{}^2 A \right) e^{i\omega t/2}$$

**Step C: Evaluate the parametric forcing term** The forcing term is $\omega_0^2 \theta f \sin(\omega t)$. Let's rewrite the sine using Euler's formula: $f\sin(\omega t) = \frac{f}{2i}(e^{i\omega t} - e^{-i\omega t})$. Multiply this by $\theta$:

$$\theta f \sin(\omega t) = (A e^{i\omega t/2} + A^* e^{-i\omega t/2}) \frac{f}{2i}(e^{i\omega t} - e^{-i\omega t})$$

Multiply it out term by term:

$$= \frac{f}{2i} \left[ A e^{i 3\omega t/2} - A e^{-i\omega t/2} + A^* e^{i\omega t/2} - A^* e^{-i 3\omega t/2} \right]$$

Again, we extract only the term proportional to $e^{i\omega t/2}$:

$$\text{Resonant forcing} = \frac{f}{2i} A^* e^{i\omega t/2} = -i\frac{f}{2} A^* e^{i\omega t/2}$$

_(Note: We do not multiply the forcing term by the cubic $-\theta^3/6$ part of the sine expansion, because $f \cdot \theta^3$ would be the product of two small parameters, which is negligible)._

**Step D: Bring it all together**

Substitute all our resonant terms back into the original differential equation, dropping the $e^{i\omega t/2}$ common factor.

- From $\ddot{\theta}$: $2i\omega_0\dot{A} - \omega_0^2 A - \omega_0\delta A$
    
- From $2\lambda\dot{\theta}$: $2\lambda(i\omega_0 A)$ _(we drop the $2\lambda\dot{A}$ term as it is order $\epsilon^2$)_
    
- From $\omega_0^2\sin\theta$: $\omega_0^2 A - \frac{\omega_0^2}{2}\vert{}A\vert{}^2 A$
    
- From forcing: $-i\frac{f\omega_0^2}{2} A^*$
    

Summing them to zero:

$$[2i\omega_0\dot{A} - \omega_0^2 A - \omega_0\delta A] + [2i\lambda\omega_0 A] + \left[\omega_0^2 A - \frac{\omega_0^2}{2}\vert{}A\vert{}^2 A\right] - \left[i\frac{f\omega_0^2}{2} A^*\right] = 0$$

Notice that the $-\omega_0^2 A$ and $+\omega_0^2 A$ terms perfectly cancel each other out:

$$2i\omega_0\dot{A} - \omega_0\delta A + 2i\lambda\omega_0 A - \frac{\omega_0^2}{2}\vert{}A\vert{}^2 A - i\frac{f\omega_0^2}{2} A^* = 0$$

Now, we isolate $\dot{A}$. Move everything else to the right side:

$$2i\omega_0\dot{A} = \omega_0\delta A - 2i\lambda\omega_0 A + \frac{\omega_0^2}{2}\vert{}A\vert{}^2 A + i\frac{f\omega_0^2}{2} A^*$$

Divide every term by $2i\omega_0$. Remember that $1/i = -i$:

$$\dot{A} = -i\frac{\delta}{2} A - \lambda A - i\frac{\omega_0}{4} \vert{}A\vert{}^2 A + \frac{f\omega_0}{4} A^*$$

Finally, factor out the common $A$ on the linear terms to get the clean, canonical amplitude equation:

$$\dot{A} = -\left( \lambda + i\frac{\delta}{2} \right)A + \frac{f\omega_0}{4} A^* - i\frac{\omega_0}{4} \vert{}A\vert{}^2 A$$
### 3. Linear stability of the null solution

To study the linear stability of the trivial solution ($A = 0$), we drop the non-linear term $-i\frac{\omega_0}{4}\vert{}A\vert{}^2A$ from our derived amplitude equation. This leaves us with the linearized system:

$$\dot{A} = -\left( \lambda + i\frac{\delta}{2} \right)A + \frac{f\omega_0}{4} A^*$$

We decompose the complex amplitude into its real and imaginary parts by substituting $A = X + iY$ and $A^* = X - iY$:

$$\dot{X} + i\dot{Y} = -\left( \lambda + i\frac{\delta}{2} \right)(X + iY) + \frac{f\omega_0}{4} (X - iY)$$

Expanding the right-hand side:

$$\dot{X} + i\dot{Y} = -\lambda X - i\lambda Y - i\frac{\delta}{2}X + \frac{\delta}{2}Y + \frac{f\omega_0}{4}X - i\frac{f\omega_0}{4}Y$$

By equating the real and imaginary parts, we construct a system of two coupled linear differential equations:

- **Real part:** $\dot{X} = \left(-\lambda + \frac{f\omega_0}{4}\right)X + \frac{\delta}{2}Y$
    
- **Imaginary part:** $\dot{Y} = -\frac{\delta}{2}X - \left(\lambda + \frac{f\omega_0}{4}\right)Y$
    

This can be written in matrix form as $\begin{pmatrix} \dot{X} \\ \dot{Y} \end{pmatrix} = J \begin{pmatrix} X \\ Y \end{pmatrix}$, where the Jacobian $J$ is:

$$J = \begin{pmatrix} -\lambda + \frac{f\omega_0}{4} & \frac{\delta}{2} \\ -\frac{\delta}{2} & -\lambda - \frac{f\omega_0}{4} \end{pmatrix}$$

The null solution $(X,Y) = (0,0)$ is stable if both eigenvalues of $J$ have negative real parts. This requires the trace of $J$ to be negative and the determinant of $J$ to be positive.

- **Trace:** $\text{Tr}(J) = \left(-\lambda + \frac{f\omega_0}{4}\right) + \left(-\lambda - \frac{f\omega_0}{4}\right) = -2\lambda$. Since the damping $\lambda$ is strictly positive, the trace is always negative.
    
- **Determinant:** $\text{Det}(J) = \left(-\lambda + \frac{f\omega_0}{4}\right)\left(-\lambda - \frac{f\omega_0}{4}\right) - \left(\frac{\delta}{2}\right)\left(-\frac{\delta}{2}\right) = \lambda^2 - \left(\frac{f\omega_0}{4}\right)^2 + \frac{\delta^2}{4}$.
    

For the system to become unstable (which physically corresponds to parametric resonance occurring), the determinant must be negative ($\text{Det}(J) < 0$):

$$\left(\frac{f\omega_0}{4}\right)^2 > \lambda^2 + \frac{\delta^2}{4}$$

In the $(\delta, f)$ parameter plane, the stability boundary is given by the hyperbola $f = \frac{4}{\omega_0}\sqrt{\lambda^2 + \frac{\delta^2}{4}}$. The null solution is stable below this curve and unstable inside the "V-shaped" region above it, which is the principal Arnold tongue for parametric resonance.

### 4. Stationary solutions of the amplitude equation

To find the steady-state amplitudes of the non-linear system, we set $\dot{A} = 0$ in the full equation:

$$0 = -\left( \lambda + i\frac{\delta}{2} \right)A + \frac{f\omega_0}{4} A^* - i\frac{\omega_0}{4} \vert{}A\vert{}^2 A$$

We write $A$ in polar form as $A = R e^{i\phi}$ (which implies $A^* = R e^{-i\phi}$ and $\vert{}A\vert{}^2 = R^2$):

$$0 = -\lambda R e^{i\phi} - i\frac{\delta}{2} R e^{i\phi} + \frac{f\omega_0}{4} R e^{-i\phi} - i\frac{\omega_0}{4} R^3 e^{i\phi}$$

Assuming a non-trivial solution ($R \neq 0$), we divide the entire equation by $R e^{i\phi}$:

$$0 = -\lambda - i\frac{\delta}{2} + \frac{f\omega_0}{4} e^{-2i\phi} - i\frac{\omega_0}{4} R^2$$

Using Euler's formula ($e^{-2i\phi} = \cos(2\phi) - i\sin(2\phi)$), we separate the real and imaginary parts:

1. **Real:** $0 = -\lambda + \frac{f\omega_0}{4} \cos(2\phi) \implies \lambda = \frac{f\omega_0}{4} \cos(2\phi)$
    
2. **Imaginary:** $0 = -\frac{\delta}{2} - \frac{f\omega_0}{4} \sin(2\phi) - \frac{\omega_0}{4} R^2 \implies \frac{\delta}{2} + \frac{\omega_0}{4} R^2 = -\frac{f\omega_0}{4} \sin(2\phi)$
    

To eliminate the phase $\phi$, we square both equations and add them together, utilizing the identity $\cos^2(2\phi) + \sin^2(2\phi) = 1$:

$$\lambda^2 + \left(\frac{\delta}{2} + \frac{\omega_0}{4} R^2\right)^2 = \left(\frac{f\omega_0}{4}\right)^2$$

Now, we solve for the stationary amplitude $R$:

$$\left(\frac{\delta}{2} + \frac{\omega_0}{4} R^2\right)^2 = \left(\frac{f\omega_0}{4}\right)^2 - \lambda^2$$

$$\frac{\delta}{2} + \frac{\omega_0}{4} R^2 = \pm \sqrt{\left(\frac{f\omega_0}{4}\right)^2 - \lambda^2}$$

$$R^2 = -\frac{2\delta}{\omega_0} \pm \frac{4}{\omega_0} \sqrt{\left(\frac{f\omega_0}{4}\right)^2 - \lambda^2}$$

This demonstrates that for a given detuning $\delta$ and forcing $f$, there can be either zero, one, or two non-trivial steady-state amplitudes depending on the discriminant.

### Part 2: High Frequency Limit ($\omega \gg \omega_0$)

#### 1. Derivation of the slow equation

Neglecting damping ($\lambda = 0$), the full equation is:

$$\ddot{\theta} + \omega_0^2\sin\theta(1 + f\sin\omega t) = 0$$

We decompose the motion into a slowly evolving angle $\theta_l$ and a rapid fluctuation $\theta_r$ such that $\theta = \theta_l + \theta_r$. Because $\theta_r$ is small, we Taylor expand the sine term: $\sin(\theta_l + \theta_r) \approx \sin\theta_l + \theta_r\cos\theta_l$. Substituting this into the equation of motion:

$$\ddot{\theta}_l + \ddot{\theta}_r + \omega_0^2(\sin\theta_l + \theta_r\cos\theta_l)(1 + f\sin\omega t) = 0$$

$$\ddot{\theta}_l + \ddot{\theta}_r + \omega_0^2\sin\theta_l + \omega_0^2 f \sin\theta_l \sin\omega t + \omega_0^2 \theta_r \cos\theta_l + \omega_0^2 f \theta_r \cos\theta_l \sin\omega t = 0$$

We separate the dynamics by time scale. The dominant terms driving the rapid, small-amplitude fluctuations ($\theta_r$) are those oscillating explicitly at the high frequency $\omega$:

$$\ddot{\theta}_r \approx -\omega_0^2 f \sin\theta_l \sin\omega t$$

To find $\theta_r$, we integrate twice with respect to the fast time $t$, treating the slow variables ($\theta_l$) as constants:

$$\dot{\theta}_r = \frac{\omega_0^2 f}{\omega} \sin\theta_l \cos\omega t$$

$$\theta_r = \frac{\omega_0^2 f}{\omega^2} \sin\theta_l \sin\omega t$$

Next, we obtain the equation for the slow evolution by taking the time average of the full equation over one fast period ($T = 2\pi/\omega$). Let $\langle \cdot \rangle$ denote this time average.

- $\langle \ddot{\theta}_r \rangle = 0$
    
- $\langle \sin\omega t \rangle = 0$
    
- $\langle \theta_r \cos\theta_l \rangle = \cos\theta_l \langle \theta_r \rangle = 0$
    

The only non-trivial fast term that does not average to zero is the product of $f\sin\omega t$ and $\theta_r$:

$$\langle \omega_0^2 f \theta_r \cos\theta_l \sin\omega t \rangle = \omega_0^2 f \cos\theta_l \langle \theta_r \sin\omega t \rangle$$

Substitute our derived expression for $\theta_r$:

$$= \omega_0^2 f \cos\theta_l \left\langle \left( \frac{\omega_0^2 f}{\omega^2} \sin\theta_l \sin\omega t \right) \sin\omega t \right\rangle$$

$$= \frac{\omega_0^4 f^2}{\omega^2} \sin\theta_l \cos\theta_l \langle \sin^2\omega t \rangle$$

Since the average of $\sin^2(\omega t)$ over one period is $1/2$, this term becomes $\frac{\omega_0^4 f^2}{2\omega^2} \sin\theta_l \cos\theta_l$.

Bringing the surviving averaged terms back together:

$$\ddot{\theta}_l + \omega_0^2\sin\theta_l + \frac{\omega_0^4 f^2}{2\omega^2} \sin\theta_l \cos\theta_l = 0$$

Factoring out $\omega_0^2 \sin\theta_l$ yields the exact requested slow equation:

$$\ddot{\theta}_l + \omega_0^2\sin\theta_l \left( 1 + \frac{f^2\omega_0^2}{2\omega^2}\cos\theta_l \right) = 0$$

#### 2. Stabilization of the inverted position

The inverted equilibrium position corresponds to $\theta_l = \pi$. We analyze its stability by introducing a small angular perturbation $\epsilon$, such that $\theta_l = \pi + \epsilon$.

Using the lowest-order Taylor expansions for the perturbation ($\sin(\pi + \epsilon) \approx -\epsilon$ and $\cos(\pi + \epsilon) \approx -1$), we linearize the slow equation around the inverted state:

$$\ddot{\epsilon} + \omega_0^2(-\epsilon) \left( 1 + \frac{f^2\omega_0^2}{2\omega^2}(-1) \right) = 0$$

$$\ddot{\epsilon} - \omega_0^2 \left( 1 - \frac{f^2\omega_0^2}{2\omega^2} \right) \epsilon = 0$$

For this position to be stable, the system must act as a simple harmonic oscillator ($\ddot{\epsilon} + \Omega_{eff}^2 \epsilon = 0$), meaning the coefficient of $\epsilon$ must provide a positive restoring force. This requires:

$$-\left( 1 - \frac{f^2\omega_0^2}{2\omega^2} \right) > 0$$

$$1 - \frac{f^2\omega_0^2}{2\omega^2} < 0$$

$$\frac{f^2\omega_0^2}{2\omega^2} > 1$$

Therefore, the inverted pendulum is dynamically stabilized when the forcing amplitude $f$ exceeds the critical threshold $f_c = \frac{\sqrt{2}\omega}{\omega_0}$. This is the classic mechanism of the Kapitza pendulum.

#### 3. Nature of the bifurcation

Let $K = \frac{f^2\omega_0^2}{2\omega^2}$ act as our dimensionless control parameter. We expand the restoring force to the next non-linear order ($\epsilon^3$) near the critical point $K \approx 1$.

Using $\sin(\pi+\epsilon) \approx -\epsilon + \epsilon^3/6$ and $\cos(\pi+\epsilon) \approx -1 + \epsilon^2/2$:

$$\ddot{\epsilon} + \omega_0^2\left(-\epsilon + \frac{\epsilon^3}{6}\right) \left[1 + K\left(-1 + \frac{\epsilon^2}{2}\right)\right] = 0$$

Retaining terms up to $\mathcal{O}(\epsilon^3)$:

$$\ddot{\epsilon} - \omega_0^2(1 - K)\epsilon + \omega_0^2\left(\frac{1 - 4K}{6}\right)\epsilon^3 = 0$$

This represents a system moving in an effective potential $V(\epsilon) = -\frac{\omega_0^2}{2}(1-K)\epsilon^2 + \frac{\omega_0^2}{4}\left(\frac{1-4K}{6}\right)\epsilon^4$.

- For $K < 1$ (weak forcing), the quadratic term is negative, making $\epsilon = 0$ an unstable local maximum.
    
- For $K > 1$ (strong forcing), the quadratic term becomes positive, making $\epsilon = 0$ a stable local minimum. However, because $\left(\frac{1-4K}{6}\right)$ is strictly negative near $K=1$, the $\epsilon^4$ coefficient is negative.
    

This indicates that as the origin becomes stable, two new _unstable_ equilibrium branches ($\epsilon \approx \pm\sqrt{2(K-1)}$) are symmetrically created and move outward. Because the stabilization of the central fixed point coincides with the shedding of unstable equilibria, this represents a **subcritical pitchfork bifurcation**.
![[Pasted image 20261005233736.png]]

# 19. — Two-dimensional jets

The problem explores the fluid dynamics of a two-dimensional laminar jet discharging into a fluid of the same density. Because the jet is thin and spreads slowly in the longitudinal direction ($z$), the transverse gradients ($\partial_x$) are much larger than the longitudinal gradients ($\partial_z$). This is the fundamental assumption of boundary layer theory. Furthermore, because the ambient fluid is stationary, the pressure is uniform everywhere, eliminating the pressure gradient term ($\partial_x P = 0$).

Here is the step-by-step derivation for the first two questions.

### 1. Conservation of Momentum Flux

We start with the longitudinal momentum balance and the continuity (incompressibility) equation provided in the problem:

- Continuity: $\partial_x u_x + \partial_z u_z = 0$
    
- Momentum: $u_x \partial_x u_z + u_z \partial_z u_z = \partial_x (\nu \partial_x u_z)$
    

First, we rewrite the left side of the momentum equation into a conservative form. By applying the product rule for derivatives, we can express the advection terms as:

$$\partial_x (u_x u_z) + \partial_z (u_z^2) = (u_z \partial_x u_x + u_x \partial_x u_z) + 2 u_z \partial_z u_z $$

$$\partial_x (u_x u_z) + \partial_z (u_z^2) = u_x \partial_x u_z + u_z \partial_z u_z + u_z (\partial_x u_x + \partial_z u_z)$$

Because the fluid is incompressible, the term inside the parentheses on the far right is exactly zero ($\partial_x u_x + \partial_z u_z = 0$). Thus, the momentum equation can be perfectly rewritten as:

$$\partial_x (u_x u_z) + \partial_z (u_z^2) = \nu \partial_{xx} u_z$$

Next, we integrate this entire equation over the transverse coordinate $x$ from $-\infty$ to $+\infty$:

$$\int_{-\infty}^{\infty} \partial_x (u_x u_z) dx + \int_{-\infty}^{\infty} \partial_z (u_z^2) dx = \int_{-\infty}^{\infty} \nu \partial_{xx} u_z dx$$

For the first term and the term on the right-hand side, the integrals of exact derivatives evaluate to the boundary values at $x \to \pm\infty$:

$$[u_x u_z]_{-\infty}^{\infty} + \int_{-\infty}^{\infty} \partial_z (u_z^2) dx = \nu [\partial_x u_z]_{-\infty}^{\infty}$$

Physically, far away from the center of the jet ($x \to \pm\infty$), the fluid is at rest. This requires that the longitudinal velocity vanishes ($u_z \to 0$) and the viscous shear stress vanishes ($\partial_x u_z \to 0$). Therefore, the boundary evaluation terms are completely zero:

$$0 + \int_{-\infty}^{\infty} \partial_z (u_z^2) dx = 0$$

Using Leibniz's rule, we can pull the $z$-derivative outside the integral since the integration limits are constants:

$$\frac{d}{dz} \int_{-\infty}^{\infty} u_z^2 dx = 0$$

Because the derivative with respect to $z$ is zero, the integral itself must be a constant. This proves that the specific kinematic momentum flux $F$ is conserved along the jet axis:

$$F = \int_{-\infty}^{\infty} u_z^2 dx = \text{constant}$$

### 2. The Stream Function and Incompressibility

In 2D incompressible flows, it is mathematically convenient to replace the two coupled velocity components ($u_x, u_z$) with a single scalar function called the **stream function**, denoted as $\psi$.

The problem defines the velocity components as the spatial derivatives of the stream function:

$$u_x = -\partial_z \psi$$

$$u_z = \partial_x \psi$$

To show that this definition automatically satisfies the incompressibility condition ($\partial_x u_x + \partial_z u_z = 0$), we substitute these expressions directly into the continuity equation:

$$\partial_x (-\partial_z \psi) + \partial_z (\partial_x \psi) = 0$$

$$-\partial_{xz} \psi + \partial_{zx} \psi = 0$$

By Clairaut's theorem (or Schwarz's theorem) of multivariable calculus, as long as the stream function $\psi$ has continuous second partial derivatives, the mixed partial derivatives commute ($\partial_{xz} \psi = \partial_{zx} \psi$).

Therefore, the two terms perfectly cancel each other out, confirming that any flow field derived from a valid stream function is mathematically guaranteed to be incompressible.

**3. Scaling Analysis for the Self-Similar Solution** We are given the stream function ansatz $\psi(x,z) = C_0 \sigma^\alpha(z) f(\xi)$, where $\xi = \frac{x}{\sigma(z)}$. First, compute the longitudinal velocity $u_z$ by taking the transverse derivative:

$$u_z = \partial_x \psi = C_0 \sigma^\alpha f'(\xi) \frac{\partial \xi}{\partial x} = C_0 \sigma^\alpha f'(\xi) \left( \frac{1}{\sigma} \right) = C_0 \sigma^{\alpha - 1} f'(\xi)$$

Substitute this velocity into the momentum flux integral $F = \int_{-\infty}^{\infty} u_z^2 \, dx$:

$$F = \int_{-\infty}^{\infty} \left( C_0 \sigma^{\alpha - 1} f'(\xi) \right)^2 dx$$

Transform the integration variable from $dx$ to $d\xi$ using $dx = \sigma \, d\xi$:

$$F = C_0^2 \sigma^{2\alpha - 2} \int_{-\infty}^{\infty} f'^2(\xi) \sigma \, d\xi = C_0^2 \sigma^{2\alpha - 1} \int_{-\infty}^{\infty} f'^2(\xi) \, d\xi$$

The problem provides the normalization condition $\int_{-\infty}^{\infty} f'^2(\xi) \, d\xi = 1$. Applying this reduces the integral to:

$$F = C_0^2 \sigma^{2\alpha - 1}$$

Because momentum flux $F$ is a strict constant independent of the longitudinal coordinate $z$, the $z$-dependent term $\sigma(z)$ must vanish from this expression. This requires the exponent of $\sigma$ to be zero:

$$2\alpha - 1 = 0 \implies \alpha = 1/2$$

Substituting $\alpha = 1/2$ leaves $F = C_0^2$, which forces the constant to be $C_0 = \sqrt{F}$.

**4. Determining the Velocity Fields $u_x$ and $u_z$**

Substituting the constants back into the ansatz yields the fully scaled stream function:

$$\psi(x,z) = \sqrt{F} \sigma^{1/2} f(\xi)$$

The longitudinal velocity $u_z$ was derived in the previous step:

$$u_z = \partial_x \psi = \sqrt{F} \sigma^{-1/2} f'(\xi)$$

The transverse velocity $u_x$ is found by taking the negative $z$-derivative of the stream function ($u_x = -\partial_z \psi$). We apply the product rule to $\sigma^{1/2} f(\xi)$ and the chain rule to the internal variable $\xi$.

Note that $\partial_z \xi = \partial_z (x \sigma^{-1}) = -x \sigma^{-2} \sigma' = -\xi \sigma^{-1} \sigma'$.

$$u_x = -\partial_z \left[ \sqrt{F} \sigma^{1/2} f(\xi) \right] = -\sqrt{F} \left[ \frac{1}{2} \sigma^{-1/2} \sigma' f(\xi) + \sigma^{1/2} f'(\xi) \partial_z \xi \right]$$

$$u_x = -\sqrt{F} \left[ \frac{1}{2} \sigma^{-1/2} \sigma' f(\xi) - \sigma^{-1/2} \sigma' \xi f'(\xi) \right]$$

Factoring out the common terms yields the final expression for transverse velocity:

$$u_x = \sqrt{F} \sigma^{-1/2} \sigma' \left[ \xi f'(\xi) - \frac{1}{2} f(\xi) \right]$$

**5. Deriving the Convective and Viscous Momentum Terms**

To build the convective acceleration term $u_x \partial_x u_z + u_z \partial_z u_z$, we first need the spatial derivatives of $u_z$:

- Transverse derivative:
    
    $$\partial_x u_z = \partial_x \left( \sqrt{F} \sigma^{-1/2} f' \right) = \sqrt{F} \sigma^{-1/2} f'' \left( \frac{1}{\sigma} \right) = \sqrt{F} \sigma^{-3/2} f''$$
    
- Longitudinal derivative:
    
    $$\partial_z u_z = \partial_z \left( \sqrt{F} \sigma^{-1/2} f' \right) = \sqrt{F} \left[ -\frac{1}{2} \sigma^{-3/2} \sigma' f' + \sigma^{-1/2} f'' (-\xi \sigma^{-1} \sigma') \right] = -\sqrt{F} \sigma^{-3/2} \sigma' \left[ \frac{1}{2} f' + \xi f'' \right]$$
    

Now, multiply these by the respective velocity fields:

$$u_x \partial_x u_z = \left( \sqrt{F} \sigma^{-1/2} \sigma' \left[ \xi f' - \frac{1}{2} f \right] \right) \left( \sqrt{F} \sigma^{-3/2} f'' \right) = F \sigma^{-2} \sigma' \left[ \xi f' f'' - \frac{1}{2} f f'' \right]$$

$$u_z \partial_z u_z = \left( \sqrt{F} \sigma^{-1/2} f' \right) \left( -\sqrt{F} \sigma^{-3/2} \sigma' \left[ \frac{1}{2} f' + \xi f'' \right] \right) = F \sigma^{-2} \sigma' \left[ -\frac{1}{2} f'^2 - \xi f' f'' \right]$$

Summing these two components causes the $\xi f' f''$ cross-terms to perfectly cancel out:

$$u_x \partial_x u_z + u_z \partial_z u_z = F \sigma^{-2} \sigma' \left[ -\frac{1}{2} f f'' - \frac{1}{2} f'^2 \right] = -F \frac{\sigma'}{\sigma^2} \frac{1}{2} \left( f f'' + f'^2 \right)$$

Using the derivative identity $(f^2)'' = \frac{d}{d\xi}(2f f') = 2f f'' + 2f'^2$, this convective term simplifies directly to the requested form:

$$-F \frac{\sigma'}{\sigma^2} \frac{1}{4} (f^2)''$$

For the viscous diffusion term on the right-hand side of the momentum equation, we take the transverse derivative of $\partial_x u_z$:

$$\partial_x (\nu \partial_x u_z) = \nu \partial_x \left( \sqrt{F} \sigma^{-3/2} f'' \right) = \nu \sqrt{F} \sigma^{-3/2} f''' \left( \frac{1}{\sigma} \right)$$

This yields the exact form required:

$$\partial_x (\nu \partial_x u_z) = \nu \sqrt{F} \frac{1}{\sigma^{5/2}} f'''$$

### 6. Defining the Jet Width $\sigma(z)$

Equating the convective acceleration term and the viscous diffusion term from equations (21) and (22) yields the full momentum balance:

$$-F \frac{\sigma'}{\sigma^2} \frac{1}{4} (f^2)'' = \nu \sqrt{F} \frac{1}{\sigma^{5/2}} f'''$$

To find a self-similar solution, the $z$-dependent terms ($\sigma$) must separate entirely from the $\xi$-dependent terms ($f$). Rearranging the equation to group the variables:

$$-\frac{\sqrt{F}}{4\nu} \sigma^{1/2} \sigma' = \frac{f'''}{(f^2)''}$$

For this equality to hold for all $z$ and $\xi$, both sides must equal a dimensionless separation constant. To obtain the simplest possible ordinary differential equation (ODE) for $f(\xi)$, namely $f''' + (f^2)'' = 0$, we must set the left side to $-1/4$. This is equivalent to setting the constant to 6 in the final $\sigma$ expression.

Let us enforce this balance:

$$\frac{\sqrt{F}}{4\nu} \sigma^{1/2} \frac{d\sigma}{dz} = \frac{1}{4}$$

$$\sigma^{1/2} d\sigma = \frac{4\nu}{\sqrt{F}} \frac{1}{4} dz \implies \sigma^{1/2} d\sigma = \frac{\nu}{\sqrt{F}} dz$$

Wait, balancing the factor of 4:

$$-\frac{\sqrt{F}}{4\nu} \sigma^{1/2} \sigma' (f^2)'' = f'''$$

To match $f''' + (f^2)'' = 0$, we need $\frac{\sqrt{F}}{4\nu} \sigma^{1/2} \sigma' = 1$, so $\sigma^{1/2} \sigma' = \frac{4\nu}{\sqrt{F}}$.

Integrating with respect to $z$, assuming the jet originates as a point source at $z=0$ ($\sigma(0) = 0$):

$$\int_{0}^{\sigma} \sigma^{1/2} d\sigma = \int_{0}^{z} \frac{4\nu}{\sqrt{F}} dz$$

$$\frac{2}{3} \sigma^{3/2} = \frac{4\nu z}{\sqrt{F}} \implies \sigma^{3/2} = \frac{6\nu z}{\sqrt{F}}$$

Solving for $\sigma(z)$ retrieves the exact scaling requested:

$$\sigma(z) = \left( \frac{6\nu z}{\sqrt{F}} \right)^{2/3}$$

### 7. Solving the Similarity Profile $f(\xi)$

With the chosen separation constant, the governing PDE collapses into the non-linear ODE:

$$f''' + (f^2)'' = 0$$

Integrate once with respect to $\xi$:

$$f'' + (f^2)' = C_1$$

Because the longitudinal velocity $u_z \propto f'$ and viscous shear $u_{zx} \propto f''$ must decay to zero at the far edges of the jet ($\xi \to \pm \infty$), the constant $C_1 = 0$. Using $(f^2)' = 2ff'$:

$$f'' + 2f f' = 0$$

Integrate a second time:

$$f' + f^2 = c^2$$

Here, $c^2$ is an integration constant. Because the stream function is antisymmetric across the center of the jet, $f(0) = 0$. As $\xi \to \pm \infty$, $f' \to 0$, which requires $f \to \pm c$.

This is a separable first-order ODE:

$$\frac{df}{c^2 - f^2} = d\xi$$

Integrating using the identity $\int \frac{dx}{1-x^2} = \text{artanh}(x)$:

$$\frac{1}{c} \text{artanh}\left(\frac{f}{c}\right) = \xi \implies f(\xi) = c \tanh(c\xi)$$

To determine the constant $c$, we apply the integral normalization condition:

$$\int_{-\infty}^{\infty} f'^2(\xi) d\xi = 1$$

Compute the derivative $f'(\xi) = c^2 \text{sech}^2(c\xi)$ and substitute:

$$\int_{-\infty}^{\infty} c^4 \text{sech}^4(c\xi) d\xi = 1$$

Let $u = c\xi$, so $d\xi = du / c$:

$$c^3 \int_{-\infty}^{\infty} \text{sech}^4(u) du = c^3 \int_{-\infty}^{\infty} (1 - \tanh^2 u) \text{sech}^2 u \, du = 1$$

$$c^3 \left[ \tanh(u) - \frac{1}{3}\tanh^3(u) \right]_{-\infty}^{\infty} = c^3 \left[ \left(1 - \frac{1}{3}\right) - \left(-1 + \frac{1}{3}\right) \right] = c^3 \left(\frac{4}{3}\right) = 1$$

$$c^3 = \frac{3}{4} \implies c = \left( \frac{3}{4} \right)^{1/3}$$

### 8. Deducing the Final Solution Fields

**The Stream Function ($\psi$):**

Using $\psi = \sqrt{F} \sigma^{1/2} f(\xi)$, we substitute $\sigma(z)$ and $f(\xi)$:

$$\sigma^{1/2} = \left( \frac{6\nu z}{\sqrt{F}} \right)^{1/3} = 6^{1/3} \nu^{1/3} z^{1/3} F^{-1/6}$$

$$\psi = F^{1/2} \left( 6^{1/3} \nu^{1/3} z^{1/3} F^{-1/6} \right) \left[ \left( \frac{3}{4} \right)^{1/3} \tanh(c\xi) \right] = F^{1/3} \nu^{1/3} z^{1/3} \left( 6 \cdot \frac{3}{4} \right)^{1/3} \tanh(c\xi)$$

$$\psi = \left( \frac{9}{2} \right)^{1/3} \left( F \nu z \right)^{1/3} \tanh\left[ \left( \frac{3}{4} \right)^{1/3} \xi \right]$$

**The Longitudinal Velocity ($u_z$):**

Using $u_z = \sqrt{F} \sigma^{-1/2} f'(\xi)$:

$$f'(\xi) = c^2 \text{sech}^2(c\xi) = \left( \frac{3}{4} \right)^{2/3} \frac{1}{\cosh^2(c\xi)}$$

$$\sigma^{-1/2} = \left( \frac{6\nu z}{\sqrt{F}} \right)^{-1/3} = 6^{-1/3} \nu^{-1/3} z^{-1/3} F^{1/6}$$

$$u_z = F^{1/2} \left( 6^{-1/3} \nu^{-1/3} z^{-1/3} F^{1/6} \right) \left( \frac{3}{4} \right)^{2/3} \frac{1}{\cosh^2(c\xi)} = F^{2/3} \nu^{-1/3} z^{-1/3} \left( \frac{1}{6} \cdot \frac{9}{16} \right)^{1/3} \frac{1}{\cosh^2(c\xi)}$$

$$u_z = \left( \frac{3}{32} \right)^{1/3} \left( \frac{F^2}{\nu z} \right)^{1/3} \frac{1}{\cosh^2\left[ \left( \frac{3}{4} \right)^{1/3} \xi \right]}$$

**The Transverse Velocity ($u_x$):**

Using $u_x = \sqrt{F} \sigma^{-1/2} \sigma' \left[ \xi f'(\xi) - \frac{1}{2} f(\xi) \right]$.

From the Step 6 derivation, we established $\sigma^{1/2} \sigma' = \frac{4\nu}{\sqrt{F}}$, meaning $\sqrt{F} \sigma^{-1/2} \sigma' = 4\nu \sigma^{-1}$.

$$4\nu \sigma^{-1} = 4\nu \left( \frac{\sqrt{F}}{6\nu z} \right)^{2/3} = 4 \nu^{1/3} F^{1/3} 6^{-2/3} z^{-2/3} = \left( \frac{F \nu}{z^2} \right)^{1/3} \frac{4}{36^{1/3}}$$

Next, rewrite the bracketed term using $c$:

$$\xi f' - \frac{1}{2} f = \xi c^2 \text{sech}^2(c\xi) - \frac{1}{2} c \tanh(c\xi) = \frac{c}{2} \left[ 2c\xi \frac{1}{\cosh^2(c\xi)} - \tanh(c\xi) \right]$$

Multiply the prefactors together:

$$\left( \frac{F \nu}{z^2} \right)^{1/3} \frac{4}{36^{1/3}} \cdot \frac{c}{2} = \left( \frac{F \nu}{z^2} \right)^{1/3} \left( \frac{8}{36} \right)^{1/3} \left( \frac{3}{4} \right)^{1/3} = \left( \frac{F \nu}{z^2} \right)^{1/3} \left( \frac{24}{144} \right)^{1/3} = \left( \frac{F \nu}{z^2} \right)^{1/3} \frac{1}{6^{1/3}}$$

Combining this prefactor with the bracketed terms yields the exact final formula:

$$u_x = \frac{1}{6^{1/3}} \left( \frac{F \nu}{z^2} \right)^{1/3} \left\{ \left( \frac{3}{4} \right)^{1/3} 2\xi \frac{1}{\cosh^2\left[ \left( \frac{3}{4} \right)^{1/3} \xi \right]} - \tanh\left[ \left( \frac{3}{4} \right)^{1/3} \xi \right] \right\}$$



# 20. — Bullard’s dynamo (1955)

### The Physics of Bullard's Homopolar Dynamo

Bullard's dynamo is a foundational electro-mechanical model designed to demonstrate how a self-exciting magnetic field can be sustained by fluid motion, serving as a simplified analog for the geodynamo in Earth's liquid iron core.

The system consists of a conducting disk rotating in a magnetic field. According to Faraday's law of induction, the radial motion of the disk's charges through the magnetic field generates an electromotive force (EMF) between the center and the rim. By connecting a wire from the rim back to the center via a coil, this EMF drives a current. The coil is arranged so that this very current produces the magnetic field the disk is rotating in, creating a powerful non-linear feedback loop.

(Note: Questions 1 & 2 are duplicated as 3 & 4 in the provided problem sheet. They are addressed together below).

### 1 & 3. Origin of the Terms in the Coupled Equations

The system is governed by two coupled non-linear differential equations representing the electrical and mechanical energy balances.

**The Electrical Equation:** $L \frac{di}{dt} = -Ri + M\omega i$

This is Kirchhoff's voltage law for the closed circuit.

- **$L \frac{di}{dt}$**: The voltage drop across the inductor, representing the self-inductance of the wire loop resisting changes in current.
    
- **$-Ri$**: The Ohmic dissipation. This is the voltage lost to the electrical resistance of the wire, disk, and sliding frictional contacts.
    
- **$M\omega i$**: The motional EMF (Faraday induction). The disk rotates at angular velocity $\omega$ through a magnetic field proportional to the current $i$. The parameter $M$ represents the mutual induction linking the mechanical rotation to the generated voltage. This is the energy source term for the electrical circuit.
    

**The Mechanical Equation:** $J \frac{d\omega}{dt} = \Gamma - Mi^2$

This is Newton's second law for rotation ($\sum \tau = J \alpha$).

- **$J \frac{d\omega}{dt}$**: The rate of change of the disk's angular momentum, where $J$ is the moment of inertia.
    
- **$\Gamma$**: The constant external driving torque applied to the shaft (representing the thermal convection forces in the Earth's core maintaining fluid motion).
    
- **$-Mi^2$**: The Lorentz braking torque (Laplace force). When current flows radially through the disk in the presence of the vertical magnetic field, it experiences a force perpendicular to both, which opposes the rotation of the disk according to Lenz's law. Because the magnetic field is proportional to $i$ and the radial current is exactly $i$, the opposing torque scales as $i \times i = i^2$.
    

### 2 & 4. Determining the Fixed Points ($\Gamma > 0$)

A fixed point (steady state) occurs when both the current and the angular velocity are constant, meaning $\frac{di}{dt} = 0$ and $\frac{d\omega}{dt} = 0$.

Setting the electrical equation to zero:

$0 = -Ri + M\omega i \implies i(M\omega - R) = 0$

This gives two possibilities: either $i = 0$ or $\omega = R/M$.

Setting the mechanical equation to zero: $0 = \Gamma - Mi^2 \implies Mi^2 = \Gamma$ Because we are given $\Gamma > 0$, the current must be non-zero to balance the external torque: $i^2 = \frac{\Gamma}{M} \implies i = \pm \sqrt{\frac{\Gamma}{M}}$

Since $i \neq 0$, the electrical equation forces $\omega = R/M$.

Therefore, the system has exactly two fixed points, differing only by the polarity of the current (the direction of the magnetic field):

$$(i^*, \omega^*) = \left( \sqrt{\frac{\Gamma}{M}}, \frac{R}{M} \right) \quad \text{and} \quad \left( -\sqrt{\frac{\Gamma}{M}}, \frac{R}{M} \right)$$

### 5. Change of Variables and the Parameter $a$

To simplify the stability analysis, we non-dimensionalize the equations to match the target form:

$$\frac{d\Omega}{d\tau} = 1 - I^2$$

$$\frac{dI}{d\tau} = a(\Omega - 1)I$$

**Step 1: Scaling the variables**

We define dimensionless variables based on the steady-state values found in the previous step:

- Normalized current: $I = \frac{i}{i^*} = i \sqrt{\frac{M}{\Gamma}}$
    
- Normalized angular velocity: $\Omega = \frac{\omega}{\omega*} = \omega \frac{M}{R}$
    
- Normalized time: $\tau = \frac{t}{t_0}$ (where $t_0$ is a characteristic time scale to be determined).
    

Substitute $i = I \sqrt{\frac{\Gamma}{M}}$ and $\omega = \Omega \frac{R}{M}$ into the mechanical equation:

$$J \frac{d}{dt} \left( \Omega \frac{R}{M} \right) = \Gamma - M \left( I \sqrt{\frac{\Gamma}{M}} \right)^2$$

$$J \frac{R}{M} \frac{d\Omega}{dt} = \Gamma - \Gamma I^2 = \Gamma(1 - I^2)$$

Introduce the dimensionless time $\tau = t / t_0$, meaning $\frac{d}{dt} = \frac{1}{t_0} \frac{d}{d\tau}$:

$$\left( \frac{J R}{M \Gamma t_0} \right) \frac{d\Omega}{d\tau} = 1 - I^2$$

For this to exactly match the target equation $\frac{d\Omega}{d\tau} = 1 - I^2$, the term in parentheses must equal 1. This uniquely determines the time scale $t_0$:

$$t_0 = \frac{J R}{M \Gamma}$$

**Step 2: Deriving the parameter $a$**

Substitute the scalings and the time derivative $\frac{d}{dt} = \frac{M \Gamma}{J R} \frac{d}{d\tau}$ into the electrical equation:

$$L \left( \frac{M \Gamma}{J R} \right) \frac{d}{d\tau} \left( I \sqrt{\frac{\Gamma}{M}} \right) = -R \left( I \sqrt{\frac{\Gamma}{M}} \right) + M \left( \Omega \frac{R}{M} \right) \left( I \sqrt{\frac{\Gamma}{M}} \right)$$

Divide the entire equation by $\sqrt{\frac{\Gamma}{M}}$:

$$\frac{L M \Gamma}{J R} \frac{dI}{d\tau} = -R I + R \Omega I = R(\Omega - 1)I$$

Isolate $\frac{dI}{d\tau}$:

$$\frac{dI}{d\tau} = \left( \frac{J R^2}{L M \Gamma} \right) (\Omega - 1)I$$

Matching this to the target form $\frac{dI}{d\tau} = a(\Omega - 1)I$, we find the explicit expression for $a$:

$$a = \frac{J R^2}{L M \Gamma}$$

**Physical Interpretation of $a$:**

The parameter $a$ is a dimensionless ratio of the two fundamental relaxation timescales of the system.

- The **mechanical timescale** ($t_m$) describes how quickly the disk accelerates under torque: $t_m = \frac{J \omega^*}{\Gamma} = \frac{J (R/M)}{\Gamma} = \frac{J R}{M \Gamma}$.
    
- The **electrical timescale** ($t_e$) describes how quickly current decays in the RL circuit: $t_e = \frac{L}{R}$.
    

By dividing the two, we recover $a$:

$$\frac{t_m}{t_e} = \frac{\frac{J R}{M \Gamma}}{\frac{L}{R}} = \frac{J R^2}{L M \Gamma} = a$$

Therefore, $a$ represents the ratio of the mechanical response time to the electromagnetic response time. If $a \gg 1$, the mechanics respond very slowly compared to the rapid fluctuations of the electrical circuit.

### 6. Stability of the Fixed Point in the Non-Linear Regime

To study the stability of the non-dimensionalized system at the fixed point $(I, \Omega) = (1, 1)$, we first evaluate the Jacobian matrix of the linearized system. The governing equations are:

$$f(I, \Omega) = \frac{d\Omega}{d\tau} = 1 - I^2$$

$$g(I, \Omega) = \frac{dI}{d\tau} = a(\Omega - 1)I$$

The Jacobian matrix $J$ is:

$$J(I, \Omega) = \begin{pmatrix} \frac{\partial f}{\partial \Omega} & \frac{\partial f}{\partial I} \\ \frac{\partial g}{\partial \Omega} & \frac{\partial g}{\partial I} \end{pmatrix} = \begin{pmatrix} 0 & -2I \\ aI & a(\Omega - 1) \end{pmatrix}$$

Evaluating this at the fixed point $(1, 1)$ gives:

$$J(1, 1) = \begin{pmatrix} 0 & -2 \\ a & 0 \end{pmatrix}$$

The eigenvalues $\lambda$ are found by solving the characteristic equation $\text{Det}(J - \lambda \mathbf{I}) = 0$:

$$\lambda^2 - \text{Tr}(J)\lambda + \text{Det}(J) = 0 \implies \lambda^2 + 2a = 0 \implies \lambda = \pm i\sqrt{2a}$$

Because the eigenvalues are purely imaginary, linear stability analysis classifies the fixed point as a **center** (marginally stable). However, according to the Hartman-Grobman theorem, a center is a borderline, structurally unstable case in linear theory; higher-order non-linear terms could cause the trajectories to weakly spiral inward (stable focus) or outward (unstable focus). Therefore, linear stability analysis alone is strictly **not valid** for concluding the true non-linear behavior. We must find a conserved quantity to prove the trajectories remain closed.

### 7. The First Integral and Trajectory Properties

To find the exact phase-space trajectories, we eliminate time by dividing the two differential equations to find $\frac{d\Omega}{dI}$:

$$\frac{d\Omega}{dI} = \frac{\frac{d\Omega}{d\tau}}{\frac{dI}{d\tau}} = \frac{1 - I^2}{a(\Omega - 1)I}$$

This is a separable differential equation. Grouping the variables yields:

$$a(\Omega - 1) d\Omega = \left(\frac{1}{I} - I\right) dI$$

Integrating both sides produces:

$$a\left(\frac{1}{2}\Omega^2 - \Omega\right) = \ln\vert{}I\vert{} - \frac{1}{2}I^2 + C$$

Rearranging this to group all state variables on one side gives the first integral (a conserved "energy" function) of the form $g(\Omega) + h(I) = K$:

$$a\left(\frac{1}{2}\Omega^2 - \Omega\right) + \frac{1}{2}I^2 - \ln\vert{}I\vert{} = K$$

The constant $K$ is determined entirely by the initial conditions $\Omega_0$ and $I_0$:

$$K = a\left(\frac{1}{2}\Omega_0^2 - \Omega_0\right) + \frac{1}{2}I_0^2 - \ln\vert{}I_0\vert{}$$

**Important Property:** Because the system possesses this exact first integral, it is a conservative system. The trajectories $(I(t), \Omega(t))$ are constrained to follow the contour lines of this constant function $K$. These contours form strictly closed, nested loops around the fixed point $(1,1)$. This mathematically proves that the non-linear system does indeed behave as a perfect center, undergoing perpetual, undamped oscillations.

### 8. Phase Plane Trajectories Sketch

In the $(I, \Omega)$ plane, the trajectories look like deformed ellipses (or kidney-bean shapes) orbiting the fixed points at $(1, 1)$ and $(-1, 1)$. Because of the non-linear $\ln\vert{}I\vert{}$ term, the loops become highly asymmetric for large initial perturbations, bulging outward while being strictly bounded away from the $I = 0$ vertical axis.

### 9. Magnetic Field Reversals

A magnetic field reversal requires the current $i$ (and thus $I$) to change sign, which physically means the trajectory must cross the $I = 0$ axis. Looking at the conserved quantity, as $I \to 0$, the term $\ln\vert{}I\vert{} \to -\infty$. For the trajectory to reach $I = 0$, the function would require an infinite amount of "energy" ($K \to \infty$). Since the initial conditions dictate a finite constant $K$, the trajectory is enclosed by an infinite potential barrier at the $I = 0$ axis. Consequently, this idealized version of Bullard's dynamo **cannot** produce a reversal of the current; a dynamo starting with a positive field will oscillate but remain strictly positive forever.

### Effect of Low Friction on the Phase Portrait

By introducing mechanical friction, the system becomes dissipative, altering the mechanics to:

$$J \frac{d\omega}{dt} = \Gamma - \alpha \omega - M i^2$$

(The electrical equation $L \frac{di}{dt} = -Ri + M\omega i$ remains unchanged).

#### 1. Threshold Torque ($\Gamma_c$) and Fixed Points

To find the steady states, we set the time derivatives to zero:

$0 = -Ri + M\omega i \implies i = 0 \text{ or } \omega = \frac{R}{M}$

- **Case 1 (Trivial State):** If $i = 0$, the mechanical equation gives $0 = \Gamma - \alpha \omega \implies \omega = \frac{\Gamma}{\alpha}$.
    
- **Case 2 (Dynamo State):** If $\omega = \frac{R}{M}$, substitute this into the mechanical equation:
    
    $$0 = \Gamma - \alpha \left(\frac{R}{M}\right) - M i^2 \implies M i^2 = \Gamma - \frac{\alpha R}{M}$$
    
    For the dynamo to operate, there must be a real, non-zero solution for the current $i$. This requires $M i^2 > 0$, meaning:
    
    $$\Gamma > \frac{\alpha R}{M}$$
    
    Therefore, the critical threshold torque required to observe a non-zero current is **$\Gamma_c = \frac{\alpha R}{M}$**.
    

Beyond $\Gamma_c$, there are **three** valid mathematical solutions:

1. The non-dynamo state: $(i, \omega) = (0, \Gamma/\alpha)$
    
2. Positive dynamo state: $(i, \omega) = \left( +\sqrt{\frac{\Gamma - \Gamma_c}{M}}, \frac{R}{M} \right)$
    
3. Negative dynamo state: $(i, \omega) = \left( -\sqrt{\frac{\Gamma - \Gamma_c}{M}}, \frac{R}{M} \right)$
    

#### 2. Stability of the Fixed Points ($\Gamma > \Gamma_c$)

We evaluate the Jacobian for the new frictional system:

$$J(i, \omega) = \begin{pmatrix} \frac{1}{L}(-R + M\omega) & \frac{M}{L}i \\ -\frac{2M}{J}i & -\frac{\alpha}{J} \end{pmatrix}$$

- **At the trivial state $(0, \Gamma/\alpha)$:**
    
    $$J = \begin{pmatrix} \frac{1}{L}(-R + M \frac{\Gamma}{\alpha}) & 0 \\ 0 & -\frac{\alpha}{J} \end{pmatrix}$$
    
    The eigenvalues are the diagonal elements. The second eigenvalue is $\lambda_2 = \frac{M}{L\alpha}(\Gamma - \Gamma_c)$. Since $\Gamma > \Gamma_c$, $\lambda_2$ is strictly positive. The trivial state is a **saddle point (unstable)**.
    
- **At the dynamo states $(\pm i^*, R/M)$:**
    
    $$J = \begin{pmatrix} 0 & \pm\frac{M}{L}i^* \\ \mp\frac{2M}{J}i^* & -\frac{\alpha}{J} \end{pmatrix}$$
    
    The trace is $\text{Tr}(J) = -\frac{\alpha}{J}$ (which is strictly negative due to friction).
    
    The determinant is $\text{Det}(J) = 0 - \left(\pm\frac{M}{L}i^*\right)\left(\mp\frac{2M}{J}i^*\right) = \frac{2M^2}{JL}(i^*)^2$ (which is strictly positive).
    
    Because the trace is negative and the determinant is positive, both dynamo fixed points are **stable** (either stable spirals or stable nodes, depending on the magnitude of the friction).
    

**Phase Portrait Sketch Description:**

The $(i, \omega)$ plane features an unstable saddle at the origin (current $i=0$). Trajectories originating near the $i=0$ axis are repelled away from it and spiral inward toward one of the two stable fixed points at $\pm i^*$. The system now exhibits damped oscillations, eventually settling into a steady, constant magnetic field.

#### 3. Bifurcation Diagram (Bullard with Friction)

Plotting the equilibrium current $i$ as a function of the motor torque $\Gamma$ reveals a classic **supercritical pitchfork bifurcation**:

- For $0 < \Gamma \le \Gamma_c$, there is a single solid horizontal line along $i = 0$ representing the only stable state.
    
- At $\Gamma = \Gamma_c$, the $i = 0$ line becomes dashed (indicating it has lost stability).
    
- Two new solid branches emerge symmetrically from the critical point, growing outwards to the right as parabolic curves defined by $i = \pm \sqrt{\frac{1}{M}(\Gamma - \Gamma_c)}$.
- ![[Pasted image 20260925011517.png]]


# 21. — Rikitake’s dynamo (1956)

### 1. Conservation of the Rotational Difference

The coupled mechanical equations for the rotational speeds of the two dynamos are:

$$\frac{d\Omega_1}{d\tau} = 1 - I_1 I_2$$

$$\frac{d\Omega_2}{d\tau} = 1 - I_1 I_2$$

Subtracting the second equation from the first eliminates the non-linear coupling term:

$$\frac{d\Omega_1}{d\tau} - \frac{d\Omega_2}{d\tau} = (1 - I_1 I_2) - (1 - I_1 I_2) = 0$$

$$\frac{d}{d\tau}(\Omega_1 - \Omega_2) = 0$$

Integrating this with respect to time shows that the difference in angular velocities is a strict constant. We define this constant as $2\Omega_0$. If we set $\Omega_1 = \Omega$, it immediately follows that $\Omega_2 = \Omega - 2\Omega_0$. This conservation law elegantly reduces the 4D dynamical system into a 3D system.

### 2. Determining the Fixed Points (Rikitake System)

To find the steady states (fixed points), we set the time derivatives of the reduced 3D system to zero:

1. **Mechanical equilibrium:**
    
    $$\frac{d\Omega}{d\tau} = 1 - I_1 I_2 = 0 \implies I_1 I_2 = 1 \implies I_2 = \frac{1}{I_1}$$
    
    Let the equilibrium current of the first dynamo be $I_1 = I^*$. Consequently, $I_2 = \frac{1}{I^*}$.
    
2. **Electrical equilibrium (Dynamo 1):**
    
    $$\frac{dI_1}{d\tau} = -a I_1 + a \Omega_1 I_2 = 0 \implies I_1 = \Omega I_2$$
    
    Substitute the values $I_1 = I^*$ and $I_2 = 1/I^*$:
    
    $$I^* = \Omega^* \left(\frac{1}{I^*}\right) \implies \Omega^* = (I^*)^2$$
    
3. **Electrical equilibrium (Dynamo 2):**
    
    $$\frac{dI_2}{d\tau} = -a I_2 + a \Omega_2 I_1 = 0 \implies I_2 = \Omega_2 I_1$$
    
    Substitute $I_2 = 1/I^*$, $I_1 = I^*$, and $\Omega_2 = \Omega^* - 2\Omega_0$:
    
    $$\frac{1}{I^*} = (\Omega^* - 2\Omega_0) I^* \implies \Omega^* - 2\Omega_0 = \frac{1}{(I^*)^2}$$
    
    Because we previously found $\Omega^* = (I^*)^2$, this translates to $\Omega^* - 2\Omega_0 = \frac{1}{\Omega^*}$.
    

This confirms the fixed point takes the exact form $X^* = \left(I^*, \frac{1}{I^*}, \Omega^*\right)$.

### 3. The Jacobian Matrix at $X^*$

The Jacobian matrix $J$ linearizes the system around the equilibrium. We take the partial derivatives of the flow vector field $f = (\dot{I_1}, \dot{I_2}, \dot{\Omega})$ with respect to the state variables $(I_1, I_2, \Omega)$:

$$J = \begin{pmatrix} \frac{\partial \dot{I_1}}{\partial I_1} & \frac{\partial \dot{I_1}}{\partial I_2} & \frac{\partial \dot{I_1}}{\partial \Omega} \\ \frac{\partial \dot{I_2}}{\partial I_1} & \frac{\partial \dot{I_2}}{\partial I_2} & \frac{\partial \dot{I_2}}{\partial \Omega} \\ \frac{\partial \dot{\Omega}}{\partial I_1} & \frac{\partial \dot{\Omega}}{\partial I_2} & \frac{\partial \dot{\Omega}}{\partial \Omega} \end{pmatrix} = \begin{pmatrix} -a & a\Omega & aI_2 \\ a(\Omega - 2\Omega_0) & -a & aI_1 \\ -I_2 & -I_1 & 0 \end{pmatrix}$$

Evaluating this matrix at the fixed point $X^*$ using our derived relations $\Omega^* = (I^*)^2$ and $(\Omega^* - 2\Omega_0) = \frac{1}{(I^*)^2}$ yields:

$$J(X^*) = \begin{pmatrix} -a & a(I^*)^2 & a/I^* \\ a/(I^*)^2 & -a & aI^* \\ -1/I^* & -I^* & 0 \end{pmatrix}$$

### 4. Eigenvalue Analysis and Dynamics

Expanding the determinant along the third row:

$$\left( -\frac{1}{k} \right) \begin{vmatrix} a k^2 & \frac{a}{k} \\ -a - \lambda & a k \end{vmatrix} - (-k) \begin{vmatrix} -a - \lambda & \frac{a}{k} \\ \frac{a}{k^2} & a k \end{vmatrix} + (-\lambda) \begin{vmatrix} -a - \lambda & a k^2 \\ \frac{a}{k^2} & -a - \lambda \end{vmatrix} = 0$$

Evaluate each 2x2 determinant:

1. **First term:**
    
    $$-\frac{1}{k} \left[ a^2 k^3 - \frac{a}{k}(-a - \lambda) \right] = -\frac{1}{k} \left[ a^2 k^3 + \frac{a^2}{k} + \frac{a\lambda}{k} \right] = -a^2 k^2 - \frac{a^2}{k^2} - \frac{a\lambda}{k^2}$$
    
2. **Second term:**
    
    $$+k \left[ -a^2 k - a k \lambda - \frac{a^2}{k^3} \right] = -a^2 k^2 - a k^2 \lambda - \frac{a^2}{k^2}$$
    
3. **Third term:**
    
    $$-\lambda \left[ (a + \lambda)^2 - a^2 \right] = -\lambda \left( a^2 + 2a\lambda + \lambda^2 - a^2 \right) = -\lambda^2 (\lambda + 2a)$$
    

Combine the first and second terms:

$$\left( -a^2 k^2 - \frac{a^2}{k^2} - \frac{a\lambda}{k^2} \right) + \left( -a^2 k^2 - a k^2 \lambda - \frac{a^2}{k^2} \right) = -2a^2 \left( k^2 + \frac{1}{k^2} \right) - a\lambda \left( k^2 + \frac{1}{k^2} \right)$$

Define $M = a \left( k^2 + \frac{1}{k^2} \right) = a \left( (I^*)^2 + \frac{1}{(I^*)^2} \right)$. Substituting $M$ back into the full characteristic equation yields:

$$-2a M - M \lambda - \lambda^2 (\lambda + 2a) = 0$$

Group terms to factor by grouping:

$$-M(\lambda + 2a) - \lambda^2(\lambda + 2a) = 0 \implies -(\lambda + 2a)(\lambda^2 + M) = 0$$

Setting each factor to zero determines the three eigenvalues:

- **Real eigenvalue:**
    
    $$\lambda_1 = -2a$$
    
- **Purely imaginary conjugate eigenvalues:**
    
    $$\lambda^2 + M = 0 \implies \lambda_{\pm} = \pm i \sqrt{M} = \pm i \sqrt{a (I^*)^2 + \frac{a}{(I^*)^2}}$$


- **Conclusion on Dynamics:** The real eigenvalue $\lambda_1 = -2a$ is strictly negative (since $a > 0$), meaning the system is strongly attracting along one specific direction in phase space, continually squashing the volume of trajectories. However, the other two eigenvalues are purely imaginary. In linear stability theory, this indicates a local "center" (marginal stability). Because centers are structurally unstable, the non-linear terms dictate the true long-term behavior, which ultimately drives the trajectories into chaotic spirals moving outward from the fixed points.
-  Dissipation and Volume Contraction in Phase Space

	The divergence of the 3D flow vector field $\mathbf{f} = (\dot{I}_1, \dot{I}_2, \dot{\Omega})$ represents the rate of change of an infinitesimal phase space volume $V(t)$:
	
	$$\nabla \cdot \mathbf{f} = \frac{\partial \dot{I}_1}{\partial I_1} + \frac{\partial \dot{I}_2}{\partial I_2} + \frac{\partial \dot{\Omega}}{\partial \Omega} = -a + (-a) + 0 = -2a$$
	
	By Liouville's theorem, any initial volume in phase space contracts exponentially over time:
	
	$$V(t) = V(0) e^{-2a \tau}$$
	
	Because $a > 0$, the system is **strictly dissipative**. Trajectories originating anywhere in 3D phase space are rapidly compressed onto a subset of zero volume (an attractor). The trace of the Jacobian matrix equals the sum of the eigenvalues:
	
	$$\text{Tr}(J) = \lambda_1 + \lambda_+ + \lambda_- = -2a + i\sqrt{M} - i\sqrt{M} = -2a$$
* Stable Manifold ($\lambda_1 = -2a < 0$)

	The negative real eigenvalue $\lambda_1 = -2a$ provides a strong directional contraction along a 1D stable manifold. Any deviation from the fixed point in this direction decays exponentially as $e^{-2a\tau}$. Physically, this rapidly squashes 3D trajectories onto a 2D center manifold passing through $X^*$.

*  Non-Hyperbolicity and Structural Instability ($\lambda_{\pm} = \pm i \sqrt{M}$)

- **Linear Center Prediction:** In linear stability analysis, purely imaginary eigenvalues ($\text{Re}(\lambda_{\pm}) = 0$) correspond to a **center**. The linear approximation predicts that small perturbations around $X^*$ will orbit indefinitely in closed, concentric, non-decaying ellipses with natural frequency $\omega_0 = \sqrt{a (I^*)^2 + a / (I^*)^2}$.
    
- **Failure of Linearization:** Because the real parts of $\lambda_{\pm}$ are identically zero, $X^*$ is a **non-hyperbolic fixed point**. The Hartman-Grobman Theorem does not apply, meaning the linear approximation cannot determine the true stability of the non-linear system.
    

*  Non-Linear Destabilization and Chaotic Spirals

- **Energy Pumping by Non-Linear Terms:** In the full non-linear Rikitake system, higher-order non-linear terms act as a subtle feedback mechanism that pumps energy into the oscillatory modes.
    
- **Outward Spiraling:** Rather than remaining on closed ellipses, trajectories on the 2D surface slowly **spiral outward** away from $X^*$, with increasing oscillation amplitude over time.
    
- **Global Flips & Strange Attractor:** Due to the mirror symmetry of the equations, a twin fixed point $-X^* = \left(-I^*, -\frac{1}{I^*}, \Omega^*\right)$ exists with identical eigenvalues. As the outward spiral around $X^*$ grows sufficiently large, the trajectory escapes the local influence of $X^*$ and is captured by the neighborhood of $-X^*$. It then begins spiraling outward around $-X^*$ until it flips back unpredictably


    

### 5. Synchronization When $\Omega_0 = 0$

If the dynamos are perfectly identical without a rotational offset ($\Omega_0 = 0$), then $\Omega_1 = \Omega_2 = \Omega$. To show that $I_1^2 - I_2^2$ decays to zero, we take its time derivative:

$$\frac{d}{d\tau}(I_1^2 - I_2^2) = 2I_1 \frac{dI_1}{d\tau} - 2I_2 \frac{dI_2}{d\tau}$$

Substitute the differential equations for $\dot{I_1}$ and $\dot{I_2}$:

$$2I_1(-aI_1 + a\Omega I_2) - 2I_2(-aI_2 + a\Omega I_1)$$

$$= -2aI_1^2 + 2a\Omega I_1 I_2 + 2aI_2^2 - 2a\Omega I_1 I_2$$

$$= -2a(I_1^2 - I_2^2)$$

This is a linear decay equation. Its solution is $(I_1^2 - I_2^2) = C e^{-2a\tau}$. Because $a > 0$, the difference strictly goes to zero as $\tau \to \infty$. This proves that eventually, $I_1^2 = I_2^2$, meaning the dynamos must lock into either the $I_2 = I_1$ or $I_2 = -I_1$ mode.

- **Mode $I_2 = -I_1$:** The mechanical equation becomes $\frac{d\Omega}{d\tau} = 1 + I_1^2$, forcing the rotational speed to accelerate to infinity. This is physically impossible for a bounded steady state.
    
- **Mode $I_2 = I_1$:** The system perfectly collapses into $\frac{dI}{d\tau} = a(\Omega - 1)I$ and $\frac{d\Omega}{d\tau} = 1 - I^2$. This is exactly the single Bullard dynamo equation.
    
- **Deduction:** The two dynamos completely synchronize. Once locked, the system executes the stable, non-reversing periodic oscillations characteristic of a single Bullard dynamo.

### Question 6: Behavior of the System ($\Omega_0 \neq 0$)

Looking at the time series in Figure 4 for parameters $a=1$ and $\Omega_0=0.5$, we can observe the physical evolution of the current $I_1(t)$ and the angular velocity $\Omega(t)$.

- **Aperiodic Oscillations:** The current $I_1(t)$ oscillates around a positive mean value with an amplitude that systematically grows over time.
    
- **Sudden Reversals:** Once the amplitude reaches a critical threshold, the current suddenly crosses the zero axis and begins oscillating around a negative mean value, again with growing amplitude, until it flips back.
    
- **Conclusion on Periodicity:** The intervals between these reversals are irregular, and the exact waveforms do not repeat. Therefore, the system is **not periodic**. It exhibits unpredictable, aperiodic behavior characteristic of deterministic chaos.

The term "size" in this context refers to the **dimension of the phase space** (the number of independent variables). By using the conservation law $\Omega_1 - \Omega_2 = 2\Omega_0$, we reduced the 4D Rikitake system to a **3D system** governed by $(I_1, I_2, \Omega)$.

Observing the time series, the system behaves chaotically: it oscillates aperiodically and undergoes unpredictable reversals. This chaotic behavior is perfectly compatible with its 3D size due to the **Poincaré-Bendixson theorem**. This foundational theorem in dynamical systems states that continuous, autonomous, deterministic chaos cannot exist in 1D or 2D systems. In 2D space, trajectories cannot cross each other without intersecting, meaning any bounded motion must eventually settle into a fixed point or a simple closed loop (a periodic limit cycle).

To achieve chaos, trajectories must be able to continuously stretch and fold without ever intersecting themselves or locking into a repeating loop. This requires a minimum of **three dimensions** to allow trajectories to pass "over" and "under" one another, forming a strange attractor. Because the Rikitake system is exactly 3D, it meets the strict minimum dimensional requirement to support the chaotic reversals observed in the time series.


### Part 2: Time Series Analysis via the 1D Map (Incorporating Figure 5)

The provided plot in Figure 5 visually maps the amplitude of successive angular velocity maxima, $\Omega_{k+1} = f(\Omega_k)$. The black curve represents the map $f(x)$, and the red line represents the identity line $y = x$.

#### 1. Stability of the Fixed Point

A fixed point occurs where the next maximum equals the current one ($\Omega_{k+1} = \Omega_k$), which is visually represented by the intersection of the black curve and the red line $y=x$. Looking at Figure 5, there is a single intersection point near $\Omega_k \approx 3.75$. At this exact point, the black curve is sloping downward, and its steepness is visibly much greater than the slope of the red line (which is $+1$). Mathematically, this confirms $\vert{}f'(\Omega_k)\vert{} > 1$ at the fixed point. Because the absolute derivative exceeds 1, any small perturbation will be amplified with each iteration. Therefore, the fixed point is **strictly unstable**.

#### 2. Instability of Limit Cycles and Periodicity

Figure 5 visually supports the assumption that the curve is incredibly steep almost everywhere, specifically $\vert{}f'(x)\vert{} > 1$ for all $x$ except perhaps exactly at the asymptotic peak (which trajectories never perfectly hit). If we assume a limit cycle of length $p$ exists—a repeating sequence $(x_1, \dots, x_p)$ where $x_{k+1} = f(x_k)$ and $x_{p+1} = x_1$—we must evaluate the stability of the $p$-th iterate map, $f^p(x)$. Using the chain rule, the derivative along the entire cycle is the product of the local derivatives:

$$\vert{}(f^p)'(x_1)\vert{} = \vert{}f'(x_p)\vert{} \cdot \vert{}f'(x_{p-1})\vert{} \cdots \vert{}f'(x_1)\vert{}$$

Because every individual $\vert{}f'(x_i)\vert{} > 1$, multiplying them together guarantees a product strictly greater than 1:

$$\vert{}(f^p)'(x_1)\vert{} > 1$$

This proves that any theoretical limit cycle is **unstable**. Since there are no stable fixed points and no stable limit cycles to trap the dynamics, the trajectory $(I_1(t), I_2(t), \Omega(t))$ **cannot be periodic**. It will never perfectly repeat its past behavior.

#### 3. Lyapunov Exponent and Trajectory Nature

The Lyapunov exponent $\lambda$ measures the average exponential rate at which two infinitesimally close trajectories separate over time. For a 1D map, it is calculated as:

$$\lambda = \lim_{n \to \infty} \frac{1}{n} \sum_{k=0}^{n-1} \ln \vert{}f'(\Omega_k)\vert{}$$

Because $\vert{}f'(x)\vert{} > 1$ universally across the map, it follows that $\ln \vert{}f'(x)\vert{} > 0$ for every single point evaluated along the trajectory. Averaging strictly positive numbers guarantees a **strictly positive Lyapunov exponent** ($\lambda > 0$). A positive Lyapunov exponent is the definitive mathematical signature of **chaos**. It dictates that the trajectory exhibits extreme sensitivity to initial conditions; any microscopic change in the initial state will grow exponentially until the long-term state of the system is completely unpredictable.

#### 4. Compatibility with Earth's Magnetic Field

This model aligns beautifully with the paleomagnetic record of the Earth's magnetic field. Geological evidence shows that Earth's field does not reverse like a perfectly tuned clock. Instead, it maintains a stable polarity for highly irregular, unpredictable durations (epochs ranging from tens of thousands to tens of millions of years) before rapidly flipping.

The Rikitake dynamo captures this essential physics: deterministic electro-mechanical equations naturally produce a chaotic strange attractor, resulting in non-periodic, structurally unpredictable magnetic reversals driven purely by internal fluid dynamics, without the need for external triggers.


![[Pasted image 20261006000554.png]]
# 22. — Epidemiology

This problem explores the mathematical foundations of epidemic modeling using a continuous-time renewal equation, often applied to model the spread of airborne viruses like SARS-CoV-2. The model bridges individual viral shedding kinetics with macroscopic population transmission dynamics.

**1. Deriving the Integral Equation for the Infection Rate $I(t)$** The infection rate $I(t)$ represents the number of newly infected individuals per unit time. To find the number of new infections at current time $t$, we must look at all individuals infected in the past and determine how much they are transmitting the virus now.

- Consider an index case infected at a past time $t - t'$. The time elapsed since their infection is $t'$.
    
- Their rate of emitting virus at this current time $t$ is proportional to the transmissibility function $\psi(t')$.
    
- The total number of people infected at that specific past time was $I(t - t')$.
    
- The parameter $R$ is the epidemic reproduction rate (total expected secondary cases per index case in a fully susceptible population). To ensure the total integral of the transmission kernel over an individual's entire infectious period equals $R$, we normalize the transmissibility profile. Since $\int_0^\infty \psi(t') dt' = T$, the transmission rate generated by one individual at infection age $t'$ is $\frac{R}{T}\psi(t')$.
    
- Integrating this over all possible past infection times $t'$ (from $0$ to $\infty$) gives the total "infectious pressure" currently exerted on the population.
    
- However, not everyone is susceptible. Only a fraction $(1 - P(t))$ of the population can actually catch the disease.
    

Multiplying the total infectious pressure by the susceptible fraction yields the required renewal equation:

$$I(t) = (1 - P(t)) \frac{R}{T} \int_0^\infty I(t - t') \psi(t') dt'$$

**2. Expressing Epidemic Prevalence $P(t)$** Assuming each infection leads to long-term immunization, the epidemic prevalence $P(t)$ is the fraction of the total population $N$ that has ever been infected up to time $t$. To find the total number of people infected, we integrate the infection rate $I(t')$ from the beginning of the epidemic ($t \to -\infty$) up to the current time $t$. Dividing this by the total population $N$ gives the fraction:

$$P(t) = \frac{1}{N} \int_{-\infty}^t I(t') dt'$$

**3. The Exponential Solution for a Fully Susceptible Population** Early in an epidemic, a very small fraction of the population is immune, so $P(t) \approx 0$. The integral equation simplifies to a linear Volterra integral equation:

$$I(t) = \frac{R}{T} \int_0^\infty I(t - t') \psi(t') dt'$$

We test the exponential ansatz $I(t) = I_0 e^{\sigma t}$ (where $\sigma$ is the exponential growth rate):

$$I_0 e^{\sigma t} = \frac{R}{T} \int_0^\infty I_0 e^{\sigma (t - t')} \psi(t') dt'$$

Factor out $I_0 e^{\sigma t}$ on the right side:

$$I_0 e^{\sigma t} = I_0 e^{\sigma t} \frac{R}{T} \int_0^\infty e^{-\sigma t'} \psi(t') dt'$$

Divide both sides by $I_0 e^{\sigma t}$:

$$1 = \frac{R}{T} \int_0^\infty e^{-\sigma t'} \psi(t') dt'$$

Rearranging this to solve for the reproduction rate $R$ yields the classic Euler-Lotka equation for epidemiology:

$$R = \frac{T}{\int_0^\infty e^{-\sigma t'} \psi(t') dt'}$$

This shows that $R$ is inversely proportional to the Laplace transform of the transmissibility profile evaluated at the growth rate $\sigma$.

**4. Expanding to a Multi-Class Society** When society is divided into $J$ classes, the transmission dynamics become heterogeneous. We deconstruct the given Equation (30) logically:

- **Total Infectious Pressure:** An individual in class $k$ has a relative onward transmissibility $\mathscr{T}_k$. The total infectious pressure generated by all past cases across all classes is the weighted sum: $\int_0^\infty \sum_k \mathscr{T}_k I_k(t - t') \frac{R}{T} \psi(t') dt'$.
    
- **Contact Distribution:** Assuming homogeneous mixing (mass action), this total infectious pressure is distributed across the population. The fraction of this pressure that targets class $j$ is simply their demographic proportion of the total population: $\frac{N_j}{N}$.
    
- **Class Susceptibility:** Individuals in class $j$ have a relative susceptibility $\mathscr{S}_j$.
    
- **Depletion of Susceptibles:** Only the fraction $(1 - P_j(t))$ of class $j$ remains capable of being infected.
    

Multiplying the total infectious pressure by the demographic fraction, the relative susceptibility, and the remaining susceptible fraction directly produces the coupled integral equation for the infection rate in class $j$:

$$I_j(t) = \mathscr{S}_j (1 - P_j(t)) \frac{R}{T} \frac{N_j}{N} \int_0^\infty \sum_k \mathscr{T}_k I_k(t - t') \psi(t') dt'$$

**5. Multi-Class Exponential Growth Relation** At the onset of the epidemic, all $P_j(t) \approx 0$. We assume that the infection rates in all classes grow exponentially at the same shared rate $\sigma$, such that $I_j(t) = I_{j,0} e^{\sigma t}$. Substitute this into the linearized multi-class equation:

$$I_{j,0} e^{\sigma t} = \mathscr{S}_j \frac{R}{T} \frac{N_j}{N} \int_0^\infty \sum_k \mathscr{T}_k I_{k,0} e^{\sigma (t - t')} \psi(t') dt'$$

Factor out $e^{\sigma t}$ and cancel it from both sides:

$$I_{j,0} = \mathscr{S}_j \frac{N_j}{N} \frac{R}{T} \left( \sum_k \mathscr{T}_k I_{k,0} \right) \int_0^\infty e^{-\sigma t'} \psi(t') dt'$$

To solve this eigenvalue problem, we must eliminate the constant amplitudes $I_{j,0}$. We do this by multiplying both sides of the equation by $\mathscr{T}_j$ and then summing over all classes $j$:

$$\sum_j \mathscr{T}_j I_{j,0} = \sum_j \left( \mathscr{T}_j \mathscr{S}_j \frac{N_j}{N} \right) \frac{R}{T} \left( \sum_k \mathscr{T}_k I_{k,0} \right) \int_0^\infty e^{-\sigma t'} \psi(t') dt'$$

Notice that the term $\left(\sum_k \mathscr{T}_k I_{k,0}\right)$ on the right side is identical to the term $\left(\sum_j \mathscr{T}_j I_{j,0}\right)$ on the left side. Assuming there is an active epidemic, this sum is non-zero, allowing us to divide it out completely:

$$1 = \left( \sum_j \mathscr{T}_j \mathscr{S}_j \frac{N_j}{N} \right) \frac{R}{T} \int_0^\infty e^{-\sigma t'} \psi(t') dt'$$

Solving for $R$, we find the generalized relationship linking the reproduction rate to the exponential growth rate and the heterogeneous population traits:

$$R = \frac{T}{\left( \sum_j \mathscr{T}_j \mathscr{S}_j \frac{N_j}{N} \right) \int_0^\infty e^{-\sigma t'} \psi(t') dt'}$$

This result beautifully demonstrates that in a heterogeneous mixing model, the effective reproduction rate is scaled by the population-weighted average of the product of relative susceptibility and relative transmissibility
.

