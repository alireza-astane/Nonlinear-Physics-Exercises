
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
    
2. Global Mass Conservation (Eq 2): The total cross-sectional area scales linearly with time:
    
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

