# Lagrange Point Analysis
This project contains an exploration of the restricted three body probelm wherein the motion of a small "massless test particle" is explored. Specifically, this project is focused on an exploration of the Lagrange points. 

The Lagrange points are points of equillibrium for small mass object under the graviational influence of two larger bodies. Under most conditions, the restricted three body problem is a chaotic system, but given particular starting conditions, it is possible for the smaller body to remain in an equillibrium point relative to the position of the two larger bodies.

## Project Overview
* Derive the position of the Lagrange points
* Asses the stability of these points
* Simulate dynamics of a small mass placed in each Lagrange point affected by a small perturbation

## Identification of the Lagrange Points
The procedure for finding the Lagrange points is fairly straightforward. The following will replicate the derivaiton of the points as conducted by Neil J. Cornish in [1].

We seek solutions which maintain a fixed distance between the three bodies. If $M_1$ and $M_2$ are the two larger masses, and $\vec{r_1}$ and $\vec{r_2}$ are the position vectors of the respective masses, then the total force exerted on a third mass $m$ at a position $\vec{r}$ will be

$$ 
\vec{F}
= -\frac{GM_1m}{|\vec{r}-\vec{r_1}|^3} (\vec{r}-\vec{r_1})
-\frac{GM_2m}{|\vec{r}-\vec{r_2}|^3} (\vec{r}-\vec{r_2})
$$

To find stationary solutions, we adopt a co-rotating frame in which $M_1$ and $M_2$ are held fixed. This frame has the origin fixed at the center of mass of the two larger bodies and has an angular frequency $\Omega$ given by Kepler's law:

$$
\Omega^2 R^3 = G(M_1 + M_2)
$$

The effective potential on $m$ is

$$
U_\Omega = U - 
\vec{v} \cdot (\vec{\Omega} \times \vec{r}) -
\frac{1}{2}(\vec{\Omega} \times \vec{r})\cdot(\vec{\Omega} \times \vec{r})
$$

With $U$ being the gravitational potential ($m$ can be set to unity without a loss of generality), 

$$
U = -G \left( \frac{M_1}{r_1}+\frac{M_2}{r_2} \right)
.$$

The velocity dependant terms in the potential do not effect the positions of the equillibrium points, but are key in understanding the stability of motion about the equillibrium points.

The potential $U_\Omega$ is plotted with $M_1=10$, $M_2=1$, $R=10$, $G=1$, and $\vec{v}=0$ in figure 1 with numerically calculated Lagrange points labled.

![Heatmap of U_eff](plots/ueff.png) <br>
**Figure 1**

Taking the negative gradint of this potential, we get the force,

$$
\vec{F} = \Omega^2 \begin{pmatrix}
x-\frac{\beta(x+\alpha R)R^3}{((x+\alpha R)^2 + y^2)^{3/2}} - \frac{\alpha(x-\beta R)R^3}{((x-\beta R)^2 + y^2)^{3/2}} \\
y-\frac{\beta y R^3}{((x+\alpha R)^2 + y^2)^{3/2}} - \frac{\alpha y R^3}{((x-\beta R)^2 + y^2)^{3/2}}
\end{pmatrix}
$$

Where

$$ \alpha = \frac{M2}{M1+M2} $$
$$ \beta = \frac{M1}{M1+M2} $$

The force corresponding to the potential is plotted in figure 2.

![Heatmap of F_eff](plots/feff.png) <br>
**Figure 2**

Since the system is symmetric about the $x$-axis, the $y$ component of the force must vanish. To find the first three Lagrange points, we can set $y=0$ and use the change of coordinate $x=R(u+\beta)$. The desired 3 coordinates are the solutions to the three fifth-order equations

$$
u^2((1-s_1) + 3u + 3u^2 +u^3) = \alpha(s_0+2s_0 u + (1+ s_0 - s_1)u^2 +2u^3 +u^4)
$$

where $s_0=\text{sign}(u)$ and $s_1=\text{sign}(u+1)$. The cases that give $L_1$, $L_2$, and $L_3$ have $(s_0, s_1)$ equal to (-1, 1), (1, 1), and (-1, -1) respectively. Each equation has only one real root to the quintic equation, which gives the Lagrange points. However, we are unable to find a general solution to these equations for arbitrary $\alpha$ so we take the approximation of $\alpha \ll 1$ which applies for relevant equations such as the Earth-Sun system $(\alpha \approx 10^{-6})$ and (less strongly) the Earth-Moon system $(\alpha \approx 10^{-2})$. To first order in $\alpha$ we find these first three Lagrange points to be positioned at:

$$
L_1: \begin{pmatrix} R[1-(\frac{\alpha}{3})^{1/3}] \\ 0 \end{pmatrix}
$$
$$
L_2: \begin{pmatrix} R[1+(\frac{\alpha}{3})^{1/3}] \\ 0 \end{pmatrix}
$$
$$
L_3: \begin{pmatrix} -R[1+\frac{5}{12}\alpha] \\ 0 \end{pmatrix}
$$

The accuracy of these approximations depends on the mass ratio $\frac{\beta}{\alpha}$. The relative magnitude of the error is shown in figure 3,

![Mass ratio accuracy](plots/mass_ratio_accuracy.png) <br>
**Figure 3**

To find the fourth and fifth Lagrange points, we must take a different approach. These points are found where the centrifugal force is balanced with the gravitional force from each of the bodies. To do this, we can resolve the force into parallel and perpecdicular components relative to $\vec{r}$. The perpendicular component yields

$$
F_{\Omega}^{\perp}=\alpha \beta y \Omega^2 R^3 \left(\frac{1}{((x-\beta R)^2 + y^2)^{3/2}}-\frac{1}{((x+\alpha R)^2 + y^2)^{3/2}}\right).
$$

Setting this to 0 for non-zero $y$ tells us that the distance between the test points and each mass must be equal. Using this fact, the parrallel component of the force becomes

$$
F_{Ω}^\parallel=\Omega^2\frac{x^2+y^2}{R}\left(\frac{1}{R^3}-\frac{1}{((x-\beta R)^2 + y^2)^{3/2}}\right),
$$

which, when set equal to 0 gives that the distance that the test points must be a distance R from each mass. This means that the location of $L_4$ and $L_5$ are given by the tip of an equillateral triangle of side length $R$ with the other two points being centered on the COM of the larger masses.

$$
L_4: \begin{pmatrix} \left(\frac{R}{2} (\beta-\alpha)\right) \\ \frac{\sqrt{3}}{2}R \end{pmatrix}
$$
$$
L_5: \begin{pmatrix} \left(\frac{R}{2} (\beta-\alpha)\right) \\ -\frac{\sqrt{3}}{2}R \end{pmatrix}
$$

Note that $L_4$ and $L_5$ are exact, not subject to the same approximation used to find the first three Lagrange points.

## Stability of the Lagrange Points
Now that we have found the 5 Lagrange points, the next step is to determine the stability of these points. While It is usually sufficient to look at the shape of the potential to find the stability of an equillibrium point, the fact that our potential is velocity dependant makes this methodology insufficient. Instead we must perform a linear stability analysis about each Lagrange point by linearizing the eqution of motion about each equillibrium point and solving for small perturbations. Writing

$$ 
x=x_i+\delta x, \hspace{.5cm}
v_x = \delta v_x,
$$
$$ 
y=y_i+\delta y, \hspace{.5cm}
v_y = \delta v_y,
$$

Where $(x_i, y_i)$ is the position of the $i$-th Lagrange point, and $v_j$ is the $j$-th spatial velocity component of the test particle in the rotating frame $\Omega$. With this, the linearized equatinos of motion become 

$$
\frac{d}{dt} 
\begin{pmatrix} \delta x \\ \delta y \\ \delta v_x \\ \delta v_y \end{pmatrix}
= \begin{pmatrix}
0 & 0 & 1 & 0 \\
0 & 0 & 0 & 1 \\
\frac{d^2 U_\Omega}{dx^2} & \frac{d^2 U_\Omega}{dxdy} & 0 & 2 \Omega \\
\frac{d^2 U_\Omega}{dydx} & \frac{d^2 U_\Omega}{dy^2} & -2 \Omega & 0 \\
\end{pmatrix}
$$

In the following sections, the second derivatives will be evaluated at $\vec{r} = (x_i, y_i)$. 

The general soluition for the perturbation $\xi (\delta x, \delta y, \delta v_x, \delta v_y )^T$ is a linear combination of eigenmodes:

$$
\xi(t) = \sum^4_{k=1}c_ke^{\lambda_kt}.
$$

Writing each eigenvalue as $\lambda = \gamma + i\omega$, the behavior is controlled by its real part $\gamma=\operatorname{Re}(\lambda)$. Because of this, any positive real component of $\lambda$ will result in unstable equillibrium points.

To measure the instability of unstable equillibrium points we can calculate the $e$-folding time, $\tau$ which corresponds to the characteristic time required for a small intial displacement to grow by a factor of $e$.

$$
\tau=\frac{1}{\gamma}, \hspace{.5cm}\text{for } \gamma>0
$$

### Stability of $L_1$ and $L_2$

The curvature of the effective potential near $L_1$ and $L_2$ shows that they are saddle points:

$$
\frac{d^2U_\Omega}{dx^2}=-9\Omega \hspace{5mm}
\frac{d^2U_\Omega}{dy^2}=3\Omega
$$
$$
\frac{d^2U_\Omega}{dxdy}=\frac{d^2U_\Omega}{dydx}=0
$$

Calculating the eigenvalues of this system gives

$$
\lambda=\pm \Omega \sqrt{1+2\sqrt{7}},\hspace{2mm} \pm i \Omega \sqrt{2\sqrt{7}-1}.
$$

The existance of a real, positive eigenvalue means that these equillibrium points are unstable. The $e$-folding time for these points is about

$$
\tau \approx \frac{2}{5\Omega}=23 \text{ days.}
$$

This means that a satellite placed exactly at one of these points will drift away within a few months if not corrected.

However, there are still many projects that use these points by exploiting the purely imaginary eigenvalues of the system. While the positive real eigenvalue drives exponential divergence, the purely imaginary pair ($\pm i\Omega\sqrt{2\sqrt{7}-1}$) corresponds to bounded, oscillatory motion that, in the full three-dimensional system, gives rise to periodic and quasi-periodic trajectories known as halo and Lissajous orbits. Missions such as the James Webb Space Telescope (JWST) at $L_2$ [2] and the Solar and Heliospheric Observatory (SOHO) at $L_1$ [3] are injected directly into these halo orbits so that their initial displacement along the unstable eigenmode is nearly zero. Because real-world perturbations—such as solar radiation pressure and gravitational tugs from other planets—inevitably reintroduce a small component along the unstable direction, these spacecraft execute small, periodic station-keeping maneuvers (every few weeks for JWST [2] and roughly every 90 days for SOHO [3]) to zero out the growing unstable mode before it can pull them off course, allowing them to maintain their orbits for years with minimal propellant.

### Stability of $L_3$

$L_3$, like it's colinear counterparts, is also located at a addle point of the effective potential with curvature

$$
\frac{d^2U_\Omega}{dx^2}=-\frac{1}{4}(12+7\alpha)\Omega^2 \hspace{5mm}
\frac{d^2U_\Omega}{dy^2}=\frac{7 \alpha \Omega^2}{8}
$$
$$
\frac{d^2U_\Omega}{dxdy}=\frac{d^2U_\Omega}{dydx}=0
$$

The linearized evolution matrix eigenvalues to first order in $\alpha$ are

$$
\lambda = \pm \frac{i}{8}(8+7\alpha)\Omega, \hspace{2mm} \pm \frac{1}{2}\sqrt{\frac{21\alpha}{2}}\Omega
$$

Once again, the existence of a real, positive eigenvalue means that this point is unstable. The $e$-folding time of this point is 

$$
\tau \approx 56.7 \text{ years.}
$$

### Stability of L4 and L5

The curvature of effective potential at L4 and L5 

$$
\frac{d^2U_\Omega}{dx^2}=-\frac{3}{4}\Omega^2 \hspace{5mm}
\frac{d^2U_\Omega}{dy^2}=-\frac{9}{4}\Omega^2
$$
$$
\frac{d^2U_\Omega}{dxdy}=\frac{d^2U_\Omega}{dydx}=\frac{3\sqrt{3}}{4}(1-2\alpha)\Omega^2=\frac{3\sqrt{3}}{4}\kappa\Omega^2
$$

where $\kappa=\frac{M_1-M_2}{M_1+M_2}$. The eigenvalues of the linearized evolution matrix are

$$
\lambda = \pm i\frac{\Omega}{2}\sqrt{2\pm\sqrt{27\kappa^2-23}}
$$

Since the points will be stable if the eigenvalues are purely imaginary, we get the stability conditions:

$$
\kappa^2 > \frac{23}{27}, \text{  and  } 
\sqrt{27\kappa^2-23} \leq2
$$

The second condition is always satisfied since $\kappa < 1$, while the first requires that 

$$
M_1 \geq \frac{1}{2}M_2(25+3\sqrt{69})
$$

When this stability condition is met, the $L_4$ and $L_5$ points yeild stable orbits and are refered to as Trojan points. The mass ratios for both the Earth-Moon and Sun-Earth system satisfy these conditions.

### Numeric Stability

Using the REBOUND N-body code [4] (integrated using IAS15, a 15th order Gauss-Radu integrator [5]) and the same parameters used to generate figures 1-3, we can simulate what happends if we place a test mass at each of the Lagrange points, and enact a small randomized perturbation $<1e-4$ on their x and y positions.

Figures 4-8 are example gifs showing motion of test particles located at each of the Lagrange points after a small perturbation in their position for a mass ratio $\frac{M_1}{M_2}=10$

![Perturbed L1 Dynamics](plots/l1.gif) <br>
**Figure 4**

![Perturbed L2 Dynamics](plots/l2.gif) <br>
**Figure 5**

![Perturbed L3 Dynamics](plots/l3.gif) <br>
**Figure 6**

![Perturbed L4 Dynamics](plots/l4.gif) <br>
**Figure 7**

![Perturbed L5 Dynamics](plots/l5.gif) <br>
**Figure 8**

We can test the stability criteria found for $L_4$ and $L_5$ by running additional simulations with a higher mass ratio of $\frac{M_1}{M_2}=25$, which just barely clears the stability requirement.

![Perturbed L4 Dynamics](plots/l4_heavy.gif) <br>
**Figure 9**

![Perturbed L5 Dynamics](plots/l5_heavy.gif) <br>
**Figure 10**

As we can see, the particles remain stably in place for the duration of the simulation which is equivalent to more than 50,000 orbits of the larger bodies (years).

## Sources:
[1]: N.J Cornish, *The Lagrange Points*, https://science.nasa.gov/wp-content/uploads/2023/07/3322_lagrange.pdf

[2]: J. Petersen, *L2 Station Keeping Maneuver Strategy for the James Webb Space Telescope* Preprint, https://ntrs.nasa.gov/api/citations/20190028877/downloads/20190028877.pdf

[3]: C.E Roberts, *Long Term Missions at the Sun-Earth Liberation Point L1: ACE, SOHO, and WIND*, AAS 11-485,https://ntrs.nasa.gov/api/citations/20110022484/downloads/20110022484.pdf 

[4]: H. Rein and S.-F. Liu, *REBOUND: an open-source multi-purpose N-body code for collisional dynamics* A&A, 537 (2012) A128, https://doi.org/10.1051/0004-6361/201118085

[5]: H. Rein, D.S. Spiegel, *ias15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits*, Monthly Notices of the Royal Astronomical Society, Volume 446, Issue 2, 11 January 2015, Pages 1424–1437, https://doi.org/10.1093/mnras/stu2164

[6]: Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.

[7]: P. Virtanen, et al., (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272. 10.1038/s41592-019-0686-2.

[8]: J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007

[9]: Waskom, M. L., (2021). seaborn: statistical data visualization. Journal of Open Source Software, 6(60), 3021, https://doi.org/10.21105/joss.03021