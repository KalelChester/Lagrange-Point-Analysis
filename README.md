# PH 621 (Classical Mechancis) Computational Final Project
This project contains an exploration of the restricted three body probelm wherein the motion of a small "massless test particle" is explored. Specifically, this project is focused on an exploration of the Lagrange points. 

The Lagrange points are points of equillibrium for small mass object under the graviational influence of two larger bodies. Under most conditions, the restricted three body problem is a chaotic system, but given particular starting conditions, it is possible for the smaller body to remain in an equillibrium point relative to the position of the two larger bodies.

## Project Overview
* Derive the position of the Lagrange points
* Asses the stability of these points
* Simulate dynamics of a small mass placed in each Lagrange point affected by a small perturbation

## Identification of the Lagrange Points
The procedure for finding the Lagrange points is fairly straightforward. The following will replicate the derivaiton of the points as conducted by Neil J. Cornish in [^1].

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
U = -G(\frac{M_1}{r_1}+\frac{M_2}{r_2})
.$$

The velocity dependant terms in the potential do not effect the positions of the equillibrium points, but are key in understanding the stability of motion about the equillibrium points.

The potential $U_\Omega$ is plotted with $M_1=10$, $M_2=1$, $R=10$, $G=1$, and $\vec{v}=0$ in Figure 1 with numerically calculated Lagrange points labled.

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

The force corresponding to the potential is plotted in Figure 2.

![Heatmap of F_eff](plots/feff.png) <br>
**Figure 2**

Since the system is symmetric about the $x$-axis, the $y$ component of the force must vanish. To find the first three Lagrange points, we can set $y=0$ and use the change of coordinate $x=R(u+\beta)$. The desired 3 coordinates are the solutions to the three fifth-order equations

$$
u^2((1-s_1) + 3u + 3u^2 +u^3) = \alpha(s_0+2s_0 u + (1+ s_0 - s_1)u^2 +2u^3 +u^4)
$$

where $s_0=\text{sign}(u)$ and $s_1=\text{sign}(u+1)$. The cases that give L1, L2, and L3 have $(s_0, s_1)$ equal to (-1, 1), (1, 1), and (-1, -1) respectively. Each equation has only one real root to the quintic equation, which gives the Lagrange points. However, we are unable to find a general solution to these equations for arbitrary $\alpha$ so we take the approximation of $\alpha \ll 1$ which applies for relevant equations such as the Earth-Sun system $(\alpha \approx 10^{-6})$ and (less strongly) the Earth-Moon system $(\alpha \approx 10^{-2})$. To first order in $\alpha$ we find these first three Lagrange points to be positioned at:

$$
L1: \begin{pmatrix} R[1-(\frac{\alpha}{3})^{1/3}] \\ 0 \end{pmatrix}
$$
$$
L2: \begin{pmatrix} R[1+(\frac{\alpha}{3})^{1/3}] \\ 0 \end{pmatrix}
$$
$$
L3: \begin{pmatrix} -R[1+\frac{5}{12}\alpha] \\ 0 \end{pmatrix}
$$

The accuracy of these approximations depends on the mass ratio $\frac{\beta}{\alpha}$. The relative magnitude of the error is shown in Figure 3,

![Mass ratio accuracy](plots/mass_ratio_accuracy.png) <br>
**Figure 3**

To find the fourth and fifth Lagrange points, we must take a different approach. These points are found where the centrifugal force is balanced with the gravitional force from each of the bodies. **CONTINUE WITH L4 and L5**

## Stability of the Lagrange Points
TODO

## Dynamics Near the Lagrange Points
TODO


Sources: <br>
[^1] https://web.archive.org/web/20150907090116/http://www.physics.montana.edu/faculty/cornish/lagrange.pdf

Simulations in this paper made use of the REBOUND N-body code \citep{rebound}. The simulations were integrated using IAS15, a 15th order Gauss-Radau integrator \citep{reboundias15}. 


@ARTICLE{rebound,
       author = {{Rein}, H. and {Liu}, S. -F.},
        title = "{REBOUND: an open-source multi-purpose N-body code for collisional dynamics}",
      journal = {\aap},
     keywords = {methods: numerical, planets and satellites: rings, protoplanetary disks, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Mathematics - Dynamical Systems, Physics - Computational Physics},
         year = 2012,
        month = jan,
       volume = {537},
          eid = {A128},
        pages = {A128},
          doi = {10.1051/0004-6361/201118085},
archivePrefix = {arXiv},
       eprint = {1110.4876},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2012A&A...537A.128R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
@ARTICLE{reboundias15,
       author = {{Rein}, Hanno and {Spiegel}, David S.},
        title = "{IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits}",
      journal = {\mnras},
     keywords = {gravitation, methods: numerical, planets and satellites: dynamical evolution and stability, Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics, Mathematics - Numerical Analysis},
         year = 2015,
        month = jan,
       volume = {446},
       number = {2},
        pages = {1424-1437},
          doi = {10.1093/mnras/stu2164},
archivePrefix = {arXiv},
       eprint = {1409.4779},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.1424R},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}