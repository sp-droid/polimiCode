





<center><big><big><b>
Politecnico di Milano
<br><br>Master in Space Engineering
<br>Academic Year 2022 / 23<br><br>Orbital mechanics course
</b></big></big></center>







<center><big><big><b>
Final assignment
<br><br>Interplanetary and planetary
<br>explorer missions
</b></big></big></center>








<right><big><b>
GROUP NUMBER: 2204
<br><br>name (name id)
<br><br>name (name id)
<br><br>name (name id)
<br><br>name (name id)
</b></big></right>




<center><big><big><b>January 2023</b></big></big></center>
<div style="page-break-after: always; break-after: page;"></div>

# Table of Contents

[TOC]

<div style="page-break-after: always; break-after: page;"></div>

---

# Assignment 1: Interplanetary Explorer

>The PoliMi Space Agency is carrying out a feasibility study for a potential Interplanetary Explorer Mission visiting an asteroid in the Solar System, with an intermediate flyby on a planet.
>
>As part of the mission analysis team, you are requested to perform the preliminary mission analysis. You have to study the transfer options from Earth to the arrival asteroid, with a powered gravity assist (flyby) at the intermediate planet, and propose a solution based on the mission cost (measured through the total Δ𝑣).
>
>The mission departs from Earth, and the flyby planet and arrival asteroid have been decided by the science team. Constraints on earliest departure and latest arrival have also been set by the launch provider, the systems engineering team, and the Agency’s leadership.



<div style="page-break-after: always; break-after: page;"></div>

---

# Assignment 2: Planetary Explorer

>The PoliMi Space Agency wants to launch a Planetary Explorer Mission, to perform Earth observation.
>
>As part of the mission analysis team, you are requested to carry out the orbit analysis and ground track estimation. You have to study the effects of orbit perturbations, and compare different propagation methods. Also, you have to characterize the ground track, and propose an orbit modification to reach a repeating ground track (for better communications with our network of ground stations).



## 2.1 Nominal orbit summary



### 2.1.1 Assigned* parameters

| **Semi-major axis**  | **Eccentricity** | **Inclination** | **RAAN** | **Argument of periapsis** | **True anomaly at t=0** |
| -------- | -------- | -------- | -------- | -------- | -------- |
| <u>26,619 km</u> | <u>0.7452</u> | <u>62.9089º</u> | 60º | 30º | 0º |

**Initial date**: November 1st, 2021 at 0h 0' 0'' GMT

*****RAAN, argument of periapsis and initial true anomaly were left to our criteria, which were chosen in order to have a good movie and after fixing the initial date.



### 2.1.2 Assigned perturbations

- **J<sub>2</sub> effect** due to Earth's oblateness

$$
\Large
\mathbf{a_{J_2}} =   \frac{3}{2}\frac{J_2 \mu R_E^2}{r^5}\left [ \left ( 5\frac{z^2}{r^2}-1 \right )x\widehat{\mathbf{i}} + \left ( 5\frac{z^2}{r^2}-1 \right )y\widehat{\mathbf{j}} + \left ( 5\frac{z^2}{r^2}-3 \right )z\widehat{\mathbf{k}}\right ]
$$

<u>**J<sub>2</sub>** =  0.0010826</u> as provided by the function *astroConstants* 


- **Atmospheric drag** due to the absorption and reflection of massive particles against the spacecraft

$$
\Large
\mathbf{\ddot{r}} = -\frac{1}{2}C_D\frac{A}{m}\rho (h,t)v_{rel}\mathbf{v_{rel}}
$$

$$
\Large
where: \hspace{4mm} \mathbf{v_{rel}} = \mathbf{v}-\mathbf{\omega_E\wedge} \mathbf{r}
$$
The density at very high altitudes significantly depends on day/night cycles and the solar activity at any given time, which by itself makes the modelling of aerodynamic drag a challenge. There are several models, like the Jacchia-Bowman 2008, that adress this issue with varying levels of success. However, we are going to employ a simple exponential approximation that depends only on the altitude:
$$
\Large
\rho = \rho_0 exp\left ( -\frac{h-h_0}{H_0} \right ) \hspace{2mm};\hspace{2mm} \rho_0 = 2\cdot 10^{-8} \hspace{2mm};\hspace{2mm} h_0 = 122 km \hspace{2mm};\hspace{2mm} H_0 = H_0|(\rho=3\cdot 10^{-12},h=400km)
$$

<u>**C<sub>D</sub>** = 2.1</u>

<u>**A/m** = 0.0095 m<sup>2</sup>/kg</u>

<div style="page-break-after: always; break-after: page;"></div>

### 2.1.3 Other perturbations

In the last point, the comparison with real data, we include an additional third model that features the following perturbations:

- **Atmospheric drag** (Explained earlier)
- **Solar radiation pressure**

$$
\Large
\mathbf{\ddot{r}} = -\upsilon \phi _EC_R\frac{A}{m}\mathbf{\widehat{r}} \hspace{2mm};\hspace{2mm} \phi _E = \phi _0\left ( \frac{R{sun}}{r_{sun}} \right )^2 \hspace{2mm};\hspace{2mm} \phi _0 = \sigma T_{sun}^4 
$$

$$
\Large
where: \hspace{4mm} \sigma = Stefan-Boltzmann\hspace{2mm}constant, \phi = Solar\hspace{1mm}flux
$$

$$
\Large
and: \hspace{4mm} \upsilon = shadow function(r,r_{sun},R_E,R_{sun})
$$

This shadow function is implemented and working and can be seen in the referenced book.

- **Earth gravitational model 1996** (And therefore substitutes the J<sub>2</sub>-only perturbation)

$$
\Large
U(r,\phi,\xi) = -\sum_{n=2}^{360}\frac{\mu R_E^n}{r^{n+1}}\sum_{m=0}^{n}\left ( C_{nm}cos(m\phi)+S_{nm}sin(m\phi) \right )\widetilde{P}_{nm}\left ( \xi  \right )
$$

$$
\Large
where: \hspace{4mm} r^2 = x^2+y^2+z^2\hspace{1mm};\hspace{1mm} \phi = tan^{-1}\left ( \frac{y}{x} \right )\hspace{1mm};\hspace{1mm}\xi = \frac{z}{r}\hspace{1mm};\hspace{1mm}\widetilde{P} = norm.\hspace{1mm}Legendre \hspace{1mm}polynom.
$$

Deriving -U(x,y,z) gives the perturbing accelerations we are seeking.

With this it's also possible to calculate the undulation potential, which shows Earth's uneven gravity field:

<figure>
  <div align="center"><img src="assets\undulation.png" alt="Undulation">
  <figcaption>Fig.1 - Undulation potential T(lat,lon) </figcaption>
</figure>

- **Sun and Moon 3rd body perturbations**

$$
\Large
\mathbf{\ddot{r}} = GM_2\left ( \frac{\mathbf{r_2-r}}{\left \| r_2-r \right \|^3}-\frac{\mathbf{r_2}}{r_2^3} \right )
$$

Where **r** and **r<sub>2</sub>** are the position vectors of the spacecraft and the planet of mass **M<sub>2</sub>** relative to the Earth, **2** being the Sun or the Moon.

- **Relativistic effects**

$$
\Large
\mathbf{\ddot{r}} = -\frac{GM_E}{r}\left [ (4\frac{GM}{c^2r}-\frac{v^2}{c^2})\mathbf{r}+\frac{4}{c^2}(\mathbf{r}\cdot \mathbf{v})\mathbf{v} \right ]
$$



## 2.2 Groundtrack

As a spacecraft orbits around a celestial body, the path drawn on it is deviated from an ellipse by the rotation of the body and the perturbations affecting the spacecraft.

### 2.2.1 Orbit GT

<figure>
  <div align="center"><img src="assets\gt1.png" alt="Undulation">
  <figcaption>Fig.1 - 1 orbit groundtrack (perturbed-right)</figcaption>
</figure>

<figure>
  <div align="center"><img src="assets\gt1d.png" alt="Undulation">
  <figcaption>Fig.1 - 1 day groundtrack (perturbed-right)</figcaption>
</figure>

<figure>
  <div align="center"><img src="assets\gt10d.png" alt="Undulation">
  <figcaption>Fig.1 - 10 days groundtrack (perturbed-right)</figcaption>
</figure>

### 2.2.2 Repeating GT

## 2.3 Orbit propagation

### 2.3.1 Analytical vs. numerical integration

### 2.3.2 Cartesian vs. Gauss's equations methods

#### 2.3.2.1 High frequency filtering

## 2.4 Evolution of the orbit

### 2.4.1 Short term

We have made two videos of the first 3 orbits. First one shows the orbit from a fixed point relative to Earth's center:

https://www.youtube.com/watch?v=gYX6A3jnGZs

And the second shows the orbit from the point of view of the spacecraft, centered on Earth:

https://www.youtube.com/watch?v=YlP3BiDUItA

### 2.4.2 Long term

## 2.5 Real data comparison

### 2.5.1 Object selection

The aim was to find the satellite or debris with the closest 3 first orbital elements. After a search through some tracking websites, this is the result:

|  | **Semi-major axis** | **Eccentricity** | **Inclination** |
| ------------------- | ---------------- | --------------- | -------- |
| Assigned | 26,619 km           | 0.7452           | 62.9089º        |
| Real | 23,716 km  | 0.712        | 63º     |

The object is designated SL-6 R/B(2), a block-ML from a Molniya rocket that is orbiting the planet since 1999. It is catalogued in NORAD as 25850U and it's described as a 900kg, cylinder shaped (2.6x2.58x2.58m) rocket body, with an average cross section of 7.8824 and radar cross section of 2 sq. meters.

### 2.5.2 Parameter estimation



### 2.5.3 Comparison vs. assigned model



### 2.5.4 Comparison vs. our model

<div style="page-break-after: always; break-after: page;"></div>

---

# Bibliography

[1]: https://www.hindawi.com/journals/tswj/2014/163949 "Translational propagation LHS"
[2]: http://yann.lecun.com/exdb/mnist/ "MNIST dataset"
[3]: https://scholar.google.es/scholar_url?url=https://www.researchgate.net/file.PostFileLoader.html%3Fid%3D56d1ad946307d916ce8b4569%26assetKey%3DAS%253A333743404404737%25401456582036953&hl=en&sa=X&ei=kTnIYqLtO_-Ty9YPo6Gi8A4&scisig=AAGBfm04DGsPU0sKp62oS_g-YQq3UCh_wA&oi=scholarr "Cost sensitive learning"
[4]: https://scholar.google.es/scholar_url?url=https://wires.onlinelibrary.wiley.com/doi/abs/10.1002/widm.1249&hl=en&sa=X&ei=BzrIYtGEBo_ymgG-lo_YBA&scisig=AAGBfm0_vKPlbNkn_FSo7O_I7qO0WOqW0Q&oi=scholarr "Ensemble learning"
[5]: https://proceedings.neurips.cc/paper/2014/file/5ca3e9b122f61f8f06494c97b1afccf3-Paper.pdf "GAN"
[6]: http://www.bio.miami.edu/dana/330/330F19_9.html "Whittaker's climograph"
[7]: https://www.mit.edu/~dbertsim/papers/Optimization/Simulated%20annealing.pdf "SA"
[8]: https://ieeexplore.ieee.org/document/488968 "PSO"
[9]: https://link.springer.com/article/10.1007/s11042-020-10139-6 "GA"
[10]: https://arxiv.org/abs/1406.2572 "Saddle point problem in high dimensional non-convex optimization"

[Translational propagation LHS][1]

[MNIST dataset][2]

[Cost sensitive learning][3]

[Ensemble learning][4]

[GAN][5]

[Whittaker][6]

[Simulated Annealing][7]

[Particle Swarm Optimization][8]

[Genetic Algorithm][9]

[Saddle points in high dimensional data][10]

<div style="page-break-after: always; break-after: page;"></div>

<div style="page-break-after: always; break-after: page;"></div>
