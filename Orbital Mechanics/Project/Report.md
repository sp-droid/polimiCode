

![logo](C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\logo.png)

<center><big><big><big><b>
Politecnico di Milano
<br><br>Master in Space Engineering
<br>Academic Year 2022 / 23<br><br>Orbital mechanics course
</b></big></big></big></center>


<center><big><big><b>
Final project: Interplanetary and planetary
<br>explorer missions
</b></big></big></center>

| <img src="C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\projectlogo.PNG" alt="projectlogo" style="zoom:67%;" /> |
| :----------------------------------------------------------: |



<right><big><b>
VASP GROUP, NUMBER 2204
<br><br>    Alessandro Michelazzi (10709804)
<br><br>    Pablo Arbelo (10904636)
<br><br>    Stefano Marinelli (10705548)
<br><br>    Veronica Cerni (10700624)
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

<u>**Repeating groundtrack k:m ratio:** 2</u>

*****RAAN, argument of periapsis and initial true anomaly were left to our criteria, which were chosen in order to have a good movie and after fixing the initial date.



### 2.1.2 Assigned perturbations

- **J<sub>2</sub> effect** due to Earth's oblateness<sup><a href="#ref2.1">[2.1]</a></sup>

$$
\Large
\mathbf{a_{J_2}} =   \frac{3}{2}\frac{J_2 \mu R_E^2}{r^5}\left [ \left ( 5\frac{z^2}{r^2}-1 \right )x\widehat{\mathbf{i}} + \left ( 5\frac{z^2}{r^2}-1 \right )y\widehat{\mathbf{j}} + \left ( 5\frac{z^2}{r^2}-3 \right )z\widehat{\mathbf{k}}\right ]
$$

<u>**J<sub>2</sub>** =  0.0010826</u> as provided by the function *astroConstants* 


- **Atmospheric drag** due to the absorption and reflection of massive particles against the spacecraft<sup><a href="#ref2.1">[2.1]</a></sup>

$$
\Large
\mathbf{\ddot{r}} = -\frac{1}{2}C_D\frac{A}{m}\rho (h,t)v_{rel}\mathbf{v_{rel}}
$$

$$
\Large
where: \hspace{4mm} \mathbf{v_{rel}} = \mathbf{v}-\mathbf{\omega_E\wedge} \mathbf{r}
$$
The density at very high altitudes significantly depends day/night cycles and the solar activity at any given time, which by itself makes the modelling of aerodynamic drag a challenge. Moreover, the atmosphere is very heterogeneous in chemical composition at high h and that has to be taken into account aswell. There are several models, like the Jacchia-Bowman 2008,<sup><a href="#ref2.2">[2.2]</a></sup> that adress these issues with varying levels of success. However, we are going to employ a simple exponential approximation that depends only on the altitude:
$$
\Large
\rho = \rho_0 exp\left ( -\frac{h-h_0}{H_0} \right ) \hspace{2mm};\hspace{2mm} \rho_0 = 2\cdot 10^{-8} \hspace{2mm};\hspace{2mm} h_0 = 122 km \hspace{2mm};\hspace{2mm} H_0 = H|(\rho=3\cdot 10^{-12},h=400km)
$$

<u>**C<sub>D</sub>** = 2.1</u> ; <u>**A/m** = 0.0095 m<sup>2</sup>/kg</u>

The reference density values for both heights were selected from usual values.<sup><a href="#ref2.3">[2.3]</a></sup>

<div style="page-break-after: always; break-after: page;"></div>

### 2.1.3 Other perturbations

In the last point, the comparison with real data, we include an additional third model that features the following perturbations:

- **Atmospheric drag** (Explained earlier)
- **Solar radiation pressure**<sup><a href="#ref2.4">[2.4]</a></sup>

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

This shadow function is implemented in the code and is explained in the referenced book.<sup><a href="#ref2.4">[2.4]</a></sup>

- **Earth gravitational model 1996** (And therefore substitutes the J<sub>2</sub>-only perturbation)<sup><a href="#ref2.4">[2.4]</a></sup><sup><a href="#ref2.5">[2.5]</a></sup><sup><a href="#ref2.6">[2.6]</a></sup><sup><a href="#ref2.7">[2.7]</a></sup>

$$
\Large
U(r,\phi,\xi) = -\sum_{n=2}^{360}\frac{\mu R_E^n}{r^{n+1}}\sum_{m=0}^{n}\left ( C_{nm}cos(m\phi)+S_{nm}sin(m\phi) \right )\widetilde{P}_{nm}\left ( \xi  \right )
$$

$$
\Large
where: \hspace{4mm} r^2 = x^2+y^2+z^2\hspace{1mm};\hspace{1mm} \phi = tan^{-1}\left ( \frac{y}{x} \right )\hspace{1mm};\hspace{1mm}\xi = \frac{z}{r}\hspace{1mm};\hspace{1mm}\widetilde{P} = norm.\hspace{1mm}Legendre \hspace{1mm}polynom.
$$

Deriving -U(x,y,z) gives the perturbing accelerations we are seeking.

With this it's also possible to calculate the geoid's undulation, which shows Earth's uneven gravity field:

| ![undulation](C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\undulation.png) |
| :----------------------------------------------------------: |
|             <b>Fig.2.1 - Geoid's undulation</b>              |


- **Sun and Moon 3rd body perturbations**<sup><a href="#ref2.4">[2.4]</a></sup>

$$
\Large
\mathbf{\ddot{r}} = GM_2\left ( \frac{\mathbf{r_2-r}}{\left \| r_2-r \right \|^3}-\frac{\mathbf{r_2}}{r_2^3} \right )
$$

Where **r** and **r<sub>2</sub>** are the position vectors of the spacecraft and the planet of mass **M<sub>2</sub>** relative to the Earth, **2** being the Sun or the Moon.

- **Relativistic effects**<sup><a href="#ref2.4">[2.4]</a></sup>

$$
\Large
\mathbf{\ddot{r}} = -\frac{GM_E}{r}\left [ (4\frac{GM}{c^2r}-\frac{v^2}{c^2})\mathbf{r}+\frac{4}{c^2}(\mathbf{r}\cdot \mathbf{v})\mathbf{v} \right ]
$$



## 2.2 Groundtrack

As a spacecraft orbits around a celestial body, the path drawn on it is deviated from an ellipse by the rotation of the body and the perturbations affecting the spacecraft.

From here, the solver used is **ode113** with **1e-13** **relative** tolerance and **1e-14** **absolute** tolerance.

### 2.2.1 Orbit GT

| ![gt1](C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\gt1.png) |
| :----------------------------------------------------------: |
| ![gt1d](C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\gt1d.png) |
| ![gt10d](C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\gt10d.png) |
| <b>Fig.2.2 - 1 orbit, 1 day & 10 day groundtrack (unperturbed-perturbed)</b> |

### 2.2.2 Repeating GT

The orbit is already very close to having the k:m ratio of 2. Solve for the new semi-major axis from the equations:<sup><a href="#ref2.8">[2.8]</a></sup>
$$
\Large
T^2 = 4\pi ^2\frac{a^3}{\mu }\hspace{1mm};\hspace{1mm}\frac{T}{T_E} = \frac{m}{k}\hspace{1mm};\hspace{1mm}T_E = 23h\hspace{1mm}56m\hspace{1mm}4.1s
$$
| **Semi-major axis a** | **Repeating GT a** | **Delta** |
| ------------------- | --------------- | -------- |
| 26,619 km           | 26,562 km           | -57.236 km        |

| <img src="C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\repgt2.png" alt="repgt2"  /> |
| :----------------------------------------------------------: |
|    <b>Fig.2.3 -  k = 2, m = 1 (unperturbed-perturbed)</b>    |

However, when propagated further into the future, the perturbed one deviates significantly:

| ![repgt20](C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\repgt20.png) |
| :----------------------------------------------------------: |
|   <b>Fig.2.4 -  k = 20, m = 10 (unperturbed-perturbed)</b>   |

The aerodynamic drag reduces e but most importantly, a (smaller period, shorter gt), while J<sub>2</sub> acts on the RAAN (shifted gt). J<sub>2</sub> can be compensated with a new a, using the secular variations that we know it causes, but drag mandates performing periodic burns to recover the right value of semi-major axis.

## 2.3 Perturbed orbit propagation

When choosing the propagation time, our criteria has been propagating until having multiple periods of the keplerian elements' long term oscillations, and in our case at they do not require many orbits because our simplified perturbations do not depend on the sun or the moon, but only on the orbit's initial state. We do expect some interaction between the two, as the semi-major axis plays a role in the secular variations due to J<sub>2</sub>, so we are propagating for 10 orbits. 100 orbits were used for the filtering to minimize the effect of the edges on the secular values.

### 2.3.1 Cartesian vs. Gauss's equations methods

The Gauss's planetary equations allow us to integrate directly on the keplerian elements, after transforming the perturbations to either RSW or TNH frames, as opposed to the simple method of integrating the equation of motion in cartesian coordinates.<sup><a href="#ref2.4">[2.4]</a></sup>

After 10 orbits, the relative error between both methods continues near the set relative tolerance.

| <img src="C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\carvskep.png" alt="carvskep"  /> |
| :----------------------------------------------------------: |
|   <b>Fig.2.5 -  Relative error (car-kep)/[a,1,pi,2pi]</b>    |

In terms of computational cost, our implementation of the Gauss's planetary eqs. is about 20% slower than plain cartesian, but this is because the drag function was developed in cartesian coordinates and we transform the perturbation directly from cartesian to RSW, when it would be faster to work on TNH. When removing the drag perturbation and leaving the J<sub>2</sub>, the cartesian method is more than 30% slower. 

We find that the computational complexity is very problem dependent and since we worked out the EGM96 model, SRP and other perturbations in cartesian we are going to proceed with them for the rest of the project.

### 2.3.2 High frequency filtering

| (rounded to 3 d) | $\mathbf{\dot{a}}$ (km/s) | $\mathbf{\dot{e}}$ (-/s) | $\mathbf{\dot{i}}$  º/s | $\mathbf{\dot{\Omega }}$  º/s | $\mathbf{\dot{\omega }}$  º/s | $\mathbf{\dot{\theta }}$  º/s |
| -- | ------------------- | ---------------- | --------------- | -------- | ------------------------- | ---------------- |
| **Secular variation** | -2.93e-7      | -3.45e-12 | 2.97e-11 | -1.79e-6 | 7.32e-8          | 8.35e-3 |
| **R<sup>2</sup>** | 0.94              | 0.94         | 0.69         | 1.00  | 1.00                   | 1.00   |

The secular variations above have been obtained by applying a low-pass filter in the form of a moving window on the perturbed orbit, and a linear regression with y-intercept to get the slope. The size of the moving window is 1 orbital period. To evaluate the goodness of fit we can calculate the R-squared statistic, which tells us how well the linear fit predicts the secular variations.

The low value of R-squared for the inclination points to a lower confidence. After testing with and without each perturbation, it seems this effect comes from the J<sub>2</sub>, even though we suspected it could come from the effect of the wind on the drag's force direction, but this is one order of magnitude smaller.

| <img src="C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\filtered.png" alt="filtered"  /> |
| :----------------------------------------------------------: |
| <b>Fig.2.6 -  Filtered data in moving windows of 100, 25 & 5% of the initial orbital period</b> |

## 2.4 Evolution of the orbit

### 2.4.1 Short term

We have made two videos of the first 3 orbits. First one shows the orbit from a fixed point relative to Earth's center:

https://www.youtube.com/watch?v=gYX6A3jnGZs

And the second shows the orbit from the point of view of the spacecraft, centered on Earth:

https://www.youtube.com/watch?v=YlP3BiDUItA

### 2.4.2 Long term

| <img src="C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\longterm.png" alt="longterm"  /> |
| :----------------------------------------------------------: |
| <b>Fig.2.7 -  Long term propagation (not showing every orbit)</b> |

Propagated for two years or around 1464 orbits from the initial assigned orbital conditions and perturbations.

<div style="page-break-after: always; break-after: page;"></div>

## 2.5 Real data comparison

### 2.5.1 Object selection and procedure

The aim was to find the tracked space object with the closest 3 first orbital elements: 

| t=0 | **Semi-major axis** | **Eccentricity** | **Inclination** | RAAN | Argument of periapsis | True anomaly |
| ------------------- | ---------------- | --------------- | -------- | -------- | -------- | -------- |
| **Assigned** | 26,619 km           | 0.7452           | 62.9089º        | -       | -       | -       |
| **Real (01-11-2021)** | 26,910 km | 0.7399      | 62.5711º | 74.4297º | 278.065º | 174.3877º |

The object is designated as SL-6 R/B(2), a block-ML from a Molniya rocket that is orbiting the planet since 1999. It is catalogued in NORAD as 25850U and it's described<sup><a href="#ref2.9">[2.9]</a></sup> as a <u>900kg</u>, cylinder shaped (2.6x2.58x2.58m) rocket body, with an <u>average cross section of 7.8824</u> and radar cross section of 2 sq. meters. 

The procedure was the following: after selecting the object in tracking websites,<sup><a href="#ref2.10">[2.10]</a></sup><sup><a href="#ref2.11">[2.11]</a></sup> we obtained its TLEs from space-track.org from 1-11-2021 to 1-11-2022,<sup><a href="#ref2.12">[2.12]</a></sup> propagated them in the NASA Horizons website (to obtain very accurate ephemerides for every timestep) for the same duration but every hour,<sup><a href="#ref2.13">[2.13]</a></sup> and finally compared them with two models; the one with assigned perturbations and the full model that was detailed in [section 2.1.3](###2.1.3 Other perturbations).

The <u>coefficient of drag was assumed to be 2.2</u> like other similar space objects.<sup><a href="#ref2.14">[2.14]</a></sup> Other assumptions include the <u>temperature of the Sun's surface to be constant and equal to 5778 K</u>, <u>Earth's axial tilt to be a constant 23.436º</u> (data from January 2021) or the <u>radiation pressure coefficient to be 1.05</u>, calculated iteratively until achieving the best results and testing in the [0.7,1.3] range.

### 2.5.2 Results and discussion
| <img src="C:\Users\a1pab\Desktop\1POLIMI\polimiCode\Orbital Mechanics\Project\assets\realvsmodels.png" alt="realvsmodels" style="zoom:70%;"  /> |
| :----------------------------------------------------------: |
| <b>Fig.2.8 -  1 year propagation of historic data vs. J<sub>2</sub>+drag model vs. Full model</b> |

We are also computing the secular values. This time the moving window size chosen is half a year to account for the longest oscillation in the simulated perturbations, which is Earth's translation around the Sun.
| (rounded to 3 d) + [R<sup>2</sup>] | $\mathbf{\dot{a}}$ (km/s) | $\mathbf{\dot{e}}$ (-/s) | $\mathbf{\dot{i}}$  º/s | $\mathbf{\dot{\Omega }}$  º/s | $\mathbf{\dot{\omega }}$  º/s | $\mathbf{\dot{\theta }}$  º/s |
| ------------------- | ---------------- | --------------- | -------- | -------- | -------- | -------- |
| **NORAD 25850U** | -2.06e-8 [0.98] | 1.51e-10 /s [0.95] | 1.29e-8 º/s [0.99] | -1.54e-6 º/s [0.99] | 9.96e-8 º/s [0.99] | 7.38e-5 º/s [0.99] |
| **Full model** | -1.46e-8 [0.90]           | 1.53e-10 /s [0.95]       | 1.48e-8 º/s [0.99]      | -1.53e-6 º/s [0.99]           | 8.84e-8 º/s [0.99]            | 7.54e-5 º/s [0.99]            |
| **Assigned model** | 9.89e-10 [0.40]           | 2.19e-13 /s [0.99]       | -1.00e-11 º/s [0.99]    | -1.42e-6 º/s [0.99]           | 9.40e-8 º/s [0.99]            | 8.59e-5 º/s [0.99]            |

The full model was able to much more accurately predict the real secular effects, always being within 1 order of magnitude in absolute terms, whereas the assigned model fell very short in relation to the shape parameters, being completely incapable of predicting them. However, the J<sub>2</sub> term alone seems able to predict the other 3 variables. Especially evident are the two week and half a year oscillations of the Moon and Sun in the eccentricity, that are lacking in the second model.

On a different note, the computational cost of the full model was much higher. What it took around 5 seconds for the simple model, took 3 hours for the complex one. The most likely reason for such disparity is probably the implementation of the EGM96, which has more than 65,000 pairs of coefficients and consequently, iterations per integration step. There is probably a better way to implement it.


<div style="page-break-after: always; break-after: page;"></div>

---

# Reference

<a id="ref2.1">[2.1]</a> - [Course material Chapter 5](https://webeep.polimi.it/)
<a id="ref2.2">[2.2]</a> - Bowman, B. R., Tobiska, W. K., Marcos, F. A., Huang, C. Y., Lin, C. S., & Burke, W. J. (2008). A new empirical thermospheric density model JB2008 using new solar and geomagnetic indices. *AIAA/AAS Astrodynamics Specialist Conference and Exhibit*. [DOI](https://doi.org/10.2514/6.2008-6438)
<a id="ref2.3">[2.3]</a> - [Usual density values in high altitude](https://www.eoas.ubc.ca/courses/atsc113/flying/met_concepts/03-met_concepts/03a-std_atmos/index.html#:~:text=At%20that%20400%20km%20altitude,%2D12%20kg%2Fm3.)
<a id="ref2.4">[2.4]</a> - Montenbruck, O., & Gill, E. (2000). Satellite Orbits. *Satellite Orbits*. [DOI](https://doi.org/10.1007/978-3-642-58351-3)
<a id="ref2.5">[2.5]</a> - [EGM96 coefficients](https://cddis.nasa.gov/926/egm96/getit.html)
<a id="ref2.6">[2.6]</a> - [EGM96 equations](https://people.sc.fsu.edu/~lb13f/projects/space_environment/egm96.php)
<a id="ref2.7">[2.7]</a> - [Fully normalized Legendre polynomials](http://mitgcm.org/~mlosch/geoidcookbook/node11.html)
<a id="ref2.8">[2.8]</a> - [Course material Chapter 2](https://webeep.polimi.it/)
<a id="ref2.9">[2.9]</a> - [Detailed tracking-discosweb.esoc.esa.int](https://discosweb.esoc.esa.int/objects?page=259&pageSize=100)
<a id="ref2.10">[2.10]</a> - [Tracking-satview.org](https://www.satview.org/?sat_id=25850U)
<a id="ref2.11">[2.11]</a> - [Tracking-2yo.com](https://www.n2yo.com/satellite/?s=25850#results)
<a id="ref2.12">[2.12]</a> - [TLEs from space-track.org](https://www.space-track.org)
<a id="ref2.13">[2.13]</a> - [TLE NASA-Horizons propagator](https://ssd.jpl.nasa.gov/horizons/app.html#/)
<a id="ref2.14">[2.14]</a> - Mostaza Prieto, D., Graziano, B. P., & Roberts, P. C. E. (2014). Spacecraft drag modelling. In *Progress in Aerospace Sciences* (Vol. 64). [DOI](https://doi.org/10.1016/j.paerosci.2013.09.001)