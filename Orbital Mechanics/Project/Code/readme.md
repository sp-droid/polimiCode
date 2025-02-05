# Specific requirements

**GROUP ID**: 2204

**a**: 2.6619 x10^4 km

**e**: 0.7452

**i**: 62.9089 deg

**Repeating GT ratio k-m**: 2:1

**Perturbation 1**: J2

**Perturbation 2**: Air drag (cD = 2.1, A/M = 0.0095 m^2/kg)

# General requirements

1. **<u>Nominal orbit</u>**, indicating its initial values and main characteristics
- <u>Data in WeBeep does not include Ω, 𝜔, and 𝑓0. You can choose them freely</u>
2. **<u>Ground track</u>**
  - **Plot the ground track of the nominal orbit** over **1 orbit, 1 day, and 10 days**, for the **unperturbed** 2BP
  - **Modify the semimajor axis to obtain a repeating ground track**, and plot it for the **unperturbed** 2BP:
    - Use the ratio for satellite and Earth revolutions given in WeBeep
  - **Plot again the ground tracks for the nominal orbit and the modified orbit obtained in the last point, adding the**
    **assigned perturbations to the orbit propagation** (𝐽2 + see table in WeBeep)
    - Does the repeating ground track solution from the modified orbit still work under the presence of perturbations? Why?
      **IMPORTANT**: The modified value of semimajor axis should only be used for the ground track analysis. **For the rest of the assignment, use the nominal value given in WeBeep**
3. **<u>Propagate the orbit with the assigned perturbations</u>** (𝐽2 + see table in WeBeep), in:
- **Cartesian coordinates**
- **Keplerian elements through Gauss’s planetary equations**
4. **<u>Plot the history of the Keplerian elements</u>**:
- Choose an **adequate and reasoned propagation time based on the assigned perturbations** (that is, <u>long</u>
  <u>enough to observe the effect of all the perturbations</u>)
- Choose **proper units for time**, use **degrees for angles**
- **Compare and analyse the evolution of each element**
- **Compare both propagation methods** (<u>in terms of relative error, computational time</u>, etc)
5. **<u>Represent the evolution of the orbit</u>** (image or movie)
6. **<u>Filtering of high frequencies</u>**:
- **Use a low-pass filter** (e.g., movmean) to remove high frequencies in the orbital elements, retrieving
  the **long-period and/or secular evolution**. **You can do more than one** filter (i.e., with different cut-off
  frequencies) for the different time scales
- **Plot** the results (you can plot together filtered
  and unfiltered evolution)
7. **<u>Comparison with real data</u>**:
- **Select** an object (i.e., satellite, debris, rocket body,
  etc.) **in the same** orbital **region**, and download its
  orbital elements (<u>for a significant time span</u>)
- **Propagate its orbit using your model**, <u>using as</u>
  <u>initial condition the orbital elements for the real</u>
  <u>object at the initial time</u>
- **Compare the downloaded** elements **with the**
  **results from your model**.