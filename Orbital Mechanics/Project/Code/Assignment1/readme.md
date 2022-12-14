# Specific requirements

**GROUP ID**: 2204

**Departure**: Earth

**Flyby planet**: Saturn

**Arrival NEO**: 85

**Earliest Departure**: 30-September-2026

**Latest Arrival**: 28-March-2061


# General requirements

- Use the method of **patched conics**
- **Do not consider planetary departure and insertion**, that is:
  - Initial heliocentric orbit is equal to that of the departure planet
  - Final heliocentric orbit is equal to that of the arrival asteroid
- Use **uplanet** (in WeBeep) to compute the ephemerides of planets, and **ephNEO** (in WeBeep) for the ephemerides of NEOs (see comments in the code for names of NEOs)
- The figure of merit for the mission is the **total cost in terms of Ξπ£tot**
  - Other criteria should be taken into account, such as altitude restrictions during the flyby

# General outputs

The mission analysis should cover the following points:

1. **<u>Design process</u>**, detailing:
	- **Initial choice for the time windows**, justifying it based on the characteristics of the mission
		- Do not just take the whole time interval provided in the mission requirements for both departure and arrival windows. Choose and justify them based on the characteristics of your mission.
	- **Additional constraints** considered (such as minimum altitude of the closest approach during the flyby)
	- **Strategy followed** to explore, analyze and compare the different transfer options
	- **Justified selection of a final solution**
	- **Plots and data** supporting your design choices (e.g., π₯π cost plots, preliminary estimates,β¦ )
2. **<u>Final solution</u>**, including:
	- **Heliocentric trajectory**
		- Departure, flyby, and arrival times
		- Characterization of the interplanetary transfer arcs
		- Plot of the heliocentric trajectory, together with the orbits of the three celestial objects and their positions at departure, flyby, and arrival
	- **Flyby** (powered gravity assist)
		- Altitude of the closest approach
		- Time duration of the flyby (considering a finite SOI)
		- Comparison of the total velocity change due to the flyby Ξπ£fb with the cost of the powered maneuver at perigee Ξπ£ga
		- Plot of the incoming and outcoming hyperbola arcs
	- **Cost of the mission in terms of π«ππ­π¨π­**
		- Detail the separate values of Ξπ£dep, Ξπ£ga, and Ξπ£arr
