####
# Orbital parameters and initial conditions configuration file
#####

import numpy as np
import astropy.units as u
import astropy.constants as const


# UNITS: kpc, Msun, Myr
# constants
G = (const.G).to(u.kpc**3 * u.Msun**-1 * u.Myr**-2).value # kpc^3/Msun/Myr

# orbit parameters
m0 = 1e12 # Msun
m1 = 5e12 # Msun
M = m0 + m1

a = 1 # (G*M*T**2/(4*np.pi**2))**(1/3)   # semi-major axis (kpc, Kepler's Law)
e = 1e-12 # 0.5
p = a*(1 - e**2)   # semi-latus rectum (kpc)
