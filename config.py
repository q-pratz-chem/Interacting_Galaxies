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
m1 = 1e12 # Msun
m2 = 1e12 # Msun
M = m1 + m2

# can get a = (G*M*T**2/(4*np.pi**2))**(1/3)
a = 100    # semi-major axis (kpc, Kepler's Law)
ecc = 1e-12 # 0.5
p = a*(1 - ecc**2)   # semi-latus rectum (kpc)

# starting position of second galaxy 
r0 = (m2/M)*a   #for circular orbit
    
# starting velocity of the second galaxy
v0 = np.sqrt(G*m2/a)*np.sqrt((m2*(1-ecc**2))/M)

# make sure r is close to analytic value over time
err_tol = 1e-1  # allowable error 
