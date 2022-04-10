import numpy as np
from config import * 

def f(r, **kwargs):
    """
    Function defining system of equations of motion for two bodies
      orbiting about a center of mass
    
    Inputs: r -- 2D vector w/ elems r (distance from center of mass),
                    phi (angle from distance of closest approach)
    Returns: 2D array of functions f(theta) and f(phi) 
             evaluated at r, phi
    """
#     G = kwargs.get('G', 4.4985e-12) # kpc^3/Msun/Myr^2
#     M = kwargs.get('M', 1e12)       # Msun
#     p = kwargs.get('p', 1)          # kpc
#     e = kwargs.get('e', 0.)         # dimensionless
    
    r, phi = r[0], r[1]
    
    fphi = np.sqrt(G*M*p)/r**2      # phi', angular velocity
    fr = e*r**2*np.sin(phi)*fphi/p  # r', radial velocity
    
    return np.array([fr, fphi], float)