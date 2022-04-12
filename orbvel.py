import numpy as np
from config import *

def f_xy(rs, **kwargs):
    """
    Calculate f(x,y,t) = x'' for a given [x,y]
    Inputs:
        rs, 2-D array [x,y]
    
    Returns 2-D array [x'', y'']
    
    """
    
    x = rs[0]
    y = rs[1]
    coeff = G*M/np.pow(x**2 + y**2, 1.5)
    return coeff*rs


def v(rs, **kwargs):
    """
    Evaluate first time derivative of r, phi for closed orbits.
    
    Inputs: r -- 2D vector w/ elems r (separation) 
                 and phi (angle)
    Returns: 2D array of functions v(theta) and v(omega) 
             evaluated at these points
    """
    r, phi = rs[0], rs[1]
    vphi = np.sqrt(G*M*p)/r**2  # phi'
    vr = (e/p)*r**2*np.sin(phi)*vphi
    return np.array([vr, vphi], float)

def f(rs, **kwargs):
    """
    Function defining system of equations of motion for two bodies
      orbiting about a center of mass
    
    Inputs: r -- 2D vector w/ elems r (distance from center of mass),
                    phi (angle from distance of closest approach)
    Returns: 2D array of functions f(theta) and f(phi) 
             evaluated at r, phi
    """
    r, phi = rs[0], rs[1]
    vr, vphi = v(rs)
    #fphi = np.sqrt(G*M*p)/r**2  # this is actually phi'
    fphi = -2*np.sqrt(G*M*p)*r**-3*vr  # f(phi) = phi''
    # fr = ecc*r**2*np.sin(phi)*fphi/p
    fr = e*np.sqrt(G*M/p)*np.cos(phi)*vphi  # f(r) = r''
    
    return np.array([fr, fphi], float)

# def f(r, **kwargs):
#     """
#     Function defining system of equations of motion for two bodies
#       orbiting about a center of mass
#     Inputs: r -- 2D vector w/ elems r (distance from center of mass),
#                     phi (angle from distance of closest approach)
#     Returns: 2D array of functions f(theta) and f(phi) 
#              evaluated at r, phi
#     """
#     G = kwargs.get('G', 4.4985e-12) # kpc^3/Msun/Myr^2
#     M = kwargs.get('M', 1e12)      # Msun
#     p = kwargs.get('p', 1)      # kpc
#     e = kwargs.get('e', 0.)        # dimensionless
#     r, phi = r[0], r[1]
#     fphi = np.sqrt(G*M*p)/r**2        # phi', angular velocity
#     fr = e*r**2*np.sin(phi)*fphi/p    # r', radial velocity
#     return np.array([fr, fphi], float)