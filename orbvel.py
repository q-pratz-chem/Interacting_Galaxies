import numpy as np
from config import *

def v_init(r, phi, **kwargs):
    """
    Evaluate first time derivative of r, phi.
    
    Inputs: r -- 2D vector w/ elems r (separation) 
                 and phi (angle)
    Returns: 2D array of functions v(theta) and v(omega) 
             evaluated at these points
    """
        
    r0 = r
    phi0 = -np.arccos((p - r0)/(r0*e))
    
    # first derivative of (r0, phi0)
    fphi0 = np.sqrt(G*M*p)/r0**2 
    fr = e*r0**2*np.sin(phi0)*fphi0/p

    # starting velocity of the second galaxy (cartesian)
    vx0 = -r0*np.sin(phi0)*fphi0 + fr*np.cos(phi0)
    vy0 = r0*np.cos(phi0)*fphi0 + fr*np.sin(phi0)

    return np.array([vx0, vy0], float)

def f(r, **kwargs):
    """
    Evaluate acceleration of a body with respect to the center of mass
      with another body, f(x), f(y)
      
    Inputs: 
        Pair of coordinates [x,y]
    
    Returns:
        Array, [f(x), f(y)]
    """
        
    x, y = rs[0], rs[1]  
    r = np.sqrt(x**2 + y**2)
    
    # x, y components of Acceleration 'a'
    ax = G*M*x/r**3    # m_rest: mass of the body at rest
    ay = G*M*y/r**3    
    
    return np.array([ax, ay], float)


# def f(rs, **kwargs):
#     """
#     Function defining system of equations of motion for two bodies
#       orbiting about a center of mass
    
#     Inputs: r -- 2D vector w/ elems r (distance from center of mass),
#                     phi (angle from distance of closest approach)
#     Returns: 2D array of functions f(theta) and f(phi) 
#              evaluated at r, phi
#     """
#     r, phi = rs[0], rs[1]
#     vr, vphi = v(rs)
#     #fphi = np.sqrt(G*M*p)/r**2  # this is actually phi'
#     fphi = -2*np.sqrt(G*M*p)*r**-3*vr  # f(phi) = phi''
#     # fr = ecc*r**2*np.sin(phi)*fphi/p
#     fr = e*np.sqrt(G*M/p)*np.cos(phi)*vphi  # f(r) = r''
    
#     return np.array([fr, fphi], float)

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