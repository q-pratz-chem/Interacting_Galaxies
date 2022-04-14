import numpy as np
from config import *

def rk4(t0, tmax, r0, v0, f, h, err_tol, r_ref):
    """
    Calculate orbit of a body about center of mass using 4th-order Runge-Kutte
    method.
    
    Inputs:
      t0, tmax: start, end time of simulation, Myr
      r0: initial position of second body, kpc
      v0: initial velocity (in +y direction) of second body, kpc/Myr
      f: function for getting vx', vy'
      h: time interval
      err_tol: amount of allowed deviation from r_ref
      r_ref: expected value for r2
    
    Returns:
      R: 2-D array of x1,y1,x2,y2 coordinates at each step
      err: relative error of path compared to circle/ellipse/etc.
      
    """
    
    # Initialize arrays for time, position vector, & velocity vector
    ti = t0
    t = np.array(ti)
    # mass 1 starts at point of furthest separation
    x1 = r0  # kpc
    y1 = 0   # kpc
    vx1 = 0
    vy1 = v0
    # use center of mass m1x1 = m2x2 to get position, vel of mass 2
    x2 = -(m1/m2)*r0
    y2 = 0
    vx2 = 0
    vy2 = -(m1/m2)*v0
     
    print(f"initial R: {[x1, y1, x2, y2]}, f(R) = {f(np.array([x1, y1, x2, y2]))}")
    # initialize arrays of positions at each step
    ri = np.array([x1, y1, x2, y2], float)
    
    # initialize arrays of velocities at each step
    vinit = np.array([vx1, vy1, vx2, vy2], float)
    vi = vinit
    v_ = np.array([vi]) # Initial vx, vy
    
    # initialize position-velocity mega-vector
    r = np.concatenate([ri,vi])
    R = np.array([r])
    
    
    debug = False
    
    # Initialize steps for error check for circular orbit
    r_check = np.sqrt(r[0]**2 + r[1]**2)  # current radius, r
    rel_err = (r_ref - r_check)/r_ref
    print(f"init err: {rel_err}")
    err = np.array([rel_err])
    
    # solve eqns of motion using R-K method
    while np.abs(rel_err) < err_tol:
        # break if t = tmax
        if ti > tmax:
            if debug: print(f"Time: {ti} > {tmax}")
            break
        else:
            # calculate r(t), r'(t)
            
            # compute k1, k2, k3, k4 for r
            k1  = h*get_f(r, f)
            k2 = h*get_f(r + k1*0.5, f) 
            k3 = h*get_f(r + k2*0.5, f) 
            k4 = h*get_f(r + k3, f)  

            # update ti, Ri
            ti += h
            r += (1/6)*(k1 + 2*k2 + 2*k3 + k4)

            # calculate new rel error
            r_check = np.sqrt(r[0]**2 + r[1]**2)  # current radius, r
            rel_err = (r_ref - r_check)/r_ref
            if debug: print(f"{ti} {r_check} {rel_err}")

            # save new ri to r
            t = np.append(t, ti)
            R = np.append(R, np.array([r]), axis=0)
            err = np.append(err, rel_err) # relative error
            
    # components of cartesian coords
    # R: [x1, x2, y1, y2]
    X = R[:,[0,2]]  # x1, x2
    Y = R[:,[1,3]]  # y1, y2
    return t, X, Y, err # v_ , err
