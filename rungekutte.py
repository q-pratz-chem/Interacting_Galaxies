import numpy as np
from config import *

def get_path_err(Ri):
    """
    Get analytical solution of galaxy path for pair of coordinates, [x, y].
    
    The analytic solution of an ellipse is:
        (x**2 / a**2) + (y**2 / b**2) = 1
    
    We calculate the left-hand side of the equation, and
      return the error = (1 - LHS).
    
    """
    
    x, y = Ri 
    
    # parabolic case
    if e == 1.:
        r_test = np.sqrt(x**2 + y**2)
        phi_test = np.tan(x / y)
        
        r_true = 2*rmin / (1 + np.cos(phi))
    
    # closed orbit case
    else:
        b = a * np.sqrt(1 - e**2)
        test_point = (x**2 / a**2) + (y**2 / b**2)
        err = 1 - test_point


    return err
    

def rk4(t0, tmax, r0, phi0, f, h, dr_tol, r_ref):
    """
    Calculate orbit of a body about center of mass using 4th-order Runge-Kutte
    method.
    
    Inputs:
      t0, tmax: start, end time of simulation, s
      r0: initial position of second body, m
      phi0: initial angular position of second body, rad
      f: function for getting r', phi'
      h: time interval
      dr_tol: amount of allowed deviation from r_ref
      r_ref: expected value for r
    
    Returns:
      R: 2-D array of x,y coordinates
      err: relative error of each (x,y) pair
      
    """
    
    # initial conditions
    x0 = r0 * np.cos(phi0)
    y0 = r0 * np.cos(phi0)
    R0 = [x0, y0]

    # initialize arrays for t, orbital elements
    t = np.array(t0)
    ti = t0
    Ri = np.array(R0)
    R = np.array([Ri])
    
    debug = True
    
    # break if r goes outside tolerance
    err = np.array([])
    err_i = get_path_err(Ri)
    
    while err_i < dr_tol:
        # break if t = tmax
        if ti > tmax:
            if debug: print(f"Time: {ti} > {tmax}")
            break
        else:
            # calculate y(t), y'(t)
            if debug: print(ti, Ri, np.abs(r_check - r_ref))
            err = np.append(err, get_path_err(Ri))
            
            # compute k1, k2, k3, k4 for Ri (r_i, phi_i, r_i', phi_i')
            k1 = h*f(Ri)
            k2 = h*f(Ri + k1*0.5) 
            k3 = h*f(Ri + k2*0.5) 
            k4 = h*f(Ri + k3)  

            # update ti, Ri
            ti += h
            Ri += (1/6)*(k1 + 2*k2 + 2*k3 + k4)

            # update r_check
            r_check = Ri

            # save new ri to r
            t = np.append(t, ti)
            R = np.append(R, np.array([Ri]), axis=0)
            
#     # convert results to cartesian coords
#     # r: [r, phi, r', phi']
#     x = R[:,0]*np.cos(R[:,1])  # r * cos phi
#     y = R[:,0]*np.sin(R[:,1])  # r * sin phi

    return R, err