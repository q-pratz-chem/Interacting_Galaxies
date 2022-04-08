import numpy as np

def rk4(m0, m1, t0, tmax, r0, phi0, f, h, dr_tol, r_ref):
    """
    Calculate orbit of a body about center of mass using 4th-order Runge-Kutte
    method.
    
    Inputs:
      m0, m1: mass of each body, kg
      t0, tmax: start, end time of simulation, s
      r0: initial position of second body, m
      phi0: initial angular position of second body, rad
      f: function for getting r', phi'
      h: time interval
      dr_tol: amount of allowed deviation from r_ref
      r_ref: expected value for r
    
    """
    
    # initial conditions
    R0 = [r0, phi0]

    # initialize arrays for t, orbital elements
    t = np.array(t0)
    ti = t0
    Ri = np.array(R0)
    R = np.array([Ri])
    
    debug = False
    
    # break if r goes outside tolerance
    r_check = Ri[0]
    err = np.array([])
    
    while np.abs(r_check - r_ref) < dr_tol:
        # break if t = tmax
        if ti > tmax:
            if debug: print(f"Time: {ti} > {tmax}")
            break
        else:
            # calculate y(t), y'(t)
            if debug: print(ti, Ri, np.abs(r_check - r_ref))
            err = np.append(err, (r_ref - r_check)/r_ref)
            
            # compute k1, k2, k3, k4 for Ri (r_i, phi_i, r_i', phi_i')
            k1 = h*f(Ri)
            k2 = h*f(Ri + k1*0.5) 
            k3 = h*f(Ri + k2*0.5) 
            k4 = h*f(Ri + k3)  

            # update ti, Ri
            ti += h
            Ri += (1/6)*(k1 + 2*k2 + 2*k3 + k4)

            # update r_check
            r_check = Ri[0]

            # save new ri to r
            t = np.append(t, ti)
            R = np.append(R, np.array([Ri]), axis=0)
            
    # convert results to cartesian coords
    # r: [r, phi, r', phi']
    x = R[:,0]*np.cos(R[:,1])  # r * cos phi
    y = R[:,0]*np.sin(R[:,1])  # r * sin phi

    return R, x, y, err