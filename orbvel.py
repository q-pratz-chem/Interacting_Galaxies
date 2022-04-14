import numpy as np
from config import *

def f(rs, **kwargs):
    """
    Evaluate the second derivative of the position vector in 2D.
    
    Inputs:  array with x-, y- components of the position vector
    Returns: array of x-,y- components of acceleration 
             evaluated at these points
    """
        
    # position of mass1, mass2    
    x1, y1 = rs[0], rs[1]  
    x2, y2 = rs[2], rs[3] 
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    # x, y components of Acceleration a(t) = x''(t), which is also v'(t)
    ax1 = G*m2*(x2 - x1)/r**3    
    ay1 = G*m2*(y2 - y1)/r**3    
    ax2 = G*m1*(x1 - x2)/r**3    
    ay2 = G*m1*(y1 - y2)/r**3    
    
    return np.array([ax1, ay1, ax2, ay2], float)

def get_f(r, f):
    """
    Calculate updated x, y, vx, vy for each particle using a function 
        f(x,y,t) = v'(t)
    
    Inputs:
        r, 1-D array [x1,y1,x2,y2,vx1, vy1, vx2, vy2]
        
    Returns:
        fxy, 1-D array v(t), f(x,y,t)
        
    """
    
    x1, y1 = r[0], r[1]
    vx1, vy1 = r[4], r[5]
    x2, y2 = r[2], r[3]
    vx2, vy2 = r[6], r[7]
    
    fxy = np.array([vx1, vy1, vx2, vy2])
    fv = f(np.array([x1, y1, x2, y2]))
    
    return np.array(np.concatenate([fxy, fv]))
    