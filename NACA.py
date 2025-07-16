# FUNCTION - GENERATE AIRFOIL COORDINATES FOR NACA 4-DIGIT SERIES AIRFOILS

import numpy as np
import math as math

def GENERATE_NACA4(NameAirfoil, c=1.0, n=100, power=1.0):
    m = int(NameAirfoil[0])*0.01
    p = int(NameAirfoil[1])*0.1
    t = int(NameAirfoil[2:4])*0.01
    beta = np.linspace(0, np.pi, n)
    x_dist = (1 - np.cos(beta)) / 2
    x = (x_dist**power) * c
    #yt = np.sqrt((0.25-(x-0.5)*(x-0.5))/9)
    yt = 5 * t * c * (0.2969 * np.sqrt(x / c) - 0.1260 * (x / c) - 0.3516 * (x / c)**2 
                      + 0.2843 * (x / c)**3 - 0.1015 * (x / c)**4)
    
    if m == 0 and p == 0:
        xu, yu = x, yt
        xl, yl = x, -yt
    else:
        yc = np.where(x < p * c, m * x / (p**2) * (2 * p - x / c),
                      m * (c - x) / ((1 - p)**2) * (1 + x / c - 2 * p))
        dyc_dx = np.where(x < p * c, 2 * m / (p**2) * (p - x / c),
                          2 * m / ((1 - p)**2) * (p - x / c))
        theta = np.arctan(dyc_dx)
        xu = x - yt * np.sin(theta)
        yu = yc + yt * np.cos(theta)
        xl = x + yt * np.sin(theta)
        yl = yc - yt * np.cos(theta)
    
    return np.append(xl[::-1], xu[1:]), np.append(yl[::-1], yu[1:])

GENERATE_NACA4('0018')