# FUNCTION - COMPUTE TANGENTIAL VELOCITY AND PRESSURE COEFFICIENTS FOR SPVP METHOD

import numpy as np
import math as math
import matplotlib.pyplot as plt

def COMPUTE_Vt_Cp(numPan,Vinf,beta,lam,gamma,J,L,rhoinf,Delta_Cp,Rs,a,n,pore_entry,pore_out,omega_in,omega_out,b,delta,XC):
    # Compute velocities
    Vt = np.zeros(numPan)                                                           # Initialize tangential velocity
    Vn = np.zeros(numPan)
    V = np.zeros(numPan) 
    Cp = np.zeros(numPan)                                                           # Initialize pressure coefficient
    for i in range(numPan):                                                         # Loop over all panels
        normal_vs_pore_in = delta[i]-np.pi*(omega_in)/180
        normal_vs_pore_out = delta[i]-np.pi*(omega_out)/180
        term1 = Vinf*np.sin(beta[i])                                                # Uniform flow term
        term2 = (1/(2*np.pi))*sum(lam*J[i,:])                                       # Source panel terms when j is not equal to i
        term3 = gamma/2                                                             # Vortex panel term when j is equal to i
        term4 = -(gamma/(2*np.pi))*sum(L[i,:])                                      # Vortex panel terms when j is not equal to i
        V_mean = 0.5*rhoinf*Vinf**2*Delta_Cp/Rs*n/a
        if abs(V_mean)>0.6*Vinf:
            V_mean = 0.6*Vinf*V_mean/abs(V_mean)
        if i in pore_entry:
            term5 = V_mean*np.sin(normal_vs_pore_in) 
        elif i in pore_out:
            term5 = V_mean*np.sin(normal_vs_pore_out)
        else:
            term5 = 0
        Vt[i] = term1 + term2 + term3 + term4 + term5                                      # Compute tangential velocity on panel i
        Vn[i] = b[i] + Vinf*2*np.pi*np.cos(beta[i])
        V[i] =math.sqrt(Vt[i]**2+Vn[i]**2)
        Cp[i] = 1-(V[i]/Vinf)**2                                                   # Compute pressure coefficient on panel i
    """print('V_mean. = ',V_mean)
    print('V = ',V[pore_out[int(len(pore_out)/2)]])
    print('Vn = ',Vn[pore_out[int(len(pore_out)/2)]])
    print('Vt = ',Vt[pore_out[int(len(pore_out)/2)]])
    print('V_mean*sin() = ',V_mean*np.sin(delta[pore_out[int(len(pore_out)/2)]]-np.pi*(omega_in)/180) )
    print(omega_out)
    print(delta[pore_out[int(len(pore_out)/2)]]*180/np.pi)"""

    return Vt, Cp
