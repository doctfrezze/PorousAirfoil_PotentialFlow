# FUNCTION - COMPUTE RHS VECTOR b FOR PANEL VORTEX/PANEL METHOD INCLUDING KUTTA CONDITION

import numpy as np
import math as math

def COMPUTE_b_SPVP(numPan,Vinf,beta,rhoinf,Delta_Cp,Rs,a,n,pore_entry,pore_out,omega_in,omega_out,delta):
    b = np.zeros(numPan)                                                            # Initialize the b array
    for i in range(numPan):                                                         # Loop over all i panels (rows)
        normal_vs_pore_in = delta[i]-np.pi*(omega_in)/180
        normal_vs_pore_out = delta[i]-np.pi*(omega_out)/180
        V_mean = 0.5*rhoinf*Vinf**2*Delta_Cp/Rs*n/a
        if abs(V_mean)>0.6*Vinf:
            V_mean = 0.6*Vinf*V_mean/abs(V_mean)
        if i in pore_entry:
            b[i] = -Vinf*2*np.pi*np.cos(beta[i]) + V_mean*np.cos(normal_vs_pore_in) # Compute RHS array
        elif i in pore_out:
            b[i] = -Vinf*2*np.pi*np.cos(beta[i]) + V_mean*np.cos(normal_vs_pore_out)
        else:
            b[i] = -Vinf*2*np.pi*np.cos(beta[i])                                   

    

    # Last element of b array (Kutta condition)
    b = np.append(b,-Vinf*2*np.pi*(np.sin(beta[0]) + np.sin(beta[numPan-1])))       # Add Kutta condition equation RHS to b array
    return b