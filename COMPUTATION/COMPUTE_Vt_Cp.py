# FUNCTION - COMPUTE TANGENTIAL VELOCITY AND PRESSURE COEFFICIENTS FOR SPVP METHOD

import numpy as np
import math as math

def COMPUTE_Vt_Cp(numPan,Vinf,beta,lam,gamma,J,L):
    # Compute velocities
    Vt = np.zeros(numPan)                                                           # Initialize tangential velocity
    Cp = np.zeros(numPan)                                                           # Initialize pressure coefficient
    for i in range(numPan):                                                         # Loop over all panels
        term1 = Vinf*np.sin(beta[i])                                                # Uniform flow term
        term2 = (1/(2*np.pi))*sum(lam*J[i,:])                                       # Source panel terms when j is not equal to i
        term3 = gamma/2                                                             # Vortex panel term when j is equal to i
        term4 = -(gamma/(2*np.pi))*sum(L[i,:])                                      # Vortex panel terms when j is not equal to i
        
        Vt[i] = term1 + term2 + term3 + term4                                       # Compute tangential velocity on panel i
        Cp[i] = 1-(Vt[i]/Vinf)**2                                                   # Compute pressure coefficient on panel i

    return Vt, Cp
