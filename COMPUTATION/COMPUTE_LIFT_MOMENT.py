# FUNCTION - COMPUTE LIFT (CL) AND MOMENT (CM) COEFFICIENTS FROM PRESSURE DISTRIBUTION IN SOURCE PANEL METHOD

import numpy as np
import math as math
def COMPUTE_LIFT_MOMENT(Cp,S,beta,phi,AoAR,XC,YC,low_point = [],high_point = []):
    if low_point == [] and high_point == []:
        CL = sum(-Cp*S*np.sin(beta))
        CM = sum(Cp*(XC-0.25)*S*np.cos(phi) + Cp*YC*S*np.sin(phi))                           # Moment coefficient []
    # Compute normal and axial force coefficients
    else:
        CL = sum(-Cp[i]*S[i]*np.sin(beta[i]) for i in low_point)                             # Normal force coefficient []
        CL += sum(-Cp[i]*S[i]*np.sin(beta[i]) for i in high_point)
        CM = sum(Cp[i]*(XC[i]-0.25)*S[i]*np.cos(phi[i]) + Cp[i]*YC[i]*S[i]*np.sin(phi[i]) for i in low_point)           
        CM += sum(Cp[i]*(XC[i]-0.25)*S[i]*np.cos(phi[i]) + Cp[i]*YC[i]*S[i]*np.sin(phi[i]) for i in high_point)                 

    return CL,CM