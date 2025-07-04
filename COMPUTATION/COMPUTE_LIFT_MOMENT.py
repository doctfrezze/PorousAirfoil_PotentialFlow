# FUNCTION - COMPUTE LIFT (CL) AND MOMENT (CM) COEFFICIENTS FROM PRESSURE DISTRIBUTION IN SOURCE PANEL METHOD

import numpy as np
import math as math
def COMPUTE_LIFT_MOMENT(Cp,S,beta,phi,AoAR,XC,YC):
    # Compute normal and axial force coefficients
    CL = sum(-Cp*S*np.sin(beta))                                                         # Normal force coefficient []
    CM = sum(Cp*(XC-0.25)*S*np.cos(phi) + Cp*YC*S*np.sin(phi))                           # Moment coefficient []

    return CL,CM