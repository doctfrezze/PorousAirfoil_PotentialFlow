# FUNCTION - COMPUTE LIFT (CL) AND MOMENT (CM) COEFFICIENTS FROM PRESSURE DISTRIBUTION IN SOURCE PANEL METHOD

import numpy as np
import math as math
def COMPUTE_LIFT_MOMENT(Cp,S,beta,phi,AoAR,XC):
    # Compute normal and axial force coefficients
    CN = -Cp*S*np.sin(beta)                                                         # Normal force coefficient []
    CA = -Cp*S*np.cos(beta)                                                         # Axial force coefficient []

    # Compute lift and moment coefficients
    CL = sum(CN*np.cos(AoAR)) - sum(CA*np.sin(AoAR))                                # Decompose axial and normal to lift coefficient []
    CM = sum(Cp*(XC-0.25)*S*np.cos(phi))                                            # Moment coefficient []

    return CL,CM