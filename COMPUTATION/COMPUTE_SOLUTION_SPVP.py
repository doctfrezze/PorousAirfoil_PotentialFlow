# FUNCTION - SOLVE LINEAR SYSTEM FOR SOURCE STRENGTHS AND VORTEX STRENGTH IN SPVP METHOD

import numpy as np
import math as math

def COMPUTE_SOLUTION_SPVP(A,b):
    # Compute result array
    resArr = np.linalg.solve(A,b)                                                   # Solve system of equation for all source strengths and single vortex strength

    # Separate lam and gamma values from result 
    lam   = resArr[0:len(resArr)-1]                                                 # All panel source strengths
    gamma = resArr[len(resArr)-1]                                                   # Constant vortex strength
    return lam, gamma