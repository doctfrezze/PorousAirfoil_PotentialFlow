# FUNCTION - COMPUTE Nx AND Ny GEOMETRIC INTEGRALS FOR VORTEX PANEL METHOD

import numpy as np
import math as math

np.seterr('raise')

def STREAMLINE_VPM(XP,YP,XB,YB,phi,S):
    
    # Number of panels
    numPan = len(XB)-1                                                          # Number of panels (control points)
    
    # Initialize arrays
    Nx = np.zeros(numPan)                                                       # Initialize Nx integral array
    Ny = np.zeros(numPan)                                                       # Initialize Ny integral array
    
    # Compute Nx and Ny
    for j in range(numPan):                                                     # Loop over all panels
        # Compute intermediate values
        A = -(XP-XB[j])*np.cos(phi[j]) - (YP-YB[j])*np.sin(phi[j])              # A term
        B  = (XP-XB[j])**2 + (YP-YB[j])**2                                      # B term
        Cx = np.sin(phi[j])                                                     # Cx term (X-direction)
        Dx = -(YP-YB[j])                                                        # Dx term (X-direction)
        Cy = -np.cos(phi[j])                                                    # Cy term (Y-direction)
        Dy = XP-XB[j]                                                           # Dy term (Y-direction)
        E  = math.sqrt(B-A**2)                                                  # E term
        if (E == 0 or np.iscomplex(E) or np.isnan(E) or np.isinf(E)):           # If E term is 0 or complex or a NAN or an INF
            Nx[j] = 0                                                           # Set Nx value equal to zero
            Ny[j] = 0                                                           # Set Ny value equal to zero
        else:
            # Compute Nx, Ref [1]
            term1 = 0.5*Cx*np.log((S[j]**2 + 2*A*S[j]+B)/B);                    # First term in Nx equation
            term2 = ((Dx-A*Cx)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));   # Second term in Nx equation
            Nx[j] = term1 + term2;                                              # Compute Nx integral
            
            # Compute Ny, Ref [1]
            term1 = 0.5*Cy*np.log((S[j]**2 + 2*A*S[j]+B)/B);                    # First term in Ny equation
            term2 = ((Dy-A*Cy)/E)*(math.atan2((S[j]+A),E) - math.atan2(A,E));   # Second term in Ny equation
            Ny[j] = term1 + term2;                                              # Compute Ny integral
            
        # Zero out any problem values
        if (np.iscomplex(Nx[j]) or np.isnan(Nx[j]) or np.isinf(Nx[j])):         # If Nx term is complex or a NAN or an INF
            Nx[j] = 0                                                           # Set Nx value equal to zero
        if (np.iscomplex(Ny[j]) or np.isnan(Ny[j]) or np.isinf(Ny[j])):         # If Ny term is complex or a NAN or an INF
            Ny[j] = 0                                                           # Set Ny value equal to zero
    
    return Nx, Ny                                                               # Return both Nx and Ny matrices
