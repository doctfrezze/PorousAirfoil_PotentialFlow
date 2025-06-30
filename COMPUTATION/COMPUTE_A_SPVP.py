# FUNCTION - BUILD INFLUENCE MATRIX A FOR PANEL VORTEX/PANEL METHOD INCLUDING KUTTA CONDITION

import numpy as np
import math as math

def COMPUTE_A_SPVP(numPan,I,K,J,L):
    # Populate A matrix
    A = np.zeros([numPan,numPan])                                                   # Initialize the A matrix
    for i in range(numPan):                                                         # Loop over all i panels
        for j in range(numPan):                                                     # Loop over all j panels
            if (i == j):                                                            # If the panels are the same
                A[i,j] = np.pi                                                      # Set A equal to pi
            else:                                                                   # If panels are not the same
                A[i,j] = I[i,j]                                                     # Set A equal to I

    # Right column of A matrix
    newAV = np.zeros((numPan,1))                                                    # Used to enlarge the A matrix to account for gamma column
    A     = np.hstack((A,newAV))                                                    # Horizontally stack the A matrix with newAV to get enlarged matrix
    for i in range(numPan):                                                         # Loop over all i panels (rows)
        A[i,numPan] = -sum(K[i,:])                                                  # Add gamma term to right-most column of A matrix

    # Bottom row of A matrix
    newAH = np.zeros((1,numPan+1))                                                  # Used to enlarge the A matrix to account for Kutta condition equation
    A     = np.vstack((A,newAH))                                                    # Vertically stack the A matrix with newAH to get enlarged matrix
    for j in range(numPan):                                                         # Loop over all j panels (columns)
        A[numPan,j] = J[0,j] + J[numPan-1,j]                                        # Source contribution of Kutta condition equation
    A[numPan,numPan] = -(sum(L[0,:] + L[numPan-1,:])) + 2*np.pi                     # Vortex contribution of Kutta condition equation 
    return A