# SOURCE/VORTEX PANEL METHOD

import numpy as np
import math as math

from COMPUTATION.COMPUTE_IJ_SPM import COMPUTE_IJ_SPM
from COMPUTATION.COMPUTE_KL_VPM import COMPUTE_KL_VPM
from COMPUTATION.COMPUTE_A_SPVP import COMPUTE_A_SPVP
from COMPUTATION.COMPUTE_b_SPVP import COMPUTE_b_SPVP
from COMPUTATION.COMPUTE_SOLUTION_SPVP import COMPUTE_SOLUTION_SPVP
from COMPUTATION.COMPUTE_Vt_Cp import COMPUTE_Vt_Cp
from COMPUTATION.COMPUTE_LIFT_MOMENT import COMPUTE_LIFT_MOMENT
from PLOT import PLOT
from NACA import GENERATE_NACA4
from GEOMETRY import GEOMETRY
from PANEL_DIRECTIONS import PANEL_DIRECTIONS

def SPVP(XB, YB,numPan,Vinf,AoA,rhoinf=1,Delta_Cp=0,Rs=1,a=1,n=1,pore_entry=[],pore_out=[]):
    # Convert angle of attack to radians
    AoAR = AoA*(np.pi/180)                                                          # Angle of attack [rad]
    

    # CHECK PANEL DIRECTIONS - FLIP IF NECESSARY
    PANEL_DIRECTIONS(numPan,XB,YB)

    # PANEL METHOD GEOMETRY - REF [1]
    XC,YC,S,phi,delta,beta = GEOMETRY(numPan,XB,YB,AoAR)

    # COMPUTE SOURCE AND VORTEX PANEL STRENGTHS - REF [10]
    I, J = COMPUTE_IJ_SPM(XC,YC,XB,YB,phi,S)                                        # Call COMPUTE_IJ_SPM function (Refs [2] and [3])
    K, L = COMPUTE_KL_VPM(XC,YC,XB,YB,phi,S)                                        # Call COMPUTE_KL_VPM function (Refs [6] and [7])

    A = COMPUTE_A_SPVP(numPan,I,K,J,L)
    b = COMPUTE_b_SPVP(numPan,Vinf,beta,rhoinf,Delta_Cp,Rs,a,n,pore_entry,pore_out)

    lam, gamma = COMPUTE_SOLUTION_SPVP(A,b)

    # COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS
    Vt, Cp = COMPUTE_Vt_Cp(numPan,Vinf,beta,lam,gamma,J,L)

    # COMPUTE LIFT AND MOMENT COEFFICIENTS
    CL,CM = COMPUTE_LIFT_MOMENT(Cp,S,beta,phi,AoAR,XC,YC)
    return XC,YC,S,delta,Cp,phi,lam,gamma,CL,CM


if __name__ == '__main__':
    # User-defined knowns
    Vinf = 1                                                                        # Freestream velocity [] (just leave this at 1)
    AoA  = 0                                                                        # Angle of attack [deg]
    # Plotting flags
    flagPlot = [0,      # Airfoil with panel normal vectors
                0,      # Geometry boundary pts, control pts, first panel, second panel
                1,      # Cp vectors at airfoil surface panels
                1,      # Pressure coefficient comparison (XFOIL vs. VPM)
                1,      # Airfoil streamlines
                1]      # Pressure coefficient contour

    # AirFoil panels
    numPan = 100
    NameAirfoil = "0018"
    XB, YB = GENERATE_NACA4(NameAirfoil,n=int(numPan/2+1),power=3)
    XC,YC,S,delta,Cp,phi,lam,gamma,CL,CM = SPVP(XB, YB,numPan,Vinf,AoA)
    # %% PLOT
    PLOT(flagPlot,XB,YB,numPan,XC,YC,S,delta,Cp,phi,Vinf,AoA,lam,gamma)


    #%% RESULTS
    # Print the results to the Console
    print("======= RESULTS =======")
    print("Lift Coefficient (CL)")
    print("  SPVP : %2.8f" % CL)                                                    # From this SPVP code
    print("Moment Coefficient (CM)")
    print("  SPVP : %2.8f" % CM)                                                    # From this SPVP code