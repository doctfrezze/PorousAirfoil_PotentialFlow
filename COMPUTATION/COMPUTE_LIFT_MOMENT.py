# FUNCTION - COMPUTE LIFT (CL) AND MOMENT (CM) COEFFICIENTS FROM PRESSURE DISTRIBUTION IN SOURCE PANEL METHOD

import numpy as np
import math as math
def COMPUTE_LIFT_MOMENT(Cp,Fluid_characteristics,Airfoil_geometry,Pore_characteristics, Delta_Cp = 0,Cp_inter_low=0,Cp_inter_high=0,A=0,is_porous = 0):
    S = Airfoil_geometry['S']
    beta = Airfoil_geometry['beta']
    phi = Airfoil_geometry['phi']
    XC = Airfoil_geometry['XC']
    YC = Airfoil_geometry['YC']
    if not is_porous:
        CL = sum(-Cp*S*np.sin(beta))
        CD = sum(-Cp*S*np.cos(beta))
        CM = sum(Cp*(XC-0.25)*S*np.cos(phi) + Cp*YC*S*np.sin(phi))                           # Moment coefficient []
    # Compute normal and axial force coefficients
    else:
        low_point = Pore_characteristics['low_point']
        high_point = Pore_characteristics['high_point']
        S_pore_low = Pore_characteristics['phi_pore_low']
        S_pore_high = Pore_characteristics['phi_pore_high']
        phi_pore_low = Pore_characteristics['phi_pore_low']
        phi_pore_high = Pore_characteristics['phi_pore_high']
        AoAR = Fluid_characteristics['AoAR']

        CL = sum(-Cp[i]*S[i]*np.sin(beta[i]) for i in low_point)
        print("CL = ",CL)                             # Normal force coefficient []
        CL += sum(-Cp[i]*S[i]*np.sin(beta[i]) for i in high_point)
        print("CL = ",CL)
        CL += sum(-Cp_inter_low*S_pore_low*np.sin(phi_pore_low+(np.pi/2)-AoAR))
        print("CL = ",CL)
        CL += sum(-Cp_inter_high*S_pore_high*np.sin(phi_pore_high+(np.pi/2)-AoAR))
        print("CL = ",CL)
        CL += -2*Delta_Cp*A*np.sin(phi_pore_high[0])
        print("CL = ",CL)
        print('.  ')

        CD = sum(-Cp[i]*S[i]*np.cos(beta[i]) for i in low_point)                             
        CD += sum(-Cp[i]*S[i]*np.cos(beta[i]) for i in high_point)
        CD += sum(-Cp_inter_low*S_pore_low*np.cos(phi_pore_low+(np.pi/2)-AoAR))
        CD += sum(-Cp_inter_high*S_pore_high*np.cos(phi_pore_high+(np.pi/2)-AoAR))
        print('CD = ',CD)
        CD += -2*Delta_Cp*A*np.cos(phi_pore_high[0])
        print('CD = ',CD)
        print('. ')

        CM = sum(Cp[i]*(XC[i]-0.25)*S[i]*np.cos(phi[i]) + Cp[i]*YC[i]*S[i]*np.sin(phi[i]) for i in low_point)           
        CM += sum(Cp[i]*(XC[i]-0.25)*S[i]*np.cos(phi[i]) + Cp[i]*YC[i]*S[i]*np.sin(phi[i]) for i in high_point)                 

    return CL,CM,CD