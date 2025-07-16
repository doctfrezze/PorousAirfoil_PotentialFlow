# FUNCTION - COMPUTE RHS VECTOR b FOR PANEL VORTEX/PANEL METHOD INCLUDING KUTTA CONDITION

import numpy as np
import math as math

def COMPUTE_b_SPVP(Airfoil_geometry,Fluid_characteristics,Delta_Cp,Pore_characteristics={},is_porous = False):
    numPan = Airfoil_geometry['numPan']
    delta = Airfoil_geometry['delta']
    beta = Airfoil_geometry['beta']

    Vinf = Fluid_characteristics['Vinf']
    rhoinf = Fluid_characteristics['rhoinf']
    mu = Fluid_characteristics['mu']
    if is_porous:

        Rs = Pore_characteristics['Rs']
        a = Pore_characteristics['a']
        n = Pore_characteristics['n']
        pore_entry = Pore_characteristics['pore_entry']
        pore_out = Pore_characteristics['pore_exit']
        omega_in = Pore_characteristics['omega_in']
        omega_out = Pore_characteristics['omega_out']
        Dh = Pore_characteristics['Dh']
        L = Pore_characteristics['L']
    


    b = np.zeros(numPan)                                                            # Initialize the b array
    for i in range(numPan):                                                         # Loop over all i panels (rows)
        if is_porous:
            normal_vs_pore_in = delta[i]-np.pi*(omega_in)/180
            normal_vs_pore_out = delta[i]-np.pi*(omega_out)/180
            V_mean = 0.5*rhoinf*Vinf**2*Delta_Cp/Rs*n/a
            Re_laminar = Dh*rhoinf*V_mean/mu
            if Re_laminar >= 2000:                                                      #Turbulence 
                Delta_P = Delta_Cp*0.5*rhoinf*Vinf**2
                V_mean=2.868*((Delta_P**4)*(Dh**5)/((L**4)*mu*(rhoinf**3)))**(1/7)
                print('V_mean = ',V_mean)
            if i in pore_entry:
                b[i] = -Vinf*2*np.pi*np.cos(beta[i]) + V_mean#*np.cos(normal_vs_pore_in) # Compute RHS array
            elif i in pore_out:
                b[i] = -Vinf*2*np.pi*np.cos(beta[i]) + V_mean#*np.cos(normal_vs_pore_out)
            else:
                b[i] = -Vinf*2*np.pi*np.cos(beta[i])        
        else:
            b[i] = -Vinf*2*np.pi*np.cos(beta[i])

    

    # Last element of b array (Kutta condition)
    b = np.append(b,-Vinf*2*np.pi*(np.sin(beta[0]) + np.sin(beta[numPan-1])))       # Add Kutta condition equation RHS to b array
    return b