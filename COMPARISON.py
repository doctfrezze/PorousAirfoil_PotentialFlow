import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib import path

from SPVP_Airfoil import SPVP
from X_FOIL.XFOIL import XFOIL_DATA
from GEOMETRY.GEOMETRY import GENERATE_GEOMETRY
from GEOMETRY.NACA import GENERATE_NACA4
from COMPUTATION.COMPUTE import COMPUTE_LIFT_MOMENT

Vinf = 1               # Freestream velocity [unitless, normalized]
rhoinf = 1             # Freestream density
Re = 160000              # Reynolds number
mu = Vinf * rhoinf / Re  # Dynamic viscosity
AoA  = 4               # Angle of attack [degrees]
AoAR = np.pi * AoA / 180  # Angle of attack [radians]


# Airfoil parameters
numPan = 100             # Number of panels (discretization segments)
NameAirfoil = "0018"     # NACA 4-digit airfoil code
power = 3                # Point spacing exponent for clustering near leading/trailing edges

#%% ================================
#          Geometry initialization
#==================================
# Generate boundary coordinates for the NACA airfoil
XB, YB = GENERATE_NACA4(NameAirfoil, n=int(numPan / 2 + 1), power=power)

# Compute panel geometry based on boundary coordinates and angle of attack
XC, YC, S, phi, delta, beta = GENERATE_GEOMETRY(numPan, XB, YB, AoAR)

#%% ================================
#       Dictionary-based storage
#==================================
# Store fluid parameters
Fluid_characteristics = {
    'Vinf' : Vinf,
    'rhoinf' : rhoinf,
    'AoA' : AoA,
    'AoAR' : AoAR,
    'mu' : mu
}

# Store airfoil geometry data
Airfoil_geometry = {
    'XB' : XB,
    'YB' : YB,
    'XC' : XC,
    'YC' : YC,
    'S' : S,
    'phi' : phi,
    'delta' : delta,
    'beta' : beta,
    'numPan' : numPan,
    'power' : power,
    'NameAirfoil' : NameAirfoil
}

#%% ================================
#        SPVP Main Calculation
#==================================
# Compute aerodynamic results using the SPVP method
Cp, lam, gamma, CL, CM, CD = SPVP(
    Fluid_characteristics,
    Airfoil_geometry,
    is_porous=0  # 0 = non-porous airfoil
)
x_xfoil,cp_xfoil = XFOIL_DATA(NameAirfoil,4)

t = 0.24
c=1
y_xfoil = 5 * t * c * (0.2969 * np.sqrt(x_xfoil / c) - 0.1260 * (x_xfoil / c) - 0.3516 * (x_xfoil / c)**2 
                      + 0.2843 * (x_xfoil / c)**3 - 0.1015 * (x_xfoil / c)**4)

XC_xfoil,YC_xfoil,S_xfoil,phi_xfoil,delta_xfoil,beta_xfoil = GENERATE_GEOMETRY(len(x_xfoil)-1,x_xfoil,y_xfoil,AoAR)



#%% PLOT
fig = plt.figure()
fig.suptitle("NACA"+NameAirfoil)
midIndS_xfoil = int(np.floor(len(cp_xfoil)/2))                                          # Airfoil middle index for XFOIL data
plt.plot(x_xfoil[midIndS_xfoil:len(x_xfoil)], cp_xfoil[midIndS_xfoil:len(x_xfoil)],
         color='r',label='XFOIL Lower')
plt.plot(x_xfoil[0:midIndS_xfoil], cp_xfoil[0:midIndS_xfoil],
         color='b',label='XFOIL Upper')

midIndS = int(np.floor(len(Cp)/2))                                          # Airfoil middle index for SPVP data
plt.plot(XC[midIndS+1:len(XC)],Cp[midIndS+1:len(XC)],                       # Plot Cp for upper surface of airfoil from panel method
            'ks',markerfacecolor='b',label='SPVP Upper')
plt.plot(XC[0:midIndS],Cp[0:midIndS],                                       # Plot Cp for lower surface of airfoil from panel method
            'ks',markerfacecolor='r',label='SPVP Lower')

plt.gca().invert_yaxis()  # Cp plots ont souvent l'axe Y inversé
plt.xlabel("x/c")
plt.ylabel("Cp")
plt.title("Cp at α = 4°")
plt.grid(True)
plt.legend()
plt.show()

