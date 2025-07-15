import numpy as np
import math as math
import matplotlib.pyplot as plt


from SPVP_Airfoil import SPVP
from Hydraulic_Resistance import Hydraulic_Resistance
from PLOT import PLOT_ALL
from PLOT import PLOT_AIRFOIL
from PLOT import PLOT_CP_COMPARISON
from PLOT import PLOT_CP_PRESSURE_SIDE
from PLOT import PLOT_CP_SUCCION_SIDE
from Hydraulic_GEOMETRY import Hydraulic_GEOMETRY
from Hydraulic_GEOMETRY import Refine_GEOMETRY
from Hydraulic_GEOMETRY import pressure_succion_side
from GEOMETRY import GEOMETRY
from NACA import GENERATE_NACA4

def INIT_POROUS_GEOMETRY(AoA,NameAirfoil,numPan,omega_in,omega_out,out_point,entry_point,power=1,is_straight=1):
    AoAR = AoA*(np.pi/180)                                                          # Angle of attack [rad]
    XB, YB = GENERATE_NACA4(NameAirfoil,n=int(numPan/2+1),power=power)
    XC,YC,S,phi,delta,beta = GEOMETRY(numPan,XB,YB,AoAR)
    if is_straight:
        omega_in = math.atan2(YC[out_point]-YC[entry_point],XC[out_point]-XC[entry_point])
        omega_out = math.atan2(YC[out_point]-YC[entry_point],XC[out_point]-XC[entry_point])
    XB,YB,XC,YC,S,phi,delta,beta,entry_point,out_point = Refine_GEOMETRY(XB,NameAirfoil,entry_point,out_point,AoAR)
    numPan = len(XB)-1

    pore_entry = Hydraulic_GEOMETRY(XC,YC,omega_in,a,entry_point)
    pore_exit = Hydraulic_GEOMETRY(XC,YC,omega_out,a,out_point)

    return XB,YB,XC,YC,S,phi,delta,beta,entry_point,out_point,numPan,pore_entry,pore_exit,omega_in,omega_out

if __name__ == "__main__":
    #%% User-defined knowns
    Vinf = 1                                                                        # Freestream velocity [] (just leave this at 1)
    rhoinf = 1                                                                      # Density [] (just leave this at 1)
    Re = 160000                                                                      # Reynolds number
    AoA  = 4                                                                        # Angle of attack [deg]
    
    numPan = 100
    power = 3
    NameAirfoil = "0018"
    
    # Plotting flags
    flagPlot = [1,      # Airfoil with panel normal vectors
                1,      # Geometry boundary pts, control pts, first panel, second panel
                1,      # Cp vectors at airfoil surface panels
                1,      # Pressure coefficient comparison (XFOIL vs. VPM)
                0,      # Airfoil streamlines
                0]      # Pressure coefficient contour
    
    #Pore geometry
    type = 'rectangle'
    pore_geometry = [0.157,0.025]
    L = 0.89
    a = 0.025                       #Height of the pores
    n = 1/0.157
    entry_point = 36
    out_point = 5
    omega_in = 0
    omega_out = 0

    #Convergence Variables
    max_iter = 100
    tol = 1e-8
    err = 100
    iter = 1

    #%% Initialisation 
    # Fluid characteristics
    AoAR = AoA*np.pi/180
    mu = Vinf*rhoinf/Re

    #Pores characteristics
    Rs,Dh = Hydraulic_Resistance(mu,L,type,pore_geometry)
    
    XB,YB,XC,YC,S,phi,delta,beta,entry_point,out_point,numPan,pore_entry,pore_exit,omega_in,omega_out = INIT_POROUS_GEOMETRY(AoA,NameAirfoil,numPan,omega_in,omega_out,out_point,entry_point,power=power,is_straight=1)


    low_point,high_point = pressure_succion_side(numPan,pore_entry,pore_exit)
    
    pore_intern_co_X_low = [XC[i] for i in low_point]
    pore_intern_co_Y_low = [YC[i] for i in low_point]
    pore_intern_co_X_high = [XC[i] for i in high_point]
    pore_intern_co_Y_high = [YC[i] for i in high_point]
    
    
    #%% Macro variable
    Fluid_characteristics = {
        "Vinf" : Vinf,
        "rhoinf" : rhoinf,
        "Re" : Re,
        "mu" : mu,
        'AoA' : AoA,
        'AoAR' : AoAR
    }

    Pore_characteristics = {
        'type' : type,
        'pore_geometry' : pore_geometry,
        'L' : L,
        'a' : a,
        'n' : n,
        'entry_point' : entry_point,
        'out_point' : out_point,
        'omega_in' : omega_in,
        'omega_out' : omega_out,
        'Rs' : Rs,
        'pore_entry' : pore_entry,
        'pore_exit' : pore_exit,
        'Dh' : Dh
    }

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

    



    #%% First round without porous

    Cp,lam,gamma,CL,CM = SPVP(Fluid_characteristics,Airfoil_geometry,is_porous = 0)
    Cp_Solid = Cp
    CL_Solid = CL
    CM_Solid = CM

    #%% Loop with porous
    while err > tol and iter < max_iter: 
        Delta_Cp = Cp[entry_point]-Cp[out_point]
        print("Delta_cp = ", Delta_Cp)
        Cp,lam,gamma,CL,CM = SPVP(Fluid_characteristics,Airfoil_geometry,Pore_characteristics,is_porous = 1,Delta_Cp=Delta_Cp, low_point= low_point, high_point = high_point)
        err = abs((Delta_Cp-(Cp[entry_point]-Cp[out_point]))/Delta_Cp)
        iter += 1
        if iter == max_iter:
            print("Maximum number of iterations reached")
        if err < tol:
            print("Convergence reached in ", iter, ' iterations')
    
    # %% CALCULATION OF FINAL CP
    Cp_inter_low = []
    for i in range(len(pore_intern_co_X_low)):
        Cp_local = Cp[entry_point] - Delta_Cp*(pore_intern_co_X_low[i]-XB[entry_point])/(XB[out_point]-XB[entry_point])
        Cp_inter_low.append(Cp_local)
    
    Cp_inter_high = []
    for i in range(len(pore_intern_co_X_high)):
        Cp_local = Cp[entry_point] - Delta_Cp*(pore_intern_co_X_high[i]-XB[entry_point])/(XB[out_point]-XB[entry_point])
        Cp_inter_high.append(Cp_local)
    
    # %% PLOT Result
    print("CL Porous = ",CL)
    print("CL Solid = ",CL_Solid)
    PLOT_AIRFOIL(XB,YB,low_point,high_point,alone=0)
    PLOT_CP_COMPARISON(XB,XC,Cp,Cp_Solid,label1='Porous',label2='Solid',alone = False)
    PLOT_CP_PRESSURE_SIDE(XC,YC, Cp, Cp_inter_low, low_point, pore_intern_co_X_low, alone = False)
    PLOT_CP_SUCCION_SIDE(XC,YC, Cp, Cp_inter_high, high_point, pore_intern_co_X_high, alone = False)
    
    PLOT_ALL(flagPlot,XB,YB,numPan,XC,YC,S,delta,Cp,phi,Vinf,AoA,lam,gamma)



