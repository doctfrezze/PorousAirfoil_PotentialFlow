import numpy as np
import math as math
import matplotlib.pyplot as plt


from SPVP_Airfoil import SPVP
from COMPUTATION.Hydraulic_Resistance import Hydraulic_Resistance
from PLOT import *
from GEOMETRY.Hydraulic_GEOMETRY import *
from GEOMETRY.GEOMETRY import GENERATE_GEOMETRY
from GEOMETRY.NACA import GENERATE_NACA4
from COMPUTATION.COMPUTE import COMPUTE_LIFT_MOMENT

def INIT_POROUS_GEOMETRY(AoA,NameAirfoil,numPan,omega_in,omega_out,out_point,entry_point,a,power=1,is_straight=1):
    AoAR = AoA*(np.pi/180)                                                          # Angle of attack [rad]
    XB, YB = GENERATE_NACA4(NameAirfoil,n=int(numPan/2+1),power=power)
    XC,YC,S,phi,delta,beta = GENERATE_GEOMETRY(numPan,XB,YB,AoAR)
    if is_straight:
        omega_in = math.atan2(YC[out_point]-YC[entry_point],XC[out_point]-XC[entry_point])
        omega_out = math.atan2(YC[out_point]-YC[entry_point],XC[out_point]-XC[entry_point])
    XB,YB,XC,YC,S,phi,delta,beta,entry_point,out_point = Refine_GEOMETRY(XB,NameAirfoil,entry_point,out_point,AoAR)
    numPan = len(XB)-1

    pore_entry = Hydraulic_GEOMETRY(XC,YC,omega_in,a,entry_point)
    pore_exit = Hydraulic_GEOMETRY(XC,YC,omega_out,a,out_point)

    low_point,high_point = pressure_succion_side(numPan,pore_entry,pore_exit)
    
    pore_intern_co_XB_low = np.linspace(XB[low_point[-1]+1],XB[low_point[0]],len(low_point)+1)
    pore_intern_co_YB_low = np.linspace(YB[low_point[-1]+1],YB[low_point[0]],len(low_point)+1)

    pore_intern_co_XB_high = np.linspace(XB[high_point[-1]+1],XB[high_point[0]],len(high_point)+1)
    pore_intern_co_YB_high = np.linspace(YB[high_point[-1]+1],YB[high_point[0]],len(high_point)+1)

    S_pore_low,phi_pore_low, pore_intern_co_XC_low, pore_intern_co_YC_low = Pore_Geometry(pore_intern_co_XB_low,pore_intern_co_YB_low)
    S_pore_high,phi_pore_high, pore_intern_co_XC_high, pore_intern_co_YC_high = Pore_Geometry(pore_intern_co_XB_high,pore_intern_co_YB_high)

    return XB,YB,XC,YC,S,phi,delta,beta,entry_point,out_point,numPan,pore_entry,pore_exit,omega_in,omega_out,low_point,high_point,pore_intern_co_XB_low,pore_intern_co_YB_low,pore_intern_co_XB_high,pore_intern_co_YB_high,S_pore_low,phi_pore_low, pore_intern_co_XC_low, pore_intern_co_YC_low,S_pore_high,phi_pore_high, pore_intern_co_XC_high, pore_intern_co_YC_high

def POROUS_SPVP(tol,max_iter,Pore_characteristics,Fluid_characteristics,Airfoil_geometry):
    
    entry_point = Pore_characteristics['entry_point']
    out_point = Pore_characteristics['out_point']
    pore_intern_co_XC_low = Pore_characteristics['pore_intern_co_XC_low']
    pore_intern_co_XC_high = Pore_characteristics['pore_intern_co_XC_high']
    A = Pore_characteristics['A']
    low_point = Pore_characteristics['low_point']
    high_point = Pore_characteristics['high_point']

    XB = Airfoil_geometry['XB']

    #%% First round without porous

    Cp,lam,gamma,CL,CM,CD = SPVP(Fluid_characteristics,Airfoil_geometry,is_porous = 0)
    Cp_Solid = Cp
    CL_Solid = CL
    CM_Solid = CM
    CD_Solid = CD
    err = 100000
    iter =0
    #%% Loop with porous
    while err > tol and iter < max_iter: 
        Delta_Cp = Cp[entry_point]-Cp[out_point]
        Cp,lam,gamma,CL,CM,CD = SPVP(Fluid_characteristics,Airfoil_geometry,Pore_characteristics,is_porous = 1,Delta_Cp=Delta_Cp, low_point= low_point, high_point = high_point)
        err = abs((Delta_Cp-(Cp[entry_point]-Cp[out_point]))/Delta_Cp)
        print("Delta_cp = ", Delta_Cp, '     err = ', err)
        iter += 1
        if iter == max_iter:
            print("Maximum number of iterations reached")
        if err < tol:
            print("Convergence reached in ", iter, ' iterations')
    
    # %% CALCULATION OF FINAL CP
    Cp_inter_low = np.zeros(len(pore_intern_co_XC_low))
    for i in range(len(pore_intern_co_XC_low)):
        Cp_inter_low[i] = Cp[entry_point] - Delta_Cp*(pore_intern_co_XC_low[i]-XB[entry_point])/(XB[out_point]-XB[entry_point])
    
    Cp_inter_high = np.zeros(len(pore_intern_co_XC_high))
    for i in range(len(pore_intern_co_XC_high)):
        Cp_inter_high[i] = Cp[entry_point] - Delta_Cp*(pore_intern_co_XC_high[i]-XB[entry_point])/(XB[out_point]-XB[entry_point])
    CL,CM,CD = COMPUTE_LIFT_MOMENT(Cp,Fluid_characteristics,Airfoil_geometry,Pore_characteristics,Cp_inter_low,Cp_inter_high,Delta_Cp,A,is_porous=1)
    return Cp, Cp_Solid, Cp_inter_low,Cp_inter_high, CL, CL_Solid, CD, CD_Solid,lam,gamma

if __name__ == "__main__":
    
    #%% User-defined knowns
    Vinf = 1                                                                        # Freestream velocity [] (just leave this at 1)
    rhoinf = 1                                                                      # Density [] (just leave this at 1)
    Re = 160000                                                                      # Reynolds number
    AoA  = 1                                                                     # Angle of attack [deg]
    
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
    n = 1/0.166
    entry_point = 36
    out_point = 5
    omega_in = 0
    omega_out = 0

    #Convergence Variables
    max_iter = 100
    tol = 1e-8
    err = 100

    #%% Initialisation 
    # Fluid characteristics
    AoAR = AoA*np.pi/180
    mu = Vinf*rhoinf/Re
    #Pores characteristics
    Rs,Dh,A = Hydraulic_Resistance(mu,L,type,pore_geometry)


    AoAR = AoA*np.pi/180
    XB,YB,XC,YC,S,phi,delta,beta,entry_point,out_point,numPan,pore_entry,pore_exit,omega_in,omega_out,low_point,high_point,pore_intern_co_XB_low,pore_intern_co_YB_low,pore_intern_co_XB_high,pore_intern_co_YB_high,S_pore_low,phi_pore_low, pore_intern_co_XC_low, pore_intern_co_YC_low,S_pore_high,phi_pore_high, pore_intern_co_XC_high, pore_intern_co_YC_high = INIT_POROUS_GEOMETRY(AoA,NameAirfoil,numPan,omega_in,omega_out,out_point,entry_point,a,power=power,is_straight=1)
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
        'Dh' : Dh,
        'A' : A,
        'pore_intern_co_XB_low' : pore_intern_co_XB_low,
        'pore_intern_co_YB_low' : pore_intern_co_YB_low,
        'pore_intern_co_XB_high' : pore_intern_co_XB_high,
        'pore_intern_co_YB_high' : pore_intern_co_YB_high,
        'S_pore_low' : S_pore_low,
        'phi_pore_low' : phi_pore_low,
        'pore_intern_co_XC_low' : pore_intern_co_XC_low,
        'pore_intern_co_YC_low' : pore_intern_co_YC_low,
        'S_pore_high' : S_pore_high,
        'phi_pore_high' : phi_pore_high,
        'pore_intern_co_XC_high' : pore_intern_co_XC_high,
        'pore_intern_co_YC_high' : pore_intern_co_YC_high,
        'low_point' : low_point,
        'high_point' : high_point
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

    # %% CALCULATION
    Cp, Cp_Solid, Cp_inter_low,Cp_inter_high, CL, CL_Solid, CD, CD_Solid,lam,gamma = POROUS_SPVP(tol,max_iter,Pore_characteristics,Fluid_characteristics,Airfoil_geometry)



    # %% PLOT Result
    PLOT_AIRFOIL(XB,YB,low_point,high_point,alone=0)
    PLOT_CP_COMPARISON(XB,XC,Cp,Cp_Solid,pore_entry,pore_exit,label1='Porous',label2='Solid',alone = False)
    PLOT_CP_PRESSURE_SIDE(XC,YC, Cp, Cp_inter_low, low_point, pore_intern_co_XC_low, alone = False)
    PLOT_CP_SUCCION_SIDE(XC,YC, Cp, Cp_inter_high, high_point, pore_intern_co_XC_high, alone = False)
    
    print('CL_Porous = ', CL)
    print('CL_solid = ', CL_Solid)
    PLOT_ALL(flagPlot,XB,YB,numPan,XC,YC,S,delta,Cp,phi,Vinf,AoA,lam,gamma)



