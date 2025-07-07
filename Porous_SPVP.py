import numpy as np
import math as math
import matplotlib.pyplot as plt


from SPVP_Airfoil import SPVP
from Hydraulic_Resistance import Hydraulic_Resistance
from PLOT import PLOT_ALL
from PLOT import PLOT_CP_COMPARISON
from Hydraulic_GEOMETRY import Hydraulic_GEOMETRY
from Hydraulic_GEOMETRY import Refine_GEOMETRY
from GEOMETRY import GEOMETRY
from NACA import GENERATE_NACA4


if __name__ == "__main__":
    # User-defined knowns
    Vinf = 1                                                                        # Freestream velocity [] (just leave this at 1)
    rhoinf = 1                                                                      # Density [] (just leave this at 1)
    Re = 2000                                                                      # Reynolds number
    AoA  = 4                                                                        # Angle of attack [deg]

    # Plotting flags
    flagPlot = [0,      # Airfoil with panel normal vectors
                0,      # Geometry boundary pts, control pts, first panel, second panel
                1,      # Cp vectors at airfoil surface panels
                1,      # Pressure coefficient comparison (XFOIL vs. VPM)
                1,      # Airfoil streamlines
                1]      # Pressure coefficient contour
    
    #Pore geometry
    type = 'rectangle'
    pore_geometry = [0.157,0.025]
    L = 0.89
    a = 0.025
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


    

    mu = Vinf*rhoinf/Re
    Rs = Hydraulic_Resistance(mu,L,type,pore_geometry)



    numPan = 100
    power = 3
    NameAirfoil = "0018"


    AoAR = AoA*(np.pi/180)                                                          # Angle of attack [rad]
    XB, YB = GENERATE_NACA4(NameAirfoil,n=int(numPan/2+1),power=power)
    XC,YC,S,phi,delta,beta = GEOMETRY(numPan,XB,YB,AoAR)
    omega_in = math.atan2(YC[out_point]-YC[entry_point],XC[out_point]-XC[entry_point])
    omega_out = math.atan2(YC[out_point]-YC[entry_point],XC[out_point]-XC[entry_point])
    XB,YB,XC,YC,S,phi,delta,beta,entry_point,out_point = Refine_GEOMETRY(XB,NameAirfoil,entry_point,out_point,AoAR)
    numPan = len(XB)-1


   
            

    pore_entry = Hydraulic_GEOMETRY(XC,YC,omega_in,a,entry_point)
    pore_exit = Hydraulic_GEOMETRY(XC,YC,omega_out,a,out_point)


    fig = plt.figure(1)                                                         # Create the figure
    plt.cla()                                                                   # Clear the axes
    XB_up = np.concatenate([XB[:pore_exit[0]],XB[pore_entry[-1]:]])
    YB_up = np.concatenate([YB[:pore_exit[0]],YB[pore_entry[-1]:]])

    XB_down = np.concatenate([XB[pore_exit[-1]:pore_entry[0]]])
    YB_down = np.concatenate([YB[pore_exit[-1]:pore_entry[0]]])
    plt.fill(XB_up,YB_up,'k')                                                         # Plot the airfoil
    plt.fill(XB_down,YB_down,'k')
    plt.axis('equal')                                                           # Set axes equal
    plt.show()


    X_entry = []
    Y_entry = []
    X_out = []
    Y_out = []
    for i in pore_entry:
        X_entry.append(XC[i])
        Y_entry.append(YC[i])
    for i in pore_exit:
        X_out.append(XC[i])
        Y_out.append(YC[i])
    plt.plot(XB,YB,'ks')
    plt.plot(X_entry,Y_entry,c='red')
    plt.plot(XC[entry_point],YC[entry_point],'rs')
    plt.plot(X_out,Y_out,c='blue')
    plt.plot(XC[out_point],YC[out_point],'gs')
    plt.plot(XB[out_point],YB[out_point],'bs')
    plt.plot(XB[out_point+1],YB[out_point+1],'bs')
    #plt.show()











    #First round without porous

    XC,YC,S,delta,Cp,phi,lam,gamma,CL,CM = SPVP(XB,YB,numPan,Vinf,AoA,rhoinf,0,Rs,a,n,pore_entry,pore_exit,omega_in,omega_out)
    Cp_Solid = Cp
    CL_Solid = CL
    CM_Solid = CM

    #First round with porous
    while err > tol and iter < max_iter: 
        Delta_Cp = Cp[entry_point]-Cp[out_point]
        print("Delta_cp = ", Delta_Cp)
        print('Cp[entry_point] = ', Cp[entry_point])
        print('Cp[out_point] = ', Cp[out_point])
        print("V_moyenne pore = ",0.5*rhoinf*Vinf**2*Delta_Cp/Rs*n/a)
        XC,YC,S,delta,Cp,phi,lam,gamma,CL,CM = SPVP(XB,YB,numPan,Vinf,AoA,rhoinf,Delta_Cp,Rs,a,n,pore_entry,pore_exit,omega_in,omega_out)
        err = abs((Delta_Cp-(Cp[entry_point]-Cp[out_point]))/Delta_Cp)
        iter += 1
        if iter == max_iter:
            print("Maximum number of iterations reached")
        if err < tol:
            print("Convergence reached in ", iter, ' iterations')
    
    # %% PLOT Result
    print("CL Porous = ",CL)
    print("CL Solid = ",CL_Solid)
    PLOT_CP_COMPARISON(XB,XC,Cp,Cp_Solid,label1='Porous',label2='Solid')
    #PLOT_ALL(flagPlot,XB,YB,numPan,XC,YC,S,delta,Cp,phi,Vinf,AoA,lam,gamma)