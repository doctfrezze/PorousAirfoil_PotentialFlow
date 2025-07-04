import numpy as np
import math as math
import matplotlib.pyplot as plt


from SPVP_Airfoil import SPVP
from Hydraulic_Resistance import Hydraulic_Resistance
from PLOT import PLOT
from Hydraulic_GEOMETRY import Hydraulic_GEOMETRY
from Hydraulic_GEOMETRY import Refine_GEOMETRY
from GEOMETRY import GEOMETRY
from NACA import GENERATE_NACA4


if __name__ == "__main__":
    # User-defined knowns
    Vinf = 1                                                                        # Freestream velocity [] (just leave this at 1)
    rhoinf = 1                                                                      # Density [] (just leave this at 1)
    Re = 1000                                                                       # Reynolds number
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
    L = 0.6
    a = 0.025
    n = 1/0.167
    entry_point = 45
    out_point = 20

    #Convergence Variables
    max_iter = 100
    tol = 1e-8
    err = 100
    iter = 1


    

    mu = Vinf*rhoinf/Re
    Rs = Hydraulic_Resistance(mu,L,type,pore_geometry)



    numPan = 100
    power = 1
    NameAirfoil = "0018"


    AoAR = AoA*(np.pi/180)                                                          # Angle of attack [rad]
    XB, YB = GENERATE_NACA4(NameAirfoil,n=int(numPan/2+1),power=power)
    XC,YC,S,phi,delta,beta = GEOMETRY(numPan,XB,YB,AoAR)
    XB,YB,XC,YC,S,phi,delta,beta,entry_point,out_point = Refine_GEOMETRY(XB,NameAirfoil,entry_point,out_point,AoAR)
    numPan = len(XB)-1



            

    pore_entry = Hydraulic_GEOMETRY(XC,YC,delta,a,entry_point)
    pore_exit = Hydraulic_GEOMETRY(XC,YC,delta,a,out_point)

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
    plt.show()











    #First round without porous

    XC,YC,S,delta,Cp,phi,lam,gamma,CL,CM = SPVP(XB,YB,numPan,Vinf,AoA,rhoinf,0,Rs,a,n,pore_entry,pore_exit)
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
        XC,YC,S,delta,Cp,phi,lam,gamma,CL,CM = SPVP(XB,YB,numPan,Vinf,AoA,rhoinf,Delta_Cp,Rs,a,n,pore_entry,pore_exit)
        err = abs((Delta_Cp-(Cp[entry_point]-Cp[out_point]))/Delta_Cp)
        iter += 1
        if iter == max_iter:
            print("Maximum number of iterations reached")
        if err < tol:
            print("Convergence reached in ", iter, ' iterations')
    
    # %% PLOT Result
    print("CL Porous = ",CL)
    print("CL Solid = ",CL_Solid)
    fig = plt.figure(4)                                                         # Create figure
    plt.cla()                                                                   # Get ready for plotting
    midIndS = int(np.floor(len(Cp)/2))                                          # Airfoil middle index for VPM data
    plt.plot(XC[midIndS+1:len(XC)],Cp[midIndS+1:len(XC)],                       # Plot Cp for upper surface of airfoil from panel method
                'ks',markerfacecolor='b',label='Porous Upper')
    plt.plot(XC[0:midIndS],Cp[0:midIndS],                                       # Plot Cp for lower surface of airfoil from panel method
                'ks',markerfacecolor='r',label='Porous Lower')
    plt.plot(XC[midIndS+1:len(XC)],Cp_Solid[midIndS+1:len(XC)],                       # Plot Cp for upper surface of airfoil from panel method
                color='b',label='Solid Upper')
    plt.plot(XC[0:midIndS],Cp_Solid[0:midIndS],                                       # Plot Cp for lower surface of airfoil from panel method
                color='r',label='Solid Lower')
    plt.plot(XC[entry_point],Cp[entry_point],'ks')
    plt.xlim(0,1)                                                               # Set X-limits
    plt.xlabel('X Coordinate')                                                  # Set X-label
    plt.ylabel('Cp')                                                            # Set Y-label
    plt.title('Pressure Coefficient')                                           # Set title
                                                                        # Display plot
    plt.legend()                                                                # Display legend
    plt.gca().invert_yaxis()                                                    # Invert Cp (Y) axis
    plt.show()
    print(len(phi))
    PLOT(flagPlot,XB,YB,numPan,XC,YC,S,delta,Cp,phi,Vinf,AoA,lam,gamma)