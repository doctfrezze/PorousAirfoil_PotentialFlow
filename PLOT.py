# FUNCTION - GENERATE VARIOUS PLOTS FOR SPVP METHOD INCLUDING GEOMETRY, PRESSURE COEFFICIENTS, STREAMLINES, AND CONTOURS

import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib import path

from COMPUTATION.STREAMLINE_SPM import STREAMLINE_SPM
from COMPUTATION.STREAMLINE_VPM import STREAMLINE_VPM



def PLOT(flagPlot,XB,YB,numPan,XC,YC,S,delta,Cp,phi,Vinf,AoA,lam,gamma):
    AoAR = AoA*(np.pi/180)  
    if (flagPlot[4] == 1 or flagPlot[5] == 1):                                      # If we are plotting streamlines or pressure coefficient contours
        # Grid parameters
        nGridX = 100                                                                # X-grid for streamlines and contours
        nGridY = 100                                                                # Y-grid for streamlines and contours
        xVals  = [min(XB)-0.5, max(XB)+0.5]                                         # X-grid extents [min, max]
        yVals  = [min(YB)-0.3, max(YB)+0.3]                                         # Y-grid extents [min, max]
        
        # Streamline parameters
        slPct  = 25                                                                 # Percentage of streamlines of the grid
        Ysl    = np.linspace(yVals[0],yVals[1],int((slPct/100)*nGridY))             # Create array of Y streamline starting points
        Xsl    = xVals[0]*np.ones(len(Ysl))                                         # Create array of X streamline starting points
        XYsl   = np.vstack((Xsl.T,Ysl.T)).T                                         # Concatenate X and Y streamline starting points
        
        # Generate the grid points
        Xgrid  = np.linspace(xVals[0],xVals[1],nGridX)                              # X-values in evenly spaced grid
        Ygrid  = np.linspace(yVals[0],yVals[1],nGridY)                              # Y-values in evenly spaced grid
        XX, YY = np.meshgrid(Xgrid,Ygrid)                                           # Create meshgrid from X and Y grid arrays
        
        # Initialize velocities
        Vx     = np.zeros([nGridX,nGridY])                                          # Initialize X velocity matrix
        Vy     = np.zeros([nGridX,nGridY])                                          # Initialize Y velocity matrix
        
        # Path to figure out if grid point is inside polygon or not
        AF     = np.vstack((XB.T,YB.T)).T                                           # Concatenate XB and YB geometry points
        afPath = path.Path(AF)                                                      # Create a path for the geometry
        
        # Solve for grid point X and Y velocities
        for m in range(nGridX):                                                     # Loop over X-grid points
            print("m: %i" % m)
            for n in range(nGridY):                                                 # Loop over Y-grid points
                XP     = XX[m,n]                                                    # Current iteration's X grid point
                YP     = YY[m,n]                                                    # Current iteration's Y grid point
                Mx, My = STREAMLINE_SPM(XP,YP,XB,YB,phi,S)                          # Compute streamline Mx and My values
                Nx, Ny = STREAMLINE_VPM(XP,YP,XB,YB,phi,S)                          # Compute streamline Nx and Ny values
                
                # Check if grid points are in object
                # - If they are, assign a velocity of zero
                if afPath.contains_points([(XP,YP)]):                               # If (XP,YP) is in the body
                    Vx[m,n] = 0                                                     # Set X-velocity equal to zero
                    Vy[m,n] = 0                                                     # Set Y-velocity equal to zero
                else:
                    Vx[m,n] = (Vinf*np.cos(AoAR) + sum(lam*Mx/(2*np.pi))            # Compute X-velocity
                                                + sum(-gamma*Nx/(2*np.pi)))
                    Vy[m,n] = (Vinf*np.sin(AoAR) + sum(lam*My/(2*np.pi))            # Compute Y-velocity
                                                + sum(-gamma*Ny/(2*np.pi)))
        
        # Compute grid point velocity magnitude and pressure coefficient
        Vxy  = np.sqrt(Vx**2 + Vy**2)                                               # Compute magnitude of velocity vector []
        CpXY = 1 - (Vxy/Vinf)**2                                                    # Pressure coefficient []

    if (flagPlot[0] == 1):
        fig = plt.figure(1)                                                         # Create the figure
        plt.cla()                                                                   # Clear the axes
        plt.fill(XB,YB,'k')                                                         # Plot the airfoil
        X = np.zeros(2)                                                             # Initialize 'X'
        Y = np.zeros(2)                                                             # Initialize 'Y'
        for i in range(numPan):                                                     # Loop over all panels
            X[0] = XC[i]                                                            # Set X start of panel orientation vector
            X[1] = XC[i] + S[i]*np.cos(delta[i])                                    # Set X end of panel orientation vector
            Y[0] = YC[i]                                                            # Set Y start of panel orientation vector
            Y[1] = YC[i] + S[i]*np.sin(delta[i])                                    # Set Y end of panel orientation vector
            if (i == 0):                                                            # If it's the first panel index
                plt.plot(X,Y,'b-',label='First Panel')                              # Plot normal vector for first panel
            elif (i == 1):                                                          # If it's the second panel index
                plt.plot(X,Y,'g-',label='Second Panel')                             # Plot normal vector for second panel
            else:                                                                   # If it's neither the first nor second panel index
                plt.plot(X,Y,'r-')                                                  # Plot normal vector for all other panels
        plt.xlabel('X Units')                                                       # Set X-label
        plt.ylabel('Y Units')                                                       # Set Y-label
        plt.title('Panel Geometry')                                                 # Set title
        plt.axis('equal')                                                           # Set axes equal
        plt.legend()                                                                # Display legend

    # FIGURE: Geometry with the following indicated:
    # - Boundary points, control points, first panel, second panel
    if (flagPlot[1] == 1):
        fig = plt.figure(2)                                                         # Create figure
        plt.cla()                                                                   # Get ready for plotting
        plt.plot(XB,YB,'k-')                                                        # Plot airfoil panels
        plt.plot([XB[0], XB[1]],[YB[0], YB[1]],'b-',label='First Panel')            # Plot first panel
        plt.plot([XB[1], XB[2]],[YB[1], YB[2]],'g-',label='Second Panel')           # Plot second panel
        plt.plot(XB,YB,'ko',markerfacecolor='k',label='Boundary Pts')               # Plot boundary points (black circles)
        plt.plot(XC,YC,'ko',markerfacecolor='r',label='Control Pts')                # Plot control points (red circles)
        plt.xlabel('X Units')                                                       # Set X-label
        plt.ylabel('Y Units')                                                       # Set Y-label
        plt.axis('equal')                                                           # Set axes equal
        plt.legend()                                                                # Display legend

    # FIGURE: Cp vectors at airfoil control points
    if (flagPlot[2] == 1):
        fig = plt.figure(3)                                                         # Create figure
        plt.cla()                                                                   # Get ready for plotting
        Cps = np.absolute(Cp*0.15)                                                  # Scale and make positive all Cp values
        X = np.zeros(2)                                                             # Initialize X values
        Y = np.zeros(2)                                                             # Initialize Y values
        for i in range(len(Cps)):                                                   # Loop over all panels
            X[0] = XC[i]                                                            # Control point X-coordinate
            X[1] = XC[i] + Cps[i]*np.cos(delta[i])                                  # Ending X-value based on Cp magnitude
            Y[0] = YC[i]                                                            # Control point Y-coordinate
            Y[1] = YC[i] + Cps[i]*np.sin(delta[i])                                  # Ending Y-value based on Cp magnitude
            
            if (Cp[i] < 0):                                                         # If pressure coefficient is negative
                plt.plot(X,Y,'r-')                                                  # Plot as a red line
            elif (Cp[i] >= 0):                                                      # If pressure coefficient is zero or positive
                plt.plot(X,Y,'b-')                                                  # Plot as a blue line
        plt.fill(XB,YB,'k')                                                         # Plot the airfoil as black polygon
        plt.xlabel('X Units')                                                       # Set X-label
        plt.ylabel('Y Units')                                                       # Set Y-label
        plt.gca().set_aspect('equal')                                               # Set aspect ratio equal

    # FIGURE: Pressure coefficient
    if (flagPlot[3] == 1):
        fig = plt.figure(4)                                                         # Create figure
        plt.cla()                                                                   # Get ready for plotting
        midIndS = int(np.floor(len(Cp)/2))                                          # Airfoil middle index for VPM data
        plt.plot(XC[midIndS+1:len(XC)],Cp[midIndS+1:len(XC)],                       # Plot Cp for upper surface of airfoil from panel method
                    'ks',markerfacecolor='b',label='VPM Upper')
        plt.plot(XC[0:midIndS],Cp[0:midIndS],                                       # Plot Cp for lower surface of airfoil from panel method
                    'ks',markerfacecolor='r',label='VPM Lower')
        plt.xlim(0,1)                                                               # Set X-limits
        plt.xlabel('X Coordinate')                                                  # Set X-label
        plt.ylabel('Cp')                                                            # Set Y-label
        plt.title('Pressure Coefficient')                                           # Set title
                                                                            # Display plot
        plt.legend()                                                                # Display legend
        plt.gca().invert_yaxis()                                                    # Invert Cp (Y) axis

    # FIGURE: Airfoil streamlines
    if (flagPlot[4] == 1):
        fig = plt.figure(5)                                                         # Create figure
        plt.cla()                                                                   # Get ready for plotting
        np.seterr(under="ignore")                                                   # Ignore underflow error message
        plt.streamplot(XX,YY,Vx,Vy, linewidth=0.5, density=40, color='r',           # Plot streamlines
                        arrowstyle='-', start_points=XYsl)
        plt.clim(vmin=0, vmax=2)
        plt.fill(XB,YB,'k')                                                         # Plot airfoil as black polygon
        plt.xlabel('X Units')                                                       # Set X-label
        plt.ylabel('Y Units')                                                       # Set Y-label
        plt.gca().set_aspect('equal')                                               # Set axes equal
        plt.xlim(xVals)                                                             # Set X-limits
        plt.ylim(yVals)                                                             # Set Y-limits

    # FIGURE: Pressure coefficient contour
    if (flagPlot[5] == 1):
        fig = plt.figure(6)                                                         # Create figure
        plt.cla()                                                                   # Get ready for plotting
        plt.contourf(XX,YY,CpXY,500,cmap='jet')                                     # Plot contour
        plt.fill(XB,YB,'k')                                                         # Plot airfoil as black polygon
        plt.xlabel('X Units')                                                       # Set X-label
        plt.ylabel('Y Units')                                                       # Set Y-label
        plt.gca().set_aspect('equal')                                               # Set axes equal
        plt.xlim(xVals)                                                             # Set X-limits
        plt.ylim(yVals)                                                             # Set Y-limits
    if (flagPlot != [0,0,0,0,0,0]):
        plt.show()                                                                  # Display plots