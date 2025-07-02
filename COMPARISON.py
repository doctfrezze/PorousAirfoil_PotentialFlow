import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib import path

from SPVP_Airfoil import SPVP
from XFOIL import XFOIL_DATA
from GEOMETRY import GEOMETRY
from COMPUTATION.COMPUTE_LIFT_MOMENT import COMPUTE_LIFT_MOMENT

#%% Parameters
NameAirfoil = "2412"
numPan = 300
Vinf = 1
AoA =4
AoAR = AoA*(np.pi/180) 

#%% COMPUTATION AND XFOIL DATA
XB,YB,XC,YC,S,delta,Cp,phi,lam,gamma,CL,CM = SPVP(NameAirfoil,numPan,Vinf,AoA,power=1)
print(CL)
print(CM)
x_xfoil,cp_xfoil = XFOIL_DATA(NameAirfoil,4)

t = 0.24
c=1
y_xfoil = 5 * t * c * (0.2969 * np.sqrt(x_xfoil / c) - 0.1260 * (x_xfoil / c) - 0.3516 * (x_xfoil / c)**2 
                      + 0.2843 * (x_xfoil / c)**3 - 0.1015 * (x_xfoil / c)**4)

XC_xfoil,YC_xfoil,S_xfoil,phi_xfoil,delta_xfoil,beta_xfoil = GEOMETRY(len(x_xfoil)-1,x_xfoil,y_xfoil,AoAR)


CL_xfoil,CM_xfoil = COMPUTE_LIFT_MOMENT(cp_xfoil[:-1],S_xfoil,beta_xfoil,phi_xfoil,AoAR,x_xfoil[:-1],y_xfoil[:-1])

print("CL = ",CL)
print("CM = ",CM)
print("CL_xfoil = ",CL_xfoil)
print("CM_xfoil = ",CM_xfoil)
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


###############  COMPARISAON OF CL AND CM OVER A RANGE OF BETWEEN ANGL########################################

#%% PARAMETERS
# Chemin vers ton fichier
filename = "XFOIL_DATA/"+NameAirfoil+"_CL_CM.dat"

#%% XFOIL DATA LOADING
# Ouvre le fichier pour voir combien de lignes d'en-tête il y a
with open(filename, 'r') as f:
    for i, line in enumerate(f):
        if line.strip().startswith('-----'):
            header_end = i + 1
            break

# Charger les données après les lignes d'en-tête
data = np.loadtxt(filename, skiprows=header_end)
alpha_xfoil = data[:, 0]
cl_xfoil = data[:, 1]
cm_xfoil = data[:, 4]


#%% COMPUTATION OF CL AND CM WITH SPVP
CL = np.zeros(21)
CM = np.zeros(21)
for i in range(0,21):
    print(i/len(CL)*100, "%")
    XB,YB,XC,YC,S,delta,Cp,phi,lam,gamma,CL[i],CM[i] = SPVP(NameAirfoil,numPan,Vinf,alpha_xfoil[i],power=3)



#%% PLOT
# Tracer CL vs Alpha
fig2 = plt.figure(figsize=(8,4))
fig2.suptitle("NACA"+NameAirfoil)
plt.subplot(1, 2, 1)
plt.plot(alpha_xfoil, cl_xfoil, marker='o', color='r',label ="XFOIL")
plt.plot(alpha_xfoil, CL, marker='o', color='b',label ="SPVP")
plt.xlabel('Alpha [deg]')
plt.ylabel('CL')
plt.title('CL vs Alpha ('+str(numPan)+' panels)')
plt.legend()
plt.grid(True)

# Tracer CM vs Alpha
plt.subplot(1, 2, 2)
plt.plot(alpha_xfoil, cm_xfoil, marker='s', color='r',label ="XFOIL")
plt.plot(alpha_xfoil, CM, marker='s', color='b',label ="SPVP")
plt.xlabel('Alpha [deg]')
plt.ylabel('CM')
plt.title('CM vs Alpha ('+str(numPan)+' panels)')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
