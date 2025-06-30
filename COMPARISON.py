import numpy as np
import math as math
import matplotlib.pyplot as plt
from matplotlib import path

from SPVP_Airfoil import SPVP
from XFOIL import XFOIL_DATA


NameAirfoil = "2412"
numPan = 500
Vinf = 1
AoA = 4


XB,YB,XC,YC,S,delta,Cp,phi,lam,gamma,CL,CM = SPVP(NameAirfoil,numPan,Vinf,AoA,power=1)
x_xfoil,cp_xfoil = XFOIL_DATA(NameAirfoil,4)


plt.figure()
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
plt.title("Coefficient de Pression à α = 4°")
plt.grid(True)
plt.legend()


#######################################################

# Chemin vers ton fichier
filename = "XFOIL_DATA/2412_CL_CM.dat"

# Ouvre le fichier pour voir combien de lignes d'en-tête il y a
with open(filename, 'r') as f:
    for i, line in enumerate(f):
        print(f"Ligne {i+1}: {line.strip()}")
        if line.strip().startswith('-----'):
            header_end = i + 1
            break

# Charger les données après les lignes d'en-tête
data = np.loadtxt(filename, skiprows=header_end)

# Extraire colonnes :
# typiquement :
# Colonne 0 → Alpha
# Colonne 1 → CL
# Colonne 4 → CM
alpha_xfoil = data[:, 0]
cl_xfoil = data[:, 1]
cm_xfoil = data[:, 4]

CL = np.zeros(21)
CM = np.zeros(21)
for i in range(0,21):
    print(i/len(CL)*100, "%")
    XB,YB,XC,YC,S,delta,Cp,phi,lam,gamma,CL[i],CM[i] = SPVP(NameAirfoil,numPan,Vinf,alpha_xfoil[i],power=1)




# Tracer CL vs Alpha
plt.figure(figsize=(8,4))
plt.subplot(1, 2, 1)
plt.plot(alpha_xfoil, cl_xfoil, marker='o', color='r',label ="XFOIL")
plt.plot(alpha_xfoil, CL, marker='o', color='b',label ="SPVP")
plt.xlabel('Alpha [deg]')
plt.ylabel('CL')
plt.title('CL vs Alpha')
plt.legend()
plt.grid(True)

# Tracer CM vs Alpha
plt.subplot(1, 2, 2)
plt.plot(alpha_xfoil, cm_xfoil, marker='s', color='r',label ="XFOIL")
plt.plot(alpha_xfoil, CM, marker='s', color='b',label ="SPVP")
plt.xlabel('Alpha [deg]')
plt.ylabel('CM')
plt.title('CM vs Alpha')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
