import numpy as np
import math as math
import matplotlib.pyplot as plt


from GEOMETRY import GEOMETRY

def Hydraulic_GEOMETRY(XC,YC,delta,a,centre_point):
    # Define the list of control points that are in the pores
    control_points_in_pores = []
    a_slope = math.tan(delta[centre_point])
    b_low = YC[centre_point]- math.tan(delta[centre_point])*(XC[centre_point]+a/2/math.sin(delta[centre_point]))
    b_high = YC[centre_point] - math.tan(delta[centre_point])*(XC[centre_point]-a/2/math.sin(delta[centre_point]))
    if b_low > b_high:
        b_low, b_high = b_high, b_low
    for i in range(len(XC)):
        if YC[i] >= a_slope*XC[i]+b_low and YC[i] <= a_slope*XC[i]+b_high:
            control_points_in_pores.append(i)
    return filter(control_points_in_pores,centre_point)

def filter(liste, A):
    """
    Retourne tous les éléments connectés à A via des pas de +1 ou -1 dans la liste.
    """
    liste_set = set(liste)
    visited = set()
    à_visiter = [A]

    while à_visiter:
        courant = à_visiter.pop()
        if courant in liste_set and courant not in visited:
            visited.add(courant)
            if courant + 1 in liste_set:
                à_visiter.append(courant + 1)
            if courant - 1 in liste_set:
                à_visiter.append(courant - 1)
    
    return sorted(list(visited), key=liste.index)



def Refine_GEOMETRY(XB,NameAirfoil, entry_point, out_point,AoAR,n_refinement=5):
    m = int(NameAirfoil[0])*0.01
    p = int(NameAirfoil[1])*0.1
    t = int(NameAirfoil[2:4])*0.01
    c = 1
    if entry_point <= 0 or entry_point >= len(XB)-1:
        raise ValueError("i doit être un indice interne (1 <= i <= len(X)-2)")

    x_prev_in = XB[entry_point]
    x_curr_in = XB[entry_point+1]

    x_prev_out = XB[out_point]
    x_curr_out = XB[out_point+1]

    # Points de transition
    trans_in = np.linspace(x_prev_in, x_curr_in, n_refinement + 2)[1:-1]
    trans_out = np.linspace(x_prev_out, x_curr_out, n_refinement + 2)[1:-1]
    if entry_point < out_point:
        XB = np.concatenate([XB[:entry_point+1],trans_in, XB[entry_point+1:out_point+1],trans_out,XB[out_point+1:]])
        entry_point += (n_refinement-1)/2
        out_point +=n_refinement+(n_refinement-1)/2
    elif out_point < entry_point:
        XB = np.concatenate([XB[:out_point+1],trans_out,XB[out_point+1:entry_point+1],trans_in,XB[entry_point+1:]])
        out_point += (n_refinement-1)/2
        entry_point +=n_refinement+(n_refinement-1)/2
    
    for i in range(len(XB)):
        if XB[i] < XB[i+1]:
            lim = i
            break
    

    YB_DOWN = -5 * t * c * (0.2969 * np.sqrt(XB[:lim+1] / c) - 0.1260 * (XB[:lim+1] / c) - 0.3516 * (XB[:lim+1] / c)**2 
                      + 0.2843 * (XB[:lim+1] / c)**3 - 0.1015 * (XB[:lim+1] / c)**4)
    YB_UP = 5 * t * c * (0.2969 * np.sqrt(XB[lim+1:] / c) - 0.1260 * (XB[lim+1:] / c) - 0.3516 * (XB[lim+1:] / c)**2 
                      + 0.2843 * (XB[lim+1:] / c)**3 - 0.1015 * (XB[lim+1:] / c)**4)
    YB = np.concatenate([YB_DOWN,YB_UP])
    numPan = len(XB)-1
    XC,YC,S,phi,delta,beta = GEOMETRY(numPan,XB,YB,AoAR)
    print(len(phi))
    return XB,YB,XC,YC,S,phi,delta,beta,int(entry_point),int(out_point)


     
def plot_line(a, b, x_min=0, x_max=1):
        """
        Trace la droite y = a*x + b entre x_min et x_max

        Paramètres :
        - a : pente de la droite
        - b : ordonnée à l'origine
        - x_min, x_max : bornes de l'axe x
        """
        x = np.linspace(x_min, x_max, 100)  # 100 points entre x_min et x_max
        y = a * x + b  # équation de la droite

        plt.plot(x, y, label=f'y = {a}x + {b}')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Représentation d\'une droite')
        plt.legend()
        plt.grid(True)
