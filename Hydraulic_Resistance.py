import numpy as np
import math as math

def Hydraulic_Resistance(mu,L,type='circular',geometry=[1]):
    # if type = 'circular' ==> geometry = [radius]
    # if type = 'rectangular' ==> geometry = [length,width]
    # if type = 'any' ==> geometry = [hydraulic diameter]
    if type == 'circular' and len(geometry) == 1:
        Dh = 2*geometry[0]
    elif type == 'rectangle' and len(geometry) == 2:
        Dh = 4*geometry[0]*geometry[1]/(2*geometry[0]+2*geometry[1])
    elif type == 'any' and len(geometry) == 1:
        Dh = geometry[0]
    else:
        print('Error: geometry is not defined correctly')
    R = 128*mu*L/np.pi/Dh**4
    print('Resistance = ',R)
    return R
    
