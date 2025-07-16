Delta_P = 0.2

geometry = [0.157,0.025]
w=geometry[0]
h=geometry[1]
Dh = 4*w*h/(2*(w+h))
L = 0.89
mu = 1/160000
rhoinf = 1

V_mean=2.868*((Delta_P**4)*(Dh**5)/((L**4)*mu*(rhoinf**3)))**(1/7)
print(V_mean)