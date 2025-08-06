import matplotlib.pyplot as plt
import numpy as np

NumPan = [20,50,100,200,300,400,500,750,1000,2000]
CD_equal = [-0.00365971, 0.00095220,  0.00008811, -0.00019255, -0.00017924, -0.00015313, -0.00013156, -0.00009572, -0.00007481, -0.00003978]
CD_power1 = [0.00535930, 0.00170207, 0.00069274, 0.00030715, 0.00019597, 0.00014364, 0.00011330, 0.00007408, 0.00005501, 0.00002709]
CD_power3 = [0.01038477, 0.00383386, 0.00181370, 0.00087362, 0.00057388, 0.00042697, 0.00033985, 0.00022495, 0.00016808, 0.00008353]
CD_power0_66 = [-0.00510953, -0.00027364, -0.00077545, -0.00107796, -0.00098743, -0.00087043, -0.00076976, -0.00059133, -0.00047912, -0.00027376]

CD_power3 = np.array(CD_power3)
CD_power1 = np.array(CD_power1)
CD_equal = np.array(CD_equal)
CD_power0_66 = np.array(CD_power0_66)


NumPan_ellipse = [20,50,100,200,300,400,500,750]
CD_equal_ellipse = [-0.00000731, 0.00000269, -0.00000189, -0.00000183, 0.00000012, -0.00000192, 0.00000017, 0.00000001]
CD_power1_ellipse = [0,0,0,0,0,0,0,0]
CD_power3_ellipse = [0.00001790, -0.00000208, 0.00000002, 0.00000182, -0.00000026,  0.00000085,  -0.00000028, -0.00000043]
CD_power0_66_ellipse = [-0.00000731, 0.00000269, -0.00000189,-0.00000183]

fig = plt.figure()
plt.title('Err of CD vs Number of panels (NACA 0018)')
plt.plot(NumPan,abs(CD_equal),c='r',marker = 'd',label='Equal size panels')
plt.plot(NumPan,abs(CD_power1),c='g',marker = '.',label='Power = 1')
plt.plot(NumPan,abs(CD_power3),c='b',marker = 'v',label='Power = 3')
plt.plot(NumPan,abs(CD_power0_66),c='y',marker = 'h',label='Power = 2/3')

plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Panels number')
plt.ylabel('Error of Drag Coefficient $C_D$')
plt.legend()
plt.ylim(0,0.015)
plt.show()