# my plot try from data from LED data LEDOutput.txt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt # for plotting stuff
outp='LEDOutput.txt' # read file is assumed
inp='LEDInput.csv' # read file imput parameters
p=np.genfromtxt(inp,delimiter=',')
p=p[1,:]
x=p[0]
dx=p[1]
t=p[2]
dt=p[3]
a=p[4]
u = np.genfromtxt(outp)
n,m=u.shape
#~ mpl.rcParams['legend.fontsize']=10
fig=plt.figure()
ax=fig.gca(projection='3d')
#~ fig, ax=plt.subplots()
#~ ax.plot(np.linspace(0,m-1,m),np.zeros((m)),u[0,:] ,linewidth=4,label='Time=0')
ax.plot(dx*np.linspace(0,m-1,m),np.zeros((m)),u[0,:] ,linewidth=4)
#~ for i in range(1,n): ax.plot(np.linspace(0,m-1,m),i*np.ones((m)),u[i,:],label=str(i))
for i in range(1,n): ax.plot(dx*np.linspace(0,m-1,m),i*dt*np.ones((m)),u[i,:])
#~ ax.legend()
ax.set_title(r'$\frac{du}{dt}+a\cdot\frac{du}{dx}=f(x,t)$')
ax.set_xlabel(r'$bar\ length$')
ax.set_ylabel(r'$time$')
ax.set_zlabel(r'$u(x,t)$')
plt.show()
