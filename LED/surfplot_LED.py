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
fig=plt.figure()
ax=fig.gca(projection='3d')
x=dx*np.linspace(0,m-1,m)
y=dt*np.linspace(0,n-1,n)
x,y = np.meshgrid(x,y)
z=u
surf=ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='jet',
        linewidth=0, antialiased=False)

ax.set_title(r'$\frac{du}{dt}+a\cdot\frac{du}{dx}=f(x)$')
ax.set_xlabel(r'$bar\ length$')
ax.set_ylabel(r'$time$')
ax.set_zlabel(r'$u(x,t)$')
fig.colorbar(surf)
plt.show()
