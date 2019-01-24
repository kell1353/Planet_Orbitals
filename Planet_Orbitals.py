from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
from matplotlib import animation

fig = plt.figure()
ax = fig.gca(projection='3d')
##ax.set_aspect("equal")

limit = 4750000
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-1, 1)

#draw a point
#ax.scatter([0], [0], [0], color="y", s=50)
#ax.scatter([57900], [0], [0], color="r", s=2)
#ax.scatter([108200], [0], [0], color="r", s=0.00171431445)
#ax.scatter([149600 ], [0], [0], color="b", s=0.3002187524)
#ax.scatter([227900], [0], [0], color="r", s=0.03228483066)
ax.scatter([778600], [0], [0], color="r", s=95.44643082)
#ax.scatter([1433500], [0], [0], color="y", s=28.56352619)
#ax.scatter([2872500], [0], [0], color="b", s=4.364989565)
#ax.scatter([4495100], [0], [0], color="b", s=5.129365618)

#ax = plt.axes(xlim=(0, 10), ylim=(0, 10))
#circ = plt.Circle((5, -5), 0.75, fc='b')

#def init():
    circ.center = (0, 0)
    ax.add_patch(circ)
    return circ,

#def animate(i):
    x, y = circ.center
    x = 5 + 3 * np.sin(np.radians(i))
    y = 5 + 3 * np.cos(np.radians(i))
    circ.center = (x, y)
    return circ,

#anim = animation.FuncAnimation(fig,animate,init_func=init,frames=360,interval=20,blit=True)

fig.set_facecolor('black')
ax.set_facecolor('black') 
ax.grid(False) 
ax.w_xaxis.pane.fill = False
ax.w_yaxis.pane.fill = False
ax.w_zaxis.pane.fill = False
plt.show()
