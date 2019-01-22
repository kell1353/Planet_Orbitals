from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations


fig = plt.figure()
ax = fig.gca(projection='3d')
##ax.set_aspect("equal")
limit = 4750000

ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-1, 1)

### draw cube
##r = [-1, 1]
##for s, e in combinations(np.array(list(product(r, r, r))), 2):
##    if np.sum(np.abs(s-e)) == r[1]-r[0]:
##        ax.plot3D(*zip(s, e), color="b")

#draw sphere
##u, v = np.mgrid[0:2*np.pi:30j, 0:np.pi:15j]
##x = np.cos(u)*np.sin(v)
##y = np.sin(u)*np.sin(v)
##z = np.cos(v)
##ax.plot_wireframe(x, y, z, color="r")

##r = 1
##pi = np.pi
##cos = np.cos
##sin = np.sin
##phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
##x = r*sin(phi)*cos(theta)
##y = r*sin(phi)*sin(theta)
##z = r*cos(phi)
##
##ax.plot_surface(
##    x, y, z,  rstride=1, cstride=1, color='c', alpha=0.6, linewidth=0)

#draw a point
ax.scatter([0], [0], [0], color="y", s=50)
ax.scatter([57900], [0], [0], color="r", s=2)
ax.scatter([108200], [0], [0], color="r", s=0.00171431445)
ax.scatter([149600 ], [0], [0], color="b", s=0.3002187524)
ax.scatter([227900], [0], [0], color="r", s=0.03228483066)
ax.scatter([778600], [0], [0], color="r", s=95.44643082)
ax.scatter([1433500], [0], [0], color="y", s=28.56352619)
ax.scatter([2872500], [0], [0], color="b", s=4.364989565)
ax.scatter([4495100], [0], [0], color="b", s=5.129365618)

# draw a vector
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


##class Arrow3D(FancyArrowPatch):
##
##    def __init__(self, xs, ys, zs, *args, **kwargs):
##        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
##        self._verts3d = xs, ys, zs
##
##    def draw(self, renderer):
##        xs3d, ys3d, zs3d = self._verts3d
##        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
##        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
##        FancyArrowPatch.draw(self, renderer)
##
##a = Arrow3D([0, 1], [0, 1], [0, 1], mutation_scale=20,
##            lw=1, arrowstyle="-|>", color="k")
##ax.add_artist(a)
fig.set_facecolor('black')
ax.set_facecolor('black') 
ax.grid(False) 
ax.w_xaxis.pane.fill = False
ax.w_yaxis.pane.fill = False
ax.w_zaxis.pane.fill = False
plt.show()
