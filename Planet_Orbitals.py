import numpy as np
from itertools import product, combinations
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from matplotlib import cm, colors

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Set the aspect ratio to 1 so our spheres looks spherical
ax.set_aspect('equal')


phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2*np.pi, 100)
phi, theta = np.meshgrid(phi, theta)

# The Cartesian coordinates of the unit sphere
def draw_sphere(r, x, y, z, c):
    x = (r* np.sin(phi) * np.cos(theta)) + x
    y = (r*np.sin(phi) * np.sin(theta)) + y
    z = (r*np.cos(phi)) + z
    ax.plot_surface(x, y, z,  rstride=1, cstride=1, color=c)
    #ax.plot_surface(x, y, z,  rstride=10, cstride=10, color=c)   (this may allow program to run faster)
    
limit = 4750000
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-limit, limit)

#https://matplotlib.org/gallery/color/named_colors.html             (color lists link)

#draw planet as spheres at specific points
sun = draw_sphere(695956, 0, 0, 0, 'cornsilk')
##mercury = draw_sphere(2440, 57900, 0, 0, 'lightgrey')
##venus = draw_sphere(6052, 108200, 0, 0, 'oldlace')
##earth = draw_sphere(6371, 149600, 0, 0, 'dodgerblue')
##mars = draw_sphere(3390, 227900, 0, 0, 'indianred')
##jupiter = draw_sphere(69911, 778600, 0, 0, 'orange')
##saturn = draw_sphere(58232, 1433500, 0, 0, 'navajowhite')
##uranus = draw_sphere(25362, 2872500, 0, 0, 'lightblue')
##neptune = draw_sphere(24622, 4495100, 0, 0, 'blue')


#ax.scatter([0], [0], [0], color="y", s=50)
#ax.scatter([57900], [0], [0], color="r", s=2)
#ax.scatter([108200], [0], [0], color="r", s=0.00171431445)
#ax.scatter([149600], [0], [0], color="b", s=0.3002187524)
#ax.scatter([227900], [0], [0], color="r", s=0.03228483066)
#ax.scatter([778600], [0], [0], color="r", s=95.44643082)
#ax.scatter([1433500], [0], [0], color="y", s=28.56352619)
#ax.scatter([2872500], [0], [0], color="b", s=4.364989565)
#ax.scatter([4495100], [0], [0], color="b", s=5.129365618)

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

# Turn off the axis planes
#ax.set_axis_off()

fig.set_facecolor('black')
ax.set_facecolor('black') 
ax.grid(False) 
ax.w_xaxis.pane.fill = False
ax.w_yaxis.pane.fill = False
ax.w_zaxis.pane.fill = False
plt.show()
