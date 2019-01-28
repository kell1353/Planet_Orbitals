import numpy as np
from itertools import product, combinations
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

#fig = plt.figure()
#Setting the size of the window for the plot
fig = plt.figure(figsize=(8,5.5))
ax = fig.add_subplot(111, projection='3d')
# Set the aspect ratio to 1 so our spheres look spherical
ax.set_aspect('equal')

phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2*np.pi, 100)
phi, theta = np.meshgrid(phi, theta)

# The Cartesian coordinates of the unit sphere
def draw_sphere(r, x, y, z, c):
    x = (r* np.sin(phi) * np.cos(theta)) + x
    y = (r*np.sin(phi) * np.sin(theta)) + y
    z = (r*np.cos(phi)) + z
    #ax.plot_surface(x, y, z,  rstride=1, cstride=1, color=c)
    ax.plot_surface(x, y, z,  rstride=4, cstride=4, color=c)   #(having higher strides allow the program to run faster)

limit = 4550000000
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-limit, limit)

#https://matplotlib.org/gallery/color/named_colors.html             (color lists link)

#draw planet as spheres at specific points
sun = draw_sphere(695956*50, 0, 0, 0, 'yellow')                     #figure out how to scale this when you zoom in to graph.
mercury = draw_sphere(2440000, 57900000, 0, 0, 'lightgrey')
venus = draw_sphere(6052000, 108200000, 0, 0, 'oldlace')
earth = draw_sphere(6371000, 149600000, 0, 0, 'dodgerblue')
mars = draw_sphere(3390000, 227900000, 0, 0, 'indianred')
jupiter = draw_sphere(69911000, 778600000, 0, 0, 'orange')
saturn = draw_sphere(58232000, 1433500000, 0, 0, 'navajowhite')
uranus = draw_sphere(25362000, 2872500000, 0, 0, 'lightblue')
neptune = draw_sphere(24622000, 4495100000, 0, 0, 'blue')



##def init():
##    circ.center = (0, 0)
##    ax.add_patch(circ)
##    return circ,
##
##def animate(i):
##    x, y = circ.center
##    x = 5 + 3 * np.sin(np.radians(i))
##    y = 5 + 3 * np.cos(np.radians(i))
##    circ.center = (x, y)
##    return circ,

##anim = animation.FuncAnimation(fig,animate,init_func=init,frames=360,interval=20,blit=True)



# Turn off the axis planes
#ax.set_axis_off()

# Turn off tick labels
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

fig.set_facecolor('black')
ax.set_facecolor('black') 
ax.grid(False) 
ax.w_xaxis.pane.fill = False
ax.w_yaxis.pane.fill = False
ax.w_zaxis.pane.fill = False


plt.show()
