import numpy as np
from itertools import product, combinations
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D


#Setting the size of the window for the plot
fig = plt.figure(figsize=(8,5.5))
ax = fig.add_subplot(111, projection='3d')
# Set the aspect ratio to 1 so our spheres look spherical
ax.set_aspect('equal')

limit = 4550000000
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-limit, limit)


# The Cartesian coordinates of the unit sphere
phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2*np.pi, 100)
phi, theta = np.meshgrid(phi, theta)

def draw_sphere(r, x, y, z, c):
    x = (r* np.sin(phi) * np.cos(theta)) + x
    y = (r*np.sin(phi) * np.sin(theta)) + y
    z = (r*np.cos(phi)) + z
    #ax.plot_surface(x, y, z,  rstride=1, cstride=1, color=c)
    ax.plot_surface(x, y, z,  rstride=4, cstride=4, color=c)                #(having higher strides allow the program to run faster)

    
#https://stackoverflow.com/questions/42264232/draw-ellipse-in-matplotlib-given-the-focii (Info for ellipses)

def elliptical_orbit(semimajor_axis, semiminor_axis, eccentricity):
    # Compute ellipse parameters
    a1 = 0                                                                                         # Foci one x-coordinate
    b1 = 0                                                                                         # Foci one y-coordinate
    b2 = 0                                                                                         # Foci two y-coordinate
    x0 = semimajor_axis * eccentricity                                                             # Center x-value
    y0 = (b1 + b2) / 2                                                                             # Center y-value
    f = np.sqrt((a1 - x0)**2 + (b1 - y0)**2)                                                       # Distance from center to focus
    a2 = (2 * x0) - b2                                                                             # Foci two x-coordinate
    phi_ell = np.arctan2((b2 - b1), (a2 - a1))                                                     # Angle between major axis and x-axis

    # Parametric plot in t
    resolution = 1000
    t = np.linspace(0, 2*np.pi, resolution)
    x = x0 + semimajor_axis * np.cos(t) * np.cos(phi_ell) - semiminor_axis * np.sin(t) * np.sin(phi_ell)
    y = y0 + semimajor_axis * np.cos(t) * np.sin(phi_ell) + semiminor_axis * np.sin(t) * np.cos(phi_ell)

    # Plot ellipse
    plt.plot(x, y, linewidth=.75)



"""Calculate Orbits"""
mercury_orbit = elliptical_orbit(57910000, 56670300, .2056)
venus_orbit = elliptical_orbit(108210000, 107997400, .0067)
earth_orbit = elliptical_orbit(149600000, 149978300, .0167)
mars_orbit = elliptical_orbit(227920000, 226990500, .0935)
jupiter_orbit = elliptical_orbit(778570000, 778064300, .0489)
saturn_orbit = elliptical_orbit(1433530000, 488114900, .0565)
uranus_orbit = elliptical_orbit(2872460000, 2866961900, .0457)
neptune_orbit = elliptical_orbit(4495060000, 4499727700, .0113)


#A way to plot a planet on a given line
N = 4
tuples = []
for i in range(N):
    theta1 = 2*np.pi*i/N
    x = np.cos(theta1)
    y = np.sin(theta1)
    draw_sphere(69911000, (x*778600000), (y*778600000), 0, 'orange')       # on a unit circle
    
##draw_sphere(69911000, (np.cos(2*np.pi*0/4)+778600000), (np.sin(2*np.pi*0/4)+0), 0, 'orange')
##draw_sphere(69911000, (np.cos(2*np.pi*1/4)+0), (np.sin(2*np.pi*1/4)+778600000), 0, 'orange')
##draw_sphere(69911000, (np.cos(2*np.pi*2/4)-778600000), (np.sin(2*np.pi*2/4)+0), 0, 'orange')
##draw_sphere(69911000, (np.cos(2*np.pi*3/4)+0), (np.sin(2*np.pi*3/4)-778600000), 0, 'orange')



#https://matplotlib.org/gallery/color/named_colors.html (color lists link)    
"""Draw planet as spheres at specific points"""
sun = draw_sphere(695956*50, 0, 0, 0, 'yellow')                     #####figure out how to scale this when you zoom in to graph.
mercury = draw_sphere(2440000, 57910000, 0, 0, 'lightgrey')
venus = draw_sphere(6052000, 108210000, 0, 0, 'oldlace')
earth = draw_sphere(6371000, 149600000, 0, 0, 'dodgerblue')
mars = draw_sphere(3390000, 227920000, 0, 0, 'indianred')
jupiter = draw_sphere(69911000, 778570000, 0, 0, 'orange')
saturn = draw_sphere(58232000, 1433530000, 0, 0, 'navajowhite')
uranus = draw_sphere(25362000, 2872460000, 0, 0, 'lightblue')
neptune = draw_sphere(24622000, 4495060000, 0, 0, 'blue')


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
