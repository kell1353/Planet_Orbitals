import numpy as np
from itertools import product, combinations
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import itertools


# Setting the size of the window for the plot
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111, projection='3d')
fig.subplots_adjust(top=1,bottom=0,left=0,right=1)

# Set the aspect ratio to 1 so our spheres look spherical
ax.set_aspect('equal')

# Turn off the axis planes
ax.set_axis_off()

# Turn off tick labels
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

##fig.set_facecolor('black')
##ax.set_facecolor('black') 
ax.grid(False) 
ax.w_xaxis.pane.fill = False
ax.w_yaxis.pane.fill = False
ax.w_zaxis.pane.fill = False

limit = 4500000000
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-limit, limit)


# The Cartesian coordinates of the unit sphere          #####possibly redesign sphere
phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2*np.pi, 100)
phi, theta = np.meshgrid(phi, theta)

def draw_sphere(r, x, y, z, c):
    x = (r* np.sin(phi) * np.cos(theta)) + x
    y = (r*np.sin(phi) * np.sin(theta)) + y
    z = (r*np.cos(phi)) + z
    s = 7
    planet = ax.plot_surface(x, y, z,  rstride=s, cstride=s, color=c)                #(having higher strides allow the program to run faster)

    
#https://stackoverflow.com/questions/42264232/draw-ellipse-in-matplotlib-given-the-focii (Info for ellipses)
def elliptical_orbit_line(semimajor_axis, semiminor_axis, eccentricity):
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


def plot_multi_planet(r, color, semimajor_axis, semiminor_axis, eccentricity):
# A way to plot a planet on a given planets orbital line
    a1 = 0                                                                                                          # Foci one x-coordinate
    b1 = 0                                                                                                          # Foci one y-coordinate
    b2 = 0                                                                                                          # Foci two y-coordinate
    x0 = semimajor_axis * eccentricity                                                # Center x-value
    y0 = (b1 + b2) / 2                                                                                     # Center y-value
    a2 = (2 * x0) - b2                                                                                     # Foci two x-coordinate
    phi_ell = np.arctan2((b2 - b1), (a2 - a1))                                      # Angle between major axis and x-axis
    x = x0 + semimajor_axis * np.cos(2*np.pi*i/N) * np.cos(phi_ell) - semiminor_axis * np.sin(2*np.pi*i/N) * np.sin(phi_ell)
    y = y0 + semimajor_axis * np.cos(2*np.pi*i/N) * np.sin(phi_ell) + semiminor_axis * np.sin(2*np.pi*i/N) * np.cos(phi_ell)
    planet = draw_sphere(r, x, y, 0, color)


"""Calculate Orbits"""
mercury_orbit = elliptical_orbit_line(57910000, 56670300, .2056)
venus_orbit = elliptical_orbit_line(108210000, 107997400, .0067)
earth_orbit = elliptical_orbit_line(149600000, 149978300, .0167)
mars_orbit = elliptical_orbit_line(227920000, 226990500, .0935)
jupiter_orbit = elliptical_orbit_line(778570000, 778064300, .0489)
saturn_orbit = elliptical_orbit_line(1433530000, 1431240077, .0565)
uranus_orbit = elliptical_orbit_line(2872460000, 2866961900, .0457)
neptune_orbit = elliptical_orbit_line(4495060000, 4499727700, .0113)


# A way to plot a planet on a given circular line
##N = 8
##for i in range(N):
##    theta1 = 2*np.pi*i/N
##    x = np.cos(theta1)
##    y = np.sin(theta1)
##    draw_sphere(69911000, (x*778600000), (y*778600000), 0, 'orange')                              # on a unit circle


#https://matplotlib.org/gallery/color/named_colors.html (color lists link)    
"""Draw planet as spheres at specific points"""
"""Note these planet radii are scaled up by a factor of 10^3. The sun is just multiplied by 50."""
sun = draw_sphere(695956*50, 0, 0, 0, 'yellow')                     #####figure out how to scale this when you zoom in to graph.
    
repeat = 0
while repeat < 2:
    repeat = repeat + 1
    # A way to plot a planet on a given planets orbital line
    N=4
    for i in range(N):
        mercury = plot_multi_planet(2440000, 'lightgrey', 57910000, 56670300, .2056)
        venus = plot_multi_planet(6052000, 'navajowhite', 108210000, 107997400, .0067)
        earth = plot_multi_planet(6371000, 'dodgerblue', 149600000, 149978300, .0167)
        mars = plot_multi_planet(3390000, 'indianred', 227920000, 226990500, .0935)
        jupiter = plot_multi_planet(69911000, 'orange', 778570000, 778064300, .0489)
        saturn = plot_multi_planet(58232000, 'navajowhite', 1433530000, 1431240077, .0565)
        uranus = plot_multi_planet(25362000, 'lightblue', 2872460000, 2866961900, .0457)
        neptune = plot_multi_planet(24622000, 'blue', 4495060000, 4499727700, .0113)
##        fig.canvas.draw()
##        fig.canvas.draw_idle()
        plt.draw()
        plt.pause(.0000000001)
       # fig.remove(planet)
##        ax.clear(jupiter)

#https://matplotlib.org/gallery/color/named_colors.html (color lists link)    
"""Draw planet as spheres at specific points"""
"""Note these planet radii are scaled up by a factor of 10^3. The sun is just multiplied by 50."""
##sun = draw_sphere(695956*50, 0, 0, 0, 'yellow')                     #####figure out how to scale this when you zoom in to graph.
##mercury = draw_sphere(2440000, 57910000, 0, 0, 'lightgrey')
##venus = draw_sphere(6052000, 108210000, 0, 0, 'oldlace')
##earth = draw_sphere(6371000, 149600000, 0, 0, 'dodgerblue')
##mars = draw_sphere(3390000, 227920000, 0, 0, 'indianred')
##jupiter = draw_sphere(69911000, 778570000, 0, 0, 'orange')
##saturn = draw_sphere(58232000, 1433530000, 0, 0, 'navajowhite')
##uranus = draw_sphere(25362000, 2872460000, 0, 0, 'lightblue')
##neptune = draw_sphere(24622000, 4495060000, 0, 0, 'blue')


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

print("done")
#plt.show()
