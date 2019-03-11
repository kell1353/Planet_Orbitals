import pylab
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
plt.title('Solar System', fontsize = 14, color = 'white')
fig.subplots_adjust(top=1,bottom=0,left=0,right=1)


# Turn off the axis planes
ax.set_axis_off()

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

limit = 4500000000
ax.set_xlim(-limit, limit)
ax.set_ylim(-limit, limit)
ax.set_zlim(-limit, limit)

# Set the aspect ratio to 1 so our spheres look spherical
ax.set_aspect('equal')

# The Cartesian coordinates of the unit sphere              #####possibly redesign sphere
points_range = 25
phi = np.linspace(0, np.pi, points_range)
theta = np.linspace(0, 2*np.pi, points_range)
phi, theta = np.meshgrid(phi, theta)



def draw_sphere(r, x, y, z, c):
    x = (r* np.sin(phi) * np.cos(theta)) + x
    y = (r*np.sin(phi) * np.sin(theta)) + y
    z = (r*np.cos(phi)) + z
    s = 1
    planet = ax.plot_surface(x, y, z,  rstride=s, cstride=s, color=c)                #(having higher strides allow the program to run faster)


    
#https://stackoverflow.com/questions/42264232/draw-ellipse-in-matplotlib-given-the-focii (Info for ellipses)
def elliptical_orbit_line(semimajor_axis, semiminor_axis, eccentricity, deg_inclination):
    # Compute ellipse parameters
    a1, b1, b2 = 0, 0, 0                                                                                               # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
    x0 = semimajor_axis * eccentricity                                                            # Center x-value
    y0 = (b1 + b2) / 2                                                                                                 # Center y-value
    f = np.sqrt((a1 - x0)**2 + (b1 - y0)**2)                                                       # Distance from center to focus
    a2 = (2 * x0) - b2                                                                                                  # Foci two x-coordinate
    phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                                   # Angle between major axis and x-axis

    z_list = []
    # Parametric plot in t
    resolution = 100
    t = np.linspace(0, 2*np.pi, resolution)
    x = ((semimajor_axis*eccentricity) + semimajor_axis * np.cos(t) * np.cos(np.radians(phi_ell)) - semiminor_axis * np.sin(t) * np.sin(np.radians(phi_ell)))*np.cos(np.radians(deg_inclination))
    y = y0 + semimajor_axis * np.cos(t) * np.sin(np.radians(phi_ell)) + semiminor_axis * np.sin(t) * np.cos(np.radians(phi_ell))

    for i in range(0, len(x)):  
        r = np.sqrt(((x[i]-0)**2))
        if x[i] > 0:
            r = -(np.sqrt(((x[i]-0)**2)))
        z = r*np.tan(np.radians(deg_inclination))
        z_list.append(z)

    # Plot ellipse
    plt.plot(x, y, z_list, linewidth=.25)



def plot_object(name, semimajor_axis, semiminor_axis, eccentricity, deg_inclination, planet_rad, color):
        global mercury; global me_text
        global venus; global v_text
        global earth; global e_text
        global mars; global ma_text
        global jupiter; global j_text
        global saturn; global s_text
        global uranus; global u_text
        global neptune; global n_text
        global pluto; global pluto_text

        # ellipse parameters
        init_theta = 0
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate
        
        # Orbital Period calulation in years (using Kepler's Laws)
        planet_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        # Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/planet_year 

        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))

        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (planet_rad* np.sin(phi) * np.cos(theta)) + x
        y0 = (planet_rad*np.sin(phi) * np.sin(theta)) + y
        z0 = (planet_rad*np.cos(phi)) + z

        if name == "Mercury":
            mercury = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color= color)
            me_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Venus":
            venus = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color= color)
            v_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Earth":
            earth= ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color= color)
            e_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Mars":
            mars = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color= color)
            ma_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Jupiter":
            jupiter = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color= color)
            j_text = ax.text(x, y, z + planet_rad*2, name, color='white', ha = "center")
        elif name == "Saturn":
            saturn = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color= color)
            s_text = ax.text(x, y, z + planet_rad*3, name, color='white', ha = "center")
        elif name == "Uranus":
            uranus = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color= color)
            u_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Neptune":
            neptune = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color= color)
            n_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Pluto":
            pluto = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color= color)
            pluto_text = ax.text(x, y, z + planet_rad*100, name, color='white', ha = "center")
    

        
"""Draw Celestial Body Orbits"""
mercury_orbit = elliptical_orbit_line(57910000, 56672817, .2056, 7.000)
venus_orbit = elliptical_orbit_line(108210000, 108207571, .0067, 3.390)
earth_orbit = elliptical_orbit_line(149600000, 149579138, .0167, 0)
mars_orbit = elliptical_orbit_line(227920000, 226921546, .0935, 1.850)
jupiter_orbit = elliptical_orbit_line(778570000, 777638581, .0489, 1.304)
saturn_orbit = elliptical_orbit_line(1433530000, 1431240078, .0565, 2.485)
uranus_orbit = elliptical_orbit_line(2872460000, 2869458880, .0457, .772)
neptune_orbit = elliptical_orbit_line(4495060000, 4494773004, .0113, 1.769)
pluto_orbit =  elliptical_orbit_line(5906380000, 5720653186, .2488, 17.16)


#https://matplotlib.org/gallery/color/named_colors.html (color lists link)    
"""Draw planet as spheres at specific points"""
"""Note the planet radii are scaled up by a factor of 10^3. The sun is just multiplied by 50."""
grav_constant =  6.67*(10**(-11))
mass_sun = 1.98855*(10**30)
earth_year = 1.0010
planet_smooth = 3

sun = draw_sphere(695956*50, 0, 0, 0, 'yellow')                     #####figure out how to scale this when you zoom in to graph.

""" Starting parameters """
#repeat_amount = 5
repeat_amount = float(input("How many years would you like to see?: "))
N = 30                            # The lower the value the faster and more jumpy the planets orbit. (25 is ideal)
min_range = 0
repeat = 0
calc_range = N


while repeat < repeat_amount:
    repeat = repeat + 1
    # A way to plot a planet on a given planets orbital line
    for i in range(min_range, calc_range):
        
        # Show the amount of years have pass ed during the orbits
        year_calc = round(i/N, 1)
        fig.suptitle(str(year_calc) + ' Earth years' , y = .8, fontsize = 9, color = 'white')


        mercury_planet =  plot_object("Mercury", 57910000, 56672817, .2056, 7.000, 2440000, 'lightgrey')
        venus_planet =  plot_object("Venus", 108210000, 108207571, .0067, 3.390, 6052000, 'navajowhite')
        earth_planet =  plot_object("Earth", 149600000, 149579138, .0167, 0.000, 6371000, 'dodgerblue')
        mars_planet =  plot_object("Mars", 227920000, 226921546, .0935, 1.850, 3390000, 'indianred')
        jupiter_planet =  plot_object("Jupiter", 778570000, 777638581, .0489, 1.304, 69911000, 'orange')
##        ##jupiter1 = ax.contour(x0, y0, z0,  20, cmap=cm.Oranges) 
        saturn_planet =  plot_object("Saturn", 1433530000, 1431240078, .0565, 2.485, 58232000, 'navajowhite')
        uranus_planet =  plot_object("Uranus", 2872460000, 2869458880, .0457, .772, 25362000, 'lightblue')
        neptune_planet =  plot_object("Neptune", 4495060000, 4494773004, .0113, 1.769, 24622000, 'blue')
        pluto_planet =  plot_object("Pluto", 5906380000, 5720653186, .2488, 17.16, 1187000,'lightgrey')
##        ##pl_line, = plt.plot([x, x], [y, y], [z + (r_pl * 2), z + r_pl*500], linewidth=.25, color = "white")  

        min_range = i + 1
        if year_calc == repeat_amount:
            plt.show()           
        else:
            # Removing the previously plotted celestial objects
            plt.pause(.00000000000000000000000000001)
            ax.collections.remove(mercury), ax.collections.remove(venus), ax.collections.remove(earth)
            ax.collections.remove(mars), ax.collections.remove(jupiter), ax.collections.remove(saturn)
            ax.collections.remove(uranus), ax.collections.remove(neptune), ax.collections.remove(pluto)
            #for npl in jupiter1.collections:
                    #npl.remove()
            # Removing the previously drawn labels and lines
            me_text.remove(), v_text.remove(), e_text.remove()
            ma_text.remove(), j_text.remove(), s_text.remove()
            u_text.remove(), n_text.remove(), pluto_text.remove()  #ax.lines.remove(pl_line)

    calc_range = (min_range+N)



"""Other Planet Calculations"""

def calc_avgVelocity(averageDistance, semiminor_axis):
    # Average Orbital Velocity calculation in kilometers per second (using Kepler's Law)
    avgVelocity = (np.sqrt((grav_constant*mass_sun)*((2/(averageDistance*1000)) - (1/(semiminor_axis*1000)))))/1000

"""Mercury """
avgDist_mercury = 57900000
semiminor_axis =  56672817
vel_mercury = calc_avgVelocity(avgDist_mercury, semiminor_axis)

""" Venus """
avgDist_venus = 108200000
semiminor_axis =  108207571
vel_venus = calc_avgVelocity(avgDist_venus, semiminor_axis)

""" Earth """
avgDist_earth = 149600000
semiminor_axis =  149579138
vel_earth = calc_avgVelocity(avgDist_earth, semiminor_axis)

""" Mars """
avgDist_mars = 227900000
semiminor_axis =  226921546
vel_mars = calc_avgVelocity(avgDist_mars, semiminor_axis)

""" Jupiter """
avgDist_jupiter = 778600000
semiminor_axis =  777638581
vel_jupiter = calc_avgVelocity(avgDist_jupiter, semiminor_axis)

""" Saturn """
avgDist_saturn = 1433500000
semiminor_axis =  1431240078
vel_saturn = calc_avgVelocity(avgDist_saturn, semiminor_axis)

""" Uranus """
avgDist_uranus = 287250000
semiminor_axis =  2869458880
vel_uranus = calc_avgVelocity(avgDist_uranus, semiminor_axis)

""" Neptune """
avgDist_neptune = 4495100000
semiminor_axis =  4494773004
vel_neptune = calc_avgVelocity(avgDist_neptune, semiminor_axis)

""" Pluto """
avgDist_pluto = 3670000000
semiminor_axis =  5720653186
vel_pluto = calc_avgVelocity(avgDist_pluto, semiminor_axis)

    
print("The program has finished after calculating/plotting the planets orbits for " + str(repeat_amount) + " Earth years")
