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
        z_list.append(r*np.tan(np.radians(deg_inclination)))

    # Plot ellipse
    plt.plot(x, y, z_list, linewidth=.25)



##def plot_planet(semimajor_axis, semiminor_axis, eccentricity, deg_inclination):
##        global pluto
##        # planet parameters
##        avgDist_pl = 3670000000
##        r_pl = 1187000
##        # orbit parameters
##        deg_inclination, eccentricity, init_theta = 17.16, .2488, 0
##        semimajor_axis, semiminor_axis = 5906380000, 5720653186
##        # ellipse parameters
##        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
##        x_center = semimajor_axis * eccentricity                                                # Center x-value
##        y_center = (b1 + b2) / 2                                                                                     # Center y-value
##        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate
##        
##        # Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
##        vel_n = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_pl*1000)) - (1/(semiminor_axis*1000)))))/1000
##        # Orbital Period calulation in years (using Kepler's Laws)
##        pl_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
##        # Scaling the orbital period in refernce to earths year
##        scaled_per = earth_year/pl_year 
##
##        # Getting individual orbital locations
##        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
##        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))
##
##        r = np.sqrt(((x - 0)**2))
##        if x > 0:
##                r = -(np.sqrt(((x - 0)**2)))
##        z = r*np.tan(np.radians(deg_inclination))
##
##        x0 = (r_pl* np.sin(phi) * np.cos(theta)) + x
##        y0 = (r_pl*np.sin(phi) * np.sin(theta)) + y
##        z0 = (r_pl*np.cos(phi)) + z
##        pluto = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='lightgrey')
##        #pl_text = ax.text(x, y, z + r_pl*100, "Pluto", color='white', ha = "center")
    
        
"""Draw Orbits"""
mercury_orbit = elliptical_orbit_line(57910000, 56670300, .2056, 7)
venus_orbit = elliptical_orbit_line(108210000, 107997400, .0067, 3.390)
earth_orbit = elliptical_orbit_line(149600000, 149978300, .0167, 0)
mars_orbit = elliptical_orbit_line(227920000, 226990500, .0935, 1.850)
jupiter_orbit = elliptical_orbit_line(778570000, 778064300, .0489, 1.304)
saturn_orbit = elliptical_orbit_line(1433530000, 1431240077, .0565, 2.485)
uranus_orbit = elliptical_orbit_line(2872460000, 2866961900, .0457, .772)
neptune_orbit = elliptical_orbit_line(4495060000, 4499727700, .0113, 1.769)
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

        
        ######################################################### Mercury #########################################################
        
        # planet parameters
        avgDist_me = 57900000
        r_me = 2440000
        # orbit parameters
        deg_inclination, eccentricity, init_theta = 7.000, .2056, 0
        semimajor_axis, semiminor_axis = 57910000, 56670300
        # ellipse parameters
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate

        #Orbital Velocity calculation in kilometers per second (using Kepler's Laws)
        vel_me = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_me*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        me_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/me_year
    
        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))
        
        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (r_me* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_me*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_me*np.cos(phi)) + z
        mercury = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='lightgrey')
        me_text = ax.text(x, y, z + r_me*10, "Mercury", color='white', ha = "center")




        ######################################################### Venus #########################################################
        
        #Initial planet information
        avgDist_v = 108200000
        r_v = 6052000
        # orbit parameters
        deg_inclination, eccentricity, init_theta = 3.390, .0067, 0
        semimajor_axis, semiminor_axis = 108210000, 107997400
        # ellipse parameters
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate
        
       #Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_v = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        v_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/v_year 
        
        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))
        
        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (r_v* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_v*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_v*np.cos(phi)) + z
        venus = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='navajowhite')
        v_text = ax.text(x, y, z + r_v*10, "Venus", color='white', ha = "center")




        ######################################################### Earth #########################################################
        
        # planet parameters
        avgDist_e = 149600000
        r_e = 6371000
        # orbit parameters
        deg_inclination, eccentricity, init_theta = 0.000, .0167, 0
        semimajor_axis, semiminor_axis = 149600000, 149978300
        # ellipse parameters
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate

        #Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_e = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        e_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/e_year
        
        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))
        
        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (r_e* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_e*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_e*np.cos(phi)) + z
        earth = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='dodgerblue')
        e_text = ax.text(x, y, z + r_e*10, "Earth", color='white', ha = "center")

        


        ######################################################### Mars #########################################################

        # planet parameters
        avgDist_ma = 227900000
        r_ma = 3390000
        # orbit parameters
        deg_inclination, eccentricity, init_theta = 1.850, .0935, 0
        semimajor_axis, semiminor_axis = 227920000, 226990500
        # ellipse parameters
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate

        # Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_ma = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        # Orbital Period calulation in years (using Kepler's Laws)
        ma_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        # Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/ma_year
             
        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))
        
        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (r_ma* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_ma*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_ma*np.cos(phi)) + z
        mars = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='indianred')
        ma_text = ax.text(x, y, z + r_ma*10, "Mars", color='white', ha = "center")
        


    
        ######################################################### Jupiter #########################################################
        
        # planet parameters
        avgDist_j = 778600000
        r_j = 69911000
        # orbit parameters
        deg_inclination, eccentricity, init_theta = 1.304, .0489, 0
        semimajor_axis, semiminor_axis  = 778570000, 778064300
        # ellipse parameters
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate

        # Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_j = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        # Orbital Period calulation in years (using Kepler's Laws)
        j_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        # Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/j_year
        
        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))
        
        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (r_j* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_j*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_j*np.cos(phi)) + z
        jupiter = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='orange')
        j_text = ax.text(x, y, z + r_j*2, "Jupiter", color='white', ha = "center")
        #jupiter1 = ax.contour(x0, y0, z0,  20, cmap=cm.Oranges) 




        ######################################################### Saturn #########################################################

        # planet parameters
        avgDist_s = 1433500000
        r_s = 58232000
        # orbit parameters
        deg_inclination, eccentricity, init_theta = 2.485, .0565, 0
        semimajor_axis, semiminor_axis = 1433530000, 1431240077
        # ellipse parameters
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate

       # Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_s = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        # Orbital Period calulation in years (using Kepler's Laws)
        s_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        # Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/s_year
        
        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))
        
        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (r_s* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_s*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_s*np.cos(phi)) + z
        saturn = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='navajowhite')
        s_text = ax.text(x, y, z + r_s*3, "Saturn", color='white', ha = "center")




        ######################################################### Uranus #########################################################

        # planet parameters
        avgDist_u = 2872500000
        r_u = 25362000
        # orbit parameters
        deg_inclination, eccentricity, init_theta = .772, .0457, 0
        semimajor_axis, semiminor_axis = 2872460000, 2866961900
        eccentricity = .0457
        # ellipse parameters
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate

       # Orbital Velocity Calculation in megameters per second (using Kepler's Laws)
        vel_u = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        # Orbital Period calulation in years (using Kepler's Laws)
        u_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        # Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/u_year
        
        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))
        
        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (r_u* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_u*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_u*np.cos(phi)) + z
        uranus = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='lightblue')
        u_text = ax.text(x, y, z + r_u*10, "Uranus", color='white', ha = "center")




        ######################################################### Neptune #########################################################
        
        # planet parameters
        avgDist_n = 4495100000
        r_n = 24622000
        # orbit parameters
        deg_inclination, eccentricity, init_theta = 1.769, .0113, 0
        semimajor_axis, semiminor_axis = 4495060000, 4499727700
        # ellipse parameters
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate

        # Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_n = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        # Orbital Period calulation in years (using Kepler's Laws)
        n_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        # Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/n_year

        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))
        
        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (r_n* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_n*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_n*np.cos(phi)) + z
        neptune = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='blue')
        n_text = ax.text(x, y, z + r_n*10, "Neptune", color='white', ha = "center")




        ######################################################### Pluto #########################################################

        # planet parameters
        avgDist_pl = 3670000000
        r_pl = 1187000
        # orbit parameters
        deg_inclination, eccentricity, init_theta = 17.16, .2488, 0
        semimajor_axis, semiminor_axis = 5906380000, 5720653186
        # ellipse parameters
        a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
        x_center = semimajor_axis * eccentricity                                                # Center x-value
        y_center = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x_center) - b2                                                                                     # Foci two x-coordinate
        
        # Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_n = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_pl*1000)) - (1/(semiminor_axis*1000)))))/1000
        # Orbital Period calulation in years (using Kepler's Laws)
        pl_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
        # Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/pl_year 

        # Getting individual orbital locations
        x = (x_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(np.radians(init_theta)))*np.cos(np.radians(deg_inclination))
        y = (y_center + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(np.radians(init_theta)) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(np.radians(init_theta)))

        r = np.sqrt(((x - 0)**2))
        if x > 0:
                r = -(np.sqrt(((x - 0)**2)))
        z = r*np.tan(np.radians(deg_inclination))

        x0 = (r_pl* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_pl*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_pl*np.cos(phi)) + z
        pluto = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color='lightgrey')
        pl_text = ax.text(x, y, z + r_pl*100, "Pluto", color='white', ha = "center")
        #pl_line, = plt.plot([x, x], [y, y], [z + (r_pl * 2), z + r_pl*500], linewidth=.25, color = "white")  
        
        #pluto_planet =  plot_planet(5906380000, 5720653186, .2488, 17.16)


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
            u_text.remove(), n_text.remove(), pl_text.remove()  #ax.lines.remove(pl_line)

    calc_range = (min_range+N)
    
print("The program has finished after calculating/plotting the planets orbits for " + str(repeat_amount) + " Earth years")
