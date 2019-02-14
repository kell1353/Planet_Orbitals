import numpy as np
from itertools import product, combinations
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import itertools


# Setting the size of the window for the plot
fig = plt.figure(figsize=(8,5))
fig.suptitle('Solar System', fontsize = 14)
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
def elliptical_orbit_line(semimajor_axis, semiminor_axis, eccentricity):
    # Compute ellipse parameters
    a1 = 0                                                                                                                      # Foci one x-coordinate
    b1 = 0                                                                                                                      # Foci one y-coordinate
    b2 = 0                                                                                                                      # Foci two y-coordinate
    x0 = semimajor_axis * eccentricity                                                            # Center x-value
    y0 = (b1 + b2) / 2                                                                                                 # Center y-value
    f = np.sqrt((a1 - x0)**2 + (b1 - y0)**2)                                                       # Distance from center to focus
    a2 = (2 * x0) - b2                                                                                                  # Foci two x-coordinate
    phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                                   # Angle between major axis and x-axis

    # Parametric plot in t
    resolution = 1000
    t = np.linspace(0, 2*np.pi, resolution)
    x = x0 + semimajor_axis * np.cos(t) * np.cos(phi_ell) - semiminor_axis * np.sin(t) * np.sin(phi_ell)
    y = y0 + semimajor_axis * np.cos(t) * np.sin(phi_ell) + semiminor_axis * np.sin(t) * np.cos(phi_ell)

    # Plot ellipse
    plt.plot(x, y, linewidth=.75)



"""Draw Orbits"""
mercury_orbit = elliptical_orbit_line(57910000, 56670300, .2056)
venus_orbit = elliptical_orbit_line(108210000, 107997400, .0067)
earth_orbit = elliptical_orbit_line(149600000, 149978300, .0167)
mars_orbit = elliptical_orbit_line(227920000, 226990500, .0935)
jupiter_orbit = elliptical_orbit_line(778570000, 778064300, .0489)
saturn_orbit = elliptical_orbit_line(1433530000, 1431240077, .0565)
uranus_orbit = elliptical_orbit_line(2872460000, 2866961900, .0457)
neptune_orbit = elliptical_orbit_line(4495060000, 4499727700, .0113)


#https://matplotlib.org/gallery/color/named_colors.html (color lists link)    
"""Draw planet as spheres at specific points"""
"""Note the planet radii are scaled up by a factor of 10^3. The sun is just multiplied by 50."""
grav_constant =  6.67*(10**(-11))
mass_sun = 1.98855*(10**30)
earth_year = 1.0010
planet_smooth = 3

sun = draw_sphere(695956*50, 0, 0, 0, 'yellow')                     #####figure out how to scale this when you zoom in to graph.

""" Starting parameters """
repeat_amount = 4
N = 50
min_range = 0
repeat = 0
calc_range = N


while repeat < repeat_amount:
    repeat = repeat + 1
    # A way to plot a planet on a given planets orbital line
    for i in range(min_range, calc_range):

        
        
        ######################################################### Mercury #########################################################
        
        #Initial planet information
        avgDist_me = 57900000
        r_me = 2440000
        semimajor_axis = 57910000
        semiminor_axis = 56670300
        color = 'lightgrey'
        eccentricity = .2056
        a1 = 0                                                                                                          # Foci one x-coordinate
        b1 = 0                                                                                                          # Foci one y-coordinate
        b2 = 0                                                                                                          # Foci two y-coordinate
        x0 = semimajor_axis * eccentricity                                                # Center x-value
        y0 = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x0) - b2                                                                                     # Foci two x-coordinate
        phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                      # Angle between major axis and x-axis


        #Orbital Velocity calculation in kilometers per second (using Kepler's Laws)
        vel_me = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_me*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        me_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)

                   #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/me_year

        
        #Getting individual orbital locations
        x = x0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per))* np.sin(phi_ell)
        y = y0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell)
        z = 0
        
        #Plotting the spheres in their location
        x0 = (r_me* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_me*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_me*np.cos(phi)) + z
        mercury = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color=color)




        ######################################################### Venus #########################################################
        
        #Initial planet information
        avgDist_v = 108200000
        r_v = 6052000
        semimajor_axis = 108210000
        semiminor_axis = 107997400
        color = 'navajowhite'
        eccentricity = .0067
        a1 = 0                                                                                                          # Foci one x-coordinate
        b1 = 0                                                                                                         # Foci one y-coordinate
        b2 = 0                                                                                                          # Foci two y-coordinate
        x0 = semimajor_axis * eccentricity                                                # Center x-value
        y0 = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x0) - b2                                                                                     # Foci two x-coordinate
        phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                      # Angle between major axis and x-axis
        
 
       #Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_v = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        v_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)

                   #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/v_year 
        
        #Getting individual orbital locations
        x = x0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell)
        y = y0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell)
        z = 0
        
        #Plotting the spheres in their location
        x0 = (r_v* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_v*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_v*np.cos(phi)) + z
        venus = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color=color)




        ######################################################### Earth #########################################################
        
        #Initial planet information
        avgDist_e = 149600000
        r_e = 6371000
        semimajor_axis = 149600000
        semiminor_axis = 149978300
        color = 'dodgerblue'
        eccentricity = .0167
        a1 = 0                                                                                                          # Foci one x-coordinate
        b1 = 0                                                                                                          # Foci one y-coordinate
        b2 = 0                                                                                                          # Foci two y-coordinate
        x0 = semimajor_axis * eccentricity                                                # Center x-value
        y0 = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x0) - b2                                                                                     # Foci two x-coordinate
        phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                      # Angle between major axis and x-axis


        #Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_e = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        e_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)

                   #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/e_year
        
        #Getting individual orbital locations      
        x = x0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell)
        y = y0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell)
        z = 0
        
        #Plotting the spheres in their location
        x0 = (r_e* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_e*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_e*np.cos(phi)) + z
        earth = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color=color)

        


        ######################################################### Mars #########################################################

        #Initial planet information
        avgDist_ma = 227900000
        r_ma = 3390000
        semimajor_axis = 227920000
        semiminor_axis = 226990500
        color = 'indianred'
        eccentricity = .0935
        a1 = 0                                                                                                          # Foci one x-coordinate
        b1 = 0                                                                                                          # Foci one y-coordinate
        b2 = 0                                                                                                          # Foci two y-coordinate
        x0 = semimajor_axis * eccentricity                                                # Center x-value
        y0 = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x0) - b2                                                                                     # Foci two x-coordinate
        phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                      # Angle between major axis and x-axis


        #Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_ma = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        ma_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)

                   #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/ma_year
        
        
        #Getting individual orbital locations                                        
        x = x0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell)
        y = y0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell)
        z = 0
        
        #Plotting the spheres in their location
        x0 = (r_ma* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_ma*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_ma*np.cos(phi)) + z
        mars = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color=color)
        


    
        ######################################################### Jupiter #########################################################
        
        #Initial planet information
        avgDist_j = 778600000
        r_j = 69911000
        semimajor_axis = 778570000
        semiminor_axis = 778064300
        color = 'orange'
        eccentricity = .0489
        a1 = 0                                                                                                          # Foci one x-coordinate
        b1 = 0                                                                                                          # Foci one y-coordinate
        b2 = 0                                                                                                          # Foci two y-coordinate
        x0 = semimajor_axis * eccentricity                                                # Center x-value
        y0 = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x0) - b2                                                                                     # Foci two x-coordinate
        phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                      # Angle between major axis and x-axis


        #Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_j = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        j_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)

                   #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/j_year
        
        #Getting individual orbital locations      
        x = x0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell)
        y = y0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell)
        z = 0
        
        #Plotting the spheres in their location
        x0 = (r_j* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_j*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_j*np.cos(phi)) + z
        jupiter = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color=color) 




        ######################################################### Saturn #########################################################

        #Initial planet information
        avgDist_s = 1433500000
        r_s = 58232000
        semimajor_axis = 1433530000
        semiminor_axis = 1431240077
        color = 'navajowhite'
        eccentricity = .0565
        a1 = 0                                                                                                          # Foci one x-coordinate
        b1 = 0                                                                                                          # Foci one y-coordinate
        b2 = 0                                                                                                          # Foci two y-coordinate
        x0 = semimajor_axis * eccentricity                                                # Center x-value
        y0 = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x0) - b2                                                                                     # Foci two x-coordinate
        phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                      # Angle between major axis and x-axis


       #Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_s = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        s_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)

                   #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/s_year
        
        #Getting individual orbital locations      
        x = x0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell)
        y = y0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell)
        z = 0

        #Plotting the spheres in their location
        x0 = (r_s* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_s*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_s*np.cos(phi)) + z
        saturn = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color=color)




        ######################################################### Uranus #########################################################

        #Initial planet information
        avgDist_u = 2872500000
        r_u = 25362000
        semimajor_axis = 2872460000
        semiminor_axis = 2866961900
        color = 'lightblue'
        eccentricity = .0457
        a1 = 0                                                                                                          # Foci one x-coordinate
        b1 = 0                                                                                                          # Foci one y-coordinate
        b2 = 0                                                                                                          # Foci two y-coordinate
        x0 = semimajor_axis * eccentricity                                                # Center x-value
        y0 = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x0) - b2                                                                                     # Foci two x-coordinate
        phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                      # Angle between major axis and x-axis


       #Orbital Velocity Calculation in megameters per second (using Kepler's Laws)
        vel_u = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        u_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)

                   #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/u_year
        
        #Getting individual orbital locations      
        x = x0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell)
        y = y0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell)
        z = 0

        #Plotting the spheres in their location
        x0 = (r_u* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_u*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_u*np.cos(phi)) + z
        uranus = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color=color)




        ######################################################### Neptune #########################################################
        
        #Initial planet information
        avgDist_n = 4495100000
        r_n = 24622000
        semimajor_axis = 4495060000
        semiminor_axis = 4499727700
        color = 'blue'
        eccentricity = .0113
        a1 = 0                                                                                                          # Foci one x-coordinate
        b1 = 0                                                                                                          # Foci one y-coordinate
        b2 = 0                                                                                                          # Foci two y-coordinate
        x0 = semimajor_axis * eccentricity                                                # Center x-value
        y0 = (b1 + b2) / 2                                                                                     # Center y-value
        a2 = (2 * x0) - b2                                                                                     # Foci two x-coordinate
        phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                      # Angle between major axis and x-axis


        #Orbital Velocity Calculation in kilometers per second (using Kepler's Laws)
        vel_n = (np.sqrt((grav_constant*mass_sun)*((2/(avgDist_v*1000)) - (1/(semiminor_axis*1000)))))/1000
        #Orbital Period calulation in years (using Kepler's Laws)
        n_year = ((2*np.pi)*np.sqrt(((semimajor_axis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)

                   #Scaling the orbital period in refernce to earths year
        scaled_per = earth_year/n_year 
        
        #Getting individual orbital locations      
        x = x0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell) - semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell)
        y = y0 + semimajor_axis * np.cos(2*np.pi*i/(N/scaled_per)) * np.sin(phi_ell) + semiminor_axis * np.sin(2*np.pi*i/(N/scaled_per)) * np.cos(phi_ell)
        z = 0

        #Plotting the spheres in their location
        x0 = (r_n* np.sin(phi) * np.cos(theta)) + x
        y0 = (r_n*np.sin(phi) * np.sin(theta)) + y
        z0 = (r_n*np.cos(phi)) + z
        neptune = ax.plot_surface(x0, y0, z0,  rstride=planet_smooth, cstride=planet_smooth, color=color) 


        
        plt.draw()
        plt.pause(.00000000000000000000000000001)
        ax.collections.remove(mercury)
        ax.collections.remove(venus)
        ax.collections.remove(earth)
        ax.collections.remove(mars)
        ax.collections.remove(jupiter)
        ax.collections.remove(saturn)
        ax.collections.remove(uranus)
        ax.collections.remove(neptune)

        min_range = i + 1
    calc_range = (min_range)*repeat_amount 

print("The Program has Finished")
