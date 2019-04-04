import numpy as np
import PIL
import os, os.path
from mayavi import mlab
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import vtk

repeat_amount = float(input("How many years would you like to see?: "))
fig = mlab.figure(bgcolor=(1,1,1),size = (600,400))

# Setting the size of the window for the plot
##fig = plt.figure(figsize=(8,5))
##ax = fig.add_subplot(111, projection='3d')
##plt.title('Solar System', fontsize = 14, color = 'white')
##fig.subplots_adjust(top=1,bottom=0,left=0,right=1)

# The Cartesian coordinates of the unit sphere     
points_range = 100
phi = np.linspace(0, 2 * np.pi, points_range)
theta = np.linspace(0, 2*np.pi, points_range)
phi, theta = np.meshgrid(phi, theta)


def draw_sphere(r, x, y, z, c):
    x = (r* np.sin(phi) * np.cos(theta)) + x
    y = (r*np.sin(phi) * np.sin(theta)) + y
    z = (r*np.cos(phi)) + z
    planet = mlab.mesh(x, y, z, color=(1, 1, 0))              #(having higher strides allow the program to run faster)

    
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
    mlab.plot3d(x, y, z_list, color=(0,0,0), tube_radius=1000000/2)


path = "/Users/Austin Keller/Desktop/planet textures"
u = np.linspace(0, 2 * np.pi, 180)
v = np.linspace(0, np.pi, 180)

def get_texture(name, texture):
    name.actor.actor.mapper.scalar_visibility=False
    name.actor.enable_texture = True
    name.actor.tcoord_generator_mode = 'sphere'
    name.actor.actor.texture = texture
    cylinder_mapper = name.actor.tcoord_generator
    cylinder_mapper.prevent_seam = 0

def get_ring_texture(name, texture):
    name.actor.actor.mapper.scalar_visibility=False
    name.actor.enable_texture = True
    name.actor.tcoord_generator_mode = 'plane'
    name.actor.actor.texture = texture


"""Mercury"""
me_texture = vtk.vtkTexture()
me_textureReader = vtk.vtkJPEGReader()
me_textureReader.SetFileName(os.path.join(path, 'mercury_scaled.jpg'))
me_texture.SetInputConnection(me_textureReader.GetOutputPort())
"""Venus"""
v_texture = vtk.vtkTexture()
v_textureReader = vtk.vtkJPEGReader()
v_textureReader.SetFileName(os.path.join(path, 'venus_atmosphere_scaled.jpg'))
v_texture.SetInputConnection(v_textureReader.GetOutputPort())
"""Earth"""
e_texture = vtk.vtkTexture()
e_textureReader = vtk.vtkJPEGReader()
e_textureReader.SetFileName(os.path.join(path, 'earth_scaled_flip.jpg'))
e_texture.SetInputConnection(e_textureReader.GetOutputPort())
"""Mars"""
ma_texture = vtk.vtkTexture()
ma_textureReader = vtk.vtkJPEGReader()
ma_textureReader.SetFileName(os.path.join(path, 'mars.jpg'))
ma_texture.SetInputConnection(ma_textureReader.GetOutputPort())
"""Jupiter"""
j_texture = vtk.vtkTexture()
j_textureReader = vtk.vtkJPEGReader()
j_textureReader.SetFileName(os.path.join(path, 'jupiter_scaled_flip.jpg'))
j_texture.SetInputConnection(j_textureReader.GetOutputPort())
"""Saturn"""
s_texture = vtk.vtkTexture()
s_textureReader = vtk.vtkJPEGReader()
s_textureReader.SetFileName(os.path.join(path, 'saturn_scaled.jpg'))
s_texture.SetInputConnection(s_textureReader.GetOutputPort())
r_texture = vtk.vtkTexture()
r_textureReader = vtk.vtkJPEGReader()
r_textureReader.SetFileName(os.path.join(path, 'saturnringcolor.jpg'))
r_texture.SetInputConnection(r_textureReader.GetOutputPort())
"""Uranus"""
u_texture = vtk.vtkTexture()
u_textureReader = vtk.vtkJPEGReader()
u_textureReader.SetFileName(os.path.join(path, 'uranus_scaled.jpg'))
u_texture.SetInputConnection(u_textureReader.GetOutputPort())
"""Neptune"""
n_texture = vtk.vtkTexture()
n_textureReader = vtk.vtkJPEGReader()
n_textureReader.SetFileName(os.path.join(path, 'neptune_scaled_rotate.jpg'))
n_texture.SetInputConnection(n_textureReader.GetOutputPort())
"""Pluto"""
pl_texture = vtk.vtkTexture()
pl_textureReader = vtk.vtkJPEGReader()
pl_textureReader.SetFileName(os.path.join(path, 'pluto_scaled.jpg'))
pl_texture.SetInputConnection(pl_textureReader.GetOutputPort())




def plot_object(name, semimajor_axis, semiminor_axis, eccentricity, deg_inclination, eq_rad, polar_rad):
        global mercury; global me_text; global mercury_source
        global venus; global v_text; global venus_source
        global earth; global e_text; global earth_source
        global mars; global ma_text; global mars_source
        global jupiter; global j_text; global jupiter_source
        global saturn; global rings; global s_text; global saturn_source
        global uranus; global u_text; global uranus_source
        global neptune; global n_text; global neptune_source
        global pluto; global pluto_text; global pluto_source
        global x0; global y0; global z0
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
        
        if name == "Mercury":
            x0 = eq_rad * np.outer(np.cos(u), np.sin(v)) + x
            y0 = eq_rad * np.outer(np.sin(u), np.sin(v)) + y
            z0 = polar_rad * np.outer(np.ones_like(u), np.cos(v)) + z
            mercury = mlab.mesh(x0, y0, z0)
            get_texture(mercury, me_texture)
            mercury_source =  mercury.mlab_source
##            me_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Venus":
            x0 = eq_rad * np.outer(np.cos(u), np.sin(v)) + x
            y0 = eq_rad * np.outer(np.sin(u), np.sin(v)) + y
            z0 = polar_rad * np.outer(np.ones_like(u), np.cos(v)) + z
            venus = mlab.mesh(x0, y0, z0)
            get_texture(venus, v_texture)
            venus_source =  venus.mlab_source
##            v_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Earth":
            x0 = eq_rad * np.outer(np.cos(u), np.sin(v)) + x
            y0 = eq_rad * np.outer(np.sin(u), np.sin(v)) + y
            z0 = polar_rad * np.outer(np.ones_like(u), np.cos(v)) + z
            earth = mlab.mesh(x0, y0, z0)
            get_texture(earth, e_texture)
            earth_source =  earth.mlab_source
##            e_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Mars":
            x0 = eq_rad * np.outer(np.cos(u), np.sin(v)) + x
            y0 = eq_rad * np.outer(np.sin(u), np.sin(v)) + y
            z0 = polar_rad * np.outer(np.ones_like(u), np.cos(v)) + z
            mars = mlab.mesh(x0, y0, z0)
            get_texture(mars, ma_texture)
            mars_source =  mars.mlab_source
##            ma_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Jupiter":
            x0 = eq_rad * np.outer(np.cos(u), np.sin(v)) + x
            y0 = eq_rad * np.outer(np.sin(u), np.sin(v)) + y
            z0 = polar_rad * np.outer(np.ones_like(u), np.cos(v)) + z
            jupiter = mlab.mesh(x0, y0, z0)
            get_texture(jupiter, j_texture)
            jupiter_source =  jupiter.mlab_source
##            j_text = ax.text(x, y, z + planet_rad*2, name, color='white', ha = "center")
        elif name == "Saturn":
            x0 = eq_rad * np.outer(np.cos(u), np.sin(v)) + x
            y0 = eq_rad * np.outer(np.sin(u), np.sin(v)) + y
            z0 = polar_rad * np.outer(np.ones_like(u), np.cos(v)) + z
            saturn = mlab.mesh(x0, y0, z0)
            get_texture(saturn, s_texture)
            # Saturn's Rings
            c, a = 1.9*polar_rad, eq_rad/2
            xr = (c + a*np.cos(theta)) * np.cos(phi) + x
            yr = (c + a*np.cos(theta)) * np.sin(phi) + y
            zr = 0 * np.sin(theta) + z
            rings = mlab.mesh(xr, yr, zr)
            get_ring_texture(rings, r_texture)
            saturn_source =  saturn.mlab_source
##            s_text = ax.text(x, y, z + planet_rad*3, name, color='white', ha = "center")
        elif name == "Uranus":
            x0 = eq_rad * np.outer(np.cos(u), np.sin(v)) + x
            y0 = eq_rad * np.outer(np.sin(u), np.sin(v)) + y
            z0 = polar_rad * np.outer(np.ones_like(u), np.cos(v)) + z
            uranus = mlab.mesh(x0, y0, z0)
            get_texture(uranus, u_texture)
            uranus_source =  uranus.mlab_source
##            u_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Neptune":
            x0 = eq_rad * np.outer(np.cos(u), np.sin(v)) + x
            y0 = eq_rad * np.outer(np.sin(u), np.sin(v)) + y
            z0 = polar_rad * np.outer(np.ones_like(u), np.cos(v)) + z
            neptune = mlab.mesh(x0, y0, z0)
            get_texture(neptune, n_texture)
            neptune_source =  neptune.mlab_source
##            n_text = ax.text(x, y, z + planet_rad*10, name, color='white', ha = "center")
        elif name == "Pluto":
            x0 = eq_rad * np.outer(np.cos(u), np.sin(v)) + x
            y0 = eq_rad * np.outer(np.sin(u), np.sin(v)) + y
            z0 = polar_rad * np.outer(np.ones_like(u), np.cos(v)) + z
            pluto = mlab.mesh(x0, y0, z0)
            get_texture(pluto, pl_texture)
            pluto_source =  pluto.mlab_source
##            pluto_text = ax.text(x, y, z + planet_rad*100, name, color='white', ha = "center")


            

######################################## End of Definitions ########################################




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


"""Draw planet as spheres at specific points"""
"""Note the planet radii are scaled up by a factor of 10^3. The sun is just multiplied by 50."""
grav_constant =  6.67*(10**(-11))
mass_sun = 1.98855*(10**30)
earth_year = 1.0010

sun = draw_sphere(695956*50, 0, 0, 0, 'yellow')                     #####figure out how to scale this when you zoom in to graph.

""" Starting parameters """
N = 25                            # The lower the value the faster and more jumpy the planets orbit. (25 is ideal)
min_range = 0
repeat = 0
calc_range = N


while repeat < repeat_amount:
    repeat = repeat + 1
    # A way to plot a planet on a given planets orbital line
    for i in range(min_range, calc_range):
        
        # Show the amount of years have passed during the orbits
        year_calc = round(i/N, 1)
        #fig.suptitle(str(year_calc) + ' Earth years' , y = .8, fontsize = 9, color = 'white')
        
        mercury_planet =  plot_object("Mercury", 57910000, 56672817, .2056, 7.000, 2439700, 2439700)
        venus_planet =  plot_object("Venus", 108210000, 108207571, .0067, 3.390, 6051800, 6051800)
        earth_planet =  plot_object("Earth", 149600000, 149579138, .0167, 0.000, 6378100, 6356100)
        mars_planet =  plot_object("Mars", 227920000, 226921546, .0935, 1.850, 3396200, 3376200)
        jupiter_planet =  plot_object("Jupiter", 778570000, 777638581, .0489, 1.304, 71492000, 66854000)
        saturn_planet =  plot_object("Saturn", 1433530000, 1431240078, .0565, 2.485, 60268000, 54364000)
        uranus_planet =  plot_object("Uranus", 2872460000, 2869458880, .0457, .772, 25559000, 24973000)
        neptune_planet =  plot_object("Neptune", 4495060000, 4494773004, .0113, 1.769, 24764000, 24341000)
        pluto_planet =  plot_object("Pluto", 5906380000, 5720653186, .2488, 17.16, 1187000, 1187000)
##        ##pl_line, = plt.plot([x, x], [y, y], [z + (r_pl * 2), z + r_pl*500], linewidth=.25, color = "white")
        min_range = i + 1
        if year_calc == repeat_amount:
            cam = mercury.scene.camera
            cam.zoom(3)
            mlab.show()           
        else:
##            mercury_source.reset(x0=x0, y0=y0, z0=z0)
##            venus_source.reset(x0=x0, y0=y0, z0=z0), earth_source.reset(x0=x0, y0=y0, z0=z0)
##            mars_source.reset(x0=x0, y0=y0, z0=z0), jupiter_source.reset(x0=x0, y0=y0, z0=z0), saturn_source.reset(x0=x0, y0=y0, z0=z0)
##            uranus_source.reset(x0=x0, y0=y0, z0=z0), neptune_source.reset(x0=x0, y0=y0, z0=z0), pluto_source.reset(x0=x0, y0=y0, z0=z0)
            #mlab.clf()
            ######################## this is where some slowdown happens ############################
##            plt.pause(.000000000000000000000000000000000000000000000000000000001)
            # Removing the previously plotted celestial objects
            mercury.remove(), venus.remove(), earth.remove()
            mars.remove(), jupiter.remove(), saturn.remove(), rings.remove()
            uranus.remove(), neptune.remove(), pluto.remove()
            # Removing the previously drawn labels and lines
##            me_text.remove(), v_text.remove(), e_text.remove()
##            ma_text.remove(), j_text.remove(), s_text.remove()
##            u_text.remove(), n_text.remove(), pluto_text.remove()  #ax.lines.remove(pl_line)
            mlab.gcf()
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

