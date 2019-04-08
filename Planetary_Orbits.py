from mayavi import mlab
import numpy as np
import math
import time
import vtk
import os, os.path

repeat_amount = int(input("How many years would you like to see?: "))
fig = mlab.figure('Solar System', bgcolor = (0,0,0), size = (600,400))

path = "/Users/Austin Keller/Desktop/planet textures"

def get_texture(name, texture):   
    name.actor.actor.mapper.scalar_visibility=False
    name.actor.enable_texture = True
    name.actor.tcoord_generator_mode = 'sphere'
    name.actor.actor.texture = texture
    cylinder_mapper = name.actor.tcoord_generator
    cylinder_mapper.prevent_seam = 0

def get_flat_texture(name, texture):
    name.actor.actor.mapper.scalar_visibility=False
    name.actor.enable_texture = True
    name.actor.tcoord_generator_mode = 'plane'
    name.actor.actor.texture = texture

def draw_sun(r, x, y, z):
    global star
    x = (r* np.sin(phi) * np.cos(theta)) + x
    y = (r*np.sin(phi) * np.sin(theta)) + y
    z = (r*np.cos(phi)) + z
    star = mlab.mesh(x, y, z)              #(having higher strides allow the program to run faster)

def draw_sphere(e_r, p_r):
    global x; global y; global z
    x = (e_r* np.sin(phi) * np.cos(theta))
    y = (e_r*np.sin(phi) * np.sin(theta))
    z = (p_r*np.cos(phi))

def elliptical_orbit_line(semimajor_axis, semiminor_axis, eccentricity, deg_inclination, rad):
    # Compute ellipse parameters
    a1, b1, b2 = 0, 0, 0                                                                                               # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
    x0 = semimajor_axis * eccentricity                                                            # Center x-value
    y0 = (b1 + b2) / 2                                                                                                 # Center y-value
    f = np.sqrt((a1 - x0)**2 + (b1 - y0)**2)                                                       # Distance from center to focus
    a2 = (2 * x0) - b2                                                                                                  # Foci two x-coordinate
    phi_ell = np.arctan2((b1 - b2), (a1 - a2))                                                   # Angle between major axis and x-axis
    # Parametric plot in t
    z_list = []
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
    mlab.plot3d(x, y, z_list, color=(1,1,1), tube_radius = rad)


""" Constant Variables """
"""Note the planet radii are scaled up by a factor of 10^3. The sun is just multiplied by 50."""
grav_constant =  6.67*(10**(-11))
mass_sun = 1.98855*(10**30)
earth_year = 1.0010
# Elliptical Constants
init_theta = 0
a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
t = np.linspace(0, 20, 200)
# Spherical constants
points_range = 50
phi = np.linspace(0, 2*np.pi, points_range)
theta = np.linspace(0, 2*np.pi, points_range)
phi, theta = np.meshgrid(phi, theta)


##################################################### Sun ##################################################### 
"""Note the planet radii are scaled up by a factor of 10^3. The sun is just multiplied by 50."""
sun = draw_sun(695956*50, 0, 0, 0)
sun_texture = vtk.vtkTexture()
sun_textureReader = vtk.vtkJPEGReader()
sun_textureReader.SetFileName(os.path.join(path, 'sun_scaled.jpg'))
sun_texture.SetInputConnection(sun_textureReader.GetOutputPort())
get_texture(star, sun_texture)


##################################################### Mercury ##################################################### 
# Planetary variables
me_degInclination = 7.000
me_semimajorAxis = 57910000; me_semiminorAxis = 56672817
me_eqRad = 2439700; me_polarRad = 2439700
# Elliptical variables
me_eccentricity = .2056
me_xCenter = me_semimajorAxis * me_eccentricity                                    # Center x-value
me_yCenter = (b1 + b2) / 2                                                                                         # Center y-value
me_a2 = (2 * me_xCenter) - b2                                                                                 # Foci two x-coordinate
me_planetYear = ((2*np.pi)*np.sqrt(((me_semimajorAxis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
me_scaledPer = earth_year/me_planetYear
# Planetary Functions
draw_sphere(me_eqRad, me_polarRad)
mercury = mlab.mesh(x, y, z)
mercury_orbit = elliptical_orbit_line(me_semimajorAxis, me_semiminorAxis, me_eccentricity, me_degInclination, 1000000/2)
# Texture Wrapping
me_texture = vtk.vtkTexture()
me_textureReader = vtk.vtkJPEGReader()
me_textureReader.SetFileName(os.path.join(path, 'mercury_scaled.jpg'))
me_texture.SetInputConnection(me_textureReader.GetOutputPort())
get_texture(mercury, me_texture)

##################################################### Venus ##################################################### 
# Planetary variables
v_degInclination = 3.390
v_semimajorAxis = 108210000; v_semiminorAxis = 108207571
v_eqRad = 6051800; v_polarRad = 6051800
# Elliptical variables
v_eccentricity = .0067
v_xCenter = v_semimajorAxis * v_eccentricity                                         # Center x-value
v_yCenter = (b1 + b2) / 2                                                                                     # Center y-value
v_a2 = (2 * v_xCenter) - b2                                                                                 # Foci two x-coordinate
v_planetYear = ((2*np.pi)*np.sqrt(((v_semimajorAxis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
v_scaledPer = earth_year/v_planetYear
# Planetary Functions
draw_sphere(v_eqRad, v_polarRad)
venus = mlab.mesh(x, y, z)
venus_orbit = elliptical_orbit_line(v_semimajorAxis, v_semiminorAxis, v_eccentricity, v_degInclination, 1000000/2)
# Texture Wrapping
v_texture = vtk.vtkTexture()
v_textureReader = vtk.vtkJPEGReader()
v_textureReader.SetFileName(os.path.join(path, 'venus_atmosphere_scaled.jpg'))
v_texture.SetInputConnection(v_textureReader.GetOutputPort())
get_texture(venus, v_texture)

##################################################### Earth ##################################################### 
# Planetary variables
e_degInclination = 0.000
e_semimajorAxis = 149600000; e_semiminorAxis = 149579138
e_eqRad = 6378100; e_polarRad = 6356100
# Elliptical variables
e_eccentricity = .0167
e_xCenter = e_semimajorAxis * e_eccentricity                                        # Center x-value
e_yCenter = (b1 + b2) / 2                                                                                     # Center y-value
e_a2 = (2 * e_xCenter) - b2                                                                                 # Foci two x-coordinate
e_planetYear = ((2*np.pi)*np.sqrt(((e_semimajorAxis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
e_scaledPer = earth_year/e_planetYear
# Planetary Functions
draw_sphere(e_eqRad, e_polarRad)
earth = mlab.mesh(x, y, z)
earth_orbit = elliptical_orbit_line(e_semimajorAxis, e_semiminorAxis, e_eccentricity, e_degInclination, 1000000/2)
# Texture Wrapping
e_texture = vtk.vtkTexture()
e_textureReader = vtk.vtkJPEGReader()
e_textureReader.SetFileName(os.path.join(path, 'earth_scaled_flip.jpg'))
e_texture.SetInputConnection(e_textureReader.GetOutputPort())
get_texture(earth, e_texture)

##################################################### Mars ##################################################### 
# Planetary variables
ma_degInclination = 1.850
ma_semimajorAxis = 227920000; ma_semiminorAxis = 226921546
ma_eqRad = 3396200; ma_polarRad = 3376200
# Elliptical variables
ma_eccentricity = .0935
ma_xCenter = ma_semimajorAxis * ma_eccentricity                                     # Center x-value
ma_yCenter = (b1 + b2) / 2                                                                                         # Center y-value
ma_a2 = (2 * ma_xCenter) - b2                                                                                 # Foci two x-coordinate
ma_planetYear = ((2*np.pi)*np.sqrt(((ma_semimajorAxis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
ma_scaledPer = earth_year/ma_planetYear
# Planetary Functions
draw_sphere(ma_eqRad, ma_polarRad)
mars = mlab.mesh(x, y, z)
mars_orbit = elliptical_orbit_line(ma_semimajorAxis, ma_semiminorAxis, ma_eccentricity, ma_degInclination, 1000000/2)
# Texture Wrapping
ma_texture = vtk.vtkTexture()
ma_textureReader = vtk.vtkJPEGReader()
ma_textureReader.SetFileName(os.path.join(path, 'mars.jpg'))
ma_texture.SetInputConnection(ma_textureReader.GetOutputPort())
get_texture(mars, ma_texture)

##################################################### Jupiter ##################################################### 
# Planetary variables
j_degInclination = 1.304
j_semimajorAxis = 778570000; j_semiminorAxis = 777638581
j_eqRad = 71492000; j_polarRad = 66854000
# Elliptical variables
j_eccentricity = .0489
j_xCenter = j_semimajorAxis * j_eccentricity                                        # Center x-value
j_yCenter = (b1 + b2) / 2                                                                                  # Center y-value
j_a2 = (2 * j_xCenter) - b2                                                                               # Foci two x-coordinate
j_planetYear = ((2*np.pi)*np.sqrt(((j_semimajorAxis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
j_scaledPer = earth_year/j_planetYear
# Planetary Functions
draw_sphere(j_eqRad, j_polarRad)
jupiter = mlab.mesh(x, y, z)
jupiter_orbit = elliptical_orbit_line(j_semimajorAxis, j_semiminorAxis, j_eccentricity, j_degInclination, 1000000/2)
# Texture Wrapping
j_texture = vtk.vtkTexture()
j_textureReader = vtk.vtkJPEGReader()
j_textureReader.SetFileName(os.path.join(path, 'jupiter_scaled_flip.jpg'))
j_texture.SetInputConnection(j_textureReader.GetOutputPort())
get_texture(jupiter, j_texture)

##################################################### Saturn ##################################################### 
# Planetary variables
s_degInclination = 2.485
s_semimajorAxis = 1433530000; s_semiminorAxis = 1431240078
s_eqRad = 60268000; s_polarRad = 54364000
# Elliptical variables
s_eccentricity = .0565
s_xCenter = s_semimajorAxis * s_eccentricity                                        # Center x-value
s_yCenter = (b1 + b2) / 2                                                                                   # Center y-value
s_a2 = (2 * s_xCenter) - b2                                                                                # Foci two x-coordinate
s_planetYear = ((2*np.pi)*np.sqrt(((s_semimajorAxis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
s_scaledPer = earth_year/s_planetYear
# Saturn's Rings
c, a = 1.9*s_polarRad, s_eqRad/2
ring_x = (c + a*np.cos(theta)) * np.cos(phi) 
ring_y = (c + a*np.cos(theta)) * np.sin(phi) 
ring_z = 0 * np.sin(theta) 
# Planetary Functions
draw_sphere(s_eqRad, s_polarRad)
saturn = mlab.mesh(x, y, z)
saturn_orbit = elliptical_orbit_line(s_semimajorAxis, s_semiminorAxis, s_eccentricity, s_degInclination, 1000000/2)
# Saturn's Rings
rings = mlab.mesh(x, y, z)
# Texture Wrapping
s_texture = vtk.vtkTexture()
s_textureReader = vtk.vtkJPEGReader()
s_textureReader.SetFileName(os.path.join(path, 'saturn_scaled.jpg'))
s_texture.SetInputConnection(s_textureReader.GetOutputPort())
get_texture(saturn, s_texture)
# Saturn's Rings Texture Wrapping
r_texture = vtk.vtkTexture()
r_textureReader = vtk.vtkJPEGReader()
r_textureReader.SetFileName(os.path.join(path, 'saturnringcolor.jpg'))
r_texture.SetInputConnection(r_textureReader.GetOutputPort())
get_flat_texture(rings, r_texture)

##################################################### Uranus ##################################################### 
# Planetary variables
u_degInclination = 0.772
u_semimajorAxis = 2872460000; u_semiminorAxis = 2869458880
u_eqRad = 25559000; u_polarRad = 24973000
# Elliptical variables
u_eccentricity = .0457
u_xCenter = u_semimajorAxis * u_eccentricity                                      # Center x-value
u_yCenter = (b1 + b2) / 2                                                                                   # Center y-value
u_a2 = (2 * u_xCenter) - b2                                                                               # Foci two x-coordinate
u_planetYear = ((2*np.pi)*np.sqrt(((u_semimajorAxis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
u_scaledPer = earth_year/u_planetYear
# Planetary Functions
draw_sphere(u_eqRad, u_polarRad)
uranus = mlab.mesh(x, y, z)
uranus_orbit = elliptical_orbit_line(u_semimajorAxis, u_semiminorAxis, u_eccentricity, u_degInclination, 1000000/2)
# Texture Wrapping
u_texture = vtk.vtkTexture()
u_textureReader = vtk.vtkJPEGReader()
u_textureReader.SetFileName(os.path.join(path, 'uranus_scaled.jpg'))
u_texture.SetInputConnection(u_textureReader.GetOutputPort())
get_texture(uranus, u_texture)

##################################################### Neptune ##################################################### 
# Planetary variables
n_degInclination = 1.769
n_semimajorAxis = 4495060000; n_semiminorAxis = 4494773004
n_eqRad = 24764000; n_polarRad = 24341000
# Elliptical variables
n_eccentricity = .0113
n_xCenter = n_semimajorAxis * n_eccentricity                                      # Center x-value
n_yCenter = (b1 + b2) / 2                                                                                   # Center y-value
n_a2 = (2 * n_xCenter) - b2                                                                               # Foci two x-coordinate
n_planetYear = ((2*np.pi)*np.sqrt(((n_semimajorAxis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
n_scaledPer = earth_year/n_planetYear
# Planetary Functions
draw_sphere(n_eqRad, n_polarRad)
neptune = mlab.mesh(x, y, z)
neptune_orbit = elliptical_orbit_line(n_semimajorAxis, n_semiminorAxis, n_eccentricity, n_degInclination, 1000000/2)
# Texture Wrapping
n_texture = vtk.vtkTexture()
n_textureReader = vtk.vtkJPEGReader()
n_textureReader.SetFileName(os.path.join(path, 'neptune_scaled_rotate.jpg'))
n_texture.SetInputConnection(n_textureReader.GetOutputPort())
get_texture(neptune, n_texture)

##################################################### Pluto ##################################################### 
# Planetary variables
pl_degInclination = 17.16
pl_semimajorAxis = 5906380000; pl_semiminorAxis = 5720653186
pl_eqRad = 1187000; pl_polarRad = 1187000
# Elliptical variables
pl_eccentricity = .2488
pl_xCenter = pl_semimajorAxis * pl_eccentricity                                     # Center x-value
pl_yCenter = (b1 + b2) / 2                                                                                     # Center y-value
pl_a2 = (2 * pl_xCenter) - b2                                                                               # Foci two x-coordinate
pl_planetYear = ((2*np.pi)*np.sqrt(((pl_semimajorAxis*1000)**3)/(grav_constant*mass_sun)))/(60*60*24*365)
pl_scaledPer = earth_year/pl_planetYear
# Planetary Functions
draw_sphere(pl_eqRad, pl_polarRad)
pluto = mlab.mesh(x, y, z)
pluto_orbit = elliptical_orbit_line(pl_semimajorAxis, pl_semiminorAxis, pl_eccentricity, pl_degInclination, 1000000/2)
# Texture Wrapping
pl_texture = vtk.vtkTexture()
pl_textureReader = vtk.vtkJPEGReader()
pl_textureReader.SetFileName(os.path.join(path, 'pluto_scaled.jpg'))
pl_texture.SetInputConnection(pl_textureReader.GetOutputPort())
get_texture(pluto, pl_texture)


""" Animation constants"""
N = 30
scaling_factor = 20000000
calc_range = ((N + 1) * repeat_amount)

"""Draw planet as spheres at specific points"""
@mlab.animate(delay = 50)
def anim():
    for i in range(0, calc_range):
        if i > 0:
                mercury_line.remove(), venus_line.remove(), earth_line.remove()
                mars_line.remove()#, jupiter_line.remove(), saturn_line.remove()
                ##uranus_line.remove(), neptune_line.remove(), pluto_line.remove()
                pluto_line.remove()
                
        time.sleep(.0000001)
        # Show the amount of years have passed during the orbits
        year_calc = round(i/N, 2)
        #sub_title = mlab.text3d(0, 0, 40000000, str(year_calc) + ' Earth years', color = (1,1,1), scale = .7)
        #print(i)
        """ Mercury """
        me_x = (me_xCenter + me_semimajorAxis * np.cos(2*np.pi*i/(N/me_scaledPer)) * np.cos(np.radians(init_theta)) - me_semiminorAxis * np.sin(2*np.pi*i/(N/me_scaledPer))* np.sin(np.radians(init_theta)))*np.cos(np.radians(me_degInclination))
        me_y = (me_yCenter + me_semimajorAxis * np.cos(2*np.pi*i/(N/me_scaledPer)) * np.sin(np.radians(init_theta)) + me_semiminorAxis * np.sin(2*np.pi*i/(N/me_scaledPer)) * np.cos(np.radians(init_theta)))
        r1 = np.sqrt(((me_x - 0)**2))
        if me_x > 0:
            r1 = -(np.sqrt(((me_x - 0)**2)))
        me_z = r1*np.tan(np.radians(me_degInclination))
        # Planent Text
        mercury_line = mlab.text3d(me_x, me_y, me_z + (me_polarRad*10), 'Mercury', scale=(scaling_factor , scaling_factor , scaling_factor))
        # Setting x, y, z values
        mercury.mlab_source.set(x = me_eqRad * np.sin(phi) * np.cos(theta) + (me_x), y = me_eqRad * np.sin(phi) * np.sin(theta) + (me_y), z = me_polarRad * np.cos(phi) + (me_z))
        """ Venus """
        v_x = (v_xCenter + v_semimajorAxis * np.cos(2*np.pi*i/(N/v_scaledPer)) * np.cos(np.radians(init_theta)) - v_semiminorAxis * np.sin(2*np.pi*i/(N/v_scaledPer))* np.sin(np.radians(init_theta)))*np.cos(np.radians(v_degInclination))
        v_y = (v_yCenter + v_semimajorAxis * np.cos(2*np.pi*i/(N/v_scaledPer)) * np.sin(np.radians(init_theta)) + v_semiminorAxis * np.sin(2*np.pi*i/(N/v_scaledPer)) * np.cos(np.radians(init_theta)))
        r1 = np.sqrt(((v_x - 0)**2))
        if v_x > 0:
            r1 = -(np.sqrt(((v_x - 0)**2)))
        v_z = r1*np.tan(np.radians(v_degInclination))
        # Planent Text
        venus_line = mlab.text3d(v_x, v_y, v_z + (v_polarRad*7), 'Venus', scale=(scaling_factor , scaling_factor , scaling_factor))
        # Setting x, y, z values
        venus.mlab_source.set(x = v_eqRad * np.sin(phi) * np.cos(theta) + (v_x), y = v_eqRad * np.sin(phi) * np.sin(theta) + (v_y), z = v_polarRad * np.cos(phi) + (v_z))
        """ Earth """
        e_x = (e_xCenter + e_semimajorAxis * np.cos(2*np.pi*i/(N/e_scaledPer)) * np.cos(np.radians(init_theta)) - e_semiminorAxis * np.sin(2*np.pi*i/(N/e_scaledPer))* np.sin(np.radians(init_theta)))*np.cos(np.radians(e_degInclination))
        e_y = (e_yCenter + e_semimajorAxis * np.cos(2*np.pi*i/(N/e_scaledPer)) * np.sin(np.radians(init_theta)) + e_semiminorAxis * np.sin(2*np.pi*i/(N/e_scaledPer)) * np.cos(np.radians(init_theta)))
        r1 = np.sqrt(((e_x - 0)**2))
        if e_x > 0:
            r1 = -(np.sqrt(((e_x - 0)**2)))
        e_z = r1*np.tan(np.radians(e_degInclination))
        # Planent Text
        earth_line = mlab.text3d(e_x, e_y, e_z + (e_polarRad*7), 'Earth', scale=(scaling_factor , scaling_factor , scaling_factor))
        # Setting x, y, z values
        earth.mlab_source.set(x = e_eqRad * np.sin(phi) * np.cos(theta) + (e_x), y = e_eqRad * np.sin(phi) * np.sin(theta) + (e_y), z = e_polarRad * np.cos(phi) + (e_z))
        """ Mars """
        ma_x = (ma_xCenter + ma_semimajorAxis * np.cos(2*np.pi*i/(N/ma_scaledPer)) * np.cos(np.radians(init_theta)) - ma_semiminorAxis * np.sin(2*np.pi*i/(N/ma_scaledPer))* np.sin(np.radians(init_theta)))*np.cos(np.radians(ma_degInclination))
        ma_y = (ma_yCenter + ma_semimajorAxis * np.cos(2*np.pi*i/(N/ma_scaledPer)) * np.sin(np.radians(init_theta)) + ma_semiminorAxis * np.sin(2*np.pi*i/(N/ma_scaledPer)) * np.cos(np.radians(init_theta)))
        r1 = np.sqrt(((ma_x - 0)**2))
        if ma_x > 0:
            r1 = -(np.sqrt(((ma_x - 0)**2)))
        ma_z = r1*np.tan(np.radians(ma_degInclination))
        # Planent Text
        mars_line = mlab.text3d(ma_x, ma_y, ma_z + (ma_polarRad*10), 'Mars', scale=(scaling_factor , scaling_factor , scaling_factor))
        # Setting x, y, z values
        mars.mlab_source.set(x = ma_eqRad * np.sin(phi) * np.cos(theta) + (ma_x), y = ma_eqRad * np.sin(phi) * np.sin(theta) + (ma_y), z = ma_polarRad * np.cos(phi) + (ma_z))
        """ Jupiter """
        j_x = (j_xCenter + j_semimajorAxis * np.cos(2*np.pi*i/(N/j_scaledPer)) * np.cos(np.radians(init_theta)) - j_semiminorAxis * np.sin(2*np.pi*i/(N/j_scaledPer))* np.sin(np.radians(init_theta)))*np.cos(np.radians(j_degInclination))
        j_y = (j_yCenter + j_semimajorAxis * np.cos(2*np.pi*i/(N/j_scaledPer)) * np.sin(np.radians(init_theta)) + j_semiminorAxis * np.sin(2*np.pi*i/(N/j_scaledPer)) * np.cos(np.radians(init_theta)))
        r1 = np.sqrt(((j_x - 0)**2))
        if j_x > 0:
            r1 = -(np.sqrt(((j_x - 0)**2)))
        j_z = r1*np.tan(np.radians(j_degInclination))
        # Planent Text
        ##jupiter_line = mlab.text3d(j_x, j_y, j_z + (j_polarRad*2), 'Jupiter', scale=(scaling_factor , scaling_factor , scaling_factor))
        # Setting x, y, z values
        jupiter.mlab_source.set(x = j_eqRad * np.sin(phi) * np.cos(theta) + (j_x), y = j_eqRad * np.sin(phi) * np.sin(theta) + (j_y), z = j_polarRad * np.cos(phi) + (j_z))
        """ Saturn """
        s_x = (s_xCenter + s_semimajorAxis * np.cos(2*np.pi*i/(N/s_scaledPer)) * np.cos(np.radians(init_theta)) - s_semiminorAxis * np.sin(2*np.pi*i/(N/s_scaledPer))* np.sin(np.radians(init_theta)))*np.cos(np.radians(s_degInclination))
        s_y = (s_yCenter + s_semimajorAxis * np.cos(2*np.pi*i/(N/s_scaledPer)) * np.sin(np.radians(init_theta)) + s_semiminorAxis * np.sin(2*np.pi*i/(N/s_scaledPer)) * np.cos(np.radians(init_theta)))
        r1 = np.sqrt(((s_x - 0)**2))
        if s_x > 0:
            r1 = -(np.sqrt(((s_x - 0)**2)))
        s_z = r1*np.tan(np.radians(s_degInclination))
        # Planent Text
        ##saturn_line = mlab.text3d(s_x, s_y, s_z + (s_polarRad*2), 'Saturn', scale=(scaling_factor , scaling_factor , scaling_factor))
        # Setting x, y, z values
        saturn.mlab_source.set(x = s_eqRad * np.sin(phi) * np.cos(theta) + (s_x), y = s_eqRad * np.sin(phi) * np.sin(theta) + (s_y), z = s_polarRad * np.cos(phi) + (s_z))
        # Saturn's Rings
        rings.mlab_source.set(x = (c + a*np.cos(theta)) * np.cos(phi)  + (s_x), y = (c + a*np.cos(theta)) * np.sin(phi)  + (s_y), z = 0 * np.sin(theta) + (s_z))        
        """ Uranus """
        u_x = (u_xCenter + u_semimajorAxis * np.cos(2*np.pi*i/(N/u_scaledPer)) * np.cos(np.radians(init_theta)) -u_semiminorAxis * np.sin(2*np.pi*i/(N/u_scaledPer))* np.sin(np.radians(init_theta)))*np.cos(np.radians(u_degInclination))
        u_y = (u_yCenter + u_semimajorAxis * np.cos(2*np.pi*i/(N/u_scaledPer)) * np.sin(np.radians(init_theta)) + u_semiminorAxis * np.sin(2*np.pi*i/(N/u_scaledPer)) * np.cos(np.radians(init_theta)))
        r1 = np.sqrt(((u_x - 0)**2))
        if u_x > 0:
            r1 = -(np.sqrt(((u_x - 0)**2)))
        u_z = r1*np.tan(np.radians(u_degInclination))
        # Planent Text
        ##uranus_line = mlab.text3d(u_x, u_y, u_z + (u_polarRad*3), 'Uranus', scale=(scaling_factor , scaling_factor , scaling_factor))
        # Setting x, y, z values
        uranus.mlab_source.set(x = u_eqRad * np.sin(phi) * np.cos(theta) + (u_x), y = u_eqRad * np.sin(phi) * np.sin(theta) + (u_y), z = u_polarRad * np.cos(phi) + (u_z))
        """ Neptune """
        n_x = (n_xCenter + n_semimajorAxis * np.cos(2*np.pi*i/(N/n_scaledPer)) * np.cos(np.radians(init_theta)) - n_semiminorAxis * np.sin(2*np.pi*i/(N/n_scaledPer))* np.sin(np.radians(init_theta)))*np.cos(np.radians(n_degInclination))
        n_y = (n_yCenter + n_semimajorAxis * np.cos(2*np.pi*i/(N/n_scaledPer)) * np.sin(np.radians(init_theta)) + n_semiminorAxis * np.sin(2*np.pi*i/(N/n_scaledPer)) * np.cos(np.radians(init_theta)))
        r1 = np.sqrt(((n_x - 0)**2))
        if n_x > 0:
            r1 = -(np.sqrt(((n_x - 0)**2)))
        n_z = r1*np.tan(np.radians(n_degInclination))
        # Planent Text
##         neptune_line = mlab.text3d(n_x, n_y, n_z + (n_polarRad*3), 'Neptune', scale=(scaling_factor , scaling_factor , scaling_factor))
        # Setting x, y, z values
        neptune.mlab_source.set(x = n_eqRad * np.sin(phi) * np.cos(theta) + (n_x), y = n_eqRad * np.sin(phi) * np.sin(theta) + (n_y), z = n_polarRad * np.cos(phi) + (n_z))
        """ Pluto """
        pl_x = (pl_xCenter + pl_semimajorAxis * np.cos(2*np.pi*i/(N/pl_scaledPer)) * np.cos(np.radians(init_theta)) - pl_semiminorAxis * np.sin(2*np.pi*i/(N/pl_scaledPer))* np.sin(np.radians(init_theta)))*np.cos(np.radians(pl_degInclination))
        pl_y = (pl_yCenter + pl_semimajorAxis * np.cos(2*np.pi*i/(N/pl_scaledPer)) * np.sin(np.radians(init_theta)) + pl_semiminorAxis * np.sin(2*np.pi*i/(N/pl_scaledPer)) * np.cos(np.radians(init_theta)))
        r1 = np.sqrt(((pl_x - 0)**2))
        if pl_x > 0:
            r1 = -(np.sqrt(((pl_x - 0)**2)))
        pl_z = r1*np.tan(np.radians(pl_degInclination))
        # Planent Text
        pluto_line = mlab.text3d(pl_x, pl_y, pl_z + (pl_polarRad*100), 'Pluto', scale=(scaling_factor , scaling_factor , scaling_factor))
        # Setting x, y, z values
        pluto.mlab_source.set(x = pl_eqRad * np.sin(phi) * np.cos(theta) + (pl_x), y = pl_eqRad * np.sin(phi) * np.sin(theta) + (pl_y), z = pl_polarRad * np.cos(phi) + (pl_z))
        yield

anim()
mlab.show()
