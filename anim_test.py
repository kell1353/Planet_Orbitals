from mayavi import mlab
import numpy as np
import math
import time
import vtk
import os, os.path

path = "/Users/Austin Keller/Desktop/planet textures"

def get_texture(name, texture):   
        name.actor.actor.mapper.scalar_visibility=False
        name.actor.enable_texture = True
        name.actor.tcoord_generator_mode = 'sphere'
        name.actor.actor.texture = texture
        cylinder_mapper = name.actor.tcoord_generator
        cylinder_mapper.prevent_seam = 0

def draw_sun(r, x, y, z):
    x = (r* np.sin(phi) * np.cos(theta)) + x
    y = (r*np.sin(phi) * np.sin(theta)) + y
    z = (r*np.cos(phi)) + z
    planet = mlab.mesh(x, y, z, color=(1, 1, 0))              #(having higher strides allow the program to run faster)

def draw_sphere(e_r, p_r):
    global x; global y; global z
    x = (e_r* np.sin(phi) * np.cos(theta))
    y = (e_r*np.sin(phi) * np.sin(theta))
    z = (p_r*np.cos(phi))


f = mlab.figure()
ellipse_line_rad =  2
N=20
points_range = 25


""" Constant Variables """
earth_year = 1.0010
# Elliptical Constants
init_theta = 0
a1, b1, b2 = 0, 0, 0                                        # Foci one x-coordinate, # Foci one y-coordinate, # Foci two y-coordinate
# Spherical constants
phi = np.linspace(0, np.pi, points_range)
theta = np.linspace(0, 2*np.pi, points_range)
phi, theta = np.meshgrid(phi, theta)

"""Note the planet radii are scaled up by a factor of 10^3. The sun is just multiplied by 50."""
sun = draw_sun(1, 0, 0, 0) 

"""Draw Celestial Body Orbits"""
t = np.linspace(0, 20, 200)
mlab.plot3d(ellipse_line_rad * np.sin(t), ellipse_line_rad * np.cos(t), 0*t, tube_radius=.025)

##################################################### Earth ##################################################### 
# Functions
e_r = .5
draw_sphere(e_r, e_r)
earth = mlab.mesh(x, y, z)
#get_object_values(semimajor_axis, semiminor_axis, deg_inclination, 0.000, 6, 6)
# Texture variables
e_texture = vtk.vtkTexture()
e_textureReader = vtk.vtkJPEGReader()
e_textureReader.SetFileName(os.path.join(path, 'earth_scaled_flip.jpg'))
e_texture.SetInputConnection(e_textureReader.GetOutputPort())
get_texture(earth, e_texture)

##################################################### Jupiter ##################################################### 
j_r = .5
draw_sphere(j_r, j_r)
jupiter = mlab.mesh(x, y, z)

j_texture = vtk.vtkTexture()
j_textureReader = vtk.vtkJPEGReader()
j_textureReader.SetFileName(os.path.join(path, 'jupiter_scaled_flip.jpg'))
j_texture.SetInputConnection(j_textureReader.GetOutputPort())
get_texture(jupiter, j_texture)


@mlab.animate(delay=100)
def anim():
    #i = 0
    #while i < 100:
    for i in range(0, 101):
        if i > 0:
                earth_line.remove()
                jupiter_line.remove()
        
        time.sleep(.00000001)
        #plot_object(14, 14, .0167, 0.000, 6, 6)
        """ Earth """
        earth.mlab_source.set(x = e_r * np.sin(phi) * np.cos(theta) + (ellipse_line_rad * np.cos(2*np.pi*i/N)), y = e_r * np.sin(phi) * np.sin(theta) + (ellipse_line_rad * np.sin(2*np.pi*i/N)), z = e_r * np.cos(phi) + 0)
        x0 = (ellipse_line_rad * np.cos(2*np.pi*i/N))
        y0 = (ellipse_line_rad * np.sin(2*np.pi*i/N))
        z0 =  0
        earth_line = mlab.text3d(x0, y0, 1, 'Earth', scale=(.1, .1, .1))
        """ Jupiter """
        jupiter.mlab_source.set(x = j_r * np.sin(phi) * np.cos(theta) + (4 * np.cos(2*np.pi*i/N)), y = j_r * np.sin(phi) * np.sin(theta) + (4 * np.sin(2*np.pi*i/N)), z = j_r * np.cos(phi) + 0)
        x0 = (4 * np.cos(2*np.pi*i/N))
        y0 = (4 * np.sin(2*np.pi*i/N))
        z0 =  0
        jupiter_line = mlab.text3d(x0, y0, 1, 'Jupiter', scale=(.1, .1, .1))
        #i += 1
        yield

anim()
mlab.show()
