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


f = mlab.figure()
ellipse_line_rad =  2
N=20
points_range = 25

phi = np.linspace(0, np.pi, points_range)
theta = np.linspace(0, 2*np.pi, points_range)
phi, theta = np.meshgrid(phi, theta)

t = np.linspace(0, 20, 200)
mlab.plot3d(ellipse_line_rad * np.sin(t), ellipse_line_rad * np.cos(t), 0*t, tube_radius=.05)

r = .5
x = r * np.sin(phi) * np.cos(theta) 
y = r * np.sin(phi) * np.sin(theta) 
z = r * np.cos(phi)
jupiter = mlab.mesh(x, y, z)

j_texture = vtk.vtkTexture()
j_textureReader = vtk.vtkJPEGReader()
j_textureReader.SetFileName(os.path.join(path, 'jupiter_scaled_flip.jpg'))
j_texture.SetInputConnection(j_textureReader.GetOutputPort())
get_texture(jupiter, j_texture)

@mlab.animate(delay=100)
def anim():
    i = 0
    while i < 100:
        time.sleep(.00000001)
        #s.mlab_source.set(x=5 * np.sin(phi) * np.cos(theta) + i, y=5 * np.sin(phi) * np.sin(theta)+ i, z=5 * np.cos(phi) + i)
        jupiter.mlab_source.set(x = r * np.sin(phi) * np.cos(theta) + (ellipse_line_rad * np.cos(2*np.pi*i/N)), y = r * np.sin(phi) * np.sin(theta) + (ellipse_line_rad * np.sin(2*np.pi*i/N)), z = r * np.cos(phi) + 0)
        i += 1
        yield

anim()
mlab.show()
