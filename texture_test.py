from mayavi import mlab
import PIL
import vtk
import numpy as np
import os, os.path


path = "/Users/Austin Keller/Desktop/planet textures"
erad    = 5 # equatorial radius (meters)
prad    = 5 # polar radius (meters)
erot    = 7.2921158553e-5 # earth rotation rate (radians/sec)


def get_texture(name, texture):   
    if name == rings:
        name.actor.actor.mapper.scalar_visibility=False
        name.actor.enable_texture = True
        name.actor.tcoord_generator_mode = 'plane'
        name.actor.actor.texture = texture
    else: 
        name.actor.actor.mapper.scalar_visibility=False
        name.actor.enable_texture = True
        name.actor.tcoord_generator_mode = 'sphere'
        name.actor.actor.texture = texture
        cylinder_mapper = name.actor.tcoord_generator
        cylinder_mapper.prevent_seam = 0
    

u = np.linspace(0, 2 * np.pi, 180)
v = np.linspace(0, np.pi, 180)

# Rings
theta = np.linspace(0, 2.*np.pi, 100)
phi = np.linspace(0, 2*np.pi, 100)
theta, phi = np.meshgrid(theta, phi)
c, a = 1.9*prad, prad/2
x = (c + a*np.cos(theta)) * np.cos(phi)
y = (c + a*np.cos(theta)) * np.sin(phi)
z = 0 * np.sin(theta)

rings = mlab.mesh(x, y, z)

r_texture = vtk.vtkTexture()
r_textureReader = vtk.vtkJPEGReader()
r_textureReader.SetFileName(os.path.join(path, 'saturnringcolor.jpg'))
r_texture.SetInputConnection(r_textureReader.GetOutputPort())
get_texture(rings, r_texture)


### First Planet
x = erad * np.outer(np.cos(u), np.sin(v))
y = erad * np.outer(np.sin(u), np.sin(v))
z = prad * np.outer(np.ones_like(u), np.cos(v))
#mlab.figure(size=(800, 800), bgcolor=(0.16, 0.28, 0.46))

earth = mlab.mesh(x,y,z)

e_texture = vtk.vtkTexture()
e_textureReader = vtk.vtkJPEGReader()
e_textureReader.SetFileName(os.path.join(path, 'saturn_scaled.jpg'))
e_texture.SetInputConnection(e_textureReader.GetOutputPort())
get_texture(earth, e_texture)




### Second Planet
x0 = erad * np.outer(np.cos(u), np.sin(v)) + 30
y0 = erad * np.outer(np.sin(u), np.sin(v))  
z0 = prad * np.outer(np.ones_like(u), np.cos(v))

jupiter = mlab.mesh(x0,y0,z0)

j_texture = vtk.vtkTexture()
j_textureReader = vtk.vtkJPEGReader()
j_textureReader.SetFileName(os.path.join(path, 'neptune_scaled_rotate.jpg'))
j_texture.SetInputConnection(j_textureReader.GetOutputPort())
get_texture(jupiter, j_texture)


#mlab.view()
mlab.show()

