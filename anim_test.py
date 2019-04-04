from mayavi import mlab
import numpy as np
import math
import time


N=20

points_range = 25
phi = np.linspace(0, np.pi, points_range)
theta = np.linspace(0, 2*np.pi, points_range)
phi, theta = np.meshgrid(phi, theta)

f = mlab.figure()

r = .5

t = np.linspace(0, 20, 200)
mlab.plot3d(np.sin(t), np.cos(t), 0*t, tube_radius=.05)

x = r * np.sin(phi) * np.cos(theta) 
y = r * np.sin(phi) * np.sin(theta) 
z = r * np.cos(phi)
s = mlab.mesh(x, y, z)

@mlab.animate(delay=100)
def anim():
    i = 0
    while i < 100:
        time.sleep(.00000001)
        #s.mlab_source.set(x=5 * np.sin(phi) * np.cos(theta) + i, y=5 * np.sin(phi) * np.sin(theta)+ i, z=5 * np.cos(phi) + i)
        s.mlab_source.set(x = r * np.sin(phi) * np.cos(theta) + np.cos(2*np.pi*i/N), y = r * np.sin(phi) * np.sin(theta) + np.sin(2*np.pi*i/N), z = r * np.cos(phi) + 0)
        i += 1
        yield

anim()
mlab.show()
