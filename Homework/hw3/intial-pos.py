import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n = 50
np.random.seed(12345)
r = (np.random.random(n))**(1./3)
theta = np.arccos(2*np.random.random(n)-1)
phi = 2*np.pi*np.random.random(n)
x = r*np.sin(theta)*np.cos(phi)
y = r*np.sin(theta)*np.sin(phi)
z = r*np.cos(theta)
pos = zip(x,y,z)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, marker="o")

plt.show()
