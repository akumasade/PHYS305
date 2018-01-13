import numpy as np
import matplotlib.pyplot as plt

n = 1000
xmin = 0.0
xmax = 10.
x = np.linspace(xmin, xmax, n+1)
y = np.cos(x)
z = np.sin(x)
w = y**2 - z**2

plt.figure()
plt.plot(x, y, 'r', label='cos x')
plt.plot(x, z, 'b', label='sin x')
plt.plot(x, w, 'g--', label='cos$^2$ x - sin$^2$ x')
plt.xlim(-1, 12)
plt.ylim(-1.75, 1.2)
plt.xlabel('x value')
plt.ylabel('y value')
plt.title('Simple plot operations')
plt.legend(loc='lower left')
plt.savefig('plot1.pdf')
plt.show()
