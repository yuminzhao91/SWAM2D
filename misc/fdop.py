#!/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt

n = 101
dx = 2.*np.pi/float(n)
f = np.zeros(n, dtype=np.float32)
g = np.zeros(n, dtype=np.float32)
h1 = np.zeros(n, dtype=np.float32)
h2 = np.zeros(n, dtype=np.float32)

for i in range(0,n):
    x = (-1.*np.pi)+float(i)/float(n)*(2.*np.pi)
    f[i] = np.sin(x)
    g[i] = np.cos(x)

for i in range(0,n-1):
    h1[i] = (f[i]-f[i-1])/dx
    h2[i] = (f[i+1]-f[i])/dx
    
print np.amax(h1), 1./np.amax(h1)
d = np.linspace(-np.pi, np.pi, n)
plt.plot(d, f)
plt.plot(d, g)
plt.plot(d, h1)
plt.plot(d, h2)
plt.show()

