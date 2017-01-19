#!/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt

cp = 600.
R = 0.1

npml = 21
dpml = 0.5

n = 5

L = float(npml-1)*dpml
d0 = 3.*cp*np.log(1./R)/(2.*L)

d = np.zeros(npml, dtype=np.float32)

for ipml in range(0, npml):
    q = float(ipml)*dpml
    d[ipml] = d0*(q/L)**n

plt.plot(d)
plt.show()
