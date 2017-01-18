#!/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt

# see PML Festa
d = 0.0092
n = 11
a = np.zeros(n, dtype=np.float32)

h = float(n-1)*d

for i in range(0, n):
    q = float(i)*d
    B = q/h
    a[i] = 1.-5*(4./(h*d))*np.power(B, 5)

plt.plot(a)
plt.show()
