#!/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt

# see PML Festa
d = 1.
n = 21
a = np.zeros(n, dtype=np.float32)

h = float(n-1)*d

for i in range(0, n):
    q = float(i)*d
    B = q/h
    a[i] = 5.*(300./(h*d))*np.power(B, 5)

plt.plot(a)
plt.show()
