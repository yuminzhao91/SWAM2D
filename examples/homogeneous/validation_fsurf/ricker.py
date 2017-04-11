#!/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt

f1 = np.fromfile('fricker_fwm', dtype=np.float32)
f2 = np.fromfile('fricker.bin', dtype=np.float32)

plt.plot(f1, color='black')
plt.plot(f2[0:1024], color='red')
plt.show()
