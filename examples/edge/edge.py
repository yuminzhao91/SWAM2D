#!/bin/python2.7

import numpy as np

n1 = 401
n2 = 401

vp = np.zeros((n1, n2), dtype=np.float32, order='F')
vs = np.zeros((n1, n2), dtype=np.float32, order='F')
ro = np.zeros((n1, n2), dtype=np.float32, order='F')

vp[:, :] = 400.
vp[200:401,200:401] = 600.

vs[:, :] = 250.
vs[200:401,200:401] = 375.

ro[:, :] = 1500.
ro[200:401,200:401] = 2250.

fvp = open('fvp.bin', 'wb')
vp2 = vp.swapaxes(0,1)
vp2.tofile(fvp)
fvp.close()

fvs = open('fvs.bin', 'wb')
vs2 = vs.swapaxes(0,1)
vs2.tofile(fvs)
fvs.close()

fro = open('fro.bin', 'wb')
ro2 = ro.swapaxes(0,1)
ro2.tofile(fro)
fro.close()
