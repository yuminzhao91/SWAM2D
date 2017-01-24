import numpy as np
import matplotlib.pyplot as plt

n1 = 121
n2 = 601
h = 0.1

vp = np.zeros((n1, n2), dtype=np.float32)
vs = np.zeros((n1, n2), dtype=np.float32)
ro = np.zeros((n1, n2), dtype=np.float32)

for i1 in range(0, n1):
    z = float(n1-1-i1)*h
    if(z > 8.):
        vs[n1-1-i1,:] = 340.
        vp[n1-1-i1,:] = 340.*np.sqrt(3.)
    elif(z <=8. and z > 4.):
        vs[n1-1-i1,:] = 300.
        vp[n1-1-i1,:] = 300.*np.sqrt(3.)
    elif(z <=4. and z > 2.):
        vs[n1-1-i1,:] = 240.
        vp[n1-1-i1,:] = 240.*np.sqrt(3.)
    elif(z <= 2.):
        vs[n1-1-i1,:] = 200.
        vp[n1-1-i1,:] = 200.*np.sqrt(3.)

ro[:,:] = 1000.

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
