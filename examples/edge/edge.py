#!/bin/python2.7

"""
Generate input parameter file, acquisition file and direct access binary
model files to run with SWAM2D.
"""

def genparam():
    """
    Generate input parameter file for SWAM2D.
    """
    fparam = open('param.input', 'w')
    fparam.write('#[run]\n')
    fparam.write('test\n')
    fparam.write('0.5 0.00025\n')
    fparam.write('#[materials]\n')
    fparam.write('fvp.bin fvs.bin fro.bin\n')
    fparam.write('401 401 0.5\n')
    fparam.write('#[boundaries]\n')
    fparam.write('0 20 0.7\n')
    fparam.write('#[source]\n')
    fparam.write('2 2 2\n')
    fparam.write('30. 0.05\n')
    fparam.write('100. 50.\n') 
    fparam.write('#[receiver]\n')
    fparam.write('facqui.ascii\n')
    fparam.write('0.00025\n')
    fparam.close()

def genacqui(x0, z0, dx, dz, n):
    """
    Genrate acquisition file for SWAM2D.
    """
    facqui = open('facqui.ascii', 'w')
    for i in range(0, n):
        xr = x0+float(i)*dx
        zr = z0+float(i)*dz
        facqui.write(str(xr)+' '+str(zr)+'\n')
    facqui.close()
    
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
