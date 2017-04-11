#!/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import subprocess as sp

# >> Postscript output parametrization using LaTeX
plt.rcParams['text.latex.preamble']=[r"\usepackage{times} \usepackage{txfonts} \usepackage[french]{babel} \RequirePackage[utf8]{inputenc}"]

# >> Parametrization options
params = {'backend': 'PS',
          'font.family' : 'serif',
          'text.usetex' : True,
          'font.size' : 11,
          'font.weight' : 'bold',
          'axes.labelsize' : 11,
          'xtick.labelsize' : 11,
          'ytick.labelsize' : 11,
          }

# >> Apply parametrization
plt.rcParams.update(params)

tmp1 = np.loadtxt('Uz_file_an')
tmp2 = np.fromfile('recz.bin', dtype=np.float32)
tmp2 = tmp2.reshape(4, 2001).swapaxes(1,0)

perc=20.
r1 = 25.
r2 = 50.
r3 = 75.
r4 = 100.

t  = tmp1[:,0]
a1 = tmp1[:,1]/np.amax(np.abs(tmp1[:,1]))*perc+r1
a2 = tmp1[:,2]/np.amax(np.abs(tmp1[:,1]))*perc+r2
a3 = tmp1[:,3]/np.amax(np.abs(tmp1[:,1]))*perc+r3
a4 = tmp1[:,4]/np.amax(np.abs(tmp1[:,1]))*perc+r4

w1 = -1.*tmp2[1:,0]/np.amax(np.abs(tmp2[:,0]))*perc+r1
w2 = -1.*tmp2[1:,1]/np.amax(np.abs(tmp2[:,0]))*perc+r2
w3 = -1.*tmp2[1:,2]/np.amax(np.abs(tmp2[:,0]))*perc+r3
w4 = -1.*tmp2[1:,3]/np.amax(np.abs(tmp2[:,0]))*perc+r4

tmp1x = np.loadtxt('Ux_file_an')
tmp2x = np.fromfile('recx.bin', dtype=np.float32)
tmp2x = tmp2x.reshape(4, 2001).swapaxes(1,0)

tx  = tmp1x[:,0]
a1x = tmp1x[:,1]/np.amax(np.abs(tmp1x[:,1]))*perc+r1
a2x = tmp1x[:,2]/np.amax(np.abs(tmp1x[:,1]))*perc+r2
a3x = tmp1x[:,3]/np.amax(np.abs(tmp1x[:,1]))*perc+r3
a4x = tmp1x[:,4]/np.amax(np.abs(tmp1x[:,1]))*perc+r4

w1x = tmp2x[1:,0]/np.amax(np.abs(tmp2x[:,0]))*perc+r1
w2x = tmp2x[1:,1]/np.amax(np.abs(tmp2x[:,0]))*perc+r2
w3x = tmp2x[1:,2]/np.amax(np.abs(tmp2x[:,0]))*perc+r3
w4x = tmp2x[1:,3]/np.amax(np.abs(tmp2x[:,0]))*perc+r4

font1 = plt.figure(figsize=(6.0,9.0))

subplotsize=[2.,2.]
figuresize=[7.,9.] 
left = 0.1*(1.-subplotsize[0]/figuresize[0])
right = 1.-left
bottom = 0.1*(1.-subplotsize[1]/figuresize[1])
top = 1.-bottom

fig1 = font1.add_subplot(3,2,1, aspect='auto')
font1.subplots_adjust(left=left,right=right,bottom=bottom,top=top)

y = ([25, 50., 75., 100.])
fig1.set_ylabel(r"Offset [m]")
fig1.set_xlabel(r"Time [s]")
fig1.set_ylim(0.,120.)
fig1.set_xlim(0.,0.5)
#fig1.set_xticks(x)
fig1.set_yticks(y)
#fig1.tick_params(axis='x')
#fig1.tick_params(axis='y')
fig1.text(0.01, 110.0, r"\bf horizontal component", fontsize=10)

fig1.plot(t, a1, color='lightgray', linewidth=4)
fig1.plot(t, a2, color='lightgray', linewidth=4)
fig1.plot(t, a3, color='lightgray', linewidth=4)
fig1.plot(t, a4, color='lightgray', linewidth=4)
fig1.plot(t, w1, color='black', linewidth=1)
fig1.plot(t, w2, color='black', linewidth=1)
fig1.plot(t, w3, color='black', linewidth=1)
fig1.plot(t, w4, color='black', linewidth=1)
#fig1.set_aspect('equal')

fig2 = font1.add_subplot(3,2,2)
#fig2.set_ylabel(r"Offset [m]")
fig2.set_xlabel(r"Time [s]")
#fig2.set_xticks(x)
fig2.set_ylim(0.,120.)
fig2.set_yticks(y)
#fig2.tick_params(axis='x')
#fig2.tick_params(axis='y')
#fig2.text(5., 195., r"\bf t=0.2s", fontsize=10)
fig2.text(0.01, 110.0, r"\bf vertical component", fontsize=10)

fig2.plot(t, a1x, color='lightgray', linewidth=4)
fig2.plot(t, a2x, color='lightgray', linewidth=4)
fig2.plot(t, a3x, color='lightgray', linewidth=4)
fig2.plot(t, a4x, color='lightgray', linewidth=4)
fig2.plot(t, w1x, color='black', linewidth=1)
fig2.plot(t, w2x, color='black', linewidth=1)
fig2.plot(t, w3x, color='black', linewidth=1)
fig2.plot(t, w4x, color='black', linewidth=1)

font1.savefig('validation_trac_pml.ps')
sp.call('ps2eps -l -B -s b0 -c -n -f validation_trac_pml.ps', shell=True)
sp.call('rm -f validation_trac_pml.ps', shell=True)
