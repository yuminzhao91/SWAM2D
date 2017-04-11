#encoding: utf-8
#!/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('PS')
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

# >> Read input binary files
dfig1 = np.fromfile('validation_pml/snapp00801', dtype=np.float32)
dfig1 = dfig1.reshape(481, 481).swapaxes(1,0)

dfig2 = np.fromfile('validation_pml/snapp01201', dtype=np.float32)
dfig2 = dfig2.reshape(481, 481).swapaxes(1,0)

dfig3 = np.fromfile('validation_pml/snapp01601', dtype=np.float32)
dfig3 = dfig3.reshape(481, 481).swapaxes(1,0)

dfig4 = np.fromfile('validation_pml/snapp02001', dtype=np.float32)
dfig4 = dfig4.reshape(481, 481).swapaxes(1,0)

font1 = plt.figure(figsize=(6.0,9.0))

subplotsize=[2.,2.]
figuresize=[7.,9.] 
left = 0.1*(1.-subplotsize[0]/figuresize[0])
right = 1.-left
bottom = 0.1*(1.-subplotsize[1]/figuresize[1])
top = 1.-bottom

fig1 = font1.add_subplot(3,2,1)
font1.subplots_adjust(left=left,right=right,bottom=bottom,top=top)

x = []
y = []
for i in range(0, 11):
    x.append(float(i-1)*40.)
    y.append(float(i-1)*40.)

pmlbox = np.zeros((5,2), dtype=np.float32)
pmlbox[0,0] = 0.
pmlbox[0,1] = 0.
pmlbox[1,0] = 200.
pmlbox[1,1] = 0.
pmlbox[2,0] = 200.
pmlbox[2,1] = 200.
pmlbox[3,0] = 0.
pmlbox[3,1] = 200.
pmlbox[4,0] = 0.
pmlbox[4,1] = 0.

fig1.set_xlabel(r" ") #Distance [m]")
fig1.set_ylabel(r"Depth [m]")
fig1.set_xticks(x)
fig1.set_yticks(y)
fig1.tick_params(axis='x')
fig1.tick_params(axis='y')
fig1.text(5., 195., r"\bf t=0.2s", fontsize=10)
fig1.imshow(dfig1, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig1.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

fig2 = font1.add_subplot(3,2,2)
fig2.set_xlabel(r" ") #Distance [m]")
fig2.set_ylabel(r" ") #Depth [m]")
fig2.set_xticks(x)
fig2.set_yticks(y)
fig2.tick_params(axis='x')
fig2.tick_params(axis='y')
fig2.text(5., 195., r"\bf t=0.3s", fontsize=10)

fig2.imshow(dfig2, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig2.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

fig3 = font1.add_subplot(3,2,3)
fig3.set_xlabel(r"Distance [m]")
fig3.set_ylabel(r"Depth [m]")
fig3.set_xticks(x)
fig3.set_yticks(y)
fig3.tick_params(axis='x')
fig3.tick_params(axis='y')
fig3.text(5., 195., r"\bf t=0.4s", fontsize=10)
fig3.imshow(dfig3, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig3.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

fig4 = font1.add_subplot(3,2,4)
fig4.set_xlabel(r"Distance [m]")
fig4.set_ylabel(r" ") #Depth [m]")
fig4.set_xticks(x)
fig4.set_yticks(y)
fig4.tick_params(axis='x')
fig4.tick_params(axis='y')
fig4.text(5., 195., r"\bf t=0.5s", fontsize=10)
fig4.imshow(dfig4, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig4.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

font1.savefig('validation_snap_pml.ps', transparent=True)
sp.call('ps2eps -l -B -s b0 -c -n -f validation_snap_pml.ps', shell=True)
sp.call('rm -f validation_snap_pml.ps', shell=True)

# >> TRACES

tmp1 = np.loadtxt('validation_pml/Uz_file_an')
tmp2 = np.fromfile('validation_pml/recz.bin', dtype=np.float32)
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

tmp1x = np.loadtxt('validation_pml/Ux_file_an')
tmp2x = np.fromfile('validation_pml/recx.bin', dtype=np.float32)
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

font1.savefig('validation_trac_pml.ps', transparent=True)
sp.call('ps2eps -l -B -s b0 -c -n -f validation_trac_pml.ps', shell=True)
sp.call('rm -f validation_trac_pml.ps', shell=True)
