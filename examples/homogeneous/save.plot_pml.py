#encoding: utf-8
#!/bin/python2.7

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

# >> Read input binary files
dfig1 = np.fromfile('validation_pml0/snapp00801', dtype=np.float32)
dfig1 = dfig1.reshape(481, 481).swapaxes(1,0)

dfig2 = np.fromfile('validation_pml1/snapp00801', dtype=np.float32)
dfig2 = dfig2.reshape(481, 481).swapaxes(1,0)

dfig3 = np.fromfile('validation_pml0/snapp01201', dtype=np.float32)
dfig3 = dfig3.reshape(481, 481).swapaxes(1,0)

dfig4 = np.fromfile('validation_pml1/snapp01201', dtype=np.float32)
dfig4 = dfig4.reshape(481, 481).swapaxes(1,0)

dfig5 = np.fromfile('validation_pml0/snapp01601', dtype=np.float32)
dfig5 = dfig5.reshape(481, 481).swapaxes(1,0)

dfig6 = np.fromfile('validation_pml1/snapp01601', dtype=np.float32)
dfig6 = dfig6.reshape(481, 481).swapaxes(1,0)

font1 = plt.figure(figsize=(7.0,9.0))

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
fig1.set_title(r"\bf without PML")
fig1.text(5., 195., r"\bf t=0.2s", fontsize=10)
fig1.imshow(dfig1, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig1.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

fig2 = font1.add_subplot(3,2,2)
fig2.set_xlabel(r" ") #Distance [m]")
fig2.set_ylabel(r"Depth [m]")
fig2.set_xticks(x)
fig2.set_yticks(y)
fig2.tick_params(axis='x')
fig2.tick_params(axis='y')
fig2.text(5., 195., r"\bf t=0.2s", fontsize=10)
fig2.set_title(r"\bf with PML")

fig2.imshow(dfig2, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig2.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

fig3 = font1.add_subplot(3,2,3)
fig3.set_xlabel(r" ") #Distance [m]")
fig3.set_ylabel(r"Depth [m]")
fig3.set_xticks(x)
fig3.set_yticks(y)
fig3.tick_params(axis='x')
fig3.tick_params(axis='y')
fig3.text(5., 195., r"\bf t=0.3s", fontsize=10)
fig3.imshow(dfig3, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig3.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

fig4 = font1.add_subplot(3,2,4)
fig4.set_xlabel(r" ") #Distance [m]")
fig4.set_ylabel(r"Depth [m]")
fig4.set_xticks(x)
fig4.set_yticks(y)
fig4.tick_params(axis='x')
fig4.tick_params(axis='y')
fig4.text(5., 195., r"\bf t=0.3s", fontsize=10)
fig4.imshow(dfig4, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig4.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

fig5 = font1.add_subplot(3,2,5)
fig5.set_xlabel(r"Distance [m]")
fig5.set_ylabel(r"Depth [m]")
fig5.set_xticks(x)
fig5.set_yticks(y)
fig5.tick_params(axis='x')
fig5.tick_params(axis='y')
fig5.text(5., 195., r"\bf t=0.4s", fontsize=10)
fig5.imshow(dfig5, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig5.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

fig6 = font1.add_subplot(3,2,6)
fig6.set_xlabel(r"Distance [m]")
fig6.set_ylabel(r"Depth [m]")
fig6.set_xticks(x)
fig6.set_yticks(y)
fig6.tick_params(axis='x')
fig6.tick_params(axis='y')
fig6.text(5., 195., r"\bf t=0.4s", fontsize=10)
fig6.imshow(dfig6, cmap='gray', aspect='equal', extent=[-20,220,220,-20], vmin=-5.0e-3, vmax=5.0e-3)
fig6.plot(pmlbox[:,0], pmlbox[:,1], linestyle=':', linewidth=1, color="black")

font1.savefig('validation_pml.ps')
sp.call('ps2eps -l -B -s b0 -c -n -f validation_pml.ps', shell=True)
sp.call('rm -f validation_pml.ps', shell=True)
