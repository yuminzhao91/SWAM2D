import numpy as np
import matplotlib.pyplot as plt

tmp1 = np.loadtxt('Uz_file_an')
tmp2 = np.fromfile('recz.bin', dtype=np.float32)
tmp2 = tmp2.reshape(4, 2001).swapaxes(1,0)

t  = tmp1[:,0]
a1 = tmp1[:,1]/np.amax(np.abs(tmp1[:,1]))
a2 = tmp1[:,2]/np.amax(np.abs(tmp1[:,1]))
a3 = tmp1[:,3]/np.amax(np.abs(tmp1[:,1]))
a4 = tmp1[:,4]/np.amax(np.abs(tmp1[:,1]))

w1 = -1.*tmp2[1:,0]/np.amax(np.abs(tmp2[:,0]))
w2 = -1.*tmp2[1:,1]/np.amax(np.abs(tmp2[:,0]))
w3 = -1.*tmp2[1:,2]/np.amax(np.abs(tmp2[:,0]))
w4 = -1.*tmp2[1:,3]/np.amax(np.abs(tmp2[:,0]))

tmp1x = np.loadtxt('Ux_file_an')
tmp2x = np.fromfile('recx.bin', dtype=np.float32)
tmp2x = tmp2x.reshape(4, 2001).swapaxes(1,0)

tx  = tmp1x[:,0]
a1x = tmp1x[:,1]/np.amax(np.abs(tmp1x[:,1]))
a2x = tmp1x[:,2]/np.amax(np.abs(tmp1x[:,1]))
a3x = tmp1x[:,3]/np.amax(np.abs(tmp1x[:,1]))
a4x = tmp1x[:,4]/np.amax(np.abs(tmp1x[:,1]))

w1x = tmp2x[1:,0]/np.amax(np.abs(tmp2x[:,0]))
w2x = tmp2x[1:,1]/np.amax(np.abs(tmp2x[:,0]))
w3x = tmp2x[1:,2]/np.amax(np.abs(tmp2x[:,0]))
w4x = tmp2x[1:,3]/np.amax(np.abs(tmp2x[:,0]))

plt.subplot(421)
plt.plot(a1, color='lightgray', linewidth=8)
plt.plot(w1, color='black', linewidth=1)

plt.subplot(423)
plt.plot(a2, color='lightblue', linewidth=6)
plt.plot(w2, color='black', linewidth=1)

plt.subplot(425)
plt.plot(a3, color='lightblue', linewidth=6)
plt.plot(w3, color='black', linewidth=1)

plt.subplot(427)
plt.plot(a4, color='lightblue', linewidth=6)
plt.plot(w4, color='black', linewidth=1)

plt.subplot(422)
plt.plot(a1x, color='lightblue', linewidth=6)
plt.plot(w1x, color='black', linewidth=1)

plt.subplot(424)
plt.plot(a2x, color='lightblue', linewidth=6)
plt.plot(w2x, color='black', linewidth=1)

plt.subplot(426)
plt.plot(a3x, color='lightblue', linewidth=6)
plt.plot(w3x, color='black', linewidth=1)

plt.subplot(428)
plt.plot(a4x, color='lightblue', linewidth=8)
plt.plot(w4x, color='black', linewidth=1)

plt.show()

print np.argmax(np.abs(a1)), np.argmax(np.abs(w1))
