import numpy as np
import matplotlib.pyplot as plt

tmp1 = np.loadtxt('Ux_file_an')
tmp2 = np.fromfile('recx.bin', dtype=np.float32)
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

plt.subplot(411)
plt.plot(a1, color='lightblue', linewidth=4)
plt.plot(w1, color='black', linewidth=1)

plt.subplot(412)
plt.plot(a2, color='lightblue', linewidth=4)
plt.plot(w2, color='black', linewidth=1)

plt.subplot(413)
plt.plot(a3, color='lightblue', linewidth=4)
plt.plot(w3, color='black', linewidth=1)

plt.subplot(414)
plt.plot(a4, color='lightblue', linewidth=4)
plt.plot(w4, color='black', linewidth=1)

plt.show()

print np.argmax(np.abs(a1)), np.argmax(np.abs(w1))
