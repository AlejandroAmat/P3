import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt("res_x.txt", dtype=float)
x_no_clip = np.loadtxt("res_x_no_clipp.txt", dtype=float)

fm = 20000
time = np.arange(0,len(x)).astype(float)
time = time/fm




plt.subplot(211)
plt.title("Original vs Clipping")
plt.plot(time, x_no_clip, linewidth =0.5)
plt.xlabel('Tiempo (s)')
plt.ylabel('Amplitud')
plt.grid(True)
plt.subplot(212)
plt.xlabel('Tiempo(s)')
plt.ylabel('Amplitud')
plt.plot(time, x, linewidth =0.5)
plt.grid(True)
plt.show()