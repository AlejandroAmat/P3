import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt("prueba.f0", dtype=float)
x_ideal = np.loadtxt("prueba.f0ref", dtype=float)


time = np.arange(0,len(x)).astype(float)
time = time*0.015




plt.subplot(211)
plt.title("Comparativa Estimación Ideal")
plt.plot(time, x_ideal, linewidth =0.5)
plt.xlabel('Tiempo (s)')
plt.ylabel('Frec(Hz)')
plt.grid(True)
plt.subplot(212)
plt.xlabel('Tiempo(s)')
plt.ylabel('Frec(Hz)')
plt.plot(time, x, linewidth =0.5)
plt.grid(True)
plt.show()