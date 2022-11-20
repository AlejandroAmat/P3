import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt("res_x.txt", dtype=float)
r = np.loadtxt("rex_r.txt", dtype=float)

fm = 20000
time = np.arange(0,len(x)).astype(float)
time = time/fm

muestras = np.arange(0,len(r)).astype(int)



plt.subplot(211)
plt.title("Segmento de Señal Sonoro")
plt.plot(time, x, linewidth =0.5)
plt.xlabel('Tiempo (s)')
plt.ylabel('Amplitud')
plt.grid(True)
plt.subplot(212)
plt.xlabel('Muestras')
plt.ylabel('Autocorrelación Calculada')
plt.plot(muestras, r, linewidth =0.5)
plt.grid(True)
plt.show()