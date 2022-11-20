import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt("hamming.txt", dtype=float)



m = np.arange(0,len(x)).astype(float)

plt.title("Filtro Hamming")
plt.plot(m ,x , linewidth =0.5)
plt.xlabel('Muestra')
plt.ylabel('Amplitud')
plt.grid(True)
plt.show()