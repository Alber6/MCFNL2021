import numpy as np
import matplotlib.pyplot as plt

s0 = 10e-6
t0 = 10*s0

t = np.linspace(0, 500e-6, num=(1000 + 1))
gauss = np.exp( - np.power(t - t0, 2) / 2 / s0**2)


plt.figure("Dominio del tiempo")
plt.plot(t,gauss)

plt.figure("Dominio de la frecuencia")
n=gauss.size
dt=t[1]-t[0]
f=np.fft.fftfreq(n, dt)
gausst=np.fft.fft(gauss)

plt.plot(f,np.absolute(gausst))


plt.show()
