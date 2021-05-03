from mesh import Mesh
from solver import FDTD, Utilities, Source, Panel
from viewer import Animator
import copy
import numpy as np





malla1=Mesh(200,0.001,4,0.04,110,140)
malla2=Mesh(200,0.001,1,0,110,140)
pulso=Source('gauss',40,12,20)


et1k1= FDTD(malla1,pulso,5e-8).FDTDLoop(110,140)[0]
e2tk1= FDTD(malla2,pulso,5e-8).FDTDLoop(110,140)[0]
et1k2= FDTD(malla1,pulso,5e-8).FDTDLoop(110,140)[1]
e2tk2= FDTD(malla2,pulso,5e-8).FDTDLoop(110,140)[1]

ex_film=FDTD(malla1,pulso,5e-8).FDTDLoop(110,140)[2]
Animator().animationex(ex_film,malla1)

r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(5e-8,et1k1)

#Analytical Result
omega = np.linspace(1, 3e12, 100001) * 2 *np.pi
x = np.abs(Panel(3e-2,  4,   0.04).R(omega))
y = np.abs(Panel(3e-2,  4,   0.04).T(omega))



Animator().fftgraph(freq,r,t)
Animator().fftexact(omega, x, y)