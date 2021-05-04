from mesh import Mesh
from solver import FDTD, Utilities, Source, Panel
from viewer import Animator
import copy
import numpy as np


sigma=0.04
permittivity=4

ddx=0.001
start_m=110
end_m=140

thickness=(end_m-start_m)*ddx
time=5e-8


malla1=Mesh(200,ddx,permittivity,sigma,start_m,end_m)
malla2=Mesh(200,ddx,1,0,start_m,end_m)
pulso=Source('gauss',40,12,20)


et1k1, et1k2, ex_film = FDTD(malla1,pulso,time).FDTDLoop(110,140)
e2tk1, e2tk2, _= FDTD(malla2,pulso,time).FDTDLoop(110,140)

Animator().animationex(ex_film,malla1)


#FFT results  

r= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[0]
t= Utilities().FFT(et1k1,e2tk1,et1k2,e2tk2)[1]
freq=Utilities().frequency(time,et1k1)

#Analytical Result

omega = np.linspace(1, 3e12, 100001) * 2 *np.pi
x = np.abs(Panel(thickness,  permittivity,   sigma).R(omega))
y = np.abs(Panel(thickness,  permittivity,   sigma).T(omega))



Animator().fftgraph(freq,r,t)
Animator().fftexact(omega, x, y)
