# Simple 3D Moecular Dynamics code (I wrote this for my term project)

```python
import numpy as np
import matplotlib.pyplot as plt
from random import *


Natom = 20
dens = 1.0  # density (1.2 for fcc)
t  = 0
itern = 1000
x = np.zeros((Natom),float) 	
y = np.zeros((Natom), float) 	
z = np.zeros((Natom),float)
vx = np.zeros((Natom), float) 	
vy = np.zeros((Natom), float) 	
vz = np.zeros((Natom),float)
fx = np.zeros((Natom,2),float) 	# x component of force
fy = np.zeros((Natom,2), float) 	# y component of force
fz = np.zeros((Natom,2),float)
PE1 = np.zeros((itern))
KE1 = np.zeros((itern))
Etot1 = np.zeros((itern))
time = np.arange((itern))

L = 8 		# side of square with atoms

def initialposvel(): 
	for i in range(0,Natom): 
		x[i] = randint(0,L) 
		y[i] = randint(0,L) 
                z[i] = randint(0,L)
		vx[i] = random() 
		vy[i] = random() 
                vz[i] = random()

def sign(a,b):
	if(b >= 0.0):
		return abs(a)
	else:
		return -abs(a)

def Forces(t,w,PE,PEorW):  
	r2cut = 6. 
	PE = 0.
	for i in range(0,Natom):
		fx[i][t] = fy[i][t] = fz[i][t] = 0.0 
	for i in range(0,Natom-1):
		for j in range(i+1,Natom):
			dx = x[i] - x[j] 
			dy = y[i] - y[j] 
                        dz = z[i] - z[j]
			if(abs(dx) > 0.50*L): 
				dx = dx - sign(L,dx) 
			if(abs(dy) > 0.50*L): 
				dy = dy - sign(L,dy)
                        if(abs(dz) > 0.50*L): 
				dz = dz - sign(L,dz)
			r2 = dx*dx + dy*dy + dz*dz
			if(r2 < r2cut): 
				if(r2 == 0.): 
					r2 = 0.0001
				invr2 = 1./r2 
				wij = 48.*(invr2**3 - 0.5)* invr2**3
				fijx = wij*invr2*dx
				fijy = wij*invr2*dy
                                fijz = wij*invr2*dz
				fx[i][t] = fx[i][t] + fijx
				fy[i][t] = fy[i][t] + fijy
				fz[i][t] = fz[i][t] + fijz
				fx[j][t] = fx[j][t] - fijx 
				fy[j][t] = fy[j][t] - fijy
				fz[j][t] = fz[j][t] - fijz 
				PE = PE + 4.*(invr2**3)*((invr2**3) - 1.)
				w = w + wij
	if(PEorW == 1):
		return PE
	else:
		return w


t1 = 0
PE = 0.0
h = 0.031 #step
hover2 = h/2.0
KE = 0.0
w = 0.0
initialposvel()
for i in range(0,Natom):
	KE = KE+ (vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i])/2.0
PE = Forces(t1,w,PE,1)
time =1

while t<itern:
	#rate(100)
	for i in range(0,Natom):
		PE = Forces(t1,w,PE,1)
		x[i] = x[i] + h*(vx[i] + hover2*fx[i][t1]) # velocity Verlet
		y[i] = y[i] + h*(vy[i] + hover2*fy[i][t1]) 
		z[i] = z[i] + h*(vz[i] + hover2*fz[i][t1]) 
		if x[i] <= 0.0:
			x[i] = x[i] + L 
		if x[i] >= L:
			x[i] = x[i] - L
		if y[i] <= 0.0:
			y[i] = y[i] + L
		if y[i] >= L:
			y[i] = y[i] - L
		if z[i] <= 0.0:
			z[i] = z[i] + L
		if z[i] >= L:
			z[i] = z[i] - L

	t2 =1
	PE = Forces(t2,w,PE,1)
        KE = 0
	for i in range(0,Natom):
		vx[i] = vx[i] + hover2*(fx[i][t1]+fx[i][t2])
		vy[i] = vy[i] + hover2*(fy[i][t1]+fy[i][t2])
		vz[i] = vz[i] + hover2*(fz[i][t1]+fz[i][t2])
		KE = KE + (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])/2

        Etot = PE+KE
        PE1[t] = PE
        KE1[t] = KE
        Etot1[t] = Etot
        t = t+1


time = np.arange(0,itern)*h
plt.figure()
plt.plot(time,PE1,label='PE')
plt.plot(time,KE1,label='KE')
plt.plot(time,Etot1,label='Total')
plt.legend()
plt.show()

```
