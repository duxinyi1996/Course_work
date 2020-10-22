import math
import matplotlib.pyplot as plt
totalx = 4
deltax = 0.0001
halfwidth = 1
num = totalx/deltax
deltanum = halfwidth/deltax
psi= [0]*int(num)
e =40
v = [0]*int(num)
for j in range (0,int(num)):
	v[j]=50*math.sin(j/num)
psi[0] = 0
psi[1] = 0.0001;
for i in range (1,int(num)-1):
	print (i)
	psi[i+1] =2*psi[i]-psi[i-1]+2*(math.pow(deltax,2))*(v[i]-e)*psi[i]
plt.plot(psi)
plt.show()