import math
import matplotlib as mpl
import matplotlib.pyplot as plt
step = 0.01
time = 10000
x =[0]
t =[0]
xx=[0]
for i in range(0,time):
	x +=[x[i]+(x[i]+1)*step]
	t +=[t[i]+ step]
	xx += [math.exp(t[i+1])-1]
	if (x[i]+1)>15:
		break

plt.plot(t,x,'+-',label="Euler method")
plt.plot(t,xx,label="Exact result")
plt.xlabel('time')
plt.ylabel('x')
plt.show()