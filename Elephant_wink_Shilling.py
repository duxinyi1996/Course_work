#!/usr/bin/python
import numpy as np 
from matplotlib import pyplot as plt 
from matplotlib import animation 

# elephant parameters
p1, p2, p3, p4 = (50 - 30j, 18 + 8j, 12 - 10j, -14 - 60j)
p5 = 1 + 20j # eyepiece


def fourier(t, C):
	f = np.zeros(t.shape)
	A, B = C.real, C.imag
	for k in range(len(C)):
		f = f + A[k]*np.cos(k*t)+ B[k]*np.sin(k*t)
	return f

# First set up the figure, the axis, and the plot element we want to animate 
fig = plt.figure() 
ax = plt.axes(xlim=(-80, 100), ylim=(-80, 100)) 
line, = ax.plot([], [], '.') 

# initialization function: plot the background of each frame 
def init(): 
	line.set_data([], []) 
	return line, 
	
# animation function. This is called sequentially 
# note: i is framenumber 
def animate(i): 
	npar = 6
	t = np.linspace(0, 2.*np.pi, 1000)
	Cx = np.zeros((npar,), dtype='complex')
	Cy = np.zeros((npar,), dtype='complex')
		
	Cx[1] = p1.real*1j
	Cx[2] = p2.real*1j
	Cx[3] = p3.real
	Cx[4] = -3*i
	Cx[5] = p4.real
	
	Cy[1] = p4.imag + p1.imag*1j
	Cy[2] = p2.imag*1j
	Cy[3] = p3.imag*1j
	x = fourier(t, Cx)
	y = fourier(t, Cy)
	if i==0:
		ttx = -np.linspace( p5.imag- p5.real, p5.imag + p5.real, 10)
		tty = p5.imag * np.ones(10)
	elif i==1:
		tty = np.linspace( p5.imag- p5.real, p5.imag + p5.real, 10)
		ttx = -p5.imag * np.ones(10)
	x = np.append(x,ttx)
	y = np.append(y,tty)
	line.set_data(y, -x) 
	return line, 
# call the animator. blit=True means only re-draw the parts that have changed. 
anim = animation.FuncAnimation(fig, animate, init_func=init, 
								frames=2, interval=100, blit=True) 
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264']) 
plt.show() 
