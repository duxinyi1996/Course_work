import numpy as np
import matplotlib.pyplot as plt
import math

def f_r(r):
    fr = (1+ math.exp((r-R)/a))**(-1)
    return fr

def df_r(r):
    dfr = - math.exp((r-R)/a)/(a * (math.exp((r-R)/a)+1)**2)
    return dfr

def second_derivative_n(matrix,l):
    for i in range(0,size):
        matrix[i,i] += -2
        if i > 0:
            matrix[i-1,i] += 1
        if i < size-1:
            matrix[i+1,i] += 1
    if (l % 2) == 0:
        matrix[0,0] += -1
    else:
        matrix[0,0] += 1
    matrix = - p_n * matrix * (step**(-2))  
    return matrix

def spherical_n(matrix,l):
    for i in range(0,size):
        matrix[i,i] += p_n * l* (l+1)/(r[i])**2
    return matrix
    
def Vn_ls(matrix,l,s,j):
    for i in range(0,size):
        matrix[i,i] += (f_r(r[i])) * V_n + V_ls * V_n * ((j*(j+1)-l*(l+1)-s*(s+1))/2) * ((r0**2)/r[i]) * (df_r(r[i]))
    return matrix

    
def second_derivative_p(matrix,l):
    for i in range(0,size):
        matrix[i,i] += -2
        if i > 0:
            matrix[i-1,i] += 1
        if i < size-1:
            matrix[i+1,i] += 1
    if (l % 2) == 0:
        matrix[0,0] += -1
    else:
        matrix[0,0] += 1
    matrix = - p_p * matrix * (step**(-2))  
    return matrix

def spherical_p(matrix,l):
    for i in range(0,size):
        matrix[i,i] += p_p * l* (l+1)/(r[i])**2
    return matrix


def Vp_ls(matrix,l,s,j):
    for i in range(0,size):
        matrix[i,i] += (f_r(r[i])) * V_p + V_ls * V_p * ((j*(j+1)-l*(l+1)-s*(s+1))/2) * ((r0**2)/r[i]) * (df_r(r[i]))
    return matrix

def Vp_column(matrix,l,s,j):
    for i in range(0,size):
        if r[i]<R:
            matrix[i,i] += 0.1602*8.9875*Z*(3-((r[i])**2)/R**2)/(2*R)
        else:
            matrix[i,i] += 0.1602*8.9875*Z/r[i]
    return matrix
#----------------------------------------------------------------------#
#   Main
#----------------------------------------------------------------------#
 
#   Constant
Z = 82
N = 126
A = 208
r0 = 1.27
R = (A**(1/3))*r0
a = 0.67
p_n = ((197.326)**2)/(2*938.27)
p_p = ((197.326)**2)/(2*939.56)
pi = 3.1415926


V_n = (-51 +33*(N-Z)/A)
V_p = (-51 -33*(N-Z)/A)
V_ls = -0.44

step = 0.01
size =1500


name = ['s','p','d','f','g','h','i','j']
s = 1/2
l = 0
j = 1/2

#   The radius array
r = []
for i in range (0,size):
    r += [(i+1/2)*step]
    
print('Neutrons')
print('-------------------------------------------------')
for l in range (8):
    for k in range (2):
        plt.figure()
        j=l-1/2*(2*k-1)
        if j>0:
            matrix_n = np.zeros((size,size))
            matrix_n = second_derivative_n(matrix_n,l)
            matrix_n = spherical_n(matrix_n,l)
            matrix_n = Vn_ls(matrix_n,l,s,j)
            eign_n, vector_n = np.linalg.eig(matrix_n)
            eign_n = np.real(eign_n)
            x=0
            for i in range (0, len(eign_n)):
                if eign_n[i]<0:
                    x+=1
                    print(x, name[l], ', J =', j, ', E_n =', eign_n[i])
                    aa = str(x) + str(name[l]) + '  J =' + str(j) + '  E_n =' + str(eign_n[i]) +'--'
                    bb = str(name[l]) + ', J =' + str(j) + ', E_n '+'.png'
                    plt.plot(r, vector_n[:,i])
                    plt.legend(aa)
        plt.xlabel('r(fm)')
        plt.ylabel('eigenvector')
        plt.title(aa)
        plt.savefig(bb)


        
print('-------------------------------------------------')
print('Protons')
print('-------------------------------------------------')
for l in range (8):
    for k in range (2):
        j=l-1/2*(2*k-1)
        plt.figure()
        if j>0:
            matrix_p = np.zeros((size,size))
            matrix_p = second_derivative_p(matrix_p,l)
            matrix_p = spherical_p(matrix_p,l)
            matrix_p = Vp_ls(matrix_p,l,s,j)
            matrix_p = Vp_column(matrix_p,l,s,j)
            eign_p, vector_p = np.linalg.eig(matrix_p)
            eign_p = np.real(eign_p)
            y=0
            for i in range (0, len(eign_p)):
                if eign_p[i]<0:
                    y+=1
                    print(y,name[l],', J =',j,', E_p =',eign_p[i])
                    aa = str(y) + str(name[l]) + '  J =' + str(j) + '  E_p =' + str(eign_p[i])+'--'
                    bb = str(name[l]) + ', J =' + str(j) + ', E_p ' +'.png'
                    plt.plot(r, vector_p[:,i])
                    plt.legend(aa)
        plt.xlabel('r(fm)')
        plt.ylabel('eigenvector')
        plt.title(bb)
        plt.savefig(bb)

print('-------------------------------------------------')

