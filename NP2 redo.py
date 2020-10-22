import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *
from sympy.physics.wigner import wigner_3j
from sympy.physics.hydrogen import R_nl
from scipy.linalg import fractional_matrix_power as matpow
from scipy import special

def C(l1,l2,L):
    return (2*l2 +1)* np.float(wigner_3j(l1,L,l2,0,0,0)**2)

def second_derivative(matrix,l):
    for i in range(0,size):
        matrix[i,i] += -2/(r[i]**2)*(- (1/2))* (step**(-2))
        if i > 0:
            matrix[i-1,i] += 1/(r[i]*r[i-1])*(- (1/2))* (step**(-2))
        if i < size-1:
            matrix[i+1,i] += 1/(r[i]*r[i+1])*(- (1/2))* (step**(-2))
    '''if (l % 2) == 0:
        matrix[0,0] += -1/(r[0]**2)
    else:
        matrix[0,0] += 1/(r[0]**2)'''

    return matrix

def spherical(matrix,l):
    for i in range(0,size):
        matrix[i,i] += (1/2)*((l+1/2)**2)/(r[i])**2
    return matrix
        
def Column(matrix,zz):
    for i in range(0,size):
        matrix[i,i] += -Z[zz]/r[i]
    return matrix
    
def V_h(matrix,h_nl):
    # book's method
    V_h = np.zeros(size)
    V = np.zeros(size)
    W = np.zeros(size)
    G = np.zeros(size)
    for i in range(size):
        for nn1, ll1 in h_nl:
            G[i] += 2*(2* ll1 +1)*(h_nl[nn1, ll1][i]**2)
    for i in range(size):
        V[i] = V[i-1]+ G[i]
        W[i] = W[i-1]+ G[i]/r[i]
    for i in range(size):
        V_h[i] = step * (V[i] / r[i] + (W[size - 1] - W[i]))
    # direct method (too slow)
    '''V_h = []
    n_h = np.zeros(size)
    for nn,ll in h_nl:
        for i in range(size):
            n_h[i] += (2 * (2 * l + 1) * (np.array(h_nl[nn,ll])[i])** 2)
    for i in range(0, size):
        integral = 0
        for j in range(0, size):
            if r[i] >= r[j]:
                integral += step * n_h[j] / r[i]
            elif r[i] < r[j]:
                integral += step * n_h[j] / r[j]
        V_h += [integral]'''
    #efficient
    '''h_sum = np.zeros(size)
    for n, l in h_nl:
        h_sum += 2 * (2 * l + 1) * (h_nl[n, l])** 2
    I1 = np.asarray([np.sum(h_sum[:i+1]) for i in range(size)]) / r * step
    I2 = np.asarray([np.sum(1 / r[i:] * h_sum[i:]) for i in range(size)]) * step
    V_h = I1 + I2'''
    
    return V_h

def get_y(h_nl,nn1,ll1):
    y = np.zeros(size)
    # book's method
    for nn2, ll2 in h_nl:
        for J_2 in range(max(ll1, ll2), ll1+ll2+1):   # J_2 means J/2
            L = 2 * J_2 - ll1 - ll2
            V = np.zeros(size)
            W = np.zeros(size)
            F = np.zeros(size)
            for i in range(size):
                V[i]= V[i-1]+ h_nl[nn1,ll1][i]* C(ll1,ll2,L)* h_nl[nn2,ll2][i]*(r[i]**L)
                W[i]= W[i-1]+ h_nl[nn1,ll1][i]* C(ll1,ll2,L)* h_nl[nn2,ll2][i]/r[i]**(L+1)
            for i in range(size):
                F[i]= step*(V[i]/r[i]**(L+1)+(W[size-1]-W[i])*r[i]**L)
                y[i] += F[i] * h_nl[nn2,ll2][i]
    # direct method (too slow)
    '''for i in range(size):
        integral = 0
        for j in range(size):
            for nn2,ll2 in h_nl:
                for J_2 in range(max(ll1, ll2), ll1+ll2+1):
                    sum_l = 2 * J_2 - ll1 - ll2                
                    if r[i] >= r[j]:
                        integral += (h_nl[nn1,ll1][j]) * step * (h_nl[nn2,ll2][j]) * (h_nl[nn2,ll2][i]) * (r[j]** sum_l) / (r[i] ** (sum_l+1)) * C(ll1,ll2,sum_l)
                    elif r[i] < r[j]:
                        integral += (h_nl[nn1,ll1][j]) * step * (h_nl[nn2,ll2][j]) * (h_nl[nn2,ll2][i]) * (r[i]** sum_l) / (r[j] ** (sum_l+1)) * C(ll1,ll2,sum_l)
        y[i]= integral'''
    # efficient
    '''for nn2, ll2 in h_nl:
        for J_2 in range(max(ll1, ll2), ll1+ll2+1):   # J_2 means J/2
            L = 2 * J_2 - ll1 - ll2
            I1 = np.asarray([np.sum(r[:i+1] ** L * h_nl[nn1, ll1][:i+1] * h_nl[nn2, ll2][:i+1]) for i in range(size)]) * step
            I2 = np.asarray([np.sum(r[i:] ** (-(L + 1)) * h_nl[nn1, ll1][i:] * h_nl[nn2, ll2][i:]) for i in range(size)]) * step
            y += r ** -(L + 1) * h_nl[nn2, ll2] * C(ll1, ll2, L) * I1
            y += r ** L * h_nl[nn2, ll2] * C(ll1, ll2, L) * I2'''

    return y
def ortho_normal(h_nl,zz):


    '''group = []
    n_list = []
    l_list = []
    for n, l in h_nl:
        n_list += [n]
        l_list += [l]
        group += [Matrix(h_nl[n,l])]
    group = GramSchmidt(group,True)
    for i in range(0,len(group)):
        for j in range(size):
            h_nl[n_list[i],l_list[i]][j] = float(group[i][j,0])
    for n,l in h_nl:
        h_nl[n,l]= h_nl[n,l]/ math.sqrt(step) / np.linalg.norm(h_nl[n, l])'''
        
    n_list = nn[:zz+1]
    l_list = ll[:zz+1]
    n_forl = dict()
    ndex_forl = dict()
    l_forn = dict()
    ldex_forn = dict()
    l_set = set(l_list)
    # n,n' read
    for l in l_set:
        ndex_forl[l] = []
        n_forl[l] = []
    for i in range(zz+1):
        ndex_forl[l_list[i]] += [n_list[i]]
        n_forl[l_list[i]] += [h_nl[n_list[i],l_list[i]]]
    
    # n,n'
    for l in l_set:
        n_forl[l] = np.asarray(n_forl[l])
        #
        integral = n_forl[l] @ np.transpose(n_forl[l]) * step
        n_forl[l] = matpow(integral,-0.5) @ n_forl[l]
        for i in range(len(ndex_forl[l])): 
            h_nl[ndex_forl[l][i],l]= n_forl[l][i]

    '''    new = np.zeros(group.shape)
        for i in range(len(group)):
            new[i] = group[i]
            for j in range(0, i):
                new[i] -= np.dot(group[i], new[j]) / (np.dot(new[j], new[j])) * new[j]
        for i in range(len(ndex_forl[l])):
            h_nl[ndex_forl[l][i], l] = new[i]'''
            
    # l,l' read
    '''for n in n_list:
        l_forn[n] =[]
        ldex_forn[n] =[]
    for n,l in h_nl:
        l_forn[n] +=[h_nl[n,l]]
        ldex_forn[n] += [l]'''
    # l,l'
    '''for n in l_forn:
        group = np.array(l_forn[n])
        #
        integral = group @ np.transpose(group) *step
        group = matpow(integral,-0.5) @ group
        for i in range(len(ldex_forn[n])):
            h_nl[n, ldex_forn[n][i]] = group[i]'''
        #
    '''new = np.zeros(group.shape)
        for i in range(len(group)):
            new[i] = group[i]
            for j in range(0, i):
                new[i] -= np.dot(group[i], new[j]) / (np.dot(new[j], new[j])) * new[j]
        for i in range(len(ldex_forn[n])):
            h_nl[n, ldex_forn[n][i]] = new[i]'''   
    '''for n,l in h_nl:
        h_nl[n,l]= h_nl[n,l]/ math.sqrt(step* np.dot(h_nl[n, l],h_nl[n, l]))'''
    return h_nl
    
def normal(h_nl):
    for n,l in h_nl:
        h_nl[n,l]= h_nl[n,l]/ math.sqrt(step* np.dot(h_nl[n, l],h_nl[n, l]))
    return h_nl




    

#----------------------------------------------------------------------#
#   Main
#----------------------------------------------------------------------#
size = 1000
r1 = 25

Z = [2,4,10,12,18]
name =['He','Be','Ne','Mg','Ar']
nn = [1,2,2,3,3]
ll = [0,0,1,0,1]


for zz in range(0,len(Z)):   # 1 <-> len(Z)
    print('for Atom',name[zz])
    # build r & u grid
    r0 = 10**(-7)/Z[zz]
    step = (1/(size-1))* math.log(r1/r0)
    u = np.zeros(size)
    r = np.zeros(size)
    for i in range (size):
        u[i]= math.log(r0)+i*step # u = ln(r0) + delta*(i-1)
        r[i]= math.exp(u[i])          # r = e^u
    h_nl = dict()
    
    for zzz in range(zz+1):
        wave = np.zeros(size)
        for i in range(size):
            wave[i] = float(R_nl(nn[zzz],ll[zzz],r[i],Z[zz])*r[i]**(3/2))
        h_nl[nn[zzz],ll[zzz]] = np.array(wave)
        '''psi_nl = np.exp(- Z[zz] * r / nn[zzz]) * (2*Z[zz]*r/nn[zzz])**ll[zzz] * special.eval_genlaguerre(nn[zzz]-ll[zzz]-1, 2*ll[zzz]+1, 2*Z[zz]*r/nn[zzz])
        h = r ** (3/2) * psi_nl
        h_nl[nn[zzz], ll[zzz]] = h'''
    h_nl = normal(h_nl)
    
    e_nl = dict()
    e_old = dict()
    h_new = dict()
    for n, l in h_nl:
        e_old[n,l] = 5
        
    iterate=0
    while 1:
        iterate += 1
        for n, l in h_nl:
            '''h_nl[n, l] = h_nl[n, l] / math.sqrt(step * np.dot(h_nl[n, l],h_nl[n, l]))'''
            matrix = np.zeros((size, size))
            matrix = second_derivative(matrix, l)
            matrix = Column(matrix, zz)
            matrix = spherical(matrix, l)
            V_nl = V_h(matrix, h_nl)
            y_nl = get_y(h_nl, n, l)
            e = np.sum( h_nl[n, l]*(matrix @ h_nl[n, l]+ V_nl * h_nl[n, l] - y_nl)) * step
            e = min( e , 0 )
            for i in range(size):
                matrix[i, i] += -e
            x = np.linalg.solve(matrix + np.diag(V_nl),y_nl)
            h_new[n, l] = x
            e_nl[n, l] = e
        
        h_nl = normal(h_new)
        h_nl = ortho_normal(h_nl,zz)
        
        h_new = dict()
        
        discrepancy = 0.0
        index = 0
        for n, l in h_nl:
            if e_nl[n,l] >= 0:
                index += 1
            for i in range(size):
                discrepancy += abs(e_old[n,l] - e_nl[n, l])
        if discrepancy < 0.0001 and index == 0 and iterate > 20:
            print(e_nl)
            break
        else:
            e_old = e_nl
            e_nl= dict()

       
        '''plt.figure()
        plt.plot(r,h_nl[n,l])
        plt.show()'''





        
        
    

