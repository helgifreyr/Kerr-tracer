# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 20:08:23 2018

@author: helgi
"""

from sympy import *
import time
import multiprocessing
from itertools import product
tstart = time.time()

t, r, theta, phi = symbols('t, r, theta, phi')
rho, Delta = symbols('rho Delta')
coords = [t,r,theta,phi]
s = Symbol('s')
n = len(coords)
a, M = symbols('a, M')
#a=0.1
#rho = r**2 + a**2 * cos(theta)**2
#Delta = r**2 - 2*M*r + a**2


rplus = M + sqrt(M**2 - a**2)

tt = 2 * M * r / rho(r,theta) - 1
rr = rho(r,theta) / Delta(r)
thth = rho(r,theta)
phph = (Delta(r) + (2*M*r*(r**2+a**2))/rho(r,theta)) * sin(theta)**2
tph = -4*a*M*r*sin(theta)**2/rho(r,theta)

g = zeros(4,4)
g[0,0] = tt
g[1,1] = rr
g[2,2] = thth
g[3,3] = phph
g[3,1] = tph
g[1,3] = -tph
ig = simplify(g.inv())

#print(diff(g[0,0],r))

def udotu():
    return simplify(sum([sum([g[k,j]*diff(coords[k](s),s)*diff(coords[j](s),s) for k in range(n)]) for j in range(n)]))

def Chr(i,j,k):
    return simplify(0.5*sum([ig[i,s]*(diff(g[s,j],coords[k])+diff(g[s,k],coords[j])-diff(g[j,k],coords[s])) for s in range(n)]))
    
def Geodesic(i):
    return -simplify(sum([sum([Chr(i,j,k)*diff(coords[j](s),s)*diff(coords[k](s),s) for k in range(n)]) for j in range(n)]))
 
#for i in range(n):
    #print('Geodesic('+str(i)+')')
    #print(Geodesic(i))
print(solve(1+udotu(),Derivative(t(s), s)))
print(solve(r**2-1,r))


 
 # if __name__ == '__main__':
    # indices = range(n)
    # cores = multiprocessing.cpu_count()-2
    # print(indices)
    # with multiprocessing.Pool(processes=cores) as pool:
        # results = pool.starmap(Geodesic, product(indices, repeat=1))
    # print(results)
    # # num_processes = 6
    # # p = multiprocessing.Pool(num_processes)
    # # mp_solutions = p.map(Geodesic, indices)
    # tend = time.time()
    # tmp = tend - tstart
    # print(" %8.3f seconds" % tmp)