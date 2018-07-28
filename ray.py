from RHS import RHS

def writeToFile(proper_time,solutions,p,endPoint):
    with open('output.dat', 'w') as f:
        r1 = solutions[0][3]
        f.write('# '+' '.join([str(i) for i in p])+' '+str(r1)+' '+str(endPoint)+'\n')
        # Print & save the solution.
        for time, solution in zip(proper_time, solutions):
            f.write(str(time)+' '.join([" %s" % i for i in solution])+'\n')    

def solveGeodesic(RightHandSide,s,p,x0):
    
    solutions = [[0,0,0,0] for i in range(numpoints)]
    solutions[0] = x0
    rs,M,a = p

    # ODE solver parameters
    abserr = 1.0e-6
    relerr = 1.0e-6

    # loop each time step, kill loop if r is too close to rs
    endPoint = 0
    for idx in range(1,int(numpoints)):
        if x0[3]>1.01*rs and x0[3]<1.1*r1:
            solutions[idx] = odeint(RightHandSide, x0, [s[idx-1],s[idx]], args=(p,), atol=abserr, rtol=relerr)[1]
            x0 = solutions[idx]
        else:
            if x0[3]<=1.05*rs:
                endPoint = 1
            if x0[3]>=5*r1:
                endPoint = 0
            writeToFile(s[:idx],solutions[:idx],p,endPoint)
            break
        #time.sleep(0.1)
    return s[:idx], solutions[:idx]

def run(RHS, s, p, x0):
    s, xsol=solveGeodesic(RHS,s,p,x0)
    return s, xsol

# Use ODEINT to solve the differential equations defined by the vector field
from scipy.integrate import odeint
from numpy import sqrt, zeros, pi, linspace, append, loadtxt, pi, abs, sin, cos, tan, sqrt
import sys
import time
from random import random, uniform, randint
tstart = time.time()

# Parameter values. these provide the initial velocities via line 11
M = 1
a = 0.5
rs = M+sqrt(M**2-a**2)
rs = 2*M

# Initial conditions
s1 = 0
r1 = 25 # in multiples of rs
th1 = pi/2.
ph1 = pi/2. # this is where "around" the BH the ray starts, 0 for east, pi/2 for north, etc
th1 = uniform(0,pi)
ph1 = uniform(0,2*pi)

xr1 = -0.1 # this is initial r velocity if you'd like
xth1 = 0.001 # up above the BH and then shoots below items
xth1 = -0.001 # down below the BH and then shoots up above it
xph1 = -0.001 # a bit up
xph1 = 0 # straight in
xph1 = 0.001 # a bit down
xt1 = 0

xth1 = randint(-1,1)*uniform(0.0005,0.001)
xph1 = randint(-1,1)*uniform(0.0005,0.001)


# affine parameter values
sf = 10**4
numpoints = sf*20+1
s = linspace(s1,sf,numpoints)
# Pack up the parameters and initial conditions:
p = [rs,M,a]
x0 = [xt1,s1,xr1,r1,xth1,th1,xph1,ph1]

# loop over b's
#for i in linspace(-3,1,31):
#    run(RHS, s, [rs, a, i], x0)
run(RHS, s, p, x0)
tend = time.time()
print(tend-tstart)
