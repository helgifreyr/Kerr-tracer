from RHS import RHS
import tarfile

def writeToFile(proper_time,solutions,p,endPoint):
    """
    Write the solutions to a file
    """
    with open('output.dat', 'w') as f:
        r1 = solutions[0][3]
        th1 = solutions[0][5]
        ph1 = solutions[0][7]
        xth1 = solutions[0][4]
        xph1 = solutions[0][6]
        f.write('# '+' '.join([str(i) for i in p])+' '+str(r1)+' '+str(endPoint)+'\n')
        # Print & save the solution.
        for time, solution in zip(proper_time, solutions):
            f.write(str(time)+' '.join([" %s" % i for i in solution])+'\n')    
    # tarName = 'data/r1='+'%4.2f'%r1 + '-th1='+'%4.3f'%th1 + '-ph1='+'%4.3f'%ph1 + '-xth1='+'%4.3f'%xth1 + '-xph1='+'%4.3f'%xph1+'.tar.gz'
    # make_tarfile(tarName, 'output.dat')

def solveGeodesic(RightHandSide,s,p,x0):
    """
    Solve the equations. We need to check at each time step whether the particle is too close to the black hole or has escaped.
    Define escaped as 10% further away from the BH than it started. This makes no sense if we start close to the BH.
    """
    solutions = [[0,0,0,0] for i in range(numpoints)]
    solutions[0] = x0
    rs,M,a = p

    # ODE solver parameters
    abserr = 1.0e-6
    relerr = 1.0e-6

    # endPoint is 0 if the ray escapes and 1 if it falls in. keep this in case I ever want to draw the shadow
    endPoint = 0
    # loop each time step, kill loop if r is too close to rs
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

def make_tarfile(output_filename, source_dir):
    """
    tar the solution and save it
    """
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

# Use ODEINT to solve the differential equations defined by the vector field
from scipy.integrate import odeint
from numpy import sqrt, zeros, pi, linspace, append, loadtxt, pi, abs, sin, cos, tan, sqrt
import os
import time
from random import random, uniform, randint
tstart = time.time()

# Parameter values. these provide the initial velocities via line 11
M = 1
a = 0.5
rs = M+sqrt(M**2-a**2)
rs = 2*M

# Initial conditions
# randomize the initial position
s1 = 0
r1 = uniform(15,25) # in multiples of rs
# th1 = pi/2.
# ph1 = pi/2. # this is where "around" the BH the ray starts, 0 for east, pi/2 for north, etc
th1 = uniform(0,pi)
ph1 = uniform(0,2*pi)

xr1 = -0.1 # this is initial r velocity if you'd like
# here below are tests for the angles
# xth1 = 0.001 # up above the BH and then shoots below items
# xth1 = -0.001 # down below the BH and then shoots up above it
# xph1 = -0.001 # a bit up
# xph1 = 0 # straight in
# xph1 = 0.001 # a bit down

# randomize the initial angles but always shoot towards the black hole (ish)
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
