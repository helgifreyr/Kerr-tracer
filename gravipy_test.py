from gravipy import *

t, r, th, ph, M = symbols('t, r, th, ph, M')
x = Coordinates('s', [t, r, th, ph])
a, M = symbols('a, M')
#a=0
rho = r**2 + a**2 * cos(th)**2
Delta = r**2 - 2*M*r + a**2

rplus = M + sqrt(M**2 - a**2)

tt = 2 * M * r / rho - 1
rr = rho / Delta
thth = rho
phph = (Delta + (2*M*r*(r**2+a**2))/rho) * sin(th)**2
tph = -4*a*M*r*sin(th)**2/rho

Metric = zeros(4,4)
Metric[0,0] = tt
Metric[1,1] = rr
Metric[2,2] = thth
Metric[3,3] = phph
Metric[3,1] = tph
Metric[1,3] = -tph
g = MetricTensor('g', x, Metric) 



s = Symbol('s')
w = Geodesic('w', g, s)
#print(w(All).transpose())
for i in range(1,len(x.c)+1):
    print('D'+str(x(-i))+'/ds = ' + str(w(i)))
    print("----------------------------------")