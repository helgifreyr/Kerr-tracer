from sympy import *

t, r, th, ph = symbols('t, r, th, ph')
xt, xr, xth, xph = symbols('xt, xr, xth, xph')
rho, Delta, a, M = symbols('rho Delta a M')

print(solve((((2*M*r - rho)*xt**2 + (2*M*r*(a**2 + r**2) + Delta*rho)*sin(th)**2*xph**2)*Delta + Delta*rho**2*xth**2 + rho**2*xr**2)/(Delta*rho), xt)[1])