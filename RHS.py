from scipy import sin, cos, tan, sqrt

def RHS(w, s, p):
    """
    Defines the differential equations for the system 

    Arguments:
        x :  vector of the state variables:
                  x = [t, xr=r', r, xth=th', th, xph=ph', ph]
        s :  affine parameter
        p :  vector of the parameters:
                  p = [rs,M,a]
    """
    xt, t, xr, r, xth, th, xph, ph = w
    rs, M, a = p
    rho = r**2 + a**2 * cos(th)**2
    Delta = r**2 - 2*M*r + a**2
    xt = sqrt(-(Delta*rho*xph**2*sin(th)**2 + 2*M*a**2*r*xph**2*sin(th)**2 + 2*M*r**3*xph**2*sin(th)**2 + rho**2*xth**2 + rho**2*xr**2/Delta)/(2*M*r - rho))
    
    geodesic = [0,0,0,0]
    geodesic[0] = 2.0*M*(a**2*r*sin(2*th)*xth + (a**2*cos(th)**2 - r**2)*xr)*xt/((a**2*cos(th)**2 + r**2)*(-2*M*r + a**2*cos(th)**2 + r**2))

    geodesic[1] = -(-1.0*r*(a**2*cos(th)**2 + r**2)**2*(2*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(-2*M*r + a**2 + r**2))*(-2*M*r + a**2 + r**2)**2*xth**2 + (16.0*M**2*a**2*r*(a**2*cos(th)**2 - r**2)*(-2*M*r + a**2 + r**2)**2*sin(th)**2 + (a**2*cos(th)**2 + r**2)**2*(r*(-2*M*r + a**2 + r**2) + (M - r)*(a**2*cos(th)**2 + r**2))*(2*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(-2*M*r + a**2 + r**2)))*xr**2 + (-2*M*r + a**2 + r**2)*(-8.0*M*a**3*r*(-2*M*r + a**2 + r**2)*(2*M*a**2*r*cos(th)**2 + 2*M*r**3 - a**4*cos(th)**2 - a**2*r**2*cos(th)**2 - a**2*r**2 - r**4)*sin(th)**3*cos(th)*xph*xth + 4.0*M*a*r*(-2*M*r + a**2 + r**2)*(2*M*r**2*(a**2 + r**2) - M*(a**2 + 3*r**2)*(a**2*cos(th)**2 + r**2) + (M - r)*(a**2*cos(th)**2 + r**2)**2)*sin(th)**2*xph*xr + 4.0*M*a*(r*(2*M*r**2*(a**2 + r**2) - M*(a**2 + 3*r**2)*(a**2*cos(th)**2 + r**2) + (M - r)*(a**2*cos(th)**2 + r**2)**2) + (a**2*cos(th)**2 - r**2)*(2*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(-2*M*r + a**2 + r**2)))*(-2*M*r + a**2 + r**2)*sin(th)**2*xph*xr - 1.0*M*(a**2*cos(th)**2 - r**2)*(2*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(-2*M*r + a**2 + r**2))*(-2*M*r + a**2 + r**2)*xt**2 + 2*a**2*(16.0*M**2*r**2*(a**2 + r**2)*(-2*M*r + a**2 + r**2) + (a**2*cos(th)**2 + r**2)**2*(-2.0*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(2.0*M*r - 1.0*a**2 - 1.0*r**2)))*sin(th)*cos(th)*xr*xth + (2*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(-2*M*r + a**2 + r**2))*(-2.0*M*r + 1.0*a**2 + 1.0*r**2)*(2*M*r**2*(a**2 + r**2) - M*(a**2 + 3*r**2)*(a**2*cos(th)**2 + r**2) + (M - r)*(a**2*cos(th)**2 + r**2)**2)*sin(th)**2*xph**2))/((a**2*cos(th)**2 + r**2)*(16*M**2*a**2*r**2*(-2*M*r + a**2 + r**2)*sin(th)**2 + (a**2*cos(th)**2 + r**2)**2*(2*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(-2*M*r + a**2 + r**2)))*(-2*M*r + a**2 + r**2))

    geodesic[2] = -(-1.0*M*a**2*r*(-2*M*r + a**2 + r**2)*sin(2*th)*xt**2 - 0.5*a**2*(a**2*cos(th)**2 + r**2)**2*(-2*M*r + a**2 + r**2)*sin(2*th)*xth**2 + 0.5*a**2*(a**2*cos(th)**2 + r**2)**2*sin(2*th)*xr**2 + 2.0*r*(a**2*cos(th)**2 + r**2)**2*(-2*M*r + a**2 + r**2)*xr*xth - (2.0*M*a**2*r*(a**2 + r**2)*sin(th)**2 + 1.0*(a**2*cos(th)**2 + r**2)*(2*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(-2*M*r + a**2 + r**2)))*(-2*M*r + a**2 + r**2)*sin(th)*cos(th)*xph**2)/((a**2*cos(th)**2 + r**2)**3*(-2*M*r + a**2 + r**2))

    geodesic[3] = -(-8.0*M*a*r*(a**2*cos(th)**2 + r**2)**2*(-2*M*r + a**2 + r**2)*(a**2*sin(th)**2 + a**2 + r**2)*xr*xth + 4.0*M*a*(a**2*cos(th)**2 + r**2)**2*(r*(r*(-2*M*r + a**2 + r**2) + (M - r)*(a**2*cos(th)**2 + r**2)) + (a**2*cos(th)**2 - r**2)*(2*M*r - a**2 - r**2))*tan(th)*xr**2 - (a**2*cos(th)**2 + r**2)**2*(4.0*M*a*r**2*(-2*M*r + a**2 + r**2)*xth**2 + (2.0*M*r**2*(a**2 + r**2) - M*(1.0*a**2 + 3.0*r**2)*(a**2*cos(th)**2 + r**2) + 1.0*(M - r)*(a**2*cos(th)**2 + r**2)**2)*xph*xr)*(-2*M*r + a**2 + r**2)*tan(th) + 2*(16.0*M**2*a**2*r**2*(a**2 + r**2)*(-2*M*r + a**2 + r**2)*sin(th)**2 + (a**2*cos(th)**2 + r**2)**2*(2.0*M*a**2*r*(a**2 + r**2)*sin(th)**2 + 1.0*(a**2*cos(th)**2 + r**2)*(2*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(-2*M*r + a**2 + r**2))))*(-2*M*r + a**2 + r**2)*xph*xth + (-2*M*r + a**2 + r**2)*(-4.0*M**2*a*r*(a**2*cos(th)**2 - r**2)*(-2*M*r + a**2 + r**2)*xt**2 + 4.0*M*a*r*(-2*M*r + a**2 + r**2)*(2*M*r**2*(a**2 + r**2) - M*(a**2 + 3*r**2)*(a**2*cos(th)**2 + r**2) + (M - r)*(a**2*cos(th)**2 + r**2)**2)*sin(th)**2*xph**2 + (16.0*M**2*a**2*r*(a**2*cos(th)**2 - r**2)*(-2*M*r + a**2 + r**2)*sin(th)**2 + (a**2*cos(th)**2 + r**2)**2*(-2.0*M*r**2*(a**2 + r**2) + M*(1.0*a**2 + 3.0*r**2)*(a**2*cos(th)**2 + r**2) + 1.0*(-M + r)*(a**2*cos(th)**2 + r**2)**2))*xph*xr)*tan(th))/((a**2*cos(th)**2 + r**2)*(16*M**2*a**2*r**2*(-2*M*r + a**2 + r**2)*sin(th)**2 + (a**2*cos(th)**2 + r**2)**2*(2*M*r*(a**2 + r**2) + (a**2*cos(th)**2 + r**2)*(-2*M*r + a**2 + r**2)))*(-2*M*r + a**2 + r**2)*tan(th))

    # Create f = w' = (xt',t',xr',r',xth',th',xph',ph'):
    f = [geodesic[0],
         xt,
         geodesic[1],
		 xr,
		 geodesic[2],
		 xth,
		 geodesic[3],
		 xph]
    return f