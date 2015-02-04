#Define the equations of motion for an entry vehicle under various assumptions
from math import sin, cos
import numpy as np
#3DOF, Non-rotating Earth (i.e. Coriolis terms are excluded)
def entry_3dof(x, t, planet, vehicle):
    r,theta,phi,v,gamma,psi = x
    
    g = planet.mu/r**2
    h = r - planet.radius
    rho,a = planet.atmosphere(h)
    M = v/a
    cD,cL = vehicle.aerodynamic_coefficients(M)
    f = 0.5*rho*vehicle.area*v**2/vehicle.mass
    L = f*cL
    D = f*cD
    
    dh = v*sin(gamma)
    dtheta = v*cos(gamma)*cos(psi)/r/cos(phi)
    dphi = v*cos(gamma)*sin(psi)/r
    dv = -D - g*sin(gamma)
    dgamma = L/V*cos(sigma) + cos(gamma)*(v/r - g)
    dpsi = -L*sin(sigma)/v/cos(gamma) - v*cos(gamma)*cos(psi)*tan(phi)/r
    return np.array([dh, dtheta, dphi, dv, dgamma, dpsi])
	
#3DOF, Rotating Earth Model
def entry_vinhs(x,t):

    return

#2DOF, Longitudinal Model
def entry_2dof(x, t, planet, vehicle):
    r,s,v,gamma = x
    
    g = planet.mu/r**2
    h = r - planet.radius
    rho,a = planet.atmosphere(h)
    M = v/a
    cD,cL = vehicle.aerodynamic_coefficients(M)
    f = 0.5*rho*vehicle.area*v**2/vehicle.mass
    L = f*cL
    D = f*cD
    
    dh = v*sin(gamma)
    ds = v*cos(gamma)
    dv = -D - g*sin(gamma)
    dgamma = L/V*cos(sigma) + cos(gamma)*(v/r - g)

return np.array([dh,ds,dv,dgamma])