#Define the equations of motion for an entry vehicle under various assumptions
from math import sin, cos
import numpy as np
#3DOF, Non-rotating Earth (i.e. Coriolis terms are excluded)
def entry_3dof(x,t):
	dh = v*sin(gamma)
	dv = -D - g*sin(gamma)
	dgamma = L/V*cos(sigma) + cos(gamma)*(v/r - g)
	return np.array([]})
	
#3DOF, Rotating Earth Model
def entry_vinhs(x,t):
	
	return

#2DOF, Longitudinal Model
def entry_2dof(x,t):

	dh = v*sin(gamma)
	ds = v*cos(gamma)
	dv = -D - g*sin(gamma)
	dgamma = L/V*cos(sigma) + cos(gamma)*(v/r - g)
	
	return np.array([dh,ds,dv,dgamma])