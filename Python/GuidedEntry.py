import numpy as np
from scipy.integrate import odeint
from scipy import linalg
from scipy.interpolate import interp1d

import sys
from os import path
sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
from Utils.perturbUtils import writeInputLog, inputLogInit
inputLogInit()

from EntryGuidance.EntryEquations import Entry

#Consider using an fsm



entry = Entry()
r0, theta0, phi0, v0, gamma0, psi0 = (3540.0e3, np.radians(-90.07), np.radians(-43.90),
                                      5505.0, np.radians(-14.15), np.radians(4.99))
x0 = np.array([r0, theta0, phi0, v0, gamma0, psi0])
time = np.linspace(0,180,500)


flat =  lambda x,t: 0

X = odeint(entry.dynamics(flat), x0, time)
r,theta,phi,v,gamma,psi = X[:,0], X[:,1], X[:,2], X[:,3], X[:,4], X[:,5]
h = [1e-3*(R-entry.planet.radius) for R in r]
#Compute other variables to log, altitude, DR, CR, bank angle, etc
#Write time histories to files for viewing via the guis
#Construct a single python dictionary and save it for use with MCP codes.

writeInputLog()
fmt = '%-19.4f'
with file('./Results/trajectory.txt','w') as result:
    
    result.write('#Trajectory Data\n')
    for state in ['time','altitude','radius','longitude','latitude','velocity','fpa','heading']:
        result.write("{0:20}".format(state)) 
    result.write('\n')
    for unit in ['s','km','m','rad','rad','m/s','rad','rad']:
        result.write("{0:20}".format(unit)) 
    result.write('\n')
    # np.savetxt(result, np.concatenate((time.T,X),axis=1), fmt = fmt)
    # np.savetxt(result, np.hstack((time.T,X)), fmt = fmt)
    np.savetxt(result, np.c_[time,h,X], fmt = fmt)