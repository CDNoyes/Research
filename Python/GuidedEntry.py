import numpy as np
from scipy.integrate import odeint,ode
from scipy import linalg
from scipy.interpolate import interp1d

import sys
from os import path
sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
from Utils.loggingUtils import loggingInit,loggingFinish
import Utils.perturbUtils as perturb
loggingInit()

from EntryGuidance.EntryEquations import Entry
from EntryGuidance.HighElevationPlanner import HEPBank

#Consider using an fsm

entry = Entry()
r0, theta0, phi0, v0, gamma0, psi0,s0 = (3540.0e3, np.radians(-90.07), np.radians(-43.90),
                                      5505.0, np.radians(-14.15), np.radians(4.99), 780e3)
x0 = np.array([r0+perturb.getVar('altitude'), theta0+perturb.getVar('longitude'), phi0+perturb.getVar('latitude'), v0, gamma0, psi0, s0])
time = np.linspace(0,300,1000)

hep = lambda x,t: HEPBank(t,5,75,133)
X = odeint(entry.dynamics(hep), x0, time, mxstep = 1000)

# print entry.tf
try:
    idx = np.where(np.diff(X[:,3])==0)[0][0]
except:
    idx = X.shape[0]

r,theta,phi,v,gamma,psi,s = X[0:idx,0], X[0:idx,1], X[0:idx,2], X[0:idx,3], X[0:idx,4], X[0:idx,5], (s0-X[0:idx,6])/1000
h = [1e-3*(R-entry.planet.radius) for R in r]
bank = [np.degrees(b) for b in hep(0,time[0:idx])]
# s = cumtrapz(entry.planet.radius/1000*v*np.cos(gamma)/r, time, initial=0)

#Compute other variables to log, altitude, DR, CR, bank angle, etc
#Construct a single python dictionary and save it for use with MCP codes.

loggingFinish(states = ['time','bank','altitude','radius','longitude','latitude','velocity','fpa', 'heading', 'range'],
              units =  ['s',   'deg', 'km',      'm',     'rad',      'rad',     'm/s',     'rad', 'rad',     'km'], 
              data = np.c_[time[0:idx],bank,h,X[0:idx,0:6],s])
