import numpy as np
from scipy.integrate import odeint, romb, cumtrapz
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
# Write an addEntry method for the trajectory log, pass tuples of state, units, data
# Pull initial state from DB w/ or w/o perturbing


entry = Entry()
r0, theta0, phi0, v0, gamma0, psi0,s0 = (3540.0e3, np.radians(-90.07), np.radians(-43.90),
                                      5505.0, np.radians(-14.15), np.radians(4.99), 780e3)
x0 = np.array([r0+perturb.getVar('altitude'), theta0+perturb.getVar('longitude'), phi0+perturb.getVar('latitude'), v0, gamma0, psi0, s0])
time = np.linspace(0,300,500)

hep = lambda x,t: HEPBank(t,5,75,133)
X = odeint(entry.dynamics(hep), x0, time)

print entry.tf
idx = np.where(np.diff(X[:,3])==0)
# print idx.size
# print idx[0].size
print time[idx[0][0]]
r,theta,phi,v,gamma,psi,s = X[:,0], X[:,1], X[:,2], X[:,3], X[:,4], X[:,5], (s0-X[:,6])/1000
h = [1e-3*(R-entry.planet.radius) for R in r]
bank = [np.degrees(b) for b in hep(0,time)]
# s = cumtrapz(entry.planet.radius/1000*v*np.cos(gamma)/r, time, initial=0)

#Compute other variables to log, altitude, DR, CR, bank angle, etc
#Construct a single python dictionary and save it for use with MCP codes.

loggingFinish(states = ['time','bank','altitude','radius','longitude','latitude','velocity','fpa', 'heading', 'range'],
              units =  ['s',   'deg', 'km',      'm',     'rad',      'rad',     'm/s',     'rad', 'rad',     'km'], 
              data = np.c_[time,bank,h,X,s])
