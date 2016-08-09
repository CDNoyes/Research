import numpy as np
from scipy.integrate import odeint, trapz
from scipy import linalg
from scipy.interpolate import interp1d

from Utils.redirect import stdout_redirected
from Utils.loggingUtils import loggingInit,loggingFinish
import Utils.perturbUtils as perturb
loggingInit() # Has to be called before any perturbations are pulled, i.e. before Entry equations

from EntryGuidance.EntryEquations import Entry
from EntryGuidance.HighElevationPlanner import Optimize
from EntryGuidance.Triggers import BiasDeployParachute,DeployParachute, findTriggerPoint
#Consider using an fsm

entry = Entry()
r0, theta0, phi0, v0, gamma0, psi0,s0 = (3540.0e3, np.radians(-90.07), np.radians(-43.90),
                                      5505.0, np.radians(-14.15), np.radians(4.99), 780e3)
x0 = np.array([r0+perturb.getVar('altitude'), theta0+perturb.getVar('longitude'), phi0+perturb.getVar('latitude'), v0, gamma0, psi0, s0])
time = np.linspace(0,500,1500)

n = 2
hep,ts = Optimize(x0,n)
ts = np.hstack((np.zeros(3-n),ts))

with stdout_redirected():    
    X = odeint(entry.dynamics(hep), x0, time, mxstep = 100)


idx = findTriggerPoint(X,time)

energy = entry.energy(X[0:idx,0],X[0:idx,3])
eSwitch = interp1d(time[0:idx],energy,'cubic')(ts)        
    
r,theta,phi,v,gamma,psi,s = X[0:idx,0], np.degrees(X[0:idx,1]), np.degrees(X[0:idx,2]), X[0:idx,3], np.degrees(X[0:idx,4]), np.degrees(X[0:idx,5]), (s0-X[0:idx,6])/1000
h = [entry.altitude(R,km=True) for R in r]
bank = [np.degrees(hep(xx,tt)) for xx,tt in zip(X[0:idx,:],time[0:idx])]
range = [entry.planet.range(*x0[[1,2,5]],lonc=lon,latc=lat,km=True) for lon,lat in zip(theta,phi)]
energy = energy[0:idx] #entry.energy(r,v)
L,D = entry.aeroforces(r,v)
Dbar = trapz(D,energy)
Lbar = trapz(L,energy)


loggingFinish(states = ['time','energy','bank','altitude','radius','longitude','latitude','velocity','fpa', 'heading', 'DR', 'CR','lift', 'drag'],
              units =  ['s',   '-',     'deg', 'km',      'm',     'deg',      'deg',     'm/s',     'deg', 'deg',     'km', 'km','m/s^2','m/s^2'], 
              data = np.c_[time[0:idx], energy, bank, h,   X[0:idx,0:6],                                                range,     L,      D],
              summary = ([('state','-','trajSummary'),
                         ('lbar','m/s^2',Lbar),
                         ('dbar','m/s^2',Dbar),
                         ('lodbar','-',Lbar/Dbar),
                         ('ts1','s',ts[0]),('ts2','s',ts[1]),('ts3','s',ts[2]),
                         ('es1','-',eSwitch[0]),('es2','-',eSwitch[1]),('es3','-',eSwitch[2])
                         ]))
