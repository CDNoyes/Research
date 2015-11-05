#Import the items we need
import numpy as np
from scipy.integrate import odeint
from scipy import linalg
from scipy.interpolate import interp1d
from Vanderpol import Vanderpol
#from EntryEquations import Entry
import matplotlib.pyplot as plt

#np has array, ones, zeros, eye, diag, linspace as well as trig functions
#linalg has eig, solve, transpose

#%% Integrate the vanderpol oscillator system for 5 seconds:
initial_state = np.array([3.0, 0.0]) #Define the system's initial state
time = np.linspace(0,5)

vdp = Vanderpol(0.25)
def no_control(x,t):
    return 0
uncontrolled = lambda x,t: 0    # Alternative anonymous definition
u = 1.5*np.sin(time)
control_interp = interp1d(np.linspace(0,5.1), u, kind = 'cubic') # Returns a function object
control = lambda x,T: control_interp(T)
	
y = odeint(vdp.dynamics(control),initial_state, time)
y2 = odeint(vdp.dynamics(uncontrolled), initial_state, time, Dfun = vdp.jacobian())
#x1 = y2[:,0]
#x2 = y2[:,1]
x1,x2 = y2[:,0],y2[:,1]
# Plotting is identical to Matlab other than the namespace :
plt.plot(y[:,0], y[:,1], label = 'without gradient')
plt.hold(True)
plt.plot(x1, x2, label = 'with gradient')
plt.axis([-3, 3, -3, 3]) 
plt.xlabel('x_1')
plt.ylabel('x_2')
plt.title('Uncontrolled Vanderpol Oscillator')
plt.legend()
#plt.legend(loc = 'best', ncol = 2) # Can optionally specify location and number of columns in the legend
#plt.legend('Without Grad','With Grad') # Can be used just like Matlab

#%% Entry Equations!

entry = Entry()
r0, theta0, phi0, v0, gamma0, psi0 = (3540.0e3, np.radians(-90.07), np.radians(-43.90),
                                      5505.0, np.radians(-14.15), np.radians(4.99))
x0 = np.array([r0, theta0, phi0, v0, gamma0, psi0])
time = np.linspace(0,180,500)


flat =  lambda x,t: 0
turn = lambda x,t: np.radians(87)


X = odeint(entry.dynamics(turn), x0, time)
r,theta,phi,v,gamma,psi = X[:,0], X[:,1], X[:,2], X[:,3], X[:,4], X[:,5]
plt.figure()
plt.plot(v, (r-entry.planet.radius)/1000, label = '\sigma = 0')
plt.xlabel('Velocity [m/s]')
plt.ylabel('Altitude [km]')