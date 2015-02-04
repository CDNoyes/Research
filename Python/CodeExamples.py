#Import the items we need
import numpy as np
from scipy.integrate import odeint
#import matplotlib as mpl
import matplotlib.pyplot as plt

#np has array, ones, zeros, eye, diag, linspace


initial_state = np.array([3.0, 0.0]) #Define the system's initial state
time = np.linspace(0,5)

def vanderpol(x, t):
    mu = 0.25
    return [x[1], -x[0]+mu*(1-x[0]**2)*x[1]]
	
def vanderpol_gradient(x,t):
    mu = 0.25
    return [[0, 1], [-1-2*mu*x[0]*x[1], mu*(1-x[0]**2)]]
	
y = odeint(vanderpol,initial_state, time)
y2 = odeint(vanderpol, initial_state, time, Dfun = vanderpol_gradient)
x1 = y2[:,0]
x2 = y2[:,1]
#Identical to Matlab other than the namespace :
plt.plot(y[:,0], y[:,1], label = 'without gradient')
plt.hold(True)
plt.plot(x1, x2, label = 'with gradient')
plt.axis([-3, 3, -3, 3]) 
plt.xlabel('x_1')
plt.ylabel('x_2')
plt.title('Uncontrolled Vanderpol Oscillator')
plt.legend()
#plt.legend(loc = 'best', ncol = 2) #Can optionally specify location and number of columns in the legend
#plt.legend('Without Grad','With Grad') #Can be used just like Matlab