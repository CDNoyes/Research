from __future__ import division
import numpy as np
import scipy as sp

import matplotlib.pyplot as plt

#Define methods for solution to control problems via differential dynamic programming or the closely related iterative LQR/LQG methods


class Problem:
    def __init__(self, dynamics, constraints, initialState, jacobian = None):
        self.dynamics = dynamics
        self.constraints = constraints
        self.initialState = initialState
        if jacobian is not None:
            self.jacobian = jacobian    
        else:
            self.jacobian = numericalDiff(dynamics).complex
        # self.bounds = bounds

def eLQR(Problem):
    return None
    
    #%%
class numericalDiff:
    def __init__(self, dynamics):
        self.dynamics = dynamics #The dynamics should take only one input, the current state and/or the current control at which to differentiate. See example below.
        
    def complex(self, State):
    #Check inputs
        if isinstance(State, (list,tuple)):
            State = np.array(State)
        if isinstance(State, (float,int)):
            State = np.array([State])
        stepSize = 1e-6
        h = np.array([1j])
        h.imag = stepSize     
        offset = np.eye(State.size)*h
        return np.squeeze(np.array([np.imag((self.dynamics(State+offset[:,iter])))/stepSize for iter in range(State.size)])).T
        
    def central(self, State):
        #Check inputs
        if isinstance(State, (list,tuple)):
            State = np.array(State)
        if isinstance(State, (float,int)):
            State = np.array([State])
        stepSize = 1e-6
        offset = np.eye(State.size)*stepSize
        return np.squeeze(np.array([(self.dynamics(State+offset[:,iter]) - self.dynamics(State-offset[:,iter]))/stepSize/2 for iter in range(State.size)])).T
    

def numDiffExample():
    f = lambda x,u: np.array([x[1],u[0]])
    n = 2
    m = 1 
    x0 = np.array([1.0, 0.0])
    u0 = np.array([-1.0])
    # dfdx = numericalDiff(lambda x: f(x,u0)).central(x0)
    dfdxc = numericalDiff(lambda x: f(x,u0)).complex(x0)
     
    dfdu = numericalDiff(lambda u: f(x0,u)).complex(u0)
    df = numericalDiff(lambda X: f(X[0:n],X[n:n+m])).complex(np.concatenate((x0,u0),axis=1)) #Both at the same time
    # print dfdx
    print dfdxc
    print dfdu 
    print df 
    return None
    
numDiffExample()