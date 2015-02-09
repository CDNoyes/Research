import numpy as np

class Vanderpol:
    def __init__(self, mu):
        self.mu = mu

    def __vdp(self,x,t,control_fun):
        x1,x2 = x
        dx1 = x2
        dx2 = -x1+self.mu*(1-x1**2)*x2+control_fun(x,t)
        return np.array([dx1, dx2])
   
    def dynamics(self, control_fun):
#        return lambda x,t: np.array([x[1], -x[0] + self.mu*(1-x[0]**2)*x[1] + control_fun(x,t)])
        return lambda x,t: self.__vdp(x,t,control_fun)
            
        
    def jacobian(self):
        return lambda x,t: np.array([[0, 1], [-1-2*self.mu*x[0]*x[1], self.mu*(1-x[0]**2)]])
        
