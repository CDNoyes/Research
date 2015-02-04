from math import exp

class Planet:
    def __init__(self, name = 'Mars'):
        self.name = name.capitalize()
        
        if self.name == 'Mercury':
            self.radius = float('nan')  # equatorial radius, m
            self.omega = float('nan')     # angular rate of planet rotation, rad/s
            self.mu = float('nan')      # gravitational parameter, m^3/s^2
            print 'Planet model not yet implemented!'
        elif self.name == 'Venus':
            self.radius = float('nan')    
            self.omega = float('nan')
            self.mu = float('nan')
            
        elif self.name == 'Earth':
            self.radius = 6378.1e3
            self.omega = 7.292115e-5
            self.mu = 3.98600e14
            
        elif self.name == 'Mars':
            self.radius = 3397e3    
            self.omega = 7.095e-5     
            self.mu = 4.2830e13
        
        elif self.name == 'Saturn':
            self.radius = float('nan')    
            self.omega = float('nan')
            self.mu = float('nan')
            
        elif self.name == 'Jupiter':
            self.radius = float('nan')    
            self.omega = float('nan')
            self.mu = float('nan')
            
        elif self.name == 'Uranus':
            self.radius = float('nan')    
            self.omega = float('nan')
            self.mu = float('nan')
            
        elif self.name == 'Neptune':
            self.radius = float('nan')    
            self.omega = float('nan')
            self.mu = float('nan')
            
        else:
            print 'Input planet name, '+ self.name +', is not valid'
        
    def atmosphere(self, h):
        if self.name == 'Mars':
        #Density computation:
            rho0 = 0.0158
            scale_height = 9354.5
            rho = rho0*exp(-h/scale_height)
        # Local speed of sound computation:
            coeff = [223.8, -0.2004e-3, -1.588e-8, 1.404e-13]
            a = 0
            i = 0
            for c in coeff:
                a = a + c*h**i
                i = i+1
            
            return rho,a
        else:
          print 'Atmosphere model not yet implemented!'
          return float('nan'), float('nan')