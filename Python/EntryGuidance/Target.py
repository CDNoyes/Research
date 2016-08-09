from Planet import Planet
from numpy import radians

class Target:

    def __init__(self, DR=780., CR=0.):
        # from functools import partial
        
        # self.coord2range = partial(planet.range, 
        # self.range2coord = partial(planet.coord, heading0=heading, lonc, latc)

        self.DR = DR
        self.CR = CR
        # self.lat,self.lon = self.range2coord(DR,CR,heading)
        
        return
    
    # def updateHeading(heading):
        
        # return
    
    # def update(DR,CR):
    
        # return