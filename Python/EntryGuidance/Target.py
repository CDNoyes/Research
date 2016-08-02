from Planet import Planet
from numpy import radians

class Target:

    def __init__(self, planet = Planet('Mars'), heading = radians(4.99), DR=780., CR=0., Lon=None, Lat=None):
    
    self.coord2range = planet.range
    self.range2coord = planet.coord
    
    if Lat is not None and Lon is not None:
        self.lat = lat
        self.lon = lon
        self.DR,self.CR = self.coord2range(Lat,Lon)
    else:
        self.DR = DR
        self.CR = CR
        self.lat,self.lon = self.range2coord(DR,CR,heading)
    
    return self
    
    def updateHeading(heading):
        
        return
    
    def update(DR,CR):
    
        return