
def Parachute(alt,vel):
    '''Altitude in km, vel in m/s'''
    from scipy.interpolate import interp1d
    v = [438.5,487]
    h = [6,7.98]
    
    val = (alt >= h[0]) and (vel >= 310) and (vel <= v[0])
    if val:
        return val
        
    elif (vel > v[0]) and (vel < v[1]) and (alt >= interp1d(v,h)(vel)):
        v2 = [476.4,v[1]]
        h2 = [16.73,h[1]]
        if vel > v2[0]:
            return (alt<=interp1d(v2,h2)(vel))
        else:
            return True
    else:
        return False
        
     
def ParachuteTest():
    import numpy as np
    import matplotlib.pyplot as plt
    h = np.linspace(5,17)
    v = np.linspace(300,500)
    
    for alt in h:
        for vel in v:
            if Parachute(alt,vel):
                plt.plot(vel,alt,marker = 'o',color='b')
            else:
                plt.plot(vel,alt,marker = 'x',color ='r')
            
    plt.show()