import numpy as np
from scipy.integrate import odeint, trapz
from scipy import linalg
from scipy.interpolate import interp1d
from scipy.io import savemat

from Utils.redirect import stdout_redirected
from Utils.loggingUtils import loggingInit,loggingFinish
loggingInit() # Has to be called before any perturbations are pulled, i.e. before Entry equations

from EntryGuidance.EntryEquations import Entry
from EntryGuidance.Planet import Planet
from EntryGuidance.EntryVehicle import EntryVehicle
from EntryGuidance.Triggers import BiasDeployParachute, DeployParachute, findTriggerPoint
from EntryGuidance.HighElevationPlanner import Optimize
#Consider using an fsm

#for sample in pdf.sample(nSamples):
def Simulate(sample):
    CD,CL,rho0,sh = sample
    entry = Entry(PlanetModel = Planet(rho0=rho0,scaleHeight=sh), VehicleModel = EntryVehicle(CD=CD,CL=CL))
    # entryNom = Entry()
    r0, theta0, phi0, v0, gamma0, psi0,s0 = (3540.0e3, np.radians(-90.07), np.radians(-43.90),
                                             5505.0,   np.radians(-14.15), np.radians(4.99),   780e3)
    x0 = np.array([r0, theta0, phi0, v0, gamma0, psi0, s0])
    # x0Nom = np.array([r0, theta0, phi0, v0, gamma0, psi0, s0])
    time = np.linspace(0,500,1500)

    n = 2
    # hep,ts = Optimize(x0,n)
    hep = lambda x,t: np.pi
    ts = np.zeros(3)
    # ts = np.hstack((np.zeros(3-n),ts))

    with stdout_redirected():    
        X = odeint(entry.dynamics(hep), x0, time, mxstep = 100)

    idx = findTriggerPoint(X,time)
    tInterp = np.linspace(0,time[idx-1],500)
    xInterp = interp1d(time[0:idx],X[0:idx,:],'cubic',axis=0)(tInterp)
    return xInterp




# energy = entry.energy(X[0:idx,0],X[0:idx,3])
# eSwitch = interp1d(time[0:idx],energy,'cubic')(ts)        
    
# r,theta,phi,v,gamma,psi,s = X[0:idx,0], np.degrees(X[0:idx,1]), np.degrees(X[0:idx,2]), X[0:idx,3], np.degrees(X[0:idx,4]), np.degrees(X[0:idx,5]), (s0-X[0:idx,6])/1000
# h = [entry.altitude(R,km=True) for R in r]
# bank = [np.degrees(hep(xx,tt)) for xx,tt in zip(X[0:idx,:],time[0:idx])]
# range = [entry.planet.range(*x0[[1,2,5]],lonc=np.radians(lon),latc=np.radians(lat),km=True) for lon,lat in zip(theta,phi)]
# energy = energy[0:idx] #entry.energy(r,v)
# L,D = entry.aeroforces(r,v)
# Dbar = trapz(D,energy)
# Lbar = trapz(L,energy)


# loggingFinish(states = ['time','energy','bank','altitude','radius','longitude','latitude','velocity','fpa', 'heading', 'DR', 'CR','lift', 'drag'],
              # units =  ['s',   '-',     'deg', 'km',      'm',     'deg',      'deg',     'm/s',     'deg', 'deg',     'km', 'km','m/s^2','m/s^2'], 
              # data = np.c_[time[0:idx], energy, bank, h,   r,      theta,       phi,      v,         gamma, psi,       range,     L,      D],
              # summary = ([('state','-','trajSummary'),
                         # ('lbar','m/s^2',Lbar),
                         # ('dbar','m/s^2',Dbar),
                         # ('lodbar','-',Lbar/Dbar),
                         # ('ts1','s',ts[0]),('ts2','s',ts[1]),('ts3','s',ts[2]),
                         # ('es1','-',eSwitch[0]),('es2','-',eSwitch[1]),('es3','-',eSwitch[2])
                         # ]))

if __name__ == '__main__':
    import multiprocessing as mp
    import chaospy as cp


    CD          = cp.Uniform(-0.10, 0.10)   # CD
    CL          = cp.Uniform(-0.10, 0.10)   # CL
    rho0        = cp.Normal(0, 0.0333)      # rho0
    scaleHeight = cp.Uniform(-0.05,0.05)    # scaleheight
    pdf = cp.J(CD,CL,rho0,scaleHeight)
    mp.freeze_support()
    samples = pdf.sample(48)    
    # pool = mp.Pool(mp.cpu_count()/2.)
    pool = mp.Pool(6)
    stateTensor = pool.map(Simulate,samples.T)

    savemat('./data/testingMPandMAT',{'states':stateTensor, 'samples':samples})