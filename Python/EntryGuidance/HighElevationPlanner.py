import numpy as np
from EntryEquations import Entry
from Triggers import BiasDeployParachute, findTriggerPoint
from Target import Target
from scipy.integrate import odeint
from Utils.redirect import stdout_redirected

def HEPBankReduced(T,t1,t2,minBank = np.radians(15.), maxBank = np.radians(85.)):
    maxRate = np.radians(20.)
    ti1 = t1 + 2*(maxBank)/maxRate
    ti2 = t2 + (maxBank+minBank)/maxRate
    isScalar = False
    bank = []
    if isinstance(T,(float,int,np.float32,np.float64)):
        T = [T]
        isScalar = True

    for t in T:

        if t <= t1:
            bank.append(maxBank)
        elif t >= t1 and t <= ti1:
            bank.append(maxBank-maxRate*(t-t1))
        elif t >= ti1 and t <= t2:
            bank.append(-maxBank)
        elif t >= t2 and t <= ti2:
            bank.append(-maxBank + maxRate*(t-t2))
        else:
            bank.append(minBank)
    if isScalar:
        try:
            return bank[0]
        except:
            print "Bank angle comp failed"
            return -1.
    else:
        return bank

def HEPBank(T,t1,t2,t3,minBank = np.radians(15.), maxBank = np.radians(85.)):

    maxRate = np.radians(20.)
    ti1 = t1 + (maxBank+minBank)/maxRate
    ti2 = t2 + 2*(maxBank)/maxRate
    ti3 = t3 + (maxBank+minBank)/maxRate
    isScalar = False
    bank = []
    if isinstance(T,(float,int,np.float32,np.float64)):
        T = [T]
        isScalar = True

    for t in T:
        if t < t1 and t >= 0:
            bank.append(-minBank) 
        elif t >= t1 and t <= ti1:
            bank.append(-minBank+maxRate*(t-t1))
        elif t >= ti1 and t <= t2:
            bank.append(maxBank)
        elif t >= t2 and t <= ti2:
            bank.append(maxBank-maxRate*(t-t2))
        elif t >= ti2 and t <= t3:
            bank.append(-maxBank)
        elif t >= t3 and t <= ti3:
            bank.append(-maxBank + maxRate*(t-t3))
        else:
            bank.append(minBank)
    if isScalar:
        try:
            return bank[0]
        except:
            print "Bank angle comp failed"
            return -1.
    else:
        return bank

def checkFeasibility(T,sign=-1):
    
    sig = sign*np.diff(T) #np.array([T[0]-T[1],T[1]-T[2]])
    cost = (sum(sig+np.abs(sig)) + sum(abs(T)-T))*1e5
    if cost:
        return max(cost,1e7)
    else:
        return 0
        
def HEPCost(T, x0, entry = Entry(Trigger = BiasDeployParachute), target = Target(), bank = HEPBank, getIV = (lambda x,t: t), check=checkFeasibility):
    
    J = check(T)
    if J > 300:
        return J
        
    hep = lambda x,t: bank(getIV(x,t),*T)

        
    time = np.linspace(0,450,1000)

    try:
        with stdout_redirected():
            X = odeint(entry.dynamics(hep), x0, time)
    except TerminateSimulation:
        pass

    idx = findTriggerPoint(X,time)-1
    h     = entry.altitude(X[idx,0], km=True)
    dr,cr = entry.planet.range(*x0[[1,2,5]], lonc = X[idx,1], latc = X[idx,2],km=True)

    Wh = 0.8*1e-6
    J = -h*Wh + ((target.DR-dr)**2 + (target.CR-cr)**2)**0.5

    return J

    
def Optimize(x0,n=3,iv='time'):
    
    from scipy.optimize import minimize, minimize_scalar
       
    if (n == 3):
        guess = [50,100,133]
        bankFun = HEPBank
    elif (n==2):
        guess = [70,100]
        bankFun = HEPBankReduced
    else:   
        print "Optimize: improper input for n."
        return None
    
    if iv.lower() == 'time':
        getIV = lambda x,t: t
        check = checkFeasibility
        
    elif 'vel' in iv.lower():
        if n == 2:
            guess = [-4700,-2100]
        elif n == 3:
            guess = [-5100,-4500,-2500]
        getIV = lambda x,t: -x[3]
        check = lambda v: checkFeasibility(-v,sign=1)
    
    entry = Entry(Trigger = BiasDeployParachute)
    sol = minimize(HEPCost,np.array(guess),args = (x0, entry, Target(), bankFun, getIV, check), method='Nelder-Mead',tol=1e-5,options={'disp':True})
    print "The {0} switching {1}s are {2}".format(n,iv,sol.x)
    return (lambda x,t: bankFun(getIV(x,t),*sol.x)), sol.x
    
