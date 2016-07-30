

def createResultsDir():
    import os
    if not os.path.exists('./Results/'):
        os.mkdir('./Results/')

def loggingInit():
    from perturbUtils import inputLogInit
    inputLogInit()
    createResultsDir()

    return

def loggingFinish(states,units,data):
    import numpy as np
    from perturbUtils import writeInputLog
    writeInputLog()
    
    fmt = '%-19.4f'
    with file('./Results/trajectory.txt','w') as result:
        
        result.write('#Trajectory Data\n')
        for state in states:
            result.write("{0:20}".format(state)) 
        result.write('\n')
        for unit in units:
            result.write("{0:20}".format(unit)) 
        result.write('\n')
        np.savetxt(result, data, fmt = fmt)