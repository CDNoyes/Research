import imuQ
from numpy.linalg import norm

db = imuQ.Project.loadDB()
unitSystem = 'SI'

global inputLog
       
def getVar(nominal, name, isVariable = True):
    import logging
    # This is a safe method for grabbing scalar values or class structures from the XML database. 
    # Set isVariable to False when pulling a class instead of a single variable.
    # If isVariable is False, in general you will need to manually call the addEntry method to add the class perturbation(s) to the input.dat file.
    modelName = 'MarsEDL'
    try:
        if isVariable:
            value = db.get('multiphysicsModel['+modelName+']::variable['+name+']').float(unitSystem=unitSystem)
            unit = db.get('multiphysicsModel['+modelName+']::variable['+name+']').units(unitSystem=unitSystem)
            addEntry(name, value, unit)
            logging.debug('Added variable ' + name + ' to inputs.dat')
        else:
            value = db.get('multiphysicsModel['+modelName+']::'+name)
        
        return value
    except:
        if isVariable:
            logging.debug("Getting variable " +name.split(']')[0] + " failed.")
            return nominal
        else:
            raise NameError(name + ' not found in XML databse.')
 
def getVec3DCartesian(name,elements = ['_X','_Y','_Z']):
    vec = getVar(0,'vector3D[{0}]'.format(name), isVariable=False)
    val = [vec.x.float(unitSystem='SI'),vec.y.float(unitSystem='SI'),vec.z.float(unitSystem='SI')]
    units = [vec.x.units(unitSystem='SI'),vec.y.units(unitSystem='SI'),vec.z.units(unitSystem='SI')]
    for text, v, unit in zip(elements, val, units):
        addEntry(name+text, v, unit)

    return val

def getVec3DRotation(name, returnSOA=False):
    from Math import SOA_Py as soa
    
    R = getVar(0,'vector3D[{0}]'.format(name), isVariable=False)
    q = soa.SOAQuaternion()
    q.identity()
    qtemp = soa.SOAQuaternion()
    for rotation in R.rotations:
        axis = [0,0,0]
        axis[rotation.axis.value-1] = 1
        addEntry(name+"_"+rotation.name, rotation.angle.float(unitSystem=unitSystem),rotation.angle.units(unitSystem=unitSystem))
        qtemp.setAngleAxis(soa.SOAVector3(axis), rotation.angle.float(unitSystem=unitSystem))
        q = qtemp*q
    if returnSOA:
        return q
    else:
        return [qi for qi in q()] # Return as list rather than tuple
    

#Input log methods -  could potentially have their own file but they integrate seamlessly with the perturbUtils.
def inputLogInit():
    from collections import OrderedDict
    global inputLog
    inputLog = {}
    inputLog['values'] = OrderedDict()
    inputLog['units'] = OrderedDict()
    return None
    
def addEntry(name, value, unit):
    global inputLog
    #Add an entry to the inputLog that gets written at the start of the sim.
    inputLog['values'][name.split(']')[0]] = value
    if unit is None or 'dimensionless' in unit:
        inputLog['units'][name.split(']')[0]] = '-'
    else:
        inputLog['units'][name.split(']')[0]] = unit
    
    return None
    
def writeInputLog():
    import logging
    
    global inputLog
    logging.info('Writing ' + str(len(inputLog['values'].keys())) + ' Dakota inputs to file ./results/inputs.dat')
    spacing = 40
    file = open('./results/inputs.dat','w') # Create summary file to write to
    file.write('#Simulation inputs: '+ str(len(inputLog.keys())) +'\n')
    file.write('Input name'+ ' '*(spacing-10)) #10 is len('Input name')
  
    
    for key in inputLog['values']:
        file.write(key+ ' '*(spacing-len(key)))
    file.write('\n')
    file.write('Value'+ ' '*(spacing-5))
    for key in inputLog['values']:
        file.write(str(inputLog['values'][key]) + ' '*(spacing-len(str(inputLog['values'][key]))))
    file.write('\n')
    file.write('Units'+ ' '*(spacing-5))        
    for key in inputLog['values']:
        file.write(str(inputLog['units'][key]) + ' '*(spacing-len(str(inputLog['units'][key]))))
    file.close()
        
        
    return None