#!/usr/bin/env python

if __name__ == "__main__":

    from numpy import pi, loadtxt
    import os
    import imuQ
    import sys
    from os import path
    sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )

    #  Load the database

    db = imuQ.Project.loadDB()

    #  Define the responses

    responses = db.get('ModelDB[MarsEDL]::pythonProtocol[TimeSimulation]:responses')
    try:

        with open('./results/trajectory.txt') as myfile:
            head = [next(myfile) for x in xrange(3)]
        
            varnames = head[1].split()
            # varunits = head[2].split()
            data = loadtxt(myfile, skiprows=3)
      
        get = lambda name,data: data[-1,varnames.index(name)]    
  
                  
        #  Compute the system response quantities
        
        for SRQName in responses.responseFunctions.keys():
                          
            responses.responseFunctions[SRQName].value = get(SRQName,data)
    
        #  Output the responses
        
        responses.printDakotaResultsFile()
        
    except Exception as e:

        #  Output the responses with a FAIL flag for DAKOTA
    
        responses.printDakotaFailureFile(e)

    print "Trajectory simulation complete."

