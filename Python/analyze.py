import sys
from os import path
sys.path.append( path.join(path.realpath('.'), '../../../MCP') )
from math import pi
import numpy as np
import matplotlib.pyplot as plt
import MonteCarlo as mc
import MonteCarloAnalysis as mca

method = 'lhs'
date = '080716'
nSwitchTimesAlpha = ['two','three']
nSwitchTimes = ['2','3']
# nSwitchTimesAlpha = ['three']
# nSwitchTimes = ['3']
data = [mc.MonteCarloData(loadFile = './data/{0}_{1}p_{2}'.format(method,n,date)) for n in nSwitchTimes]

# states = ['altitude','velocity','longitude','latitude','DR','CR']
inputs = ['CD','CL','densitySH','density0']
states = ['time','altitude','velocity','DR','CR']
events = ['Init','Final']

pairs = {state : (state,'Final') for state in states}
plist = [(state,'Final') for state in states]
ds = {}
# mc.compareConvergence(data, nSwitchTimes, nBins = 10, nSamples = 100, percentiles = [1,5,10,50,90,95,99], title=None, maxSize = 2000):


# prop = Fuel.derive([a+b for a,b in zip(O2.data,Fuel.data)])
# rcs_prop = lhsData.getDataSet('N2_remaining','OS_released')
# rcs_prop = rcs_prop.derive([p-0.379 for p in rcs_prop.data]) #correction

# Some basic plots
for nSwitch,d in zip(nSwitchTimesAlpha,data):
    # mca.jointPlot(d,pairs['velocity'],pairs['altitude'],title='{0} Switching Times'.format(nSwitch.capitalize()))
    # plt.savefig('./data/results/{0}-AltVel'.format(nSwitch))
    # mca.jointPlot(d,pairs['CR'],pairs['DR'],draw_marginal_histogram=False,title='{0} Switching Times'.format(nSwitch.capitalize()))
    # plt.savefig('./data/results/{0}-Range'.format(nSwitch))
    # mca.generateDataSetReport(d,plist, savePath = './data/results/', MCDataName = nSwitch)
    
    print nSwitch
    for s in states:
        dataset = d.getDataSet(s,'Final')
        print s
        print dataset.getMean()
        print dataset.getPercentile(1)
        print dataset.getPercentile(99)
        ds[s] = dataset
    
    cd = d.getInputSet('CD')
    cl = d.getInputSet('CL')
    
    dbar = [case.summaryData['dbar']['value'] for case in d.data]
    lbar = [case.summaryData['lbar']['value'] for case in d.data]
    lodbar = [case.summaryData['lodbar']['value'] for case in d.data]
    ts1 = [case.summaryData['ts1']['value'] for case in d.data]
    ts2 = [case.summaryData['ts2']['value'] for case in d.data]
    ts3 = [case.summaryData['ts3']['value'] for case in d.data]
    es = [[case.summaryData['es{0}'.format(switch)]['value'] for case in d.data] for switch in [1,2,3]]

    
    plt.scatter(dbar,lbar,c=ds['time'].data,s=100,alpha=0.5)
    plt.title('Final Time as Function of L,D for {0} switches'.format(nSwitch))
    plt.colorbar()
    plt.show()
    
    
    
    # plt.scatter(dbar,lbar,c=ds['altitude'].data,s=100,alpha=0.5)
    # plt.colorbar()
    # plt.grid(True)
    # plt.xlabel('Dbar',fontsize=16)
    # plt.ylabel('Lbar',fontsize=16)
    # plt.title('Colored by Final Altitude [km]')
    # plt.savefig('./{0}-DLTS1'.format(nSwitch))
    # plt.show()
    if False:
        if nSwitch == 'three':
            plt.scatter(dbar,lbar,c=ts1,s=100,alpha=0.5)
            plt.colorbar()
            plt.grid(True)
            plt.xlabel('Dbar',fontsize=16)
            plt.ylabel('Lbar',fontsize=16)
            # plt.title('Colored by Final Altitude [km]')
            plt.title('Colored by Switch Time 1 ({0} switches total)'.format(nSwitch))
            # plt.axis([-15, 375,0,95])
            plt.savefig('./{0}-DLTS1'.format(nSwitch))
            plt.show()
        
        plt.scatter(dbar,lbar,c=ts2,s=100,alpha=0.5)
        plt.colorbar()
        plt.grid(True)
        plt.xlabel('Dbar',fontsize=16)
        plt.ylabel('Lbar',fontsize=16)
        plt.title('Colored by Switch Time 2 ({0} switches total)'.format(nSwitch))
        # plt.axis([-15, 375,0,95])
        # plt.savefig('./{0}-DLTS2'.format(nSwitch))
        plt.show()
        
        plt.scatter(dbar,lbar,c=ts3,s=100,alpha=0.5)
        plt.colorbar()
        plt.grid(True)
        plt.xlabel('Dbar',fontsize=16)
        plt.ylabel('Lbar',fontsize=16)
        plt.title('Colored by Switch Time 3 ({0} switches total)'.format(nSwitch))
        # plt.axis([-15, 375,0,95])
        # plt.savefig('./{0}-DLTS3'.format(nSwitch))
        plt.show() 

        # plt.scatter(ts2,ts3,c=lodbar,s=100,alpha=0.5)
        # plt.colorbar()
        # plt.grid(True)
        # plt.xlabel('Switch 2',fontsize=16)
        # plt.ylabel('Switch 3',fontsize=16)
        # plt.title('Colored by lodbar')
        # plt.axis([-15, 375,0,95])
        # plt.savefig('./{0}-ts1ts2lod'.format(nSwitch))
        # plt.show()        
        
    for energy in es:    
        plt.scatter(dbar,lbar,c=energy,s=100,alpha=0.5)
        plt.colorbar()
        plt.grid(True)
        plt.xlabel('Dbar',fontsize=16)
        plt.ylabel('Lbar',fontsize=16)
        plt.title('Colored by Switch Energy {0} ({1} switches total)'.format(es.index(energy),nSwitch))
        # plt.axis([-15, 375,0,95])
        # plt.savefig('./{0}-DLTS3'.format(nSwitch))
        plt.show() 