# Adapted from efficiency.py by M.Jacquart : https://github.com/marcjacquart/scifi_analysis

# To launch this script:
# Default : python3 scifiAnalysis.py -r <RUN_NUMBER> 
# Custom : python3 scifiAnalysis.py -pin <PATH TO INPUT FILE> -pg <PATH TO GEO FILE> -f <INPUTFILE> -g <GEOFILE>

# Folder with H6 runs : /eos/experiment/sndlhc/convertedData/commissioning-h6/
# Geometry file with H6 geometry : /eos/experiment/sndlhc/convertedData/commissioning-h6/geofile_sndlhc_H6.root
# TO DO : Create a "figures" folder where the script is executed to save them all here.

# This script:
#    - Imports the raw data and geomery files
#    - Uses the geometry to compute cluster positions
#    - Plots hits map for raw data, selected events
#    - Plots cluster size distribution for clustered data, selected events
#    - Computes and plots several scifi misalignments (station, within plane, rotations between planes, rotations between mats)
# Uses functions from analysisFunctions.py and plotFunctions.py

# General
from argparse import ArgumentParser
import os
import numpy as np
import matplotlib.pyplot as plt
import gc
# Root:
import ROOT
import rootUtils as ut
# SND:
import SndlhcGeo
import SndlhcTracking
# Custom functions defined in external file:
from analysisFunctions import (goodEvent, zPlaneArr, extendHits, distFit,
    customFitStatus, CreateClusters, ClusterperPlane, indexChanHit, separateMats)
from plotFunctions import (chi2Hist, planesHist, diffHist, allPlanesGauss, 
    diffPosHist, rotationAngle, hitMap, CluSizeDistrib)

# Script options:
# Indicate the steps at which the analysis is wanted : Raw, EvtSelect, Tracking
Steps = ["Tracking"] #["Raw", "EvtSelect", "Tracking"]

#Plotting options
HitMap = False #true to plot hit map
ClustSize = False #true to plot cluster size distribution

# Default : EventSelect and Tracking set to False, they are automatically set to True in the script if present in Steps
EventSelect = False
Tracking = False
if "Tracking" in Steps :
    fitReducedStations = True # True to loop 4-fit, 1-test.
    needXYInfo = True # Both vertical and horizontal planes must be hit (for rotation)
else :
    fitReducedStations = False # False to fit with all 5 stations
    needXYInfo = False


# Paths+name of the data and geometry root files as script arguments:
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int, required=False)
parser.add_argument("-pg", "--pathgeo", dest = "pathgeo", help = "path to geo file", required = False, default="/eos/experiment/sndlhc/testbeam/scifi/sndsw/") 
parser.add_argument("-pin", "--pathin", dest = "pathin", help = "path to input files", required = False, default="/eos/experiment/sndlhc/convertedData/commissioning-h6/")
parser.add_argument("-f", "--inputFile", dest = "inputFile", help = "input root file", required=False, default = "")
parser.add_argument("-g", "--geoFile", dest = "geoFile", help = "geometry file", required = False, default="geofile_full.Ntuple-TGeant4.root")
options = parser.parse_args()

# Open geometry and input files:
geo = SndlhcGeo.GeoInterface(options.pathgeo+options.geoFile)
if options.inputFile=="":
  f=ROOT.TFile.Open(options.pathin+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
  eventTree = f.rawConv #Extract the TTree: rawConv is the name of the tree in the root file
else:
  f=ROOT.TFile.Open(options.pathin+options.inputFile)
  eventTree = f.rawConv

# Global variables definitions:
lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
lsOfGlobals.Add(geo.modules['Scifi']) # only look at SciFi events

# GetCurrentNavigator() to browse trough the geometry files (TGeo class).
nav = ROOT.gGeoManager.GetCurrentNavigator()

# Setup tracking:
trackTask = SndlhcTracking.Tracking()
trackTask.InitTask(eventTree)
# For rotation plot to fill with each test station iteration:
if needXYInfo:
    #planes rotation
    xPosyOff_Slope = []
    xPosyOff_SlopeErr = []
    yPosxOff_Slope = []
    yPosxOff_SlopeErr = []
    #mats rotation
    xPosyOff_Slope_mat = []
    xPosyOff_SlopeErr_mat = []
    yPosxOff_Slope_mat = []
    yPosxOff_SlopeErr_mat = []

# Fit on 4 stations and use the 5th one as a test:
if fitReducedStations:
    gaussFitArr = []
    fitStationsArr = [[2,3,4,5],[1,3,4,5],[1,2,4,5],[1,2,3,5],[1,2,3,4]]
    testStationArr = [1,2,3,4,5]
else: # Only one loop iteration with all the stations to fit.
    fitStationsArr = [[1,2,3,4,5]]
    testStationArr =  [0] # Doesn't exist, but still need for plot names.

# z positions of the planes, array of 10 float
zArr = zPlaneArr(eventTree = eventTree, geo = geo)

for testStationNum in testStationArr:
    print(f'Test station: {testStationNum}') # Info of script speed
    chi2_nDfArr = [] # To fill in the loop for histogram
    nPlanesHit = [] # Same

    horDiffArr = []
    verDiffArr = []
    horPosArr = []
    verPosArr = []

    for step in Steps :
        print(f'Step : {step}')
        NbEvents = eventTree.GetEntries()
        print(f'Number of events : {NbEvents}')
        Events_const = NbEvents
        Event = 0
        AllClusters = []
        AllhitChanX = []
        AllhitChanY = []
        if step == "EvtSelect" : EventSelect = True
        if step == "Tracking" : 
            EventSelect = True
            Tracking = True
        
        print("Begining analysis")
        for sTree in eventTree :
            Event +=1
            #Clustering without tracking
            ClusterArr = CreateClusters(sTree = sTree)

            if EventSelect : # event selection
                ClupPlane = ClusterperPlane(ClusterArr)
                if not (goodEvent(eventTree = sTree, nStations = 5, allowMore = False) and len(ClupPlane) == 10 and all(cpp == 1 for cpp in ClupPlane)) : 
                    NbEvents -= 1
                    continue
            
            if Tracking : # tracking + quality selection
                fit, fitStatus = customFitStatus(trackTask = trackTask, FitStations = fitStationsArr[testStationNum-1])
                fitHitsExt = extendHits(fittedTrack = fit, zArr = zArr)
                if fitStatus.getNdf() == 0.0: # Fit has no degree of freedom
                    chi2_nDf = -1  # Impossible value to put them aside
                else: 
                    chi2_nDf = fitStatus.getChi2()/fitStatus.getNdf()
                #print(chi2_nDf)
                if not chi2_nDf<200 and chi2_nDf>=0 : # very large chi2 due to station 3 misalignment
                    NbEvents -= 1
                    continue

                chi2_nDfArr.append(chi2_nDf)
                if fitReducedStations: # Compute difference between hits and fit:
                    horDiff, verDiff, horPos, verPos , horCluSize, verCluSize, horEnergy, verEnergy = distFit(fitHits = fitHitsExt, clusterArr = trackTask.clusters, testStationNum = testStationNum)
                    if needXYInfo: # Select event if both planes are hit
                        if abs(horDiff) <1 and abs(verDiff) <1: # Search within 1cm from the fit
                            horDiffArr.append(horDiff)
                            horPosArr.append(horPos)
                            verDiffArr.append(verDiff)
                            verPosArr.append(verPos)
                    else: # Or at least one plane
                        if horDiff <1: # Don't append missing hits
                            horDiffArr.append(horDiff)
                            horPosArr.append(horPos)
                        if verDiff <1:
                            verDiffArr.append(verDiff)
                            verPosArr.append(verPos)
                
                # Delete objects from memory (don't know if really useful)
                del fit
                del fitStatus
                gc.collect()
            
            AllClusters.append(ClusterArr) #AllClusters = array of ClusterArr of each event
            
            hitChanX, hitChanY = indexChanHit(eventTree = sTree)
            AllhitChanX.append(hitChanX)
            AllhitChanY.append(hitChanY)

        # Hit Map    
        if HitMap : 
            print("Plotting Hit Maps")
            hitMap(hitChanX = AllhitChanX, hitChanY = AllhitChanY, run = str(options.runNumber), option = step, NbEvents = NbEvents)
        
        # Cluster Size distrib
        if ClustSize : 
            print("Plotting Cluster Size distribution")
            CluSizeDistrib(Clusters = AllClusters, run = str(options.runNumber), option = step, NbEvents = NbEvents)
        
        # Tracking related plots
        if Tracking : 
            print("Begining tracking plots")
            if fitReducedStations:
                print('Difference between hit and fit')
                resultFit = diffHist(horDiffArr = horDiffArr, verDiffArr = verDiffArr, stationNum = testStationNum, run = str(options.runNumber))
                gaussFitArr.append(resultFit) # in mm since diffHist() output changed
            
                # Difference between hit and fit vs position within the same plane:
                print('Vertical and horizontal shifts')
                diffPosHist(posArr = verPosArr, diffArr = verDiffArr, binsPos = np.linspace(-50,-5,90), run = str(options.runNumber), fileName = f'_OFFSET_ver_testSta{testStationNum}', labels = ['X residual [mm]', 'X cluster position [mm]'], title = f'Residual vs position in vertical plane - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = False)
                diffPosHist(posArr = horPosArr, diffArr = horDiffArr, binsPos = np.linspace(15,60,90), run = str(options.runNumber), fileName = f'_OFFSET_hor_testSta{testStationNum}', labels = ['Y residual [mm]', 'Y cluster position [mm]'], title = f'Residual vs position in horizontal plane - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = False,)
           
                if needXYInfo: # Cross terms to check rotation misalignments
                    print('Rotation between planes')
                    slope, slopeErr = diffPosHist(posArr = verPosArr, diffArr = horDiffArr, binsPos = np.linspace(-50,-5,90), run = str(options.runNumber), fileName = f'_ROT_horDiff_verPos_testSta{testStationNum}', labels = ['X cluster position [mm]', 'Y residual [mm]'], title = f'Rotation between horizontal and vertical planes (ver pos, hor residual) - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = True)
                    xPosyOff_Slope.append(slope)
                    xPosyOff_SlopeErr.append(slopeErr)
                    slope, slopeErr = diffPosHist(posArr = horPosArr, diffArr = verDiffArr, binsPos = np.linspace(15,60,90), run = str(options.runNumber), fileName = f'_ROT_verDiff_horPos_testSta{testStationNum}', labels = ['Y cluster position [mm]', 'X residual [mm]'], title = f'Rotation between horizontal and vertical planes (hor pos, ver residual) - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = True)
                    yPosxOff_Slope.append(slope)
                    yPosxOff_SlopeErr.append(slopeErr)
                    
                    #check mats rotations misalignment
                    print('Rotation between mats')
                    verPosArr1y, verPosArr2y, verPosArr3y, horPosArr1x, horPosArr2x, horPosArr3x, verDiffArr1x, verDiffArr2x, verDiffArr3x, horDiffArr1y, horDiffArr2y, horDiffArr3y = separateMats(verPosArr = verPosArr, horPosArr = horPosArr, verDiffArr = verDiffArr, horDiffArr = horDiffArr)
                    # X pos, Y residual : one plot for each y mat, slopes and err filled for each mat then for each test station
                    slope, slopeErr = diffPosHist(posArr = verPosArr1y, diffArr = horDiffArr1y, binsPos = np.linspace(-50,-5,90), run = str(options.runNumber), fileName = f'_matROT_horDiffmat1_verPos_testSta{testStationNum}', labels = ['X cluster position [mm]', 'Y residual [mm]'], title = f'Rotation between horizontal mat 1 and vertical plane - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = True)
                    xPosyOff_Slope_mat.append(slope)
                    xPosyOff_SlopeErr_mat.append(slopeErr)
                    slope, slopeErr = diffPosHist(posArr = verPosArr2y, diffArr = horDiffArr2y, binsPos = np.linspace(-50,-5,90), run = str(options.runNumber), fileName = f'_matROT_horDiffmat2_verPos_testSta{testStationNum}', labels = ['X cluster position [mm]', 'Y residual [mm]'], title = f'Rotation between horizontal mat 2 and vertical plane - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = True)
                    xPosyOff_Slope_mat.append(slope)
                    xPosyOff_SlopeErr_mat.append(slopeErr)
                    slope, slopeErr = diffPosHist(posArr = verPosArr3y, diffArr = horDiffArr3y, binsPos = np.linspace(-50,-5,90), run = str(options.runNumber), fileName = f'_matROT_horDiffmat3_verPos_testSta{testStationNum}', labels = ['X cluster position [mm]', 'Y residual [mm]'], title = f'Rotation between horizontal mat 3 and vertical plane - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = True)
                    xPosyOff_Slope_mat.append(slope)
                    xPosyOff_SlopeErr_mat.append(slopeErr)
                    # Y pos, X residual
                    slope, slopeErr = diffPosHist(posArr = horPosArr1x, diffArr = verDiffArr1x, binsPos = np.linspace(15,60,90), run = str(options.runNumber), fileName = f'_matROT_verDiffmat1_horPos_testSta{testStationNum}', labels = ['Y cluster position [mm]', 'X residual [mm]'], title = f'Rotation between vertical mat 1 and horizontal plane - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = True)
                    yPosxOff_Slope_mat.append(slope)
                    yPosxOff_SlopeErr_mat.append(slopeErr)
                    slope, slopeErr = diffPosHist(posArr = horPosArr2x, diffArr = verDiffArr2x, binsPos = np.linspace(15,60,90), run = str(options.runNumber), fileName = f'_matROT_verDiffmat2_horPos_testSta{testStationNum}', labels = ['Y cluster position [mm]', 'X residual [mm]'], title = f'Rotation between vertical mat 2 and horizontal plane - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = True)
                    yPosxOff_Slope_mat.append(slope)
                    yPosxOff_SlopeErr_mat.append(slopeErr)
                    slope, slopeErr = diffPosHist(posArr = horPosArr3x, diffArr = verDiffArr3x, binsPos = np.linspace(15,60,90), run = str(options.runNumber), fileName = f'_matROT_verDiffmat3_horPos_testSta{testStationNum}', labels = ['Y cluster position [mm]', 'X residual [mm]'], title = f'Rotation between vertical mat 3 and horizontal plane - TestStation {testStationNum} - Run' + str(options.runNumber), isCrossed = True)
                    yPosxOff_Slope_mat.append(slope)
                    yPosxOff_SlopeErr_mat.append(slopeErr)
                    
            print('Chi2 histogram')
            chi2Hist(chi2_nDfArr=chi2_nDfArr, stationNum=testStationNum, run = str(options.runNumber))
        
if Tracking :
    if fitReducedStations:
        # Translation misalignment of each station:
        allPlanesGauss(fitArr=gaussFitArr, run = str(options.runNumber))

    if needXYInfo:
        
        # Rotation misalignment of each station:
        rotationAngle(xPosyOff_Slope = xPosyOff_Slope, xPosyOff_SlopeErr = xPosyOff_SlopeErr, yPosxOff_Slope = yPosxOff_Slope, yPosxOff_SlopeErr = yPosxOff_SlopeErr, mats = False, run = str(options.runNumber))

        # Rotation misalignment of mat for each station:
        rotationAngle(xPosyOff_Slope = xPosyOff_Slope_mat, xPosyOff_SlopeErr = xPosyOff_SlopeErr_mat, yPosxOff_Slope = yPosxOff_Slope_mat, yPosxOff_SlopeErr = yPosxOff_SlopeErr_mat, mats = True, run = str(options.runNumber))
       
