# To launch this script:
# Default : python3 scifiEfficiency.py -r <RUN_NUMBER> 
# Custom : python3 scifiEfficiency.py -pin <PATH TO INPUT FILE> -pg <PATH TO GEO FILE> -f <INPUTFILE> -g <GEOFILE>

# Folder with H6 runs : /eos/experiment/sndlhc/convertedData/commissioning-h6/
# Geometry file with H6 geometry : /eos/experiment/sndlhc/convertedData/commissioning-h6/geofile_sndlhc_H6.root
# TO DO : Create a "figures" folder where the script is executed to save them all here.

# This script:
#  if DoAnalysis == True :
#    - Imports the raw data and geomery files
#    - Uses the geometry to compute cluster positions
#    - Saves quantities of interest in dataframe for given run and testStation
#  if DoAnalysis == False :
#    - Loads dataframe for given run and testStation
#  in both cases :
#    - Plots efficiency as function of position
#    - Plots efficiency with zoom onf physical gaps between and in SiPMs
#    - Plots histograms of tracks reconstructed
#    - Plots distributions of slopes, cluster sizes, energy deposit
# Uses functions from analysisFunctions.py and plotFunctions.py

# General
from argparse import ArgumentParser
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import gc
# Root:
import ROOT
import rootUtils as ut
# SND:
import SndlhcGeo
import SndlhcTracking
# Custom functions defined in external file:
from analysisFunctions import (zPlaneArr, extendHits, distFit,
    customFitStatus, CreateClusters, ClusterperPlane, trackSlope)
from plotFunctions import (Efficiency, SIPMgap, trackHistograms, slopeEfficiency)

# Paths+name of the data and geometry root files as script arguments:
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int, required=False)
parser.add_argument("-pg", "--pathgeo", dest = "pathgeo", help = "path to geo file", required = False, default="/eos/experiment/sndlhc/testbeam/scifi/sndsw/")
parser.add_argument("-pin", "--pathin", dest = "pathin", help = "path to input files", required = False, default="/eos/experiment/sndlhc/convertedData/commissioning-h6/")
parser.add_argument("-f", "--inputFile", dest = "inputFile", help = "input root file", required=False, default = "")
parser.add_argument("-g", "--geoFile", dest = "geoFile", help = "geometry file", required = False, default="geofile_full.Ntuple-TGeant4.root")
options = parser.parse_args()

# Script options
run = str(options.runNumber)
DoAnalysis = False #set to False to just load file with analysis and plot
TestStation = 3
FitStations = [1,2,4,5]

if DoAnalysis : 
    #Load and setup everything ============================================================================================ 
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

    # z positions of the planes, array of 10 float
    zArr = zPlaneArr(eventTree = eventTree, geo = geo)

    #===========================================================================================================================
    # Begin analysis

    TestPlanes = [2*(TestStation-1), 2*(TestStation-1)+1]
    FitPlanes = []
    for FitStat in FitStations :
        FitPlanes.append(2*(FitStat-1))
        FitPlanes.append(2*(FitStat-1)+1)
    keepEvent = False

    #Creates empty DataFrame
    cols = ['FitPos_hor', 'FitPos_ver',
            'FitSlope_y', 'FitSlope_x', 
            'CluPos_hor', 'CluPos_ver', 
            'CluSize_hor', 'CluSize_ver', 
            'CluEnergy_hor', 'CluEnergy_ver',
            'chi2_Ndof', 
            'Clu_match_hor', 'Clu_match_ver', 
            'Residual_hor', 'Residual_ver']
    df_Events = pd.DataFrame(columns = cols)

    # Loop over all events
    for sTree in eventTree :

        ClusterArr = CreateClusters(sTree = sTree)
        ClupPlane = ClusterperPlane(ClusterArr = ClusterArr) 
        # Event selection : 1 cluster per plane per event for all fit stations
        if not (all(ClupPlane[plane] == 1 for plane in FitPlanes)) : continue 
        if (ClupPlane[TestPlanes[0]] == 0) : match_hor = 0 # no horizontal match
        else : match_hor = 1
        if (ClupPlane[TestPlanes[1]] == 0) : match_ver = 0 # no vertical match
        else : match_ver = 1

        fit, fitStatus = customFitStatus(trackTask = trackTask, FitStations = FitStations)
        fitHitsExt = extendHits(fittedTrack = fit, zArr = zArr) #array of 10 Vector3D with coordinates of fit on each plane

        if fitStatus.getNdf() == 0.0: # Fit has no degree of freedom
            chi2_nDf = -1  # Impossible value to put them aside
        else: 
            chi2_nDf = fitStatus.getChi2()/fitStatus.getNdf()
    
        horFit = fitHitsExt[TestPlanes[0]][1]
        verFit = fitHitsExt[TestPlanes[1]][0] 
        xslope, yslope = trackSlope(fitHits = fitHitsExt)
        horDiff, verDiff, horPos, verPos, horCluSize, verCluSize, horEnergy, verEnergy = distFit(fitHits = fitHitsExt, clusterArr = trackTask.clusters, testStationNum = TestStation)

        #fill dataframe
        newEvent = {'FitPos_hor' : horFit, 'FitPos_ver': verFit,
                    'FitSlope_y': yslope, 'FitSlope_x': xslope,
                    'CluPos_hor': horPos, 'CluPos_ver': verPos,
                    'CluSize_hor': horCluSize, 'CluSize_ver': verCluSize,
                    'CluEnergy_hor': horEnergy, 'CluEnergy_ver': verEnergy,
	                'chi2_Ndof': chi2_nDf,
	                'Clu_match_hor': match_hor, 'Clu_match_ver': match_ver,
                    'Residual_hor': horDiff, 'Residual_ver': verDiff}
        df_Events = df_Events.append(newEvent, ignore_index=True)

        # Delete objects from memory (don't know if really useful)
        del fit
        del fitStatus

        gc.collect()

    # Save DataFrame in csv file    
    df_Events.to_csv(f"Run{run}_Events_testStat{TestStation}.csv", index=True)


else :
    #load file with DataFrame
    df_Events = pd.read_csv(f"Run{run}_Events_testStat{TestStation}.csv")

# Plots
#Cut on chi2/Ndof (different for station 3 because of station misalignment)
    if TestStation == 3 :
        df_Events = df_Events[df_Events['chi2_Ndof'] < 30]
    else : 
        df_Events = df_Events[df_Events['chi2_Ndof'] < 200]

#Efficiency as function of position for hor, ver and bot matches + rescale between 0.8 and 1
Efficiency(df_Events = df_Events, testStation = TestStation, run = run)

#Reconstructed tracks histogram for 2 SiPMs + corresponding efficency
SIPMgap(df_Events = df_Events, testStation = TestStation, run = run)

#Distributions of reconstructed tracks (distrib=True), energy deposit(signal=True), cluster size (clusize=True), slopes (slopes=True)
# All bool == True by default
trackHistograms(df_Events = df_Events, testStation = TestStation, run = run, distrib = True, signal = False, clusize = False)

#Efficiency plots with cut on slopes to select beam
slopeEfficiency(df_Events = df_Events, testStation = TestStation, run = run)

