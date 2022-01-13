# Modified from Marc Jacquart : https://github.com/marcjacquart/scifi_analysis
# Some functions added or modified

import ROOT
import numpy as np

def goodEvent(eventTree, nStations, allowMore):
    '''Return True if nStations (or more if allowMore) are hits in the event'''

    stations = {}
    for detectedHit in eventTree.Digi_ScifiHits:
        stations[detectedHit.GetDetectorID()//1000000] = 1

    # for detectedHit in eventTree.Digi_MuFilterHit:
    #     plane = 100*(detectedHit.GetDetectorID()//1000)
    #     stations[plane] = 1
    if allowMore:
        if len(stations) >= nStations: return True
    else:
        if len(stations) == nStations: return True
    return False

def indexStationsHit(eventTree):
    '''
    return array [X_1,X_2,X_3,...] of planes hit.
    X: 0-9 number of the plane hit. With same convention as zArr:
    2 * (plane number-1) + 1 for vertical plane
    i.e. plane 2 vertical would have number 3 here. 
    '''

    # Filling a dictionary allows overwritting to have one entery per plane.
    stations = {} 
    for detectedHit in eventTree.Digi_ScifiHits:
        detID = detectedHit.GetDetectorID()

        planeNumber = (detID // 1000000 - 1) * 2  + (detID // 100000) % 2
        stations[detID//100000] = planeNumber
    
    # Fill array to return from dict values:
    indexHitArr = [] 
    for planeID in stations.values():
        indexHitArr.append(planeID)
    return indexHitArr

def indexChanHit(eventTree): #added
    '''
    returns arrays [X_1,X_2,X_3,...] and [Y_1, Y_2, Y_3,...] of channel hit.
    X,Y: 0-1535 number of the channel hit. Convention : 
    0-127 : Station , mat 0  SIPM 0
    128-255 : Station , mat 0, SIPM 0 etc.
    1408-1535 : Station , mat 2, SIPM 3 
    '''
    stat1x, stat1y, stat2x, stat2y, stat3x, stat3y, stat4x, stat4y, stat5x, stat5y = [],[],[],[],[],[],[],[],[],[]
    chanHitx = [stat1x, stat2x, stat3x, stat4x, stat5x] #elem N-1 corresponds to station N
    chanHity = [stat1y, stat2y, stat3y, stat4y, stat5y]
    Stations = [1,2,3,4,5]

    for Hit in eventTree.Digi_ScifiHits:
        hitID = Hit.GetDetectorID()
        chanIndex = hitID - (hitID//1000)*1000
        SIPMIndex = ((hitID//1000) - (hitID//10000)*10)
        matIndex = ((hitID//10000) - (hitID//100000)*10)
        chanHitidx = chanIndex + 128*SIPMIndex + 128*4*matIndex
        if (hitID //100000)%2 == 0 : #horizontal fibers -> Y index
            for stat in Stations :
                if (hitID//1000000) == stat :
                    chanHity[stat-1].append(chanHitidx)
        if (hitID //100000)%2 == 1 : #vertical fibers -> X index
            for stat in Stations :
                if (hitID//1000000) == stat :
                    chanHitx[stat-1].append(chanHitidx)

    return chanHitx, chanHity

def extendHits(fittedTrack, zArr):
    '''
    Extend the fit position to include missing planes hits.
    fittedTrack: ROOT.genfit.Track, obtained with customFitStatus()
    zArr: array containing z position (float) of each plane. 
    Return array of 10 TVector3, aech containing the fit intersection
    with one plane.
    '''
    fitHits =[ROOT.TVector3()]*10 # Points where fit crosses the 10 planes.
    # for i in range(fittedTrack.getNumPointsWithMeasurement()):
    # Only the first can be used, other hits give same line
    state = fittedTrack.getFittedState(0)
    pos = state.getPos()
    mom = state.getMom()
    # linear fit: pos + lambda * mom
    for planeIndex in range(len(zArr)):
        lambdaPlane = (zArr[planeIndex] - pos[2]) / mom[2]
        fitHits[planeIndex] = pos + lambdaPlane * mom
    return fitHits

def crossAllPlanes(fitHitsArr,geo, verbose=False):
    '''
    Return true if fit is within the detector boundaries.
    Tell where the fit exit the detector if verbose activated
    fitHitsArr: coordinates of fit-plane intersection from extendHits().
    geo: SndlhcGeo.GeoInterface(<geometry file>) for detector boundaries.
    '''
    isInside = True # Set to False if only once out of bounds
    A,B = ROOT.TVector3(),ROOT.TVector3()

    #First plane horizontal:
    geo.modules['Scifi'].GetSiPMPosition(1000000,A,B)
    if fitHitsArr[0][0]<B[0] or A[0]<fitHitsArr[0][0]:
        isInside = False
        if verbose:
            print('first x border exceeded:')
            print(f'{fitHitsArr[0][0]}<{B[0]} or {A[0]}<{fitHitsArr[0][0]}')

    #First plane vertical:
    geo.modules['Scifi'].GetSiPMPosition(1100000,A,B)
    if fitHitsArr[1][1]<A[1] or B[1]<fitHitsArr[1][1]:
        isInside = False
        if verbose:
            print('first y border exceeded')
            print(f'{fitHitsArr[1][1]}<{A[1]} or {B[1]}<{fitHitsArr[1][1]}')

    #Last plane horizontal:
    geo.modules['Scifi'].GetSiPMPosition(5000000,A,B)
    if fitHitsArr[8][0]<B[0] or A[0]<fitHitsArr[8][0]:
        isInside = False
        if verbose:
            print('last x border exceeded')
            print(f'{fitHitsArr[8][0]}<{B[0]} or {A[0]}<{fitHitsArr[8][0]}')

    #Last plane vertical:
    geo.modules['Scifi'].GetSiPMPosition(5100000,A,B)
    if fitHitsArr[9][1]<A[1] or B[1]<fitHitsArr[9][1]:
        isInside = False
        if verbose:
            print('last y border exceeded')
            print(f'{fitHitsArr[9][1]}<{A[1]} or {B[1]}<{fitHitsArr[9][1]}')

    return isInside

def zPlaneArr(eventTree,geo):
    '''
    Return array with the 10 z coordinates of the planes.
    format: [1hor, 1ver, 2hor, ...,5hor, 5ver]
    /!\ THIS ASSUME Z FIX FOR ALL THE PLANE!
    geo: SndlhcGeo.GeoInterface(<geometry file>) for planes positions.
    eventTree: event list to get one event crossing all planes 
    to get each plane position.
    '''
    zArr = [0,0,0,0,0,0,0,0,0,0] # Fill Z coordinate of the planes
    # Warning message if z values > epsilon for different SiPM of same plane
    epsilon = 0.0001 

    for sTree in eventTree: # Need the first event with all planes hit:
        if goodEvent(eventTree = eventTree, nStations = 5, allowMore = False):
            A,B = ROOT.TVector3(),ROOT.TVector3()
            for hit in sTree.Digi_ScifiHits:
                hitID = hit.GetChannelID()
                geo.modules['Scifi'].GetSiPMPosition(hitID,A,B)
                indexArr = (hitID // 1000000-1) * 2  + (hitID // 100000) % 2
                #           add 2 for each plane        add 1 if vertical

                zVal = 0.5 * (A[2] + B[2])
                #Will overwrite if event of same plane, good to check if same z.
                if zArr[indexArr]!=0 and abs(zArr[indexArr]-zVal)>epsilon:
                    print(f'WARNING: SciFi planes {indexArr} not aligned in Z direction!')
                zArr[indexArr] = zVal # Fill z array
            break
    return zArr

def sortHitStation(clusterArr,stationArr):
    '''
    return the array of clusters from clusterArr that belong
    to one of the selected stations of stationArr
    stationArr: Array with selected stations number, ex: [1,3,4,5]
    '''
    
    clusFit = []
    clusTest = []
    for cluster in clusterArr:
        # If the cluster station is inside our list to fit:
        if (cluster.GetFirst()//1000000) in stationArr:
            clusFit.append(cluster)
        else:
            clusTest.append(cluster)
    return clusFit, clusTest

def distFit(fitHits, clusterArr, testStationNum): #modified
    '''
    Takes the cluster list of an event (clusterArr), the fit
    (constructed with the 4 stations) hits on all the the plane
    fitHits which is a 10 elements array and the 5th test station
    (testStationNum) to compute the distance between the test station's
    hits and the predicted position by the fit on the same plane.
    Missing hits are set to 1000, easy to remove afterward.

    fitHits: array of TVector3 plane-fit intersection from extendHits()
    Return horizontal fit-hit difference, vertical fit-hit difference,
    horizontal position of the selected hit, vertical position of the selected hit,
    horizontal and vertical cluster sizes
    horizontal and vertical cluster energy
    '''
    horIndex = 2 * (testStationNum - 1)
    verIndex = horIndex +1
    horClusters = [1000] # Large numer so it will be the minimum
    verClusters = [1000] # only if the array is empty
    horCluSizeArr = [0] # if min(Clusters)==1000 (ie, no hit)-> size = 0
    verCluSizeArr = [0]
    horEnergyArr = [0] # if min(CLusters)==1000 (ie, no hit) -> energy = 0
    verEnergyArr = [0]
    #Fill horizontal and vertical cluster lists of the test plane:
    for cluster in clusterArr:
        if (cluster.GetFirst()//1000000) == testStationNum:
            A,B = ROOT.TVector3(),ROOT.TVector3()
            cluster.GetPosition(A,B)
            if ((cluster.GetFirst()//100000) % 2 ) == 0: # horizontal
                horClusters.append(A[1])
                horCluSizeArr.append(cluster.GetN())
                horEnergyArr.append(cluster.GetEnergy())
            else:
                verClusters.append(A[0])
                verCluSizeArr.append(cluster.GetN())
                verEnergyArr.append(cluster.GetEnergy())

    horFit = fitHits[horIndex][1]
    verFit = fitHits[verIndex][0]
    # We want the differece (no abs) of the closest point(in absolute diff)
    horDiffIndex = np.argmin([abs(x - horFit) for x in horClusters]) # Minimal distance
    verDiffIndex = np.argmin([abs(y - verFit) for y in verClusters]) # argmin() returns index of min

    horPos = horClusters[horDiffIndex]
    verPos = verClusters[verDiffIndex]
    horCluSize = horCluSizeArr[horDiffIndex]
    verCluSize = verCluSizeArr[verDiffIndex]
    horEnergy = horEnergyArr[horDiffIndex]
    verEnergy = verEnergyArr[verDiffIndex]
    horDiff = horPos - horFit
    verDiff = verPos - verFit
    return horDiff, verDiff, horPos , verPos , horCluSize, verCluSize, horEnergy, verEnergy

def testClusterProblem(eventTree):
    ''' 
    Test for cluster separated by only one unactivated SiPM channel.
    Print the result for easily see the neighboring clusters.
    eventTree: Event list from root file.
    return: None
    '''
    print('##################')
    IDArr = []
    for hit in eventTree.Digi_ScifiHits:
        detID = int(hit.GetDetectorID())
        IDArr.append(detID)
    IDArr.sort()
    prevElement=0
    for element in IDArr:
        diff = element - prevElement
        if not (diff ==1):
            if diff < 10:
                print(f'!!!---{element - prevElement}---!!!')
            else:
                print('---')
        print(element)
        prevElement = element
        
def customFitStatus(trackTask, FitStations):
    '''
    Do manually the trackTask.ExecuteTask() so it can use arbitrary number
    of stations for the fit.
    Return the fit and fit status object containing the fitted track.
    stationArr: Array with selected stations number for the fit, ex: [1,3,4,5]
    '''
    trackTask.clusters = trackTask.scifiCluster() # Build cluster list
    clusFit, clusTest = sortHitStation(clusterArr=trackTask.clusters,stationArr=FitStations)
    # Do the manual fit with the reduced number of stations:
    # Need manual dictionnary for the fitTrack(hitlist):
    clusDict = {}
    for x in clusFit:
        A,B = ROOT.TVector3(),ROOT.TVector3()
        clusterID = x.GetFirst()
        dictEntery = {clusterID: x}
        clusDict.update(dictEntery)
    fit = trackTask.fitTrack(hitlist=clusDict) # Type ROOT.genfit.Track
    fitStatus= fit.getFitStatus()
    trackTask.event.fittedTracks = [fit] # Array to keep sndsw format
    return fit, fitStatus

def CreateClusters(sTree) : #added
    '''
    Creates the clusters for event sTree.
    Returns array of clusters, a cluster is an array of channels
    Ex : [[array of hor channels plane1],[array of ver channels plane1],[array of hor channels plane2],[array of ver channels plane2]]
    Does not need tracking contrary to scificluster()
    '''
    IDArr = []
    Cluster = []
    ClusterArr = []
    for hit in sTree.Digi_ScifiHits :
        detID = int(hit.GetDetectorID())
        IDArr.append(detID)
    IDArr.sort()
    prevElem=0
    for elem in IDArr :
        diff = elem - prevElem
        if diff ==1 or elem == IDArr[0]:
            Cluster.append(elem) 
            if elem == IDArr[-1] : ClusterArr.append(Cluster)
        else :
            ClusterArr.append(Cluster) #See what happens for first element
            Cluster = []
            Cluster.append(elem)
        prevElem = elem
    return ClusterArr

def ClusterperPlane(ClusterArr) : #added
    '''
    ClusterArr : sorted array of clusters for one event 
    Returns the number of cluster hits per plane per event in array of length 10
    '''
    clu10, clu11, clu20, clu21, clu30, clu31, clu40, clu41, clu50, clu51 = 0,0,0,0,0,0,0,0,0,0
    CluperPlane = [clu10, clu11, clu20, clu21, clu30, clu31, clu40, clu41, clu50, clu51]
    #count = 0
    #prevClu = 10
    PlaneArr = [10,11,20,21,30,31,40,41,50,51]
    for cluster in ClusterArr :
        for plane in PlaneArr :
            if (cluster[0]//100000) == plane :
                CluperPlane[PlaneArr.index(plane)] += 1

    return CluperPlane

def separateMats(verPosArr, horPosArr, verDiffArr, horDiffArr) : #added
    '''
    Separates position arrays into 3 (for each mat) and assign the Diffs to the corresponding mat
    Recall : X pos ~ (-47;-7), Y pos ~ (15;55) [cm]
    Cuts a bit inside the mats to prevent from confusion at the edges
    '''
    verPosArr1y, verPosArr2y, verPosArr3y = [],[],[] # verPosArr of events in 1st/2nd/3rd hor(==y) mat
    horPosArr1x, horPosArr2x, horPosArr3x = [],[],[] # horPosArr of events in 1st/2nd/3rd ver(==x) mat
    verDiffArr1x, verDiffArr2x, verDiffArr3x = [],[],[]
    horDiffArr1y, horDiffArr2y, horDiffArr3y = [],[],[]
    # Separation into the 3 vertical mats
    for X in verPosArr : 
        if X > -47 and X < -34.5 :
            horPosArr1x.append(horPosArr[verPosArr.index(X)])
            verDiffArr1x.append(verDiffArr[verPosArr.index(X)])
        if X > -33.25 and X < -20.75 :
            horPosArr2x.append(horPosArr[verPosArr.index(X)])
            verDiffArr2x.append(verDiffArr[verPosArr.index(X)])
        if X > -19.5 and X < -7 :
            horPosArr3x.append(horPosArr[verPosArr.index(X)])
            verDiffArr3x.append(verDiffArr[verPosArr.index(X)])

    # Separation into the 3 horizontal mats
    for Y in horPosArr : 
        if Y > 15 and Y < 27.7 :
            verPosArr1y.append(verPosArr[horPosArr.index(Y)])
            horDiffArr1y.append(horDiffArr[horPosArr.index(Y)])
        if Y > 28.75 and Y < 41.25 :
            verPosArr2y.append(verPosArr[horPosArr.index(Y)])
            horDiffArr2y.append(horDiffArr[horPosArr.index(Y)])
        if Y > 42.5 and Y < 55 :
            verPosArr3y.append(verPosArr[horPosArr.index(Y)])
            horDiffArr3y.append(horDiffArr[horPosArr.index(Y)])

    return verPosArr1y, verPosArr2y, verPosArr3y, horPosArr1x, horPosArr2x, horPosArr3x, verDiffArr1x, verDiffArr2x, verDiffArr3x, horDiffArr1y, horDiffArr2y, horDiffArr3y

def trackSlope(fitHits) : #added
    '''
    fitHits : array of 10 vector3D with track position on each plane
    fitHits ~ [hor1, ver1, hor2, ver2 ...] with 1,2... station number
    returns x slope, y slope of track
    '''
    xslope = (fitHits[9][0] - fitHits[1][0])/(fitHits[9][2] - fitHits[1][2])
    yslope = (fitHits[8][1] - fitHits[0][1])/(fitHits[8][2] - fitHits[0][2])

    return xslope, yslope

