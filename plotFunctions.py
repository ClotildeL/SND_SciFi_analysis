# Modified from Marc Jacquart : https://github.com/marcjacquart/scifi_analysis
# Some functions added or modified


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit

# Custom legend:
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerBase

class HandlerColormap(HandlerBase):
    '''
    Custom handles for cluster legend,
    Adapted from stackoverflow.com/questions/55501860/
    '''
    def __init__(self, cmap, num_stripes=10, **kw):
        HandlerBase.__init__(self, **kw)
        self.cmap = cmap
        self.num_stripes = num_stripes
    def create_artists(self, legend, orig_handle, 
                       xdescent, ydescent, width, height, fontsize, trans):
        stripes = []
        for i in range(self.num_stripes):
            s = Rectangle([xdescent + i * width / self.num_stripes, ydescent], 
                          width / self.num_stripes, 
                          height, 
                          fc=self.cmap((2 * i + 1) / (2 * self.num_stripes)), 
                          transform=trans)
            stripes.append(s)
        return stripes


def valueToColor(value, cmap_name='nipy_spectral', vmin=-171, vmax=220):
    '''
    Colormap from z coordinate to distinguish between SciFi planes.
    v_min, v_max: color range value in mm.
    value: z plane value to set the color.
    '''
    
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = matplotlib.cm.get_cmap(cmap_name)
    rgb = cmap(norm(abs(value)))[:3]  # rgba[:3] -> rgb
    color = matplotlib.colors.rgb2hex(rgb)
    return color

def extendedFitArr(trackTask, fitHits):
    '''
    Put fit intersection with planes in an array.
    Add by linear extrapolation the missing hits.
    return the full array of 10 ROOT.TVector3().
    '''

    # Fit infos
    fitArr = []
    # trackTask.event.fittedTracks are filled with the fitTrack() function
    for aTrack in trackTask.event.fittedTracks:
        for i in range(aTrack.getNumPointsWithMeasurement()):

            state = aTrack.getFittedState(i)
            pos = state.getPos()
            fitArr.append(pos)

    # Extend the fit display to missed hits points:
    for pos in fitHits:
        fitArr.append(pos) 
    return fitArr

def display3dTrack(arrPosStart, arrPosStop, trackTask, fitHits):
    '''
    arrPosStart/stop: position of the activated fibers, A&B from the getPos() function
    trackTask is the fit object from SndlhcTracking.Tracking()
    fitHits to display individually the missed hits of the fit on the planes.

    Uses matplotlib to display a 3d plot of the track, the fit and missing hits.
    '''
    print(fitHits)
    fig= plt.figure(figsize = (7.5, 6),dpi=500, tight_layout=True)
    ax = plt.axes(projection="3d")

    fitArr = extendedFitArr(trackTask=trackTask, fitHits=fitHits)
    
    # cm to mm conversion *10:
    arrPosStart = [10 * element for element in arrPosStart]
    arrPosStop = [10 * element for element in arrPosStop]
    fitArr = [10 * element for element in fitArr]
    fitHits = [10 * element for element in fitHits]

    for hitNumber in range(len(arrPosStart)):
        ax.plot(
            xs = [arrPosStart[hitNumber][0], arrPosStop[hitNumber][0]], 
            ys = [arrPosStart[hitNumber][1], arrPosStop[hitNumber][1]],
            zs = [arrPosStart[hitNumber][2], arrPosStop[hitNumber][2]],
            ls = '-',
            # RGB format to color different Scifi planes
            color = valueToColor(abs(arrPosStart[hitNumber][2])) )

    ax.plot(
        xs = [element[0] for element in fitArr], 
        ys = [element[1] for element in fitArr],
        zs = [element[2] for element in fitArr],
        color = 'k',
        label = 'Fit')

    ax.scatter3D(
        xs = [hit[0] for hit in fitHits], 
        ys = [hit[1] for hit in fitHits],
        zs = [hit[2] for hit in fitHits],
        color = 'b',
        marker = '^',
        label = 'Missed hits')

    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('z [mm]')
    
    # Legend fields before adding custom entery.
    handles, labels = ax.get_legend_handles_labels()
    #handler_map = matplotlib.legend.Legend.get_legend_handler_map()
    handler_map = matplotlib.legend.Legend.get_default_handler_map()
    # Define custom legend for cluster hits.
    cmap = plt.cm.nipy_spectral_r
    cmap_handles = [Rectangle((0, 0), 1, 1)]
    handler_rainbow = dict(zip(cmap_handles, [HandlerColormap(cmap)]))
    label_rainbow = ['Channel clusters']

    legend1 = plt.legend(
        loc = 'upper right',
        handles = handles, 
        labels = labels, 
        handler_map = handler_map)
    legend2 = plt.legend(
        loc = 'upper left',
        handles = cmap_handles, 
        labels = label_rainbow, 
        handler_map = handler_rainbow)
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)
    plt.savefig('figures/3d.png')
    #plt.show()
    plt.close()

def display2dTrack(arrPosStart, arrPosStop, trackTask, fitHits):
    '''
    x-z and y-z 2d projections of the 3d track.
    arrPosStart, arrPosStop: Array each containing one end of the event clusters.
    trackTask: SndlhcTracking.Tracking() object containing the fit infos.
    fitHits: Array of TVector3 of hits to display.
    '''
    verStart = []
    verStop = []
    horStart = []
    horStop = []
    
    epsilon = 0.00001
    for i in range(len(arrPosStart)):
        delta = arrPosStart[i] - arrPosStop[i]
        if delta[0] < epsilon: #/!\ Change the condition if geometry not aligned anymore.
            verStart.append(arrPosStart[i])
            verStop.append(arrPosStop[i])
        elif delta[1] < epsilon:
            horStart.append(arrPosStart[i])
            horStop.append(arrPosStop[i])
        else:
            print('ERROR: fiber neither vertical nor horizontal!')
            print('Check geometry alignment or change hor/ver conditions.')


    fitArr = extendedFitArr(trackTask = trackTask, fitHits = fitHits)

    fig, (ax1, ax2) = plt.subplots(
        2,
        figsize = (5,7),
        dpi = 500,
        tight_layout = True)

    ax1.grid(axis = 'y')
    ax2.grid(axis = 'y')

    # cm to mm *10 conversion:
    horStart = [10 * element for element in horStart]
    horStop = [10 * element for element in horStop]
    verStart = [10 * element for element in verStart]
    fitArr = [10 * element for element in fitArr]
    verStop = [10 * element for element in verStop]    

    # z-x plane:
    # Horizontal fibers are lines in this plane
    for i in range(len(horStart)):
        ax1.vlines(
            x= horStart[i][2],
            ymin = min(horStart[i][0],horStop[i][0]),
            ymax = max(horStart[i][0],horStop[i][0]),
            colors = 'b')

    # Vertical lines are only points in this plane
    ax1.scatter(
        x=[point[2] for point in verStart],
        y=[point[0] for point in verStart],
        color = 'b',
        marker = '.',
        label = 'Clusters')
    # Add fit:
    # Can use sort() only because it is a straight line
    ax1.plot(
        [vect[2] for vect in fitArr],
        [vect[0] for vect in fitArr],
        color = 'r',
        label = 'Fit')
    ax1.set_xlabel('z [mm]')
    ax1.set_ylabel('x [mm]')
   
    ax1.legend()
    # y-z plane:
    # Vertical fibers are lines in this plane
    for i in range(len(verStart)):
        ax2.vlines(
            x= verStart[i][2],
            ymin = min(verStart[i][1],verStop[i][1]),
            ymax = max(verStart[i][1],verStop[i][1]),
            colors = 'b')
    # Horizontal lines are only points in this plane
    ax2.scatter(
        x=[point[2] for point in horStart],
        y=[point[1] for point in horStart],
        color = 'b',
        marker = '.')
    ax2.plot(
        [vect[2] for vect in fitArr],
        [vect[1] for vect in fitArr],
        color = 'r')
    ax2.set_xlabel('z [mm]')
    ax2.set_ylabel('y [mm]')

    plt.savefig('figures/2d.png')
    #plt.show()
    plt.close()

def chi2Hist(chi2_nDfArr, run, stationNum=0):
    '''
    Chi2/nDOF histogram.
    chi2_nDfArr: Array filled with the value for each event.
    '''
    binsArr = np.linspace(0,200,400)
    fig, ax = plt.subplots(figsize=(8,6), dpi=500)
    fig.suptitle(f'Chi2/Ndof histogram - Test Station {stationNum} - Run {run}')
    ax.hist(chi2_nDfArr, bins=binsArr)
    ax.set_xlim(left=0.0,right=200)
    plt.xlabel('chi2/dof')
    plt.ylabel('Number of events')
    plt.savefig(f'figures/Run{run}_chi2Hist_Stat{stationNum}.png')
    #plt.show()
    plt.close()

def planesHist(nPlanesHit, run): #modified
    '''
    Historam of number of planes hit.
    nPlanesHit: Array filled with number of planes hit (0-10) for each event.
    '''
    fig, ax = plt.subplots(figsize=(5,5), dpi=500)
    fig.suptitle(f'Number of planes hit histogram - Run {run}')
    binsArr = np.linspace(2,11,10)
    ax.hist(nPlanesHit,bins=binsArr)
    plt.xlabel('Number of planes hit')
    plt.ylabel('Number of events')
    plt.savefig(f'figures/Run{run}_nPlanesHist.png')
    #plt.show()
    plt.close()

def gauss(x, A, x0, sigma):
    '''
    Gaussian function f(x).
    Magnitude A, offset x0 and width sigma.
    '''
    return  A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def diffHist(horDiffArr, verDiffArr, stationNum, run): #modified
    '''
    Histogram of position difference between the fit on 4 planes and
    the hit on th 5th one. Gaussian fit.
    horDiffArr, verDiffArr: vertical/horizontal difference between hit and fit,
    each event gives one element of the array.
    stationNum: The test station where the difference hit-fit is measured.
    Return [x0_x, err_x0_x, sigma_x, err_sigma_x, x0_y, err_x0_y, sigma_y, err_sigma_y]
    of gaussian fit for the full-stations plot.
    '''

    # conversion cm in mm *10:
    horDiffArr = [10 * element for element in horDiffArr]
    verDiffArr = [10 * element for element in verDiffArr]

    fig, (ax1, ax2) = plt.subplots(2)
    ax1.grid(axis = 'y')
    ax2.grid(axis = 'y')
    fig.set_size_inches(8, 8)
    fig.suptitle(f'Difference between hit and fit - Test Station {stationNum} - Run {run}',fontsize='x-large',fontweight='bold')
    binsArr = np.linspace(-2.5,2.5,500)
    hist1n, hist1Bins, hist1Patches = ax1.hist(verDiffArr, bins=binsArr)
    ax1.set_xlabel(r'$\Delta$ x [mm]')
    ax1.set_ylabel('Number of events')

    hist2n, hist2Bins, hist2Patches = ax2.hist(horDiffArr, bins=binsArr)
    ax2.set_xlabel(r'$\Delta$ y [mm]')
    ax2.set_ylabel('Number of events')

    diff = hist1Bins[1] - hist1Bins[0]
    binCenters = [hist1Bins[i] + diff for i in range(len(hist1Bins)-1)]
    
    # Gaussian fits:
    param1, cov1 = curve_fit(f=gauss, xdata=binCenters, ydata = hist1n)
    errA1=np.sqrt(cov1[0][0])
    errX01=np.sqrt(cov1[1][1])
    errSigma1=np.sqrt(cov1[2][2])
    ax1.plot(
        binCenters, 
        [gauss(
            x = binCenters[i],
            A = param1[0],
            x0 = param1[1],
            sigma = param1[2]) 
            for i in range(len(binCenters))],
        color = 'r',
        label = (f'Gaussian fit: x0 = ({param1[1]:.2} ± {errX01:.2}) mm'
               + f'\n                     sigma = ({abs(param1[2]):.2} ± {errSigma1:.2}) mm'))
    ax1.legend()

    param2, cov2 = curve_fit(f=gauss, xdata=binCenters, ydata = hist2n)
    errA2 = np.sqrt(cov2[0][0])
    errX02 = np.sqrt(cov2[1][1])
    errSigma2 = np.sqrt(cov2[2][2])
    ax2.plot(
        binCenters, 
        [gauss(
            x = binCenters[i],
            A = param2[0],
            x0 = param2[1],
            sigma = param2[2]) 
            for i in range(len(binCenters))],
        color = 'r',
        label = (f'Gaussian fit: x0 = ({param2[1]:.2} ± {errX02:.2}) mm'
               + f'\n                     sigma = ({abs(param2[2]):.2} ± {errSigma2:.2}) mm'))
    ax2.legend()
    plt.savefig(f'figures/Run{run}_diffHistGauss_testStation{stationNum}.png')
    #plt.show()

    plt.close()
    # Sigma can be fitted negative, but return abs() by convention.
    resultFit = [param1[1], errX01, abs(param1[2]), errSigma1, param2[1], errX02, abs(param2[2]), errSigma2]
    return resultFit

def allPlanesGauss(fitArr, run): #modified
    '''
    Global plot for the 5 test stations.
    fitArr: array of (return format from diffHist()) for each station
    '''

    # fitArr already in mm
    stationsArr = [1,2,3,4,5]
    zeroArr = [0,0,0,0,0]
    # Extract data from array: See order of return of diffHist()
    x0_x = [x[0] for x in fitArr]
    err_x0_x = [x[1] for x in fitArr]
    sigma_x = [x[2] for x in fitArr]
    err_sigma_x = [x[3] for x in fitArr]
    x0_y = [x[4] for x in fitArr]
    err_x0_y = [x[5] for x in fitArr]
    sigma_y = [x[6] for x in fitArr]
    err_sigma_y = [x[7] for x in fitArr]

    fig, ((ax1, ax3),(ax2, ax4)) = plt.subplots(
        nrows = 2,
        ncols = 2,
        sharex = 'col')
    fig.set_size_inches(8, 8)
    fig.suptitle(f'Fit result for all stations - Run {run}')
    ax1.errorbar(
        x = stationsArr,
        y = x0_x,
        xerr = zeroArr,
        yerr = err_x0_x,
        ls = '',
        marker = 'x',
        markeredgecolor = 'k')
    ax2.errorbar(
        x = stationsArr,
        y = sigma_x,
        xerr = zeroArr,
        yerr = err_sigma_x,
        ls = '',
        marker = 'x',
        markeredgecolor = 'k')
    ax3.errorbar(
        x = stationsArr,
        y = x0_y,
        xerr = zeroArr,
        yerr = err_x0_y,
        ls = '',
        marker = 'x',
        markeredgecolor = 'k')
    ax4.errorbar(
        x = stationsArr,
        y = sigma_y,
        xerr = zeroArr,
        yerr = err_sigma_y,
        ls = '',
        marker = 'x',
        markeredgecolor = 'k')
    ax1.set_ylabel('X offset [mm]')
    ax2.set_ylabel(r'$\sigma_x$ [mm]')
    ax3.set_ylabel('Y offset [mm]')
    ax4.set_ylabel(r'$\sigma_y$ [mm]')
    ax2.set_xlabel('Test station')
    ax4.set_xlabel('Test station')
    ax1.set_ylim(bottom = -3, top = 3)
    ax3.set_ylim(bottom = -3, top = 3)
    ax2.set_ylim(bottom = 0.05, top = 0.35)
    ax4.set_ylim(bottom = 0.05, top = 0.35)
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()
    plt.savefig(f'figures/Run{run}_FullStationsDiff.png')
    #plt.show()
    plt.close()

def diffPosHist(posArr, diffArr, binsPos, labels, run, fileName, title, isCrossed):
    '''
    Difference fit-hit versus position on plane histogram.
    diffArr, posArr, binsPos must be given in cm, conversion in mm inside function.

    posArr: cluster position, one element of the array for each event.
    diffArr: difference fit-hit, one element of the array for each event.
    binsPos: histogram bins for the cluster position. Must set the limits.
    according to the plane geometry.
    labels: [xlabel, ylabel] for the axes.
    fileName: file name without extension, to save the plot.
    isCrossed: True if X-Y axis are mixed in the same plot to check rotations.
    Return: slope and its uncertainty for rotation global plot, else None.
    '''

    fig, ax = plt.subplots(figsize=(8,6),dpi=500)
    fig.suptitle(title)

    plt.rcParams.update({'font.size': 10})

    # cm to mm conversion *10:
    diffArr = [10 * element for element in diffArr]
    posArr = [10 * element for element in posArr]
    binsPos = [10 * element for element in binsPos]

    # Histogram limits:
    yMin = -2.5
    yMax = 2.5
    xMin = -400
    xMax = 400
    nBins = [60,30]
    if isCrossed: # Difference in y axis for better horizontal fit.
        binx = binsPos
        biny = np.linspace(-2.5,2.5,50)
        xHist = posArr
        yHist = diffArr
    else:
        binx = np.linspace(-2.5,2.5,50)
        biny = binsPos
        xHist = diffArr
        yHist = posArr
    
    plt.hist2d(
        xHist,
        yHist,
        bins = [binx,biny],
        cmap = plt.get_cmap('gist_heat_r'))

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Number of events')

    # Linear fit for rotation offset:
    if isCrossed:
        linFitModel, cov = np.polyfit(posArr, diffArr, 1,cov = True) # Compute slope and intercept
        linFitFunc = np.poly1d(linFitModel) # Make the line equation
        #print(linFitModel)
        ax.plot(
            binx,
            linFitFunc(binx),
            color = 'b',
            label = f'Slope: {np.format_float_scientific(linFitModel[0], precision = 2,)}±{np.format_float_scientific(np.sqrt(cov[0][0]), precision = 2,)}')
        plt.legend()
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    
    plt.savefig(f'figures/Run{run}_{fileName}.png')
    #plt.show()
    plt.close()
    if isCrossed:
        return linFitModel[0], np.sqrt(cov[0][0]) # Return slope and its uncertainty.

def rotationAngle(xPosyOff_Slope, xPosyOff_SlopeErr, yPosxOff_Slope, yPosxOff_SlopeErr, mats, run): #modified
    '''
    Global rotation plot for each stations, 
    Use what is returned by diffPosHist().
    xPosyOff: Cross the X position and the Y offset
    yPosxOff: Cross the Y position and the X offset
    mats : True if rotation angle measured for each mat -> separation into the three mats for each station
    '''
    x1 = [0.95, 1.95, 2.95, 3.95, 4.95]
    x11 = [0.97, 1.97, 2.97, 3.976, 4.97]
    x12 = [0.99, 1.99, 2.99, 3.99, 4.99]
    x2 = [1.05, 2.05, 3.05, 4.05, 5.05]
    x21 = [1.07, 2.07, 3.07, 4.07, 5.07]
    x22 = [1.09, 2.09, 3.09, 4.09, 5.09]
    zeroes = [0, 0, 0, 0, 0]

    # Both plane slopes must be the same in absolute value.
    xPosyOff_Slope = [abs(element) for element in xPosyOff_Slope]
    yPosxOff_Slope = [abs(element) for element in yPosxOff_Slope]

    if mats :
        xPosyOff_Slope1, xPosyOff_SlopeErr1, xPosyOff_Slope2, xPosyOff_SlopeErr2, xPosyOff_Slope3, xPosyOff_SlopeErr3 = [],[],[],[],[],[]
        yPosxOff_Slope1, yPosxOff_SlopeErr1, yPosxOff_Slope2, yPosxOff_SlopeErr2, yPosxOff_Slope3, yPosxOff_SlopeErr3 = [],[],[],[],[],[]
        for i in range(15) : # 15 == 3 mats x 5 stations
            if i%3 == 0 :
                xPosyOff_Slope1.append(xPosyOff_Slope[i])
                xPosyOff_SlopeErr1.append(xPosyOff_SlopeErr[i])
                yPosxOff_Slope1.append(yPosxOff_Slope[i])
                yPosxOff_SlopeErr1.append(yPosxOff_SlopeErr[i])
            if i%3 == 1 :
                xPosyOff_Slope2.append(xPosyOff_Slope[i])
                xPosyOff_SlopeErr2.append(xPosyOff_SlopeErr[i])
                yPosxOff_Slope2.append(yPosxOff_Slope[i])
                yPosxOff_SlopeErr2.append(yPosxOff_SlopeErr[i])
            if i%3 == 2 :
                xPosyOff_Slope3.append(xPosyOff_Slope[i])
                xPosyOff_SlopeErr3.append(xPosyOff_SlopeErr[i])
                yPosxOff_Slope3.append(yPosxOff_Slope[i])
                yPosxOff_SlopeErr3.append(yPosxOff_SlopeErr[i])

        fig, ax = plt.subplots(figsize=(8,6),dpi=500)
        fig.suptitle(f'Rotation angles between mats - Run {run}', fontsize=16)
        ax.errorbar(x = x1, y = xPosyOff_Slope1, xerr = zeroes, yerr = xPosyOff_SlopeErr1, ls = '', ecolor = 'b', marker = 'x', markeredgecolor = 'b', label = 'X position, Y offset-mat1')
        ax.errorbar(x = x11, y = xPosyOff_Slope2, xerr = zeroes, yerr = xPosyOff_SlopeErr2, ls = '', ecolor = 'royalblue', marker = 'x', markeredgecolor = 'royalblue', label = 'X position, Y offset-mat2')
        ax.errorbar(x = x12, y = xPosyOff_Slope3, xerr = zeroes, yerr = xPosyOff_SlopeErr3, ls = '', ecolor = 'navy', marker = 'x', markeredgecolor = 'navy', label = 'X position, Y offset-mat3')
        
        ax.errorbar(x = x2, y = yPosxOff_Slope1, xerr = zeroes, yerr = yPosxOff_SlopeErr1, ls = '', ecolor = 'r', marker = 'x', markeredgecolor = 'r', label = 'Y position, X offset-mat1')
        ax.errorbar(x = x21, y = yPosxOff_Slope2, xerr = zeroes, yerr = yPosxOff_SlopeErr2, ls = '', ecolor = 'coral', marker = 'x', markeredgecolor = 'coral', label = 'Y position, X offset-mat2')
        ax.errorbar(x = x22, y = yPosxOff_Slope3, xerr = zeroes, yerr = yPosxOff_SlopeErr3, ls = '', ecolor = 'orange', marker = 'x', markeredgecolor = 'orange', label = 'Y position, X offset-mat3')

        ax.grid(axis = 'y')
        plt.tick_params(bottom=False)
        plt.xlabel('Station number')
        plt.ylabel('Relative rotation angle [rad]')
        plt.legend()
        plt.savefig(f'figures/Run{run}_angles_mats.png')
        plt.show()
        plt.close()

    else : 
        fig, ax = plt.subplots(figsize=(8,6),dpi=500)
        fig.suptitle(f'Rotation angles between planes - Run {run}', fontsize=16)
        ax.errorbar(x = x1, y = xPosyOff_Slope, xerr = zeroes, yerr = xPosyOff_SlopeErr, ls = '', marker = 'x', markeredgecolor = 'b', label = 'X position, Y offset')
        ax.errorbar( x = x2, y = yPosxOff_Slope, xerr = zeroes, yerr = yPosxOff_SlopeErr, ls = '', marker = 'x', markeredgecolor = 'r', label = 'Y position, X offset')

        ax.grid(axis = 'y')
        plt.tick_params(bottom=False)
        plt.xlabel('Station number')
        plt.ylabel('Relative rotation angle [rad]')
        plt.legend()
        plt.savefig(f'figures/Run{run}_angles.png')
        plt.close()

def hitMap(hitChanX, hitChanY, option, run, NbEvents) : #added
    '''
    Histogram of hit channels
    title as parameter to be changed depending on run number
    hitChanX : array(array(arrays(channel indices for each station) for each event)), channel index = 0-1535
    ex: hitChanX = [[[chan indices of x hits on stat1 for event1], [chan indices of x hits on stat2 for event1]], etc.]
    option : string, ex : 'raw', 'eventselection' ...
    run : run number
    '''  
    hitChanX1,hitChanX2,hitChanX3,hitChanX4,hitChanX5 = [],[],[],[],[]
    hitChanXtot = [hitChanX1,hitChanX2,hitChanX3,hitChanX4,hitChanX5]
    for i in range(len(hitChanX)) : #len(hitChanX) == nb of events
        for j in range(len(hitChanX[i])) : #len(hitChanX[i]) == nb of station hit in event i
            for k in range(len(hitChanX[i][j])) : #len(hitChanX[i][j]) == nb of hits in station j in event i
                hitChanXtot[j].append(hitChanX[i][j][k])
       
    hitChanY1,hitChanY2,hitChanY3,hitChanY4,hitChanY5 = [],[],[],[],[]
    hitChanYtot = [hitChanY1,hitChanY2,hitChanY3,hitChanY4,hitChanY5]
    for i in range(len(hitChanY)) : #len(hitChanY) == nb of events
        for j in range(len(hitChanY[i])) :
            for k in range(len(hitChanY[i][j])) :
                hitChanYtot[j].append(hitChanY[i][j][k])
        
    for j in range(5) : # loop on stations 1->5
        fig, (ax1, ax2) = plt.subplots(2)
        fig.set_size_inches(8, 8)
        fig.suptitle(f'Hit map - Station {j+1} - Run {run} - {option} : Nb of events = {NbEvents}', fontsize=16)
        binsArr = np.linspace(0,1535,1536)
        ax1.hist(hitChanXtot[j],bins=binsArr)
        ax1.set_xlabel('Vertical channel hits', fontsize=14)
        ax1.set_ylabel('Number of events', fontsize=14)
        ax2.hist(hitChanYtot[j],bins=binsArr)
        ax2.set_xlabel('Horizontal channel hits', fontsize=14)
        ax2.set_ylabel('Number of events', fontsize=14)
        plt.savefig(f'figures/Run{run}_HitMap_Stat{j+1}_{option}.png')
        #plt.show()
        plt.close()
       
        fig, (ax1, ax2) = plt.subplots(2)
        fig.set_size_inches(8, 8)
        fig.suptitle(f'Hit map Zoom - Station {j+1} - Run {run} - {option} : Nb of events = {NbEvents}', fontsize=16)
        binsArr = np.linspace(0,600,601)
        ax1.hist(hitChanXtot[j],bins=binsArr)
        ax1.set_xlabel('Vertical channel hits', fontsize=14)
        ax1.set_ylabel('Number of events', fontsize=14)
        ax2.hist(hitChanYtot[j],bins=binsArr)
        ax2.set_xlabel('Horizontal channel hits', fontsize=14)
        ax2.set_ylabel('Number of events', fontsize=14)
        plt.savefig(f'figures/Run{run}_HitMap_Stat{j+1}-Zoom_{option}.png')
        #plt.show()
        plt.close()

def CluSizeDistrib(Clusters, option, run, NbEvents) : #added
    '''
    Histogram of cluster size distribution
    Clusters = array of ClusterArr, ClusterArr = array of cluster for each event, cluster = array of hit chan
    -> Clusters == [[[clu1 evt1],[clu2 evt1],[clu3 evt1]],[[clu1 evt2],[clu2 evt2]], ... , [[clu1 evtN], [clu2 evtN]]]
    '''
    CluSize1, CluSize2, CluSize3, CluSize4, CluSize5 = [], [], [], [], []
    CluSize = [CluSize1, CluSize2, CluSize3, CluSize4, CluSize5]
    Stations = [1, 2, 3, 4, 5]
    Boards = [0,1,2]
    cluB11h, cluB12h, cluB13h, cluB21h, cluB22h, cluB23h, cluB31h, cluB32h, cluB33h, cluB41h, cluB42h, cluB43h, cluB51h, cluB52h, cluB53h = [], [], [], [], [],[], [], [], [], [], [], [], [], [], []
    cluB11v, cluB12v, cluB13v, cluB21v, cluB22v, cluB23v, cluB31v, cluB32v, cluB33v, cluB41v, cluB42v, cluB43v, cluB51v, cluB52v, cluB53v = [], [], [], [], [],[], [], [], [], [], [], [], [], [], []
    CluS1h = [cluB11h, cluB12h, cluB13h]
    CluS2h = [cluB21h, cluB22h, cluB23h]
    CluS3h = [cluB31h, cluB32h, cluB33h]
    CluS4h = [cluB41h, cluB42h, cluB43h]
    CluS5h = [cluB51h, cluB52h, cluB53h]
    CluS1v = [cluB11v, cluB12v, cluB13v]
    CluS2v = [cluB21v, cluB22v, cluB23v]
    CluS3v = [cluB31v, cluB32v, cluB33v]
    CluS4v = [cluB41v, cluB42v, cluB43v]
    CluS5v = [cluB51v, cluB52v, cluB53v]
    CluS1 = [CluS1h, CluS1v]
    CluS2 = [CluS2h, CluS2v]
    CluS3 = [CluS3h, CluS3v]
    CluS4 = [CluS4h, CluS4v]
    CluS5 = [CluS5h, CluS5v]
    CluSizeB = [CluS1, CluS2, CluS3, CluS4, CluS5]
    for cluArr in Clusters :
        for clu in cluArr :
            for stat in Stations :
                if (clu[0]//1000000) == stat :
                    CluSize[stat-1].append(len(clu))
                    for board in Boards :
                        if ((clu[0]//10000) - (clu[0]//100000)*10) == board and (clu[0]//100000)%2 == 0:
                            CluSizeB[stat-1][0][board].append(len(clu))
                        if ((clu[0]//10000) - (clu[0]//100000)*10) == board and (clu[0]//100000)%2 == 1:
                            CluSizeB[stat-1][1][board].append(len(clu))
    
    # Cluster size histograms for each station
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5)
    fig.set_size_inches(18, 12)
    fig.suptitle(f'Cluster size distribution - Run {run} - {option} : Nb of events = {NbEvents}', fontsize=20)
    ax1 = plt.subplot(231)
    ax2 = plt.subplot(232)
    ax3 = plt.subplot(233)
    ax4 = plt.subplot(234)
    ax5 = plt.subplot(235)
    binsArr = np.linspace(0,10,11)
    ax1.hist(CluSize[0],bins=binsArr)
    ax1.set_xlabel('Number of channel per cluster, station 1', fontsize=16)
    ax1.set_ylabel('Number of events', fontsize=16)
    ax2.hist(CluSize[1],bins=binsArr)
    ax2.set_xlabel('Number of channel per cluster, station 2', fontsize=16)
    ax2.set_ylabel('Number of events', fontsize=16)
    ax3.hist(CluSize[2],bins=binsArr)
    ax3.set_xlabel('Number of channel per cluster, station 3', fontsize=16)
    ax3.set_ylabel('Number of events', fontsize=16)
    ax4.hist(CluSize[3],bins=binsArr)
    ax4.set_xlabel('Number of channel per cluster, station 4', fontsize=16)
    ax4.set_ylabel('Number of events', fontsize=16)
    ax5.hist(CluSize[4],bins=binsArr)
    ax5.set_xlabel('Number of channel per cluster, station 5', fontsize=16)
    ax5.set_ylabel('Number of events', fontsize=16)
    plt.savefig(f'figures/Run{run}_CluSize_{option}.png')
    #plt.show()
    plt.close()

    # Cluster size histograms for each board 
    # vertical clusters
    for s in range(1,6) : #loop on station numbers
        fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, constrained_layout=True)
        fig.set_size_inches(18, 5)
        fig.suptitle(f'Cluster size distribution per vertical board, station {s} - Run {run} - {option}', fontsize=20)
        binsArr = np.linspace(0,10,11)
        ax1.hist(CluSizeB[s-1][1][0],bins=binsArr)
        ax1.set_xlabel('Number of channel per cluster, board 1', fontsize=16)
        ax1.set_ylabel('Number of events', fontsize=16)
        ax2.hist(CluSizeB[s-1][1][1],bins=binsArr)
        ax2.set_xlabel('Number of channel per cluster, board 2', fontsize=16)
        ax2.set_ylabel('Number of events', fontsize=16)
        ax3.hist(CluSizeB[s-1][1][2],bins=binsArr)
        ax3.set_xlabel('Number of channel per cluster, board 3', fontsize=16)
        ax3.set_ylabel('Number of events', fontsize=16)
        plt.savefig(f'figures/Run{run}_CluSizeBoard-Stat{s}-Ver_{option}.png')
        #plt.show()
        plt.close()

        # horizontal clusters
    for s in range(1,6) : #loop on station numbers
        fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, constrained_layout=True)
        fig.set_size_inches(18, 5)
        fig.suptitle(f'Cluster size distribution per horizontal board, station {s} - Run {run} - {option}', fontsize=20)
        binsArr = np.linspace(0,10,11)
        ax1.hist(CluSizeB[s-1][0][0],bins=binsArr)
        ax1.set_xlabel('Number of channel per cluster, board 1', fontsize=16)
        ax1.set_ylabel('Number of events', fontsize=16)
        ax2.hist(CluSizeB[s-1][0][1],bins=binsArr)
        ax2.set_xlabel('Number of channel per cluster, board 2', fontsize=16)
        ax2.set_ylabel('Number of events', fontsize=16)
        ax3.hist(CluSizeB[s-1][0][2],bins=binsArr)
        ax3.set_xlabel('Number of channel per cluster, board 3', fontsize=16)
        ax3.set_ylabel('Number of events', fontsize=16)
        plt.savefig(f'figures/Run{run}_CluSizeBoard-Stat{s}-Hor_{option}.png')
        #plt.show()
        plt.close()

def Efficiency(df_Events, testStation, run) : #added
    '''
    df_Events : DataFrame with events informations on track and cluster match
    Plots : 2D hist of tracks position
            2D hist of efficiency as function of position with hor matches, ver maches, hor+ver matches
            2D hist of same but rescaled on efficiencies between 0.8 and 1
    '''
    # Dataframes with conditions on matches
    #events with hor and ver matches
    Matches = df_Events[(df_Events['Clu_match_hor'] == 1) & (df_Events['Clu_match_ver'] == 1) & (df_Events['Residual_hor'] < 1) & (df_Events['Residual_ver'] < 1)] 
    # events with hor matches
    horMatches = df_Events[(df_Events['Clu_match_hor'] == 1) & (df_Events['Residual_hor'] < 1)] 
    # events with ver matches
    verMatches = df_Events[(df_Events['Clu_match_ver'] == 1) & (df_Events['Residual_ver'] < 1)] 

    # Efficiency vs position
    nbins = 100
    xbins = np.linspace(-500, -50, nbins)
    ybins = np.linspace(130, 570, nbins)

    h_ref, edges_x_ref, edges_y_ref = np.histogram2d(df_Events.FitPos_ver*10, df_Events.FitPos_hor*10, bins=[xbins, ybins])
    h_hor, edges_x_hor, edges_y_hor = np.histogram2d(horMatches.FitPos_ver*10, horMatches.FitPos_hor*10, bins=[xbins, ybins])
    h_ver, edges_x_ver, edges_y_ver = np.histogram2d(verMatches.FitPos_ver*10, verMatches.FitPos_hor*10, bins=[xbins, ybins])
    h_both, edges_x_both, edges_y_both = np.histogram2d(Matches.FitPos_ver*10, Matches.FitPos_hor*10, bins=[xbins, ybins])

    # Transpose histograms since histogram doesn't follow cartesian convention
    h_ref = h_ref.T
    h_hor = h_hor.T
    h_ver = h_ver.T
    h_both = h_both.T

    # Mask values = 0 in tracks hist
    h_ref = np.ma.masked_equal(h_ref, 0)
    mask = np.ma.getmaskarray(h_ref)
    h_hor = np.ma.masked_array(h_hor, mask)
    h_ver = np.ma.masked_array(h_ver, mask)
    h_both = np.ma.masked_array(h_both, mask)

    # Plot tracks histogram
    fig, ax1  = plt.subplots()
    fig.set_size_inches(8,8)
    fig.suptitle(f'Reconstructed track histogram - Test Station {testStation} - Run {run}')
    plot1 = ax1.imshow(h_ref, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_ref[0], edges_x_ref[-1], edges_y_ref[0], edges_y_ref[-1]])
    ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=14)
    ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=14)
    cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Number of events', size=14)
    plt.savefig(f'figures/Run{run}_TrackHist_Stat{testStation}_{nbins}bins.png')
    #plt.show()
    plt.close()

    # Efficiencies
    Effhor = h_hor/h_ref
    Effver = h_ver/h_ref
    Effboth = h_both/h_ref

    # Plot Efficiency vs x-y position
    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, constrained_layout=True)
    fig.set_size_inches(18, 5)
    fig.suptitle(f'Efficiency - Test Station {testStation} - Run {run}')
    plot1 = ax1.imshow(Effhor, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edges_x_hor[0], edges_x_hor[-1], edges_y_hor[0], edges_y_hor[-1]]) 
    plot2 = ax2.imshow(Effver, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edges_x_ver[0], edges_x_ver[-1], edges_y_ver[0], edges_y_ver[-1]])
    plot3 = ax3.imshow(Effboth, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edges_x_both[0], edges_x_both[-1], edges_y_both[0], edges_y_both[-1]])
    ax1.title.set_text('Horizontal match')
    ax2.title.set_text('Vertical match')
    ax3.title.set_text('Horizontal and vertical match')
    ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Efficiency', size=12)
    cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='Efficiency', size=12)
    cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='Efficiency', size=12)
    plt.savefig(f'figures/Run{run}_Efficiency-XYcorrel_Stat{testStation}_{nbins}bins.png')
    #plt.show()
    plt.close()

    # Zoom on Eff between 0.8 and 1
    Effhor = np.ma.masked_less(Effhor, 0.8)
    Effver = np.ma.masked_less(Effver, 0.8)
    Effboth = np.ma.masked_less(Effboth, 0.8)

    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, constrained_layout=True)
    fig.set_size_inches(18, 5)
    fig.suptitle(f'Efficiency : 0.8 to 1 - Test Station {testStation} - Run {run}')
    plot1 = ax1.imshow(Effhor, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edges_x_hor[0], edges_x_hor[-1], edges_y_hor[0], edges_y_hor[-1]]) 
    plot2 = ax2.imshow(Effver, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edges_x_ver[0], edges_x_ver[-1], edges_y_ver[0], edges_y_ver[-1]])
    plot3 = ax3.imshow(Effboth, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edges_x_both[0], edges_x_both[-1], edges_y_both[0], edges_y_both[-1]])
    ax1.title.set_text('Horizontal match')
    ax2.title.set_text('Vertical match')
    ax3.title.set_text('Horizontal and vertical match')
    ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Efficiency', size=12)
    cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='Efficiency', size=12)
    cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='Efficiency', size=12)
    plt.savefig(f'figures/Run{run}_Efficiency-XYcorrel_Zoom_Stat{testStation}_{nbins}bins.png')
    #plt.show()
    plt.close()

def slopeEfficiency(df_Events, testStation, run) : #added
    '''
    df_Events : DataFrame with events informations on track and cluster match
    Plots : 2D hist of reconstructed tracks with cut on slopes
            2D hist of efficiency as fct of pos with hor matches, ver maches, hor+ver matches with cut on slopes
            
    '''
    # Cuts definition : determined thanks to slope plots in trackHistograms
    cut0x = -1.5
    cutx = -1.3
    cut0y = -0.1
    cuty = 0.1
    # Cut on slopes
    slopex = df_Events.FitSlope_x * 180/np.pi
    slopey = df_Events.FitSlope_y * 180/np.pi
    df_Events_parallel = df_Events[(slopex >= cut0x) & (slopex < cutx) & (slopey >= cut0y) & (slopey < cuty)]
    #df_Events_slope = df_Events[(slopex < cut0x) | (slopex >= cutx) | (slopey < cut0y) | (slopey >= cuty)]
    df_Events_slope = df_Events[((slopex < cut0x) | (slopex >= cutx)) & ((slopey < cut0y) | (slopey >= cuty))]
    # dataframes with matches
    Matches_par = df_Events_parallel[(df_Events_parallel['Clu_match_hor'] == 1) & (df_Events_parallel['Clu_match_ver'] == 1) & (df_Events_parallel['Residual_hor'] < 1) & (df_Events_parallel['Residual_ver'] < 1)] #events with hor and ver matches
    horMatches_par = df_Events_parallel[(df_Events_parallel['Clu_match_hor'] == 1) & (df_Events_parallel['Residual_hor'] < 1)] # events with hor matches
    verMatches_par = df_Events_parallel[(df_Events_parallel['Clu_match_ver'] == 1) & (df_Events_parallel['Residual_ver'] < 1)] # events with ver matches
    Matches_sl = df_Events_slope[(df_Events_slope['Clu_match_hor'] == 1) & (df_Events_slope['Clu_match_ver'] == 1) & (df_Events_slope['Residual_hor'] < 1) & (df_Events_slope['Residual_ver'] < 1)] #events with hor and ver matches
    horMatches_sl = df_Events_slope[(df_Events_slope['Clu_match_hor'] == 1) & (df_Events_slope['Residual_hor'] < 1)] # events with hor matches
    verMatches_sl = df_Events_slope[(df_Events_slope['Clu_match_ver'] == 1) & (df_Events_slope['Residual_ver'] < 1)] # events with ver matches

    # Efficiency vs position
    nbins = 100
    xbins = np.linspace(-500, -50, nbins)
    ybins = np.linspace(130, 570, nbins)

    hpar_ref, edgespar_x_ref, edgespar_y_ref = np.histogram2d(df_Events_parallel.FitPos_ver*10, df_Events_parallel.FitPos_hor*10, bins=[xbins, ybins])
    hpar_hor, edgespar_x_hor, edgespar_y_hor = np.histogram2d(horMatches_par.FitPos_ver*10, horMatches_par.FitPos_hor*10, bins=[xbins, ybins])
    hpar_ver, edgespar_x_ver, edgespar_y_ver = np.histogram2d(verMatches_par.FitPos_ver*10, verMatches_par.FitPos_hor*10, bins=[xbins, ybins])
    hpar_both, edgespar_x_both, edgespar_y_both = np.histogram2d(Matches_par.FitPos_ver*10, Matches_par.FitPos_hor*10, bins=[xbins, ybins])
    hsl_ref, edgessl_x_ref, edgessl_y_ref = np.histogram2d(df_Events_slope.FitPos_ver*10, df_Events_slope.FitPos_hor*10, bins=[xbins, ybins])
    hsl_hor, edgessl_x_hor, edgessl_y_hor = np.histogram2d(horMatches_sl.FitPos_ver*10, horMatches_sl.FitPos_hor*10, bins=[xbins, ybins])
    hsl_ver, edgessl_x_ver, edgessl_y_ver = np.histogram2d(verMatches_sl.FitPos_ver*10, verMatches_sl.FitPos_hor*10, bins=[xbins, ybins])
    hsl_both, edgessl_x_both, edgessl_y_both = np.histogram2d(Matches_sl.FitPos_ver*10, Matches_sl.FitPos_hor*10, bins=[xbins, ybins])

    # Transpose histograms since histogram doesn't follow cartesian convention
    hpar_ref = hpar_ref.T
    hpar_hor = hpar_hor.T
    hpar_ver = hpar_ver.T
    hpar_both = hpar_both.T
    hsl_ref = hsl_ref.T
    hsl_hor = hsl_hor.T
    hsl_ver = hsl_ver.T
    hsl_both = hsl_both.T

    # Mask values = 0 in tracks hist
    hpar_ref = np.ma.masked_equal(hpar_ref, 0)
    maskpar = np.ma.getmaskarray(hpar_ref)
    hpar_hor = np.ma.masked_array(hpar_hor, maskpar)
    hpar_ver = np.ma.masked_array(hpar_ver, maskpar)
    hpar_both = np.ma.masked_array(hpar_both, maskpar)
    hsl_ref = np.ma.masked_equal(hsl_ref, 0)
    masksl = np.ma.getmaskarray(hsl_ref)
    hsl_hor = np.ma.masked_array(hsl_hor, masksl)
    hsl_ver = np.ma.masked_array(hsl_ver, masksl)
    hsl_both = np.ma.masked_array(hsl_both, masksl)

    # Plot tracks histogram
    fig, (ax1,ax2)  = plt.subplots(nrows = 1, ncols = 2, constrained_layout = True)
    fig.set_size_inches(15,8)
    fig.suptitle(f'Reconstructed track histogram with cut on slopes - Test Station {testStation} - Run {run}')
    plot1 = ax1.imshow(hpar_ref, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edgespar_x_ref[0], edgespar_x_ref[-1], edgespar_y_ref[0], edgespar_y_ref[-1]])
    plot2 = ax2.imshow(hsl_ref, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edgessl_x_ref[0], edgessl_x_ref[-1], edgessl_y_ref[0], edgessl_y_ref[-1]])
    ax1.title.set_text(f'{cut0x}° <= x slope < {cutx}° and {cut0y}° <= y slope < {cuty}°')
    ax2.title.set_text(f'x and/or y slope out of ranges')
    ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=14)
    ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=14)
    ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=14)
    ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=14)
    cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Number of events', size=14)
    cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='Number of events', size=14)
    plt.savefig(f'figures/Run{run}_TrackHist-cutslopes-X{cut0x}-{cutx}-Y{cut0y}-{cuty}_Stat{testStation}_{nbins}bins.png')
    #plt.show()
    plt.close()

    # Efficiencies
    Effhor_par = hpar_hor/hpar_ref
    Effver_par = hpar_ver/hpar_ref
    Effboth_par = hpar_both/hpar_ref
    Effhor_sl = hsl_hor/hsl_ref
    Effver_sl = hsl_ver/hsl_ref
    Effboth_sl = hsl_both/hsl_ref

    # Plots
    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, constrained_layout=True)
    fig.set_size_inches(18, 5)
    fig.suptitle(f'Efficiency, {cut0x}° <= x slope < {cutx}° and {cut0y}° <= y slope < {cuty}° - Test Station {testStation} - Run {run}')
    plot1 = ax1.imshow(Effhor_par, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edgespar_x_hor[0], edgespar_x_hor[-1], edgespar_y_hor[0], edgespar_y_hor[-1]]) 
    plot2 = ax2.imshow(Effver_par, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edgespar_x_ver[0], edgespar_x_ver[-1], edgespar_y_ver[0], edgespar_y_ver[-1]])
    plot3 = ax3.imshow(Effboth_par, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edgespar_x_both[0], edgespar_x_both[-1], edgespar_y_both[0], edgespar_y_both[-1]])
    ax1.title.set_text('Horizontal match')
    ax2.title.set_text('Vertical match')
    ax3.title.set_text('Horizontal and vertical match')
    ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Efficiency', size=12)
    cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='Efficiency', size=12)
    cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='Efficiency', size=12)
    plt.savefig(f'figures/Run{run}_Efficiency-slopein-X{cut0x}-{cutx}-Y{cut0y}-{cuty}_Stat{testStation}_{nbins}bins.png')
    #plt.show()
    plt.close()

    fig, (ax1, ax2, ax3) = plt.subplots(nrows = 1, ncols = 3, constrained_layout=True)
    fig.set_size_inches(18, 5)
    fig.suptitle(f'Efficiency, x and/or y slope out of ranges - Test Station {testStation} - Run {run}')
    plot1 = ax1.imshow(Effhor_sl, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edgessl_x_hor[0], edgessl_x_hor[-1], edgessl_y_hor[0], edgessl_y_hor[-1]]) 
    plot2 = ax2.imshow(Effver_sl, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edgessl_x_ver[0], edgessl_x_ver[-1], edgessl_y_ver[0], edgessl_y_ver[-1]])
    plot3 = ax3.imshow(Effboth_sl, cmap= plt.get_cmap('RdYlGn'), origin='lower', extent=[edgessl_x_both[0], edgessl_x_both[-1], edgessl_y_both[0], edgessl_y_both[-1]])
    ax1.title.set_text('Horizontal match')
    ax2.title.set_text('Vertical match')
    ax3.title.set_text('Horizontal and vertical match')
    ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
    ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
    cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Efficiency', size=12)
    cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='Efficiency', size=12)
    cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='Efficiency', size=12)
    plt.savefig(f'figures/Run{run}_Efficiency-slopeout-X{cut0x}-{cutx}-Y{cut0y}-{cuty}_Stat{testStation}_{nbins}bins.png')
    #plt.show()
    plt.close()

def SIPMgap(df_Events, testStation, run, superpose = False) : #added
    '''
    df_Events : DataFrame with events informations on track and cluster match
    superpose : set to True to plot also with all SiPM of one direction superposed 2 by 2
    Plots : 1D histogram of track repartition for all SiPMs of one direction superposed 2 by 2
            Efficiency plot as function of x (y) for all SiPMs of one direction superposed 2 by 2
            1D histogram of track repartition for two SiPMs in beam region
            Efficiency plot as function of x (y) for two SiPMs in beam region
            Efficiency zoomed on gap between two SiPMs in beam region
            Efficiency zoomed on gap between silicon dice in SiPMs in beam region
    '''
    # Histogram of tracks - zoom in beam region
    nbins = 300
    width = 65.3 #width of 2 SiPMs : 256*0.25 (channels) + 2*0.25 (tofpets gaps) + 2*0.4 (SiPMs gaps)
    xmin = -473.1
    xmax = xmin+width
    ymin = 158.6
    ymax = ymin+width
    xbin = np.linspace(xmin, xmax, nbins)
    ybin = np.linspace(ymin, ymax, nbins)
    # Location of 2 SiPMs in beam region
    xbin_best = np.linspace(xmin + width, xmax + width, nbins)
    ybin_best = np.linspace(ymin , ymax , nbins)
    # Location of gaps 1) between two SiPMS, 2) between two tofpets
    ybin_zoom = np.linspace(188, 192, 100)
    ybin_zoom2 = np.linspace(204, 208, 100)

    # Create histograms
    h_x, edges_x = np.histogram(df_Events.FitPos_ver*10, bins = xbin)
    h_y, edges_y = np.histogram(df_Events.FitPos_hor*10, bins = ybin)
    h_x_simple, edges_x_simple = np.histogram(df_Events.FitPos_ver*10, bins = xbin_best)
    h_y_simple, edges_y_simple = np.histogram(df_Events.FitPos_hor*10, bins = ybin_best)
    h_y_zoom, edges_y_zoom = np.histogram(df_Events.FitPos_hor*10, bins = ybin_zoom)
    h_y_zoom2, edges_y_zoom2 = np.histogram(df_Events.FitPos_hor*10, bins = ybin_zoom2)

    #Efficiency
    horMatches = df_Events[(df_Events['Clu_match_hor'] == 1) & (df_Events['Residual_hor'] < 1)]
    verMatches = df_Events[(df_Events['Clu_match_ver'] == 1) & (df_Events['Residual_ver'] < 1)]
    h_x_match, edges_x_match = np.histogram(verMatches.FitPos_ver*10, bins = xbin)
    h_y_match, edges_y_match = np.histogram(horMatches.FitPos_hor*10, bins = ybin)
    h_x_match_simple, edges_x_match_simple = np.histogram(verMatches.FitPos_ver*10, bins = xbin_best)
    h_y_match_simple, edges_y_match_simple = np.histogram(horMatches.FitPos_hor*10, bins = ybin_best)
    h_y_match_zoom, edges_y_match_zoom = np.histogram(horMatches.FitPos_hor*10, bins = ybin_zoom)
    h_y_match_zoom2, edges_y_match_zoom2 = np.histogram(horMatches.FitPos_hor*10, bins = ybin_zoom2)
    
    xbin = np.delete(xbin, -1)
    ybin = np.delete(ybin, -1)
    xbin_best = np.delete(xbin_best, -1)
    ybin_best = np.delete(ybin_best, -1)
    ybin_zoom = np.delete(ybin_zoom, -1)
    ybin_zoom2 = np.delete(ybin_zoom2, -1)

    # Superpose 
    for i in range(1,6) :
        xbins = np.linspace(xmin + width*i, xmax + width*i, nbins)
        ybins = np.linspace(ymin + width*i, ymax + width*i, nbins)
        h_xtmp, edges_xtmp = np.histogram(df_Events.FitPos_ver*10, bins = xbins)
        h_ytmp, edges_ytmp = np.histogram(df_Events.FitPos_hor*10, bins = ybins)
        h_xtmp_match, edges_xtmp_match = np.histogram(verMatches.FitPos_ver*10, bins = xbins)
        h_ytmp_match, edges_ytmp_match = np.histogram(horMatches.FitPos_hor*10, bins = ybins)
        h_x += h_xtmp
        h_y += h_ytmp
        h_x_match += h_xtmp_match
        h_y_match += h_ytmp_match

    h_x_ref = np.ma.masked_equal(h_x, 0)
    h_y_ref = np.ma.masked_equal(h_y, 0)
    mask_x = np.ma.getmaskarray(h_x_ref)
    mask_y = np.ma.getmaskarray(h_y_ref)
    h_x_match = np.ma.masked_array(h_x_match, mask_x)
    h_y_match = np.ma.masked_array(h_y_match, mask_y)
    h_y_zoom = np.ma.masked_equal(h_y_zoom, 0)
    mask_zoom = np.ma.getmaskarray(h_y_zoom)
    h_y_match_zoom = np.ma.masked_array(h_y_match_zoom, mask_zoom)
    h_y_zoom2 = np.ma.masked_equal(h_y_zoom2, 0)
    mask_zoom2 = np.ma.getmaskarray(h_y_zoom2)
    h_y_match_zoom2 = np.ma.masked_array(h_y_match_zoom2, mask_zoom2)

    # Plot tracks histogram
    # Two SiPMs
    fig, (ax1, ax2) = plt.subplots(2)
    fig.set_size_inches(8,8)
    fig.suptitle(f'Reconstructed track histogram for two SiPMs - Test Station {testStation} - Run {run}')
    ax1.plot(xbin_best, h_x_simple, 'b', drawstyle = 'steps-post')
    ax1.set_xlabel('X position [mm]')
    ax1.set_ylabel('Number of events')
    ax1.set_ylim(bottom=0.0, top=250.0)
    ax1.grid()
    ax2.plot(ybin_best, h_y_simple, 'b', drawstyle = 'steps-post')
    ax2.set_xlabel('Y position[mm]')
    ax2.set_ylabel('Number of events')
    ax2.set_ylim(bottom=0.0, top=250.0)
    ax2.grid()
    plt.savefig(f'figures/Run{run}_simpleSiPMTrackHist_Stat{testStation}_{nbins}bins.png')
    #plt.show()
    plt.close()

    if superpose : 
        # Superposition
        fig, (ax1, ax2) = plt.subplots(2)
        fig.set_size_inches(8,8)
        fig.suptitle(f'Reconstructed track histogram for two SiPMs - Test Station {testStation} - Run {run}')
        ax1.plot(xbin, h_x, 'b', drawstyle = 'steps-post')
        ax1.set_xlabel('X position of first ver SiPM [mm]')
        ax1.set_ylabel('Number of events (cumulated for all SiPMs)')
        ax1.set_ylim(bottom=0.0, top=400.0)
        ax1.grid()
        ax2.plot(ybin, h_y, 'b', drawstyle = 'steps-post')
        ax2.set_xlabel('Y position of first hor SiPM [mm]')
        ax2.set_ylabel('Number of events (cumulated for all SiPMs)')
        ax2.set_ylim(bottom=0.0, top=400.0)
        ax2.grid()
        plt.savefig(f'figures/Run{run}_SiPMTrackHist_Stat{testStation}_{nbins}bins.png')
        #plt.show()
        plt.close()
    

    # Plot Efficiency
    # Two SiPMs
    fig, (ax1, ax2) = plt.subplots(2)
    fig.set_size_inches(8,8)
    fig.suptitle(f'Efficiency for two SiPMs - Test Station {testStation} - Run {run}')
    ax1.plot(xbin_best, h_x_match_simple/h_x_simple, 'b', drawstyle = 'steps-post')
    ax1.set_xlabel('X position [mm]')
    ax1.set_ylabel('Efficiency')
    ax1.set_ylim(bottom=0.0)
    ax1.grid()
    ax2.plot(ybin_best, h_y_match_simple/h_y_simple, 'b', drawstyle = 'steps-post')
    ax2.set_xlabel('Y position [mm]')
    ax2.set_ylabel('Efficiency')
    ax2.set_ylim(bottom=0.0)
    ax2.grid()
    plt.savefig(f'figures/Run{run}_simpleSiPMEfficiency_Stat{testStation}_{nbins}bins.png')
    #plt.show()
    plt.close()
    # Zoom between two SiPMs
    fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols= 2, constrained_layout=True)
    fig.set_size_inches(10,5)
    fig.suptitle(f'Zoom on efficiency in SiPMs gaps - Test Station {testStation} - Run {run}')
    ax1.plot(ybin_zoom, h_y_match_zoom/h_y_zoom, 'b', drawstyle = 'steps-post')
    ax2.plot(ybin_zoom2, h_y_match_zoom2/h_y_zoom2, 'b', drawstyle = 'steps-post')
    ax1.title.set_text('Between two SiPMs')
    ax2.title.set_text('Between two silicon pieces in SiPM')
    ax1.set_xlabel('Y position [mm]', fontsize=12)
    ax1.set_ylabel('Efficiency', fontsize=12)
    ax1.set_ylim(bottom=0.0)
    ax1.grid()
    ax2.set_xlabel('Y position [mm]', fontsize=12)
    ax2.set_ylabel('Efficiency', fontsize=12)
    ax2.set_ylim(bottom=0.0)
    ax2.grid()
    plt.savefig(f'figures/Run{run}_ZoomSiPMEfficiency_Stat{testStation}_{nbins}bins.png')
    #plt.show()
    plt.close()

    if superpose : 
        # Superposition
        fig, (ax1, ax2) = plt.subplots(2)
        fig.set_size_inches(8,8)
        fig.suptitle(f'Mean efficiency for two SiPMs - Test Station {testStation} - Run {run}')
        ax1.plot(xbin, h_x_match/h_x_ref, 'b', drawstyle = 'steps-post')
        ax1.set_xlabel('X position of first ver SiPM [mm]')
        ax1.set_ylabel('Efficiency (mean for all SiPMs)')
        ax1.set_ylim(bottom=0.0)
        ax1.grid()
        ax2.plot(ybin, h_y_match/h_y_ref, 'b', drawstyle = 'steps-post')
        ax2.set_xlabel('Y position of first hor SiPM [mm]')
        ax2.set_ylabel('Efficiency (mean for all SiPMs)')
        ax2.set_ylim(bottom=0.0)
        ax2.grid()
        plt.savefig(f'figures/Run{run}_SiPMEfficiency_Stat{testStation}_{nbins}bins.png')
        #plt.show()
        plt.close()
    
def trackHistograms(df_Events, testStation, run, distrib = True, signal = True, clusize = True, slopes = True) : #added
    '''
    distrib : True to plot track histogram
    signal : True to plot histograms with signal intensity norm
    clusize : True to plot histograms with cluster size norm
    slopes : True to plot histograms with x,y absolute slope (sqrt(x^2+y^2)) and histograms with signed x,y slopes
    '''
    # Selection of matches
    Matches = df_Events[(df_Events['Clu_match_hor'] == 1) & (df_Events['Clu_match_ver'] == 1) & (df_Events['Residual_hor'] < 1) & (df_Events['Residual_ver'] < 1)] #events with hor and ver matches
    horMatches = df_Events[(df_Events['Clu_match_hor'] == 1) & (df_Events['Residual_hor'] < 1)] # events with hor matches
    verMatches = df_Events[(df_Events['Clu_match_ver'] == 1) & (df_Events['Residual_ver'] < 1)] # events with ver matches

    # Histograms
    nbins = 100
    xbins = np.linspace(-500, -50, nbins)
    ybins = np.linspace(130, 570, nbins)
    ones_ref = np.ones(df_Events.shape[0])
    ones_hor = np.ones(horMatches.shape[0])
    ones_ver = np.ones(verMatches.shape[0])
    ones_both = np.ones(Matches.shape[0])

    h_ref, edges_x_ref, edges_y_ref = np.histogram2d(df_Events.FitPos_ver*10, df_Events.FitPos_hor*10, bins=[xbins, ybins], weights = ones_ref)
    h_hor, edges_x_hor, edges_y_hor = np.histogram2d(horMatches.FitPos_ver*10, horMatches.FitPos_hor*10, bins=[xbins, ybins], weights = ones_hor)
    h_ver, edges_x_ver, edges_y_ver = np.histogram2d(verMatches.FitPos_ver*10, verMatches.FitPos_hor*10, bins=[xbins, ybins], weights = ones_ver)
    h_both, edges_x_both, edges_y_both = np.histogram2d(Matches.FitPos_ver*10, Matches.FitPos_hor*10, bins=[xbins, ybins], weights = ones_both)

    # Transpose histograms since histogram doesn't follow cartesian convention
    h_ref = h_ref.T
    h_hor = h_hor.T
    h_ver = h_ver.T
    h_both = h_both.T

    # Mask values = 0
    h_ref = np.ma.masked_equal(h_ref, 0)
    mask_ref = np.ma.getmaskarray(h_ref)
    h_hor = np.ma.masked_equal(h_hor, 0)
    mask_hor = np.ma.getmaskarray(h_hor)
    h_ver = np.ma.masked_equal(h_ver, 0)
    mask_ver = np.ma.getmaskarray(h_ver)
    h_both = np.ma.masked_equal(h_both, 0)
    mask_both = np.ma.getmaskarray(h_both)

    # Signal
    if signal : 
        # Weights
        sig_ref = np.sqrt(np.square(df_Events.CluEnergy_hor) + np.square(df_Events.CluEnergy_ver))
        sig_hor = np.sqrt(np.square(horMatches.CluEnergy_hor) + np.square(horMatches.CluEnergy_ver))
        sig_ver = np.sqrt(np.square(verMatches.CluEnergy_hor) + np.square(verMatches.CluEnergy_ver))
        sig_both = np.sqrt(np.square(Matches.CluEnergy_hor) + np.square(Matches.CluEnergy_ver))
        # Histograms
        hSig_ref, edgesSig_x_ref, edgesSig_y_ref = np.histogram2d(df_Events.FitPos_ver*10, df_Events.FitPos_hor*10, bins=[xbins, ybins], weights = sig_ref)
        hSig_hor, edgesSig_x_hor, edgesSig_y_hor = np.histogram2d(horMatches.FitPos_ver*10, horMatches.FitPos_hor*10, bins=[xbins, ybins], weights = sig_hor)
        hSig_ver, edgesSig_x_ver, edgesSig_y_ver = np.histogram2d(verMatches.FitPos_ver*10, verMatches.FitPos_hor*10, bins=[xbins, ybins], weights = sig_ver)
        hSig_both, edgesSig_x_both, edgesSig_y_both = np.histogram2d(Matches.FitPos_ver*10, Matches.FitPos_hor*10, bins=[xbins, ybins], weights = sig_both)
        # Transpose histograms since histogram doesn't follow cartesian convention
        hSig_ref = hSig_ref.T
        hSig_hor = hSig_hor.T
        hSig_ver = hSig_ver.T
        hSig_both = hSig_both.T
        # Mask values = 0 in hists
        hSig_ref = np.ma.masked_array(hSig_ref, mask_ref)
        hSig_hor = np.ma.masked_array(hSig_hor, mask_hor)
        hSig_ver = np.ma.masked_array(hSig_ver, mask_ver)
        hSig_both = np.ma.masked_array(hSig_both, mask_both)
        # Histograms to plot + rescaling
        Sig_ref = hSig_ref/h_ref
        Sig_hor = hSig_hor/h_hor
        Sig_ver = hSig_ver/h_ver
        Sig_both = hSig_both/h_both
        Sig_ref = np.ma.masked_greater(Sig_ref, 12)
        Sig_hor = np.ma.masked_greater(Sig_hor, 12)
        Sig_ver = np.ma.masked_greater(Sig_ver, 12)
        Sig_both = np.ma.masked_greater(Sig_both, 12)

    # Cluster size
    if clusize : 
        # Weights
        size_ref = np.sqrt(np.square(df_Events.CluSize_hor) + np.square(df_Events.CluSize_ver))
        size_hor = np.sqrt(np.square(horMatches.CluSize_hor) + np.square(horMatches.CluSize_ver))
        size_ver = np.sqrt(np.square(verMatches.CluSize_hor) + np.square(verMatches.CluSize_ver))
        size_both = np.sqrt(np.square(Matches.CluSize_hor) + np.square(Matches.CluSize_ver))
        # Histograms
        hSize_ref, edgesSize_x_ref, edgesSize_y_ref = np.histogram2d(df_Events.FitPos_ver*10, df_Events.FitPos_hor*10, bins=[xbins, ybins], weights = size_ref)
        hSize_hor, edgesSize_x_hor, edgesSize_y_hor = np.histogram2d(horMatches.FitPos_ver*10, horMatches.FitPos_hor*10, bins=[xbins, ybins], weights = size_hor)
        hSize_ver, edgesSize_x_ver, edgesSize_y_ver = np.histogram2d(verMatches.FitPos_ver*10, verMatches.FitPos_hor*10, bins=[xbins, ybins], weights = size_ver)
        hSize_both, edgesSize_x_both, edgesSize_y_both = np.histogram2d(Matches.FitPos_ver*10, Matches.FitPos_hor*10, bins=[xbins, ybins], weights = size_both)
        # Transpose histograms since histogram doesn't follow cartesian convention
        hSize_ref = hSize_ref.T
        hSize_hor = hSize_hor.T
        hSize_ver = hSize_ver.T
        hSize_both = hSize_both.T
        # Mask values = 0 in hists
        hSize_ref = np.ma.masked_array(hSize_ref, mask_ref)
        hSize_hor = np.ma.masked_array(hSize_hor, mask_hor)
        hSize_ver = np.ma.masked_array(hSize_ver, mask_ver)
        hSize_both = np.ma.masked_array(hSize_both, mask_both)
        # Histograms to plot + rescaling
        Size_ref = hSize_ref/h_ref
        Size_hor = hSize_hor/h_hor
        Size_ver = hSize_ver/h_ver
        Size_both = hSize_both/h_both
        Size_ref = np.ma.masked_greater(Size_ref, 3)
        Size_hor = np.ma.masked_greater(Size_hor, 3)
        Size_ver = np.ma.masked_greater(Size_ver, 3)
        Size_both = np.ma.masked_greater(Size_both, 3)

    # Slopes
    if slopes : 
        # Weights + conversion into degrees
        # Norm
        slope_ref = np.sqrt(np.square(df_Events.FitSlope_y) + np.square(df_Events.FitSlope_x))*180/np.pi
        slope_hor = np.sqrt(np.square(horMatches.FitSlope_y) + np.square(horMatches.FitSlope_x))*180/np.pi
        slope_ver = np.sqrt(np.square(verMatches.FitSlope_y) + np.square(verMatches.FitSlope_x))*180/np.pi
        slope_both = np.sqrt(np.square(Matches.FitSlope_y) + np.square(Matches.FitSlope_x))*180/np.pi
        # x slope 
        slopex_ref = df_Events.FitSlope_x*180/np.pi
        slopex_hor = horMatches.FitSlope_x*180/np.pi
        slopex_ver = verMatches.FitSlope_x*180/np.pi
        slopex_both = Matches.FitSlope_x*180/np.pi
        # y slope
        slopey_ref = df_Events.FitSlope_y*180/np.pi
        slopey_hor = horMatches.FitSlope_y*180/np.pi
        slopey_ver = verMatches.FitSlope_y*180/np.pi
        slopey_both = Matches.FitSlope_y*180/np.pi


        # Histograms
        # Norm
        hSlope_ref, edgesSlope_x_ref, edgesSlope_y_ref = np.histogram2d(df_Events.FitPos_ver*10, df_Events.FitPos_hor*10, bins=[xbins, ybins], weights = slope_ref)
        hSlope_hor, edgesSlope_x_hor, edgesSlope_y_hor = np.histogram2d(horMatches.FitPos_ver*10, horMatches.FitPos_hor*10, bins=[xbins, ybins], weights = slope_hor)
        hSlope_ver, edgesSlope_x_ver, edgesSlope_y_ver = np.histogram2d(verMatches.FitPos_ver*10, verMatches.FitPos_hor*10, bins=[xbins, ybins], weights = slope_ver)
        hSlope_both, edgesSlope_x_both, edgesSlope_y_both = np.histogram2d(Matches.FitPos_ver*10, Matches.FitPos_hor*10, bins=[xbins, ybins], weights = slope_both)
        # x slope
        hSlopex_ref, edgesSlopex_x_ref, edgesSlopex_y_ref = np.histogram2d(df_Events.FitPos_ver*10, df_Events.FitPos_hor*10, bins=[xbins, ybins], weights = slopex_ref)
        hSlopex_hor, edgesSlopex_x_hor, edgesSlopex_y_hor = np.histogram2d(horMatches.FitPos_ver*10, horMatches.FitPos_hor*10, bins=[xbins, ybins], weights = slopex_hor)
        hSlopex_ver, edgesSlopex_x_ver, edgesSlopex_y_ver = np.histogram2d(verMatches.FitPos_ver*10, verMatches.FitPos_hor*10, bins=[xbins, ybins], weights = slopex_ver)
        hSlopex_both, edgesSlopex_x_both, edgesSlopex_y_both = np.histogram2d(Matches.FitPos_ver*10, Matches.FitPos_hor*10, bins=[xbins, ybins], weights = slopex_both)
        # y slope
        hSlopey_ref, edgesSlopey_x_ref, edgesSlopey_y_ref = np.histogram2d(df_Events.FitPos_ver*10, df_Events.FitPos_hor*10, bins=[xbins, ybins], weights = slopey_ref)
        hSlopey_hor, edgesSlopey_x_hor, edgesSlopey_y_hor = np.histogram2d(horMatches.FitPos_ver*10, horMatches.FitPos_hor*10, bins=[xbins, ybins], weights = slopey_hor)
        hSlopey_ver, edgesSlopey_x_ver, edgesSlopey_y_ver = np.histogram2d(verMatches.FitPos_ver*10, verMatches.FitPos_hor*10, bins=[xbins, ybins], weights = slopey_ver)
        hSlopey_both, edgesSlopey_x_both, edgesSlopey_y_both = np.histogram2d(Matches.FitPos_ver*10, Matches.FitPos_hor*10, bins=[xbins, ybins], weights = slopey_both)
        
        # Transpose histograms since histogram doesn't follow cartesian convention
        # Norm
        hSlope_ref = hSlope_ref.T
        hSlope_hor = hSlope_hor.T
        hSlope_ver = hSlope_ver.T
        hSlope_both = hSlope_both.T
        # x slope
        hSlopex_ref = hSlopex_ref.T
        hSlopex_hor = hSlopex_hor.T
        hSlopex_ver = hSlopex_ver.T
        hSlopex_both = hSlopex_both.T
        # y slope
        hSlopey_ref = hSlopey_ref.T
        hSlopey_hor = hSlopey_hor.T
        hSlopey_ver = hSlopey_ver.T
        hSlopey_both = hSlopey_both.T

        # Mask values = 0 in hists
        # Norm
        hSlope_ref = np.ma.masked_array(hSlope_ref, mask_ref)
        hSlope_hor = np.ma.masked_array(hSlope_hor, mask_hor)
        hSlope_ver = np.ma.masked_array(hSlope_ver, mask_ver)
        hSlope_both = np.ma.masked_array(hSlope_both, mask_both)
        # x slope
        hSlopex_ref = np.ma.masked_array(hSlopex_ref, mask_ref)
        hSlopex_hor = np.ma.masked_array(hSlopex_hor, mask_hor)
        hSlopex_ver = np.ma.masked_array(hSlopex_ver, mask_ver)
        hSlopex_both = np.ma.masked_array(hSlopex_both, mask_both)
        # y slope
        hSlopey_ref = np.ma.masked_array(hSlopey_ref, mask_ref)
        hSlopey_hor = np.ma.masked_array(hSlopey_hor, mask_hor)
        hSlopey_ver = np.ma.masked_array(hSlopey_ver, mask_ver)
        hSlopey_both = np.ma.masked_array(hSlopey_both, mask_both)

        # Histograms to plot + rescaling
        # Norm
        Slope_ref = hSlope_ref/h_ref
        Slope_hor = hSlope_hor/h_hor
        Slope_ver = hSlope_ver/h_ver
        Slope_both = hSlope_both/h_both
        Slope_ref = np.ma.masked_greater(Slope_ref, 6)
        Slope_hor = np.ma.masked_greater(Slope_hor, 6)
        Slope_ver = np.ma.masked_greater(Slope_ver, 6)
        Slope_both = np.ma.masked_greater(Slope_both, 6)
        # x slope
        scale = 3
        Slopex_ref = hSlopex_ref/h_ref
        Slopex_hor = hSlopex_hor/h_hor
        Slopex_ver = hSlopex_ver/h_ver
        Slopex_both = hSlopex_both/h_both
        Slopex_ref = np.ma.masked_greater(Slopex_ref, scale)
        Slopex_hor = np.ma.masked_greater(Slopex_hor, scale)
        Slopex_ver = np.ma.masked_greater(Slopex_ver, scale)
        Slopex_both = np.ma.masked_greater(Slopex_both, scale)
        Slopex_ref = np.ma.masked_less(Slopex_ref, - scale)
        Slopex_hor = np.ma.masked_less(Slopex_hor, - scale)
        Slopex_ver = np.ma.masked_less(Slopex_ver, - scale)
        Slopex_both = np.ma.masked_less(Slopex_both, - scale)
        # y slope 
        Slopey_ref = hSlopey_ref/h_ref
        Slopey_hor = hSlopey_hor/h_hor
        Slopey_ver = hSlopey_ver/h_ver
        Slopey_both = hSlopey_both/h_both
        Slopey_ref = np.ma.masked_greater(Slopey_ref, scale)
        Slopey_hor = np.ma.masked_greater(Slopey_hor, scale)
        Slopey_ver = np.ma.masked_greater(Slopey_ver, scale)
        Slopey_both = np.ma.masked_greater(Slopey_both, scale)
        Slopey_ref = np.ma.masked_less(Slopey_ref, - scale)
        Slopey_hor = np.ma.masked_less(Slopey_hor, - scale)
        Slopey_ver = np.ma.masked_less(Slopey_ver, - scale)
        Slopey_both = np.ma.masked_less(Slopey_both, - scale)

    # Plots
    # Distribution
    if distrib : 
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, constrained_layout=True)
        fig.set_size_inches(8, 8)
        fig.suptitle(f'Track reconstruction histograms - Test Station {testStation} - Run {run}')
        plot1 = ax1.imshow(h_ref, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_ref[0], edges_x_ref[-1], edges_y_ref[0], edges_y_ref[-1]]) 
        plot2 = ax2.imshow(h_both, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_both[0], edges_x_both[-1], edges_y_both[0], edges_y_both[-1]])
        plot3 = ax3.imshow(h_hor, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_hor[0], edges_x_hor[-1], edges_y_hor[0], edges_y_hor[-1]])
        plot4 = ax4.imshow(h_ver, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_ver[0], edges_x_ver[-1], edges_y_ver[0], edges_y_ver[-1]])
        ax1.title.set_text('No selection on match')
        ax2.title.set_text('Horizontal and vertical matches')
        ax3.title.set_text('Horizontal match')
        ax4.title.set_text('Vertical match')
        ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax4.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax4.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Number of events', size=12)
        cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='Number of events', size=12)
        cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='Number of events', size=12)
        cbar4 = plt.colorbar(plot4, ax = ax4).set_label(label='Number of events', size=12)
        plt.savefig(f'figures/Run{run}_trackhist-distrib_Stat{testStation}_{nbins}bins.png')
        #plt.show()
        plt.close()
    
    # Signal
    if signal : 
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, constrained_layout=True)
        fig.set_size_inches(8,8)
        fig.suptitle(f'Energy deposit histograms - Test Station {testStation} - Run {run}')
        plot1 = ax1.imshow(Sig_ref, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_ref[0], edges_x_ref[-1], edges_y_ref[0], edges_y_ref[-1]]) 
        plot2 = ax2.imshow(Sig_both, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_both[0], edges_x_both[-1], edges_y_both[0], edges_y_both[-1]])
        plot3 = ax3.imshow(Sig_hor, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_hor[0], edges_x_hor[-1], edges_y_hor[0], edges_y_hor[-1]])
        plot4 = ax4.imshow(Sig_ver, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_ver[0], edges_x_ver[-1], edges_y_ver[0], edges_y_ver[-1]])
        ax1.title.set_text('No selection on match')
        ax2.title.set_text('Horizontal and vertical matches')
        ax3.title.set_text('Horizontal match')
        ax4.title.set_text('Vertical match')
        ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax4.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax4.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Energy deposit', size=12)
        cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='Energy deposit', size=12)
        cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='Energy deposit', size=12)
        cbar4 = plt.colorbar(plot4, ax = ax4).set_label(label='Energy deposit', size=12)
        plt.savefig(f'figures/Run{run}_trackhist-signal_Stat{testStation}_{nbins}bins.png')
        #plt.show()
        plt.close()

    # Cluster size
    if clusize : 
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, constrained_layout=True)
        fig.set_size_inches(8,8)
        fig.suptitle(f'Cluster size histograms - Test Station {testStation} - Run {run}')
        plot1 = ax1.imshow(Size_ref, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_ref[0], edges_x_ref[-1], edges_y_ref[0], edges_y_ref[-1]]) 
        plot2 = ax2.imshow(Size_both, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_both[0], edges_x_both[-1], edges_y_both[0], edges_y_both[-1]])
        plot3 = ax3.imshow(Size_hor, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_hor[0], edges_x_hor[-1], edges_y_hor[0], edges_y_hor[-1]])
        plot4 = ax4.imshow(Size_ver, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_ver[0], edges_x_ver[-1], edges_y_ver[0], edges_y_ver[-1]])
        ax1.title.set_text('No selection on match')
        ax2.title.set_text('Horizontal and vertical matches')
        ax3.title.set_text('Horizontal match')
        ax4.title.set_text('Vertical match')
        ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax4.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax4.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Cluster size', size=12)
        cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='Cluster size', size=12)
        cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='Cluster size', size=12)
        cbar4 = plt.colorbar(plot4, ax = ax4).set_label(label='Cluster size', size=12)
        plt.savefig(f'figures/Run{run}_trackhist-clusize_Stat{testStation}_{nbins}bins.png')
        #plt.show()
        plt.close()

    # Slope
    if slopes : 
        # Norm
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, constrained_layout=True)
        fig.set_size_inches(8,8)
        fig.suptitle(f'Slopes histograms - Test Station {testStation} - Run {run}')
        plot1 = ax1.imshow(Slope_ref, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_ref[0], edges_x_ref[-1], edges_y_ref[0], edges_y_ref[-1]]) 
        plot2 = ax2.imshow(Slope_both, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_both[0], edges_x_both[-1], edges_y_both[0], edges_y_both[-1]])
        plot3 = ax3.imshow(Slope_hor, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_hor[0], edges_x_hor[-1], edges_y_hor[0], edges_y_hor[-1]])
        plot4 = ax4.imshow(Slope_ver, cmap= plt.get_cmap('viridis'), origin='lower', extent=[edges_x_ver[0], edges_x_ver[-1], edges_y_ver[0], edges_y_ver[-1]])
        ax1.title.set_text('No selection on match')
        ax2.title.set_text('Horizontal and vertical matches')
        ax3.title.set_text('Horizontal match')
        ax4.title.set_text('Vertical match')
        ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax4.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax4.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='Absolute slope [deg]', size=12)
        cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='Absolute slope [deg]', size=12)
        cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='Absolute slope [deg]', size=12)
        cbar4 = plt.colorbar(plot4, ax = ax4).set_label(label='Absolute slope [deg]', size=12)
        plt.savefig(f'figures/Run{run}_trackhist-slopes_Stat{testStation}_{nbins}bins.png')
        #plt.show()
        plt.close()

        # x slope
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, constrained_layout=True)
        fig.set_size_inches(8,8)
        fig.suptitle(f'x slopes histograms - Test Station {testStation} - Run {run}')
        plot1 = ax1.imshow(Slopex_ref, cmap= plt.get_cmap('Spectral'), origin='lower', extent=[edges_x_ref[0], edges_x_ref[-1], edges_y_ref[0], edges_y_ref[-1]]) 
        plot2 = ax2.imshow(Slopex_both, cmap= plt.get_cmap('Spectral'), origin='lower', extent=[edges_x_both[0], edges_x_both[-1], edges_y_both[0], edges_y_both[-1]])
        plot3 = ax3.imshow(Slopex_hor, cmap= plt.get_cmap('Spectral'), origin='lower', extent=[edges_x_hor[0], edges_x_hor[-1], edges_y_hor[0], edges_y_hor[-1]])
        plot4 = ax4.imshow(Slopex_ver, cmap= plt.get_cmap('Spectral'), origin='lower', extent=[edges_x_ver[0], edges_x_ver[-1], edges_y_ver[0], edges_y_ver[-1]])
        ax1.title.set_text('No selection on match')
        ax2.title.set_text('Horizontal and vertical matches')
        ax3.title.set_text('Horizontal match')
        ax4.title.set_text('Vertical match')
        ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax4.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax4.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='x slope [deg]', size=12)
        cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='x slope [deg]', size=12)
        cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='x slope [deg]', size=12)
        cbar4 = plt.colorbar(plot4, ax = ax4).set_label(label='x slope [deg]', size=12)
        plt.savefig(f'figures/Run{run}_trackhist-Xslopes_Stat{testStation}_{nbins}bins.png')
        #plt.show()
        plt.close()

        # y slope
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2, constrained_layout=True)
        fig.set_size_inches(8,8)
        fig.suptitle(f'y slopes histograms - Test Station {testStation} - Run {run}')
        plot1 = ax1.imshow(Slopey_ref, cmap= plt.get_cmap('Spectral'), origin='lower', extent=[edges_x_ref[0], edges_x_ref[-1], edges_y_ref[0], edges_y_ref[-1]]) 
        plot2 = ax2.imshow(Slopey_both, cmap= plt.get_cmap('Spectral'), origin='lower', extent=[edges_x_both[0], edges_x_both[-1], edges_y_both[0], edges_y_both[-1]])
        plot3 = ax3.imshow(Slopey_hor, cmap= plt.get_cmap('Spectral'), origin='lower', extent=[edges_x_hor[0], edges_x_hor[-1], edges_y_hor[0], edges_y_hor[-1]])
        plot4 = ax4.imshow(Slopey_ver, cmap= plt.get_cmap('Spectral'), origin='lower', extent=[edges_x_ver[0], edges_x_ver[-1], edges_y_ver[0], edges_y_ver[-1]])
        ax1.title.set_text('No selection on match')
        ax2.title.set_text('Horizontal and vertical matches')
        ax3.title.set_text('Horizontal match')
        ax4.title.set_text('Vertical match')
        ax1.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax1.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax2.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax2.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax3.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax3.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        ax4.set_xlabel('X position [mm] (vertical fibers)', fontsize=12)
        ax4.set_ylabel('Y position [mm] (horizontal fibers)', fontsize=12)
        cbar1 = plt.colorbar(plot1, ax = ax1).set_label(label='y slope [deg]', size=12)
        cbar2 = plt.colorbar(plot2, ax = ax2).set_label(label='y slope [deg]', size=12)
        cbar3 = plt.colorbar(plot3, ax = ax3).set_label(label='y slope [deg]', size=12)
        cbar4 = plt.colorbar(plot4, ax = ax4).set_label(label='y slope [deg]', size=12)
        plt.savefig(f'figures/Run{run}_trackhist-Yslopes_Stat{testStation}_{nbins}bins.png')
        #plt.show()
        plt.close()
    