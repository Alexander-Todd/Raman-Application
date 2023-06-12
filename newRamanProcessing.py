import pandas as pd
import glob
import numpy as np
import copy
from tqdm.notebook import tqdm
from scipy import interpolate
from scipy.signal import savgol_filter
import pywt # wavelet smoother
# from scipy.sparse import csc_matrix, eye, diags
from scipy import sparse  # For whittaker smoother
from scipy.sparse.linalg import splu  # For whittaker smoother
from BaselineRemoval import BaselineRemoval
from statsmodels.robust import mad  # For wavelet smooth
from scipy.sparse import csc_matrix, eye, diags  # For WhittakerSmooth
from scipy.sparse.linalg import spsolve  # For WhittakerSmooth
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit  # for lorentzian fitting
import scipy.integrate  # For area under peak
import matplotlib.pyplot as mpl
from numpy.polynomial.polynomial import polyfit

# ---------------------------------------------------------------------------------------------------------------------------------------
# Unchanged process defs
# ---------------------------------------------------------------------------------------------------------------------------------------

# def smooth(array, method = 'Savitzky–Golay', window = 3, polynomial = 0, axis = 1, fourier_values = 3, wavelet = 'db29',
#            wavelet_level = 1, lambda_val = 50000, d = 2):
def smooth(wns, array, method_dict):
    axis = 1
    [[method, parameters]] = method_dict.items()
    variables = list(parameters.values())
    # variables = [eval(i) for i in vars_strings]
        
    array = np.transpose(np.stack(array))
    smoothed_array = np.zeros(np.shape(array))

    if method == 'Savitzky–Golay':
        window, polynomial = variables
        for spectra in range(np.shape(array)[axis]):
            smoothed_array[:,spectra] = savgol_filter(array[:,spectra],window,polynomial)

    elif method == 'FFT':
        fourier_values = int(variables[0])
        for spectra in range(np.shape(array)[axis]):
            padded_array = np.pad(array[:,spectra],
                                 (100, 100), 'constant',
                                 constant_values=((array[:,spectra][0],array[:,spectra][-1])))
            rft = np.fft.rfft(padded_array)
            rft[fourier_values:] = 0
            # error created if array has odd length as rft halves and rounds down
            # and when doubled again the sizes don't match
            if np.shape(array)[0] % 2 == 0:
                theproblem = np.fft.irfft(rft)[100:-100]
            else:
                theproblem = np.fft.irfft(rft)[99:-100]
            smoothed_array[:,spectra] = theproblem

    elif method == 'Wavelet':
        wavelet, wavelet_level = variables
        for spectra in range(np.shape(array)[axis]):
            smoothed_array[:,spectra] = waveletSmooth(array[:,spectra], wavelet=wavelet, level=wavelet_level)

    elif method == 'Whittaker':
        lambda_val, d = variables
        for spectra in range(np.shape(array)[axis]):
            smoothed_array[:,spectra] = whittaker_smooth(array[:,spectra], lambda_val, d = d)
    
    else:
        raise Exception("Smoothing method not found!")
        
    return np.transpose(smoothed_array)

def waveletSmooth( x, wavelet="db4", level=1):
    # calculate the wavelet coefficients
    coeff = pywt.wavedec( x, wavelet, mode="per" )
    # calculate a threshold
    sigma = mad( coeff[-level] )
    # changing this threshold also changes the behavior,
    # but I have not played with this very much
    uthresh = sigma * np.sqrt( 2*np.log( len( x ) ) )
    coeff[1:] = ( pywt.threshold( i, value=uthresh, mode="soft" ) for i in coeff[1:] )
    # reconstruct the signal using the thresholded coefficients
    y = pywt.waverec( coeff, wavelet, mode="per" )
    return y


def speyediff(N, d, format='csc'): # function for whittaker smoother
    """
    (utility function)
    Construct a d-th order sparse difference matrix based on
    an initial N x N identity matrix

    Final matrix (N-d) x N
    """
    assert not (d < 0), "d must be non negative"
    shape     = (N-d, N)
    diagonals = np.zeros(2*d + 1)
    diagonals[d] = 1.
    for i in range(d):
        diff = diagonals[:-1] - diagonals[1:]
        diagonals = diff
    offsets = np.arange(d+1)
    spmat = sparse.diags(diagonals, offsets, shape, format=format)
    return spmat

def whittaker_smooth(y, lmbd, d = 2):
    """
    Implementation of the Whittaker smoothing algorithm,
    based on the work by Eilers [1].
    [1] P. H. C. Eilers, "A perfect smoother", Anal. Chem. 2003, (75), 3631-3636

    The larger 'lmbd', the smoother the data.
    For smoothing of a complete data series, sampled at equal intervals
    This implementation uses sparse matrices enabling high-speed processing
    of large input vectors

    ---------
    Arguments :
    y       : vector containing raw data
    lmbd    : parameter for the smoothing algorithm (roughness penalty)
    d       : order of the smoothing
    ---------

    Returns :
    z       : vector of the smoothed data.
    """

    m = len(y)
    E = sparse.eye(m, format='csc')
    D = speyediff(m, d, format='csc') # defined above
    coefmat = E + lmbd * D.conj().T.dot(D)
    z = splu(coefmat).solve(y)
    return z

# Why are there two Whittaker smoothing functions??
def WhittakerSmooth(x,w,lambda_,differences=1):
    '''
    Penalized least squares algorithm for background fitting

    input
        x: input data (i.e. chromatogram of spectrum)
        w: binary masks (value of the mask is zero if a point belongs to peaks and one otherwise)
        lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background
        differences: integer indicating the order of the difference of penalties

    output
        the fitted background vector

    https://github.com/zmzhang/airPLS/blob/master/airPLS.py
    '''
    X=np.matrix(x)
    m=X.size
    i=np.arange(0,m)
    E=eye(m,format='csc')
    D=E[1:]-E[:-1] # numpy.diff() does not work with sparse matrix. This is a workaround.
    W=diags(w,0,shape=(m,m))
    A=csc_matrix(W+(lambda_*D.T*D))
    B=csc_matrix(W*X.T)
    background=spsolve(A,B)
    return np.array(background)

# def removeCosmicRaySpikes(array, method = 'wavenumber', threshold_diff=5, threshold_wn=5, return_CRS_positions=False):

def removeCosmicRaySpikes(wns, array, method_dict):

    [[method, parameters]] = method_dict.items()
    parameters = list(parameters.values())

    if method == 'Modified Z-Score':
        score_threshold = int(parameters[0])
        order = int(parameters[1])
        width = int(parameters[2])

        output = []
        for spectrum in array:
            processed = fixer(spectrum, width, score_threshold, order)
            output.append(processed)
        return output

    return_CRS_positions = False
    array = np.transpose(np.stack(array))
    removed_spikes = []
    despiked_raman_array = copy.deepcopy(array)

    if method == 'Wavenumber':
        threshold_wn = int(parameters[0])
        for row, x in enumerate(array):
            for index, y in enumerate(x):
                if y > np.mean(x) + np.std(x):
                    new_array = np.delete(x, index)
                    if  y > np.mean(new_array) + threshold_wn * np.std(new_array):
                        despiked_raman_array[row,index] = np.mean(new_array)
                        removed_spikes.append([index,row,y])

        if return_CRS_positions == True:
            return despiked_raman_array, removed_spikes
        else:
            return np.transpose(despiked_raman_array)

    elif method == 'Differential':
        threshold_diff = int(parameters[0])
        intensitites = array
        spikes = 1
        while spikes != 0:
            spikes_desc = spikes
            spikes = 0
            # the mean of the standard deviation of the difference of the differences
            min_threshold = 0 - np.mean(np.std(np.diff(intensitites,2),axis=1))*threshold_diff
            for spectrum in tqdm(range(np.shape(intensitites)[0]),desc='Despike (' + str(spikes_desc) + ' Spikes Remaing)' , leave=False):
                second_diff_spectrum = np.diff(intensitites[spectrum,:],2)
                for wavenumber in range(np.shape(intensitites)[1]-2):
                    if min_threshold > second_diff_spectrum[wavenumber] or abs(min_threshold) < second_diff_spectrum[wavenumber]:
                        removed_spikes.append([wavenumber+1,spectrum,intensitites[spectrum,wavenumber+1]])
                        intensitites[spectrum,wavenumber+1] = intensitites[spectrum,wavenumber+1]*0.95
                        #spikes += 1
        if return_CRS_positions == True:
            return intensitites, removed_spikes
        else:
            return np.transpose(intensitites)

    elif method == 'Consensus':
        threshold_wn, threshold_diff = parameters
        removed_spikes_wn = []
        removed_spikes_diff = []
        despiked_raman_array = copy.deepcopy(array)
        row = 0
        for x in tqdm(array, desc='Despike', leave=False):
            index = 0
            for y in x:
                if y > np.mean(x) + np.std(x):
                    new_array = np.delete(x, index)
                    if  y > np.mean(new_array) + threshold_wn * np.std(new_array):
                        removed_spikes_wn.append([index,row,y,np.mean(new_array)])
                        #print('spike')
                index += 1
            row += 1

        intensitites = array
        spikes = 1
        while spikes != 0:
            spikes_desc = spikes
            spikes = 0
            min_threshold = 0 - np.mean(np.std(np.diff(intensitites,2),axis=1))*threshold_diff
            for spectrum in tqdm(range(np.shape(intensitites)[0]),desc='Despike (' + str(spikes_desc) + ' Spikes Remaing)' , leave=False):
                second_diff_spectrum = np.diff(intensitites[spectrum,:],2)
                for wavenumber in range(np.shape(intensitites)[1]-2):
                    if min_threshold > second_diff_spectrum[wavenumber] or abs(min_threshold) < second_diff_spectrum[wavenumber]:
                        removed_spikes_diff.append([wavenumber+1,spectrum,intensitites[spectrum,wavenumber+1]])

        #print(np.shape(np.array(removed_spikes_wn)))
        #print(np.shape(np.array(removed_spikes_diff)))
        try:
            mask = np.isin(np.array(removed_spikes_wn)[:,1:3],np.array(removed_spikes_diff)[:,1:3])
            detected_spikes = np.array(removed_spikes_wn)[mask[:,0],0:2].astype(int)
            new_values = np.array(removed_spikes_wn)[mask[:,0],3]
            index = 0
            for x in detected_spikes[:,0]:
                despiked_raman_array[detected_spikes[index,1],x] = new_values[index]
                index += 1
        # Catch the times that the first method detects no peaks
        except IndexError:
            pass
        return np.transpose(despiked_raman_array)
    else:
        print('Error: Method not found = ' + str(method))

# def z_score(intensity):
#  mean_int = np.mean(intensity)
#  std_int = np.std(intensity)
#  z_scores = (intensity - mean_int) / std_int
#  return z_scores

def modified_z_score(intensity):
 median_int = np.median(intensity)
 mad_int = np.median([np.abs(intensity - median_int)])
 modified_z_scores = 0.6745 * (intensity - median_int) / mad_int
 return modified_z_scores

def fixer(y,m, threshold, order): 
    spikes = abs(np.array(modified_z_score(np.diff(y, order)))) > threshold
    y_out = y.copy() # So we don’t overwrite y
    for i in np.arange(len(spikes)):
        if spikes[i] != 0: # If we have an spike in position i
            init_w = np.arange(i-m,i+1+m) # we select 2 m + 1 points around our spike
            n = np.arange(i-2, i+2) # we remove points ajdacent to the spike
            w = np.setdiff1d(init_w, n)
            w2 = w[spikes[w] == 0] # From such interval, we choose the ones which are not spikes
            value = np.mean(y[w2])
            y_out[i-2] = value
            y_out[i-1] = value # and we average their values
            y_out[i] = value
            y_out[i+1] = value
            y_out[i+2] = value
    
    return y_out


def new_fixer(wns, intensities, width=1):
    # Get the spikes wns in the brackets
    spikes = []
    spikes.sort()  # arrange in ascending order
    spikes_left=spikes  # one array for changing, other for reference
    intensities_out = intensities.copy()

    # Since arranged we know position 0 is always the lowest
    # Remove spikes from array until none left
    while len(spikes_left) > 0:
        spikes_subtracted = []
        min_wn = wns.index(spikes_left[0]) - width  # index of wn at lower end
        max_wn = wns.index(spikes_left[0]) + width  # index of wn at higher end
        spikes_subtracted.append(spikes_left[0])  # wn val of removed spike
        spikes_left = np.delete(spikes_left, 0)  # remove the value we used
        # In case there are other spikes in the range we have selected:
        while True:
            if spikes_left[0] < wns[max_wn]:
                max_wn = wns.index(spikes_left[0]) + width
                spikes_subtracted.append(spikes_left[0])
                spikes_left = np.delete(spikes_left, 0)
            else:
                break
        # Make x and y arrays for the polyfit
        wns_range = wns[min_wn:max_wn]
        ints_range = intensities_out[min_wn:max_wn]
        # Remove the spike values from the ranges
        for i in spikes_subtracted:
            index = wns_range.index(i)
            wns_range = np.delete(wns_range, index)
            ints_range = np.delete(ints_range, index)
        degree = len(wns_range) - 1
        polynomial = np.polyfit(wns_range, ints_range, degree)
        coefficients = len(polynomial)
        for spike in spikes_subtracted:
            val = 0
            for i in range(coefficients):
                val = val + polynomial[i]*(spike**i)
            intensities_out[wns.index(spike)] = val
    
    return intensities_out
    

def new_remove_crs(wns, intensities, method_dict):
    # threshold should be parameter
    # Want diff, threshold, width?, degree?
    [[method, parameters]] = method_dict.items()
    parameters = list(parameters.values())

    if method == "Modified Z-Score":
        threshold = int(parameters[0])
        order = int(parameters[1])
        width = int(parameters[2])
        try:
            graphing = bool(parameters[3])
        except TypeError:
            print("""Invalid boolean value given for modified Z score graphing.
            Value should be 0 (False) or 1 (True).
            Will assume graphing to be False.""")
            graphing = False
        output = []
        for count, spectrum in enumerate(intensities):
            processed = modifiedZscore(wns[count], spectrum, threshold, order, width, graphing)
            output.append(processed)
        return output
    
    else:
        print("Only Modified Z score has been implemented")
        return intensities

def modifiedZscore(wns, spectrum, threshold, diff, width, graphing):
    """Takes one spectrum as input, including x and y values. Calculates
    the modified Z score according to the differential order given and
    removes the spikes at the points where the score is above the threshold.
    A cubic least squares polynomial is fitted to the non-spike points
    within the window (spike position +/- the width).
    """
    # Modified Z score calculation:
    delta_spec = np.diff(spectrum, diff)
    median_int = np.median(delta_spec)
    mad_int = np.median([np.abs(delta_spec - median_int)])
    modZscore = 0.6745 * (delta_spec - median_int) / mad_int
    modZscore = np.array(np.abs(modZscore))
    # wns shortens by 1 for each diff, is like a pyramid so
    # is lost from both ends working inward
    accessibleWns = wns
    for i in range(diff):
        if (i+1) % 2 == 0:
            accessibleWns = np.delete(accessibleWns, 0)
        else:   
            accessibleWns = np.delete(accessibleWns, -1)
    spikes = []  # storage for the spikes found
    # for singleArrayScores in intensity_modified_z_score: 
    for i in range(len(modZscore)):
        if modZscore[i] > threshold:
            # Need to correct for the offset from diff, not sure why -2 though
            spikes_wn = accessibleWns[i]
            print("Threshold exceeded at", spikes_wn)
            spikes.append(spikes_wn)

    # Graphing useful to understand if there is a threshold that will work
    if graphing:
        plt.plot(np.transpose(accessibleWns), np.transpose(modZscore))
        plt.plot(np.transpose(accessibleWns), threshold*np.ones(len(accessibleWns)), label = 'threshold')
        plt.title('Modified Z-Score of ∇x(i)', fontsize = 20)
        plt.xticks(fontsize = 15)
        plt.yticks(fontsize = 15)
        plt.xlabel('Wavenumbers ($cm^{-1}$)', fontsize = 20)
        plt.ylabel('Score', fontsize = 20)
        plt.show()

    spikes.sort()  # arrange in ascending order
    spikes_left = spikes  # one array for changing, other for reference
    spectrum_out = spectrum.copy() # dunno why this is needed
    wns = list(wns)  # Need list for .index method
    # Spectrum is often in highest to lowest, though don't want to assume
    if wns[0] > wns[1]:
        inverted = True  # i.e. highest to lowest wns
    else:
        inverted = False  # lowest to highest wns

    # Since sorting spikes we know position 0 is always the lowest
    # Remove spikes from array until none left
    while len(spikes_left) > 0:
        spikes_subtracted = []
        if inverted:
            min_wn = wns.index(spikes_left[0]) + width  # index of wn at lower end
            max_wn = wns.index(spikes_left[0]) - width  # index of wn at higher end
        else:
            min_wn = wns.index(spikes_left[0]) - width  # index of wn at lower end
            max_wn = wns.index(spikes_left[0]) + width  # index of wn at higher end

        spikes_subtracted.append(spikes_left[0])  # wn val of removed spike
        spikes_left = np.delete(spikes_left, 0)  # remove the value we used
        # In case there are other spikes in the range we have selected:
        while True:
            try:
                if spikes_left[0] <= wns[max_wn]:
                    if inverted:
                        max_wn = wns.index(spikes_left[0]) - width
                    else:
                        max_wn = wns.index(spikes_left[0]) + width
                    spikes_subtracted.append(spikes_left[0])
                    spikes_left = np.delete(spikes_left, 0)
                else:
                    break
            except IndexError:
                break  # If there are no more spikes left in array
        # Make x and y arrays for the polyfit
        if inverted:
            wns_range = wns[max_wn:min_wn+1]
            spec_range = spectrum_out[max_wn:min_wn+1]
        else:
            wns_range = wns[min_wn:max_wn+1]
            spec_range = spectrum_out[min_wn:max_wn+1]
        # Remove the spike values from the ranges
        for i in spikes_subtracted:
            wns_range = list(wns_range)  # need lists for .index
            spec_range = list(spec_range)
            index = wns_range.index(i)
            wns_range = np.delete(wns_range, index)
            spec_range = np.delete(spec_range, index)
        # Fixed degree can help cut out noise with larger windows.
        # Will cause errors for smaller window (e.g. width=1)
        # degree = len(wns_range) - 1
        degree = 3
        polynomial = np.polyfit(wns_range, spec_range, degree)
        coefficients = len(polynomial)
        for spike in spikes_subtracted:
            val = 0
            for i in range(coefficients):
                val = val + polynomial[-(i+1)]*(spike**i)
            spectrum_out[wns.index(spike)] = val
    
    return spectrum_out


def baselineALS(y, lam, p, niter=10):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w) # Do not create a new matrix, just update diagonal values
        Z = W + D
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z

def oval(h,w,y):
    return np.sqrt(((1-((y**2)/(h**2)))*(w**2)))

def rollingBall(spectrum,ball_H,ball_W):
    baseline = []
    oval_shape = []
    index = -ball_W
    for o in range(ball_W*2):
        oval_shape.append(oval(ball_W,ball_H,abs(index)))
        index += 1
    for x in range(len(spectrum)):
        ball = spectrum[x]-ball_H
        for y in range(ball_W*2):
            if (x-ball_W+y > 0) & (x-ball_W+y < len(spectrum)) & (index != 0):
                if spectrum[x-ball_W+y] < ball + oval_shape[y]:
                    ball = spectrum[x-ball_W+y] - oval_shape[y]
        baseline.append(ball)
    return np.array(baseline) + ball_H

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def airPLS(x, lambda_=30, porder=1, itermax=15):
    '''
    Adaptive iteratively reweighted penalized least squares for baseline fitting

    input
        x: input data (i.e. chromatogram of spectrum)
        lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background, z
        porder: adaptive iteratively reweighted penalized least squares for baseline fitting

    output
        the fitted background vector

    https://github.com/zmzhang/airPLS/blob/master/airPLS.py
    '''
    m=x.shape[0]
    w=np.ones(m)
    for i in range(1,itermax+1):
        z=WhittakerSmooth(x,w,lambda_, porder)
        d=x-z
        dssn=np.abs(d[d<0].sum())
        if(dssn<0.001*(abs(x)).sum() or i==itermax):
            if(i==itermax):
                print('WARNING max iteration reached!')
            break
        w[d>=0]=0 # d>0 means that this point is part of a peak, so its weight is set to 0 in order to ignore it
        w[d<0]=np.exp(i*np.abs(d[d<0])/dssn)
        w[0]=np.exp(i*(d[d<0]).max()/dssn)
        w[-1]=w[0]
    return z

def manualBaseline(array,mask,polyorder=20):
    spectra_copy = copy.deepcopy(array)
    spectra_copy[mask] = np.nan
    idx = np.isfinite(spectra_copy)
    ab = np.polyfit(np.arange(len(array))[idx], spectra_copy[idx], 20)
    p = np.polyval(ab,range(len(array)))
    return p

# def baselineCorrection(array, method = 'ALS', lam = 10**4, p = 0.01, niter = 10, fourier_values = 3,
#                        fourier_type = 'RDFT', polynomial = 3, iterations = 10, ball_H = 0.1,
#                        ball_W = 25, window_size = 100, lambda_airPLS = 30, porder = 1, itermax = 15,
#                        return_baseline = False, mask=False, manualPolyOrder=20):

def baselineCorrection(wns, array, method_dict):
    [[method, parameters]] = method_dict.items()
    parameters = list(parameters.values())
    return_baseline = False

    array = np.stack(array)
    baselined_array = np.zeros(np.shape(array))

    if method == 'ALS':
        lam, p, niter = parameters
        for index, spectra in enumerate(array):
            baselined_array[index,:] = spectra - baselineALS(spectra, lam=lam, p=p, niter=niter)

    elif method == 'airPLS':
        airPls, porder, itermax = parameters
        for index, spectra in enumerate(array):
            baselined_array[index,:] = spectra - airPLS(spectra, lambda_=airPls, porder=porder, itermax=itermax)

    elif method == 'FFT_DFT' or method == 'FFT_RDFT' or method == 'FFT_Hermitian':
        fourier_values = int(parameters[0])
        for index, spectra in enumerate(array):
            if method == 'FFT_DFT':
                rft = np.fft.fft(spectra)
                rft[:fourier_values] = 0
                baselined_array[index,:] = np.fft.ifft(rft)
            elif method == 'FFT_RDFT':
                rft = np.fft.rfft(spectra)
                rft[:fourier_values] = 0
                baselined_array[index,:] = np.fft.irfft(rft)
            elif method == 'FFT_Hermitian':
                rft = np.fft.hfft(spectra)
                rft[:fourier_values] = 0
                baselined_array[index,:] = np.fft.ihfft(rft)

    elif method == 'ModPoly':
        polynomial = int(parameters[0])
        for index, spectra in enumerate(array):
            baseObj=BaselineRemoval(spectra)
            baselined_array[index,:] = baseObj.ModPoly(polynomial)

    elif method == 'IModPoly':
        polynomial = int(parameters[0])
        for index, spectra in enumerate(array):
            baseObj=BaselineRemoval(spectra)
            baselined_array[index,:] = baseObj.IModPoly(polynomial)

    elif method == 'Zhang':
        polynomial = int(parameters[0])
        for index, spectra in enumerate(array):
            baseObj=BaselineRemoval(spectra)
            baselined_array[index,:] = baseObj.ZhangFit(polynomial)

    elif method == 'SWiMA':
        iterations = int(parameters[0])
        for index, spectra in enumerate(array):
            window = 3
            working_spectra = spectra
            for repeat in range(iterations):
                working_spectra = np.pad(working_spectra, 2, mode='reflect')
                smoothed_array = savgol_filter(working_spectra,window,0)
                a = working_spectra-smoothed_array
                a[a > 0] = 0
                working_spectra = a + smoothed_array
                window += 2
            baselined_array[index,:] = spectra - working_spectra[(window-3):-(window-3)]

    elif method == 'RollingBall':
        ball_H, ball_W = parameters
        for index, spectra in enumerate(array):
            baselined_array[index,:] = spectra - rollingBall(spectra,ball_H,ball_W)

    elif method == 'Average':
        window_size = int(parameters[0])
        for index, spectra in enumerate(array):
            f = np.pad(spectra, (int(window_size/2)), 'constant', constant_values=(spectra[0], spectra[-1]))
            baselined_array[index,:] = spectra - moving_average(f, int(window_size))[0:len(spectra)]

    elif method == 'manual':
        mask, manualPolyOrder = parameters
        for index, spectra in enumerate(array):
            baselined_array[index,:] = spectra - manualBaseline(spectra,mask,polyorder=manualPolyOrder)

    if return_baseline == True:
        return np.transpose(baselined_array), np.transpose(array-baselined_array)
    else:
        return baselined_array

# def normalise(array, axis = 1, method = 'max_within_range', normalisation_indexes = [890,910], wavenumbers=False,
#               zero_min=False, custom_values = False, return_normalisation_values = False):
    
def normalise(wns, array, method_dict):
    # Methods that don't require parameters first to prevent errors




    [[method, parameters]] = method_dict.items()
    normalisation_indexes = list(parameters.values())[0]
    
    if method == "Ip Normalisation":
        return normalise_ip(array, int(normalisation_indexes))

    # What is happening below here?? No idea, causes errors.
    # -----------------------------------------------------------------
    array = np.transpose(np.stack(array))
    axis = 1
    return_normalisation_values = False
    custom_values = False

    # Is this functionality needed?
    # if zero_min == True:
    #     array = array + abs(np.min(array))

    normalised_array = np.zeros(np.shape(array))

    # why did I make this like this??
    normalisation_indexes_2 = normalisation_indexes

    normalisation_indexes_2 = sorted(normalisation_indexes_2)
    normalisation_values = []

    if method == 'Scale':
        max_value = np.max(array)
        normalisation_values.append(max_value)
        normalised_array = array / max_value
    else:
        for spectra in range(np.shape(array)[axis]):
            if method == 'Maximum in range':
                max_value = max(array[normalisation_indexes_2[0]:normalisation_indexes_2[1],spectra])
                normalised_array[:,spectra] = array[:,spectra] / max_value
                normalisation_values.append(max_value)
            elif method == 'Maximum whole array':
                max_value = max(array[:,spectra])
                normalised_array[:,spectra] = array[:,spectra] / max_value
                normalisation_values.append(max_value)
            elif method == 'Whole array':
                max_value = sum(array[:,spectra])
                normalised_array[:,spectra] = array[:,spectra] / max_value
                normalisation_values.append(max_value)
            elif method == 'Single point':  
                normalised_array[:,spectra] = array[:,spectra] / array[normalisation_indexes_2[0],spectra]
                normalisation_values.append(array[normalisation_indexes_2[0],spectra])
            elif method == 'Area':
                areaunderstudy = array[:,spectra]
                max_value = sum(array[normalisation_indexes_2[0]:normalisation_indexes_2[1],spectra])
                normalised_array[:,spectra] = array[:,spectra] / max_value
                normalisation_values.append(max_value)
            elif method == 'Interpolate area':
                f = interpolate.interp1d(range((normalisation_indexes_2[1]-normalisation_indexes_2[0])),
                                         array[normalisation_indexes_2[0]:normalisation_indexes_2[1],spectra],
                                         kind='quadratic')
                normalised_array[:,spectra] = array[:,spectra] / sum(f(np.arange(0,
                                                                    (normalisation_indexes_2[1]-normalisation_indexes_2[0])-1,
                                                                     0.1)))
                normalisation_values.append(sum(f(np.arange(0,(normalisation_indexes_2[1]-normalisation_indexes_2[0])-1,0.1))))
            elif method == 'Max in interpolation range':
                f = interpolate.interp1d(range((normalisation_indexes_2[1]-normalisation_indexes_2[0])),
                                         array[normalisation_indexes_2[0]:normalisation_indexes_2[1],spectra],
                                         kind='quadratic')
                normalised_array[:,spectra] = array[:,spectra] / np.max(f(np.arange(0,
                                                                    (normalisation_indexes_2[1]-normalisation_indexes_2[0])-1,
                                                                     0.1)))
                normalisation_values.append(np.max(f(np.arange(0,(normalisation_indexes_2[1]-normalisation_indexes_2[0])-1,0.1))))

            elif method == 'custom_values':
                normalised_array[:,spectra] = array[:,spectra] / custom_values[spectra]
                normalisation_values.append(custom_values[spectra])

    if method == 'Area':
        max_value = np.max(np.mean(normalised_array[normalisation_indexes_2[0]:normalisation_indexes_2[1],:],axis=1))
        normalised_array = normalised_array / max_value
        normalisation_values = np.array(normalisation_values) * max_value
    elif method == 'Interpolate area':
        max_value = np.max(np.mean(normalised_array[normalisation_indexes_2[0]:normalisation_indexes_2[1],:],axis=1))
        normalised_array = normalised_array / max_value
        normalisation_values = np.array(normalisation_values) * max_value
    elif method == 'custom_values':
        max_value = np.max(np.mean(normalised_array,axis=1))
        normalised_array = normalised_array / max_value
        normalisation_values = np.array(normalisation_values) * max_value
    elif method == 'Whole array':
        max_value = np.max(np.mean(normalised_array,axis=1))
        normalised_array = normalised_array / max_value
        normalisation_values = np.array(normalisation_values) * max_value
    else:
        pass
    if return_normalisation_values == True:
        return normalised_array, np.array(normalisation_values)
    else:
        return np.transpose(normalised_array)
    
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

def normalise_ip(spectraldata, p):
    # spectraldata = np.transpose(spectraldata)
    normalised = []
    # Assuming "spectraldata" is an array of spectra, is this reasonable?
    for spectrum in spectraldata:
        thesum = 0
        for intensity in spectrum:
            thesum += intensity**p
        normalised.append(spectrum / (thesum**(1/p)))
    # Please let's change the setup so I don't have to keep transposing!
    # I don't understand why this keeps being necessary
    return normalised

def signal_to_noise(intensities_array, wns_array, wavenumber):
    '''Return the signal to noise of a given wavenumber.'''
    datapoints = []
    for count, array in enumerate(wns_array):
        # Find the index that is closest to given wn
        index = np.argmin(np.abs(array - wavenumber)) 
        # Store the intensity value in array
        datapoints.append(intensities_array[count][index])
    # Find average value and standard deviation of this point
    m = np.mean(datapoints)
    sd = np.std(datapoints)
    snr = m/sd
    return snr

def get_sample_names(filepaths, separator="%"):
   storage = []
   for filepath in filepaths:
       tail = filepath
       while True:
           head, sep, tail = tail.partition("\\")
           if len(sep) > 0:
               pass
           else:
               storage.append(head)
               break
        
   sample_names = []
   for sample_name in storage:
       head, tail, sep = sample_name.partition(separator)
       while len(head) > 0:
           try:
               head = float(head)
               sample_names.append(abs(head))  # Can't have -ve concs!
               break
           except ValueError:
               head = head[1:]
   return sample_names


def new_read_df_from_files(files):
    # initialise the top level lists
    x_arrays = []
    y_arrays = []
    
    # loop through the files
    for count, file_name in enumerate(files):
        # initialise the arrays to hold the data from each spectrum
        x_values = []
        y_values = []
        # Open and read the file
        with open(file_name, 'r') as file:
            # loop through the lines in the file
            for line in file:
                x, y = line.strip().split('\t')
                x_values.append(float(x))
                y_values.append(float(y))

        # append the spectral data
        x_arrays.append(x_values)
        y_arrays.append(y_values)

    # First values (highest wn) has anomalous value, delete
    x_arrays = np.delete(x_arrays, 0, 1)
    y_arrays = np.delete(y_arrays, 0, 1)

    # return the top level lists
    return np.array(x_arrays), np.array(y_arrays)

function_storage = {
    "Smooth": smooth,
    "Normalise": normalise,
    "Baseline Correction": baselineCorrection,
    "Remove Cosmic Ray Spikes": new_remove_crs,
}

def lorentzian(x, a, b, c, d):
    return d+c*((1/math.pi)*(b / ((x-a)**2 + b**2)))

def integrationfunction(intensities, wavenumbers, wn_width, popts=None, graphing=False):
    # intensities = np.transpose(intensities)
    # wavenumbers = np.transpose(wavenumbers)

    # find maximum point in the corrected arrays, take this as representative peak
    max = 0
    for i, spectrum in enumerate(intensities):
        for j, intensity in enumerate(spectrum):
            if intensity > max:
                max = intensity
                index_of_max = [i,j]

    wn_of_peak = wavenumbers[index_of_max[0]][index_of_max[1]]
    lowest_wn = wn_of_peak-(wn_width/2)
    highest_wn = wn_of_peak+(wn_width/2)

    # make lists of just the wn/intensity pairs that are within the
    # range of wn desired
    list_of_x = []
    list_of_y = []
    for i, spectrum in enumerate(wavenumbers):
        upper_list_x = []  # clearing/prepping a new list that can be appended
        upper_list_y = []
        for j, wn in enumerate(spectrum):
            if wn >= lowest_wn and wn <= highest_wn:
                upper_list_x.append(wn)
                if intensities[i][j] < 0:
                    upper_list_y.append(0)
                else:
                    upper_list_y.append(intensities[i][j])
        list_of_x.append(upper_list_x)
        list_of_y.append(upper_list_y)
       
    x = np.linspace(np.floor(lowest_wn), np.ceil(highest_wn), 100)
    integration = []

    # If not given optimal values to start with then the spectrum with
    # the clearest curve i.e. highest peak is used for estimates
    if popts is not None:
        a_guess = popts[0]
        b_guess = popts[1]
        c_guess = popts[2]
    else:
        a_guess = wn_of_peak
        b_guess = np.std(list_of_x[index_of_max[0]])
        c_guess = np.max(list_of_y[index_of_max[0]]) - np.min(list_of_y[index_of_max[0]])
    
    # For each spectrum we fit a lorentzian function
    for i, spectrum in enumerate(list_of_x):
        try:
            popt, pcov = curve_fit(lorentzian, list_of_x[i], list_of_y[i], p0=[a_guess, b_guess, c_guess, 0], maxfev=10000)

            # Can display the graphs of the data alongside the fitted curve
            # if graphing is set to true at function input
            if graphing: 
                ym = lorentzian(x, popt[0], popt[1], popt[2], popt[3])
                fig = mpl.figure()
                ax = fig.add_subplot(111)
                ax.scatter(list_of_x[i], list_of_y[i])
                ax.plot(x, ym, c='r', label='Best fit')
                ax.legend()
                mpl.show()

            # Find the area under the lornetz curve
            integral_area, certainty = scipy.integrate.quad(lorentzian, spectrum[0], spectrum[-1], args=(popt[0], popt[1], popt[2], popt[3]))
            integration.append(abs(integral_area))

        except RuntimeError:  # In case not finding optimal params
            print("Could not find optimal parameters for wn width ", wn_width, " and spectrum number ", i)
            integration.append(0)
        
    return integration, np.sqrt(np.diag(pcov)), popt


def wn_variation(intensities, wavenumbers, width_values):
    areas_list = []
    confidence_list = []

    # Use 50 wn_width to find initial optimal values
    area_under, fit_confidence, popts = integrationfunction(intensities, wavenumbers, 50)

    # For each value in the list of width_values perform the curve
    # fitting and integration
    for value in width_values:
        area_under, fit_confidence, popts = integrationfunction(intensities, wavenumbers, value, popts)
        areas_list.append(area_under)
        confidence_list.append(fit_confidence)
    return areas_list, confidence_list


def quickProcess(file_paths, sample_names, method_dict):
    wavenumbers, intensities = new_read_df_from_files(file_paths)
    data = {"Sample type": sample_names,
            "Raw data": list(intensities)}
    big_storage_df = pd.DataFrame(data) 

    for process in method_dict:
        # The process would take parameters from the passed dict
        # The functions inside function_storage don't currently follow this syntax
        if method_dict[process] is not False:
            intensities = function_storage[process](wavenumbers, intensities, method_dict[process])
            # store intensities in df 
            big_storage_df[process] = list(intensities)
        
    # integration_of_lorentzian, error = integrationfunction(intensities, wavenumbers, 40)
    
    return big_storage_df, wavenumbers, intensities


if __name__ == "__main__":
# ---------------------------------------------------------------------------------------------------------------------------------------
# For debugging purposes, to give the file somthing to work with. This will be done in the importing file.

    SpectraDir = 'EthanolData\\60water lens 10 3 1 0.3 0.1 0.03 0.01 60slit 300s 5um     below'
    SpectraDir = SpectraDir + '\*.txt'

    file_paths = glob.glob(SpectraDir, 
                    recursive = True)
    
    # file_paths = ['EthanolData\\60water lens 10 3 1 0.3 0.1 0.03 0.01 60slit 300s 5um below\\Hide\\60water lens-60slit-10%ETOH-30s 10 accumulations-2.txt']

    sample_names = get_sample_names(file_paths)

    normTestDict = {
        "Remove Cosmic Ray Spikes": False,
        "Smooth": {"Whittaker":{"Lambda":5, "d": 4}},
        "Normalise": {"Ip Normalisation":{"p": 5}},
        "Baseline Correction": {"airPLS": {"Lambda": 10, "Polynomial order": 1, "Maximum iterations": 50}}        
        }

    theDFofDFs, wns, intensities = quickProcess(file_paths, sample_names, normTestDict)



    # widthTest = np.linspace(9, 50, 43)

    # areas, confidence_of_curve = wn_variation(intensities, wns, widthTest)
    # polyfit_residuals = []

    # for integral in areas:
    #     try:
    #         coefficients, residiuals, rank, something, somethingelse = np.polyfit(concentrations, integral, 1, full=True)
    #     except TypeError:
    #         residiuals = 0
    #     polyfit_residuals.append(float(residiuals))

    # RSE = []
    # for sum_of_squares in polyfit_residuals:
    #     residual_standard_error = (sum_of_squares/(len(polyfit_residuals)-2))**(1/2)
    #     RSE.append(residual_standard_error)

    # # This new array can be plotted.
    # plt.rcParams['figure.figsize'] = [18,10]
    # font = {'family' : 'DejaVu Sans',
    #         'weight' : 'normal',
    #         'size'   : 24}
    # plt.rc('font', **font)

    # plt.axhline(y=0,color='0.2')
    # plt.scatter(concentrations, integration)
    # plt.scatter(widthTest, confidence_of_curve, marker="o", label="Curve confidence")
    # plt.scatter(widthTest, RSE, marker="v", label = "Best fit residual")

    # plt.plot(np.unique(concentrations), np.poly1d(coefficients)(np.unique(concentrations)))
    # print(residiuals)

    # plt.autoscale(enable=True, axis='x', tight=True)
    # plt.title("Calibration curve, 5-FU, TIRF, 100um")
    # plt.xlabel("Wavenumber width")
    # plt.ylabel("Best fit residual")

    # un/comment lines depending on whether standard or log plot is wanted:
    # plt.xlabel('Log10(100*Concentration(mM))')
    # plt.xlabel('Concentration (mM)')
    # plt.ylabel('Normalised intensity')
    # plt.legend(sample_names)

    # plt.show()

    # print(theDFofDFs)

    # plt.rcParams['figure.figsize'] = [18,10]
    # font = {'family' : 'DejaVu Sans',
    #         'weight' : 'normal',
    #         'size'   : 24}
    # plt.rc('font', **font)

    # plt.axhline(y=0,color='0.2')
    # plt.plot(wns, intensities, label=sample_names)
    # plt.autoscale(enable=True, axis='x', tight=True)
    # plt.title("graph title")
    # plt.xlabel("Wavenumbers (cm$^{-1}$)")
    # plt.ylabel("AU")
    # plt.legend()
    # plt.show()