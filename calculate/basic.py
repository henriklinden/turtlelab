from __future__ import division
from scipy import signal
import numpy as np
from pylab import rms_flat

def nerve_kernel(nerve, sr):
    bandwidth = 0.100 #sec
    #nerve_rms = sqrt(mean((nerve.*kernel').^2));
    
    res = np.convolve(np.absolute(nerve),kernel(bandwidth,sr),'same')
    #subtract noise level
    #estimate of noise level: lowest point in all data - careful with edges
    noise = np.min(res[10*(bandwidth*sr):-10*(bandwidth*sr)])
    return res-noise

def nerve_top(sig, sr):
    mx = signal.argrelmax(sig,order=int(np.ceil(sr*0.20)))[0] #100 ms on each side
    idx = sig[mx]>0.0007
    return mx[idx]
    
def nerve_cycle_width(sig, maxs):
    res = np.zeros(maxs.shape)
    for i,mx in enumerate(maxs):
        h_height = sig[mx]/2
        idx = np.argwhere(sig<h_height)
        after = idx[idx>mx][0]
        before = idx[idx<mx][-1]
        res[i] = after-before
    return res
    
def spike_fr(intra, sr, start, stop):
    #calculate underlying firing rate from a spike train
    bandwidth = 0.100 #sec
    length = np.ceil((stop-start)*sr)
    
    hist,bin_edges = np.histogram(intra, bins=length, range=(start,stop))
    
    kern = kernel(bandwidth, sr)*sr
    return signal.fftconvolve(hist,kern,'same')
    
def kernel(bandwidth, sr):
    #kernel_t = -5*bandwidth*sr:1:5*bandwidth*sr;
    #kernel = (1/(sqrt(2*pi)*bandwidth*sr))*exp(-(kernel_t.^2)/(2*(bandwidth*sr).^2));
    kernel_t = np.arange(-5*bandwidth*sr,5*bandwidth*sr+1,1)
    kernel = (1/(np.sqrt(2*np.pi)*bandwidth*sr))*np.exp(-(np.square(kernel_t))/(2*np.square(bandwidth*sr)))
    return kernel
        
def filter_nerve(nerve, sr):
    N = 3
    Fc = 10
    b, a = signal.butter(N, 2*Fc/sr, 'highpass', output='ba')
    return signal.filtfilt(b,a, nerve-np.mean(nerve))

def isi_unit(spiketrain):
    isi = np.diff(spiketrain)
    hist, bin_edges = np.histogram(isi,50)
    return hist, bin_edges

def xcorr(sig1, sig2, maxlags=None, norm='none'):
    #cross correlation
    #TODO: Tapering!!!
    a_len = sig1.shape[0]
    b_len = sig2.shape[0]
    
    a = sig1-np.mean(sig1)
    b = sig2-np.mean(sig2)
    if a_len>b_len:
        b = np.zeros(a_len)
        tmp = np.floor((a_len-b_len)/2)
        b[tmp:tmp+b_len] = sig2-np.mean(sig2)
    elif b_len>a_len:
        a = np.zeros(b_len)
        tmp = np.floor((b_len-a_len)/2)
        a[tmp:tmp+a_len] = sig1-np.mean(sig1)
    N = len(a)

    #Array flipped convolution: correlation.
    res = signal.fftconvolve(a,b[::-1], mode='full')
    if maxlags == None:
        maxlags = N-1
        lags = np.arange(0, 2*N-1)
    else:
        assert maxlags < N
        lags = np.arange(N-maxlags-1, N+maxlags)
    #print N
    #print np.shape(res)
    if norm == 'biased':
        res = res[lags] / float(N)    # do not use /= !!
    elif norm == 'unbiased':
        res = res[lags] / (float(N)-abs(np.arange(-N+1, N)))[lags]
    elif norm == 'coeff':
        rms = np.sqrt(np.mean(np.absolute(sig1)**2)) * np.sqrt(np.mean(np.absolute(sig2)**2))
        res = res[lags] / rms / float(N)
    else:
        res = res[lags]
    lags = np.arange(-maxlags, maxlags+1)
    return res,lags

def spike_correlogram(spk_tr1, spk_tr2, window=0.005, sr=40, maxlags=0.1, normalize=False):
    """
    :param ref: Spike times of the reference train. It is 
        assumed that the spike times are ordered sequentially.
    
    :param comp: Spike times of the comparing train that shifts in time.
        It is assumed that the spike times are ordered sequentially.
    
    :param window: Time window to say that two spikes are synchronized.
        This has a default value of 5.
        
    :param dt: The binwidth.
    
    :param timewindow: Correlogram range between [-timewindow, timewindow].
        
    :param normalize: If True, values are normalized by the rate expected by
        chance at each lag, where chance is defined as 
        ``2 * frequency(comp) * window * len(ref)``.
        
    :return: A vector of length ``1 + 2*timewindow/dt``.
    """
    lags = np.arange(0,(2*maxlags)+window,window)
    lags = lags-lags[-1]/2
    res = np.zeros(lags.shape)
    bin_edges = np.append((lags-(window/2)),(lags[-1]+(window/2)))
    #bin_edges = bin_edges*sr
    #
    for idx1 in spk_tr1:
        tmp = spk_tr2-idx1
        hist, ttt = np.histogram(tmp, bins=bin_edges)
        res += hist
    #
    if normalize:
        res = res/np.sqrt(len(spk_tr1)*len(spk_tr2))
    else:
        res = res
    #        
    return res, lags

def spike_correlation(spk_tr1, spk_tr2, window=0.005, normalize=False):
    """
    :param ref: Spike times of the reference train. It is 
        assumed that the spike times are ordered sequentially.
    
    :param comp: Spike times of the comparing train that shifts in time.
        It is assumed that the spike times are ordered sequentially.
    
    :param window: Time window to say that two spikes are synchronized.
        This has a default value of 5.
        
    :param dt: The binwidth.
    
    :param timewindow: Correlogram range between [-timewindow, timewindow].
        
    :param normalize: If True, values are normalized by the rate expected by
        chance at each lag, where chance is defined as 
        ``2 * frequency(comp) * window * len(ref)``.
        
    :return: A vector of length ``1 + 2*timewindow/dt``.
    """
    res = 0
    for idx1 in spk_tr1:
        tmp = spk_tr2-idx1
        #print np.sum([(tmp>window)&(tmp<-window)])
        res += np.sum([(tmp>-window)&(tmp<window)])
        
    #
    if normalize:
        res = res/np.sqrt(len(spk_tr1)*len(spk_tr2))
    else:
        res = res
    #        
    return res