"""
General zero-phase filter.  Thanks to Ed Gross for most
of this.  It's mostly a python translation of matlab code he
provided. (this is a generalization of the code in tidal_filter.py)
"""

from numpy import *
from scipy.signal.filter_design import butter

from filtfilt import filtfilt
# seems to have issues:
# from scipy.signal import filtfilt

def lowpass(data,in_t,cutoff,order=4):
    """
    data: vector of data
    in_t: sample times
    cutoff: cutoff period in the same units as in_t

    returns vector same as data, but with high frequencies removed
    """
    
    # Step 1: Determine dt from data (complicated due to TRIM time series precision legacy) 
    mean_dt=mean(diff(in_t))

    # print 'Data Time Step: %8.4f'%mean_dt

    # test size of data for vectors only

    Wn = mean_dt / cutoff  # 90 = half # of minutes in 3 hours.

    # print 'Using %1ith Order Butterworth Filter, with %g time-unit cutoff (Wn=%12.6f).'%(order,2*cutoff,Wn)

    B,A = butter(order, Wn)
    
    data_filtered = filtfilt(B,A,data)

    return data_filtered

def lowpass_gotin(data,in_t_days,*args,**kwargs):
    """ Approximate Gotin's tidal filter
    Note that in preserving the length of the dataset, the ends aren't really
    valid
    """
    mean_dt_h = 24*mean(diff(in_t_days))

    # how many samples are in 24 hours?
    N24 = round(24. / mean_dt_h)
    # and in 25 hours?
    N25 = round(25. / mean_dt_h)

    A24 = ones(N24) / float(N24)
    A25 = ones(N25) / float(N25)

    data = convolve(data,A24,'same')
    data = convolve(data,A24,'same')
    data = convolve(data,A25,'same')

    return data
    


# def bode(G):
#     figure()
#     w=arange(1e-4j, 1e-1j, 1e-6j)
#     y = polyval(G.num, w) / polyval(G.den, w)
#     mag = 20.0*log10(abs(y))
#     pha = arctan2(y.getimag(), y.getreal())*180.0/pi
#     for i in arange(1, pha.shape[0]):
#         if abs(pha[i]-pha[i-1]) > 179:
#             pha[i:] -= 360.
#  
#     subplot(211)
#     semilogx(w.imag, mag)
#     grid()
#     gca().xaxis.grid(True, which='minor')
#     
#     gcf().text(0.5, 0.91, r'Bode diagram: %s' %G.desc, 
#                horizontalalignment='center', fontsize=16)    
#     ylabel(r'Magnitude (db)')
#     
#     subplot(212)
#     semilogx(w.imag, pha)
#     grid()
#     gca().xaxis.grid(True, which='minor')
#     ylabel(r'Phase (deg)')
#     yticks(arange(0, pha.min()-30, -30))
#  
#     gcf().canvas.mpl_connect('button_press_event', line_properties)
#  
#     return mag, pha
