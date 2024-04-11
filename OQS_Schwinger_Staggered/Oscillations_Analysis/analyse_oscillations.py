import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

def get_peaks_of_peaks(peaks):
    
    indices = []
    peaks_of_peaks = []
    upward = peaks[1] > peaks[0]
    
    for (idx, element) in enumerate(peaks):
        
        if idx == 0:
            continue
        
        if upward and (element < peaks[idx-1]):
            
            indices.append(idx-1)
            peaks_of_peaks.append(peaks[idx-1])
            upward = False
            
        elif element > peaks[idx-1]:
            
            upward = True

    return indices, peaks_of_peaks

def exponential_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

def exponential_decay_two(x, a, b):
    return a * np.exp(-b * x)

def inverse_linear(x, a, b, c):
    return a*(b/x) + c

plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.figsize'] = (9, 7)

N = 6
e = 1
x = 1/(e^2)
ma = 0.1
l_0 = 10
env_corr_type = "delta"
sigma_over_a = 3.0
aD_0 = 0.15
# beta = 0.1
# aT = 1/beta
dt = 0.04 
steps = 3000
aD_0 = 0.1

aT_list = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

b_list = []

for aT in aT_list:

    f = h5py.File(f'at_vs_pnd_{aT}.h5', 'r')

    at, pnd = f['at_list'][:], f['pnd'][:]
    
    coefficients = np.polyfit(at, pnd, deg = 1)
    slope, intercept = coefficients
    
    ax1.plot(at, at*slope + intercept, color = 'k', zorder = 10)
    
    ax1.plot(at, pnd, label = f'aT = {aT}')

    peaks_indices, _ = find_peaks(pnd)

    at_peaks = at[peaks_indices]
    pnd_peaks = pnd[peaks_indices]

    peaks_of_peaks_indices, pnd_peaks_of_peaks = get_peaks_of_peaks(pnd_peaks)
    at_peaks_of_peaks = at_peaks[peaks_of_peaks_indices]
                      
    ax1.scatter(at_peaks_of_peaks, pnd_peaks_of_peaks, marker = 'x', color = 'red')
    ax2.scatter(at_peaks_of_peaks, pnd_peaks_of_peaks - (at_peaks_of_peaks*slope + intercept), marker = 'x', label = f'aT = {aT}')
    
    popt, pcov = curve_fit(exponential_decay_two, at_peaks_of_peaks, pnd_peaks_of_peaks - (at_peaks_of_peaks*slope + intercept))
    a_opt, b_opt = popt
    ax2.plot(at_peaks_of_peaks, exponential_decay_two(at_peaks_of_peaks, *popt), label = f'Fit: {b_opt}')
    b_list.append(b_opt)
    
popt, pcov = curve_fit(exponential_decay, aT_list, b_list)
a_opt, b_opt, c_opt = popt
ax3.scatter(aT_list, b_list)
ax3.plot(aT_list, exponential_decay(np.array(aT_list), *popt), label = f'Fit: {b_opt}')

ax1.legend()
ax2.legend()
ax3.legend()
ax1.set_title(f'N={N},ma={ma},l_0={l_0},aD_0={aD_0},tau={dt},e={e}')
ax2.set_title(f'N={N},ma={ma},l_0={l_0},aD_0={aD_0},tau={dt},e={e}')
plt.show()

# fft_spectrum = np.fft.rfft(pnd)
# freq = np.fft.rfftfreq(len(pnd), d=1./(at[1]-at[0]))

# fft_spectrum_abs = np.abs(fft_spectrum)

# plt.plot(freq, fft_spectrum_abs)
# plt.xlabel("frequency, Hz")
# plt.ylabel("Amplitude, units")
# plt.show()
