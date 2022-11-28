"""
Author: Joppe Swart
Description: Script to calculate the setjy coefficients manualy for polarization


"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import binom

#freq = np.array([628,665,700,730,765,800,836,   885,930,990,1045,1100,1150,1315,1365,1415,1470,1510,1630,1685,   1777,1835,2100])
freq = np.array([885.0,930.0,990.0,1045.0,1100.0,1150.0,1315.0,1365.0,1415.0,1470.0,1510.0,1630.0,1685.0])
reffreq = 1280.0
#polindex = np.array([3.21,3.71,4.47,5.10,5.77,6.36,6.90,   7.36,7.99,8.45,8.74,8.95,9.25,9.67,9.71,9.75,9.83,9.96,10.08,10.13])
polindex = np.array([7.36,7.99,8.45,8.74,8.95,9.25,9.67,9.71,9.75,9.83,9.96,10.08,10.13])
polangle = np.array([21.86,22.22,22.94,23.67,24.03,24.39,25.84,26.20,26.20,26.38,26.56,27.10,27.10])

def polindex_coff(freq, a, b, c, d):
    x = (freq-reffreq)
    return a + b*x + c*x**2 + d*x**3

def polang_coff(freq, a, b, c, d):
    x = (freq-reffreq)
    return a + b*x + c*x**2 + d*x**3
    
def logS(freq):
    a0 = 1.2480
    a1 = -0.4507
    a2 = -0.1798
    a3 = 0.0357
    return 10**(a0 + a1*np.log10(freq*1e-3) + a2*(np.log10(freq*1e-3))**2 + a3*(np.log10(freq*1e-3))**3)

def spix_coff(freq, a, b, c):
    S_0 = 15.77
    x = freq/reffreq
    return S_0*(x)**(a + b*np.log10(x)+ c*np.log10(x)**2)


spix, _ = curve_fit(spix_coff, freq, logS(freq))
pi, _ = curve_fit(polindex_coff, freq, polindex)
a, _ = curve_fit(polang_coff, freq, polangle)

print('pi = ', pi)
print('a = ', a)
print('spix = ',spix)

frequency = np.linspace(885.0, 1685.0, 200)

#print(frequency)
plt.figure()
plt.title('Polindex')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Polarization [%]')
plt.plot(freq, polindex, linestyle = '--', label='MeerKAT Memo')
plt.plot(frequency, polindex_coff(frequency, *pi), label='Curvefit Polindex')
plt.legend(loc='best')
plt.savefig(filename='polindex.png')
plt.close()

plt.figure()
plt.title('Polangle')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Polarization angle [deg]')
plt.plot(freq, polangle, linestyle = '--', label='MeerKAT Memo')
plt.plot(frequency, polang_coff(frequency, *a), label='Curvefit Polangle')
plt.legend(loc='best')
plt.savefig(filename='polangle.png')
plt.close()

plt.figure()
plt.title('Spectral index')
plt.xlabel('Flux density [Jy]')
plt.ylabel('Polarization angle [deg]')
plt.plot(frequency, logS(frequency), linestyle = '--', label='Perley Butler 2017')
x = freq/reffreq
plt.plot(frequency, spix_coff(frequency, *spix), alpha = 0.7, label='Curvefit Spix')
plt.legend(loc = 'best')
plt.savefig(filename='spix_fit.png')
plt.close()

