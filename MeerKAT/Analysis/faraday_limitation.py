import numpy as np

def deltaLambdasquared(deltanu,nu_central):
    """
    Channel width in wavelength squared (units of m^2)
    Assuming a tophat channel bandpass (which is fine for deltanu<<nu_central).


    Inputs:
    deltanu    -- float -- channel width in Hertz (usually 8MHz)
    nu_central -- float -- channel central frequency in Hertz

    Returns:
    deltaLambdasquared -- float -- channel width in wavelength squared (m^2)

    """
    c = 299792458 # m / s
    return 2*c**2*deltanu/(nu_central**3)*(1+0.5*(deltanu/nu_central)**2)

def phimax(dlambdasq):
    """
    The maximum Faraday depth to which one has more than 50% sensitivity

    Given in rad/m^2 (if dlambdasq is given in m^2)
    """

    return np.sqrt(3)/dlambdasq

def DeltaLambdasquared(nu_min,nu_max):
    """
    Not to be confused with deltaLambdaSquared
    This is the total bandwidth in wavelength squared, i.e.,
    lambda_max^2 - lambda_min^2

    Inputs
    nu_min -- float -- minimum frequency in Hz
    nu_max -- float -- maximum frequency in Hz

    Returns
    DeltaLambdasquared -- float -- Total bandwidth in wavelength squared (m^2)

    """
    c = 299792458 # m / s

    lambda_max = c/nu_min
    lambda_min = c/nu_max

    return lambda_max**2 - lambda_min**2

print('resolution ', 2*np.sqrt(3)/DeltaLambdasquared(0.9e9, 1.6e9))
print('max scale ', np.pi/(299792458/1.6e9)**2)
print('phi max ', phimax(deltaLambdasquared(3e6, 1.2e9)))
