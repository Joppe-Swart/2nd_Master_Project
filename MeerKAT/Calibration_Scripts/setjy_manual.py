import numpy as np
def casa_flux_model(lnunu0, iref, *args):
    """
    Compute model:
    iref * 10**lnunu0 ** (args[0] + args[1] * lnunu0 + args[1] * lnunu0 ** 2 + args[0] * lnunu0 ** 3)
    """
    exponent = np.sum([arg * (lnunu0 ** (power  ))
                       for power, arg in enumerate(args)], axis=0)
    return iref * (10**lnunu0) **(exponent)
    
def fit_flux_model(nu, s, nu0, sigma, sref, order=5):  
    from scipy.optimize import curve_fit
    from scipy.special import binom
    """
    Fit a flux model of given order from :
    S = fluxdensity *(freq/reffreq)**(spix[0]+spix[1]*log(freq/reffreq)+..)
    Very rarely, the requested fit fails, in which case fall
    back to a lower order, iterating until zeroth order. If all
    else fails return the weighted mean of the components.
    Finally convert the fitted parameters to a
    katpoint FluxDensityModel:
    log10(S) = a + b*log10(nu) + c*log10(nu)**2 + ...
    Parameters
    ----------
    nu : np.ndarray
        Frequencies to fit in Hz
    s : np.ndarray
        Flux densities to fit in Jy
    nu0 : float
        Reference frequency in Hz
    sigma : np.ndarray
        Errors of s
    sref : float
        Initial guess for the value of s at nu0
    order : int (optional)
        The desired order of the fitted flux model (1: SI, 2: SI + Curvature ...)
    """

    init = [sref, -0.7] + [0] * (order - 1)
    lnunu0 = np.log10(nu/nu0)
    for fitorder in range(order, -1, -1):
        try:
            popt, _ = curve_fit(casa_flux_model, lnunu0, s, p0=init[:fitorder + 1], sigma=sigma)
        except RuntimeError:
            log.warn("Fitting flux model of order %d to CC failed. Trying lower order fit." %
                     (fitorder,))
        else:
            coeffs = np.pad(popt, ((0, order - fitorder),), "constant")
            return [nu0] +  coeffs.tolist()
    # Give up and return the weighted mean
    coeffs = [np.average(s, weights=1./(sigma**2))] + [0] * order
    return [nu0]+  coeffs.tolist()

def convert_flux_model(nu=np.linspace(0.9,2,200)*1e9 , a=1,b=0,c=0,d=0,Reffreq= 1.0e9) :
    """
    Convert a flux model from the form:
    log10(S) = a + b*log10(nu) + c*log10(nu)**2 + ...
    to an ASA style flux model in the form:
    S = fluxdensity *(freq/reffreq)**(spix[0]+spix[1]*log(freq/reffreq)+..)
    Parameters
    ----------
    nu : np.ndarray
        Frequencies to fit in Hz
    a,b,c,d : float
        parameters of a log flux model.
    Reffreq : float
        Reference frequency in Hz
    returns :
    reffreq,fluxdensity,spix[0],spix[1],spix[2]
    """
    MHz = 1e6
    S = 10**(a + b*np.log10(nu/MHz) +c*np.log10(nu/MHz)**2 + d*np.log10(nu/MHz)**3)
    return fit_flux_model(nu, S , Reffreq,np.ones_like(nu),sref=1 ,order=3)

#name=0408-65 epoch=2016 ra=04h08m20.4s dec=-65d45m09s a=-0.9790 b=3.3662 c=-1.1216  
a=-0.9790
b=3.3662
c=-1.1216
d=0.0861
reffreq,fluxdensity,spix0,spix1,spix2 =  convert_flux_model(np.linspace(0.9,2,200)*1e9,a,b,c,d)
f_cal_alt = '0408-658'

setjy(vis='Bullet_CLuster_Hanningsmooth.ms',
      field=f_cal_alt,
      spix=[spix0, spix1, spix2, 0],
      fluxdensity = fluxdensity,
      reffreq='%f Hz'%(reffreq),
      standard='manual')
