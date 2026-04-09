"""
start = time.time()

count = 0
    for gal_id,data in spectra_dict.items():
    w = data["wavelength"]
    f = data["flux"]  
    z = data["z"] 

end = time.perf_counter()

print("Average Time taken:", (end - start) / # of galaxies)
"""

import numpy as np
import time 

def lineflux(wave_rest, spec, cont_mask=None, spec_unc=None, line_mask=None, deg=1):

    """
    SUMMARY:
    Measuring the flux of a spectral line by fitting nearby cont


    VARIABLES AND DEFINTIONS:
    """

    wave_cont=wave_rest[cont_mask] #Restframe continuum wavelength 
    line_wave=wave_rest[line_mask] #Line wavelength within the line mask
    cont_spec=spec[cont_mask] #Spectrum within the continuum mask
    line_spec=spec[line_mask] #Line spectra within the line mask
    
    
    if spec_unc is not None:
        cont_unc=spec_unc[cont_mask] #uncertainty array in continuum windows
        coeffs_fit, cov=np.polyfit(wave_cont, cont_spec,deg, w=1/cont_unc, cov=True)#fitting the continuum with low noise points
        coeff_unc= np.sqrt(np.diag(cov)) #ask again
        cont_fit=np.polyval(coeffs_fit, line_wave) #continuum under the line
        V=np.vander(line_wave, deg+1) #matrix 
        cont_fit_var=np.sum(V@cov*V, axis=1) #continuum fit uncertainty at each wavelength
        
        #uncertainty in line region data
        line_unc = spec_unc[line_mask] #line uncertainty 
        line_flux = np.trapz(line_spec - cont_fit, x=line_wave) #measured flux above the continuum line
        total_var = line_unc**2 + cont_fit_var #total variance per wavelength point
        
#         x = line_wave
#         dx = np.diff(x)
#         weights = np.zeros_like(x)
#         weights[1:-1] = (dx[:-1] + dx[1:]) / 2
#         weights[0] = dx[0] / 2
#         weights[-1] = dx[-1] / 2
#         total_var = ((weights * line_unc)**2)
#         line_flux_unc = np.sqrt(np.sum(total_var + cont_fit_var))
#         print(weights.shape, line_unc.shape, line_flux.shape, total_var.shape, cont_fit_var.shape)
        line_flux_unc = np.sqrt(np.sum(total_var))
        return line_flux, line_flux_unc
        
    else:
        coeffs_fit=np.polyfit(wave_cont, cont_spec, deg=deg) #Coeffecients of continuum polynomial fit 
        cont_fit=np.polyval(coeffs_fit, line_wave) #Evaluating Continuum fit at line wavelength (line_wave)
        line_flux=np.trapz(line_spec - cont_fit, x=line_wave) #Integrating Line flux within the line mask 
        return line_flux
    
def equiv_width(wave_rest, spec, cont_mask, line_mask, spec_unc=None, deg=1):
    wave_cont = wave_rest[cont_mask]
    line_wave = wave_rest[line_mask]
    cont_spec = spec[cont_mask]
    line_spec = spec[line_mask]

    coeffs_fit = np.polyfit(wave_cont, cont_spec, deg=deg)
    cont_fit = np.polyval(coeffs_fit, line_wave)
    ew = np.trapz(1 - (line_spec / cont_fit), x=line_wave)

    if spec_unc is not None:
        line_unc = spec_unc[line_mask]
        ew_unc = np.sqrt(np.sum(((line_unc / cont_fit) ** 2)))
        return ew, ew_unc

    return ew

#equivelent width function 



def breakstrength(wave_rest, spec, blue_mask, red_mask, spec_unc=None, use_median=True): #defining what variables are needed

    """
    SUMMARY:
    Balmer and D4000 Breaks, measuring spectral break strength as a ratio of continuums. 
    Wavelength restframe arrays calculates the ratio flux of provided red and blue masks.
    Uncertainty of wavelength restframes calculates the ratio of propagated errors of provided red and blue masks.
    If the uncertainty of the error is provided, then both break strength and uncertainty are returned.  


    VARIABLES AND DEFINTIONS:
    Assumes provided variables falls within assumed #'s 
    Spec (array, length N): Spectrum corresponding to provided rest frame wavelength array 
    wave_rest (array): Restframe wavelength array. Assumes that all red/blue masks defines fall within wavelength array 
    """
    
    #wavlengths instead of mask, blue_wave = x red and blue = y, etc etc 
    
    blue_wave = wave_rest[blue_mask] #restframe wavelength in blue window
    red_wave=wave_rest[red_mask] #restframe wavelength in red window
    blue_spec=spec[blue_mask] #flux in blue window
    red_spec=spec[red_mask] #flux in red window 
    #if use median = true calc median if else caluclate mean
    #use averages instead 
    blue_med = np.median(blue_spec) #medians of blue continuum level
    red_med = np.median(red_spec) #medians of red continuum level
    
    break_strength = red_med / blue_med #spectral break stregnth ratio of red to blue continuum
    
    if spec_unc is not None: #uncertainty / uncertain data within the spectrum?
        blue_unc = spec_unc[blue_mask] #uncertainty in blue window
        red_unc = spec_unc[red_mask] #uncertainty in red window
        
        blue_err = 1.253 * np.sqrt(np.mean(blue_unc**2)) / np.sqrt(len(blue_unc)) #error value in blue window
        red_err = 1.253 * np.sqrt(np.mean(red_unc**2)) / np.sqrt(len(red_unc)) #error vlaue in red window
        
        break_unc = break_strength * np.sqrt((red_err / red_med)**2 + (blue_err / blue_med)**2) #Uncertainty of flux within each window 
        
        return break_strength, break_unc
    
    else:  # if no uncertainty 
        return break_strength

    
"""
Summary:
Achieving rest-frame colors from user provided spectrum. 
Using individual filters to determine how bright a spectrum would appear if measured in a given wavelength range. 

Variables:
Filter = sedpy instance of filters that has been loaded in
wave_rest = restframe wavelength in angstroms
spec = spectrum in F_lam

colors = {} defining colors 
"""

def restframe_color(wave_rest, spec, filters, spec_unc=None):
    
    assert len(filters) ==2, f"Only two filters are allowed, recived {len(filter)}." #make sure statement is true , lenght of filters assuming list must be two 
    #f string attacthed 
    
    mags = [] #mags in a list 
    mags_err = []
    
    for filt in filters:
        mag = filt.ab_mag(wave_rest, spec)
        mags.append(mag)
        if spec_unc is not None:
            weights = filt.wavelength * filt.transmission * filt.dwave
            weights = np.interp(wave_rest,filt.wavelength,weights) # interpolate onto the observed wavelength grid
            mag_err = np.sqrt(np.sum((weights * spec_err)**2)) # flux_err = spec_err
            mags_err.append(mag_err) 
        
    color = mags[0] - mags[1]
    
    if spec_unc is not None:
        color_err = np.sqrt(mags_err[0]*2 + mags_err[1]**2)
        return color, color_err
    
    else:
        return color 

  



