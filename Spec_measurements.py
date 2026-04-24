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

from astropy.modeling import models, fitting
from scipy.interpolate import UnivariateSpline
import numpy as np
import matplotlib.pyplot as plt


def fit_line_emission_absorption(
    wave, flux, flux_err=None,
    central=4861, outer_window=30, inner_window=5, middle_window=None,
    cont_type='spline', poly_degree=1, center_jitter=3,
    n_sigma_int=5, int_type='gaussian', plot=False):
    """
    Fit a spectral line with combined emission and absorption components and
    compute line flux and equivalent width (EW) with uncertainties.

    This function models a spectral feature using a combination of:
    (1) a smooth continuum (spline or polynomial), and
    (2) two Gaussian components representing emission and absorption.
    It then integrates the resulting profile to estimate the total line flux
    and equivalent width, including multiple sources of uncertainty.

    Parameters
    ----------
    wave : array-like
        Wavelength array.
    flux : array-like
        Flux array corresponding to `wave`.
    flux_err : array-like, optional
        Flux uncertainties. If None, zeros are assumed (no observational errors).
    central : float, optional
       Central wavelength of the spectral line (default: 4861 Å for Hbeta).
    outer_window : float, optional
        Half-width (in wavelength units) of the region used for continuum fitting.
    inner_window : float, optional
        Half-width around `central` used for fitting the Gaussian components.
    middle_window : float, optional
        Exclusion region around the line center for continuum fitting.
        Defaults to `inner_window` if not provided.
    cont_type : {'spline', 'poly'}, optional
        Type of continuum model:
        - 'spline': Univariate spline fit
        - 'poly': Polynomial fit of degree `poly_degree`
    poly_degree : int, optional
        Degree of polynomial used if `cont_type='poly'`.
    center_jitter : float, optional
        Allowed deviation of Gaussian centers from `central`.
    n_sigma_int : float, optional
        Number of effective sigma used to define the integration window.
    int_type : {'native', 'gaussian'}, optional
        Integration method:
        - 'native': integrate directly on observed data grid
        - 'gaussian': integrate model on fine grid
    plot : bool, optional
        If True, plot the data, continuum, and fitted model.

    Returns
    -------
    line_flux : float
        Integrated flux of the line.
    line_flux_err : float
        Total uncertainty on line flux, combining observational,
        continuum, and model uncertainties.
    EW : float
        Equivalent width of the line.
    EW_err : float
        Uncertainty on the equivalent width.

    Notes
    -----
    - The emission component is constrained to have positive amplitude,
      while the absorption component is constrained to be negative.
    - The integration window is determined dynamically using an
      area-weighted combination of the fitted Gaussian widths.
    - Uncertainties include:
        * observational noise (from `flux_err`)
        * continuum modeling uncertainty (RMS of residuals)
        * Gaussian parameter covariance (if using model integration)
    - The equivalent width is defined as:
        EW = ∫ (1 - F_line / F_cont) dλ

    Raises
    ------
    ValueError
        If `cont_type` or `int_type` is not recognized.
    """

    if middle_window is None:
        middle_window = inner_window

    # -------------------------------
    # Region selection
    # -------------------------------
    mask_region = (wave >= central - outer_window) & (wave <= central + outer_window)
    w_reg = wave[mask_region]
    f_reg = flux[mask_region]
    f_err_reg = flux_err[mask_region] if flux_err is not None else np.zeros_like(f_reg)

    mask_line = (wave >= central - inner_window) & (wave <= central + inner_window)
    w_fit = wave[mask_line]
    f_fit = flux[mask_line]
    f_err_fit = flux_err[mask_line]

    # -------------------------------
    # Continuum fit
    # -------------------------------
    cont_mask = (w_reg < central - middle_window) | (w_reg > central + middle_window)
    w_cont, f_cont, f_cont_err = w_reg[cont_mask], f_reg[cont_mask], f_err_reg[cont_mask]

    if cont_type == 'spline':
        s_factor = len(w_cont) * np.median(np.diff(f_cont))**2
        cont_model = UnivariateSpline(w_cont, f_cont, s=s_factor)
    elif cont_type == 'poly':
        poly = models.Polynomial1D(degree=poly_degree)
        poly.c1.fixed = True
        poly.c1.value = 0.0
        cont_model = fitting.LinearLSQFitter(calc_uncertainties=True)(poly, w_cont, f_cont,weights=1/f_cont_err)
    else:
        raise ValueError("cont_type must be 'spline' or 'poly'")

    cont_vals = cont_model(w_reg)
    cont_rms = np.std(f_cont - cont_model(w_cont))

    # -------------------------------
    # Fit Gaussians (continuum-subtracted)
    # -------------------------------
    flux_sub = f_fit - cont_model(w_fit)

    g_em = models.Gaussian1D(
        amplitude=np.max(flux_sub), mean=central, stddev=2.0
    )
    g_em.amplitude.bounds = (0, 1e6)
    g_em.mean.bounds = (central - center_jitter, central + center_jitter)
    g_em.stddev.bounds = (0.1, 5)

    g_abs = models.Gaussian1D(
        amplitude=np.min(flux_sub), mean=central, stddev=5.0
    )
    g_abs.amplitude.bounds = (-1e6, 0)
    g_abs.mean.bounds = (central - center_jitter, central + center_jitter)
    g_abs.stddev.bounds = (5, 15)

    fitter = fitting.LevMarLSQFitter(calc_uncertainties=True)
    g_fit = fitter(g_em + g_abs, w_fit, flux_sub, maxiter=50000,weights=1/f_err_fit)

    # -------------------------------
    # Extract parameters
    # -------------------------------
    amp_em, mean_em, std_em = g_fit[0].parameters
    amp_abs, mean_abs, std_abs = g_fit[1].parameters

    # -------------------------------
    # Area-weighted window
    # -------------------------------
    w_em = np.abs(amp_em * std_em)
    w_abs = np.abs(amp_abs * std_abs)
    denom = w_em + w_abs

    if denom > 0:
        sigma_eff = (w_em * std_em + w_abs * std_abs) / denom
        central_eff = (w_em * mean_em + w_abs * mean_abs) / denom
    else:
        sigma_eff, central_eff = std_em, mean_em

    sigma_window = n_sigma_int * sigma_eff
    w_min, w_max = central_eff - sigma_window, central_eff + sigma_window

    # -------------------------------
    # Native grid (for uncertainties ALWAYS)
    # -------------------------------
    line_mask = (w_reg >= w_min) & (w_reg <= w_max)
    w_native = w_reg[line_mask]
    cont_native = cont_vals[line_mask]
    f_native = f_reg[line_mask]
    f_err_native = f_err_reg[line_mask]

    delta_w = np.gradient(w_native)

    # -------------------------------
    # Choose evaluation grid
    # -------------------------------
    if int_type == 'native':
        w_eval = w_native
        cont_eval = cont_native
        f_eval = f_native

    elif int_type == 'gaussian':
        w_eval = np.linspace(w_min, w_max, 1000)
        cont_eval = np.interp(w_eval, w_reg, cont_vals)
        f_eval = cont_eval + g_fit(w_eval)

    else:
        raise ValueError("int_type must be 'native' or 'gaussian'")

    # -------------------------------
    # Flux + EW
    # -------------------------------
    line_flux = np.trapz(f_eval, x=w_eval)
    EW = np.trapz(1 - (f_eval / cont_eval), x=w_eval)

    # -------------------------------
    # Observational + continuum errors
    # -------------------------------
    line_flux_err_obs = np.sqrt(np.sum((f_err_native * delta_w)**2))
    line_flux_err_cont = np.sqrt(np.sum((cont_rms * delta_w)**2))

    EW_err_obs = np.sqrt(
        np.sum((f_err_native / cont_native * delta_w)**2) +
        np.sum(((f_native - cont_native) * cont_rms / cont_native**2 * delta_w)**2)
    )

    # -------------------------------
    # Gaussian model uncertainty (analytic)
    # -------------------------------
    line_flux_err_model = 0.0
    EW_err_model = 0.0

    cov = fitter.fit_info['param_cov']

    if int_type == 'gaussian' and cov is not None:

        sqrt2pi = np.sqrt(2 * np.pi)
        J = np.zeros(len(g_fit.parameters))

        # parameter order: [A_em, mean_em, std_em, A_abs, mean_abs, std_abs]
        J[0] = std_em * sqrt2pi
        J[2] = amp_em * sqrt2pi
        J[3] = std_abs * sqrt2pi
        J[5] = amp_abs * sqrt2pi

        line_flux_err_model = np.sqrt(J @ cov @ J)

        cont_mean = np.mean(cont_native)
        EW_err_model = line_flux_err_model / cont_mean

    # -------------------------------
    # Final uncertainties
    # -------------------------------
    line_flux_err = np.sqrt(
        line_flux_err_obs**2 +
        line_flux_err_cont**2 +
        line_flux_err_model**2
    )

    EW_err = np.sqrt(EW_err_obs**2 + EW_err_model**2)

    # -------------------------------
    # Plot
    # -------------------------------
    if plot:
        x = np.linspace(np.min(w_reg), np.max(w_reg), 1000)
        plt.figure(figsize=(8, 5))
        plt.step(w_reg, f_reg, where='mid', color='k', label='Data')
        plt.plot(w_reg, cont_vals, '--b', label='Continuum')
        plt.plot(x, cont_model(x) + g_fit(x), '--g', label='Fit')
        plt.axvspan(w_min, w_max, color='orange', alpha=0.2)
        plt.legend()
        plt.xlabel("Wavelength")
        plt.ylabel("Flux")
        plt.show()

    return line_flux, line_flux_err, EW, EW_err

# def lineflux(wave_rest, spec, cont_mask=None, spec_unc=None, line_mask=None, deg=1):

#     """
#     SUMMARY:
#     Measuring the flux of a spectral line by fitting nearby cont


#     VARIABLES AND DEFINTIONS:
#     """

#     wave_cont=wave_rest[cont_mask] #Restframe continuum wavelength 
#     line_wave=wave_rest[line_mask] #Line wavelength within the line mask
#     cont_spec=spec[cont_mask] #Spectrum within the continuum mask
#     line_spec=spec[line_mask] #Line spectra within the line mask
    
    
#     if spec_unc is not None:
#         cont_unc=spec_unc[cont_mask] #uncertainty array in continuum windows
#         coeffs_fit, cov=np.polyfit(wave_cont, cont_spec,deg, w=1/cont_unc, cov=True)#fitting the continuum with low noise points
#         coeff_unc= np.sqrt(np.diag(cov)) #ask again
#         cont_fit=np.polyval(coeffs_fit, line_wave) #continuum under the line
#         V=np.vander(line_wave, deg+1) #matrix 
#         cont_fit_var=np.sum(V@cov*V, axis=1) #continuum fit uncertainty at each wavelength
        
#         #uncertainty in line region data
#         line_unc = spec_unc[line_mask] #line uncertainty 
#         line_flux = np.trapz(line_spec - cont_fit, x=line_wave) #measured flux above the continuum line
#         total_var = line_unc**2 + cont_fit_var #total variance per wavelength point
        
# #         x = line_wave
# #         dx = np.diff(x)
# #         weights = np.zeros_like(x)
# #         weights[1:-1] = (dx[:-1] + dx[1:]) / 2
# #         weights[0] = dx[0] / 2
# #         weights[-1] = dx[-1] / 2
# #         total_var = ((weights * line_unc)**2)
# #         line_flux_unc = np.sqrt(np.sum(total_var + cont_fit_var))
# #         print(weights.shape, line_unc.shape, line_flux.shape, total_var.shape, cont_fit_var.shape)
#         line_flux_unc = np.sqrt(np.sum(total_var))
#         return line_flux, line_flux_unc
        
#     else:
#         coeffs_fit=np.polyfit(wave_cont, cont_spec, deg=deg) #Coeffecients of continuum polynomial fit 
#         cont_fit=np.polyval(coeffs_fit, line_wave) #Evaluating Continuum fit at line wavelength (line_wave)
#         line_flux=np.trapz(line_spec - cont_fit, x=line_wave) #Integrating Line flux within the line mask 
#         return line_flux
    
# def equiv_width(wave_rest, spec, cont_mask, line_mask, spec_unc=None, deg=1):
#     wave_cont = wave_rest[cont_mask]
#     line_wave = wave_rest[line_mask]
#     cont_spec = spec[cont_mask]
#     line_spec = spec[line_mask]

#     coeffs_fit = np.polyfit(wave_cont, cont_spec, deg=deg)
#     cont_fit = np.polyval(coeffs_fit, line_wave)
#     ew = np.trapz(1 - (line_spec / cont_fit), x=line_wave)

#     if spec_unc is not None:
#         line_unc = spec_unc[line_mask]
#         ew_unc = np.sqrt(np.sum(((line_unc / cont_fit) ** 2)))
#         return ew, ew_unc

#     return ew

# #equivelent width function 



def breakstrength(wave_rest, spec, blue_wave_edges, red_wave_edges, spec_unc=None): #defining what variables are needed
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
     
    blue_mask = (wave_rest>=blue_wave_edges[0])&(wave_rest<=blue_wave_edges[1])
    red_mask = (wave_rest>=red_wave_edges[0])&(wave_rest<=red_wave_edges[1])
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
        
        blue_err = np.sqrt(np.sum(blue_unc**2))/len(blue_spec)
        red_err = np.sqrt(np.sum(red_unc**2))/len(red_spec)
        
        break_unc = break_strength * np.sqrt((red_err / red_med)**2 + (blue_err / blue_med)**2) #Uncertainty of flux within each window 
        
        return break_strength, break_unc
    
    else:  
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
    
    assert len(filters) == 2, f"Only two filters are allowed, received {len(filters)}."
    
    mags = [] 
    mags_err = []
    
    for filt in filters:
        mag = filt.ab_mag(wave_rest, spec)
        mags.append(mag)

        if spec_unc is not None:
            filt_trans = np.interp(wave_rest, filt.wavelength, filt.transmission)
            dwave = np.gradient(wave_rest)
            flux = np.nansum(spec * filt_trans * dwave) / np.nansum(filt_trans * dwave)
            flux_err = np.sqrt(np.nansum((filt_trans * spec_unc * dwave)**2)) / np.nansum(filt_trans * dwave)
            mag_err = 1.086 * flux_err / flux
            mags_err.append(mag_err)
    
    color = mags[0] - mags[1]
    
    if spec_unc is not None:
        color_err = np.sqrt(mags_err[0]**2 + mags_err[1]**2)
        return color, color_err
    else:
        return color


  



