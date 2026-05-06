# SUMMARY
Within the local universe, star formation is generally regulated by young, massive stars and supernovae in which cold gas is accreted. Over time, this balance leads to relatively smooth star formation fluctuations when averaged over long timescales. However, recent studies of these distant galaxies’ SEDs provide strong evidence of star formation occurring very rapidly and intensely over very short timescales, described as "burstiness." Galaxy brightness has been closely linked to this burstiness, as these events contribute to a galaxy’s observed luminosity (Wang et al. 2025). Because of the intensity of burstiness of the birth of young stellar populations, analyses of individual galaxies have become increasingly complicated by the effect of outshining. This effect has led to uncertainty in key galaxy properties, such as understanding galaxy stellar masses and ages. 

Observations of high-redshift galaxies allow us to piece together the history of our universe and how those galaxies came to be. The spectral energy distribution (SED) of galaxies describe the distribution of flux density across the electromagnetic spectrum, capture light and information about a galaxy’s current and past star formation activity, and contain variables such as stellar continuum, gas emission, dust attenuation and re-emission, and active galactic nucleus information.

Within this project, data from the SDSS and ManGA surveys are analyzed and calculated to find the equivlent line widths, continuum breaks, and restframe colors of spectra. Measuremnts are meant to serve as a universal tool for all those accsessing this code, and to analyze data, as we are trying to avoid significant differnces in potential measurements by anyone else .. **(re word this to we dont want significaint measurements so this code is like a universal tool for anyone to use so we dont have any biases or discrepencies)

-------------------------------------------------------------------------------------
# Codes imported and functions used/Important notes before starting
This project includes some of the following libraries: 

Numpy
Scipy
Matplotlib
Astropy
Time
Astroquery

Ensure that all appropriate libraries are downloaded before **replicating?**

-------------------------------------------------------------------------------------

# Equivelent Width Lines
Definitions **fix this**
Flux of Spectrum  - 
Restframe Wavelegnth of Spectrum - 
Uncertainty of flux -
Restframe wavelength - 
etc - 


Equivelenet Width Lines measure how strong a spectral line is in comparsion to the rest of its surrounding continuum. The line flux will measure how much brighter or fainter a feature is to the rest of the continuum, such as emission and absoprtion lines. Large lines above the continuum, alongside a measurement of a negative equivlent width represent emission lines. Large lines below the continuum, alongside a measurement of a positive equivelent width represent absoprtion lines. 

Equivelent Width Lines are measured using the formula of an integral of 1 subtracted by the emission line flux over the continuum flux:

∫ (1 - F_line / F_cont) dλ

Blue and red wavelegnths are on either side of the line continuum that captures the line continuum itself and its wings and surrounding continuum. The continuum is then esitmated by fitting a linear line to the continuum, and masks are created to have a spefific array of the continuum flux within the windows defined. 

