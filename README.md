# SUMMARY
Within the local universe, star formation is generally regulated by young, massive stars and supernovae in which cold gas is accreted. Over time, this balance leads to relatively smooth star formation fluctuations when averaged over long timescales. However, recent studies of these distant galaxies’ SEDs provide strong evidence of star formation occurring very rapidly and intensely over very short timescales, described as "burstiness." Galaxy brightness has been closely linked to this burstiness, as these events contribute to a galaxy’s observed luminosity (Wang et al. 2025). Because of the intensity of burstiness of the birth of young stellar populations, analyses of individual galaxies have become increasingly complicated by the effect of outshining. This effect has led to uncertainty in key galaxy properties, such as understanding galaxy stellar masses and ages. 

Observations of high-redshift galaxies allow us to piece together the history of our universe and how those galaxies came to be. The spectral energy distribution (SED) of galaxies describes the distribution of flux density across the electromagnetic spectrum, captures light and information about a galaxy’s current and past star formation activity, and contain variables such as stellar continuum, gas emission, dust attenuation and re-emission, and active galactic nucleus information.

Within this project, data from the Sloan Digital Sky Survey (SDSS) and Mapping Nearby Galaxies at APO (MaNGA) surveys are analyzed and calculated to find the equivlent line widths, continuum breaks, and restframe colors of spectra. The goal of these measurements is to create a consistent and universal analysis tool that can be used by anyone accessing the codebase. By ensuring consistency within our program, we aim to reduce discrepancies, biases, and inconsistencies that may arise when spectra are analyzed differently and or independently. 

-------------------------------------------------------------------------------------
# Codes imported and functions used/Important notes before starting
This project includes some of the following libraries: 

Numpy
Scipy
Matplotlib
Astropy
Time
Astroquery

Ensure that all required libraries are downloaded in your environment, and it is highly recommended to briefly review the imported packages beforehand to verify if any additional installations are needed.

-------------------------------------------------------------------------------------
# Equivelent Width Lines
Equivelenet Width Lines measure how strong a spectral line is in comparsion to the rest of its surrounding continuum. The line flux will measure how much brighter or fainter a feature is to the rest of the continuum, such as emission and absoprtion lines. Large lines above the continuum, alongside a measurement of a negative equivlent width represent emission lines. Large lines below the continuum, alongside a measurement of a positive equivelent width represent absoprtion lines. 

Equivelent Width Lines are measured using the formula of an integral of 1 subtracted by the emission line flux over the continuum flux:

∫ (1 - F_line / F_cont) dλ

Blue and red wavelegnths are on either side of the line continuum that captures the line continuum itself and its wings and surrounding continuum. The continuum is then estimated by fitting a linear line to the continuum, and masks are created to have a specific array of the continuum flux within the windows defined. 

-------------------------------------------------------------------------------------
# Continuum Breaks

Continuum breaks are tracers of older stellar populations, calculated as the ratio between the median flux in two windows on either side of the break (often termed the "red" and "blue" sides.) The Balmer break traces hydrogen absorption from cooler, intermediate-aged stars, while the DN4000 break traces metal abopsrtion features in  significantly older stellar populations.

-------------------------------------------------------------------------------------
# RestFrame Wavlegenths
Restframe Wavelengths are the intrisic colors a galaxy emits within its restframes, assuming the galaxy is not affected by redshift. 

Colors are calculated by calculating the magnitude of a spectrum within a specific bandpass, and using a library named astro-sedpy, we assume a photometric sensitivty curve, which is the amount of light a filter is allowed to let through at a specific wavelength.

Colors are calculated by subtracting the rest from magnitudes of a spectrum in a filter, e.g:

Color = mag_A - mag_B
