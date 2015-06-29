# D-DRAMA

One day this collection of scripts will fit exoplanet and brown dwarf spectra in the best way possible.



Degrade_spectra.py decreases a spectral resolution of synthetic spectra.
Since I want to be very philosophical, I use a variable sigma for a gaussina filter.
For this I calculate a gauss filter for each wavelength. This is too expensive for the original 
synthetic spectra (that have more than 100,000 wavelength points).
That's why I do the follwoing:
1. "Smooth" spectra using gaussian filter with constant sigma and interpolation. This reduces the number of wavelength points.
2. Apply gaussin filter witn variable sigma (basically, matrix multiplication).



Normalize_spectra.py normalizes (obviously) spectra in a straightforward and (probably) 
wrong way - by dividing spectra to a flux at a given wavelength (following  Bonnefoy et al. 2010).


