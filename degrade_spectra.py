#This script degrades spectral resolution of synthetic spectra xs

import numpy as np
import asciidata
from scipy import ndimage, interpolate
import os

def gauss(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) / np.sqrt(2 * np.pi * np.power(sig, 2))     


def degrade_spectra_Jband(Teff, logg, met):
    #Each model has its own directory
    root = '/Users/duisiya/astro/ABPic'
    dirname = str(Teff) + '_' + 'logg' + str(logg) + 'met' + str(met)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname, 'J_band_mod.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()

    #read observed spectrum
    spec_obs = asciidata.open(os.path.join(dirname, 'J_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    
    
    #BINNING SPECTRUM
    #
    #the synthetic spectrum needs to be binned to have fewer point in the wavelength dimension (~1700 instead of >100000) 
    dwave = (max(wave_obs) - min(wave_obs)) / len (wave_obs) #"fwhm" of the observed spectrum
    dwave_sigma = dwave / (2 * np.sqrt(2 * np.log(2)))
    
    #before interpolating the spectrum to a new wavelength grid, smooth it with a Gaussian
    flux_conv = ndimage.filters.gaussian_filter(flux_synth, dwave_sigma)
    
    #interpolating the smoothed synthetic spectrum on the observed wavelength grid
    func_interp = interpolate.splrep(wave_synth, flux_conv) #searching for interpolation function
    flux_synth_smooth = interpolate.splev(wave_obs, func_interp) #interpolating on a new grid
    
    
    #DEGRADING SPECTRAL RESOLUTION
    #
    R = 2000. #resolution of observed spectra in J band
    FWHM = wave_synth / R #FWHM of a line profile at a given wavelength
    G_sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
           
    #Convolution of the synthetic spectrum with a Gaussian function that has a variable sigma
    gauss_matr = np.vstack([np.hstack([gauss(x, mu, sig) for x in wave_obs]) for mu, sig in zip(wave_obs, G_sigma)])
    flux_synth_conv = np.dot(flux_synth_smooth * dwave, gauss_matr)  

    #writing the degraded spectrum to file
    flux_smooth_file = open(os.path.join(dirname, 'low_res_synth_spectrum_Jband.dat'), 'w')
    for n, m in zip(wave_obs, flux_synth_conv):   
        flux_smooth_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_smooth_file.close()



def degrade_spectra_Hband(Teff, logg, met):
    #Each model has its own directory
    root = '/Users/duisiya/astro/ABPic'
    dirname = str(Teff) + '_' + 'logg' + str(logg) + 'met' + str(met)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname, 'H_band_mod.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    #read observed spectrum
    spec_obs = asciidata.open(os.path.join(dirname, 'H_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    
    #BINNING SPECTRUM
    #
    #the synthetic spectrum needs to be binned to have fewer point in the wavelength dimension (~1700 instead of >100000) 
    dwave = (max(wave_obs) - min(wave_obs)) / len (wave_obs) #"fwhm" of the observed spectrum
    dwave_sigma = dwave / (2 * np.sqrt(2 * np.log(2)))
    
    #before interpolating the spectrum to a new wavelength grid, smooth it with a Gaussian
    flux_conv = ndimage.filters.gaussian_filter(flux_synth, dwave_sigma)
    
    #interpolating the smoothed synthetic spectrum on the observed wavelength grid
    func_interp = interpolate.splrep(wave_synth, flux_conv) #searching for interpolation function
    flux_synth_smooth = interpolate.splev(wave_obs, func_interp) #interpolating on a new grid
    
    #DEGRADING SPECTRAL RESOLUTION
    #
    R = 1500. #resolution of observed spectra in J band
    FWHM = wave_synth / R #FWHM of a line profile at a given wavelength
    G_sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    
    #Convolution of the synthetic spectrum with a Gaussian function that has a variable sigma
    gauss_matr = np.vstack([np.hstack([gauss(x, mu, sig) for x in wave_obs]) for mu, sig in zip(wave_obs, G_sigma)])
    flux_synth_conv = np.dot(flux_synth_smooth * dwave, gauss_matr)  
    
    #writing the degraded spectrum to file
    flux_smooth_file = open(os.path.join(dirname, 'low_res_synth_spectrum_Hband.dat'), 'w')
    for n, m in zip(wave_obs, flux_synth_conv):   
        flux_smooth_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_smooth_file.close()



def degrade_spectra_Kband(Teff, logg, met):
    #Each model has its own directory
    root = '/Users/duisiya/astro/ABPic'
    dirname = str(Teff) + '_' + 'logg' + str(logg) + 'met' + str(met)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname, 'K_band_mod.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    #read observed spectrum
    spec_obs = asciidata.open(os.path.join(dirname, 'K_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    
    #BINNING SPECTRUM
    #
    #the synthetic spectrum needs to be binned to have fewer point in the wavelength dimension (~1700 instead of >100000) 
    dwave = (max(wave_obs) - min(wave_obs)) / len (wave_obs) #"fwhm" of the observed spectrum
    dwave_sigma = dwave / (2 * np.sqrt(2 * np.log(2)))
    
    #before interpolating the spectrum to a new wavelength grid, smooth it with a Gaussian
    flux_conv = ndimage.filters.gaussian_filter(flux_synth, dwave_sigma)
    
    #interpolating the smoothed synthetic spectrum on the observed wavelength grid
    func_interp = interpolate.splrep(wave_synth, flux_conv) #searching for interpolation function
    flux_synth_smooth = interpolate.splev(wave_obs, func_interp) #interpolating on a new grid

    #DEGRADING SPECTRAL RESOLUTION
    #    
    R = 1500. #resolution of observed spectra in J band
    FWHM = wave_synth / R #FWHM of a line profile at a given wavelength
    G_sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    
    #Convolution of the synthetic spectrum with a Gaussian function that has a variable sigma
    gauss_matr = np.vstack([np.hstack([gauss(x, mu, sig) for x in wave_obs]) for mu, sig in zip(wave_obs, G_sigma)])
    flux_synth_conv = np.dot(flux_synth_smooth * dwave, gauss_matr)  
    
    #writing the degraded spectrum to file
    flux_smooth_file = open(os.path.join(dirname, 'low_res_synth_spectrum_Kband.dat'), 'w')
    for n, m in zip(wave_obs, flux_synth_conv):   
        flux_smooth_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_smooth_file.close()



def degrade_spectra_C_O_Jband(Teff, logg, met, C_O):
    #Each model has its own directory
    root = '/Users/duisiya/astro/ABPic'
    dirname = str(Teff) + '_' + 'logg' + str(logg) + 'met' + str(met) + 'C_O' + str(C_O)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname, 'J_band_mod.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    #read observed spectrum
    spec_obs = asciidata.open(os.path.join(dirname, 'J_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    
    #BINNING SPECTRUM
    #
    #the synthetic spectrum needs to be binned to have fewer point in the wavelength dimension (~1700 instead of >100000) 
    dwave = (max(wave_obs) - min(wave_obs)) / len (wave_obs) #"fwhm" of the observed spectrum
    dwave_sigma = dwave / (2 * np.sqrt(2 * np.log(2)))
    
    #before interpolating the spectrum to a new wavelength grid, smooth it with a Gaussian
    flux_conv = ndimage.filters.gaussian_filter(flux_synth, dwave_sigma)
    
    #interpolating the smoothed synthetic spectrum on the observed wavelength grid
    func_interp = interpolate.splrep(wave_synth, flux_conv) #searching for interpolation function
    flux_synth_smooth = interpolate.splev(wave_obs, func_interp) #interpolating on a new grid

    #DEGRADING SPECTRAL RESOLUTION
    #     
    R = 2000. #resolution of observed spectra in J band
    FWHM = wave_synth / R #FWHM of a line profile at a given wavelength
    G_sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
      
    #Convolution of the synthetic spectrum with a Gaussian function that has a variable sigma
    gauss_matr = np.vstack([np.hstack([gauss(x, mu, sig) for x in wave_obs]) for mu, sig in zip(wave_obs, G_sigma)])
    flux_synth_conv = np.dot(flux_synth_smooth * dwave, gauss_matr)  
    
    #writing the degraded spectrum to file
    flux_smooth_file = open(os.path.join(dirname, 'low_res_synth_spectrum_Jband.dat'), 'w')
    for n, m in zip(wave_obs, flux_synth_conv):   
        flux_smooth_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_smooth_file.close()



def degrade_spectra_C_O_Hband(Teff, logg, met, C_O):
    #Each model has its own directory
    root = '/Users/duisiya/astro/ABPic'
    dirname = str(Teff) + '_' + 'logg' + str(logg) + 'met' + str(met) + 'C_O' + str(C_O)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname, 'H_band_mod.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    #read observed spectrum
    spec_obs = asciidata.open(os.path.join(dirname, 'H_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    
    #BINNING SPECTRUM
    #
    #the synthetic spectrum needs to be binned to have fewer point in the wavelength dimension (~1700 instead of >100000) 
    dwave = (max(wave_obs) - min(wave_obs)) / len (wave_obs) #"fwhm" of the observed spectrum
    dwave_sigma = dwave / (2 * np.sqrt(2 * np.log(2)))
    
    #before interpolating the spectrum to a new wavelength grid, smooth it with a Gaussian
    flux_conv = ndimage.filters.gaussian_filter(flux_synth, dwave_sigma)
    
    #interpolating the smoothed synthetic spectrum on the observed wavelength grid
    func_interp = interpolate.splrep(wave_synth, flux_conv) #searching for interpolation function
    flux_synth_smooth = interpolate.splev(wave_obs, func_interp) #interpolating on a new grid

    #DEGRADING SPECTRAL RESOLUTION
    #     
    R = 1500. #resolution of observed spectra in J band
    FWHM = wave_synth / R #FWHM of a line profile at a given wavelength
    G_sigma = FWHM / (2 * np.sqrt(2 * np.log(2))) 

    #Convolution of the synthetic spectrum with a Gaussian function that has a variable sigma
    gauss_matr = np.vstack([np.hstack([gauss(x, mu, sig) for x in wave_obs]) for mu, sig in zip(wave_obs, G_sigma)])
    flux_synth_conv = np.dot(flux_synth_smooth * dwave, gauss_matr)  
    
    #writing the degraded spectrum to file
    flux_smooth_file = open(os.path.join(dirname, 'low_res_synth_spectrum_Hband.dat'), 'w')
    for n, m in zip(wave_obs, flux_synth_conv):   
        flux_smooth_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_smooth_file.close()



def degrade_spectra_C_O_Kband(Teff, logg, met, C_O):
    #Each model has its own directory
    root = '/Users/duisiya/astro/ABPic'
    dirname = str(Teff) + '_' + 'logg' + str(logg) + 'met' + str(met) + 'C_O' + str(C_O)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname, 'K_band_mod.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    #read observed spectrum
    spec_obs = asciidata.open(os.path.join(dirname, 'K_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    
    #BINNING SPECTRUM
    #
    #the synthetic spectrum needs to be binned to have fewer point in the wavelength dimension (~1700 instead of >100000) 
    dwave = (max(wave_obs) - min(wave_obs)) / len (wave_obs) #"fwhm" of the observed spectrum
    dwave_sigma = dwave / (2 * np.sqrt(2 * np.log(2)))
    
    #before interpolating the spectrum to a new wavelength grid, smooth it with a Gaussian
    flux_conv = ndimage.filters.gaussian_filter(flux_synth, dwave_sigma)
    
    #interpolating the smoothed synthetic spectrum on the observed wavelength grid
    func_interp = interpolate.splrep(wave_synth, flux_conv) #searching for interpolation function
    flux_synth_smooth = interpolate.splev(wave_obs, func_interp) #interpolating on a new grid

    #DEGRADING SPECTRAL RESOLUTION
    #     
    R = 1500. #resolution of observed spectra in J band
    FWHM = wave_synth / R #FWHM of a line profile at a given wavelength
    G_sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    
    #Convolution of the synthetic spectrum with a Gaussian function that has a variable sigma
    gauss_matr = np.vstack([np.hstack([gauss(x, mu, sig) for x in wave_obs]) for mu, sig in zip(wave_obs, G_sigma)])
    flux_synth_conv = np.dot(flux_synth_smooth * dwave, gauss_matr)  
    
    #writing the degraded spectrum to file
    flux_smooth_file = open(os.path.join(dirname, 'low_res_synth_spectrum_Kband.dat'), 'w')
    for n, m in zip(wave_obs, flux_synth_conv):   
        flux_smooth_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_smooth_file.close()


root = '/Users/duisiya/astro/ABPic/models/'
l = os.listdir(root)

for f in l:
    if not f.startswith('.') and os.path.isfile(os.path.join(root, f)): #ignore hidden files
        Teff = str(f[3:7])
        logg = str(f[8:12])
        met = str(f[13:16])
        degrade_spectra_Jband(Teff, logg, met)
        degrade_spectra_Hband(Teff, logg, met)
        degrade_spectra_Kband(Teff, logg, met)

root = '/Users/duisiya/astro/ABPic/models_C_O/'
l = os.listdir(root)



for f in l:
    if not f.startswith('.') and os.path.isfile(os.path.join(root, f)): #ignore hidden files
        Teff=str(f[3:7])
        logg=str(f[8:12])
        met=str(f[13:16])
        C_O=str(f[22:26])
        degrade_spectra_C_O_Jband(Teff, logg, met, C_O)
        degrade_spectra_C_O_Hband(Teff, logg, met, C_O)
        degrade_spectra_C_O_Kband(Teff, logg, met, C_O)