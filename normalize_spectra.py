import numpy as np
import asciidata
import os
import asciidata

def normalize_spectra_Jband(Teff, logg, met):
    #Each model has its own directory
    root='/Users/duisiya/astro/ABPic'
    dirname=str(Teff)+'_'+'logg'+str(logg)+'met'+str(met)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname,'low_res_synth_spectrum_Jband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()


    #Continuum normalization (assuming continuum at 1.226 microns as in Bonnefoy et al. 2010)
    index = min(range(len(wave_synth)), key = lambda i: abs(wave_synth[i] - 1.226))

    flux_synth_norm = flux_synth / flux_synth[index]
    flux_norm_file = open(os.path.join(dirname, 'low_res_synth_norm_Jband.dat'), 'w')
    for n, m in zip(wave_synth, flux_synth_norm):   
        flux_norm_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_norm_file.close()

def normalize_spectra_Hband(Teff, logg, met):
    #Each model has its own directory
    root='/Users/duisiya/astro/ABPic'
    dirname=str(Teff)+'_'+'logg'+str(logg)+'met'+str(met)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname,'low_res_synth_spectrum_Hband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    
    #Continuum normalization (assuming continuum at 1.65 microns as in Bonnefoy et al. 2010)
    index = min(range(len(wave_synth)), key = lambda i: abs(wave_synth[i] - 1.65))
    
    flux_synth_norm = flux_synth / flux_synth[index]
    flux_norm_file = open(os.path.join(dirname, 'low_res_synth_norm_Hband.dat'), 'w')
    for n, m in zip(wave_synth, flux_synth_norm):   
        flux_norm_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_norm_file.close()

def normalize_spectra_Kband(Teff, logg, met):
    #Each model has its own directory
    root='/Users/duisiya/astro/ABPic'
    dirname=str(Teff)+'_'+'logg'+str(logg)+'met'+str(met)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname,'low_res_synth_spectrum_Kband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    
    #Continuum normalization (assuming continuum at 2.2 microns as in Bonnefoy et al. 2010)
    index = min(range(len(wave_synth)), key = lambda i: abs(wave_synth[i] - 2.2))
    
    flux_synth_norm = flux_synth / flux_synth[index]
    flux_norm_file = open(os.path.join(dirname, 'low_res_synth_norm_Kband.dat'), 'w')
    for n, m in zip(wave_synth, flux_synth_norm):   
        flux_norm_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_norm_file.close()

def normalize_spectra_C_O_Jband(Teff, logg, met, C_O):
    #Each model has its own directory
    root='/Users/duisiya/astro/ABPic'
    dirname=str(Teff)+'_'+'logg'+str(logg)+'met'+str(met)+'C_O'+str(C_O)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname,'low_res_synth_spectrum_Jband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    
    #Continuum normalization (assuming continuum at 1.226 microns as in Bonnefoy et al. 2010)
    index = min(range(len(wave_synth)), key = lambda i: abs(wave_synth[i] - 1.226))
    
    flux_synth_norm = flux_synth / flux_synth[index]
    flux_norm_file = open(os.path.join(dirname, 'low_res_synth_norm_Jband.dat'), 'w')
    for n, m in zip(wave_synth, flux_synth_norm):   
        flux_norm_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_norm_file.close()

def normalize_spectra_C_O_Hband(Teff, logg, met, C_O):
    #Each model has its own directory
    root='/Users/duisiya/astro/ABPic'
    dirname=str(Teff)+'_'+'logg'+str(logg)+'met'+str(met)+'C_O'+str(C_O)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname,'low_res_synth_spectrum_Hband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    
    #Continuum normalization (assuming continuum at 1.65 microns as in Bonnefoy et al. 2010)
    index = min(range(len(wave_synth)), key = lambda i: abs(wave_synth[i] - 1.65))
    
    flux_synth_norm = flux_synth / flux_synth[index]
    flux_norm_file = open(os.path.join(dirname, 'low_res_synth_norm_Hband.dat'), 'w')
    for n, m in zip(wave_synth, flux_synth_norm):   
        flux_norm_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_norm_file.close()

def normalize_spectra_C_O_Kband(Teff, logg, met, C_O):
    #Each model has its own directory
    root='/Users/duisiya/astro/ABPic'
    dirname=str(Teff)+'_'+'logg'+str(logg)+'met'+str(met)+'C_O'+str(C_O)
    
    #read synthetic spectrum
    spec_synth = asciidata.open(os.path.join(dirname,'low_res_synth_spectrum_Kband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    
    #Continuum normalization (assuming continuum at 2.2 microns as in Bonnefoy et al. 2010)
    index = min(range(len(wave_synth)), key = lambda i: abs(wave_synth[i] - 2.2))
    
    flux_synth_norm = flux_synth / flux_synth[index]
    flux_norm_file = open(os.path.join(dirname, 'low_res_synth_norm_Kband.dat'), 'w')
    for n, m in zip(wave_synth, flux_synth_norm):   
        flux_norm_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    flux_norm_file.close()


root = '/Users/duisiya/astro/ABPic/models/'
l = os.listdir(root)

for f in l:
    if not f.startswith('.') and os.path.isfile(os.path.join(root, f)): #ignore hidden files
        Teff = str(f[3:7])
        logg = str(f[8:12])
        met = str(f[13:16])
        normalize_spectra_Jband(Teff, logg, met)
        normalize_spectra_Hband(Teff, logg, met)
        normalize_spectra_Kband(Teff, logg, met)

root = '/Users/duisiya/astro/ABPic/models_C_O/'
l = os.listdir(root)

for f in l:
    if not f.startswith('.') and os.path.isfile(os.path.join(root, f)): #ignore hidden files
        Teff=str(f[3:7])
        logg=str(f[8:12])
        met=str(f[13:16])
        C_O=str(f[22:26])
        normalize_spectra_C_O_Jband(Teff, logg, met, C_O)
        normalize_spectra_C_O_Hband(Teff, logg, met, C_O)
        normalize_spectra_C_O_Kband(Teff, logg, met, C_O)

#THIS OPTION IS TO BE DISCUSSED
#
#The continuum is presented as a combination of sin and cos functions.
#The matrix equation F = AX, where F is flux, A is coefficient matrix (values for sin and cos ar each wavelenght) and X is solutution matrix

#K = 2
#L = (max(wave_synth) - min(wave_synth)) *2
#A = np.vstack([np.vstack([np.cos(2*np.pi*k*wave_synth/L) for k in range(K)]),
#               np.vstack([np.sin(2*np.pi*k*wave_synth/L) for k in range(1,K)])]) #coefficient matrix A
#
#ATA = np.dot(A, A.T)
#ATy = np.dot(A, flux_synth)
#X = np.linalg.solve(ATA, ATy)#solution
#F = np.dot(X, A) #continuum



