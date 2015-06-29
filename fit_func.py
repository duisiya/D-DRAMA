import numpy as np
import asciidata
import os


    
def fit_spectra_Jband(Teff, logg, met):
    dirname = str(Teff) + '_' + 'logg' + str(logg) + 'met' + str(met)
    
    spec_obs = asciidata.open(os.path.join(dirname, 'J_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    sigma_obs = spec_obs[2].tonumpy()

    spec_synth = asciidata.open(os.path.join(dirname,'low_res_synth_norm_Jband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()

    ivar = 1 / np.power(sigma_obs, 2)
    L = (max(wave_synth) - min(wave_synth)) * 2

    K = 10
    hs = np.vstack([np.vstack([np.cos(2 * np.pi * k * wave_synth/L) for k in range(K)]), np.vstack([np.sin(2 * np.pi * k * wave_synth/L) for k in range(1,K)])]) #create base matrix

    A = hs * flux_synth
    ATA = np.dot(A, (ivar * A).T)
    ATy = np.dot(A, ivar * flux_obs)
    x = np.linalg.solve(ATA, ATy)
    fit = np.dot(A.T, x)
    chi_sq = np.sum(ivar * (flux_obs - fit) ** 2)
    func = np.dot(x, hs)
    
    best_fit_file = open(os.path.join(dirname, 'best_fit_Jband.dat'), 'w')
    for n, m in zip(wave_synth, fit):
        best_fit_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    best_fit_file.close()

    fitting_func_file = open(os.path.join(dirname, 'fitting_func_Jband.dat'), 'w')
    for n, m in zip(wave_synth, func):
        fitting_func_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    fitting_func_file.close()
    return chi_sq

def fit_spectra_Hband(Teff, logg, met):
    dirname = str(Teff) + '_' + 'logg' + str(logg) + 'met' + str(met)
        
    spec_obs = asciidata.open(os.path.join(dirname, 'H_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    sigma_obs = spec_obs[2].tonumpy()
    
    spec_synth = asciidata.open(os.path.join(dirname, 'low_res_synth_norm_Hband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    ivar = 1 / np.power(sigma_obs, 2)
    L = (max(wave_synth) - min(wave_synth)) * 2
    
    K = 10
    hs = np.vstack([np.vstack([np.cos(2 * np.pi * k * wave_synth/L) for k in range(K)]), np.vstack([np.sin(2 * np.pi * k * wave_synth/L) for k in range(1,K)])]) #create base matrix
    
    A = hs * flux_synth
    ATA = np.dot(A, (ivar * A).T)
    ATy = np.dot(A, ivar * flux_obs)
    x = np.linalg.solve(ATA, ATy)
    fit = np.dot(A.T, x)
    chi_sq = np.sum(ivar * (flux_obs - fit) ** 2)
    func = np.dot(x, hs)

    best_fit_file = open(os.path.join(dirname, 'best_fit_Hband.dat'), 'w')
    for n, m in zip(wave_synth, fit):
        best_fit_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    best_fit_file.close()

    fitting_func_file = open(os.path.join(dirname, 'fitting_func_Hband.dat'), 'w')
    for n, m in zip(wave_synth, func):
        fitting_func_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    fitting_func_file.close()
    return chi_sq


def fit_spectra_Kband(Teff, logg, met):
    dirname = str(Teff) + '_' + 'logg' + str(logg) + 'met' + str(met)
        
    spec_obs = asciidata.open(os.path.join(dirname, 'K_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    sigma_obs = spec_obs[2].tonumpy()
    
    spec_synth = asciidata.open(os.path.join(dirname, 'low_res_synth_norm_Kband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    ivar = 1 / np.power(sigma_obs, 2)
    L = (max(wave_synth) - min(wave_synth)) * 2
    
    K = 10
    hs = np.vstack([np.vstack([np.cos(2 * np.pi * k * wave_synth/L) for k in range(K)]), np.vstack([np.sin(2 * np.pi * k * wave_synth/L) for k in range(1,K)])]) #create base matrix
    
    A = hs * flux_synth
    ATA = np.dot(A, (ivar * A).T)
    ATy = np.dot(A, ivar * flux_obs)
    x = np.linalg.solve(ATA, ATy)
    fit = np.dot(A.T, x)
    chi_sq = np.sum(ivar * (flux_obs - fit) ** 2)
    func = np.dot(x, hs)
    
    best_fit_file = open(os.path.join(dirname, 'best_fit_Kband.dat'), 'w')
    for n, m in zip(wave_synth, fit):
        best_fit_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    best_fit_file.close()

    fitting_func_file = open(os.path.join(dirname, 'fitting_func_Kband.dat'), 'w')
    for n, m in zip(wave_synth, func):
        fitting_func_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    fitting_func_file.close()
    return chi_sq

def fit_spectra_C_O_Jband(Teff_C_O, logg_C_O, met_C_O, C_O):
    dirname = str(Teff_C_O) + '_' + 'logg' + str(logg_C_O) + 'met' + str(met_C_O) + 'C_O' + str(C_O)
    
    spec_obs = asciidata.open(os.path.join(dirname, 'J_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    sigma_obs = spec_obs[2].tonumpy()
    
    spec_synth = asciidata.open(os.path.join(dirname, 'low_res_synth_norm_Jband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    ivar = 1 / np.power(sigma_obs, 2)
    L = (max(wave_synth) - min(wave_synth)) * 2
    
    K = 10
    hs = np.vstack([np.vstack([np.cos(2 * np.pi * k * wave_synth/L) for k in range(K)]), np.vstack([np.sin(2 * np.pi * k * wave_synth/L) for k in range(1,K)])]) #create base matrix
    
    A = hs * flux_synth
    ATA = np.dot(A, (ivar * A).T)
    ATy = np.dot(A, ivar * flux_obs)
    x = np.linalg.solve(ATA, ATy)
    fit = np.dot(A.T, x)
    chi_sq = np.sum(ivar * (flux_obs - fit) ** 2)
    func = np.dot(x, hs)
    
    best_fit_file = open(os.path.join(dirname, 'best_fit_Jband.dat'), 'w')
    for n, m in zip(wave_synth, fit):
        best_fit_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    best_fit_file.close()

    fitting_func_file = open(os.path.join(dirname, 'fitting_func_Jband.dat'), 'w')
    for n, m in zip(wave_synth, func):
        fitting_func_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    fitting_func_file.close()
    return chi_sq


def fit_spectra_C_O_Hband(Teff_C_O, logg_C_O, met_C_O, C_O):
    dirname = str(Teff_C_O) + '_' + 'logg' + str(logg_C_O) + 'met' + str(met_C_O) + 'C_O' + str(C_O)
    
    spec_obs = asciidata.open(os.path.join(dirname, 'H_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    sigma_obs = spec_obs[2].tonumpy()
    
    spec_synth = asciidata.open(os.path.join(dirname, 'low_res_synth_norm_Hband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    ivar = 1 / np.power(sigma_obs, 2)
    L = (max(wave_synth) - min(wave_synth)) * 2
    
    K = 10
    hs = np.vstack([np.vstack([np.cos(2 * np.pi * k * wave_synth/L) for k in range(K)]), np.vstack([np.sin(2 * np.pi * k * wave_synth/L) for k in range(1,K)])]) #create base matrix
    
    A = hs * flux_synth
    ATA = np.dot(A, (ivar * A).T)
    ATy = np.dot(A, ivar * flux_obs)
    x = np.linalg.solve(ATA, ATy)
    fit = np.dot(A.T, x)
    chi_sq = np.sum(ivar * (flux_obs - fit) ** 2)
    func = np.dot(x, hs)
    
    best_fit_file = open(os.path.join(dirname, 'best_fit_Hband.dat'), 'w')
    for n, m in zip(wave_synth, fit):
        best_fit_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    best_fit_file.close()

    fitting_func_file = open(os.path.join(dirname, 'fitting_func_Hband.dat'), 'w')
    for n, m in zip(wave_synth, func):
        fitting_func_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    fitting_func_file.close()
    return chi_sq


def fit_spectra_C_O_Kband(Teff_C_O, logg_C_O, met_C_O, C_O):
    dirname = str(Teff_C_O) + '_' + 'logg' + str(logg_C_O) + 'met' + str(met_C_O) + 'C_O' + str(C_O)
    
    spec_obs = asciidata.open(os.path.join(dirname, 'K_band_obs.dat'))
    wave_obs = spec_obs[0].tonumpy()
    flux_obs = spec_obs[1].tonumpy()
    sigma_obs = spec_obs[2].tonumpy()
    
    spec_synth = asciidata.open(os.path.join(dirname, 'low_res_synth_norm_Kband.dat'))
    wave_synth = spec_synth[0].tonumpy()
    flux_synth = spec_synth[1].tonumpy()
    
    ivar = 1 / np.power(sigma_obs, 2)
    L = (max(wave_synth) - min(wave_synth)) * 2
    
    K = 10
    hs = np.vstack([np.vstack([np.cos(2 * np.pi * k * wave_synth/L) for k in range(K)]), np.vstack([np.sin(2 * np.pi * k * wave_synth/L) for k in range(1,K)])]) #create base matrix
    
    A = hs * flux_synth
    ATA = np.dot(A, (ivar * A).T)
    ATy = np.dot(A, ivar * flux_obs)
    x = np.linalg.solve(ATA, ATy)
    fit = np.dot(A.T, x)
    chi_sq = np.sum(ivar * (flux_obs - fit) ** 2)
    func = np.dot(x, hs)
        
    best_fit_file = open(os.path.join(dirname, 'best_fit_Kband.dat'), 'w')
    for n, m in zip(wave_synth, fit):
        best_fit_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    best_fit_file.close()

    fitting_func_file = open(os.path.join(dirname, 'fitting_func_Kband.dat'), 'w')
    for n, m in zip(wave_synth, func):
        fitting_func_file.write(str(n) + ' ' + str(m) + ' ' + '\n')
    fitting_func_file.close()
    return chi_sq



root = '/Users/duisiya/astro/ABPic/models/'
l = os.listdir(root)
root2 = '/Users/duisiya/astro/ABPic/'



#OK, I don't like this business with arrays myself but solar and super-solar C/O models have they names written in slighltly different way that why I read them separetly. Arrays is the way to merge parameters for solor and super-solar C/O models after model fitting is done. I need to think about something more sophisticated here.

Teff_arr = []
logg_arr = []
met_arr = []
chi_sq_Jband_arr = []
chi_sq_Hband_arr = []
chi_sq_Kband_arr = []

for f in l:
    if not f.startswith('.') and os.path.isfile(os.path.join(root, f)): #ignore hidden files
        Teff = str(f[3:7])
        logg = str(f[8:12])
        met = str(f[13:16])
        Teff_arr.append(Teff)
        logg_arr.append(logg)
        met_arr.append(met)
        fit_spectra_Jband(Teff, logg, met)
        fit_spectra_Hband(Teff, logg, met)
        fit_spectra_Kband(Teff, logg, met)
        chi_sq_Jband = fit_spectra_Jband(Teff, logg, met)
        chi_sq_Hband = fit_spectra_Hband(Teff, logg, met)
        chi_sq_Kband = fit_spectra_Kband(Teff, logg, met)

        chi_sq_Jband_arr.append(chi_sq_Jband)
        chi_sq_Hband_arr.append(chi_sq_Hband)
        chi_sq_Kband_arr.append(chi_sq_Kband)


root = '/Users/duisiya/astro/ABPic/models_C_O/'
l = os.listdir(root)

Teff_C_O_arr = []
logg_C_O_arr = []
met_C_O_arr = []
C_O_arr = []
chi_sq_C_O_Jband_arr = []
chi_sq_C_O_Hband_arr = []
chi_sq_C_O_Kband_arr = []

for f in l:
    if not f.startswith('.') and os.path.isfile(os.path.join(root, f)): #ignore hidden files
        Teff_C_O = str(f[3:7])
        logg_C_O = str(f[8:12])
        met_C_O = str(f[13:16])
        C_O = str(f[22:26])
        Teff_C_O_arr.append(Teff_C_O)
        logg_C_O_arr.append(logg_C_O)
        met_C_O_arr.append(met_C_O)
        C_O_arr.append(C_O)
        fit_spectra_C_O_Jband(Teff_C_O, logg_C_O, met_C_O, C_O)
        fit_spectra_C_O_Hband(Teff_C_O, logg_C_O, met_C_O, C_O)
        fit_spectra_C_O_Kband(Teff_C_O, logg_C_O, met_C_O, C_O)
        chi_sq_C_O_Jband = fit_spectra_C_O_Jband(Teff_C_O, logg_C_O, met_C_O, C_O)
        chi_sq_C_O_Hband = fit_spectra_C_O_Hband(Teff_C_O, logg_C_O, met_C_O, C_O)
        chi_sq_C_O_Kband = fit_spectra_C_O_Kband(Teff_C_O, logg_C_O, met_C_O, C_O)
        chi_sq_C_O_Jband_arr.append(chi_sq_C_O_Jband)
        chi_sq_C_O_Hband_arr.append(chi_sq_C_O_Hband)
        chi_sq_C_O_Kband_arr.append(chi_sq_C_O_Kband)

Teff_comb = np.vstack([np.vstack(Teff_arr), np.vstack(Teff_C_O_arr)])
logg_comb = np.vstack([np.vstack(logg_arr), np.vstack(logg_C_O_arr)])
met_comb = np.vstack([np.vstack(met_arr), np.vstack(met_C_O_arr)])
C_O_comb = np.vstack([np.vstack([0.0] * len(Teff_arr)), np.vstack(C_O_arr)])
chi_sq_comb_Jband = np.vstack([np.vstack(chi_sq_Jband_arr), np.vstack(chi_sq_C_O_Jband_arr)])
chi_sq_comb_Hband = np.vstack([np.vstack(chi_sq_Hband_arr), np.vstack(chi_sq_C_O_Hband_arr)])
chi_sq_comb_Kband = np.vstack([np.vstack(chi_sq_Kband_arr), np.vstack(chi_sq_C_O_Kband_arr)])

chi_sq_comb_all = chi_sq_comb_Jband + chi_sq_comb_Hband + chi_sq_comb_Kband

results = np.hstack([Teff_comb, logg_comb, met_comb, C_O_comb, chi_sq_comb_Jband, chi_sq_comb_Hband, chi_sq_comb_Kband, chi_sq_comb_all])

chi_sq_file = open(os.path.join(root2, 'chi_sq_values.dat'), 'w')
for line in results:
    chi_sq_file.write(' '.join(line) + '\n')
chi_sq_file.close()

#Find indices for array elements with lowest chi-squares
index_Jband = np.argmin(chi_sq_comb_Jband)
index_Hband = np.argmin(chi_sq_comb_Hband)
index_Kband = np.argmin(chi_sq_comb_Kband)
index_all = np.argmin(chi_sq_comb_all)

chi_sq_best_file = open(os.path.join(root2, 'best_fits.dat'), 'w')
chi_sq_best_file.write(' '.join(results[index_Jband]) + ' ' + 'best Jband fit' + '\n')
chi_sq_best_file.write(' '.join(results[index_Hband]) + ' ' + 'best Hband fit' + '\n')
chi_sq_best_file.write(' '.join(results[index_Kband]) + ' ' + 'best Kband fit' + '\n')
chi_sq_best_file.write(' '.join(results[index_all]) + ' ' + 'best fit all bands' + '\n')
chi_sq_best_file.close()







