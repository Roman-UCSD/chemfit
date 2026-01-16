############################################################
#                                                          #
#                Goettingen PHOENIX PRESET                 #
#                                                          #
#                                                          #
############################################################

import os, pickle
import glob
from astropy.io import fits
import scipy as scp
import re
from ftplib import FTP

settings = {
    ### Where to store downloaded models ###
    'griddir': original_settings['griddir'],

    ### List of models available for download ###
    'model_list': original_settings['model_list'],

    ### FTP server with models ###
    'ftp_server': 'phoenix.astro.physik.uni-goettingen.de',

    ### Directory on FTP server with models ###
    'ftp_dir': '/HiResFITS/',

    ### Goettingen grid wavelength file ###
    'wl_filename': original_settings['wl_filename'],

    ### Which parameters to fit? ###
    'fit_dof': ['zscale', 'teff', 'logg'],

    ### Virtual grid dimensions ###
    'virtual_dof': {'redshift': [-200, 200]},
}

def read_grid_dimensions(flush_cache = False):
    global __model_parameters

    f = open(settings['model_list'], 'r')
    models = f.read().strip().split('\n')
    f.close()

    __model_parameters = {}
    for model in models:
        params = {}
        file = model.split('/')[-1]
        if file.endswith('.txt') or file.startswith('WAVE_PHOENIX'):
            continue
        parsed = re.findall('lte([0-9]+)([+-][0-9.]+)([+-][0-9.]+).PHOENIX-ACES-AGSS-COND-2011-HiRes.fits', file)
        if len(parsed) != 1:
            parsed = re.findall('lte([0-9]+)([+-][0-9.]+)([+-][0-9.]+)\.Alpha=([+-][0-9.]+)\.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits', file)
            if len(parsed) != 1:
                print(file)
            params['alpha'] = float(parsed[0][3])
        else:
            params['alpha'] = 0.0
        params['teff'] = int(parsed[0][0])
        params['logg'] = -float(parsed[0][1])
        params['zscale'] = float(parsed[0][2])
        model_id = 't{teff:.3f}l{logg:.3f}z{zscale:.3f}a{alpha:.3f}'.format(**params).replace('-0.000', '0.000')
        __model_parameters[model_id] = {**params, 'filename': model}

    grid = {}
    for param in ['teff', 'logg', 'zscale', 'alpha']:
        grid[param] = np.unique([__model_parameters[model][param] for model in __model_parameters])

    return grid

def read_grid_model(params, grid):
    global __model_parameters, wl_grid

    try:
        __model_parameters
    except:
        read_grid_dimensions()

    try:
        wl_grid
    except:
        h = fits.open(settings['wl_filename'])
        wl_grid = h[0].data
        h.close()

    model_id = 't{teff:.3f}l{logg:.3f}z{zscale:.3f}a{alpha:.3f}'.format(**params).replace('-0.000', '0.000')
    if model_id not in __model_parameters:
        raise FileNotFoundError('Model with parameters {} not available'.format(params))
    filename = '{}/{}.fits'.format(settings['griddir'], model_id)
    if not os.path.isfile(filename):
        ftp = FTP(settings['ftp_server'])
        ftp.login(user = 'anonymous')
        ftp.cwd(settings['ftp_dir'])
        f = open(filename, 'wb')
        ftp.retrbinary('RETR ' + __model_parameters[model_id]['filename'], f.write)
        f.close()
    wl = wl_grid * 1.0
    h = fits.open(filename)
    flux = h[0].data
    h.close()

    # Trim the spectrum on both sides to make sure we can do redshift corrections
    wl_range = [np.min(wl * (1 + settings['virtual_dof']['redshift'][1] * 1e3 / scp.constants.c)), np.max(wl * (1 + settings['virtual_dof']['redshift'][0] * 1e3 / scp.constants.c))]
    mask_left = wl < wl_range[0]; mask_right = wl > wl_range[1]; mask_in = (~mask_left) & (~mask_right)
    meta = {'left': [wl[mask_left], flux[mask_left]], 'right': [wl[mask_right], flux[mask_right]]}
    return wl[mask_in], flux[mask_in], meta

def preprocess_grid_model(wl, flux, params, meta):
    # Restore the full (untrimmed) spectrum
    wl_full = np.concatenate([meta['left'][0], wl, meta['right'][0]])
    flux_full = np.concatenate([meta['left'][1], flux, meta['right'][1]])

    # Apply the redshift
    wl_redshifted = wl_full * (1 + params['redshift'] * 1e3 / scp.constants.c)

    # Re-interpolate back into the original wavelength grid
    flux = np.interp(wl, wl_redshifted, flux_full)
    return flux

