############################################################
#                                                          #
#                      GRIDFIT PRESET                      #
#                                                          #
#  This preset is used to measure the stellar parameters   #
#  of the spectrum (Teff, log(g), [M/H], [a/M] and [C/M])  #
#  by fitting a model to the spectrum obtained by          #
#  interpolating the model grid                            #
#                                                          #
############################################################

import h5py
import pickle
import scipy as scp

if 'grid_filename' not in original_settings:
    raise ValueError('You must create a local settings preset (settings/local/gridfit.py) that defines the "grid_filename" key to point to the model grid HDF5 file')



# This is the "defective" mask that removes parts of the spectrum where the models do not match the spectra of 
# Arcturus and the Sun well
defective_mask = [[3800, 4000], [4006, 4012], [4065, 4075], [4093, 4110], [4140, 4165], [4170, 4180], [4205, 4220],
                 [4285, 4300], [4335, 4345], [4375, 4387], [4700, 4715], [4775, 4790], [4855, 4865], [5055, 5065],
                 [5145, 5160], [5203, 5213], [5885, 5900], [6355, 6365], [6555, 6570], [7175, 7195], [7890, 7900],
                 [8320, 8330], [8490, 8505], [8530, 8555], [8650, 8672]]

settings = {
    ### Model grid settings ###
    'grid_filename': original_settings['grid_filename'],    # Model directory must be specified in local settings

    ### Spacing between trial redshifts when redshift search is requested in public__gridfit() ###
    'redshift_search_spacing': 50,

    ### Maximum absolute shifts in best-fit parameters between iterations at convergence ###
    'conv_shifts': {'teff': 1, 'logg': 0.01, 'zscale': 0.01, 'alpha': 0.01, 'carbon': 0.01},

    ### Which parameters to fit for? All dimensions of the grid + redshift ###
    'fit_dof': ['zscale', 'alpha', 'teff', 'logg', 'carbon', 'redshift'],

    ### Bounds for virtual parameters, i.e. parameters that are not dimensions of the grid ###
    'virtual_dof': {'redshift': [-500, 500]},

    ### Default initial guesses ###
    'default_initial': {'redshift': 0.0},

    ### Update the fitting mask ###
    'masks': apply_standard_mask(defective_mask, original_settings['masks']),

    # Print notifications?
    'silent': False,
}

# Reference to chemfit.chemfit() which will be populated when this preset is initialized
def main__chemfit():
    pass

_global = {}

def notify(message, color = 'k'):
    """Print a status notification
    
    This is a wrapper around `print()`, which makes it suppressible by `settings['silent']`. We also allow color printing, with the
    default color ('k') meaning "uncolored"
    
    Parameters
    ----------
    message : str
        Status message to display
    color : str, optional
        Color of the message. Accepted values are 'k' for uncolored, 'r' for red, 'g' for green, 'y' for yellow, 'b' for blue, 'm' for
        magenta and 'c' for cyan
    """
    prefix = {'k': '', 'r': '\033[31m', 'g': '\033[32m', 'y': '\033[33m', 'b': '\033[34m', 'm': '\033[35m', 'c': '\033[36m'}[color]
    suffix = ['', '\033[0m'][color != 'k']
    if not settings['silent']:
        print(prefix + message + suffix, flush = True)

def read_grid_dimensions():
    """Determine the available dimensions in the model grid and the grid points
    available in those dimensions
    
    For HDF5 model grids, the dimensions are stored in the pickled header dataset
    
    Returns
    -------
    dict
        Dictionary of lists, keyed be grid axis names. The lists contain unique
        values along the corresponding axis that are available in the grid
    """
    global _global

    # Load the grid from the HDF5 header
    with h5py.File(settings['grid_filename'], 'r') as f:
        header = pickle.loads(bytes(f['header'][()]))
    grid = {'teff': header['teff'], 'logg': header['logg'], 'zscale': header['zscale'], 'alpha': header['alpha'], 'carbon': header['carbon']}
    for key in grid:
        grid[key] = np.array(grid[key])

    # Prepare a map of grid points to indices for quick lookup
    _global['index_map'] = {param: {value: i for i, value in enumerate(grid[param])} for param in grid}

    # Load the wavelength grid
    _global['wl'] = header['wl']

    return grid

def read_grid_model(params, grid):
    """Load a specific model spectrum from the model grid

    In order to handle redshift, the function will trim the wavelength range of the model spectrum
    on both sides to make sure that the resulting wavelength range remains within the model coverage
    at all redshifts between the bounds in `settings['virtual_dof']['redshift']`. The trimmed parts
    of the spectrum are provided as additional model data in `meta` and may be used by the
    preprocessor to apply the redshift correction
    
    Parameters
    ----------
    params : dict
        Dictionary of model parameters. A value must be provided for each grid
        axis, keyed by the axis name
    grid   : dict
        Model grid dimensions, previously obtained with `read_grid_dimensions()`
    
    Returns
    -------
    wl : array_like
        Grid of model wavelengths in A
    flux : array_like
        Corresponding flux densities
    meta : dict
        Dictionary with trimmed parts of the spectrum for redshift calculations
    """
    global _global

    indices = tuple([np.where(grid[param] == params[param])[0][0] for param in ['teff', 'logg', 'zscale', 'alpha', 'carbon']])

    wl = _global['wl'] * 1.0
    with h5py.File(settings['grid_filename'], 'r') as f:
        flux = f['cont'][indices] * f['line'][indices]

    # Trim the spectrum on both sides to make sure we can do redshift corrections
    wl_range = [np.min(wl * (1 + settings['virtual_dof']['redshift'][1] * 1e3 / scp.constants.c)), np.max(wl * (1 + settings['virtual_dof']['redshift'][0] * 1e3 / scp.constants.c))]
    mask_left = wl < wl_range[0]; mask_right = wl > wl_range[1]; mask_in = (~mask_left) & (~mask_right)
    meta = {'left': [wl[mask_left], flux[mask_left]], 'right': [wl[mask_right], flux[mask_right]]}

    return wl[mask_in], flux[mask_in], meta

def preprocess_grid_model(wl, flux, params, meta):
    """Apply redshift correction

    Parameters
    ----------
    wl : array_like
        Grid of model wavelengths in A (trimmed to accommodate all redshifts)
    flux : array_like
        Corresponding flux densities
    params : dict
        Parameters of the model, including desired redshift
    meta : dict
        Trimmed parts of the spectrum
    
    Returns
    -------
    array_like
        Redshifted flux density
    """
    # Restore the full (untrimmed) spectrum
    wl_full = np.concatenate([meta['left'][0], wl, meta['right'][0]])
    flux_full = np.concatenate([meta['left'][1], flux, meta['right'][1]])

    # Apply the redshift
    wl_redshifted = wl_full * (1 + params['redshift'] * 1e3 / scp.constants.c)

    # Re-interpolate back into the original wavelength grid
    flux = np.interp(wl, wl_redshifted, flux_full)
    return flux

def public__get_carbon_wrt_solar(params):
    """Convert 'carbon' in the best-fit parameters dictionary calculated with `chemfit.chemfit()`
    to the solar scale

    Parameters
    ----------
    params : dict
        Dictionary of best-fit parameter values
    
    Returns
    -------
    number
        [C/M] re-expressed on the scale, such that [C/M] represents solar abundance
    """
    # Load the grid from the HDF5 header
    with h5py.File(settings['grid_filename'], 'r') as f:
        header = pickle.loads(bytes(f['header'][()]))

    # Convert carbon_map from individual points to regular grid interpolator
    zscale_grid = sorted(set(x for x, y in header['carbon_map']))
    logg_grid = sorted(set(y for x, y in header['carbon_map']))
    carbon_map = np.zeros([len(zscale_grid), len(logg_grid)])
    for i, x in enumerate(zscale_grid):
        for j, y in enumerate(logg_grid):
            carbon_map[i,j] = header['carbon_map'][(x, y)]
    carbon_map = scp.interpolate.RegularGridInterpolator([zscale_grid, logg_grid], carbon_map)
    
    # Run the conversion
    return params['carbon'] + carbon_map([params['zscale'], params['logg']])[0]

def public__gridfit(wl, flux, ivar, initial = {}, phot = {}, separate_redshift = True, search_redshift = False, max_iter = 20):
    """Determine the stellar parameters and redshift (radial velocity) of a star from its spectrum and
    (optionally) its photometry

    When photometry is provided, the fit is carried out iteratively such that the effective temperature
    is determined purely from photometry, and all other parameters are determined purely from the
    spectrum. The function also provides the ability to assign distinct redshifts to individual arms of
    the spectrograph to handle inconsistent wavelength calibration between the arms

    Spectral fitting tends to be sensitive to the initial redshift guess. The function provides the
    ability to run the preliminary fit multiple times with different initial redshifts in order to
    ensure convergence

    Parameters
    ----------
    wl : dict
        Spectrum wavelengths keyed by spectrograph arm
    flux : dict
        Spectrum flux densities keyed by spectrograph arm
    ivar : dict
        Spectrum weights (inverted variances) keyed by spectrograph arm
    initial : dict
        Initial guesses for the stellar parameters, keyed by parameter. Default initial guesses will
        be used for unspecified parameters (see `default`). If `search_redshift` is `False`, then an
        initial guess of redshift is mandatory. Note that statistical priors are not supported by `gridfit()`
        (use `chemfit()` instead)
    phot : dict, optional
        Photometric colors of the star. Each color is keyed by `BAND1#BAND2`, where `BAND1` and `BAND2`
        are the transmission profile filenames of the filters, as required by `synphot()`. Each element
        is a 2-element tuple, where the first element is the measured color, and the second element is the
        uncertainty in the measurement. The dictionary may also include optional elements `reddening`
        (E(B-V), single numerical value), and `mag_system` (one of the magnitude systems supported by
        `synphot()`, single string)
    separate_redshift : bool, optional
        If enabled, separate redshifts will be determined for each arm of the spectrograph. Otherwise,
        all arms will be forced to have the same redshift
    search_redshift : bool, optional
        If enabled, the function will try a multitude of initial redshift guesses across the entire allowed
        range (given by `settings['virtual_dof']['redshift']`) in steps given by
        `settings['redshift_search_spacing']`. While this approach is more likely to produce the best final
        result, it can be computationally demanding. Alternatively, this option may be disabled, in which
        case a good redshift guess must be provided in `initial`
    max_iter : number, optional
        Maximum number of refinement iterations, which are used to process photometric priors and/or separate
        arm redshifts. The fitter may carry out fewer iterations than this maximum value if good convergence
        of parameters is attained

    Returns
    -------
    dict
        The format of the output is similar to `chemfit.chemfit()`, except the best-fit redshift values
        (`['fit']['redshift']`) are provided as a dictionary keyed by spectrograph arms instead of a singular
        value. If `separate_redshift` is `True`, the priliminary best-fit parameters obtained from the spectra
        of individual arms are stored in `['extra']['fit_individual']`. The final rest frame wavelength vectors
        (keyed by spectrograph arms) are provided in `['extra']['rest_wl']`
    """
    extra = {}

    # Process initial guesses
    default = {'teff': 4000, 'logg': 1.5, 'zscale': -1.0, 'alpha': 0.0, 'carbon': 0.0}
    initial = {**default, **initial}
    if not search_redshift:
        if 'redshift' not in initial:
            raise ValueError('Initial redshift must be either specified or the redshift search must be enabled (gridfit(..., search_redshift=True))')
        redshift_trials = [initial['redshift']]
    else:
        if 'redshift' not in settings['fit_dof']:
            raise ValueError('redshift must be included in settings[\'fit_dof\'] to use gridfit(..., search_redshift=True)')
        redshift_trials = np.arange(*settings['virtual_dof']['redshift'], settings['redshift_search_spacing'])[1:]
    if np.any([np.ndim(initial[param]) > 0 for param in initial]):
        raise ValueError('Statistical priors are not supported in gridfit()')

    if len(phot) > 0 and 'teff' not in settings['fit_dof']:
        raise ValueError('teff must be included in settings[\'fit_dof\'] to use photometric priors')

    # Get preliminary redshifts and stellar parameters as a starting point for iterative refinements
    if separate_redshift:
        if 'redshift' not in settings['fit_dof']:
            raise ValueError('redshift must be included in settings[\'fit_dof\'] to use gridfit(..., separate_redshift=True)')
        notify('Calculating preliminary redshifts for individual arms...')
        extra['fit_individual'] = {}
        redshift = {}
        for arm in wl:
            redshift[arm] = (0.0, np.inf)
            for redshift_trial in redshift_trials:
                fit = main__chemfit({arm: wl[arm]}, {arm: flux[arm]}, {arm: ivar[arm]}, initial = {**initial, 'redshift': redshift_trial}, phot = phot)
                if 'redshift' not in fit['errors']:
                    redshift[arm] = (fit['fit']['redshift'], np.nan)
                    extra['fit_individual'][arm] = {'fit': fit['fit'], 'errors': fit['errors']}
                elif fit['errors']['redshift'] < redshift[arm][1]:
                    redshift[arm] = (fit['fit']['redshift'], fit['errors']['redshift'])
                    extra['fit_individual'][arm] = {'fit': fit['fit'], 'errors': fit['errors']}
            notify('Redshift in {} = {:.3f} ± {:.3f}'.format(arm, *redshift[arm]), color = 'g')
        notify('Calculating preliminary combined fit...')
        rest_wl = {arm: wl[arm] / (1 + redshift[arm][0] * 1e3 / scp.constants.c) for arm in wl}
        fit = main__chemfit(rest_wl, flux, ivar, initial = {**initial, 'redshift': 0.0}, phot = phot, dof = list(set(settings['fit_dof']) - {'redshift'}))
        result = {param: fit['fit'][param] for param in default}
    else:
        notify('Calculating preliminary fit...')
        redshift = (0.0, np.inf)
        for redshift_trial in redshift_trials:
            fit = main__chemfit(wl, flux, ivar, initial = {**initial, 'redshift': redshift_trial}, phot = phot)
            if 'redshift' not in fit['errors']:
                redshift = (fit['fit']['redshift'], np.nan)
                result = {param: fit['fit'][param] for param in default}
            elif fit['errors']['redshift'] < redshift[1]:
                redshift = (fit['fit']['redshift'], fit['errors']['redshift'])
                result = {param: fit['fit'][param] for param in default}
        notify('Redshift = {:.3f} ± {:.3f}'.format(*redshift), color = 'g')
        redshift = {arm: redshift for arm in wl}
    notify('Preliminary fit computed: {}'.format({key: np.round(result[key], 3) for key in result}), color = 'g')

    # Blank flux for purely photometric fits
    blank_flux = {arm: np.full(len(flux[arm]), np.nan) for arm in flux}

    is_phot_step = len(phot) > 0 and 'teff' in settings['fit_dof'] # Do we need photometric steps
    # Which parameters to fit for in the spectroscopic step
    spec_step_params = set(settings['fit_dof'])
    if is_phot_step:
        spec_step_params -= {'teff'}
    if separate_redshift:
        spec_step_params -= {'redshift'}
    spec_step_params = list(spec_step_params)

    # Carry out iterative refinements
    for n_iter in range(max_iter):
        notify('Starting iteration {}'.format(n_iter + 1))

        rest_wl = {arm: wl[arm] / (1 + redshift[arm][0] * 1e3 / scp.constants.c) for arm in wl}

        prev = copy.deepcopy(result)

        # Photometric step (teff only)
        if is_phot_step:
            fit = main__chemfit(rest_wl, blank_flux, ivar, initial = {**result, 'redshift': 0}, phot = phot, dof = ['teff'])
            result = {param: fit['fit'][param] for param in default}
            notify('Photometric step: {}, redshift: {}'.format({key: np.round(result[key], 3) for key in result}, {key: np.round(redshift[key][0], 3) for key in redshift}), color = 'c')

        # Spectroscopic step
        fit = main__chemfit(rest_wl, flux, ivar, initial = {**result, 'redshift': 0}, dof = spec_step_params)
        result = {param: fit['fit'][param] for param in default}
        if 'redshift' in fit['errors']:
            redshift = {arm: (redshift[arm][0] + fit['fit']['redshift'], fit['errors']['redshift']) for arm in wl}
        notify('Spectroscopic step: {}, redshift: {}'.format({key: np.round(result[key], 3) for key in result}, {key: np.round(redshift[key][0], 3) for key in redshift}), color = 'c')

        # Redshift step
        if separate_redshift:
            for arm in wl:
                fit = main__chemfit({arm: rest_wl[arm]}, {arm: flux[arm]}, {arm: ivar[arm]}, initial = {**result, 'redshift': 0}, dof = ['redshift'])
                redshift[arm] = (redshift[arm][0] + fit['fit']['redshift'], fit['errors']['redshift'])
            notify('Redshift step: {}, redshift: {}'.format({key: np.round(result[key], 3) for key in result}, {key: np.round(redshift[key][0], 3) for key in redshift}), color = 'c')

        # Evaluate convergence
        if np.all([np.abs(prev[param] - result[param]) <= settings['conv_shifts'][param] for param in settings['conv_shifts']]):
            notify('Stellar parameters have converged', color = 'g')
            break
        elif n_iter == max_iter - 1:
            notify('Maximum number of iterations reached', color = 'r')

    # Get the final uncertainties and the fit object
    notify('Calculating covariance matrix...')
    rest_wl = {arm: wl[arm] / (1 + redshift[arm][0] * 1e3 / scp.constants.c) for arm in wl}
    extra['rest_wl'] = rest_wl
    fit = main__chemfit(rest_wl, flux, ivar, initial = {**result, 'redshift': 0}, phot = phot, method = 'cov')
    fit['fit']['redshift'] = {arm: redshift[arm][0] for arm in redshift}
    fit['extra'].update(extra)

    return fit
