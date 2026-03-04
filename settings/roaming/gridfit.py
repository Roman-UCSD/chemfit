############################################################
#                                                          #
#                      GRIDFIT PRESET                      #
#                                                          #
#  This preset is used to measure the stellar parameters   #
#  of the spectrum (Teff, log(g), [M/H], [a/M] and [C/M])  #
#  by fitting a model obtained by linear interpolation of  #
#  a precomputed model grid                                #
#                                                          #
#  This preset also supports photometric priors and        #
#  transforms the input spectrum into the rest frame of    #
#  reference                                               #
#                                                          #
############################################################

import h5py
import pickle
import scipy as scp
import gc, ctypes

if 'grid_filename' not in original_settings:
    raise ValueError('You must create a local settings preset (settings/local/gridfit.py) that defines the "grid_filename" key to point to the model grid HDF5 file')

settings = {
    ### Model grid settings ###
    'grid_filename': original_settings['grid_filename'],    # Model directory must be specified in local settings

    ### Spacing between trial redshifts when redshift search is requested in gridfit() ###
    'redshift_search_spacing': 50,

    ### Maximum absolute shifts in best-fit parameters between iterations at convergence ###
    'conv_shifts': {'teff': 1, 'logg': 0.01, 'zscale': 0.01, 'alpha': 0.01, 'carbon': 0.01},

    ### Which parameters to fit for? All dimensions of the grid + redshift ###
    'fit_dof': ['zscale', 'alpha', 'teff', 'logg', 'carbon', 'redshift'],

    ### Bounds for virtual parameters, i.e. parameters that are not dimensions of the model grid ###
    'virtual_dof': {'redshift': [-500, 500]},

    ### Default initial guesses. Note that this parameter only applies to `chemfit.chemfit()`. GRIDFIT has its ###
    ### own handling of initial values ###
    'default_initial': {'redshift': 0.0},

    ### Models are loaded in absolute flux density units if `True`, or as continuum-normalized spectra if `False` ###
    'use_model_continuum': True,

    ### Fitting masks in the rest and laboratory frames of reference ###
    'masks': {
        # In this rest frame, we use the "defective" mask that removes parts of the spectrum where the models do not
        # match the spectra of Arcturus and the Sun well
        'rest': [[100, 4000], [4006, 4012], [4065, 4075], [4093, 4110], [4140, 4165], [4170, 4180], [4205, 4220],
                 [4285, 4300], [4335, 4345], [4375, 4387], [4700, 4715], [4775, 4790], [4855, 4865], [5055, 5065],
                 [5145, 5160], [5203, 5213], [5885, 5900], [6355, 6365], [6555, 6570], [7175, 7195], [7890, 7900],
                 [8320, 8330], [8490, 8505], [8530, 8555], [8650, 8672]],

        # In the lab frame, we use the "aggressive" telluric mask that attempts to completely remove all regions affected
        # by telluric absoprption
        'lab': [[6270, 6330], [6860, 6970], [7150, 7400], [7590, 7715], [8100, 8380], [8915, 9910], [10730, 12300], [12450, 12900]],
    },

    ### Print notifications? ###
    'silent': False,
}

# References to chemfit.chemfit(), chemfit.flush_convolution_weights_cache(), chemfit.compute_covariance(),
# chemfit.estimate_continuum(), chemfit.ranges_to_mask() and chemfit.ModelGridInterpolator() which will be populated when
# this preset is initialized
def main__chemfit():
    pass
def main__flush_convolution_weights_cache():
    pass
def main__compute_covariance():
    pass
def main__estimate_continuum():
    pass
def main__ranges_to_mask():
    pass
class main__ModelGridInterpolator:
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

def free_memory():
    """Helper function to de-allocate memory in between chemfit runs

    This function both dispatches the Python garbage collector, and also runs glibc's `malloc_trim()` routine, as it is sometimes
    necessary to fully release garbage-collected memory
    """
    gc.collect()
    libc = ctypes.CDLL('libc.so.6')
    libc.malloc_trim(0)
    return

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

    # Check if the grid has already been loaded
    if ('grid_filename' in _global) and (_global['grid_filename'] == settings['grid_filename']):
        return _global['grid']

    # Load the grid from the HDF5 header
    with h5py.File(settings['grid_filename'], 'r') as f:
        header = pickle.loads(bytes(f['header'][()]))
    _global['grid_filename'] = settings['grid_filename']
    _global['grid'] = {'teff': header['teff'], 'logg': header['logg'], 'zscale': header['zscale'], 'alpha': header['alpha'], 'carbon': header['carbon']}
    for key in _global['grid']:
       _global['grid'][key] = np.array(_global['grid'][key])

    # Prepare a map of grid points to indices for quick lookup
    _global['index_map'] = {param: {value: i for i, value in enumerate(_global['grid'][param])} for param in _global['grid']}

    # Load the wavelength grid
    _global['wl'] = header['wl']

    return _global['grid']

def read_grid_model(params, grid):
    """Load a specific model spectrum from the model grid

    In order to handle redshift, the function will trim the wavelength range of the model spectrum
    on both sides to make sure that the resulting wavelength range remains within the model coverage
    at all redshifts between the bounds in `settings['virtual_dof']['redshift']`. The full spectrum
    will instead be provided as the meta data. When `preprocess_grid_model() preprocesses the spectrum,
    it will interpolate the full spectrum onto the appropriately redshifted wavelength grid

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
        Grid of model wavelengths in A, trimmed on both sides to accommodate all allowed redshifts
    flux : array_like
        Normally this array would contain the corresponding flux densities; however, the fitter will
        never use this flux array (it will instead use the interpolated flux from `preprocess_grid_model()`).
        We therefore return a blank array of appropriate size to save memory
    meta : dict
        Dictionary with the wavelength and flux density of the full (untrimmed) spectrum
    """
    global _global

    indices = tuple([np.where(grid[param] == params[param])[0][0] for param in ['teff', 'logg', 'zscale', 'alpha', 'carbon']])

    wl = _global['wl'] * 1.0
    with h5py.File(settings['grid_filename'], 'r') as f:
        if settings['use_model_continuum']:
            flux = f['cont'][indices] * f['line'][indices]
        else:
            flux = f['line'][indices]

    # Trim the spectrum on both sides to make sure we can do redshift corrections
    wl_range = [np.min(wl * (1 + settings['virtual_dof']['redshift'][1] * 1e3 / scp.constants.c)), np.max(wl * (1 + settings['virtual_dof']['redshift'][0] * 1e3 / scp.constants.c))]
    mask_in = (wl >= wl_range[0]) & (wl <= wl_range[1])
    meta = {'wl_full': wl, 'flux_full': flux}

    return wl[mask_in], np.broadcast_to(np.nan, (np.count_nonzero(mask_in),)), meta

def delete_grid_model(meta):
    """This function does nothing

    Normally this function is used to safely delete model data; however, in the GRIDFIT
    mode, the models do not contain any data that will not be deleted by the garbage
    collector, so we do not do anything here

    Parameters
    ----------
    meta : dict
        Dictionary with additional model data, as loaded by `read_grid_model()`
    """
    return

def preprocess_grid_model(wl, flux, params, meta):
    """Apply redshift correction

    Parameters
    ----------
    wl : array_like
        Grid of model wavelengths in A (trimmed to accommodate all redshifts)
    flux : array_like
        Corresponding flux densities (not used)
    params : dict
        Parameters of the model, including desired redshift
    meta : dict
        Trimmed parts of the spectrum

    Returns
    -------
    array_like
        Redshifted flux density
    """
    # Apply the redshift
    wl_redshifted = meta['wl_full'] * (1 + params['redshift'] * 1e3 / scp.constants.c)

    # Re-interpolate back into the original wavelength grid
    return np.interp(wl, wl_redshifted, meta['flux_full'])

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
    (optionally) photometry

    When photometry is provided, the fit is carried out iteratively such that the effective temperature
    is determined purely from photometry, and all other parameters are determined purely from the
    spectrum. This strategy reduces the sensitivity of the Teff estimate (and all other stellar parameters
    that depend on it) to the abundances of individual elements and NLTE effects. The function also provides
    the ability to assign distinct redshifts to individual arms of the spectrograph to handle inconsistent
    wavelength calibration between the arms

    Spectral fitting tends to be sensitive to the initial redshift guess. The function provides the
    ability to run the preliminary fit multiple times with different initial redshifts in order to
    find the best initial guess that maximizes the likelihood of good convergence

    Parameters
    ----------
    wl : dict of array_like
        Observed spectrum wavelengths keyed by spectrograph arm. Must be provided in A, air, laboratory
        (i.e., observed, telluric) frame of reference
    flux : dict of array_like
        Spectrum flux densities per unit wavelength, keyed by spectrograph arm. The continuum level in the
        spectrum will be ignored by the fitter, so the flux vectors can be provided in any units and with any
        normalization
    ivar : dict of array_like
        Spectrum weights (inverse variances) keyed by spectrograph arm. The fitter will attempt to rescale
        the weights to match the observed scatter in the data with respect to the best-fit model, so the
        exact normalization of the inverse variance arrays is less important (see `compute_covariance()`)
    initial : dict
        Initial guesses for the stellar parameters, keyed by parameter. Default initial guesses will
        be used for unspecified parameters (see `default`). If `search_redshift` is `False`, then an
        initial guess of redshift is mandatory. Note that statistical priors are not supported by `gridfit()`
        (use `chemfit()` instead)
    phot : dict, optional
        Photometric color(s) of the star. Each color is keyed by `BAND1#BAND2`, where `BAND1` and `BAND2`
        are the transmission profile filenames of the filters, stored in `settings['filter_dir']`
        (e.g. 'PAN-STARRS_PS1_g.dat#PAN-STARRS_PS1_i.dat'). Each element is a 2-element tuple, where the
        first element is the measured color, and the second element is the 1-sigma uncertainty in the
        measurement. The dictionary may also include optional keys `reddening` (E(B-V), single numerical
        value), and `mag_system` (one of the magnitude systems supported by `synphot()`, single string). If
        omitted, the default values will be assumed (`settings['default_reddening']` and 
        `settings['default_mag_system']`)
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
        of parameters is attained (where "good" is defined according to `settings['conv_shifts']`)

    Returns
    -------
    dict
        Fitting results. This is a multi-level dictionary with the following keys:

        - [fit]: Inferred values of best-fit stellar parameters and redshift. The best-fit redshift value is
                 separated by spectrograph arms (they would be identical if `separate_redshift` is `False`)
        - [errors]: Errors in the best-fit parameters
        - [extra][observed][wl]: Observed spectrum wavelengths in A, in air, transformed to rest frame
        - [extra][observed][flux]: Observed flux densities, as provided in `flux`
        - [extra][observed][ivar]: Inverse variances of the pixels of the observed spectrum, as in `ivar`
        - [extra][model][wl]: Best-fit model wavelengths, identical to [extra][observed][wl]
        - [extra][model][flux]: Best-fit model flux in erg / s / cm^2 / A / sr
        - [extra][model][cont]: Continuum correction between the observed flux vector and the model
        - [extra][model][line]: Line absorption coefficient of each pixel in the model spectrum
        - [extra][fit_individual]: Best-fit stellar parameters and their errors from the preliminary fits to
                                   individual arms of the spectrograph. Returned only if `separate_redshift` is
                                   `True`
        - [localfit]: Recommended input to LOCALFIT (run LOCALFIT as `localfit(**fit['extra']['localfit'])`).
                      The recommended input will propagate the flux and inverse variance arrays to LOCALFIT
                      unchanged, transform the wavelength arrays to rest frame (which is necessary, as LOCALFIT
                      does not support redshift), propagate the best-fit stellar parameters, propagate the inferred
                      photometric error on Teff (if photometry is provided), define the fitting mask with the lab
                      component transformed to rest frame, and make GRIDFIT's Jacobian available to LOCALFIT, so
                      that LOCALFIT could propagate the errors in stellar parameters into abundances
        - [extra][mask]: Fitting mask, which evaluates to `True` for all spectral pixels that were used in the fit
        - [extra][cost]: Chi-squared of the spectral pixels, (observed - model)^2 / error^2
        - [extra][arm_index]: The arm index of each pixel in [extra][observed] and [extra][model]. The arm index is an
                              integer that corresponds to the serial number of the arm that the pixel belongs to. The
                              arms are indexed in the alphabetical order starting with 0
        - [extra][lab_mask]: Mask that evaluates to `True` for all spectral pixels excluded due to the laboratory
                             frame mask, `settings['masks']['lab']` (usually telluric absorption)
        - [extra][gridfit_n_iter]: Total number of refinement iterations carried out by GRIDFIT
        - [extra][dof]: List of degrees of freedom included in the covariance matrix and the Jacobian
        - [extra][wl_cov]: Wavelengths of the spectral pixels included in the Jacobian ([extra][observed][wl] but with
                           bad and masked out pixels removed)
        - [extra][cov]: Covariance matrix, which is an array of shape NxN where N is the number of degrees of freedom
                        listed in [extra][dof]
        - [extra][jac]: Jacobian matrix, which is an array of shape MxN where N is the number of degrees of freedom
                        listed in [extra][dof], and M is the number of spectral pixels listed in [extra][wl_cov]
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

    # Set up the initial fitting mask
    settings['mask'] = settings['masks']['rest'] + settings['masks']['lab']

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
        main__flush_convolution_weights_cache() # Remove the calculated convolution weights as we no longer need them since the wavelength grid changed
        fit = main__chemfit(rest_wl, flux, ivar, initial = {**initial, 'redshift': 0.0}, phot = phot, dof = sorted(list(set(settings['fit_dof']) - {'redshift'})))
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
    spec_step_params = sorted(list(spec_step_params))

    # Carry out iterative refinements
    for n_iter in range(max_iter):
        notify('Starting iteration {}'.format(n_iter + 1))

        best_redshift = [redshift[arm][0] for arm in wl][np.argmin([redshift[arm][1] for arm in wl])]
        settings['mask'] = settings['masks']['rest'] + (np.array(settings['masks']['lab']) / (1 + best_redshift * 1e3 / scp.constants.c)).tolist()
        rest_wl = {arm: wl[arm] / (1 + redshift[arm][0] * 1e3 / scp.constants.c) for arm in wl}
        main__flush_convolution_weights_cache()

        prev = copy.deepcopy(result)

        # Photometric step (teff only)
        if is_phot_step:
            fit = main__chemfit(wl, blank_flux, ivar, initial = {**result, 'redshift': best_redshift}, phot = phot, dof = ['teff'])
            result = {param: fit['fit'][param] for param in default}
            e_teff = fit['errors']['teff']
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

    # Evaluate the final uncertainties in best-fit parameters
    notify('Calculating covariance matrix...')
    best_redshift = [redshift[arm][0] for arm in wl][np.argmin([redshift[arm][1] for arm in wl])]
    settings['mask'] = settings['masks']['rest'] + (np.array(settings['masks']['lab']) / (1 + best_redshift * 1e3 / scp.constants.c)).tolist()
    rest_wl = {arm: wl[arm] / (1 + redshift[arm][0] * 1e3 / scp.constants.c) for arm in wl}
    main__flush_convolution_weights_cache()
    fit = main__chemfit(rest_wl, flux, ivar, initial = {**result, 'redshift': 0}, method = 'jac')
    # If photometric teff is used, propagate the uncertainty in teff correctly into the final covariance matrix
    if is_phot_step:
        best_fit = [fit['fit'][param] for param in fit['extra']['fit']['dof']]
        fixed_param_errors = np.full(len(best_fit), np.nan)
        fixed_param_errors[fit['extra']['fit']['dof'].tolist().index('teff')] = e_teff
        fit['extra']['cov'] = main__compute_covariance(fit['extra']['jac'], fit['extra']['fit']['x'], fit['extra']['fit']['y'], fit['extra']['fit']['sigma'], fit['extra']['fit']['f'], best_fit, fixed_param_errors = fixed_param_errors)
        fit['errors'] = {param: np.sqrt(fit['extra']['cov'][i,i]) for i, param in enumerate(fit['extra']['fit']['dof'])}

    # Re-organize the output for easy analysis of data quality
    use_model_continuum = settings['use_model_continuum']
    mask = copy.deepcopy(settings['mask'])
    try:
        settings['use_model_continuum'] = False
        settings['mask'] = []

        # Refit the continuum without the fitting mask to fill in the masked gaps
        cont_full = main__estimate_continuum(fit['extra']['observed']['wl'], fit['extra']['observed']['flux'] / fit['extra']['model']['flux'], fit['extra']['observed']['ivar'] * fit['extra']['model']['flux'] ** 2, arm_index = fit['extra']['arm_index'])
        fit['extra']['model']['cont'][~fit['extra']['mask']] = cont_full[~fit['extra']['mask']]

        # Retrieve the continuum-normalized model spectrum
        interpolator = main__ModelGridInterpolator(detector_wl = rest_wl)
        fit['extra']['model']['line'] = interpolator({**fit['fit'], 'redshift': 0})[1]

        # Compute chi-squared
        fit['extra']['cost'] = (fit['extra']['observed']['flux'] - fit['extra']['model']['flux'] * fit['extra']['model']['cont']) ** 2.0 * fit['extra']['observed']['ivar']

        # Make the telluric mask available to the user
        fit['extra']['lab_mask'] = main__ranges_to_mask(fit['extra']['observed']['wl'], np.array(settings['masks']['lab']) / (1 + best_redshift * 1e3 / scp.constants.c), True)

    finally:
        settings['use_model_continuum'] = use_model_continuum
        settings['mask'] = mask

    # Finalize the fit object
    extra['gridfit_n_iter'] = n_iter + 1
    fit['fit']['redshift'] = {arm: redshift[arm][0] for arm in redshift}
    fit['extra'].update(extra)
    del fit['interpolator_statistics'] # We do not keep interpolator statistics because they only represent the last call of chemfit() and not this entire analysis run
    fit['extra']['dof'] = fit['extra']['fit']['dof']
    fit['extra']['wl_cov'] = fit['extra']['fit']['x']
    del fit['extra']['fit']

    # Add LOCALFIT inputs to the fit object
    fit['localfit'] = {'wl': rest_wl, 'flux': flux, 'ivar': ivar, 'mask': copy.deepcopy(settings['mask']), 'jacobian': [fit['extra']['dof'], fit['extra']['jac']]}
    fit['localfit']['gridfit'] = {param: fit['fit'][param] for param in read_grid_dimensions()}
    if is_phot_step:
        fit['localfit']['gridfit']['teff'] = (fit['localfit']['gridfit']['teff'], fit['errors']['teff'])

    free_memory()

    return fit
