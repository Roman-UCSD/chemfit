import numpy as np
import scipy as scp
import pickle
import gzip
import os
import re
import hashlib
import itertools
import warnings
import copy
import inspect
import importlib.util

# scipy.optimize.curve_fit() throws an exception when the optimizer exhausts the maximum allowed number of function
# evaluations. Since we are demanding fairly narrow tolerance, reaching the maximum number of evaluations is not
# necessarily fatal (it merely means that the actual tolerance is somewhat less than our default standard). So we would
# like to deescalate the error to a warning and proceed. Unfortunately, I am unable to find a more elegant way of doing
# it except by overriding the optimizer function that curve_fit() calls
original_least_squares = scp.optimize._minpack_py.least_squares
def least_squares_wrapper(*args, **kwargs):
    res = original_least_squares(*args, **kwargs)
    if (not res.success) and (res.status == 0):
        warn('Optimizer failed to reach desired convergence after the maximum number ({}) of function evaluations'.format(res.nfev))
        res.success = True
    return res
scp.optimize._minpack_py.least_squares = least_squares_wrapper

script_dir = os.path.dirname(os.path.realpath(__file__))
settings = {}

warnings_stack = []
warnings_messages = {}

def warn(message):
    """Issue a warning. Wrapper for `warnings.warn()`
    
    Add the warning message to the stack. If required, also throw the normal Python warning

    Instead of storing warning messages in the stack in full, we associate unique numerical
    identifiers with each distinct message and store them in a dictionary. The stack then
    only contains the identifiers to save memory
    
    Parameters
    ----------
    message : str
        Warning message
    """
    global warnings_stack, warnings_messages

    if message not in warnings_messages:
        warning_id = len(warnings_messages)
        warnings_messages[message] = warning_id
    else:
        warning_id = warnings_messages[message]

    warnings_stack += [warning_id]
    if settings['throw_python_warnings']:
        warnings.warn(message)

def read_grid_model(params, grid):
    """Load a specific model spectrum from the model grid

    All models within the same model grid are expected to have the same wavelength
    sampling
    
    This function definition is a template. The actual implementation is deferred to the
    settings preset files
    
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
        Corresponding flux densities in wavelength space in arbitrary units
    meta : dict
        Dictionary with additional model data that will be made available to the model
        preprocessor
    """
    raise NotImplementedError()

def preprocess_grid_model(wl, flux, params, meta):
    """Preprocess a model spectrum

    The purpose of the preprocessor is to apply the effect of virtual degrees of freedom
    to the model spectrum. If additional model data are required to do that, they may be
    loaded by `read_grid_model()` and returned as the `meta` argument, which is made
    available to the preprocessor
    
    This function definition is a template. The actual implementation is deferred to the
    settings preset files
    
    Parameters
    ----------
    wl : array_like
        Grid of model wavelengths in A, as loaded by `read_grid_model()`
    flux : array_like
        Corresponding flux densities in wavelength space in arbitrary units, as loaded by
        `read_grid_model()`
    params : dict
        Parameters of the model, including both real and virtual degrees of freedom
    meta : dict
        Dictionary with additional model data, as loaded by `read_grid_model()`
    
    Returns
    -------
    flux : array_like
        Processed flux. The flux array must be sampled over the same wavelengths as the
        original model (i.e. the wavelength array in `wl`)
    """
    raise NotImplementedError()

def read_grid_dimensions(flush_cache = False):
    """Determine the available dimensions in the model grid and the grid points
    available in those dimensions
    
    This function definition is a template. The actual implementation is deferred to the
    settings preset files
    
    Parameters
    ----------
    flush_cache : bool, optional
        If True, discard cache and scan the grid afresh
    
    Returns
    -------
    dict
        Dictionary of lists, keyed be grid axis names. The lists contain unique
        values along the corresponding axis that are available in the grid
    """
    raise NotImplementedError()

def apply_standard_mask(exclude, original_mask = False):
    """Helper function to extend standard masks. Primarily intended to be used in settings files
    
    A "standard" fitting mask does not have any arm or parameter-specific components (i.e. only
    `settings['masks']['all']['all']` and `settings['masks']['continuum']` are populated, and the
    continuum mask is the inverse of the `['all']['all']` mask. We also require the first included
    interval to start at 100 A and the last included interval to end at 100000 A

    This function adds additional excluded regions to an existing standard mask, or creates a new
    standard mask from scratch if the existing (i.e. original) standard mask is not provided
    
    Parameters
    ----------
    exclude : list
        List of intervals to exclude from the mask
    original_mask : dict, optional
        Original fitting mask (`settings['masks']`) or `False` to create a new mask
    
    Returns
    -------
    dict
        Updated (or newly created) fitting mask that can be placed directly into `settings['masks']`
    """
    if np.min(exclude) <= 100 or np.max(exclude) >= 100000:
        raise ValueError('Mask wavelengths must be 100 < wl < 100000')

    # If an original mask is provided, check that it is standard
    if (type(original_mask) is not bool) and ((len(original_mask['all']['all']) > 0) or (len(original_mask['continuum']) > 0)):
        if len(original_mask['all']['all']) < 1 or original_mask['all']['all'][0][0] != 100 or original_mask['all']['all'][-1][-1] != 100000:
            raise ValueError('Cannot extend non-standard mask')
        original_exclude = np.array(original_mask['all']['all']).flatten()[1:-1]
        original_exclude = np.reshape(original_exclude, [len(original_exclude) // 2, 2]).tolist()
        if len(original_exclude) != len(original_mask['continuum']) or (not np.all(original_exclude == original_mask['continuum'])):
            raise ValueError('Cannot extend non-standard mask')

    masks = {'all': {'all': []}, 'continuum': []}
    if type(original_mask) is not bool:
        masks['continuum'] += original_mask['continuum']
    masks['continuum'] += list(copy.deepcopy(exclude))
    masks['continuum'].sort(key = lambda x: x[0])
    merged = []
    for interval in masks['continuum']:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            merged[-1][1] = max(merged[-1][1], interval[1])
    masks['continuum'] = merged
    masks['all']['all'] = np.concatenate([[100], np.array(masks['continuum']).flatten(), [100000]])
    assert np.all(np.diff(masks['all']['all']) > 0)
    masks['all']['all'] = np.reshape(masks['all']['all'], [len(masks['continuum']) + 1, 2]).tolist()
    return masks

def initialize(*presets):
    """Load the settings presets

    Instrument specifications, fitting parameters, model grid handling and other required data and procedures
    are stored in a collection of Python scripts referred to as presets. The scripts may be found in the
    "settings" directory. Each script is expected to define a `settings` dictionary, the entries of which will
    then be used to update the global `settings` dictionary. Optionally, the scripts may also define the
    `read_grid_model()` and `read_grid_dimensions()` functions, following the templates in this file. These
    functions will be used to load model spectra during the fitting process

    Note that the "default" preset will be loaded automatically

    Parameters
    ----------
    presets : tuple
        Settings presets to load. For each preset, the function will attempt to load both
        "settings/local/<preset>.py" and "settings/roaming/<preset>.py" in that order, if available. It is
        recommended to store global parameters in the former file (versioned) and machine-specific parameters in
        the latter file (not versioned, i.e. ignored)
    """
    global read_grid_model, read_grid_dimensions, preprocess_grid_model, settings

    # Reset settings
    settings = {}

    # Environment variables provided to the preset scripts
    env = {'script_dir': script_dir, 'np': np, 'original_settings': copy.deepcopy(settings), 'copy': copy, 'warn': warn, 'apply_standard_mask': apply_standard_mask}

    index = 0 # Load counter to ensure unique module names for all loaded files
    for preset in ['default'] + list(presets):
        scripts = ['{}/settings/local/{}.py'.format(script_dir, preset), '{}/settings/roaming/{}.py'.format(script_dir, preset)]
        scripts = [script for script in scripts if os.path.isfile(script)]
        if len(scripts) == 0:
            raise ValueError('Settings preset {} not found'.format(preset))
        for script in scripts:
            # Load and run the file
            index += 1
            spec = importlib.util.spec_from_file_location('settings_{}'.format(index), script)
            module = importlib.util.module_from_spec(spec)
            module.__dict__.update(env)
            spec.loader.exec_module(module)

            # Update the global settings variable and make it available to the module in case it defines model read
            # functions and they need it
            settings.update(module.settings)
            env['original_settings'] = copy.deepcopy(settings)
            module.settings = settings

            # Try to overwrite the model read functions with the ones stored in the module
            try:
                read_grid_model = module.read_grid_model
            except:
                pass
            try:
                read_grid_dimensions = module.read_grid_dimensions
            except:
                pass
            try:
                preprocess_grid_model = module.preprocess_grid_model
            except:
                pass

            # See if the configuration file has public functions that need to be exposed to the user
            public_functions = [name for name in dir(module) if name.startswith(prefix := 'public__')]
            for name in public_functions:
                globals()[name[len(prefix):]] = getattr(module, name)

            # See if the module has references for base chemfit functions that need to be populated
            references = [name for name in dir(module) if name.startswith(prefix := 'main__')]
            for name in references:
                setattr(module, name, globals()[name[len(prefix):]])

# Load the default settings preset
initialize()

def convolution_integral(sigma, segment_left, segment_right, bin_left, bin_right):
    """Calculate weights of the convolution for an arbitrary flux density spectrum
    and a Gaussian kernel with fixed or slowly-varying standard deviation
    
    A source with a given (true) flux density spectrum is observed with a detector
    of given resolution that records the average flux density in a set of wavelength
    bins. This function estimates the linear coefficients `C1` and `C2` such that

    ``C1 * A + C2 * B``

    gives the contribution of the spectrum segment between the wavelengths
    `segment_left` and `segment_right`, to the detector bin with edge wavelengths
    `bin_left` and `bin_right`. It is assumed that the resolution of the detector
    is constant for each detector bin / spectrum segment pair, and that the flux
    density within the segment varies linearly with wavelength:

    ``flux density = A + B * wavelength``

    This function is primarily designed to handle model spectra, for which the flux
    density has been calculated at given wavelength points, and can be linearly
    interpolated between those points. The weight coefficients may then be evaluated
    for each segment between adjacent wavelength points, and the contributions
    of each segment can be added together to obtain the total flux density in the
    detector bin
    
    Parameters
    ----------
    sigma : array_like
        Detector resolution for the flux arriving from the spectrum segment of
        interest at the detector bin of interest. Adopted as the standard deviation
        of the Gaussian convolution kernel
    segment_left : array_like
        Lower wavelength bound of the spectrum segment of interest
    segment_right : array_like
        Upper wavelength bound of the spectrum segment of interest
    bin_left : array_like
        Lower wavelength bound of the detector bin of interest
    bin_right : array_like
        Upper wavelength bound of the detector bin of interest
    
    Returns
    -------
    C1 : array_like
        Weight coefficient of the vertical offset in the flux density spectrum segment
    C2 : array_like
        Weight coefficient of the linear slope in the flux density spectrum segment
    """
    sqrt_sigma = np.sqrt(2) * sigma
    sqpi = np.sqrt(np.pi)
    sq2pi = np.sqrt(2) * sqpi

    blsl = (segment_left - bin_left) / sqrt_sigma; blsl_sigma = segment_left ** 2 - bin_left ** 2 + sigma ** 2
    blsr = (segment_left - bin_right) / sqrt_sigma; blsr_sigma = segment_left ** 2 - bin_right ** 2 + sigma ** 2
    brsl = (segment_right - bin_left) / sqrt_sigma; brsl_sigma = segment_right ** 2 - bin_left ** 2 + sigma ** 2
    brsr = (segment_right - bin_right) / sqrt_sigma; brsr_sigma = segment_right ** 2 - bin_right ** 2 + sigma ** 2

    erfblsl = scp.special.erf(blsl); expblsl = np.exp(-blsl ** 2)
    erfblsr = scp.special.erf(blsr); expblsr = np.exp(-blsr ** 2)
    erfbrsl = scp.special.erf(brsl); expbrsl = np.exp(-brsl ** 2)
    erfbrsr = scp.special.erf(brsr); expbrsr = np.exp(-brsr ** 2)

    x = lambda x, erfxsl, erfxsr, expxsl, expxsr: (expxsl - expxsr) * sqrt_sigma / sqpi + (x - bin_left) * erfxsl + (bin_right - x) * erfxsr
    x1 = x(segment_left, erfblsl, erfblsr, expblsl, expblsr)
    x2 = x(segment_right, erfbrsl, erfbrsr, expbrsl, expbrsr)

    x = lambda erfxsl, erfxsr, expxsl, expxsr: (expxsl - expxsr) * sigma / sq2pi + 0.5 * (erfxsl - erfxsr)
    x3 = x(blsl_sigma * erfblsl, blsr_sigma * erfblsr, expblsl * (segment_left + bin_left), expblsr * (segment_left + bin_right))
    x4 = x(brsl_sigma * erfbrsl, brsr_sigma * erfbrsr, expbrsl * (segment_right + bin_left), expbrsr * (segment_right + bin_right))

    C1 = 0.5 * (x2 - x1) / (bin_right - bin_left)
    C2 = 0.5 * (sigma ** 2 * (erfblsl - erfblsr - erfbrsl + erfbrsr) - x3 + x4) / (bin_right - bin_left)
    return C1, C2

def get_bin_edges(bins):
    """Convert reference wavelengths of detector bins to edge wavelengths
    
    `bins` is the array of reference wavelengths of length M. The M-1 inner bin edges are
    taken as midpoints between adjacent values of this array. The 2 outer bin edges are
    estimated by assuming that the first and the last bins are symmetric with respect to
    their reference wavelengths. The input must be strictly ascending
    
    Parameters
    ----------
    bins : array of length M
        Reference wavelengths of all detector bins
    
    Returns
    -------
    array of length M+1
        Bin edge wavelengths
    """
    if not np.all(bins[1:] > bins[:-1]):
        raise ValueError('Bin wavelengths must be strictly ascending')

    return np.concatenate([[bins[0] - (bins[1] - bins[0]) / 2], (bins[1:] + bins[:-1]) / 2, [bins[-1] + (bins[-1] - bins[-2]) / 2]])

def convolution_weights(bins, x, sigma, clip = 5.0, mode = 'window', max_size = 25e6, flush_cache = False):
    """Calculate the complete convolution matrix for an arbitrary flux density spectrum
    and an arbitrary set of detector bins
    
    The final observed spectrum is then given by `C * flux`, where `flux` is the (true) flux
    density spectrum sampled at wavelengths `x`, and `C` is the convolution matrix calculated
    by this function. Here `*` represents matrix multiplication (dot product)
    
    It is assumed that the flux density between the sampling wavelengths is given by
    linear interpolation, while beyond the range of `x` it is zero

    All of `bins`, `x` and `sigma` must be given in the same units, e.g. A

    Since the number of detector bins and flux wavelengths may be very large, the output is
    returned as a sparse matrix in the COOrdinate format. The function also checks that the
    memory usage does not exceed the limit provided in the optional argument `max_size`

    The function implements memory caching
    
    Parameters
    ----------
    bins : array of length M
        Reference wavelengths of all detector bins. See `get_bin_edges()`
    x : array of length N
        Array of wavelengths where the flux density spectrum is sampled. The segments are
        defined as the N-1 intervals between adjacent values of this array. Must be
        strictly ascending
    sigma : float or array of length M or array of length N
        Resolution of the detector, defined as the standard deviation of the Gaussian
        convolution kernel. If scalar, constant resolution is adopted for all bins and all
        segments. If `mode` is 'window', the resolution is defined at each detector bin, and
        the length of this array must be M. If `mode` is 'dispersion', the resolution is
        defined at each wavelength in the spectrum and the length of this array must be N
    clip : float, optional
        Assume that the weights for the detector bin / spectrum segment pair are zero if
        the wavelength ranges of the bin and the segment are separated by this many values of
        `sigma` for this pair. The argument removes negligibly small weights from the result,
        allowing the final sparse matrix to take up less memory
    mode : str, optional
        If 'window', the resolution determines the width of the window, over which each
        detector bin is sampling the spectrum flux density. If 'dispersion', the resolution
        determines the dispersion range, into which each segment of the spectrum
        spreads out before it is binned by the detector. There is no distinction between
        these two modes if the resolution is constant. This argument is therefore ignored
        when `sigma` is scalar
    max_size : int, optional
        Maximum number of non-zero elements in the convolution matrix. An exception is thrown
        if the predicted number of non-zero elements exceeds this argument. The number of
        non-zero elements directly correlates with memory usage
    flush_cache : bool, optional
        If True, discard cache and calculate the convolution matrix afresh
    
    Returns
    -------
    C : scipy.sparse._coo.coo_matrix
        Convolution weight matrix, such that the dot product of this matrix and the flux vector
        gives the vector of average flux densities in the detector bins
    """
    global _convolution_weights_cache
    try:
        _convolution_weights_cache
    except:
        _convolution_weights_cache = {}

    # Since input data can be very large, we will not attempt sorting it here. Leave it up to the user
    if not np.all(bins[1:] > bins[:-1]):
        raise ValueError('Bin wavelengths must be strictly ascending')
    if not np.all(x[1:] > x[:-1]):
        raise ValueError('x must be strictly ascending')

    # If `sigma` is scalar, we can choose any mode and populate the entire array with constant values
    try:
        sigma[0]
    except:
        sigma = np.full(len(bins), sigma)
        mode = 'window'

    # Check for cached output
    if not flush_cache:
        hash_string = ''.join(list(map(lambda arg: hashlib.sha256(bytes(arg)).hexdigest(), [bins, x, sigma]))) + str(clip) + mode
        if hash_string in _convolution_weights_cache:
            return _convolution_weights_cache[hash_string]

    # Check dimensions of `sigma` depending on the mode
    if mode == 'window':
        if len(bins) != len(sigma):
            raise ValueError('In "window" mode, must provide sigma for each bin')
    elif mode == 'dispersion':
        if len(x) != len(sigma):
            raise ValueError('In "dispersion" mode, must provide sigma for each x')
        # For dispersion, we want `sigma` at spectrum segments, not wavelength points
        sigma = (sigma[1:] + sigma[:-1]) / 2

    # Calculate bin edges
    bin_edges = get_bin_edges(bins)

    # Estimate the range of wavelengths that fall within the clipping range (clip * sigma) of each bin
    clip_start = np.zeros(len(bins), dtype = int)
    clip_end = np.zeros(len(bins), dtype = int)
    for i in range(len(bins)):
        if mode == 'dispersion':
            clip_indices = np.where(((x[:-1] - bin_edges[i + 1]) < (clip * sigma)) & ((bin_edges[i] - x[1:]) < (clip * sigma)))[0]
            if len(clip_indices) == 0:
                clip_start[i] = 0; clip_end[i] = 0
            else:
                clip_start[i] = np.min(clip_indices); clip_end[i] = np.max(clip_indices) + 1
        else:
            # If sigma is constant for the bin, we can find the clipping range more efficiently with np.searchsorted()
            clip_start[i] = np.searchsorted(x[1:], bin_edges[i] - clip * sigma[i])
            clip_end[i] = np.searchsorted(x[:-1], bin_edges[i + 1] + clip * sigma[i])

    # Get the row and column indices of each non-zero element of the (to be calculated) convolution matrix. Columns correspond to
    # detector bins, and rows correspond to spectrum segments
    row_lengths = clip_end - clip_start
    nonzero_count = np.sum(row_lengths)
    if nonzero_count > max_size:
        raise ValueError('Maximum memory usage exceeded. Requested: {}, allowed: {}. See the optional `max_size` argument'.format(nonzero_count, max_size))
    row_end = np.cumsum(row_lengths)
    row_start = row_end - row_lengths
    row_indices = np.zeros(nonzero_count, dtype = int)
    col_indices = np.zeros(nonzero_count, dtype = int)
    for i in range(len(bins)):
        row_indices[row_start[i] : row_end[i]] = np.arange(clip_start[i], clip_end[i])
        col_indices[row_start[i] : row_end[i]] = i

    # Obtain the values of every argument required by convolution_integral() for each non-zero element of the convolution matrix
    segment_left = x[:-1][row_indices]
    segment_right = x[1:][row_indices]
    bin_left = bin_edges[col_indices]
    bin_right = bin_edges[col_indices + 1]
    if mode == 'dispersion':
        sigma_all = sigma[row_indices]
    else:
        sigma_all = sigma[col_indices]

    # Run the convolution integral calculator
    C1, C2 = convolution_integral(sigma_all, segment_left, segment_right, bin_left, bin_right)

    # The way convolution_integral() works, we now have "observed = C1 * offset + C2 * slope". Instead, we want to convert these
    # two matrices into a single matrix C that gives "observed = C * flux". We can do this with SciPy's sparse matrix operations.
    # First, pad both C1 and C2 with empty rows and get their C1[:-1] and C1[1:] slices (low and high as they are called here)
    C1_low  = scp.sparse.coo_matrix((C1, (col_indices, row_indices + 1)), shape = [len(bins), len(x)])
    C2_low  = scp.sparse.coo_matrix((C2, (col_indices, row_indices + 1)), shape = [len(bins), len(x)])
    C1_high = scp.sparse.coo_matrix((C1, (col_indices, row_indices)), shape = [len(bins), len(x)])
    C2_high = scp.sparse.coo_matrix((C2, (col_indices, row_indices)), shape = [len(bins), len(x)])
    # Now do the same with the wavelength vector and do the conversion
    padded_x = np.insert(x, [0, len(x)], [0, 0])
    C = C1_high + (C2_low - C1_low.multiply(padded_x[:-2])).multiply(1 / (padded_x[1:-1] - padded_x[:-2])) - (C2_high - C1_high.multiply(padded_x[1:-1])).multiply(1 / (padded_x[2:] - padded_x[1:-1]))

    _convolution_weights_cache[hash_string] = C
    return C

def combine_arms(wl = None, flux = None, return_arm_index = False):
    """Combine wavelengths and fluxes recorded by individual spectrograph arms into
    a single spectrum

    The function returns the combined spectrum with wavelengths sorted in ascending
    order
    
    Parameters
    ----------
    wl : None or list or dict, optional
        If `wl` is a dictionary, it must be keyed by the identifiers of the arms that
        need to be combined. The values are 1D arrays of reference wavelengths of
        detector bins (see `get_bin_edges()`). Alternatively, provide a list of arm
        identifiers and the "typical" wavelength sampling for each arm (as defined in
        `settings['arms']`) will be assumed. Alternatively, set to `None` to use all
        arms defined in `settings['arms']` with "typical" wavelength sampling
    flux : None or dict, optional
        Dictionary of fluxes corresponding to the wavelength bins in `wl`. Must have
        the same keys and array lengths as `wl`. Alternatively, set to `None` if only
        wavelengths are required in the output
    return_arm_index : bool, optional
        If `True`, return the `arm_index` array that allows each value in the output to
        be traced to its arm of origin
    
    Returns
    -------
    wl : array_like
        Combined and sorted array of reference wavelengths for all arms
    flux : array_like
        Corresponding array of fluxes. Only returned if the optional argument `flux` is
        not `None`
    arm_index : array_like
        Corresponding array of spectrograph arm indices. Only returned if the optional
        argument `return_arm_index` is set to `True`. Each value in the array is an
        integer that determines the serial number of the arm of origin. Serial numbers
        are assigned starting with 0 in the alphabetical order of arm names
    """
    # Populate wl if not given
    if wl is None:
        wl = {arm: settings['arms'][arm]['wl'] for arm in settings['arms']}
    elif type(wl) is list:
        wl = {arm: settings['arms'][arm]['wl'] for arm in wl}
    else:
        wl = copy.deepcopy(wl) # We will be transforming wavelengths in place, so get a copy
    flux = copy.deepcopy(flux) # Same for flux

    # If flux is given, make sure its arms and dimensionality match wl
    if flux is not None:
        if set(wl.keys()) != set(flux.keys()):
            raise ValueError('The spectrograph arms in the `wl` and `flux` dictionaries do not match')
        if not np.all([len(wl[key]) == len(flux[key]) for key in wl]):
            raise ValueError('The dimensions of `wl` and `flux` do not match')

    # Combine the arms into a single spectrum
    keys = sorted(list(wl.keys()))
    arm_index = np.concatenate([np.full(len(wl[key]), i) for i, key in enumerate(keys)])
    wl = np.concatenate([wl[key] for key in keys])
    if flux is not None:
        flux = np.concatenate([flux[key] for key in keys])
    sort = np.argsort(wl)
    wl = wl[sort]
    arm_index = arm_index[sort]
    if flux is not None:
        flux = flux[sort]

    if return_arm_index:
        if flux is not None:
            return wl, flux, arm_index
        else:
            return wl, arm_index
    else:
        if flux is not None:
            return wl, flux
        else:
            return wl

def simulate_observation(wl, flux, detector_wl = None, mask_unmodelled = True, clip = 5, combine = True):
    """Simulate observation of the model spectrum by a spectrograph
    
    If the model does not fully cover the range of a spectrograph arm, the output flux density
    in the affected detector wavelength bins will be set to `np.nan`. This behavior can be
    disabled by setting `mask_unmodelled` to False, in which case the edge effects of the
    convolution will be left in the output spectrum
    
    Parameters
    ----------
    wl : array_like
        Model wavelengths
    flux : array_like
        Model flux densities corresponding to `wl`. It is assumed that the flux density between
        the values of `wl` can be obtained by linear interpolation, and the flux density beyond
        the range of `wl` is zero
    detector_wl : None or list or dict, optional
        Reference wavelengths of the detector bins in each arm of the spectrograph. If
        `detector_wl` is a dictionary, it must be keyed by the identifiers of the arms that
        are used in this observation. The values are 1D arrays of reference wavelengths of
        detector bins (see `get_bin_edges()`). Alternatively, provide a list of arm
        identifiers and the "typical" wavelength sampling for each arm (as defined in
        `settings['arms']`) will be assumed. Alternatively, set to `None` to use all
        arms defined in `settings['arms']` with "typical" wavelength sampling
    mask_unmodelled : bool, optional
        If True, set to `np.nan` the flux density in all bins that are affected by the finite
        wavelength range of the model spectrum. Otherwise, the bins near the edges of the
        model will suffer from convolution edge effects, and the bins far beyond the wavelength
        range of the model will receive zero flux
    clip : float, optional
        Sigma-clipping parameter for the convolution calculator (see `convolution_weights()`).
        This parameter also determines the range within which the lack of model coverage at a
        particular wavelength can affect the detector bins
    combine : bool, optional
        If True, return a single wavelength and a single flux array that represents the combined
        spectrum across all arms of the spectrograph (see `combine_arms()`). Otherwise, return
        the spectra in individual arms
    
    Returns
    -------
    wl : dict or array_like
        Wavelengths of the observed spectrum with the spectrograph (keyed by spectrograph arm
        if `combine` is True)
    flux : dict or array_like
        Corresponding flux densities (keyed by spectrograph arm if `combine` is True)
    """
    # Populate detector_wl if not given
    if detector_wl is None:
        detector_wl = {arm: settings['arms'][arm]['wl'] for arm in settings['arms']}
    elif type(detector_wl) is list:
        detector_wl = {arm: settings['arms'][arm]['wl'] for arm in detector_wl}

    # Resample the spectrum onto the detector bins of each arm
    detector_flux = {}
    for arm in detector_wl:
        if 'sigma' not in settings['arms'][arm]:
            sigma = settings['arms'][arm]['FWHM'] / (2 * np.sqrt(2 * np.log(2)))
        else:
            if 'FWHM' in settings['arms'][arm]:
                raise ValueError('Both FWHM and sigma provided for arm {}'.format(arm))
            sigma = settings['arms'][arm]['sigma']
        C = convolution_weights(detector_wl[arm], wl, sigma, clip = clip)
        detector_flux[arm] = C * flux
        # Remove wavelengths that exceed the modelled range
        if mask_unmodelled:
            message = 'In spectrograph arm {} the model does not cover the full wavelength range of the detector. Affected bins were set to np.nan'.format(arm)
            first = C.getcol(0).nonzero()[0]; last = C.getcol(-1).nonzero()[0]
            if len(first) == 0 and len(last) == 0:
                # If neither the first nor the last spectrum segment contribute to detected flux, the model either
                # covers the entire range of the detector or completely misses it
                if (np.min(wl) < np.min(detector_wl[arm])) and (np.max(wl) > np.max(detector_wl[arm])):
                    continue
                else:
                    detector_flux[arm] = np.full(len(detector_flux[arm]), np.nan)
                    warn(message)
                    continue
            warn(message)
            mask = np.full(len(detector_flux[arm]), True)
            if len(first) != 0:
                mask[detector_wl[arm] <= detector_wl[arm][np.max(first)]] = False
            if len(last) != 0:
                mask[detector_wl[arm] >= detector_wl[arm][np.min(last)]] = False
            detector_flux[arm][~mask] = np.nan


    # Combine the results into a single spectrum
    if combine:
        return combine_arms(detector_wl, detector_flux)
    else:
        return detector_wl, detector_flux

class ModelGridInterpolator:
    """Handler class for interpolating the model grid to arbitrary stellar parameters
    
    The class provides the `interpolate()` method to carry out the interpolation. The
    interpolation is linear with dynamic fetching of models from disk and caching for
    loaded models
    
    Attributes
    ----------
    statistics : dict
        Statistical data to track the interpolator's performance. Includes the following
        keys:
            'num_models_used': Total number of models read from disk *if* no caching was
                               used
            'num_models_loaded': Actual number of models read from disk
            'num_interpolations': Total number of interpolator calls
            'num_interpolators_built': Total number of interpolator objects constructed
                                       (scipy.interpolate.RegularGridInterpolator)
    """
    def __init__(self, resample = True, detector_wl = None, synphot_bands = [], mag_system = settings['default_mag_system'], reddening = settings['default_reddening'], max_models = 10000):
        """
        Parameters
        ----------
        resample : bool, optional
            If True, resample the models to the detector bins before interpolation. The
            resampling is carried out with `simulate_observation()`. Otherwise,
            interpolate without resampling over the original wavelength grid of the models
        detector_wl : None or list or dict, optional
            Reference wavelengths of the detector bins in each arm of the spectrograph. The
            argument is only relevant if `resample` is set to True. If `detector_wl` is
            a dictionary, it must be keyed by the identifiers of the arms that
            are used in this observation. The values are 1D arrays of reference wavelengths of
            detector bins (see `get_bin_edges()`). Alternatively, provide a list of arm
            identifiers and the "typical" wavelength sampling for each arm (as defined in
            `settings['arms']`) will be assumed. Alternatively, set to `None` to use all
            arms defined in `settings['arms']` with "typical" wavelength sampling
        synphot_bands : list, optional
            If synthetic photometry for the interpolated models is required, provide the desired
            colors as elements of this list. Each color must be a 2-element tuple with the
            filenames of the bluer and redder filters in the `settings['filter_dir']` directory
        mag_system : str, optional
            Magnitude system to use for synthetic photometry. Must be supported by `synphot()`
        reddening : float, optional
            Optical reddening parameter (E(B-V)) to use for synthetic photometry
        max_models : int, optional
            Maximum number of models to keep in the loaded models cache. If exceeded, the
            models loaded earliest will be removed from the cache. Higher numbers lead to
            less frequent disk access (and hence faster performance) but higher memory usage
        """
        # Populate detector_wl if not given
        if resample:
            if detector_wl is None:
                detector_wl = {arm: settings['arms'][arm]['wl'] for arm in settings['arms']}
            elif type(detector_wl) is list:
                detector_wl = {arm: settings['arms'][arm]['wl'] for arm in detector_wl}
            setattr(self, '_detector_wl', detector_wl)
        setattr(self, '_resample', resample)

        # Load grid dimensions
        setattr(self, '_grid', read_grid_dimensions())

        # Loaded models
        setattr(self, '_loaded', {})
        setattr(self, '_loaded_ordered', [])
        setattr(self, '_max_models', max_models)

        # Models embedded into the current interpolator
        setattr(self, '_interpolator_models', set())

        # Synthetic photometry parameters
        setattr(self, '_synphot_bands', synphot_bands)
        setattr(self, '_mag_system', mag_system)
        setattr(self, '_reddening', reddening)

        # Holders of statistical information
        setattr(self, 'statistics', {'num_models_used': 0, 'num_models_loaded': 0, 'num_interpolations': 0, 'num_interpolators_built': 0})

    def _build_interpolator(self, x):
        # Separate out real and virtual grid dimensions
        virtual_x = {}; real_x = {}
        for key in x:
            if key in settings['virtual_dof']:
                virtual_x[key] = x[key]
            else:
                real_x[key] = x[key]

        # Make sure real_x has the right dimensions
        if set(real_x.keys()) != set(list(self._grid.keys())):
            raise ValueError('Model grid dimensions and requested interpolation target dimensions do not match')

        # Make sure we are not exceeding the bounds of virtual dimensions
        for key in virtual_x:
            if (virtual_x[key] < settings['virtual_dof'][key][0]) or (virtual_x[key] > settings['virtual_dof'][key][1]):
                raise ValueError('Virtual dimension {} bounds exceeded'.format(key))

        # Which models are required to interpolate to real_x?
        subgrid = {}
        for key in self._grid.keys():
            if (real_x[key] > np.max(self._grid[key])) or (real_x[key] < np.min(self._grid[key])):
                raise ValueError('Model grid dimensions exceeded along {} axis'.format(key))
            if real_x[key] in self._grid[key]:
                subgrid[key] = np.array([real_x[key]])
            else:
                subgrid[key] = np.array([np.max(self._grid[key][self._grid[key] < real_x[key]]), np.min(self._grid[key][self._grid[key] > real_x[key]])])
        required_models = set(['|'.join((np.array(model) + 0.0).astype(str)) for model in itertools.product(*[subgrid[key] for key in sorted(list(subgrid.keys()))])])
        self.statistics['num_models_used'] += len(required_models)

        # Determine which of the required models have not been loaded yet
        new_models = [model for model in required_models if model not in self._loaded]

        # If the total number of models exceeds max_models, delete the ones loaded earlier
        if len(new_models) + len(self._loaded) > self._max_models:
            to_delete = len(new_models) + len(self._loaded) - self._max_models
            remaining = []
            for model in self._loaded_ordered:
                if model in required_models or to_delete == 0:
                    remaining += [model]
                else:
                    del self._loaded[model]
                    to_delete -= 1
            self._loaded_ordered = remaining

        # Helper function to run the model preprocessor and resample the model (if necessary)
        def preprocess(model, x, virtual_x, do_preprocess):
            if not do_preprocess:
                return model
            params = {key: x[i] for i, key in enumerate(sorted(list(subgrid.keys())))}
            params.update(virtual_x)
            wl = self._model_wl
            flux = preprocess_grid_model(wl, model[0] * 1.0, params, model[1])

            # Synthetic photometry
            colors = np.zeros(len(self._synphot_bands))
            if len(self._synphot_bands) != 0:
                if 'teff' not in params:
                    raise ValueError('The model grid must have "teff" as one of the axes to compute synthetic photometry')
                phot = synphot(wl, flux, params['teff'], np.unique(self._synphot_bands), mag_system = self._mag_system, reddening = self._reddening)
                for i, color in enumerate(self._synphot_bands):
                    if np.isnan(phot[color[0]]) or np.isnan(phot[color[1]]):
                        raise ValueError('Could not calculate synthetic color {} for model {}'.format(color, params))
                    # Note that the order of bands is reversed in the color calculation, since synphot() returns bolometric corrections, not magnitudes
                    colors[i] = phot[color[1]] - phot[color[0]]

            if self._resample:
                wl, flux = simulate_observation(wl, flux, self._detector_wl)
            try:
                self._wl
            except:
                setattr(self, '_wl', wl)

            if len(colors) > 0:
                flux = np.concatenate([colors, flux])

            return flux

        # Decide if we want to run preprocessing / resampling after the model is read or after it is accessed
        # The former choice makes sense if no virtual dimensions are being fit. Since virtual parameters can
        # change their values between accesses, we need to run preprocessing / resampling after each access
        resample_after_read = len(list(set(settings['fit_dof']) & set(settings['virtual_dof']))) == 0

        # Parallelization
        if ('parallel' in settings) and (settings['parallel'] == True):
            request = []
            for model in new_models:
                params = np.array(model.split('|')).astype(float)
                keys = sorted(list(real_x.keys()))
                params = {keys[i]: params[i] for i in range(len(keys))}
                params_ordered = [params[key] for key in sorted(list(subgrid.keys()))]
                request += [params]
            if len(request) > 1:
                read_grid_model(request, self._grid, batch = True)

        # Load the new models
        for model in new_models:
            params = np.array(model.split('|')).astype(float)
            keys = sorted(list(real_x.keys()))
            params = {keys[i]: params[i] for i in range(len(keys))}
            params_ordered = [params[key] for key in sorted(list(subgrid.keys()))]
            wl, flux, meta = read_grid_model(params, self._grid)
            if not np.all(np.diff(wl) > 0.0):
                raise ValueError('Model wavelengths must be strictly ascending')

            try:
                self._model_wl
            except:
                setattr(self, '_model_wl', wl)
            self._loaded[model] = (flux, meta)
            self._loaded[model] = preprocess(self._loaded[model], params_ordered, virtual_x, resample_after_read)
            self._loaded_ordered += [model]
        self.statistics['num_models_loaded'] += len(new_models)

        # Build the interpolator
        subgrid_ordered = [subgrid[key] for key in sorted(list(subgrid.keys()))]
        meshgrid = np.meshgrid(*subgrid_ordered, indexing = 'ij')
        spectra = np.vectorize(lambda *x: preprocess(self._loaded['|'.join((np.array(x) + 0.0).astype(str))], x, virtual_x, not resample_after_read), signature = ','.join(['()'] * len(self._grid)) + '->(n)')(*meshgrid)
        setattr(self, '_interpolator', scp.interpolate.RegularGridInterpolator(subgrid_ordered, spectra))
        self._interpolator_models = required_models
        self.statistics['num_interpolators_built'] += 1
        return self._interpolator

    def interpolate(self, params):
        """Carry out model grid interpolation
        
        Parameters
        ----------
        params : dict
            Dictionary of target stellar parameters. A value must be provided for each grid
            axis, keyed by the axis name
        
        Returns
        -------
        wl : array_like
            Grid of model wavelengths in A
        flux : array_like
            Interpolated flux densities in wavelength space in arbitrary units corresponding to the
            wavelengths in `wl`
        phot : array_like
            If requested (i.e. if `synphot_bands` is not empty), synthetic colors for the bands
            listed in `synphot_bands` with the desired reddening in the desired magnitude system
        """
        # Separate out real and virtual grid dimensions
        virtual_params = {}; real_params = {}
        for key in params:
            if key in settings['virtual_dof']:
                virtual_params[key] = params[key]
            else:
                real_params[key] = params[key]

        self.statistics['num_interpolations'] += 1
        interpolator = self._build_interpolator(params)
        result = interpolator([real_params[key] for key in sorted(list(real_params.keys()))])[0]
        if len(self._synphot_bands) != 0:
            return self._wl, result[len(self._synphot_bands):], result[:len(self._synphot_bands)]
        else:
            return self._wl, result

    def __call__(self, params):
        return self.interpolate(params)

def ranges_to_mask(arr, ranges, in_range_value = True, strict = False):
    """Convert a list of value ranges into a boolean mask, such that all values in `arr` that
    fall in any of the ranges correspond to `in_range_value`, and the rest correspond to
    `not in_range_value`
    
    Parameters
    ----------
    arr : array_like
        Array of values
    ranges : list of tuple
        List of two-element tuples, corresponding to the lower and upper bounds of each range
    in_range_value : bool, optional
        If True, in-range elements are set to True, and the rest are set to False in the mask.
        Otherwise, in-range values are set to False, and the rest are set to True
    strict : bool, optional
        If True, use strict comparison (`lower_bound < value < upper_bound`). Otherwise, use
        non-strict comparison (`lower_bound <= value <= upper_bound`)
    
    Returns
    -------
    array_like
        Boolean mask array of the same shape as `arr`
    """
    mask = np.full(np.shape(arr), not in_range_value)
    for window in ranges:
        if strict:
            mask[(arr > window[0]) & (arr < window[1])] = in_range_value
        else:
            mask[(arr >= window[0]) & (arr <= window[1])] = in_range_value
    return mask

def estimate_continuum(wl, flux, ivar, npix = 100, k = 3, masks = None, arm_index = None):
    """Estimate continuum correction in the spectrum using a spline fit
    
    The function carries out a weighted spline fit to a spectrum given by wavelengths in `wl`,
    flux densities in `flux` and using the weights in `ivar` (usually inverted variances)

    The wavelength regions given in `settings['masks']['continuum']` are excluded from the fit.
    If one of the excluded regions overlaps with the edge of the spectrum, the spline fit near
    the edge may be poorly conditioned (the spline will be extrapolated in that region, leading
    to potentially very large edge effects). If that part of the continuum is then used in
    stellar parameter determination, extremely poor convergence is likely. As such, when the
    spectral masks used by the fitter are passed in the optional `masks` parameter, it will be
    updated to mask out the affected region of the spectrum
    
    Parameters
    ----------
    wl : array_like
        Spectrum wavelengths
    flux : array_like
        Spectrum flux densities
    ivar : array_like
        Spectrum weights (inverted variances)
    npix : int, optional
        Desired interval between spline knots in pixels. The actual interval will be adjusted
        to keep the number of pixels in each spline segment identical
    k : int, optional
        Spline degree. Defaults to cubic
    masks : dict, optional
        Dictionary of boolean masks, keyed by stellar parameters. If given, this argument will be
        modified to exclude the regions of the spectrum potentially affected by spline extrapolation
        from the main fitter
    arm_index : array_like
        Indices of spectrograph arms corresponding to the provided wavelengths and fluxes. See
        `combine_arms()`
    
    Returns
    -------
    array_like
        Estimated continuum correction multiplier at each wavelength in `wl`
    """
    # Separate the spectrum into individual spectrograph arms if `uninterrupted_cont` is `False`, or merge
    # all arms together otherwise
    if arm_index is None or settings['uninterrupted_cont']:
        spline_index = np.zeros(len(wl))
    else:
        spline_index = arm_index

    # Build a mask of values to be included in the continuum estimation. If parameter masks are provided,
    # update them to avoid edge effects when necessary
    mask = (ivar > 0) & (~np.isnan(ivar)) & (~np.isnan(flux))
    for bad_continuum_range in settings['masks']['continuum']:
        for index in np.unique(spline_index):
            include = spline_index == index
            bad_continuum_range_mask = ranges_to_mask(wl[include], [bad_continuum_range], False)
            # Check for potential edge effects and remove the affected region from the fit
            if masks is not None:
                if (len(bad_continuum_range_mask[mask[include]]) > 0) and ((not bad_continuum_range_mask[mask[include]][-1]) or (not bad_continuum_range_mask[mask[include]][0])):
                    for param in masks:
                        masks[param][include] &= ranges_to_mask(wl[include], [bad_continuum_range], False)
            mask[include] &= bad_continuum_range_mask

    # Fit the spline
    result = np.full(len(wl), np.nan)
    for index in np.unique(spline_index):
        include = mask & (spline_index == index)
        t = wl[include][np.round(np.linspace(0, len(wl[include]), int(len(wl[include]) / npix))).astype(int)[1:-1]]
        if k == 0:
            result[spline_index == index] = np.full(np.count_nonzero(spline_index == index), np.sum(flux[include] * ivar[include]) / np.sum(ivar[include]))
        elif len(wl[include]) > k:
            spline = scp.interpolate.splrep(wl[include], flux[include], w = ivar[include], t = t, k = k)
            result[spline_index == index] = scp.interpolate.splev(wl[spline_index == index], spline)
        else:
            result[spline_index == index] = 1.0
    return result

def fit_model(wl, flux, ivar, initial, priors, dof, errors, masks, interpolator, arm_index, phot, method):
    """Fit the model to the spectrum
    
    Helper function to `chemfit()`. It sets up a model callback for `scp.optimize.curve_fit()` with
    the appropriate signature, defines the initial guesses and bounds for all free parameters, applies
    parameter masks and initiates the optimization routine

    The results of the fit and the associated errors are placed in the `initial` and `errors` arguments
    
    Parameters
    ----------
    wl : array_like
        Spectrum wavelengths
    flux : array_like
        Spectrum flux densities
    ivar : array_like
        Spectrum weights (inverted variances)
    initial : dict
        Initial guesses for the stellar parameters, keyed by parameter. The updated parameters after the
        optimization are stored in this dictionary as well
    priors : dict of 2-element tuples
        Prior estimates of the stellar parameters, keyed by parameter. Each element is a tuple with the
        first element storing the best estimate and the second element storing its uncertainty. All tuples
        of length other than 2 are ignored
    dof : list
        List of parameters to be optimized. The rest are treated as fixed to their initial values
    errors : dict
        Estimate errors in the best-fit parameter values are placed in this dictionary
    masks : dict
        Dictionary of boolean masks, keyed by stellar parameters. The masks determine which wavelengths are
        included in the fit for each parameter. If multiple parameters are fit simultaneously, the masks are
        logically added (or)
    interpolator : ModelGridInterpolator
        Model grid interpolator object that will be used to construct models during optimization
    arm_index : array_like
        Indices of spectrograph arms corresponding to the provided wavelengths and fluxes. See
        `combine_arms()`
    phot : dict
        Photometric colors of the star. Each color is keyed by `BAND1#BAND2`, where `BAND1` and `BAND2` are
        the transmission profile filenames of the filters, as required by `synphot()`. Each element is a
        2-element tuple, where the first element is the measured color, and the second element is the
        uncertainty in the measurement. The dictionary may also include optional elements `reddening`
        (E(B-V), single numerical value,), and `mag_system` (one of the magnitude systems supported by
        `synphot()`, single string)
    method : str
        Method to determine the best-fit stellar parameters. Must correspond to a callable function of form
        `fit_{method}(f, x, y, p0, sigma, bounds)`, where the arguments have the same meaning as those
        used by `scipy.optimize.curve_fit()`. The callable must return 3 values: best-fit parameter values
        in the same order as accepted by `f()`, errors in the best-fit parameter values, and a dictionary of
        additional data from the fit (e.g. covariance matrices, MCMC chains etc)
    Returns
    -------
    dict
        Additional data from the fit as returned by the fitting method callable function. If
        `settings['return_diagnostics']`, will also return detailed diagnostic data, including the observed
        spectrum, fitting masks and the best-fit model
    """
    # Create a dictionary to store diagnostic data in
    global _fitter_diagnostic_storage
    _fitter_diagnostic_storage = {}
    diagnostic = _fitter_diagnostic_storage

    # This function will be passed to curve_fit() as the model callback. The <signature> comment is a placeholder
    # to be replaced with the interpretation of the function signature later
    def f(x, params, mask, data_wl = wl, data_flux = flux, data_ivar = ivar, priors = priors, interpolator = interpolator, arm_index = arm_index, phot = phot, diagnostic = _fitter_diagnostic_storage):
        # <signature>

        # Load the requested model
        if len(interpolator._synphot_bands) != 0:
            model_wl, model_flux, model_phot = interpolator(params)
        else:
            model_wl, model_flux = interpolator(params)
        cont = estimate_continuum(data_wl, data_flux / model_flux, data_ivar * model_flux ** 2, npix = settings['cont_pix'], k = settings['spline_order'], arm_index = arm_index)
        diagnostic['model_wl'] = model_wl; diagnostic['model_flux'] = model_flux; diagnostic['model_cont'] = cont
        model_wl = model_wl[mask]; model_flux = (cont * model_flux)[mask]

        # Add priors
        index = 1
        for param in sorted(list(params.keys())):
            if len(np.atleast_1d(priors[param])) == 2:
                model_wl = np.concatenate([np.array([-index]), model_wl])
                model_flux = np.concatenate([np.array([params[param]]), model_flux])
                index += 1

        # Add photometric colors
        for i, color in enumerate(interpolator._synphot_bands):
            model_wl = np.concatenate([np.array([-index * 100]), model_wl])
            model_flux = np.concatenate([np.array([model_phot[i]]), model_flux])
            index += 1

        return model_flux

    # Define p0 and bounds
    p0 = [initial[param] for param in dof]
    def get_bounds(axis):
        if axis in interpolator._grid:
            return [np.min(interpolator._grid[axis]), np.max(interpolator._grid[axis])]
        elif axis in settings['virtual_dof']:
            return settings['virtual_dof'][axis]
        else:
            raise ValueError('Unknown axis {}'.format(axis))
    bounds = np.array([get_bounds(axis) for axis in dof]).T

    # Construct  and apply the fitting mask by superimposing the masks of individual parameters and removing bad pixels
    mask = np.full(len(wl), False)
    for param in dof:
        mask |= masks[param]
    mask &= (ivar > 0) & (~np.isnan(ivar)) & (~np.isnan(flux))
    mask &= ~np.isnan(interpolator(initial)[1])
    x = wl[mask]; y = flux[mask]; sigma = ivar[mask] ** -0.5

    # Since we do not a priori know the number of parameters being fit, we need to dynamically update the signature of the
    # model callback, f(). Unfortunately, there appears to be no better way to do that than by retrieving the source code
    # of the function with inspect, updating it, and reevaluating with exec()
    f = inspect.getsource(f).split('\n')
    f[0] = f[0].replace('params', ', '.join(dof))
    f[0] = f[0].replace('mask', 'mask = mask')
    f[1] = f[1].replace('# <signature>', 'params = {' + ', '.join(['\'{}\': {}'.format(param, [param, initial[param]][param not in dof]) for param in initial]) + '}')
    scope = {'priors': priors, 'phot': phot, 'interpolator': interpolator, 'arm_index': arm_index, 'mask': mask, 'np': np, 'wl': wl, 'flux': flux, 'ivar': ivar, 'estimate_continuum': estimate_continuum, 'settings': settings, 'synphot': synphot, '_fitter_diagnostic_storage': _fitter_diagnostic_storage}
    exec('\n'.join(f)[f[0].find('def'):], scope)
    f = scope['f']

    # As we add statistical and photometric priors, keep track of which one is which so we can provide their
    # labels to the user in diagnostic data
    prior_labels = []

    # Add priors. Each prior is just an extra pixel in the spectrum
    index = 1
    for param in sorted(list(priors.keys())):
        if len(np.atleast_1d(priors[param])) == 2:
            x = np.concatenate([np.array([-index]), x])
            y = np.concatenate([np.array([priors[param][0]]), y])
            sigma = np.concatenate([np.array([priors[param][1]]), sigma])
            prior_labels += [param]
            index += 1

    # Add photometric colors
    for color in interpolator._synphot_bands:
        color = '#'.join(color)
        x = np.concatenate([np.array([-index * 100]), x])
        y = np.concatenate([np.array([phot[color][0]]), y])
        sigma = np.concatenate([np.array([phot[color][1]]), sigma])
        prior_labels += [color]
        index += 1

    if method == 'cov':
        # If only the covariance is requested, we assume that "p0" has the best-fit parameters
        # and estimate the covariance matrix from the Jacobian using the standard formula
        #       COV = (J^T x diag(IVAR) x J)^-1 * CHI^2_REDUCED
        # where
        #       CHI^2_REDUCED = (SUM((OBSERVED - MODEL)^2 * IVAR) / (N_points - N_params))
        jacobian = scp.optimize._numdiff.approx_derivative(lambda args: f(x, *args), p0, bounds = bounds)
        residuals = (y - f(x, *p0)) / sigma
        chi2_red = np.sum(residuals ** 2) / (np.shape(jacobian)[0] - np.shape(jacobian)[1])
        weighted_jacobian = jacobian * (1 / sigma)[:, np.newaxis]
        cov = np.linalg.inv(weighted_jacobian.T @ weighted_jacobian) * chi2_red
        fit = [p0, np.diag(cov) ** 0.5, {'cov': cov}]

    else:
        # Run the optimizer
        fit = globals()['fit_{}'.format(method)](f, x, y, p0, sigma, bounds)

    # Save the results in "initial" and "errors"
    for i, param in enumerate(dof):
        initial[param] = fit[0][i]
        errors[param] = fit[1][i]

    # Provide diagnostic data if requested
    if settings['return_diagnostics']:
        fit[2]['observed'] = {'wl': wl, 'flux': flux, 'ivar': ivar}
        fit[2]['mask'] = mask
        fit[2]['fit'] = {'x': x, 'y': y, 'sigma': sigma, 'f': f(x, *fit[0]), 'p0': p0, 'bounds': bounds, 'dof': dof, 'priors': prior_labels[::-1]}
        fit[2]['model'] = {'wl': diagnostic['model_wl'], 'flux': diagnostic['model_flux'], 'cont': diagnostic['model_cont']}
        fit[2]['cost'] = (fit[2]['fit']['f'] - fit[2]['fit']['y']) ** 2.0 / fit[2]['fit']['sigma'] ** 2.0
        fit[2]['arm_index'] = arm_index

    return fit[2]

def fit_gradient_descent(f, x, y, p0, sigma, bounds):
    """Fit a 2D data series to a model using the Trust Region Reflective gradient descent algorithm
    
    The fit is carried out using `scipy.optimize.curve_fit()`. The covariance matrix is returned as
    additional data
    
    Parameters
    ----------
    See the parameters of `scipy.optimize.curve_fit()`
    
    Returns
    -------
    best : array_like
        Best-fit model parameters
    errors : array_like
        Errors in the best-fit model parameters
    extra : dict
        Dictionary with a single key, 'cov', that contains the covariance matrix of the fit
    """
    options = copy.deepcopy(settings['gradient_descent']['curve_fit'])
    # If the total number of data points is 1, it is impossible to normalize chi^2 so we need to
    # treat the provided errors as absolute regardless of the settings
    if len(x) == 1:
        options['absolute_sigma'] = True
    fit = scp.optimize.curve_fit(f, x, y, p0 = p0, sigma = sigma, bounds = bounds, **options)
    best = fit[0]
    errors = np.sqrt(np.diagonal(fit[1]))
    return best, errors, {'cov': fit[1]}

def fit_mcmc(f, x, y, p0, sigma, bounds):
    """Fit a 2D data series to a model using Markov Chain Monte Carlo (MCMC) sampling
    
    The fit is carried out using `emcee.EnsembleSampler()`. The MCMC chains for individual walkers are
    returned as additional data. The initial positions of the walkers are drawn from a random uniform
    distribution, within the prescribed bounds. The best-fit parameters and their errors are calculated
    as the median and the standard deviation of the chains after the removal of the burn-in steps
    
    Parameters
    ----------
    See the parameters of `scipy.optimize.curve_fit()`
    
    Returns
    -------
    best : array_like
        Best-fit model parameters
    errors : array_like
        Errors in the best-fit model parameters
    extra : dict
        Dictionary with a single key, 'chain', that contains the full MCMC chains for each walker
    """
    try:
        import emcee
    except:
        raise ImportError('emcee not installed')

    def mcmc_convergence(chain, c = 5):
        """Calculate the convergence parameters of an MCMC chain
        
        This function takes the chain output by emcee, and computes the autocorrelation length and the
        Geweke drift for each dimension of the parameter space
        
        Parameters
        ----------
        chain : array_like
            MCMC chain as returned by `emcee.EnsembleSampler().get_chain()`
        c : number, optional
            Step size for the autocorrelation window search (see `emcee.autocorr.integrated_time()`)
        
        Returns
        -------
        autocorr : array_like
            2D array with autocorrelation lengths for each parameter (first dimension) and each walker
            (second dimension), expressed as the number of autocorrelation lengths contained within the
            provided chain. Larger values indicate better convergence (emcee documentaion recommends
            requiring the minimum value of this array to exceed 50)
        geweke : array_like
            Array of Geweke drifts (z-scores) for each parameter. Larger values indicate poor convergence
        """
        nsteps, nwalkers, ndim = np.shape(chain)
        autocorr = np.zeros([ndim, nwalkers])
        geweke = np.zeros(ndim)
        for i in range(ndim):
            for w in range(nwalkers):
                f = emcee.autocorr.function_1d(chain[:, w, i])
                taus = 2.0 * np.cumsum(f) - 1.0
                windows = emcee.autocorr.auto_window(taus, c)
                tau = taus[windows]
                autocorr[i,w] = nsteps / tau

            flatchain = chain.reshape(chain.shape[0] * chain.shape[1], -1)
            a = flatchain[:len(flatchain) // 4, i]
            b = flatchain[-len(flatchain) // 4:, i]
            geweke[i] = (np.mean(a) - np.mean(b)) / (np.var(a) + np.var(b)) ** 0.5

        return autocorr, geweke

    # Choose initial walker positions
    if settings['mcmc']['initial'] == 'gradient_descent':
        best, errors, extra_gd = fit_gradient_descent(f, x, y, p0, sigma, bounds)
    initial = []
    for i in range(len(p0)):
        # Uniformly random initial walker positions
        if settings['mcmc']['initial'] == 'uniform':
            initial += [np.random.uniform(bounds[0][i], bounds[1][i], settings['mcmc']['nwalkers'])]
        # Gaussian initial positions based on gradient descent
        elif settings['mcmc']['initial'] == 'gradient_descent':
            initial += [scp.stats.truncnorm.rvs(loc = best[i], scale = errors[i], a = (bounds[0][i] - best[i]) / errors[i], b = (bounds[1][i] - best[i]) / errors[i], size = settings['mcmc']['nwalkers'])]
        else:
            raise ValueError('Unrecognized initial walker distribution {}'.format(settings['mcmc']['initial']))

    initial = np.vstack(initial).T

    def log_likelihood(p0, x, y, sigma, bounds):
        for i in range(len(p0)):
            if p0[i] <= bounds[0][i] or p0[i] >= bounds[1][i]:
                return -np.inf
        model = f(x, *p0)
        return np.sum(scp.stats.norm.logpdf(y, model, sigma))

    # Run the MCMC sampler
    sampler = emcee.EnsembleSampler(settings['mcmc']['nwalkers'], np.shape(initial)[1], log_likelihood, args = [x, y, sigma, bounds])
    sampler.run_mcmc(initial, settings['mcmc']['nsteps'], progress = settings['mcmc']['progress'])
    chain = sampler.get_chain(flat = False)
    autocorr, geweke = mcmc_convergence(chain)
    flatchain = chain[settings['mcmc']['discard']:,:,:].reshape((chain.shape[0] - settings['mcmc']['discard']) * chain.shape[1], -1)
    extra = {'chain': chain, 'log_prob': sampler.get_log_prob(), 'initial': initial, 'autocorr': autocorr, 'geweke': geweke}
    if settings['mcmc']['initial'] == 'gradient_descent':
        extra['gradient_descent'] = extra_gd
        extra['gradient_descent']['fit'] = best
        extra['gradient_descent']['errors'] = errors
    return np.median(flatchain, axis = 0), np.std(flatchain, axis = 0), extra

def chemfit(wl, flux, ivar, initial, phot = {}, method = 'gradient_descent', dof = False):
    """Determine the stellar parameters of a star given its spectrum
    
    Parameters
    ----------
    wl : dict
        Spectrum wavelengths keyed by spectrograph arm
    flux : dict
        Spectrum flux densities keyed by spectrograph arm
    ivar : dict
        Spectrum weights (inverted variances) keyed by spectrograph arm
    initial : dict
        Initial guesses for the stellar parameters, keyed by parameter. Each parameter supported
        by the model grid must be listed, except those for which default initial guesses are defined
        in `settings['default_initial']`. The value of each element is either a float or a
        2-element tuple. In the former case, the value is treated as the initial guess to the
        fitter. Otherwise, the first element is treated as an initial guess and the second value
        is treated as the prior uncertainty in the parameter
    phot : dict, optional
        Photometric colors of the star (if available). The colors and the spectrum will be fit to
        the models simultaneously to attain stricter constraints on the stellar parameters. Each
        color is keyed by `BAND1#BAND2`, where `BAND1` and `BAND2` are the transmission profile
        filenames of the filters, as required by `synphot()`. Each element is a 2-element tuple,
        where the first element is the measured color, and the second element is the uncertainty
        in the measurement. The dictionary may also include optional elements `reddening` (E(B-V),
        single numerical value), and `mag_system` (one of the magnitude systems supported by
        `synphot()`, single string)
    method : str
        Method to use for model fitting. Currently supported methods are 'gradient_descent' that
        employs the Trust Region Reflective algorithm implemented in `scipy`, and 'mcmc' that uses
        the MCMC sampler implemented in `emcee`. This parameter can also be set to 'cov' in which
        case no fitting will be carried out. Instead, it will be assumed that `initial` is already
        the best-fit values, and the fitter will compute the covariance matrix from the local
        Jacobian using the standard error propagation formula
    dof : list
        List of degrees of freedom to fit for. Defaults to `False`, in which case the list is adopted
        from `settings['fit_dof']`
    
    Returns
    -------
    dict
        Fitting results. Dictionary with the following keys:
            'fit': Best-fit stellar parameters
            'errors': Standard errors in the best-fit parameters from the diagonal of covariance matrix
            'interpolator_statistics': Interpolator statistics (see `ModelGridInterpolator().statistics`)
            'warnings': Warnings issued during the fitting process
            'extra': Additional fit data and diagnostic data depending on the chosen fitting method
    """
    # Remember the length of pre-existing warnings stack
    warnings_stack_length = len(warnings_stack)

    # Set default initial guesses
    initial = copy.deepcopy(initial)
    if 'default_initial' in settings:
        for param in settings['default_initial']:
            if (param not in initial):
                initial[param] = settings['default_initial'][param]

    # Combine arms
    wl_combined, flux_combined, arm_index = combine_arms(wl, flux, return_arm_index = True)
    ivar_combined = combine_arms(wl, ivar)[1]

    # Build the interpolator and resampler
    synphot_bands = [color.split('#') for color in phot if len(color.split('#')) == 2]
    if 'reddening' in phot:
        reddening = phot['reddening']
    else:
        reddening = settings['default_reddening']
    if 'mag_system' in phot:
        mag_system = phot['mag_system']
    else:
        mag_system = settings['default_mag_system']
    interpolator = ModelGridInterpolator(detector_wl = wl, synphot_bands = synphot_bands, reddening = reddening, mag_system = mag_system, max_models = settings['max_model_cache'])

    # Get the fitting masks for each parameter
    masks = {}
    for param in initial:
        mask = {}
        ranges_specific = []
        ranges_general = []
        for arm in wl:
            if (arm in settings['masks']) and (param in settings['masks'][arm]):
                ranges_specific += list(settings['masks'][arm][param])
            if ('all' in settings['masks']) and (param in settings['masks']['all']):
                ranges_specific += list(settings['masks']['all'][param])
            if (arm in settings['masks']) and ('all' in settings['masks'][arm]):
                ranges_general += list(settings['masks'][arm]['all'])
            if ('all' in settings['masks']) and ('all' in settings['masks']['all']):
                ranges_general += list(settings['masks']['all']['all'])
        masks[param] = ranges_to_mask(wl_combined, ranges_general)
        if len(ranges_specific) != 0:
            masks[param] &= ranges_to_mask(wl_combined, ranges_specific)
    # Run the continuum fitter once to give it a chance to update the fitting masks if it needs to
    cont = estimate_continuum(wl_combined, flux_combined, ivar_combined, npix = settings['cont_pix'], k = settings['spline_order'], masks = masks, arm_index = arm_index)

    # Preliminary setup
    fit = {param: np.atleast_1d(initial[param])[0] for param in initial}       # Initial guesses for the fitter
    errors = {}                                                                # Placeholder for fitting errors

    # Which degrees of freedom to fit for?
    if type(dof) is bool:
        dof = settings['fit_dof']

    # Run the main fitter
    extra = fit_model(wl_combined, flux_combined, ivar_combined, fit, initial, np.atleast_1d(dof), errors, masks, interpolator, arm_index, phot, method = method)

    # Get the texts of unique issued warnings
    warnings = np.unique(warnings_stack[warnings_stack_length:])
    inv_warnings_messages = {warnings_messages[key]: key for key in warnings_messages}
    warnings = [inv_warnings_messages[warning_id] for warning_id in warnings]
    interpolator_statistics = copy.deepcopy(interpolator.statistics)
    del interpolator

    return {'fit': fit, 'errors': errors, 'extra': extra, 'interpolator_statistics': interpolator_statistics, 'warnings': warnings}

def synphot(wl, flux, teff, bands, mag_system = settings['default_mag_system'], reddening = settings['default_reddening']):
    """Calculate synthetic photometry for a given spectrum

    The function supports reddening corrections. Reddening is parameterized by a single degree of freedom
    that is the optical reddening, E(B-V). The function applies extinction to each pixel of the spectrum,
    using the FM07 reddening law implemented in the `extinction` Python module. Rv=3.1 is assumed

    The output of the function is bolometric corrections. Photometric colors can then be calculated as the
    inverted differences of bolometric corrections. E.g. to compute B-V, the bolometric correction of B is
    subtracted from the bolometric correction of V. Absolute magnitudes can be estimated if the bolometric
    luminosity of the star is known, e.g. from an evolutionary model

    Parameters
    ----------
    wl : array_like
        Wavelength array of the spectrum in A
    flux : array_like
        Corresponding spectral intensity array in erg / s / cm^2 / A / sr
    teff : number
        Effective temperature of the star in K. Note that the effective temperature does not affect
        photometric colors, i.e. it offsets all bolometric corrections by the same amount
    bands : list
        Filenames that store the transmission profiles of the photometric bands to calculate synthetic
        photometry in. Only the basenames must be provided (i.e. no paths). The files must be stored in
        `settings['filter_dir']`. Each file should be a whitespace-separated ASCII table with two columns:
        wavelength in A, and transmission (between 0 and 1)
    mag_system : str, optional
        Magnitude system for synthetic photometry. This must be either `ABMAG` or one of the keys of the
        `settings['mag_systems']` dictionary. If it is the latter, the value of the dictionary must be the
        full file path to the reference spectrum of the photometric system (likely Vega). The spectrum must
        be a comma-separated CSV with two columns: wavelength in A and flux density in erg / s / cm^2 / A.
        The argument defaults to the value in `settings['default_mag_system']`
    reddening : number, optional
        Optical reddening to incorporate in the synthetic photometry, E(B-V) in mag. The argument defaults to
        the value in `settings['default_reddening']`

    Returns
    -------
    dict
        Bolometric corrections (in mag) in each of the requested photometric bands
    """
    global _synphot_cache
    try:
        _synphot_cache
    except:
        _synphot_cache = {'bands': {}, 'refs': {}}

    if not np.all(np.diff(wl) > 0.0):
        raise ValueError('Wavelengths in spectrum received for synthetic photometry are out of order')

    # Equation 5 from Gerasimov+2022. Note that 71.197425 is the IAU constant approximately equaling M_bol_Sun + 2.5 * np.log10(Lsun in Watts)
    C = 71.197425 - 2.5 * np.log10(4 * np.pi * scp.constants.sigma * (10 * scp.constants.parsec) ** 2.0 * 1000 ** 4.0)
    c = scp.constants.c * 1e+10

    result = {}
    for band in np.atleast_1d(bands):
        # Load the transmission profile of the photometric band
        if band not in _synphot_cache['bands']:
            try:
                _synphot_cache['bands'][band] = np.loadtxt(settings['filter_dir'] + '/' + band, unpack = True)
            except:
                raise ValueError('Unable to load transmission of photometric band {}'.format(band))
            if not np.all(np.diff(_synphot_cache['bands'][band][0]) > 0.0):
                raise ValueError('Wavelengths in transmission of photometric band {} are out of order'.format(band))
        filter_wl, filter_t = _synphot_cache['bands'][band]

        # Load the reference spectrum of the chosen magnitude system
        if mag_system == 'ABMAG':
            # Equation 7 from Gerasimov+2022
            f_ref_wl = wl
            f_ref_flux = (3631 * 1e-23 * c) / (wl ** 2)
        elif mag_system in settings['mag_systems']:
            if mag_system not in _synphot_cache['refs']:
                try:
                    _synphot_cache['refs'][mag_system] = np.loadtxt(settings['mag_systems'][mag_system], unpack = True, delimiter = ',')
                except:
                    raise ValueError('Unable to load reference spectrum for {}'.format(mag_system))
                if not np.all(np.diff(_synphot_cache['refs'][mag_system][0]) > 0.0):
                    raise ValueError('Wavelengths in reference spectrum {} are out of order'.format(mag_system))
            f_ref_wl, f_ref_flux = _synphot_cache['refs'][mag_system]
            if filter_wl.max() > f_ref_wl.max() or filter_wl.min() < f_ref_wl.min():
                result[band] = np.nan
                continue
        else:
            raise ValueError('Unknown magnitude system {}'.format(mag_system))
        filter_t_ref = np.interp(f_ref_wl, filter_wl, filter_t, left = 0, right = 0)

        if filter_wl.max() > wl.max() or filter_wl.min() < wl.min():
            result[band] = np.nan
            continue
        filter_t = np.interp(wl, filter_wl, filter_t, left = 0, right = 0)

        if reddening != 0.0:
            try:
                import extinction
            except:
                raise ValueError('Could not import the extinction module. For extinguished photometry, make sure the module is installed ( https://github.com/kbarbary/extinction ).')
            Rv = 3.1 # Rv must be 3.1 for the FM07 law
            A = extinction.fm07(wl, reddening * Rv)
        else:
            A = 0.0

        # Equation 6 from Gerasimov+2022
        f1 = wl * flux * filter_t * 10 ** (-A / 2.5) * np.pi
        f2 = f_ref_wl * f_ref_flux * filter_t_ref
        # Equation 4 from Gerasimov+2022
        result[band] = 2.5 * np.log10(np.trapz(f1, wl) / np.trapz(f2, f_ref_wl)) - 10 * np.log10(teff / 1000) + C

    return result
