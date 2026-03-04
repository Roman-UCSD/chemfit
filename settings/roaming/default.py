############################################################
#                                                          #
#             DEFAULT CHEMFIT SETTINGS PRESET              #
#                                                          #
#   This preset sets the default synthetic photometry      #
#   parameters, creates a test spectrograph arm, defines   #
#   standard convergence and output parameters             #
#                                                          #
############################################################

settings = {
    'griddir': None,

    ### Synthetic photometry ###
    'filter_dir': script_dir + '/bands/',                                    # Path to the transmission profile directory
    'mag_systems': {'VEGAMAG': script_dir + '/misc/vega_bohlin_2004.dat'},   # Reference spectra for magnitude systems (ABMAG is added automatically)
    'default_mag_system': 'VEGAMAG',                                         # Default magnitude system
    'default_reddening': 0.0,                                                # Default E(B-V)

    ### Spectrograph settings ###
    'arms': {                                      # Parameters of individual spectrograph arms
        'default_arm': {
            'FWHM': 2.07,                              # FWHM of the line spread function in the arm in A
            'wl': np.linspace(3800, 6500, 4096),       # Default bin wavelengths in the arm in A
        },
    },

    ### Blank fitting mask ###
    'mask': [],

    'max_model_cache': 1000, # Maximum number of models allowed in memory

    ### Optimization parameters ###
    'gradient_descent': {
        'curve_fit': {
            'ftol': 1e-10,
            'gtol': 1e-10,
            'xtol': 1e-10,
        },
    },
    'mcmc': {
        'nwalkers': 32,
        'nsteps': 5000,
        'initial': 'gradient_descent',
        'progress': True,
    },
    'cont_pix': 165,
    'spline_order': 3,
    'uninterrupted_cont': False,
}

def read_grid_model(params, grid):
    """Load a specific model spectrum from the model grid

    All models within the same model grid are expected to have the same wavelength
    sampling

    This function definition is a template. The actual implementation is deferred to another
    settings preset file that will configure chemfit to use a specific model grid

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

    This function definition is a template. The actual implementation is deferred to another
    settings preset file that will configure chemfit to use a specific model grid

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

def delete_grid_model(meta):
    """Safely delete loaded model data

    This function is called whenever a model loaded by the ModelGridInterpolator() object
    is being deleted either because the maximum allowed number of loaded models
    (`settings['max_model_cache']`) is exceeded, or because the entire interpolator
    object is being garbage collected. The function can be used to safely delete any data
    that may be missed by the garbage collector, such as disk cache

    This function definition is a template. The actual implementation is deferred to another
    settings preset file that will configure chemfit to use a specific model grid

    Parameters
    ----------
    meta : dict
        Dictionary with additional model data, as loaded by `read_grid_model()`
    """
    raise NotImplementedError()

def read_grid_dimensions(flush_cache = False):
    """Determine the available dimensions in the model grid and the grid points
    available in those dimensions

    This function definition is a template. The actual implementation is deferred to another
    settings preset file that will configure chemfit to use a specific model grid

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
