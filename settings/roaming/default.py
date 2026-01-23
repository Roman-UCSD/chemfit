############################################################
#                                                          #
#             DEFAULT CHEMFIT SETTINGS PRESET              #
#                                                          #
#   This preset sets the default synthetic photometry      #
#   parameters, creates a test spectrograph arm, provides  #
#   basic fitting masks and defines standard convergence   #
#   and output parameters                                  #
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

    ### Fitting masks ###
    'masks': {
        'all': {
            'all': [[100, 100000]],
        },
        'continuum': [],
    },

    'max_model_cache': 1000, # Maximum number of models allowed in memory

    ### Optimization parameters ###
    'return_diagnostics': True, # Return best-fit model, continuum correction and fitting masks in the chemfit.chemfit() output
    'gradient_descent': {
        'curve_fit': {
            'absolute_sigma': False,
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

    ### Warnings ###
    'throw_python_warnings': True,
}
