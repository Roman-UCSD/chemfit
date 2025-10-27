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

    ### Which parameters to fit for? All dimensions of the grid + redshift ###
    'fit_dof': ['zscale', 'alpha', 'teff', 'logg', 'carbon', 'redshift'],

    ### Bounds for virtual parameters, i.e. parameters that are not dimensions of the grid ###
    'virtual_dof': {'redshift': [-300, 300]},

    ### Default initial guesses ###
    'default_initial': {'redshift': 0.0},

    ### Update the fitting mask ###
    'masks': apply_standard_mask(defective_mask, original_settings['masks']),
}

_global = {}

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
