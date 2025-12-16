############################################################
#                                                          #
#                     LOCALFIT PRESET                      #
#                                                          #
#  This preset is designed to measure the abundances of    #
#  individual chemical elements using response functions   #
#  calculated at runtime. In order to use this preset, the #
#  stellar parameters (Teff, log(g), [M/H], [a/M], [C/M])  #
#  must be provided as input (e.g. from a prior grid       #
#  search with the GRIDFIT preset); however, the fitter    #
#  can be configured to allow some of the stellar          #
#  parameters to deviate from their input values for       #
#  better accuracy                                         #
#                                                          #
#  This preset requires BasicATLAS                         #
#                                                          #
############################################################


try:
    import atlas
except:
    raise ValueError('Could not import BasicATLAS. Make sure BasicATLAS is installed and available in your Python environment')

if 'scratch' not in original_settings:
    raise ValueError('You must create a local settings preset (settings/local/localfit.py) that defines the "scratch" key to point to the intended scratch directory (temporary storage for new models)')


import pickle
import scipy as scp
import hashlib
import os
import subprocess
import shutil
import copy
import sys

# Import PyTLAS. PyTLAS should check itself if it is properly compiled
sys.path.append('{}/PyTLAS/'.format(atlas.python_path))
import PyTLAS


# This is the "defective" mask that removes parts of the spectrum where the models do not match the spectra of 
# Arcturus and the Sun well
defective_mask = [[3800, 4000], [4006, 4012], [4065, 4075], [4093, 4110], [4140, 4165], [4170, 4180], [4205, 4220],
                 [4285, 4300], [4335, 4345], [4375, 4387], [4700, 4715], [4775, 4790], [4855, 4865], [5055, 5065],
                 [5145, 5160], [5203, 5213], [5885, 5900], [6355, 6365], [6555, 6570], [7175, 7195], [7890, 7900],
                 [8320, 8330], [8490, 8505], [8530, 8555], [8650, 8672]]

settings = {
    # Stellar parameters (teff, logg, zscale, alpha, carbon) of the observed spectrum. A subset of these parameters
    # can be allowed to deviate from the provided values using the 'gridfit_offsets' setting. All 5 parameters must be
    # specified for each fit, which is done automatically by `public__localfit()`
    'gridfit_params': False,

    # If True, a new set of ODFs and model atmospheres will be calculated for the fit. Otherwise, the precomputed grid
    # of atmosphere structures will be interpolated to the required stellar parameters. This setting is configured
    # automatically by `public__localfit()` depending on the chosen level of analysis
    'compute_new_structure': False,

    # Lowest convergence class for newly calculated model atmospheres which is considered acceptable. If a newly calculated
    # model atmosphere does not meet this criterion (has a worse convergence class), an exception is raised. This setting
    # is only relevant if 'compute_new_structure' is `True`
    'min_convergence': 'SILVER',

    # Offsets in the abundances of individual elements (from the chemistry of stellar parameters), at which the null
    # spectra and the model atmospheres (if 'compute_new_structure' and 'use_initial_abundance_offsets_in_structure' are
    # `True`) are calculated. This setting is configured automatically by `public__localfit()`
    'initial_abundance_offsets': {},

    # If True, the abundance offsets in 'initial_abundance_offsets' are incorporated in the model atmospheres to be used
    # by the fitter. This flag cannot be applied to interpolated structures, so if 'compute_new_structure' is set to `False`,
    # and this flag is set to `True`, an exception is raised. This setting is configured automatically by `public__localfit()`
    'use_initial_abundance_offsets_in_structure': True,

    # For stellar parameters that should be allowed to deviate from their 'gridfit_params' values, this setting specifies
    # the offsets from the 'gridfit_param' value that will define the interpolation grid used by the fitter
    'gridfit_offsets': {'logg': [-1.0, -0.5, 0.0, 0.5, 1.0]},

    # Which individual elements to fit, and the lines of which species (atomic, molecular and ionic) should be considered
    # when calculating the response functions for each element. All species are specified as Kurucz codes
    'elements': {
       'Ba': [56, 56.01, 56.02, 56.03],
       'Eu': [63, 63.01, 63.02, 63.03],
       'Mg': [12, 12.01, 12.02, 12.03, 12.04, 12.05, 112.0, 812.0, 112.01, 10812.0],
       'Si': [14, 14.01, 14.02, 14.03, 14.04, 14.05, 114.0, 814.0, 114.01, 80814.0, 101010114.0, 614.0, 60614.0],
       'Ca': [20, 20.01, 20.02, 20.03, 20.04, 20.05, 20.06, 20.07, 20.08, 20.09, 120.0, 820.0, 120.01, 10820.0],
       'Ti': [22, 22.01, 22.02, 22.03, 22.04, 22.05, 22.06, 22.07, 22.08, 22.09, 122.0], #, 822.0],
       'Ni': [28, 28.01, 28.02, 28.03, 28.04, 28.05, 28.06, 28.07, 28.08, 28.09, 128.0, 828.0],
       'Na': [11, 11.01, 11.02, 11.03, 11.04, 11.05, 111.0, 10811.0],
       'Fe': [26, 26.01, 26.02, 26.03, 26.04, 26.05, 26.06, 26.07, 26.08, 26.09, 126.0, 826.0],
       'Mn': [25, 25.01, 25.02, 25.03, 25.04, 25.05, 25.06, 25.07, 25.08, 25.09, 125.0, 825.0],
       'Co': [27, 27.01, 27.02, 27.03, 27.04, 27.05, 27.06, 27.07, 27.08, 27.09, 127.0, 827.0],
       'Li': [3, 3.01, 3.02, 3.03, 3.04, 3.05, 103.0],
       'Al': [13, 13.01, 13.02, 13.03, 13.04, 13.05, 113.0, 813.0, 113.01],
       'K': [19, 19.01, 19.02, 19.03, 19.04, 19.05, 119.0],
       'Cr': [24, 24.01, 24.02, 24.03, 24.04, 24.05, 24.06, 24.07, 24.08, 24.09, 124.0, 824.0],
       'La': [57, 57.01, 57.02, 57.03, 857.0],
       'Sc': [21, 21.01, 21.02, 21.03, 21.04, 21.05, 21.06, 21.07, 21.08, 21.09, 121.0, 821.0],
       'V': [23, 23.01, 23.02, 23.03, 23.04, 23.05, 23.06, 23.07, 23.08, 23.09, 123.0, 823.0],
       'Y': [39, 39.01, 39.02, 39.03, 39.04, 39.05, 839.0],
    },

    # Spectral synthesis parameters
    'wl_start': 350,           # Start wavelength in nm
    'wl_end': 1300,            # End wavelength in nm
    'res': 300000,             # Native resolution of the spectrum
    'air_wl': True,            # Use air wavelengths (vacuum if `False`)
    'atoms': 'BasicATLAS',     # Atomic linelist to use (BasicATLAS is recommended)

    # Abundance offsets at which response functions are calculated (the minimum and maximum values will determine the range of
    # offsets considered by the fitter, and for abundances between specified values linear interpolation will be used)
    'abun': np.round(np.arange(-2.0, 2.01, 0.1), 1),

    # Full list of alpha elements. The list must match the definition of alpha in the stellar parameters
    'alpha_elements': ['O', 'Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti'],

    # Print notifications?
    'silent': False,

    # Allow parallel calculation of ODFs (only relevant if 'compute_new_structure' is `True`)
    'ODF_multithreading': False,

    # We want to cache all calculated models so that the null-spectra and opacities can be reused once calculated without storing
    # them on disk
    'max_model_cache': 99999,

    # Scratch directory to store atmospheric structures and linelists. Must be specified in the local settings preset
    'scratch': original_settings['scratch'],

    # Apply the defective mask to the fits
    'masks': apply_standard_mask(defective_mask, original_settings['masks']),
}

# Define the ranges of "virtual" parameters (i.e. parameters that are not dimensions of the model grid). For this preset, the virtual parameters
# are element abundances and redshift
settings['virtual_dof'] = {'redshift': [-300, 300], **{element: [np.min(settings['abun']), np.max(settings['abun'])] for element in settings['elements']}}

# Which degrees of freedom (parameters) to consider in the fit? Here, we include redshift, element abundances and whichever stellar parameters
# are allowed to vary by the 'gridfit_offsets' setting
settings['fit_dof'] = [param for param in settings['gridfit_offsets'] if len(settings['gridfit_offsets'][param]) > 1] + ['redshift'] + [element for element in settings['elements']]

# Default initial guesses for the parameters
settings['default_initial'] = {'redshift': 0.0, **{element: 0.0 for element in settings['elements']}}

# Relationship between C12/C13 isotope ratio and surface gravity from Kirby+2015
C12C13_kirby_2015 = lambda logg: np.where(logg > 2.7, 50, np.where(logg <= 2.0, 6, 63 * logg - 120))

# Relationship between VTURB and LOGG based on DEIMOS spectra from Gerasimov+2026
VTURB_LOGG = lambda logg: 2.792 * np.exp(-0.241 * logg -0.118)

# Global variable to store the linelists
linelists = {}

# Global variable to store the header of the ATLAS restart database
ATLAS_restarts_header = atlas.restarts.load_header()
if os.path.basename(ATLAS_restarts_header['grid']) != 'light.h5':
    raise ValueError('BasicATLAS must be configured to use the full grid of restart models')

# PyTLAS library instances
xnfpelsyn = PyTLAS.init_xnfpelsyn()
synthe = PyTLAS.init_synthe()
spectrv = PyTLAS.init_spectrv()

# References to chemfit.chemfit(), chemfit.simulate_observation() and chemfit.estimate_continuum() which will be populated when this
# preset is initialized
def main__chemfit():
    pass
def main__simulate_observation():
    pass
def main__estimate_continuum():
    pass

def notify(message, color = 'k'):
    """Print a status notification
    
    This is a wrapper around `print()`, which makes it suppressible by `settings['silent']`. We also allow color printing, with the
    default color ('k') meaning "uncolored". Printed messages are flushed immediately to provide real time updates when running
    the fitter on a cluster
    
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

def generate_hash(data):
    """Generate an MD5 hash for an arbitrary data structure

    The function serializes the data structure using the `pickle` module before calculating the MD5 checksum. This function is primarily
    intended to generate unique filenames for model atmospheres, ODFs and linelists based on their parameters
    
    Parameters
    ----------
    data : object
        Any data structure that can be serialized with `pickle`
    
    Returns
    -------
    str
        MD5 checksum of the provided data structure
    """
    return hashlib.md5(pickle.dumps(data, protocol = pickle.HIGHEST_PROTOCOL)).hexdigest()

def code_to_nelion(code):
    """Convert a Kurucz species code into the corresponding NELION index, i.e. the index of that species in the XNFPEL and DOPPLE arrays
    defined by XNFPELSYN
    
    Kurucz species codes are used by ATLAS to uniquely identify atomic, ionic, and molecular species:
        - Atomic species are given by their atomic (proton) number. E.g., neutral carbon is 6.00
        - Ionic species encode the ionization state in the hundredths place. E.g., doubly ionized carbon (C III) is 6.02
        - Molecular species encode the constituent atoms as odd powers of ten, ordered by decreasing atomic mass. E.g., neutral SiC2 is
          60614.00, corresponding to Si (14 x 1 = 14), the first C atom (6 x 100 = 600) and the second C atom (6 x 10000 = 60000)
        - Molecular ions are represented by adding the ionization state in the hundredths place, in direct analogy with atomic ions.
          E.g., MgH+ is 112.01 and OH+ is 108.01
        - Negative ions are in principle represented as molecules whose constituents include electrons with the proton number of zero.
          E.g., the H- ion is represented as 100.00. However, the lines of negative ions are not included in SYNTHE, and there are no
          NELION codes that correspond to them

    The SYNTHE suite stores species-dependent quantities (e.g. number densities) in flattened arrays of shape 139x6 (Fortran ordering).
    Each species occupies a fixed position in these arrays, identified by a serial index referred to as NELION. This function maps
    a given Kurucz species code to its corresponding NELION index. An exception is raised if the input Kurucz code does not correspond
    to any species allocated in these arrays

    Parameters
    ----------
    code : number
        Input Kurucz code of the atomic, ionic or molecular species
    
    Returns
    -------
    int
        NELION index of the input species
    """
    # For supported molecular species, we define a direct mapping of Kurucz codes to NELION indices
    code_mol = [101., 106., 107., 108., 606., 607., 608., 707., 708., 808., 112., 113., 114., 812., 813., 814., 116., 120., 816., 820., 821., 822., 823., 103., 104., 105., 109., 115., 117., 121., 122., 123., 124., 125., 126.,106.01, 107.01,108.01,112.01,113.01,114.01,120.01, 111., 119.,10101.01, 817., 824., 825., 826.,10108.,60808.,10106.,60606., 127., 128., 129., 827., 828., 829.,608.01, 408., 508., 815., 10808.,10811.,10812.,10820.,10106.,10107.,10116.,10606.,10607., 10608.,10708.,60816.,61616.,70708.,70808.,80814.,80816.,1010106.,1010107.,1010606.,101010106.,101010114., 614.,60614.,60607., 6060707.,6060607., 839., 840., 857.]
    nelion_mol = [240, 246, 252, 258, 264, 270, 276, 282, 288, 294, 300, 306, 312, 318, 324, 330, 336, 342, 348, 354, 360, 366, 372, 378, 384, 390, 396, 402, 408, 414, 420, 426, 432, 438, 444, 450, 456, 462, 468, 474, 480, 486, 492, 498, 504, 510, 516, 522, 528, 534, 540, 546, 552, 558, 564, 570, 576, 582, 588, 594, 600, 606, 612, 618, 624, 630, 636, 642, 648, 654, 660, 666, 672, 678, 684, 690, 696, 702, 708, 714, 720, 726, 732, 738, 744, 750, 756, 762, 768, 774, 780, 786, 792]

    # Some NELION indices are reserved for "special ions", i.e. ions with charges larger than 5 which are supported for some elements
    special_ions = [299, 305, 311, 317, 323, 329, 335, 341, 347, 359, 365, 371, 377, 383, 389, 395, 401, 407, 419, 425, 431, 437, 443, 449, 455, 461, 467, 479, 485, 491, 497, 503, 509, 515, 521, 527]

    # If we are dealing with a molecule, apply the direct mapping defined above directly
    if code > 100:
        return nelion_mol[np.where(np.array(code_mol) == code)[0][0]]

    # If we are dealing with atomic/ionic species, first compute NELION according to the general conversion formula
    Z = int(np.floor(code))
    charge = int(np.round((code - Z) * 100.0))
    nelion = Z * 6 - 6 + (charge + 1)
    if charge > 5:
        # If the charge is between 6 and 10, a different ("special ion") conversion formula must be used instead which is supported
        # only for elements with 9 < Z < 29
        if Z > 19 and Z < 29 and charge < 10:
            nelion = 6 * (Z + charge * 10 - 30) - 1
            assert nelion in special_ions # Double check that the "special ion" formula does indeed land us in one of the "special ion" NELIONs
        else:
            # Ions with charges over 10, or ions with charges over 5 that do not satisfy the 9 < Z < 29 condition are not supported
            raise ValueError('Species {} is not supported by XNFPELSYN/SYNTHE'.format(code))
    else:
        # If the general conversion formula lands us on a NELION which is reserved for molecules or special ions, this species is not supported
        if nelion in special_ions + nelion_mol:
            raise ValueError('Species {} is not supported by XNFPELSYN/SYNTHE'.format(code))
    return nelion

def build_BasicATLAS_linelist(output_dir, wl_start, wl_end, res, C12C13, atoms, air_wl):
    """Compile a SYNTHE-formatted linelist for spectral synthesis using the internal machinery of BasicATLAS. The metadata of the linelist also
    fixes the wavelength sampling to be used in spectral synthesis. Isotope abundance ratios are fixed at this stage as well, as SYNTHE itself
    does not recognize lines by isotopes and, instead, the compiled linelist simply rescales the oscillator strengths according to the isotope
    ratios

    BasicATLAS uses the SYNBEG, RGFALLLINESNEW, RMOLECASC and RH2OFAST codes to produce three files which will be saved in `output_dir`:
        - fort.12 is the main (LTE) linelist
        - fort.19 is the NLTE linelist, which typically only includes hydrogen lines. Note that SYNTHE does not support NLTE synthesis, so both
          linelists are treated identically
        - fort.93 is the meta data, which include the total number of lines, wavelength sampling parameters, air/vacuum setting and the line
          strength cutoff in SYNTHE
    
    Parameters
    ----------
    output_dir : str
        Output directory for the compiled linelist (fort.12, fort.19 and fort.93). The directory must not already exist
    wl_start : number
        Starting wavelength for spectral synthesis (in nm)
    wl_end : number
        Final wavelength for spectral synthesis (in nm)
    res : number
        Native resolution of spectral synthesis (wavelength / step in wavelength)
    C12C13 : number
        C12-C13 isotope abundance ratio which will be "baked" into CH lines
    atoms : str
        Atomic linelist to use. Must be one of the linelists supported by BasicATLAS. The default linelist, "BasicATLAS", is recommended
    air_wl : bool
        Vacuum (if `False`) or air (if `True`) wavelengths?
    """
    # Prepare the output directory
    if os.path.isdir(output_dir):
        raise ValueError('{} already exists'.format(output_dir))
    os.mkdir(output_dir)

    # Extract the line list calculation script template from BasicATLAS
    script = atlas.templates.synthe_control
    start = script.find('# synbeg.exe initializes the computation')
    end = script.find('# synthe.exe computes line opacities')
    if start == -1 or end == -1:
        raise ValueError('LOCALFIT no longer compatible with BasicATLAS')
    script = 'cd {output}\n' + script[start:end] + '\nrm *.out *.dat fort.14 fort.20'

    # Fill out the template
    C13 = 1 / (C12C13 + 1)
    C12 = 1 - C13
    C12C13_cmd = 'echo "{} {}" > c12c13.dat'.format(np.log10(C12), np.log10(C13))
    atoms_linelist = os.path.realpath(atlas.python_path + '/data/synthe_files/{}.dat'.format(atoms))
    if not os.path.isfile(atoms_linelist):
        raise ValueError('Linelist {} not found'.format(atoms_linelist))
    cards = {
        's_files': atlas.python_path + '/data/synthe_files/',
        'd_files': atlas.python_path + '/data/dfsynthe_files/',
        'synthe_suite': atlas.python_path + '/bin/',
        'airorvac': ['VAC', 'AIR'][air_wl],
        'wlbeg': wl_start,
        'wlend': wl_end,
        'resolu': res,
        'turbv': 2.0,
        'linelist': atoms_linelist,
        'C12C13': C12C13_cmd,
        'ifnlte': 0,
        'linout': -1,
        'cutoff': 0.0001,
        'ifpred': 1,
        'nread': 0,
        'output': os.path.realpath(output_dir),
    }
    script = script.format(**cards)

    # Run the script and make sure it ran without errors
    result = subprocess.run(['bash', '-c', script], text = True, capture_output = True)
    if result.stderr != '':
        f = open(script_dump := ('{}/script.sh'.format(output_dir)), 'w')
        f.write(script)
        f.close()
        print(result.stderr, file = sys.stderr)
        raise ValueError('BasicATLAS produced an error when calculating the linelist. The bash script was saved in {}'.format(script_dump))
    if result.stdout != '':
        print(result.stdout, file = sys.stderr)

def read_grid_dimensions():
    """Normally, this function returns the available dimensions in the model grid and the grid points available in those dimensions;
    however, LOCALFIT does not use a grid of pre-computed models. Instead, this function returns a placeholder grid centered on the
    input stellar parameters (`settings['gridfit_params']`) with grid points determined based on parameter offsets in
    `settings['gridfit_offsets']`
    
    This function also pre-computes the linelist for the required combination of `settings['gridfit_params']` and
    `settings['gridfit_offsets']` values and loads it into memory. Since the C12-C13 isotope abundance ratio is "baked in" the linelist,
    we do not account for variations in the isotope ratio with stellar parameters (instead, the C12-C13 ratio corresponding to the
    central stellar parameters in `settings['gridfit_params']` will be used in all spectral synthesis)
    
    Returns
    -------
    dict
        Placeholder grid points. Each element corresponds to a grid dimension (e.g. teff, logg etc), and the values are lists that
        contain the central value and the offset values
    """
    global ATLAS_restarts_header, linelists

    # Check that the GRIDFIT stellar parameters are fully specified
    if (type(settings['gridfit_params']) is bool) or (settings['gridfit_params'].keys() != set(['teff', 'logg', 'zscale', 'alpha', 'carbon'])):
        raise ValueError('Provide stellar parameters inferred from GRIDFIT (teff,logg,zscale,alpha,carbon) in settings[\'gridfit_params\']')

    # Set up the grid points
    grid = {}
    for param in settings['gridfit_params']:
        if param not in settings['gridfit_offsets']:
            grid[param] = np.array([settings['gridfit_params'][param]])
        else:
            grid[param] = settings['gridfit_params'][param] + np.array(settings['gridfit_offsets'][param])
            grid[param][grid[param] < np.min(ATLAS_restarts_header[param])] = np.min(ATLAS_restarts_header[param])
            grid[param][grid[param] > np.max(ATLAS_restarts_header[param])] = np.max(ATLAS_restarts_header[param])
            grid[param] = np.unique(grid[param])
            if np.min(grid[param]) > settings['gridfit_params'][param] or np.max(grid[param]) < settings['gridfit_params'][param]:
                raise ValueError('The configured gridfit offsets for {} do not accommodate the initial value'.format(param))

    # Generate the master linelist on disk
    grid_id = generate_hash(grid)
    C12C13 = C12C13_kirby_2015(settings['gridfit_params']['logg'])
    if os.path.isdir(linelist_dir := ('{}/linelist_{}'.format(settings['scratch'], grid_id))):
        notify('Linelist already exists in {}'.format(linelist_dir))
    else:
        notify('Compiling linelist in {}'.format(linelist_dir), color = 'y')
        try:
            build_BasicATLAS_linelist(linelist_dir, settings['wl_start'], settings['wl_end'], settings['res'], C12C13, settings['atoms'], settings['air_wl'])
        except:
            shutil.rmtree(linelist_dir)
            raise
        notify('Linelist compiled!', color = 'g')

    # Load the master linelist into memory and sort it by element
    notify('Loading linelist...'.format(linelist_dir), color = 'y')
    linelist = PyTLAS.load_linelist(linelist_dir)
    linelists = {}
    null_mask = np.full(len(linelist[0]), True)
    for element in settings['elements']:
        element_mask = np.isin(linelist[0]['f3'], [code_to_nelion(code) for code in settings['elements'][element]])
        element_meta = copy.deepcopy(linelist[2])
        element_meta['n_lines'] = np.count_nonzero(element_mask)
        element_meta['n_lines_f19'] = 0
        null_mask &= ~element_mask
        linelists[element] = (linelist[0][element_mask], linelist[1][:0], element_meta)
    linelists['null'] = (linelist[0][null_mask], linelist[1], linelist[2])
    linelists['null'][2]['n_lines'] = np.count_nonzero(null_mask)
    notify('Linelist loaded!'.format(linelist_dir), color = 'g')

    return grid

def read_grid_model(params, grid):
    """Prepare a given grid point for abundance fitting by pre-computing the model atmosphere, null-opacities and null-spectra

    The null-spectrum corresponds to the synthetic spectrum of the model with no response-function corrections. The chemistry
    of the null-spectrum corresponds to the stellar parameters at the grid point (zscale, alpha, carbon) as well as abundance
    offsets in `settings['initial_abundance_offsets']`. The null-opacities are the opacities due to individual elements (as
    defined in `settings['elements']` in the null-spectrum)

    In order to handle redshift, the function will trim the wavelength range of the model spectrum on both sides to make sure
    that the resulting wavelength range remains within the model coverage at all redshifts between the bounds in
    `settings['virtual_dof']['redshift']`

    Parameters
    ----------
    params : dict
        Dictionary of stellar parameters that define the model grid point of interest
    grid   : dict
        Model grid dimensions, previously obtained with `read_grid_dimensions()`

    Returns
    -------
    wl : array_like
        Grid of model wavelengths in A
    flux : array_like
        Corresponding flux densities
    meta : dict
        Dictionary of model meta data with the following keys:
            asynth: 2D null-opacities of individual elements, as well as the total opacity of the null-spectrum in 'null' sub-key
            element_masks: boolean masks that select the wavelength points affected by each element
            response: placeholder array to store response functions which will be computed on demand by `preprocess_grid_model()`
            model_indices: list of integer indices that uniquely identify this model grid point
            gridfit_abun: abundance offsets of individual elements that correspond to the chemistry of stellar parameters alone
            null_abun: abundance offsets of the null-spectrum, which are `gridfit_abun` + `settings['initial_abundance_offsets']`
            structure: model atmosphere at this grid point as a Kurucz-formatted text file loaded as byte array
            null_xnfpelsyn: XNFPELSYN output of the null-spectrum, which most importantly includes continuum opacities of the
                            null-spectrum. We do not update continuum opacities when fitting response functions
            null_spectrum: Full (without redshift trimming) null-spectrum, including the wavelength grid, the corresponding flux
                           densities, continuum flux densities and the continuum-normalized flux
    """
    global ATLAS_restarts_header, linelists, xnfpelsyn, spectrv, synthe
    model_indices = tuple([list(grid[param]).index(params[param]) for param in sorted(list(params.keys()))])

    # Generate the ATLAS settings object for the structure model
    ATLAS_settings = atlas.Settings()
    ATLAS_settings.teff = np.round(params['teff'], 0)
    ATLAS_settings.logg = np.round(params['logg'], 3)
    ATLAS_settings.zscale = np.round(params['zscale'], 3)
    ATLAS_settings.Y = 0.245
    ATLAS_settings.abun = {element: np.round(params['alpha'], 2) for element in settings['alpha_elements']}
    ATLAS_settings.abun['C'] = np.round(params['carbon'] + ATLAS_restarts_header['carbon_map']([params['zscale'], params['logg']])[0], 2)
    if settings['use_initial_abundance_offsets_in_structure']:
        for element in settings['initial_abundance_offsets']:
            if element not in ATLAS_settings.abun:
                ATLAS_settings.abun[element] = 0.0
            ATLAS_settings.abun[element] = np.round(ATLAS_settings.abun[element] + settings['initial_abundance_offsets'][element], 2)

    # Define a unique name for this structure
    structure_id = '{}_{}'.format(generate_hash([ATLAS_settings.teff, ATLAS_settings.logg, ATLAS_settings.zscale, ATLAS_settings.abun]), ['interpolated', 'original'][settings['compute_new_structure']])
    structure_dir = '{}/structure_{}'.format(settings['scratch'], structure_id)
    notify('({}) Will use structure {}'.format(model_indices, structure_id))

    # If the structure already exists, do nothing
    if os.path.isdir(structure_dir):
        notify('({}) Structure already exists in {}'.format(model_indices, structure_dir))

    # Compute new structure
    elif settings['compute_new_structure']:
        ODF_id = generate_hash([ATLAS_settings.zscale, ATLAS_settings.abun])
        ODF_dir = '{}/ODF_{}'.format(settings['scratch'], ODF_id)
        if os.path.isdir(ODF_dir):
            notify('({}) ODF already exists in {}'.format(model_indices, ODF_dir))
        else:
            notify('({}) Calculating ODF in {}'.format(model_indices, ODF_dir), color = 'y')
            try:
                atlas.dfsynthe(ODF_dir, settings = ATLAS_settings, parallel = settings['ODF_multithreading'], silent = True)
            except:
                shutil.rmtree(ODF_dir)
                raise
            notify('({}) ODF calculated!'.format(model_indices), color = 'g')
        notify('({}) Calculating model atmosphere in {}'.format(model_indices, structure_dir), color = 'y')
        try:
            atlas.atlas(structure_dir, settings = ATLAS_settings, ODF = ODF_dir, silent = True)
        except:
            shutil.rmtree(structure_dir)
            raise
        structure = atlas.read_structure(structure_dir)
        max_err = np.max(np.abs(structure[0]['flux_error']))
        max_de = np.max(np.abs(structure[0]['flux_error_derivative']))
        conv_classes = {0: 'GOLD', 1: 'SILVER', 2: 'BRONZE', 3: 'UNCONVERGED'}
        if max_err < 1.0 and max_de < 10.0:
            conv = 0
        elif max_err < 10.0 and max_de < 100.0:
            conv = 1
        elif max_err < 1000.0:
            conv = 2
        else:
            conv = 3
        notify('({}) Model atmosphere calculated! Convergence: {} '.format(model_indices, conv_classes[conv]), color = 'g')
        if conv > dict((v,k) for k,v in conv_classes.items())[settings['min_convergence']]:
            shutil.rmtree(structure_dir)
            raise ValueError('Stopping because the model convergence class {} is worse than the minimum acceptable class in the settings ({})'.format(conv_classes[conv], settings['min_convergence']))

    # Interpolate the structure from the grid
    else:
        if settings['use_initial_abundance_offsets_in_structure']:
            raise ValueError('Cannot use the use_initial_abundance_offsets_in_structure flag with an interpolated structure')
        notify('({}) Building an interpolated structure in {}'.format(model_indices, structure_dir), color = 'y')
        args = {'teff': ATLAS_settings.teff, 'logg': ATLAS_settings.logg, 'zscale': ATLAS_settings.zscale, 'alpha': ATLAS_settings.abun['Mg'], 'carbon': np.round(params['carbon'], 2)}
        structure = atlas.restarts.interpolate_structure(args, header = ATLAS_restarts_header)
        model = atlas.restarts.generate_model(*structure)
        os.mkdir(structure_dir)
        f = open('{}/output_summary.out'.format(structure_dir), 'w'); f.write(model); f.close()
        f = open('{}/output_last_iteration.out'.format(structure_dir), 'w'); f.close()
        f = open('{}/output_main.out'.format(structure_dir), 'w'); f.close()
        notify('({}) Interpolated structure calculated!'.format(model_indices), color = 'g')

    # Generate a structure for model meta data
    meta = {'asynth': {'null': np.zeros([linelists['null'][2]['n_wl'], 72])}, 'element_masks': {}, 'response': {}, 'model_indices': model_indices}
    for element in settings['elements']:
        meta['response'][element] = {value: False for value in settings['abun']}

    # Generate null-opacities
    notify('({}) Computing null-opacities'.format(model_indices), color = 'y')
    meta['gridfit_abun'] = {element: params['alpha'] for element in settings['alpha_elements']}
    meta['gridfit_abun']['C'] = params['carbon'] + ATLAS_restarts_header['carbon_map']([params['zscale'], params['logg']])[0]
    meta['null_abun'] = copy.deepcopy(meta['gridfit_abun'])
    for element in settings['elements']:
        if element not in meta['null_abun']:
            meta['null_abun'][element] = 0.0
        if element in settings['initial_abundance_offsets']:
            meta['null_abun'][element] += settings['initial_abundance_offsets'][element]
    xnfpelsyn.load_structure('{}/output_summary.out'.format(structure_dir))
    meta['structure'] = xnfpelsyn.f5.copy()
    xnfpelsyn.update_abun(params['zscale'], meta['null_abun'], Y = 0.245, std_round = False)
    xnfpelsyn.run()
    meta['null_xnfpelsyn'] = copy.deepcopy(xnfpelsyn.xnfpelsyn_output)
    synthe.load_xnfpelsyn(xnfpelsyn)
    for element in list(settings['elements']) + ['null']:
        synthe.load_linelist(*linelists[element], VTURB_LOGG(params['logg']))
        synthe.run()
        meta['asynth']['null'] += synthe.asynth
        if element != 'null':
            meta['element_masks'][element] = np.any(synthe.asynth, axis = 1)
            meta['asynth'][element] = synthe.asynth[meta['element_masks'][element]]
    notify('({}) Null-opacities computed!'.format(model_indices), color = 'g')

    # Generate the null-spectrum
    notify('({}) Computing null-spectrum'.format(model_indices), color = 'y')
    synthe.load_linelist(*linelists['null'], VTURB_LOGG(params['logg']))
    synthe.asynth[:,:] = meta['asynth']['null']
    spectrv.load_xnfpelsyn(xnfpelsyn)
    spectrv.load_synthe(synthe)
    spectrv.run()
    meta['null_spectrum'] = spectrv.get_spectrum()
    notify('({}) Null-spectrum computed!'.format(model_indices), color = 'g')

    # Trim the spectrum on both sides to make sure we can do redshift corrections
    wl_range = [np.min(meta['null_spectrum'][0] * (1 + settings['virtual_dof']['redshift'][1] * 1e3 / scp.constants.c)), np.max(meta['null_spectrum'][0] * (1 + settings['virtual_dof']['redshift'][0] * 1e3 / scp.constants.c))]
    mask_left = meta['null_spectrum'][0] < wl_range[0]; mask_right = meta['null_spectrum'][0] > wl_range[1]; mask_in = (~mask_left) & (~mask_right)

    return meta['null_spectrum'][0][mask_in], meta['null_spectrum'][3][mask_in], meta

def preprocess_grid_model(wl, flux, params, meta):
    """Apply response functions and redshift corrections to a synthetic spectrum

    Response functions are calculated dynamically in real time. When the null-spectrum is first generated in `read_grid_model()`,
    a placeholder dictionary is defined in `meta` to store response functions. This function determines which response functions
    are necessary for the chemical abundances requested by the fitter, calculates any response functions that are missing and then
    assembles the final spectrum

    Parameters
    ----------
    wl : array_like
        Grid of model wavelengths in A (trimmed to accommodate all redshifts)
    flux : array_like
        Corresponding flux densities. This flux array is however unused, as it is immediately replaced by the full (without
        redshift trimming) flux array stored in `meta`
    params : dict
        Parameters of the model, including desired chemical abundances to be reproduced with response functions, and redshift
    meta : dict
        Null-spectrum, null-opacities, already computed response functions and other data, as defined in `read_grid_model()`

    Returns
    -------
    array_like
        Redshifted flux density with applied response functions
    """
    global ATLAS_restarts_header, linelists, xnfpelsyn, spectrv, synthe

    # Recover the full flux
    flux = meta['null_spectrum'][3] * 1.0

    # Determine which response functions are needed to generate this model
    required = []
    sides = {}
    for element in settings['elements']:
        assert (params[element] <= np.max(settings['abun'])) and (params[element] >= np.min(settings['abun']))
        if params[element] in settings['abun']:
            required += [[element, params[element]]]
            sides[element] = [params[element]]
        else:
            required += [[element, np.max(settings['abun'][settings['abun'] < params[element]])]]
            required += [[element, np.min(settings['abun'][settings['abun'] > params[element]])]]
            sides[element] = [required[-2][1], required[-1][1]]

    # Determine which of those are yet to be calculated
    to_compute = []
    for response in required:
        if type(meta['response'][response[0]][response[1]]) is bool:
            to_compute += [response]

    # Calculate the missing response functions
    if len(to_compute) != 0:
        notify('({}) Computing {} response function{} ({})'.format(meta['model_indices'], len(to_compute), ['', 's'][int(len(to_compute) > 1)], ','.join(['{}={}'.format(*response) for response in to_compute])), color = 'y')
        for response in to_compute:
            abun = copy.deepcopy(meta['null_abun'])
            abun[response[0]] = response[1]
            xnfpelsyn.f5[:,:] = meta['structure']
            xnfpelsyn.update_abun(params['zscale'], abun, Y = 0.245, std_round = False)
            xnfpelsyn.run()
            synthe.load_xnfpelsyn(xnfpelsyn)
            synthe.load_linelist(*linelists[response[0]], VTURB_LOGG(params['logg']))
            synthe.run()
            synthe.asynth[meta['element_masks'][response[0]]] += meta['asynth']['null'][meta['element_masks'][response[0]]] - meta['asynth'][response[0]]
            xnfpelsyn.xnfpelsyn_output = meta['null_xnfpelsyn']
            spectrv.load_xnfpelsyn(xnfpelsyn)
            spectrv.load_synthe(synthe)
            spectrv.mask[:] = meta['element_masks'][response[0]]
            spectrv.run()
            spectrum = spectrv.get_spectrum()
            meta['response'][response[0]][response[1]] = spectrum[3] - flux[meta['element_masks'][response[0]]]
        notify('({}) Response functions computed!'.format(meta['model_indices']), color = 'g')

    # Assemble the spectrum
    for element in settings['elements']:
        if len(sides[element]) == 1:
            response = meta['response'][element][sides[element][0]]
        else:
            response_high = meta['response'][element][sides[element][1]]
            response_low = meta['response'][element][sides[element][0]]
            response = response_low + (response_high - response_low) * (params[element] - sides[element][0]) / (sides[element][1] - sides[element][0])
        flux[meta['element_masks'][element]] += response

    # Apply continuum
    flux *= meta['null_spectrum'][2]

    # Apply redshift
    wl_full = meta['null_spectrum'][0]
    wl_redshifted = wl_full * (1 + params['redshift'] * 1e3 / scp.constants.c)
    flux = np.interp(wl, wl_redshifted, flux)

    return flux

def extract_results(fit, propagate_gridfit = False, detector_wl = False):
    """Extract best-fit abundances and errors from the `chemfit.chemfit()` output
    
    This function re-expresses all abundances with respect to hydrogen ([X/H]), such that we do not need to worry about
    uncertainties in best-fit metallicity and alpha-enhancement

    If `propagate_gridfit` is `True`, the function will also update the errors in the best-fit abundances to include
    the contributions due to uncertainties in other fixed stellar parameters (i.e. teff, logg and carbon if any of
    those were not allowed to vary in `settings['gridfit_offsets']`)

    The error propagation is only approximate. It assumes local linearity of the parameter space, ignores statistical
    and photometric priors used in GRIDFIT and assumes that the derivatives of the objective function with repsect to
    stellar parameters do not depend on the abundances of individual elements. It also assumes that the exact same
    spectrum and masking were used by GRIDFIT to obtain the stellar parameters, as the ones used in LOCALFIT
    
    Parameters
    ----------
    fit : dict
        Output of `chemfit.chemfit()`
    propagate_gridfit : bool, optional
        Set to `True` to propagate the uncertainties in Teff, log(g) and [C/M] (if any of those were fixed during the fit)
        into the errors in best-fit abundances. Defaults to `False`
    detector_wl : dict, optional
        If `propagate_gridfit` is set to `True`, this argument must be set to the observed wavelengths, i.e. the `wl`
        argument given to `chemfit.chemfit()`
    
    Returns
    -------
    dict
        Modified `fit` dictionary with the additional `['abun']` key that stores the following data in
        its sub-keys:
            abund   : dictionary of best-fit abundances with respect to hydrogen
            errors  : dictionary of errors in best-fit abundances
            dof     : order of parameters in the covariance matrix. Only provided if `propagate_gridfit`
                      is `True`
            cov     : covariance matrix from the Jacobian. Only provided if `propagate_gridfit` is `True`
    """
    global ATLAS_restarts_header, linelists, xnfpelsyn, spectrv, synthe

    # Re-express all abundances with respect to hydrogen
    fit['extra']['abun'] = {'abun': {}, 'errors': {}}
    for element in settings['elements']:
        ratio = '[{}/H]'.format(element)
        fit['extra']['abun']['abun'][ratio] = fit['fit'][element]
        fit['extra']['abun']['errors'][ratio] = fit['errors'][element]
        fit['extra']['abun']['abun'][ratio] += settings['gridfit_params']['zscale']
        if element in settings['alpha_elements']:
            fit['extra']['abun']['abun'][ratio] += settings['gridfit_params']['alpha']

    # Propagate the errors in teff, logg and carbon into the abundances
    if propagate_gridfit:
        notify('*** Propagating gridfit uncertainties into best-fit abundances ***', color = 'm')

        # Helper function to adapt model spectra to observations (downsampling, rebinning, masking, continuum-correction)
        def observe(wl, flux, detector_wl):
            ds_wl, ds_flux = main__simulate_observation(wl, flux, detector_wl = detector_wl)
            cont = main__estimate_continuum(ds_wl, fit['extra']['observed']['flux'] / ds_flux, fit['extra']['observed']['ivar'] * ds_flux ** 2, npix = settings['cont_pix'], k = settings['spline_order'], arm_index = fit['extra']['arm_index'])
            return (cont * ds_flux)[fit['extra']['mask']]

        # Which degrees of freedom do we need to add to the Jacobian?
        added_dof = [dof for dof in ['teff', 'logg', 'carbon'] if dof not in fit['extra']['fit']['dof']]

        # Compute derivatives for those degrees of freedom
        if len(added_dof) != 0:
            notify('Computing derivative{} in {}'.format(['', 's'][int(len(added_dof) > 1)], ','.join(added_dof)), color = 'y')

            # Generate the all-inclusive line list
            f12 = np.concatenate([linelists[element][0] for element in linelists])
            f19 = np.concatenate([linelists[element][1] for element in linelists])
            meta = copy.deepcopy(linelists['null'][2])
            meta['n_lines'] = len(f12)
            meta['n_lines_f19'] = len(f19)

            for dof in ['nominal'] + added_dof:
                params = {param: fit['fit'][param] for param in settings['gridfit_params']}
                steps = {'teff': 10, 'logg': 0.1, 'carbon': 0.1}
                if dof != 'nominal':
                    if params[dof] + steps[dof] <= np.max(ATLAS_restarts_header[dof]):
                        step = steps[dof]
                    else:
                        step = -steps[dof]
                    params[dof] += step
                ATLAS_settings = atlas.Settings()
                ATLAS_settings.teff = np.round(params['teff'], 0)
                ATLAS_settings.logg = np.round(params['logg'], 3)
                ATLAS_settings.zscale = np.round(params['zscale'], 3)
                ATLAS_settings.Y = 0.245
                ATLAS_settings.abun = {element: np.round(params['alpha'], 2) for element in settings['alpha_elements']}
                ATLAS_settings.abun['C'] = np.round(params['carbon'] + ATLAS_restarts_header['carbon_map']([params['zscale'], params['logg']])[0], 2)
                structure_id = '{}_interpolated'.format(generate_hash([ATLAS_settings.teff, ATLAS_settings.logg, ATLAS_settings.zscale, ATLAS_settings.abun]))
                structure_dir = '{}/structure_{}'.format(settings['scratch'], structure_id)
                args = {'teff': ATLAS_settings.teff, 'logg': ATLAS_settings.logg, 'zscale': ATLAS_settings.zscale, 'alpha': ATLAS_settings.abun['Mg'], 'carbon': np.round(params['carbon'], 2)}
                structure = atlas.restarts.interpolate_structure(args, header = ATLAS_restarts_header)
                if not os.path.isdir(structure_dir):
                    model = atlas.restarts.generate_model(*structure)
                    os.mkdir(structure_dir)
                    f = open('{}/output_summary.out'.format(structure_dir), 'w'); f.write(model); f.close()
                    f = open('{}/output_last_iteration.out'.format(structure_dir), 'w'); f.close()
                    f = open('{}/output_main.out'.format(structure_dir), 'w'); f.close()
                xnfpelsyn.load_structure('{}/output_summary.out'.format(structure_dir))
                xnfpelsyn.run()
                synthe.load_linelist(f12, f19, meta, VTURB_LOGG(params['logg']))
                synthe.load_xnfpelsyn(xnfpelsyn)
                synthe.run()
                spectrv.load_xnfpelsyn(xnfpelsyn)
                spectrv.load_synthe(synthe)
                spectrv.run()
                spectrum = observe(*spectrv.get_spectrum()[:2], detector_wl)
                if dof == 'nominal':
                    nominal = spectrum
                else:
                    fit['extra']['jac'] = np.hstack([fit['extra']['jac'], np.array([(spectrum - nominal) / step]).T])

            notify('Derivatives computed!', color = 'g')

        all_dof = fit['extra']['fit']['dof'].tolist() + added_dof

        # Now that we have the Jacobian matrix fully populated, we can compute the covariance matrix using the standard
        # error propagation formula:
        #            COV = (J^T x diag(IVAR) x J)^-1
        # We also scale the covariance matrix by the reduced chi-squared (SUM((OBSERVED - MODEL)^2 * IVAR) / (N_points - N_params))
        # to mimick scp.optimize.curve_fit()'s `absolute_sigma = False` setting

        # Compute reduced chi squared for the final fit
        residuals = (fit['extra']['observed']['flux'] - (fit['extra']['model']['cont'] * fit['extra']['model']['flux']))[fit['extra']['mask']] * fit['extra']['observed']['ivar'][fit['extra']['mask']] ** 0.5
        chi2_red = np.sum(residuals ** 2) / (np.shape(fit['extra']['jac'])[0] - np.shape(fit['extra']['jac'])[1])

        # Handle degrees of freedom with zero Jacobians
        singular = np.array([np.all(fit['extra']['jac'][:,i] == 0.0) for i in range(len(all_dof))])
        if np.any(singular):
            notify('The spectrum is insensitive to {}'.format(','.join(np.array(all_dof)[singular])), color = 'r')

        # Compute the covariance matrix and update the errors in abundances
        weighted_jacobian = fit['extra']['jac'][:,~singular] * (fit['extra']['observed']['ivar'] ** 0.5)[fit['extra']['mask']][:, np.newaxis]
        fit['extra']['abun']['cov'] = np.linalg.pinv(weighted_jacobian.T @ weighted_jacobian) * chi2_red
        fit['extra']['abun']['dof'] = np.array(all_dof)[~singular].tolist()
        for element in settings['elements']:
            if element in fit['extra']['abun']['dof']:
                fit['extra']['abun']['errors']['[{}/H]'.format(element)] = np.sqrt(fit['extra']['abun']['cov'][tuple([fit['extra']['abun']['dof'].index(element)] * 2)])
            else:
                fit['extra']['abun']['errors']['[{}/H]'.format(element)] = np.inf

    return fit

def public__localfit(wl, flux, ivar, gridfit, initial_abundance_offsets = {}, level = 1, niter = 5):
    """Determine the abundances of individual elements using response functions computed at runtime
    
    This is a wrapper around `chemfit.chemfit()`, which should be used instead of `chemfit.chemfit()`
    for the LOCALFIT preset
    
    Parameters
    ----------
    wl : dict
        Spectrum wavelengths keyed by spectrograph arm
    flux : dict
        Spectrum flux densities keyed by spectrograph arm
    ivar : dict
        Spectrum weights (inverted variances) keyed by spectrograph arm
    gridfit : dict
        Stellar parameters of the spectrum (required keys are 'teff', 'logg', 'zscale', 'alpha' and 'carbon') and
        the initial values for all virtual degrees of freedom that are not the abundances of individual elements
        (for now, it is just redshift)
    initial_abundance_offsets : dict, optional
        Dictionary of initial abundance offsets from the zscale+alpha+carbon nominal chemistry to be used in the
        fitter. The keys must only include elements from `settings['elements']`. Defaults to no offsets
    level : number, optional
        Level of analysis. In the order of increasing accuracy, the available levels are:
            1 : Adopt interpolated model atmospheres. Calculate the null-spectra of each model
                with no individual element enhancements. Fit for element abundances with response functions
            2 : Same as level 1, but compute new model atmospheres using DFSYNTHE/ATLAS
                instead of using interpolated structures
            3 : Same as level 1, but repeat the abundance determination `niter` times,
                updating the abundances of the null-spectrum to the output of the previous
                iteration
            4 : Same as level 3, but compute new model atmospheres instead of interpolating the grid
            5 : Same as level 4, but recompute the model atmospheres at the end of each
                iteration with the determined abundances
    niter : number, optional
        Number of iterations to carry out. Only relevant if the analysis level is 3, 4 or 5
    Returns
    -------
    dict
        Fitting results. The structure of the dictionary is the same as in the output of `chemfit.chemfit()` with
        one additional key: ['extra']['abun']. This key contains the best-fit element abundances ([X/H]) and their
        errors, as well as the full covariance matrix estimated with `extract_results(..., propagate_gridfit = False, ...)`.
        The intermediate abundances determined at the end of each iteration are also provided
    """
    global _global

    # Check that all basegrid parameters and all virtual parameters are present in `gridfit`
    basegrid_dof = ['teff', 'logg', 'zscale', 'alpha', 'carbon']
    required = basegrid_dof + [param for param in settings['virtual_dof'] if param not in settings['elements']]
    if set(required) != gridfit.keys():
        raise ValueError('Provided gridfit parameters {} do not match the requirement {}'.format(list(gridfit.keys()), required))

    if level not in [1, 2, 3, 4, 5]:
        raise ValueError('Unknown analysis level {}'.format(level))

    settings['gridfit_params'] = {param: gridfit[param] for param in basegrid_dof}
    settings['compute_new_structure'] = level in [2, 4, 5]
    settings['initial_abundance_offsets'] = {}
    settings['use_initial_abundance_offsets_in_structure'] = level in [2, 4, 5]

    # Apply initial offsets
    for element in initial_abundance_offsets:
        if element not in settings['elements']:
            raise ValueError('Unknown element {} in initial abundance offsets'.format(element))
        if initial_abundance_offsets[element] > np.max(settings['abun']) or initial_abundance_offsets[element] < np.min(settings['abun']):
            raise ValueError('Element {} in initial abundance offsets exceeds the allowed range'.format(element))
        settings['initial_abundance_offsets'][element] = initial_abundance_offsets[element]

    params = copy.deepcopy(gridfit)
    intermediate = [] # Storage for intermediate abundances at the end of each iteration

    actual_niter = [1, niter][level in [3, 4, 5]]
    for iteration in range(actual_niter):
        notify('*** Starting iteration {} ***'.format(iteration + 1), color = 'm')
        notify('gridfit_params={}\ninitial_abundance_offsets={}\ncompute_new_structure={} | use_initial_abundance_offsets_in_structure={}'.format(settings['gridfit_params'], settings['initial_abundance_offsets'], settings['compute_new_structure'], settings['use_initial_abundance_offsets_in_structure']), color = 'm')
        fit = main__chemfit(wl, flux, ivar, initial = params, method = ['gradient_descent', 'gradient_descent+jac'][int(iteration == actual_niter - 1)])
        if iteration != actual_niter - 1:
            for element in settings['elements']:
                settings['initial_abundance_offsets'][element] = fit['fit'][element]
            for param in settings['fit_dof']:
                params[param] = fit['fit'][param]
            settings['use_initial_abundance_offsets_in_structure'] = level == 5
        fit = extract_results(fit, propagate_gridfit = False)
        intermediate += [copy.deepcopy(fit['extra']['abun'])]
        notify('Completed iteration {}: {}\n'.format(iteration + 1, fit['extra']['abun']), color = 'm')
    notify('*** Completed all iterations ***\n')

    fit = extract_results(fit, propagate_gridfit = True, detector_wl = wl)
    fit['extra']['abun']['intermediate'] = intermediate
    del fit['extra']['jac']

    return fit

def public__get_global():
    global _global
    return _global