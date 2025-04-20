############################################################
#                                                          #
#                 LIVE SYNTHESIS PRESET                    #
#                                                          #
#  This preset configures chemfit to use runtime spectral  #
#  synthesis instead of interpolating an existing model    #
#  grid. This mode of operation incurs high computational  #
#  demand. It is recommended to use a reduced molecular    #
#  line list, e.g. by excluding TiO and water. Spectral    #
#  synthesis is carried out with BasicATLAS which must be  #
#  installed.                                              #
#                                                          #
############################################################

import pickle
import os, shutil
import zipfile
import sys
import random, string
import copy
import scipy as scp
import concurrent.futures

# This is the "defective" mask that removes parts of the spectrum where the models do not match the spectra of 
# Arcturus and the Sun well
defective_mask = [[3800, 4000], [4006, 4012], [4065, 4075], [4093, 4110], [4140, 4165], [4170, 4180], [4205, 4220],
                 [4285, 4300], [4335, 4345], [4375, 4387], [4700, 4715], [4775, 4790], [4855, 4865], [5055, 5065],
                 [5145, 5160], [5203, 5213], [5885, 5900], [6355, 6365], [6555, 6570], [7175, 7195], [7890, 7900],
                 [8320, 8330], [8490, 8505], [8530, 8555], [8650, 8672]]

settings = {
    ### All directories are specified in local settings ###
    'structures_dir': original_settings['structures_dir'],
    'structures_index': original_settings['structures_index'],
    'rescalc_path': original_settings['rescalc_path'],
    'scratch': original_settings['scratch'],

    ### Which elements to fit using live synthesis? [min, max, step] ###
    'elements': {'O': [-1.0, 1.0, 0.1], 'C': [-1.0, 1.0, 0.1], 'Fe': [-1.0, 1.0, 0.1]},

    ### Binning factor for new spectra ###
    'binning': 10,

    ### Suppress status messages? ###
    'silent': False,

    ### RESCALC settings ###
    'rescalc': {
        'atoms': 'BasicATLAS',
        'air_wl': True,
        'wl_start': 300,
        'wl_end': 1300,
        'res': 300000,
    },

    ### Work directory name for this session. `False` to randomize ###
    'workdir': False,

    ### Enable parallel computation of models (inherited from local settings) ###
    'parallel': original_settings['parallel'],

    ### Which parameters to fit? The abundances of individual elements enabled through live synthesis ###
    ### are not listed here as they will be added automatically ###
    'fit_dof': ['zscale', 'alpha', 'teff', 'logg', 'redshift'],

    ### Virtual grid bounds ###
    'virtual_dof': {'redshift': [-300, 300]},

    ### Default initial guesses ###
    'default_initial': {'redshift': 0.0},

    ### Fitting masks ###
    'masks': apply_standard_mask(defective_mask, original_settings['masks']),
}

# Import RESCALC
sys.path.append(settings['rescalc_path'])
import rescalc

# Relationship between C12/C13 isotope ratio and surface gravity from Kirby+2015
C12C13_kirby_2015 = lambda logg: np.where(logg > 2.7, 50, np.where(logg <= 2.0, 6, 63 * logg - 120))
# Relationship between VTURB and LOGG based on Deimos spectra from Gerasimov+2025
VTURB_LOGG = lambda logg: 2.792 * np.exp(-0.241 * logg -0.118)

# Work directory for live synthesis
livesyn_workdir = False

def notify(message):
    """Issue a status notification
    
    Wrapper for the print() function that respects the `silent` setting
    
    Parameters
    ----------
    message : str
        Message to print
    """
    if not settings['silent']:
        print(message, flush = True)

def init_livesyn():
    """Initialize the work directory for the session
    
    The work directory is where the model spectra created with real-time spectral synthesis are
    stored. Since there is likely little overlap in the models required to fit different spectra,
    the contents of the work directory are typically discarded at the end of the session. For instance,
    when running live synthesis fitting on a compute cluster, the work directory may be created in the
    local storage of the node that will be automatically discarded once the job is complete. Note that
    deletion of the work directory must be carried out by the system or the user: it will not be done
    by chemfit

    Since the contents of the work directory are often disposable, the work directory can have a random
    name. If `settings['workdir']` is set to `False`, the work directory will be created by this function
    in the scratch directory (`settings['scratch']`) with a random name

    If it is instead preferred to reuse the same models in multiple sessions, a previously used work
    directory (or an empty directory) can be provided in `settings['workdir']`, and all new models will be
    stored there instead, while the models that already exist in the directory will be made available to
    the fitter
    """
    global livesyn_workdir

    # Generate the work directory for this session
    if type(settings['workdir']) == bool:
        livesyn_workdir = '{}/livesyn_{}/'.format(settings['scratch'], ''.join(random.choices(string.ascii_letters + string.digits, k = 20)))
    else:
        livesyn_workdir = settings['workdir']
    if not os.path.isdir(livesyn_workdir):
        os.mkdir(livesyn_workdir)
        notify('LIVESYN session started in {}'.format(livesyn_workdir))
    else:
        notify('LIVESYN session resumed in {}'.format(livesyn_workdir))


def build_linelists():
    """Build the line lists for the fitting session

    The line lists are going to be the same for every fitting session so long as the RESCALC configuration
    is unchanged. For this reason, it is only necessary to precompute them once. This function will check
    for an existing line list in `settings['scratch'] + '/linelist'`, and use RESCALC to generate it if it
    does not exist. The line lists are generated for each gravity in the structures grid due to distinct
    C12/C13 isotope ratios and turbulent velocities (turbulent velocity strictly speaking does not change
    the line list; however, it changes the SYNTHE configuration file, fort.93, and that file is created
    at the same time as the line list)
    """
    global livesyn_grid

    # Check if linelist exists
    linelist_dir = settings['scratch'] + '/linelist'
    if os.path.isdir(linelist_dir):
        return
    os.mkdir(linelist_dir)
    notify('Generating linelist...')

    for logg in livesyn_grid['logg']:
        rescalc.settings = {
            **settings['rescalc'],
            'C12C13': np.round(C12C13_kirby_2015(logg), 2),
            'vturb': np.round(VTURB_LOGG(logg), 2),
        }
        linelist = '{}/logg_{}'.format(linelist_dir, logg)
        rescalc.build_linelist([], linelist, invert = True)
        notify('Linelist for logg={} done'.format(logg))


def structure_penalty(zscale_diff, alpha_diff, CO_diff, carbon_diff):
    """Calculate the flux error penalty for a selected structure
    
    When the chemical abundances are inconsistent between the model of atmospheric structure and spectral
    synthesis, the calculated spectrum may be inaccurate. In live synthesis mode, this effect is inevitable,
    since we draw structures from a pre-computed grid, which may not contain a structure with the exact abundances
    as those required by the fitter. Ideally, we would like to choose the closest structure to the desired
    abundances such that this error is minimized

    One clear symptom of the mismatch in abundances between the structure and the spectrum is the accumulation
    of flux errors in the radiative transfer equation, such that the integral of the calculated spectrum no
    longer matches the Stefan-Boltzmann prediction for the given effective temperature. In order to choose the
    best structure, we may seek to minimize this departure from the Stefan-Boltzmann law

    This function evaluates the expected errors in effective temperature for a given set of abundance offsets
    between the structure and the spectrum. These errors are drawn from two pre-computed tables: one for errors
    due to offsets in [C/H]/[O/H]; and another one for offsets in metallicity/alpha. These tables are interpolated,
    and then the two errors are added in quadrature

    The pre-computed tables were evaluated with BasicATLAS/ATLAS/SYNTHE for a star with Teff=4000 K, log(g)=1.5,
    [M/H]=-1.5. The synthetic spectra were integrated between 100 nm and 12 micron
    
    Parameters
    ----------
    zscale_diff : float
        Offset in metallicity excluding carbon and oxygen (spectrum - structure)
    alpha_diff : float
        Offset in alpha-enhancement excluding oxygen (spectrum - structure)
    CO_diff : float
        Offset in [C/H]-[O/H] (spectrum - structure)
    carbon_diff : float
        Offset in [C/H] (spectrum - structure)
    
    Returns
    -------
    float
        Expected offset in the integral of the spectrum (converted to temperature in K using Stefan-Boltzmann law) due
        to discrepant abundances. If the requested offsets exceed the pre-tabulated range, a penalty of infinity is returned
    """
    CO_penalty = [
        [9.6, 9.1, 8.5, 7.8, 7.1, 6.4, 5.6, 4.7, 3.8, 2.9, 2.0, 1.0, 0.1, -0.8, -1.7, -2.5, -3.2, -3.7, 1.8, 1.8, 2.1, 2.9, 4.3, 6.5, 9.8, 14.4, 20.7, 29.2, 40.5, 59.9, 78.6, 102.5, 136.4, 177.1, 65.9, 96.9, 137.9, 179.0, 230.8, 289.1, 359.8],
        [10.1, 9.6, 9.0, 8.4, 7.7, 7.0, 6.2, 5.4, 4.5, 3.6, 2.6, 1.7, 0.7, -0.2, -1.2, -2.1, -2.9, -3.6, -4.2, 1.3, 1.2, 1.5, 2.2, 3.6, 5.7, 9.0, 13.5, 19.7, 28.2, 39.3, 58.7, 77.3, 101.1, 135.0, 175.6, 63.4, 94.6, 135.5, 176.4, 228.0, 286.0],
        [10.5, 10.0, 9.5, 8.9, 8.3, 7.6, 6.8, 6.0, 5.2, 4.3, 3.3, 2.4, 1.4, 0.4, -0.6, -1.6, -2.5, -3.4, -4.1, -4.7, 0.7, 0.6, 0.8, 1.5, 2.8, 4.9, 8.0, 12.5, 18.7, 27.0, 38.1, 57.4, 75.9, 99.6, 133.5, 174.0, 60.6, 92.0, 132.8, 173.6, 224.9],
        [10.8, 10.4, 9.9, 9.4, 8.8, 8.1, 7.4, 6.6, 5.8, 4.9, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0, -2.0, -3.0, -3.9, -4.6, -5.3, 0.1, -0.1, 0.0, 0.7, 1.9, 3.9, 7.0, 11.4, 17.5, 25.7, 36.7, 56.0, 74.4, 98.1, 131.8, 172.3, 57.4, 89.1, 129.9, 170.5],
        [11.1, 10.7, 10.3, 9.8, 9.2, 8.6, 7.9, 7.2, 6.4, 5.5, 4.7, 3.7, 2.7, 1.7, 0.7, -0.4, -1.4, -2.5, -3.5, -4.4, -5.2, -5.9, -0.6, -0.9, -0.8, -0.2, 0.9, 2.9, 5.9, 10.2, 16.2, 24.4, 35.3, 54.5, 72.8, 96.4, 130.1, 170.5, 53.8, 85.9, 126.6],
        [11.3, 11.0, 10.6, 10.1, 9.6, 9.1, 8.4, 7.7, 7.0, 6.1, 5.3, 4.3, 3.4, 2.3, 1.3, 0.2, -0.9, -1.9, -3.0, -4.0, -5.0, -5.9, -6.7, -1.4, -1.7, -1.7, -1.2, -0.1, 1.7, 4.6, 8.8, 14.8, 22.9, 33.7, 52.9, 71.1, 98.8, 128.3, 168.7, 49.6, 82.3],
        [11.6, 11.2, 10.9, 10.5, 10.0, 9.5, 8.9, 8.2, 7.5, 6.7, 5.9, 5.0, 4.0, 3.0, 1.9, 0.9, -0.2, -1.4, -2.5, -3.6, -4.7, -5.7, -6.6, -7.5, -2.3, -2.7, -2.7, -2.3, -1.3, 0.4, 3.3, 7.4, 13.3, 21.3, 32.1, 51.1, 69.4, 97.0, 126.4, 166.7, 48.5],
        [11.7, 11.5, 11.1, 10.7, 10.3, 9.8, 9.3, 8.7, 8.0, 7.2, 6.4, 5.5, 4.6, 3.6, 2.6, 1.5, 0.4, -0.7, -1.9, -3.1, -4.2, -5.4, -6.4, -7.4, -8.3, -3.2, -3.7, -3.8, -3.5, -2.6, -0.9, 1.8, 5.9, 11.6, 19.6, 30.3, 49.3, 67.5, 95.1, 124.4, 164.6],
        [11.9, 11.6, 11.3, 11.0, 10.6, 10.1, 9.6, 9.1, 8.4, 7.7, 6.9, 6.1, 5.2, 4.2, 3.2, 2.1, 1.0, -0.1, -1.3, -2.5, -3.7, -4.9, -6.1, -7.3, -8.4, -9.3, -4.3, -4.8, -5.1, -4.9, -4.1, -2.4, 0.2, 4.2, 9.9, 17.8, 28.4, 47.4, 65.5, 93.0, 122.3],
        [12.0, 11.8, 11.5, 11.2, 10.8, 10.4, 9.9, 9.4, 8.8, 8.1, 7.4, 6.6, 5.7, 4.8, 3.8, 2.7, 1.6, 0.5, -0.7, -1.9, -3.2, -4.5, -5.7, -7.0, -8.2, -9.4, -10.4, -5.5, -6.1, -6.5, -6.3, -5.6, -4.1, -1.5, 2.4, 8.0, 15.9, 26.5, 45.4, 63.4, 90.9],
        [12.1, 11.9, 11.6, 11.4, 11.0, 10.6, 10.2, 9.7, 9.1, 8.5, 7.8, 7.0, 6.2, 5.3, 4.3, 3.3, 2.2, 1.1, -0.1, -1.3, -2.6, -3.9, -5.3, -6.6, -7.9, -9.3, -10.5, -5.8, -6.8, -7.6, -8.0, -8.0, -7.3, -5.9, -3.4, 0.4, 6.0, 13.8, 29.4, 43.2, 61.2],
        [12.2, 12.0, 11.8, 11.5, 11.2, 10.8, 10.4, 9.9, 9.4, 8.8, 8.2, 7.4, 6.6, 5.8, 4.8, 3.8, 2.8, 1.6, 0.5, -0.8, -2.1, -3.4, -4.8, -6.2, -7.6, -9.0, -10.4, -11.8, -7.2, -8.3, -9.1, -9.7, -9.7, -9.2, -7.9, -5.4, -1.7, 3.9, 11.6, 27.1, 40.9],
        [12.2, 12.1, 11.8, 11.6, 11.3, 11.0, 10.6, 10.1, 9.6, 9.1, 8.5, 7.8, 7.0, 6.2, 5.3, 4.3, 3.3, 2.2, 1.0, -0.2, -1.5, -2.9, -4.3, -5.7, -7.2, -8.7, -10.3, -11.8, -13.2, -8.7, -9.9, -10.9, -11.5, -11.7, -11.3, -10.0, -7.7, -3.9, 1.5, 9.2, 24.7],
        [12.3, 12.1, 11.9, 11.7, 11.4, 11.1, 10.7, 10.3, 9.8, 9.3, 8.7, 8.1, 7.3, 6.5, 5.7, 4.7, 3.7, 2.6, 1.5, 0.3, -1.0, -2.4, -3.8, -5.3, -6.8, -8.4, -10.0, -11.6, -13.3, -14.9, -10.5, -11.8, -12.9, -13.6, -13.9, -13.5, -12.4, -10.1, -6.4, -1.0, 6.6],
        [12.3, 12.1, 11.9, 11.7, 11.4, 11.1, 10.8, 10.4, 10.0, 9.5, 8.9, 8.3, 7.6, 6.8, 6.0, 5.1, 4.1, 3.1, 1.9, 0.7, -0.5, -1.9, -3.3, -4.8, -6.4, -8.1, -9.8, -11.5, -13.2, -15.0, -16.7, -12.5, -13.9, -15.1, -15.9, -16.3, -16.1, -15.0, -12.8, -9.2, -3.8],
        [12.3, 12.1, 11.9, 11.7, 11.5, 11.2, 10.9, 10.5, 10.1, 9.6, 9.0, 8.4, 7.8, 7.1, 6.3, 5.4, 4.4, 3.4, 2.3, 1.1, -0.1, -1.5, -2.9, -4.4, -6.1, -7.7, -9.5, -11.3, -13.2, -15.1, -16.9, -18.8, -14.7, -16.3, -17.6, -18.6, -19.1, -18.9, -17.9, -15.7, -12.2],
        [12.3, 12.1, 11.9, 11.7, 11.5, 11.2, 10.9, 10.5, 10.1, 9.6, 9.1, 8.5, 7.9, 7.2, 6.4, 5.6, 4.7, 3.7, 2.6, 1.4, 0.2, -1.2, -2.6, -4.1, -5.7, -7.5, -9.3, -11.2, -13.1, -15.1, -17.2, -19.2, -21.2, -17.3, -19.0, -20.5, -21.6, -22.2, -22.1, -21.2, -19.1],
        [12.3, 12.1, 11.9, 11.7, 11.4, 11.2, 10.8, 10.5, 10.1, 9.6, 9.1, 8.6, 8.0, 7.3, 6.5, 5.7, 4.8, 3.9, 2.8, 1.7, 0.4, -0.9, -2.4, -3.9, -5.5, -7.3, -9.1, -11.1, -13.1, -15.3, -17.4, -19.7, -21.9, -24.1, -20.3, -22.2, -23.9, -25.1, -25.8, -25.9, -25.0],
        [12.3, 12.1, 11.9, 11.6, 11.4, 11.1, 10.8, 10.4, 10.0, 9.6, 9.1, 8.5, 7.9, 7.2, 6.5, 5.7, 4.8, 3.9, 2.9, 1.7, 0.5, -0.8, -2.2, -3.8, -5.4, -7.2, -9.1, -11.1, -13.3, -15.5, -17.8, -20.2, -22.7, -25.1, -27.5, -24.0, -26.1, -27.9, -29.3, -30.2, -30.3],
        [12.2, 12.0, 11.8, 11.6, 11.3, 11.0, 10.7, 10.3, 9.9, 9.4, 8.9, 8.4, 7.8, 7.1, 6.4, 5.6, 4.7, 3.8, 2.7, 1.6, 0.4, -0.9, -2.3, -3.9, -5.6, -7.4, -9.3, -11.4, -13.7, -16.0, -18.5, -21.1, -23.7, -26.4, -29.1, -25.9, -28.5, -30.8, -32.8, -34.5, -35.5],
        [12.2, 12.0, 11.8, 11.5, 11.2, 10.9, 10.5, 10.1, 9.7, 9.2, 8.7, 8.1, 7.5, 6.8, 6.0, 5.2, 4.3, 3.4, 2.3, 1.2, 0.0, -1.3, -2.8, -4.4, -6.1, -7.9, -10.0, -12.1, -14.5, -17.0, -19.6, -22.4, -25.3, -28.3, -31.3, -34.3, -31.5, -34.4, -37.0, -39.4, -41.3],
        [12.1, 11.9, 11.7, 11.4, 11.1, 10.7, 10.3, 9.9, 9.4, 8.9, 8.3, 7.7, 7.0, 6.2, 5.4, 4.5, 3.6, 2.6, 1.5, 0.3, -1.0, -2.4, -3.9, -5.6, -7.4, -9.3, -11.5, -13.8, -16.3, -19.0, -21.9, -24.9, -28.2, -31.5, -35.0, -38.4, -42.0, -39.7, -43.1, -46.2, -49.0],
        [12.1, 11.8, 11.6, 11.2, 10.9, 10.5, 10.1, 9.6, 9.0, 8.4, 7.7, 7.0, 6.2, 5.3, 4.3, 3.2, 2.1, 0.9, -0.4, -1.9, -3.4, -5.0, -6.8, -8.8, -10.9, -13.2, -15.6, -18.3, -21.3, -24.4, -27.8, -31.4, -35.3, -39.3, -43.5, -47.9, -52.3, -51.2, -55.7, -60.2, -64.5],
        [12.0, 11.7, 11.4, 11.1, 10.7, 10.2, 9.7, 9.1, 8.4, 7.6, 6.7, 5.8, 4.6, 3.4, 1.9, 0.3, -1.5, -3.5, -5.8, -8.4, -11.4, -14.8, -18.6, -23.1, -28.1, -34.0, -40.7, -48.5, -57.5, -68.0, -80.0, -93.8, -109.6, -127.7, -148.4, -171.8, -198.4, -228.2, -253.5, -290.4, -331.3],
        [11.9, 11.6, 11.3, 10.9, 10.4, 9.9, 9.2, 8.5, 7.7, 6.7, 5.6, 4.3, 2.7, 1.0, -1.1, -3.5, -6.2, -9.4, -13.1, -17.3, -22.3, -28.2, -35.0, -43.1, -52.5, -63.6, -76.5, -91.4, -108.7, -128.7, -151.5, -177.6, -207.2, -240.7, -278.1, -319.9, -365.9, -416.4, -461.7, -520.4, -582.8],
        [11.9, 11.5, 11.2, 10.7, 10.2, 9.6, 8.9, 8.0, 7.0, 5.9, 4.5, 2.9, 1.1, -1.1, -3.6, -6.5, -10.0, -14.0, -18.7, -24.3, -30.8, -38.5, -47.5, -58.1, -70.5, -85.0, -101.7, -121.0, -143.3, -168.7, -197.7, -230.4, -267.2, -308.4, -353.9, -404.0, -458.5, -517.2, -579.8, -635.9, -704.0],
        [11.8, 11.5, 11.0, 10.6, 10.0, 9.3, 8.5, 7.6, 6.5, 5.2, 3.7, 1.9, -0.2, -2.6, -5.5, -8.9, -12.9, -17.6, -23.1, -29.6, -37.2, -46.3, -56.9, -69.3, -83.8, -100.5, -119.9, -142.2, -167.7, -196.7, -229.6, -266.6, -307.8, -353.6, -404.0, -458.8, -517.9, -580.9, -647.1, -706.1, -775.9],
        [11.8, 11.4, 11.0, 10.4, 9.8, 9.1, 8.3, 7.3, 6.1, 4.7, 3.0, 1.1, -1.2, -3.9, -7.0, -10.8, -15.2, -20.3, -26.4, -33.6, -42.2, -52.2, -63.9, -77.7, -93.6, -112.1, -133.3, -157.7, -185.5, -217.1, -252.6, -292.4, -336.7, -385.6, -439.0, -496.8, -558.8, -624.2, -692.3, -752.5, -822.8],
        [11.7, 11.3, 10.9, 10.4, 9.7, 9.0, 8.1, 7.0, 5.8, 4.3, 2.5, 0.5, -2.0, -4.8, -8.2, -12.2, -16.9, -22.4, -29.0, -36.8, -45.9, -56.7, -69.3, -84.0, -101.1, -120.8, -143.4, -169.3, -198.8, -232.1, -269.6, -311.5, -357.9, -408.9, -464.4, -524.2, -587.9, -654.8, -723.9, -794.1, -855.0],
        [11.7, 11.3, 10.8, 10.3, 9.6, 8.8, 7.9, 6.8, 5.5, 4.0, 2.1, -0.0, -2.6, -5.6, -9.1, -13.3, -18.2, -24.0, -31.0, -39.1, -48.8, -60.1, -73.4, -88.9, -106.8, -127.4, -151.1, -178.1, -208.8, -243.5, -282.4, -325.7, -373.6, -426.2, -483.2, -544.4, -609.3, -677.1, -746.8, -817.1, -878.0],
        [11.7, 11.3, 10.8, 10.2, 9.6, 8.7, 7.8, 6.7, 5.3, 3.7, 1.8, -0.4, -3.0, -6.1, -9.8, -14.1, -19.2, -25.3, -32.5, -41.0, -51.0, -62.8, -76.6, -92.6, -111.1, -132.5, -156.9, -184.8, -216.5, -252.1, -292.1, -336.5, -385.6, -439.2, -497.2, -559.4, -625.2, -693.6, -763.7, -834.0, -894.6],
        [11.7, 11.2, 10.8, 10.2, 9.5, 8.7, 7.7, 6.5, 5.2, 3.5, 1.6, -0.7, -3.4, -6.6, -10.3, -14.8, -20.0, -26.3, -33.6, -42.4, -52.7, -64.9, -79.0, -95.5, -114.4, -136.3, -161.4, -190.0, -222.3, -258.7, -299.5, -344.7, -394.6, -449.0, -507.9, -570.8, -637.1, -706.0, -776.3, -846.6, -907.0],
        [11.6, 11.2, 10.7, 10.2, 9.5, 8.6, 7.6, 6.4, 5.0, 3.4, 1.4, -0.9, -3.7, -6.9, -10.7, -15.3, -20.6, -27.0, -34.5, -43.5, -54.0, -66.4, -80.9, -97.6, -117.0, -139.3, -164.8, -193.9, -226.8, -263.8, -305.2, -351.0, -401.5, -456.6, -516.0, -579.5, -646.2, -715.4, -785.8, -856.0, -916.2],
        [11.6, 11.2, 10.7, 10.1, 9.4, 8.6, 7.6, 6.4, 4.9, 3.3, 1.3, -1.1, -3.9, -7.2, -11.0, -15.6, -21.1, -27.6, -35.2, -44.3, -55.0, -67.6, -82.3, -99.3, -119.0, -141.6, -167.5, -196.9, -230.2, -267.7, -309.5, -355.9, -406.8, -462.4, -522.2, -586.1, -653.1, -722.5, -793.0, -863.1, -923.2],
        [11.6, 11.2, 10.7, 10.1, 9.4, 8.5, 7.5, 6.3, 4.9, 3.2, 1.1, -1.2, -4.1, -7.4, -11.3, -15.9, -21.5, -28.0, -35.8, -45.0, -55.8, -68.5, -83.4, -100.6, -120.5, -143.3, -169.5, -199.2, -232.9, -270.7, -312.9, -359.6, -410.9, -466.8, -527.0, -591.1, -658.4, -727.9, -798.4, -868.5, -936.7],
        [11.6, 11.2, 10.7, 10.1, 9.4, 8.5, 7.5, 6.3, 4.8, 3.1, 1.1, -1.3, -4.2, -7.5, -11.5, -16.2, -21.7, -28.3, -36.2, -45.4, -56.4, -69.2, -84.2, -101.6, -121.6, -144.7, -171.0, -201.0, -234.9, -273.0, -315.4, -362.4, -414.1, -470.2, -530.7, -595.0, -662.5, -732.1, -802.6, -872.7, -940.7],
        [11.6, 11.2, 10.7, 10.1, 9.4, 8.5, 7.5, 6.2, 4.8, 3.0, 1.0, -1.4, -4.3, -7.6, -11.6, -16.3, -21.9, -28.6, -36.5, -45.8, -56.8, -69.8, -84.8, -102.3, -122.5, -145.7, -172.2, -202.4, -236.5, -274.7, -317.4, -364.6, -416.4, -472.8, -533.5, -598.0, -665.5, -735.2, -805.8, -875.8, -943.7],
        [11.6, 11.2, 10.7, 10.1, 9.4, 8.5, 7.4, 6.2, 4.7, 3.0, 1.0, -1.5, -4.3, -7.7, -11.7, -16.5, -22.1, -28.8, -36.7, -46.1, -57.2, -70.2, -85.3, -102.9, -123.2, -146.5, -173.1, -203.4, -237.6, -276.0, -318.9, -366.3, -418.3, -474.8, -535.6, -600.2, -667.9, -737.7, -808.2, -878.2, -946.0],
        [11.6, 11.2, 10.7, 10.1, 9.3, 8.5, 7.4, 6.2, 4.7, 3.0, 0.9, -1.5, -4.4, -7.8, -11.8, -16.6, -22.2, -28.9, -36.9, -46.3, -57.4, -70.5, -85.7, -103.3, -123.7, -147.0, -173.8, -204.2, -238.5, -277.1, -320.0, -367.5, -419.7, -476.3, -537.2, -601.9, -669.7, -739.5, -810.1, -880.0, -947.8],
        [11.6, 11.2, 10.7, 10.1, 9.3, 8.5, 7.4, 6.2, 4.7, 3.0, 0.9, -1.6, -4.4, -7.9, -11.9, -16.6, -22.3, -29.0, -37.0, -46.5, -57.6, -70.7, -86.0, -103.6, -124.0, -147.5, -174.3, -204.8, -239.2, -277.8, -320.9, -368.5, -420.7, -477.4, -538.5, -603.2, -671.0, -740.9, -811.4, -881.3, -949.1],
        [11.6, 11.2, 10.7, 10.1, 9.3, 8.5, 7.4, 6.2, 4.7, 2.9, 0.9, -1.6, -4.5, -7.9, -11.9, -16.7, -22.4, -29.1, -37.1, -46.6, -57.8, -70.9, -86.2, -103.9, -124.3, -147.8, -174.7, -205.2, -239.7, -278.4, -321.5, -369.2, -421.5, -478.3, -539.4, -604.2, -672.1, -741.9, -812.5, -882.4, -950.0]
    ]

    za_penalty = [
        [503.8, 503.7, 503.7, 503.7, 503.6, 503.6, 503.5, 503.4, 503.3, 503.2, 503.1, 502.9, 502.7, 502.4, 502.1, 501.7, 501.2, 500.6, 499.8, 498.8, 497.5, 496.0, 494.0, 491.6, 488.5, 484.6, 479.8, 473.7, 466.3, 457.0, 445.6, 431.7, 414.7, 394.3, 370.0, 341.3, 307.9, 269.6, 226.4, 178.6, 126.8],
        [502.6, 502.5, 502.5, 502.5, 502.4, 502.3, 502.3, 502.2, 502.1, 501.9, 501.7, 501.5, 501.3, 500.9, 500.5, 500.0, 499.4, 498.6, 497.6, 496.3, 494.8, 492.8, 490.3, 487.2, 483.4, 478.5, 472.5, 465.0, 455.7, 444.4, 430.4, 413.5, 393.2, 368.9, 340.3, 306.9, 268.7, 225.5, 177.8, 126.1, 71.5],
        [501.1, 501.1, 501.0, 501.0, 500.9, 500.8, 500.7, 500.6, 500.5, 500.3, 500.1, 499.8, 499.5, 499.1, 498.6, 497.9, 497.1, 496.1, 494.9, 493.3, 491.3, 488.8, 485.7, 481.8, 477.0, 471.0, 463.5, 454.2, 442.9, 429.0, 412.1, 391.8, 367.6, 339.0, 305.7, 267.5, 224.5, 176.8, 125.2, 70.7, 14.6],
        [499.3, 499.3, 499.2, 499.2, 499.1, 499.0, 498.9, 498.7, 498.5, 498.3, 498.1, 497.7, 497.3, 496.8, 496.1, 495.3, 494.3, 493.0, 491.4, 489.5, 487.0, 483.9, 480.0, 475.1, 469.1, 461.6, 452.4, 441.0, 427.2, 410.3, 390.1, 365.9, 337.4, 304.2, 266.2, 223.2, 175.7, 124.2, 69.9, 13.8, -42.4],
        [497.1, 497.0, 497.0, 496.9, 496.8, 496.7, 496.5, 496.4, 496.1, 495.9, 495.5, 495.1, 494.6, 493.9, 493.1, 492.1, 490.8, 489.2, 487.2, 484.7, 481.6, 477.7, 472.9, 466.8, 459.4, 450.2, 438.9, 425.0, 408.2, 388.0, 364.0, 335.6, 302.5, 264.5, 221.7, 174.3, 123.0, 68.8, 12.9, -43.1, -97.9],
        [494.3, 494.3, 494.2, 494.1, 494.0, 493.8, 493.6, 493.4, 493.1, 492.8, 492.4, 491.8, 491.2, 490.3, 489.3, 488.0, 486.4, 484.4, 481.9, 478.8, 474.9, 470.1, 464.1, 456.6, 447.5, 436.2, 422.4, 405.7, 385.6, 361.6, 333.3, 300.4, 262.6, 219.9, 172.7, 121.6, 67.6, 11.9, -44.0, -98.7, -151.0],
        [490.9, 490.8, 490.7, 490.6, 490.4, 490.3, 490.0, 489.7, 489.4, 489.0, 488.4, 487.8, 486.9, 485.9, 484.6, 483.0, 481.0, 478.5, 475.4, 471.5, 466.7, 460.7, 453.3, 444.1, 432.9, 419.2, 402.6, 382.6, 358.8, 330.6, 297.9, 260.2, 217.8, 170.8, 119.9, 66.1, 10.7, -45.0, -99.5, -151.7, -200.7],
        [486.6, 486.5, 486.4, 486.2, 486.1, 485.8, 485.5, 485.2, 484.8, 484.2, 483.6, 482.7, 481.7, 480.4, 478.8, 476.8, 474.3, 471.2, 467.3, 462.5, 456.6, 449.2, 440.1, 429.0, 415.3, 398.8, 379.0, 355.3, 327.3, 294.8, 257.4, 215.2, 168.5, 117.9, 64.4, 9.2, -46.3, -100.5, -152.5, -201.3, -246.8],
        [481.3, 481.2, 481.0, 480.8, 480.6, 480.3, 480.0, 479.6, 479.0, 478.4, 477.5, 476.5, 475.2, 473.6, 471.6, 469.1, 466.0, 462.2, 457.4, 451.5, 444.2, 435.2, 424.1, 410.6, 394.2, 374.6, 351.1, 323.4, 291.1, 254.0, 212.1, 165.7, 115.5, 62.3, 7.4, -47.7, -101.7, -153.4, -202.1, -247.4, -289.3],
        [474.7, 474.6, 474.4, 474.2, 473.9, 473.5, 473.1, 472.6, 471.9, 471.1, 470.1, 468.8, 467.2, 465.2, 462.8, 459.7, 455.9, 451.2, 445.3, 438.1, 429.1, 418.2, 404.8, 388.7, 369.2, 346.0, 318.5, 286.6, 249.8, 208.3, 162.4, 112.5, 59.8, 5.3, -49.4, -103.1, -154.6, -203.0, -248.1, -289.9, -328.6],
        [466.6, 466.4, 466.2, 465.9, 465.6, 465.2, 464.6, 464.0, 463.2, 462.2, 460.9, 459.3, 457.4, 454.9, 451.9, 448.2, 443.5, 437.7, 430.6, 421.8, 411.0, 397.8, 381.9, 362.7, 339.7, 312.6, 281.1, 244.8, 203.8, 158.3, 109.0, 56.8, 2.8, -51.5, -104.8, -155.9, -204.1, -249.0, -290.6, -329.2, -364.8],
        [456.7, 456.5, 456.2, 455.8, 455.4, 454.9, 454.3, 453.5, 452.5, 451.2, 449.7, 447.8, 445.4, 442.4, 438.7, 434.1, 428.4, 421.4, 412.8, 402.2, 389.3, 373.6, 354.7, 332.2, 305.5, 274.4, 238.7, 198.3, 153.4, 104.7, 53.1, -0.3, -54.0, -106.8, -157.5, -205.3, -250.0, -291.4, -329.8, -365.4, -398.2],
        [444.6, 444.3, 444.0, 443.5, 443.0, 442.4, 441.6, 440.7, 439.4, 437.9, 436.1, 433.7, 430.8, 427.2, 422.7, 417.2, 410.3, 401.9, 391.5, 378.9, 363.5, 345.1, 323.0, 296.9, 266.4, 231.3, 191.6, 147.5, 99.6, 48.7, -3.9, -57.0, -109.2, -159.4, -206.9, -251.2, -292.4, -330.6, -366.0, -398.8, -429.1],
        [429.9, 429.5, 429.1, 428.6, 428.0, 427.3, 426.3, 425.2, 423.7, 421.9, 419.6, 416.8, 413.3, 408.9, 403.5, 396.9, 388.7, 378.6, 366.3, 351.4, 333.4, 311.9, 286.4, 256.7, 222.4, 183.6, 140.4, 93.4, 43.5, -8.3, -60.6, -112.1, -161.7, -208.7, -252.6, -293.5, -331.5, -366.7, -399.4, -429.6, -457.6],
        [412.2, 411.8, 411.3, 410.7, 410.0, 409.1, 407.9, 406.5, 404.8, 402.6, 399.8, 396.5, 392.3, 387.1, 380.6, 372.7, 363.0, 351.2, 336.7, 319.3, 298.5, 273.8, 245.0, 211.7, 173.9, 131.8, 85.9, 37.1, -13.6, -64.9, -115.6, -164.4, -210.8, -254.3, -294.9, -332.6, -367.6, -400.1, -430.2, -458.1, -483.9],
        [391.0, 390.6, 390.0, 389.3, 388.4, 387.3, 386.0, 384.3, 382.2, 379.6, 376.4, 372.4, 367.4, 361.3, 353.7, 344.4, 333.0, 319.2, 302.5, 282.5, 258.7, 230.9, 198.8, 162.2, 121.5, 77.0, 29.4, -20.0, -70.1, -119.7, -167.7, -213.4, -256.4, -296.5, -333.9, -368.6, -400.9, -430.9, -458.7, -484.4, -508.3],
        [366.0, 365.5, 364.8, 364.0, 363.0, 361.7, 360.1, 358.1, 355.7, 352.6, 348.8, 344.1, 338.3, 331.1, 322.3, 311.5, 298.3, 282.4, 263.4, 240.7, 214.1, 183.3, 148.3, 109.1, 66.2, 20.3, -27.6, -76.3, -124.7, -171.6, -216.5, -258.8, -298.4, -335.4, -369.9, -401.9, -431.7, -459.3, -485.0, -508.8, -531.1],
        [336.8, 336.1, 335.4, 334.4, 333.2, 331.7, 329.9, 327.5, 324.7, 321.1, 316.7, 311.3, 304.5, 296.3, 286.1, 273.8, 258.8, 240.8, 219.4, 194.3, 165.1, 131.8, 94.4, 53.4, 9.4, -36.7, -83.8, -130.6, -176.4, -220.2, -261.7, -300.7, -337.2, -371.3, -403.1, -432.7, -460.1, -485.6, -509.3, -531.5, -552.4],
        [302.9, 302.2, 301.3, 300.2, 298.8, 297.1, 295.0, 292.3, 289.0, 284.9, 279.9, 273.7, 266.0, 256.6, 245.1, 231.2, 214.5, 194.5, 171.0, 143.6, 112.3, 77.1, 38.3, -3.5, -47.5, -92.6, -137.7, -182.0, -224.6, -265.1, -303.4, -339.4, -373.0, -404.5, -433.8, -461.1, -486.4, -510.0, -532.1, -552.9, -572.7],
        [264.3, 263.4, 262.4, 261.2, 259.6, 257.7, 255.3, 252.3, 248.5, 243.9, 238.3, 231.3, 222.7, 212.2, 199.4, 184.0, 165.7, 144.0, 118.7, 89.6, 56.8, 20.5, -18.7, -60.3, -103.1, -146.2, -188.7, -229.9, -269.3, -306.7, -341.9, -375.1, -406.2, -435.2, -462.2, -487.3, -510.7, -532.7, -553.5, -573.1, -592.0],
        [220.9, 220.0, 218.8, 217.4, 215.7, 213.5, 210.9, 207.5, 203.4, 198.3, 192.0, 184.3, 174.8, 163.3, 149.5, 132.8, 113.1, 90.1, 63.6, 33.5, 0.0, -36.4, -75.1, -115.3, -156.1, -196.6, -236.1, -274.2, -310.6, -345.0, -377.6, -408.1, -436.8, -463.5, -488.4, -511.6, -533.5, -554.1, -573.7, -592.4, -610.5],
        [173.0, 172.0, 170.8, 169.3, 167.4, 165.0, 162.1, 158.4, 154.0, 148.4, 141.6, 133.2, 123.1, 110.8, 96.0, 78.5, 57.9, 34.1, 6.9, -23.5, -56.7, -92.4, -129.6, -167.7, -206.0, -243.6, -280.1, -315.2, -348.7, -380.5, -410.5, -438.7, -465.0, -489.6, -512.7, -534.4, -554.8, -574.3, -593.0, -611.0, -626.5],
        [121.4, 120.3, 119.0, 117.3, 115.3, 112.8, 109.6, 105.8, 101.0, 95.1, 87.9, 79.1, 68.4, 55.6, 40.3, 22.2, 1.2, -22.8, -49.9, -79.7, -112.0, -146.1, -181.2, -216.9, -252.3, -287.1, -320.7, -353.1, -384.0, -413.3, -440.9, -466.8, -491.1, -513.9, -535.4, -555.7, -575.1, -593.6, -611.5, -627.0, -643.9],
        [66.8, 65.7, 64.3, 62.6, 60.5, 57.8, 54.6, 50.5, 45.6, 39.5, 32.0, 23.0, 12.0, -1.0, -16.5, -34.6, -55.5, -79.1, -105.4, -134.1, -164.7, -196.7, -229.5, -262.5, -295.2, -327.3, -358.3, -388.2, -416.6, -443.6, -469.0, -492.9, -515.4, -536.6, -556.8, -576.0, -594.4, -612.2, -627.6, -644.4, -660.8],
        [10.6, 9.5, 8.0, 6.3, 4.1, 1.4, -1.9, -6.0, -11.0, -17.2, -24.7, -33.8, -44.7, -57.6, -72.8, -90.5, -110.8, -133.5, -158.5, -185.5, -214.2, -243.9, -274.3, -304.8, -334.9, -364.5, -393.1, -420.6, -446.8, -471.6, -495.0, -517.1, -538.1, -558.0, -577.0, -595.3, -613.0, -628.3, -645.0, -661.3, -677.2],
        [-45.9, -47.0, -48.4, -50.2, -52.3, -55.0, -58.3, -62.4, -67.3, -73.4, -80.8, -89.6, -100.3, -112.8, -127.5, -144.4, -163.6, -185.0, -208.5, -233.7, -260.3, -287.8, -315.8, -343.9, -371.7, -398.9, -425.3, -450.6, -474.7, -497.5, -519.2, -539.8, -559.4, -578.2, -596.3, -613.9, -629.1, -645.7, -661.9, -677.7, -693.1],
        [-101.3, -102.4, -103.8, -105.5, -107.6, -110.2, -113.4, -117.3, -122.1, -128.0, -135.0, -143.5, -153.6, -165.6, -179.5, -195.4, -213.4, -233.4, -255.2, -278.5, -303.0, -328.4, -354.3, -380.2, -405.8, -430.8, -455.1, -478.3, -500.5, -521.7, -541.8, -561.1, -579.7, -597.6, -614.9, -630.0, -646.5, -662.6, -678.3, -693.6, -708.6],
        [-154.5, -155.6, -156.9, -158.5, -160.5, -163.0, -166.1, -169.8, -174.4, -179.9, -186.6, -194.6, -204.1, -215.3, -228.3, -243.1, -259.8, -278.4, -298.5, -320.1, -342.7, -366.1, -390.0, -413.8, -437.4, -460.4, -482.6, -504.0, -524.6, -544.2, -563.1, -581.4, -599.0, -616.1, -631.0, -647.4, -663.4, -679.0, -694.2, -709.1, -723.5],
        [-204.9, -205.9, -207.1, -208.7, -210.6, -212.9, -215.8, -219.3, -223.6, -228.8, -235.0, -242.5, -251.4, -261.8, -273.8, -287.6, -303.1, -320.2, -338.8, -358.7, -379.6, -401.2, -423.1, -445.0, -466.6, -487.7, -508.2, -528.0, -547.1, -565.5, -583.4, -600.7, -617.6, -632.3, -648.5, -664.3, -679.8, -695.0, -709.7, -724.1, -738.0],
        [-252.2, -253.1, -254.3, -255.7, -257.5, -259.7, -262.4, -265.7, -269.7, -274.5, -280.3, -287.3, -295.5, -305.2, -316.3, -329.0, -343.3, -359.2, -376.3, -394.6, -413.9, -433.7, -453.7, -473.8, -493.7, -513.1, -532.0, -550.5, -568.3, -585.8, -602.7, -619.3, -633.8, -649.7, -665.4, -680.8, -695.8, -710.5, -724.7, -738.6, -751.8],
        [-296.4, -297.3, -298.4, -299.8, -301.4, -303.5, -306.0, -309.0, -312.7, -317.2, -322.6, -329.1, -336.7, -345.6, -356.0, -367.7, -380.9, -395.4, -411.2, -428.1, -445.7, -463.8, -482.2, -500.6, -518.9, -536.8, -554.4, -571.7, -588.6, -605.1, -621.4, -635.5, -651.2, -666.7, -681.9, -696.8, -711.3, -725.5, -739.2, -752.4, -765.0],
        [-337.8, -338.7, -339.7, -341.0, -342.5, -344.5, -346.8, -349.6, -353.1, -357.2, -362.2, -368.2, -375.2, -383.5, -393.0, -403.8, -415.9, -429.3, -443.7, -459.1, -475.2, -491.8, -508.6, -525.6, -542.4, -559.1, -575.6, -591.9, -607.9, -623.8, -637.5, -653.0, -668.2, -683.2, -697.9, -712.3, -726.3, -740.0, -753.1, -765.6, -777.4],
        [-376.7, -377.5, -378.5, -379.7, -381.1, -382.9, -385.0, -387.7, -390.9, -394.7, -399.3, -404.8, -411.3, -418.9, -427.6, -437.5, -448.5, -460.7, -473.9, -488.0, -502.7, -517.9, -533.3, -549.0, -564.6, -580.3, -595.8, -611.3, -626.6, -640.0, -655.0, -670.0, -684.7, -699.2, -713.4, -727.4, -740.9, -753.9, -766.3, -778.0, -787.7],
        [-413.2, -413.9, -414.8, -415.9, -417.3, -418.9, -420.9, -423.3, -426.3, -429.8, -434.0, -439.1, -445.0, -451.9, -459.8, -468.8, -478.9, -490.0, -502.0, -514.8, -528.3, -542.2, -556.5, -571.0, -585.7, -600.4, -615.2, -629.9, -642.8, -657.5, -672.1, -686.5, -700.8, -714.8, -728.5, -741.9, -754.8, -767.1, -778.7, -788.3, -798.4],
        [-447.4, -448.1, -448.9, -449.9, -451.2, -452.7, -454.5, -456.7, -459.4, -462.6, -466.5, -471.0, -476.4, -482.7, -489.9, -498.1, -507.2, -517.3, -528.3, -540.0, -552.3, -565.2, -578.4, -592.0, -605.8, -619.8, -633.9, -646.2, -660.4, -674.5, -688.6, -702.6, -716.4, -729.9, -743.1, -755.8, -768.0, -779.5, -789.1, -799.0, -807.9],
        [-479.4, -480.0, -480.8, -481.7, -482.8, -484.2, -485.9, -487.9, -490.3, -493.2, -496.7, -500.9, -505.8, -511.4, -518.0, -525.4, -533.7, -542.9, -552.9, -563.6, -575.0, -586.9, -599.3, -612.1, -625.2, -636.6, -650.1, -663.7, -677.4, -691.1, -704.7, -718.2, -731.5, -744.4, -757.0, -769.1, -780.5, -789.9, -799.8, -808.6, -815.4],
        [-509.2, -509.8, -510.5, -511.3, -512.4, -513.6, -515.2, -517.0, -519.2, -521.8, -525.0, -528.7, -533.2, -538.3, -544.2, -550.9, -558.5, -566.9, -576.1, -586.0, -596.5, -607.7, -619.3, -631.4, -642.0, -654.7, -667.7, -680.8, -694.0, -707.2, -720.4, -733.3, -746.1, -758.4, -770.3, -781.5, -790.9, -800.6, -809.4, -816.1, -822.6],
        [-537.1, -537.6, -538.3, -539.1, -540.0, -541.2, -542.5, -544.2, -546.2, -548.6, -551.5, -554.9, -558.8, -563.5, -568.9, -575.0, -581.9, -589.6, -598.1, -607.3, -617.1, -627.6, -638.6, -648.2, -660.1, -672.3, -684.8, -697.4, -710.2, -722.9, -735.5, -747.9, -760.0, -771.7, -782.8, -792.0, -801.6, -810.3, -816.9, -823.4, -827.8],
        [-563.2, -563.7, -564.3, -565.0, -565.9, -566.9, -568.2, -569.7, -571.5, -573.7, -576.3, -579.4, -583.0, -587.3, -592.2, -597.8, -604.1, -611.3, -619.1, -627.7, -636.9, -644.9, -655.4, -666.3, -677.7, -689.5, -701.5, -713.6, -725.8, -738.1, -750.1, -761.9, -773.4, -784.2, -793.2, -802.7, -811.3, -817.8, -824.2, -828.5, -831.7],
        [-587.7, -588.2, -588.8, -589.4, -590.2, -591.2, -592.4, -593.7, -595.4, -597.4, -599.8, -602.6, -605.9, -609.8, -614.3, -619.5, -625.4, -632.0, -639.3, -647.3, -654.1, -663.5, -673.5, -684.0, -694.9, -706.1, -717.6, -729.3, -741.0, -752.7, -764.1, -775.3, -785.9, -794.7, -804.0, -812.4, -818.8, -825.1, -829.4, -832.5, -834.4],
        [-610.9, -611.4, -611.9, -612.5, -613.3, -614.1, -615.2, -616.5, -618.0, -619.9, -622.1, -624.7, -627.7, -631.3, -635.4, -640.2, -645.7, -651.9, -656.8, -664.4, -672.6, -681.6, -691.1, -701.1, -711.5, -722.3, -733.3, -744.5, -755.7, -766.7, -777.5, -787.8, -796.4, -805.5, -813.7, -820.0, -826.2, -830.4, -833.4, -835.2, -835.8]
    ]

    axis = np.round(np.arange(-2, 2.001, 0.1), 1)
    CO_penalty = scp.interpolate.RegularGridInterpolator([axis, axis], CO_penalty, bounds_error = False, fill_value = np.inf)(np.array([CO_diff, carbon_diff]).T)
    za_penalty = scp.interpolate.RegularGridInterpolator([axis, axis], za_penalty, bounds_error = False, fill_value = np.inf)(np.array([zscale_diff, alpha_diff]).T)
    return np.sqrt(CO_penalty ** 2.0 + za_penalty ** 2.0)


def select_structure(teff, logg, zscale, eheu):
    """Choose the most appropriate structure from the structures grid for a given chemical composition

    The structure is chosen by minimizing the flux error penalties calculated with `structure_penalty()`.
    The function also returns the lowest penalty and the abundance offsets between the spectrum and the
    structure that may be passed to RESCALC for spectral synthesis
    
    Parameters
    ----------
    teff : float
        Desired effective temperature
    logg : float
        Desired surface gravity
    zscale : float
        Desired metallicity, [M/H]
    eheu : dict
        Dictionary of desired abundances, [A/M]
    
    Returns
    -------
    structure : str
        Unique identifier of the chosen structure
    penalty : float
        Associated flux error penalty in Kelvin
    offsets: dict
        Dictionary of abundance offsets
    """
    global livesyn_grid, livesyn_structures_index

    # Determine the requested alpha-enhancement, ["a"/"M"]. We assume ["a"/"M"]=[Mg/Fe]
    if 'Mg' not in eheu:
        eheu['Mg'] = 0
    alpha_required = eheu['Mg']

    # teff and logg must match the structure exactly
    if (teff not in livesyn_grid['teff']) or (logg not in livesyn_grid['logg']):
        raise ValueError('No structure available for (TEFF,LOGG)=({},{})'.format(teff, logg))

    # Determine the requested [C/H]-[O/H] and [C/H]
    if 'O' not in eheu:
        eheu['O'] = 0
    if 'C' not in eheu:
        eheu['C'] = 0
    CO_required = eheu['C'] - eheu['O']
    carbon_required = zscale + eheu['C']

    # Choose the best structure
    candidates = np.array([structure for structure in livesyn_structures_index if (livesyn_structures_index[structure]['teff'] == teff) and (livesyn_structures_index[structure]['logg'] == logg)])
    zscale_available, alpha_available, carbon_available = np.array([[livesyn_structures_index[structure]['zscale'], livesyn_structures_index[structure]['alpha'], livesyn_structures_index[structure]['carbon']] for structure in candidates]).T
    zscale_offsets = zscale - zscale_available
    alpha_offsets = alpha_required - alpha_available
    CO_offsets = CO_required - (carbon_available - alpha_available)
    carbon_offsets = carbon_required - (zscale_available + carbon_available)
    penalties = structure_penalty(zscale_offsets, alpha_offsets, CO_offsets, carbon_offsets)
    best = np.argmin(penalties)

    # Construct the required abundance shifts for this structure
    zscale_elements = ['C', 'O', 'Li', 'Be', 'B', 'N', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es']
    alpha_elements = ['Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti']
    offsets = {}
    for element in zscale_elements:
        model_value = zscale_available[best]
        if element in ['O'] + alpha_elements:
            model_value += alpha_available[best]
        if element == 'C':
            model_value += carbon_available[best]
        required_value = zscale
        if element in eheu:
            required_value += eheu[element]
        offsets[element] = np.round(required_value - model_value, 2)

        # Some assertions to make absolutely sure we used the right penalties
        if element == 'Am': # Stand in for metallicity
            assert offsets[element] == np.round(zscale_offsets[best], 2)
        if element == 'Mg': # Stand in for alpha
            assert np.round(offsets[element] - zscale_offsets[best], 2) == np.round(alpha_offsets[best], 2)
        if element == 'C':
            assert offsets[element] == np.round(carbon_offsets[best], 2)
        if element == 'O':
            carbon = 0.0
            if 'C' in offsets: carbon = offsets['C']
            assert np.round(carbon - offsets[element], 2) == np.round(CO_offsets[best], 2)

        if offsets[element] == 0:
            del offsets[element]

    return candidates[best], penalties[best], offsets


def synthesize(output, structure, linelist, eheu):
    global livesyn_structures_index

    # Prepare the model directory
    if os.path.isdir(output):
        return
    os.mkdir(output)

    # Extract the structure from the grid
    if structure not in livesyn_structures_index:
        raise ValueError('Unknown structure {}'.format(structure))
    z = zipfile.ZipFile(settings['structures_dir'] + '/' + structure[:structure.find('/')], 'r')
    structure_files = [filename for filename in z.namelist() if filename.startswith(structure[structure.find('/') + 1:]) and filename.endswith('.out')]
    for filename in structure_files:
        data = z.read(filename)
        f = open('{}/{}'.format(output, filename.split('/')[-1]), 'wb')
        f.write(data)
        f.close()
    z.close()

    # Run spectral synthesis
    rescalc.synthesize(output, eheu, linelist, output + '/model')

def generate_spectrum(model):
    global livesyn_workdir

    # Do not include O in live synthesis alpha
    alpha_elements = ['Ne', 'Mg', 'Si', 'S', 'Ar', 'Ca', 'Ti']

    # Select the structure and determine abundance offsets
    eheu = {}
    for element in alpha_elements:
        eheu[element] = model['alpha']
    for element in settings['elements']:
        if element not in eheu:
            eheu[element] = 0
        eheu[element] += model['live_{}'.format(element)]
    eheu = {element: np.round(eheu[element], 2) for element in eheu}
    structure, penalty, offsets = select_structure(model['teff'], model['logg'], model['zscale'], eheu)
    model_name = 't{teff:.2f}_g{logg:.2f}_z{zscale:.2f}_a{alpha:.2f}'.format(**model)
    for element in settings['elements']:
        model_name += '_{}{:.2f}'.format(element, model['live_{}'.format(element)])
    model_name = model_name.replace('-0.00', '0.00')

    # Run the synthesis
    if type(livesyn_workdir) is bool:
        init_livesyn()
    model_dir = '{}/{}'.format(livesyn_workdir, model_name)
    notify('Calculating spectrum for {}\nStructure: {} (penalty: {:.1f})'.format(model_name, structure, penalty))
    synthesize(model_dir, structure, settings['scratch'] + '/linelist/logg_{}'.format(model['logg']), offsets)

    # Parse results and bin them down
    wl, flux = np.loadtxt(model_dir + '/model/synthe_1/spectrum.asc', skiprows = 2, unpack = True, usecols = [0, 1])
    wl = wl[:(len(wl) // settings['binning']) * settings['binning']].reshape(-1, settings['binning']).mean(axis = 1)
    flux = flux[:(len(flux) // settings['binning']) * settings['binning']].reshape(-1, settings['binning']).mean(axis = 1)
    np.save(model_dir + '.npy', np.array([wl, flux]))
    shutil.rmtree(model_dir)

def generate_spectra(params):
    build_linelists()
    if len(params) == 1:
        return generate_spectrum(params[0])
    else:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = executor.map(generate_spectrum, params)

def read_grid_dimensions():
    global livesyn_grid

    grid = {param: livesyn_grid[param] for param in ['zscale', 'alpha', 'teff', 'logg']}
    for element in settings['elements']:
        grid['live_{}'.format(element)] = np.round(np.arange(settings['elements'][element][0], settings['elements'][element][1] + 0.0001, settings['elements'][element][2]), 2)
    return grid

def read_grid_model(params, grid, batch = False):
    global livesyn_workdir

    if type(livesyn_workdir) is bool:
        init_livesyn()

    # Parallel mode
    if batch:
        request = []
        for instance in params:
            model_name = 't{teff:.2f}_g{logg:.2f}_z{zscale:.2f}_a{alpha:.2f}'.format(**instance)
            for element in settings['elements']:
                model_name += '_{}{:.2f}'.format(element, instance['live_{}'.format(element)])
            model_name = model_name.replace('-0.00', '0.00')
            model = '{}/{}.npy'.format(livesyn_workdir, model_name)
            if not os.path.isfile(model):
                request += [instance]
        notify('Pre-computing {} new models'.format(len(request)))
        generate_spectra(request)
        return

    # Regular mode
    model_name = 't{teff:.2f}_g{logg:.2f}_z{zscale:.2f}_a{alpha:.2f}'.format(**params)
    for element in settings['elements']:
        model_name += '_{}{:.2f}'.format(element, params['live_{}'.format(element)])
    model_name = model_name.replace('-0.00', '0.00')
    model = '{}/{}.npy'.format(livesyn_workdir, model_name)
    if os.path.isfile(model):
        wl, flux = np.load(model)
    else:
        generate_spectra([params])
        wl, flux = np.load(model)

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

# Load the structure index
f = open(settings['structures_index'], 'rb')
livesyn_structures_index = pickle.load(f)
f.close()

# Extract the grid axes
livesyn_grid = {}
for structure in livesyn_structures_index:
    for param in livesyn_structures_index[structure]:
        if param not in livesyn_grid:
            livesyn_grid[param] = []
        if livesyn_structures_index[structure][param] not in livesyn_grid[param]:
            livesyn_grid[param] += [livesyn_structures_index[structure][param]]
for param in livesyn_grid:
    livesyn_grid[param] = np.unique(livesyn_grid[param])

# Add live synthesis elements to the list of degrees of freedom
settings['fit_dof'] = settings['fit_dof'] + ['live_{}'.format(element) for element in settings['elements']]

# Provide default initial guesses for live synthesis elements
settings['default_initial'] = {**settings['default_initial'], **{'live_{}'.format(element): np.round((settings['elements'][element][0] + settings['elements'][element][1]) / 2.0, 2) for element in settings['elements']}}
