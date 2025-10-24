############################################################
#                                                          #
#                     LOCALFIT PRESET                      #
#                                                          #
#  This preset is used to measure the abundances of        #
#  individual chemical elements using response functions   #
#  calculated at runtime. In order to use this preset, the #
#  stellar parameters (Teff, log(g), [M/H], [a/M], [C/M])  #
#  must be known already, e.g. from a prior grid search    #
#  fitting using the GRIDFIT preset                        #
#                                                          #
#  This preset requires BasicATLAS                         #
#                                                          #
############################################################


import atlas
import pickle
import scipy as scp
import hashlib
import os
import zipfile
import subprocess
import shutil
import copy


# This is the "defective" mask that removes parts of the spectrum where the models do not match the spectra of 
# Arcturus and the Sun well
defective_mask = [[3800, 4000], [4006, 4012], [4065, 4075], [4093, 4110], [4140, 4165], [4170, 4180], [4205, 4220],
                 [4285, 4300], [4335, 4345], [4375, 4387], [4700, 4715], [4775, 4790], [4855, 4865], [5055, 5065],
                 [5145, 5160], [5203, 5213], [5885, 5900], [6355, 6365], [6555, 6570], [7175, 7195], [7890, 7900],
                 [8320, 8330], [8490, 8505], [8530, 8555], [8650, 8672]]

settings = {
    # Stellar parameters (teff, logg, zscale, alpha, carbon) of the observed spectrum. The atmosphere structure will
    # be based on these parameters. The parameters must be specified by for each fit, which is done automatically
    # by `public__localfit()`
    'gridfit_params': False,

    # If True, a new set of ODFs and a new atmosphere structure will be calculated for the fit. Otherwise, the
    # precomputed structures grid is interpolated to the desired stellar parameters. This setting is configured
    # automatically by `public__localfit()` depending on the required level of analysis
    'compute_new_structure': True,

    # Minimum convergence class for newly calculated models that is considered acceptable. If a newly calculated
    # model atmosphere does not meet this criterion (has worse convergence class), an exception is raised
    'min_convergence': 'SILVER',

    # Offsets in the abundances of individual elements (from the chemistry of stellar parameters), at which the null
    # spectrum is calculated. This setting is configured automatically by `public__localfit()`
    'initial_abundance_offsets': {},

    # If True, the initial abundance offsets are incorporated in the atmosphere structure that will be computed for the
    # fit. Otherwise, only the chemistry of stellar parameters is considered. This flag is cannot be applied to
    # interpolated structures, so if `compute_new_structure` is set to False, and this flag is set to True, an exception
    # is raised. This setting is configured automatically by `public__localfit()`
    'use_initial_abundance_offsets_in_structure': False,

    # Which individual elements to fit, and the lines of which species (atomic, molecular and ionic) should be considered
    # when calculating the abundance of each element. All species are specified as Kurucz codes
    'elements': {
       # 'Ba': [56, 56.01, 56.02, 56.03, 56.04, 56.05],
       # 'Eu': [63, 63.01, 63.02, 63.03, 63.04, 63.05],
       'Mg': [12, 12.01, 12.02, 12.03, 12.04, 12.05, 112.0, 812.0, 112.01, 10812.0],
       # 'Si': [14, 14.01, 14.02, 14.03, 14.04, 14.05, 114.0, 814.0, 114.01, 80814.0, 101010114.0, 614.0, 60614.0],
       # 'Ca': [20, 20.01, 20.02, 20.03, 20.04, 20.05, 120.0, 820.0, 120.01, 10820.0],
       # 'Ti': [22, 22.01, 22.02, 22.03, 22.04, 22.05, 122.0], #, 822.0],
       # 'Ni': [28, 28.01, 28.02, 28.03, 28.04, 28.05, 128.0, 828.0],
       # 'Na': [11, 11.01, 11.02, 11.03, 11.04, 11.05, 111.0, 10811.0],
       'Fe': [26, 26.01, 26.02, 26.03, 26.04, 26.05, 126.0, 826.0],
       # 'Mn': [25, 25.01, 25.02, 25.03, 25.04, 25.05, 125.0, 825.0],
       # 'Co': [27, 27.01, 27.02, 27.03, 27.04, 27.05, 127.0, 827.0],
       # 'Li': [3, 3.01, 3.02, 3.03, 3.04, 3.05, 103.0],
       # 'Al': [13, 13.01, 13.02, 13.03, 13.04, 13.05, 113.0, 813.0, 113.01],
       # 'K': [19, 19.01, 19.02, 19.03, 19.04, 19.05, 119.0],
       # 'Cr': [24, 24.01, 24.02, 24.03, 24.04, 24.05, 124.0, 824.0],
       # 'La': [57, 57.01, 57.02, 57.03, 57.04, 57.05, 857.0],
       # 'Sc': [21, 21.01, 21.02, 21.03, 21.04, 21.05, 121.0, 821.0],
       # 'V': [23, 23.01, 23.02, 23.03, 23.04, 23.05, 123.0, 823.0],
       # 'Y': [39, 39.01, 39.02, 39.03, 39.04, 39.05, 839.0],
    },

    # Minimum change in normalized spectrum between the lowest and highest abundances of an element in order to be included
    # in that element's fitting mask
    'threshold': 0.01,

    # Some elements have significant effect not just on their lines, but also on the lines of other elements through their
    # impact on chemical equilibrium. For these elements, we recompute all lines when the abundances are varied
    'higher_order_impact': ['Mg', 'Si', 'Na', 'Ca', 'Al'],

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

    # Scratch directory to store newly calculated structures, spectra, response functions and linelists. Must be specified in the
    # local settings preset
    'scratch': original_settings['scratch'],

    # Apply the defective mask to the fits
    'masks': apply_standard_mask(defective_mask, original_settings['masks']),
}

# Define the ranges of "virtual" parameters (i.e. parameters that are not dimensions of the model grid). For this preset, we do not use a model
# grid at all, so all parameters are virtual. They include redshift and the abundances of all individual elements defined above
settings['virtual_dof'] = {'redshift': [-300, 300], **{element: [np.min(settings['abun']), np.max(settings['abun'])] for element in settings['elements']}}

# Which degrees of freedom (parameters) to consider in the fit? By default, we include all of them
settings['fit_dof'] = ['redshift'] + [element for element in settings['elements']]

# Default initial guesses for the parameters. We start at 0 for all of them
settings['default_initial'] = {'redshift': 0.0, **{element: 0.0 for element in settings['elements']}}

# Relationship between C12/C13 isotope ratio and surface gravity from Kirby+2015
C12C13_kirby_2015 = lambda logg: np.where(logg > 2.7, 50, np.where(logg <= 2.0, 6, 63 * logg - 120))

# Relationship between VTURB and LOGG based on Deimos spectra from Gerasimov+2026
VTURB_LOGG = lambda logg: 2.792 * np.exp(-0.241 * logg -0.118)

# Global variable storage used by the functions
_global = {}

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

def build_linelist(species, output_dir, C12C13, air_wl, wl_start, wl_end, res, vturb, atoms, invert = False):
    """Build a SYNTHE linelist (fort.12) for a given wavelength range and list of species to include

    This function first dispatches the SYNBEG code from the SYNTHE suite that prepares the header of the calculation (fort.93). The
    header stores the synthesis parameters, including turbulent velocity. Even though the chosen turbulent velocity does not
    affect the linelist itself, it must be specified at this stage for SYNBEG to work. We then dispatch all of the codes that convert
    the available line databases into the SYNTHE format, including RGFALLLINESNEW (atomic lines), RMOLECASC (diatomics), RH2OFAST (water).
    The isotope ratio of carbon-12 to carbon-13 must also be specified at this stage, as it is treated as a correction to the oscillator
    strength by RMOLECASC
    
    Parameters
    ----------
    species : list
        Which species to include the lines of? Must be a list of Kurucz codes
    output_dir : str
        Directory to store the output. Must not exist beforehand
    C12C13 : number
        Carbon-12 to carbon-13 isotope ratio
    air_wl : bool
        If `True`, use air wavelengths. Otherwise, vacuum wavelengths
    wl_start : number
        Start wavelength in nm
    wl_end : number
        End wavelength in nm
    res : number
        Native resolution of the spectrum, R = lambda / delta-lambda
    vturb : number
        Turbulent velocity in km/s
    atoms : str
        Atomic line list to use. 'BasicATLAS' is recommended
    invert : bool, optional
        If `True`, include all species except for the ones specified in `species`. Otherwise (default), include only the species specified in
        `species`
    """
    # Helper function to convert Kurucz species codes into NELION index. The NELION index is how they are stored in the linelist
    def code_to_nelion(code):
        code_mol = [101., 106., 107., 108., 606., 607., 608., 707., 708., 808., 112., 113., 114., 812., 813., 814., 116., 120., 816., 820., 821., 822., 823., 103., 104., 105., 109., 115., 117., 121., 122., 123., 124., 125., 126.,106.01, 107.01,108.01,112.01,113.01,114.01,120.01, 111., 119.,10101.01, 817., 824., 825., 826.,10108.,60808.,10106.,60606., 127., 128., 129., 827., 828., 829.,608.01, 408., 508., 815., 10808.,10811.,10812.,10820.,10106.,10107.,10116.,10606.,10607., 10608.,10708.,60816.,61616.,70708.,70808.,80814.,80816.,1010106.,1010107.,1010606.,101010106.,101010114., 614.,60614.,60607., 6060707.,6060607., 839., 840., 857.]
        nelion_mol = [240, 246, 252, 258, 264, 270, 276, 282, 288, 294, 300, 306, 312, 318, 324, 330, 336, 342, 348, 354, 360, 366, 372, 378, 384, 390, 396, 402, 408, 414, 420, 426, 432, 438, 444, 450, 456, 462, 468, 474, 480, 486, 492, 498, 504, 510, 516, 522, 528, 534, 540, 546, 552, 558, 564, 570, 576, 582, 588, 594, 600, 606, 612, 618, 624, 630, 636, 642, 648, 654, 660, 666, 672, 678, 684, 690, 696, 702, 708, 714, 720, 726, 732, 738, 744, 750, 756, 762, 768, 774, 780, 786, 792]

        if code > 100:
            return nelion_mol[np.where(np.array(code_mol) == code)[0][0]]

        Z = int(np.floor(code))
        charge = int(np.round((code - Z) * 100.0))
        nelion = Z * 6 - 6 + (charge + 1)
        if Z > 19 and Z < 29 and charge > 5:
            nelion = 6 * (Z + charge * 10 - 30) - 1
        return nelion

    # Raise an exception if the output directory already exists, create it otherwise
    if os.path.isdir(output_dir):
        raise ValueError('{} already exists'.format(output_dir))
    os.mkdir(output_dir)

    # Get the part of the SYNTHE launcher script in BasicATLAS that prepares the line list
    script = atlas.templates.synthe_control
    start = script.find('# synberg.exe initializes the computation')
    end = script.find('# synthe.exe computes line opacities')
    if start == -1 or end == -1:
        raise ValueError('LOCALFIT no longer compatible with BasicATLAS')
    script = 'cd {output}\n' + script[start:end]

    # Run the extracted part of the launcher script
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
        'turbv': vturb,
        'linelist': atoms_linelist,
        'C12C13': C12C13_cmd,
        'ifnlte': 0,
        'linout': -1,
        'cutoff': 0.0001,
        'ifpred': 1,
        'nread': 0,
        'output': os.path.realpath(output_dir),
    }
    f = open(output_dir + '/packager.com', 'w')
    f.write(script.format(**cards))
    f.close()
    command = 'bash {}/packager.com'.format(output_dir)
    session = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    stdout, stderr = session.communicate()
    stdout = stdout.decode().strip()
    stderr = stderr.decode().strip()
    if stderr != '':
        raise ValueError('Command {} returned an error: {}'.format(command, stderr))
    if stdout != '':
        print(stdout, file = sys.stderr)

    if len(species) > 0 or (not invert):
        # Filter lines to include the species we want
        f = open(output_dir + '/fort.12', 'rb')
        dt = np.dtype('i4,i4,f4,i4,f4,f4,f4,f4,f4,i4')
        linelist = np.fromfile(f, dtype = dt, count = -1)
        f.close()
        mask = np.full(len(linelist), False)
        for element in species:
            mask |= linelist['f3'] == code_to_nelion(element)
        if invert:
            mask = ~mask
        f = open(output_dir + '/fort.12', 'wb')
        linelist[mask].tofile(f)
        f.close()

        # Update the total number of lines in fort.93
        f = open(output_dir + '/fort.93', 'rb')
        dt = np.dtype('i4,i4')
        headers = np.fromfile(f, dtype = dt, count = -1)
        f.close()
        headers[0][1] = np.count_nonzero(mask)
        f = open(output_dir + '/fort.93', 'wb')
        headers.tofile(f)
        f.close()

    # Remove unnecessary output files
    for filename in os.listdir(output_dir):
        if not filename.startswith('fort.'):
            os.remove('{}/{}'.format(output_dir, filename))

def build_all_linelists(output_dir, elements, higher_order_impact, C12C13, air_wl, wl_start, wl_end, res, vturb, atoms):
    """Build all linelists required for calculating response functions
    
    For a set of response functions, we require the null spectrum linelist which includes all species, and a line list
    for each considered element that only includes the lines of that element (and its ions/molecules)

    For the linelists of selected elements listed in `higher_order_impact`, we include not just the lines of those elements,
    but the lines of all considered elements. This is because the elements listed in `higher_order_impact` may have significant
    impact on all lines in the spectrum through offsets in the chemical equilibrium (as well as the corresponding offsets in
    the continuum opacity)
    
    Parameters
    ----------
    output_dir : str
        Directory to store the output
    elements : list
        List of elements to compute linelists for
    higher_order_impact : list
        List of elements whose linelists must include not just their lines, but also the lines of all other considered elements
    C12C13 : number
        Carbon-12 to carbon-13 isotope ratio. See `build_linelist()`
    air_wl : bool
        If `True`, use air wavelengths. Otherwise, vacuum wavelengths
    wl_start : number
        Start wavelength in nm
    wl_end : number
        End wavelength in nm
    res : number
        Native resolution of the spectrum, R = lambda / delta-lambda
    vturb : number
        Turbulent velocity in km/s. See `build_linelist()`
    atoms : str
        Atomic line list to use. 'BasicATLAS' is recommended
    """
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    for element in elements:
        if element not in higher_order_impact:
            species = elements[element]
        else:
            species = sum(list(elements.values()), [])
        build_linelist(species, output_dir + '/{}'.format(element), C12C13, air_wl, wl_start, wl_end, res, vturb, atoms)

    # Build the null (all-inclusive) line list as well
    build_linelist([], output_dir + '/null', C12C13, air_wl, wl_start, wl_end, res, vturb, atoms, invert = True)

def update_abundances(filename, abun_offsets):
    """Update chemical abundances in the output_summary.out file, which BasicATLAS uses to specify the abundance vector
    used in spectral synthesis
    
    Parameters
    ----------
    filename : str
        File name of the output_summary.out file to be updated
    abun_offsets : dict
        Dictionary of required abundance offsets
    """
    # Load the existing abundances and convert them to the standard format
    elements, params = atlas.parse_atlas_abundances(filename, lookbehind = 1, params = ['ABUNDANCE SCALE'])
    abun = atlas.Settings().abun_atlas_to_std(elements, np.log10(params['ABUNDANCE SCALE']))

    # Apply the required abundance offsets
    for element in abun_offsets:
        if element not in abun['abun']:
            abun['abun'][element] = 0
        abun['abun'][element] += abun_offsets[element]

    # Use the templates in BasicATLAS to generate a new output_summary.out file with the offset abundances and
    # save it
    elements = atlas.Settings().abun_std_to_atlas(**abun)
    template = atlas.templates.atlas_control
    template = template[template.find('ABUNDANCE SCALE'):template.find('\n', template.find('ABUNDANCE CHANGE 99'))]
    sub = {'element_{}'.format(i): elements[i] for i in range(1, 100)}
    template = template.format(abundance_scale = 10 ** abun['zscale'], **sub)
    f = open(filename, 'r')
    content = f.read()
    f.close()
    start = content.find('ABUNDANCE SCALE')
    end = content.find('\n', content.find('ABUNDANCE CHANGE 99'))
    content = content[:start] + template + content[end:]
    f = open(filename, 'w')
    f.write(content)
    f.close()

def get_opacity_header_size(filename):
    """Get the size of the header in a SYNTHE opacity file (fort.9)
    
    SYNTHE opacity maps contain a header and a data table. When doing arithmetic operations on opacity maps, we only
    change the data table and not the header, so we need to know how long the header is. This function reads a fort.9
    file in 4-byte chunks until it runs into the end signature of the header
    
    Parameters
    ----------
    filename : str
        File name of the opacity map
    
    Returns
    -------
    number
        Size of the header in the opacity file in 4-byte chunks
    """
    size = 1
    f = open(filename, 'rb')
    dt = np.dtype('i4')
    data = np.fromfile(f, dtype = dt, count = 1)
    for i in range(2):
        while data[0] != 9860:
            size += 1
            data = np.fromfile(f, dtype = dt, count = 1)
        size += 1
        data = np.fromfile(f, dtype = dt, count = 1)
    f.close()
    return size

def load_opacity_table(filename):
    """Load a SYNTHE opacity file (fort.9)
    
    The contents of the file are loaded in a format on which arithmetic operations can be performed
    
    Parameters
    ----------
    filename : str
        File name of the opacity map
    
    Returns
    -------
    header : array_like
        File header that precedes the opacity map
    data : array_like
        Opacity map
    footer: array_like
        File footer that follows the opacity map
    """
    size = get_opacity_header_size(filename)
    f = open(filename, 'rb')
    dt = np.dtype('i4')
    header = np.fromfile(f, dtype = dt, count = size)
    dt = np.dtype(','.join(['f4'] * 99 + ['i4','i4']))
    data = np.fromfile(f, dtype = dt, count = -1)
    dt = np.dtype('i4')
    footer = np.fromfile(f, dtype = dt, count = -1)
    f.close()
    return header, data, footer

def save_opacity_table(filename, header, data, footer):
    """Save a SYNTHE opacity file (fort.9)
    
    Parameters
    ----------
    filename : str
        File name of the output file
    header : array_like
        File header that precedes the opacity map in the format returned by `load_opacity_table()`
    data : array_like
        Opacity map  in the format returned by `load_opacity_table()`
    footer: array_like
        File footer that follows the opacity map  in the format returned by `load_opacity_table()`
    """
    f = open(filename, 'wb')
    header.tofile(f)
    data.tofile(f)
    footer.tofile(f)
    f.close()

def synthesize(structure, abun_offsets, linelist, output_dir, update_opacity = [], update_continuum = False, run_synthesis = True):
    """Compute a synthetic spectrum
    
    This function computes a synthetic spectrum using SYNTHE (opacity calculator) and SPECTRV (radiative transfer) for an ATLAS model
    in the BasicATLAS format. Unlike `atlas.synthe()` in BasicATLAS, this implementation is more general, and allows for custom line
    lists and opacity tables.

    We begin with an ATLAS atmosphere structure (`structure`) and a linelist previously calculated with `build_linelist()` (`linelist`).
    The linelist must also contain the header (fort.93) which stores the synthesis parameters (wavelength range, resolution etc). The
    ATLAS model is then copied into a new directory and the abundance offsets (`abun_offsets`) are applied if the synthetic spectrum
    needs to have a different chemistry compared to the structure. Note that we use XNFPELSYN to re-evaluate the chemical equilibrium,
    so the effect of abundance offsets on the molecular number densities and ionization balance will be accounted for.

    The function will then dispatch XNFPELSYN and SYNTHE to compute the line and continuum opacities (fort.9 and fort.10). At this point,
    the user has the option of replacing the continuum opacity with their own fort.10 file (this is useful if the continuum opacity is
    expected to be very similar to that of another model, so we can ignore the difference and get more numerically consistent results).
    The user can also update the line opacity by adding or subtracting other fort.9 opacity tables.

    Finally, we run spectral synthesis with SPECTRV and convert the result to ASCII with CONVERFSYNNMTOA. Spectral synthesis can be
    skipped if the user only needs opacity maps
    
    Parameters
    ----------
    structure : str
        File name of the directory that contains the ATLAS atmosphere structure in the BasicATLAS format. Since we may be updating
        the abundances of the model, it will be copied to the output directory first, such that the original model is not changed
    abun_offsets : dict
        Offsets in element abundances, see `update_abundances()`
    linelist : str
        File name of the linelist directory produced by `build_linelist()`
    output_dir : str
        Directory to store the output. Must not exist beforehand
    update_opacity : list, optional
        List of other line opacity maps (fort.9) to add or subtract from the opacity of the model. Each element of the list is a
        2-tuple with the first element storing the file name of the opacity map, and the second element storing the coefficient,
        such that `+1` represents addition, and `-1` represents subtraction
    update_continuum : str, optional
        File path of the continuum opacity file (fort.10) to be used instead of the calculated continuum opacity. If `False`, the
        calculated continuum opacity will be used
    run_synthesis : bool, optional
        If `False`, do not run spectral synthesis and exit once the opacity maps are calculated
    """
    # Check that the output directory does not already exist and create it
    if os.path.isdir(output_dir):
        raise ValueError('{} already exists'.format(output_dir))
    os.mkdir(output_dir)

    # Copy the ATLAS model into the output, so we can safely update its abundances
    shutil.copy(structure + '/output_summary.out', output_dir + '/')

    # Update abundances
    update_abundances(output_dir + '/output_summary.out', abun_offsets)

    # Convert the ATLAS model into SYNTHE input
    f = open(output_dir + '/output_summary.out', 'r')
    model = f.read()
    f.close()
    # Remove ATLAS turbulent velocity from the output
    model = model.split('\n')
    in_output = False
    for i, line in enumerate(model):
        if line.find('FLXRAD,VCONV,VELSND') != -1:
            in_output = True
            continue
        if line.find('PRADK') != -1:
            in_output = False
            continue
        if in_output:
            model[i] = line[:-40] + ' {:9.3E}'.format(0) + line[-30:]
    model = '\n'.join(model)
    f = open(output_dir + '/output_synthe.out', 'w')
    f.write(atlas.templates.synthe_prependix + model)
    f.close()

    # Retrieve the SYNTHE launcher script, except replace the SYNBEG calls and linelist creation with the already computed linelist
    script = atlas.templates.synthe_control
    start = script.find('# synberg.exe initializes the computation')
    mid = script.find('# synthe.exe computes line opacities')
    end = script.find('# spectrv.exe computes the synthetic spectrum')
    if start == -1 or end == -1 or mid == -1:
        raise ValueError('rescalc no longer compatible with BasicATLAS')
    script = script[:start] + '\nln -s {linelist}/fort.12 ./\nln -s {linelist}/fort.19 ./\ncp {linelist}/fort.14 ./\ncp {linelist}/fort.20 ./\ncp {linelist}/fort.93 ./\n' + script[mid:end]

    # Run SYNTHE
    cards = {
      's_files': atlas.python_path + '/data/synthe_files/',
      'd_files': atlas.python_path + '/data/dfsynthe_files/',
      'synthe_suite': atlas.python_path + '/bin/',
      'synthe_solar': os.path.realpath(output_dir + '/output_synthe.out'),
      'output_dir': os.path.realpath(output_dir),
      'synthe_num': 1,
      'linelist': os.path.realpath(linelist),
    }
    f = open(output_dir + '/launcher.com', 'w')
    f.write(script.format(**cards))
    f.close()
    command = 'bash {}/launcher.com'.format(output_dir)
    session = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    stdout, stderr = session.communicate()
    stdout = stdout.decode().strip()
    stderr = stderr.decode().strip()
    if stderr != '':
        raise ValueError('Command {} returned an error: {}'.format(command, stderr))
    if stdout != '':
        print(stdout, file = sys.stderr)

    # Save the opacity table for future use
    os.rename(original := os.path.realpath(output_dir + '/synthe_1/fort.9'), new := os.path.realpath(output_dir + '/opacity.9'))
    os.symlink(new, original)
    # Save the continuum
    os.rename(original_cont := os.path.realpath(output_dir + '/synthe_1/fort.10'), new_cont := os.path.realpath(output_dir + '/continuum.10'))
    os.symlink(new_cont, original_cont)

    # If spectral synthesis is not required, we can clean up the SYNTHE run directory and quit at this point
    if not run_synthesis:
        shutil.rmtree(output_dir + '/synthe_1')
        return

    # Update the opacity table if necessary
    if len(update_opacity) > 0:
        header, data, footer = load_opacity_table(new)
        for opacity in update_opacity:
            update = load_opacity_table(opacity[0])[1]
            for i in range(99):
                data['f{}'.format(i)] += opacity[1] * update['f{}'.format(i)]
        os.remove(original)
        save_opacity_table(original, header, data, footer)

    # Update continuum if necessary
    if type(update_continuum) != bool:
        os.remove(original_cont)
        shutil.copy(update_continuum, original_cont)

    # Run SPECTRV and CONVERFSYNNMTOA
    script = atlas.templates.synthe_control
    script = 'cd {output_dir}\ncd synthe_{synthe_num}/\n\n' + script[end:]
    f = open(output_dir + '/launcher.com', 'w')
    f.write(script.format(**cards))
    f.close()
    command = 'bash {}/launcher.com'.format(output_dir)
    session = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    stdout, stderr = session.communicate()
    stdout = stdout.decode().strip()
    stderr = stderr.decode().strip()
    if stderr != '':
        raise ValueError('Command {} returned an error: {}'.format(command, stderr))
    if stdout != '':
        print(stdout, file = sys.stderr)

def combine_abundances(abun_1, abun_2):
    """Helper function to add two abundance dictionaries
    
    Parameters
    ----------
    abun_1 : dict
        First abundance dictionary
    abun_2 : dict
        Second abundance dictionary
    
    Returns
    -------
    dict
        Combined abundance dictionary
    """
    abun = copy.deepcopy(abun_1)
    for element in abun_2:
        abun[element] = np.round(abun[element] + abun_2[element], 2)
    return abun

def clean_response(response_dir):
    """Helper function to remove all output of `synthesize()` from disk except for the synthetic spectrum itself
    
    This function is useful when the opacity tables are no longer needed and can be safely removed to conserve
    disk space
    
    Parameters
    ----------
    response_dir : str
        Output directory of `synthesize()`
    """
    for filename in os.listdir(response_dir):
        if filename != 'synthe_1':
            os.remove('{}/{}'.format(response_dir, filename))
    for filename in os.listdir('{}/synthe_1'.format(response_dir)):
        if filename != 'spectrum.asc':
            os.remove('{}/synthe_1/{}'.format(response_dir, filename))

def read_grid_dimensions():
    """Carry out the necessary setup for model fitting using response functions calculated at runtime
    
    This function normally communicates to chemfit which spectral models are available in the model grid; however,
    in the LOCALFIT preset we do not use a model grid. Instead, we use this function to prepare everything we need
    for response function fitting which includes the following:

        (1) Atmosphere structure model (newly computed with DFSYNTHE/ATLAS or interpolated)
        (2) Linelists
        (3) Null spectrum and its opacity map
        (4) Contributions of individual elements to the null opacity map
        (5) Element masks
        (6) The `_global` dictionary that will store all of the response functions computed at runtime as well as
            the data necessary to carry out those calculations
    
    Returns
    -------
    dict
        Stellar parameters, for which the response functions will be calculated
    """
    global _global

    if type(settings['gridfit_params']) is bool:
        raise ValueError('Provide stellar parameters in settings[\'gridfit_params\']')

    # Get the header of the structures grid and make sure it is complete
    header = atlas.restarts.load_header()
    if os.path.basename(header['grid']) != 'light.h5':
        raise ValueError('BasicATLAS must be configured to use the full grid of restart models')

    # Generate the ATLAS settings object for the structure model
    ATLAS_settings = atlas.Settings()
    ATLAS_settings.teff = np.round(settings['gridfit_params']['teff'], 0)
    ATLAS_settings.logg = np.round(settings['gridfit_params']['logg'], 3)
    ATLAS_settings.zscale = np.round(settings['gridfit_params']['zscale'], 3)
    ATLAS_settings.Y = 0.245
    ATLAS_settings.abun = {element: np.round(settings['gridfit_params']['alpha'], 2) for element in settings['alpha_elements']}
    ATLAS_settings.abun['C'] = np.round(settings['gridfit_params']['carbon'] + header['carbon_map']([settings['gridfit_params']['zscale'], settings['gridfit_params']['logg']])[0], 2)
    if settings['use_initial_abundance_offsets_in_structure']:
        for element in settings['initial_abundance_offsets']:
            if element not in ATLAS_settings.abun:
                ATLAS_settings.abun[element] = 0.0
            ATLAS_settings.abun[element] = np.round(ATLAS_settings.abun[element] + settings['initial_abundance_offsets'][element], 2)

    # Define a unique name for this structure
    structure_name = '{}_{}'.format(hashlib.md5(pickle.dumps([ATLAS_settings.teff, ATLAS_settings.logg, ATLAS_settings.zscale, ATLAS_settings.abun])).hexdigest(), ['interpolated', 'original'][settings['compute_new_structure']])
    notify('Will use structure {}'.format(structure_name))

    # Create a directory for structures if it does not already exist
    if not os.path.isdir(structures_dir := ('{}/localfit_structures'.format(settings['scratch']))):
        os.mkdir(structures_dir)

    # If the structure already exists, do nothing
    structure_dir = '{}/{}'.format(structures_dir, structure_name)
    _global['structure'] = structure_dir
    if os.path.isdir(structure_dir):
        notify('Structure {} already exists'.format(structure_name), color = 'g')

    # Compute new structure
    elif settings['compute_new_structure']:
        notify('Calculating ODF for {}'.format(structure_name), color = 'y')
        ODF_dir = '{}_ODF'.format(structure_dir)
        atlas.dfsynthe(ODF_dir, settings = ATLAS_settings, parallel = True, silent = True)
        notify('ODF calculated!', color = 'g')
        notify('Calculating model atmosphere for {}'.format(structure_name), color = 'y')
        atlas.atlas(structure_dir, settings = ATLAS_settings, ODF = ODF_dir, silent = True)
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
        notify('Model atmosphere calculated! Convergence: {} '.format(conv_classes[conv]), color = 'g')
        if conv > dict((v,k) for k,v in conv_classes.items())[settings['min_convergence']]:
            raise ValueError('Stopping because the model convergence class {} is worse than the minimum acceptable class in the settings ({})'.format(conv_classes[conv], settings['min_convergence']))

    # Interpolate the structure from the grid
    else:
        if settings['use_initial_abundance_offsets_in_structure']:
            raise ValueError('Cannot use the use_initial_abundance_offsets_in_structure flag with an interpolated structure')
        notify('Building an interpolated structure', color = 'y')
        params = {'teff': ATLAS_settings.teff, 'logg': ATLAS_settings.logg, 'zscale': ATLAS_settings.zscale, 'alpha': ATLAS_settings.abun['Mg'], 'carbon': np.round(settings['gridfit_params']['carbon'], 2)}
        structure = atlas.restarts.interpolate_structure(params, header = header)
        model = atlas.restarts.generate_model(*structure)
        os.mkdir(structure_dir)
        f = open('{}/output_summary.out'.format(structure_dir), 'w'); f.write(model); f.close()
        f = open('{}/output_last_iteration.out'.format(structure_dir), 'w'); f.close()
        f = open('{}/output_main.out'.format(structure_dir), 'w'); f.close()
        notify('Saved interpolated structure in {}'.format(structure_dir), color = 'g')

    # Prepare the linelist
    if not os.path.isdir(linelists_dir := ('{}/localfit_linelists'.format(settings['scratch']))):
        os.mkdir(linelists_dir)
    C12C13 = C12C13_kirby_2015(ATLAS_settings.logg)
    vturb = VTURB_LOGG(ATLAS_settings.logg)
    args = [settings['elements'], settings['higher_order_impact'], C12C13, settings['air_wl'], settings['wl_start'], settings['wl_end'], settings['res'], vturb, settings['atoms']]
    linelist_name = hashlib.md5(pickle.dumps(args)).hexdigest()
    notify('Will use linelist {}'.format(linelist_name))
    linelist_dir = '{}/{}'.format(linelists_dir, linelist_name)
    _global['linelist'] = linelist_dir
    if os.path.isdir(linelist_dir):
        notify('Linelist {} already exists'.format(linelist_name), color = 'g')
    else:
        notify('Building linelist {}'.format(linelist_name), color = 'y')
        build_all_linelists(linelist_dir, *args)
        notify('Finished building linelist {}'.format(linelist_name), color = 'g')

    # Determine the sampling of abundance offsets
    abun_offsets = {element: list(settings['abun']) for element in settings['elements']}
    for element in settings['initial_abundance_offsets']:
        abun_offsets[element] += [settings['initial_abundance_offsets'][element]]
        abun_offsets[element] = np.unique(np.round(np.array(abun_offsets[element]) - settings['initial_abundance_offsets'][element], 2)).tolist()
    null_offsets = {element: 0.0 for element in settings['elements']}
    if not settings['use_initial_abundance_offsets_in_structure']:
        for element in settings['initial_abundance_offsets']:
            null_offsets[element] = np.round(null_offsets[element] + settings['initial_abundance_offsets'][element], 2)

    # Set up the work directory
    if not os.path.isdir(synthesis_dir := ('{}/localfit_synth'.format(settings['scratch']))):
        os.mkdir(synthesis_dir)
    synth_name = hashlib.md5(pickle.dumps([structure_name, linelist_name, abun_offsets, null_offsets])).hexdigest()
    workdir = '{}/{}'.format(synthesis_dir, synth_name)
    if not os.path.isdir(workdir):
        os.mkdir(workdir)
    _global['workdir'] = workdir
    notify('Will use response functions {}'.format(synth_name))

    # Compute null spectrum
    if not os.path.isdir(workdir + '/null'):
        notify('Calculating the null spectrum', color = 'y')
        synthesize(structure_dir, null_offsets, linelist_dir + '/null', workdir + '/null')

    # Compute null opacities for each element
    for element in settings['elements']:
        if not os.path.isdir('{}/{}_null'.format(workdir, element)):
            notify('Calculating the null opacity for {}'.format(element), color = 'y')
            synthesize(structure_dir, null_offsets, '{}/{}'.format(linelist_dir, element), '{}/{}_null'.format(workdir, element), run_synthesis = False)

    _global['response'] = {element: {abun: False for abun in abun_offsets[element]} for element in settings['elements']}
    _global['masks'] = {}
    null_cont, null_line = np.loadtxt('{}/null'.format(workdir) + '/synthe_1/spectrum.asc', skiprows = 2, usecols = [2, 3], unpack = True)
    _global['null'] = {'cont': null_cont, 'line': null_line}
    _global['null_offsets'] = null_offsets

    # Compute the masks for each element
    for element in settings['elements']:
        response_max = '{}/{}_max'.format(workdir, element); response_min = '{}/{}_min'.format(workdir, element)
        if not os.path.isdir(response_max) or not os.path.isdir(response_min):
            notify('Calculating the mask for {}'.format(element), color = 'y')
        update_continuum = workdir + '/null/continuum.10'
        update_opacity = [[workdir + '/{}_null/opacity.9'.format(element), -1.0], [workdir + '/null/opacity.9', +1.0]]
        if not os.path.isdir(response_dir := response_max):
            synthesize(structure_dir, combine_abundances(null_offsets, {element: np.max(abun_offsets[element])}), '{}/{}'.format(linelist_dir, element), response_dir, update_opacity = update_opacity, update_continuum = update_continuum)
            clean_response(response_dir)
        if not os.path.isdir(response_dir := response_min):
            synthesize(structure_dir, combine_abundances(null_offsets, {element: np.min(abun_offsets[element])}), '{}/{}'.format(linelist_dir, element), response_dir, update_opacity = update_opacity, update_continuum = update_continuum)
            clean_response(response_dir)
        spectrum_high = np.loadtxt('{}/{}_max'.format(workdir, element) + '/synthe_1/spectrum.asc', skiprows = 2, usecols = [3])
        spectrum_low  = np.loadtxt('{}/{}_min'.format(workdir, element) + '/synthe_1/spectrum.asc', skiprows = 2, usecols = [3])
        _global['masks'][element] = np.abs(spectrum_high - spectrum_low) > settings['threshold']
        # If this is not one of the higher_order_impact, we can consider the max and min spectra to be valid response functions
        # Otherwise, we want to recompute them later with updated continua
        if element not in settings['higher_order_impact']:
            _global['response'][element][np.max(abun_offsets[element])] = spectrum_high[_global['masks'][element]]
            _global['response'][element][np.min(abun_offsets[element])] = spectrum_low[_global['masks'][element]]

    for element in settings['elements']:
        _global['response'][element][0.0] = _global['null']['line'][_global['masks'][element]]

    return {param: [settings['gridfit_params'][param]] for param in settings['gridfit_params']}

def wl_grid():
    """Helper function to get the standard SYNTHE wavelength grid based on the synthesis parameters in the settings
    
    Returns
    -------
    array_like
        Standard wavelength grid
    """
    return np.exp(np.arange(np.ceil(np.log(settings['wl_start']) / (lgr := np.log(1.0 + 1.0 / settings['res']))), np.floor(np.log(settings['wl_end']) / lgr) + 1) * lgr) * 10

def read_grid_model(params, grid):
    """Retrieve the wavelength grid of the spectral models and trim it to accommodate redshift corrections
    
    Normally, this function loads the spectral model from the model grid. Since we do not use a model grid,
    we use this function only to retrieve the wavelength grid on which the response functions will be calculated
    (the standard SYNTHE grid in `wl_grid()`). The flux vector is filled with ones, and it will be completely
    replaced with the actual spectrum in `preprocess_grid_model()`, assembled from response functions

    We trim the wavelength grid on both sides by whatever amount necessary to make sure that the model remains
    within that wavelength range at all allowed redshift values
    
    Parameters
    ----------
    params : dict
        Not used
    grid : dict
        Not used
    
    Returns
    -------
    wl : array_like
        Trimmed wavelength grid
    flux : array_like
        Placeholder for a flux array
    meta : dict
        Not used
    """
    wl = wl_grid()
    flux = np.ones(len(wl))

    # Trim the spectrum on both sides to make sure we can do redshift corrections
    wl_range = [np.min(wl * (1 + settings['virtual_dof']['redshift'][1] * 1e3 / scp.constants.c)), np.max(wl * (1 + settings['virtual_dof']['redshift'][0] * 1e3 / scp.constants.c))]
    mask_left = wl < wl_range[0]; mask_right = wl > wl_range[1]; mask_in = (~mask_left) & (~mask_right)
    return wl[mask_in], flux[mask_in], {}

def preprocess_grid_model(wl, flux, params, meta):
    """Assemble a spectral model for a given chemical composition using response functions
    
    This function calculates and caches all of the response functions it needs at runtime
    
    Parameters
    ----------
    wl : array_like
        Wavelength grid to calculate the spectral model at
    flux : array_like
        Placeholder for a flux array, not used
    params : dict
        Model parameters, including the abundances of individual elements
    meta : dict
        Not used
    
    Returns
    -------
    array_like
        The flux array of the calculated spectral model
    """
    global _global

    # The input is provided in terms of abundance offsets from the gridfit fit. We want to convert those into offsets from the null spectrum
    requested = {element: params[element] for element in settings['elements']}
    for element in settings['initial_abundance_offsets']:
        requested[element] -= np.round(settings['initial_abundance_offsets'][element], 2)

    # Determine which response functions are needed to generate this model
    required = []
    sides = {}
    for element in settings['elements']:
        abun = np.array(list(_global['response'][element].keys()))
        if requested[element] in abun:
            required += [[element, requested[element]]]
            sides[element] = [requested[element]]
        else:
            required += [[element, np.max(abun[abun < requested[element]])]]
            required += [[element, np.min(abun[abun > requested[element]])]]
            sides[element] = [required[-2][1], required[-1][1]]

    # Determine which of those are yet to be calculated
    batch = []
    for response in required:
        if type(_global['response'][response[0]][response[1]]) is bool:
            batch += [response]

    # Calculate the missing response functions
    if len(batch) != 0:
        for response in batch:
            if response[0] not in settings['higher_order_impact']:
                update_continuum = _global['workdir'] + '/null/continuum.10'
            else:
                update_continuum = False
            update_opacity = [[_global['workdir'] + '/{}_null/opacity.9'.format(response[0]), -1.0], [_global['workdir'] + '/null/opacity.9', +1.0]]
            if not os.path.isdir(response_dir := (_global['workdir'] + '/{}_{}'.format(*response))):
                notify('Calculating response function for {}={}'.format(*response), color = 'y')
                synthesize(_global['structure'], combine_abundances(_global['null_offsets'], {response[0]: response[1]}), '{}/{}'.format(_global['linelist'], response[0]), response_dir, update_opacity = update_opacity, update_continuum = update_continuum)
                clean_response(response_dir)
            spectrum = np.loadtxt(_global['workdir'] + '/{}_{}'.format(*response) + '/synthe_1/spectrum.asc', skiprows = 2, usecols = [3])
            _global['response'][response[0]][response[1]] = spectrum[_global['masks'][response[0]]]

    # Assemble the spectrum
    flux = _global['null']['line'] * 1.0
    for element in settings['elements']:
        if len(sides[element]) == 1:
            response = _global['response'][element][sides[element][0]]
        else:
            response_high = _global['response'][element][sides[element][1]]
            response_low = _global['response'][element][sides[element][0]]
            response = response_low + (response_high - response_low) * (requested[element] - sides[element][0]) / (sides[element][1] - sides[element][0])
        flux[_global['masks'][element]] += response - _global['null']['line'][_global['masks'][element]]
    flux *= _global['null']['cont']

    # Apply redshift
    wl_full = wl_grid()
    wl_redshifted = wl_full * (1 + params['redshift'] * 1e3 / scp.constants.c)
    flux = np.interp(wl, wl_redshifted, flux)

    return flux

def extract_results(fit, propagate_gridfit = False, detector_wl = False):
    """Extract best-fit abundances and errors from the `chemfit.chemfit()` output
    
    The function re-expresses all abundances with respect to hydrogen ([X/H]), such that we do
    not need to worry about uncertainties in best-fit metallicity and alpha-enhancement.

    If `propagate_gridfit` is `True`, the function will also update the errors in the best-fit
    abundances to include the contributions due to uncertainties in Teff, log(g) and [C/M]

    The error propagation is only approximate. It assumes local linearity of the parameter space,
    ignores statistical and photometric priors and assumes that the derivatives of the objective
    function with repsect to stellar parameters do not depend on the abundances of individual elements.
    It also assumes that the exact same spectrum and masking were used to obtain the stellar parameters,
    as the ones used in localfit
    
    Parameters
    ----------
    fit : dict
        Output of `chemfit.chemfit()`
    propagate_gridfit : bool, optional
        Set to `True` to propagate the uncertainties in Teff, log(g) and [C/M] into the errors
        in best-fit abundances. Defaults to `False`
    detector_wl : dict, optional
        If `propagate_gridfit` is set to `True`, this argument must be set to the observed wavelength
        array, i.e. the `wl` argument given to `chemfit.chemfit()`
    
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
    global _global

    # Re-express all abundances with respect to hydrogen
    fit['extra']['abun'] = {'abun': {}, 'errors': {}}
    for element in settings['elements']:
        ratio = '[{}/H]'.format(element)
        fit['extra']['abun']['abun'][ratio] = fit['fit'][element]
        fit['extra']['abun']['errors'][ratio] = fit['errors'][element]
        fit['extra']['abun']['abun'][ratio] += np.round(settings['gridfit_params']['zscale'], 3)
        if element in settings['alpha_elements']:
            fit['extra']['abun']['abun'][ratio] += np.round(settings['gridfit_params']['alpha'], 2)

    # Propagate the errors in teff, logg and carbon into the abundances
    if propagate_gridfit:
        notify('Propagating gridfit uncertainties into best-fit abundances', color = 'y')

        # Helper function to adapt model spectra to observations (downsampling, rebinning, masking, continuum-correction)
        def observe(wl, flux, detector_wl):
            ds_wl, ds_flux = main__simulate_observation(wl, flux, detector_wl = detector_wl)
            cont = main__estimate_continuum(ds_wl, fit['extra']['observed']['flux'] / ds_flux, fit['extra']['observed']['ivar'] * ds_flux ** 2, npix = settings['cont_pix'], k = settings['spline_order'], arm_index = fit['extra']['arm_index'])
            return (cont * ds_flux)[fit['extra']['mask']]

        # Set up the Jacobian matrix. Note that we exclude zscale and alpha here because the abundances are expressed
        # on the absolute scale ([X/H])
        dof = ['teff', 'logg', 'carbon'] + list(settings['virtual_dof'])
        jacobian = np.zeros([np.count_nonzero(fit['extra']['mask']), len(dof)])

        # First, evaluate the Jacobian for the gridfit parameters: teff, logg and carbon. To do this, we need to
        # synthesize spectra with small offsets in these parameters. We use interpolated structures
        cases = {
            'nominal': {'teff_offset': 0.0, 'logg_offset': 0.0, 'carbon_offset': 0.0},
            'teff': {'teff_offset': 1.0, 'logg_offset': 0.0, 'carbon_offset': 0.0},
            'logg': {'teff_offset': 0.0, 'logg_offset': 0.01, 'carbon_offset': 0.0},
            'carbon': {'teff_offset': 0.0, 'logg_offset': 0.0, 'carbon_offset': 0.1},
        }
        for case in cases:
            header = atlas.restarts.load_header()
            params = copy.deepcopy(settings['gridfit_params'])
            # Apply the offsets. If the offsets take us outside the structure grid, flip the sign of the offset
            for param in ['teff', 'logg', 'carbon']:
                params[param] += cases[case]['{}_offset'.format(param)]
                if params[param] > np.max(header[param]):
                    params[param] -= 2 * cases[case]['{}_offset'.format(param)]
                    cases[case]['{}_offset'.format(param)] = -cases[case]['{}_offset'.format(param)]
            # Produce an interpolated structure for the offset stellar parameters
            structure_name = '{}_interpolated'.format(hashlib.md5(pickle.dumps(params)).hexdigest())
            structure_dir = '{}/localfit_structures/{}'.format(settings['scratch'], structure_name)
            if not os.path.isdir(structure_dir):
                structure = atlas.restarts.interpolate_structure(params, header = header)
                model = atlas.restarts.generate_model(*structure)
                os.mkdir(structure_dir)
                f = open('{}/output_summary.out'.format(structure_dir), 'w'); f.write(model); f.close()
                f = open('{}/output_last_iteration.out'.format(structure_dir), 'w'); f.close()
                f = open('{}/output_main.out'.format(structure_dir), 'w'); f.close()
            # Synthesize the spectrum
            if not os.path.isdir(spectrum_dir := ('{}/spectrum'.format(structure_dir))):
                synthesize(structure_dir, {}, _global['linelist'] + '/null', spectrum_dir)
            cases[case]['spectrum'] = observe(*np.loadtxt(spectrum_dir + '/synthe_1/spectrum.asc', skiprows = 2, usecols = [0, 1], unpack = True), detector_wl)
        # Calculate the derivatives and put them in the Jacobian
        jacobian[:,0] = (cases['teff']['spectrum'] - cases['nominal']['spectrum']) / cases['teff']['teff_offset']
        jacobian[:,1] = (cases['logg']['spectrum'] - cases['nominal']['spectrum']) / cases['logg']['logg_offset']
        jacobian[:,2] = (cases['carbon']['spectrum'] - cases['nominal']['spectrum']) / cases['carbon']['carbon_offset']

        # Now evaluate the Jacobian for all of the other parameters, which we can do using finite differences and
        # `read_grid_model` / `preprocess_grid_model()`
        wl_grid, flux, meta = read_grid_model(None, None)
        params = copy.deepcopy(fit['fit'])
        nominal = observe(wl_grid, preprocess_grid_model(wl_grid, flux, params, None), detector_wl = detector_wl)
        step = 0.01
        for i in range(3, len(dof)):
            params = copy.deepcopy(fit['fit'])
            # As before, step in the other direction if the bounds are exceeded
            if params[dof[i]] + step > settings['virtual_dof'][dof[i]][1]:
                diff = -step
            else:
                diff = step
            params[dof[i]] += diff
            # Compute the synthetic spectrum and put it in the Jacobian
            spectrum = observe(wl_grid, preprocess_grid_model(wl_grid, flux, params, None), detector_wl = detector_wl)
            jacobian[:,i] = (spectrum - nominal) / diff

        # Now that we have the Jacobian matrix fully populated, we can compute the covariance matrix using the standard
        # error propagation formula:
        #            COV = (J^T x diag(IVAR) x J)^-1
        # We also scale the covariance matrix by the reduced chi-squared (SUM((OBSERVED - MODEL)^2 * IVAR) / (N_points - N_params))
        # to mimick scp.optimize.curve_fit()'s `absolute_sigma = False` setting

        # Compute reduced chi squared for the final fit
        residuals = (fit['extra']['observed']['flux'] - (fit['extra']['model']['cont'] * fit['extra']['model']['flux']))[fit['extra']['mask']] * fit['extra']['observed']['ivar'][fit['extra']['mask']] ** 0.5
        chi2_red = np.sum(residuals ** 2) / (np.shape(jacobian)[0] - np.shape(jacobian)[1])

        # Compute the covariance matrix and update the errors in abundances
        weighted_jacobian = jacobian * (fit['extra']['observed']['ivar'] ** 0.5)[fit['extra']['mask']][:, np.newaxis]
        fit['extra']['abun']['cov'] = np.linalg.inv(weighted_jacobian.T @ weighted_jacobian) * chi2_red
        fit['extra']['abun']['dof'] = dof
        for element in settings['elements']:
            fit['extra']['abun']['errors']['[{}/H]'.format(element)] = np.sqrt(fit['extra']['abun']['cov'][tuple([dof.index(element)] * 2)])

    return fit

def public__localfit(wl, flux, ivar, gridfit, level = 1, niter = 5):
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
    level : number, optional
        Level of analysis. In the order of increasing accuracy, the available levels are:
            1 : Adopt an interpolated model atmosphere at `gridfit` parameters. Calculate
                the null-spectrum at the same `gridfit` parameters using that model
                atmosphere. Fit for element abundances with response functions
            2 : Same as level 1, but compute a new model atmosphere at `gridfit` parameters
                instead of using an interpolated model
            3 : Same as level 1, but repeat the abundance determination `niter` times,
                updating the abundances of the null-spectrum to the output of the previous
                iteration
            4 : Same as level 3, but compute a new model atmosphere at `gridfit` parameters
                instead of using an interpolated model
            5 : Same as level 4, but recompute the model atmosphere at the end of each
                iteration to the determined abundances
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

    params = copy.deepcopy(gridfit)
    intermediate = [] # Storage for intermediate abundances at the end of each iteration

    for iteration in range([1, niter][level in [3, 4, 5]]):
        notify('*** Starting iteration {} ***'.format(iteration + 1), color = 'm')
        fit = main__chemfit(wl, flux, ivar, initial = params)
        for element in settings['elements']:
            params[element] = fit['fit'][element]
            settings['initial_abundance_offsets'][element] = fit['fit'][element]
        fit = extract_results(fit, propagate_gridfit = False)
        settings['use_initial_abundance_offsets_in_structure'] = level == 5
        intermediate += [copy.deepcopy(fit['extra']['abun'])]
        notify('Completed iteration {}: {}\n'.format(iteration + 1, fit['extra']['abun']), color = 'm')
    notify('*** Completed all iterations ***')

    fit = extract_results(fit, propagate_gridfit = True, detector_wl = wl)
    fit['extra']['abun']['intermediate'] = intermediate

    return fit
