import numpy as np

def C12C13_kirby_2015(logg):
    """
    Relationship between C12/C13 isotope ratio and surface gravity from Kirby+2015
    """
    return np.where(logg > 2.7, 50, np.where(logg <= 2.0, 6, 63 * logg - 120))

def VTURB_LOGG(logg):
    """
    Relationship between VTURB and LOGG based on Deimos spectra of M13
    """
    return 2.792 * np.exp(-0.241 * logg -0.118)

settings = {
    'scripts_dir': '/expanse/lustre/projects/csd835/rgerasim/pfsgrid/scripts/',
    'runs_dir': '/expanse/lustre/projects/csd835/rgerasim/pfsgrid/runs/',
    'results_dir': '/expanse/lustre/projects/csd835/rgerasim/pfsgrid/results/',

    'ODF': {
        'output_dir': '/expanse/lustre/projects/csd835/rgerasim/pfsgrid/ODF/',
    },

    'ATLAS': {
        'output_dir': '/expanse/lustre/projects/csd835/rgerasim/pfsgrid/ATLAS/',
        'max_flux_error': 1.0,
        'max_flux_error_derivative': 10.0,
    },

    'SYNTHE': {
        'wl_start': 350,
        'wl_end': 1300,
        'res': 300000,
        'air_wl': True,
        'vturb': {np.round(logg, 2): np.round(VTURB_LOGG(logg), 2) for logg in np.arange(0.0, 5.001, 0.25)},
        'atoms': 'BasicATLAS',
        'C12C13': {np.round(logg, 2): np.round(C12C13_kirby_2015(logg), 2) for logg in np.arange(0.0, 5.001, 0.25)},
    },

    'rescalc': {
        'linelist': '/expanse/lustre/projects/csd835/rgerasim/pfsgrid/linelist/',
        'threshold': 0.01,
        'abun': list(np.round(np.arange(-1.0, 1.01, 0.1), 1)),
        'elements': {
           # 'Ba': [56, 56.01, 56.02, 56.03, 56.04, 56.05],
           # 'Eu': [63, 63.01, 63.02, 63.03, 63.04, 63.05],
           # 'Mg': [12, 12.01, 12.02, 12.03, 12.04, 12.05, 112.0, 812.0, 112.01, 10812.0],
           # 'Si': [14, 14.01, 14.02, 14.03, 14.04, 14.05, 114.0, 814.0, 114.01, 80814.0, 101010114.0, 614.0, 60614.0],
           # 'Ca': [20, 20.01, 20.02, 20.03, 20.04, 20.05, 120.0, 820.0, 120.01, 10820.0],
           # 'Ti': [22, 22.01, 22.02, 22.03, 22.04, 22.05, 122.0], #, 822.0],
           # 'Ni': [28, 28.01, 28.02, 28.03, 28.04, 28.05, 128.0, 828.0],
           # 'Na': [11, 11.01, 11.02, 11.03, 11.04, 11.05, 111.0, 10811.0],
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
        # For which elements do we want to model the effect on the lines of other elements and continuum?
        'higher_order_impact': ['Mg', 'Si', 'Na', 'Ca', 'Al'],

        'conserve_space': True,
    },

    'slurm': {
        'queue': 'shared',
        'account': 'csd958',
        'nodes': 1,
        'mem_per_cpu': 2000,
    }
}


