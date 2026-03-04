############################################################
#                                                          #
#                         PFS PRESET                       #
#                                                          #
#   This preset configures chemfit to analyze PFS spectra. #
#   It defines typical wavelength coverage and resultion   #
#   of the arms of PFS                                     #
#                                                          #
############################################################


settings = {
    ### Spectrograph settings ###
    'arms': {
        'blue': {
            'FWHM': 2.07,
            'wl': np.linspace(3800, 6500, 4096),
        },
        'red_lr': {
            'FWHM': 2.63,
            'wl': np.linspace(6300, 9700, 4096),
        },
        'red_mr': {
            'FWHM': 1.368,
            'wl': np.linspace(7100, 8850, 4096),
        },
        'ir': {
            'FWHM': 2.4,
            'wl': np.linspace(9400, 12600, 4096),
        }
    },
}
