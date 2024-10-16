"""
@author: I.Escala (iescala@carnegiescience.edu)

A Python implentation of feh.pro (E.N. Kirby), built for compatability with
PARSEC isochrones
CMD v3.3 defaults
filter: Subaru/Hyper-Suprime Cam (AB mags)

Linear age: 1 - 12 Gyr in steps of 1 Gyr
[M/H]: -2 to 0.3 dex in steps of 0.01 dex
"""

from __future__ import absolute_import
from scipy.interpolate import griddata
import numpy as np
import glob
import sys

# MNI --BEGIN--
from scipy import interpolate
import os

# MNI --END--


Class phot():
    def __init__(self, mags, mags_err=None, filters=None):
        self.mag = mags[:,1]
        self.clr = mags[:,0]-mags[:,1]
        if mags_err is not None:
            self.magerr = mags_err[:,1]
            self.clrerr = np.sqrt(mags_err[:,0]**2 + mags_err[:,1]**2)
        
    def measure_parameters(self, dm=None, ddm=None, extrap=True, enforce_bounds=False):

        """
        PURPOSE: Calculate photometric metallicity with PARSEC isochrones.

        CALLING SEQUENCE:
        -----------------
        phot = from_phot(mag_in, clr_in, err_mag_in=err_mag_in, err_color_in=err_color_in,
                dm=dm, ddm=ddm)

        INPUTS:
        -------
        mag_in: An array of extinction-corrected magnitudes.
        clr_in: An array of dereddened colors.  mag_in is always the
        redder filter. For example, if your color is V-I, then mag_in
        must be I.  clr_in must have the same number of elements as
        mag_in.
        dm: The absolute distance modulus (m-M)_0.  It should have
            one element for all the input stars, or it should have as many
            elements as the input stars.  In other words, each star can
            have its own distance modulus.
        ddm: The error on the distance modulus.

        err_mag_in: An array of photometric errors on the magnitudes.  It
            must have the same number of elements as mag_in.
        err_clr_in: An array of photometric errors on the colors.  It
            must have the same number of elements as mag_in.
        extrap: boolean: Whether to extrapolate beyond the most metal-poor isochrone and beyond
                the tip of the red giant branch. Default is true.
        enforce_bounds: boolean: If a value of [Fe/H], Teff, or logg that is interpolated
                        based on the isochrones is outside of the range of the isochrones,
                        set the values to NaN. Default is False.

        OUTPUTS:
        --------
        A structure with as many elements as mag_in.  The structure
        contains the following tags:
        AGES:    An array of stellar ages in Gyr.  The number of ages
            depends on the choice of isochrone.
        FEH:     An array of photometric metallicities corresponding to
            each age.
        TEFF:    An array of effective temperatures (K) corresponding to
            each age.
        LOGG:    An array of surface gravities (log cm/s^2) corresponding
            to each age.
        ERR_FEH: An array of errors on the photometric metallicity
            corresponding to each element in FEH.  This tag exists only if
            the inputs include photometric errors.
        ERR_TEFF: An array of errors on the effective temperatures
            corresponding to each element in TEFF.  This tag exists only
            if the inputs include photometric errors.
        ERR_TEFF_ARR: A 2D array containing the full MCMC distribution for the photometric
            errors, with dimensions equal to NAGES (number of stellar ages for a given isochrone)
            and NMC (the number of iterations in the MCMC that determines the photometric errors).
        ERR_LOGG: An array of errors on the surface gravities
            corresponding to each element in LOGG.  This tag exists only
            if the inputs include photometric errors.
        """

        kind = 'parsec'; filter = 'hsc'; Zsun = 0.0152

        if err_mag_in is None and err_color_in is None:
            err_mag_in = np.zeros(len(mag_in)); err_color_in = np.zeros(len(color_in))

        if isinstance(mag_in, list) or isinstance(mag_in, np.ndarray):
            pass
        else:
            mag_in = [mag_in]; color_in = [color_in]
            if err_mag_in is not None and err_color_in is not None:
                err_mag_in = [err_mag_in]; err_color_in = [err_color_in]



        mag_in = np.array(mag_in); color_in = np.array(color_in)
        err_mag_in = np.array(err_mag_in); err_color_in = np.array(err_color_in)

        n = len(mag_in)
        dm_uniq = np.unique(dm)
        ndm = len(dm_uniq)
        if ndm == 1 and n > 1: dm = np.full(n, dm)

        if ddm is not None:
            err_mag_in = np.sqrt(err_mag_in**2. + ddm**2.)


        # MNI

        ## ORIGINAL
        #iso_files = glob.glob(f'isochrones/feh_{filter}_*_{kind}.dat')
        ## ORIGINAL

        iso_files = glob.glob(os.path.dirname(__file__) + f'/../isochrones/feh_{filter}_*_{kind}.dat')

        # MNI

        nages = len(iso_files)
        nmc = 1000

        phots = {}
        for i in range(n):
            phot = {}
            phot['mag'] = mag_in[i]
            phot['color'] = color_in[i]
            phot['err_mag'] = err_mag_in[i]
            phot['err_color'] = err_color_in[i]
            phot['ages'] = np.full(nages, -999.)
            phot['feh'] = np.full(nages, -999.)
            phot['teff'] = np.full(nages, -999.)
            phot['logg'] = np.full(nages, -999.)
            phot['err_feh'] = np.full(nages, -999.)
            phot['err_teff'] = np.full(nages, -999.)
            phot['err_logg'] = np.full(nages, -999.)
            phot['err_teff_arr'] = np.full((nages, nmc), -999.)
            phots[i] = phot

        for kk,file in enumerate(iso_files):

            usecols = (0,2,7,8,9,25,27)

            Zi, age, Te, logg, label, mag_blue, mag_red = np.loadtxt(file, skiprows=12, unpack=True,
                                                                            usecols=usecols)
            feh_val = np.unique(Zi)

            # MNI -- BEGIN --

            age = 10**age

            # MNI -- END --


            print(f'{kind.upper()} isochrone set using {filter.upper()} filter, assuming age = '+\
            f'{str(int(np.unique(age)[0]/1.e9))} Gyr')



            mag_grid = []; clr_grid = []; feh_grid = []
            logt_grid = []; logg_grid = []

            for ii in range(len(feh_val)):

                wrgb = np.where((label == 3) & (Zi == feh_val[ii]) )[0]

                mag = mag_red[wrgb]
                color = mag_blue[wrgb]-mag_red[wrgb]
                feh = np.log10(Zi[wrgb]/Zsun)
                grav = logg[wrgb]
                teff = Te[wrgb]

                if extrap:

                    #Extrapolate beyond the tip of the red giant branch
                    m, b = np.polyfit(color[-2:], mag[-2:], 1)

                    #Extend the colors redward by 1.2 mags
                    xx = np.arange(color[-1], color[-1]+1.2, 0.1)
                    yy = xx*m + b

                    #Fill in this extension with the values for Teff, Logg,
                    #and [Fe/H] from the TRGB
                    zz = np.full(len(xx), feh[-1])
                    gg = np.full(len(xx), grav[-1])
                    tt = np.full(len(xx), teff[-1])

                    #Extend the grid with these values
                    color = np.append(color, xx)
                    mag = np.append(mag, yy)
                    feh = np.append(feh, zz)
                    grav = np.append(grav, gg)
                    teff = np.append(teff, tt)

                mag_grid.extend(mag)
                clr_grid.extend(color)
                feh_grid.extend(feh)
                logg_grid.extend(grav)
                logt_grid.extend(teff)

            mag_grid = np.array(mag_grid); clr_grid = np.array(clr_grid)
            logg_grid = np.array(logg_grid); logt_grid = np.array(logt_grid)
            feh_grid = np.array(feh_grid)




            for jj in range(ndm):

                mag_grid += dm_uniq[jj]
                wdm = np.where( dm == dm_uniq[jj] )
                if len(wdm) == 1: wdm = wdm[0]

                feh_i = griddata((clr_grid, mag_grid), feh_grid, (color_in[wdm], mag_in[wdm]))
                logt_i = griddata((clr_grid, mag_grid), logt_grid, (color_in[wdm], mag_in[wdm]))
                logg_i = griddata((clr_grid, mag_grid), logg_grid, (color_in[wdm], mag_in[wdm]))
                t_i = 10**logt_i

                out_of_bounds = np.isnan(feh_i)

                #Extrapolate blueward of the bluest isochrone ONLY for Teff and Logg
                #Use interp2d for these points ONLY since it is less robust than griddata
                if len(feh_i[out_of_bounds]) > 0:

                    #Identify indices for points that are out of bounds
                    arg = np.arange(0, len(feh_i))[out_of_bounds]

                    #Construct an interpolator that can extrapolate
                    #Note that griddata cannot extrapolate, but it is much more robust than *Spline functions
                    #*Spline functions very dependent on smoothing parameter s and can generate wild swings

                    # ORIGINAL

                    #f_t = SmoothBivariateSpline(clr_grid, mag_grid, logt_grid, kx=2, ky=2)
                    #f_g = SmoothBivariateSpline(clr_grid, mag_grid, logg_grid, kx=2, ky=2)
                    #f_f = SmoothBivariateSpline(clr_grid, mag_grid, feh_grid, kx=2, ky=2)
                    # ORIGINAL

                    # MNI
                    f_t = interpolate.SmoothBivariateSpline(clr_grid, mag_grid, logt_grid, kx=2, ky=2)
                    f_g = interpolate.SmoothBivariateSpline(clr_grid, mag_grid, logg_grid, kx=2, ky=2)
                    f_f = interpolate.SmoothBivariateSpline(clr_grid, mag_grid, feh_grid, kx=2, ky=2)

                    # MNI


                    for ll in arg:

                        logt_l = f_t(color_in[wdm][ll], mag_in[wdm][ll])[0][0]
                        logg_l = f_g(color_in[wdm][ll], mag_in[wdm][ll])[0][0]
                        feh_l = f_f(color_in[wdm][ll], mag_in[wdm][ll])[0][0]

                        logt_i[ll] = logt_l
                        logg_i[ll] = logg_l
                        feh_i[ll] = feh_l
                        t_i[ll] = 10**logt_l




                if enforce_bounds:
                    w = np.where( (feh_i < feh_grid.min()) | (feh_i > feh_grid.max()) )[0]
                    if len(w) > 0:
                        feh_i[w] = np.nan

                    w = np.where( (logt_i < logt_grid.min()) | (logt_i > logt_grid.max()) )[0]

                    if len(w) > 0:
                        t_i[w] = np.nan

                    w = np.where( (logg_i < logg_grid.min()) | (logg_i > logg_grid.max()) )[0]
                    if len(w) > 0:
                        logg_i[w] = np.nan

                #LOOK AT THIS AS A POTENTIAL SOURCE OF THE ISSUE: WORKS FOR APOGEE, BUT DOES IT WORK FOR
                #PANDAS?


                for iw in range(len(wdm)):

                    phots[wdm[iw]]['feh'][kk] = feh_i[iw]
                    phots[wdm[iw]]['teff'][kk] = t_i[iw]
                    phots[wdm[iw]]['ages'][kk] = np.unique(age)
                    phots[wdm[iw]]['logg'][kk] = logg_i[iw]

                if err_mag_in is not None and err_color_in is not None:
                    for ii in range(len(wdm)):

                        if np.isnan(feh_i[ii]):
                            phots[wdm[ii]]['err_feh'][kk] = 0.

                        else:
                            rand = np.random.normal(size=(nmc,2))
                            clr_mc = color_in[wdm[ii]] + rand[:,0]*err_color_in[wdm[ii]]
                            mag_mc = mag_in[wdm[ii]] + rand[:,1]*err_mag_in[wdm[ii]]

                            feh_mc = griddata((clr_grid, mag_grid), feh_grid, (clr_mc, mag_mc), fill_value=-999.)
                            wgood = np.where( (feh_mc >= -5.) & (feh_mc <= 2.) )[0]

                            if len(wgood) > 2:
                                phots[wdm[ii]]['err_feh'][kk] = np.std(feh_mc[wgood])
                            else:
                                phots[wdm[ii]]['err_feh'][kk] = np.nan

                        if np.isnan(t_i[ii]):
                            phots[wdm[ii]]['err_teff'][kk] = 0.

                        else:
                            rand = np.random.normal(size=(nmc,2))
                            clr_mc = color_in[wdm[ii]] + rand[:,0]*err_color_in[wdm[ii]]
                            mag_mc = mag_in[wdm[ii]] + rand[:,1]*err_mag_in[wdm[ii]]

                            logt_mc = griddata((clr_grid, mag_grid), logt_grid, (clr_mc, mag_mc), fill_value=-999.)
                            wgood = np.where( (10**logt_mc >= 3000.) & (10**logt_mc <= 10000) )[0]


                            if len(wgood) > 2:
                                phots[wdm[ii]]['err_teff'][kk] = np.std(10**logt_mc[wgood])
                            else:
                                phots[wdm[ii]]['err_teff'][kk] = np.nan

                            for mm in range(nmc):
                                phots[wdm[ii]]['err_teff_arr'][kk][mm] = 10**logt_mc[mm]

                        if np.isnan(logg_i[ii]):
                            phots[wdm[ii]]['err_logg'][kk] = 0.

                        else:
                            rand = np.random.normal(size=(nmc,2))
                            clr_mc = color_in[wdm[ii]] + rand[:,0]*err_color_in[wdm[ii]]
                            mag_mc = mag_in[wdm[ii]] + rand[:,1]*err_mag_in[wdm[ii]]

                            logg_mc = griddata((clr_grid, mag_grid), logg_grid, (clr_mc, mag_mc), fill_value=-999.)
                            wgood = np.where( (logg_mc >= -2.) & (logg_mc <= 7.) )[0]

                            if len(wgood) > 2:
                                phots[wdm[ii]]['err_logg'][kk] = np.std(logg_mc[wgood])
                            else:
                                phots[wdm[ii]]['err_logg'][kk] = np.nan

                mag_grid -= dm_uniq[jj]

        #Sort the dictionary according to stellar age
        wsort = np.argsort(phots[0]['ages'])




        for i in range(len(phots)):

            phots[i]['ages'] = phots[i]['ages'][wsort]
            phots[i]['feh'] = phots[i]['feh'][wsort]
            phots[i]['teff'] = phots[i]['teff'][wsort]
            phots[i]['logg'] = phots[i]['logg'][wsort]

            if err_mag_in is not None and err_color_in is not None:
                phots[i]['err_feh'] = phots[i]['err_feh'][wsort]
                phots[i]['err_teff'] = phots[i]['err_teff'][wsort]
                phots[i]['err_logg'] = phots[i]['err_logg'][wsort]
                phots[i]['err_teff_arr'] = phots[i]['err_teff_arr'][wsort]

        return phots
