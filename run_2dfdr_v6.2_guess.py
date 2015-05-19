#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

""" 
Cross-correlate individual spectra against the AMBRE grid and provide initial
guesses. This is effectively what GUESS will do, but Ireland & Lin have yet to
respond.
"""

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import cPickle as pickle
import matplotlib
import os
from collections import OrderedDict
from glob import glob
from warnings import simplefilter
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.constants import c as speed_of_light

import specutils
from convert import motions

# Suppress "polyfit may be poorly conditioned" messages
simplefilter("ignore", np.RankWarning)
c = speed_of_light.to("km/s").value

# CONTROLS
UPDATE_FILES = False
SAVE_BEST_FIT_MODELS = False
continuum_degree = {
    1: 4,
    2: 4,
    3: 4,
    4: 3
}
ccf_mask = [
    [7500, 7700]
]
chi_sq_mask = [
    [7590, 7617],
    [7620, 7700]
]

vrad_limits = (-1000, 1000)
standardised_file_folder = "/home/acasey/miner3/galah/standardised_spectra/2dfdr_v6.2/"
output_figures_folder = "/home/acasey/miner3/galah/standardised_spectra/2dfdr_v6.2_figures/"
output_models_folder = ""


# CODE
def get_unique_id(filename):
    base = filename.split("/")[-1]
    date_fco_ccd, fib_galah_id = base.split("_")
    date_fco_ccd = date_fco_ccd[:-1] # strip the ccd number from the right
    fib_galah_id = fib_galah_id[:-5] # strip the .fits text from the right
    return "?_".join([date_fco_ccd, fib_galah_id]) + ".fits"

# Identify unique stars
standardised_files = glob(os.path.join(standardised_file_folder, "*.fits"))
id_filename_masks =  map(get_unique_id, standardised_files)

unique_filename_masks = set(id_filename_masks)

# Load in the synthetic spectra from the pickled fluxes.
# Available from: https://zenodo.org/record/14977/files/galah-ambre-grid.pkl
with open("galah-ambre-grid.pkl", "rb") as fp:
    model_grid, model_wavelengths, model_intensities, pixels = pickle.load(fp)
    model_intensities = model_intensities.reshape(model_grid.size, sum(pixels))

# Create a set of rows to store tabular data in.
N, rows = len(unique_filename_masks), []
for i, filename_mask in enumerate(unique_filename_masks):
    print("At star {0}/{1}: {2}".format(i + 1, N, filename_mask))

    # get the spectra
    filenames = glob(os.path.join(standardised_file_folder, filename_mask))

    print("Found {0} spectra matching mask {1}".format(len(filenames),
        filename_mask))

    images = map(fits.open, filenames)

    # Check that MOON_DEG exists, otherwise put it in.
    if "MOON_DEG" not in images[0][0].header:
        moon_deg = motions.moon_distance(images[0][0].header)
        print("Inserting MOON_DEG ({0:.1f}) into image headers".format(moon_deg))
        for image in images:
            image[0].header["MOON_DEG"] = moon_deg
            image[0].header.comments["MOON_DEG"] \
                = "Distance to the moon during observing (deg)"

    # Get a list of all the CCDs we actually have.
    ccds = [image[0].header["CCD"] for image in images]

    # Do CCF against all of the arms.
    ccf_results = {}
    for filename, image in zip(filenames, images):

        obs_dispersion = image[0].header["CRVAL1"] + image[0].header["CDELT1"] \
            * (np.arange(image[0].data.size) - image[0].header["CRPIX1"])
        ccd = image[0].header["CCD"]
        si, ei = map(sum, (pixels[:ccd - 1], pixels[:ccd]))

        # Apply the CCF mask.
        ccf_flux_data = image[0].data.copy()
        for mask_region in ccf_mask:
            sj, ej = np.searchsorted(obs_dispersion, mask_region)
            ccf_flux_data = np.nan

        v, v_err, R = specutils.ccf.cross_correlate(
            obs_dispersion, ccf_flux_data, model_wavelengths[si:ei],
            model_intensities[:, si:ei], rebin="observed", rebin_method="fast",
            ccf_method="fast", continuum_degree=continuum_degree[ccd])

        # Limit the available solutions.
        if isinstance(vrad_limits, (tuple, list)) \
        and vrad_limits[0] is not None and vrad_limits[1] is not None:
            R[~((vrad_limits[1] > v) * (v > vrad_limits[0]))] = np.nan

        # Get the best value by the CCF peak.
        index = np.nanargmax(R)
        ccf_results[filename] = (v[index], v_err[index], R[index], index)

    best_model_indices = [each[-1] for each in ccf_results.values()]
    v_median = np.nanmedian([each[0] for each in ccf_results.values()])
    v_err_median = np.nanmedian([each[1] for each in ccf_results.values()])
    chi_sqs = np.zeros(len(best_model_indices))
    dofs = np.zeros(len(best_model_indices))

    all_coefficients = []
    for i, index in enumerate(best_model_indices):

        # For each arm, calculate the best coefficients.
        index_coefficients = []
        chi_sq, dof = 0, -len(model_grid.dtype.names)
        for j, image in enumerate(images):

            disp = image[0].header["CRVAL1"] + image[0].header["CDELT1"] \
                * (np.arange(image[0].data.size) - image[0].header["CRPIX1"])
            ccd = image[0].header["CCD"]
            si, ei = map(sum, (pixels[:ccd - 1], pixels[:ccd]))

            closest_intensities = np.interp(disp,
                model_wavelengths[si:ei] * (1. + v_median/c),
                model_intensities[index, si:ei], left=1, right=1)
            continuum = image[0].data/closest_intensities

            # Apply chi-sq mask prior to continuum determination.
            for start, end in chi_sq_mask:
                sj, ej = np.searchsorted(disp, [start, end])
                continuum[sj:ej] = np.nan

            finite = np.isfinite(continuum)
            coefficients = np.polyfit(disp[finite], continuum[finite],
                deg=continuum_degree[ccd], w=image[1].data[finite])
            closest_flux = closest_intensities * np.polyval(coefficients, disp)

            # Apply chi-sq mask to the closest fluxes.
            for start, end in chi_sq_mask:
                sj, ej = np.searchsorted(disp, [start, end])
                closest_flux[sj:ej] = np.nan

            chi = (closest_flux - image[0].data)/image[1].data
            finite = np.isfinite(chi)

            # Update the chi^2 and finite values.
            chi_sq += np.sum(chi[finite]**2)
            dof += finite.sum()

            index_coefficients.append(coefficients)

        all_coefficients.append(index_coefficients)
        chi_sqs[i], dofs[i] = chi_sq, dof

    index = np.nanargmin(chi_sqs/dofs)
    best_model = model_grid[best_model_indices[index]]
    coefficients = all_coefficients[index]

    # Update the images with the CCF information.
    for image in images:
        image[4].header["VRAD"] = v_median
        image[4].header["U_VRAD"] = v_err_median
        image[4].header["TEFF"] = best_model[0]
        image[4].header["LOGG"] = best_model[1]
        image[4].header["FE_H"] = best_model[2]

        # Update the checksum.
        image[4].add_checksum()

    # Update the images with the normalised spectra.
    for coeffs, image in zip(coefficients, images):
        disp = image[0].header["CRVAL1"] + image[0].header["CDELT1"] \
            * (np.arange(image[0].data.size) - image[0].header["CRPIX1"])
        continuum = np.polyval(coeffs, disp)

        image[2].data = image[0].data/continuum
        image[3].data = image[1].data/continuum

        # Update checksums.
        image[2].add_checksum()
        image[3].add_checksum()

    # Create an image.
    fig, axes = plt.subplots(4)
    fig.subplots_adjust()
    model_hdus = []
    for ax, coeff, image in zip(axes, coefficients, images):

        disp = image[0].header["CRVAL1"] + image[0].header["CDELT1"] \
                * (np.arange(image[0].data.size) - image[0].header["CRPIX1"])
        continuum = np.polyval(coeff, disp)

        ax.fill_between(disp, (image[0].data - image[1].data)/continuum,
            (image[0].data + image[1].data)/continuum, facecolor="#BBBBBB",
            edgecolor="#BBBBBB")
        ax.plot(disp, image[0].data/continuum, c='k')

        ccd = image[0].header["CCD"]
        si, ei = map(sum, (pixels[:ccd - 1], pixels[:ccd]))

        mod_intensities = model_intensities[best_model_indices[index], si:ei]
        # Apply chi-sq mask.
        for mask_region in chi_sq_mask:
            sj, ej = np.searchsorted(
                model_wavelengths[si:ei] * (1. + v_median/c), mask_region)
            mod_intensities[sj:ej] = np.nan

        resampled_mod_intensities = np.interp(disp,
            model_wavelengths[si:ei] * (1. + v_median/c),
            mod_intensities, left=np.nan, right=np.nan) * continuum
            
        if SAVE_BEST_FIT_MODELS:
            hdu = fits.PrimaryHDU if ax == axes[0] else fits.ImageHDU
            model_hdus.append(hdu(data=resampled_mod_intensities,
                header=image[0].header, do_not_scale_image_data=True))

        ax.plot(disp, resampled_mod_intensities, c='r')

        ax.set_xlim(disp[0], disp[-1])
        ax.set_ylim(0, 1.2)
        ax.set_yticks([0, 0.25, 0.50, 0.75, 1.0])

    axes[0].set_title("TEFF / LOGG / [FE/H] = {0:.0f} / {1:.2f} / {2:.2f} with"\
        " chi^2/d.o.f = {3:.1f}/{4:.0f} = {5:.1f}".format(
            best_model[0], best_model[1], best_model[2],
            chi_sqs[index], dofs[index], chi_sqs[index]/dofs[index]),
        size=10)
    axes[-1].set_xlabel("Wavelength")
    fig.tight_layout()
    figure_filename = os.path.join(output_figures_folder,
        filename_mask.replace("?", "X").replace(".fits", ".png"))
    fig.savefig(figure_filename)
    print("Created figure {}".format(figure_filename))
    plt.close("all")

    if UPDATE_FILES:
        for image, filename in zip(images, filenames):
            print("Updated {}".format(filename))
            image.writeto(filename, clobber=True)

    if SAVE_BEST_FIT_MODELS:
        # Save the best fit model of each pixels to a new file.
        [each.add_checksum() for each in model_hdus]
        model_filename = os.path.join(output_models_folder,
            filename_mask.replace("?", "X").replace(".fits", "_ccf_model.fits"))
        hdulist = fits.HDUList(model_hdus)
        hdulist.writeto(model_filename)
        print("Created file {}".format(model_filename))

    # Create a table row with heaps of relevant information.
    # FILENAME, RA, DEC, TOTALEXP, NAME, GALAH_ID, UTSTART, UTEND, UTDATE, RUN,
    # OBSNUM, FIBRE, PMRA, PMDEC, MAG, DESCR, V_HELIO, V_RAD, E_V_RAD,
    # BLUE_ARM, GREEN_ARM, RED_ARM, IR_ARM, CCF_TEFF, CCF_LOGG, CCF_FEH, CCF_R_CHI_SQ
    galah_id = image[0].header["NAME"].split("_")[1] \
        if image[0].header["NAME"].lower().startswith("galahic_") else -1

    rows.append(OrderedDict([
        ("RA", image[0].header["RA"]),
        ("DEC", image[0].header["DEC"]),
        ("NAME", image[0].header["NAME"]),
        ("GALAH_ID", galah_id),
        ("TEMP_UID", filename_mask[:-5]),
        ("FILENAME", "|".join(sorted(map(os.path.basename, filenames)))),
        ("FIBRE", image[0].header["FIBRE"]),
        ("PMRA", image[0].header["PMRA"]),
        ("PMDEC", image[0].header["PMDEC"]),
        ("MAG", image[0].header["MAG"]),
        ("DESCR", image[0].header["DESCR"]),
        ("V_HELIO", image[0].header["V_HELIO"]),
        ("V_RAD", v_median),
        ("E_V_RAD", v_err_median),
        ("BLUE_ARM", 1 in ccds),
        ("GREEN_ARM", 2 in ccds),
        ("RED_ARM", 3 in ccds),
        ("IR_ARM", 4 in ccds),
        ("CCF_TEFF", best_model[0]),
        ("CCF_LOGG", best_model[1]),
        ("CCF_FEH", best_model[2]),
        ("CCF_R_CHI_SQ", np.round(chi_sqs[index]/dofs[index], 2))
        ]))

# Create table and save to disk
table = Table(rows=rows, names=rows[0].keys())
table.write("GALAH_iDR1_2DFDR_Summary.fits")

