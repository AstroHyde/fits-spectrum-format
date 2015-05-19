#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

""" Utilities for converting data into the GALAH standard FITS format. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["verify_hermes_origin", "dummy_normalisation_hdus", "dummy_ccf_hdu",
    "dummy_no_sky_hdu", "WG6_GIT_HASH"]

from subprocess import check_output
from astropy.io import fits

WG6_GIT_HASH = check_output("git rev-parse --short HEAD".split()).strip()


def verify_hermes_origin(image):
    assert image[0].header["INSTRUME"].strip() == "HERMES-2dF"


def dummy_normalisation_hdus(hdulist=None, flux_data=None, sigma_data=None):
    """
    Return dummy HDUs for the normalised flux and associated uncertainty.
    """

    # Return a dummy normalisation HDU
    hdu_normed_flux = fits.ImageHDU(data=flux_data, header=None,
        do_not_scale_image_data=True)

    # Add header information
    if hdulist is not None:
        for k in ("CRVAL1", "CDELT1", "CRPIX1", "CTYPE1", "CUNIT1"):
            hdu_normed_flux.header[k] = hdulist[0].header[k]
            hdu_normed_flux.header.comments[k] = hdulist[0].header.comments[k]

    hdu_normed_flux.header["EXTNAME"] = "normalised_spectrum"
    hdu_normed_flux.header.comments["EXTNAME"] = "Normalised spectrum flux"
    hdu_normed_flux.add_checksum()

    hdu_normed_sigma = fits.ImageHDU(data=sigma_data, header=None,
        do_not_scale_image_data=True)

    # Add header information
    if hdulist is not None:
        for k in ("CRVAL1", "CDELT1", "CRPIX1", "CTYPE1", "CUNIT1"):
            hdu_normed_sigma.header[k] = hdulist[0].header[k]
            hdu_normed_sigma.header.comments[k] = hdulist[0].header.comments[k]

    hdu_normed_sigma.header["EXTNAME"] = "normalised_sigma"
    hdu_normed_sigma.header.comments["EXTNAME"] = "Normalised flux sigma"
    hdu_normed_sigma.add_checksum()

    return [hdu_normed_flux, hdu_normed_sigma]


def dummy_ccf_hdu(hdulist=None, data=None):
    """
    Return a dummy HDU for the best-fitting cross-correlation function.
    """

    hdu_ccf = fits.ImageHDU(data=data, header=None,
        do_not_scale_image_data=True)

    hdu_ccf.header["EXTNAME"] = "CCF"
    hdu_ccf.header.comments["EXTNAME"] = "CCF from best-fitting template"

    # Add empty default values.
    default_headers = [
        ("VRAD", "Radial velocity (km/s)"),
        ("U_VRAD", "Uncertainty on radial velocity (km/s)"),
        ("TEFF", "Effective temperature"),
        ("LOGG", "Surface gravity"),
        ("FE_H", "Metallicity ([Fe/H])"),
        ("ALPHA_FE", "Alpha-enhancement ([alpha/Fe])")
    ]
    for key, comment in default_headers:
        hdu_ccf.header[key] = "NaN"
        hdu_ccf.header.comments[key] = comment
    
    hdu_ccf.add_checksum()

    return hdu_ccf


def dummy_no_sky_hdu(hdulist=None, data=None):
    """
    Return a dummy HDU for the spectrum before sky subtraction.
    """

    hdu_sky = fits.ImageHDU(data=data, header=None,
        do_not_scale_image_data=True)
    hdu_sky.header["EXTNAME"] = "no_sky"
    hdu_sky.header.comments["EXTNAME"] = "Spectrum before sky subtraction"
    
    hdu_sky.add_checksum()
    return hdu_sky