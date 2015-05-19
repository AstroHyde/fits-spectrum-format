#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

""" Convert IRAF reduced images to extracted spectra in GALAH FITS format. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["from_iraf"]

from subprocess import check_output

import numpy as np

from astropy.io import fits
from astropy.constants import c as speed_of_light
import utils

WG6_GIT_HASH = check_output("git rev-parse --short HEAD".split()).strip()


def from_iraf(reduced_filename, dummy_hdus=True):
    """
    Returns a list of `astropy.io.fits.hdu.hdulist.HDUList` objects
    (1 HDU List per IRAF program object).

    :param reduced_filename:
        The path of the IRAF-reduced combined filename.

    :type reduced_filename:
        str
    """

    image = fits.open(reduced_filename)
    utils.verify_hermes_origin(image)

    # Header modifications.
    keywords = ["NAME", "RA", "DEC", "PMRA", "PMDEC", ("MAGNITUDE", "MAG"),
        ("COMMENT", "DESCR")]
    keyword_comments = {
        "NAME": "Input source name",
        "DESCR": "Input source comment",
        "MAG": "Input source magnitude",
        "RA": "Right Ascension (degrees)",
        "DEC": "Declination (degrees)",
        "PMRA": "Proper motion in RA (mas/yr)",
        "PMDEC": "Proper motion in DEC (mas/yr)",
        "FIBRE": "Fibre number"
    }
    header_template = image[ext["data"]].header.copy()
    for column in ("CDELT2", "CRPIX2", "CRVAL2", "CTYPE2", "CUNIT2"):
        del header_template[column]

    ccd_number, ccd_name = get_ccd_number(image)
    header_template["CCD"] = ccd_number
    header_template.comments["CCD"] = "{0} camera".format(ccd_name)
    header_template["WG6_HASH"] = utils.WG6_GIT_HASH
    header_template.comments["WG6_HASH"] = "WG6 standardisation commit hash"

    extracted_sources = []
    for apid in [_ for _ in image[0].header.keys() if _.startswith("APID")]:

        # APID1 == index 0
        program_index = int(apid[4:].strip()) - 1

        flux = image[0].data[program_index, :]
        header = header_template.copy()

        # Add header information from the fibre table.
        # The program_index references from 0-391 (fibres 1-392) because the
        # guide fibres are removed. So in order to reference back to the real
        # fibres we have to add values.
        guide_fibre_indices = np.array([40, 90, 140, 190, 240, 293, 340, 390])
        _ = (program_index + 1 >= guide_fibre_indices).sum()
        fibre = (program_index + 1 + _ >= guide_fibre_indices).sum()

        header["FIBRE"] = fibre

        """
        fibre_header = image[0].data[program_index]
        for keyword in keywords:
            if isinstance(keyword, (str, unicode)):
                header[keyword] = fibre_header[keyword]
            else:
                from_keyword, to_keyword = keyword
                header[to_keyword] = fibre_header[from_keyword]

        # The 2dfdr fibre table has the RA and DEC in radians. Who does that?!
        if "RA" in keywords and "DEC" in keywords:
            header["RA"] *= 180./np.pi
            header["DEC"] *= 180./np.pi

        # Add associated comments for the fibre table information
        for keyword, comment in keyword_comments.items():
            header.comments[keyword] = comment

        """

        # Create the HDUList.
        hdu_flux = fits.PrimaryHDU(data=flux, header=header,
            do_not_scale_image_data=True)
        hdu_flux.header["EXTNAME"] = "input_spectrum"
        hdu_flux.header.comments["EXTNAME"] = "Spectrum flux"

        hdu_sigma = fits.ImageHDU(data=None, header=None,
            do_not_scale_image_data=True)
        
        # Get the wavelength for this spectrum.
        for key in ("CRVAL1", "CDELT1", "CRPIX1", "CTYPE1", "CUNIT1"):
            hdu_sigma.header[key] = hdu_flux.header[key]
            hdu_sigma.header.comments[key] = hdu_flux.header.comments[key]
        hdu_sigma.header["EXTNAME"] = "input_sigma"
        hdu_sigma.header.comments["EXTNAME"] = "Flux sigma"



        # Calculate the barycentric motion.
        v_bary, v_helio = motions.from_header(header)

        hdu_flux.header["V_BARY"] = v_bary.to("km/s").value
        hdu_flux.header["V_HELIO"] = v_helio.to("km/s").value
        hdu_flux.header.comments["V_BARY"] = "Barycentric motion (km/s)"
        hdu_flux.header.comments["V_HELIO"] = "Heliocentric motion (km/s)"
        hdu_flux.header["HISTORY"] \
            = "Corrected for heliocentric motion (V_HELIO)"

        # Add motion correction information.
        for hdu in (hdu_flux, hdu_sigma):
            hdu.header["CRVAL1"] *= 1. + (v_helio/speed_of_light).value
            hdu.header["CDELT1"] *= 1. + (v_helio/speed_of_light).value

        # Create HDUList and add checksums.
        hdulist = fits.HDUList([hdu_flux, hdu_sigma])
        [hdu.add_checksum() for hdu in hdulist]

        # Add dummy extensions
        if dummy_hdus:
            hdulist.extend(utils.dummy_normalisation_hdus(hdulist))
            hdulist.append(utils.dummy_ccf_hdu(hdulist))
            hdulist.append(utils.dummy_no_sky_hdu(hdulist))

        extracted_sources.append(hdulist)