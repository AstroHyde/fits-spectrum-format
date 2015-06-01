#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

""" Convert 2dfdr reduced images to extracted spectra in GALAH FITS format. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["from_2dfdr"]

import numpy as np

from astropy.io import fits
from astropy.constants import c as speed_of_light
import motions
import utils

def get_ccd_number(image):
    crval1 = image[0].header["CRVAL1"]
    #4801, 5760, 6610, 7741
    if crval1 > 7000:
        ccd_number, ccd_name = 4, "IR"
    elif crval1 > 6000:
        ccd_number, ccd_name = 3, "RED",
    elif crval1 > 5000:
        ccd_number, ccd_name = 2, "GREEN",
    elif crval1 > 4000:
        ccd_number, ccd_name = 1, "BLUE"
    else:
        raise ValueError("CRVAL1 value is weird: {0}".format(crval1))
    return (ccd_number, ccd_name)


def read_2dfdr_extensions(image):
    """
    Parse the image extensions from a 2dfdr-reduced image.

    :param image:
        The 2dfdr-reduced image.

    :type image:
        :class:`astropy.io.fits.hdu.hdulist.HDUList`
    """

    extensions = {
        "data": 0,
    }
    extnames = map(str.strip, [hdu.header.get("EXTNAME", "") for hdu in image])
    extensions["variance"] = extnames.index("VARIANCE")
    extensions["fibres"] = extnames.index("FIBRES")
    return extensions


def from_2dfdr(reduced_filename, dummy_hdus=True):
    """
    Returns a list of `astropy.io.fits.hdu.hdulist.HDUList` objects
    (1 HDU List per 2dfdr program object).

    :param reduced_filename:
        The path of the 2dfdr-reduced combined filename.

    :type reduced_filename:
        str
    """

    image = fits.open(reduced_filename)
    ext = read_2dfdr_extensions(image)

    # Verifications
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
    for program_index in np.where(image[ext["fibres"]].data["TYPE"] == "P")[0]:

        flux = image[ext["data"]].data[program_index, :]
        variance = image[ext["variance"]].data[program_index, :]
        header = header_template.copy()

        # Add header information from the fibre table.
        fibre_header = image[ext["fibres"]].data[program_index]
        for keyword in keywords:
            if isinstance(keyword, (str, unicode)):
                header[keyword] = fibre_header[keyword]
            else:
                from_keyword, to_keyword = keyword
                header[to_keyword] = fibre_header[from_keyword]

        header["FIBRE"] = program_index + 1
        
        # The 2dfdr fibre table has the RA and DEC in radians. Who does that?!
        if "RA" in keywords and "DEC" in keywords:
            header["RA"] *= 180./np.pi
            header["DEC"] *= 180./np.pi

        # Add associated comments for the fibre table information
        for keyword, comment in keyword_comments.items():
            header.comments[keyword] = comment

        # Create the HDUList.
        hdu_flux = fits.PrimaryHDU(data=flux, header=header,
            do_not_scale_image_data=True)
        hdu_flux.header["EXTNAME"] = "input_spectrum"
        hdu_flux.header.comments["EXTNAME"] = "Spectrum flux"
        # For some reason 2DFDR says 'CUNIT1' references units for axis 2
        hdu_flux.header.comments["CUNIT1"] = "Units for axis 1"

        hdu_sigma = fits.ImageHDU(data=variance**0.5, header=None,
            do_not_scale_image_data=True)
        
        for key in ("CRVAL1", "CDELT1", "CRPIX1", "CTYPE1", "CUNIT1"):
            hdu_sigma.header[key] = hdu_flux.header[key]
            hdu_sigma.header.comments[key] = hdu_flux.header.comments[key]
        hdu_sigma.header["EXTNAME"] = "input_sigma"
        hdu_sigma.header.comments["EXTNAME"] = "Flux sigma"

        # Calculate the distance to the moon.
        hdu_flux.header["MOON_DEG"] = motions.moon_distance(hdu_flux.header)
        hdu_flux.header.comments["MOON_DEG"] \
            = "Distance to the moon during observing (deg)"

        # Calculate the barycentric motion.
        v_bary, v_helio = motions.sol_corrections(header)

        hdu_flux.header["V_BARY"] = np.round(v_bary.to("km/s").value, 2)
        hdu_flux.header["V_HELIO"] = np.round(v_helio.to("km/s").value, 2)
        hdu_flux.header.comments["V_BARY"] = "Barycentric motion (km/s)"
        hdu_flux.header.comments["V_HELIO"] = "Heliocentric motion (km/s)"

        # Just to be consistent with the IRAF pipeline
        hdu_flux.header["WG6_VHEL"] = True
        hdu_flux.header.comments["WG6_VHEL"] \
            = "Was the V_HELIO correction calculated by WG6?"

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

    image.close()

    return extracted_sources



