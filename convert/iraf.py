#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

""" Convert IRAF reduced images to extracted spectra in GALAH FITS format. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["from_iraf"]

import os
import numpy as np

from astropy.io import fits
from astropy.constants import c as speed_of_light

import motions
import utils

def from_iraf(reduced_filename, fibre_table, dummy_hdus=True,
    helio_tolerance=0.1):
    """
    Returns a list of `astropy.io.fits.hdu.hdulist.HDUList` objects
    (1 HDU List per IRAF program object).

    :param reduced_filename:
        The path of the IRAF-reduced combined filename.

    :type reduced_filename:
        str

    :param fibre_table:
        A table containing all of the fibre information from all reduced files.

    :type fibre_table:
        :class:`astropy.table.Table`

    :param dummy_hdus: [optional]
        Add dummy HDUs (normalisation, normalised sigma, no sky) to the HDU list

    :type dummy_hdus:
        bool

    :param helio_tolerance: [optional]
        Acceptable absolute level of tolerance (in km/s) between Janez Koz's
        calculations and those made in this script.

    :type helio_tolerance:
        float
    """

    image = fits.open(reduced_filename)
    utils.verify_hermes_origin(image)

    # Header modifications.
    keywords = ["NAME", "RA", "DEC", "PMRA", "PMDEC", "MAG", "DESCR", "FIBRE"]
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
    header_template = image[0].header.copy()
    delete_past = [_[0] for _ in header_template.cards].index("WCSDIM")
    del header_template[delete_past:]

    # 2DFDR puts in a header called 'TDFDRVER'
    # Let's use IRAFVER
    header_template["IRAFVER"] = 0.0 # This gets updated from the fibre table.
    header_template.comments["IRAFVER"] = "IRAF reduction pipeline version"

    ccd_number, ccd_name = get_ccd_number(image)
    header_template["CCD"] = ccd_number
    header_template.comments["CCD"] = "{0} camera".format(ccd_name)
    header_template["WG6_HASH"] = utils.WG6_GIT_HASH
    header_template.comments["WG6_HASH"] = "WG6 standardisation commit hash"

    bad_pixel_filename = os.path.join(os.path.dirname(reduced_filename),
        "../../cosmicmask/badpix{}".format(os.path.basename(
            image[0].header["FILEORIG"])).replace(".fits", ".ms.fits"))

    bad_pixel_image = fits.open(bad_pixel_filename)

    extracted_sources = []
    for apid in [_ for _ in image[0].header.keys() if _.startswith("APID")]:

        fibre_name = image[0].header[apid]
        # Skip sky or dead fibres
        if fibre_name.startswith("FIBRE NOT IN USE") \
        or fibre_name.startswith("skyv1_"):
            continue

        # APID1 == index 0
        apid = int(apid[4:].strip())
        program_index = apid - 1

        flux = image[0].data[program_index, :]
        header = header_template.copy()

        # Add header information from the fibre table.
        matched_fibres \
            = (fibre_table["combined"] == os.path.basename(reduced_filename)) \
            * (fibre_table["apt_num"] == apid)

        if not np.any(matched_fibres):
            raise ValueError("could not find aperture number {0} in {1} in the"\
                " fibre table".format(
                    apid, os.path.basename(reduced_filename)))

        # Since the RA/DEC does not change, we just need the first match/
        row = fibre_table[matched_fibres][0]
        fibre_header = {
            "NAME": row["obj_name"],
            "DESCR": "", # No description passed.
            "MAG": row["mag"],
            "RA": row["ra"] * 180./np.pi,
            "DEC": row["de"] * 180./np.pi,
            "PMRA": "NaN",  # IRAF v4.1 does not output PMRA to the fibre table
            "PMDEC": "NaN", # IRAF v4.1 does not output PMDEC to the fibre table
            "FIBRE": row["fib_num"]
        }
        # Update the IRAFVER with the information from the row.
        header["IRAFVER"] = row["version"]

        # Input the RA, DEC and other keywords from the database.
        for keyword in keywords:
            header[keyword] = fibre_header[keyword]

        # Add associated comments for the fibre table information
        for keyword, comment in keyword_comments.items():
            header.comments[keyword] = comment

        # Create the HDUList.
        hdu_flux = fits.PrimaryHDU(data=flux, header=header,
            do_not_scale_image_data=True)
        hdu_flux.header["CRVAL1"] = float(image[0].header["WS_{}".format(apid)])
        hdu_flux.header.comments["CRVAL1"] = "Co-ordinate value of axis 1"
        hdu_flux.header["CDELT1"] = float(image[0].header["WD_{}".format(apid)])
        hdu_flux.header.comments["CDELT1"] = "Co-ordinate increment along axis 1"

        # 2DFDR (and the subsequent standard format) uses 2048 as their CRPIX1
        # value, but I have already processed the 2DFDR stuff, so we will just
        # 'homogenise' this a little bit in the next processing version.
        hdu_flux.header["CRPIX1"] = 0.0
        hdu_flux.header.comments["CRPIX1"] = "Reference pixel along axis 1"
        hdu_flux.header["CTYPE1"] = "Wavelength"
        hdu_flux.header.comments["CTYPE1"] = "Label for axis 1"
        hdu_flux.header["CUNIT1"] = "Angstroms"
        hdu_flux.header.comments["CUNIT1"] = "Units for axis 1"

        hdu_flux.header["EXTNAME"] = "input_spectrum"
        hdu_flux.header.comments["EXTNAME"] = "Spectrum flux"

        hdu_sigma = fits.ImageHDU(data=bad_pixel_image[0].data[program_index, :],
            header=None, do_not_scale_image_data=True)
        
        # Get the wavelength for this spectrum.
        for key in ("CRVAL1", "CDELT1", "CRPIX1", "CTYPE1", "CUNIT1"):
            hdu_sigma.header[key] = hdu_flux.header[key]
            hdu_sigma.header.comments[key] = hdu_flux.header.comments[key]
        hdu_sigma.header["EXTNAME"] = "input_sigma"
        hdu_sigma.header.comments["EXTNAME"] = "Flux sigma"

        hdu_flux.header["MOON_DEG"] = motions.moon_distance(hdu_flux.header)
        hdu_flux.header.comments["MOON_DEG"] \
            = "Distance to the moon during observing (deg)"

        # Calculate the barycentric motion.
        v_bary, v_helio = motions.sol_corrections(header)

        hdu_flux.header["V_BARY"] = np.round(v_bary.to("km/s").value, 2)
        hdu_flux.header["V_HELIO"] = np.round(v_helio.to("km/s").value, 2)
        hdu_flux.header.comments["V_BARY"] = "Barycentric motion (km/s)"
        hdu_flux.header.comments["V_HELIO"] = "Heliocentric motion (km/s)"
        hdu_flux.header["HISTORY"] \
            = "Corrected for heliocentric motion (V_HELIO)"

        # For provenance..
        hdu_flux.header["WG6_VHEL"] = True
        hdu_flux.header.comments["WG6_VHEL"] \
            = "Was the V_HELIO correction re-calculated by WG6"

        if "DOPCOR{:02d}".format(apid) in image[0].header:
            koz_dopcor = image[0].header["DOPCOR{:02d}".format(apid)]
            koz_dopcor = float(koz_dopcor.split()[0])

            # DOPCOR refers to the correction, not the measured value.
            helio_difference = koz_dopcor + v_helio.to("km/s").value
            print("Helio difference is {0:.2f} km/s (Koz - WG6) (tolerance: |"\
                "{1:.2f}|)".format(helio_difference, helio_tolerance))
            assert abs(helio_difference) < helio_tolerance

            # Since we have this value, update the headers accordingly.
            hdu_flux.header["V_HELIO"] = np.round(-koz_dopcor, 2)
            hdu_flux.header["V_HEL_RC"] = False

        # The IRAF reductions come corrected for heliocentric motion, so we do
        # not need to apply any correction.

        # Create HDUList and add checksums.
        hdulist = fits.HDUList([hdu_flux, hdu_sigma])
        [hdu.add_checksum() for hdu in hdulist]

        # Add dummy extensions
        if dummy_hdus:
            hdulist.extend(utils.dummy_normalisation_hdus(hdulist))
            hdulist.append(utils.dummy_ccf_hdu(hdulist))
            hdulist.append(utils.dummy_no_sky_hdu(hdulist))

        # Should we update the normalisation HDUs with that provided by Koz?
        # [TODO]

        extracted_sources.append(hdulist)


def get_ccd_number(image):
    lambdac = image[0].header["LAMBDAC"]
    if lambdac > 7000:
        ccd_number, ccd_name = 4, "IR"
    elif lambdac > 6000:
        ccd_number, ccd_name = 3, "RED",
    elif lambdac > 5000:
        ccd_number, ccd_name = 2, "GREEN",
    elif lambdac > 4000:
        ccd_number, ccd_name = 1, "BLUE"
    else:
        raise ValueError("LAMBDAC value is weird: {0}".format(lambdac))
    return (ccd_number, ccd_name)