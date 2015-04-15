#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

""" Standardise the GALAH iDR1 globular and open cluster calibrator dataset. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import os
from glob import glob

import convert

INPUT_PATH = "/Users/arc/research/galah/data/calibrators/2dfdr_input_data/"
OUTPUT_PATH = "/Users/arc/research/galah/data/calibrators/2dfdr_standardised/"

filenames = glob(os.path.join(INPUT_PATH, "*.fits"))
for filename in filenames:
    standardised_spectra = convert.from_2dfdr(filename)

    for spectrum in standardised_spectra:

        fibre = spectrum[0].header["FIBRE"]
        year, month, day = spectrum[0].header["UTDATE"].split(":")
        output_filename = os.path.join(OUTPUT_PATH, 
            "cal_{year}{month}{day}{image_number}{fibre:03d}_{ccd}.fits.gz".format(
                year=year, month=month, day=day, fibre=fibre,
                ccd=spectrum[0].header["CCD"],
                image_number=os.path.basename(filename)[7:10]))

        print("{filename} extracted fibre {fibre} to {output_filename}".format(
            **locals()))
        spectrum.writeto(output_filename, clobber=True)