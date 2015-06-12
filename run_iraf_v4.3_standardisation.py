#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

""" Convert IRAF reduced images to extracted spectra in GALAH FITS format. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import multiprocessing as mp
import os
from glob import glob

from astropy.table import Table

import convert

threads = 10

input_folder = "/priv/miner3/skymap/acasey/galah/standardised_spectra/temp_iraf/"
output_path = "/priv/miner3/skymap/acasey/galah/standardised_spectra/iraf_v4.3/"
fibre_table = Table.read(os.path.join(input_folder, "dr43.csv"), format="csv")



def convert_and_save(fits_filename, ccd):

    print("\t\t\tLoading {}".format(fits_filename))

    spectra = convert.from_iraf(fits_filename, fibre_table)
    print("\t\t\tFound {} program fibres:".format(len(spectra)))

    for spectrum in spectra:

        #YYMMDD*1000 + night_fco_id*10 + ccd + "_" + galah_id
        date = int(spectrum[0].header["FILEORIG"].split("/")[4])
        name = spectrum[0].header.get("NAME", None)

        if isinstance(name, (str, unicode)) \
        and name.lower().startswith("galahic_"):
            name = name.split("_")[1]
        else:
            name = "FIB{0}".format(spectrum[0].header["FIBRE"])

        temp_id = 0
        while os.path.exists(output_path, "{0}_{1}_{2}_{3}.fits".format(date,
            temp_id, ccd, name)):
            temp_id += 1

        output_filename = os.path.join(output_path,
            "{0}_{1}_{2}_{3}.fits".format(date, temp_id, ccd, name))

        print("\t\t\t\tCreated file {}".format(output_filename))
        spectrum.writeto(output_filename)

    return None



pool = mp.Pool(threads)
folders = glob(os.path.join(input_folder, "??????/combined/fluxed/"))
for i, folder in enumerate(folders):
    print("Working on folder {}".format(folder))

    fits_filenames = glob("{0}/*.ccd_1.fits".format(folder))

    for filename in fits_filenames:
        print("Distributing {}".format(filename))
        ccd = int(filename.split("_")[-1].split(".")[0])
        pool.apply_async(convert_and_save, args=(filename, ccd))

pool.close()
pool.join()
