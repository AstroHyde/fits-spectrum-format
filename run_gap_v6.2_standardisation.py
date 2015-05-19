#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

""" Convert 2dfdr reduced images to extracted spectra in GALAH FITS format. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import multiprocessing as mp
import os
from glob import glob

import convert

threads = 10
output_path = "/priv/miner3/skymap/acasey/galah/standardised_spectra/gap_v6.2/"
folders = sorted(glob("/priv/miner3/galah/galahsk/GAP/tdfdr_data/*/"))


def convert_and_save(fits_filename, fco, ccd):

    print("\t\t\tLoading {}".format(fits_filename))

    spectra = convert.from_2dfdr(fits_filename)
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

        output_filename = os.path.join(output_path, 
            "{0}{1}_{2}.fits".format(date * 1000 + fco, ccd, name))

        print("\t\t\t\tCreated file {}".format(output_filename))
        spectrum.writeto(output_filename)

    return None



pool = mp.Pool(threads)

for i, folder in enumerate(folders):
    print("Working on folder {}".format(folder))

    for ccd in range(1, 5):
        print("\tChecking CCD {}".format(ccd))
        fco_folders = sorted(glob("{0}ccd_{1}/*/".format(folder, ccd)))

        for fco, fco_folder in enumerate(fco_folders, start=1):
            if fco_folder.split("/")[-2].lower() == "bias": continue

            fits_filename = glob("{0}/*_combined.fits".format(fco_folder))
            
            print("\t\tIn FCO {0} (path {1}, size {2})".format(fco, fco_folder,
                len(fits_filename)))
            
            if len(fits_filename) == 0: continue
            elif len(fits_filename) > 1:
                raise WTFError()
            else:
                # Extract the 1D spectra.
                pool.apply_async(convert_and_save,
                    args=(fits_filename[0], fco, ccd))

pool.close()
pool.join()
