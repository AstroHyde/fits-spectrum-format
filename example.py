
from astropy.io import fits

import numpy as np


def load_GALAH(filename, normalised=False, rest=True, **kwargs):
    """
    Load a Spectrum1D object from the GALAH standardised FITS format.

    :param filename:
        The path of the filename to load.

    :type filename:
        str
    """

    image = fits.open(filename, **kwargs)

    data_ext = 0 if not normalised else 2

    flux = image[data_ext].data
    if flux.size == 0 and normalised:
        raise ValueError("no normalised spectrum found")

    variance = image[data_ext + 1].data

    disp = image[data_ext].header["CRVAL1"] \
        + (np.arange(flux.size) - image[data_ext].header.get("CRPIX1", 0)) \
        * image[data_ext].header["CDELT1"]

    header_columns = ["CCD", "WG6_HASH", "NAME", "RA", "DEC",
        "PMRA", "PMDEC", "MAG", "DESCR", "FIBRE", "MOON_DEG", "V_HELIO"]

    headers = dict(zip(header_columns, [image[0].header.get(k, None) \
        for k in header_columns]))

    if rest:
        vrad = image[4].header.get("VRAD", np.nan)
        if np.isfinite(vrad):
            disp *= (1 - vrad/299792.458)
        else:
            logger.warn("No velocity information found! Spectrum not at rest!")

    return (disp, flux, variance, headers)