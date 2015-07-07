# GALAH Spectrum FITS Standard
This repository describes the GALAH spectrum FITS format and provides tools to deal with it.

# Processing Steps

Stacked spectra from a single epoch (e.g., 3 sequential GALAH observations) of the same field will be reduced by  multiple reduction pipelines. Once reduced, they should be converted to the FITS format below so that we have 1 file per star per setup (e.g., `STAR_X_blue.fits`). At this stage the FITS files will have the following extensions (`EXTNAME`-header: description):

- `input_spectrum`: Flux (heliocentric-corrected wavelengths)
- `input_sigma`: Sigma
- `normalised_spectrum`: Normalised flux (heliocentric-corrected wavelengths)
- `normalised_sigma`: Normalised sigma
- `CCF`: Cross-Correlation Function (CCF) from best-fitting template, including associated information
- `no_sky_spectrum`: Flux (heliocentric-corrected wavelengths) before sky subtraction
- `no_sky_sigma` (TBC): Sigma before sky subtraction

So **if you want flux, you should always take the data from the first extension** (and associated uncertainties from the second extension). Similarly **if you want normalised flux you should always take the data from the third extension**. This is the data format that will go the analysis nodes.

# Load Example

Below (and in `example.py`) you will find an example Python function that will load GALAH spectra. You can select to return normalised or unnormalised spectra, and select whether to correct the spectra for the measured radial velocity.

````python

import numpy as np
from astropy.io import fits

def load_GALAH(filename, normalised=False, rest=True, **kwargs):
    """
    Load a Spectrum1D object from the GALAH standardised FITS format.

    :param filename:
        The path of the filename to load.

    :type filename:
        str
        
    :param normalised: [optional]
        Return normalised spectra, if it exists.
    
    :type normalised:
        bool
    
    :param rest: [optional]
        If the radial velocity has been measured and put into the headers,
        return rest frame spectra.
    
    :type rest:
        bool
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
````

# Conversion Example

The code below will extract individual program spectra from a `2DFDR`-reduced image and write them to file.

````python

import convert

# 19jan20015red.fits is a 2dfdr-reduced combined file
standardised_spectra = convert.from_2dfdr("19jan20015red.fits")

for spectrum in standardised_spectra:
    
    # Information like FIBRE, NAME, etc are available through the header.
    spectrum.writeto("19jan20015red_{FIBRE}.fits".format(**spectrum[0].header))
````

Here is an example of what one spectrum looks like:
````python
In [2]: spectrum
Out[2]: 
[<astropy.io.fits.hdu.image.PrimaryHDU at 0x105794990>,
 <astropy.io.fits.hdu.image.ImageHDU at 0x105794a50>,
 <astropy.io.fits.hdu.image.ImageHDU at 0x1057a73d0>,
 <astropy.io.fits.hdu.image.ImageHDU at 0x1057af650>,
 <astropy.io.fits.hdu.image.ImageHDU at 0x103e5b610>,
 <astropy.io.fits.hdu.image.ImageHDU at 0x10342ed10>,
 ]
````

Here are the names for each extension (header keyword `EXTNAME`):

- `input_spectrum`: The observed spectrum flux
- `input_sigma`: Sigma on the observed spectrum flux
- `normalised_spectrum`: Normalised spectrum flux
- `normalised_sigma`: Sigma on normalised spectrum flux
- `CCF`: Cross-correlation function from best-fitting template
- `no_sky_spectrum`: The observed spectrum flux before sky subtraction
- `no_sky_sigma` (TBC): Sigma on the observed spectrum flux before sky subtraction

Note that we have `normalised_flux`, `normalised_sigma`, and `CCF` extensions, but because we haven't normalised or cross-correlated the spectrum yet, there are currently no data in those extensions. You can check to see if there is any data in an extension by checking the `DATASUM` header keyword. When `DATASUM` is zero, there are no data for that extension. 

We place dummy extensions here so that the analysis groups can be sure they are always referencing the correct data extension. This is a lesson learned from the Gaia-ESO Survey inserting extensions over time, ruffling feathers with the analysis groups. To future-proof your code we recommend you reference extensions by their `EXTNAME` keyword, but if you just reference extensions by their index you can be assured any additional extensions will only be appended;  **no re-ordering of extensions will happen between major data releases**.

Once ``GUESS`` (Lin & Ireland) or ``oracle`` (Casey) has run over the standardised FITS files, these `normalised_flux`, `normalised_sigma`, and `CCF` extensions will have data present.
