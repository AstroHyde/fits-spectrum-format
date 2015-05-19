# GALAH FITS Standard
This repository describes the GALAH FITS-file format and provides tools to deal with it.

# Processing Steps

Stacked spectra from a single epoch (e.g., 3 sequential GALAH observations) of the same field will be reduced by  multiple reduction pipelines. Once reduced, they should be converted to the FITS format below so that we have 1 file per star per setup (e.g., `STAR_X_blue.fits`). At this stage the FITS files will have the following extensions:

- Flux (heliocentric-corrected wavelengths)
- Sigma
- Normalised flux (heliocentric-corrected wavelengths)
- Normalised sigma
- Cross-Correlation Function (CCF) from best-fitting template, including associated information
- Flux before sky subtraction

So **if you want flux, you should always take the data from the first extension** (and associated uncertainties from the second extension). Similarly **if you want normalised flux you should always take the data from the third extension**. This is the data format that will go the analysis nodes.

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
 <astropy.io.fits.hdu.image.ImageHDU at 0x103e5b610>]
````

Here are the names for each extension (header keyword `EXTNAME`):

- `input_spectrum`: The observed spectrum flux
- `input_sigma`: Sigma on the observed spectrum flux
- `normalised_spectrum`: Normalised spectrum flux
- `normalised_sigma`: Sigma on normalised spectrum flux
- `CCF`: Cross-correlation function from best-fitting template

Note that we have `normalised_flux`, `normalised_sigma`, and `CCF` extensions, but because we haven't normalised or cross-correlated the spectrum yet, there are currently no data in those extensions. You can check to see if there is any data in an extension by checking the `DATASUM` header keyword. When `DATASUM` is zero, there are no data for that extension. We place dummy extensions here so that the analysis groups can be sure they are always referencing the correct data extension. (This is a lesson learned from the Gaia-ESO Survey inserting extensions over time, ruffling feathers with the analysis groups)

Once ``GUESS`` (Lin & Ireland) or ``oracle`` (Casey) has run over the standardised FITS files, these `normalised_flux`, `normalised_sigma`, and `CCF` extensions will have data present.
