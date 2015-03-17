# GALAH FITS Standard
This repository describes the GALAH FITS-file format and provides tools to deal with it.

# Processing Steps

Stacked spectra from a single epoch (e.g., 3 sequential GALAH observations) of the same field will be reduced by  multiple reduction pipelines. Once reduced, they should be converted to the FITS format below so that we have 1 file per star per setup (e.g., `STAR_X_blue.fits`). At this stage the format has just two extensions:

- Flux (barycentric-corrected wavelengths)
- Sigma

After standardisation, we can stack multiple exposures by matching their header information. When data are stacked from multiple epochs, the format would have more extensions:

- Stacked flux (barycentric-corrected wavelengths)
- Stacked sigma
- Flux from epoch 1 (barycentric-corrected wavelengths)
- Sigma from epoch 1
- Flux from epoch 2 (barycentric-corrected wavelengths)
- Sigma from epoch 2

After stacking, we can have Jane Lin's code run over all the data to provide a normalised spectrum and a first pass of stellar parameters. If Jane isn't keen on doing this then I can make Oracle do it. This step will insert a few extra extensions:

- Stacked flux (rest-frame wavelengths)
- Stacked sigma
- Normalised flux (rest-frame wavelengths)
- Normalised sigma
- CCF from best-fitting synthetic template?
- Flux from epoch 1 (barycentric-corrected wavelengths)
- Sigma from epoch 1
- Flux from epoch 2 (barycentric-corrected wavelengths)
- Sigma from epoch 2

This is the data format that will go the analysis nodes.

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
