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

Here is an example of what one spectrum looks like:
````python
In [5]: spectrum
Out[5]: 
[<astropy.io.fits.hdu.image.PrimaryHDU at 0x7fb256c05a10>,
 <astropy.io.fits.hdu.image.ImageHDU at 0x7fb256c05ad0>]
````

Here are some example headers from an extracted standardised spectra (first extension):
````
SIMPLE  =                    T / conforms to FITS standard                      
BITPIX  =                  -32 / array data type                                
NAXIS   =                    1 / number of array dimensions                     
NAXIS1  =                 4096                                                  
EXTEND  =                    T                                                  
COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H 
LONGSTRN= 'OGIP 1.0'           / The HEASARC Long String Convention may be used.
COMMENT   This FITS file may contain long string keyword values that are        
COMMENT   continued over multiple keywords.  The HEASARC convention uses the &  
COMMENT   character at the end of each substring which is then continued        
COMMENT   on the next keyword which has the name CONTINUE.                      
DCT_DATE= 'Oct  9 2013'        / DCT release date                               
DCT_VER = 'r3_138  '           / DCT version number                             
DETECXE =                 4096 / Last column of detector                        
DETECXS =                    1 / First column of detector                       
DETECYE =                 4112 / Last row of detector                           
DETECYS =                    1 / First row of detector                          
FIRMVSYS= 'System: AAO2 CCD Controller System Controller V1.28 130312' / System 
FIRMVSEQ= 'Sequencer: AAO2 CCD Controller Sequencer V1.3 110321' / Sequencer fir
DETECTOR= 'E2V6-G  '           / Detector name                                  
XPIXSIZE=                  15. / X Pixel size in microns                        
YPIXSIZE=                  15. / Y Pixel size in microns                        
CONFIGID=                   14 / Controller configuration Id                    
DETECTID=                   14 / Controller detector Id                         
ABCKPLID=                    6 / Analog backplane Id                            
VIDPBID =                    7 / Video personality board Id                     
CLKPBID =                   21 / Clock personality board Id                     
BRDID_1 =                65537 / Controller board #1 Id                         
BRDSN_1 =                  124 / Controller board #1 serial #                   
BRDID_2 =               327681 / Controller board #2 Id                         
BRDSN_2 =                  125 / Controller board #2 serial #                   
BRDID_3 =               131074 / Controller board #3 Id                         
BRDSN_3 =                  113 / Controller board #3 serial #                   
BRDID_4 =               131074 / Controller board #4 Id                         
BRDSN_4 =                  112 / Controller board #4 serial #                   
BRDID_5 =               131074 / Controller board #5 Id                         
BRDSN_5 =                  111 / Controller board #5 serial #                   
BRDID_6 =               131074 / Controller board #6 Id                         
BRDSN_6 =                  108 / Controller board #6 serial #                   
BRDID_7 =               196610 / Controller board #7 Id                         
BRDSN_7 =                   99 / Controller board #7 serial #                   
BRDID_8 =               262146 / Controller board #8 Id                         
BRDSN_8 =                  127 / Controller board #8 serial #                   
BRDID_9 =               458754 / Controller board #9 Id                         
BRDSN_9 =                  106 / Controller board #9 serial #                   
BRDID_10=               393217 / Controller board #10 Id                        
BRDSN_10=                  128 / Controller board #10 serial #                  
METHOD  = 'Normal ccd control method' / Observing method                        
SPEED   = 'FAST    '           / Readout speed                                  
READAMP = 'LEFT_TWO'           / Readout amplifier                              
EXPOSED =                 790. / Exposure time (seconds)                        
ELAPSED =               1096.5 / Elapsed time (seconds)                         
TOTALEXP=               1886.4 / Total exposure (seconds)                       
RO_GAIN =                 3.01 / Readout amplifier (inverse) gain (e-/ADU)      
RO_NOISE=                 4.36 / Readout noise (electrons)                      
DETTEMP =             167.5141 / Detector temperature (K)                       
HTRVOLTS=             10.57418 / Detector heater voltage                        
SPECTID = 'GN      '           / Spectrograph ID                                
SFILT1  = '        '           / Lower filter wheel position                    
SFILT2  = '        '           / Upper filter wheel position                    
ORIGIN  = 'AAO     '           / Originating Institution                        
TELESCOP= 'Anglo-Australian Telescope' / Telescope Name                         
ALT_OBS =                 1164 / Altitude of observatory in metres              
LAT_OBS =            -31.27704 / Observatory latitude in degrees                
LONG_OBS=             149.0661 / Observatory longitude in degrees               
RCT_VER = 'r3_71HB '           / Run Control Task version number                
RCT_DATE= '19-Oct-2013'        / Run Control Task version date                  
OBJECT  = 'LRa03_1P_2 CoRoT HERMES' / Object name                               
RUNCMD  = 'RUN     '           / Run command                                    
RADECSYS= 'FK5     '           / FK5 reference system                           
EQUINOX =                2000. / J2000 equinox                                  
INSTRUME= 'HERMES-2dF'         / Instrument in use                              
TDFCTVER= 'r11_20  '           / 2dF Control Task Version                       
TDFCTDAT= '06-Jan-2014'        / 2dF Control Task Version Date                  
AAOMOSFT=                    1 / Indicates file generated by AAOmega software.  
NDFCLASS= 'MFOBJECT'           / Data Reduction class name (NDFCLASS)           
CFG_FILE= 'LRa03_1P_2_FINAL_upd_p1.sds' / 2dF fibre configuration filename      
PARKMASK=                    F / 2dF park position mask NOT installed           
TDFCTEVT= '10      '           / 2dF Control Task Event                         
DUMMY   =                    0                                                  
HARTMANA= 'Open    '           / Hartmann shutter A position                    
HARTMANB= 'Open    '           / Hartmann shutter B position                    
DCHROIC1= 'GA Blue '           / Dichroic1 name                                 
DCHROIC2= 'GA Green'           / Dichroic2 name                                 
DCHROIC3= 'GA Red  '           / Dichroic3 name                                 
ELECTEMP=                  25. / Computer or electronics rack temperature in deg
SPEC1TMP=                 19.2 / Spectrograph room 1 temperature in degrees     
SPEC2TMP=                 19.3 / Spectrograph room 2 temperature in degrees     
FRM1TEMP=                 18.5 / Spectrograph table or frame 1 temperature in de
FRM2TEMP=                 18.8 / Spectrograph frame 2 temperature in degrees    
SLITTEMP=                 19.3 / Slit temperature in degrees                    
SLITMASK= 'OUT     '           / Indicates the state of the slit mask           
HPOSCONF= 'YES     '           / TRUST value of the spectrograph                
SLITDITH= 'OFF     '           / Is slit dithering enabled                      
SOURCE  = 'Plate 1 '                                                            
FLAPS   = 'OPEN    '           / 2df Flap Status                                
RUN     =                   15 / Run number                                     
OBSNUM  =                   15 / Observation number                             
GRPNUM  =                   15 / Group Number                                   
GRPMEM  =                    1 / Group member                                   
GRPMAX  =                    0 / Group maximum                                  
OBSTYPE = 'OBJECT  '           / Observation type                               
UTDATE  = '2014:01:19'         / UT date                                        
EPOCH   =     2014.05062235796 / Current Epoch, Years A.D.                      
UTSTART = '11:45:20'           / UT start                                       
UTEND   = '12:03:37'           / UT end                                         
STSTART = '05:36:47'           / ST start                                       
STEND   = '05:55:07'           / ST end                                         
UTMJD   =     56676.4898162442 / Modified Julian Date (UTC)                     
TOPEND  = '2DF     '           / Telescope top-end                              
AXIS    = 'B       '           / Current optical axis                           
AXIS_X  = -2.8216156240575E-05 / Optical axis x (mm)                            
AXIS_Y  = -9.91443977869001E-05 / Optical axis y (mm)                           
TRACKING= 'TRACKING'           / Telescope is tracking.                         
TELFOC  =                36.83 / Telescope Focus (mm)                           
MEANRA  =     92.8225278531304 / 06 11 17.41                                    
MEANDEC =     4.92237729090512 / +04 55 20.6                                    
TEL_PA  =      89.929896985001 / Rotator position angle                         
HASTART =    -8.82219841776634 / HA at start of run                             
ZDSTART =     37.1586642310223 / ZD at start of run                             
APPRA   =     1.62345877736346 / Current apparent place position right ascension
APPDEC  =   0.0857865057989631 / Current apparent place position declination    
WINDOW  = 'HERMES_2ROs_LR'     / Observing window (file name)                   
PISTON  =     2661.97619047619 / Camera focus piston 111803 (encoder)           
TILTSPEC=     1456.27904353667 / Camera focus spectral 17554 (encoder)          
TILTSPAT=                    0 / Camera focus spatial (encoder)                 
GRATID  = '8       '           / VPH grating Id ?                               
GRATTILT=                    0 / Grating tilt to be symmetric (degree)          
GRATLPMM=                3196. / Disperser ruling (lines/mm)                    
ORDER   =                    1 / Dispersion order                               
DISPENC =                    0 / Disperser angle reading (encoder unit)         
CAMENC  =                    0 / Camera angle reading (encoder unit)            
GRATANGL=                    0 / Disperser angle reading (degree)               
GRATANGR=                    0 / Requested disperser angle reading (degree)     
CAMANGL =                    0 / Camera angle reading (degree)                  
CAMANGR =                    0 / Requested camera angle reading (degree)        
LAMBDAC =                5788. / Central wavelength in Angstrom                 
LAMBDAB =                5788. / Optimal blaze wavelength in Angstrom           
LAMBDCR =                5788. / Requested central wavelength in Angstrom       
LAMBDBR =                5788. / Requested optimal wavelength in Angstrom       
DISPERS =                0.055 / Central dispersion (Angstrom/pixel)            
SPLYTEMP=                 16.6 / GREEN camera supply line temperature in degrees
RTRNTEMP=                  4.8 / GREEN camera return line temperature in degrees
CAMRTEMP=                 19.2 / GREEN camera temperature in degrees            
HAEND   =    -4.23990715626662 / HA at end of run                               
ZDEND   =     36.4179661996652 / ZD at end of run                               
WINDOXS1=                    1 / First column of window 1                       
WINDOXE1=                 4096 / Last column of window 1                        
WINDOYS1=                    1 / First row of window 1                          
WINDOYE1=                 2056 / Last row of window 1                           
WINDOXS2=                 4100 / First column of window 2                       
WINDOXE2=                 4149 / Last column of window 2                        
WINDOYS2=                    1 / First row of window 2                          
WINDOYE2=                 2056 / Last row of window 2                           
XEFFSIZE=                  15. / Effective X pixel size in microns              
YEFFSIZE=                  15. / Effective Y pixel size in microns              
RO_GAIN1=                 2.99 / Readout amplifier (inverse) gain (e-/ADU)      
RO_NOIS1=                 4.45 / Readout noise (electrons)                      
FILEORIG= '/data/aatobs/OptDet_data/140119/ccd_2/19jan20015.fits' / The filename
TDFDRVER= '1.40    '           / 2dfdr version                                  
HISTORY Subtracted bias frame                                                   
CRVAL1  =    5760.573037003231 / Co-ordinate value of axis 1                    
CDELT1  =  0.05467527599566951 / Co-ordinate increment along axis 1             
CRPIX1  =   2.048000000000E+03 / Reference pixel along axis 1                   
CTYPE1  = 'Wavelength'         / Label for axis 1                               
CUNIT1  = 'Angstroms'          / Units for axis 2                               
HISTORY GAUSS extracting using tram map 19jan20013tlm.fits                      
HISTORY Divided by fibre flat field 19jan20013red.fits                          
SCRUNCH =                    T / Has data been scrunched?                       
GAP_VER =                  6.4 / GALAH Analysis Pipeline version                
GAP_HASH= '7ae04f7 '           / GALAH Analysis Pipeline commit hash            
FIBRE   =                    2 / Fibre number                                   
NAME    = '102375900'          / Input source name                              
RA      =    93.32441666666666 / Right Ascension (degrees)                      
DEC     =               4.8276 / Declination (degrees)                          
PMRA    =               0.0021 / Proper motion in RA (mas/yr)                   
PMDEC   =              -0.0036 / Proper motion in DEC (mas/yr)                  
MAG     =                 11.7 / Input source magnitude                         
DESCR   = 'COROT_STAR'         / Input source comment                           
EXTNAME = 'input_spectrum'     / Rebinned spectrum                              
V_BARY  =   -12.10617809601693 / Barycentric motion (km/s)                      
V_HELIO =   -12.10974532765815 / Heliocentric motion (km/s)                     
CHECKSUM= 'Q4jjS2hhQ2hhQ2hh'   / HDU checksum updated 2015-03-17T12:54:46       
DATASUM = '4270491391'         / data unit checksum updated 2015-03-17T12:54:46 
HISTORY No throughput calibration performed                                     
HISTORY Sky subtracting using sky fibres                                        
HISTORY No telluric correction performed                                        
HISTORY Corrected for barycentric motion (V_BARY)   
````

The headers are largely the same in the second extension, except for the `EXTNAME`, `CHECKSUM` and `DATASUM` properties:

````
EXTNAME = 'input_sigma'        / Rebinned sigma                                 
CHECKSUM= 'Y0rAa0o4T0o9Y0o9'   / HDU checksum updated 2015-03-17T12:54:46       
DATASUM = '2955136283'         / data unit checksum updated 2015-03-17T12:54:46 
````
