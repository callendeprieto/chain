## HORuS chain

This an IDL/GDL pipeline for reducing the data from the HORuS spectrograph on the 10.4-m Gran Telescopio Canarias (GTC). HORuS is a high-resolution (R=25,000) wide-coverage (380-700 nm, with gaps in the red), fiber-fed, echelle spectrograph. HORuS had a previous life as the Utrech Echelle Spectrograph (UES) on the 4.2-m William Herschel Telescope in La Palma. UES was decommissioned in 1997, and later donated to the IAC, where was transformed into HORuS in 2010-2016. The instrument is available on GTC since January 2019. 

This chain is a project to produce a fully automated pipeline for reducing echelle spectra, applied first to HORuS. It performs bias removal, cosmic ray cleaning, order tracing, extraction, and wavelength calibration. The list of things to improve is not short, but slow and continuous progress is happening:

- sky subtraction: tough since the standard binning is fairly agressive and there is only one sky fiber.
- wavelength calibration: currently a 4th or 5th order polynomial is used for each order, but it would be better to do a 2D calibration for all orders at once.
- optimal extraction: so far the signal is added (but working on this as of January 2020).
- RV determination: Yeisson Osorio is starting on a separate package for this
- continuum normalization and order merging: not clear the latter is a good idea, but so far the former is achieved by simply dividing the extracted spectra by the extracted flatfield. 
- 1x1 (unbinned) data: coming soon ... but the focus lately has been on 8x2 binniing data, which is the vast majority.

This software was created by Carlos Allende Prieto and it is distributed under the MIT license. It uses a bunch of open software tools: GDL, la_cosmic implementing Pieter van Dokkum's algorithm (written by Josh Bloom), the IDL astro library, the Coyote library, and other things I'm probably forgetting.

### Installing the chain

You need to first have  IDL or GDL installed. GDL is open source and can be installed easily in linux machines: in ubuntu/debian just type

  `sudo apt-get install gnudatalanguage`

(fedora users would need replace 'apt-get' by 'yum', and macports users by 'ports')

Then you need to install the IDL astro library (see https://idlastro.gsfc.nasa.gov/), and the chain itself. For that, create a directory in your home for installing idl libraries

  `cd $HOME`

  `mkdir idl`

  `cd idl`

and download the software inside that directory using, for example, git

  `git clone https://github.com/wlandsman/idlastro astro`

  `git clone https://github.com/idl-coyote/coyote coyote`

  `git clone https://github.com/callendeprieto/chain chain`

To get the libraries in your idl/gdl path you need to create a startupfile

  `echo 'device,decompose=0' > $HOME/idl/.idl_startup`

  `echo '!PATH= "'$HOME'/idl/chain:'$HOME'/idl/astro/pro:'$HOME'/idl/coyote:'$HOME'/idl/coyote/public:"  + !PATH'  >> .idl_startup`


and set an environmental variable that points to that file (assuming you are using bash)

  `echo 'export IDL_STARTUP=$HOME/idl/.idl_startup' >> $HOME/.bashrc`

(cshell users will instead do
  `echo 'setenv IDL_STARTUP $HOME/idl/.idl_startup' >> $HOME/.cshrc`  
 )

Then source your rc file

  `source $HOME/.bashrc`

If the install has succeded you should be able to 'find' the chain code from any directory in your computer. Try (replacing 'gdl' by 'idl' for IDL users) the following, and if the answers you get are similar to those following '-->' you're all set to run the chain

  `which gdl`
   --> /usr/bin/gdl 

  `echo $IDL_STARTUP`
   --> /home/callende/idl/.idl_startup

  `cat $IDL_STARTUP`

   --> device,decompose=0

   --> !PATH= "/home/callende/idl/chain:/home/callende/idl/astro/pro:/home/callende/idl/coyote:/home/callende/idl/coyote/public:"  + !PATH


  `cd $HOME`

  `gdl` 

  GDL> `.r chain`
  --> % Compiled module: CHAIN

  GDL> `.r readfits`
  --> % Compiled module: READFITS.

  GDL> `print,!PATH`
  --> /home/callende/idl/chain:/home/callende/idl/astro/pro:/home/callende/idl/coyote:/home/callende/idl/coyote/public:/usr/share/gnudatalanguage/lib/dicom:/usr/share/gnudatalanguage/lib/envi:/usr/share/gnudatalanguage/lib


  GDL> `exit`


### Reducing the data

Try reducing your data (2x8 binned observations) by entering IDL/GDL in the directory where the all the input FITS files are and typing

  GDL> `chain`

Among the output you will find a 'logfile' describing what data have been found and what processing has been done to them, x*fits files with the extracted 2D (27 orders x 2048 wavelengths) wavelength-calibrated spectra, and n*fits files, which are identical to the x*fits but with the edges trimmed and divided by the flatfield to approximately remove the instrumental response (blaze function).


### Reading the data

Included in the chain you'll find the utility rs.pro which allows reading the x*fits or n*fits data products from the pipeline

  GDL> `rs,'x0002326040-20191011-HORS-Spectroscopy.fits',y,w=x,/pl`

you can easily read the data in python using the attached horus.py script

  `from horus import rs,hs`

  `x,y,v = rs('x0002326040-20191011-HORS-Spectroscopy.fits')`



Tenerife, April 2020
Carlos

