'''
example: 

from horus import rs,hs
(x,y,v) = rs('xspectrum.fits')

'''

from astropy.io import fits
import numpy as np

def rs(file):
  d = fits.open(file)
  data = d[0].data
  print(data.shape)
  
  x = data[:,:,0]
  y = data[:,:,1]
  v = data[:,:,2]

  return(x,y,v)

def hs(file):
  d = fits.open(file)
  return d[0].header


