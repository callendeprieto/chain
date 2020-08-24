'''
example: 

from horus import rs,hs
(x,y,v) = rs('xspectrum.fits')

'''

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob


def rs(file):
  d = fits.open(file)
  data = d[0].data
  print(file,data.shape)
  x = data[:,:,0]
  y = data[:,:,1]
  v = data[:,:,2]
  return(x,y,v)

def hs(file):
  d = fits.open(file)
  return(d[0].header)

def mkpng(files):
  for entry in files:
    x,y,v = rs(entry)
    hd = hs(entry)
    nx,ny = x.shape
    plt.clf()
    cmap=plt.cm.jet(np.linspace(0,0.9,nx))
    plt.figure(num=None, figsize=(8*1.5, 6), dpi=100, facecolor='w', edgecolor='k')
    for i in range(nx):
      plt.plot(y[i,:]/np.median(y[i,:])+i,color=cmap[i],linewidth=0.5)
      plt.annotate("%7.2f" % x[i,0], (-150,i+0.6)  )
      plt.annotate("%7.2f" % x[i,ny-1], (ny,i+0.6) )
    plt.title(entry+'   '+hd['OBJECT'])
    plt.xlabel('wavelength (nm)')
    plt.ylabel('aperture')
    ax = plt.gca()
    ax.xaxis.set_ticklabels([])
    plt.ylim([0,i+2.])
    plt.xlim([-200,ny+200])
    plt.savefig(entry[:-4]+'png')
  return()

if __name__ == "__main__":

  print('creating png images ...')
  files = glob.glob("[xn]*fits")
  mkpng(files)


