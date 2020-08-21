from binaryio import *
from astropy.io import fits
import numpy as np
import scipy.integrate as integrate


c = 29979245800.

# Set cosmology
# h100 = Hubble Constant, units of 100 km/s/Mpc
# omega_m = Matter density ratio to critical density at z=0
h100 = 0.704
omega_m = 0.272


# What scale factors are you observing (need not be the same as output)
aexps = [0.41]

# Set rest wavelength in microns. Here we use Halpha
rest_wave = .65628

# Enter list of VelocityCube Filenames, without the .bin
names = [
"Velcube_x_152_l17_sfr_Ha_fd0.01",
"Velcube_y_152_l17_sfr_Ha_fd0.01",
"Velcube_z_152_l17_sfr_Ha_fd0.01",
]



# Here we go! Nothing more needed.
 

# Calculate z values for filenames
z_strs = []
for a in aexps:
  z_str = "{0:.2f}".format(1./a - 1.)
  z_strs.append(z_str)


# Loop over velcubes
for name in names:
  ifile = open(name+".bin","rb")
  # this means: read from binary 3 integers
  head1 = bread_int(ifile,3)
  ix,iy,nv = head1

  # read data
  toto_raw = bread_dbl(ifile,ix*iy*nv)
  toto = np.array(toto_raw)
  toto_raw = np.array(toto_raw)
  # a copy is made to correct velocity space to wavelength space

  toto = np.reshape(toto,head1,order='F')
  toto_raw = np.reshape(toto_raw,head1,order='F') 

  xrange = bread_dbl(ifile,2)
  yrange = bread_dbl(ifile,2)
  zrange = bread_dbl(ifile,2)

  print "Found an xrange of ", xrange
  print "Found an yrange of ", yrange
  print "Found an zrange of ", zrange

  xc = xrange[0]
  yc = yrange[0]

  unitl = bread_dbl(ifile,1)
  print "Found an unitl of ", unitl

  varray = bread_float(ifile,nv)
  print "Found an varray range of ", varray[0], varray[-1]

  # For each scale factor to observe, do unique velocity correction
  for a,z_str in zip(aexps,z_strs):
    dp = c/(1.e7*h100)*1000.*integrate.quad(lambda x: 1./np.sqrt(omega_m*x + (1.-omega_m)*x**4),a,1.)[0] # kpc
    da = dp*a
    dl = dp/a*3.08e21 # cm
    plate_scale = da/206264.80624709636
    print "At a=", a, " found a plate_scale of (kpc/'') ", plate_scale

    lam_c = rest_wave/a
    delv = (varray[-1] - varray[0])/(nv-1)
    beta = (varray[-1] - varray[0])/2/(c/1.e5)
    del_lam = (np.sqrt( (1.+beta)/(1.-beta) ) - 1.)*lam_c/((nv+1)/2)

    # v - lambda conversion
    for i in range(nv):
      iv = i - (nv-1)/2
      q = 2.*iv*del_lam/lam_c + iv**2*del_lam**2/lam_c**2
      fv = (nv-1.0)/2.0 + c/1.e5*(q/(2.+q))/delv

      im1 = int(fv)
      ip1 = im1 + 1
      fp1 = fv - im1
      fm1 = 1.0 - fp1
      if ip1 >= nv:
        ip1 = im1
        fm1 = 1.0
        fp1 = 0.0
      toto[:,:,i] = fm1*toto_raw[:,:,im1] + fp1*toto_raw[:,:,ip1]

    lam0 = lam_c - del_lam*nv/2

    # convert units
    toto = toto*(plate_scale)**2/(4.*np.pi*dl**2)
    toto = toto/del_lam

    xc0 = xc*unitl*1000./plate_scale*1000.
    yc0 = yc*unitl*1000./plate_scale*1000.

    delx = (xrange[1]-xrange[0])/(ix-1)*unitl*1000./plate_scale*1000.
    dely = (yrange[1]-yrange[0])/(iy-1)*unitl*1000./plate_scale*1000.

    print "values on ", np.min(toto), np.max(toto)

    # initialise the Fits file
    hdu = fits.PrimaryHDU(toto.transpose())

#   print "Image has pixel dimensions (mas x mas x microns) " , delx, dely, delv
#   See http://jsoc.stanford.edu/doc/keywords/STEREO/STEREO_site_standard_fits_keywords.txt for possible keywords 

    hdu.header['NAXIS1'] = ix
    hdu.header['NAXIS2'] = iy
    hdu.header['NAXIS3'] = nv

    hdu.header['CRPIX1'] = 1
    hdu.header['CRPIX2'] = 1
    hdu.header['CRPIX3'] = 1

    hdu.header['CTYPE1'] = 'x'
    hdu.header['CTYPE2'] = 'y'
    hdu.header['CTYPE3'] = 'WAVELENGTH'

    hdu.header['CUNIT1'] = 'mas'
    hdu.header['CUNIT2'] = 'mas'
    hdu.header['CUNIT3'] = 'microns'

    hdu.header['CDELT1'] = delx
    hdu.header['CDELT2'] = dely
    hdu.header['CDELT3'] = del_lam

    hdu.header['CRVAL1'] = xc0
    hdu.header['CRVAL2'] = yc0
    hdu.header['CRVAL3'] = lam0

    hdu.header['FUNITS'] = 'erg/s/cm2/um/arcsec2'
    hdu.header['SPECRES'] = del_lam

    hdu.writeto(name+'_z'+z_str+'.fits')

