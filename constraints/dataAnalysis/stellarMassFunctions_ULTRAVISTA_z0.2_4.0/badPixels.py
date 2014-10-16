# Output a list of coorindates of bad pixels in the UltraVISTA survey.
# Andrew Benson (28-August-2014)

import sys
import numpy
from astropy import wcs
from astropy.io import fits

# Generate a set of mangle polygons which describe the geometry of bad pixels in the ULTRAVISTA survey. These correspond to the
# "cut-out" in the lower left corner of the survey mask as shown in Fig. 1 of Muzzin et al. (2013;
# http://adsabs.harvard.edu/abs/2013ApJS..206....8M). According to A. Muzzin (private communication) these regions are defined as
# those with values less than 0.02 in the Ks-band weight map.
# Andrew Benson (03-September-2014)

# Get work directory path name.
workDirectoryName=sys.argv[1]

# Open the weight map FITS file. http://ultravista.org/release1/data/UVISTA_Ks_15_12_10_skysub_015_v1.weight.fits
hdulist = fits.open(workDirectoryName+'UVISTA_Ks_15_12_10_skysub_015_v1.weight.fits')
# Parse the WCS.
w = wcs.WCS(hdulist[0].header)

# Determine bounding region.
pixcrd = numpy.array([[5000, 6000], [15000, 15750], ], numpy.float_)
world = w.wcs_pix2world(pixcrd, 1)
print(world)

# Find pixels with weight less than 0.02.
scidata = hdulist[0].data
shape   = scidata.shape
f1=open(workDirectoryName+'UVISTA_Ks_15_12_10_bad.ply', 'w')
for j in range(5000, 15000):
    print(5000,j,15000)
    iHigh=6000
    for i in range(6000, 15750):
        if ( scidata[i, j] > 0.0 and scidata[i, j] <= 0.02 ):
            iHigh=i
    pixcrd = numpy.array([[j, 6000], [j+1, iHigh+1], [j+1, 6000]], numpy.float_)
    world = w.wcs_pix2world(pixcrd, 1)
    f1.write(str(world[2,0])+" "+str(world[0,0])+" "+str(world[0,1])+" "+str(world[1,1])+"\n")
f1.close()
# Close the weight map file.
hdulist.close
