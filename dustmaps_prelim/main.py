from __future__ import print_function
from astropy.coordinates import SkyCoord
from astropy import units as u
from dustmaps.sfd import SFDQuery

from dustmaps.config import config
config['data_dir'] = '.'			#location where sfd maps has been stored after downloading


 # ---- taking input from the user -----------------------

ra=float(input("enter Right Ascension (in degrees) - "))
dec=float(input("enter Declination (in degrees) - "))
dis=float(input("enter Distance - "))

#---------------------------------------------------------

#here, we are using 'ICRS' representation of coordinates

coords = SkyCoord(ra*u.deg, dec*u.deg, dis*u.pc, frame='icrs')			# assembling all the inputs together

#print(coords)				//to see the current coordinates

sfd = SFDQuery()
ebv = sfd(coords)
print()
print('E(B-V) = {:.3f} mag'.format(ebv))
print()