'''
Function to generate the value of extinction, E(B-V) from sky coordinates in Right Ascension (RA) and Declination (Dec) from SFD Dustmap 
Must have downloaded SFD Dustmap
'''

from __future__ import print_function
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery 
import astropy.units as u

def EBV_from_Coord(rightAscension,declination):
    coords = SkyCoord(ra=rightAscension*u.degree, dec=declination*u.degree, frame='icrs')
    sfd = SFDQuery()
    ebv = sfd(coords)
    return (ebv)



'''
Testing
x = EBV_from_Coord(2,5)
print("E(B-V) = %.3f" % x)
'''
