'''
Function to generate the value of extinction E(B-V) from sky coordinates in Right Ascension (RA) and Declination (Dec) from SFD Dustmap. 
Must have downloaded SFD Dustmap
'''

from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery 
import astropy.units as u

'''
PARAMETERS:
rightAcsension: (int/float) angular distance of a particular point measured eastward along equator
declination: (int/float) angular distance of a point north or south of the celestial equator
-----------
RETURNS:
ebv: (magnitudes) color excess - degree of interstellar extinction/reddening
'''
def EBV_from_Coord(rightAscension,declination):
    coords = SkyCoord(ra=rightAscension*u.degree, dec=declination*u.degree, frame='icrs')
    sfd = SFDQuery()
    ebv = sfd(coords)
    return (ebv)
