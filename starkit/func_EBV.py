
from __future__ import print_function
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
import astropy.units as au

def getEBV_effect(rtAs_n,dec_n):
     '''
     Generates E(B-V) i.e magnitude of Extinction effect with sky coordinates.
     parameters
     ----------
     rightAscension: angular distance of a particular point measured eastward along equator
     declination: angular distance of a point north or south of the celestial equator
     

     source: https://en.wikipedia.org/wiki/International_Celestial_Reference_System
     '''
     rt_As=rightAscension*au.degree
     dec=declination*au.degree
     sk_loc = SkyCoord(rt_As, dec, frame='icrs')
     sfd = SFDQuery()
     ebv = sfd(sk_loc)
     '''
     sdf():returns reddening in a unit that is similar to magnitudes of E(B-V).
     '''
     return (ebv)

'''
NOTE: SFD Dustmap is a must, to use getEBV_effect()
'''
