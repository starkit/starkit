import numpy as np
from astropy.modeling import models 
from dust_extinction.parameter_averages import F99
from astropy import units as au

# -- Values taken from docs.astropy.org for examples --

# Guassian1D model
Gau_mod=models.Gaussian1D

def Getfrom_Gmodel(ampli, men, st_dev):
	amplitude=ampli
	mean=men
	stddev=st_dev
	return Gau_mod(amplitude,mean,stddev)     

G_result=Getfrom_Gmodel(1.2,0.9,0.5) 				  
print(G_result)



# Blackbody1D model
Blb_mod=models.BlackBody1D 	

def Getfrom_BBmodel(t,bf):
	temp=t                                                            # Temperature (convert to Kelvin)
	bflux=bf                                                          #Bolometric flux from the  star
	return Blb_mod(temp,bflux)

B_result=Getfrom_BBmodel(5000*au.AA)
print(B_result)
					
wavlens = np.logspace(np.log10(1000), np.log10(3e4), num=1000)*au.AA      # 1D array parameter for spectrum & wavelenghts 
spectrum = blackbody_lambda(wavlens, 10000*au.K)

# defining F99 to model the Extinction effect 
ext = F99(Rv=3.1)
spectrum_ext = spectrum*ext.extinguish(wavlens, Ebv=0.5)                    # With reference from docs.astropy.org







