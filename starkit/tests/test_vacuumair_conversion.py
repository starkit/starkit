from starkit.utils import vacuumair_conversion

#test wavelength
#this program is for single test case
#that's why I used wavelength variable so that in future anyone can test by assigning value to this variable
wavelength=10

def test_convert_vacuum2air():
    a=vacuumair_conversion.convert_vacuum2air(wavelength)   
    a=str(a)                                           
    assert a=="9.999444010085922 Angstrom"
def test_convert_air2vacuum():
    b=vacuumair_conversion.convert_air2vacuum(wavelength)
    b=str(b)    
    assert b=="10.0005550758127 Angstrom"
