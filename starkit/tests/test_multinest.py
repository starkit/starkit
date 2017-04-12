from starkit.fitkit import multinest
def test_result_obj():
    # test the results object
    r = multinest.MultiNestResult.from_multinest_basename('sun10000',['teff','logg','feh','vz'])

    # calculate the 1 sigma levels for each parameter
    sigmas = r.calculate_sigmas(1)
    print sigmas
    print r.median
    r.plot_triangle()
    r2 = multinest.MultiNestResult.from_multinest_basename('sun10000',['teff','logg','feh','vz'],equal_weights=True)
    sigmas2 = r.calculate_sigmas(1)
    print sigmas2
