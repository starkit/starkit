"""Regression test for starkit.evolutionkit.base, which raised TypeError on import."""


def test_grid_connector_inputs_are_the_stellar_parameters():
    from starkit.evolutionkit.base import GridConnector

    assert GridConnector.inputs == ("teff", "logg", "mh", "radius")
