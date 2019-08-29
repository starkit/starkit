from starkit.base.operations import (DoubleSpectrum, SpectrographOperationModel,
                                     StellarOperationModel)

from starkit.base.operations.spectrograph import (Interpolate, Normalize,
                                                  NormalizeParts)
from starkit.base.operations.imager import Photometry
from starkit.base.operations.stellar import CCM89Extinction
from starkit.fitkit.likelihoods import SpectralChi2Likelihood


def assemble_model(spectral_grid=None, spectrum=None, normalize_parts=None,
                   normalize_npol=None, filter_set=None, mag_type='vega',
                   **kwargs):
    """

    Parameters
    ----------

    spectral_grid: ~specgrid.SpectralGrid
        spectral grid to be used in observation
    spectrum: ~specutils.Spectrum1D
        spectrum to be used for interpolation, if None neither interpolation nor
            will be performed [default None]

    normalize_npol: int
        degree of polynomial to be used for interpolation, only if not None and
        spectrum is not None will the normalization plugin be used [default None]

    plugin_names: ~li st of ~str
        select between the following available plugin choices:
        {stellar_operations}
        {instrument_operations}

    Returns
    -------
        : ~Model
        Model with the requested operations


    """

    parameters = kwargs.copy()

    def assemble_model_part(operations):
        observation_model = None
        for operation in operations:
            if ((filter_set is not None)
                and (operation.__name__.startswith('CCM89'))):
                continue

            param_values = {}
            for param_name in operation.param_names:
                if param_name in parameters:
                    param_values[param_name] = parameters.pop(param_name)

            if param_values != {}:
                if getattr(operation, 'requires_observed_spectrum', False):
                    param_values['observed'] = spectrum

                if hasattr(operation, 'from_grid'):
                    current_stellar_operation = operation.from_grid(
                            spectral_grid, **param_values)

                else:
                    current_stellar_operation = operation(**param_values)
                if observation_model is None:
                    observation_model = current_stellar_operation
                else:
                    observation_model = (observation_model |
                                         current_stellar_operation)

        return observation_model

    stellar_operations = assemble_model_part(
        StellarOperationModel.__subclasses__())
    spectrograph_operations = assemble_model_part(
        SpectrographOperationModel.__subclasses__())

    if spectrum is not None:
        if spectrograph_operations is None:
            spectrograph_operations = Interpolate(spectrum)
        else:
            spectrograph_operations = (spectrograph_operations |
                                       Interpolate(spectrum))
        if normalize_npol is not None:
            if normalize_parts:
                spectrograph_operations = (spectrograph_operations |
                                       NormalizeParts(spectrum,
                                                      normalize_parts,
                                                      normalize_npol))
            else:
                spectrograph_operations = (spectrograph_operations |
                                       Normalize(spectrum, normalize_npol))

    if filter_set is not None:
        if 'a_v' in kwargs:
            imager_operations = (CCM89Extinction(a_v=parameters.pop('a_v'))
                                 | Photometry.from_grid(
                spectral_grid, filter_set=filter_set, mag_type=mag_type))
        else:
            imager_operations = Photometry.from_grid(
                spectral_grid, filter_set=filter_set, mag_type=mag_type)
    else:
        imager_operations = None

    if parameters != {}:
        raise ValueError('Given parameters {0} not understood - please do not '
                         'add grid parameters as these are already given by '
                         'the grid'.format(
            ', '.join(list(parameters.keys()))))


    if imager_operations is not None and spectrograph_operations is not None:
        if spectral_grid is not None:
            starkit_model = spectral_grid | stellar_operations
        else:
            starkit_model = stellar_operations
        starkit_model = starkit_model | DoubleSpectrum()
        starkit_model = starkit_model | (spectrograph_operations &
                                         imager_operations)
    elif imager_operations is not None:
        if spectral_grid is not None:
            starkit_model = (spectral_grid | stellar_operations |
                             imager_operations)
        else:
            starkit_model = (stellar_operations | imager_operations)

    elif spectrograph_operations is not None:
        if spectral_grid is not None:
            starkit_model = (spectral_grid | stellar_operations
                             | spectrograph_operations)
        else:
            starkit_model = (stellar_operations | spectrograph_operations)

    else:
        if spectral_grid is not None:
            starkit_model = spectral_grid | stellar_operations
        else:
            starkit_model = stellar_operations

    return starkit_model



    ObservationModel.__class__.fit_parameters = property(fit_parameters_property)
    ObservationModel.rename('ObservationModel')

    return ObservationModel


def assemble_likelihood(spectral_grid=None, spectrum=None, normalize_parts=None,
                   normalize_npol=None, filter_set=None, mag_type='vega',
                   **kwargs):
    pass