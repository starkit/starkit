"""Tests that the operations are exported by their modules."""


def test_stellar_exports_every_operation():
    namespace = {}
    exec('from starkit.base.operations.stellar import *', namespace)
    for name in ['RotationalBroadening', 'DopplerShift', 'RadialVelocity',
                 'CCM89Extinction', 'Distance']:
        assert name in namespace


def test_spectrograph_exports_every_operation():
    namespace = {}
    exec('from starkit.base.operations.spectrograph import *', namespace)
    for name in ['InstrumentRConstant', 'InstrumentDeltaLambdaConstant',
                 'Interpolate', 'Normalize', 'NormalizeParts']:
        assert name in namespace
