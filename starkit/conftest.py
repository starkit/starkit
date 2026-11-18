"""Fixtures shared by the StarKit tests."""

from pathlib import Path

import numpy as np
import pytest

SVO_DATA = Path(__file__).parent / 'tests' / 'data'


@pytest.fixture
def svo(monkeypatch):
    """
    Serve the SVO Filter Profile Service files saved in tests/data instead of
    downloading them.

    Returns a dict that records the requested URLs ('urls') and the URLs
    removed from the download cache ('cleared'). A test can set 'error_file'
    to the file served for filters that are not saved.
    """
    from starkit import filters

    requests = {'urls': [], 'cleared': []}

    def download_file(url, cache=False, **kwargs):
        requests['urls'].append(url)
        filter_id = url.split('ID=')[1]
        path = SVO_DATA / 'svo_{0}.xml'.format(filter_id.replace('/', '_'))
        return str(path) if path.exists() else requests['error_file']

    monkeypatch.setattr(filters, 'download_file', download_file)
    monkeypatch.setattr(filters, 'clear_download_cache',
                        lambda url: requests['cleared'].append(url))
    return requests


@pytest.fixture
def flat_vega(monkeypatch):
    """Replace the Vega spectrum by f_lambda = 1 from 1000 to 30000 Angstrom."""
    from starkit import filters

    wavelength = np.linspace(1000.0, 30000.0, 29001)
    monkeypatch.setattr(filters, 'vega_spectrum',
                        lambda: (wavelength, np.ones(wavelength.size)))
