#!/usr/bin/env python

from pathlib import Path

from setuptools import setup

ROOT = Path(__file__).resolve().parent

version_ns = {}
exec((ROOT / "starkit" / "version.py").read_text(encoding="utf-8"), version_ns)


setup(version=version_ns["version"])
