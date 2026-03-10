#!/usr/bin/env python

from pathlib import Path
from setuptools import setup, find_packages

ROOT = Path(__file__).resolve().parent
PACKAGE_NAME = "starkit"

version_ns = {}
exec((ROOT / "starkit" / "version.py").read_text(encoding="utf-8"), version_ns)

scripts = []
scripts_dir = ROOT / "scripts"
if scripts_dir.exists():
    scripts = [
        str(path)
        for path in scripts_dir.iterdir()
        if path.is_file() and path.name != "README.rst"
    ]

package_data = {
    PACKAGE_NAME: ["data/*"],
}

# Preserve the old behavior of shipping any generated .c files under starkit/
c_files = []
package_root = ROOT / PACKAGE_NAME
if package_root.exists():
    for path in package_root.rglob("*.c"):
        c_files.append(str(path.relative_to(package_root)))
package_data[PACKAGE_NAME].extend(c_files)

setup(
    name="starkit",
    version=version_ns["version"],
    description="Collection of tools for working with stellar data",
    long_description="Collection of tools for working with stellar data",
    author="Wolfgang Kerzendorf, Tuan Do",
    author_email="wkerzendorf@gmail.com",
    license="BSD",
    url="https://github.com/thorsbro/starkit",
    packages=find_packages(),
    include_package_data=True,
    package_data=package_data,
    install_requires=["astropy"],
    scripts=scripts,
    zip_safe=False,
)
