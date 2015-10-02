************
Installation
************

We recommend you use `Anaconda <http://continuum.io/downloads>`_ to install
the necessary requirements for Starkit to work.

Once you have anaconda installed please make a new environment with the prerequisites
for starkit in the following way::

    conda create -n starkit --file https://raw.githubusercontent.com/wkerzendorf/starkit/master/conda_requirements.txt python=2
    source activate starkit
    pip install specutils

For now until a new version of Astropy comes out that fixes these problems::

    conda install cython
    pip install git+https://github.com/astropy/astropy

Once this is installed, there are two ways to install starkit. For simple use::

    pip install git+https://github.com/wkerzendorf/starkit

For to download a full development version of starkit please do::

    git clone https://github.com/wkerzendorf/starkit
    cd starkit
    python setup.py develop
