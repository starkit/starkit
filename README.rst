StarKit
=======

************
Installation
************

We recommend you use `Anaconda <http://continuum.io/downloads>`_ to install
the necessary requirements for Starkit to work.

Once you have anaconda installed please make a new environment with the prerequisites
for starkit in the following way. This will create an environment called `starkit`::

    conda env create --file https://raw.githubusercontent.com/starkit/starkit/master/starkit_env27.yml python=2
    source activate starkit
    pip install specutils

For now until a new version of Astropy comes out that fixes these problems::

    conda install cython
    pip install astropy==1.1.2

If you are using ipython, also install::
  
    conda install ipython 

Once this is installed, there are two ways to install starkit. For simple use::

    pip install git+https://github.com/starkit/starkit

For to download a full development version of starkit please do::

    git clone https://github.com/starkit/starkit
    cd starkit
    python setup.py develop

