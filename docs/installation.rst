************
Installation
************

We recommend you use `Anaconda <http://continuum.io/downloads>`_ to install
the necessary requirements for Starkit to work.

Once you have anaconda installed please make a new environment with the prerequisites
for starkit in the following way::

    curl -O https://raw.githubusercontent.com/starkit/starkit/master/starkit_env27.yml
    conda create -n starkit --file starkit_env27.yml
    source activate starkit

Once this is installed, there are two ways to install starkit. For simple use::

    pip install git+https://github.com/starkit/starkit

For to download a full development version of starkit please do::

    git clone https://github.com/starkit/starkit
    cd starkit
    python setup.py develop
