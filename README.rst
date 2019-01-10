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


Then you can install any other packages you like with
  
    conda install <your package>

Once this is installed, there are two ways to install starkit. For simple use::

    pip install git+https://github.com/starkit/starkit

For to download a full development version of starkit please do::

    git clone https://github.com/starkit/starkit
    cd starkit
    python setup.py develop


## Example publications that use StarKit:

- Do, Tuan; Kerzendorf, Wolfgang; Konopacky, Quinn; Marcinik, Joseph M.; Ghez, Andrea; Lu, Jessica R.; Morris, Mark R., [Super-solar Metallicity Stars in the Galactic Center Nuclear Star Cluster: Unusual Sc, V, and Y Abundances](https://ui.adsabs.harvard.edu/#abs/2018ApJ...855L...5D/abstract)
- Feldmeier-Krause, A.; Kerzendorf, W.; Neumayer, N.; Schödel, R.; Nogueras-Lara, F.; Do, T.; de Zeeuw, P. T.; Kuntschner, H., [KMOS view of the Galactic Centre - II. Metallicity distribution of late-type stars](https://ui.adsabs.harvard.edu/#abs/2017MNRAS.464..194F/abstract)


