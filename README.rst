StarKit
=======

.. image:: https://dev.azure.com/starkit/starkit/_apis/build/status/starkit-CI?branchName=master
   :target: https://dev.azure.com/starkit/starkit/_build/latest?definitionId=2&branchName=master


Installation
************

We recommend you use `Anaconda <https://www.anaconda.com/distribution/>`_ to install
the necessary requirements for Starkit to work.

Once you have anaconda installed please make a new environment with the prerequisites
for starkit in the following way. This will create an environment called `starkit`::

    curl -O https://raw.githubusercontent.com/starkit/starkit/master/starkit_env3.yml
    
    # install using yml file into an environment called starkit. 
    # If you want to call it something else, or already have a starkit enivornment, 
    # you can change the -n argument
    
    conda env create --file starkit_env3.yml -n starkit
    source activate starkit

For now until a new version of Astropy comes out that fixes some problems of memory leak, you need to do::

    conda install cython
    pip install specutils

Then you can additionally install any other packages you like with::
  
    conda install <your package>

Once this is installed, there are two ways to install starkit. For simple use::

    pip install git+https://github.com/starkit/starkit

To download a full development version of starkit, please do::

    git clone https://github.com/starkit/starkit
    cd starkit
    python setup.py develop

Example publications that use StarKit
**************************************

- Do, Tuan; Kerzendorf, Wolfgang; Konopacky, Quinn; Marcinik, Joseph M.; Ghez, Andrea; Lu, Jessica R.; Morris, Mark R., `Super-solar Metallicity Stars in the Galactic Center Nuclear Star Cluster: Unusual Sc, V, and Y Abundances <https://ui.adsabs.harvard.edu/#abs/2018ApJ...855L...5D/abstract>`_
- Feldmeier-Krause, A.; Kerzendorf, W.; Neumayer, N.; Schödel, R.; Nogueras-Lara, F.; Do, T.; de Zeeuw, P. T.; Kuntschner, H., `KMOS view of the Galactic Centre - II. Metallicity distribution of late-type stars <https://ui.adsabs.harvard.edu/#abs/2017MNRAS.464..194F/abstract>`_

A few test grids can be found at https://starkit.github.io/starkit/io/available_grids.html. If you use any of these grids, please make sure that the grid from which it is created (like Phoenix grid) is also cited.
