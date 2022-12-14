Installation
============

To install kosmatau3d, download the git repo, then simply install using,

.. code:: console
   
   $ pip install ./kosmatau3d

This will install kosmatau3d and most of its dependencies. There are a couple of packages that will not be installed, such as :code:`astrokit` and :code:`healpy`.
These are only used in the :code:`comparison` submodule, which is used to process observation files and compare with the synthetic observations from kosmatau3d.
You will need to install these or similar packages to use this submodule.
:code:`astrokit` is a package developed by Slawa Kabanovic during his PhD, and it can be used to work with observational data, create slices, and calculate observational error.
:code:`healpy` is the python port of the well-kown and maintained HEALPix code.
It only installs well in a Linux environment, so that is why it is not a strict dependency.
Ideally the user is in a linux environment, so downloading this dependency is not an issue.
