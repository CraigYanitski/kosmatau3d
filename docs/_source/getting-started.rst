:html_theme.sidebar_primary.remove:

***************
Getting Started
***************

The best way to get started with all of the finer details of how to use the 
code and how it works is to review the 
`notebook <https://github.com/CraigYanitski/kosmatau3d/blob/main/notebooks/single-voxel/voxel.ipynb>`_.
Since the code is still in early development, all of the examples show the 
features of a voxel_ in the framework of this code.
The three-dimensional uses of this code have thus-far been limited to models 
of the Milky Way as they appear in Yanitski (2023).
Some of the results are shown in `Three-dimensional PDR Models`_, but potential 
users are encouraged to contact yanitski@ph1.uni-koeln.de to iterate on their 
use-case.

The Single-Voxel model
======================

.. _voxel:

.. sidebar:: Voxel representation


   .. figure:: _static/uniform_RT-small_comp.png
      :alt: voxel diagram
      :width: 500

      A diagram representing how the voxels work.
      The region in the red box is the ISM we want to model, with the clumps 
      shown as scatter points.
      The blue boxes are the voxels.
      In order for the probabilistic calculation to be accurate, the clumps are 
      randomly positioned in the voxel.

This is the base functionality of :code:`kosmatau3d`.
A *voxel*, or volumetric pixel, is one cell in a three-dimensional spatial
grid.
It contains one or more *ensembles* of KOSMA-:math:`\tau` clumps (see 
:doc:`Theory <theory>` for more information about these), which each contribute 
to the emissivity and absorption in the voxel.
It performs a probabilistic calculation of the emission by considering all 
possible combinations of clumps in the line-of-sight (see the schematic below).

Below is a minimal example of how to initialise and set the properties in a voxel.
Note that there are plenty of other parameters that can be adjusted to modify 
the physics of the ensembles contained within the voxel, and most of them are 
explained in :class:`kosmatau3d.models.voxel.Voxel` and in the jupyter
`notebook <https://github.com/CraigYanitski/kosmatau3d/blob/main/notebooks/single-voxel/voxel.ipynb>`_.

:doc:`generic <generic>`

.. sidebar::
    
    something

.. code-block:: python
    :caption:
    :linenos:

    >>> from kosmatau3d import models
    >>> vox = models.Voxel()
    >>> vox.set_properties(
    ...     voxel_size=1,
    ...     clump_mass_number=[3],
    ...     clump_mass_range=[[0,2]],
    ...     ensemble_mass=1e1,
    ...     ensemble_density=1e5,
    ...     fuv=1e5,
    ... )

This will create a voxel with side length 
:math:`\ell_\mathrm{vox}=1\,\mathrm{pc}` (so the volume is 
:math:`1\, \mathrm{pc}^3`) containing an ensemble with three distinct clump 
masses: :math:`1`, :math:`10`, and :math:`100\, M_\odot`.
The number of each of these clumps in the ensemble is calculated from the 
`clump-mass distribution` as explained in :doc:`Theory <theory>`.
Since the mass of the ensemble in this example is lower than the maximum 
mass of the clumps in the ensemble, this example will only contain a fraction 
of the largest clump.
This showcases a key property of voxels in :code:`kosmatau3d`: the 
scale-invariance of the calculations.

Now that a voxel has been initialised with an ensemble, it is possible to 
obtain intrinsic properties such as the fractional abundance by,

..  code:: python

    >>> vox.get_abundances(species=["H2", "H", "C+", "CO"])

where the abundances are normalised by the total hydrogen abundance 
(:math:`n_\mathrm{H} = n_\mathrm{H^0} + 2\, n_\mathrm{H_2}`).
Extrinsic properties such as column density are also available at this 
point, but for the emission it is necessary to first execute the computation:

..  code:: python

    >>> vox.calculate_emission()

This will compute the emissivity and absorption for the line and continuum emission
from the base KOSMA-:math:`\tau` models.
One can then access the result using,

..  code:: python

    >>> vox.get_species_emissivity()

Read the documentation of :class:`kosmatau3d.models.voxel.Voxel` to find out 
the different arguments for this method.
In the current procedure for calculating the emission, we use the peak intensity 
and optical depth as well as post-processing the HI :math:`21\, \mathrm{cm}` line 
emission by making an isothermal approximation.
for that reason we have split the accessing methods for the continuum and line 
emission into :code:`get_dust_...` and :code:`get_species_...`, respectively, 
and there is a kwarg :code:`hi` that can be set to :code:`True` to get the 
HI line emission.
The emission values that are available are,

* emissivity :math:`\epsilon_\nu` in :math:`\frac{K}{pc}`
* absorption :math:`\kappa_\nu` in :math:`\frac{1}{pc}`
* intensity :math:`I_\nu` in :math:`K`
* optical depth :math:`\tau_\nu` (dimensionless)

The intensity and optical depth require integrating over the length-scale of the 
voxel, so they should not be used in three-dimensional models.

Three-dimensional PDR Models
============================

.. sidebar:: Galactic model

   .. figure:: _static/integrated_C+1.png
      :alt: model showing integrated C+

      One Galactic model, where the voxels are coloured according to the 
      \[CII\] :math:`158\, \mu\mathrm{m}` intensity integrated over the spectrum
      (:math:`-350` to :math:`350\, \mathrm{km\, s^{-1}}`).

The most complex functionality of kosmatau3d, and the reason for its 
development since its conception 
(`Andree-Labsch et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017A%26A...598A...2A/abstract>`_). 
This has so-far been used to model the Orion Bar and the Milky Way, though more 
models will soon be developed.
The figure below depicts the \[CII\] :math:`158\, \mu\mathrm{m}` integrated 
intensity in each voxel of one of the galactic models, which is then used to 
compute the synthetic observation.

The benefit of using :code:`kosmatau3d` voxels for the model is two-fold: 
it uniquely accounts for the inhomogeneity and shadowing in PDRs and has 
a velocity component to its emission.
Because of this, it is important to first get the velocity information for each 
voxel in the model correct before integrating the radiative transfer equation.
The procedure is described in detail in Ch. IV of 
`Yanitski (2023) <https://kups.ub.uni-koeln.de/71850/>`_, and it results in a 
longitude-velocity diagram like below.
We focus on galactic latitude :math:`b\! =\! 0` to avoid the complications 
regarding partially-filled voxels.

We use the :class:`kosmatau3d.models.model.Model()` to initialise all voxels, 
compute their radiative properties, and save the relevant data in FITS files.
It shares many kwargs with :meth:`set_properties`, though some are renamed 
to fit the context of a three-dimensional model.
A minimal working example to create a galactic model with voxel size 
:math:`\ell_\mathrm{vox}=400\,\mathrm{pc}` is,

..  code:: python

    >>> from kosmatau3d import models
    >>> kwargs = { ... }
    >>> galaxy = models.Model(resolution=400, 
    ...     history_path='.', 
    ...     folder='temp', 
    ...     **kwargs)
    >>> galaxy.calculate_model()

Here :code:`kwargs` can be used to specify any of the model parameters.
A distinct difference in making the full model is that the kwargs are given 
when initialising the object instance rather than through a separate method.
While all of the physical and emissive properties are calculated at this stage, 
a synthetic observation requires the :mod:`kosmatau.radiative_transfer`:

..  code:: python
    
    >>> import numpy as np
    >>> models.radiative_transfer.calculateObservation(directory='temp/', 
    ...     slRange=[(-np.pi, np.pi), (-np.pi/2, np.pi/2)], 
    ...     nsl=[180, 90])

This will result in a synthetic datacube of the region for all of the included 
transitions (by default all of them) and a subset of the dust continuum (where 
22 wavelengths are used; enough to span the FIR emission).
From the synthetic intensity datacube, it is possible to get the 
position-velocity diagram as below.

.. figure:: _static/model_C+1.png
   :alt: model synthetic C+ 1
   :width: 500

   The synthetic emission resulting from the model above.
   Note the large-scale velocity structure of the Milky Way is replicated.

It should be noted, though, that the procedure described thus-far in this 
section is for **one** model, but for scientific modelling it is likely useful 
to analyse the sythetic emission from a grid of models to constrain some 
parameters.
There is a convenient method to do this with the github repository.
From the root directory of the repo, we can run a grid of models using,

..  code:: bash

    $ mkdir ../kt3_models
    $ python terminal_scripts/run_model_grid.py -f ../kt3_models -m 0

By default, this will run a grid of three models of varying resolution 
(specifically :math:`400\,\mathrm{pc}`, :math:`400\,\mathrm{pc}`, and 
:math:`400\,\mathrm{pc}`), though you may notice that it takes a long time to 
finish.
For that reason, it might be better to set :code:`-m 8` for example to 
multiprocess the radiative transfer calculation.
At the moment, :meth:`calculateModel()` does not have the ability to utilise 
multiprocessing.
