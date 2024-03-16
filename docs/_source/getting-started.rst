***************
Getting Started
***************

The Single-Voxel model
======================

This is the base functionality of :code:`kosmatau3d`.
A *voxel*, or volumetric pixel, is one cell in a three-dimensional spatial
grid.
It contains one or more *ensembles* of KOSMA-:math:`\tau` clumps (see 
:doc:`Theory <theory>` for more information about these), which each contribute 
to the emissivity and absorption in the voxel.
It performs a probabilistic calculation of the emission by considering all 
possible combinations of clumps in the line-of-sight (see the schematic below).


.. figure:: _static/uniform_RT-small_comp.png
   :alt: voxel diagram
   :width: 500

   A diagram representing how the voxels work.
   The region in the red box is the ISM we want to model, with the clumps 
   shown as scatter points.
   The blue boxes are the voxels.
   In order for the probabilistic calculation to be accurate, the clumps are 
   randomly positioned in the voxel.

Three-dimensional PDR Models
============================

The most complex functionality of kosmatau3d, and the reason for its 
development since its conception 
(`Andree-Labsch et al. 2017 <https://ui.adsabs.harvard.edu/abs/2017A%26A...598A...2A/abstract>`_). 
This has so-far been used to model the Orion Bar and the Milky Way, though more 
models will soon be developed.
The figure below depicts the \[CII\] :math:`158\, \mu\mathrm{m}` integrated 
intensity in each voxel of one of the galactic models, which is then used to 
compute the synthetic observation.

.. figure:: _static/integrated_C+1.png
   :alt: model showing integrated C+

   One Galactic model, where the voxels are coloured according to the 
   \[CII\] :math:`158\, \mu\mathrm{m}` intensity integrated over the spectrum
   (:math:`-350` to :math:`350\, \mathrm{km\, s^{-1}}`).

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

.. figure:: _static/model_C+1.png
   :alt: model synthetic C+ 1

   The synthetic emission resulting from the model above.
   Note the large-scale velocity structure of the Milky Way is replicated.

