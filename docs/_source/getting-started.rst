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
It performs a probabilistic calculation of the the emission by considering all 
possible combinations of clumps in the line-of-sight.

.. figure:: _static/uniform_RT-small_cropped.png
   :alt: voxel diagram
   :width: 500

   A diagram representing how the voxels work.
   The region in the red box is the ISM we want to model, with the clumps 
   shown as scatter points.
   The blue boxes are the voxels.

Three-dimensional PDR Models
============================

The most complex functionality of kosmatau3d, and the reason for its 
development since its conception (Andree-Labsch et al. 2017). 
This has so-far been used to model the Orion Bar and the Milky Way, though more 
models will soon be developed.

.. figure:: _static/integrated_C+1.png
   :alt: model showing integrated C+

   One Galactic model, where the voxels are coloured according to the 
   \[CII\] :math:`158\, \mu\mathrm{m}` intensity integrated over the spectrum
   (:math:`-350` to :math:`350\, \mathrm{km\, s^{-1}}`).