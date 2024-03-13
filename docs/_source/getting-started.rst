***************
Getting Started
***************

The Single-Voxel model
======================

This is the base functionality of :code:`kosmatau3d`.
A *voxel*, or volumetric pixel, is one cell in a three-dimensional spatial
grid.
It contains one or more *ensembles* of KOSMA-:math:`\tau` clumps (see 
`Theory <theory>` for more information about these), which each contribute to 
the emissivity and absorption in the voxel.
It performs a probabilistic calculation of the the emission by considering all 
possible combinations of clumps in the line-of-sight.

.. image:: _static/uniform_RT-small.pdf

Three-dimensional PDR Models
============================

The most complex functionality of kosmatau3d, and the reason for its 
development since its conception (Andree-Labsch et al. 2017). 
This has so-far been used to model the Orion Bar and the Milky Way, though more 
models will soon be developed.

.. image:: _static/integrated_C+1.png

.. rubric:: Footnotes

.. [#f1]

   The spherical KOSMA-:math:`\tau` PDR models are referred to as clumps in the 
   context of :code:`kosmatau3d`.
   The distinction must be made that while the *clumps* referenced in this 
   documentation are the astronomical clumps that will eventually collapse into
   stars, they are approximately in hydrostatic equilibrium.
   Thus we are able to compute an instantaneous synthetic  observation.
   :code:`kosmatau3d` is a *fractal* model to simulate the inhomogeneous 
   structure of the ISM using a multitude of these smaller clumps (see eg. 
   Stutzki et al. 1998).
