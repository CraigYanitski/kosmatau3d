***************
Getting Started
***************

The Single-Voxel model
======================

This is the base functionality of :code:`kosmatau3d`, and has been successfully 
used in PDR research.

KOSMA-:math:`\tau` Properties
=============================

It is imperative to understand the functionality and assumptions used in 
KOSMA-:math:`\tau` when using :code:`kosmatau3d` [#f1]_. 
KOSMA-:math:`\tau` is an isotropically-radiated, spherically symmetric PDR code.
It is maintained by Markus Röllig at the Universität zu Köln, and utilises 
chemical models from Amiel Sternberg at Tel Aviv University.

It assumes a piece-wise density profile approximating that of a Bonnor-Ebert 
sphere, that is, a gaseous sphere in pressure equilibrium with its environment.
For the total hydrogen surface density of the PDR :math:`n_s`, we define
the density profile as,

.. math::
   
   n(r) = 
   \begin{dcases}
   \left( \frac{r}{r_\mathrm{cl}} \right)^{-\gamma} & r > r_\mathrm{core} \\
   \left( \frac{r_\mathrm{core}}{r_\mathrm{cl}} \right)^{-\gamma} & r \leq r_\mathrm{core}
   \end{dcases},

for a given core radius :math:`r_\mathrm{core}` and powerlaw index 
:math:`\gamma`.
The for the KOSMA-:math:`\tau` grids in :code:`kosmatau3d`, we use the values
:math:`0.2 r_\mathrm{cl}` and :math:`1.5`, respectively.
With this definition, we can also write the total hydrogen number as,

.. math::
   N_H &= \int_0^{r_\mathrm{cl}} \mathrm{d}r 4 \pi r^2 n(r), \\
   &= \frac{4 \pi}{3} n_s r_\mathrm{cl}^3 \left( 1 - \frac{\gamma}{3} \left( 
   \frac{r_\mathrm{core}}{r_\mathrm{cl}} \right)^{3-\gamma} \right).



Three-dimensional PDR Models
============================

The most complex functionality of kosmatau3d, and the reason for its 
development since its conception (Andree-Labsch et al. 2017). 
This has so-far been used to model the Orion Bar and the Milky Way, though more 
models will soon be developed.


.. rubric:: Footnotes

.. [#f1]

   The spherical KOSMA-:math:`\tau` PDR models are referred to as clumps in the 
   context of :code:`kosmatau3d`.
   The distinction must be made that the `clumps` referenced in this 
   documentation are not by any means the astronomical clumps that will 
   eventually collapse into stellar clusters.
   :code:`kosmatau3d` is a `fractal` model to simulate the inhomogeneous 
   structure of the ISM using a multitude of smaller clumps (see eg. Stutzki 
   et al. 1998).
