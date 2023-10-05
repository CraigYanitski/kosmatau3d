KOSMA-:math:`\tau` Properties
=============================

It is imperative to understand the functionality and assumptions used in KOSMA-:math:`\tau` when using :code:`kosmatau3d` [#f1]_. KOSMA-:math:`\tau` is an isotropically-radiated, spherically symmetric PDR code.
It is maintained by Markus Röllig at the Universität zu Köln, and utilises chemical models from Amiel Sternberg at Tel Aviv University.
It assumes the density profile of a Bonnor-Ebert sphere, that is, a gaseous sphere in pressure equilibrium with its environment.

.. rubric:: Footnotes

.. [#f1]

   The spherical KOSMA-:math:`\tau` PDR models are referred to as clumps in the context of :code:`kosmatau3d`.
   The distinction must be made that the `clumps` referenced in this documentation are not by any means the astronomical clumps that will eventually collapse into stellar clusters.
   :code:`kosmatau3d` is a `fractal` model to simulate the inhomogeneous structure of the ISM using a multitude of smaller clumps (see eg. Stutzki et al. 1998).
