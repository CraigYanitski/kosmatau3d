:html_theme.sidebar_primary.remove:

**************
Fractal Theory
**************

KOSMA-:math:`\tau` Properties
=============================

It is imperative to understand the functionality and assumptions used in 
KOSMA-:math:`\tau` when using :code:`kosmatau3d` [#f1]_. 
KOSMA-:math:`\tau` is an isotropically-radiated, spherically symmetric PDR code.
It is maintained by Markus Röllig at the Universität zu Köln, and utilises 
chemical models from Amiel Sternberg at Tel Aviv University.

Density
-------

It assumes a piece-wise density profile approximating that of a Bonnor-Ebert 
sphere, that is, a gaseous sphere in pressure equilibrium with its environment.
For the total hydrogen surface density of the PDR :math:`n_s`, we define
the density profile as,

.. math::
   
   n_\mathrm{H, cl}(r) = n_\mathrm{s}
   \left\{
      \begin{aligned}
         \left( \frac{r}{r_\mathrm{cl}} \right)^{-\gamma} 
         & \hspace{0.5cm} & r > r_\mathrm{core} \\
         \left( \frac{r_\mathrm{core}}{r_\mathrm{cl}} \right)^{-\gamma} 
         & \hspace{0.5cm} & r \leq r_\mathrm{core}
      \end{aligned}
   \right. \hspace{1cm} ,

for a given core radius :math:`r_\mathrm{core}` and powerlaw index 
:math:`\gamma`.
The for the KOSMA-:math:`\tau` grids in :code:`kosmatau3d`, we use the values
:math:`0.2 r_\mathrm{cl}` and :math:`1.5`, respectively.
With this definition, we can also write the total hydrogen number as,

.. math::
   N_\mathrm{H, cl} &= \int_0^{r_\mathrm{cl}} \mathrm{d}r\, 4 \pi\, r^2\, n(r), \\
   &= \frac{4 \pi}{3}\, n_\mathrm{s}\, r_\mathrm{cl}^3 \left( 1 - \frac{\gamma}{3} \left( 
   \frac{r_\mathrm{core}}{r_\mathrm{cl}} \right)^{3-\gamma} \right)
   \hspace{1cm} .

For :math:`\gamma\! >\! 0`, the number density increases towards the core, 
and is constant and maximal for :math:`r\! <\! r_\mathrm{core}` with 
corresponding density :math:`n_\mathrm{core}`.
It is also possible to derive the total hydrogen number as a function of 
density, but the piece-wise nature of our density profile means we will need
to separate our equation for :math:`n\! <\! n_\mathrm{core}` and 
:math:`n\! =\! n_\mathrm{core}`.
Momentarily neglecting the core, so 
:math:`n\! \in\! \left[ n_\mathrm{s}, n_\mathrm{core} \right)`, we can write 
the dependence of the total number of hydrogen atoms as a function of density:

.. math::
   \frac{\mathrm{d}N_\mathrm{H, cl}}{\mathrm{d}n} &= 
   \frac{\mathrm{d}N_\mathrm{H, cl}}{\mathrm{d}r}
   \left( \frac{\mathrm{d}n_\mathrm{H, cl}}{\mathrm{d}r} \right)^{-1} 
   \hspace{1cm} , \\
   &= \left( -\gamma\, r_\mathrm{cl}\, n_\mathrm{s} 
   \left( \frac{n}{n_\mathrm{s}} \right)^{\frac{\gamma + 1}{\gamma}} \right)^{-1} 
   4\pi\, r_\mathrm{s}^2 \left( \frac{n}{n_\mathrm{s}} \right)^{- \frac{2}{\gamma}} 
   n_\mathrm{s} \frac{n}{n_\mathrm{s}} \hspace{1cm} , \\
   &= - \frac{4\pi\, r_\mathrm{cl}}{\gamma} 
   \left( \frac{n}{n_\mathrm{s}} \right)^{-\frac{3}{\gamma}} \hspace{1cm} ,

where we have expressed the radius as a function of density. 
For :math:`n_\mathrm{H, cl}\! =\! n_\mathrm{core}`, we can simply perform a 
spherical integration with a constant density to derive the total number of 
hydrogen atoms. 
Since the core has constant density, we can write the final form of the density
dependence of the total number of hydrogen atoms:

.. math::
   \frac{\mathrm{d}N_\mathrm{H, cl} (n)}{\mathrm{d}n} = 
   \left\{
      \begin{aligned}
         - \frac{4\pi\, r_\mathrm{cl}}{\gamma} 
         \left( \frac{n_\mathrm{cl}}{n_\mathrm{s}} \right)^{-\frac{3}{\gamma}} 
         & \hspace{0.5cm} & n_\mathrm{s} < n < n_\mathrm{core} \\
         0 & \hspace{0.5cm} & n = n_\mathrm{core}
      \end{aligned}
   \right. \hspace{1cm} ,

The density probability distribution function (PDF) for the spherical clump 
can be defined as,

.. math::
   \mathcal{P}_\mathrm{cl}(n) \equiv N_\mathrm{H, cl}^{-1} 
   \frac{\mathrm{d}N_\mathrm{H, cl} (n)}{\mathrm{d}n} 
   \hspace{1cm} .

Using the density profile of the KOSMA-:math:`\tau` clumps, as well as manually 
integrating the core to derive its probability, we obtain,

.. math::
   \mathcal{P}_\mathrm{cl}(n) = N_\mathrm{H, cl}^{-1} 4\pi r_\mathrm{cl}^3
   \left\{
      \begin{aligned}
         - \frac{1}{\gamma} \left( \frac{n}{n_\mathrm{s}} \right)^{-\frac{3}{\gamma}} 
         & \hspace{0.5cm} & n_\mathrm{s} < n < n_\mathrm{core} \\
         \frac{1}{3} \left( \frac{n}{n_\mathrm{s}} \right)^{\frac{\gamma - 3}{\gamma}} 
         & \hspace{0.5cm} & n = n_\mathrm{core}
      \end{aligned}
   \right. \hspace{1cm} .

How this is utilised in :code:`kosmatau3d` for the fractal ISM will soon be 
explained in the Ensembles_ section.

Far-UV radiation
----------------

Currently [#f2]_ KOSMA-:math:`\tau` uses a modified Draine spectrum for the spectral 
energy distribution (SED) incident on the clump surface, which can be scaled 
by a given factor.
We therefore denote the far-UV radiation as :math:`\chi` in units of 
:math:`\chi_\mathrm{D}`.

More information about the various far-UV spectra will be provided here upon 
the initial publication for :code:`kosmatau3d`.

Ensembles
=========

Here I will soon explain the mathematics and statistics associated with 
groups of clumps in an ensemble.
Most of this information will be placed here after the publication of Yanitski 
et al. (2024).

The basis of the fractal approximation of the ISM is the clump mass 
distribution:

.. math::
   \frac{\mathrm{d} N_\mathrm{cl}}{\mathrm{d} m_\mathrm{cl}} 
   \propto m_\mathrm{cl}^{-\alpha} \hspace{1cm} ,

where :math:`N_\mathrm{cl}` is the number of clumps of a given mass 
:math:`m_\mathrm{cl}` and :math:`\alpha` is the index of the clump mass 
distribution.
This distribution is bourne out of turbulence, and may vary between star 
forming regions.
The value we used in the Milky Way is :math:`1.84`, which was determined from 
the Polaris Flare 
(`Heithausen et al. 1998 <https://ui.adsabs.harvard.edu/abs/1998A%26A...331L..65H/abstract>`_).

We combine the clump mass distribution with a mass-size relation,

.. math::
   m_\mathrm{cl} \propto r_\mathrm{cl}^\varpi \hspace{1cm} ,

where :math:`r_\mathrm{cl}` is the radius of the clump and :math:`\varpi` is 
mass-size index (determined to be :math:`2.31` in 
`Heithausen et al. 1998 <https://ui.adsabs.harvard.edu/abs/1998A%26A...331L..65H/abstract>`_), 
in order to approximately derive the clump surface density distribution.

.. todo::
   This derivation is being completed at the moment.
   Please be patient :-)

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

.. [#f2]

   At some point in the future this will be extended to utilise a user-defined 
   SED, but currently there is nobody developing this.
   It is particularily important in order to use KOSMA-:math:`\tau` to model 
   the X-ray dominated regions (XDRs) around active galactic nuclei (AGNs).