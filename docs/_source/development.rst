:html_theme.sidebar_primary.remove:

***********
Development
***********

As seen on the home page, this code is under active development.
There are many planned extensions/improvements, so we offer this page to 
describe them.

Developers
==========

Active contributors
-------------------

.. Here the columns are for the developer and their affiliation.
   We will later include information about what they develop.

=======================================================   ===================
name                                                      affiliation
=======================================================   ===================
`Craig N. Yanitski <https://github.com/CraigYanitski>`_   Universität zu Köln
`Aditi Sinha <https://github.com/aditi0009>`_             Universität zu Köln
=======================================================   ===================

Previous contributors
---------------------

.. Here the columns are for the previous contributor and the version(s) of the
   code they worked on.

===================   ===========================
name                  version(s) [#f1]_
===================   ===========================
Silke Andree-Labsch   (deprecated version)
Christoph Bruckmann   (deprecated version)
===================   ===========================

Planned development
===================

Below we list the developments planned for upcoming versions of 
:code:`kosmatau3d`.

+---------+-------------------------------------------------------------------+
| version | changes                                                           |
+=========+===================================================================+
| v1.0.8/ | - new method in :class:`models.Voxel()` to return species column  |
|         |   density                                                         |
| v1.0.10 |                                                                   |
|         |   - should have similar inputs to                                 |
|         |     :meth:`models.Voxel.get_abundances()` and                     |
|         |     :attr:`models.constants.voxel_size` for calculation.          |
+---------+-------------------------------------------------------------------+
| v1.0.11 | - column density method new features                              |
|         |                                                                   |
|         |   - will include a voxel size-invariant option that adds the      |
|         |     column density of each clump in the ensemble                  |
|         |                                                                   |
|         |   - will have beam size and possible filling factor as optional   |
|         |     arguments (to mimic the work in §7 of Röllig &                |
|         |     Ossenkopf-Okada 2022)                                         |
+---------+-------------------------------------------------------------------+
| ...     | ...                                                               |
+---------+-------------------------------------------------------------------+
| v1.1.0  | - restructuring of code                                           |
|         |                                                                   |
|         |   - :mod:`models.masspoints` -> :mod:`models.clumps`              |
|         |                                                                   |
|         |   - :mod:`models.observations` -> :mod:`models.model_data`        |
|         |                                                                   |
|         | - removal of circular imports                                     |
+---------+-------------------------------------------------------------------+
| v1.1.x  | - New features                                                    |
|         |                                                                   |
|         |   - python implementation of the Mathematica routines from        |
|         |     Markus Röllig                                                 |
|         |                                                                   |
|         |     - will be placed in :mod:`kosmatau3d.kosmatau`                |
|         |                                                                   |
|         |     - should help streamline the processing of the                |
|         |       KOSMA-:math:`\tau` output and compare to :code:`kosmatau3d` |
|         |       results                                                     |
|         |                                                                   |
|         |   - unify loading of data files to use :code:`pandas`             |
|         |                                                                   |
|         |   - parsing of KOSMA-:math:`\tau` grid parameters directly from   |
|         |     header                                                        |
|         |                                                                   |
|         |     - should be simple to implement using :code:`pandas` to open  |
|         |       the file                                                    |
|         |                                                                   |
|         |     - one should also implement a failsafe to ignore parameters   |
|         |       that do not change, since this will cause an error          |
+---------+-------------------------------------------------------------------+
| ...     | ...                                                               |
+---------+-------------------------------------------------------------------+


Potential development
=====================

The following features have been identified as useful, but currently nobody
is implementing them.
The maintainers of :code:`kosmatau3d` are happy to work with anybody who wishes
to develop these features.

+-------------+---------------------------------------------------------------+
| feature     | explanation                                                   |
+=============+===============================================================+
| cython      | From the beginning it was known that compiling the code in    |
|             | cython would dramatically improve its efficiency.             |
|             | The issue with this is that the code needed to be developed   |
|             | before we would understand the most-efficient implementation  |
|             | (you can blame Craig Yanitski for that).                      |
|             | Now that the code is mostly working, it should not be too     |
|             | much work to rewrite in cython.                               |
+-------------+---------------------------------------------------------------+
| GUI         | This existed in the first iteration of KOSMA-:math:`\tau` 3D  |
|             | developed by Silke Andree-Labsch (see Andree-Labsch et al.    |
|             | 2017).                                                        |
|             | While it existed to give some order to the series of scripts  |
|             | that existed before implementing an object-orientated         |
|             | approach, it is still a good idea to implement this and give  |
|             | some order to the slew of properties offered by the model.    |
|             | Likewise it should be possible to view an interactive         |
|             | representation of the model as well as compute the synthetic  |
|             | observation.                                                  |
+-------------+---------------------------------------------------------------+
| recursive   | :mod:`radiativeTransfer` is currently functioning by          |
| RT          | looping through the voxels in a sightline to integrate the    |
|             | radiative transfer equation.                                  |
|             | This can and should be cleaned-up by rewriting this as a      |
|             | recursive function that can be called for one sightline.      |
|             | It might involve improving the current linear approximation   |
|             | used in the computation (see                                  |
|             | `Yanitski 2023 <https://kups.ub.uni-koeln.de/71850/>`_).      |
+-------------+---------------------------------------------------------------+

.. rubric:: Footnotes

.. [#f1]

   The *deprecated* version of :code:`kosmatau3d` is the version written by 
   Silke Andree-Labsch c. 2015 (KOSMA-:math:`\tau` 3D) that predates the 
   publicly-available code
