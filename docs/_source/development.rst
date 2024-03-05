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

===================   ============================
name                  version(s)
===================   ============================
Silke Andree-Labsch   (depreciated version) [#f1]_
Christoph Bruckmann   (depreciated version) [#f1]_
===================   ============================

Planned development
===================

Below we list the developments planned for upcoming versions of 
:code:`kosmatau3d`.

+---------+-------------------------------------------------------------------+
| version | changes                                                           |
+=========+===================================================================+
| v1.1.0  | - restructuring of code                                           |
|         |                                                                   |
|         |   - :mod:`masspoints` -> :mod:`clumps`                      |
|         |                                                                   |
|         |   - :mod:`observations` -> :mod:`model_data`                |
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
| recursive   | :mod:`radiativeTransfer` is currently functioning by       |
| radiative   | looping through the voxels in a sightline to integrate the    |
| transfer    | radiative transfer equation.                                  |
|             | This can and should be cleaned-up by rewriting this as a      |
|             | recursive function that can be called for one sightline.      |
|             | It might involve improving the current linear approximation   |
|             | used in the computation (see                                  |
|             | `Yanitski 2023 <https://kups.ub.uni-koeln.de/71850/>`_).      |
+-------------+---------------------------------------------------------------+

.. rubric:: Footnotes

.. [#f1]

   The *depreciated* version of :code:`kosmatau3d` is the version written by 
   Silke Andree-Labsch c. 2015 (KOSMA-:math:`\tau` 3D) that predates the 
   publically-available code

.. toctree::
   :maxdepth: 2
   :caption: Content
