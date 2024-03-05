***********
Development
***********

As seen on the home page, this code is under active development.
There are many planned extensions/improvements, so we offer this page to 
describe them.

Developers
==========

Active contributers
-------------------

.. Here the columns are for the developer and their affiliation.
   We will later include information about what they develop.

=======================================================   ===================
=======================================================   ===================
`Craig N. Yanitski <https://github.com/CraigYanitski>`_   Universität zu Köln
`Aditi Sinha <https://github.com/aditi0009>`_             Universität zu Köln
=======================================================   ===================

Previous contributers
---------------------

.. Here the columns are for the previous contributer and the version(s) of the
   code they worked on

===================   =====================
===================   =====================
Silke Andree-Labsch   (depreciated version)
Christoph Bruckmann   (depreciated version)
===================   =====================

Planned development
===================

Below we list the developments planned for upcoming versions of 
:code:`kosmatau3d`.

+---------+-------------------------------------------------------------------+
| version | changes                                                           |
+=========+===================================================================+
| v1.1.0  | - restructuring of code                                           |
|         |   - :code:`masspoints` -> :code:`clumps`                          |
+---------+-------------------------------------------------------------------+
| test    | ...                                                               |
| multi-  |                                                                   |
| lines   |                                                                   |
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

.. toctree::
   :maxdepth: 2
   :caption: Content
