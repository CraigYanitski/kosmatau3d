.. kosmatau3d documentation master file, created by
   sphinx-quickstart on Fri Mar 11 09:04:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

****************************************
The documentation of a clumpy PDR model.
****************************************

:code:`kosmatau3d` is a clumpy photo-dissociation region (PDR) model built on 
the KOSMA-:math:`\tau` PDR model from the Universität zu Köln and Tel Aviv 
University.
It's functionality and use in modelling the Galactic cooling lines in the 
Milky Way are explained in 
`Yanitski (2023)`_.
It will soon be able to download its required data from the 
`InterStellar Medium DataBase (ISMDB)`_ hosted at the Observatoire de 
Paris-Meudon when their API is completed.

.. _ISMDB: https://ismdb.obspm.fr/
.. _InterStellar Medium DataBase (ISMDB): ISMDB_

.. caution::
   This documentation and the corresponding code are under active development.
   Please send all questions and issues to yanitski@ph1.uni-koeln.de.

.. rubric:: References

Yanitski 2023, :emphasis:`The Milky Way with kosmatau3d: Modelling the Galactic 
cooling lines using clumpy PDRs`, PhD thesis, (Universität zu Köln)

.. _yanitski2023url: https://kups.ub.uni-koeln.de/71850/
.. _Yanitski 2023: yanitski2023url_
.. _Yanitski (2023): yanitski2023url_

.. toctree::
   :maxdepth: 2
   :caption: Sections:

   Installation <installation>
   Getting Started <getting-started>
   API <api>
   Development <development>
   Acknowledgements <acknowledgements>

.. Clumps <properties>
   Theory <theory>
   Voxels <voxel>
   3D Models <model>

******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
