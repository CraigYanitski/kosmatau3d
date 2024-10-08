**************************
Comparison to observations
**************************

Of course every astrophysicist has their own method of comparing observed and 
synthetic observations, but a few tools are included in :code:`kosmatau3d`.
What is crutial to remember is what is being calculated and how is it stored.
:code:`kosmatau3d` does not offer the standard datacubes for a particular 
wavelength/transition, but rather it saves a sythetic map for each transition 
specified in the model (a four-dimensional dataset).
This applies for both the species transition and the dust continuum emissions, 
each of which are stored in an :code:`astropy` HDU.
There is also an HDU containing the sightline positions, which is a remanent 
from a previous iteration of the code.
Thus the file :code:`channel_intensity.fits` contains three HDUs: one for 
position, one for species transition, and one for dust emission.

The simplest way to view the data is by using the built-in 
:code:`kosmatau3d.models.plotting.view_map()` method.
This is built to load the synthetic observations and enable one view all 
channel maps, integrated intensity maps, and position-velocity diagrams.

.. toctree::
   :maxdepth: 2
   :caption: Contents
