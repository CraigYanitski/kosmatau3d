This is the development branch for KOSMA-tau^3. My intention here is to increase encapsulation
in an effort to improve the code's efficiency and make it easier to develop overall. A side-
effect of this is that I will be intimately familiar with all of the nuances of setting up a
PDR model.

This development branch is very much a work-in-progress and will not yet compile. I am slowly
linking all of these classes to reproduce the galaxy model by Christof Bruckmann.

Among other things, I believe I have successfully taken out the incredibly time-consuming use
of the scipy function griddata() for each species for each masspoint in each voxel. I have replaced
it with an initial call to either the scipy function LinearNDInterpolation() or Rbf() to create
a list of functions that will interpolate the grid of data. This should save on time enormously
for large, high-resolution models.
