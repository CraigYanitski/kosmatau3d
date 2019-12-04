This is the development branch for KOSMA-tau^3. My intention here is to increase encapsulation
in an effort to improve the code's efficiency and make it easier to develop overall. A side-
effect of this is that I will be intimately familiar with all of the nuances of setting up a
PDR model.

This development branch is finally functional. It has not been converted to the master branch
since I am still debugging and ensuring it can reproduce the same result as the old code. The
most recent version of KOSMA-tau 3D is the galaxy model developed by Christof Bruckmann. When this
branch can create the same model, I will convert it to the master branch. Comparison of the two
branches is rather difficult since I am also debugging the old code. The major changes are
described in the document treatise.pdf.

Among other things, I believe I have successfully taken out the incredibly time-consuming use
of the scipy function griddata() for each species for each masspoint in each voxel. I have replaced
it with an initial call to either the scipy function LinearNDInterpolation() or Rbf() to create
a list of functions that will interpolate the grid of data. This has greatly reduced the execution
time of large, high-resolution models.
