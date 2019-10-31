This is a stable(ish) version of the kosma-tau_3D code I slightly modified to
run in python3.7.

Basically all I have done is make a global dictionary truely global,
split the code into blocks, and created a nice terminal prompt to run all of this.
I also added parallelisation, but this is very much a work-in-progress.
It should complete the necessary function in parallel, but then it starts another
loop afterwards to save these values to the global dictionary. This is
one of the many reasons I loathe global variables.

As for the way the code is segmented, it begins by determining the
number and position of voxels in the model in setupModel(),
then it has two short segments to determine which molecules/elements it
will consider (initVelocities()) and to calculate the FUV background of the
voxel grid (calculateFUV(); the FUV is taken from a model in Wolfire et al. 2002).
The parallelised function occurs in the lineEmissionAbsorption() segment, where
it needs to loop over all of the voxels, their ensembles,  combinations,
and all of the KOSMA-tau clumps required for the superposition of clumpy
PDRs. Finally, the radiative transfer along the line-of-sight is calculated
in writeFiles(), in addition to the creation of the necessary .fits and
.dat files needed for logging.

Now that this is on a git repository, I can properly develop it and correct the
aforementioned issues. . .