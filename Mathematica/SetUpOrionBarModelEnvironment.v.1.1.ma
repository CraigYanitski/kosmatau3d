#!/usr/local/bin/MathKernel -script

SetDirectory["./Mathematica/"]
modelname = "ORIONBARgrid";
Get[modelname <> ".modelDefinitions", "`*"];

Print["Reading Ensemble Parameters into Mathematica notebook..."]
SetDirectory["../data/"]
par = Import["parameters.dat", "Table"];
l = Length[par]-3;

emissivity = Table[0,{i,1,l}];
For[i = 4, i <= l + 3, i++, {emissivity[[i-3]] = clumpyFluxErg[{par[[2]][[2]],par[[3]][[2]]},\
{par[[i]][[2]],par[[i]][[3]],par[[i]][[4]],par[[i]][[5]],1(*1pc*),par[[i]][[6]],par[[i]][[7]]}]\
*(1^2/par[[1]][[2]]^2)/(30.656776*10^(17)*par[[1]][[2]])}]

(*volume emissivity in erg s^-1 cm^-3 sr^-1*)
(*solid angle of one pixel is pixel_edge_length[pc]^2/1pc^2*)
(*Intensity divided by pixel length in cm gives volume emissivity in erg s^-1 cm^-3 sr^-1*)

Print["Exporting volume emissivities to ./data/emissivity.dat..."]
Export["emissivity.dat", emissivity]

