Tb_linecentre.dat contains the line centre intensities [K] for single
clumps. Calulated assuming a linewidth of FWHM = 1.67 km/s.

surface density, mass, and uv field (10Log(..) scale) 
are given in the first three columns. 

species afterwarts in the order

In the input files the species need to be given in the following order 
(Careful, new onion output may be different)!!:
species:                 {C+, C, O, CO, 13CO, 13C+, 13C, HCO+, H13CO+, H3O+}
transitions per species: {1,  3, 3, 49, 49,   1,    3,   15,   30,     17}

##########################################################################
tau_linecentre.dat clump averaged optical depth. 
Otherwise it's the same as Tb_linecentre.dat.

clump averaged tau's calculated via
tauAv = -Log[2/(Rcl^2)*NIntegrate[Exp[-tau[r]]*r, {r, 0, Rcl}]];

##########################################################################
RhoMassAFUV.dat
contains table with clump surface densities (10Log...scale), 
masses (10Log...scale), and clump averaged FUV extinction 
A_FUV (not tau, but A). 

#########################################################################





with open('tau.dat') as ta:
    for line in ta:
        tautemp = line.split()[114:115][0]
        tautemp = float(tautemp)
        tauAv.append(tautemp)

rhoMassUV_tau=[]
with open('tau.dat') as ta:
    for line in ta:
        rhoMassUV_tau.append(line.split()[0:3])
