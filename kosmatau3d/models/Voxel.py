import importlib as il
import numpy as np
from scipy.interpolate import interp1d
from copy import copy

from . import constants
from . import interpolations
from . import observations
from . import species

from . import ensemble
from . import combinations
from . import masspoints

class Voxel(object):
  '''
  This is a class to handle each voxel in KOSMA-tau^3. It contains ensembles
  of spherical PDR simulations (to mimic the fractal structure of PDRs), a FUV
  field element, and element absorption, separated for dust and molecules. It
  should have one clump ensemble and one interclump ensemble. This is to account
  for the diffuse surroundings of a clump.
  '''
  # PRIVATE
  def __init__(self, index=0, debugging=False):
    
    self.__index = index         #index of voxel in VoxelGrid, sort of like its ID in the overall model
    
    self.__debugging = debugging
    
    self.__mass = 0
    self.__clumpMass = 0
    self.__velocity = 0     #velocity of mass at voxel point
    self.__ensembleDispersion = 0   #dispersion of velocity at voxel point
    
    self.__clumpVelocityIndeces = []
    self.__modelMass = []

    self.__intensity = 0      #intensity of emissions at voxel point
    self.__opticalDepth = 0   #optical depth at voxel point
    
    self.__FUV = 0.
    self.__Afuv = 0.
    
    self.__x = 0
    self.__y = 0
    self.__z = 0
    
    self.__r = 0
    self.__phi = 0
    
    return

  def __setMass(self):
    self.__mass = self.__clumpMass.sum()

  def __setClumpMass(self, r):
    mass = [interpolations.interpolateClumpMass(r), interpolations.interpolateInterclumpMass(r)]
    self.__clumpMass = constants.clumpMassFactor*np.asarray(mass).mean(1)
    return

  def __setVelocity(self, r):
    velocity = interpolations.interpolateRotationalVelocity(r)
    
    if constants.fromEarth:

      # Calculate the correction to the voxel velocity vectors
      relativeRpol = np.sqrt((self.__x-constants.rGalEarth)**2+self.__y**2)
      relativePhi = np.arctan2(self.__y, self.__x-constants.rGalEarth)
      relativeSigma = np.arccos((self.__r**2+relativeRpol**2-constants.rGalEarth**2)/(2*self.__r*relativeRpol))
      sigma = np.arctan2(self.__z, abs(self.__x-constants.rGalEarth))

      # Correct the relative velocity of the voxel
      velocityEarth = interpolations.interpolateRotationalVelocity(constants.rGalEarth)
      velocityCirc = velocity.mean() - velocityEarth*self.__r/constants.rGalEarth

      self.__velocity = np.sign(relativePhi) * velocityCirc * np.sin(relativeSigma) * np.cos(sigma)

      if self.__r==0: self.__velocity = 0
      #self.__velocity = (velocity.mean()) * np.sin(self.__phi)
    
    else:
      self.__velocity = np.array(velocity)

    # Use this to check the evaluation of the velocity field. It is still not working correctly...
    #print(self.__velocity)

    ensembleDispersion = interpolations.interpolateVelocityDispersion(r)
    
    self.__ensembleDispersion = ensembleDispersion.mean()
    
    return

  def __setEnsembleDensity(self, r):
    density = interpolations.interpolateDensity(r)
    self.__ensembleDensity = constants.densityFactor*density.mean()
    return

  def __setExtinction(self):
    self.__UVextinction = interpolations.interpolateFUVextinction(self.__ensembleDensity, self.__clumpMass+self.__interclumpMass)
    return

  def __setFUV(self, r, z):
    # This is in units of the Draine field
    fuv = interpolations.interpolateFUVfield(r, z)/constants.normUV*constants.globalUV
    self.__FUV = np.clip(fuv, 1, None)
    #self.__FUV = FUVfield(fuv)
    return

  def __str__(self):
    return 'Voxel {}'.format(self.__index)

  # PUBLIC

  # def reloadModules(self):
  #   il.reload(Ensemble)
  #   il.reload(FUVfield)
  #   self.__clump.reloadModules()
  #   self.__interclump.reloadModules()
  #   return

  def setIndex(self, index):
    self.__index = index
    return

  def getIndex(self):
    return self.__index

  def setPosition(self, x, y, z, r, phi):
    self.__x = x
    self.__y = y
    self.__z = z
    self.__r = r
    self.__phi = phi
    return

  def getFUV(self):
    return self.__FUV

  def setProperties(self, voxel_size=1, molecules='all', dust='PAH', alpha=1.84, gamma=2.31, \
                    clumpMassNumber=[3,1], clumpMassRange=[[0,2],[-2]], \
                    clumpNmax=[1, 100], velocityRange=[-10,10], velocityNumber=51, \
                    clumpMass=100, velocity=0., \
                    ensembleDispersion=1, volumeFactor=None, ensembleDensity=[15000, 1911], FUV=[20000, 1], fromFile=False):
    '''
      This method calculates the radii assuming an origin of (0,0). It then averages
     over a subgrid of 3x3. It might be improved later by having functionality to
     change the subgrid dimensions.

     There are a few kwargs that can be used to initialise a Voxel instance:

     MODEL PARAMETERS

     alpha: The power law distribution required for initial mass function dN(M)/dM = M^-alpha.
            The default is 1.84 as determined in Heithausen et al. (1998).
     gamma: The power law mass-size relation satisfying M = C R^gamma. The default is 2.31 as
            determined in Heithausen et al. (1998).
     velocityRange: The range (list) of the observed radial velocities. By default it
                    is [-10, 10].
     velocityNumber: The number of observed velocities in the specified range (including
                     endpoints). The default is 51.
     clumpMassNumber: The number of clump masses for each set. This should be a list of
                      integers with the same length as clumpMassRange. The default is
                      [3, 1].
     clumpMassRange: The range (list) of clump masses for each set. This is a list of
                     ranges with the same length as clumpMassNumber. The ranges can have
                     length 1, in which case it only evaluates that clump mass. The
                     masses are in units of dex(Msol). By default it is [[0, 2], [-2]].
     clumpDensity: This changes the default clump density in the constants module. This
                   overrides the voxel density kwarg, and is useful for creating interclumps
                   medium that is not scaled according to the clumps. the default is
                   [None, 1911], so the first clump set takes the voxel density. This is in
                   units of cm^-3.
     clumpFUV: This changes the default clump FUV in the constants module. This overrides
               the voxel FUV kwarg. It is a list of length of the number of clump sets.
               Use a value of None to use the voxel FUV field. The default is [None, 1], in
               units of the Draine field.
     clumpNmax: The maximum number of the largest clump in the ensemble. This helps to limit
                the calculation of a large combination of clumps by rescaling the voxel (thus
                eliminating any change to the brightness temperature). Larger combinations are used
                if the rescaled voxel is smaller than the largest clump. The default is [1, 100]. 
     voxel_size: The voxel size in parsecs. The default is 1 pc.
     molecules: The molecules included in the model. This is a list of strings, where each
                string has the element name (as in the grid) followed by the transition
                number (Which is taken from KOSMA-tau). It is set to the first transition
                of C+, C, and CO by default.
     dust: This is a string to limit how much of the dust continuum is included. The
           available continuum ranges from 1 nm to 3.1 mm. This can be limited to the range
           of molecular lines by setting dust='molecular' (22 dust lines), or to the range of
           PAH dust features by setting dust='PAH' (201 dust lines). All dust lines (333)
           will be calculated if dust=''

     VOXEL PROPERTIES

     clumpMass: This can be either a float or list of floats, depending on how many clump
                sets you have. The default is 100 Msol.
     mass: This is mostly useless. It is used for analysing a grid of Voxel instances. The
           defalt is 100 Msol.
     velocity: The radial velocity of the observed Voxel instance. This is used with the
               model velocities to determine the number of clumps seen at a given velocity.
               The default is 0 km/s.
     ensembleDispersion: The velocity dispersion of the ensemble. This is used with the
                         model velocities to determine the number of clumps at a given
                         velocity. The default is 1 km/s.
     ensembleDensity: The observed hydrogen density in the voxel. This is different from and overridden by
              any default density specified in the constants module. the default is
              15000 cm^-3.
     FUV: The FUV field of the voxel. All clumps are isotropically radiated with this
          field. It is different from and overridden by any default FUV field specified
          in the constants module. The default is 20000 Draine units.

     GENERAL FLAGS

     fromFile: This flag can be set to retrieve the voxel properties from a file. It is
               False by default.
     verbose: This is mainly used for debugging purposes. It prints various statements
              as the code is evaluated

    '''
    #print('Voxel instance initialised')

    x,y = np.meshgrid(np.linspace(self.__x-.5*constants.voxel_size, self.__x+.5*constants.voxel_size,3), \
                      np.linspace(self.__y-.5*constants.voxel_size, self.__y+.5*constants.voxel_size,3))
    r = np.array([x.flatten(), y.flatten()]).T
    r = np.linalg.norm(r, axis=1)

    if fromFile:

      self.__setClumpMass(r)
      self.__setMass()
      self.__setVelocity(r)
      self.__setEnsembleDensity(r)
      #self.__setExtinction()
      self.__setFUV(self.__r, self.__z)

      # This can convert the velocity from Cartesian to radial. It is assumed a scalar value is a radial velocity.
      if isinstance(self.__velocity, float):
        velocity = self.__velocity
      else:
        velocity = np.linalg.norm(self.__velocity)

    else:

      constants.voxel_size = voxel_size
      constants.changeMassFunctionParameters(alpha=alpha, gamma=gamma)
      constants.changeVelocityRange(velocityRange)
      constants.changeVelocityNumber(velocityNumber)
      constants.addClumps(massRange=clumpMassRange, num=clumpMassNumber, density=ensembleDensity, fuv=FUV, Nmax=clumpNmax, reset=True)
      constants.changeDustWavelengths(dust)

      observations.methods.initialise()
      species.addMolecules(molecules)
      interpolations.initialise()

      masspoints.reinitialise()
      combinations.reinitialise()
      ensemble.reinitialise()
      
      self.__velocity = velocity

      if isinstance(clumpMass, list) or isinstance(clumpMass, np.ndarray):
        self.__clumpMass = clumpMass
      else:
        self.__clumpMass = [clumpMass] * len(constants.clumpMassNumber)
      
      if isinstance(ensembleDispersion, list) or isinstance(ensembleDispersion, np.ndarray):
        self.__ensembleDispersion = ensembleDispersion
      else:
        self.__ensembleDispersion = [ensembleDispersion] * len(constants.clumpMassNumber)

      if volumeFactor:
        ensembleDensity = clumpMass*constants.massSolar/constants.massH/volumeFactor/constants.voxel_size**3/constants.pc**3/100**3
      if isinstance(ensembleDensity, list) or isinstance(ensembleDensity, np.ndarray):
        self.__ensembleDensity = ensembleDensity
      else:
        self.__ensembleDensity = [ensembleDensity] * len(constants.clumpMassNumber)

      if isinstance(FUV, list) or isinstance(FUV, np.ndarray):
        self.__FUV = FUV
      else:
        self.__FUV = [FUV] * len(clumpMassNumber)

      velocity = self.__velocity

    # This will allow the code to reuse the standard clump density constants for voxel sets (ie. do not alter the model constants)
    # density = copy(constants.clumpDensity)
    # fuv = copy(constants.clumpFUV)
    # for i in range(len(density)):
      # if density[i]=='auto': density[i] = self.__ensembleDensity
      # if fuv[i]==None: fuv[i] = self.__FUV

    masspoints.setMasspointData(density=self.__ensembleDensity, FUV=self.__FUV)
    ensemble.initialise(velocity=velocity, ensembleDispersion=self.__ensembleDispersion, clumpMass=self.__clumpMass)
    combinations.initialise(clumpCombination=[ensemble.clumpCombinations[ens][ensemble.clumpLargestIndex[ens]] for ens in range(len(constants.clumpMassNumber))])

    for ens in range(len(constants.clumpMassNumber)):
      self.__modelMass.append((ensemble.clumpDeltaNji[ens].sum(1)*10**constants.clumpLogMass[ens]).sum())
      if abs(self.__modelMass[ens]-self.__clumpMass[ens])>0.1*self.__clumpMass[ens]:
        print('ERROR: Voxel {} mass difference for clump set {} greater than 10%'.format(self.__index, ens+1))

    # clumpAfuv,interclumpAfuv = combinations.getAfuv()
    # clumpAfuv *= ensemble.CLmaxProbability
    # interclumpAfuv *= ensemble.ICmaxProbability

    self.__clumpVelocities = copy(ensemble.clumpVelocities)
    self.__clumpVelocityIndeces = copy(ensemble.clumpIndeces)   #all list entries should be the same

    # This gives an error if there are too many clumps in a line-of-sight; tau_FUV is too large for this equation...
    Afuv = combinations.getAfuv()
    self.__tauFUV = [-np.log((ensemble.CLmaxProbability[ens].prod(1)*Afuv[ens]).sum()) for ens in range(len(constants.clumpMassNumber))]
    
    return

  def getPosition(self):
    return np.array([self.__x, self.__y, self.__z])

  def getDensity(self):
    return self.__ensembleDensity

  def getMass(self):
    return self.__clumpMass

  def getClumpMass(self):
    return self.__clumpMass

  def getModelMass(self):
    return self.__modelMass

  def getVelocity(self):
    return self.__velocity, self.__ensembleDispersion

  def getClumpVelocity(self):
    return self.__clumpVelocities, self.__clumpVelocityIndeces

  def getFUVabsorption(self):
    return self.__tauFUV

  def calculateEmission(self, verbose=False):
    
    masspoints.calculateEmission(tauFUV=self.__tauFUV)
    combinations.calculateEmission()

    clumpIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
    clumpOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]
    clumpIntensityDust = [[] for _ in range(len(constants.clumpMassNumber))]
    clumpOpticalDepthDust = [[] for _ in range(len(constants.clumpMassNumber))]

    interclumpIntensity = []
    interclumpOpticalDepth = []

    iDust = constants.wavelengths[constants.nDust].size

    # Clump
    for ens in range(len(constants.clumpMassNumber)):

      vel = constants.velocityRange   #v_obs
      clumpVel = ensemble.clumpVelocities[ens]   #v_sys
      factor = np.exp(-(vel.reshape(1,-1)-clumpVel.reshape(-1,1)-self.__velocity)**2/2/constants.clumpDispersion**2)
      
      for i,probability in enumerate(ensemble.clumpProbability[ens]):

        intensity = copy(combinations.clumpIntensity[ens][:,iDust:])
        intensity = np.array([intensity*factor[ensemble.clumpIndeces[ens][i],j] for j in range(vel.size)])  #shape (v_obs, combination, wavelength)
        clumpIntensity[ens].append(np.array([(probability.prod(1)*intensity[j].T).T for j in range(vel.size)]))  #shape (v_sys, v_obs, combination, wavelength)
        intensityDust = copy(combinations.clumpIntensity[ens][:,:iDust])  #shape (combination, wavelength)
        clumpIntensityDust[ens].append((probability.prod(1)*intensityDust.T).T)  #shape (v_obs, combination, wavelength)

        opticalDepth = copy(combinations.clumpOpticalDepth[ens][:,iDust:])
        opticalDepth = np.array([opticalDepth*factor[ensemble.clumpIndeces[ens][i],j] for j in range(vel.size)])
        clumpOpticalDepth[ens].append(np.array([(probability.prod(1)*np.exp(-opticalDepth[j].T)).T for j in range(vel.size)]))
        opticalDepthDust = copy(combinations.clumpOpticalDepth[ens][:,:iDust])
        clumpOpticalDepthDust[ens].append((probability.prod(1)*np.exp(-opticalDepthDust.T)).T)

      # All of these have shape (velocity, wavelength)
      clumpIntensity[ens] = (np.array(clumpIntensity[ens]).sum(2)).sum(0).astype(np.float64)
      clumpOpticalDepth[ens] = (-np.log(np.array(clumpOpticalDepth[ens]).sum(2))).sum(0).astype(np.float64)
      clumpIntensityDust[ens] = np.array([np.array(clumpIntensityDust[ens]).sum(1).sum(0) for _ in range(vel.size)]).astype(np.float64)
      clumpOpticalDepthDust[ens] = np.array([-np.log(np.array(clumpOpticalDepthDust[ens]).sum(1)).sum(0) for _ in range(vel.size)]).astype(np.float64)

    # Sum over ensembles
    self.__intensity = np.append(clumpIntensityDust, clumpIntensity, axis=2)
    self.__opticalDepth = np.append(clumpOpticalDepthDust, clumpOpticalDepth, axis=2)

    if verbose:
      print('Voxel emission calculated.')
    
    return

  def getEmissivity(self):
    return self.__intensity/constants.voxel_size

  def getAbsorption(self):
    return self.__opticalDepth/constants.voxel_size

  def getIntensity(self):
    epsilon = self.getEmissivity()
    kappa = self.getAbsorption()
    intensity = epsilon/kappa*(1-np.exp(-kappa*constants.voxel_size))
    i_nan = np.isnan(intensity)
    intensity[i_nan] = epsilon[i_nan]
    return intensity

  def plotMolecule(self, molecule='', quantity='intensity', title='', logscale=False):
    # Plot the molecule emission with respect to the voxel velocity structure

    import matplotlib.pyplot as plt

    nDust = constants.wavelengths[constants.nDust].size

    if molecule in species.molecules and isinstance(molecule, str):
      molecule = [molecule]

    elif isinstance(molecule, list):
      pass

    else:
      molecule = species.molecules

    if quantity=='emissivity':
      value = self.getEmissivity()
      ylabel = r'$\epsilon_{\lambda} \ \left( \frac{K}{pc} \right)$'

    elif quantity=='absorption':
      value = self.getAbsorption()
      ylabel = r'$\kappa_{\lambda} \ \left( \frac{1}{pc} \right)$'

    elif quantity=='intensity':
      value = self.getIntensity()
      ylabel = r'$I_{\lambda} \ \left( K \right)$'

    fig,axes = plt.subplots(len(value), figsize=(10, 5*len(value)))

    for ens in range(len(value)):

      if isinstance(axes, np.ndarray):
        ax = axes[ens]

      else:
        ax = axes

      vel = constants.velocityRange#[self.__clumpVelocityIndeces[ens]]
      labels = []

      for mol in molecule:

        if not mol in species.molecules:
          print('Species {} not in model.'.format(mol))
          pass

        i = nDust + np.where(mol==np.asarray(species.molecules))[0][0]
        
        if logscale:
          ax.semilogy(vel, value[ens][:,i], marker='o', ms=2, ls='-', lw=1)
        else:
          ax.plot(vel, value[ens][:,i], marker='o', ms=2, ls='-', lw=1)
        labels.append(mol)

      ax.legend(labels=labels, fontsize=14, bbox_to_anchor=(1.05, 1))
      ax.set_xlabel(r'$V_r \ (\frac{km}{s})$', fontsize=20)
      ax.set_ylabel(ylabel, fontsize=20)
      plt.setp(ax.get_xticklabels(), fontsize=14)
      plt.setp(ax.get_yticklabels(), fontsize=14)

      if title:
        ax.set_title(title, fontsize=20)
      else:
        ax.set_title('Clump set {} {}'.format(ens+1, quantity), fontsize=20)
    
    fig.tight_layout()

    plt.show()

    return

  def plotSpectrum(self, quantity='intensity', title=''):
    # Plot the either the intesity or optical depth spectrum at the voxel velocity (the velocity with the largest ensemble).

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    moleculeColour = { 'C+'    : 'xkcd:orange', \
                      'C'      : 'xkcd:red', \
                      'O'      : 'xkcd:violet', \
                      'CO'     : 'xkcd:maroon', \
                      '13C+'   : 'xkcd:light blue', \
                      '13C'    : 'xkcd:blue', \
                      '13CO'   : 'xkcd:sapphire', \
                      'HCO+'   : 'xkcd:neon green', \
                      'H13CO+' : 'xkcd:forrest green', \
                      'H3O+'   : 'xkcd:vomit'
                    }

    constants.resortWavelengths()
    nDust = constants.wavelengths[constants.nDust].size

    if quantity=='emissivity':
      value = self.__intensity/constants.voxel_size
      ylabel = r'$\epsilon_{\lambda} \ \left( \frac{K}{pc} \right)$'

    elif quantity=='absorption':
      value = self.__opticalDepth/constants.voxel_size
      ylabel = r'$\kappa_{\lambda} \ \left( \frac{1}{pc} \right)$'

    elif quantity=='intensity':
      value = self.__intensity/self.__opticalDepth * (1-np.exp(-self.__opticalDepth))
      ylabel = r'$I_{\lambda} \ \left( K \right)$'

    fig,axes = plt.subplots(len(value), figsize=(10, 5*len(value)))

    for ens in range(len(value)):

      if isinstance(axes, np.ndarray):
        ax = axes[ens]

      else:
        ax = axes

      vel = int(value[ens].shape[0]/2)
      dustInterp = interp1d(constants.wavelengths[constants.nDust], value[ens][vel,:nDust], fill_value='extrapolate')

      ax.loglog(1e6*constants.wavelengths[constants.nDust], value[ens][vel,:nDust], ls='--', lw=0.5, c='xkcd:black')
      
      molecules = []
      colour = []

      for i,molecule in enumerate(species.molecules):
        valTemp = dustInterp(species.moleculeWavelengths[i])
        if not molecule.split()[0] in molecules:
          colour.append(moleculeColour[molecule.split()[0]])
          molecules.append(molecule.split()[0])
        ax.loglog([1e6*species.moleculeWavelengths[i]]*2, [valTemp,valTemp+value[ens][vel,nDust+i]], ls='-', lw=1, c=colour[-1])

      lines = [Line2D([0], [0], color=moleculeColour[molecule], lw=1) for molecule in molecules]
      # mols  = ['dust', ]
      leg = ax.legend(labels=molecules, handles=lines, fontsize=16, bbox_to_anchor=(1.05, 1))
      
      ax.set_xlabel(r'$\lambda \ (\mu m)$', fontsize=20)
      ax.set_ylabel(ylabel, fontsize=20)
      ax.set_ylim([value[ens][vel,:nDust].min()*10e-3, value[ens][vel,:nDust].max()*10e3])
      plt.setp(ax.get_xticklabels(), fontsize=14)
      plt.setp(ax.get_yticklabels(), fontsize=14)

      if title:
        ax.set_title(title, fontsize=20)
      else:
        ax.set_title('Clump Set {} {} spectrum'.format(ens+1, quantity), fontsize=20)
    
    fig.tight_layout()

    plt.show()

    return
    
  def printVoxel(self):
    # OUTDATED
    if self.__debugging:
      iClump,tauClump,FUVclump = self.__clump.getEnsembleEmission()
      iInterclump,tauInterclump,FUVinterclump = self.__interclump.getEnsembleEmission()
    names = self.__species[0].getMolecules() + self.__species[1].getDust()
    transitions = self.__species[0].getTransitions() + self.__species[1].getTransitions()
    totalIntensity = np.array([names, transitions, self.__intensity.sum(2).max(1)]).T
    clumpIntensity = np.array([names, transitions, iClump.sum(2).max(1)]).T
    interclumpIntensity = np.array([names, transitions, iInterclump.sum(2).max(1)]).T
    totalTau = np.array([names, transitions, self.__opticalDepth.sum(2).max(1)]).T
    clumpTau = np.array([names, transitions, tauClump.sum(2).max(1)]).T
    interclumpTau = np.array([names, transitions, tauInterclump.sum(2).max(1)]).T
    print('\nVoxel {}\n  ->Cartesian position: ({}, {}, {})'.format(self.__index, self.__x, self.__y, self.__z))
    print('  ->mass {} M_sol'.format(self.__mass))
    print('  ->intensity')
    for i in range(len(names)):
      print('    {} {}: {}'.format(totalIntensity[i][0],totalIntensity[i][1],totalIntensity[i][2]))
    if self.__debugging:
      print('    -clump')
      for i in range(len(names)):
        print('      {} {}: {}'.format(clumpIntensity[i][0],clumpIntensity[i][1],clumpIntensity[i][2]))
      print('    -interclump')
      for i in range(len(names)):
        print('      {} {}: {}'.format(interclumpIntensity[i][0],interclumpIntensity[i][1],interclumpIntensity[i][2]))
    print('  ->optical depth')
    for i in range(len(names)):
      print('    {} {}: {}'.format(totalTau[i][0],totalTau[i][1],totalTau[i][2]))
    if self.__debugging:
      print('    -clump')
      for i in range(len(names)):
        print('      {} {}: {}'.format(clumpTau[i][0],clumpTau[i][1],clumpTau[i][2]))
      print('    -interclump')
      for i in range(len(names)):
        print('      {} {}: {}'.format(interclumpTau[i][0],interclumpTau[i][1],interclumpTau[i][2]))
    print('  ->FUV field {}'.format(self.__FUV.getFUV()))
    return