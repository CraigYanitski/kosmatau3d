from kosmatau3d import models
import numpy as np

# Use these to change the voxel properties, which is necessary in the default mode.
parameters = { \
                # Model parameters
                'voxel_size' : 1, \
                 'molecules' : 'all', \
                      'dust' : 'PAH', \
            'clumpMassRange' : [[0, 2], [-2]], \
           'clumpMassNumber' : [3, 1], \
                 'clumpNmax' : [1, 100], \
             'velocityRange' : [-10, 10], \
            'velocityNumber' : 201, \

                # Voxel properties
                  'velocity' : 0., \
        'ensembleDispersion' : [2, 5], \
                 'clumpMass' : [1.5e5, 1e2], \
              'volumeFactor' : [0.9, 0.1], \
           'ensembleDensity' : 1e6, \
                       'FUV' : 1e4 \
  }

vox = models.Voxel()
vox.setProperties(**parameters, timed=True)
vox.calculateEmission(verbose=True, timed=True)
vox.plotMolecule(molecule=['C+ 1', 'C 1', 'CO 1', '13C+ 1', '13C 1', '13CO 1'], quantity='intensity', logscale=False)