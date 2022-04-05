from kosmatau3d import models
import numpy as np
params = {
        'voxel_size': 500,
        'molecules': 'all',
        'dust': ['240um', '550um'],
        'clumpMassRange': [[0, 2], [-2]],
        'clumpMassNumber': [3, 1],
        'clumpNmax': [1, 100],
        'velocityRange': [-350, 350],
        'velocityNumber': 701,
        'velocity': 108.8,
        'ensembleMass': [13058, 13058],
        'ensembleDensity': [2354, 1911],
        'ensembleDispersion': 4.75,
        'FUV': [1., 1.],
        'crir': 2e-16,
        'dilled': True
        }

vox = models.Voxel()
vox.setProperties(**params)
vox.calculateEmission()
print([np.shape(models.ensemble.clumpCombinations[ens]) for ens in range(models.constants.ensembles)])
vox.plotMolecule(molecule='CO 1')
