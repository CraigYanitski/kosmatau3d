#python3.7

from kosmatau3d import models

timed = False
debug = False

print('KOSMA-tau^3')

# Edit these parameters according to the model you want to produce.
parameters = {
                # Model information
                'history_path': r'\Users\cyani\projects\pdr\KT3_history',
                'directory': r'\MilkyWay',
                'x': 36,
                'y': 36,
                'z': 2,
                'modelType': 'disk',
                
                # Model parameters
                'resolution': 500,
                'molecules': 'all',
                # 'molecules' : ['C+ 1', 'C 1', 'CO 1'],
               'dust' : 'PAH',
                # 'dust': 'molecular',
                'clumpMassRange': [[0, 2], [-2]],
                'clumpMassNumber': [3, 1],
                'clumpNmax': [1, 100],
                'velocityRange': [-300, 300],
                'velocityNumber': 1000,
                
                # Property factors
                'clumpMassFactor': [1, 1],
                'FUVfactor': 1,
                'densityFactor': 1,
                'globalUV': 10
              }

kosma = models.Model(**parameters)

kosma.calculateModel(timed=timed, debug=debug)
# kosma.writeEmission(debug=debug)
