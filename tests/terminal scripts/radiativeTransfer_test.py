from kosmatau3d import models

directory = r'C:\Users\cyani\projects\KT3_history\MilkyWay\r1000_cm1-1_d1_uv10'

# Calculate integrated intensity maps
models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical', sl=[100,50], terminal=True, debug=False)