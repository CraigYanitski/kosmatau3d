from kosmatau3d import models

# Edit this directory to point to the model you want to use.
directory = r'C:\Users\cyani\projects\pdr\KT3_history\MilkyWay\r250_cm1-1_d1_uv10'

# Calculate integrated intensity maps (adjust sl to change the number of sightlines you want in the map; [longitude, latitude]).
models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical', sl=[50,25], terminal=True, debug=False)