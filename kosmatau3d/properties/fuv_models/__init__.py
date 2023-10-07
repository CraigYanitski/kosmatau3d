'''
A subpackage to compare the far-UV spectral shapes used in
the literature to the one used in KOSMA-tau.
'''

from .models import *


def help():
    print('This package contains various models for the FUV energy density spectrum. ' +
          'All methods take a wavelength argument in Angstroms and return the energy density (lambda u_lambda).')
    print('  - u_habing: The original and oft-cited FUV description for interstellar radiation (Habing 1968).')
    print('  - u_draine: A refined FUV description equal roughly to 1.71 times the Habing field (Draine 1978). ' +
          'This method has an optional chi kwarg which can be used to scale the spectrum')
    print('  - u_mezger: An approximation using a broken power-law spectrum (Mezger et al. 1983).')
    print('  - u_zucconi: Typically this should be a superposition of damped blackbodies, though the method ' +
          'features just one blackbody (Zucconi et al. 2003). This requires a ´temp´ kwarg in Kelvin and a ' +
          'dimensionless ´scale´ kwarg')
    print('  - u_kosma: This is the spectrum utilised in KOSMA-tau, which uses a base Draine spectrum with an ' +
          'enhanced energy density for wavelengths >2000 Angstrom (Mezger et al. 1983). This may also use the ' +
          '´chi´ kwarg')
    print('\n\nComparison\n')
    compare_models()
    plot_comparison((912, 2066))

    return
