from setuptools import setup, find_packages

setup(
    name='kosmatau3d',
    version='0.1',
    description='package for the modelling of photo-dissociation regions',
    author='Craig Yanitski',
    author_email='yanitski@ph1.uni-koeln.de',
    packages=find_packages(),  # same as name
    # external packages as dependencies
    install_requires=['numpy>=1.18.5',
                      'numba>=0.49.1',
                      'astropy>=4.0',
                      'scipy>=1.4.1',
                      'sympy',
                      'matplotlib>=3.1.2',
                      'tqdm>=4.36.1',
                      'jupyterlab'],
    include_package_data=True,
    package_data={'': ['grid/*.dat', 'input/*/*.dat', 'molecular_data/*', 'preamble/*']}
)
