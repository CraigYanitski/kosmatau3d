from setuptools import setup, find_packages

with open('README.md', 'r') as f:
  long_description = f.read()

setup(name = 'kosmatau3d',
      version = '1.0.3',
      description = 'package for the subgrid modelling of photo-dissociation regions using KOSMA-tau',
      url = 'https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d',
      author = 'C.N. Yanitski',
      author_email = 'yanitski@ph1.uni-koeln.de',
      long_description = long_description,
      long_description_content_type = "text/markdown",
      classifiers = ["Programming Language :: Python :: 3",
                     "License :: OSI Approved :: MIT License",
                     "Operating System :: OS Independent",
                     ],
      packages = find_packages(exclude=('docs','graphics','history','tests')),
      # external packages as dependencies
      install_requires = ['astropy',
                          'dill',
                          'jupyterlab',
                          'matplotlib',
                          'numba>=0.49.1',
                          'numpy',
                          'scipy',
                          'sphinx',
                          'sympy',
                          'tqdm',
                          ],
      include_package_data = True,
      package_data = {'': ['grid/*.dat', 'input/*/*.dat', 'molecular_data/*', 'preamble/*']},
      python_requires = ">=3.6"
)
