#setup script for DIT software

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import glob

#version number
version = '0.01dev'

# Treat everything in scripts except README.rst as a script to be installed
scripts = glob.glob('scripts/*')
try:
	scripts.remove('scripts/README.rst')
except ValueError:
	pass


setup(name='specgrid',
	  description='Interpolation of stellar spectral grids',
	  author='Wolfgang E. Kerzendorf',
      version=version,
      packages=['specgrid'],
      package_data={'specgrid': ['data/*']},
      scripts=scripts
      )
      
