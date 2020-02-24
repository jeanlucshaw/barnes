from distutils.core import setup, Extension

barneslib = Extension('barneslib', sources=['barneslib.cc'])

setup(name='barneslib',
      version= '1.0',
      description='Sources for Barnes interpolation adapted from Gri plotting software',
      ext_modules=[barneslib])
