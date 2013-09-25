
from distutils.core import setup

setup(name='papi',
      version='1.0.0',
      description='PANIC data reduction pipeline',
      author='Jose M. Ibanez',
      author_email='jmiguel@iaa.es',
      url='https://github.com/ppmim/papi',
      license='GPLv3',
      packages=['papi'],
      requires=['numpy','scipy','pyfits','pyraf'],
      classifiers=[
       "Programming Language :: Python",
       'Development Status :: 3 - Alpha',
       "Environment :: Other Environment",
       "Intended Audience :: Science/Research",
       "License :: OSI Approved :: GNU General Public License (GPL)",
       "Operating System :: OS Independent",
       "Topic :: Scientific/Engineering :: Astronomy",
      ],
)
