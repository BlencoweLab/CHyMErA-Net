import os.path
from setuptools import find_packages
from distutils.core import setup

here = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(here, 'chymeranet/version.py')).read())

setup(name='chymeranet',
      version=__version__,
      description='Deep learning framework for predicting efficient Cas12a guides',
      url='',
      author='Kevin Ha',
      author_email='k.ha@mail.utoronto.ca',
      license='',
      packages=find_packages(),
      package_data={'chymera': ['data/aux_scaler.human_mouse.full.pkl',
                                'data/CNN.human_mouse.full.h5']},
      install_requires=['numpy>=1.14.5',
                        'pandas>=0.23.4',
                        'biopython>=1.69',
                        'scikit-learn>=0.21.1',
                        'tensorflow>=1.13.1, <2',
                        'keras>=2.2.4',
                        'tqdm',
                        'joblib'],
      setup_requires=['pytest-runner'],
      tests_require=['pytest', 'pytest-mock'],
      python_requires='~=3.5',
      entry_points={
            'console_scripts': [
                'chymeranet = chymeranet:main',
            ]
      }
)
