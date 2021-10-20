from setuptools import setup
from setuptools import find_packages

# To use a consistent encoding

from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))


# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

name = 'genome_integration'
version = '1.0'
release = '1.0'

setup(
    name=name,

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=version,
    python_requires='!=3.0.*, !=3.1.*, !=3.2.*,!=3.3.*, !=3.4.*, !=3.5.*, >=3.6.*, <4',

    description='Genome integration, a personal library',
    long_description=long_description,

    # The project's main homepage.
    url='none',

    # Author details
    author='Adriaan van der Graaf',
    author_email='adriaan.vd.graaf@gmail.com',

    # Choose your license
    license='Copyright 2019 Adriaan van der Graaf',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Investigators interested in causal inference',
        'Topic :: Genome analysis :: multiomics integration',

        # Pick your license as you wish (should match "license" above)
        'License :: MIT licence 2020',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: >=3.7',
        'Programming Language :: Python :: 3.7'
    ],

    # What does your project relate to?
    keywords='multiomics',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    install_requires=['numpy', 'scipy', 'sklearn', 'statsmodels', 'plinkio', 'requests', 'bitarray', 'bgen_reader'],

    extras_require={
        'dev': ['check-manifest'],
        'test': ['coverage'],
    },

    setup_requires=['pytest-runner'],
    tests_require=['pytest'],


    entry_points={
        'console_scripts': [
            'genome_integration=genome_integration:main',
        ],
    },
    include_package_data=True,
)
