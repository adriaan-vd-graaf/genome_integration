"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup
from setuptools import find_packages
from distutils.extension import Extension
# To use a consistent encoding



from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='genome_integration',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.1.0.dev0',
    python_requires='!=3.0.*, !=3.1.*, !=3.2.*,!=3.3.*, !=3.4.*, !=3.5.*, >=3.6.*, <4',

    description='Genome integration, a personal library',
    long_description=long_description,

    # The project's main homepage.
    url='none',

    # Author details
    author='Adriaan van der Graaf',
    author_email='adriaan.vd.graaf@gmail.com',

    # Choose your license
    license='Copyright 2018 Adriaan van der Graaf',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 1 - Work in progress',

        # Indicate who your project is intended for
        'Intended Audience :: Author',
        'Topic :: Genome analysis :: multiomics integration',

        # Pick your license as you wish (should match "license" above)
        'License :: Copyright 2019 Adriaan van der Graaf',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: >=3.6',
        'Programming Language :: Python :: 3.6'
    ],

    # What does your project relate to?
    keywords='multiomics',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    install_requires=['numpy', 'scipy', 'sklearn', 'statsmodels', 'plinkio', 'requests'],

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
