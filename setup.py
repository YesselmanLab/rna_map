#!/usr/bin/env python

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.md').read()
doclink = """
Documentation
-------------

The full documentation is at http://seq_tools.rtfd.org."""
history = ""

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='dreem',
    version='0.1.0',
    description='simple functions for manipulating sequences and secondary structures in pandas dataframe format',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    url='https://github.com/jyesselm/dreem',
    packages=[
        'dreem',
    ],
    package_dir={'dreem': 'dreem'},
    py_modules=[
        'dreem/bit_vector',
        'dreem/fastq',
        'dreem/logger',
        'dreem/mapping',
        'dreem/parameters',
        'dreem/run',
        'dreem/settings',
        'dreem/util'
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='seq_tools',
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    entry_points = {
        'console_scripts' : [
            'dreem = dreem.run:main',
        ]
    }
)
