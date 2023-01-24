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
"""
history = ""

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# TODO make sure to string format for markdown
setup(
    name='rna_map',
    version='0.3.0',
    description='',
    long_description=readme + '\n\n' + doclink + '\n\n' + history,
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    url='https://github.com/YesselmanLab/rna_map',
    packages=[
        'rna_map',
    ],
    package_dir={'rna_map': 'rna_map'},
    py_modules=[
        'rna_map/bit_vector',
        'rna_map/cli',
        'rna_map/exception',
        'rna_map/external_cmd',
        'rna_map/logger',
        'rna_map/mapping',
        'rna_map/parameters',
        'rna_map/run',
        'rna_map/sam',
        'rna_map/settings',
        'rna_map/util'
    ],
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='rna_map',
    classifiers=[
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: PyPy',
    ],
    entry_points = {
        'console_scripts' : [
            'rna-map = rna_map.cli : cli',
        ]
    }
)
