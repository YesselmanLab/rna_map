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

setup(
    name='dreem',
    version='0.2.0',
    description='',
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
            'dreem = dreem.run : main',
            'dreem-docker = dreem.run_docker:main',
            'dreem-multi = dreeem.run_multi:main'
        ]
    }
)
