#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = [
    'wheel>=0.23.0',
    'Cython>=0.20.2',
    'argparse>=1.2.1',
    'networkx>=2.6.3',
    'numpy>=1.21.2',
    'pandas>=1.3.5',
    'scikit_learn>=1.0.2',
    'scipy>=1.7.3',
    'sympy>=1.9',
    'tqdm>=4.65.0'
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='pcrp',
    version='0.1.0-dev0',
    description='pcrp - network based core circadian clock associated protein predictor.',
    long_description=readme + '\n\n' + history,
    author=['Vikram Singh', 'Vikram Singh'],
    author_email=['vikram.singh7571@gmail.com', 'vikramsingh@cuhimachal.ac.in'],
    url='https://github.com/vikramsinghlab/Pcrp',
    packages=[
        'pcrp',
    ],
    entry_points={'console_scripts': ['pcrp = pcrp.__main__:main']},
    package_dir={'pcrp':
                 'pcrp'},
    include_package_data=True,
    install_requires=requirements,
    license="GPLv3",
    zip_safe=False,
    keywords='pcrp',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',

        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
