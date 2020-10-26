#!/usr/bin/env python
# -*- coding: utf-8 -*-

import importlib
import csv
import os
import sys
from setuptools import setup, find_packages
import io
import versioneer

# Constants
DISTNAME = 'DLSuR'
DESCRIPTION = 'Dynamic light scattering microrheology data analysis package'
URL = 'https://github.com/PamCai/DLSuR'
LICENSE = 'BSD-3-Clause'


# Parse description
with io.open('README.md', encoding='utf8') as f:
    README = f.read().split('\n')
    LONG_DESCRIPTION = '\n'.join([x for x in README if not x[:3]=='[!['])


VERSION = versioneer.get_version()
CMDCLASS = versioneer.get_cmdclass()


# Parse requirements.txt
with io.open('requirements.txt', 'r') as f:
    REQUIREMENTS = f.read().split('\n')


setup_info = dict(
    # Metadata
    name=DISTNAME,
    version=VERSION,
    author=('Brad Krajina, Pam Cai,'),
    author_email=('bradkrajina@gmail.com, pamela.5.cai@gmail.com, '),
    url=URL,
    classifiers=['Intended Audience :: Education',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: BSD License',
                 'Natural Language :: English',
                 'Operating System :: MacOS',
                 'Operating System :: POSIX :: Linux',
                 'Programming Language :: Python :: 3.5',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.7',
                 'Programming Language :: Python :: 3.8',
                 'Programming Language :: Python :: 3.9',
                 'Topic :: Scientific/Engineering :: Information Analysis',
                 'Topic :: Scientific/Engineering :: Mathematics',
                 'Topic :: Scientific/Engineering :: Physics',
                 ],
    cmdclass=CMDCLASS,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    python_requires='>=3.5',
    license=LICENSE,
    packages=find_packages(exclude=('test',)),
    install_requires=REQUIREMENTS,
    zip_safe=True,
)

setup(**setup_info)