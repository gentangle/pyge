#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Setup dot py."""

from __future__ import absolute_import, print_function

import os
from glob import glob
from os.path import dirname, join, relpath, splitext

from Cython.Build import cythonize
from setuptools import find_packages, setup
from setuptools.extension import Extension


def read(*names, **kwargs):
    """Read description files."""
    path = join(dirname(__file__), *names)
    with open(path, encoding=kwargs.get("encoding", "utf8")) as fh:
        return fh.read()


long_description = "{}\n{}".format(
    read("README.rst"),
    read("CHANGELOG.rst"),
)

setup(
    name="pyge",
    version="0.8.9",
    description="Python library to compute the Gaussian Entanglement",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    license="MIT License",
    author="Leonardo Salicari",
    author_email="leonardo.salicari@gmail.com",
    url="https://github.com/loscati/pyge",
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list:
        # http://pypi.python.org/pypi?%3Aaction=list_classifiers
        # 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: POSIX",
        "Operating System :: MacOS",
        "Operating System :: Microsoft",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    project_urls={
        "webpage": "https://github.com/loscati/pyge",
        "Documentation": "https://pyge.readthedocs.io/en/latest/",
        "Changelog": "https://github.com/loscati/pyge/blob/master/CHANGELOG.rst",
        "Issue Tracker": "https://github.com/loscati/pyge/issues",
        "Discussion Forum": "https://github.com/loscati/pyge/discussions",
    },
    keywords=["bioinformatics", "topology", "entanglement"],
    python_requires=">=3.9",
    install_requires=[
        # https://stackoverflow.com/questions/14399534
        "cython",
        "numpy",
        "numba",
        "mdanalysis",
        "biopython",
        "loguru",
    ],
    ext_modules=cythonize(
        [
            Extension(splitext(relpath(path).replace(os.sep, "."))[0], sources=[path])
            for root, _, _ in os.walk("pyge")
            for path in glob(join(root, "*.pyx"))
        ],
        # Important because it permits the correct compilation of GE cython routine
        language_level=3,
        annotate=True,
    ),
)
