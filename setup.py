#!/usr/bin/python
# -*- coding: utf-8 -*-

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gwidecodeml",  # Replace with your own username
    version="1.1",
    author="Laura G. Macias",
    author_email="laugmacias@gmail.com",
    description="Testing positive selection in a genome-wide framework",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lauguma/GWideCodeML",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=['ete3', 'biopython', 'scipy'],
    include_package_data=True,
    python_requires='>=3.',
    entry_points={
        "console_scripts": [
            "gwidecodeml = gwidecodeml.__main__:main",
            "bonferroni_codeml = gwidecodeml.scripts.bonferroni_correction:main",
            "cds2concat = gwidecodeml.scripts.cds2concat:main",
            "fas2msa = gwidecodeml.scripts.fas2msa:main"
        ]
    }
)