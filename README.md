# dnastorage
[![License: LGPL v3](https://img.shields.io/badge/License-LGPLv3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Build Status](https://travis-ci.com/dna-storage/dnastorage.svg?branch=master)](https://travis-ci.com/dna-storage/dnastorage)

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [License](#license)
- [Issues](https://github.com/dna-storage/dnastorage/issues)

# Overview

Core encoding, decoding, and file manipulation support for modeling DNA-based information storage systems.

# Documentation

As documentation for the software becomes available, it will be placed under the docs folder.

# System Requirements

## Hardware Requirements
The dnastorage module requires only a standard computer with enough RAM and compute power to support the needed operations. However, encoding or decoding large files may perform poorly, necessitating a more capable system.

## Software Requirements
### OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:

+ macOS: Catalina 10.15.3
+ Linux: Ubuntu 18.04.3

Note that most OSes will support our software by using Docker.

### Python Dependences

Our code has been tested on python versions 3.6 to 3.8. It has the following dependences:

```
nose
sphinx
editdistance
statistics
biopython
matplotlib
numpy
```

# Installation Guide

If you already have python 3 installed on your system, the simplest thing to do is download or checkout the code from GitHub.  Then, run the following commands:

    git clone https://github.com/dna-storage/dnastorage
    cd dnastorage
    
I recommend making a virtual environment first, but this is optional:

    python -m venv venv
    source venv/bin/activate
    
Then, use pip to install the requirements and packages:

    pip install -r requirements.txt
    pip install .

   
# License

This software is released under the LGPLv3 license.

# Acknowledgment

This work was supported by the National Science Foundation (Grants CNS-1650148, CNS-1901324, ECCS 2027655) and a North Carolina State University Research and Innovation Seed Funding Award (Grant 1402-2018-2509).
