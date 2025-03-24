
# Introduction

The Ultraviolet Imaging Telescope (UVIT) is one of the five payloads onboard the AstroSat space observatory.
This document provides instructions on how to run the UVIT Level2 pipeline (hereafter referred to as UVPipe) on an input UVIT Level1 dataset.  UVPipe is developed in C++ and Python at the UVIT Payload Operations Centre (POC).

Please note that the UVIT POC will process all UVIT Level1 datasets using the UVPipe version available at the time. The processed Level2 datasets will be available in the AstroSat Archive.
> **IMPORTANT:** UVPipe was developed for execution at the UVIT POC. If you are running it independently, please verify the obtained results. Report any errors or failures to the UVIT POC.

# Requirements

UVPipe has been tested exclusively on GNU/Linux operating systems. It requires a minimum of 64 GB of random-access memory (RAM). Storage requirements vary based on the size of the Level1 dataset.

## Required Python packages

UVPipe requires Python version 3.9 or above. The required Python packages can be installed using `pip`. For example: ```pip install aafitrans```. All required Python packages are listed below.
* Aafitrans
* Astroalign
* Astropy
* Curvit
* Joblib
* Matplotlib
* NumPy
* Openpyxl
* Pandas
* Photutils
* Reproject
* Scipy
* Scikit-image
* XlsxWriter

## Other software requirements

* Astrometry.net
* Podman

# UVPipe setup

Ensure that all the required dependencies listed above are installed.   Once done, download the UVPipe source code from the [UVPipe GitHub Repository](https://github.com/prajwel/UVPipe)

To download the latest version of UVPipe, use the following link:
📥 **[Download UVPipe (Latest Version)](https://github.com/prajwel/UVPipe/archive/refs/heads/main.zip)**
