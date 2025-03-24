# Introduction

The Ultraviolet Imaging Telescope (UVIT) is one of the five payloads onboard the AstroSat space observatory.
This document provides instructions on how to run the UVIT Level 2 pipeline (hereafter referred to as UVPipe) on an input UVIT Level1 dataset.  UVPipe is developed in C++ and Python at the UVIT Payload Operations Centre (POC).

Please note that the UVIT POC will process all UVIT Level1 datasets using the UVPipe version available at the time. The processed Level2 datasets will be available in the AstroSat Archive.
> **IMPORTANT:** UVPipe was developed for execution at the UVIT POC. If you are running it independently, please verify the obtained results. Report any errors, warnings, or failures to the UVIT POC.

# Requirements

UVPipe has been tested exclusively on GNU/Linux operating systems. It requires a minimum of 64 GB of random-access memory (RAM). Storage requirements vary based on the size of the Level1 dataset.

