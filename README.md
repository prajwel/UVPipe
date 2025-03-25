
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

# Level1 dataset format

UVPipe will exclusively run on Level1 files in gzipped tarball format. If the Level1 dataset is in ZIP format, first convert it to gzipped tarball format using the `Level1_zip_to_gzipped_tarball.sh` script (found inside the `misc/` directory). For example:

``` bash
bash Level1_zip_to_gzipped_tarball.sh LEVL1AS1UVT20160322T01_051T01_9000000390_02613.zip
```

# UVPipe setup

Ensure that all the required dependencies listed above are installed. Next, download the UVPipe source code from the [UVPipe GitHub Repository](https://github.com/prajwel/UVPipe).

To download the latest version, use the following link:  
📥 [Download UVPipe (Latest Version)](https://github.com/prajwel/UVPipe/archive/refs/heads/main.zip).

### Extracting the files
Unzip the downloaded file. This will create a directory named `UVPipe-main`, which has the following structure:

```bash
UVPipe-main
├── LICENSE.txt
├── misc/
├── README.md
├── UVPipe_Driver/
└── UVPipe_Pilot/
```

- `UVPipe_Driver/` → Contains the C++ scripts, which should be run **first**.
- `UVPipe_Pilot/` → Contains the Python scripts, which should be run **after** the C++ scripts.
- `misc/` → Contains miscellaneous useful scripts.

### How to run the UVPipe C++ scripts

> **IMPORTANT:** The UVPipe C++ scripts are executed within a container using Podman.

Ensure that the Level1 dataset is placed in a working directory. For example, assume the working directory is: `/home/prajwel/where_Level1_data_is_kept/`. Follow the below steps to build and run the container image using Podman.

1. Open a terminal and navigate to the `UVPipe_Driver` directory. Then, set a temporary directory (TMPDIR). For example:
``` bash
export TMPDIR=/home/prajwel/podman_tmp
```

2. Build the container image using the following command:
``` bash
podman build -t uvpipe_driver .
```

3. Keep the appropriate `UVIT_DriverModule.par` file alongside the Level1 dataset at `/home/prajwel/where_Level1_data_is_kept/`. Commonly used variations of the `UVIT_DriverModule.par` file can be found at `UVPipe_Driver/Driver_module_param_files`.

4. Run the `uvpipe_driver` container with the following command:
``` bash
podman run --rm -it -v /home/prajwel/where_Level1_data_is_kept:/app/work_area:Z uvpipe_driver
```

> Note 1: Unless you make any changes to the contents of the `UVPipe_Driver` directory, steps 1 and 2 needs to be executed only once. The image, once built, will be available until it is removed.

> Note 2: To reset the Podman environment and clean up all builds, use the following command:
> ``` bash
> podman system reset
> ```

### How to run the UVPipe Python scripts

> **IMPORTANT:** Before running the UVPipe Python scripts, ensure that the UVPipe C++ scripts have already been executed on the Level1 dataset, as described in the previous section.

Follow the steps below to set up and execute the Python scripts:

1. Install and Configure Astrometry.net
    - Ensure that the `astrometry.net` package is installed.
    - Download the Astrometry.net index files specifically prepared for UVIT Astrometry from:  
    📥 [Download UVIT Astrometry Index Files (ZIP)](https://zenodo.org/records/12684908/files/UVIT_astrometry.zip?download=1)
    - Extract the downloaded ZIP file into a directory of your choice.
    - Open the `astrometry.cfg` file located in the `UVPipe_Pilot` directory and configure `astrometry.net` to use the downloaded index files by specifying the `add_path` parameter. For example:
        ``` bash
        ...
        # In which directories should we search for indices?
        add_path /home/prajwel/UVIT_astrometry/GAIA_epoch_2000_index_files_for_UVIT_astrometry
        add_path /home/prajwel/UVIT_astrometry/5200
        add_path /home/prajwel/UVIT_astrometry/6100
        add_path /home/prajwel/UVIT_astrometry/6000
        ...
        ```

2. Prepare the Working Directory
    - Copy all files from the `UVPipe_Pilot` directory to the working directory containing the Level1 dataset.
    - Ensure that the products generated by the UVPipe C++ scripts are present in this directory.

3. Run the UVPipe Python Scripts
    - Navigate to the working directory and execute the following command:
      ``` bash
      bash UVIT_pilot.sh
      ```

# Checking the UVPipe run

UVPipe includes multiple built-in checks at various stages of execution. However, users should primarily focus on the following key checks:

1. **`combining.log`**
   - This file is generated after each run. Open it and check for any failures or errors.
   - In the final sections of the `combining.log` file, `R` values are listed for two bright sources for each channel-filter-window combination.
   - Ideally, the `R` value should be less than 0.3.
   - Refer to `UVPipe_Pilot/README.txt` for further details.

2. **`_animation.gif` Files**
   - Two `_animation.gif` files will be generated, corresponding to two bright sources for each channel-filter-window combination.
   - Verify that the source centre consistently falls inside the black circle displayed in the animation.
