
Astrosat/UVIT Level2 data products (README version: 1.0.1)
==========================================================

Maintained by UVIT Payload Operations Center (UVIT POC), IIA
Email: uvit.atc@iiap.res.in
Website: https://www.iiap.res.in/projects/uvit/


CONTENTS
--------
1. Definitions of "combine" and "coadd"
2. Level2 data archive
3. Directory structure and files in an extracted Level2 archive
4. Additional information


1. Definitions of "combine" and "coadd"
---------------------------------------
The terms "combine" and "coadd" have special meanings in the UVIT Level2 data. Their definitions are given below.

* Frame refers to ~1 second integration mode (IM) frames in VIS, and photon counting mode (PC) frames in UV channels (PC frames are at ~29 fps for the full window mode).

* Coadded image refers to episode-wise images that are created by combining frames in an episode. The term “coadd” refers to this process.

* Combined image refers to multi-episode images that are created by combining multiple coadded images or by combining frames from multiple episodes. The term “combine” refers to this process.


2. Level2 data archive
----------------------
An example of a UVIT Level2 output file would be "LEVL2AS1UVT20160101A01_123T01_0123456789.tar_V2.7". The file naming convention is summarized below:
char 01-05 - LEVL2, denotes Level2 data product.
char 06-08 - AS1, for Astrosat 1.
char 09-11 - 3 character instrument code, UVT for UV Imaging Telescope.
char 12-19 - Date of observation (YYYYMMDD), as listed in Level1 Data.
char 20-26 - Proposal ID, a unique identification string for each proposal. Char 20-22 is the code for the observing cycle.
char 20-40 - Observation ID, a unique identification string for each pointing.
char 41-44 - File extension, UNIX/Linux tar format.
char 45-49 - File Version information. The major version number changes with the version of the Level 1 set. The minor version will be ".7".


3. Directory structure and files in an extracted Level2 archive
---------------------------------------------------------------
The extracted archive will have the following directory structure (using the example filename from Section 1). Please note that this is an example structure, and the number of files may vary depending on the Level2 data generated from Level1 datasets.

20160101_A01_123T01_0123456789_level2/
└── uvit
    ├── data_products
    │   ├── FUV_data
    │   │   ├── AS1A01_123T01_0123456789uvtFIIPC00F4A_l2err.fits
    │   │   ├── AS1A01_123T01_0123456789uvtFIIPC00F4A_l2exp.fits
    │   │   ├── AS1A01_123T01_0123456789uvtFIIPC00F4I_l2err.fits
    │   │   ├── AS1A01_123T01_0123456789uvtFIIPC00F4I_l2exp.fits
    │   │   ├── AS1A01_123T02_0123456789uvtFIIPC00F4A_l2img.fits
    │   │   ├── AS1A01_123T02_0123456789uvtFIIPC00F4I_l2img.fits
    │   │   └── AS1A01_123T02_0123456789uvtFIIPC00F4_l2ce.fits
    │   ├── NUV_data
    │   │   ├── AS1A01_123T01_0123456789uvtNIIPC00F4A_l2err.fits
    │   │   ├── AS1A01_123T01_0123456789uvtNIIPC00F4A_l2exp.fits
    │   │   ├── AS1A01_123T01_0123456789uvtNIIPC00F4I_l2err.fits
    │   │   ├── AS1A01_123T01_0123456789uvtNIIPC00F4I_l2exp.fits
    │   │   ├── AS1A01_123T02_0123456789uvtNIIPC00F4A_l2img.fits
    │   │   ├── AS1A01_123T02_0123456789uvtNIIPC00F4I_l2img.fits
    │   │   └── AS1A01_123T02_0123456789uvtNIIPC00F4_l2ce.fits
    │   └── VIS_data
    │       ├── AS1A01_123T01_0123456789uvtVIIIM00F2A_l2ql.fits
    │       └── AS1A01_123T01_0123456789uvtVIIIM00F2I_l2ql.fits
    ├── DISCLAIMER.txt
    ├── LEVL1AS1UVT20160101A01_123T01_0123456789_03478_V2.2_dqr.xml
    ├── README.txt
    └── subsidiary_data_and_information
        ├── combining_process_information
        │   ├── astrometry_error_vector_plots
        │   ├── curvit_plots
        │   └── logs
        └── episode_data
            ├── driver_module
            ├── uvt_01
            │   ├── F_01
            │   ├── N_01
            │   └── V_01
            └── uvt_02
                ├── F_02
                ├── N_02
                └── V_02

The explanation of the Level2 data structure is given below:

* The parent 20160101_A01_123T01_0123456789_level2 directory contains a single subdirectory called "uvit".

* The uvit directory contains "data_products" and "subsidiary_data_and_information" subdirectories.

* The data_products directory contains "FUV_data", "NUV_data", and "VIS_data" subdirectories. They contain the UVIT Level2 products.

* The UV combined data ([F/N]UV_data) directories contain count-rate, count-rate error, and exposure maps in Instrument and astronomical coordinate systems. A table listing time-ordered UV events from all contributing Episodes (events list) is also provided.
	- The count-rate images have the suffix *l2img.fits.
	- The count-rate error images have the suffix *l2err.fits.
	- The exposure maps have the suffix *l2exp.fits.
	- The events lists have the suffix *l2ce.fits.
	- The instrument coordinate system data have "I_" in their filenames.
	- The astronomical coordinate system data have "A_" in their filenames.

* The VIS combined data (VIS_data) directory contains the quick-look images in Instrument and astronomical coordinate systems.
	- The instrument coordinate system data have "I_" in their filenames.
	- The astronomical coordinate system data have "A_" in their filenames.
	- No photometric calibration exists for the VIS channel data.

* All Level2 products have prefixes following the same convention (example: AS1A01_123T01_0123456789uvtFIIPC00F4). The convention is described below:
char 01-03 - AS1, for Astrosat 1.
char 04-10 - Proposal ID, a unique identification string for each proposal. Char 04-06 is the code for the observing cycle.
char 04-24 - Observation ID, a unique identification string for each pointing.
char 25-27 - 3 character instrument code, uvt for UV Imaging Telescope.
char 28    - Channel (F for FUV, N for NUV or V for VIS).
char 29-30 - Data Source.
char 31-32 - Mode (PC for Photon Counting or IM for Integrating Mode).
char 33-34 - Window
			------------------------------------------
			W W     	field size (pixel x pixel)
			------------------------------------------
			0 0 for 	512 X 512 -- full field
			2 2 for 	100 X 100
			3 3 for 	150 X 150
			4 4 for 	200 X 200
			5 5 for 	250 X 250
			6 6 for		300 X 300
			7 7 for		350 X 350
			------------------------------------------
char 35-36 - Filter wheel setting.

* The subsidiary_data_and_information directory contains "combining_process_information" and "episode_data" subdirectories.

* The combining_process_information directory contains "astrometry_error_vector_plots", "curvit_plots", and "logs" subdirectories.

* The astrometry_error_vector_plots directory contains astrometry error vector plots.

* The curvit_plots directory contains Curvit-generated plots for two bright sources. The plots contain source light curves and are part of the UV episode combining goodness inspection.

* The logs directory contains the combining process log as a text file.

* The episode_data directory contains the episode data. It has directories called "driver_module" and "uvt_01", "uvt_02", etc.

* The driver_module directory contains one text file with the UVIT driver module parameters set for that particular run.

* The uvt_01 directory, for example, contains F_01, N_01, and V_01 directories. They contain FUV, NUV, and VIS data, respectively.

* The UV episode data directory contains the count-rate and exposure maps in instrument coordinates, the count-rate image in astronomical coordinates, and the events list.

* The VIS episode data directory contains the VIS drift series and the coadded VIS image.


4. Additional information
-------------------------
A) Combining methods:
The Level2 pipeline provides three methods of combining UV episode data. The method used is denoted in the header with the COMBMETH keyword. The COMBMETH keyword can have three values:
	i) DEFAULT_METHOD - the default method. The method uses a correlation technique. All data is used.
	ii) ALT_METHOD1 - first alternate method. The method uses a correlation technique. In this method, the background pixels are first removed from the data before the correlation technique is applied.
	iii) ALT_METHOD2 - second alternate method. The method uses a correlation technique. Only the data of a bright source is used when the correlation technique is applied.

B) Curvit-generated plots:
Curvit is an open-source Python package that facilitates the creation of light curves from UVIT data. The software documentation is available at https://curvit.readthedocs.io/en/latest/. We have used Curvit's capabilities to check the goodness of the UV episodes' combining process. Light curves are generated for two bright sources in each channel-filter-window combination for goodness checks. The light curve plots are visually inspected to identify poorly combined episodes. Please note that the errors given for each data point in the light-curve plots are one standard deviation error.

Apart from visual checks of the Curvit plots, a variability measure (called "R") is used to identify bad combining. Let C be an array containing the episode-wise count-rate values, c_n. Then
R = (MAX(C) – MIN(C)) /  MEAN(C)

The measure is sensitive to even a single misaligned episode. While all datasets are inspected for quality, UVIT POC specially inspects those datasets with R values above or equal to 0.3 (corresponds to a 30% fall in count-rate from the maximum to the minimum count-rate with respect to the mean count-rate). A 1.9-pixel shift gives a flux offset of ~30%. Please note that, even if the combining was good, a >30% fall is possible due to:
	i) intrinsic variation of the astrophysical sources,
	ii) variations in the goodness of drift correction and
	iii) edge effects in episode-wise data.
Therefore, the UVIT POC checks two bright sources in each channel-filter-window combination.

C) Source encircled energy ratio:
The source encircled energy ratio (SEER) is the ratio between fluxes measured with 2.5 and 5 subpixel aperture radii centred on a source. The SEER values are estimated for two bright sources in the field of view. The estimated SEER values (SEER_1 and SEER_2) will be given in the header for each UV product. The SEER values are expected to be 0.70 for the NUV channel and 0.68 for the FUV channel (Tandon et al., 2020).
