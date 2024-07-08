/* 
 * File:   macro_def.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on Novemeber 19, 2014, 5:05 PM
 */

#ifndef MACRO_DEF_H
#define	MACRO_DEF_H


#define MKF_HDU_NAME "MKF"
#define HK_HDU_NAME "UVIT_HK"
#define PACKET_SIZE 2048                          //number of ytes in a packet
//#define DATA_INGEST_COLUMNS 8
#define BYTES_PER_EVENT 6                     //bytes for each event in data packet
#define CENTROID_BYTES 2016
#define GTI_BYTES 9
#define BUFFERSIZE 100000
#define PIXELNO 1008
#define PIXEL_DARK 48
#define MAX_EVENTS_PER_PACKET 336
#define EVENTDATA "EVENT_DATA"
//#define RAWFRAMESIZE  512
#define SCIENCEDATA_HDUNAME_FIRST   "DETECTOR_SETTING"
#define SCIENCEDATA_HDUNAME   "DETECTOR_DATA"
#define TCT_FUV_COLNO 3
#define TCT_NUV_COLNO 4
#define TCT_VIS_COLNO 5
#define NO_OF_BITS_PER_BYTE 8
#define TTYPE_KEYWORDNAME_SKIP 5
#define MONITOR2_BYTES 48
#define GTI_FILE_NAME  "GTIParams.txt"
//#define FINALFRAMESIZE_REGAVG 4800
//#define FINALFRAMESIZE_REGAVG 600
#define PIX_PADD_SIZE  600
#define SIZE_FOR_BUFFER 512
#define DARK_OUTPUT_DIR "Dark"
#define ARC_PERSEC_TO_DEGREE 3600
#define PERCENTGE_ERR_ALLOWED 30
#define SIGMA_FAC_DARK  1
#define EXE_FOR_PCMODE uvt_pc_ra
#define EXE_FOR_IMMODE  uvt_im_ra
#define CALDB_DARK_DIR DARK
#define FOCAL_LENGTH 4.6378
#define DEG_TO_RAD M_PI 180.0
#define PIX_PER_DEG_NUV 30/(512*60)
#define PIX_PER_DEG_FUV 30/(512*60)
#define PIX_PER_DEG_VIS 30/(512*60)
#define darkSize 512                    
//#define IMG_DIM_FI 600
//#define IMG_DIM_FI 9600
//#define IMG_DIM_FI 4800
#define PADD_DEFAULT 44
#define INSIDE_TEMP_NUV  11
#define INSIDE_TEMP_FUV 3
#define INSIDE_TEMP_VIS 19
#define OUTSIDE_TEMP_NUV 17
#define OUTSIDE_TEMP_FUV 9
#define OUTSIDE_TEMP_VIS  25
#define MAX_FIRSTCUTPIX_CMPR  500
#define MAX_REFINED_PIX_SIZE 20
//#define COLOR_RED  \033[1;34m

//typedef struct FrameIntegration_Arr {
//    float img_pixels_sig[IMG_DIM_FI * IMG_DIM_FI];
//    float Frame_pixels_exp[IMG_DIM_FI * IMG_DIM_FI];
//    int frameNo_fi;
//    double frameTime_fi;
//};
#endif	/* MACRO_DEF_H */

