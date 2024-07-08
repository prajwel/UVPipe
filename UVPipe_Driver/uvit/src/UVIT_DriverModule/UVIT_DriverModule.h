/* 
 * File:   UVIT_Drivermodule.h
 * Author: dhruv
 *
 * Created on October 3, 2016, 4:20 PM
 */

#ifndef UVIT_DRIVERMODULE_H
#define	UVIT_DRIVERMODULE_H
#include<string>
#include<stdlib.h>
#include<iostream>
#include<pil.h>
#include<fitsio.h>
#include<uvtUtils.h>
#include<algorithm>
#include<spMatrix.h>
#include<macro_def.h>
#include <uvtFullFrameAst.h>
using namespace std;
typedef struct Star 
{
    float intensity;
    int x, y;

};
class UVIT_DriverModule
{
public:
    int zipFlag;
    int UTC_flag;
 int readPILParameters(int argc,char **argv);
 int RunRelativeAspectChain();  
 int RunL2PCChain(string channel_l2,int cnt_num);
 int RunRAPCChain();
 int calculateShift(float *array);
 void doCentroiding (vector<int> &X , vector<int> &Y , int centroidwindow , float *arr , int h , int w);
 int findStar_algo1 (float *inputArray,int channel_id);
 int findShiftsNtheta (int totalelements , vector<float> &Xone , vector<float> &Yone , vector<float> &Xtwo , vector<float> &Ytwo ,
        vector<float> &DiffOfX , vector<float> &DiffOfY , bool flag_theta_computation,double &Xdx , double &Ydy , double &Theta);
// int  matchStars (int numrowsFirstfile , int numrowsSecfile , float divFact ,
//        float *xlocFirst , float *ylocFirst , float *xlocSec , float *ylocSec , vector<float> &matchPixelXone ,
//        vector<float> &matchPixelYone , vector<float> &matchPixelXtwo , vector<float> &matchPixelYtwo ,
//        vector<float> &matchPixelDiffX , vector<float> &matchPixelDiffY);
 int matchStars(int numrowsFirstfile, int numrowsSecfile, float divFact,
        float *xlocFirst, float *ylocFirst, float *xlocSec, float *ylocSec,double &DiffX,double &DiffY,bool &flag_TestTwo ,bool &flag_TestOne );
 int copyAllheaderKeys (char* infile);
 int removeRecords(vector<float>  &Xone , vector<float> &Yone ,vector<float> &Xtwo , vector<float> &Ytwo,vector<float> &DiffOfX , vector<float> &DiffOfY ,float *ints1,float*ints2,
    vector<float> &newXone,vector<float> &newYone,vector<float> &newXtwo,vector<float> &newYtwo,vector<float> &newDiffX,vector<float> &newDiffY,vector<float> &new_one_ints,vector<float> &new_two_ints);
 
 int ApplySubSampling (float* inputarray , int in_xsize , int in_ysize , float* outputarray , int out_xsize , int out_ysize);
 int crc_flag,crc_flagrapc,crc_flagpc,crc_flagpcfuv;
 
 vector<string> header_info;
 float *InputArray; //
    string modulename; //
    string level1indir; //Level-1 directory path from parameter file uvtRelativeAspectIM.par
    string level1indirpc;
    string caldbindir; //CalDB directory (CALDB/uvit) path from parameter file uvtRelativeAspectIM.par
    string level2outdir; //Output directory path from parameter file uvtRelativeAspectIM.par
    string channel,channelRAPC; //FUV / NUV / VIS
    string prev_Output_DirL2;
    vector<string> NUV_FAILS_track;
    int dropframe,darkframe_flag;
    int NUVonNUVflag;
    int FUVonFUVflag;
    int manualMode;
    int junkFrameFlag,gti_flag;
    int all_or_custom,valid_bit;
    int flatfieldFlag,padding_dim;
    int no_ofFramesToAcc,qe_mcpFlag,subdivisionFlag,subDivision_size,star_detect_algo_flag,refine_Winsize,
    centroid_Winsize,search_win_size,match_stars_file_flag,frames_toDiscard,nFrameToAverage,freqDomainFilter_Flag,type_Filtering,
    fitting_flag,orderpitch,orderyaw,orderroll,shift_rotation_algo,flag_thetaComp;
    float min_num_stars;
    double threshold,freqvalue,delta_time;
    float nbhd_dist;
    float backgrd_fact,rms_mul_factor_default,thrJunkFrame;
    int wtd_uc, wtd_bp, wtd_ff, wtd_pp, wtd_sd, wtd_cr, wtd_ac, wtd_fsc, wtd_dd, wtd_od, wtd_de, wtd_rfc, wtd_qemcp; //flags for writing outputs to disk
    bool last_FileFlag;
    int clobber, history; //for overwrite and history respectively in FITS file from parameter file 


    char mode[5]; //for PIL (query-learn(ql), auto(a), hidden(h))
    //for pipeline use
    bool flag_tarExtracted;
//for level2 PC 
    vector<float> Cx,Cy,Ci;//storing the centroid terms.
    vector<int> Fx,Fy,Rx,Ry;
    vector<float> Fval,Rval;
    vector<Star> star_track;
   
    char modepc[5];
    char rasfilepc[PIL_LINESIZE];
    int clobberpc, historypc;
    
    string level2outdirc; //Output directory path from parameter file uvtRelativeAspectPC.par
    
    //caldb_Handler caldb_handler; //Object for caldb_handler with methods to access different calibfile.fits file
  //  DataInfo datainfo; //object for DataInfo   
  //  float *frame_sig_Data, *frame_exp_Data; //stores frame pixel values for signal and exposure
  //  vector<string> header_info; //Stores Level-1 header information 
    int wtd_ucpc, wtd_bppc, wtd_ffpc, wtd_pppc, wtd_sdpc, wtd_crpc, wtd_fipc, wtd_fscp, wtd_ddpc, wtd_odpc, wtd_depc, wtd_rfcpc, wtd_qemcppc, wtd_centBiaspc, wtd_centCorrpc, wtd_snrpc, wtd_wmpc; //modulewise flags to decide whether to write to disk or not
    //char moduleoutdir_dd[NAMESIZE]; // Output directory path for detector distortion
    int nFrameIntegrate_fipc; //Number of frames to be integrated
    int nFrameDiscard_fipc; //Number of frames to be discarded in frame integration
    double star_time_ForDatasepc,end_time_ForDatasepc;
    
     long no_of_records ;    
//for dataIngest
    int parity_flagpc; //Parity flag (1-check parity for all 3 words, 2-check parity for x and y only) 
    int dropframepc; //yes for dropping frame if CRC fails, no for dropping packet if CRC fails
    
    int gti_flagpc;
    int all_or_custompc, valid_bitpc;
    
    //for Mask Bad pix
 
    float thr_multiphpc;
    int  flag_Thetapc;
    //for pix padding 


    //for sub division
   
    char rad_searchpc[PIL_LINESIZE];
    int fi_flagpc;
    //for cosmic Ray corrections              
   float sec_thr_Crpc; //for cosmic ray correction    
    float first_thr_Crpc; //Threshold for cosmic ray correction

   

    //for centroid Correction
    
    int centCorrflagpc,centBiasflagpc,DetectDistflagpc,OpticDistflagpc;
    ; //CalDB Effective Area Energy Filename
    int cent_corr_win_size; //Centroid correction  window sizes

    //for centroid bias
  
   
     int centroidlimitpc;
    int shift_N_Rotate_algopc;
    int star_detect_algo_flagpc, centroid_Winsizepc, refine_Winsizepc;
    float sd_mul_factorpc, sd_multi_factor_defaultpc;
    
    float diff_Distpc;
    int minimum_No_of_Starspc;
  
    int flag_thetaComppc;

    //for full frame Astorometry
    char catalogpathpc[PIL_LINESIZE];
    char att_timecolpc[NAMESIZE];
    char att_qcolpc[NAMESIZE];
    char databasenamepc[NAMESIZE];
    char outtarpathpc[NAMESIZE];
    int IMG_DIM_FIpc;
    int FINALFRAMESIZE_REGAVGpc;
    int search_algo_ctlgpc;     
    char len_apc[PIL_LINESIZE] ,len_bpc[PIL_LINESIZE];
    

   //for FUV
     char modepcfuv[5];
    char rasfilepcfuv[PIL_LINESIZE];
    int clobberpcfuv, historypcfuv;
    
    string level2outdircfuv; //Output directory path from parameter file uvtRelativeAspectPC.par
    vector<string> OutputDir_l2pcNUV,OutputDir_l2pcFUV;
    //caldb_Handler caldb_handler; //Object for caldb_handler with methods to access different calibfile.fits file
  //  DataInfo datainfo; //object for DataInfo   
  //  float *frame_sig_Data, *frame_exp_Data; //stores frame pixel values for signal and exposure
  //  vector<string> header_info; //Stores Level-1 header information 
    int wtd_ucpcfuv, wtd_bppcfuv, wtd_ffpcfuv, wtd_pppcfuv, wtd_sdpcfuv, wtd_crpcfuv, wtd_fipcfuv, wtd_fscpfuv, wtd_ddpcfuv, wtd_odpcfuv, wtd_depcfuv, wtd_rfcpcfuv, wtd_qemcppcfuv, wtd_centBiaspcfuv, wtd_centCorrpcfuv, wtd_snrpcfuv, wtd_wmpcfuv; //modulewise flags to decide whether to write to disk or not
    //char moduleoutdir_dd[NAMESIZE]; // Output directory path for detector distortion
    int nFrameIntegrate_fipcfuv; //Number of frames to be integrated
    int nFrameDiscard_fipcfuv; //Number of frames to be discarded in frame integration
    double star_time_ForDatasepcfuv,end_time_ForDatasepcfuv;
    
     long no_of_recordsfuv ;    
//for dataIngest
    int parity_flagpcfuv; //Parity flag (1-check parity for all 3 words, 2-check parity for x and y only) 
    int dropframepcfuv; //yes for dropping frame if CRC fails, no for dropping packet if CRC fails
    
    int gti_flagpcfuv;
    int all_or_custompcfuv, valid_bitpcfuv;
    
    //for Mask Bad pix
 
    float thr_multiphpcfuv;
    int  flag_Thetapcfuv;
    //for pix padding 


    //for sub division
   
    char rad_searchpcfuv[PIL_LINESIZE];
    int fi_flagpcfuv;
    //for cosmic Ray corrections              
    float sec_thr_For_Crfuv; //for cosmic ray correction    
    float first_thr_For_Crfuv; //Threshold for cosmic ray correction

   

    //for centroid Correction
    
    int centCorrflagpcfuv,centBiasflagpcfuv,DetectDistflagpcfuv,OpticDistflagpcfuv,qemcpflagpcfuv,qemcpflagpc;
    ; //CalDB Effective Area Energy Filename
    int cent_corr_win_sizefuv; //Centroid correction  window sizes

    //for centroid bias
  
   
     int centroidlimitpcfuv;
    int shift_N_Rotate_algopcfuv;
    int star_detect_algo_flagpcfuv, centroid_Winsizepcfuv, refine_Winsizepcfuv;
    float sd_mul_factorpcfuv, sd_multi_factor_defaultpcfuv;
    

    
    float diff_Distpcfuv;
    int minimum_No_of_Starspcfuv;
  
    int flag_thetaComppcfuv;

    //for full frame Astorometry
    char catalogpathpcfuv[PIL_LINESIZE];
    char att_timecolpcfuv[NAMESIZE];
    char att_qcolpcfuv[NAMESIZE];
    char databasenamepcfuv[NAMESIZE];
    char outtarpathpcfuv[NAMESIZE];
    int IMG_DIM_FIpcfuv;
    int FINALFRAMESIZE_REGAVGpcfuv;
    int search_algo_ctlgpcfuv;     
    char len_apcfuv[PIL_LINESIZE] ,len_bpcfuv[PIL_LINESIZE];
    
    //relativeaspect pc parameter 
   
    float err_perrapc;  //May be removed after checking
    int centroidlimitrapc; //May be removed after checking
    int no_ofFramesToAccrapc; //May be removed after checking
     
   float rms_mul_factor_defaultrapc; //May be removed after checking
      float *badExparrayrapc;//May be removed after checking
    char moderapc[5];       //for PIL (query-learn(ql), auto(a), hidden(h))
    int clobberrapc, historyrapc; //for overwrite and history respectively in FITS file from parameter file 
    
    string level1indirrapc; //Level-1 directory path from parameter file uvtRelativeAspectPC.par
    string level2outdirrapc; //Output directory path from parameter file uvtRelativeAspectPC.par
    string caldbindirrapc; //CalDB directory (CALDB/uvit) path from parameter file uvtRelativeAspectPC.par
    string channelrapc; //FUV / NUV / VIS
    
    char stardirrapc[10]; //Star directory path to store framewise star (first cut and refined peaks)
    char centroidDirrapc[10]; //Centroid directory path to store framewise star centroids
    int star_detect_algo_flagrapc; //Flag for SAC/PI or Joe's algorithm in star detection
    int centroid_Winsizerapc; //Window size for Star centroid computation
    int footprint_win_star_detectrapc; //Footprint window size in star detection for Joe's algorithm
    float primary_threshold_Valrapc; //Threshold in star detection for Joe's algorithm
    float secondary_threshold_Valrapc; //Threshold in star detection for Joe's algorithm
    
   int refine_Winsizerapc; //Window size for refined star peak determination in SAC/PI algorithm
   float rms_mul_factorrapc; //RMS factor for star detection in SAC/PI algorithm
   float sd_multi_factor_defaultrapc;  //Standard deviation factor for star detection in SAC/PI algorithm
    int subDivision_sizerapc; //Frame size after subdivision (e.g 9600)

    int frames_toDiscardrapc; //Number of frames to be discarded in reference frame calculation
    int nFrameToAveragerapc; //Number of frames to be averaged in reference frame calculation
    int FreqDomainFilterFlagrapc; //Frequency domain filter flag
    int fitting_flagrapc; //Polynomial fitting flag for roll, pitch and yaw
    int win_sizerapc;
    int orderrollrapc, orderpitchrapc, orderyawrapc; // Polynomial fitting order for roll, pitch and yaw
    
    int type_Filteringrapc; //1-Polynomial fitting & 0-Spatial domain based on cutt-off frequency
    double cr_thresholdrapc; //Threshold for cosmic ray correction
    double freqvaluerapc; //Cut-off frequency
    double poly_fit_intervalrapc; //Piecewise polynomial fitting interval (in seconds) for drift series
   

//    double t_darkframestart;	//Time for dark frame at end
//Refactor this to pair_nbhd_distance
    float pair_nbhd_distancerapc;  //Neighbourhood for star pairing/searching
//    double t_curr;  //Time of current frame
//    int xsize, ysize; //Width and Height of the frame
//    char RelAspFilename[FLEN_FILENAME]; //Path of relative aspect series file
////Refactor this to search_win_size
//    int win_size;	//Window size for searching stars in consecutive frames
//    int win_xsize,win_ysize; //Width and Height of the window (e.g (100, 100), (512, 512))
//Refactor dark_subtraction_flag
    int dark_subtraction_flagrapc; //Dark Subtraction flag
//Refactor to following
    //vector<int> refined_peaks_x_vect, refined_peaks_y_vect, firstcut_peaks_x_vect, firstcut_peaks_y_vect;
    //vector<int> refined_peaks_x_vect, refined_peaks_y_vect, firstcut_peaks_x_vect, firstcut_peaks_y_vect; //Stores x, y location of refined and firstcut peaks
//Refactor to following
    //vector<float> centroid_x_vect, centroid_y_vect, temp_centroid_x_vect, temp_centroid_y_vect; //Stores x, y location of final and temporary centroid 
    //vector<float> centroid_x_vect, centroid_y_vect, centroid_x_temp_vect, centroid_y_temo_vect;
//Refactor to following
    //vector<float> c_int, R_int, F_int, ci_tempi; //Stores intensity for centroid, refined peaks, firstcut peaks and temporary centroid
    //vector<float> c_int, R_int, F_int, ci_tempi;

//Flag for shift and rotation computation algorithm (1-SAC, 2-SAC, 3-IUCAA)
   int shift_rotation_algorapc;

    //vector<string> ref_frame_module_filename_track;
    int parity_flagrapc; //Parity flag (1-check parity for all 3 words, 2-check parity for x and y only)
    //vector<float> x_shift, y_shift, theta_shift; //Stores x-y shifts and rotations

//Refactor to following
    //vector<double> ref_frame_time_data; //Stores frame time for frames in reference frame calculation
    //vector<double> ref_frame_time_data;

    //float darkFrameend_data[darkSize*darkSize]; //Stores pixels of dark frame for end frame
    //float darkFramestart_data[darkSize*darkSize]; //Stores pixels of dark frame for start frame

    //Stores x,y and intensity values for all frames used for reference frame calculation
    //vector<float> x_Of_refFrame, y_Of_refFrame, int_Of_refFrame;
    //Stores average x,y and intensity values for all frames used for reference frame calculation
    //vector<float> x_ref_cumm, y_ref_cumm, int_ref_cumm; 
    //vector<float> roll_vect, pitch_vect, yaw_vect; 

//To be refined
    int gti_flagrapc; 
    int all_or_customrapc, valid_bitrapc;
    //char in_Gtifile[FLEN_FILENAME];


    //int cent_corr_win_size, cent_bias_win_size; 

    //for dataingest
    int dropframerapc; //yes for dropping frame if CRC fails, no for dropping packet if CRC fails
    
    //for unit conversion
    int unitConversionFlagrapc; //YES if unit conversion needs to run else NO

    //for flat field correction module
    int flatfieldFlagrapc; //YES if flat field correction is to be done, No if not required

    //for pixel padding
    int padding_dimrapc; //required dimension after padding

    //for Accumulate
    int nFrameIntegrate_firapc; //number of frames to integrate
    int nFrameDiscard_firapc; //Number of frames to be discarded in frame integration

    //for corrrection of temperature effects
    int qemcpFlagrapc; //YES to apply
     
    //for pixel subdivision      
    int subdivisionFlagrapc; //YES to apply
    
    int wtd_ucrapc, wtd_bprapc, wtd_ffrapc, wtd_pprapc, wtd_sdrapc, wtd_crrapc, wtd_firapc, wtd_fscrapc, wtd_ddrapc, wtd_odrapc, wtd_derapc, wtd_rfcrapc, wtd_qemcprapc, wtd_centBiasrapc, wtd_centCorrrapc;
    //float *darkCompute_array;
    int cmpr_framesrapc;
//    char eventfile[PIL_LINESIZE]; //name of input event file
//    char imgfile[PIL_LINESIZE]; //name 
//    float *darksubtr_array;
    int nCompareFramesrapc;
    float Threshold_crrapc;
    int IMG_DIM_FIrapc;
 int  flag_thetaComprapc;
 
    int minimum_No_of_Starsrapc;
    int match_stars_file_flagrapc;
    float thr_multiphrapc;
    int biasRowsrapc;
    float backgrd_factrapc;
    double *fraction_biasrapc, *x_corrrapc, *y_corrrapc;
    
     string orbnumrapc;
      bool paramfile_Varaskrapc,tar_extracted_flag_PCrapc;
       vector<string> sciencedata_filerapc,status_inforapc;
    
   
    
    

    
    
    
};


#endif	/* UVIT_DRIVERMODULE_H */

