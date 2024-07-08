/* 
 * File:   uvtPcRaCommanArr.h
 * Author: uvit
 *
 * Created on April 9, 2015, 2:40 PM
 */

#ifndef  UVTLEVEL2PC_H
#define UVTLEVEL2PC_H

#include<string>
#include<iostream>
#include<stdlib.h>
#include<uvtUtils.h>
#include<caldb_Handler.h>
#include<DataInfo.h>
#include<macro_def.h>
#include<uvtFullFrameAst.h>
using namespace std;
#define darkSize 512
//#define IMG_DIM_FI 600

//class uvtLevel2PC;

typedef struct FrameIntegration_Arr_strct {
    
    //vector<float> img_pixels_sig;//[IMG_DIM_FI * IMG_DIM_FI];
    vector<float>  Frame_pixels_exp ;//[IMG_DIM_FI * IMG_DIM_FI];
    int frameNo_fi; 
    double frameTime_fi;
};
//typedef struct FrameIntegration_Arr_strct {
//    
//    float *img_pixels_sig;
//    float *Frame_pixels_exp;
//    int frameNo_fi; 
//    double frameTime_fi;
//};
typedef struct New_location{
    
    float diff_Distance,X_locations,y_locations,value_intensity;
    
};
typedef struct Star 
{
    float intensity;
    int x, y;

};

class uvtLevel2PC 
{
    int zipFlag;
    float *InputArray;
    int UTC_flag;
    double integrationTime_frmUTC,bzero_mjd,bscale_mjd;
double integrationTime_frm_curr_option;
long total_InputPacket,start_pktNo,usedpkts;
    string modulename;
int crc_flag;
int cntTotalFrameInput;
long crcfailedpckts,parityfailedEvents,pcktremained,eventsremained,pktsFromL1;
double CrcFactor,ParityFactor;
    char mode[5];
    char rasfile[PIL_LINESIZE];
    int clobber, history;
    string level1indir; //Level-1 directory path from parameter file uvtRelativeAspectPC.par
    string caldbindir; //CalDB directory (CALDB/uvit) path from parameter file uvtRelativeAspectPC.par
    string level2outdir; //Output directory path from parameter file uvtRelativeAspectPC.par
    string channel; //FUV / NUV / VIS
    caldb_Handler caldb_handler; //Object for caldb_handler with methods to access different calibfile.fits file
    DataInfo datainfo; //object for DataInfo   
    float *frame_sig_Data, *frame_exp_Data; //stores frame pixel values for signal and exposure
    vector<string> header_info; //Stores Level-1 header information 
    int wtd_uc, wtd_bp, wtd_ff, wtd_pp, wtd_sd, wtd_cr, wtd_fi, wtd_fsc, wtd_dd, wtd_od, wtd_de, wtd_rfc, wtd_qemcp, wtd_centBias, wtd_centCorr, wtd_snr, wtd_wm; //modulewise flags to decide whether to write to disk or not
    char moduleoutdir_dd[NAMESIZE]; // Output directory path for detector distortion
    char moduleoutdir_bp[NAMESIZE]; // Output directory path for mask bad pixel
    char moduleoutdir_uc[NAMESIZE]; // Output directory path for unit conversion
    char moduleoutdir_ff[NAMESIZE]; // Output directory path for flatfield correction
    char moduleoutdir_qemcp[NAMESIZE]; // Output directory path for quantum efficiency
    char moduleoutdir_pp[NAMESIZE]; // Output directory path for pixel padding
    char moduleoutdir_sd[NAMESIZE]; // Output directory path for star detection
    char moduleoutdir_rav[NAMESIZE]; // Output directory path for registration and averaging
    char moduleoutdir_ravFlipped[NAMESIZE];//output directory for storing flipped image
    char moduleoutdir_radecimg[NAMESIZE]; // Output directory path for storing the RA-DEC converted image
    char moduleoutdir_cr[NAMESIZE]; // Output directory path for cosmic ray correction
    char moduleoutdir_fi[NAMESIZE]; // Output directory path for frame integration
    char moduleoutdir_od[NAMESIZE]; // Output directory path for optical distortion
    char moduleoutdir_snr[NAMESIZE]; // Output directory path for shift and rotate        
    char moduleoutdir_centroidCorr[NAMESIZE]; // Output directory path for centroid correction
    char moduleoutdir_centroidBias[NAMESIZE]; // Output directory path for centroid bias
    char moduleoutdir_expo[NAMESIZE];//output directory for storing exposure frames
    vector<string> name_track; //stores the  filename generated at output
    string snrimageFilename;//stores the shift and rotation  filename
    string regAvgFilename;
    string radecFilename,radecNoiceFilename,radecexpFilename;
    double star_time_ForDatase,end_time_ForDatase;
    vector<string> frameIntegrationfrmname;
     long no_of_records ;   
     int lastFileFlag;
bool flag_Roll_Angleinvalid;
double peakValueExp;
long ntotal,nzerocentroid;
double multFactor_ZeroCentroid;
//for dataIngest
    int parity_flag; //Parity flag (1-check parity for all 3 words, 2-check parity for x and y only) 
    int dropframe; //yes for dropping frame if CRC fails, no for dropping packet if CRC fails
    char imgfile[PIL_LINESIZE]; //name 
    char darkdir[FLEN_FILENAME]; //Stores Output Dark directory name
    char dataIngest_out_dir[FLEN_FILENAME]; //Stores DataIngest output directory path
    float darkFrameend_data[darkSize*darkSize]; //Stores pixels of dark frame for end frame
    float darkFramestart_data[darkSize*darkSize]; //Stores pixels of dark frame for start frame
    int gti_flag;
    int all_or_custom, valid_bit;
    char in_Gtifile[FLEN_FILENAME];
    int win_xsize, win_ysize; //Width and Height of the window (e.g (100, 100), (512, 512))       
    char eventfile[PIL_LINESIZE]; //name of input event file from DataIngest
    int xsize, ysize; //Width and Height of the frame
    vector<string> scienceDataSuccess;
    //vector<string> Total_ScienceDataFiles;
    //for unit  conversion
    char dstartpath[PIL_LINESIZE]; // Stores path for dark frame at begining
    char dendpath[PIL_LINESIZE]; // Stores path for dark frame at end
    double t_darkframeend; //Time for dark frame at begining
    double t_darkframestart; //Time for dark frame at end
    double t_curr; //Time of current frame

    //for Mask Bad pix
    char badpixfile[PIL_LINESIZE]; //CalDB bad pixel Filename
    char ExpTemplateFile[PIL_LINESIZE];
    float *badpixdata; //stores caldb data for bad pixel correction
    float thr_multiph;
    int  flag_Theta;
    //for pix padding 
    vector<long> ListOf_OutsideEvnts;

    //for sub division
    int subDivision_size;
    char rad_search[PIL_LINESIZE];
    int fi_flag;
    //for cosmic Ray corrections              
    float second_mult_Factor_CR; //for cosmic ray correction    
    float first_mult_factor_CR; //Threshold for cosmic ray correction

    //for flat field correction 
    char flatfieldfile[PIL_LINESIZE]; //CalDB flat field Filename
    float *flatfielddata; //stores caldb data for flat field correction

    //for QE and MCP correction
    char filter[FLEN_FILENAME];
    char qe_factor_file[PIL_LINESIZE]; //CalDB FilePath for QE_MCP correction
    float *cal_temp, *cal_f0, *cal_f1, *cal_f2, *cal_f3, *cal_f4, *cal_f5, *cal_f6, *cal_f7;
    double temperature; //Temperature from LBT file for current frame (in QE & MCP correction)
    long nCalDBTempValues; //Number of Temperature values (currently 20 in CalDB)
    long nrows_lbt;//Number of rows in LBT file
    double *time_lbt; //Time in LBT for each row
    float *insideTemp;//Inside temperature
    float *outsideTemp; //Outside temperature
    float *qe_mg_factor;//QE & MCP factor 
    char lbtfile[FLEN_FILENAME];//path of LBT file
       //for centroid Correction
    char centroidEAfile[PIL_LINESIZE];
    int centCorrflag,centBiasflag,DetectDistflag,OpticDistflag,qemcpflag;
    ; //CalDB Effective Area Energy Filename
    int cent_corr_win_size; //Centroid correction  window sizes

    //for centroid bias
    char centroidbiasfile[PIL_LINESIZE]; //CalDB Centroid bias Filename
    int biasRows;
    double *fraction_bias, *x_corr, *y_corr;
    int cent_bias_win_size; //Centroid bias  window sizes  

    //for detector distortion correction
    char detector_distortion_corr_file[PIL_LINESIZE]; //CalDB detector distortion Filename
    float *X_detect_distortion, *Y_detector_distortion; //stores detector distortion correction values from CalDB


    //for optical distortion correction
    char optical_distortion_corr_file[PIL_LINESIZE]; //CalDB optical distortion Filename
    float *X_optical_distortion, *Y_optical_distortion; //stores optical distortion correction values from CalDB

    //for frame integration
    int nFrameIntegrate_fi; //Number of frames to be integrated
    int nFrameDiscard_fi; //Number of frames to be discarded in frame integration
   int att_flag_val;

    //for shift and rotate
    double *time_drifts, *roll_ras, *pitch_ras, *yaw_ras;

    double *delta_x, *delta_y, *delta_theta;
    vector<long> track_rown;

    // for  registration and averaging
    int centroidlimit;
    int shift_N_Rotate_algo;
    int star_detect_algo_flag, centroid_Winsize, refine_Winsize,refine_Winsize_track,centroid_Winsize_track;
    float sd_mul_factor, sd_multi_factor_default;
    vector<Star> star_track;
    float diff_Dist;
    int minimum_No_of_Stars;
    vector<int> Fx, Fy, Rx, Ry;
    vector<float> Fval, Rval, Cx, Cy, Ci;
    float mult_fact;
    int flag_thetaComp;
    float roll_angle_applied;
    string filenamesnr;
    //for full frame Astorometry
    char catalogpath[PIL_LINESIZE];
    char att_timecol[NAMESIZE];
    char att_qcol[NAMESIZE];
    char databasename[NAMESIZE];
    char outtarpath[NAMESIZE];
    int IMG_DIM_FI;
    int FINALFRAMESIZE_REGAVG;
     int search_algo_ctlg;     
     char len_a[PIL_LINESIZE] ,len_b[PIL_LINESIZE];
     float Total_exp_time ;
     char dateOfObs[FLEN_FILENAME],timeobs[FLEN_FILENAME];\
     double RAPNT_FRMATT,DECPNT_FRMATT;
public:
  
    uvtLevel2PC();
    ~uvtLevel2PC();
    int uvtLevel2PCprocess();

    int readPILParameters(int argc, char** argv);
    vector<string> getFailedScienceDataFilename(){
        return this->scienceDataSuccess;
    }
//    vector<string> getTotalScienceDataFilename(){
//        return this->Total_ScienceDataFiles;
//    } 
   int copyAllheaderKeys(char *infile);
    int readImage(char * caldb_file, int hduno, float *caldb_data);
    int readcaldbtable(char * caldb_file, int hduno, vector<float> &caldb_data);
    int readcentroidbiasFile();
    int setDirectoryStructure(char *Dir, const char *subdir);
    int writeOutputImageToDisk(char *id, char *outDir, char *dir, char *subscript, float *Array, char *namepre, double ftime, unsigned short fno, int sizex, int sizey);
    int writeOutputTblToDisk(char *id, char *outDir, char *dir, char *subscript, char *namepre, double ftime, unsigned short fno, char **type1, char**tform1, int tfields1, char *extname, vector<float> &X, vector<float> &Y, vector<float> &val);

    double readDarkFrame(char * path, float*Array);
   // int darkFrameComputation(float *Array);
  //  int darkFrameSubtraction(float *Array, float *frame_data);
    int takeDarkinfo();



    int performUnitConversion(float *frmdata, float *expdata, double intgrntime, int sizex, int sizey);
   // int performCosmicRayCorr(float *frmsigdata, float *frmexpdata, int sizex, int sizey, float threshold_cr);
    int performFlatFieldCorr(float* frmsigdata, float* frmflatfield, float *Xfract, float *Yfract, int size);
    int performSubDivision(float *frmsigdata, int sizex, int sizey, float * subdivideddata, int size_subdivx, int size_subdivy);


    int performDistortionCorr(vector<float> &X, vector<float> &Y, float *Xdistr, float *Ydistr, int sizex, long caldbsize);
    int performMaskBadPix(unsigned short* int_x, unsigned short* int_y, float* badpixArr, unsigned char * Max_Min, unsigned char* badflag, unsigned char* multflag, int nrows, int x_size, int y_size, float thr_multphn);
    int performPixPadding(unsigned short *X_int, unsigned short *Y_int, int nrows);
    int performCentroidCorr(double *t, unsigned short *xint, unsigned short *yint, float *xf, float *yf, unsigned char *mc, float *newXfract, float *newYfract, float *darkbeginData, float *darkenddata, double integrtn_time, long EA, int nrows);

    int performCentroidBias(long nrows, float *xFrac, float *yFrac, float *new_xFrac, float *new_yFrac);
    int performDistCorrection(long nrows, float *xFrac, float *yFrac, float *x_Distortion, float *y_Distortion, int caldbdim);
//    int performFrameIntegration(long nrows, unsigned short*frame_no, int xsize, int Ndiscard, int Nacc, float *xFrac, float *yFrac, unsigned short *mult_phn, unsigned short *bad_Flag, double *t, float *ENP, float *one_dim_img, float *one_dim_exp, vector<FrameIntegration_Arr_strct> &vect);
//    int  performFrameIntegration (long nrows , long start_row,unsigned short*frame_no , int xsize , int Ndiscard , int Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr_strct> &vect,
//        float *outputSigArr,float *outputExpArray,float &outputFrmtime,long &outputFrmNo,long &last_index_for_frame);
    int performFrameIntegration (long nrows , long start_row,unsigned short*frame_no , int xsize , int Ndiscard , int Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr_strct> &vect
        ,float *noice_map_Sig_Array,float *outputSigArr,float *outputExpArray,float &outputFrmtime,long &outputFrmNo,long &last_index_for_frame);
    int performUnitConversion(float *frmdata, double intgrntime, int size);
    //int performShiftNRot (string sciencedataFile,long nrows , double *t , float *xf , float *yf , float *xi_fi , float *yi_fi);
    //int performShiftNRot (string sciencedataFile,long nrows , double *t , float *xf , float *yf , float *xi_fi , float *yi_fi,vector<float> &Shifts_x,vector<float> &Shifts_y,
      //  vector<float> &Shifts_theta);
    int performShiftNRot (string sciencedataFile,long nrows ,unsigned short *frmno, double *t , float *xf , float *yf , float *xi_fi , float *yi_fi,vector<float> &Shifts_x,vector<float> &Shifts_y,
        vector<float> &Shifts_theta,vector<unsigned short> &frm_no_vect,unsigned short *multphotonFlagIn,vector<unsigned short> &multphotonFlagOut);
//    int matchStars(int numrowsFirstfile, int numrowsSecfile, float divFact,
//            float *xlocFirst, float *ylocFirst, float *xlocSec, float *ylocSec, vector<float> &matchPixelXone,
//            vector<float> &matchPixelYone, vector<float> &matchPixelXtwo, vector<float> &matchPixelYtwo,
//            vector<float> &matchPixelDiffX, vector<float> &matchPixelDiffY);
    int matchStars (int numrowsFirstfile , int numrowsSecfile , float divFact ,
        float *xlocFirst , float *ylocFirst , float *xlocSec , float *ylocSec ,float * int1,float *int2, vector<float> &matchPixelXone ,
        vector<float> &matchPixelYone , vector<float> &matchPixelXtwo , vector<float> &matchPixelYtwo ,vector<float> &int_one,vector<float> &int_two,
        vector<float> &matchPixelDiffX , vector<float> &matchPixelDiffY);
   // int findShiftsNtheta(int totalelements, vector<float> &Xone, vector<float> &Yone, vector<float> &Xtwo, vector<float> &Ytwo,
     //       vector<float> &DiffOfX, vector<float> &DiffOfY, double &Xdx, double &Ydy, double &Theta);
    int findShiftsNtheta (int totalelements , vector<float> &Xone , vector<float> &Yone ,vector<float> &int1, vector<float> &Xtwo , vector<float> &Ytwo ,
       vector<float> &int2, vector<float> &DiffOfX , vector<float> &DiffOfY , bool flag_theta_computation,double &Xdx , double &Ydy , double &Theta);
    
//    
//    int findShiftsNtheta (int totalelements , vector<float> &Xone , vector<float> &Yone , vector<float> &Xtwo , vector<float> &Ytwo ,
//        vector<float> &DiffOfX , vector<float> &DiffOfY , bool flag_theta_computation,double &Xdx , double &Ydy , double &Theta);
    int  removeRecords(vector<float>  &Xone , vector<float> &Yone ,vector<float> &Xtwo , vector<float> &Ytwo,vector<float> &DiffOfX , vector<float> &DiffOfY ,float *ints1,float*ints2,
   vector<float> &newXone,vector<float> &newYone,vector<float> &newXtwo,vector<float> &newYtwo,vector<float> &newDiffX,vector<float> &newDiffY,vector<float> &new_one_ints,vector<float> &new_two_ints);

    int performQEMCPcorrection(float*sigdata, int size, double fact);
    int getHistory (vector<string> &vhistory);


    int findStar_algo1(float *inputArray);
    void doCentroiding(vector<int> &X, vector<int> &Y, int centroidwindow, float *arr, int h, int w);
    int ApplySubSampling(float* inputarray, int in_xsize, int in_ysize, float* outputarray, int out_xsize, int out_ysize);

//for driver module
     vector<string> sciencedata_file,status_info;
     bool paramfile_Varask_iml2;
     int tar_extracted_flag_PC;
     string orbnum;
     char * setRASfile(char * namefile )
    {
       sprintf( this->rasfile,"%s",namefile);
    }
  vector<string>  getSciencedatafile_info(){
            
            return sciencedata_file;
        };
     vector<string> getStatus_info(){
         return status_info;
     };
    int uvt_pc_l2_process();
    
      void setLevel1Dir(string  infile)
     {
       this->level1indir=  infile;
     }
        void setCaldbDir(string  calfile)
     {
       this->caldbindir=  calfile;
     }
     
     void setLevel2Dir(string  outfile)
     {
       this->level2outdir=  outfile;
     }


};

#endif	/* UVTLEVEL2PC_H */
