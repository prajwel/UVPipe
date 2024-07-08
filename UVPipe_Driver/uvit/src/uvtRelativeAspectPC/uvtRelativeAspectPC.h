///* 
// * File:   uvtRelativeAspectPC.h
// * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh

/* 
 * File:   uvtRelativeAspectPC.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on April 9, 2015, 2:40 PM
 */

//Header program for computing relative aspect series in photon-count mode
#ifndef UVTRELATIVEASPECTPC_H
#define	UVTRELATIVEASPECTPC_H
#include<string>
#include<iostream>
#include<stdlib.h>
#include<uvtUtils.h> //Common utilities
#include<caldb_Handler.h> //Handling CalDB filepath
#include<DataInfo.h>  //Interface for subsequent process, carries forward keywords
#include<macro_def.h> //Macros


//Structure for Frame integration process, stores Signal and Exposure image pixels, Frame number and Frame Time 
typedef struct FrameIntegration_Arr_strct {
   vector<float> pixels_sig_fi;//[IMG_DIM_FI * IMG_DIM_FI]; //Stores Signal image pixel values
   // float pixels_exp_fi[IMG_DIM_FI * IMG_DIM_FI]; //Stores Exposure image pixel values
    int frame_number_fi; //frame number
    double frame_time_fi; //frame time
};

//Class for uvtRelativeAspectPC
class uvtRelativeAspectPC 
{
    
    
    string modulename;
    int UTC_flag;
int crc_flag;
    float err_per;  //May be removed after checking
    int centroidlimit; //May be removed after checking
    int no_ofFramesToAcc; //May be removed after checking
     char caldbDir[PIL_LINESIZE]; //May be removed after checking
   float rms_mul_factor_default; //May be removed after checking
      float *badExparray;//May be removed after checking
    char mode[5];       //for PIL (query-learn(ql), auto(a), hidden(h))
    int clobber, history; //for overwrite and history respectively in FITS file from parameter file 
    int zipFlag;
    string level1indir; //Level-1 directory path from parameter file uvtRelativeAspectPC.par
    string level2outdir; //Output directory path from parameter file uvtRelativeAspectPC.par
    string caldbindir; //CalDB directory (CALDB/uvit) path from parameter file uvtRelativeAspectPC.par
    string channel; //FUV / NUV / VIS
    
    char stardir[10]; //Star directory path to store framewise star (first cut and refined peaks)
    char centroidDir[10]; //Centroid directory path to store framewise star centroids
    int star_detect_algo_flag; //Flag for SAC/PI or Joe's algorithm in star detection
    int centroid_Winsize; //Window size for Star centroid computation
    int footprint_win_star_detect; //Footprint window size in star detection for Joe's algorithm
    float primary_threshold_Val; //Threshold in star detection for Joe's algorithm
    float secondary_threshold_Val; //Threshold in star detection for Joe's algorithm
    
   int refine_Winsize; //Window size for refined star peak determination in SAC/PI algorithm
   float rms_mul_factor; //RMS factor for star detection in SAC/PI algorithm
   float sd_multi_factor_default;  //Standard deviation factor for star detection in SAC/PI algorithm
    int subDivision_size; //Frame size after subdivision (e.g 9600)

    int frames_toDiscard; //Number of frames to be discarded in reference frame calculation
    int nFrameToAverage; //Number of frames to be averaged in reference frame calculation
    int FreqDomainFilterFlag; //Frequency domain filter flag
    int fitting_flag; //Polynomial fitting flag for roll, pitch and yaw
    int orderroll, orderpitch, orderyaw; // Polynomial fitting order for roll, pitch and yaw
    
    int type_Filtering; //1-Polynomial fitting & 0-Spatial domain based on cutt-off frequency
    double cr_threshold; //Threshold for cosmic ray correction
    double freqvalue; //Cut-off frequency
    double poly_fit_interval; //Piecewise polynomial fitting interval (in seconds) for drift series
   
    DataInfo datainfo; //object for DataInfo
    
    caldb_Handler caldb_handler; //Object for caldb_handler with methods to access different calibfile.fits file

    char badpixfile[PIL_LINESIZE]; //CalDB bad pixel Filename
    float *badpixdata; //stores caldb data for bad pixel correction
    char flatfieldfile[PIL_LINESIZE];//CalDB flat field Filename
    float *flatfielddata; //stores caldb data for flat field correction
    
    char centroidEAfile[PIL_LINESIZE]; //CalDB Effective Area Energy Filename
    char centroidbiasfile[PIL_LINESIZE]; //CalDB Centroid bias Filename
    char detector_distortion_corr_file[PIL_LINESIZE];//CalDB detector distortion Filename
     float *X_detect_distortion, *Y_detector_distortion; //stores detector distortion correction values from CalDB
    char optical_distortion_corr_file[PIL_LINESIZE];//CalDB optical distortion Filename
     float  *X_optical_distortion, *Y_optical_distortion; //stores optical distortion correction values from CalDB
    char qe_factor_file[PIL_LINESIZE];//CalDB quantum efficiency Filename
    
    char moduleoutdir_dd[NAMESIZE]; // Output directory path for detector distortion
    char moduleoutdir_bp[NAMESIZE];// Output directory path for mask bad pixel
    char moduleoutdir_uc[NAMESIZE];// Output directory path for unit conversion
    char moduleoutdir_ff[NAMESIZE];// Output directory path for flatfield correction
    char moduleoutdir_qemcp[NAMESIZE];// Output directory path for quantum efficiency
    char moduleoutdir_pp[NAMESIZE];// Output directory path for pixel padding
    char moduleoutdir_sd[NAMESIZE];// Output directory path for star detection
    char moduleoutdir_cr[NAMESIZE];// Output directory path for cosmic ray correction
    char moduleoutdir_fi[NAMESIZE];// Output directory path for frame integration
    char moduleoutdir_sc[NAMESIZE];// Output directory path for star centroids
    char moduleoutdir_od[NAMESIZE];// Output directory path for optical distortion
    char moduleoutdir_de[NAMESIZE];// Output directory path for drift computation
    char moduleoutdir_rfc[NAMESIZE];// Output directory path for reference frame calculation
    char moduleoutdir_centroidCorr[NAMESIZE];// Output directory path for centroid correction
    char moduleoutdir_centroidBias[NAMESIZE];// Output directory path for centroid bias correction
    char dstartpath[PIL_LINESIZE]; // Stores path for dark frame at begining
    char dendpath[PIL_LINESIZE]; // Stores path for dark frame at end
    double t_darkframeend;	//Time for dark frame at begining
    double t_darkframestart;	//Time for dark frame at end
//Refactor this to pair_nbhd_distance
    float pair_nbhd_distance;  //Neighbourhood for star pairing/searching
    double t_curr;  //Time of current frame
    int xsize, ysize; //Width and Height of the frame
    char RelAspFilename[FLEN_FILENAME]; //Path of relative aspect series file
//Refactor this to search_win_size
    int win_size;	//Window size for searching stars in consecutive frames
    int win_xsize,win_ysize; //Width and Height of the window (e.g (100, 100), (512, 512))
//Refactor dark_subtraction_flag
    int dark_subtraction_flag; //Dark Subtraction flag
//Refactor to following
    //vector<int> refined_peaks_x_vect, refined_peaks_y_vect, firstcut_peaks_x_vect, firstcut_peaks_y_vect;
    vector<int> refined_peaks_x_vect, refined_peaks_y_vect, firstcut_peaks_x_vect, firstcut_peaks_y_vect; //Stores x, y location of refined and firstcut peaks
//Refactor to following
    //vector<float> centroid_x_vect, centroid_y_vect, temp_centroid_x_vect, temp_centroid_y_vect; //Stores x, y location of final and temporary centroid 
    vector<float> centroid_x_vect, centroid_y_vect, centroid_x_temp_vect, centroid_y_temo_vect;
//Refactor to following
    //vector<float> c_int, R_int, F_int, ci_tempi; //Stores intensity for centroid, refined peaks, firstcut peaks and temporary centroid
    vector<float> c_int, R_int, F_int, ci_tempi;

//Flag for shift and rotation computation algorithm (1-SAC, 2-SAC, 3-IUCAA)
   int shift_rotation_algo;

    vector<string> ref_frame_module_filename_track;
    int parity_flag; //Parity flag (1-check parity for all 3 words, 2-check parity for x and y only)
    vector<float> x_shift, y_shift, theta_shift; //Stores x-y shifts and rotations

//Refactor to following
    //vector<double> ref_frame_time_data; //Stores frame time for frames in reference frame calculation
    vector<double> ref_frame_time_data;

    float darkFrameend_data[darkSize*darkSize]; //Stores pixels of dark frame for end frame
    float darkFramestart_data[darkSize*darkSize]; //Stores pixels of dark frame for start frame

    //Stores x,y and intensity values for all frames used for reference frame calculation
    vector<float> x_Of_refFrame, y_Of_refFrame, int_Of_refFrame;
    //Stores average x,y and intensity values for all frames used for reference frame calculation
    vector<float> x_ref_cumm, y_ref_cumm, int_ref_cumm; 
    vector<float> roll_vect, pitch_vect, yaw_vect; 

//To be refined
    int gti_flag; 
    int all_or_custom, valid_bit;
    char in_Gtifile[FLEN_FILENAME];


    int cent_corr_win_size, cent_bias_win_size; 

    //for dataingest
    int dropframe; //yes for dropping frame if CRC fails, no for dropping packet if CRC fails
    vector<int> size_rows;
    char darkdir[FLEN_FILENAME];
    char Indir_dataIngest[FLEN_FILENAME];
    int att_flag_val;
    //for unit conversion
    int unitConversionFlag; //YES if unit conversion needs to run else NO

    //for flat field correction module
    int flatfieldFlag; //YES if flat field correction is to be done, No if not required

    //for pixel padding
    int padding_dim; //required dimension after padding

    //for Accumulate
    int nFrameIntegrate_fi; //number of frames to integrate
    int nFrameDiscard_fi; //Number of frames to be discarded in frame integration

    //for corrrection of temperature effects
    int qemcpFlag; //YES to apply
     float *temp,*f0,*f1,*f2,*f3,*f4,*f5,*f6,*f7;
      char qeFile[PIL_LINESIZE];
  char lbtfile[FLEN_FILENAME];
  double temperature;
    //for pixel subdivision      
    int subdivisionFlag; //YES to apply
    vector<string> header_info;
    bool lastFrame_flag;
    float *frame_Data, *frame_ExpData;
    float *frame_Data_Padded, *frame_ExpoData_padded, *frame_fc_data;
    float *frame_Data_subdivided, *frame_ExpData_subdivided;
    int wtd_uc, wtd_bp, wtd_ff, wtd_pp, wtd_sd, wtd_cr, wtd_fi, wtd_fsc, wtd_dd, wtd_od, wtd_de, wtd_rfc, wtd_qemcp, wtd_centBias, wtd_centCorr;
    float *darkCompute_array;
    int cmpr_frames;
    char eventfile[PIL_LINESIZE]; //name of input event file
    char imgfile[PIL_LINESIZE]; //name 
    float *darksubtr_array;
    int nCompareFrames;
    float Threshold_cr;
    int IMG_DIM_FI;
 int  flag_thetaComp;
 float  first_mult_factor_CR,second_mult_Factor_CR;
public:
    int getHistory (vector<string> &vhistory);
    int setDirectoryStructure(char *Dir, const char *subdir);
    int writeOutputImageToDisk(char *id, char *outDir, char *dir, char *subscript, float *Array, char *namepre, double ftime, unsigned short fno, int sizex, int sizey);
    double* DataFitting(double *t, double *X, int order, int nRecords);
    void doFft(double data[], long nn, int isign);
    int transformToUVITFrame(double *t, double *r, double *p, double *y);
    int transformToSpacecraftFrame(double *t, double *r, double *p, double *y);
    double readDarkFrame(char * path, float*Array);
    int darkFrameComputation(float *Array);
    int darkFrameSubtraction(float *Array, float *frame_data);
    int takeDarkinfo();
    uvtRelativeAspectPC();
    ~uvtRelativeAspectPC();
     long  nCalDBTempValues;
         long nrows_lbt ;
        double *time_lbt;
          float *insideTemp;
          float *outsideTemp;
        float *qe_mg_factor;
        //  int  readQEMCPFile();
//           int getTemp ();
    char filter[FLEN_FILENAME];
    int performQEMCPcorrection(float*sigdata, int size, double fact);
    int uvtRelativeAspectPCProcess();
    int readPILParameters(int argc, char** argv);
    int minimum_No_of_Stars;
    int match_stars_file_flag;
    float thr_multiph;
    int biasRows;
    float backgrd_fact;
    double *fraction_bias, *x_corr, *y_corr;
    int copyAllheaderKeys(char *infile);
    int readImage(char * caldb_file, int hduno, float *caldb_data);
    int readcaldbtable(char * caldb_file, int hduno, vector<float> &caldb_data);
    int performUnitConversion(float *frmdata, float *expdata, double intgrntime, int sizex, int sizey);
    //int performCorrection(float *frmsigdata, float *frmexpdata, float *badpixarry,int sizex,int sizey,double  intgrntime);
    int performCosmicRayCorr(float *frmsigdata, float *frmexpdata, int sizex, int sizey, float threshold_cr);
    int performFlatFieldCorr(float* frmsigdata, float* frmflatfield, float *Xfract, float *Yfract, int size);
    int performSubDivision(float *frmsigdata, int sizex, int sizey, float * subdivideddata, int size_subdivx, int size_subdivy);
    int writeOutputTblToDisk(char *id, char *outDir, char *dir, char *subscript, char *namepre, double ftime, unsigned short fno, char **type1, char**tform1, int tfields1, char *extname, vector<float> &X, vector<float> &Y, vector<float> &val);
    //void     getKeywordValue (char *file , int type,char *keyname , int hdunum , char *val);
    int performAccOFFrame(float *frmsigdata, float *sumdata, int sizex, int sizey, int numoffrmAcc);
    int performDistortionCorr(vector<float> &X, vector<float> &Y, float *Xdistr, float *Ydistr, int sizex, long caldbsize);
    int calculateRefFrame(vector<int> &sizevect, vector<float> &xref, vector<float> &yref, vector<float> &intref, int avgfact, vector<float> &xrefcumm, vector<float> &yrefcumm, vector<float> &intrefcumm);
    int performMaskBadPix(unsigned short* int_x, unsigned short* int_y, float* badpixArr, unsigned char * Max_Min, unsigned char* badflag, unsigned char* multflag, int nrows, int x_size, int y_size, float thr_multphn);
    int performPixPadding(unsigned short *X_int, unsigned short *Y_int, int nrows);
    int performCentroidCorr(double *t, unsigned short *xint, unsigned short *yint, float *xf, float *yf, unsigned char *mc, float *newXfract, float *newYfract, float *darkbeginData, float *darkenddata, double integrtn_time, long EA, int nrows);
    int readcentroidbiasFile();
    int performCentroidBias(long nrows, float *xFrac, float *yFrac, float *new_xFrac, float *new_yFrac);
    int performDistCorrection(long nrows, float *xFrac, float *yFrac, float *x_Distortion, float *y_Distortion, int caldbdim);
//  int  performFrameIntegration (long nrows , long start_row,unsigned short*frame_no , int xsize , int Ndiscard , int Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr_strct> &vect,
//        float *outputSigArr,float *outputExpArray,double &outputFrmtime,long &outputFrmNo,long &last_index_for_frame);
    int  performFrameIntegration (long nrows , long start_row,unsigned short*frame_no , int xsize , int Ndiscard , int Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr_strct> &vect,
        float *outputSigArr,float *outputExpArray,double &outputFrmtime,long &outputFrmNo,long &last_index_for_frame);
    int performUnitConversion(float *frmdata, double intgrntime, int size);

    //for driver module
    string orbnum;
      bool paramfile_Varask,tar_extracted_flag_PC;
       vector<string> sciencedata_file,status_info;
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
       const char * getRelAspFilename()
    {
        return this->RelAspFilename;
    }
       vector<string> getStatus_info(){
         return status_info;
     }
       vector<string>  getSciencedatafile_info(){
            
            return sciencedata_file;
        };
};

#endif	/* UVTRLEATIVEASPECTPC_H */






// * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
// *
// * Created on April 9, 2015, 2:40 PM
// */
//
//#ifndef UVTRELATIVEASPECTPC_H
//#define	UVTRELATIVEASPECTPC_H
//#include<string>
//#include<iostream>
//#include<stdlib.h>
//#include<uvtUtils.h> //Common utilities
//#include<caldb_Handler.h> //Handling CalDB filepath
//#include<DataInfo.h>  //Interface for subsequent process, carries forward keywords
//#include<macro_def.h> //Macros
//
//
////Structure for Frame integration process, stores Signal and Exposure image pixels, Frame number and Frame Time 
//typedef struct FrameIntegration_Arr {
//    float pixels_sig_fi[IMG_DIM_FI * IMG_DIM_FI]; //Stores Signal image pixel values
//    float pixels_exp_fi[IMG_DIM_FI * IMG_DIM_FI]; //Stores Exposure image pixel values
//    int frame_number_fi; //frame number
//    double frame_time_fi; //frame time
//};
//
////Class for uvtRelativeAspectPC
//class uvtRelativeAspectPC {
//    float err_per;  //May be removed after checking
//    int centroidlimit; //May be removed after checking
//    int no_ofFramesToAcc; //May be removed after checking
//     char caldbDir[PIL_LINESIZE]; //May be removed after checking
//      float *badExparray;//May be removed after checking
//    char mode[5];       //for PIL
//    int clobber, history; //for overwrite and history respectively in FITS file from parameter file 
//    
//    string level1indir; //Level-1 directory path from parameter file uvtRelativeAspectPC.par
//    string level2outdir; //Output directory path from parameter file uvtRelativeAspectPC.par
//    string caldbindir; //CalDB directory (CALDB/uvit) path from parameter file uvtRelativeAspectPC.par
//    string channel; //FUV / NUV / VIS
//    
//    char stardir[10]; //Star directory path to store framewise star (first cut and refined peaks)
//    char centroidDir[10]; //Centroid directory path to store framewise star centroids
//    int star_detect_algo_flag; //Flag for SAC/PI or Joe's algorithm in star detection
//    int centroid_Winsize; //Window size for Star centroid computation
//    int footprint_win_star_detect; //Footprint window size in star detection for Joe's algorithm
//    float primary_threshold_Val; //Threshold in star detection for Joe's algorithm
//    float secondary_threshold_Val; //Threshold in star detection for Joe's algorithm
//    
//   int refine_Winsize; //Window size for refined star peak determination in SAC/PI algorithm
//   float rms_mul_factor; //RMS factor for star detection in SAC/PI algorithm
//   float sd_multi_factor_default;  //Standard deviation factor for star detection in SAC/PI algorithm
//   float rms_mul_factor_default;
//    int subDivision_size; //Frame size after subdivision (e.g 9600)
//    
//
//
//    int frames_toDiscard; //Number of frames to be discarded in reference frame calculation
//    int nFrameToAverage; //Number of frames to be averaged in reference frame calculation
//    int FreqDomainFilterFlag; //Frequency domain filter flag
//    int fitting_flag; //Polynomial fitting for roll, pitch and yaw
//    int orderroll, orderpitch, orderyaw; // Polynomial fitting order for roll, pitch and yaw
//    
//    int type_Filtering; //1-Polynomial fitting & 0-Spatial domain based on cutt-off frequency
//    double cr_threshold; //Threshold for cosmic ray correction
//    double freqvalue; //Cut-off frequency
//    double poly_fit_interval; //Piecewise polynomial fitting interval (in seconds) for drift series
//   
//    DataInfo datainfo; //object for DataInfo
//    
//    caldb_Handler caldb_handler; //Object for caldb_handler with methods to access different calibfile.fits file
//
//    char badpixfile[PIL_LINESIZE]; //CalDB bad pixel Filename
//    float *badpixdata; //stores caldb data for bad pixel correction
//    char flatfieldfile[PIL_LINESIZE];//CalDB flat field Filename
//    float *flatfielddata; //stores caldb data for flat field correction
//    
//    char centroidEAfile[PIL_LINESIZE]; //CalDB Effective Area Energy Filename
//    char centroidbiasfile[PIL_LINESIZE]; //CalDB Centroid bias Filename
//    char detector_distortion_corr_file[PIL_LINESIZE];//CalDB detector distortion Filename
//     float *X_detect_distortion, *Y_detector_distortion; //stores detector distortion correction values from CalDB
//    char optical_distortion_corr_file[PIL_LINESIZE];//CalDB optical distortion Filename
//     float  *X_optical_distortion, *Y_optical_distortion; //stores optical distortion correction values from CalDB
//    char qe_factor_file[PIL_LINESIZE];//CalDB quantum efficiency Filename
//    
//    char moduleoutdir_dd[NAMESIZE]; // Output directory path for detector distortion
//    char moduleoutdir_bp[NAMESIZE];// Output directory path for mask bad pixel
//    char moduleoutdir_uc[NAMESIZE];// Output directory path for unit conversion
//    char moduleoutdir_ff[NAMESIZE];// Output directory path for flatfield correction
//    char moduleoutdir_qemcp[NAMESIZE];// Output directory path for quantum efficiency
//    char moduleoutdir_pp[NAMESIZE];// Output directory path for pixel padding
//    char moduleoutdir_sd[NAMESIZE];// Output directory path for star detection
//    char moduleoutdir_cr[NAMESIZE];// Output directory path for cosmic ray correction
//    char moduleoutdir_fi[NAMESIZE];// Output directory path for frame integration
//    char moduleoutdir_sc[NAMESIZE];// Output directory path for star centroids
//    char moduleoutdir_od[NAMESIZE];// Output directory path for optical distortion
//    char moduleoutdir_de[NAMESIZE];// Output directory path for drift computation
//    char moduleoutdir_rfc[NAMESIZE];// Output directory path for reference frame calculation
//    char moduleoutdir_centroidCorr[NAMESIZE];// Output directory path for centroid correction
//    char moduleoutdir_centroidBias[NAMESIZE];// Output directory path for centroid bias correction
//    char dstartpath[PIL_LINESIZE];
//    char dendpath[PIL_LINESIZE];
//    double t_darkframeend;
//    double t_darkframestart;
//    float diff_Dist;
//    double t_curr;
//    int xsize, ysize;
//    int centroid_row;
//    char RelAspFilename[FLEN_FILENAME];
//    int Win_NB;
//    int win_xsize,win_ysize;
//    long caldb_dim;
//    int darkFrame_Flag;
//    vector<int> Rx_vect, Ry_vect, Fx_vect, Fy_vect;
//    vector<float> cx_vect, cy_vect, c_int, R_int, F_int, cx_tempx, cy_tempy, ci_tempi;
//    int option_LeastSquare;
//    vector<string> name_track;
//    int parity_flag;
//    vector<float> x_shift, y_shift, theta_shift;
//    vector<double> frm_time_data;
//    int DATA_ARRAYSIZE;
//    float darkFrameend_data[darkSize*darkSize];
//    float darkFramestart_data[darkSize*darkSize];
//    vector<float> x_ref_cumm, y_ref_cumm, int_ref_cumm;
//    int gti_flag, windowsize, win_centBias;
//    int all_or_custom, valid_bit;
//    char in_Gtifile[FLEN_FILENAME];
//    vector<float> roll_vect, pitch_vect, yaw_vect, x_Of_refFrame, y_Of_refFrame, int_Of_refFrame;
//    //        float darkFrameend_data[darkSize*darkSize];
//    //        float darkFramestart_data[darkSize*darkSize];
//    //for dataingest
//    int dropframe; //yes for dropping frame if CRC fails, no for dropping packet if CRC fails
//    vector<int> size_rows;
//    char darkdir[FLEN_FILENAME];
//    char Indir_dataIngest[FLEN_FILENAME];
//    //for unit conversion
//    int unitConversionFlag; //YES if unit conversion needs to run else NO
//
//    //for flat field correction module
//    int flatfieldFlag; //YES if flat field correction is to be done, No if not required
//
//    //for pixel padding
//    int padding_dim; //required dimension after padding
//
//    //for Accumulate
//    int nFrameIntegrate_fi; //number of frames to integrate
//    int nFrameDiscard_fi; //Number of frames to be discarded in frame integration
//    //for corrrection of temperature effects
//    int qemcpFlag; //YES to apply
//     float *temp,*f0,*f1,*f2,*f3,*f4,*f5,*f6,*f7;
//      char qeFile[PIL_LINESIZE];
//  char lbtfile[FLEN_FILENAME];
//  double temperature;
//    //for pixel subdivision      
//    int subdivisionFlag; //YES to apply
//    vector<string> header_info;
//
//    float *frame_Data, *frame_ExpData;
//    float *frame_Data_Padded, *frame_ExpoData_padded, *frame_fc_data;
//    float *frame_Data_subdivided, *frame_ExpData_subdivided;
//    int wtd_uc, wtd_bp, wtd_ff, wtd_pp, wtd_sd, wtd_cr, wtd_fi, wtd_fsc, wtd_dd, wtd_od, wtd_de, wtd_rfc, wtd_qemcp, wtd_centBias, wtd_centCorr;
//    float *darkCompute_array;
//    int cmpr_frames;
//    char eventfile[PIL_LINESIZE]; //name of input event file
//    char imgfile[PIL_LINESIZE]; //name 
//    float *darksubtr_array;
//    int nCompareFrames;
//    float Threshold_cr;
//public:
//
//    int setDirectoryStructure(char *Dir, const char *subdir);
//    int writeOutputImageToDisk(char *id, char *outDir, char *dir, char *subscript, float *Array, char *namepre, double ftime, unsigned short fno, int sizex, int sizey);
//    double* DataFitting(double *t, double *X, int order, int nRecords);
//    void doFft(double data[], long nn, int isign);
//    int transformToUVITFrame(double *t, double *r, double *p, double *y);
//    int transformToSpacecraftFrame(double *t, double *r, double *p, double *y);
//    double readDarkFrame(char * path, float*Array);
//    int darkFrameComputation(float *Array);
//    int darkFrameSubtraction(float *Array, float *frame_data);
//    int takeDarkinfo();
//    uvtRelativeAspectPC();
//    ~uvtRelativeAspectPC();
//     int  nCalDBTempValues;
//         long nrows_lbt ;
//        double *time_lbt;
//          float *insideTemp;
//          float *outsideTemp;
//        float *qe_mg_factor;
//          int  readQEMCPFile();
//           int getTemp ();
//    char filter[FLEN_FILENAME];
//    int performQEMCPcorrection(float*sigdata, int size, double fact);
//    int uvtRelativeAspectPCProcess();
//    int read(int argc, char** argv);
//    int minimum_No_of_Stars;
//    int match_stars_file_flag;
//    float thr_multiph;
//    int biasRows;
//    double *fraction_bias, *x_corr, *y_corr;
//    int copyAllheaderKeys(char *infile);
//    int readImage(char * caldb_file, int hduno, float *caldb_data);
//    int readcaldbtable(char * caldb_file, int hduno, vector<float> &caldb_data);
//    int performUnitConversion(float *frmdata, float *expdata, double intgrntime, int sizex, int sizey);
//    //int performCorrection(float *frmsigdata, float *frmexpdata, float *badpixarry,int sizex,int sizey,double  intgrntime);
//    int performCosmicRayCorr(float *frmsigdata, float *frmexpdata, int sizex, int sizey, float threshold_cr);
//    int performFlatFieldCorr(float* frmsigdata, float* frmflatfield, float *Xfract, float *Yfract, int size);
//    int performSubDivision(float *frmsigdata, int sizex, int sizey, float * subdivideddata, int size_subdivx, int size_subdivy);
//    int writeOutputTblToDisk(char *id, char *outDir, char *dir, char *subscript, char *namepre, double ftime, unsigned short fno, char **type1, char**tform1, int tfields1, char *extname, vector<float> &X, vector<float> &Y, vector<float> &val);
//    //void     getKeywordValue (char *file , int type,char *keyname , int hdunum , char *val);
//    int performAccOFFrame(float *frmsigdata, float *sumdata, int sizex, int sizey, int numoffrmAcc);
//    int performDistortionCorr(vector<float> &X, vector<float> &Y, float *Xdistr, float *Ydistr, int sizex, long caldbsize);
//    int calculateRefFrame(vector<int> &sizevect, vector<float> &xref, vector<float> &yref, vector<float> &intref, int avgfact, vector<float> &xrefcumm, vector<float> &yrefcumm, vector<float> &intrefcumm);
//    int performMaskBadPix(unsigned short* int_x, unsigned short* int_y, float* badpixArr, unsigned char * Max_Min, unsigned char* badflag, unsigned char* multflag, int nrows, int x_size, int y_size, float thr_multphn);
//    int performPixPadding(unsigned short *X_int, unsigned short *Y_int, int nrows);
//    int performCentroidCorr(double *t, unsigned short *xint, unsigned short *yint, float *xf, float *yf, unsigned char *mc, float *newXfract, float *newYfract, float *darkbeginData, float *darkenddata, double integrtn_time, long EA, int nrows);
//    int readcentroidbiasFile();
//    int performCentroidBias(long nrows, float *xFrac, float *yFrac, float *new_xFrac, float *new_yFrac);
//    int performDistCorrection(long nrows, float *xFrac, float *yFrac, float *x_Distortion, float *y_Distortion, int caldbdim);
//    int performFrameIntegration(long nrows, unsigned short*frame_no, int xsize, int Ndiscard, int Nacc, float *xFrac, float *yFrac, unsigned short *mult_phn, unsigned short *bad_Flag, double *t, float *ENP, float *one_dim_img, float *one_dim_exp, vector<FrameIntegration_Arr> &vect);
//    int performUnitConversion(float *frmdata, double intgrntime, int size);
//
//};
//
//#endif	/* UVTPCRACOMMANARR_H */
