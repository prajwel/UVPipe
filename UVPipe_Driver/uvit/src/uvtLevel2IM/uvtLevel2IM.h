/* 
 * File:   uvtImL2_commonArray.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on December 9, 2013, 4:11 PM
 */

#ifndef  UVTLEVEL2IM_H
#define UVTLEVEL2IM_H
#include<string>
#include<iostream>
#include<stdlib.h>
#include<uvtUtils.h>
#include<caldb_Handler.h>
#include<DataInfo.h>
#include<uvtFullFrameAst.h>
#include<Database.h>
#include<Attitude.h>
#define darkSize 512

typedef struct Star_l2
{
    float intensity;
    int x, y;
};

class uvtLevel2IM 
{
    float *InputArray;
    int UTC_flag,crc_flag;
    string modulename;
    string level1indir; //Level-1 directory path from parameter file uvtLevel2IM.par
    string caldbindir; //CalDB directory (CALDB/uvit) path from parameter file uvtLevel2IM.par
    string level2outdir; //Output directory path from parameter file uvtLevel2IM.par
    string channel; //FUV / NUV / VIS
    caldb_Handler caldb_handler; //object for  caldb_handler
    DataInfo datainfo; //object for DataInfo
    vector<string> header_info; //Stores Level-1 header information  
    int wtd_uc, wtd_bp, wtd_ff, wtd_pp, wtd_sd, wtd_cr, wtd_ac, wtd_fsc, wtd_dd, wtd_od, wtd_de, wtd_rfc, wtd_qemcp, wtd_snr, wtd_wm, wtd_ravg; //flags for writing outputs to disk
    char moduleoutdir_dd[NAMESIZE]; // Output directory path for detector distortion
    char moduleoutdir_bp[NAMESIZE]; //Output directory path for masking bad pixel
    char moduleoutdir_uc[NAMESIZE]; //Output directory path for unit conversion
    char moduleoutdir_ff[NAMESIZE]; //output directory path for flatfield corr
    char moduleoutdir_qemcp[NAMESIZE]; //Output directoru path for qe and MCP
    char moduleoutdir_pp[NAMESIZE]; //Output directory path for pix padding 
    char moduleoutdir_sd[NAMESIZE]; //Output directory path for Sub Division
    char moduleoutdir_cr[NAMESIZE]; //Output directory path for Cosmic Ray correction
    char moduleoutdir_ac[NAMESIZE]; //Output directory path for accumulation
    char moduleoutdir_sc[NAMESIZE]; //Output directory path for star detection
    char moduleoutdir_od[NAMESIZE]; //Output directory path for Optical Distortion
    char moduleoutdir_de[NAMESIZE]; //Output directory path for drift series
    char moduleoutdir_rfc[NAMESIZE]; //Output directory path for reference frame calculation
    char moduleoutdir_snr[NAMESIZE]; //output directory path for shift and rotate
    char moduleoutdir_wm[NAMESIZE]; //output directory path for weighted mean
    char moduleoutdir_ravg[NAMESIZE]; //output directory path for registration and averaging
    int clobber, history;

    char mode[5];
    float *frame_Data, *frame_ExpData;
    float *frame_Data_Padded, *frame_ExpoData_padded, *frame_fc_data;
    float *frame_Data_subdivided, *frame_ExpData_subdivided;
    int xsize, ysize;
    vector<string> ref_frame_module_filename_track; //Keeps track of filename for output of reference frame calculation(no need for LEVEL2 IM)
    vector<int> track_invalidPix;


    //for dataingest
    int gti_flag;
    int all_or_custom, valid_bit;
    char in_Gtifile[FLEN_FILENAME];
    char lbtfile[FLEN_FILENAME];
    int dropframe; //yes for dropping frame if CRC fails, no for dropping packet if CRC fails
    char dataIngest_out_dir[FLEN_FILENAME];
    int win_xsize, win_ysize; //Width and Height of the window (e.g (100, 100), (512, 512))
    int parity_flag;
    char darkdir[FLEN_FILENAME];
    int att_flag_val;

    //for unitconversion                     
    int darkframe_flag;
    float *darkCompute_array;
    int unitConversionFlag; //YES if unit conversion needs to run else NO
    float darkFrameend_data[darkSize*darkSize]; //Stores pixels of dark frame for end frame
    float darkFramestart_data[darkSize*darkSize]; //Stores pixels of dark frame for start frame
    char dstartpath[PIL_LINESIZE]; //stores  path for dark start path.
    char dendpath[PIL_LINESIZE]; //stores path for dark end path
    double t_curr;//stores time of current frame in processing
    double t_darkframeend; //stores dark start path
    double t_darkframestart; //stores dark end path

    //for bad pixel 
    char badpixfile[PIL_LINESIZE]; //CalDB bad pixel Filename
    float *badpixdata; //stores caldb data for bad pixel correction

    //for cosmic ray correction
      double cr_threshold;
    
    //for flat field correction 
    int flatfieldFlag; //YES if flat field correction is to be done, No if not required
    char flatfieldfile[PIL_LINESIZE]; //CalDB flat field Filename
    float *flatfielddata; //stores caldb data for flat field correction
char* rad_search;
    
    //for QE and MCP corrections  
    char filter[FLEN_FILENAME];
    int qe_mcpFlag; //YES to apply
    double temperature;
    long nCalDBTempValues; //Number of Temperature values (currently 20 in CalDB)
    long nrows_lbt; //Number of rows in LBT file
    double *time_lbt; //Time in LBT for each row
    float *insideTemp; //Inside temperature
    float *outsideTemp; //Outside temperature
    float *qe_mg_factor; //QE & MCP factor 
    char qe_factor_file[PIL_LINESIZE];
    float *cal_temp, *cal_f0, *cal_f1, *cal_f2, *cal_f3, *cal_f4, *cal_f5, *cal_f6, *cal_f7;

    //for pixel padding 
    int padding_dim; //required dimension after padding

    //for sub division         
    int subdivisionFlag; //YES to apply
    int subDivision_size;

    //for detetctor distortion correction
    char detector_distortion_corr_file[PIL_LINESIZE];
    float *X_detect_distortion, *Y_detect_distortion;

    //for optical distortion correction

    char optical_distortion_corr_file[PIL_LINESIZE];
    float *X_optical_distortion, *Y_optical_distortion;

    //for shift and rotate
    float mult_fact;
    char rasfile[PIL_LINESIZE];

    //for finding Weighted mean
    int no_ofWeigh;


    //for registration and averaging
    int centroidlimit;
    int star_detect_algo_flag, centroid_Winsize, refine_Winsize;
   // int star_detect_algo_flag;
    vector<Star_l2> star_track;
    float diff_Dist;
    int minimum_No_of_Stars;
    vector<int> Fx, Fy, Rx, Ry;
    float sd_mul_factor, sd_multi_factor_default;
    vector<float> Fval, Rval, Cx, Cy, Ci;
    int flag_thetaComp;
    
    //for Full frame Astrometry
     int search_algo_ctlg;     
    char len_a[PIL_LINESIZE] ,len_b[PIL_LINESIZE];
    char databasename[NAMESIZE];//path of catalog
    char catalogpath[PIL_LINESIZE];
    char att_timecol[NAMESIZE];
    char att_qcol[NAMESIZE];
    int FINALFRAMESIZE_REGAVG;

public:
    uvtLevel2IM(); //Constructor
    ~uvtLevel2IM(); //Destructor

    int uvtLevel2IMProcess();//performs LEVEL2 operations

    int readPILParameters(int argc, char** argv);//Reads PIL parameters
    double readDarkFrame(char * path, float*Array);//Reads output dark frames from DataIngest
    int darkFrameComputation(float *Array);
    int takeDarkinfo();

    int copyAllheaderKeys(char *infile);
    int setDirectoryStructure(char *Dir, const char *subdir);//Makes directory structure
    int writeOutputImageToDisk(char *id, char *outDir, char *dir, char *subscript, float *Array, char *namepre, double ftime, unsigned short fno, int sizex, int sizey);
    int writeOutputTblToDisk(char *id, char *outDir, char *dir, char *subscript, char *namepre, double ftime, unsigned short fno, char **type1, char**tform1, int tfields1, char *extname, vector<float> &X, vector<float> &Y, vector<float> &val);
    //  int  transformToUVITFrame(double *t,double *r,double *p,double *y);
    // int  transformToSpacecraftFrame(double *t,double *r,double *p,double *y);

    int readRASfile(char *ras_file, float *time_ras, float *roll_ras, float *pitch_ras, float *yaw_ras);//read relative aspect file  generated  after relative aspect chain.

int getHistory (vector<string> &vhistory);
    int findStar_algo1(float *inputArray);//fnding stars and centroids of image frame for registration and averaging
    void doCentroiding(vector<int> &X, vector<int> &Y, int centroidwindow, float *arr, int h, int w);//perform centroiding
    int matchStars(int numrowsFirstfile, int numrowsSecfile, float divFact,
            float *xlocFirst, float *ylocFirst, float *xlocSec, float *ylocSec, vector<float> &matchPixelXone,
            vector<float> &matchPixelYone, vector<float> &matchPixelXtwo, vector<float> &matchPixelYtwo,
            vector<float> &matchPixelDiffX, vector<float> &matchPixelDiffY);//matching  stars from 2 consecutive frames
//    int findShiftsNtheta(int totalelements, vector<float> &Xone, vector<float> &Yone, vector<float> &Xtwo, vector<float> &Ytwo,
//            vector<float> &DiffOfX, vector<float> &DiffOfY, double &Xdx, double &Ydy, double &Theta);//finding shifts and theta between two consecutive frame in registration and averaging
    int findShiftsNtheta (int totalelements , vector<float> &Xone , vector<float> &Yone , vector<float> &Xtwo , vector<float> &Ytwo ,
        vector<float> &DiffOfX , vector<float> &DiffOfY ,bool flag_theta_computation, double &Xdx , double &Ydy , double &Theta);//finding shifts and theta between two consecutive frame in registration and averaging
    int ApplySubSampling(float* inputarray, int in_xsize, int in_ysize, float* outputarray, int out_xsize, int out_ysize);//subsampling to convert 9600 frame to 600 frame
    int performDistortionCorr(float *frmdata, float *frmtempdata, float *Xdistr, float *Ydistr, int sizex, int sizey);//perform distortion correction 
  
    //for driver module
      vector<string> sciencedata_file,status_info;
     bool paramfile_Varask_iml2,tar_extracted_flag_IM;
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

#endif	/*UVTLEVEL2IM_H*/

