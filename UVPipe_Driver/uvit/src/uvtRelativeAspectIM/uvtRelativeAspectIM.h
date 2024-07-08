/* 
 * File:   uvtImRa_commonArray.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on December 9, 2013, 4:11 PM
 */

#ifndef UVTRELATIVEASPECTIM_H
#define	UVTRELATIVEASPECTIM_H
#include<string>
#include<iostream>
#include<stdlib.h>
#include<uvtUtils.h>
#include<caldb_Handler.h>
#include<DataInfo.h>
#include<map>

#define darkSize 512

class uvtRelativeAspectIM {
    int  zipFlag;
    float *InputArray; //
int crc_flag;
    string modulename; //
    string level1indir; //Level-1 directory path from parameter file uvtRelativeAspectIM.par
    string caldbindir; //CalDB directory (CALDB/uvit) path from parameter file uvtRelativeAspectIM.par
    string level2outdir; //Output directory path from parameter file uvtRelativeAspectIM.par
    string channel; //FUV / NUV / VIS
    caldb_Handler caldb_handler; //object for  caldb_handler
    vector<int> track_invalidPix; //vector for storing the invalid pixel 
    DataInfo datainfo; //object for DataInfo
    vector<string> header_info; //Stores Level-1 header information  
    long counter_JunkFrame;
  
    int wtd_uc, wtd_bp, wtd_ff, wtd_pp, wtd_sd, wtd_cr, wtd_ac, wtd_fsc, wtd_dd, wtd_od, wtd_de, wtd_rfc, wtd_qemcp; //flags for writing outputs to disk
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
    int clobber, history; //for overwrite and history respectively in FITS file from parameter file 


    char mode[5]; //for PIL (query-learn(ql), auto(a), hidden(h))
    //for pipeline use
    int UTC_flag;
    float *frame_Data, *frame_ExpData;
    float *frame_Data_Padded, *frame_ExpoData_padded, *frame_fc_data;
    float *frame_Data_subdivided, *frame_ExpData_subdivided;
    int xsize, ysize; //Width and Height of the window (e.g (100, 100), (512, 512))

    //for DataIngest
    int gti_flag;
    int junkFrameFlag;
    int all_or_custom, valid_bit;
    char in_Gtifile[FLEN_FILENAME];
    int dropframe; //yes for dropping frame if CRC fails, no for dropping packet if CRC fails
    char Indir_dataIngest[FLEN_FILENAME];
    int win_xsize, win_ysize; //Width and Height of the window (e.g (100, 100), (512, 512))
    int parity_flag; //Parity flag(no need for IM mode)    
    char darkdir[FLEN_FILENAME];
    int att_flag_val;
    int  flag_thetaComp;

    //for unitconversion
    int darkframe_flag;
    float *darkCompute_array;
    int unitConversionFlag; //YES if unit conversion needs to run else NO
    float darkFrameend_data[darkSize*darkSize]; //Stores pixels of dark frame for end frame
    float darkFramestart_data[darkSize*darkSize]; //Stores pixels of dark frame for start frame
    char dstartpath[PIL_LINESIZE]; //stores  path for dark start path.
    char dendpath[PIL_LINESIZE]; //stores path for dark end path.
    double t_curr; //stores time of current frame in processing
    double t_darkframeend; //stores dark start path
    double t_darkframestart; //stores dark end path
     

    //for bad pixel 
    char badpixfile[PIL_LINESIZE]; //CalDB bad pixel Filename
    float *badpixdata; //stores caldb data for bad pixel correction

    //for flat field correction 
    int flatfieldFlag; //YES if flat field correction is to be done, No if not required
    char flatfieldfile[PIL_LINESIZE]; //CalDB flat field Filename
    float *flatfielddata; //stores caldb data for flat field correction

    //for QE  and MCP correction
    char filter[FLEN_FILENAME]; //Stores filter name (F0, F1 ....)  
    int qe_mcpFlag; //YES to apply
    char qe_factor_file[PIL_LINESIZE];
    double temperature;
    float *cal_temperature, *cal_f0, *cal_f1, *cal_f2, *cal_f3, *cal_f4, *cal_f5, *cal_f6, *cal_f7;//caldb temperaure and filter values 
    long nCalDBTempValues; //Number of Temperature values (currently 20 in CalDB)
    long nrows_lbt; //Number of rows in LBT file
    double *time_lbt; //Time in LBT for each row
    float *insideTemp; //Inside temperature
    float *outsideTemp; //Outside temperature
    float *qe_mg_factor; //QE & MCP factor 
    char lbtfile[FLEN_FILENAME];//path of LBT file
    
    //for pixel padding
    int padding_dim; //required dimension after padding

    //for Accumulation 
    int no_ofFramesToAcc; //number of frames to accumulate
    // int no_ofFramesToAcc; //Number of framed to be accumulated in case of IM mode

    //for subdivision
    int subdivisionFlag; //YES to apply
    int subDivision_size;

    //for star detection
    int algo_Square_Size;
    float primary_threshold_Val; //Threshold in star detection for Joe's algorithm
    float secondary_threshold_Val; //Threshold in star detection for Joe's algorithm
    float backgrd_fact; //value of the background 
    float diff_dist; //Neighbourhood for star pairing/searching
    char stardir[10]; //Star directory path to store framewise star (first cut and refined peaks)
    char centroidDir[10]; //Star directory path to store framewise star (first cut and refined peaks)
    int star_detect_algo_flag, centroid_Winsize, refine_Winsize;
    float rms_mul_factor, rms_mul_factor_default,thrJunkFrame;
    vector<int> Rx_vect, Ry_vect, Fx_vect, Fy_vect;
    vector<float> cx_vect, cy_vect, c_int, R_int, F_int, cx_tempx, cy_tempy, ci_tempi;
    int search_win_size; //Window size for searching stars in consecutive frames
    int centroid_row; //for storing the number of stars found after centroiding in star detection
    float min_num_stars; //minimum number of stars for Star detection


    //for reference frame calculation
    int frames_toDiscard; //Number of frames to be discarded in reference frame calculation
    int nFrameToAverage; //Number of frames to be averaged in reference frame calculation
    //Stores x,y and intensity values for all frames used for reference frame calculation
    vector<float> x_Of_refFrame, y_Of_refFrame, int_Of_refFrame;
    //Stores average x,y and intensity values for all frames used for reference frame calculation
    vector<float> x_ref_cumm, y_ref_cumm, int_ref_cumm;
    vector<int> size_rows; //Stores the number of stars in each frame considered for averaging in reference frame calculation
    vector<string> ref_frame_module_filename_track; //Keeps track of filename for output of reference frame calculation
    vector<double> ref_frame_time_data; //Stores frame time for frames in reference frame calculation

    //for detector distortion  correction 
    char detector_distortion_corr_file[PIL_LINESIZE]; //CalDB detector distortion Filename
    float *X_detect_distortion, *Y_detect_distortion; //stores detector distortion correction values from CalDB

    //for optical distortion correction file
    char optical_distortion_corr_file[PIL_LINESIZE]; //CalDB optical distortion Filename
    float *X_optical_distortion, *Y_optical_distortion; //stores optical distortion correction values from CalDB

    //for drift computation
    int freqDomainFilter_Flag; ////Frequency domain filter flag
    int fitting_flag; //Polynomial fitting flag for roll, pitch and yaw
    int orderpitch, orderyaw, orderroll; // Polynomial fitting order for roll, pitch and yaw  
    int type_Filtering; //1-Polynomial fitting & 0-Spatial domain based on cutt-off frequency
    double threshold, freqvalue, delta_time;
    int shift_rotation_algo; //Flag for shift and rotation computation algorithm (1-SAC, 2-SAC, 3-IUCAA)
    vector<float> x_shift, y_shift, theta_shift; //Stores x-y shifts and rotations    
    int match_stars_file_flag; //flag to decide whether match star txt file to be generated or not in drift computation
float err_per;  //May be removed after checking
float nbhd_dist;
    //for relative aspect calculation chain
    char RelAspFilename[FLEN_FILENAME];


public:

    int readPILParameters(int argc, char** argv); //Reads PIL parameters
    int setDirectoryStructure(char *Dir, const char *subdir); //Makes directory structure
    int writeOutputImageToDisk(char *id, char *outDir, char *dir, char *subscript, float *Array, char *namepre, double ftime, unsigned short fno, int sizex, int sizey);
    int copyAllheaderKeys(char *infile); //Copies header from Level-1 file and stores in header info
    int writeOutputTblToDisk(char *id, char *outDir, char *dir, char *subscript, char *namepre, double ftime, unsigned short fno, char **type1, char**tform1, int tfields1, char *extname, vector<float> &X, vector<float> &Y, vector<float> &val);
    int transformToUVITFrame(double *t, double *r, double *p, double *y);
    int transformToSpacecraftFrame(double *t, double *r, double *p, double *y);
    double readDarkFrame(char * path, float*Array); //Reads output dark frames from DataIngest
    int darkFrameComputation(float *Array);
    int takeDarkinfo();

    uvtRelativeAspectIM(); //constructor
   ~ uvtRelativeAspectIM(); //destructor
    int uvtRelativeAspectIMProcess(); //Performs Relative Aspect Calculation
    int performAccOFFrame(float *frmsigdata, float *sumdata, int sizex, int sizey, int numoffrmAcc); //performing Accumulation of number of frames
    int performDistortionCorr(vector<float> &X, vector<float> &Y, float *Xdistr, float *Ydistr, int sizex, long caldbsize); //perform distortion correction
    int calculateRefFrame(vector<int> &sizevect, vector<float> &xref, vector<float> &yref, vector<float> &intref, int avgfact, vector<float> &xrefcumm, vector<float> &yrefcumm, vector<float> &intrefcumm); //calculate reference frame
    int getHistory (vector<string> &vhistory);
    //for driver module
     string orbnum;//stores orbit number information.
    bool paramfile_Varask,tar_extracted_flag_IM;
     vector<string> sciencedata_file, status_info;
        vector<string>  getSciencedatafile_info(){
                        return sciencedata_file;
    };
        
     vector<string> getStatus_info(){
         return status_info;
     };
      
     const char * getRelAspFilename()
    {
        return this->RelAspFilename;
    }
   
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

#endif	/* UVTRELATIVEASPECTIM_H */

