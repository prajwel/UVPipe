/* 
 * File:   uvtImRa_commonArray.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on December 9, 2013, 4:11 PM
 */

#ifndef UVTIMRA_COMMONARRAY_H
#define	UVTIMRA_COMMONARRAY_H
#include<string>
#include<iostream>
#include<stdlib.h>
#include<uvtUtils.h>
#include<caldb_Handler.h>
#include<DataInfo.h>
#include<map>
#include<DataIngest.h>

#define darkSize 512
#define IMG_DIM_FI 600
typedef struct FrameIntegration_Arr {
   vector<float> pixels_sig_fi; //Stores Signal image pixel values
    vector<float> pixels_exp_fi; //Stores Exposure image pixel values
    int frame_number_fi; //frame number
    double frame_time_fi; //frame time
};
class uvtImRa_commonArray{
    
    float *InputArray;
    string modulename;
    string level1indir;
    string caldbindir;
    string channel; 
    //char caldbindir[PIL_LINESIZE];
       int darkframe_flag;
         int algo_Square_Size;
    float primary_threshold_Val, secondary_threshold_Val;
     float err_per;
    string level2outdir;
     char mode[5];
    char stardir[10];
    char centroidDir[10];
    int clobber, history;
    int subDivision_size;
    int centroidlimit;
    int algo_flag,centroid_Winsize,refine_Winsize;
    float rms_mul_factor,rms_mul_factor_default;
    int frames_toDiscard,avg_Factor;
    int freqDomainFilter_Flag,fitting_flag, orderpitch,orderyaw,orderroll,no_ofFramesToAcc,type_Filtering;
    double threshold,freqvalue,delta_time;
    caldb_Handler  caldb_handler;
       DataInfo datainfo;
         char caldbDir[PIL_LINESIZE];
             float *badpixdata,*flatfielddata;
             float *badExparray;
             char badpixfile[PIL_LINESIZE];
             char flatfieldfile[PIL_LINESIZE];
              char distrotionCorrfile[PIL_LINESIZE];
              char opticalDistCorrfile[PIL_LINESIZE];
              char op_distrotionCorrfile[PIL_LINESIZE];
              char qeFile[PIL_LINESIZE];
                        char moduleoutdir_dd[NAMESIZE];
                  char moduleoutdir_bp[NAMESIZE];
                char moduleoutdir_uc[NAMESIZE];
                char moduleoutdir_ff[NAMESIZE];
                char moduleoutdir_qemcp[NAMESIZE];
            char  moduleoutdir_pp[NAMESIZE];
             char  moduleoutdir_sd[NAMESIZE];
            char  moduleoutdir_cr[NAMESIZE];
       char    moduleoutdir_ac[NAMESIZE];
      char  moduleoutdir_sc[NAMESIZE];
       char  moduleoutdir_od[NAMESIZE];
        char  moduleoutdir_de[NAMESIZE];
           char  moduleoutdir_rfc[NAMESIZE];
           char moduleoutdir_fi[NAMESIZE];
        char dstartpath[PIL_LINESIZE];
        
        char dendpath[PIL_LINESIZE];
        double t_darkframeend;
        double t_darkframestart;
           float diff_Dist;
         double t_curr;
                int xsize,ysize;
                int centroid_row;
                float *X_dist,*Y_dist,*X_odist,*Y_odist;
             long caldb_dim;
              int darkFrame_Flag;
                vector<int> Rx_vect,Ry_vect,Fx_vect,Fy_vect;
                vector<float> cx_vect,cy_vect,c_int,R_int,F_int,cx_tempx,cy_tempy,ci_tempi;
  int option_LeastSquare;
  vector<string> name_track ;
  int parity_flag;
 vector<float> x_shift,y_shift,theta_shift;
 vector<double> frm_time_data;
 int DATA_ARRAYSIZE;
   float darkFrameend_data[darkSize*darkSize];
        float darkFramestart_data[darkSize*darkSize];
        vector<float> x_ref_cumm,y_ref_cumm,int_ref_cumm;
          int gti_flag;
          int  all_or_custom,valid_bit;
         char in_Gtifile[FLEN_FILENAME];
        vector<float> roll_vect,pitch_vect,yaw_vect,x_Of_refFrame,y_Of_refFrame,int_Of_refFrame;
//        float darkFrameend_data[darkSize*darkSize];
//        float darkFramestart_data[darkSize*darkSize];
    //for dataingest
    int dropframe;                //yes for dropping frame if CRC fails, no for dropping packet if CRC fails
     double temperature;
     char lbtfile[FLEN_FILENAME];
   vector<int> size_rows;
   char darkdir[FLEN_FILENAME];
   char Indir_dataIngest[FLEN_FILENAME];
    vector<int> track_invalidPix ; 
    //for unit conversion
    int unitConversionFlag;    //YES if unit conversion needs to run else NO
    
    //for flat field correction module
    int flatfieldFlag; //YES if flat field correction is to be done, No if not required
    
    //for pixel padding
    int padding_dim;    //required dimension after padding
    
    //for Accumulate
    int Nacc;           //number of frames to accumulate
    
    //for corrrection of temperature effects
     int qemcpFlag;                           //YES to apply
     
    //for pixel subdivision      
    int subdivisionFlag;                    //YES to apply
    int Win_NB;
    int win_xsize,win_ysize;
    vector<string> header_info;
    map<int ,int> invalid_pix_duplicates_indices;
       float *temp,*f0,*f1,*f2,*f3,*f4,*f5,*f6,*f7;
    
    float *frame_Data,*frame_ExpData;
    float *frame_Data_Padded,*frame_ExpoData_padded,*frame_fc_data;
    float *frame_Data_subdivided, *frame_ExpData_subdivided;
    int wtd_uc,wtd_bp,wtd_ff,wtd_pp,wtd_sd,wtd_cr,wtd_ac,wtd_fi,wtd_fsc,wtd_dd,wtd_od,wtd_de,wtd_rfc,wtd_qemcp;
    float *darkCompute_array;
    char eventfile[PIL_LINESIZE]; //name of input event file
        float *darksubtr_array;
        char RelAspFilename[FLEN_FILENAME];
public:
   
    int setDirectoryStructure(char *Dir, const char *subdir) ;
    int writeOutputImageToDisk(char *id, char *outDir, char *dir , char *subscript,float *Array,char *namepre, double ftime, unsigned short fno,int sizex,int sizey);
   double*   DataFitting(double *t,double *X,int order,int nRecords);
        void  doFft(double data[],  long nn, int isign);
          int  transformToUVITFrame(double *t,double *r,double *p,double *y);
            int  transformToSpacecraftFrame(double *t,double *r,double *p,double *y);
             double  readDarkFrame(char * path, float*Array);
      int  darkFrameComputation(float *Array);
    // int darkFrameSubtraction(float *Array, float *frame_data);
     int takeDarkinfo ();
    uvtImRa_commonArray();
    ~uvtImRa_commonArray();
     
     char filter[FLEN_FILENAME];
   // int  performQEMCPcorrection(float *sigdata,int xsize,int ysize,double fact);
        int uvtImRacomArrProcess();
        int nFrameDiscard_fi,nFrameIntegrate_fi;
         int read(int argc, char** argv);
         int minimum_No_of_Stars;
           int match_stars_file_flag;
      long  nCalDBTempValues;
         long nrows_lbt ;
        double *time_lbt;
          float *insideTemp;
          float *outsideTemp;
        float *qe_mg_factor;
          //int  readQEMCPFile();
       //    int getTemp ();
     int  copyAllheaderKeys(char *infile);
   //  int   readImage(char * caldb_file,int hduno,float *caldb_data);
     int  readcaldbtable(char * caldb_file,int hduno,vector<float> &caldb_data);
     int performFrameIntegration (long nrows , unsigned short*frame_no , int xsize , int Ndiscard , int  Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float  *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr> &vect);
    // int performUnitConversion(float *frmdata, float *expdata,double  intgrntime,int sizex,int sizey);
    // int performCorrection(float *frmsigdata, float *frmexpdata, float *badpixarry,int sizex,int sizey,double  intgrntime);
    // int performCosmicRayCorr(float *frmsigdata, float *frmexpdata,int sizex,int sizey,float threshold_cr);
    // int performFlatFieldCorr(float *frmsigdata, float *frmflatfield,int sizex,int sizey);
    // int performSubDivision(float *frmsigdata,int sizex,int sizey,float * subdivideddata,int size_subdivx,int size_subdivy);
 int writeOutputTblToDisk (char *id , char *outDir , char *dir , char *subscript  , char *namepre , double ftime , unsigned short fno ,char **type1,char**tform1,int tfields1,char *extname,vector<float> &X ,vector<float> &Y,vector<float> &val);
//void     getKeywordValue (char *file , int type,char *keyname , int hdunum , char *val);
int performAccOFFrame(float *frmsigdata,float *sumdata,int sizex,int sizey,int numoffrmAcc);
int performDistortionCorr(vector<float> &X,vector<float> &Y,float *Xdistr,float *Ydistr,int sizex,long caldbsize);
int calculateRefFrame(vector<int> &sizevect,vector<float> &xref,vector<float> &yref,vector<float> &intref, int avgfact,vector<float> &xrefcumm,vector<float> &yrefcumm,vector<float> &intrefcumm);

};

#endif	/* UVTIMRA_COMMONARRAY_H */

