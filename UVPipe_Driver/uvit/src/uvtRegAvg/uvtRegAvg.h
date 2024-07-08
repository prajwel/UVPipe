/* 
 * File:   uvtRegAvg.h
 * Author: uvit
 *
 * Created on November 27, 2013, 3:31 PM
 */

#ifndef UVTREGAVG_H
#define	UVTREGAVG_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>

#define MODULENAME "uvtRegAvg"
#define caldbfinalsize 600 
 #define  TEMPLATE_SIZE  1024
#define NO_OF_ELEMENTS_TO_CMPR 2
typedef struct Star
{

    float intensity ;
    int x , y ;


    
} ;
class uvtRegAvg{
    private:
         char modulename[NAMESIZE];
        char inputdatadir[PIL_LINESIZE];
        char caldbDir[PIL_LINESIZE];
       char outdir[PIL_LINESIZE];   
        int clobber;
        int history;
        char mode[10];            //pil mode
         int nframes,numFrames2Gen;
         char infofile_in[PIL_LINESIZE];   //full path
        char infofile_out[PIL_LINESIZE];
        caldb_Handler caldb_handler;
        DataInfo datainfo;
        int xsize,ysize;
        int caldbsize;
          long  caldbdim;
           float mult_term;
            int min_stars_match;
        char sigframedir[PIL_LINESIZE];     //only signal frame directory name , not full path
        char expframedir[PIL_LINESIZE]; 
      //  char centroidframedir[PIL_LINESIZE];
      char **sigframelist;   //for input , used for IM only
      char **expoframelist;   //for input , used for IM only
      char distrotionCorrfile[PIL_LINESIZE];
      char nameprefix[FLEN_FILENAME]; 
      char moduleoutdir[PIL_LINESIZE];
      float caldb_xdist_final[caldbfinalsize*caldbfinalsize];
      float caldb_ydist_final[caldbfinalsize*caldbfinalsize];
      char eventfile[PIL_LINESIZE]; 
      float x_corr,y_corr;
      vector<string> key_records;
      float sd_mul_factor,sd_multi_factor_default;
      int FINALFRAMESIZE_REGAVG;
         /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
                //char starDir[PIL_LINESIZE];
                   vector<Star> star_track;
                char centroidDir[PIL_LINESIZE];
                char **centroidframelist;   //for input , used for IM only
                int algo_flag,minimum_No_of_Stars;
                int centroidlimit;
                float diff_dist;
                int option_LeastSquare;
                  int centroid_Winsize,refine_Winsize;
         vector<int> Fx,Fy,Rx,Ry;
         vector<float> Fval,Rval,Cx,Cy, Ci;
          int algo_Square_Size;
         float primary_threshold_Val,secondary_threshold_Val;
          int getHistory(vector<string> &vhistory);
        
       
        int reg_Avg();
        int readDestortionFile();              //reads bad pixel file and sets data in badpixdata;
      int  Applypadding(float  *inputArray,float *outputArray);
      int correlate(float *array1,int h1,int w1, float *array2,int h2,int w2,float xshift,float yshift,int bp_p1,int bp_p2);
      void StarDetectionAndCentroids(char infile[],char outfile[],int sqrSize,double primaryThreshold,double secondaryThreshold);
    int check_MaxCenterPixVal(double center,double **sqr,int sqr_size);
      public:
        uvtRegAvg();
        ~uvtRegAvg();
        int read(int argc,char **argv);
        void display();
        int uvtRegAvgProcess();
      int  read (char *inputdatadir  , char *outdir , int clobber , int history ,int algoflag, float sdMutltifactor ,int centlimit,int min_nostars,int refinewinsize,int centroidwinsize,int algosqrsize,float prithr,float secthr,float opLeastsqr,float diffDist);
        template<class T>
int match(T *searchwindow,int h1,int w1,T *templatewindow,int h2,int w2,float *x,float *y);
        int restart();
    int     findStar_algo1 (float *inputArray);
    void doCentroiding (vector<int> &X , vector<int> &Y , int centroidwindow , float *arr , int h , int w);
        const char *getModuleOutdir() const { return moduleoutdir; } 
        int ApplySubSampling(float *inputarray,int in_xsize,int in_ysize,float *outputarray,int out_xsize,int out_ysize);
          int findShiftsNtheta (int totalelements , vector<float> &Xone , vector<float> &Yone , vector<float> &Xtwo , vector<float> &Ytwo , 
                                                                        vector<float> &DiffOfX , vector<float> &DiffOfY , double &Xdx , double &Ydy , double &Theta);
    int matchStars (int numrowsFirstfile , int numrowsSecfile , float divFact , float *xlocFirst , float *ylocFirst ,
        float *xlocSec , float *ylocSec , vector<float> &matchPixelXone , vector<float> &matchPixelYone , vector<float> &matchPixelXtwo , vector<float> &matchPixelYtwo , vector<float> &matchPixelDiffX , vector<float> &matchPixelDiffY);
        
    };

#endif	/* UVTREGAVG_H */


