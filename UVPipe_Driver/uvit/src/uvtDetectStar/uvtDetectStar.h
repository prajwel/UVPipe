/* 
 * File:   uvtDetectStar.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */


#ifndef UVTDETECTSTAR_H
#define	UVTDETECTSTAR_H
//#define starDir star
//#define centroidDir centroid
#include<pil.h>
#include<fitsio.h>
#include<uvtUtils.h>
#include<DataInfo.h>
#include<spMatrix.h>

using namespace std;
typedef struct Star
{
    float intensity ;
   float  x , y ;
};
class uvtDetectStar{
   private:
       
       char modulename[NAMESIZE];
       
        char inputdatadir[PIL_LINESIZE];
        char outdir[PIL_LINESIZE];   
        int algo_flag;
        float thr_intensity_refinement;                          //flag=1 for algorithm 1 and flag=2 for algorithm 2, algorithm 1 uses rms and windowing based star detection
        int clobber;
        int history;
        float sd_mul_factor,sd_multi_factor_default,sd_mul_factor_second;
        int centroidlimit;     //number of maximum centroids to be accepted for further processing
        char mode[10];
        vector<string> key_records;
        int cnt_det;
        char moduleoutdir[NAMESIZE];
         char infofile_in[PIL_LINESIZE];
         char infofile_out[PIL_LINESIZE];
         DataInfo datainfo;
         int xsize,ysize;
         char nameprefix[FLEN_VALUE];
        char sigframedir[PIL_LINESIZE];        //name of signal frame directory
        char expoframedir[PIL_LINESIZE];
        long nframes;     //used in case of IM
         char **sigframelist;   //for input , used for IM only
          char **expoframelist;   //for input , used for IM only
          int win_search;
         char starDir[PIL_LINESIZE];
         char centroidDir[PIL_LINESIZE];
                  
         int centroid_Winsize,refine_Winsize,Nacc;
         float backgrnd_fact;
         
         vector<int> Fx,Fy,Rx,Ry;
         vector<float> Fval,Rval,Cx,Cy, Ci,cx_ref_vect,cy_ref_vect,ci_ref_vect;
          int algo_Square_Size;
         float primary_threshold_Val,secondary_threshold_Val;
        int getHistory(vector<string> &vhistory);
        int findStar_algo1(float *inputArray); // background  based 
        int findStar_algo2(float *inputArray); //Joe's algorithm
        int findStar_algo3 (float *inputArray );//SAC algorithm
        int findStar_algo4 (float *inputArray );
        int findStar_algo_temp (float *inputArray );
     //   findStar_algo2();
//        int detectStarPC();
        int detectStars();
        bool first_frameFlag;
             vector<Star> star_track;
        
     public:
        uvtDetectStar();
        ~uvtDetectStar();
          
        // char *starDir[10]="Star";
        //char *centroidDir[10]="Centroid";
   //     string starDir="star";
     //   string centroidDir="centroid";
        int matchStars_TorefFrame (vector<float> &x,vector<float> &y,vector<float> &intensity,float *inputArray);
        int read(int argc,char **argv);
        int read (char* inputdatadir , char* outdir , int algo_flag , float threshold , int refine_window , int centroid_window , float  num_min_stars , float prm_thr , float sec_thr , int algo_square , int windw_search,int clobber , int history);
        void display();
        int uvtDetectStarProcess();
        int check_MaxCenterPixVal(double center,double **sqr,int sqr_size);
        void doCentroiding(vector<int> &X, vector<int> &Y, int centroidwindow, float *arr, int h,int w);
         //int  getStarCentroid(float* inputArray, int nacc,float bckgrd_fct, float rms_factor, int size_x, int size_y, float winsize, float cent_winsize, vector<int> &rx, vector<int> &ry, vector<float> &cx, vector<float> &cy, vector<int>  &fx,vector<int> &fy ,vector<float> &rint, vector<float> &cint,vector<float> &fint,char *mode, int min_No_of_stars,bool flag_check,int Window_nh);
       int getStarCentroid (float* inputArray  , int nacc,float bckgrd_fct,float rms_factor , int size_x ,
        int size_y , float winsize , float cent_winsize  , vector<int> &rx , vector<int> &ry ,
        vector<float> &cx , vector<float> &cy , vector<int> &fx , vector<int> &fy , vector<float> &rint ,
        vector<float> &cint , vector<float> &fint,char *mode, float min_No_of_stars,bool flag_check,int Window_nh,int algoflg);
       
        const char *getModuleOutdir() const { return moduleoutdir; } 
        static bool  glob_static_var;
         
    };


#endif	/* UVTDETECTSTAR_H */

