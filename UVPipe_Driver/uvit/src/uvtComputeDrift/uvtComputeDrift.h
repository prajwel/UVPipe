/* 
 * File:   uvtComputeDrift.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#ifndef UVTDRIFTCOMPUTATION_H
#define	UVTDRIFTCOMPUTATION_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>

#define MODULENAME "uvtComputeDrift"
#define caldbfinalsize 600 
//#define DATA_ARRAYSIZE 8
#define IMAGE_ARRAYSIZE 600
#define pi 4*atan(1.0)

class uvtDriftComputation {
private:
    char modulename[NAMESIZE];
    char inputdatadir[PIL_LINESIZE];
    char caldbDir[PIL_LINESIZE];
     char attFile[PIL_LINESIZE];
    //char eventfile[PIL_LINESIZE];
    double freqvalue;
    char outdir[PIL_LINESIZE];
    int clobber;
    int history;
    char mode[10]; //pil mode
    int nframes;
    char infofile_in[PIL_LINESIZE]; //full path
    char infofile_out[PIL_LINESIZE];
    int FreqDomainFilter_Flag, type_Filtering;
    caldb_Handler caldb_handler;
    DataInfo datainfo;
    int xsize, ysize;
    float scaleXCH1, scaleYCH1,scaleXVIS,scaleYVIS, Shift_X, Shift_Y, angle_NUVtoVIS;
    int orderpitch, orderyaw, orderroll;
    int fittingflag;
    int caldbsize;
    int option_LeastSquare;
    float diff_dist,diff_dist_cp;
    char nameprefix[FLEN_FILENAME];
    char infile_drift[FLEN_FILENAME];
    char moduleoutdir[PIL_LINESIZE];
    float caldb_xdist_final[caldbfinalsize*caldbfinalsize];
    float caldb_ydist_final[caldbfinalsize*caldbfinalsize];
    double delta_time;
    char eventfile[PIL_LINESIZE];
    int nframes_power2;
    vector<double> x_arr, y_arr, theta_arr;
    double *Xshift_arr, *Yshift_arr, *theta_arrfinal;
    double *time;
    int minimum_noTargetStars;
    float err_per;
    int algo_flag_value;
    float div_fact;
    vector<string> key_records;
    vector<string> vhistorystr;
    int match_starsfile_gen;
int  flag_theta;
    /**
     * Function to generate history for the module, to be written to the output file
     * @param vhistory
     * @return 
     */
    //char starDir[PIL_LINESIZE];
    char centroidDir[PIL_LINESIZE];
    char **centroidframelist; //for input , used for IM only


    float *dx, *dy;

    vector<double> x_shift_vect, y_shift_vect, theta_vect;
    int getHistory(vector<string> &vhistory);

    int findDrift_WithAlgo_1();
    int findDrift_WithAlgo_3();
    int writeDrift();
int cntglobal;

public:
    uvtDriftComputation();
    ~uvtDriftComputation();
    int read(int argc, char **argv);
    
//    int read (char *inputdatadir  ,char * attFile_in,float percent , float diff_Dist , int freqDomainFilter_Flag , int type_filter ,
//        double d_time , double freqvalue , int fitting_flag , int orderpitch , int orderyaw , int orderroll ,
//        char *outdir , int op_leastsquare , int flag_matchstars, int algoFlag,int clobber ,int history);
    
    
    int read (char *inputdatadir  ,char * attFile_in,float percent , float diff_Dist , int freqDomainFilter_Flag , int type_filter ,
        double d_time , double freqvalue , int fitting_flag , int orderpitch , int orderyaw , int orderroll ,
        char *outdir , int op_leastsquare , int flag_matchstars, int algoFlag,int flag_thetacomp,int clobber ,int history);
    int read_pc(char *inputdatadir, float diff_Dist, int freqDomainFilter_Flag, int type_filter, double d_time, double freqvalue, int fitting_flag, int orderpitch, int orderyaw, int orderroll, char *outdir, int op_leastsquare, int clobber, int history);
    int local_Smoothening(vector<double> &time_data, vector<double> &delta_term, int order, double delta_time, vector<double> &delta_term_final, vector<double> &time_data_final);
    void display();
    int uvtDriftComputationProcess();
        

    const char *getModuleOutdir() const {
        return moduleoutdir;
    }
//  int removeRecords(vector<float>  &Xone , vector<float> &Yone ,vector<float> &Xtwo , vector<float> &Ytwo,vector<float> &DiffOfX , vector<float> &DiffOfY ,float *ints1,float*ints2,
//   vector<float> &newXone,vector<float> &newYone,vector<float> &newXtwo,vector<float> &newYtwo,vector<float> &newDiffX,vector<float> &newDiffY,vector<float> &new_one_ints,vector<float> &new_two_ints);
    int removeRecords(vector<float>  &Xone , vector<float> &Yone ,vector<float> &Xtwo , vector<float> &Ytwo,vector<float> &DiffOfX , vector<float> &DiffOfY ,vector<float> &ints1,vector<float> &ints2,
   vector<float> &newXone,vector<float> &newYone,vector<float> &newXtwo,vector<float> &newYtwo,vector<float> &newDiffX,vector<float> &newDiffY,vector<float> &new_one_ints,vector<float> &new_two_ints);
   
    int doFiltering(vector<double> &xshiftvect, vector<double> &yshiftvect, vector<double> &thetavect);
//    int findShiftsNtheta (int totalelements ,vector<float>  &Xone , vector<float> &Yone , vector<float> &Xtwo , vector<float> &Ytwo , 
//                                                                        vector<float> &DiffOfX , vector<float> &DiffOfY , double &Xdx , double &Ydy , double &Theta);
//     int findShiftsNtheta (int totalelements , vector<float>  &Xone , vector<float> &Yone , float *int1, vector<float> &Xtwo , vector<float> &Ytwo , float *ints2,
//                                                                        vector<float> &DiffOfX , vector<float> &DiffOfY ,bool flag_theta_computation, double &Xdx , double &Ydy , double &Theta,double &minEleX,double &maxeEleX,double &minEleY,double &maxEleY   );
    //int matchStars(int numrowsFirstfile, int numrowsSecfile, float divFact, float *xlocFirst, float *ylocFirst, float *xlocSec, float *ylocSec, vector<float> &Xref, vector<float> &Yref, vector<float> &Xarr, vector<float> &Yarr, vector<float> &tempXarr, vector<float> &tempYarr);
     
  int findShiftsNtheta (int totalelements , vector<float>  &Xone , vector<float> &Yone , vector<float> int1, vector<float> &Xtwo , vector<float> &Ytwo ,vector<float> ints2,
                                                                        vector<float> &DiffOfX , vector<float> &DiffOfY ,bool flag_theta_computation, double &Xdx , double &Ydy , double &Theta,double &minEleX,double &maxeEleX,double &minEleY,double &maxEleY);
    
    int matchStars (int numrowsFirstfile , int numrowsSecfile , float divFact , float *xlocFirst , float *ylocFirst ,
        float *xlocSec , float *ylocSec , float * int1,float *int2,vector<float> &matchPixelXone , vector<float>  &matchPixelYone ,
        vector<float> &matchPixelXtwo , vector<float>  &matchPixelYtwo , vector<float> &int_one,vector<float> &int_two, 
        vector<float> &matchPixelDiffX , vector<float> &matchPixelDiffY);
     int local_Smoothening (vector<double> &time , vector<double> &data, int order,double delta_time);
      int local_Smoothening_new (vector<double> &time , vector<double> &data, int order,double delta_time);
    
};

#endif	/* UVTDRIFTEXERSISE_H */

