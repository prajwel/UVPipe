/* 
 * File:   uvtComputeJitter.h
 * Author: uvit
 *
 * Created on November 26, 2013, 1:32 PM
 */

#ifndef UVTCOMPUTEJITTER_H
#define	UVTCOMPUTEJITTER_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>

#define MODULENAME "uvtComputeJitter"
#define caldbfinalsize 600 
#define IMAGE_ARRAYSIZE 512


class uvtComputeJitter 
{
private:
    char modulename[NAMESIZE];
    char inputdatadir[PIL_LINESIZE];
    char caldbDir[PIL_LINESIZE];
    char outdir[PIL_LINESIZE];
    int clobber;
    int history;
    char mode[10]; //pil mode
   
    char infofile_in[PIL_LINESIZE]; //full path
    char infofile_out[PIL_LINESIZE];
    char input_driftFile[PIL_LINESIZE];
    int FreqDomainFilter_Flag, type_Filtering;
    caldb_Handler caldb_handler;
  
    int orderpitch, orderyaw, orderroll, fittingflag;
    DataInfo datainfo;
    int xsize, ysize;
    int caldbsize;
    double freqvalue;
    long caldbdim;
    float *x_Distortion, *y_Distortion;
    char distrotionCorrfile[PIL_LINESIZE];
    char nameprefix[FLEN_FILENAME];
    char moduleoutdir[PIL_LINESIZE];
    char gyrofile[PIL_LINESIZE];
    int gyroData_power2;
    long gyro_rows;
    int nrows_gyro;
    double *time;
    
    double *G1_R, *G1_P, *G2_Y, *G2_P, *G3_Y, *G3_R;
    vector<double> G1_R_final, G1_P_final, G2_Y_final, G2_P_final, G3_Y_final, G3_R_final, time_final_gyro;
    vector<double> time_gyro, integrated_timeGyro, integrated_r, integrated_y, integrated_p;
    vector<double> delta_x_final, delta_y_final, delta_theta_final, time_final;
    float pix_perdeg_X,pix_perdeg_Y;
  
    float scaleChannel, eulerangle_one, eulerangle_two, eulerangle_three;

    /**
     * Function to generate history for the module, to be written to the output file
     * @param vhistory
     * @return 
     */
    

    char centroidDir[PIL_LINESIZE];
    char **centroidframelist; //for input , used for IM only
    float *dx, *dy;
   
    int getHistory(vector<string> &vhistory);
    int computeJitter();
   
    double *time_data;
    int readGyroFile(); //reads bad pixel file and sets data in badpixdata;
    int extractGyro(long nrows);
    //   int  Applypadding(float  *inputArray,float *outputArray);
    int transformToUVITFrame(double *t, double*r, double *p, double *y);
    int transformToSpacecraftFrame(double *t, double *r, double *p, double *y);
    vector<double> vect_beforeExtr_x, vect_beforeExtr_y, vect_beforeExtr_theta, vect_beforeExtr_time;
    vector<string> key_records;
    
public:
    uvtComputeJitter();
    ~uvtComputeJitter();
    int read(int argc, char **argv);
    int read (char *inputdatadir , char *caldbDir,char *gyro_file , int freqDomainFilter_Flag , double freqvalue , int fitting_flag , int orderpitch , int orderyaw , int orderroll , char *outdir , int Type_Filtr,int clobber , int history);
    void display();
    int uvtComputeJitterProcess();
    int doIntegration();
    int doFiltering();
 
    int writeJitter();
    const char *getModuleOutdir() const {
        return moduleoutdir;
    }

    
    };

#endif	/* UVTCOMPUTEJITTER_H */

