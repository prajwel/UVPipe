/* 
 * File:   uvtShiftRot.h
 * Author: uvit
 *
 * Created on November 27, 2013, 9:36 AM
 */

#ifndef UVTSHIFTROT_H
#define	UVTSHIFTROT_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
using namespace std;

#define MODULENAME "uvtShiftRot"
//#define IMG_DIM 512
#define pi 4*atan(1.0)
#define IMAGE_ARRAYSIZE 600
#define STARTFRAME 1
#define ENDFRAME 10
#define centroidframedir Centroid
class uvtShiftRot{
   private:
        char modulename[NAMESIZE];
        
        char inputdatadir[PIL_LINESIZE];
     
        char outdir[PIL_LINESIZE];   
        int clobber;
        int history;
        char mode[10];            //pil mode
       // int padding_dimension;
         char **centroidframelist; 
        char moduleoutdir[PIL_LINESIZE];   //full path of module output directory
        char sigframedir[PIL_LINESIZE];     //only signal frame directory name , not full path
        char expframedir[PIL_LINESIZE];  //only exposure frame directory name
        double threshold_value,Threshold;
        int nCompareFrames ;
        char rasfile[PIL_LINESIZE];
        char eventfile[PIL_LINESIZE];
        char infofile_in[PIL_LINESIZE];    //full path
        char infofile_out[PIL_LINESIZE];
        long no_of_records;
        int nframes;
        vector<string> key_records;
        char **sigframelist;   //for input , used for IM only
          char **expoframelist;   //for input , used for IM only
        int xsize,ysize;
         double *time , *roll_ras , *pitch_ras , *yaw_ras ;
    double *roll , *pitch , *yaw ;
    double *delta_x , *delta_y , *delta_theta ;
        char nameprefix[FLEN_FILENAME]; 
        DataInfo datainfo;
        
        /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
        int getHistory(vector<string> &vhistory);
     int    readRASfile();
        int uvtshiftrotatePC();
        int uvtshiftrotateIM();
      //  int readFlatFieldFile();               //reads bad pixel file and sets data in badpixdata;
        
     public:
        uvtShiftRot();
        ~uvtShiftRot();
        int read(int argc,char **argv);
        int read(char *inputdatadir,char *rasFile,char *outdir,int clobber,int history);
        void display();
        int uvtShiftRotProcess();
          const char *getModuleOutdir() const { return moduleoutdir; } 
  
    };


#endif	/* UVTSHIFTROT_H */

