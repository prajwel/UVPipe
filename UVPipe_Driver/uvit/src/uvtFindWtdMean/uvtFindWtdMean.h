/* 
 * File:   uvtFindWtdMean.h
 * Author: uvit
 *
 * Created on November 27, 2013, 11:45 AM
 */

#ifndef UVTFINDWTDMEAN_H
#define	UVTFINDWTDMEAN_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>

#define MODULENAME "uvtFindWtdMean"
#define caldbfinalsize 600 
class uvtFindWtdMean{
    private:
        char modulename[NAMESIZE];
        char inputdatadir[PIL_LINESIZE];
        char caldbDir[PIL_LINESIZE];
        char outdir[PIL_LINESIZE];   
        int clobber;
        int history;
        char mode[10];            //pil mode
        int nframes;
        char infofile_in[PIL_LINESIZE];   //full path
        char infofile_out[PIL_LINESIZE];
        caldb_Handler caldb_handler;
        DataInfo datainfo;
        int xsize,ysize;
        int caldbsize;
        long  caldbdim;
        float *x_Distortion,*y_Distortion;
        char distrotionCorrfile[PIL_LINESIZE];
        char nameprefix[FLEN_FILENAME]; 
        char moduleoutdir[PIL_LINESIZE];
        float caldb_xdist_final[caldbfinalsize*caldbfinalsize];
        float caldb_ydist_final[caldbfinalsize*caldbfinalsize];
        char eventfile[PIL_LINESIZE]; 
        int  no_ofWeigh;
         /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
                //char starDir[PIL_LINESIZE];
                char centroidDir[PIL_LINESIZE];
                char **centroidframelist;   //for input , used for IM only
                char **signalframelist;
                char **exposureframelist;
                char sigframedir[PIL_LINESIZE];        //name of signal frame directory
        char expoframedir[PIL_LINESIZE];
                 int index1;
                float *dx,*dy;
             //   int *X,*Y;
                int getHistory(vector<string> &vhistory);
        
       
        int findWtdMeanIM();
          int findWtdMeanPC();
        
      //  int readDestortionFile();              //reads bad pixel file and sets data in badpixdata;
      int  Applypadding(float  *inputArray,float *outputArray);
     public:
       
         uvtFindWtdMean();
         
        ~uvtFindWtdMean();
        vector<string> key_records;
        int read(int argc,char **argv);
     int   read(char *inputdatadir, char *outdir, int num_ofweight,int clobber, int history);
        void display();
        int uvtFindWeightedMeanProcess();
        
        const char *getModuleOutdir() const { return moduleoutdir; } 
        
    };

#endif	/* UVTFINDWTDMEAN_H */

