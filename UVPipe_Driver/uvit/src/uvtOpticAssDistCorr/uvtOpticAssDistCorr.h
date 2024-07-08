/* 
 * File:   uvtOpticAssDistCorr.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#ifndef UVTOPTICALDESTORTION_H
#define	UVTOPTICALDESTORTION_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>

#define MODULENAME "uvtOpticAssDistCorr"
#define caldbfinalsize 600 
class uvtOpticAssDistCorr{
    private:
        char modulename[NAMESIZE];
        char inputdatadir[PIL_LINESIZE];
        char caldbDir[PIL_LINESIZE];
        //char eventfile[PIL_LINESIZE];
        char filter[3];
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
        char eventfile[PIL_LINESIZE];
        char distortionCorrfile[PIL_LINESIZE];
        char nameprefix[FLEN_FILENAME]; 
        char moduleoutdir[PIL_LINESIZE];
        vector<string> key_records;
        /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
                //char starDir[PIL_LINESIZE];
                char centroidDir[PIL_LINESIZE];
                char **centroidframelist;   //for input , used for IM only
            
                 float x_Distortion[CALDB_DIST_SIZE*CALDB_DIST_SIZE], y_Distortion[CALDB_DIST_SIZE*CALDB_DIST_SIZE];
                float *dx,*dy;
                int *X,*Y;
                 long caldbdim;
                 float caldb_xdist_final[caldbfinalsize*caldbfinalsize];
         float caldb_ydist_final[caldbfinalsize*caldbfinalsize];
                int getHistory(vector<string> &vhistory);
        
       vector<int> vex,vey;
       vector<float>vedx,vedy;
        int opticalDistortionIM();
          int opticalDistortionPC();
        int readDistortionFile();               //reads bad pixel file and sets data in badpixdata;
        
     public:
        uvtOpticAssDistCorr();
        ~uvtOpticAssDistCorr();
        int read(int argc,char **argv);
        int read(char *inputdatadir,char *caldbDir,char *outdir,int clobber,int history);
        void display();
        int uvtOpticalDistcorrProcess();
      //  int  Applypadding(float  *inputArray,float *outputArray);
        const char *getModuleOutdir() const { return moduleoutdir; } 
        
    };

#endif	/* UVTOPTICALDESTORTION_H */

