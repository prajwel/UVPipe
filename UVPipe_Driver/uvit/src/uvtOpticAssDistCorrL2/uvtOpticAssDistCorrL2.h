/* 
 * File:   uvtOpticAssDistCorrL2.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on November 27, 2013, 2:15 PM
 */

#ifndef UVTOPTICALDESTORTION_L2_H
#define	UVTOPTICALDESTORTION_L2_H


#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>

#define MODULENAME "uvtOpticAssDistCorr_l2"
#define caldbfinalsize 600 
class uvtOpticalDestortion_l2{
    private:
        char modulename[NAMESIZE];
        char inputdatadir[PIL_LINESIZE];
        char caldbDir[PIL_LINESIZE];
        //char eventfile[PIL_LINESIZE];
       
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
        char sigframedir[PIL_LINESIZE];     //only signal frame directory name , not full path
        char expframedir[PIL_LINESIZE];  
        char **sigframelist;   //for input , used for IM only
          char **expoframelist;   //for input , used for IM only
        char distortionCorrfile[PIL_LINESIZE];
        char nameprefix[FLEN_FILENAME]; 
        char moduleoutdir[PIL_LINESIZE];
        float *caldb_xdist_final;
         float *caldb_ydist_final;
     char eventfile[PIL_LINESIZE]; 
         /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
                //char starDir[PIL_LINESIZE];
                char centroidDir[PIL_LINESIZE];
                char **centroidframelist;   //for input , used for IM only
            
               
                float *dx,*dy;
             //   int *X,*Y;
                int getHistory(vector<string> &vhistory);
        
       
        int opticalDistortionIM();
          int opticalDestrotionPC();
        
    //    int readDestortionFile();              //reads bad pixel file and sets data in badpixdata;
    //  int  Applypadding(float  *inputArray,float *outputArray);
     public:
        uvtOpticalDestortion_l2();
        ~uvtOpticalDestortion_l2();
        int read(int argc,char **argv);
        int read(char *inputdatadir,char *caldbDir,char *outdir,int clobber,int history);
        void display();
        int uvtOpticaDistoProcess();
        
        const char *getModuleOutdir() const { return moduleoutdir; } 
        
    };

#endif	/* UVTOPTICALDESTORTION_L2_H */

