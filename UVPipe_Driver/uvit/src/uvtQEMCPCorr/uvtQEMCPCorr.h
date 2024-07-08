/* 
 * File:   uvtQEMCPCorr.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#ifndef UVTQEMCPCORR_H
#define	UVTQEMCPCORR_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>
using namespace std;

#define MODULENAME "uvtQEMCPCorr"
#define INSIDE_TEMP_NUV  11
#define INSIDE_TEMP_FUV 3
#define INSIDE_TEMP_VIS 19
#define OUTSIDE_TEMP_NUV 17
#define OUTSIDE_TEMP_FUV 9
#define OUTSIDE_TEMP_VIS  25

class uvtQEMCPCorr{
    private:
        char modulename[NAMESIZE];
        char indir[PIL_LINESIZE];             //input data directory/directory containing signal framelist and exposure framelistfiles
        char tempvsfilterfile[PIL_LINESIZE];
         char caldbDir[PIL_LINESIZE]; 
        char lbtfile[PIL_LINESIZE];                      //for finding temperature
        double temperature;                       //will  be taken from LBT file
        //char eventfile[PIL_LINESIZE];
        char outdir[PIL_LINESIZE];             //directory where output directory needs to be created
        char moduleoutdir[PIL_LINESIZE];
        char signalframedir[PIL_LINESIZE];
           char infofile_in[PIL_LINESIZE];    //full path
        char infofile_out[PIL_LINESIZE];
        char expframedir[PIL_LINESIZE];
        char sigframelistfile[PIL_LINESIZE];           //for output
        char expframelistfile[PIL_LINESIZE];
        char inframelistfile[PIL_LINESIZE];
        char inexpframelistfile[PIL_LINESIZE];
      //  char obsid[];
        double *temperaturearr;
        DataInfo datainfo;
        vector<string> key_record;
        
        caldb_Handler caldb_handler;
       
         int xsize,ysize;
        char nameprefix[FLEN_FILENAME]; 
        char filter[FLEN_FILENAME];
        double **qemcpdata;               //buffer to store qemcp data values for each filter in a row wrt to temperature
        int numfilter;
        int nCalDBTempValues;
        int clobber;
        int history;
        char eventfile[PIL_LINESIZE];   //name of input event file
        char imgfile[PIL_LINESIZE];        //name 
        char mode[10];            //pil mode
        int nframes; double factor ;
        char qeFile[PIL_LINESIZE];
        char **sigframelist,**expframelist;   //for input , used for IM only
         long nrows_lbt ;
        double *time_lbt;
          float *insideTemp;
          float *outsideTemp;
           float *qe_mg_factor;
        /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
        int getHistory(vector<string> &vhistory);
        int qemcpCorrectionPC();
        int qemcpCorrectionIM();
        int readQEMCPFile();               //reads temperature vs filter file
        int getTemp();                   //reads lbt file to get temperature for data for QE MCP gain
        
     public:
          float *temp,*f0,*f1,*f2,*f3,*f4,*f5,*f6,*f7;
        uvtQEMCPCorr();
        ~uvtQEMCPCorr();
        int read(int argc,char **argv);
        int read(char *indir,char *caldir,char *lbtfile,char*outdir,
                                int clobber,int history);
        void display();
        int uvtQEMCPCorrectionProcess();
       int readQEMCPFileforCmnArray (char* filename,char *lbtfilename,int &caldb_counts,double &temperature,char *detector, vector<float> &temp,vector<float> &f0,vector<float> &f1,vector<float> &f2,vector<float> &f3,vector<float> &f4,vector<float> &f5,vector<float> &f6);
       const char *getModuleOutdir() const { return moduleoutdir; } 
        
    };
    
 #endif	/* UVTQEMCPCORRECTION_H */

