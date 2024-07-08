/* 
 * File:   uvtCentroidBias.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on October 24, 2013, 4:46 PM
 */

#ifndef UVTCENTROIDBIAS_H
#define	UVTCENTROIDBIAS_H


#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include  <caldb_Handler.h>

using namespace std;

#define MODULENAME "uvtCentroidBias"

class uvtCentroidBias{
    private:
        char modulename[NAMESIZE];
        char inputdatadir[PIL_LINESIZE];
        char caldbDir[PIL_LINESIZE];
        //char eventfile[PIL_LINESIZE];
        char outdir[PIL_LINESIZE];   
       // float  threshold_multph;             //threshold for multiple photon event in case of PC mode
        int clobber;
        int history;
        char mode[10];            //pil mode
        int windowsize;
        char centroid_algo[FLEN_FILENAME];
        char moduleoutdir[PIL_LINESIZE];   //full path of module output directory
       // char sigframedir[PIL_LINESIZE];     //only signal frame directory name , not full path
        //char expframedir[PIL_LINESIZE];  //only exposure frame directory name
        char centroidbiasfile[PIL_LINESIZE];
      
        char eventfile[PIL_LINESIZE];   //name of input event file
        char imgfile[PIL_LINESIZE];        //name 
        char infofile_in[PIL_LINESIZE];    //full path
        char infofile_out[PIL_LINESIZE];
        //float *badpixdata;
        int nframes;
        
      //  char **sigframelist;   //for input , used for IM only
        int xsize,ysize;
        char nameprefix[FLEN_FILENAME]; 
        DataInfo datainfo;
        caldb_Handler caldb_handler;
         int biasRows;
         double *fraction_bias,*x_corr,*y_corr;
        /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
         vector<string> key_record;         
        int getHistory(vector<string> &vhistory);
        
      //  int filterBadpixPC();
      //  int filterBadpixIM();
        int readcentroidbiasFile();               //reads bad pixel file and sets data in badpixdata;
        
     public:
        uvtCentroidBias();
        ~uvtCentroidBias();
        int read(int argc,char **argv);
        int read(char *inputdatadir,char *caldbDir,char *outdir,int clobber,int history);
        void display();
        int uvtCentroidBiasProcess();
        
        const char *getModuleOutdir() const { return moduleoutdir; } 
        
    };
#endif	/* UVTCENTROIDBIAS_H */

