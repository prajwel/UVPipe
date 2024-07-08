/* 
 * File:   uvtCentroidCorr.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on October 28, 2013, 11:10 AM
 */

#ifndef UVTCENTROIDCORRECTION_H
#define	UVTCENTROIDCORRECTION_H

//#define  DARKFRAME_START  "/raid-data1/data/astrowrk/uvit/Dhruv/DarkFrame_start.fits"
//#define  DARKFRAME_END  "/raid-data1/data/astrowrk/uvit/Dhruv/DarkFrame_end.fits"

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include  <caldb_Handler.h>

using namespace std;

#define MODULENAME "uvtCentroidCorrection"
#define DARKFRAME_SIZE 512
class uvtCentroidCorr{
    private:
        char modulename[NAMESIZE];
        char inputdatadir[PIL_LINESIZE];
        char caldbDir[PIL_LINESIZE];
        //char eventfile[PIL_LINESIZE];
        char outdir[PIL_LINESIZE]; 
        char dataIngestDir[PIL_LINESIZE]; 
       // float  threshold_multph;             //threshold for multiple photon event in case of PC mode
        int clobber;
        int history;
        char mode[10];            //pil mode
vector<string> key_record;       
        char moduleoutdir[PIL_LINESIZE];   //full path of module output directory
       // char sigframedir[PIL_LINESIZE];     //only signal frame directory name , not full path
        //char expframedir[PIL_LINESIZE];  //only exposure frame directory name
        char centroidEAfile[PIL_LINESIZE];
      
        char eventfile[PIL_LINESIZE];   //name of input event file
        char imgfile[PIL_LINESIZE];        //name 
        char infofile_in[PIL_LINESIZE];    //full path
        char infofile_out[PIL_LINESIZE];
        //float *badpixdata;
        int nframes;
        
      //  char **sigframelist;   //for input , used for IM only
        int xsize,ysize;
        int windowsize;
        char nameprefix[FLEN_FILENAME]; 
        char darkdir[FLEN_FILENAME]; 
           char centroid_algo[FLEN_FILENAME]; 
        DataInfo datainfo;
        caldb_Handler caldb_handler;
         float *darkFramestart_data;//[DARKFRAME_SIZE*DARKFRAME_SIZE];
         float *darkFrameend_data;//[DARKFRAME_SIZE*DARKFRAME_SIZE];
            int EA;
//            int centroid_algo;
         double *cent_x,*cent_y,*x_corr,*y_corr;
      char  dstartpath[FLEN_FILENAME];
      char  dendpath[FLEN_FILENAME];
        /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
        int getHistory(vector<string> &vhistory);
        
      //  int filterBadpixPC();
      //  int filterBadpixIM();
        int readcentroidEAFile();               //reads bad pixel file and sets data in badpixdata;
        int readDarkFrame(char * path,float * Array);
     
public:
        uvtCentroidCorr();
        ~uvtCentroidCorr();
        int read(int argc,char **argv);
        int read (char *inputdatadir , char *caldbDir , char *inputdataIngestDir,char *outdir , int clobber , int history);
        void display();
        int uvtCentroidCorrProcess();
        int takeDarkinfo ();
        const char *getModuleOutdir() const { return moduleoutdir; } 
        
    };


#endif	/* UVTCENTROIDCORRECTION_H */

