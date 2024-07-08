/* 
 * File:   uvtCosmicRayCorr.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#ifndef UVTCOSMICRAYCORR_H
#define	UVTCOSMICRAYCORR_H
#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
using namespace std;

#define MODULENAME "uvtCosmicRayCorr"
#define IMG_DIM 512


class uvtCosmicRayCorr{
   private:
        char modulename[NAMESIZE];
        
        char inputdatadir[PIL_LINESIZE];
     
        char outdir[PIL_LINESIZE];   
        int clobber;
        int history;
        char mode[10];            //pil mode
       // int padding_dimension;
        
        char moduleoutdir[PIL_LINESIZE];   //full path of module output directory
        char sigframedir[PIL_LINESIZE];     //only signal frame directory name , not full path
        char expframedir[PIL_LINESIZE];  //only exposure frame directory name
        double Threshold;
        int nCompareFrames ;
        vector<string> key_records;
         vector<string> vhistorystr ;
        
        char eventfile[PIL_LINESIZE];
        char infofile_in[PIL_LINESIZE];    //full path
        char infofile_out[PIL_LINESIZE];
     
        int nframes;
        
        char **sigframelist;   //for input , used for IM only
          char **expoframelist;   //for input , used for IM only
        int xsize,ysize;
        char nameprefix[FLEN_FILENAME]; 
        DataInfo datainfo;
        
        /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
        int getHistory(vector<string> &vhistory);
        
        int cosmicRayCorrPC();
        int cosmicRayCorrIM();
      //  int readFlatFieldFile();               //reads bad pixel file and sets data in badpixdata;
        
     public:
        uvtCosmicRayCorr();
        ~uvtCosmicRayCorr();
        int read(int argc,char **argv);
  int  read(char *input_datadir,float threshold,int compare_frames,char *out_dir,int clobber,int history);
//  int read_pc(char *input_datadir , float threshold_pc , int compare_frames , char *out_dir , int clobber , int history);
        void display();
        int uvtCosmicRayCorrProcess();
          const char *getModuleOutdir()  { return moduleoutdir; }
  
    };




#endif	/* UVTCOSMICRAYCORR_H */

