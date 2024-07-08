/* 
 * File:   uvtPixPadding.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#ifndef UVTPIXPADDING_H
#define	UVTPIXPADDING_H
#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
using namespace std;

#define MODULENAME "uvtPixPadding"

class uvtPixPadding{
    private:
        char modulename[NAMESIZE];
        
        char inputdatadir[PIL_LINESIZE];
     
        char eventfile[PIL_LINESIZE];
        char outdir[PIL_LINESIZE];   
        int clobber;
        int history;
        char mode[10];            //pil mode
        int padding_dimension;
        
        char moduleoutdir[PIL_LINESIZE];   //full path of module output directory
        char sigframedir[PIL_LINESIZE];     //only signal frame directory name , not full path
        char expframedir[PIL_LINESIZE];  //only exposure frame directory name
      
   
        char infofile_in[PIL_LINESIZE];    //full path
        char infofile_out[PIL_LINESIZE];
    
        int nframes;
        vector<string> key_record;
        
        char **sigframelist;   //for input , used for IM only
          char **expoframelist;   //for input , used for IM only
        int xsize,ysize;//dimentions for input and output image
        char nameprefix[FLEN_FILENAME]; 
        DataInfo datainfo;
        
        /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
        int getHistory(vector<string> &vhistory);
        
        int pixPaddingPC();
        int pixPaddingIM();
        
   //     int  Applypadding(float  *inputArray,float *outputArray);
     
     public:
        uvtPixPadding();
        ~uvtPixPadding();
        int read(int argc,char **argv);
        int read(char *inputdatadir,int outputdimension,char *outdir,int clobber,int history);
        void display();
        int uvtPixPaddingProcess();
          int Applypadding_commonArray(float* inputArray, float* outputArray,int padd_dim,int size_x);
        const char *getModuleOutdir()  { return moduleoutdir; }
    };





#endif	/* UVTPIXPADDING_H */

