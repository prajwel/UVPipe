/* 
 * File:   uvtUnitConversion.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#ifndef UVTUNITCONVERSION_H
#define	UVTUNITCONVERSION_H

 #include "fitsio.h"
 #include "uvtUtils.h"
 #include<pil.h>
#include<DataInfo.h>

using namespace std;

#define MODULENAME "uvtUnitConversion"
#define MODULEDIR "UnitConversion"
#define darkSize 512
#define ENPCOLNAME "EFFECTIVE_NUM_PHOTONS"

class uvtUnitConversion{
    private:
        char modulename[NAMESIZE];
        char dataIngestdir[NAMESIZE];
        char inputdatadir[PIL_LINESIZE];
        char dstartpath[PIL_LINESIZE];
        char dendpath[PIL_LINESIZE];
        char outdir[PIL_LINESIZE];   
         int clobber;
        int history;
        char mode[10];
        int  darkFrame_Flag;
        double t_darkframeend;
        double t_darkframestart;
        char sigframedir[PIL_LINESIZE];        //name of signal frame directory
        char moduleoutdir[PIL_LINESIZE];   //for output
        char eventfile[PIL_LINESIZE];
        char imgfile[PIL_LINESIZE];
        char infofile_in[PIL_LINESIZE];
        char infofile_out[PIL_LINESIZE];
        char nameprefix[FLEN_VALUE];
        char darkdir[FLEN_VALUE];
        
        long nframes;     //used in case of IM
        char **sigframelist;   //for input , used for IM only
        char **outframelist;   
        int xsize,ysize;
        float *darkCompute_array;
        float *darksubtr_array;
        DataInfo datainfo;
       double t_curr;
        float darkFrameend_data[darkSize*darkSize];
        float darkFramestart_data[darkSize*darkSize];
        vector<string> key_record;
    
        /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
        int getHistory(vector<string> &vhistory);
        
        int unitConvertPC();
        int unitConvertIM();
        
     public:
        uvtUnitConversion();
        int read(int argc,char **argv);
        int read(char *inputdatadir,char *outdir,int clobber,int history,int  dark_frame_flag,char *dstart_path,char *dend_path);
        int read_pc_param(char *inputdatadir , char *outdir , int clobber , int history);
        void display();
        int uvtUnitConversionProcess();
      double  readDarkFrame(char * path, float*Array);
      int  darkFrameComputation(float *Array);
     void  darkFrameSubtraction(float *Array, float *frame_data);
     int   takeDarkinfo();
     const char *getModuleOutdir() const { return moduleoutdir; } 
    };

    
    
//struct uvtUnitConversion_param{
//    //char sTimeFrame[FLEN_FILENAME];
//    //char eTimeFrame[FLEN_FILENAME];
//    char iFile[FLEN_FILENAME];
//    char darkFrame[FLEN_FILENAME];
//    char oFile[FLEN_FILENAME];
//    double time;
//    int TimeFlag;
//    double tOneTime;
//    int _switch;
//    int clobber;
//};
//
//void readParam(struct uvtUnitConversion_param *UUC_par, int argc, char **argv);
//void displayParam(struct uvtUnitConversion_param *FF_par);
//
//
//unsigned short **allocate_Memory_short(long height,long width);
//void printError(int);


#endif	/* UVTUNITCONVERSION_H */

