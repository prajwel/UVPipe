/* 
 * File:   uvtRefFrameCal.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#ifndef UVTREFFRAMECAL_H
#define   UVTREFFRAMECAL_H
#include<pil.h>
#include<fitsio.h>
#include<uvtUtils.h>
#include<DataInfo.h>

using namespace std;


class uvtRefFrameCal
{
   private:
       
       char modulename[NAMESIZE];
       
        char inputdatadir[PIL_LINESIZE];
        char outdir[PIL_LINESIZE];   
        int framesToDiscard;                         
        int clobber;
        int history;
        int  averageFactor;
   
        char mode[10];       
          int   number;
        char moduleoutdir[NAMESIZE];
         char infofile_in[PIL_LINESIZE];
         char infofile_out[PIL_LINESIZE];
         DataInfo datainfo;
         int xsize,ysize;
         char nameprefix[FLEN_VALUE];
        char sigframedir[PIL_LINESIZE];        //name of signal frame directory
        char expoframedir[PIL_LINESIZE];
      int nframes;     //used in case of IM
         char **centroidframelist;   //for input , used for IM only
         // char **starframelist;   //for input , used for IM only
        
      vector<string> key_records;
         char centroidDir[PIL_LINESIZE];
                  
        int getHistory(vector<string> &vhistory);
  //        int refFrameCalPC();
        int refFrameCal();
            
     public:
        uvtRefFrameCal();
        ~uvtRefFrameCal();
          
  
        int read(int argc,char **argv);
    int  read(char* inputdatadir, char* outdir, int frames_toDiscard, int avg_Factor, int clobber, int history) ;
    void display();
        int uvtRefFrameCalProcess();
        int starPairing(char *inputfile_first,char *inputfile_second,float * xfirst,float * yfirst,float * xsecond,float * ysecond,int neighbour);
      
        const char *getModuleOutdir() const { return moduleoutdir; } 
     
         
    };
    
    #endif   /*UVTREFFRAMECAL_H*/
