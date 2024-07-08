/* 
 * File:   uvtRelAspCal.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on March 10, 2014, 11:41 AM
 */
#ifndef UVTRELASPCAL_H
#define	UVTRELASPCAL_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>

#define MODULENAME "uvtRelAspCal"
#define caldbfinalsize 600 
//#define DATA_ARRAYSIZE 8

#define pi 4*atan(1.0)
class uvtRelAspCal
{
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
     int FreqDomainFilter_Flag;
        caldb_Handler caldb_handler;
         int orderpitch,orderyaw,orderroll,fittingflag;
        DataInfo datainfo;
        int xsize,ysize;
        int caldbsize;
        double freqvalue;
          long  caldbdim;
          char nameprefix[FLEN_FILENAME]; 
        char moduleoutdir[PIL_LINESIZE];
          char eventfile[PIL_LINESIZE]; 
    
     double *time;
                 
                              char jitterfile[PIL_LINESIZE];
                char driftfile[PIL_LINESIZE];
                char thermalfile[PIL_LINESIZE];
                 char driftdir[PIL_LINESIZE];
                  char jitterdir[PIL_LINESIZE];
                  char thermaldir[PIL_LINESIZE];
                        
                   
                float *dx,*dy;
                
             //   int *X,*Y;
                int getHistory(vector<string> &vhistory);
                vector<string> key_record;
      
        int computeRelAsp();
     char *RASfile;
     public:
        uvtRelAspCal();
        ~uvtRelAspCal();
         
        int read(int argc,char **argv);
        int read(char *drift_dir, char *jitter_dir,char *thermal_dir, char *outdir, int clobber, int history) ;
        void display();
        int uvtRelAspCalProcess();
        int readFitsFile (char* filename, double *timedata , double* Xdata, double* Ydata, double* thetadata,vector<double> &timeterm,vector<double> &xterm,vector<double> &yterm,vector<double> &thetaterm);
        const char *getModuleOutdir() const { return moduleoutdir; } 
      //  void   setRelAspFilename(char *relfilename)
      // {
    //    this->RASfile=relfilename;
       // cout<<"B$ "<<RASfile<<endl;
       // }
    char * getRelAspFilename() 
    const
    {
        cout<<"RAS fie is "<<RASfile<<endl;
        return RASfile;
    }
            
    };


#endif	/* UVTRELASPCAL_H */

