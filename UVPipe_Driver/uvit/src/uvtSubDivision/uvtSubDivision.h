/* 
 * File:   main.cpp
 * Authors:: Dhruv,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#ifndef UVTSUBDIVISION_H
#define	UVTSUBDIVISION_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>

using namespace std;

#define MODULENAME "uvtSubDivision"

class uvtSubDivision{
      private:
        char modulename[NAMESIZE];
        char inputdatadir[PIL_LINESIZE];
     
        //char eventfile[PIL_LINESIZE];
        char outdir[PIL_LINESIZE];   //output Directory For the Fits Files
        int clobber;
        int history;
        char mode[10];            //pil mode
        int subdivision_dimension;//dimention for the output image
        
        char moduleoutdir[PIL_LINESIZE];   //full path of module output directory
        char sigframedir[PIL_LINESIZE];     //only signal frame directory name , not full path
        char expframedir[PIL_LINESIZE];  //only exposure frame directory name
      
        char eventfile[PIL_LINESIZE];//for  PhotonCounting Mode
        char infofile_in[PIL_LINESIZE];    //path  for the input List  fit File
        char infofile_out[PIL_LINESIZE]; //path for the output List fit File
  
        int nframes;//number of Fremes
        
        char **sigframelist;   //for input , used for IM only
          char **expoframelist;   //for input , used for IM only
        int xsize,ysize; //dimention of the input image
        char nameprefix[FLEN_FILENAME]; 
        DataInfo datainfo;
        vector<string> key_records;
        
        /**
        * Function to generate history for the module, to be written to the output file
        * @param vhistory
        * @return 
        */
        int getHistory(vector<string> &vhistory);
         int  applySubDivision(float  *inputArray,float *outputArray);  //For Appling SubDivision
        int subDivisionPC();
        int subDivisionIM();
      //  int readFlatFieldFile();               //reads bad pixel file and sets data in badpixdata;
        
     public:
        uvtSubDivision();
        ~uvtSubDivision();
        int read(int argc,char **argv);
        int read(char *inputdatadir,int outputdimension,char *outdir,int clobber,int history);
        void display();
        int uvtSubDivisionProcess();
  
      
  const char *getModuleOutdir() const { return moduleoutdir; } 
};



#endif	/* UVTSUBDIVISION_H */

