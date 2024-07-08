/* 
 * File:   main.cpp
 * Authors::Dhruv,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#ifndef UVTFLATFIELDCORRECTION_H
#define	UVTFLATFIELDCORRECTION_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>


using namespace std;

#define MODULENAME "uvtFlatFieldCorr"
#define paddingsize 600
//#define IMG_DIM_PC  512

class uvtFlatFieldCorr{
private:
    char modulename[NAMESIZE];
    char inputdatadir[PIL_LINESIZE];
    char outdir[PIL_LINESIZE];
    int clobber;
    int history;
    char mode[10]; //pil mode
    char calDir[FLEN_FILENAME];
    char flatfieldfile[FLEN_FILENAME];
    char moduleoutdir[FLEN_FILENAME]; //full path of module output directory
    char sigframedir[PIL_LINESIZE]; //only signal frame directory name , not full path
    char expframedir[PIL_LINESIZE]; //only exposure frame directory name
    char eventfile[PIL_LINESIZE];
    char imgfile[PIL_LINESIZE];
    char infofile_in[PIL_LINESIZE]; //full path
    char infofile_out[PIL_LINESIZE];
    float *flatfielddata;
    int nframes;
    char **sigframelist; //for input , used for IM only
    char **expoframelist; //for input , used for IM only
    int xsize, ysize;
    char nameprefix[FLEN_FILENAME];
    int flatfield_dim1,flatfield_dim2;                 // xsize and ysize for CALDB flatfield image size  ;
    DataInfo datainfo;
    caldb_Handler caldb_handler;
    vector<string> key_records;
    /**
     * Function to generate history for the module, to be written to the output file
     * @param vhistory
     * @return 
     */
    int getHistory(vector<string> &vhistory);
    float *finalFlatFielddata ;
    int flatFieldPC();
    int flatFieldIM();
    int readFlatFieldFile(); //reads bad pixel file and sets data in badpixdata;
  //  int Applypadding(float *inputArray, float *outputArray);
public:
    uvtFlatFieldCorr();
    ~uvtFlatFieldCorr();
    int read(int argc, char **argv);
    int read(char *inputdatadir, char *caldb_dir, char *out_dir, int clobber, int history);
    void display();
    int uvtFlatFieldCorrProcess();
    const char *getModuleOutdir() {
        return moduleoutdir;
    }
};



#endif	/* UVTFLATFIELDCORRECTION_H */

