/* 
 * File:   uvtDetectDistCorr.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#ifndef UVTDETECTDESTCORR_H
#define	UVTDETECTDESTCORR_H
#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>

#define MODULENAME "uvtDetectDistCorr"
#define CALDBFINALSIZE 600 
#define PADDED_FRAMESIZE 600

class uvtDetectDistCorr {
private:
    char modulename[NAMESIZE];
    char inputdatadir[PIL_LINESIZE];
    char caldbDir[PIL_LINESIZE];
    char outdir[PIL_LINESIZE];
    int clobber;
    int history;
    char mode[10]; //pil mode
    int nframes;
    char infofile_in[PIL_LINESIZE]; //full path
    char infofile_out[PIL_LINESIZE];
    caldb_Handler caldb_handler;
    DataInfo datainfo;
    vector<string> key_record;
    int xsize, ysize;
    int caldbsize;
    long caldbdim;
    float x_Distortion[CALDB_DIST_SIZE*CALDB_DIST_SIZE], y_Distortion[CALDB_DIST_SIZE*CALDB_DIST_SIZE];
    char distortionCorrfile[PIL_LINESIZE];
    char nameprefix[FLEN_FILENAME];
    char moduleoutdir[PIL_LINESIZE];
    float caldb_xdist_final[CALDBFINALSIZE*CALDBFINALSIZE];
    float caldb_ydist_final[CALDBFINALSIZE*CALDBFINALSIZE];
    char eventfile[PIL_LINESIZE];
    vector<string> key_records;
    /**
     * Function to generate history for the module, to be written to the output file
     * @param vhistory
     * @return 
     */
    //char starDir[PIL_LINESIZE];
    char centroidDir[PIL_LINESIZE];
    char **centroidframelist; //for input , used for IM only

    //float *dx, *dy;
 
    int getHistory(vector<string> &vhistory);

    int detectDistortionIM();
    int detectDistortionPC();
   //int readDistortionFile(); //reads bad pixel file and sets data in badpixdata;
    //int Applypadding(float *inputArray, float *outputArray);
public:
    uvtDetectDistCorr();
    ~uvtDetectDistCorr();
    int read(int argc, char **argv);
    int read(char *inputdatadir, char *caldbDir, char *outdir, int clobber, int history);
    void display();
    int uvtDetectDistCorrProcess();
 //   int readCaldbDistoFile(float *Xdist, float *Ydist, char *distrotionCorrfile, long &dim);

    const char *getModuleOutdir() const {
        return moduleoutdir;
    }

};

#endif	/* UVTDETECTDESTORTION_H */

