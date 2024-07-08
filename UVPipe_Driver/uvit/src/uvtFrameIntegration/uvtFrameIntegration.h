/* 
 * File:   uvtFrameIntegration.h
 * Author: uvit
 *
 * Created on October 29, 2013, 4:02 PM
 */

#ifndef UVTFRAMEINTEGRATION_H
#define	UVTFRAMEINTEGRATION_H

#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>

#define MODULENAME "uvtFrameIntegration"
#define IMG_DIM_FI 9600
//#define FRAMESIZEFACTOR 16
//typedef struct Star
//{
//    float intensity ;
//    int x , y ;
//};
class uvtFrameIntegration {
private:
    char modulename[NAMESIZE];                      //Name of this module/program
    char inputdatadir[PIL_LINESIZE];                 //Input data directory path
    char outdir[PIL_LINESIZE];                           //Output directory path
    int clobber;                                                        //Overwrite
    int history;                                                         //History write or not
    char mode[10];                                                 //pil mode
    int nframes;                                                       //Number of input frames
    char infofile_in[PIL_LINESIZE];                      //Input Info file path 
    char infofile_out[PIL_LINESIZE];                   //Output info file path
    DataInfo datainfo;                                          //To read and write basic data information
    int xsize, ysize;                                                //Size of frame
    int Nacc, Ndiscard;                                        //Number of frames to discard at beginning, Number of frames to accumulate at a time
    char nameprefix[FLEN_FILENAME];          //prefix for input filenames
    char moduleoutdir[PIL_LINESIZE];            //Name of output module directory which will be created at outdir path
    char eventfile[PIL_LINESIZE];                    //Path for input eventfile
    vector<string> key_records;
     /**
     * Function to generate history for the module, to be written to the output file
     * @param vhistory
     * @return 
     */
      int getHistory(vector<string> &vhistory);
//      vector<Star> star_track_exp,star_track_sig;
public:
    uvtFrameIntegration();
    ~uvtFrameIntegration();
    int read(int argc, char **argv);
    int read(char *inputdata_dir, int frm_discard, int frm_cmpt, char* out_dir, int clobber, int history);
    void display();
    int uvtFrameIntProcess();
    void doCentroiding (vector<float> &X , vector<float> &Y , vector<float> Intensity,int centroidwindow , float *arr , int h , int w);
    const char *getModuleOutdir() const                          //return the name of the output directory created
    { return moduleoutdir;      }
};


#endif	/* UVTFRAMEINTEGRATION_H */

