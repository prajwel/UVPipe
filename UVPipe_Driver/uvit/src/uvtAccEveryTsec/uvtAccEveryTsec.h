/* 
 * File:   uvtAccEveryTsec
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */


#ifndef UVTACCEVERYTSEC_H
#define	UVTACCEVERYTSEC_H

#include<uvtUtils.h>
#include<DataInfo.h>
#include<fitsio.h>
#include<pil.h>
#include<map>

#define MODULENAME "uvtAccEveryTsec"

class uvtAccEveryTsec {
private:
    char modulename[NAMESIZE];

    char inputdatadir[PIL_LINESIZE];
    char outdir[PIL_LINESIZE];
    int numOfFramesToAcc; //number of frames to accumulate
    int clobber;
    int history;
    char mode[5]; //pil mode

    char moduleoutdir[PIL_LINESIZE];
    int nframes; //for input number of frames
    int nframes_out; //for number of output frames;
    int xsize, ysize;
    bool flag_Acc;
    char sigframedir[NAMESIZE];
    char expframedir[NAMESIZE];
    char nameprefix[NAMESIZE];
    char infofile_in[PIL_LINESIZE];
    char infofile_out[PIL_LINESIZE];
    char frameTimeFile[PIL_LINESIZE];
    char **sigframelist;
    char **expframelist;
    vector<string> key_records;

    map<unsigned short, double> FrameTimeMap; //map to store frame times for accumulated frames

    DataInfo datainfo;

    /**
     * Function to accumulate frames for exposure and signal frames 
     * @param dir - Name of directory containing frames 
     * @param namelist - Name list of frames
     * @param frameIdentifier - type of frame "sig" or "exp"
     * @return  - returns 0 on success, non-zero on failure
     */
    int accumulateFrames(char *dir, char **namelist, char *frameIdentifier); //dir will either exposure or signal frame dir



    /**
     * Function to generate history for the module, to be written to the output file
     * @param vhistory
     * @return 
     */
    int getHistory(vector<string> &vhistory);

public:
    uvtAccEveryTsec();
    ~uvtAccEveryTsec();
    void display();
    int read(int argc, char **argv);
    int read(char *inputdatadir, char *outdir, int Nacc, int clobber, int history);
    int uvtAccEveryTsecProcess();

    const char *getModuleOutdir() const {
        return moduleoutdir;
    }

};



#endif	/* UVTACCEVERYTSEC_H */

