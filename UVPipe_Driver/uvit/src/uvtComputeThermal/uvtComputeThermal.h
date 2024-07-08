/* 
 * File:   uvtComputeThermal.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on February 24, 2014, 11:41 AM
 */

#ifndef UVTCOMPUTETHERMAL_H
#define	UVTCOMPUTETHERMAL_H
#include<fitsio.h>
#include<uvtUtils.h>
#include<pil.h>
#include<DataInfo.h>
#include<caldb_Handler.h>

using namespace std;

#define MODULENAME "uvtComputeThermal"              //this modulename is used for creation of output directory

class uvtComputeThermal {
private:
    char modulename[NAMESIZE];

    char caldb_path[PIL_LINESIZE];
    char LBTpath[PIL_LINESIZE];
    char indir[PIL_LINESIZE];
 
    char outdir[PIL_LINESIZE];
    int clobber;
    int history;
    char mode[10]; //pil mode
    char moduleoutdir[PIL_LINESIZE]; //full path of module output directory
     int nCalDBTempValues;
    float *temp_One_caldb, *temp_Two_caldb, *temp_Three_caldb, *temp_Four_caldb, *yaw_corr_caldb, *pitch_corr_caldb;
      float scaleXCH1, scaleYCH1,scaleXVIS,scaleYVIS, Shift_X, Shift_Y, angle_NUVtoVIS;
   char nameprefix[FLEN_FILENAME];
    DataInfo datainfo;
    char filter[FLEN_FILENAME];
    char caldb_orig_path[FLEN_FILENAME];
    char InframeDir[FLEN_FILENAME];
    char frameIntegrationDir[FLEN_FILENAME];
 char infofile_in[PIL_LINESIZE];    //full path
 char dataIngestinfofile_in[PIL_LINESIZE];    //full path
 char frameIntegration_infofile_in[PIL_LINESIZE];    //full path
    caldb_Handler caldb_handler;
    char thermalFile[PIL_LINESIZE];
    int nframes;
    char sigdir[FLEN_FILENAME];
    unsigned short *MB_top_temp,*TT_top_temp,*TT_bot_temp,*BOT_ring_temp;
    /**
     * Function to generate history for the module, to be written to the output file
     * @param vhistory
     * @return 
     */
    int getHistory(vector<string> &vhistory);


    int thermalCalc();
   
   
    int readThermalFile();
    int getTempFromLBT();
    int cnt_Fortemp1_repeat, cnt_Fortemp2_repeat, cnt_Fortemp3_repeat, cnt_Fortemp4_repeat;
    long num_record_LBT;
    //float *temp_One_HK, *temp_Two_HK, *temp_Three_HK, *temp_Four_HK;
    

public:
    vector<string> key_records;
    uvtComputeThermal();
    ~uvtComputeThermal();
    int read(int argc, char **argv);
    int read (char *indir , char *dataIngestindir,char*caldbdir , char* hkdir , char *outdir , int clobber , int history);
    void display();
    int uvtThermalCalcProcess();
     int writeThermalcorrTotable (vector<double> &timeData, vector<double> &yawcorr , vector<double> &pitchcorr , vector<double> &rollcorr);
    const char *getModuleOutdir() {
        return moduleoutdir;
    }
};







#endif	/* UVTCOMPUTETHERMAL_H */

