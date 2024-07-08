
/* 
 * File:   DataIngest.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */



#ifndef DATAINGEST_H
#define DATAINGEST_H

#include <fitsio.h>
#include <uvtUtils.h>
#include "DataInfo.h"

#include <map>


using namespace std;

class DataIngest{
    
     
   char modulename[NAMESIZE];                   //Variable to store name of module with version number "DataIngest_1.0"
    char inFile[FLEN_FILENAME];                   //Input level-1 science data file
    char inFile_GTI[FLEN_FILENAME];           //Input level-1 GTI file
    char inFile_TCT[FLEN_FILENAME];          //Input level-1 TCT file
    char caldbindir[FLEN_FILENAME];    
    char inFile_LBT[FLEN_FILENAME];           //Input level-1 LBT file
    char outdir[FLEN_FILENAME];              //output directory path
    char inGtifile[FLEN_FILENAME];
    int dropFrame;                                //if yes, frame will be dropped else packet will be dropped where CRC fails
    int history;                                        //if yes, HITSORY will be added to output files
    int clobber;                                       //If yes, output directory will be overwritten if already exists
    char mode[10];                                //PIL mode           
    int parity_flag;                                  //flag to check parity of event data
    long nrow_l1;                                   //Number of rows in level-1 science data file
    double errorRate;             
    char filePrefix[FLEN_FILENAME];           //Variable to store filename prefix containing satellite id, observation id
    char outFile[FLEN_FILENAME];                         //data ingest output file
    char timeCorrectionFile[FLEN_FILENAME];             //data ingest output file containing original and corrected time
    char crcFailedList[FLEN_FILENAME];                  //Output CRC log file
    char missingFrameList[FLEN_FILENAME];               //Output missing framelist file
   // char imFramesDir[FLEN_FILENAME];                    //path for IM mode frames
    //char eventFile[FLEN_FILENAME];                      //Event File in case of PC mode
    char *eventFile;
    char outpath[FLEN_FILENAME];                        //Output Path for Data Ingest output
     char darkFramesoutPath[FLEN_FILENAME];                        //Output Path for Data Ingest output
    char IMframespath[FLEN_FILENAME];             //Output path for IM Frames
    char Darkframes_path[FLEN_FILENAME];
  // char infoFile[FLEN_FILENAME];                       //Output information file path
  char tempinfo[FLEN_FILENAME];
    char dark_infoFile[FLEN_FILENAME];   //file will contain IM frame path list
 //  char PCImageFile[PIL_LINESIZE];                               //image file for PC Mode
    char *infoFilename,*PCImageFile;
    char darkFramesPath[FLEN_FILENAME];                    
   vector<unsigned short> frm_count;
   vector<unsigned int> frm_time_int;
   vector<double>  frm_time_double;
   vector<unsigned char>  GTI_vect;
   char  param_all_or_cust[FLEN_FILENAME];
   char dark_outputDir[FLEN_FILENAME];
    int all_Or_custom,valid_gtiflag,gti_flag;
   long nframes;                                       //number of frames in case of IM
   vector<float> cpu_VIS_temp,cpu_NUV_temp,cpu_FUV_temp;
    vector<string> framelistvector,dark_framelist;       //contains framename list to be added to info file
                                                    //will be populated  in createFrames function                                                            
    DataInfo datainfo;                 //for writing data  info file
     vector<long> list ;
    map<unsigned short, double>  FrameTimeMap;
      map<int, string>  gti_map;
     int track_frmNo;
    long nrow_lbt;
    
    int setFilePaths();                              //function to set filepaths
    int doUTCCorrection();                    //function takes dataingest output file and TCT file and does UTC correction to it
    void doCRCCorrection();
   
    int createFrames();
    int readDarkFrames();
    int createVISFrames();
    int getEvents();
    int createInfoFile();
   
    int create_darkframe_InfoFile();
   int doGTIfiltering (vector<long> &rowno);
   int genL2gti();   //function to generate level2 gti 
    
    template<class T> 
    unsigned char getParity(T value);
  
     vector<string> vhistorystr ;
     vector<string> key_record;
//     float *imageframe_begin,*imageframe_end;
//     
public:
  
   // float findMedianValue(float *input, int noOfValues);
//    void checkAcceptCriatria (int bitPosm,vector<GTIParams> &param_vect);
    float findMedianValue(vector<float> &input, int noOfValues);
    DataIngest();
    ~DataIngest();
    int DataIngestProcess();
    int read(int argc, char **argv);
    int read (char* scidatafile , char *caldb_dir,char* tctfile , char* gtifile , char* lbtfile , char *darkDir , int gtiflag,int validbit,int allOrcust,char* outdir , int dropframe , int parityflag , int clobber , int history);
    void display();
    int getHistory (vector<string> &vhistory);
    const char *getModuleOutdir() const { return outpath; }
    int decodeDarkFrames ( char * filename, vector<float> &Arr,double &timeMedian);
    int readGTIparameters();
   
};

template<class T>
unsigned char DataIngest:: getParity(T value){
    int nbits=sizeof(T)*8;
    unsigned char parity=0;
    unsigned char bit;
    //cout<<endl<<"Number : "<<value<<endl;
    for(int i=0;i<nbits;i++){
        bit=value & 0x1;
        parity=parity^bit;
        value=value>>1;
        //cout<<endl<<(int)bit<<"-----"<<(int)parity;
        //cout<<(int)bit;
    }
    //cout<<endl;
    return parity;
}
template <class T>
bool is_in_range(T val, string minValue, string maxValue) {
    bool x = true;
    double minVal;
    double maxVal;
    if (minValue == "-" && maxValue != "-") {
        maxVal = atof((char*)maxValue.c_str());
        x = (val <= (T) maxVal) ? true : false;
    } else if (minValue != "-" && maxValue == "-") {
        minVal = atof((char*)minValue.c_str());
        x = (val >= (T) minVal) ? true : false;
    } else if (minValue != "-" && maxValue != "-") {
        maxVal = atof((char*)maxValue.c_str());
        minVal = atof((char*)minValue.c_str());
        x = (val >= (T) minVal && val <= (T) maxVal) ? true : false;
    } else {
        x = true;
    }
    return x;
}



struct GTIParams{
    string name;
    unsigned  short bitpos;
   string  min;
   string   max;
   string incl_Or_Excl;
};


#endif	/* DATAINGEST_H */

