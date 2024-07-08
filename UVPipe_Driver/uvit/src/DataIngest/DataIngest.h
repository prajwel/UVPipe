
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
    char inFile_MKF[FLEN_FILENAME];
    char attFilename[FLEN_FILENAME];
    char caldbindir[FLEN_FILENAME];    
    char inFile_LBT[FLEN_FILENAME];           //Input level-1 LBT file
    char outdir[FLEN_FILENAME];              //output directory path
    char inGtifile[FLEN_FILENAME];
    int dropFrame;   				//if yes, frame will be dropped else packet will be dropped where CRC fails
    int crcFlag;                             
    int history;                                        //if yes, HITSORY will be added to output files
    int clobber;                                       //If yes, output directory will be overwritten if already exists
    char mode[10];                                //PIL mode           
    int parity_flag;                                  //flag to check parity of event data
    long nrow_l1,nrows_l1_Original,totalpackateUsed;                                     //Number of rows in level-1 science data file
    string clock_master;                        //Clock master channel
    double errorRate;             
    char filePrefix[FLEN_FILENAME];           //Variable to store filename prefix containing satellite id, observation id
    char outFile[FLEN_FILENAME];                         //data ingest output file
    char timeCorrectionFile[FLEN_FILENAME];             //data ingest output file containing original and corrected time
    char crcFailedList[FLEN_FILENAME];                  //Output CRC log file
    char missingFrameList[FLEN_FILENAME];               //Output missing framelist file
   // char imFramesDir[FLEN_FILENAME];                    //path for IM mode frames
    char eventFile[FLEN_FILENAME];                      //Event File in case of PC mode
    char outpath[FLEN_FILENAME];                        //Output Path for Data Ingest output
     char darkFramesoutPath[FLEN_FILENAME];                        //Output Path for Data Ingest output
    char IMframespath[FLEN_FILENAME];             //Output path for IM Frames
    char IMDiscardedframespath[FLEN_FILENAME];
    char Darkframes_path[FLEN_FILENAME];
    char infoFile[FLEN_FILENAME];                       //Output information file path
    char dark_infoFile[FLEN_FILENAME];   //file will contain IM frame path list
    char PCImageFile[PIL_LINESIZE];                               //image file for PC Mode
    char darkFramesPath[FLEN_FILENAME];  
    bool flag_tct_Continues;
 vector<double> frmDataTime_backup;
   vector<unsigned short> frm_count;
long nFrames_InputTomakeEvents,nframes_outAsEvents,nframes_zeroCentroid;
  // vector<unsigned int> frm_time_int;
   vector<long> frm_time_int;
   vector<double>  frm_time_double;
   vector<unsigned char>  GTI_vect;
   char  param_all_or_cust[FLEN_FILENAME];
   char dark_outputDir[FLEN_FILENAME];
   int no_frms_stacked;
    int all_Or_custom,valid_gtiflag,gti_flag,UTC_Flag;
   long nframes;                                       //number of frames in case of IM
   long nstartIndex_level1file;
double Bzero_MJD,Bscale_MJD;
   vector<float> cpu_VIS_temp,cpu_NUV_temp,cpu_FUV_temp,FUV_PRI_PLATE,FUV_MB_TOP_LOC,FUV_TT_MIDDLE,FUV_TT_TOP,FUV_TT_BOTTOM,
   FUV_BOT_RING,NUV_PRI_PLATE,NUV_MB_TOP_LOC,NUV_TT_MIDDLE,NUV_TT_TOP,NUV_TT_BOTTOM,NUV_BOT_RING,SUN_ANGLE,MOON_ANGLE,
   TIME_SINCE_SAA,TEMP_HVU_VIS,TEMP_HVU_FUV,TEMP_HVU_NUV,EU_BASETEMP_TH,MCP_VOLT_VIS,MCP_VOLT_NUV,MCP_VOLT_FUV,CATHOD_VOLT_VIS,
   CATHOD_VOLT_NUV,CATHOD_VOLT_FUV,ANODE_VOLT_VIS,ANODE_VOLT_FUV,ANODE_VOLT_NUV,ACQUIRE_BIT_FW_VIS,ACQUIRE_BIT_FW_NUV,ACQUIRE_BIT_FW_FUV,
   SSM_MOTION,BR_EARTH,TARGET_IN_FOV,SAA_FLAG,RSERR_COUNTF,RSERR_COUNTN,RSERR_COUNTV;
    vector<string> framelistvector,dark_framelist;       //contains framename list to be added to info file
                                                    //will be populated  in createFrames function                                                            
    DataInfo datainfo;                 //for writing data  info file
     vector<long> list ,list_track,listbackup;
long parityfailedEvents;
    map<unsigned short, double>  FrameTimeMap;
    map<int ,vector<float> > GTI_bitPos_vect;
      vector<double> listof_attTime;
long TotalpktsFrmL1,totalPktused,TotalEventsGenerated;
      map<int, string>  gti_map;
     int track_frmNo;
    long nrow_lbt,nrows_mkf;
    vector<long> Rows_bad_Framepackets;
    int Att_flag_val;
    
    int cnt_VISForMasterClock,cnt_NUVForMasterClock,cnt_FUVForMasterClock;
    long cnt_TotForMasterClock;
    double integrationTime_frm_UTC;
    double integrationTime_WithCurrentOption;
    int setFilePaths();                              //function to set filepaths
    int doUTCCorrection();                    //function takes dataingest output file and TCT file and does UTC correction to it
    void doCRCCorrection();
   
    int createFrames();
    int createDiscardedFrames();
    int readDarkFrames();
    int createVISFrames();
    int getEvents();
    int createInfoFile();
   
    int create_darkframe_InfoFile();
   int doGTIfiltering (vector<long> &rowno);
   int genL2gti();   //function to generate level2 gti 
    
    template<class T> 
    unsigned char getParity(T value);
    //static  int  meri(int g,int k);
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
   //int read (char* scidatafile , char *caldb_dir , char* tctfile , char *mkfFile , char* gtifile , char* lbtfile ,char *attfile ,char *darkDir ,int Att_flag, int gtiflag , int validbit , int allOrcust , char* outdir , int dropframe , int parityflag , int clobber , int history);
  int read (char* scidatafile , char *caldb_dir , char* tctfile , char *mkfFile , char* gtifile , char* lbtfile ,char *attfile ,char *darkDir ,int Att_flag, int gtiflag , int validbit , int allOrcust , char* outdir , int dropframe , int parityflag , int UTCFlag,int crcflag,int clobber , int history);
    void display();
    int getHistory (vector<string> &vhistory);
    const char *getModuleOutdir() const { return outpath; }
    int decodeDarkFrames ( char * filename, vector<float> &Arr,double &timeMedian);
    int readGTIparameters();
	 int FindCoeffientFromTCT( char *inFile_TCT,string clock_master, double &BZERO_MJD ,double &BSCALE_MJD );
    
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
bool is_in_range(T val, string minValue, string maxValue,string  inc_exclude) {
    bool x = true;
    double minVal;
    double maxVal;
    if(strcmp((char*)inc_exclude.c_str(),"0")==0 ){
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
         //  if(val==17.5){ LOG(INFO)<<"GEtit"<<" "<<minVal<< " "<<maxVal<<" "<<x<<endl;}
    } 
    }
    else if (strcmp((char*)inc_exclude.c_str(),"1")==0)
    {
          if (minValue == "-" && maxValue != "-") {
        maxVal = atof((char*)maxValue.c_str());
        x = (val >= (T) maxVal) ? true : false;
    } else if (minValue != "-" && maxValue == "-") {
        minVal = atof((char*)minValue.c_str());
        x = (val <= (T) minVal) ? true : false;
    }
          else if (minValue != "-" && maxValue != "-") {
        maxVal = atof((char*)maxValue.c_str());
        minVal = atof((char*)minValue.c_str());
        x = (val < (T) minVal || val > (T) maxVal) ? true : false;
          }
        
    }
    else if (strcmp((char*)inc_exclude.c_str(),"2")==0){//for checking equal values
         minVal = atof((char*)minValue.c_str());
          x = (val == (T) minVal ) ? true : false;
     }
    else {
       
        x = true;
    }
   
    
    return x;   
   
}
//bool is_in_range(T val, string minValue, string maxValue) {
//    bool x = true;
//    double minVal;
//    double maxVal;
//   
//    if (minValue == "-" && maxValue != "-") {
//        maxVal = atof((char*)maxValue.c_str());
//        x = (val <= (T) maxVal) ? true : false;
//    } else if (minValue != "-" && maxValue == "-") {
//        minVal = atof((char*)minValue.c_str());
//        x = (val >= (T) minVal) ? true : false;
//    } else if (minValue != "-" && maxValue != "-") {
//        maxVal = atof((char*)maxValue.c_str());
//        minVal = atof((char*)minValue.c_str());
//        x = (val >= (T) minVal && val <= (T) maxVal) ? true : false;
//    } else {
//        x = true;
//    }
//    return x;
//    
//    
//}



struct GTIParams{
    string name;
    unsigned  short bitpos;
   string  min;
   string   max;
   string incl_Or_Excl;
};


#endif	/* DATAINGEST_H */

