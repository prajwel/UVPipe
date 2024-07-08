/* 
 * File:   DataInfo.h
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#ifndef  DIRECTORY_H
#define	DIRECTORY_H

#include<map>
#include<string>
#include<dirent.h>
#include<cstdlib>
#include<vector>
#include<cstring>
#include<iostream>
#include<uvtUtils.h>
#include<glog/logging.h>

using namespace std;

#define FUV 1
#define NUV 2
#define VIS 3

#define FUV_DIR "uvtF"
#define NUV_DIR "uvtN"
#define VIS_DIR "uvtV"

#define UVIT_DIR "uvit"
#define FILEPATHSEPARATOR "/"
#define DARKDIRN   "DarkN"
#define DARKDIRF   "DarkF"
#define DARKDIRV   "DarkV"
#define SCIENCEDATA_EXT "_level1.fits"
#define ATTITUDE_EXT     ".att"
#define ORBIT_EXT  ".orb"
#define GYRO_EXT   ".gyr"
#define TCT_EXT    ".tct"
#define GTI_EXT     ".gti"
#define BTI_EXT    ".bti"
#define LBT_EXT  ".lbt"
#define HK_EXT      ".hk"
#define MKF_EXT   ".mkf"
#define TILDE "~"

//class to hold the level1 filepaths using level1 data directory and also create level2 directory

typedef struct Directory{
    
    vector<string> sciencedatafile;
    vector<string> gtifile;
    vector<string>  btifile;
    
    vector<string> level2path;       //contains output path for corresponding files
       
    string hkfile;
    string lbtfile;
    string orbitfile;
    string gyrofile;
    string mkffile;
    string attfile;
    string tctfile;
    string darkDirectory;
    string orbit_no;
     //function to set the files using path of level1 directory    channel= FUV, NUV or VIS
    //returns SUCCESS if files found else returns FAILURE
    int setup(string level1dir, string level2dir, int channel);  
    int CopyAuxDir(string level1dir, string level2dir, int channel);
    void setOrbitNo(string orbno);
        
} Directory;



  int get_firstLevelSubDirs (string dirName , vector<string> &subDirs);
  int getFiles(string dirName, vector<string> &filelist);
  int getDirs(string dirName, vector<string> &subDirs);
  int uvitF_filter(const struct dirent *dptr);
  int uvitN_filter(const struct dirent *dptr);
  int uvitV_filter(const struct dirent *dptr);
 int ExpandMergeData(string level1tar,string l2path_temp,vector<string> &individual_Tarnames);
 int readScienceDataFilesPath (string level1dir  , int channel,string orbit_no,string darkDirectory,vector<string> &sciencedatafile);

#endif	/* LEVEL1_H */

