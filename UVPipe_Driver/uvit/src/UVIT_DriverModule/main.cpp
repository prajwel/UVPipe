
/* 
 * File:   main.cpp
 * Author: uvit
 *
 * Created on December 9, 2013, 4:04 PM
 */

#include <cstdlib>

#include "UVIT_DriverModule.h"
#include<glog/logging.h>
#include<csignal>
#include<Directory.h>
#include <fstream>
#include<sstream>
#define EXE_DRIFTCOMPUTE_IM "uvtComputeDrift_Manualmode"

bool compare1(struct Star star1, struct Star star2) {
    return (star1.intensity > star2.intensity);
}

void signalHandler(int signum) {
    char temp;
    cout << "\033[1;31mDo you Want to close  the Execution??[y/n].\n\033[0m"; //Do your handling here
    cin>>temp;
    if (temp == 'y') {
        exit(signum);
    }
}
int writeUsrkeywordsFrmvectDriver(char *filename, vector<string> &strvect);

int main(int argc, char** argv) {
    signal(SIGINT, signalHandler); //Interrupt will be captured here.
    while (1) {

        //pid_t parentPid = getppid();
        //    const char* cmdline=get_process_name_by_pid(parentPid);
        // string cmdline = get_process_name_by_pid(parentPid);
        //        LOG(INFO)<<cmdline;exit(1);
        //        if (strstr(cmdline.c_str(), "-") == NULL) {
        //            LOG(ERROR) << "\033[1;31m***SHELL used is NOT C-SHELL!!,Please change the shell to C-shell (just type 'csh' on terminal ) and try again!!***\033[0m";
        //            return (EXIT_FAILURE);
        //        }
        google::InitGoogleLogging(argv[0]); //Initialization of google logging library
        google::SetStderrLogging(google::INFO); //Setting path to create log files. It reads the environment variable 'GLOG_log_dir'
        UVIT_DriverModule uvitObj;
        int status = uvitObj.readPILParameters(argc, argv);
        if (status) {
            LOG(ERROR) << "Error in reading parameter for Driver module";
            return (EXIT_FAILURE);

        }
        //
////added
string tempstr,level1dirname;
tempstr.assign(uvitObj.level1indir);
string orbitNum,cmdName;

  string tarname_s;
        //tarname = uvitObj.level1indir;
        tarname_s = uvitObj.level1indir;
        string strtemp_s(tarname_s);
        int pos_new = strtemp_s.find_last_of("/");
        string tempFilename_s = "";
        //string cmd;
        string tar_date_data_s = tarname_s.substr(pos_new + 1 + 11, 8);
        string tar_OBS_id_data_s = tarname_s.substr(pos_new + 1 + 19, 21);
        string outputNUV_s, outputFUV_s;
        string orbitnum_s = tarname_s.substr(pos_new + 1 + 41, 5);
        LOG(INFO) << "The orbit num is " << orbitnum_s;

        string path_l1dir_s;
        if (pos_new > 0 && pos_new < strtemp_s.length()) {
            tempFilename_s = strtemp_s.substr(0, pos_new);
            path_l1dir_s = tempFilename_s + "/" + tar_date_data_s + "_" + tar_OBS_id_data_s + "_level1";

        } else {
            path_l1dir_s = tar_date_data_s + "_" + tar_OBS_id_data_s + "_level1";
        }

 if (DirExists((char*) path_l1dir_s.c_str())) {
                            cmdName = "rm -r " + path_l1dir_s;
LOG(INFO)<<"Executing "<<cmdName;
                            system(cmdName.c_str());
                        }
if (uvitObj.zipFlag == FALSE) {

                
            status = extractTars(tempstr, level1dirname, orbitNum); //extract level-1 data  tar file

                if (status) {
                    LOG(INFO) << "Error in extracting tar";
                    return (EXIT_FAILURE);
                }
           } else {
             status = extractZip(tempstr, level1dirname, orbitNum); //extract level-1 data  tar file

           if (status) {
             LOG(INFO) << "Error in extracting tar";
           return (EXIT_FAILURE);
          }
          }

//level1dirname="20161102_G06_157T01_9000000772_level1";
vector<string> AllFilesFromL1;
vector<string> ScienceDatafiles_FrmL1;
getFiles((char*)level1dirname.c_str(),AllFilesFromL1);

for(int i =0;i<AllFilesFromL1.size();i++)
{
if (strstr(AllFilesFromL1[i].c_str(), ".fits") != NULL) {
ScienceDatafiles_FrmL1.push_back(AllFilesFromL1[i]);
//LOG(INFO)<<"GG "<<AllFilesFromL1[i];

}
}
fitsfile *fptr_sci;
status=0;
long tot_Num_Rows;
 unsigned short *frmCount_scienceData;
vector<string> ScienceDatFile_TobeEdited;
vector<long> JumpOccureAtRowNUMTrack,TotRowsFROML1;
vector<long> StartAtRowNUMTrack;
long start_Row_FrmScienceData;
vector<long> Totlpackets;
for(int i=0;i<ScienceDatafiles_FrmL1.size();i++){

fits_open_file(&fptr_sci , ScienceDatafiles_FrmL1[i].c_str () , READONLY , &status) ;
     
   fits_movabs_hdu (fptr_sci , 3 , NULL , &status) ;
 fits_get_num_rows (fptr_sci , &tot_Num_Rows, &status) ;
   
   frmCount_scienceData = new unsigned short[tot_Num_Rows];
     
    fits_read_col (fptr_sci , TUSHORT , 11 , 1 , 1 , tot_Num_Rows , NULL , frmCount_scienceData , NULL , &status) ;
start_Row_FrmScienceData=1;
for(int index=1;index<tot_Num_Rows;index++)
{
if(frmCount_scienceData[index]<frmCount_scienceData[index-1] && frmCount_scienceData[index]==0 && frmCount_scienceData[index-1]==65535)
{
ScienceDatFile_TobeEdited.push_back(ScienceDatafiles_FrmL1[i]);
StartAtRowNUMTrack.push_back(start_Row_FrmScienceData);
JumpOccureAtRowNUMTrack.push_back(index);
start_Row_FrmScienceData=index+2;
Totlpackets.push_back(tot_Num_Rows);
}
if(start_Row_FrmScienceData!=1 && index==tot_Num_Rows-1){
ScienceDatFile_TobeEdited.push_back(ScienceDatafiles_FrmL1[i]);
StartAtRowNUMTrack.push_back(start_Row_FrmScienceData);
JumpOccureAtRowNUMTrack.push_back(index+1);
Totlpackets.push_back(tot_Num_Rows);
}

}
    
    
  fits_close_file(fptr_sci,&status);


delete[] frmCount_scienceData;

}

vector<string> scienceDataToEdited_backup;
scienceDataToEdited_backup=ScienceDatFile_TobeEdited;
//searching Maximum number (i.e N.01,N.02..etc in NUV and FUV)
int pos_s;
string numNUV,numFUV;
int numNUV_int,numFUV_int;
for(int i=0;i<ScienceDatafiles_FrmL1.size();i++)
{
pos_s=ScienceDatafiles_FrmL1[i].find("/uvtN.");
if(pos_s>=0 && pos_s<ScienceDatafiles_FrmL1[i].size()){

numNUV=ScienceDatafiles_FrmL1[i].substr(pos_s+6,2);
LOG(INFO)<<"Number is "<<numNUV;
break;
}

}
for(int i=0;i<ScienceDatafiles_FrmL1.size();i++)
{
pos_s=ScienceDatafiles_FrmL1[i].find("/uvtF.");
if(pos_s>=0 && pos_s<ScienceDatafiles_FrmL1[i].size()){

numFUV=ScienceDatafiles_FrmL1[i].substr(pos_s+6,2);
LOG(INFO)<<"Number is "<<numFUV;
break;
}
}
numNUV_int=atoi(numNUV.c_str());
numFUV_int=atoi(numFUV.c_str());





//till this
vector<string> Copied_location_scienceDataFile;
int start_RowScienceData=0;
string TomakeFileString,cmdTogive;
string nameOfgtisubstr,nameOfbtisubstr;
string intToString,btiFile,gtiFile;
string FullPathNuvNewString,FullPathFuvNewString;
string tempStringNUV,tempStringFUV;
LOG(INFO)<<ScienceDatFile_TobeEdited.size();
for(int i=0;i<ScienceDatFile_TobeEdited.size();i++){
if(strstr(ScienceDatFile_TobeEdited[i].c_str(),"/uvtN/")!=NULL){
LOG(INFO)<<"III "<<ScienceDatFile_TobeEdited[i];

start_RowScienceData=0;
do{
LOG(INFO)<<i<<" "<<ScienceDatFile_TobeEdited[i];
numNUV_int++;
intToString=convertIntToStr(numNUV_int);
if(numNUV_int<10){
tempStringNUV="0"+intToString;
intToString=tempStringNUV;
}
TomakeFileString=level1dirname+"/uvit/"+orbitNum+"/uvtN/uvtN."+intToString+"/";
FullPathNuvNewString=TomakeFileString+basename(ScienceDatFile_TobeEdited[i].c_str());
Copied_location_scienceDataFile.push_back(FullPathNuvNewString);
LOG(INFO)<<TomakeFileString<<" "<<FullPathNuvNewString;
cmdTogive="mkdir -p "+TomakeFileString;
system(cmdTogive.c_str());
cmdTogive="cp  "+ScienceDatFile_TobeEdited[i]+" "+TomakeFileString;
system(cmdTogive.c_str());
gtiFile=ScienceDatFile_TobeEdited[i].replace(ScienceDatFile_TobeEdited[i].size()-4,3,"gti");
nameOfgtisubstr=gtiFile.substr(0,gtiFile.size()-1);
LOG(INFO)<<"GTI file is "<<nameOfgtisubstr;
btiFile=ScienceDatFile_TobeEdited[i].replace(ScienceDatFile_TobeEdited[i].size()-4,3,"bti");
nameOfbtisubstr=btiFile.substr(0,btiFile.size()-1);
LOG(INFO)<<"BTI file is"<<nameOfbtisubstr;
cmdTogive="cp  "+nameOfgtisubstr+" "+TomakeFileString;
system(cmdTogive.c_str());
cmdTogive="cp  "+nameOfbtisubstr+" "+TomakeFileString;
system(cmdTogive.c_str());
//for(int j=start_RowScienceData;j<JumpOccureAtRowNUMTrack[i];j++){


//}
//start_RowScienceData=JumpOccureAtRowNUMTrack[i]+1;
i++;
if(i>ScienceDatFile_TobeEdited.size()-1) break;
//LOG(INFO)<<"The next I"<<i;
}while(ScienceDatFile_TobeEdited[i-1]==ScienceDatFile_TobeEdited[i]);
i--;
}
}
//LOG(INFO)<<"cff";

//for FUV
for(int i=0;i<ScienceDatFile_TobeEdited.size();i++){
if(strstr(ScienceDatFile_TobeEdited[i].c_str(),"/uvtF/")!=NULL){
//LOG(INFO)<<"III "<<ScienceDatFile_TobeEdited[i];

start_RowScienceData=0;
do{
numFUV_int++;
intToString=convertIntToStr(numFUV_int);
if(numFUV_int<10){
tempStringFUV="0"+intToString;
intToString=tempStringFUV;
}
TomakeFileString=level1dirname+"/uvit/"+orbitNum+"/uvtF/uvtF."+intToString+"/";
FullPathFuvNewString=TomakeFileString+basename(ScienceDatFile_TobeEdited[i].c_str());
Copied_location_scienceDataFile.push_back(FullPathFuvNewString);
LOG(INFO)<<TomakeFileString<<" "<<FullPathFuvNewString;
cmdTogive="mkdir -p "+TomakeFileString;
system(cmdTogive.c_str());
cmdTogive="cp  "+ScienceDatFile_TobeEdited[i]+" "+TomakeFileString;
system(cmdTogive.c_str());
gtiFile=ScienceDatFile_TobeEdited[i].replace(ScienceDatFile_TobeEdited[i].size()-4,3,"gti");
nameOfgtisubstr=gtiFile.substr(0,gtiFile.size()-1);
LOG(INFO)<<"GTI file is "<<nameOfgtisubstr;
btiFile=ScienceDatFile_TobeEdited[i].replace(ScienceDatFile_TobeEdited[i].size()-4,3,"bti");
nameOfbtisubstr=btiFile.substr(0,btiFile.size()-1);
LOG(INFO)<<"BTI file is"<<nameOfbtisubstr;
cmdTogive="cp  "+nameOfgtisubstr+" "+TomakeFileString;
system(cmdTogive.c_str());
cmdTogive="cp  "+nameOfbtisubstr+" "+TomakeFileString;
system(cmdTogive.c_str());
//for(int j=start_RowScienceData;j<JumpOccureAtRowNUMTrack[i];j++){


//}
//start_RowScienceData=JumpOccureAtRowNUMTrack[i]+1;
i++;
if(i==ScienceDatFile_TobeEdited.size()) break;
}while(ScienceDatFile_TobeEdited[i-1]==ScienceDatFile_TobeEdited[i]);
i--;
}
}

LOG(INFO)<<Copied_location_scienceDataFile.size()<<" "<<ScienceDatFile_TobeEdited.size();


//now remove the rows from L1;
fitsfile *ffile;
vector<long> nRowsToDelete;
for (int i=0;i<Copied_location_scienceDataFile.size();i++){

nRowsToDelete.clear();

for(int k=1;k<=Totlpackets[i];k++){
	if(k <StartAtRowNUMTrack[i] || k >JumpOccureAtRowNUMTrack[i]){

nRowsToDelete.push_back(k);


}
}//for loop finish




LOG(INFO)<<nRowsToDelete.size();

 fits_open_file(&ffile, (char*) Copied_location_scienceDataFile[i].c_str(), READWRITE, &status);
 printError(status, "***Error in opening file***", (char*) Copied_location_scienceDataFile[i].c_str());
 fits_movabs_hdu(ffile, 3, NULL, &status);
 printError(status, "***Error in reading HDU 3***", (char*) Copied_location_scienceDataFile[i].c_str());
  fits_delete_rowlist(ffile, nRowsToDelete.data(), nRowsToDelete.size(), &status);
  printError(status, "Error in Deleting the row list ", (char*) Copied_location_scienceDataFile[i].c_str());
  fits_close_file(ffile,&status);
  printError(status, "Error in closing the file ", (char*) Copied_location_scienceDataFile[i].c_str());

}

//now remove original Science Data File
scienceDataToEdited_backup.erase(unique(scienceDataToEdited_backup.begin(),scienceDataToEdited_backup.end()),scienceDataToEdited_backup.end());
for(int i=0;i<scienceDataToEdited_backup.size();i++){

cmdTogive="rm "+scienceDataToEdited_backup[i];//removing original science data file.
system(cmdTogive.c_str());

}


//making tar file from directory.
string newl1tar=tarname_s+"_EDITED_L2";
cmdTogive="tar -cvf "+newl1tar+" --directory="+path_l1dir_s+"/..// "+basename(path_l1dir_s.c_str());
LOG(INFO)<<"Executing "<<cmdTogive<<" ...";
system(cmdTogive.c_str());

uvitObj.zipFlag=FALSE;
uvitObj.level1indir=newl1tar;





////till this
        string name_OfTar = uvitObj.level1indir;
        if (uvitObj.manualMode == TRUE) {
            if (!DirExists((char*) uvitObj.prev_Output_DirL2.c_str())) {
                LOG(ERROR) << "***Auto mode run output directory is not available***";
                return (EXIT_FAILURE);
            }
            string uvit_dir_ofAutomode = uvitObj.prev_Output_DirL2 + "/uvit/";
            vector<string> orbitDirName;
            get_firstLevelSubDirs(uvit_dir_ofAutomode, orbitDirName);
            string OrbitNum = basename(orbitDirName[0].c_str());
            if (strstr(name_OfTar.c_str(), OrbitNum.c_str()) == NULL) {
                LOG(ERROR) << "***Level1 tar filepath and output directory for AUTO mode mismatch.***";
                return (EXIT_FAILURE);
            }

            stringstream cmd1;
            int ContinueFlag = 0;
            LOG(WARNING) << "BEFORE RUNNNING MANUAL MODE.Please take backup of your Auto mode previous output directory i.e " << uvitObj.prev_Output_DirL2 << " Present code will delete some stars from star list";
            LOG(INFO) << "If you Want to take backup than Press 1 for exit";
            cin>>ContinueFlag;
            if (ContinueFlag == 1) {
                LOG(INFO) << "EXITING....";
                return (EXIT_FAILURE);
            }
            bool refFrameFoundFlag = FALSE;
            fitsfile *frefframe;
            vector<long> numRows_Topreserve;
            float *Xloc_refFrame, *Yloc_refFrame, *Int_refFrame;
            long numRows_refframe = 0;
            vector<string> AllFilesFromDir, refFrameFilenames;
            int numStars_toPreserve, star_id;

            status = getFiles(uvitObj.prev_Output_DirL2, AllFilesFromDir);
            int lastindex_OfSlash;
            string parentDirectoryPath, FinalPath_parentDir;
            string l2outputPath_DriftSeries;
            vector<string> ReframeFilename_track;
            for (int i = 0; i < AllFilesFromDir.size(); i++) {
                if (strstr(AllFilesFromDir[i].c_str(), "_f1_rf_centroid") != NULL) {
                    LOG(INFO) << "Reference Frame Name->" << AllFilesFromDir[i].c_str();
                    fits_open_file(&frefframe, (char*) AllFilesFromDir[i].c_str(), READWRITE, &status);
                    printError(status, "***Error in opening file***", (char*) AllFilesFromDir[i].c_str());
                    numRows_refframe = 0;
                    fits_movabs_hdu(frefframe, 2, NULL, &status);
                    printError(status, "***Error in reading HDU 2***", (char*) AllFilesFromDir[i].c_str());
                    fits_get_num_rows(frefframe, &numRows_refframe, &status);
                    printError(status, "***Error in getting number of rows***", (char*) AllFilesFromDir[i].c_str());
                    Xloc_refFrame = new float[numRows_refframe];
                    Yloc_refFrame = new float[numRows_refframe];
                    Int_refFrame = new float[numRows_refframe];
                    fits_read_col(frefframe, TFLOAT, 1, 1, 1, numRows_refframe, NULL, (void*) Xloc_refFrame, NULL, &status);
                    printError(status, "Reading a column Fails in Ref-frame", (char*) AllFilesFromDir[i].c_str());
                    fits_read_col(frefframe, TFLOAT, 2, 1, 1, numRows_refframe, NULL, (void*) Yloc_refFrame, NULL, &status);
                    printError(status, "Reading a column Fails in Ref-frame", (char*) AllFilesFromDir[i].c_str());
                    fits_read_col(frefframe, TFLOAT, 3, 1, 1, numRows_refframe, NULL, (void*) Int_refFrame, NULL, &status);
                    printError(status, "Reading a column Fails in Ref-frame", (char*) AllFilesFromDir[i].c_str());


                    LOG(INFO) << "====================================================";
                    LOG(INFO) << "==========Displaying reference frame stars==========";
                    LOG(INFO) << "====================================================";
                    if (numRows_refframe == 0) continue;
                    cout << "STAR_ID" << setw(15) << "XLOCATION" << setw(15) << "YLOCATION" << setw(15) << "INTENSITY" << endl;
                    for (int j = 0; j < numRows_refframe; j++) {
                        cout << j + 1 << setw(15) << Xloc_refFrame[j] << setw(15) << Yloc_refFrame[j] << setw(15) << Int_refFrame[j] << endl;
                    }
                    LOG(INFO) << "How many stars to retain for tracking?";
                    cin>>numStars_toPreserve;
                    if (numStars_toPreserve == numRows_refframe) {
                        LOG(INFO) << "ALL stars retained";
                        continue;
                    }
                    if (numStars_toPreserve > numRows_refframe) {
                        LOG(ERROR) << "***Total number of stars are " << numRows_refframe << ", No of stars to preserve are greater than total available stars***";
                        return (EXIT_FAILURE);
                    }
                    if (numStars_toPreserve < 1) {
                        LOG(ERROR) << " **Cant remove total " << numRows_refframe - numStars_toPreserve << " star(s) are available ***";
                        return (EXIT_FAILURE);
                    }

                    //                   if(numStars_toRemove>0){
                    LOG(INFO) << "***Please Enter the stars_ids you want to retain.NOTE:After writing each id, press ENTER ";
                    //                   }
                    for (int indexStar = 0; indexStar < numStars_toPreserve; indexStar++) {
                        cout << "Star " << indexStar + 1 << " to be retained-> ";
                        cin>>star_id;
                        cout << endl;
                        numRows_Topreserve.push_back(star_id);
                    }
                    vector<long> ToRemoveRows;
                    bool flag_Tobepreserve = FALSE;
                    sort(numRows_Topreserve.begin(), numRows_Topreserve.end());
                    for (int indexTotStars = 0; indexTotStars < numRows_refframe; indexTotStars++) {
                        flag_Tobepreserve = FALSE;
                        for (int indexStar = 0; indexStar < numRows_Topreserve.size(); indexStar++) {
                            if (indexTotStars + 1 == numRows_Topreserve[indexStar]) {
                                flag_Tobepreserve = TRUE;
                                break;
                            }

                        }
                        if (flag_Tobepreserve == FALSE) {
                            ToRemoveRows.push_back(indexTotStars + 1);
                        }
                    }
                    fits_delete_rowlist(frefframe, ToRemoveRows.data(), ToRemoveRows.size(), &status);
                    printError(status, "Error in Deleting the row list ", (char*) AllFilesFromDir[i].c_str());
                    fits_close_file(frefframe, &status);
                    printError(status, "Error in closing the file ", (char*) AllFilesFromDir[i].c_str());

                    //taking directory path 
                    // break;
                    ReframeFilename_track.push_back(AllFilesFromDir[i]);
                }

            }
            string RefFrameDir;
            for (int index_OnTotalOrbitdata = 0; index_OnTotalOrbitdata < ReframeFilename_track.size(); index_OnTotalOrbitdata++) {
                //put loop finish here.
                lastindex_OfSlash = ReframeFilename_track[index_OnTotalOrbitdata].find_last_of("/");
                parentDirectoryPath = ReframeFilename_track[index_OnTotalOrbitdata].substr(0, lastindex_OfSlash);
                lastindex_OfSlash = parentDirectoryPath.find_last_of("/");
                FinalPath_parentDir = parentDirectoryPath.substr(0, lastindex_OfSlash);
                RefFrameDir.assign(FinalPath_parentDir);
                lastindex_OfSlash = FinalPath_parentDir.find_last_of("/");
                parentDirectoryPath = FinalPath_parentDir.substr(0, lastindex_OfSlash);
                LOG(INFO) << RefFrameDir;
                LOG(INFO) << parentDirectoryPath;
                l2outputPath_DriftSeries = parentDirectoryPath + "/" + "uvtComputeDrift_manual";
                LOG(INFO) << l2outputPath_DriftSeries;

                cmd1 << EXE_DRIFTCOMPUTE_IM << " inputdatadir = " << RefFrameDir << " GenMatchStarsFile_flag =" << ReturnYorN(uvitObj.match_stars_file_flag) << " diff_Dist =" << uvitObj.nbhd_dist << " freqDomainFilter_Flag =" << uvitObj.freqDomainFilter_Flag <<
                        " type_Filtering= " << uvitObj.type_Filtering << " fitting_flag =" << ReturnYorN(uvitObj.fitting_flag) << " order_pitch=" << uvitObj.orderpitch << " order_yaw =" << uvitObj.orderyaw <<
                        " order_roll =" << uvitObj.orderroll << " delta_Time=" << uvitObj.delta_time << " freq_value=" << uvitObj.freqvalue << " shiftRotDetAlgoFlag=" << uvitObj.shift_N_Rotate_algopc << " flagTheta =" << ReturnYorN(uvitObj.flag_thetaComp) <<
                        " algoFlag =" << uvitObj.star_detect_algo_flag << " outdir =" << l2outputPath_DriftSeries << " clobber =" << ReturnYorN(uvitObj.clobber) << " history =" << ReturnYorN(uvitObj.history) <<
                        " mode =" << uvitObj.mode;


                LOG(INFO) << cmd1.str();
                try {
                    system(cmd1.str().c_str());
                } catch (...) {
                    LOG(ERROR) << "***Error in Drift series generation";
                    return (EXIT_FAILURE);
                }

                refFrameFoundFlag = TRUE;



                if (refFrameFoundFlag == FALSE) {
                    LOG(ERROR) << "MANUAL MODE cant run!!! NO reference found on " << uvitObj.prev_Output_DirL2 << " directory";
                    return (EXIT_FAILURE);
                }
                cmd1.str(std::string()); //for clearing 
            }

            uvitObj.level2outdir = uvitObj.prev_Output_DirL2; //+"_ManualMode";
            //need to extract because RelativeAspect will not run now.
            string temp_str;
            string orbnum;
            temp_str.assign(uvitObj.level1indir);
            if (uvitObj.zipFlag == FALSE) {


                status = extractTars(temp_str, uvitObj.level1indir, orbnum); //extract level-1 data  tar file

                if (status) {
                    LOG(INFO) << "Error in extracting tar";
                    return (EXIT_FAILURE);
                }
            } else {
                status = extractZip(temp_str, uvitObj.level1indir, orbnum); //extract level-1 data  tar file

                if (status) {
                    LOG(INFO) << "Error in extracting tar";
                    return (EXIT_FAILURE);
                }
            }
        }
   
        string filename = TIME_MISMATCH_FILENAME;
        string cmd = "rm " + filename;
        system(cmd.c_str());

        string tarname;
        //tarname = uvitObj.level1indir;
        tarname = name_OfTar;
        string strtemp(tarname);
        int pos = strtemp.find_last_of("/");
        string tempFilename = "";
        //string cmd;
        string tar_date_data = tarname.substr(pos + 1 + 11, 8);
        string tar_OBS_id_data = tarname.substr(pos + 1 + 19, 21);
        string outputNUV, outputFUV;
        string orbitnum = tarname.substr(pos + 1 + 41, 5);
        LOG(INFO) << "The orbit num is " << orbitnum;

        string path_l1dir;
        if (pos > 0 && pos < strtemp.length()) {
            tempFilename = strtemp.substr(0, pos);
            path_l1dir = tempFilename + "/" + tar_date_data + "_" + tar_OBS_id_data + "_level1";

        } else {
            path_l1dir = tar_date_data + "_" + tar_OBS_id_data + "_level1";
        }
        if (uvitObj.manualMode == FALSE && uvitObj.NUVonNUVflag == FALSE  && uvitObj.FUVonFUVflag==FALSE) {

                       if (DirExists((char*) path_l1dir.c_str())) {
                           cmd = "rm -r " + path_l1dir;
                            system(cmd.c_str());
                       }
                        status = uvitObj.RunRelativeAspectChain();
                        if (status) {
                          LOG(ERROR) << "Error in Running the Relative Aspect pipeline for the dataset";
                         return (EXIT_FAILURE);
                       }


        }

        if (uvitObj.manualMode == FALSE && uvitObj.NUVonNUVflag == TRUE) {
            if (DirExists((char*) path_l1dir.c_str())) {
                cmd = "rm -r " + path_l1dir;
                system(cmd.c_str());
            }
            uvitObj.level1indirrapc = name_OfTar;
            uvitObj.level2outdirrapc = uvitObj.level2outdir;
            uvitObj.channelRAPC="NUV";
             status = uvitObj.RunRAPCChain();
            if (status) {
                LOG(ERROR) << "Error in Running the Relative Aspect pipeline for the dataset";
                return (EXIT_FAILURE);
            }
        }
        else if(uvitObj.manualMode == FALSE && uvitObj.FUVonFUVflag==TRUE ){
            if (DirExists((char*) path_l1dir.c_str())) {
                cmd = "rm -r " + path_l1dir;
                system(cmd.c_str());
            }
            uvitObj.level1indirrapc = name_OfTar;
            uvitObj.level2outdirrapc = uvitObj.level2outdir;
             uvitObj.channelRAPC="FUV";
             status = uvitObj.RunRAPCChain();
            if (status) {
                LOG(ERROR) << "Error in Running the Relative Aspect pipeline for the dataset";
                return (EXIT_FAILURE);
            }
        }
        string temp_str_l2; //=basename(uvitObj.level1indir.c_str ());

        //        searching Attitude file path in level1 directory
        vector<string> l1files;
        LOG(INFO) << path_l1dir;
        char ATTfile[NAMESIZE];
        status = getFiles(path_l1dir, l1files);
        bool fileAvailcheck = FALSE;
        for (int i = 0; i < l1files.size(); i++) {
            if (strstr(l1files[i].c_str(), ".att") != NULL) {
                strcpy(ATTfile, l1files[i].c_str());
                fileAvailcheck = TRUE;
                break;
            }

        }
        if (!fileAvailcheck) {
            LOG(ERROR) << "Attitude file is not available in the level1"; //This error is not needed actually.
            return (EXIT_FAILURE);
        }

        if (FileExists(TOTAL_NUV_FILENAME)) {
            cmd = "rm " + (string) TOTAL_NUV_FILENAME;
            system(cmd.c_str());
        }


        //Reading RAS file from Directory
        vector<string> AllFiles, RASfiles;
        LOG(INFO) << uvitObj.level2outdir;
        status = getFiles(uvitObj.level2outdir, AllFiles);
        LOG(INFO) << AllFiles.size();
        bool flag_rasFound = FALSE;
        for (int i = 0; i < AllFiles.size(); i++) {
            if (strstr(AllFiles[i].c_str(), "dr.fits") != NULL) {
                flag_rasFound = TRUE;
                strcpy(uvitObj.rasfilepc, AllFiles[i].c_str()); //RAS file path read completed.
                LOG(INFO) << uvitObj.rasfilepc;
                RASfiles.push_back(uvitObj.rasfilepc);
            }
        }
        uvitObj.level1indirpc = path_l1dir;

        //Loop for total number of RAS file genereated after uvtRelativeAspectIM
        string finalDir;
        string nuv_subDir, fuv_subDir;
        if (!flag_rasFound) {
            if (uvitObj.NUVonNUVflag == FALSE && uvitObj.FUVonFUVflag ==FALSE) {
                LOG(ERROR) << "***NO RAS file found for particular Dataset!!!,CHECKING NUV DATA...";
                LOG(ERROR) << "\033[1;31mUVT RELATIVE ASPECT PC WILL RUN ON NUV DATA....\033[0m";
                vector<string> Allfilenames_level1;
                string filename_nuvs = TOTAL_NUV_FILENAME;
                ofstream if_nuv;
                if_nuv.open((char*) filename_nuvs.c_str());
                LOG(INFO) << uvitObj.level1indir << endl;
                status = getFiles(uvitObj.level1indirpc, Allfilenames_level1);

                for (int i = 0; i < Allfilenames_level1.size(); i++) {

                    if ((strstr(Allfilenames_level1[i].c_str(), "uvtN") != NULL || strstr(Allfilenames_level1[i].c_str(), "uvtF") != NULL) && strstr(Allfilenames_level1[i].c_str(), ".fits") != NULL) {

                        if_nuv << Allfilenames_level1[i] << endl;
                        LOG(INFO) << Allfilenames_level1[i];
                    }
                }

                if_nuv.close();

                string matched_Filenames = TIME_MISMATCH_FILENAME;
                if_nuv.open((char*) matched_Filenames.c_str());
                if_nuv.close();
            }
            //return (EXIT_FAILURE);
        }



        if (flag_rasFound) {
            //setting level1 directory path
            uvitObj.last_FileFlag = FALSE;

            for (int i = 0; i < RASfiles.size(); i++) {
                //              for (int i = 0; i < 1; i++) {   
                if (i == RASfiles.size() - 1) {
                    uvitObj.last_FileFlag = TRUE;
                }

                //Running L2 pipelinr for NUV channel.
                LOG(INFO) << "L2PC runnnig for RAS file path " << RASfiles[i];
                strcpy(uvitObj.rasfilepc, RASfiles[i].c_str());
                 status = uvitObj.RunL2PCChain("NUV", i + 1);
                if (status) {
                    LOG(ERROR) << "Error in running level2 PC chain";
                    continue;
                }
                outputNUV = uvitObj.level2outdirc;

                temp_str_l2 = basename(uvitObj.level1indirpc.c_str());
                temp_str_l2.replace(30, 7, "_level2");
                finalDir = temp_str_l2;
                cmd = "rm -r " + finalDir;
                LOG(INFO) << cmd;
                system(cmd.c_str());
                cmd = "mkdir -p " + finalDir;
                LOG(INFO) << cmd;
                system(cmd.c_str());
                LOG(INFO) << uvitObj.level2outdirc;
                cmd = "cp -r  " + (string) uvitObj.level2outdirc + "/* " + finalDir;
                system(cmd.c_str());
                LOG(INFO) << cmd;
                nuv_subDir = uvitObj.outtarpathpc + (string) "/N" + convertIntToStr(i + 1);
                cmd = "mkdir -p " + nuv_subDir;
                system(cmd.c_str());

                if (uvitObj.zipFlag == FALSE) {
                    status = compressTars(uvitObj.level1indir, finalDir, nuv_subDir);
                    //       status=compressTars (uvitObj.level1indir,finalDir,uvitObj.outtarpathpc);
                    if (status) {
                        LOG(ERROR) << "***Error in compressing the tar***";

                    }
                } else if (uvitObj.zipFlag == TRUE) {
                    status = compressZip(uvitObj.level1indir, finalDir, nuv_subDir);
                    if (status) {
                        LOG(ERROR) << "***Error in compressing the zip***";

                    }
                }
                //    string cmd;
                // Running L2 pipeline for FUV channel;
                 status = uvitObj.RunL2PCChain("FUV", i + 1);
                if (status) {
                    LOG(ERROR) << "Error in running level2 PC chain";
                    continue;
                }
                //added
                temp_str_l2 = basename(uvitObj.level1indirpc.c_str());
                temp_str_l2.replace(30, 7, "_level2");
                finalDir = temp_str_l2;
                cmd = "rm -r " + finalDir;
                system(cmd.c_str());
                cmd = "mkdir -p " + finalDir;
                system(cmd.c_str());

                cmd = "cp -r  " + (string) uvitObj.level2outdirc + "/* " + finalDir;
                system(cmd.c_str());
                fuv_subDir = uvitObj.outtarpathpcfuv + (string) "/F" + convertIntToStr(i + 1);
                cmd = "mkdir -p " + fuv_subDir;
                system(cmd.c_str());
                if (uvitObj.zipFlag == FALSE) {
                    status = compressTars(uvitObj.level1indir, finalDir, fuv_subDir);
                    if (status) {
                        LOG(ERROR) << "***Error in compressing the tar***";

                    }
                } else {
                    status = compressZip(uvitObj.level1indir, finalDir, fuv_subDir);
                    if (status) {
                        LOG(ERROR) << "***Error in compressing the zip***";

                    }
                }


            }
            //till this
        }
        vector<string> Total_pc_Filename;
        vector<string> files_matched_With_VIS;
        vector<string> remaining_Files;
        if (uvitObj.NUVonNUVflag == FALSE && uvitObj.FUVonFUVflag==FALSE) {
            filename = TIME_MISMATCH_FILENAME;
            ifstream if1;
            if1.open((char*) filename.c_str());
            string filename_pc;

            while (!if1.eof()) {
                if1>>filename_pc;

                if (if1.eof()) {
                    break;
                }
                files_matched_With_VIS.push_back(filename_pc);
            }
            if1.close();


            string filename2 = TOTAL_NUV_FILENAME;
            ifstream if2;
            if2.open((char*) filename2.c_str());
            while (!if2.eof()) {
                if2>>filename_pc;

                if (if2.eof()) {
                    break;
                }
                Total_pc_Filename.push_back(filename_pc);
            }

            if2.close();

            RemoveDuplicateStringsFromVECT(Total_pc_Filename);


            //comparison
            for (int i = 0; i < files_matched_With_VIS.size(); i++) {
                LOG(INFO) << files_matched_With_VIS[i];

            }

            bool flag_found = FALSE;
            for (int i = 0; i < Total_pc_Filename.size(); i++) {
                flag_found = FALSE;
                for (int j = 0; j < files_matched_With_VIS.size(); j++) {

                    if (strcmp(Total_pc_Filename[i].c_str(), files_matched_With_VIS[j].c_str()) == 0) {
                        flag_found = TRUE;
                        break;
                    }
                }
                if (flag_found == FALSE) {

                    remaining_Files.push_back(Total_pc_Filename[i]); //this contains the science data file for which Relative aspect PC to be run;
                }


            }
        }
        //Added
        //AND TO BE REMOVED
        //remaining_Files.clear();
        //now making a temporary directory
        temp_str_l2 = basename(uvitObj.level1indirpc.c_str());
        temp_str_l2.replace(30, 7, "_level1_temp1");
        string finalDir_temp = temp_str_l2;


        //making the directory for temp PC science data files & copying the files 
        string temp_Dir;
        string level2_Origname = uvitObj.level2outdir;
        cmd = "rm -r  " + finalDir_temp;
        system(cmd.c_str());

        char temp_char[FLEN_FILENAME];
        if (remaining_Files.size() > 0) {
            LOG(INFO) << "\033[1;31mNO VIS OVERLAP FOUND FOR FOLLOWING FILES...\033[0m";
            LOG(INFO) << "\033[1;31m=====================================================================================\033[0m";
            for (int i = 0; i < remaining_Files.size(); i++) {
                LOG(INFO) << "\033[1;31m" << remaining_Files[i] << "\033[0m";
            }
            LOG(INFO) << "\033[1;31m=====================================================================================\033[0m";
            LOG(INFO) << "\033[1;31mNOW RUNNIG RELATIVE ASPECT PC CHAIN ON THESE FILES..\033[0m";
        }

        string tempstring;
        int sizestr;
        vector<string> prefix_ofGTI;

        for (int i = 0; i < remaining_Files.size(); i++) {
            tempstring = remaining_Files[i].c_str();

            sizestr = strlen(tempstring.c_str());
            LOG(INFO) << sizestr;
            tempstring.replace(sizestr - 4, 4, "gti");
            LOG(INFO) << tempstring;
            prefix_ofGTI.push_back(tempstring);
        }

        //now searching GTI file.
        vector<string> AllFiles_l1;
        vector<string> GTI_FILEPATH;
        getFiles(uvitObj.level1indirpc, AllFiles_l1);
        LOG(INFO) << AllFiles_l1.size() << " " << prefix_ofGTI.size();
        for (int i = 0; i < AllFiles_l1.size(); i++) {

            if ((strstr(AllFiles_l1[i].c_str(), "uvtN") != NULL || strstr(AllFiles_l1[i].c_str(), "uvtF") != NULL)&& (strstr(AllFiles_l1[i].c_str(), ".gti") != NULL)) {
                LOG(INFO) << AllFiles_l1[i];
                for (int j = 0; j < prefix_ofGTI.size(); j++) {
                    if (strstr(AllFiles_l1[i].c_str(), prefix_ofGTI[j].c_str()) != NULL) {
                        GTI_FILEPATH.push_back(AllFiles_l1[i]);
                        LOG(INFO) << GTI_FILEPATH[j];
                    }
                }

            }

        }


        LOG(INFO) << remaining_Files.size() << " " << GTI_FILEPATH.size();
        int cnt_n = 0;
        int cnt_f = 0;
        for (int i = 0; i < remaining_Files.size(); i++) {

            if (strstr(remaining_Files[i].c_str(), "uvtN") != NULL) {
                cnt_n++;
                sprintf(temp_char, "%02d", cnt_n);
                temp_Dir = finalDir_temp + "/uvit/" + orbitnum + "/uvtN/uvtN." + (string) temp_char + "/";
                cmd = "mkdir -p " + temp_Dir;
                system(cmd.c_str());
                cmd = "cp -r " + remaining_Files[i] + " " + temp_Dir;
                system(cmd.c_str());
                cmd = "cp -r " + GTI_FILEPATH[i] + " " + temp_Dir;
                system(cmd.c_str());
            } else if (strstr(remaining_Files[i].c_str(), "uvtF") != NULL) {
                cnt_f++;
                sprintf(temp_char, "%02d", cnt_f);
                temp_Dir = finalDir_temp + "/uvit/" + orbitnum + "/uvtF/uvtF." + (string) temp_char + "/";
                cmd = "mkdir -p " + temp_Dir;
                system(cmd.c_str());
                cmd = "cp -r " + remaining_Files[i] + " " + temp_Dir;
                system(cmd.c_str());
                cmd = "cp -r " + GTI_FILEPATH[i] + " " + temp_Dir;
                system(cmd.c_str());
            }

        }




        if (remaining_Files.size() > 0) {

            //Aux directory to be copied
            temp_Dir = finalDir_temp + "/uvit/" + orbitnum + "/aux/";
            string temp_aux_path = uvitObj.level1indirpc + "/uvit/" + orbitnum + "/aux/*";
            cmd = "mkdir -p " + temp_Dir;
            system(cmd.c_str());
            cmd = "cp -r " + temp_aux_path + " " + temp_Dir;
            system(cmd.c_str());
            temp_Dir = uvitObj.level1indirpc + "/uvit/" + orbitnum + "/uvtN/DarkN/*";
            string namePath = finalDir_temp + "/uvit/" + orbitnum + "/uvtN/DarkN/";
            cmd = "mkdir -p " + namePath;
            system(cmd.c_str());
            cmd = "cp -r " + temp_Dir + " " + namePath;
            system(cmd.c_str());
            temp_Dir = uvitObj.level1indirpc + "/uvit/" + orbitnum + "/uvtF/DarkF/*";
            namePath = finalDir_temp + "/uvit/" + orbitnum + "/uvtF/DarkF/";
            cmd = "mkdir -p " + namePath;
            system(cmd.c_str());
            cmd = "cp -r " + temp_Dir + " " + namePath;
            system(cmd.c_str());

            //copy .mkf file.
            string mkfpathout = finalDir_temp + "/uvit/" + orbitnum + "/";
            string mkfpathin = uvitObj.level1indirpc + "/uvit/" + orbitnum + "/*.mkf";
            cmd = "cp -r " + mkfpathin + " " + mkfpathout;
            system(cmd.c_str());


            string l2tar_pcOnly;
            if (uvitObj.zipFlag == FALSE) {
                l2tar_pcOnly = uvitObj.level1indir + "_PCONLY.tar.tgz";
            } else {
                l2tar_pcOnly = uvitObj.level1indir + "_PCONLY.zip";
            }
            cmd = "rm -r " + uvitObj.level1indirpc + "_original";
            system(cmd.c_str());
            cmd = "mv " + uvitObj.level1indirpc + " " + uvitObj.level1indirpc + "_original";
            system(cmd.c_str());
            cmd = "mv " + finalDir_temp + " " + uvitObj.level1indirpc;
            system(cmd.c_str());
            LOG(INFO) << cmd;
            if (uvitObj.zipFlag == TRUE) {
                cmd = "zip -r " + l2tar_pcOnly + " " + uvitObj.level1indirpc;
            } else {
                cmd = "tar -cvf " + l2tar_pcOnly + " " + uvitObj.level1indirpc;
            }

            LOG(INFO) << cmd;
            system(cmd.c_str());



            //now run relativeAspect PC 


            uvitObj.level1indirrapc = l2tar_pcOnly;
            uvitObj.level2outdirrapc = uvitObj.level2outdir + "_RAPC";

             status = uvitObj.RunRAPCChain();

            string temp_str_l2rapc; //=basename(uvitObj.level1indir.c_str ());

            //Reading RAS file from Directory
            vector<string> AllFilesrapc, RASfilesrapc;

            status = getFiles(uvitObj.level2outdirrapc, AllFilesrapc);

            bool flag_rasFoundrapc = FALSE;
            for (int i = 0; i < AllFilesrapc.size(); i++) {
                if (strstr(AllFilesrapc[i].c_str(), "dr.fits") != NULL) {
                    flag_rasFoundrapc = TRUE;
                    strcpy(uvitObj.rasfilepc, AllFilesrapc[i].c_str()); //RAS file path read completed.
                    LOG(INFO) << uvitObj.rasfilepc;
                    RASfilesrapc.push_back(uvitObj.rasfilepc);
                }
            }

            if (!flag_rasFoundrapc) {
                LOG(ERROR) << "WARNING->***NO RAS file found for particular Dataset!!!,CANT go further..EXITING***";

            }
            //setting level1 directory path
            string tarnamerapc = uvitObj.level1indirrapc;
            string strtemprapc(tarnamerapc);
            int posrapc = strtemprapc.find_last_of("/");
            string tempFilenamerapc = "";
            //string cmd;
            string tar_date_datarapc = tarnamerapc.substr(posrapc + 1 + 11, 8);
            string tar_OBS_id_datarapc = tarnamerapc.substr(posrapc + 1 + 19, 21);
            string outputNUVrapc, outputFUVrapc;
            string orbitnumrapc = tarnamerapc.substr(posrapc + 1 + 41, 5);
            //LOG(INFO)<<"The orbit num is "<<orbitnumrapc;
            if (posrapc > 0 && posrapc < strtemprapc.length()) {
                tempFilenamerapc = strtemprapc.substr(0, posrapc);
                uvitObj.level1indirrapc = tempFilenamerapc + "/" + tar_date_datarapc + "_" + tar_OBS_id_datarapc + "_level1_temp";

            } else {
                uvitObj.level1indirrapc = tar_date_datarapc + "_" + tar_OBS_id_datarapc + "_level1_temp";
            }

            //Loop for total number of RAS file genereated after uvtRelativeAspectIM
            string finalDirrapc;
            string nuv_subDirrapc, fuv_subDirrapc;
            uvitObj.level2outdir = uvitObj.level2outdir + "_RAPC";
            for (int i = 0; i < RASfilesrapc.size(); i++) {
                //Running L2 pipelinr for NUV channel.
                LOG(INFO) << "L2PC runnnig for RAS file path " << RASfilesrapc[i];

                strcpy(uvitObj.rasfilepc, RASfilesrapc[i].c_str());
                status = uvitObj.RunL2PCChain("NUV", i + 1);
                if (status) {
                    LOG(ERROR) << "Error in running level2 PC chain";
                    continue;
                }

                outputNUVrapc = uvitObj.level2outdirrapc;

                temp_str_l2rapc = basename(uvitObj.level1indirrapc.c_str());
                //                temp_str_l2rapc.replace(30, 7, "_level2_temp");
                temp_str_l2rapc.replace(30, 7, "_level2_temp");
                finalDirrapc = temp_str_l2rapc;
                cmd = "rm -r " + finalDirrapc;
                LOG(INFO) << cmd;
                system(cmd.c_str());
                cmd = "mkdir -p " + finalDirrapc;
                LOG(INFO) << cmd;
                system(cmd.c_str());
                LOG(INFO) << uvitObj.level2outdirrapc;
                cmd = "cp -r  " + (string) uvitObj.level2outdirrapc + "/* " + finalDirrapc;
                system(cmd.c_str());
                LOG(INFO) << cmd;
                nuv_subDirrapc = uvitObj.outtarpathpc + (string) "/N" + convertIntToStr(i + 1);
                cmd = "mkdir -p " + nuv_subDirrapc;
                system(cmd.c_str());
                if (uvitObj.zipFlag == FALSE) {
                    status = compressTars(uvitObj.level1indir, finalDirrapc, nuv_subDirrapc);
                    //       status=compressTars (uvitObj.level1indir,finalDir,uvitObj.outtarpathpc);
                    if (status) {
                        LOG(ERROR) << "***Error in compressing the tar***";

                    }
                } else {
                    status = compressZip(uvitObj.level1indir, finalDirrapc, nuv_subDirrapc);
                    //       status=compressTars (uvitObj.level1indir,finalDir,uvitObj.outtarpathpc);
                    if (status) {
                        LOG(ERROR) << "***Error in compressing the zip***";

                    }
                }

                //    string cmd;
                //Running L2 pipeline for FUV channel;
                 status = uvitObj.RunL2PCChain("FUV", i + 1);
                if (status) {
                    LOG(ERROR) << "Error in running level2 PC chain";
                    continue;
                }
                //added
                temp_str_l2rapc = basename(uvitObj.level1indirrapc.c_str());
                temp_str_l2rapc.replace(30, 7, "_level2_temp");
                finalDirrapc = temp_str_l2rapc;
                cmd = "rm -r " + finalDirrapc;
                system(cmd.c_str());
                cmd = "mkdir -p " + finalDirrapc;
                system(cmd.c_str());

                cmd = "cp -r  " + (string) uvitObj.level2outdirrapc + "/* " + finalDirrapc;
                system(cmd.c_str());
                fuv_subDirrapc = uvitObj.outtarpathpcfuv + (string) "/F" + convertIntToStr(i + 1);
                cmd = "mkdir -p " + fuv_subDirrapc;
                system(cmd.c_str());
                if (uvitObj.zipFlag == FALSE) {
                    status = compressTars(uvitObj.level1indir, finalDirrapc, fuv_subDirrapc);
                    if (status) {
                        LOG(ERROR) << "***Error in compressing the tar***";

                    }
                } else {
                    status = compressZip(uvitObj.level1indir, finalDirrapc, fuv_subDirrapc);
                    if (status) {
                        LOG(ERROR) << "***Error in compressing the tar***";

                    }

                }
            }







        }
        //        
        //        
        //        //NOw combine All of output RADEC files to one file.

        vector<string> files_OfDirNUV;
//        float * FinalArray = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
//        initArray(FinalArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);

        float * FinalArray = new float[uvitObj.IMG_DIM_FIpc  * uvitObj.IMG_DIM_FIpc ];
        initArray(FinalArray, uvitObj.IMG_DIM_FIpc  * uvitObj.IMG_DIM_FIpc , 0.0f);
      
        float * counterArray = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
        initArray(counterArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);

        float * counterArray_exp = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
        initArray(counterArray_exp, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);

        float * Final_ExpArray = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
        initArray(Final_ExpArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);

        float *SubdividedImage = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
        initArray(SubdividedImage, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);

        float *Subdivided_expImage = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
        initArray(Subdivided_expImage, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);


        float * Frame_data = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
        initArray(Frame_data, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);

        float *Frame_Exposure_data = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
        initArray(Frame_Exposure_data, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);



        float * Frame_FinalArray = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
        initArray(Frame_FinalArray, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);

        float * Frame_FinalExpArray = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
        initArray(Frame_FinalExpArray, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);


        float * Frame_data_cum = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
        initArray(Frame_data_cum, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);

        float * Frame_Expdata_cum = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
        initArray(Frame_Expdata_cum, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);


        LOG(INFO) << uvitObj.OutputDir_l2pcNUV.size();
        vector<float> NonZeroInt;
        float sd_temp, mean_Ofimage;
        float thr_peaks;
        vector<string> Full_framePath, snrFilePath, ExpPath,eventList_RADEC;
        vector<float> numStars;
        long cnt_stars = 0;
        //this loop finds the reference frame which should be taken as a reference.
        numStars.clear();
        double Rollangle;
        vector<float> RollAngles_track;
        float sd_multfactor = uvitObj.sd_multi_factor_defaultpc;
        int refineWinSizedefault = uvitObj.refine_Winsizepc;
        int centWinSizedefault = uvitObj.centroid_Winsizepc;
        vector<float> cx_ref, cy_ref, ci_ref; //stores the centroid values for reference frame.
        int min_stars_match, cnt, matching_points;
        vector<float> x_ref_arr, y_ref_arr, x_arr, y_arr, temp_x_arr, temp_y_arr;
        double x_dx = 0.0, y_dy = 0.0, theta_dt = 0.0;
        vector<float> New_X_ref, New_Y_ref, New_x_arr, New_y_arr, New_Xdiff, New_Ydiff;
        vector<float> int_new_one, int_new_two;
        double ctheta, stheta;
        int x_index, y_index;
        double new_index_x, new_index_y;
        float original_diff_Distpc = uvitObj.diff_Distpc;
        float *RADECimage = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
        float *RADECExpimage = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];

        float *X_loc_snr;
        float *Y_loc_snr;
        float *bad_flag_snr;
        float *mult_phtn_snr;
        float *enp_snr;
        long nrows_snr;
        double x1, x2;
        float k1, k2, mid_X_new, mid_Y_new;


        long fpixel[2];
        fpixel[0] = fpixel[1] = 1;
        char outfile[NAMESIZE];
        long naxes[2];
        int naxis;
        int bitpix;
        string finalImagePath;
        char finalfilename[NAMESIZE];
        char finalExpfilename[NAMESIZE];
        char finalNoiseMapfilename[NAMESIZE];
        string outputFullFrameAst;
        string filenamesnr; //this is not need here because we have an image.
        fitsfile *fout;
        fitsfile *fptr_snr; //pointer to open shift and rotation file.
        char filter_val[NAMESIZE];
        char curr_filtercal_val[NAMESIZE];
    //       uvitObj.OutputDir_l2pcNUV.clear();
      //                  uvitObj.OutputDir_l2pcNUV.push_back("_NUV_3");
        //                uvitObj.OutputDir_l2pcNUV.push_back("_NUV_4");
          //              uvitObj.OutputDir_l2pcNUV.push_back("_NUV_2");
//                        uvitObj.OutputDir_l2pcNUV.push_back("_NUV_13");
//                        uvitObj.OutputDir_l2pcNUV.push_back("out5962_NUV_7");
                        
                        
/*
                        uvitObj.OutputDir_l2pcNUV.clear();
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_1");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_2");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_3");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_4");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_5");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_6");
 uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_7");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_8");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_9");
 uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_10");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_11");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_12");
 uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_13");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_14");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_15");
 uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_16");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_17");
                        uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_18");
*/
                       // uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_11");
                       // uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_12");
		//	uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_13");
                  //      uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_14");
                    //    uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_15");
                      //  uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_16");
		//	uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_17");
                  //      uvitObj.OutputDir_l2pcNUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_NUV_18");
                        
        //bool flag_filtervalid=TRUE;
        vector<string> window_index_track;
        window_index_track.push_back("99");
        window_index_track.push_back("149");
        window_index_track.push_back("199");
        window_index_track.push_back("249");
        window_index_track.push_back("299");
        window_index_track.push_back("349");
        window_index_track.push_back("511");
        int Win_Size = 0;
        string Win_Size_str;

        //vector<string> header_infoL1;

        vector< vector<string> > track_headerInfoL1;
        vector<string> First_HDU_headerTrack;
        float Exp_Time_Keyword;
         bool flag_TestOnecheck=FALSE;
         bool flag_TestTwocheck=FALSE;
         float Sum_Exp_Time=0.0f;
         vector<string> IndividualOrbit_Dir;
           string current_ExtName;
bool flagRefFrameNotFound=FALSE;
        if (uvitObj.OutputDir_l2pcNUV.size() > 0) {
            for (int window_index = 0; window_index < window_index_track.size(); window_index++) {
                for (int i = 1; i <= 7; i++) {
			flagRefFrameNotFound=FALSE;
                    initArray(Frame_data_cum, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                    //initArray(FinalArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                    initArray(FinalArray, uvitObj.IMG_DIM_FIpc  * uvitObj.IMG_DIM_FIpc , 0.0f);
                    initArray(counterArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);                    
                    initArray(counterArray_exp, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                    initArray(Final_ExpArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                    initArray(Frame_Expdata_cum, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                    IndividualOrbit_Dir.size();
                    Full_framePath.clear();
                    eventList_RADEC.clear();
                    ExpPath.clear();
                    RollAngles_track.clear();
                    snrFilePath.clear();
                    numStars.clear();
                    cx_ref.clear();
                    cy_ref.clear();
                    ci_ref.clear();
                    // flag_filtervalid=TRUE;
                    sprintf(curr_filtercal_val, "%s%d", "F", (int) i);
                    //LOG(INFO)<<" 1->"<<curr_filtercal_val;
                    int cnt_integratinTime = 0;
                    double integrationTime = 0.0; //storing individual integration time for individual frame.
                    double Final_integrationTime = 0.0;

                  
                    for (int i = 0; i < uvitObj.OutputDir_l2pcNUV.size(); i++) {
                        //flag_filtervalid=TRUE;
                        
                        files_OfDirNUV.clear();
                        getFiles(uvitObj.OutputDir_l2pcNUV[i], files_OfDirNUV);
                        if (files_OfDirNUV.size() == 0 && window_index == 0 && i == 1) {
                            cmd = "rm -r " + uvitObj.OutputDir_l2pcNUV[i];
                            LOG(INFO) << cmd;
                            system(cmd.c_str());
                        }

                        for (int j = 0; j < files_OfDirNUV.size(); j++) {
                            // LOG(INFO) << files_OfDirNUV[j];



                            //                if (strstr(files_OfDirNUV[j].c_str(), "ra-dec") != NULL) {
                            //if (strstr(files_OfDirNUV[j].c_str(), "sig_regAvg") != NULL) {
				 if (strstr(files_OfDirNUV[j].c_str(), "_sig_ra-dec.") != NULL) {
//if (strstr(files_OfDirNUV[j].c_str(), "ra-dec.fits") != NULL) {                                
	                  	LOG(INFO) << "Reading " << files_OfDirNUV[j];
                                readImage((char*) files_OfDirNUV[j].c_str(), 1, Frame_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc);
                                readKeywords((char*) files_OfDirNUV[j].c_str(), 1, 5, TDOUBLE, "INT_TIME", &integrationTime, TDOUBLE, "ROLLAPPLIED", &Rollangle, TSTRING, "FILTER", &filter_val, TINT, "WIN_X_SZ", &Win_Size,TFLOAT,"EXP_TIME",&Exp_Time_Keyword);
                                //LOG(INFO)<<"333 "<<filter_val<<" "<<atoi(window_index_track[window_index].c_str())<<" "<<Win_Size;
                                
                                if (strcmp(filter_val, curr_filtercal_val) != 0 || atoi(window_index_track[window_index].c_str()) != Win_Size) {
                                    LOG(INFO) << "NO filter match found ...Checking next file..";
                                    //                            Full_framePath.clear();
                                    //                            snrFilePath.clear();

                                    break;
                                }


                                Final_integrationTime = Final_integrationTime + integrationTime;
                                cnt_integratinTime++;

                                Full_framePath.push_back(files_OfDirNUV[j]);
                                RollAngles_track.push_back(Rollangle);


//                                NonZeroInt.clear();
//                                for (int k = 0; k < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; k++) {
//                                    if (Frame_data[i] != 0) {
//                                        NonZeroInt.push_back(Frame_data[k]);
//                                    }
//                                }
//
//                                sd_temp = getSD(NonZeroInt.data(), NonZeroInt.size());
//                                mean_Ofimage = getmean(NonZeroInt.data(), NonZeroInt.size());
//                                //    thr =  sd_temp* sd_mul_factor ;
//                                thr_peaks = mean_Ofimage + sd_temp * 3.5;
//                                cnt_stars = 0;
//                                for (int index = 0; index < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; index++) {
//                                    if (Frame_data[index] > thr_peaks && Frame_data[index] > 0) {
//                                        cnt_stars++; //stars calculation
//
//                                    }
//                                }
                                numStars.push_back(Exp_Time_Keyword);


                            }
                            //if (strstr(files_OfDirNUV[j].c_str(), "exp_regAvg") != NULL) {
if (strstr(files_OfDirNUV[j].c_str(), "sig_ra-decexp.") != NULL) {                              
//if (strstr(files_OfDirNUV[j].c_str(), "ra-decexp.fits") != NULL) {    
LOG(INFO) << "Reading " << files_OfDirNUV[j];
                                readKeywords((char*) files_OfDirNUV[j].c_str(), 1, 2, TSTRING, "FILTER", &filter_val, TINT, "WIN_X_SZ", &Win_Size);
                                if (strcmp(filter_val, curr_filtercal_val) != 0 || atoi(window_index_track[window_index].c_str()) != Win_Size) {
                                    LOG(INFO) << "NO filter match found ...Checking next file..";
                                    //                            Full_framePath.clear();
                                    //                            snrFilePath.clear();

                                    break;
                                }
                                ExpPath.push_back(files_OfDirNUV[j]);
                            }
                                 if (strstr(files_OfDirNUV[j].c_str(), "l2_radec.") != NULL) {                              
//if (strstr(files_OfDirNUV[j].c_str(), "ra-decexp.fits") != NULL) {    
LOG(INFO) << "Reading " << files_OfDirNUV[j];
                                readKeywords((char*) files_OfDirNUV[j].c_str(), 2, 2, TSTRING, "FILTER", &filter_val, TINT, "WIN_X_SZ", &Win_Size);
                                if (strcmp(filter_val, curr_filtercal_val) != 0 || atoi(window_index_track[window_index].c_str()) != Win_Size) {
                                    LOG(INFO) << "NO filter match found ...Checking next file..";
                                    //                            Full_framePath.clear();
                                    //                            snrFilePath.clear();

                                    break;
                                }
                                eventList_RADEC.push_back(files_OfDirNUV[j]);
                            }

                            if (strstr(files_OfDirNUV[j].c_str(), "snr") != NULL) {

                                snrFilePath.push_back(files_OfDirNUV[j]);

                                //LOG(INFO) << files_OfDirNUV[j];
                            }

                        }



                    }
                    if (cnt_integratinTime != 0)
                        Final_integrationTime = Final_integrationTime / cnt_integratinTime;

                    if (Full_framePath.size() > 0) {

                        string temp_str, tempsnrstr, tempexpstr;
                        float temp_rollangle;
                        float temp_numStars;

                        for (int i = 0; i < numStars.size(); i++) {
                            for (int j = numStars.size() - 1; j > i; j--) {
                                if (numStars[j - 1] < numStars[j]) {
                                    swap1(numStars[j], numStars[j - 1]);
                                    swap1(Full_framePath[j], Full_framePath[j - 1]);
                                    swap1(RollAngles_track[j], RollAngles_track[j - 1]);
                                    swap1(snrFilePath[j], snrFilePath[j - 1]);
                                    swap1(ExpPath[j], ExpPath[j - 1]);
                                    swap1(eventList_RADEC[j], eventList_RADEC[j - 1]);
                                }
                            }
                        }
//                        for (int i = 1; i < numStars.size(); i++) {
//                            if (numStars[i] > numStars[i - 1]) {
//                                temp_numStars=numStars[i - 1];
//                                numStars[i - 1]=numStars[i];
//                                numStars[i]=temp_numStars;
//                                temp_str = Full_framePath[i - 1];
//                                Full_framePath[i - 1] = Full_framePath[i];
//                                Full_framePath[i] = temp_str;
//                                temp_rollangle = RollAngles_track[i - 1];
//                                RollAngles_track[i - 1] = RollAngles_track[i];
//                                RollAngles_track[i] = temp_rollangle;
//                                tempsnrstr = snrFilePath[i - 1];
//                                snrFilePath[i - 1] = snrFilePath[i];
//                                snrFilePath[i] = tempsnrstr;
//                                tempexpstr = ExpPath[i - 1];
//                                ExpPath[i - 1] = ExpPath[i];
//                                ExpPath[i] = tempexpstr;
//                            }
//                        }

                        for (int i = 0; i < numStars.size(); i++) {
                            LOG(INFO) << Full_framePath[i];
                            LOG(INFO) << snrFilePath[i];
                            LOG(INFO) << ExpPath[i];
                            LOG(INFO)<<eventList_RADEC[i];
                        }
                       
                        //Now we identified which frame to be used as a reference frame;

                        //            for (int i = 0; i < numStars.size(); i++) {
                        //                LOG(INFO) << Full_framePath[i];
                        //                LOG(INFO) << snrFilePath[i];
                        //
                        //
                        //
                        //            }
                        // exit(1);
                        //opening the frmaes.




                        char nameprfx[NAMESIZE];
                        int xsizeimage, ysizeimage;
                        double RA_pnt, DEC_pnt;
                        char date_obs[NAMESIZE], roll_angle[NAMESIZE];
                        //readKeywords((char*) Full_framePath[0].c_str(),1,7,TSTRING,"NAMEPRFX",nameprfx,TINT,"XSIZE",&xsizeimage,TINT,"YSIZE",&ysizeimage,TFLOAT,"RA_PNT",&RA_pnt
                        //     ,TFLOAT,"DEC_PNT",&DEC_pnt,TSTRING,"DATE",date_obs,TFLOAT,"ROLLAPPLIED",&roll_angle);
                        uvitObj.copyAllheaderKeys((char*) Full_framePath[0].c_str());
                        bool flag_shiftFound = TRUE;
                        int cnt_nuv = 0;
                        float *NewArraFlipped = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
                        fitsfile *fimgFile;


                        // First_HDU_headerTrack = uvitObj.header_info;
                        track_headerInfoL1.clear();
                        fitsfile *fin;
                        vector<string> L1headerInfo;
                        status = 0;
                        uvitObj.copyAllheaderKeys((char*) Full_framePath[0].c_str());
                        vector<string> FirstHDU_headertrack;
                        FirstHDU_headertrack = uvitObj.header_info;
                        string FinalStr;
                        int pos,pos1,pos2;
                        string ExtName,Orbit_Dir,temp_strForSubstr,tempsec_str;
                         Sum_Exp_Time=0.0;
                           vector<string> IndividualOrbit_Dir;
 float Max_value_Exp;
int cnter_elements=0;
vector<float> cent_X,cent_Y,cent_I;
double peakofexp=0.0f;
int indexx,indexy,cntpixels;
float Sum_ofpixels,valuetoCmpr;
RA_pnt=0.0f;DEC_pnt=0.0f;
                        for (int i = 0; i < Full_framePath.size(); i++) {
                           // 
LOG(INFO)<<"**FILE ->"<<Full_framePath[i];
                            Sum_ofpixels=0.0f;
                            flag_shiftFound = TRUE;
string Exppathvalue = ExpPath[i];
                            readImage((char*) Full_framePath[i].c_str(), 1, Frame_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc); //reading the image.
                            readImage((char*) ExpPath[i].c_str(), 1, Frame_Exposure_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc); //reading the  exposure image
                            readKeywords((char*)  Full_framePath[i].c_str(), 1, 2, TFLOAT,"EXP_TIME",&Exp_Time_Keyword,TDOUBLE,"PEAK_OF_EXP",&peakofexp);
                            Max_value_Exp=peakofexp;
                            if(i==0){
                               readKeywords((char*)  Full_framePath[i].c_str(), 1, 2, TDOUBLE,"CRVAL1",&RA_pnt,TDOUBLE,"CRVAL2",&DEC_pnt); 
                             }


//for (int i =0;i<uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc;i++)
//{
//if(Frame_Exposure_data[i]<(20*Max_value_Exp/100)){
//Frame_Exposure_data[i]=0.0f;
//}


//}



                            
                            pos = Full_framePath[i].find("/uvtN.");
                            pos1=Full_framePath[i].find("_NUV_");
                             tempsec_str=Full_framePath[i].substr(0,pos1);
                           
                             
                            ExtName = Full_framePath[i].substr(pos + 1, 7);
                            temp_strForSubstr=Full_framePath[i].substr(pos1, Full_framePath[i].length()-1);
                            pos2=temp_strForSubstr.find("/");
                          
                            Orbit_Dir=tempsec_str+Full_framePath[i].substr(pos1, pos2);
                            
                            fits_open_file(&fin, (char*) Full_framePath[i].c_str(), READONLY, &status);
                            printError(status, "Error in Opening the file", (char*) Full_framePath[i].c_str());
                            fits_movabs_hdu(fin, 2, NULL, &status);
                            printError(status, "Error in moving to 2nd HDU", (char*) Full_framePath[i].c_str());
                            copyUsrkeywrdsTovect(fin, L1headerInfo);
                            FinalStr = "EXTNAME = " + ExtName;
                              if(i==0){
                                current_ExtName="Reference_Orbit= "+Orbit_Dir+", "+ExtName;
                            }else{
                                current_ExtName="Matched_Orbit= "+Orbit_Dir+", "+ExtName;
                            }
                            L1headerInfo.push_back(FinalStr);
                            uvitObj.copyAllheaderKeys((char*) Full_framePath[i].c_str());
                            //                    track_headerInfoL1.push_back(uvitObj.header_info);
                            track_headerInfoL1.push_back(L1headerInfo);
                            fits_close_file(fin, &status);
                            //               initArray(NewArraFlipped,uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc,0.0f);               
                            //                
                            //                for (int i=0;i<uvitObj.IMG_DIM_FIpc;i++){
                            //                    for (int j=0;j<uvitObj.IMG_DIM_FIpc;j++){
                            //                        NewArraFlipped[i*uvitObj.IMG_DIM_FIpc+j]=Frame_data[i*uvitObj.IMG_DIM_FIpc+(uvitObj.IMG_DIM_FIpc-1-j)];
                            //                    }
                            //                }
                            //                initArray(Frame_data,uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc,0.0f); 
                            //                for (int i=0;i<uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc;i++){
                            //                    Frame_data[i]=NewArraFlipped[i];
                            //                    
                            //                }


                            uvitObj.sd_multi_factor_defaultpc = sd_multfactor;
                            uvitObj.refine_Winsizepc = refineWinSizedefault * 8;
                            uvitObj.centroid_Winsizepc = centWinSizedefault * 8;

				cntpixels=0;

                            if (i == 0) {//This is for reference frame.
                                status = uvitObj.findStar_algo1(Frame_data, PC); // for first frame taken as  reference
                                if (status) {
                                    // LOG (ERROR) << endl << "***Error in finding star algorithm 1 for frame  " << infile << "  ***" << endl ;
                                    flagRefFrameNotFound=TRUE;
				//return (EXIT_FAILURE);
				break;
                                }
                                //LOG(INFO) << uvitObj.Cx.size() << endl;

for (int i = 0; i < uvitObj.Cx.size(); i++) {
//LOG(INFO)<<uvitObj.Cx[i]<<" "<<uvitObj.Cy[i]<<" "<<Frame_Exposure_data[(int)(round(uvitObj.Cy[i])*uvitObj.IMG_DIM_FIpc+round(uvitObj.Cx[i]))];
cntpixels=0;Sum_ofpixels=0.0f;
indexx=round(uvitObj.Cx[i]);
indexy=round(uvitObj.Cy[i]);
	for (int j=indexx-12;j<=indexx+12;j++){
for (int k =indexy-12;k<=indexy+12;k++){
if((int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j)) >=0 && (int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))<uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc && Frame_Exposure_data[(int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))]!= INVALID_PIX_VALUE ){
Sum_ofpixels=Sum_ofpixels+Frame_Exposure_data[(int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))];
cntpixels++;
}
}
}
valuetoCmpr=Sum_ofpixels/cntpixels;
//LOG(INFO)<<"TT:"<<Max_value_Exp<<" "<<valuetoCmpr<<" "<<cntpixels;
if(valuetoCmpr>(20*Max_value_Exp/100)){
cx_ref.push_back(uvitObj.Cx[i]);
                                    cy_ref.push_back(uvitObj.Cy[i]);
                                    ci_ref.push_back(uvitObj.Ci[i]);
}

}

if(cx_ref.size()==0){
LOG(INFO)<<"***Reference frame star not found***";
flagRefFrameNotFound=TRUE;
break;
}
                              //  for (int i = 0; i < uvitObj.Cx.size(); i++) {
				//if(Frame_Exposure_data[(int)(round(uvitObj.Cy[i])*uvitObj.IMG_DIM_FIpc+round(uvitObj.Cx[i]))]!=0.0f){
                                 //   cx_ref.push_back(uvitObj.Cx[i]);
                                   // cy_ref.push_back(uvitObj.Cy[i]);
                                //    ci_ref.push_back(uvitObj.Ci[i]);
				//}

                                //}
                                LOG(INFO) << uvitObj.Cx.size() << endl;
                                uvitObj.Cx.clear();
                                uvitObj.Cy.clear();
                                uvitObj.Ci.clear();
                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                                    if (Frame_data[i] != INVALID_PIX_VALUE && Frame_Exposure_data[i] != INVALID_PIX_VALUE) {
                                        Frame_data_cum[i] = Frame_data_cum[i] + Frame_data[i] * Frame_Exposure_data[i];
                                        Frame_Expdata_cum[i] = Frame_Expdata_cum[i] + Frame_Exposure_data[i];


                                    }
                                }
				 LOG(INFO)<<"THE DATA->"<<Frame_Exposure_data[(uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc/2)+uvitObj.IMG_DIM_FIpc/2]<<" "<<Frame_Expdata_cum[(uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc/2)+uvitObj.IMG_DIM_FIpc/2]<<" "<<Exppathvalue;

                                LOG(INFO) << uvitObj.Cx.size() << endl;
                            } else {


				if(abs(RollAngles_track[i]-RollAngles_track[0])>2){
					LOG(ERROR)<<"***Roll Angle ISSUE***<<Ref Frame Roll->"<<RollAngles_track[0]<<", Current frame Roll->"<<RollAngles_track[i];
					continue;
				}


                                status = uvitObj.findStar_algo1(Frame_data, PC); // for first frame taken as  reference
                                if (status) {
                                     LOG (ERROR) << endl << "***Error in finding star     ***" << endl ;
                                    continue;
                                }
                                uvitObj.Rx.clear();
                                uvitObj.Ry.clear();
                                uvitObj.Rval.clear();
                                uvitObj.Fx.clear();
                                uvitObj.Fy.clear();
                                uvitObj.Fval.clear();
cent_X.clear();
cent_Y.clear();
cent_I.clear();
//LOG(INFO)<<uvitObj.Cx.size();
for (int i = 0; i < uvitObj.Cx.size(); i++) {
cntpixels=0;
Sum_ofpixels=0.0f;
indexx=round(uvitObj.Cx[i]);
indexy=round(uvitObj.Cy[i]);
	for (int j=indexx-12;j<=indexx+12;j++){
for (int k =indexy-12;k<=indexy+12;k++){
if((int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j)) >=0 && (int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))<uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc && Frame_Exposure_data[(int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))]!= INVALID_PIX_VALUE ){
Sum_ofpixels=Sum_ofpixels+Frame_Exposure_data[(int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))];
cntpixels++;
}
}
}
valuetoCmpr=Sum_ofpixels/cntpixels;
if(valuetoCmpr>(20*Max_value_Exp/100)){
cent_X.push_back(uvitObj.Cx[i]);
cent_Y.push_back(uvitObj.Cy[i]);
cent_I.push_back(uvitObj.Ci[i]);
}

}




				//	for(int i =0;i<uvitObj.Cx.size();i++){
//if(Frame_Exposure_data[(int)(round(uvitObj.Cy[i])*uvitObj.IMG_DIM_FIpc+round(uvitObj.Cx[i]))]!=0.0f){
    
 //cent_X.push_back(uvitObj.Cx[i]);
   //                                 cent_Y.push_back(uvitObj.Cy[i]);
     //                               cent_I.push_back(uvitObj.Ci[i]);

//}

//				}

uvitObj.Cx=cent_X;
uvitObj.Cy=cent_Y;
uvitObj.Ci=cent_I;
if(uvitObj.Cx.size()==0){
continue;
}

//LOG(INFO)<<uvitObj.Cx.size()<<" "<<Max_value_Exp;

                                //                    for (int i = 0; i < uvitObj.Cx.size(); i++) {
                                //                        LOG(INFO) << uvitObj.Cx[i] << " " << uvitObj.Cy[i];
                                //                    }
                                // LOG(INFO) << "=================";
                                //                    for (int i = 0; i < cx_ref.size(); i++) {
                                //                        LOG(INFO) << cx_ref[i] << " " << cy_ref[i];
                                //                    }
                                if (uvitObj.Cx.size() >= cx_ref.size()) {
                                    matching_points = cx_ref.size();
                                } else {
                                    matching_points = uvitObj.Cx.size();
                                }

                                flag_TestOnecheck=FALSE;
                                flag_TestTwocheck=FALSE;
                                x_ref_arr.clear();
                                y_ref_arr.clear();
                                x_arr.clear();
                                y_arr.clear();
                                temp_x_arr.clear();
                                temp_y_arr.clear();
                                //min_stars_match = (int) ((100 - PERCENTGE_ERR_ALLOWED) * matching_points / 100);
                                min_stars_match = 2;
                                cnt = 0;
                               // while (cnt < min_stars_match) {

//                                    if (uvitObj.diff_Distpc > 256) {
//                                        //LOG(ERROR) << "No matches found...";
//                                        //return(EXIT_FAILURE);
//                                        flag_shiftFound = FALSE;
//                                        break;
//                                    }
//                                    x_ref_arr.clear();
//                                    y_ref_arr.clear();
//                                    x_arr.clear();
//                                    y_arr.clear();
//                                    temp_x_arr.clear();
//                                    temp_y_arr.clear();

                                    uvitObj.matchStars(cx_ref.size(), uvitObj.Cx.size(), 8, cx_ref.data(), cy_ref.data(), uvitObj.Cx.data(), uvitObj.Cy.data(),x_dx,y_dy,flag_TestOnecheck,flag_TestTwocheck);
                                      if(flag_TestOnecheck==FALSE && flag_TestTwocheck==FALSE){
                                         LOG(INFO)<<"NOT FOUND....";
                                         continue;
                                     }
//else{


//} 
                                    
                                    
                                  //  uvitObj.diff_Distpc = uvitObj.diff_Distpc * 2;

                                //}
//                                if (flag_shiftFound == FALSE) {
//                                    continue;
//                                }
//                                uvitObj.diff_Distpc = original_diff_Distpc;
//                                New_X_ref.clear(), New_Y_ref.clear(), New_x_arr.clear(), New_y_arr.clear(), New_Xdiff.clear(), New_Ydiff.clear();
//                                int_new_one.clear(), int_new_two.clear();
//                                if (cnt > 2) {
//                                    status = uvitObj.removeRecords(x_ref_arr, y_ref_arr, x_arr, y_arr, temp_x_arr, temp_y_arr, ci_ref.data(), uvitObj.Ci.data(), New_X_ref, New_Y_ref, New_x_arr, New_y_arr, New_Xdiff, New_Ydiff, int_new_one, int_new_two);
//                                    if (status) {
//                                        LOG(ERROR) << "Error in removing the records above mean+sigma";
//                                        return (EXIT_FAILURE);
//                                    }
//                                } else {
//                                    New_X_ref = x_ref_arr;
//                                    New_Y_ref = y_ref_arr;
//                                    New_x_arr = x_arr;
//                                    New_y_arr = y_arr;
//                                    New_Xdiff = temp_x_arr;
//                                    New_Ydiff = temp_y_arr;
//                                }
//
//                                for (int i = 0; i < New_X_ref.size(); i++) {
//
//                                    LOG(INFO) << New_X_ref[i] << " " << New_Y_ref[i] << " " << New_x_arr[i] << " " << New_y_arr[i] << " " << New_Xdiff[i] << " " << New_Ydiff[i] << endl;
//
//                                }
//
//
//                                //                 status = uvitObj.findShiftsNtheta (x_ref_arr.size () , x_ref_arr , y_ref_arr  ,x_arr , y_arr , New_Xdiff , New_Ydiff , 1,x_dx , y_dy , theta_dt) ;
//                                //        if (status)
//                                //        {
//                                //            LOG (INFO) << "Error in finding shifts n theta " << endl ;
//                                //            return (EXIT_FAILURE) ;
//                                //        } 
//                                //
//
//                                status = uvitObj.findShiftsNtheta(New_X_ref.size(), New_X_ref, New_Y_ref, New_x_arr, New_y_arr, New_Xdiff, New_Ydiff, 1, x_dx, y_dy, theta_dt);
//                                if (status) {
//                                    LOG(INFO) << "Error in finding shifts n theta " << endl;
//                                    return (EXIT_FAILURE);
//                                }



                                //LOG(INFO) << x_dx << " " << y_dy << " " << theta_dt;

                                
                                
                                //NOW OPENING EVENT LIST.
                                
                                 fits_open_file (&fptr_snr , eventList_RADEC[i].c_str() , READWRITE , &status) ;
        printError (status , " Error in opening file " , (char*)eventList_RADEC[i].c_str()) ;
        fits_movabs_hdu (fptr_snr , 2 , NULL , &status) ;
        printError (status , "Error in moving to HDU " ,(char*) eventList_RADEC[i].c_str()) ;
        fits_get_num_rows (fptr_snr , &nrows_snr , &status) ;
        printError (status , "Error in reading the number of rows" ,(char*) eventList_RADEC[i].c_str()) ;
        X_loc_snr= new  float[nrows_snr];
        Y_loc_snr= new float[nrows_snr];
        bad_flag_snr = new float[nrows_snr];
        mult_phtn_snr= new float[nrows_snr];
        enp_snr= new float[nrows_snr];
        fits_read_col (fptr_snr , TFLOAT , 13 , 1 , 1 , nrows_snr  , NULL , Y_loc_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" ,(char*) eventList_RADEC[i].c_str()) ;
        fits_read_col (fptr_snr , TFLOAT , 12 , 1 , 1 , nrows_snr  , NULL , X_loc_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" ,(char*) eventList_RADEC[i].c_str()) ;
        fits_read_col (fptr_snr , TFLOAT , 9 , 1 , 1 , nrows_snr  , NULL , bad_flag_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" ,(char*) eventList_RADEC[i].c_str()) ;
        fits_read_col (fptr_snr , TFLOAT , 10 , 1 , 1 , nrows_snr  , NULL , mult_phtn_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , (char*)eventList_RADEC[i].c_str()) ;
         fits_read_col (fptr_snr , TFLOAT , 11 , 1 , 1 , nrows_snr  , NULL , enp_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , (char*)eventList_RADEC[i].c_str()) ;
        fits_close_file (fptr_snr , &status) ;
        printError (status , "Error closing the file " , (char*)eventList_RADEC[i].c_str()) ;
                                
                                //TILL THIS

                                ctheta = 0.0f, stheta = 0.0f;
                                ctheta = cos(-1.0 * theta_dt);
                                stheta = sin(-1.0 * theta_dt);
                                //   ctheta = cos( 1.0* theta_dt);
                                //                    stheta = sin(1.0 * theta_dt);
                                LOG(INFO) << "Loop Started for assigning the correction to the frames..";

                               // initArray(FinalArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                  initArray(FinalArray, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc , 0.0f);
                                initArray(counterArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                initArray(counterArray_exp, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                initArray(Final_ExpArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                //rescale image from 4800 to 9600 for bi-linear interpolation .
                                initArray(SubdividedImage, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                initArray(Subdivided_expImage, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                //scale image to 9600.
                                //  LOG(INFO)<<"111";
                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                                    //if (Frame_data[i] != INVALID_PIX_VALUE)
                                      //  Frame_data[i] = Frame_data[i] / 4;
                                    if (Frame_Exposure_data[i] != INVALID_PIX_VALUE)
                                        Frame_Exposure_data[i] = Frame_Exposure_data[i] / 4;
                                }
//                                performSubDivisionIM(Frame_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, SubdividedImage, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2);
                                 performSubDivisionIM(Frame_Exposure_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, Subdivided_expImage, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2);
//                                //LOG(INFO)<<"222";
//                                int tempk1, tempk2;
//                                //                                double tempX1,tempX2;
//                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * 2; i++) {
//                                    mid_X_new = i - uvitObj.IMG_DIM_FIpc;
//                                    for (int j = 0; j < uvitObj.IMG_DIM_FIpc * 2; j++) {
//                                        mid_Y_new = j - uvitObj.IMG_DIM_FIpc;
//
//                                        //              x1=((mid_X_new * (ctheta)) - (mid_Y_new * (stheta)))+uvitObj.IMG_DIM_FIpc- (x_dx * 8);
//                                        //              x2=((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) +uvitObj.IMG_DIM_FIpc- (y_dy * 8);  
//                                        x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc - (x_dx * 8 * 2);
//                                        x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc - (y_dy * 8 * 2);
//                                        //              if(x1>0 && x1<FINALFRAMESIZE_REGAVG && x2>0 && x2<FINALFRAMESIZE_REGAVG){
//                                        //                  imageLocation_Array_X[(int)(round(x2)*FINALFRAMESIZE_REGAVG+round(x1))]=Regavg_subSampled_Array_Sig[j*FINALFRAMESIZE_REGAVG+i];
//                                        //              }
//                                        //                                        if((int)(round(x2)*uvitObj.IMG_DIM_FIpc*2+round(x1))>0 && (int)(round(x2)*uvitObj.IMG_DIM_FIpc*2+round(x1))<uvitObj.IMG_DIM_FIpc*2*uvitObj.IMG_DIM_FIpc*2)
//                                        if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc * 2 && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc * 2) {
//                                            if (SubdividedImage[j * uvitObj.IMG_DIM_FIpc * 2 + i] != INVALID_PIX_VALUE && FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] != INVALID_PIX_VALUE) {
//                                                FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + SubdividedImage[j * uvitObj.IMG_DIM_FIpc * 2 + i];
//                                                counterArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = counterArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + 1;
//                                            } else {
//                                                FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = INVALID_PIX_VALUE;
//                                            }
//                                        }
//                                        //                                        
//
//                                    }
//
//
//
//                                }
                                
                                
                                //added
                                for (int i =0;i <nrows_snr;i++){
                                    mid_X_new=X_loc_snr[i]-uvitObj.IMG_DIM_FIpc/2;
                                    mid_Y_new=Y_loc_snr[i]-uvitObj.IMG_DIM_FIpc/2;
                                     x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc/2 - (x_dx * 8);
                                     x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc/2 - (y_dy * 8 );
                                 
                                     if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc  && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc )
                                     FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc  + round(x1))]=FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc  + round(x1))]+bad_flag_snr[i]*mult_phtn_snr[i]*enp_snr[i];
                                     
                                     
                                }
                                

                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * 2; i++) {
                                    mid_X_new = i - uvitObj.IMG_DIM_FIpc;
                                    for (int j = 0; j < uvitObj.IMG_DIM_FIpc * 2; j++) {
                                        mid_Y_new = j - uvitObj.IMG_DIM_FIpc;

                                        x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc - (x_dx * 8 * 2);
                                        x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc - (y_dy * 8 * 2);
                                        //              if(x1>0 && x1<FINALFRAMESIZE_REGAVG && x2>0 && x2<FINALFRAMESIZE_REGAVG){
                                        //                  imageLocation_Array_X[(int)(round(x2)*FINALFRAMESIZE_REGAVG+round(x1))]=Regavg_subSampled_Array_Sig[j*FINALFRAMESIZE_REGAVG+i];
                                        //              }
                                        //if((int)(round(x2)*uvitObj.IMG_DIM_FIpc*2+round(x1))>0 && (int)(round(x2)*uvitObj.IMG_DIM_FIpc*2+round(x1))<uvitObj.IMG_DIM_FIpc*2*uvitObj.IMG_DIM_FIpc*2)
                                        if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc * 2 && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc * 2) {
                                            if (Subdivided_expImage[j * uvitObj.IMG_DIM_FIpc * 2 + i] != INVALID_PIX_VALUE && Final_ExpArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] != INVALID_PIX_VALUE) {
                                                Final_ExpArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = Final_ExpArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + Subdivided_expImage[j * uvitObj.IMG_DIM_FIpc * 2 + i];
                                                counterArray_exp[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = counterArray_exp[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + 1;
                                            } else {
                                                Final_ExpArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = INVALID_PIX_VALUE;
                                            }
                                        }
                                    }
                                }
                                //  LOG(INFO)<<"333";
                                initArray(Frame_FinalArray, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                                // LOG(INFO)<<"444";
                           //     ApplyBinning(FinalArray, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2, Frame_FinalArray, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, counterArray);
                                // LOG(INFO)<<"555";

                                initArray(Frame_FinalExpArray, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                                ApplyBinning(Final_ExpArray, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2, Frame_FinalExpArray, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, counterArray_exp);
                                //                    for (int index = 0; index < uvitObj.IMG_DIM_FIpc; index++) {
                                //                        x_index = index - uvitObj.IMG_DIM_FIpc / 2;
                                //                        for (int jindex = 0; jindex < uvitObj.IMG_DIM_FIpc; jindex++) {
                                //                            if (Frame_data[jindex * uvitObj.IMG_DIM_FIpc + index] != INVALID_PIX_VALUE) {
                                //
                                //                                y_index = jindex - uvitObj.IMG_DIM_FIpc / 2;
                                //                                //ROUNDING:Correction needed.previously rounding was only on the subset of full equation.
                                //                                new_index_x = round((x_index) * ctheta - (y_index) * stheta + uvitObj.IMG_DIM_FIpc / 2 - (x_dx * 8)); //new index x
                                //                                new_index_y = round((x_index) * stheta + (y_index) * ctheta + uvitObj.IMG_DIM_FIpc / 2 - (y_dy * 8)); //new index y
                                //
                                //                                if (round(new_index_x) < uvitObj.IMG_DIM_FIpc && round(new_index_x) > 0 && round(new_index_y) > 0 && round(new_index_y) < uvitObj.IMG_DIM_FIpc) {
                                //                                    // cnt_loop ++ ;
                                //                                    //Rounding:No correction needed.As rounding is applied on each direction i.e X and Y.
                                //                                    FinalArray[(int) (round(new_index_y) * uvitObj.IMG_DIM_FIpc + round(new_index_x))] = Frame_data[jindex * uvitObj.IMG_DIM_FIpc + index];
                                //
                                //                                }
                                //                            }
                                //
                                //                        }
                                //                    }

//                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
//                                    if (Frame_FinalArray[i] != INVALID_PIX_VALUE && Frame_FinalExpArray[i] != INVALID_PIX_VALUE && Frame_data_cum[i] != INVALID_PIX_VALUE && Frame_Expdata_cum[i] != INVALID_PIX_VALUE) {
//                                        Frame_data_cum[i] = Frame_data_cum[i] + Frame_FinalArray[i] * Frame_FinalExpArray[i];
//                                        Frame_Expdata_cum[i] = Frame_Expdata_cum[i] + Frame_FinalExpArray[i];
//                                    } else {
//                                        Frame_data_cum[i] = INVALID_PIX_VALUE;
//                                        Frame_Expdata_cum[i] = INVALID_PIX_VALUE;
//                                    }
//                                    //                                    if(Frame_data_cum[i]<0){
//                                    //                                        LOG(INFO)<< Frame_FinalArray[i]<<" "<<Frame_FinalExpArray[i]<<" "<<Frame_data_cum[i];
//                                    //                                        exit(1);
//                                    //                                    }
//                                }
                                  for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                                    if (FinalArray[i] != INVALID_PIX_VALUE && Frame_FinalExpArray[i] != INVALID_PIX_VALUE && Frame_data_cum[i] != INVALID_PIX_VALUE && Frame_Expdata_cum[i] != INVALID_PIX_VALUE) {
                                        Frame_data_cum[i] = Frame_data_cum[i] + FinalArray[i];// * Frame_FinalExpArray[i];
                                        Frame_Expdata_cum[i] = Frame_Expdata_cum[i] + Frame_FinalExpArray[i];
//if(i==(uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc/2)+uvitObj.IMG_DIM_FIpc/2) LOG(INFO)<<"THE DATA->"<<Frame_FinalExpArray[i]<<" "<<Frame_data_cum[i]<<" "<<Exppathvalue;
                                    } else {
                                        Frame_data_cum[i] = INVALID_PIX_VALUE;
                                        Frame_Expdata_cum[i] = INVALID_PIX_VALUE;
                                    }
                                    //                                    if(Frame_data_cum[i]<0){
                                    //                                        LOG(INFO)<< Frame_FinalArray[i]<<" "<<Frame_FinalExpArray[i]<<" "<<Frame_data_cum[i];
                                    //                                        exit(1);
                                    //                                    }
                                }
 LOG(INFO)<<"THE DATA->"<<Frame_FinalExpArray[(uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc/2)+uvitObj.IMG_DIM_FIpc/2]<<" "<<Frame_Expdata_cum[(uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc/2)+uvitObj.IMG_DIM_FIpc/2]<<" "<<Exppathvalue;
                                
                                
                                uvitObj.Cx.clear();
                                uvitObj.Cy.clear();
                                uvitObj.Ci.clear();
                            }
                            cnt_nuv++;
				LOG(INFO)<<"File name"<<Full_framePath[i];
                                 Sum_Exp_Time=Sum_Exp_Time+Exp_Time_Keyword;
                                  IndividualOrbit_Dir.push_back(current_ExtName);
                                
                        }
if(flagRefFrameNotFound==TRUE){
LOG(INFO)<<"***NOt going Further ,reference frame Star not found***";
continue;
}
                       
                        //float tempFrameCumm = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
                        //           for (int i = 0; i < uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc; i++) {
                        //               tempFrameCumm[i]=Frame_data_cum[i];
                        //           }
                        float *tempArrayXinvert = new float [uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];

                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc; i++) {
                            for (int j = 0; j < uvitObj.IMG_DIM_FIpc; j++) {
					tempArrayXinvert[i * uvitObj.IMG_DIM_FIpc + j]=0.0f;
                                if (Frame_Expdata_cum[i * uvitObj.IMG_DIM_FIpc + (j)] != 0) {
                                    tempArrayXinvert[i * uvitObj.IMG_DIM_FIpc + j] = Frame_data_cum[i * uvitObj.IMG_DIM_FIpc + ( j)] / (Frame_Expdata_cum[i * uvitObj.IMG_DIM_FIpc + (j)]);
                                    //tempArrayXinvert[i * uvitObj.IMG_DIM_FIpc + j] = Frame_data_cum[i * uvitObj.IMG_DIM_FIpc + (j)] / Frame_Expdata_cum[i * uvitObj.IMG_DIM_FIpc + (j)];
//                                    if (tempArrayXinvert[i * uvitObj.IMG_DIM_FIpc + j] < 0) {
  //                                      LOG(INFO) << Frame_data_cum[i * uvitObj.IMG_DIM_FIpc + (uvitObj.IMG_DIM_FIpc - 1 - j)] << " " << Frame_Expdata_cum[i * uvitObj.IMG_DIM_FIpc + (uvitObj.IMG_DIM_FIpc - 1 - j)];
  //                                  }
                                }
                            }


                            //tempArrayXinvert[i * uvitObj.IMG_DIM_FIpc + j] = Frame_data_cum[i * uvitObj.IMG_DIM_FIpc + (j)]/cnt_nuv;}
                        }

                        //            vector<float> Xarr, Yarr, intarr;
                        //            for (int i = 0; i < uvitObj.IMG_DIM_FIpc; i++) {
                        //                for (int j = 0; j < uvitObj.IMG_DIM_FIpc; j++) {
                        //                    if (Frame_data_cum[j * uvitObj.IMG_DIM_FIpc + i] > 0) {
                        //                        Xarr.push_back(i);
                        //                        Yarr.push_back(j);
                        //                        intarr.push_back(Frame_data_cum[j * uvitObj.IMG_DIM_FIpc + i]);
                        //                        //                        LOG(INFO)<<i<<" "<<j<<" "<<Frame_data_cum[j*uvitObj.IMG_DIM_FIpc+i];
                        //
                        //                    }
                        //
                        //                }
                        //            }


                        //now convert X-Y image to RA_DEC by applying to Roll angle.

                        initArray(RADECimage, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                        //            ctheta=cos(RollAngles_track[0]*M_PI/180);
                        //            stheta=sin(RollAngles_track[0]*M_PI/180);
                        //reading shift And rotationFile.

                        LOG(INFO) << "Snr file->" << snrFilePath[0];
                        //            fits_open_file(&fptr_snr, snrFilePath[0].c_str(), READWRITE, &status);
                        //            printError(status, " Error in opening file ", (char*) snrFilePath[0].c_str());
                        //            fits_movabs_hdu(fptr_snr, 2, NULL, &status);
                        //            printError(status, "Error in moving to HDU ", (char*) snrFilePath[0].c_str());
                        //            fits_get_num_rows(fptr_snr, &nrows_snr, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            X_loc_snr = new float[nrows_snr];
                        //            Y_loc_snr = new float[nrows_snr];
                        //            bad_flag_snr = new float[nrows_snr];
                        //            mult_phtn_snr = new float[nrows_snr];
                        //            enp_snr = new float[nrows_snr];
                        //            fits_read_col(fptr_snr, TFLOAT, 5, 1, 1, nrows_snr, NULL, Y_loc_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_read_col(fptr_snr, TFLOAT, 4, 1, 1, nrows_snr, NULL, X_loc_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_read_col(fptr_snr, TFLOAT, 8, 1, 1, nrows_snr, NULL, bad_flag_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_read_col(fptr_snr, TFLOAT, 9, 1, 1, nrows_snr, NULL, mult_phtn_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_read_col(fptr_snr, TFLOAT, 10, 1, 1, nrows_snr, NULL, enp_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_close_file(fptr_snr, &status);
                        //            printError(status, "Error closing the file ", (char*) snrFilePath[0].c_str());



                        //             ctheta=cos((360-RollAngles_track[0])*M_PI/180);
                        //            stheta=sin((360-RollAngles_track[0])*M_PI/180);
                        //                        ctheta = cos((RollAngles_track[0]) * M_PI / 180);
                        //                        stheta = sin((RollAngles_track[0]) * M_PI / 180);
                        ctheta = cos(-(RollAngles_track[0]) * M_PI / 180);
                        stheta = sin(-(RollAngles_track[0]) * M_PI / 180);
                        //            for (int i = 0; i < nrows_snr; i++) 
                        //            {
                        //              x1=(((X_loc_snr[i]-uvitObj.IMG_DIM_FIpc/2) * (ctheta)) - ((Y_loc_snr[i]-uvitObj.IMG_DIM_FIpc/2) * (stheta)))+uvitObj.IMG_DIM_FIpc/2;
                        //              x2=(((X_loc_snr[i]-uvitObj.IMG_DIM_FIpc/2) * (stheta)) +((Y_loc_snr[i]-uvitObj.IMG_DIM_FIpc/2)* (ctheta))) +uvitObj.IMG_DIM_FIpc/2;  
                        //
                        //              if((int)((int)(round(x2))*uvitObj.IMG_DIM_FIpc+(int)(round(x1)))<uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc){ 
                        //              //Rotated_image[(int)(FINALFRAMESIZE_REGAVG-round(x2))*FINALFRAMESIZE_REGAVG+(int)(round(x1))]=Rotated_image[(int)(FINALFRAMESIZE_REGAVG-round(x2))*FINALFRAMESIZE_REGAVG+(int)(round(x1))]+effective_NumPhotons[i]*mult_temp[i]*badFlag_temp[i];
                        //                 //Rotated_image[(int)(round(x1))*FINALFRAMESIZE_REGAVG+(int)(FINALFRAMESIZE_REGAVG-round(x2))]=Rotated_image[(int)(round(x1))*FINALFRAMESIZE_REGAVG+(int)(FINALFRAMESIZE_REGAVG-round(x2))]+effective_NumPhotons[i]*mult_temp[i]*badFlag_temp[i];
                        //            RADECimage[(int)(round(x2))*uvitObj.IMG_DIM_FIpc+(int)(round(x1))]=RADECimage[(int)(round(x2))*uvitObj.IMG_DIM_FIpc+(int)(round(x1))]+enp_snr[i]*mult_phtn_snr[i]*bad_flag_snr[i];
                        //              }
                        //            }


                        //bilinear interpolation
                        initArray(counterArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        initArray(counterArray_exp, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        float *tempArrSubDivided = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
                        initArray(tempArrSubDivided, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        float *RADECimageSubdivided = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
                        initArray(RADECimageSubdivided, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        float *tempExpArrSubDivided = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
                        initArray(tempExpArrSubDivided, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        float *RADECExpimageSubdivided = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
                        initArray(RADECExpimageSubdivided, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                       // for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                         //   if (tempArrayXinvert[i] != INVALID_PIX_VALUE) {
                           //     tempArrayXinvert[i] = tempArrayXinvert[i] / 4;
                           // }
                           // if (Frame_Expdata_cum[i] != INVALID_PIX_VALUE) {
                             //   Frame_Expdata_cum[i] = Frame_Expdata_cum[i] / 4;
                           // }
                       // }
                        performSubDivisionIM(tempArrayXinvert, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, tempArrSubDivided, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2);
                        performSubDivisionIM(Frame_Expdata_cum, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, tempExpArrSubDivided, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2); //for exposure Array

                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc * 2; i++) {
                            mid_X_new = i - uvitObj.IMG_DIM_FIpc;
                            for (int j = 0; j < uvitObj.IMG_DIM_FIpc * 2; j++) {
                                mid_Y_new = j - uvitObj.IMG_DIM_FIpc;

                                //                                x1 = (-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc;
                                //                                x2 = (-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc;
                                x1 = (-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc;
                                x2 = (-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc;
                                //              if(x1>0 && x1<FINALFRAMESIZE_REGAVG && x2>0 && x2<FINALFRAMESIZE_REGAVG){
                                //                  imageLocation_Array_X[(int)(round(x2)*FINALFRAMESIZE_REGAVG+round(x1))]=Regavg_subSampled_Array_Sig[j*FINALFRAMESIZE_REGAVG+i];
                                //              }
                                //                                if(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1)>0 && round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1)<uvitObj.IMG_DIM_FIpc * 2*uvitObj.IMG_DIM_FIpc * 2)
                                if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc * 2 && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc * 2) {
                                    if (tempArrSubDivided[j * uvitObj.IMG_DIM_FIpc * 2 + i] != INVALID_PIX_VALUE && RADECimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] != INVALID_PIX_VALUE) {
                                        RADECimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = RADECimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + tempArrSubDivided[j * uvitObj.IMG_DIM_FIpc * 2 + i];
                                        counterArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = counterArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + 1;
                                    } else {
                                        RADECimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = INVALID_PIX_VALUE;
                                    }
                                }
                            }
                        }

                        //for exposure Array 

                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc * 2; i++) {
                            mid_X_new = i - uvitObj.IMG_DIM_FIpc;
                            for (int j = 0; j < uvitObj.IMG_DIM_FIpc * 2; j++) {
                                mid_Y_new = j - uvitObj.IMG_DIM_FIpc;

                                //                                x1 = (-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc;
                                //                                x2 = (-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc;
                                x1 = (-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc;
                                x2 = (-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc;
                                //              if(x1>0 && x1<FINALFRAMESIZE_REGAVG && x2>0 && x2<FINALFRAMESIZE_REGAVG){
                                //                  imageLocation_Array_X[(int)(round(x2)*FINALFRAMESIZE_REGAVG+round(x1))]=Regavg_subSampled_Array_Sig[j*FINALFRAMESIZE_REGAVG+i];
                                //              }

                                //                             if(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1)>0 && round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1)<uvitObj.IMG_DIM_FIpc * 2*uvitObj.IMG_DIM_FIpc * 2)
                                if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc * 2 && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc * 2) {
                                    if (tempExpArrSubDivided[j * uvitObj.IMG_DIM_FIpc * 2 + i] != INVALID_PIX_VALUE && RADECExpimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] != INVALID_PIX_VALUE) {
                                        RADECExpimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = RADECExpimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + tempExpArrSubDivided[j * uvitObj.IMG_DIM_FIpc * 2 + i];
                                        counterArray_exp[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = counterArray_exp[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + 1;
                                    } else {
                                        RADECExpimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = INVALID_PIX_VALUE;
                                    }
                                }


                            }



                        }

                        //float FinalArrauSubSampled = new float[uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc];
                        ApplyBinning(RADECimageSubdivided, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2, RADECimage, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, counterArray);
                        ApplyBinning(RADECExpimageSubdivided, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2, RADECExpimage, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, counterArray_exp);

                       
                        delete[] tempArrSubDivided;
                        delete[] tempExpArrSubDivided;
                        delete[] RADECimageSubdivided;
                        delete[] RADECExpimageSubdivided;

                        float *noise_Map = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
                        initArray(noise_Map, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                            if (Frame_Expdata_cum[i] * Final_integrationTime != 0.0f && tempArrayXinvert[i] != INVALID_PIX_VALUE && Frame_Expdata_cum[i] != INVALID_PIX_VALUE)
                         //                               noise_Map[i] = sqrt(RADECimage[i]) / (RADECExpimage[i] * Final_integrationTime);
                                noise_Map[i] = sqrt(tempArrayXinvert[i] * Frame_Expdata_cum[i] * Final_integrationTime) / (Frame_Expdata_cum[i] * Final_integrationTime);
                            else {
                                noise_Map[i] = 0.0f;
                            }
                        }

                        status = 0;

                        naxes[0] = naxes[1] = uvitObj.IMG_DIM_FIpc;
                        naxis = 2;
                        bitpix = FLOAT_IMG;


                        finalImagePath = level2_Origname + "NUV_Final_" + curr_filtercal_val + "_W" + (string) window_index_track[window_index];


                        LOG(INFO) << finalImagePath;
                        if (DirExists((char*) finalImagePath.c_str())) {
                            LOG(ERROR) << "Directory exists and clobber=yes";
                            cmd = (string) "rm -rf " + (string) finalImagePath;
                            system(cmd.c_str());
                        }

                        cmd = "mkdir -p " + (string) finalImagePath;
                        system(cmd.c_str());

                        sprintf(finalfilename, "%s/%s_W%s_%s", (char*) finalImagePath.c_str(),curr_filtercal_val,(char*)window_index_track[window_index].c_str(),"FinalImage_Sig.fits");
                        sprintf(finalExpfilename, "%s/%s_W%s_%s", (char*) finalImagePath.c_str(),curr_filtercal_val,(char*)window_index_track[window_index].c_str(), "FinalImage_Exp.fits");
                        sprintf(finalNoiseMapfilename, "%s/%s_W%s_%s", (char*) finalImagePath.c_str(),curr_filtercal_val,(char*)window_index_track[window_index].c_str(), "FinalImage_NoiseMap.fits");

                        fits_create_file(&fout, finalfilename, &status);
                        printError(status, "Error in creating the output Signal File", outfile);
                        fits_create_img(fout, bitpix, naxis, naxes, &status);
                        printError(status, "Error in Creating the image for Signal Fie", outfile);
                        fits_write_pix(fout, TFLOAT, fpixel, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, tempArrayXinvert, &status);
                        printError(status, "***Error in writing the pixels to output***", outfile);
                       // LOG(INFO) << track_headerInfoL1.size();
			string minNoStars_str,threshodpc_str;
			threshodpc_str="Threshold Used for DriverModule-multiplier to Sigma ="+(string)convertFloatToStr(uvitObj.sd_multi_factor_defaultpc);			
			minNoStars_str="Minimum number of star requirement in finding star in DriverModule="+(string)convertIntToStr(uvitObj.min_num_stars);
			 fits_write_history (fout, threshodpc_str.c_str(), &status);
 			fits_write_history (fout, minNoStars_str.c_str(), &status);
                        for (int i = 0; i < track_headerInfoL1.size(); i++) {
                            naxes[0] = naxes[1] = 0;
                            if (i > 0) {
                                fits_open_file(&fout, finalfilename, READWRITE, &status);
                                printError(status, "Error in Creating the image for Signal Fie", outfile);
                                fits_movabs_hdu(fout, i + 1, NULL, &status); //previous HDU
                                printError(status, "Error in Creating the image for Signal Fie1", outfile);
                            }

                            fits_create_img(fout, bitpix, naxis, naxes, &status);
                            printError(status, "Error in Creating the image for Signal Fie", outfile);
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                            fits_open_file(&fout, finalfilename, READWRITE, &status);
                            printError(status, "Error in Creating the image for Signal Fie3", outfile);
                            fits_movabs_hdu(fout, i + 2, NULL, &status);
                            for (int j = 0; j < track_headerInfoL1[i].size(); j++) {
                                if (strstr(track_headerInfoL1[i][j].c_str(), "NAXIS") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "BITPIX") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "NAXIS1") == NULL
                                        && strstr(track_headerInfoL1[i][j].c_str(), "NAXES2") == NULL) {
                                    fits_write_record(fout, track_headerInfoL1[i][j].c_str(), &status);
                                    printError(status, "Error in Creating the image for Signal Fie4", outfile);
                                }
                            }
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                        }


                        //            fits_close_file(fout, &status);
                        //            printError(status, "Error in closing the  output Signal fits file", outfile);
                        naxes[0] = naxes[1] = uvitObj.IMG_DIM_FIpc;
                        fits_create_file(&fout, finalExpfilename, &status);
                        printError(status, "Error in creating the output Signal File", outfile);
                        fits_create_img(fout, bitpix, naxis, naxes, &status);
                        printError(status, "Error in Creating the image for Signal Fie", outfile);
                        fits_write_pix(fout, TFLOAT, fpixel, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, Frame_Expdata_cum, &status);
                        printError(status, "***Error in writing the pixels to output***", outfile);
			 fits_write_history (fout, threshodpc_str.c_str(), &status);
 			fits_write_history (fout, minNoStars_str.c_str(), &status);
                        for (int i = 0; i < track_headerInfoL1.size(); i++) {
                            naxes[0] = naxes[1] = 0;
                            if (i > 0) {
                                fits_open_file(&fout, finalExpfilename, READWRITE, &status);
                                fits_movabs_hdu(fout, i + 1, NULL, &status); //previous HDU
                            }
                            fits_create_img(fout, bitpix, naxis, naxes, &status);
                            printError(status, "Error in Creating the image for Signal Fie", outfile);
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                            fits_open_file(&fout, finalExpfilename, READWRITE, &status);
                            fits_movabs_hdu(fout, i + 2, NULL, &status);
                            for (int j = 0; j < track_headerInfoL1[i].size(); j++) {
                                if (strstr(track_headerInfoL1[i][j].c_str(), "NAXIS") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "BITPIX") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "NAXIS1") == NULL
                                        && strstr(track_headerInfoL1[i][j].c_str(), "NAXES2") == NULL) {
                                    fits_write_record(fout, track_headerInfoL1[i][j].c_str(), &status);
                                }
                            }
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                        }
                        //fits_close_file(fout, &status);
                        //printError(status, "Error in closing the  output Signal fits file", outfile);
                        naxes[0] = naxes[1] = uvitObj.IMG_DIM_FIpc;
                        fits_create_file(&fout, finalNoiseMapfilename, &status);
                        printError(status, "Error in creating the output Signal File", outfile);
                        fits_create_img(fout, bitpix, naxis, naxes, &status);
                        printError(status, "Error in Creating the image for Signal Fie", outfile);
                        fits_write_pix(fout, TFLOAT, fpixel, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, noise_Map, &status);
                        printError(status, "***Error in writing the pixels to output***", outfile);
			 fits_write_history (fout, threshodpc_str.c_str(), &status);
 			fits_write_history (fout, minNoStars_str.c_str(), &status);
                        for (int i = 0; i < track_headerInfoL1.size(); i++) {
                            naxes[0] = naxes[1] = 0;
                            if (i > 0) {
                                fits_open_file(&fout, finalNoiseMapfilename, READWRITE, &status);
                                fits_movabs_hdu(fout, i + 1, NULL, &status); //previous HDU
                            }
                            fits_create_img(fout, bitpix, naxis, naxes, &status);
                            printError(status, "Error in Creating the image for Signal Fie", outfile);
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                            fits_open_file(&fout, finalNoiseMapfilename, READWRITE, &status);
                            fits_movabs_hdu(fout, i + 2, NULL, &status);
                            for (int j = 0; j < track_headerInfoL1[i].size(); j++) {
                                if (strstr(track_headerInfoL1[i][j].c_str(), "NAXIS") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "BITPIX") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "NAXIS1") == NULL
                                        && strstr(track_headerInfoL1[i][j].c_str(), "NAXES2") == NULL) {
                                    fits_write_record(fout, track_headerInfoL1[i][j].c_str(), &status);
                                }
                            }
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                        }

                        //            fits_close_file(fout, &status);
                        //            printError(status, "Error in closing the  output Signal fits file", outfile); 



                        delete[] RADECimage;
                        delete[] noise_Map;
if(Sum_Exp_Time <0) Sum_Exp_Time=INVALID_PIX_VALUE;
                      //  LOG(INFO)<<"SUMM of EXP"<<Sum_Exp_Time;
                        writeUsrkeywordsFrmvectDriver(finalfilename, FirstHDU_headertrack);
                       
                             
                       float crpix1=uvitObj.IMG_DIM_FIpc/2;
float crpix2=uvitObj.IMG_DIM_FIpc/2;
float crota1=0.0f,cdelt1=0.0f,cdelt2=0.0f;
float crota2=0.0f;
 int factor_delta=uvitObj.IMG_DIM_FIpc/600;  
                  cdelt1=cdelt2=(3.3307/3600)/factor_delta;      
                      cdelt1=-cdelt1;///cos(center_dec*M_PI/180);    
                        
     float rapnt,decpnt;                   
                        fits_open_file(&fout, finalfilename, READWRITE, &status);
                        printError(status, "Error in Creating the image for Signal Fie", finalfilename);
                         fits_update_key (fout , TFLOAT , "EXP_TIME" , &Sum_Exp_Time , NULL , &status) ;
                         printError(status, "Error in updating keyword value of EXP_TIME", finalfilename);
                       //  fits_read_key (fout , TFLOAT , "RA_PNT" , &rapnt , NULL , &status) ;
                        // fits_read_key (fout , TFLOAT , "DEC_PNT" , &decpnt , NULL , &status) ;
                        fits_update_key(fout, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(fout, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL1", &RA_pnt, "", &status);                                              printError(status,"");  
    fits_update_key(fout, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(fout, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL2", &DEC_pnt, "", &status);                                    printError(status,"");   
    fits_update_key(fout, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(fout, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
                        for (int i=0;i<IndividualOrbit_Dir.size();i++){
                             fits_write_history(fout, IndividualOrbit_Dir[i].c_str(), &status);
                             
                         }
			fits_write_history(fout,"Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV",&status);
                         fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", finalfilename);
                     
                         writeUsrkeywordsFrmvectDriver(finalExpfilename, FirstHDU_headertrack);
                        fits_open_file(&fout, finalExpfilename, READWRITE, &status);
                        printError(status, "Error in Creating the image for Signal Fie", finalExpfilename);
                         fits_update_key (fout , TFLOAT , "EXP_TIME" , &Sum_Exp_Time , NULL , &status) ;
                         printError(status, "Error in updating keyword value of EXP_TIME", finalExpfilename);
                 fits_update_key(fout, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(fout, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL1", &RA_pnt, "", &status);                                              printError(status,"");  
    fits_update_key(fout, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(fout, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL2", &DEC_pnt, "", &status);                                    printError(status,"");   
    fits_update_key(fout, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(fout, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");    
                         for (int i=0;i<IndividualOrbit_Dir.size();i++){
                             fits_write_history(fout, IndividualOrbit_Dir[i].c_str(), &status);
                             
                         }
	fits_write_history(fout,"Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV",&status);
                         fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", finalExpfilename);
                        
                        writeUsrkeywordsFrmvectDriver(finalNoiseMapfilename, FirstHDU_headertrack);
                        fits_open_file(&fout, finalNoiseMapfilename, READWRITE, &status);
                        printError(status, "Error in Creating the image for Signal Fie", finalNoiseMapfilename);
                        fits_update_key (fout , TFLOAT , "EXP_TIME" , &Sum_Exp_Time , NULL , &status) ;
                        printError(status, "Error in updating keyword value of EXP_TIME", finalNoiseMapfilename);
                        fits_update_key(fout, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(fout, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL1", &RA_pnt, "", &status);                                              printError(status,"");  
    fits_update_key(fout, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(fout, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL2", &DEC_pnt, "", &status);                                    printError(status,"");   
    fits_update_key(fout, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(fout, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
                        
                        for (int i=0;i<IndividualOrbit_Dir.size();i++){
                             fits_write_history(fout, IndividualOrbit_Dir[i].c_str(), &status);
                             
                         }
                        fits_write_history(fout,"Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV",&status);
                      
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", finalNoiseMapfilename);
                        
                        LOG(INFO) << "NOW FULL FRAMEASTROMETRY FOR  COMBINED ORBIT  IMAGES FOR NUV";


                        outputFullFrameAst = level2_Origname + "NUV_FullFrameAst_" + curr_filtercal_val+"_W" + (string) window_index_track[window_index];
                        LOG(INFO) << outputFullFrameAst;
                        if (DirExists((char*) outputFullFrameAst.c_str())) {
                            LOG(ERROR) << "Directory exists and clobber=yes";
                            cmd = (string) "rm -rf " + (string) outputFullFrameAst;
                            system(cmd.c_str());
                        }

                        cmd = "mkdir -p " + (string) outputFullFrameAst;
                        system(cmd.c_str());
 delete[] tempArrayXinvert;

                        uvtFullFrameAst ast_obj;
                        ast_obj.read((char*) finalImagePath.c_str(),finalfilename, (char*) uvitObj.caldbindir.c_str(), (char*) filenamesnr.c_str(),finalNoiseMapfilename,finalExpfilename, (char*) outputFullFrameAst.c_str(), (char*) ATTfile, (char*) uvitObj.att_timecolpc, (char*) uvitObj.att_qcolpc, (char*) uvitObj.caldbindir.c_str(),
                                sd_multfactor, uvitObj.minimum_No_of_Starspc, uvitObj.refine_Winsizepc / 8, uvitObj.centroid_Winsizepc / 8, uvitObj.databasenamepc, uvitObj.search_algo_ctlgpc,
                                uvitObj.len_apc, uvitObj.len_bpc, uvitObj.rad_searchpc, uvitObj.clobberpc, uvitObj.historypc, 1);
                        ast_obj.display();


                        status = ast_obj.uvtFullFrmAstProcess();
                        if (status) {
                            LOG(ERROR) << endl << "Error in full frame astrometry  module";
                           continue;

                        }
                        float diffAddra = ast_obj.getDiffRAval();
                        float diffAdddec = ast_obj.getDiffDECval();
                        LOG(INFO) << "Diff RA " << diffAddra << " " << diffAdddec;
                        delete[] X_loc_snr, Y_loc_snr, enp_snr, bad_flag_snr, mult_phtn_snr;
                    }
                }
            }//added

        }//if condition

        //FOR FUV
        LOG(INFO) << "FUV CHAIN START..";
        ExpPath.clear();
      
        RollAngles_track.clear();
        Full_framePath.clear();
        snrFilePath.clear();
        numStars.clear();
        uvitObj.header_info.clear();

        initArray(Frame_data_cum, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
        initArray(Frame_Expdata_cum, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
     //   initArray(FinalArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
        initArray(FinalArray, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc , 0.0f);
        initArray(Final_ExpArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
        initArray(counterArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
        initArray(counterArray_exp, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);

/*
  uvitObj.OutputDir_l2pcFUV.clear();
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_1");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_2");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_3");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_4");
uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_5");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_6");
uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_7");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_8");
uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_9");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_10");
uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_11");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_12");
uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_13");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_14");
uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_15");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_16");
uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_17");
                        uvitObj.OutputDir_l2pcFUV.push_back("/data2/swarna/uvit/udata1/m5962_v54/driver_out/_FUV_18");
*/



   //     uvitObj.OutputDir_l2pcFUV.clear();
       // uvitObj.OutputDir_l2pcFUV.push_back("output5962new_FUV_4");
      //  uvitObj.OutputDir_l2pcFUV.push_back("output5962new_FUV_5");
      //  uvitObj.OutputDir_l2pcFUV.push_back("output5962new_FUV_7");
      //  uvitObj.OutputDir_l2pcFUV.push_back("output5962new_FUV_8");
//uvitObj.OutputDir_l2pcFUV.push_back("output5962new_FUV_9");
//uvitObj.OutputDir_l2pcFUV.push_back("output5962new_FUV_10");



/*
        uvitObj.OutputDir_l2pcFUV.clear();
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_1");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_2");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_3");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_4");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_5");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_6");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_7");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_8");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_9");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_10");
       uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_11");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_12");
       uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_13");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_14");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_15");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_16");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_17");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_18");
       uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_19");
       uvitObj.OutputDir_l2pcFUV.push_back("output10945_FUV_20");
           
 uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_1");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_2");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_3");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_4");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_5");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_6");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_7");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_8");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_9");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_10");
       uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_11");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_12");
       uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_13");
        uvitObj.OutputDir_l2pcFUV.push_back("output10945_RAPC_FUV_14");
    */   
                       
                             float Exp_time_keyword;
          flagRefFrameNotFound=FALSE;
        if (uvitObj.OutputDir_l2pcFUV.size() > 0) {
            for (int window_index = 0; window_index < window_index_track.size(); window_index++) {
                for (int i = 1; i <= 7; i++) {
                    flagRefFrameNotFound=FALSE;
                    initArray(Frame_data_cum, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                    initArray(Frame_Expdata_cum, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
//                    initArray(FinalArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                  initArray(FinalArray, uvitObj.IMG_DIM_FIpc  * uvitObj.IMG_DIM_FIpc , 0.0f);
                    initArray(Final_ExpArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                    initArray(counterArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                    initArray(counterArray_exp, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                    IndividualOrbit_Dir.clear();
                    Full_framePath.clear();
                    eventList_RADEC.clear();
                    ExpPath.clear();
                    RollAngles_track.clear();
                    snrFilePath.clear();
                    numStars.clear();
                    cx_ref.clear();
                    cy_ref.clear();
                    ci_ref.clear();
                    sprintf(curr_filtercal_val, "%s%d", "F", (int) i);


                    int cnt_integratinTime = 0;
                    double integrationTime = 0.0f; //storing individual integration time for individual frame.
                    double Final_integrationTime = 0.0f;
                    for (int i = 0; i < uvitObj.OutputDir_l2pcFUV.size(); i++) {
                        files_OfDirNUV.clear();
                        getFiles(uvitObj.OutputDir_l2pcFUV[i], files_OfDirNUV);
                        if (files_OfDirNUV.size() == 0 && window_index == 0 && i == 1) {
                            cmd = "rm -r " + uvitObj.OutputDir_l2pcFUV[i];
                            LOG(INFO) << cmd;
                            system(cmd.c_str());

                        }


                        for (int j = 0; j < files_OfDirNUV.size(); j++) {


                            if (strstr(files_OfDirNUV[j].c_str(), "_sig_ra-dec.") != NULL) {
                                LOG(INFO) << "Reading " << files_OfDirNUV[j];

                                readImage((char*) files_OfDirNUV[j].c_str(), 1, Frame_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc);
                                readKeywords((char*) files_OfDirNUV[j].c_str(), 1, 5, TDOUBLE, "INT_TIME", &integrationTime, TDOUBLE, "ROLLAPPLIED", &Rollangle, TSTRING, "FILTER", &filter_val, TINT, "WIN_X_SZ", &Win_Size, TFLOAT, "EXP_TIME", &Exp_time_keyword);

                                if (strcmp(filter_val, curr_filtercal_val) != 0 || atoi(window_index_track[window_index].c_str()) != Win_Size) {
                                    LOG(INFO) << "NO filter match found ...Checking next file..";
                                    //                            Full_framePath.clear();
                                    //                            snrFilePath.clear();

                                    break;
                                }
                                //Read header data from image file.
                                Final_integrationTime = Final_integrationTime + integrationTime;
                                cnt_integratinTime++;

                                Full_framePath.push_back(files_OfDirNUV[j]);
                                RollAngles_track.push_back(Rollangle);

                                //                                NonZeroInt.clear();
                                //                                for (int k = 0; k < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; k++) {
                                //                                    if (Frame_data[i] != 0) {
                                //                                        NonZeroInt.push_back(Frame_data[k]);
                                //                                    }
                                //                                }
                                //
                                //                                sd_temp = getSD(NonZeroInt.data(), NonZeroInt.size());
                                //                                mean_Ofimage = getmean(NonZeroInt.data(), NonZeroInt.size());
                                //                                //    thr =  sd_temp* sd_mul_factor ;
                                //                                thr_peaks = mean_Ofimage + sd_temp * 3.5;
                                //                                cnt_stars = 0;
                                //                                for (int index = 0; index < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; index++) {
                                //                                    if (Frame_data[index] > thr_peaks && Frame_data[index] > 0) {
                                //                                        cnt_stars++; //stars calculation
                                //
                                //                                    }
                                //                                }
                                numStars.push_back(Exp_time_keyword);

                               // LOG(INFO)<<"Exposure Time->"<<Exp_time_keyword;
                            }
                            if (strstr(files_OfDirNUV[j].c_str(), "_sig_ra-decexp.") != NULL) {
                                LOG(INFO) << "Reading " << files_OfDirNUV[j];
                                readKeywords((char*) files_OfDirNUV[j].c_str(), 1, 2, TSTRING, "FILTER", &filter_val, TINT, "WIN_X_SZ", &Win_Size);
                                if (strcmp(filter_val, curr_filtercal_val) != 0 || atoi(window_index_track[window_index].c_str()) != Win_Size) {
                                    LOG(INFO) << "NO filter match found ...Checking next file..";
                                    //                            Full_framePath.clear();
                                    //                            snrFilePath.clear();

                                    break;
                                }
                                ExpPath.push_back(files_OfDirNUV[j]);
                            }
                               if (strstr(files_OfDirNUV[j].c_str(), "l2_radec.") != NULL) {                              
//if (strstr(files_OfDirNUV[j].c_str(), "ra-decexp.fits") != NULL) {    
LOG(INFO) << "Reading " << files_OfDirNUV[j];
                                readKeywords((char*) files_OfDirNUV[j].c_str(), 2, 2, TSTRING, "FILTER", &filter_val, TINT, "WIN_X_SZ", &Win_Size);
                                if (strcmp(filter_val, curr_filtercal_val) != 0 || atoi(window_index_track[window_index].c_str()) != Win_Size) {
                                    LOG(INFO) << "NO filter match found ...Checking next file..";
                                    //                            Full_framePath.clear();
                                    //                            snrFilePath.clear();

                                    break;
                                }
                                eventList_RADEC.push_back(files_OfDirNUV[j]);
                            }

                            if (strstr(files_OfDirNUV[j].c_str(), "snr") != NULL) {
                                snrFilePath.push_back(files_OfDirNUV[j]);
                            }

                        }


                    }
                    if (cnt_integratinTime != 0)
                        Final_integrationTime = Final_integrationTime / cnt_integratinTime;
//
//                    for (int i = 0; i < numStars.size(); i++) {
//                            LOG(INFO) <<"111"<< Full_framePath[i]<<" "<<numStars[i];
//                         
//                        }
                    
                    if (Full_framePath.size() > 0) {
                       // LOG(INFO) << numStars.size();
                        string temp_str, tempsnrstr, tempexpstr;
                        float temp_numStars;
                        float temp_rollangle;
                         for (int i = 0; i < numStars.size(); i++) {
                            for (int j = numStars.size() - 1; j > i; j--) {
                                if (numStars[j - 1] < numStars[j]) {
                                    swap1(numStars[j], numStars[j - 1]);
                                    swap1(Full_framePath[j], Full_framePath[j - 1]);
                                    swap1(RollAngles_track[j], RollAngles_track[j - 1]);
                                    swap1(snrFilePath[j], snrFilePath[j - 1]);
                                    swap1(ExpPath[j], ExpPath[j - 1]);
                                    swap1(eventList_RADEC[j], eventList_RADEC[j - 1]);
                                }
                            }
                        }
//                        for (int i = 1; i < numStars.size(); i++) {
//                            if (numStars[i] > numStars[i - 1]) {
//                                temp_numStars=numStars[i - 1];
//                                numStars[i - 1]=numStars[i];
//                                numStars[i]=temp_numStars;
//                                temp_str = Full_framePath[i - 1];
//                                Full_framePath[i - 1] = Full_framePath[i];
//                                Full_framePath[i] = temp_str;
//                                temp_rollangle = RollAngles_track[i - 1];
//                                RollAngles_track[i - 1] = RollAngles_track[i];
//                                RollAngles_track[i] = temp_rollangle;
//                                tempsnrstr = snrFilePath[i - 1];
//                                snrFilePath[i - 1] = snrFilePath[i];
//                                snrFilePath[i] = tempsnrstr;
//                                tempexpstr = ExpPath[i - 1];
//                                ExpPath[i - 1] = ExpPath[i];
//                                ExpPath[i] = tempexpstr;
//                            }
//                        }
                        //Now we identified which frame to be used as a reference frame;

//                        for (int i = 0; i < numStars.size(); i++) {
//                            LOG(INFO) << Full_framePath[i];
////                            LOG(INFO) << snrFilePath[i];
////                            LOG(INFO) << ExpPath[i];
//                        }
                        // exit(1);
                        //opening the frmaes.
                        //   vector<float> cx_ref,cy_ref,ci_ref;//stores the centroid values for reference frame.
                        // int  min_stars_match,cnt,matching_points;       
                        //  vector<float> x_ref_arr , y_ref_arr , x_arr , y_arr , temp_x_arr , temp_y_arr ;
                        //  double x_dx = 0.0 , y_dy = 0.0 , theta_dt = 0.0 ;
                        // vector<float> New_X_ref,New_Y_ref,New_x_arr,New_y_arr,New_Xdiff,New_Ydiff;
                        //   vector<float> int_new_one,int_new_two;
                        //  double ctheta,stheta;
                        //  int x_index,y_index;
                        //  double new_index_x,new_index_y;

                        //sd_multfactor=uvitObj.sd_multi_factor_defaultpc;
                        //refineWinSizedefault=uvitObj.refine_Winsizepc;
                        //centWinSizedefault=uvitObj.centroid_Winsizepc;
                        //                        uvitObj.sd_multi_factor_defaultpc = sd_multfactor;
                        //                        uvitObj.refine_Winsizepc = refineWinSizedefault;
                        //                        uvitObj.centroid_Winsizepc = centWinSizedefault;
                        //  char nameprfx[NAMESIZE];
                        //int xsizeimage,ysizeimage;
                        double RA_pnt,DEC_pnt;
                        // char date_obs[NAMESIZE],roll_angle[NAMESIZE];
                        //readKeywords((char*) Full_framePath[0].c_str(),1,7,TSTRING,"NAMEPRFX",nameprfx,TINT,"XSIZE",&xsizeimage,TINT,"YSIZE",&ysizeimage,TFLOAT,"RA_PNT",&RA_pnt
                        //     ,TFLOAT,"DEC_PNT",&DEC_pnt,TSTRING,"DATE",date_obs,TFLOAT,"ROLLAPPLIED",&roll_angle);

                        uvitObj.copyAllheaderKeys((char*) Full_framePath[0].c_str());
                        vector<string> FirstHDU_headertrack;
                        FirstHDU_headertrack = uvitObj.header_info;
                        track_headerInfoL1.clear();
                        bool flag_shiftFound = TRUE;
                        int cnt_fuv = 0;
                        fitsfile *fin;
                        vector<string> L1headerInfo;
                        status = 0;
                        string FinalStr;
                        int pos,pos1,pos2;
                        string ExtName,Orbit_Dir,temp_strForSubstr,tempsec_str;
                       
                        Sum_Exp_Time=0.0f;
vector<float> cent_X,cent_Y,cent_I;
                       float Max_value_Exp;
int cnter_elements=0;
double peakofexp=0.0f;
int indexx,indexy,cntpixels;
float Sum_ofpixels,valuetoCmpr;
   RA_pnt=0.0f;DEC_pnt=0.0f;
                        for (int i = 0; i < Full_framePath.size(); i++) {
                         
				Sum_ofpixels=0.0f;
                            flag_shiftFound = TRUE;
                            readImage((char*) Full_framePath[i].c_str(), 1, Frame_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc); //reading the image.
                            readImage((char*) ExpPath[i].c_str(), 1, Frame_Exposure_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc); //reading the  exposure image
  readKeywords((char*)  Full_framePath[i].c_str(), 1, 2, TFLOAT,"EXP_TIME",&Exp_Time_Keyword,TDOUBLE,"PEAK_OF_EXP",&peakofexp);

   Max_value_Exp=peakofexp;
 if (i==0) {
                                readKeywords((char*) Full_framePath[i].c_str(), 1, 2, TDOUBLE, "CRVAL1", &RA_pnt, TDOUBLE, "CRVAL2", &DEC_pnt);

                            }
//LOG(INFO)<<Max_value_Exp;exit(1);
//for (int i =0;i<uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc;i++)
//{
//if(Frame_Exposure_data[i]<(20*Max_value_Exp/100)){
//Frame_Exposure_data[i]=0.0f;
//}
//}                              


   




                            pos = Full_framePath[i].find("/uvtF.");
                            pos1=Full_framePath[i].find("_FUV_");
                            tempsec_str=Full_framePath[i].substr(0,pos1);
                            ExtName = Full_framePath[i].substr(pos + 1, 7);
                             temp_strForSubstr=Full_framePath[i].substr(pos1, Full_framePath[i].length()-1);
                             pos2=temp_strForSubstr.find("/");
                              Orbit_Dir=tempsec_str+Full_framePath[i].substr(pos1, pos2);
                            fits_open_file(&fin, Full_framePath[i].c_str(), READONLY, &status);
                            printError(status, "Error in Opening the file", (char*) Full_framePath[i].c_str());
                            fits_movabs_hdu(fin, 2, NULL, &status);
                            printError(status, "Error in moving to 2nd HDU", (char*) Full_framePath[i].c_str());
                            copyUsrkeywrdsTovect(fin, L1headerInfo);
                            FinalStr = "EXTNAME = " + ExtName;
                             if(i==0){
                                current_ExtName="Reference_Orbit= "+Orbit_Dir+", "+ExtName;
                            }else{
                                current_ExtName="Matched_Orbit= "+Orbit_Dir+", "+ExtName;
                            }
                            
                            L1headerInfo.push_back(FinalStr);
                            uvitObj.copyAllheaderKeys((char*) Full_framePath[i].c_str());
                            track_headerInfoL1.push_back(L1headerInfo);
                            uvitObj.sd_multi_factor_defaultpc = sd_multfactor;
                            uvitObj.refine_Winsizepc = refineWinSizedefault * 8;
                            uvitObj.centroid_Winsizepc = centWinSizedefault * 8;
                            fits_close_file(fin, &status);
                            //copy header info from the image.
                            uvitObj.copyAllheaderKeys((char*) Full_framePath[i].c_str());
cntpixels=0;   
                            if (i == 0) {//This is for reference frame.
                                status = uvitObj.findStar_algo1(Frame_data, PC); // for first frame taken as  reference
                                if (status) {
                                    // LOG (ERROR) << endl << "***Error in finding star algorithm 1 for frame  " << infile << "  ***" << endl ;
				flagRefFrameNotFound=TRUE;
				//return (EXIT_FAILURE);
				break;
                                  
                                }
                              //  LOG(INFO)<<uvitObj.Cx.size ()<<endl;

for (int i = 0; i < uvitObj.Cx.size(); i++) {
cntpixels=0;Sum_ofpixels=0.0f;
indexx=round(uvitObj.Cx[i]);
indexy=round(uvitObj.Cy[i]);
	for (int j=indexx-12;j<=indexx+12;j++){
for (int k =indexy-12;k<=indexy+12;k++){
if((int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j)) >=0 && (int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))<uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc && Frame_Exposure_data[(int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))]!= INVALID_PIX_VALUE ){
Sum_ofpixels=Sum_ofpixels+Frame_Exposure_data[(int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))];
cntpixels++;
}
}
}
valuetoCmpr=Sum_ofpixels/cntpixels;
if(valuetoCmpr>(20*Max_value_Exp/100)){
cx_ref.push_back(uvitObj.Cx[i]);
                                    cy_ref.push_back(uvitObj.Cy[i]);
                                    ci_ref.push_back(uvitObj.Ci[i]);
}

}
if(cx_ref.size()==0){
LOG(INFO)<<"***Reference frame star not found***";
flagRefFrameNotFound=TRUE;
break;
}

                        //        for (int i = 0; i < uvitObj.Cx.size(); i++) {


//if(Frame_Exposure_data[(int)(round(uvitObj.Cy[i])*uvitObj.IMG_DIM_FIpc+round(uvitObj.Cx[i]))]!=0.0f){
  //                                  cx_ref.push_back(uvitObj.Cx[i]);
    //                                cy_ref.push_back(uvitObj.Cy[i]);
      //                              ci_ref.push_back(uvitObj.Ci[i]);

//}

//                                }
//LOG(INFO)<<uvitObj.Cx.size()<<" "<<cx_ref.size();
                                uvitObj.Cx.clear();
                                uvitObj.Cy.clear();
                                uvitObj.Ci.clear();
                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                                    if (Frame_data[i] != INVALID_PIX_VALUE && Frame_Exposure_data[i] != INVALID_PIX_VALUE) {
                                        Frame_data_cum[i] = Frame_data_cum[i] + Frame_data[i] * Frame_Exposure_data[i];
                                        Frame_Expdata_cum[i] = Frame_Expdata_cum[i] + Frame_Exposure_data[i];
                                    }
                                }

                            } else {
				
				if(abs(RollAngles_track[i]-RollAngles_track[0])>2){
					LOG(ERROR)<<"***Roll Angle ISSUE***<<Ref Frame Roll->"<<RollAngles_track[0]<<", Current frame Roll->"<<RollAngles_track[i];
					continue;
				}


                                status = uvitObj.findStar_algo1(Frame_data, PC); // for first frame taken as  reference
                                if (status) {
                                     LOG (ERROR) << endl << "***Error in finding star  ***" << endl ;
                                   // return (EXIT_FAILURE);
					continue;
                                }
                                uvitObj.Rx.clear();
                                uvitObj.Ry.clear();
                                uvitObj.Rval.clear();
                                uvitObj.Fx.clear();
                                uvitObj.Fy.clear();
                                uvitObj.Fval.clear();
cent_X.clear();
cent_Y.clear();
cent_I.clear();
			for (int i = 0; i < uvitObj.Cx.size(); i++) {
cntpixels=0;Sum_ofpixels=0.0f;
indexx=round(uvitObj.Cx[i]);
indexy=round(uvitObj.Cy[i]);
	for (int j=indexx-12;j<=indexx+12;j++){
for (int k =indexy-12;k<=indexy+12;k++){
if((int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j)) >=0 && (int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))<uvitObj.IMG_DIM_FIpc*uvitObj.IMG_DIM_FIpc && Frame_Exposure_data[(int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))]!= INVALID_PIX_VALUE ){
Sum_ofpixels=Sum_ofpixels+Frame_Exposure_data[(int)(round(k)*uvitObj.IMG_DIM_FIpc+round(j))];
cntpixels++;
}
}
}
valuetoCmpr=Sum_ofpixels/cntpixels;
if(valuetoCmpr>(20*Max_value_Exp/100)){
cent_X.push_back(uvitObj.Cx[i]);
cent_Y.push_back(uvitObj.Cy[i]);
cent_I.push_back(uvitObj.Ci[i]);
}

}


uvitObj.Cx=cent_X;
uvitObj.Cy=cent_Y;
uvitObj.Ci=cent_I;
if(uvitObj.Cx.size()==0){
continue;
}

                                //                    for (int i = 0; i < uvitObj.Cx.size(); i++) {
                                //                        LOG(INFO) << uvitObj.Cx[i] << " " << uvitObj.Cy[i];
                                //                    }
                                //                    LOG(INFO) << "=================";
                                //                    for (int i = 0; i < cx_ref.size(); i++) {
                                //                        LOG(INFO) << cx_ref[i] << " " << cy_ref[i];
                                //                    }
                                if (uvitObj.Cx.size() >= cx_ref.size()) {
                                    matching_points = cx_ref.size();
                                } else {
                                    matching_points = uvitObj.Cx.size();
                                }

                                x_ref_arr.clear();
                                y_ref_arr.clear();
                                x_arr.clear();
                                y_arr.clear();
                                temp_x_arr.clear();
                                temp_y_arr.clear();
                                //                    min_stars_match = (int) ((100 -PERCENTGE_ERR_ALLOWED) * matching_points / 100);
                                min_stars_match = 1;
                                cnt = 0;
                               
                               
////                                 flag_TestOnecheck=FALSE;
                                 flag_TestTwocheck=FALSE;
                               // while (cnt < min_stars_match) {
//                                    if (uvitObj.diff_Distpc > 256) {
//                                        LOG(ERROR) << "No matches found...";
//                                        flag_shiftFound = FALSE;
//                                    }
//                                    x_ref_arr.clear();
//                                    y_ref_arr.clear();
//                                    x_arr.clear();
//                                    y_arr.clear();
//                                    temp_x_arr.clear();
//                                    temp_y_arr.clear();

                                     uvitObj.matchStars(cx_ref.size(), uvitObj.Cx.size(), 8, cx_ref.data(), cy_ref.data(), uvitObj.Cx.data(), uvitObj.Cy.data(), x_dx,y_dy,flag_TestOnecheck,flag_TestTwocheck);

                                   // uvitObj.diff_Distpc = uvitObj.diff_Distpc * 2;

                               // }
                                     if(flag_TestOnecheck==FALSE && flag_TestTwocheck==FALSE){
                                         LOG(INFO)<<"NOT FOUND....";
                                         continue;
                                     }
//                                LOG(INFO) << "TOTAL STAR FOUND-> " << cnt << " " << uvitObj.diff_Distpc << " " << uvitObj.sd_multi_factor_defaultpcfuv;
//                                if (flag_shiftFound == FALSE) {
//                                    continue;
//                                }
//
//                                uvitObj.diff_Distpc = original_diff_Distpc;
//
//                                New_X_ref.clear(), New_Y_ref.clear(), New_x_arr.clear(), New_y_arr.clear(), New_Xdiff.clear(), New_Ydiff.clear();
//                                int_new_one.clear(), int_new_two.clear();
//
//                                if (cnt > 2) {
//
//                                    status = uvitObj.removeRecords(x_ref_arr, y_ref_arr, x_arr, y_arr, temp_x_arr, temp_y_arr, ci_ref.data(), uvitObj.Ci.data(), New_X_ref, New_Y_ref, New_x_arr, New_y_arr, New_Xdiff, New_Ydiff, int_new_one, int_new_two);
//                                    if (status) {
//                                        LOG(ERROR) << "Error in removing the records above mean+sigma";
//                                        return (EXIT_FAILURE);
//                                    }
//                                } else {
//                                    
//                                    New_X_ref = x_ref_arr;
//                                    New_Y_ref = y_ref_arr;
//                                    New_x_arr = x_arr;
//                                    New_y_arr = y_arr;
//                                    New_Xdiff = temp_x_arr;
//                                    New_Ydiff = temp_y_arr;
//                                    for(int i=0;i<New_X_ref.size();i++){
//                                        LOG(INFO)<<i+1<<"->"<<New_X_ref[i]<<" "<<New_Y_ref[i]<<" "<<New_x_arr[i]<<" "<<New_y_arr[i];
//                                    }
//                                    
//                                }
//
//                                //                    for (int i = 0; i < New_X_ref.size(); i++) {
//                                //
//                                //                        LOG(INFO) << New_X_ref[i] << " " << New_Y_ref[i] << " " << New_x_arr[i] << " " << New_y_arr[i] << " " << New_Xdiff[i] << " " << New_Ydiff[i] << endl;
//                                //
//                                //                    }
//
//
//                                //                 status = uvitObj.findShiftsNtheta (x_ref_arr.size () , x_ref_arr , y_ref_arr  ,x_arr , y_arr , New_Xdiff , New_Ydiff , 1,x_dx , y_dy , theta_dt) ;
//                                //        if (status)
//                                //        {
//                                //            LOG (INFO) << "Error in finding shifts n theta " << endl ;
//                                //            return (EXIT_FAILURE) ;
//                                //        } 
//                                //
//                                status = uvitObj.findShiftsNtheta(New_X_ref.size(), New_X_ref, New_Y_ref, New_x_arr, New_y_arr, New_Xdiff, New_Ydiff, 1, x_dx, y_dy, theta_dt);
//                                if (status) {
//                                    LOG(INFO) << "Error in finding shifts n theta " << endl;
//                                    return (EXIT_FAILURE);
//                                }
//                                  LOG(INFO) <<"The DX DY and DTHETA->"<<x_dx << " " << y_dy << " " << theta_dt;

                                
                                 fits_open_file (&fptr_snr , eventList_RADEC[i].c_str() , READWRITE , &status) ;
        printError (status , " Error in opening file " , (char*)eventList_RADEC[i].c_str()) ;
        fits_movabs_hdu (fptr_snr , 2 , NULL , &status) ;
        printError (status , "Error in moving to HDU " , (char*)eventList_RADEC[i].c_str()) ;
        fits_get_num_rows (fptr_snr , &nrows_snr , &status) ;
        printError (status , "Error in reading the number of rows" , (char*)eventList_RADEC[i].c_str()) ;
        X_loc_snr= new  float[nrows_snr];
        Y_loc_snr= new float[nrows_snr];
        bad_flag_snr = new float[nrows_snr];
        mult_phtn_snr= new float[nrows_snr];
        enp_snr= new float[nrows_snr];
        fits_read_col (fptr_snr , TFLOAT , 13 , 1 , 1 , nrows_snr  , NULL , Y_loc_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , (char*)eventList_RADEC[i].c_str()) ;
        fits_read_col (fptr_snr , TFLOAT , 12 , 1 , 1 , nrows_snr  , NULL , X_loc_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , (char*)eventList_RADEC[i].c_str()) ;
        fits_read_col (fptr_snr , TFLOAT , 9 , 1 , 1 , nrows_snr  , NULL , bad_flag_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , (char*)eventList_RADEC[i].c_str()) ;
        fits_read_col (fptr_snr , TFLOAT , 10 , 1 , 1 , nrows_snr  , NULL , mult_phtn_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , (char*)eventList_RADEC[i].c_str()) ;
         fits_read_col (fptr_snr , TFLOAT , 11 , 1 , 1 , nrows_snr  , NULL , enp_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , (char*)eventList_RADEC[i].c_str()) ;
        fits_close_file (fptr_snr , &status) ;
        printError (status , "Error closing the file " , (char*)eventList_RADEC[i].c_str()) ;         
                                     
                                
                                ctheta = 0.0f, stheta = 0.0f;
                                //                                ctheta = cos(-1.0 * theta_dt);
                                //                                stheta = sin(-1.0 * theta_dt);
                                ctheta = cos(-1.0 * theta_dt);
                                stheta = sin(-1.0 * theta_dt);
                                LOG(INFO) << "Loop Started for assigning the correction to the frames..";

                             //   initArray(FinalArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                initArray(FinalArray, uvitObj.IMG_DIM_FIpc  * uvitObj.IMG_DIM_FIpc , 0.0f);
                                initArray(Final_ExpArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);

                                initArray(SubdividedImage, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                initArray(Subdivided_expImage, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                //scale image to 9600.
                                initArray(counterArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                                initArray(counterArray_exp, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);

                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
//                                    if (Frame_data[i] != INVALID_PIX_VALUE)
//                                        Frame_data[i] = Frame_data[i] / 4;
                                    if (Frame_Exposure_data[i] != INVALID_PIX_VALUE)
                                        Frame_Exposure_data[i] = Frame_Exposure_data[i] / 4;
                                }
//                                performSubDivisionIM(Frame_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, SubdividedImage, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2);
                                performSubDivisionIM(Frame_Exposure_data, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, Subdivided_expImage, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2);

//                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * 2; i++) {
//                                    mid_X_new = i - uvitObj.IMG_DIM_FIpc;
//                                    for (int j = 0; j < uvitObj.IMG_DIM_FIpc * 2; j++) {
//                                        mid_Y_new = j - uvitObj.IMG_DIM_FIpc;
//
//                                        x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc - (x_dx * 8 * 2);
//                                        x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc - (y_dy * 8 * 2);
//                                        //              if(x1>0 && x1<FINALFRAMESIZE_REGAVG && x2>0 && x2<FINALFRAMESIZE_REGAVG){
//                                        //                  imageLocation_Array_X[(int)(round(x2)*FINALFRAMESIZE_REGAVG+round(x1))]=Regavg_subSampled_Array_Sig[j*FINALFRAMESIZE_REGAVG+i];
//                                        //              }
//                                        //if((int)(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1))>0 && (int)(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1))<uvitObj.IMG_DIM_FIpc * 2*uvitObj.IMG_DIM_FIpc * 2)
//                                        if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc * 2 && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc * 2) {
//                                            if (SubdividedImage[j * uvitObj.IMG_DIM_FIpc * 2 + i] != INVALID_PIX_VALUE && FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] != INVALID_PIX_VALUE) {
//                                                FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + SubdividedImage[j * uvitObj.IMG_DIM_FIpc * 2 + i];
//                                                counterArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = counterArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + 1;
//                                            } else {
//                                                FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = INVALID_PIX_VALUE;
//                                            }
//                                        }
//
//                                    }
//                                }
 for (int i =0;i <nrows_snr;i++){
                                    mid_X_new=X_loc_snr[i]-uvitObj.IMG_DIM_FIpc/2;
                                    mid_Y_new=Y_loc_snr[i]-uvitObj.IMG_DIM_FIpc/2;
                                     x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc/2 - (x_dx * 8);
                                     x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc/2 - (y_dy * 8 );
                                 
                                     if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc  && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc )
                                     FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc  + round(x1))]=FinalArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc  + round(x1))]+bad_flag_snr[i]*mult_phtn_snr[i]*enp_snr[i];
                                     
                                     
                                }

                                //now for exposure Array
                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * 2; i++) {
                                    mid_X_new = i - uvitObj.IMG_DIM_FIpc;
                                    for (int j = 0; j < uvitObj.IMG_DIM_FIpc * 2; j++) {
                                        mid_Y_new = j - uvitObj.IMG_DIM_FIpc;

                                        x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc - (x_dx * 8 * 2);
                                        x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc - (y_dy * 8 * 2);
                                        //              if(x1>0 && x1<FINALFRAMESIZE_REGAVG && x2>0 && x2<FINALFRAMESIZE_REGAVG){
                                        //                  imageLocation_Array_X[(int)(round(x2)*FINALFRAMESIZE_REGAVG+round(x1))]=Regavg_subSampled_Array_Sig[j*FINALFRAMESIZE_REGAVG+i];
                                        //              }
                                        //if((int)(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1))>0 && (int)(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1))<uvitObj.IMG_DIM_FIpc * 2*uvitObj.IMG_DIM_FIpc * 2)
                                        if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc * 2 && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc * 2) {
                                            if (Subdivided_expImage[j * uvitObj.IMG_DIM_FIpc * 2 + i] != INVALID_PIX_VALUE && Final_ExpArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] != INVALID_PIX_VALUE) {
                                                Final_ExpArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = Final_ExpArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + Subdivided_expImage[j * uvitObj.IMG_DIM_FIpc * 2 + i];
                                                counterArray_exp[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = counterArray_exp[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + 1;
                                            } else {
                                                Final_ExpArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = INVALID_PIX_VALUE;
                                            }
                                        }

                                    }
                                }


                                initArray(Frame_FinalArray, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                          //      ApplyBinning(FinalArray, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2, Frame_FinalArray, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, counterArray);

                                initArray(Frame_FinalExpArray, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                                ApplyBinning(Final_ExpArray, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2, Frame_FinalExpArray, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, counterArray_exp);



                                //                    for (int index = 0; index < uvitObj.IMG_DIM_FIpc; index++) {
                                //                        x_index = index - uvitObj.IMG_DIM_FIpc / 2;
                                //                        for (int jindex = 0; jindex < uvitObj.IMG_DIM_FIpc; jindex++) {
                                //                            if (Frame_data[jindex * uvitObj.IMG_DIM_FIpc + index] != INVALID_PIX_VALUE) {
                                //
                                //                                y_index = jindex - uvitObj.IMG_DIM_FIpc / 2;
                                //                                //ROUNDING:Correction needed.previously rounding was only on the subset of full equation.
                                //                                new_index_x = round((x_index) * ctheta - (y_index) * stheta + uvitObj.IMG_DIM_FIpc / 2 - (x_dx * 8)); //new index x
                                //                                new_index_y = round((x_index) * stheta + (y_index) * ctheta + uvitObj.IMG_DIM_FIpc / 2 - (y_dy * 8)); //new index y
                                //
                                //                                if (round(new_index_x) < uvitObj.IMG_DIM_FIpc && round(new_index_x) > 0 && round(new_index_y) > 0 && round(new_index_y) < uvitObj.IMG_DIM_FIpc) {
                                //                                    // cnt_loop ++ ;
                                //                                    //Rounding:No correction needed.As rounding is applied on each direction i.e X and Y.
                                //                                    FinalArray[(int) (round(new_index_y) * uvitObj.IMG_DIM_FIpc + round(new_index_x))] = Frame_data[jindex * uvitObj.IMG_DIM_FIpc + index];
                                //
                                //                                }
                                //                            }
                                //
                                //                        }
                                //                    }

                                
//                                for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
//                                    if (Frame_FinalArray[i] != INVALID_PIX_VALUE && Frame_FinalExpArray[i] != INVALID_PIX_VALUE && Frame_data_cum[i] != INVALID_PIX_VALUE && Frame_Expdata_cum[i] != INVALID_PIX_VALUE) {
//                                        Frame_data_cum[i] = Frame_data_cum[i] + Frame_FinalArray[i] * Frame_FinalExpArray[i];
//                                        Frame_Expdata_cum[i] = Frame_Expdata_cum[i] + Frame_FinalExpArray[i];
//                                    }
//                                }
                                 for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                                    if (FinalArray[i] != INVALID_PIX_VALUE && Frame_FinalExpArray[i] != INVALID_PIX_VALUE && Frame_data_cum[i] != INVALID_PIX_VALUE && Frame_Expdata_cum[i] != INVALID_PIX_VALUE) {
                                        Frame_data_cum[i] = Frame_data_cum[i] + FinalArray[i];// * Frame_FinalExpArray[i];
                                        Frame_Expdata_cum[i] = Frame_Expdata_cum[i] + Frame_FinalExpArray[i];
                                    }
                                }
                                
                                
//                                LOG(INFO)<<"*****INSIDE ELSE"<<" "<<Exp_Time_Keyword;
                                uvitObj.Cx.clear();
                                uvitObj.Cy.clear();
                                uvitObj.Ci.clear();
                            }

                            cnt_fuv++;
                              Sum_Exp_Time=Sum_Exp_Time+Exp_Time_Keyword;
                              IndividualOrbit_Dir.push_back(current_ExtName);
                        }
                        //noise map calculation
if(flagRefFrameNotFound==TRUE){
LOG(INFO)<<"***NOt going Further ,reference frame Star not found***";
continue;
}



                        float * tempArrayXinvert = new float [uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc; i++) {
                            for (int j = 0; j < uvitObj.IMG_DIM_FIpc; j++) {
				tempArrayXinvert[i * uvitObj.IMG_DIM_FIpc + j]=0.0f;
                                if (Frame_Expdata_cum[i * uvitObj.IMG_DIM_FIpc + (j)] != 0)
                                    tempArrayXinvert[i * uvitObj.IMG_DIM_FIpc + j] = Frame_data_cum[i * uvitObj.IMG_DIM_FIpc + (j)] / (Frame_Expdata_cum[i * uvitObj.IMG_DIM_FIpc + (j)]);

                            }
                        }





                        RADECimage = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
                        initArray(RADECimage, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                        initArray(RADECExpimage, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);

                        initArray(counterArray, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        initArray(counterArray_exp, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        // ctheta = cos(RollAngles_track[0] * M_PI / 180);
                        //stheta = sin(RollAngles_track[0] * M_PI / 180);

                        //bilinear interpolation
                        //            fits_open_file(&fptr_snr, snrFilePath[0].c_str(), READWRITE, &status);
                        //            printError(status, " Error in opening file ", (char*) snrFilePath[0].c_str());
                        //            fits_movabs_hdu(fptr_snr, 2, NULL, &status);
                        //            printError(status, "Error in moving to HDU ", (char*) snrFilePath[0].c_str());
                        //            fits_get_num_rows(fptr_snr, &nrows_snr, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            X_loc_snr = new float[nrows_snr];
                        //            Y_loc_snr = new float[nrows_snr];
                        //            bad_flag_snr = new float[nrows_snr];
                        //            mult_phtn_snr = new float[nrows_snr];
                        //            enp_snr = new float[nrows_snr];
                        //            fits_read_col(fptr_snr, TFLOAT, 5, 1, 1, nrows_snr, NULL, Y_loc_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_read_col(fptr_snr, TFLOAT, 4, 1, 1, nrows_snr, NULL, X_loc_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_read_col(fptr_snr, TFLOAT, 8, 1, 1, nrows_snr, NULL, bad_flag_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_read_col(fptr_snr, TFLOAT, 9, 1, 1, nrows_snr, NULL, mult_phtn_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_read_col(fptr_snr, TFLOAT, 10, 1, 1, nrows_snr, NULL, enp_snr, NULL, &status);
                        //            printError(status, "Error in reading the number of rows", (char*) snrFilePath[0].c_str());
                        //            fits_close_file(fptr_snr, &status);
                        //            printError(status, "Error closing the file ", (char*) snrFilePath[0].c_str());



                        //            ctheta = cos((360 - RollAngles_track[0]) * M_PI / 180);
                        //            stheta = sin((360 - RollAngles_track[0]) * M_PI / 180);
                        //                        ctheta = cos(-(RollAngles_track[0]) * M_PI / 180);
                        //                        stheta = sin(-(RollAngles_track[0]) * M_PI / 180);
                        
                        ctheta = cos(-(RollAngles_track[0]) * M_PI / 180);
                        stheta = sin(-(RollAngles_track[0]) * M_PI / 180);
                        //            for (int i = 0; i < nrows_snr; i++) {
                        //                x1 = (((X_loc_snr[i] - uvitObj.IMG_DIM_FIpc / 2) * (ctheta)) - ((Y_loc_snr[i] - uvitObj.IMG_DIM_FIpc / 2) * (stheta))) + uvitObj.IMG_DIM_FIpc / 2;
                        //                x2 = (((X_loc_snr[i] - uvitObj.IMG_DIM_FIpc / 2) * (stheta)) +((Y_loc_snr[i] - uvitObj.IMG_DIM_FIpc / 2)* (ctheta))) + uvitObj.IMG_DIM_FIpc / 2;
                        //
                        //                if ((int) ((int) (round(x2)) * uvitObj.IMG_DIM_FIpc + (int) (round(x1))) < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc) {
                        //                    //Rotated_image[(int)(FINALFRAMESIZE_REGAVG-round(x2))*FINALFRAMESIZE_REGAVG+(int)(round(x1))]=Rotated_image[(int)(FINALFRAMESIZE_REGAVG-round(x2))*FINALFRAMESIZE_REGAVG+(int)(round(x1))]+effective_NumPhotons[i]*mult_temp[i]*badFlag_temp[i];
                        //                    //Rotated_image[(int)(round(x1))*FINALFRAMESIZE_REGAVG+(int)(FINALFRAMESIZE_REGAVG-round(x2))]=Rotated_image[(int)(round(x1))*FINALFRAMESIZE_REGAVG+(int)(FINALFRAMESIZE_REGAVG-round(x2))]+effective_NumPhotons[i]*mult_temp[i]*badFlag_temp[i];
                        //                    RADECimage[(int) (round(x2)) * uvitObj.IMG_DIM_FIpc + (int) (round(x1))] = RADECimage[(int) (round(x2)) * uvitObj.IMG_DIM_FIpc + (int) (round(x1))] + enp_snr[i] * mult_phtn_snr[i] * bad_flag_snr[i];
                        //                }
                        //            }
                        float *tempArrSubDivided = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
                        float *tempExpArrSubDivided = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
                        float *RADECimageSubdivided = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
                        float *RADECExpimageSubdivided = new float[uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2];
                        initArray(tempArrSubDivided, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        initArray(tempExpArrSubDivided, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        initArray(RADECimageSubdivided, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);
                        initArray(RADECExpimageSubdivided, uvitObj.IMG_DIM_FIpc * 2 * uvitObj.IMG_DIM_FIpc * 2, 0.0f);

                    //    for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                        //    if (tempArrayXinvert[i] != INVALID_PIX_VALUE)
                          //      tempArrayXinvert[i] = tempArrayXinvert[i] / 4;
                     //       if (Frame_Expdata_cum[i] != INVALID_PIX_VALUE)
                     //           Frame_Expdata_cum[i] = Frame_Expdata_cum[i] / 4;
                     //   }
                        performSubDivisionIM(tempArrayXinvert, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, tempArrSubDivided, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2);
                        performSubDivisionIM(Frame_Expdata_cum, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, tempExpArrSubDivided, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2);


                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc * 2; i++) {
                            mid_X_new = i - uvitObj.IMG_DIM_FIpc;
                            for (int j = 0; j < uvitObj.IMG_DIM_FIpc * 2; j++) {
                                mid_Y_new = j - uvitObj.IMG_DIM_FIpc;

                                x1 = (-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc;
                                x2 = (-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc;
                                //              if(x1>0 && x1<FINALFRAMESIZE_REGAVG && x2>0 && x2<FINALFRAMESIZE_REGAVG){
                                //                  imageLocation_Array_X[(int)(round(x2)*FINALFRAMESIZE_REGAVG+round(x1))]=Regavg_subSampled_Array_Sig[j*FINALFRAMESIZE_REGAVG+i];
                                //              }
                                // if((int)(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1))> 0 && (int)(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1))<uvitObj.IMG_DIM_FIpc * 2*uvitObj.IMG_DIM_FIpc * 2)
                                if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc * 2 && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc * 2) {
                                    if (RADECimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] != INVALID_PIX_VALUE && tempArrSubDivided[j * uvitObj.IMG_DIM_FIpc * 2 + i] != INVALID_PIX_VALUE) {
                                        RADECimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = RADECimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + tempArrSubDivided[j * uvitObj.IMG_DIM_FIpc * 2 + i];
                                        counterArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = counterArray[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + 1;
                                    } else {
                                        RADECimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = INVALID_PIX_VALUE;
                                    }
                                }
                            }
                        }

                        //for exposure Array
                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc * 2; i++) {
                            mid_X_new = i - uvitObj.IMG_DIM_FIpc;
                            for (int j = 0; j < uvitObj.IMG_DIM_FIpc * 2; j++) {
                                mid_Y_new = j - uvitObj.IMG_DIM_FIpc;

                                x1 = (-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + uvitObj.IMG_DIM_FIpc;
                                x2 = (-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + uvitObj.IMG_DIM_FIpc;
                                //              if(x1>0 && x1<FINALFRAMESIZE_REGAVG && x2>0 && x2<FINALFRAMESIZE_REGAVG){
                                //                  imageLocation_Array_X[(int)(round(x2)*FINALFRAMESIZE_REGAVG+round(x1))]=Regavg_subSampled_Array_Sig[j*FINALFRAMESIZE_REGAVG+i];
                                //              }
                                //if((int)(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1))> 0 && (int)(round(x2)*uvitObj.IMG_DIM_FIpc * 2+round(x1))<uvitObj.IMG_DIM_FIpc * 2*uvitObj.IMG_DIM_FIpc * 2)
                                if (round(x1) > 0 && round(x1) < uvitObj.IMG_DIM_FIpc * 2 && round(x2) > 0 && round(x2) < uvitObj.IMG_DIM_FIpc * 2) {
                                    if (RADECExpimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] != INVALID_PIX_VALUE && tempExpArrSubDivided[j * uvitObj.IMG_DIM_FIpc * 2 + i] != INVALID_PIX_VALUE) {
                                        RADECExpimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = RADECExpimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + tempExpArrSubDivided[j * uvitObj.IMG_DIM_FIpc * 2 + i];
                                        counterArray_exp[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = counterArray_exp[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] + 1;
                                    } else {
                                        RADECExpimageSubdivided[(int) (round(x2) * uvitObj.IMG_DIM_FIpc * 2 + round(x1))] = INVALID_PIX_VALUE;
                                    }
                                }

                            }



                        }



                        ApplyBinning(RADECimageSubdivided, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2, RADECimage, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, counterArray);
                        ApplyBinning(RADECExpimageSubdivided, uvitObj.IMG_DIM_FIpc * 2, uvitObj.IMG_DIM_FIpc * 2, RADECExpimage, uvitObj.IMG_DIM_FIpc, uvitObj.IMG_DIM_FIpc, counterArray_exp);

                       // initArray(tempArrayXinvert, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                        //initArray(Frame_Expdata_cum, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc; i++) {
                            for (int j = 0; j < uvitObj.IMG_DIM_FIpc; j++) {
                         //       tempArrayXinvert[i * uvitObj.IMG_DIM_FIpc + j] = RADECimage[i * uvitObj.IMG_DIM_FIpc + ((uvitObj.IMG_DIM_FIpc - 1) - j)];
                           //     Frame_Expdata_cum[i * uvitObj.IMG_DIM_FIpc + j] = RADECExpimage[i * uvitObj.IMG_DIM_FIpc + ((uvitObj.IMG_DIM_FIpc - 1) - j)];
                            }
                        }

                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                            RADECimage[i] = tempArrayXinvert[i];
                            RADECExpimage[i] = Frame_Expdata_cum[i];
                        }

                       
                        delete[] tempArrSubDivided;
                        delete[] tempExpArrSubDivided;
                        delete[] RADECimageSubdivided;
                        delete[] RADECExpimageSubdivided;

                        float *noise_Map = new float[uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc];
                        initArray(noise_Map, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, 0.0f);
                        for (int i = 0; i < uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc; i++) {
                            if (Frame_Expdata_cum[i] * Final_integrationTime != 0.0f && tempArrayXinvert[i] != INVALID_PIX_VALUE && Frame_Expdata_cum[i] != INVALID_PIX_VALUE)
                                noise_Map[i] = sqrt(tempArrayXinvert[i] * Frame_Expdata_cum[i] * Final_integrationTime) / (Frame_Expdata_cum[i] * Final_integrationTime);
			    else{
				noise_Map[i]=0.0f;
			}	

                        }

                        status = 0;
                        // long fpixel[2] ;
                        fpixel[0] = fpixel[1] = 1;
                        // char outfile[NAMESIZE] ;
                        // long naxes[2] ;
                        naxes[0] = naxes[1] = uvitObj.IMG_DIM_FIpc;
                        naxis = 2;
                        bitpix = FLOAT_IMG;
                        //string  finalImagePath;
                        finalImagePath = level2_Origname + "FUV_Final_" + curr_filtercal_val + "_W" + (string) window_index_track[window_index];


                        LOG(INFO) << finalImagePath;
                        if (DirExists((char*) finalImagePath.c_str())) {
                            LOG(ERROR) << "Directory exists and clobber=yes";
                            cmd = (string) "rm -rf " + (string) finalImagePath;
                            system(cmd.c_str());
                        }

                        cmd = "mkdir -p " + (string) finalImagePath;
                        system(cmd.c_str());

                        //    char finalfilename[NAMESIZE];
                        sprintf(finalfilename, "%s/%s_W%s_%s", (char*) finalImagePath.c_str(),curr_filtercal_val,(char*)window_index_track[window_index].c_str(), "FinalImage_Sig.fits");
                        sprintf(finalExpfilename, "%s/%s_W%s_%s", (char*) finalImagePath.c_str(),curr_filtercal_val,(char*)window_index_track[window_index].c_str(), "FinalImage_Exp.fits");
                        sprintf(finalNoiseMapfilename, "%s/%s_W%s_%s", (char*) finalImagePath.c_str(),curr_filtercal_val,(char*)window_index_track[window_index].c_str(), "FinalImage_NoiseMap.fits");





                        //fitsfile *fout ;
                        fits_create_file(&fout, finalfilename, &status);
                        printError(status, "Error in creating the output Signal File", outfile);
                        fits_create_img(fout, bitpix, naxis, naxes, &status);
                        printError(status, "Error in Creating the image for Signal Fie", outfile);
                        fits_write_pix(fout, TFLOAT, fpixel, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, tempArrayXinvert, &status);
                        printError(status, "***Error in writing the pixels to output***", outfile);
string minNoStars_str,threshodpc_str;
			threshodpc_str="Threshold Used for DriverModule-multiplier to Sigma ="+(string)convertFloatToStr(uvitObj.sd_multi_factor_defaultpc);			
			minNoStars_str="Minimum number of star requirement in finding star in DriverModule="+(string)convertIntToStr(uvitObj.min_num_stars);
 			 fits_write_history (fout, threshodpc_str.c_str(), &status);
 			fits_write_history (fout, minNoStars_str.c_str(), &status);
                        for (int i = 0; i < track_headerInfoL1.size(); i++) {
                            naxes[0] = naxes[1] = 0;
                            if (i > 0) {
                                fits_open_file(&fout, finalfilename, READWRITE, &status);
                                fits_movabs_hdu(fout, i + 1, NULL, &status); //previous HDU
                            }
                            fits_create_img(fout, bitpix, naxis, naxes, &status);
                            printError(status, "Error in Creating the image for Signal Fie", outfile);
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                            fits_open_file(&fout, finalfilename, READWRITE, &status);
                            fits_movabs_hdu(fout, i + 2, NULL, &status);
                            for (int j = 0; j < track_headerInfoL1[i].size(); j++) {
                                if (strstr(track_headerInfoL1[i][j].c_str(), "NAXIS") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "BITPIX") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "NAXIS1") == NULL
                                        && strstr(track_headerInfoL1[i][j].c_str(), "NAXES2") == NULL) {
                                    fits_write_record(fout, track_headerInfoL1[i][j].c_str(), &status);
                                }
                            }
			   // fits_write_record(fout, " Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV ", &status);
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                        }


                        //fits_close_file(fout, &status);
                        //printError(status, "Error in closing the  output Signal fits file", outfile);
                        naxes[0] = naxes[1] = uvitObj.IMG_DIM_FIpc;
                        fits_create_file(&fout, finalExpfilename, &status);
                        printError(status, "Error in creating the output Signal File", outfile);
                        fits_create_img(fout, bitpix, naxis, naxes, &status);
                        printError(status, "Error in Creating the image for Signal Fie", outfile);
                        fits_write_pix(fout, TFLOAT, fpixel, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, Frame_Expdata_cum, &status);
                        printError(status, "***Error in writing the pixels to output***", outfile);
			 fits_write_history (fout, threshodpc_str.c_str(), &status);
 			fits_write_history (fout, minNoStars_str.c_str(), &status);
                        for (int i = 0; i < track_headerInfoL1.size(); i++) {
                            naxes[0] = naxes[1] = 0;
                            if (i > 0) {
                                fits_open_file(&fout, finalExpfilename, READWRITE, &status);
                                fits_movabs_hdu(fout, i + 1, NULL, &status); //previous HDU
                            }
                            fits_create_img(fout, bitpix, naxis, naxes, &status);
                            printError(status, "Error in Creating the image for Signal Fie", outfile);
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                            fits_open_file(&fout, finalExpfilename, READWRITE, &status);
                            fits_movabs_hdu(fout, i + 2, NULL, &status);
                            for (int j = 0; j < track_headerInfoL1[i].size(); j++) {
                                if (strstr(track_headerInfoL1[i][j].c_str(), "NAXIS") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "BITPIX") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "NAXIS1") == NULL
                                        && strstr(track_headerInfoL1[i][j].c_str(), "NAXES2") == NULL) {
                                    fits_write_record(fout, track_headerInfoL1[i][j].c_str(), &status);
                                }
                            }
 			//fits_write_record(fout, " Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV ", &status);
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                        }

                        // fits_close_file(fout, &status);
                        //printError(status, "Error in closing the  output Signal fits file", outfile); 
                        naxes[0] = naxes[1] = uvitObj.IMG_DIM_FIpc;
                        fits_create_file(&fout, finalNoiseMapfilename, &status);
                        printError(status, "Error in creating the output Signal File", outfile);
                        fits_create_img(fout, bitpix, naxis, naxes, &status);
                        printError(status, "Error in Creating the image for Signal Fie", outfile);
                        fits_write_pix(fout, TFLOAT, fpixel, uvitObj.IMG_DIM_FIpc * uvitObj.IMG_DIM_FIpc, noise_Map, &status);
                        printError(status, "***Error in writing the pixels to output***", outfile);
			 fits_write_history (fout, threshodpc_str.c_str(), &status);
 			fits_write_history (fout, minNoStars_str.c_str(), &status);
                        for (int i = 0; i < track_headerInfoL1.size(); i++) {
                            naxes[0] = naxes[1] = 0;
                            if (i > 0) {
                                fits_open_file(&fout, finalNoiseMapfilename, READWRITE, &status);
                                fits_movabs_hdu(fout, i + 1, NULL, &status); //previous HDU
                            }
                            fits_create_img(fout, bitpix, naxis, naxes, &status);
                            printError(status, "Error in Creating the image for Signal Fie", outfile);
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                            fits_open_file(&fout, finalNoiseMapfilename, READWRITE, &status);
                            fits_movabs_hdu(fout, i + 2, NULL, &status);
                            for (int j = 0; j < track_headerInfoL1[i].size(); j++) {
                                if (strstr(track_headerInfoL1[i][j].c_str(), "NAXIS") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "BITPIX") == NULL && strstr(track_headerInfoL1[i][j].c_str(), "NAXIS1") == NULL
                                        && strstr(track_headerInfoL1[i][j].c_str(), "NAXES2") == NULL) {
                                    fits_write_record(fout, track_headerInfoL1[i][j].c_str(), &status);
                                }
                            }
			 //fits_write_record(fout, " Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV ", &status);
                            fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", outfile);
                        }

                        //fits_close_file(fout, &status);
                        //printError(status, "Error in closing the  output Signal fits file", outfile); 

                        delete[] RADECimage;
                        delete[] noise_Map;
if(Sum_Exp_Time <0) Sum_Exp_Time=INVALID_PIX_VALUE;
                        writeUsrkeywordsFrmvectDriver(finalfilename, FirstHDU_headertrack);
                       
                       float crpix1=uvitObj.IMG_DIM_FIpc/2;
float crpix2=uvitObj.IMG_DIM_FIpc/2;
float crota1=0.0f,cdelt1=0.0f,cdelt2=0.0f;
float crota2=0.0f;
 int factor_delta=uvitObj.IMG_DIM_FIpc/600;
 
    
        //cdelt1=3.357/3600;
       //  cdelt2=3.311/3600;
       cdelt1=cdelt2=(3.3373/3600)/factor_delta;
  
    
   // cdelt1=-cdelt1/cos(center_dec*M_PI/180);
     cdelt1=-cdelt1;///cos(center_dec*M_PI/180); 
                        
     float rapnt,decpnt;
                        fits_open_file(&fout, finalfilename, READWRITE, &status);
                        printError(status, "Error in Creating the image for Signal Fie", finalfilename);
                         fits_update_key (fout , TFLOAT , "EXP_TIME" , &Sum_Exp_Time , NULL , &status) ;
                         printError(status, "Error in updating keyword value of EXP_TIME", finalfilename);
                         //fits_read_key (fout , TFLOAT , "RA_PNT" , &rapnt , NULL , &status) ;
                         //fits_read_key (fout , TFLOAT , "DEC_PNT" , &decpnt , NULL , &status) ;
    fits_update_key(fout, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(fout, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL1", &RA_pnt, "", &status);                                              printError(status,"");  
    fits_update_key(fout, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(fout, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL2", &DEC_pnt, "", &status);                                    printError(status,"");   
    fits_update_key(fout, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(fout, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
                         for (int i=0;i<IndividualOrbit_Dir.size();i++){
                             fits_write_history(fout, IndividualOrbit_Dir[i].c_str(), &status);
                             
                         }
                         fits_write_history(fout,"Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV",&status);
                         fits_close_file(fout, &status);
                          printError(status, "Error in closing the  output Signal fits file", finalfilename);
                        
                         writeUsrkeywordsFrmvectDriver(finalExpfilename, FirstHDU_headertrack);
                         fits_open_file(&fout, finalExpfilename, READWRITE, &status);
                         printError(status, "Error in Creating the image for Signal Fie", finalExpfilename);
                         fits_update_key (fout , TFLOAT , "EXP_TIME" , &Sum_Exp_Time , NULL , &status) ;
                         printError(status, "Error in updating keyword value of EXP_TIME", finalExpfilename);
                         fits_update_key(fout, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(fout, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL1", &RA_pnt, "", &status);                                              printError(status,"");  
    fits_update_key(fout, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(fout, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL2", &DEC_pnt, "", &status);                                    printError(status,"");   
    fits_update_key(fout, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(fout, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
    
                        for (int i=0;i<IndividualOrbit_Dir.size();i++){
                             fits_write_history(fout, IndividualOrbit_Dir[i].c_str(), &status);
                             
                         }
			fits_write_history(fout,"Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV",&status);
                         fits_close_file(fout, &status);
                         printError(status, "Error in closing the  output Signal fits file", finalExpfilename);
                        
                        
                        writeUsrkeywordsFrmvectDriver(finalNoiseMapfilename, FirstHDU_headertrack);
                        fits_open_file(&fout, finalNoiseMapfilename, READWRITE, &status);
                        printError(status, "Error in Creating the image for Signal Fie", finalNoiseMapfilename);
                         fits_update_key (fout , TFLOAT , "EXP_TIME" , &Sum_Exp_Time , NULL , &status) ;
                         printError(status, "Error in updating keyword value of EXP_TIME", finalNoiseMapfilename);
    fits_update_key(fout, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(fout, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL1", &RA_pnt, "", &status);                                              printError(status,"");  
    fits_update_key(fout, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(fout, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(fout, TDOUBLE, "CRVAL2", &DEC_pnt, "", &status);                                    printError(status,"");   
    fits_update_key(fout, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(fout, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
                        for (int i=0;i<IndividualOrbit_Dir.size();i++){
                             fits_write_history(fout, IndividualOrbit_Dir[i].c_str(), &status);
                             
                         }
			fits_write_history(fout,"Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV",&status);
                         fits_close_file(fout, &status);
                            printError(status, "Error in closing the  output Signal fits file", finalNoiseMapfilename);
                        
                        LOG(INFO) << "NOW FULL FRAMEASTROMETRY FOR  COMBINED ORBIT  IMAGES FOR NUV";

                        //string outputFullFrameAst;
                        outputFullFrameAst = level2_Origname + "FUV_FullFrameAst_" + curr_filtercal_val+"_W" + (string) window_index_track[window_index];;
                        LOG(INFO) << outputFullFrameAst;
                        if (DirExists((char*) outputFullFrameAst.c_str())) {
                            LOG(ERROR) << "Directory exists and clobber=yes";
                            cmd = (string) "rm -rf " + (string) outputFullFrameAst;
                            system(cmd.c_str());
                        }

                        cmd = "mkdir -p " + (string) outputFullFrameAst;
                        system(cmd.c_str());
 delete[] tempArrayXinvert;
                        //string filenamesnr;//this is not need here because we have an image.
                        uvtFullFrameAst ast_obj1;
                        ast_obj1.read((char*) finalImagePath.c_str(),finalfilename, (char*) uvitObj.caldbindir.c_str(), (char*) filenamesnr.c_str(),finalNoiseMapfilename,finalExpfilename, (char*) outputFullFrameAst.c_str(), (char*) ATTfile, (char*) uvitObj.att_timecolpc, (char*) uvitObj.att_qcolpc, (char*) uvitObj.caldbindir.c_str(),
                                sd_multfactor, uvitObj.minimum_No_of_Starspc, uvitObj.refine_Winsizepc / 8, uvitObj.centroid_Winsizepc / 8, uvitObj.databasenamepc, uvitObj.search_algo_ctlgpc,
                                uvitObj.len_apc, uvitObj.len_bpc, uvitObj.rad_searchpc, uvitObj.clobberpc, uvitObj.historypc, 1);
                        ast_obj1.display();


                        status = ast_obj1.uvtFullFrmAstProcess();
                        if (status) {
                            LOG(ERROR) << endl << "Error in full frame astrometry  module";
                           continue;
                        }
                    }
                }
            }
        }

        google::ShutdownGoogleLogging();
string cmdToRemoveTar="rm -r "+newl1tar;
LOG(INFO)<<"Executing command "<<cmdToRemoveTar;
system(cmdToRemoveTar.c_str());
        return (EXIT_SUCCESS);
    }


}

int writeOutputImageToDisk(char *id, char *outDir, char *dir, char *subscript, float *Array, char *namepre, double ftime, unsigned short fno, int sizex, int sizey) {
    int status = 0;
    long fpixel[2];
    fpixel[0] = fpixel[1] = 1;
    char outfile[NAMESIZE];
    long naxes[2];
    naxes[0] = naxes[1] = sizex;
    int naxis = 2;
    int bitpix = FLOAT_IMG;
    fitsfile *fout;
    sprintf(outfile, "%s/%s/%s_t%.4f_f%d_%s_%s.fits", outDir, dir, namepre, ftime, fno, subscript, id);
    int numhdu = 0;


    fits_create_file(&fout, outfile, &status);
    printError(status, "Error in creating the output Signal File", outfile);
    fits_create_img(fout, bitpix, naxis, naxes, &status);
    printError(status, "Error in Creating the image for Signal Fie", outfile);
    fits_write_pix(fout, TFLOAT, fpixel, sizex*sizey, Array, &status);
    printError(status, "***Error in writing the pixels to output***", outfile);


    fits_close_file(fout, &status);
    printError(status, "Error in closing the  output Signal fits file", outfile);

    return (EXIT_SUCCESS);
}

int UVIT_DriverModule::calculateShift(float* array) {




    return (EXIT_SUCCESS);
}

int UVIT_DriverModule::findStar_algo1(float *inputArray, int channel_id) //algorithm for finding the peaks
{
    Fx.clear();
    Fy.clear();
    Fval.clear();
    Rx.clear();
    Ry.clear();
    Rval.clear();
    Cx.clear();
    Cy.clear();
    Ci.clear();
    int r, c;
    float *temp_array;
    vector<float> array_temp;
    //LOG(INFO) << "INSIDE " << IMG_DIM_FIpc;
    long subDivision_sizerapc = this->IMG_DIM_FIpc;
    if (channel_id == PC) {

        // LOG(INFO) << "INSIDE ";

        array_temp.clear();

        for (int i = 0; i < subDivision_sizerapc * subDivision_sizerapc; i++) {

            if (inputArray[i] != 0.0f) {
                array_temp.push_back(inputArray[i]);

            }
        }
        temp_array = new float[array_temp.size()];
        for (int in = 0; in < array_temp.size(); in++) {
            temp_array[in] = array_temp[in];

        }
    }

label:
    Fval.clear();
    Fx.clear();
    Fy.clear();
    Rx.clear();
    Ry.clear();
    Rval.clear();


    // LOG(INFO)<<sd_mul_factor<<endl;exit(1);
    if (sd_multi_factor_defaultpc < 0) {
        LOG(ERROR) << "***SD_MULTI_FACTOR is <0***";
        return (EXIT_FAILURE);
    }
    double thr = 0;
    float mean_Ofimage, sd_temp;
    ;
    if (channel_id == PC) {
        // LOG(INFO) << "INSIDE ";
        sd_temp = getSD(temp_array, array_temp.size());
        mean_Ofimage = getmean(temp_array, array_temp.size());
        //    thr =  sd_temp* sd_mul_factor ;
        thr = mean_Ofimage + sd_temp* sd_multi_factor_defaultpc;
    } else {
        sd_temp = getSD(inputArray, subDivision_sizerapc * subDivision_sizerapc);
        mean_Ofimage = getmean(inputArray, subDivision_sizerapc * subDivision_sizerapc);
        //    thr =  sd_temp* sd_mul_factor ;
        thr = mean_Ofimage + sd_temp* sd_multi_factor_defaultpc;

    }

    //LOG(ERROR) << "Threshold for first cut peaks is   " << mean_Ofimage << " + " << sd_temp << " X " << sd_multi_factor_defaultpc << " = " << thr;
    //  thr=0;//for time being.
    for (int i = 0; i < subDivision_sizerapc * subDivision_sizerapc; i++) {
        r = (i / subDivision_sizerapc);
        c = (i % subDivision_sizerapc);

        if (inputArray[i] > thr) {
            Fval.push_back(inputArray[i]);
            Fx.push_back(c + 1); //x is for column
            Fy.push_back(r + 1); //y is for row
        }

    }

    //LOG(INFO) << "SIGMA  Factor::" << sd_multi_factor_defaultpc;
    // LOG(INFO) << " Size of First cut Peaks  " << Fy.size();

    if (Fy.size() < min_num_stars) {
        sd_multi_factor_defaultpc = sd_multi_factor_defaultpc - 0.25;
        if (sd_mul_factorpc < 0) {
            LOG(ERROR) << sd_multi_factor_defaultpc << " less than 0!!!! ";
            return (EXIT_FAILURE);
        }
        goto label;

    }




    //if winsize is even, make it odd
    if (refine_Winsizepc % 2 == 0)
        refine_Winsizepc = refine_Winsizepc - 1;

    //LOG(INFO) << "Using window size : " << refine_Winsizepc << " for refining peaks ";

    //refined peaks
    vector<int> Tx, Ty;
    vector<float> Tval;

    Tx.reserve(Fx.size());
    Ty.reserve(Fy.size());
    Tval.reserve(Fval.size());


    Tx = Fx;
    Ty = Fy;
    Tval = Fval;
    Star star1;
    star_track.clear();
    star_track.reserve(Tx.size());
    for (int i = 0; i < Tx.size(); i++) {
        star1.x = Tx[i];
        star1.y = Ty[i];
        star1.intensity = Tval[i];
        star_track.push_back(star1);
        //LOG(INFO)<<Tx[i]<<" "<<Ty[i]<<" "<<Tval[i]<<endl;
    }


    sort(star_track.begin(), star_track.end(), compare1);
    //LOG(INFO)<<"sorting  finished";
    /*refining peaks logic
    refined Window size is for the refined  peaks.
   Refined peaks are found by  making window around each of the star(i.e first cut peaks)  and  finding brightest star among that window.*/
    int start_r, end_r, start_c, end_c;
    //to be removed 
    Tx.clear();
    Ty.clear();
    Tval.clear();
    //  bool flag_unique=FALSE;

    // vector<Star> ::iterator itr =star_track.begin ();
    for (int i = 0; i < star_track.size(); i++) {
        start_r = star_track[i].y - refine_Winsizepc / 2;
        end_r = star_track[i].y + refine_Winsizepc / 2;
        start_c = star_track[i].x - refine_Winsizepc / 2;
        end_c = star_track[i].x + refine_Winsizepc / 2;
        if (start_r < 0) start_r = 0;
        if (end_r >= subDivision_sizerapc) end_r = subDivision_sizerapc - 1;
        if (start_c < 0) start_c = 0;
        if (end_c >= subDivision_sizerapc) end_c = subDivision_sizerapc - 1;
        for (int fcpeak = i + 1; fcpeak < star_track.size(); fcpeak++) {
            if (star_track[fcpeak].x > start_c && star_track[fcpeak].x < end_c && star_track[fcpeak].y > start_r && star_track[fcpeak].y < end_r) {

                star_track.erase(star_track.begin() + fcpeak);
                fcpeak--;
            }
        }
        Tx.push_back(star_track[i].x);
        Ty.push_back(star_track[i].y);
        Tval.push_back(star_track[i].intensity);
    }


    /*--------------Refining peaks completed----------------*/

    float *arr_refine = new float[subDivision_sizerapc * subDivision_sizerapc]; //to store refined peaks

    for (int i = 0; i < subDivision_sizerapc * subDivision_sizerapc; i++)
        arr_refine[i] = 0;

    if (subDivision_sizerapc == 0 || subDivision_sizerapc == 0) {
        LOG(ERROR) << "***Divide by Zero***" << endl;
        return (EXIT_FAILURE);
    }
    for (int i = 0; i < Ty.size(); i++) {
        if ((Ty[i] * subDivision_sizerapc + Tx[i]) > 0 && (Ty[i] * subDivision_sizerapc + Tx[i]) < subDivision_sizerapc * subDivision_sizerapc)
            arr_refine[Ty[i] * subDivision_sizerapc + Tx[i]] = Tval[i]; //overwriting the same place..
    }


    Tx.clear();
    Ty.clear();
    Tval.clear();

    for (int i = 0; i < subDivision_sizerapc * subDivision_sizerapc; i++) {
        if (arr_refine[i] != 0) {
            Rx.push_back((i % subDivision_sizerapc));
            Ry.push_back((i / subDivision_sizerapc));
            Rval.push_back(arr_refine[i]);
        }
    }
    // LOG(INFO) << "Number of final peaks is " << Rval.size();


    if (Ry.size() < min_num_stars) {
        sd_multi_factor_defaultpc = sd_multi_factor_defaultpc - 0.25;
        if (sd_multi_factor_defaultpc < 0) {
            LOG(ERROR) << sd_multi_factor_defaultpc << " less than 0!!!! ";
            return (EXIT_FAILURE);
        }
        delete[] arr_refine;
        goto label;
    }
    /**method for  find Centroid within the Stars**/
    doCentroiding(Rx, Ry, centroid_Winsizepc, inputArray, subDivision_sizerapc, subDivision_sizerapc);

    delete[] arr_refine;
    return (EXIT_SUCCESS);
}

void UVIT_DriverModule::doCentroiding(vector<int> &X, vector<int> &Y, int centroidwindow, float *arr, int h, int w) {

    Cx.clear();
    Cy.clear();
    Ci.clear();

    float x, y, val = 0;
    //   LOG(INFO) << "THe Centroid Window" << centroidwindow << endl;
    int num = centroidwindow*centroidwindow;
    double sum_x = 0, sum_y = 0, sum = 0;

    //LOG(INFO) << centroidwindow;
    //  LOG(INFO) << "The X.size is " << X.size() << endl;
    for (int i = 0; i < X.size(); i++) {
        sum_x = 0;
        sum_y = 0, sum = 0;
        for (int j = -1 * (centroidwindow / 2); j <= (centroidwindow / 2); j++) {
            for (int k = -1 * (centroidwindow / 2); k <= (centroidwindow / 2); k++) {
                x = X[i];
                y = Y[i];
                // if(j>0 && j<h && k>0 && k<w )
//%#Added ON 20July17#% 
 if((Y[i] + j)>0 && (Y[i] + j)<h && (X[i] + k)>0 && (X[i] + k)<w)
//%#-Till this-20July17#%
{
                val = arr[(Y[i] + j) * w + (X[i] + k)];
                //   else
                //        continue;
                if (val == INVALID_PIX_VALUE) {
                    continue;
                }
                sum_x = sum_x + (x + k) * val;
                sum_y = sum_y + (y + j) * val;
                sum = sum + val;
}
            }
        }
        if (sum <= 0) {
            LOG(INFO) << "Sum of intensites for (" << X[i] << " , " << Y[i] << ")  is <=0" << val;
            LOG(INFO) << "Divide by zero error";
            continue;
        }
        Cx.push_back((float) sum_x / (float) sum);
        Cy.push_back((float) sum_y / (float) sum);
        Ci.push_back((float) sum);

    }
    float temp, tx, ty;
    for (int i = 0; i < Ci.size(); i++) {

        //LOG(INFO)<<endl<<i;
        for (int j = Ci.size() - 1; j > i; j--) {
            if (Ci[j - 1] < Ci[j]) {
                swap1(Ci[j], Ci[j - 1]);
                swap1(Cx[j], Cx[j - 1]);
                swap1(Cy[j], Cy[j - 1]);
            }
        }
    }

}

//int UVIT_DriverModule::matchStars(int numrowsFirstfile, int numrowsSecfile, float divFact,
//        float *xlocFirst, float *ylocFirst, float *xlocSec, float *ylocSec, vector<float> &matchPixelXone,
//        vector<float> &matchPixelYone, vector<float> &matchPixelXtwo, vector<float> &matchPixelYtwo,
//        vector<float> &matchPixelDiffX, vector<float> &matchPixelDiffY) {
//
//    int cnt = 0;
//    vector<float> xtrack, ytrack;
//    float temp_x1 = 0, temp_x2 = 0, temp_y1 = 0, temp_y2 = 0;
//
//
//    // cout<<numrowsFirstfile<<" "<<numrowsSecfile<<" "<<divFact<<endl;
//    for (long int index = 0; index < numrowsFirstfile; index++) //loop for finding the similar x-coordinates and y-coordinates  from the both frames
//    {
//        //cout<<index<<endl;
//        if (xlocFirst[index] != -9999 && ylocFirst[index] != -9999) {
//            temp_x1 = xlocFirst[index] / divFact;
//            temp_y1 = ylocFirst[index] / divFact;
//            for (long int i = 0; i < numrowsSecfile; i++) {
//                if (xlocSec[i] != -9999 && ylocSec[i] != -9999) {
//                    float diff_x = 0.0, diff_y = 0.0;
//                    temp_x2 = xlocSec[i] / divFact;
//                    temp_y2 = ylocSec[i] / divFact;
//
//                    vector<float>::iterator xdupli = find(xtrack.begin(), xtrack.end(), temp_x2);
//                    vector<float>::iterator ydupli = find(ytrack.begin(), ytrack.end(), temp_y2);
//
//
//                    if (xdupli != xtrack.end() && ydupli != ytrack.end()) {
//                        continue;
//                    }
//                    diff_x = temp_x2 - temp_x1;
//                    diff_y = temp_y2 - temp_y1;
//                    //                    LOG(INFO) << diff_Distpc;
//                    if ((diff_x * diff_x) < diff_Distpc && (diff_y * diff_y) < diff_Distpc) // finding the similar points
//                    {
//                        xtrack.push_back(temp_x2);
//                        ytrack.push_back(temp_y2);
//                        matchPixelXone.push_back(temp_x1);
//                        matchPixelYone.push_back(temp_y1);
//                        matchPixelXtwo.push_back(temp_x2);
//                        matchPixelYtwo.push_back(temp_y2);
//                        matchPixelDiffX.push_back(diff_x);
//                        matchPixelDiffY.push_back(diff_y);
//
//                        cnt++;
//                        break;
//
//                    }
//                }
//            }
//        }
//    }
//    return cnt;
//}


int UVIT_DriverModule::matchStars(int numrowsFirstfile, int numrowsSecfile, float divFact,
        float *xlocFirst, float *ylocFirst, float *xlocSec, float *ylocSec,double &DiffX,double &DiffY,bool &flag_TestTwo ,bool &flag_TestOne ) {

    DiffX=0.0f;
    DiffY=0.0f;
    int cnt = 0;
    vector<float> xtrack, ytrack;
    float temp_x1 = 0, temp_x2 = 0, temp_y1 = 0, temp_y2 = 0;
    float X_predicted_loc,Y_predicted_loc;
    float nextFrame_X,nextFrame_Y;
   flag_TestTwo=FALSE;
for (int i=0;i<numrowsFirstfile;i++){
LOG(INFO)<<xlocFirst[i]<<" "<<ylocFirst[i]<<" "<<xlocSec[i]<<" "<<ylocSec[i];
}
LOG(INFO)<<"===================================";
            temp_x1 = xlocFirst[0] / divFact;
            temp_y1 = ylocFirst[0] / divFact;
            
            temp_x2 = xlocSec[0] / divFact;
            temp_y2 = ylocSec[0] / divFact;
            
            float X_diff_Brightest=temp_x2-temp_x1;
            float Y_diff_Brightest=temp_y2-temp_y1;
             flag_TestOne=FALSE;
LOG(INFO)<<X_diff_Brightest<<" "<<Y_diff_Brightest<<" "<< xlocFirst[0]<<" "<<numrowsFirstfile<<" "<<numrowsSecfile;
             for (long int index = 1; index < numrowsFirstfile; index++){
                 X_predicted_loc=(xlocFirst[index] / divFact)+X_diff_Brightest;
                 Y_predicted_loc=(ylocFirst[index] / divFact)+Y_diff_Brightest;
                 for (long int i = 1; i < numrowsSecfile; i++){
                    nextFrame_X= xlocSec[i] / divFact;
                    nextFrame_Y= ylocSec[i] / divFact;
LOG(INFO)<<X_predicted_loc<<" "<<nextFrame_X<<" "<<Y_predicted_loc<<" "<<nextFrame_Y<<" ";
                    if(abs(nextFrame_X-X_predicted_loc)<2 && abs(nextFrame_Y-Y_predicted_loc)<2 ){
                        flag_TestOne=TRUE;
                        break;
                    }
                     
                 }
                 if(flag_TestOne==TRUE){
                     break;
                 }
                 
                 
             }
            
            //test2
            
            if(flag_TestOne==FALSE){
                
                LOG(INFO)<<"TEST 2 is starting";
                temp_x1 = xlocFirst[1] / divFact;
                temp_y1 = ylocFirst[1] / divFact;
            
            temp_x2 = xlocSec[0] / divFact;
            temp_y2 = ylocSec[0] / divFact;
                
            X_diff_Brightest=temp_x2-temp_x1;
            Y_diff_Brightest=temp_y2-temp_y1;
          flag_TestTwo=FALSE;
            
             for (long int index = 0; index < numrowsFirstfile; index++){
                 if(index!=1){
                 X_predicted_loc=(xlocFirst[index] / divFact)+X_diff_Brightest;
                 Y_predicted_loc=(ylocFirst[index] / divFact)+Y_diff_Brightest;
                 for (long int i = 1; i < numrowsSecfile; i++){
                    nextFrame_X= xlocSec[i] / divFact;
                    nextFrame_Y= ylocSec[i] / divFact;
                    if(abs(nextFrame_X-X_predicted_loc)<2 && abs(nextFrame_Y-Y_predicted_loc)<2 ){
                        flag_TestTwo=TRUE;
                        break;
                    }
                     
                 }
                 if(flag_TestTwo==TRUE){
                     break;
                 }
                 
                 }
             }
            
            
                
            }
             if(flag_TestOne==TRUE || flag_TestTwo==TRUE){
                 DiffX=X_diff_Brightest;
                 DiffY=Y_diff_Brightest;
                 
             }
             LOG(INFO)<<"Diff X "<<DiffX<<" "<<"Diff_Y"<<DiffY;
            
            

    // cout<<numrowsFirstfile<<" "<<numrowsSecfile<<" "<<divFact<<endl;
//    for (long int index = 0; index < numrowsFirstfile; index++) //loop for finding the similar x-coordinates and y-coordinates  from the both frames
//    {
//        //cout<<index<<endl;
//        if (xlocFirst[index] != -9999 && ylocFirst[index] != -9999) {
//            temp_x1 = xlocFirst[index] / divFact;
//            temp_y1 = ylocFirst[index] / divFact;
//            for (long int i = 0; i < numrowsSecfile; i++) {
//                if (xlocSec[i] != -9999 && ylocSec[i] != -9999) {
//                    float diff_x = 0.0, diff_y = 0.0;
//                    temp_x2 = xlocSec[i] / divFact;
//                    temp_y2 = ylocSec[i] / divFact;
//
//                    vector<float>::iterator xdupli = find(xtrack.begin(), xtrack.end(), temp_x2);
//                    vector<float>::iterator ydupli = find(ytrack.begin(), ytrack.end(), temp_y2);
//
//
//                    if (xdupli != xtrack.end() && ydupli != ytrack.end()) {
//                        continue;
//                    }
//                    diff_x = temp_x2 - temp_x1;
//                    diff_y = temp_y2 - temp_y1;
//                    //                    LOG(INFO) << diff_Distpc;
//                    if ((diff_x * diff_x) < diff_Distpc && (diff_y * diff_y) < diff_Distpc) // finding the similar points
//                    {
//                        xtrack.push_back(temp_x2);
//                        ytrack.push_back(temp_y2);
//                        matchPixelXone.push_back(temp_x1);
//                        matchPixelYone.push_back(temp_y1);
//                        matchPixelXtwo.push_back(temp_x2);
//                        matchPixelYtwo.push_back(temp_y2);
//                        matchPixelDiffX.push_back(diff_x);
//                        matchPixelDiffY.push_back(diff_y);
//
//                        cnt++;
//                        break;
//
//                    }
//                }
//            }
//        }
//    }
    //return cnt;
      return EXIT_SUCCESS;
}

int UVIT_DriverModule::findShiftsNtheta(int totalelements, vector<float> &Xone, vector<float> &Yone, vector<float> &Xtwo, vector<float> &Ytwo,
        vector<float> &DiffOfX, vector<float> &DiffOfY, bool flag_theta_computation, double &Xdx, double &Ydy, double &Theta) {

    if (totalelements <= 2) {
        flag_theta_computation = 0;
    }
    if (flag_theta_computation == 0) {
        LOG(INFO) << "NO theta value will be calculated";
        vector<double> diff_x_cumm, diff_y_cumm;
        for (int t = 0; t < totalelements; t++) {

            diff_x_cumm.push_back(Xtwo[t] - Xone[t]);
            diff_y_cumm.push_back(Ytwo[t] - Yone[t]);
        }
        Xdx = getmean(diff_x_cumm.data(), diff_x_cumm.size());
        Ydy = getmean(diff_y_cumm.data(), diff_y_cumm.size());
        Theta = 0.0f;
    } else if (shift_N_Rotate_algopc == 1) {
        spMatrix B((totalelements) * 2, 1);
        spMatrix A((totalelements) * 2, 3);
        spMatrix X(3, 1);
        int temp = 0;

        for (int t = 0; t < totalelements * 2; t = t + 2) {
            B(t, 0) = (Xtwo[temp] - Xone[temp]); //+IMAGE_ARRAYSIZE*0.5;
            B(t + 1, 0) = (Ytwo[temp] - Yone[temp]); //+IMAGE_ARRAYSIZE*0.5;

            A(t, 0) = -1.0 * (Yone[temp] - (IMG_DIM_FIpc / 8) / 2);
            A(t, 1) = 1;
            A(t, 2) = 0;

            A(t + 1, 0) = (Xone[temp] - (IMG_DIM_FIpc / 8) / 2);
            A(t + 1, 1) = 0;
            A(t + 1, 2) = 1;

            temp++;
        }

        X.ApplyLeastSquare(A, B);

        Xdx = X(1, 0);
        Ydy = X(2, 0);
        Theta = X(0, 0);
        //X (0 , 0) = X (0 , 0) ;
        // ofptr << p << setw (15) << X (0 , 0)*180 / M_PI << setw (25) << X (1 , 0) << setw (25) << X (2 , 0) << endl ;
    } else if (shift_N_Rotate_algopc == 2) {
        spMatrix A(totalelements * 2, 4);
        spMatrix B(totalelements * 2, 1);
        spMatrix X(4, 1);
        for (int aindex = 0; aindex < totalelements; aindex++) {
            A(2 * aindex, 0) = Xone[aindex] - (IMG_DIM_FIpc / 8) * 0.5;
            A(2 * aindex, 1) = -1.0 * (Yone[aindex] - (IMG_DIM_FIpc / 8) * 0.5);
            A(2 * aindex, 2) = 1.0;
            A(2 * aindex, 3) = 0.0;
            A(2 * aindex + 1, 0) = (Yone[aindex] - (IMG_DIM_FIpc / 8) * 0.5);
            A(2 * aindex + 1, 1) = Xone[aindex] - (IMG_DIM_FIpc / 8) * 0.5;
            A(2 * aindex + 1, 2) = 0.0;
            A(2 * aindex + 1, 3) = 1.0;
            B(2 * aindex, 0) = Xtwo[aindex] - (IMG_DIM_FIpc / 8) * 0.5;
            B(2 * aindex + 1, 0) = Ytwo[aindex] - (IMG_DIM_FIpc / 8) * 0.5;
        }
        X.ApplyLeastSquare(A, B);

        double theta = atan2(X(1, 0), X(0, 0));

        // ofptr << theta * 180 / M_PI << setw (15) << X (0 , 0) << setw (25) << X (1 , 0) << endl ;
        Xdx = X(2, 0);
        Ydy = X(3, 0);
        Theta = theta;
    } else if (shift_N_Rotate_algopc == 3) {

        double a11 = 0.0, a12 = 0.0, a13 = 0.0, a21 = 0.0, a22 = 0.0, a23 = 0.0, a31 = 0.0, a32 = 0.0, a33 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0;

        spMatrix A(3, 3);
        spMatrix B(3, 1);
        spMatrix X(3, 1);

        for (int k = 0; k < totalelements; k++) {
            /* weight used for a star*/
            double w = 1.0;
            //     w=pow( (25.0*photonbk_per_pixel+starmag[k]), 2.0 )/
            //        (  (25.0*photonbk_per_pixel)+(0.09*starmag[k])  );


            /* row 1 */
            double y_mod = Ytwo[k]-((IMG_DIM_FIpc / 8) * 0.5);
            double x_mod = Xtwo[k]-((IMG_DIM_FIpc / 8) * 0.5);

            y_mod = -1.0 * y_mod;
            a11 = a11 + 2.0 * w;
            a12 = 0.0;
            a13 = a13 + (2.0 * w * (y_mod));
            b1 = b1 + 2.0 * (w * DiffOfX[k]);

            /* row 2*/
            a21 = 0.0;
            a22 = a22 + 2.0 * w;
            a23 = a23 + (2.0 * w * (x_mod));
            b2 = b2 + 2.0 * (w * DiffOfY[k]);

            /* row 3 */
            a31 = a31 + (2.0 * w * (y_mod));
            a32 = a32 + (2.0 * w * (x_mod));
            a33 = a33 + (2.0 * w * (pow(x_mod, 2.0) + pow(y_mod, 2.0)));
            b3 = b3 + (2.0 * w * ((x_mod) * DiffOfY[k] + (y_mod) * DiffOfX[k]));

        }

        A(0, 0) = a11;
        A(0, 1) = a12;
        A(0, 2) = a13;
        A(1, 0) = a21;
        A(1, 1) = a22;
        A(1, 2) = a23;
        A(2, 0) = a31;
        A(2, 1) = a32;
        A(2, 2) = a33;
        B(0, 0) = b1;
        B(1, 0) = b2;
        B(2, 0) = b3;

        X.ApplyLeastSquare(A, B);

        // ofptr << X (2 , 0)*180 / M_PI << setw (15) << X (0 , 0) << setw (25) << X (1 , 0) << endl ;
        Xdx = X(0, 0);
        Ydy = X(1, 0);
        Theta = X(2, 0);
    }
    return (EXIT_SUCCESS);
}

int UVIT_DriverModule::removeRecords(vector<float> &Xone, vector<float> &Yone, vector<float> &Xtwo, vector<float> &Ytwo, vector<float> &DiffOfX, vector<float> &DiffOfY, float *ints1, float*ints2,
        vector<float> &newXone, vector<float> &newYone, vector<float> &newXtwo, vector<float> &newYtwo, vector<float> &newDiffX, vector<float> &newDiffY, vector<float> &new_one_ints, vector<float> &new_two_ints) {
    float mean_diff_X = 0.0f, mean_diff_Y = 0.0f;
    float sd_diff_X = 0.0f, sd_diff_Y = 0.0f;
    mean_diff_X = getmean(DiffOfX.data(), DiffOfX.size());
    sd_diff_X = getSD(DiffOfX.data(), DiffOfX.size());
    mean_diff_Y = getmean(DiffOfY.data(), DiffOfY.size());
    sd_diff_Y = getSD(DiffOfY.data(), DiffOfY.size());
    newDiffX.clear(), newDiffY.clear(), newXone.clear();
    newXtwo.clear();
    newYone.clear();
    newYtwo.clear();
    new_one_ints.clear(), new_two_ints.clear();
    vector<float> Distance_vect;
    for (int i = 0; i < DiffOfX.size(); i++) {
        Distance_vect.push_back(sqrt(DiffOfX[i] * DiffOfX[i] + DiffOfY[i] * DiffOfY[i]));

    }
    float mean_dist, sd_dist;
    mean_dist = getmean(Distance_vect.data(), Distance_vect.size());
    sd_dist = getSD(Distance_vect.data(), Distance_vect.size());

    for (int i = 0; i < DiffOfX.size(); i++) {

        if (Distance_vect[i] <= (mean_dist + 3 * sd_dist) && Distance_vect[i]> (mean_dist - 3 * sd_dist)) {
            newXone.push_back(Xone[i]);
            newYone.push_back(Yone[i]);
            newXtwo.push_back(Xtwo[i]);
            newYtwo.push_back(Ytwo[i]);
            newDiffX.push_back(DiffOfX[i]);
            newDiffY.push_back(DiffOfY[i]);
            new_one_ints.push_back(ints1[i]);
            new_two_ints.push_back(ints2[i]);
        }

    }


    return (EXIT_SUCCESS);

}

int UVIT_DriverModule::ApplySubSampling(float* inputarray, int in_xsize, int in_ysize, float* outputarray, int out_xsize, int out_ysize) {

    if (in_xsize % out_xsize != 0) {
        LOG(ERROR) << "Can't sub sampling with  this input array ,input Array size is not matching";
        return (EXIT_FAILURE);
    }
    float devide_term_x = in_xsize / out_xsize;
    float devide_term_y = in_ysize / out_ysize;
    int cnt_win = 0;
    float sum_win = 0.0f;
    int index_finalArray = 0;
    for (int temp_x = 0; temp_x < in_xsize; temp_x = temp_x + devide_term_x) {
        for (int temp_y = 0; temp_y < in_ysize; temp_y = temp_y + devide_term_y) {
            if (index_finalArray > out_xsize * out_ysize) {
                LOG(ERROR) << "Array is out of bound, EXEED to " << out_xsize << " !!!";
                return (EXIT_FAILURE);
            }
            cnt_win = 0;
            sum_win = 0.0f;
            for (int i = temp_y; i < temp_y + devide_term_y; i++) {
                for (int j = temp_x; j < temp_x + devide_term_x; j++) {
                    if (inputarray[j * in_xsize + i] != INVALID_PIX_VALUE) {
                        sum_win = sum_win + inputarray[j * in_xsize + i];
                        cnt_win++;
                    }

                }
            }
            if (cnt_win != 0) {
                outputarray[index_finalArray++] = (float) (sum_win / cnt_win);
            } else {
                outputarray[index_finalArray++] = INVALID_PIX_VALUE;
            }


        }

    }

    return (EXIT_SUCCESS);

}

int writeUsrkeywordsFrmvectDriver(char *filename, vector<string> &strvect) {

    fitsfile *fout;
    int status = 0;
    int numhdu = 0;
    fits_open_file(&fout, filename, READWRITE, &status);
    printError(status, "Error in opening the input File", filename);
    // LOG(INFO)<<"Successfully opened"<<endl;
    fits_get_num_hdus(fout, &numhdu, &status);

    for (int hno = 1; hno <= 1; hno++) {
        fits_movabs_hdu(fout, hno, NULL, &status);

        for (int i = 0; i < strvect.size(); i++) {
            // if (strstr(strvect[i].c_str(), "NAXIS") == NULL && strstr(strvect[i].c_str () ,"BITPIX")==NULL && strstr(strvect[i].c_str () ,"NAXIS1")==NULL 
            //        && strstr(strvect[i].c_str () ,"NAXES2")==NULL   && strstr(strvect[i].c_str(),"SIMPLE")==NULL  && strstr(strvect[i].c_str(),"EXTEND")==NULL   && strstr(strvect[i].c_str(),"DATASUM")==NULL   && strstr(strvect[i].c_str(),"CREATOR")==NULL && strstr(strvect[i].c_str(),"XSIZE")==NULL && strstr(strvect[i].c_str(),"YSIZE")==NULL && strstr(strvect[i].c_str(),"EXP_TIME")==NULL && strstr(strvect[i].c_str(),"NAMEPRFX")==NULL && strstr(strvect[i].c_str(),"CHECKSUM")==NULL) {
            if (strstr(strvect[i].c_str(), "NAXIS") == NULL && strstr(strvect[i].c_str(), "BITPIX") == NULL && strstr(strvect[i].c_str(), "NAXIS1") == NULL
                    && strstr(strvect[i].c_str(), "NAXES2") == NULL && strstr(strvect[i].c_str(), "SIMPLE") == NULL && strstr(strvect[i].c_str(), "EXTEND") == NULL && strstr(strvect[i].c_str(), "DATASUM") == NULL && strstr(strvect[i].c_str(), "CREATOR") == NULL && strstr(strvect[i].c_str(), "CHECKSUM") == NULL) {
                fits_write_record(fout, strvect[i].c_str(), &status);


                //LOG(INFO)<<strvect[i].c_str ()<<endl;
                if (status) {
                    LOG(ERROR) << endl << "***Error in writing keyword " << strvect[i].c_str() << "***";
                    fits_report_error(stderr, status);
                    fits_close_file(fout, &status);
                    printError(status, "Error in closing the input File", filename);
                    //     LOG(INFO)<<"Successfully closed with error"<<endl;
                    return (EXIT_FAILURE);
                }
            }//
        }

    }

    fits_close_file(fout, &status);
    printError(status, "Error in closing the input File", filename);
    // LOG(INFO)<<"Successfully closed"<<endl;
    //            fits_open_file (&fout , filename , READWRITE , &status) ;
    //              printError (status , "Error in opening the input File" , filename) ;
    //                fits_close_file(fout,&status);
    //            printError (status , "Error in closing the input File" , filename) ;

    return (EXIT_SUCCESS);
}

