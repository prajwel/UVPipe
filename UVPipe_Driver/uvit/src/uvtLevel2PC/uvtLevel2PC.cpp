 /* 
 + 27-Oct-2022 : testing a fix to possible BUG in passing arrays
 + while writing FITS tables post-ShiftNRot ...
 +..............................................................
 * File:   uvtPcL2CommanArr.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <stdlib.h>

#include "uvtLevel2PC.h"

#include<Directory.h>
#include<DataInfo.h>
#include<DataIngest.h>
#include <algorithm>
#include<transform.h>
#include<uvtUtils.h>
#include<set>

#include <vector>
#include<memory.h>
#include<spMatrix.h>
#include<glog/logging.h>
#include <bits/basic_string.h>
#include <iomanip>


#define  NBHD_RADIUS 5
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr   //used in doFft() function
#define moduleoutdir_unit "uvtUnitConvertion"
#define moduleoutdir_badpix "uvtMaskBadPix"
#define moduleoutdir_fltfield "uvtFlatFieldCorr"
#define moduleoutdir_qe "uvtQEMCPCorr"
#define moduleoutdir_pixpad "uvtPixPadding"
#define moduleoutdir_cosmicray "uvtCosmicRayCorrection"
#define moduleoutdir_FrameIntegration "uvtFrameIntegration"

#define moduleoutdir_detectordistortion "uvtDetectDistCorr"
#define moduleoutdir_opticaldistortion "uvtOpticDistCorr"


#define moduleoutdir_centCorr "uvtCentroidCorr"
#define moduleoutdir_centBias "uvtCentroidBias"
#define moduleoutdir_shiftNRot "uvtShiftRot"
#define moduleoutdir_findWtdmean "uvtFindWtdMean"
#define moduleoutdir_ravg  "uvtFlippedRegImage"
#define moduleoutdir_ravgflipped "uvtRegAvg"
#define moduleoutdir_RADECIMAGE "uvtRADECImage"
#define moduleoutdir_expoframes  "uvtExposureFrames"
#define moduleoutdir_subDivision "uvtSubDivision"
//#define moduleoutdir_shiftNrot "uvtShiftRot"
//#define IMAGE_ARRAYSIZE  4800


using namespace std ;
//int ApplySubSampling_Addition (float* inputarray , int in_xsize , int in_ysize , float* outputarray , int out_xsize , int out_ysize);
bool Diff_calc (struct New_location star1 , struct New_location star2) ;


bool Diff_calc(struct New_location star1 , struct New_location star2)
{
    return (star1.diff_Distance < star2.diff_Distance) ;
}



bool compare1 (struct Star star1 , struct Star star2) ;


bool compare1 (struct Star star1 , struct Star star2)
{
    return (star1.intensity > star2.intensity) ;
}


uvtLevel2PC::uvtLevel2PC () {
//    tar_extracted_flag_PC=FALSE;
}


uvtLevel2PC::~ uvtLevel2PC() {

 
}


int uvtLevel2PC::readPILParameters (int argc , char** argv)
{
    int status = 0 ;
    char temp[PIL_LINESIZE] ;

    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG (INFO) << "***Error Initializing PIL***" ;
        return status ;
    }
    // if(this->paramfile_Varask_iml2==FALSE){
       
  
    if (PIL_OK != (status = PILGetFname ("level1indir" , temp)))
    {
        LOG (INFO) << endl << "***Error reading input level1 directory name***" ;
        return status ;
    }
    level1indir.assign (temp) ;
    if (PIL_OK != (status = PILGetFname ("caldbdir" , temp)))
    {
        LOG (INFO) << endl << "***Error reading caldb directory name***" ;
        return status ;
    }
    caldbindir.assign (temp) ;
    if (PIL_OK != (status = PILGetFname ("level2outdir" , temp)))
    {
        LOG (INFO) << endl << "***Error reading output level2 directory name***" ;
        return status ;
    }
    level2outdir.assign (temp) ;
   //  }
    if (PIL_OK != (status = PILGetString ("channel" , temp)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    channel.assign (temp) ;
     if (PIL_OK != (status = PILGetBool ("utcFlag" , &UTC_flag)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }
if (PIL_OK != (status = PILGetBool ("crcflag" , &crc_flag)))
    {
        LOG (INFO) << "Error reading parity Flag" << endl ;
        return (status) ;
    }
    
 if (PIL_OK != (status = PILGetBool ("unTar_flag" , &tar_extracted_flag_PC)))
    {
        LOG (INFO) << endl << "***Error Reading GTI flag :" << gti_flag << "***" ;
        return status ;
    }
    if(tar_extracted_flag_PC==0){
          if (PIL_OK != (status = PILGetBool("zipFlag" , &zipFlag)))
    {
        LOG (INFO) << endl << "***Error reading caldb directory name***" ;
        return status ;
    }
    }
    //dataingest param
    if (PIL_OK != (status = PILGetBool ("dropframe" , &dropframe)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("parityFlag" , &parity_flag)))
    {
        LOG (INFO) << "Error reading parity Flag" << endl ;
        return (status) ;
    }
    
//    if (PIL_OK != (status = PILGetInt ("att_flagVal" , &att_flag_val)))
//    {
//        LOG (INFO) << endl << "***Error Reading att  flag value :***" ;
//        return status ;
//    }
    
 if (PIL_OK != (status = PILGetReal("startTime" , &star_time_ForDatase)))
    {
        LOG (INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    
     if (PIL_OK != (status = PILGetReal ("endTime" , &end_time_ForDatase)))
    {
        LOG (INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thr_One_cr" , &first_mult_factor_CR)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thr_Two_cr" , &second_mult_Factor_CR)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    

    if (PIL_OK != (status = PILGetBool ("GTI_FLAG" , &gti_flag)))
    {
        LOG (INFO) << endl << "***Error Reading GTI flag :" << gti_flag << "***" ;
        return status ;
    }
    if (gti_flag)
    {
        all_or_custom = 1 ;
        valid_bit = 1 ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thresholdMultph" , &thr_multiph)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    

//    if (PIL_OK != (status = PILGetInt ("algoFlag" , &star_detect_algo_flag)))
//    {
//        LOG (INFO) << endl << "***Error reading algo_flag ***" ;
//        return status ;
//    }

    if (PIL_OK != (status = PILGetReal4 ("threshold" , &sd_multi_factor_default)))
    {
        LOG (INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetInt ("minimumTargetedstars" , &minimum_No_of_Stars)))
    {
        LOG (INFO) << endl << "***Error reading algo_flag ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("refinedWinSize" , &refine_Winsize)))
    {
        LOG (INFO) << endl << "***Error reading refine window size ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("centroidWinSize" , &centroid_Winsize)))
    {
        LOG (INFO) << endl << "***Error reading centroid  window size ***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetReal4 ("diffDist" , (float*) &diff_Dist)))
    {
        LOG (INFO) << endl << "***Error reading diffDist parameter***" ;
        return status ;
    }
 
   if (PIL_OK != (status = PILGetFname ("catalogpath" , catalogpath)))
            {
                LOG(ERROR) << endl << "\033[1;31m***Error reading output directory name***" ;
                return status ;
            }
     if (PIL_OK != (status = PILGetString ("att_timecol" , att_timecol)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("database_name" , databasename)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("att_qcol" , att_qcol)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("shiftRotDetAlgoFlag" , &shift_N_Rotate_algo)))
    {
        LOG (INFO) << endl << "***Error reading shiftRotDetAlgoFlag***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("frameIntFlag" , &fi_flag)))
    {
        LOG (INFO) << endl << "***Error reading frameIntegration flag***" ;
        return status ;
    }
   
            if (PIL_OK != (status = PILGetInt ("framesDiscard" , &nFrameDiscard_fi)))
            {
                LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
                return status ;
            }
            if (PIL_OK != (status = PILGetInt ("framesCompute" , &nFrameIntegrate_fi)))
            {
                LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
                return status ;
            }
    //}
    if (PIL_OK != (status = PILGetInt ("FrameIntDim" , &IMG_DIM_FI)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
      if (PIL_OK != (status = PILGetInt ("RegAvgfrmsize" , &FINALFRAMESIZE_REGAVG)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
//Added
      if (PIL_OK != (status = PILGetInt ("QEMCP_tobedone" , &qemcpflag)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
    
    if (PIL_OK != (status = PILGetInt ("CentCorr_tobedone" , &centCorrflag)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetInt ("CentBias_tobedone" , &centBiasflag)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
   if (PIL_OK != (status = PILGetInt ("DetectDist_tobedone" , &DetectDistflag)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("OpticDist_tobedone" , &OpticDistflag)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("pathToOutputTar" , outtarpath)))
    {
        LOG (INFO) << endl << "***Error reading output tar path***" ;
        return status ;
    }
    if(!DirExists(outtarpath)){
        LOG(ERROR)<<"**Tar file path not exists!!!***";
        return(EXIT_FAILURE);
    }

//till this    
if(FINALFRAMESIZE_REGAVG>IMG_DIM_FI){
        LOG(ERROR)<<"***Subsampling size can not be less than frame size of frame integration..!!!!***";
                return(EXIT_FAILURE);
    }
//    if (PIL_OK != (status = PILGetInt ("framesCompute" , &nFrameIntegrate_fi)))
//    {
//        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
//        return status ;
//    }
    
//    if (PIL_OK != (status = PILInit (argc , argv)))
//    {
//        LOG (INFO) << "***Error Initializing PIL***" ;
//        return status ;
//    }
      if (PIL_OK != (status = PILGetInt ("search_algo_forFullFrameAst" , &search_algo_ctlg)))
        {
            LOG (INFO) << endl << "***Error Reading search method" << search_algo_ctlg << "***" ;
            return status ;
        }
 if (PIL_OK != (status = PILGetString("Radi_search" , (char*)&rad_search)))
    {
        LOG (INFO) << endl << "***Error reading the radius value***" ;
        return status ;
    }

    if(search_algo_ctlg==1 || search_algo_ctlg==3 || search_algo_ctlg==5)
    {
      if (PIL_OK != (status = PILGetString ("len_rect_a" , (char*)&len_a)))
        {
            LOG (INFO) << endl << "***Error Reading length of rectangle :" <<len_a << "***" ;
            return status ;
        }
    
      if (PIL_OK != (status = PILGetString ("len_rect_b" , (char*)&len_b)))
        {
            LOG (INFO) << endl << "***Error Reading width of rectangle:" << len_b << "***" ;
            return status ;
        }
    
    } 
    
     
     if (PIL_OK != (status = PILGetBool ("flag_thetaComp" , &flag_thetaComp)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobber << "***" ;
        return status ;
    }
    //if(this->paramfile_Varask_iml2==FALSE){
    if (PIL_OK != (status = PILGetFname ("RASfile" , rasfile)))
    {
        LOG (ERROR) << endl << "\033[1;31m***Error reading output directory name***" ;
        return status ;
    }
//     if (PIL_OK != (status = PILGetBool ("LastFileFlag" , &lastFileFlag)))
//    {
//        LOG (INFO) << endl << "***Error Reading history parameter:" << history << "***" ;
//        return status ;
//    }
   // }
    if (PIL_OK != (status = PILGetBool ("clobber" , &clobber)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobber << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("history" , &history)))
    {
        LOG (INFO) << endl << "***Error Reading history parameter:" << history << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("mode" , mode)))
    {
        LOG (INFO) << "***Error Reading mode parameter:" << history << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskuc" , &wtd_uc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing to Disk Flag For unit Conversion:" << "***" ;
        return status ;
    }
    //    if (PIL_OK != (status = PILGetInt ("Write_todiskuc" , &wtd_sd)))
    //    {
    //        LOG(ERROR) << endl << "***Error Reading Writing to Disk Flag For unit Conversion:" << "***" ;
    //        return status ;
    //    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskbp" , &wtd_bp)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the bad pixel module:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskff" , &wtd_ff)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the FlatFieldCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskpp" , &wtd_pp)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentCorr" , &wtd_centCorr)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentBias" , &wtd_centBias)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetInt ("Write_todisksd" , &wtd_sd)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the subDivision:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskcr" , &wtd_cr)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskqe" , &wtd_qemcp)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskfi" , &wtd_fi)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the Accumulated module:" << "***" ;
        return status ;
    }
//    if (PIL_OK != (status = PILGetInt ("Write_todisksc" , &wtd_fsc)))
//    {
//        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the find Star Centroids:" << "***" ;
//        return status ;
//    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskdd" , &wtd_dd)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the Detect Distortion:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskod" , &wtd_od)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the optical Distortion:" << "***" ;
        return status ;
    }
//    if (PIL_OK != (status = PILGetInt ("Write_todisksnr" , &wtd_snr)))
//    {
//        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the optical Distortion:" << "***" ;
//        return status ;
//    }
     PILClose (status) ;
  
    return (EXIT_SUCCESS) ;
}



int uvtLevel2PC::uvtLevel2PCprocess ()
{
    
    
    
     int status = 0 ; //status flag to check return status of functions
    string temp_str;
    
    
    temp_str.assign (level1indir);
   
    if(tar_extracted_flag_PC==FALSE)
    {
    level1indir="";
     if(zipFlag==FALSE)
    {
    status= extractTars (temp_str,level1indir,orbnum);//extract level-1 data  tar file
    if(status)
    {
       LOG(INFO)<<"Error in extracting tar";
LOG(ERROR)<<"CRASH TAR EXTRACTION FAILED (uvtLevel2PC.cpp)";
       return(EXIT_FAILURE);
    }
     }
     else{
           status = extractZip (temp_str , level1indir , orbnum) ;//extract level-1 data  tar file
    
    if (status)
    {
        LOG (INFO) << "Error in extracting tar" ;
LOG(ERROR)<<"CRASH ZIP EXTRACTION FAILED (uvtLevel2PC.cpp)";
        return (EXIT_FAILURE) ;
    } 
         
     }
    }
    else {
        string temp_Dir_name=level1indir+"/uvit/";
        vector<string> orbitNum;
        get_firstLevelSubDirs(temp_Dir_name,orbitNum);
        if(orbitNum.size()<=0){
            LOG(ERROR)<<"NO orbit number named directory  found after scaning";
            return(EXIT_FAILURE);
        }
        orbnum=(string)basename(orbitNum[0].c_str());
        LOG(INFO)<<orbnum;
    }
       
   
    string temp_str_l2=basename(level1indir.c_str ());
    temp_str_l2.replace (30,7,"_level2");
    
    fitsfile *fras;
    fits_open_file (&fras ,rasfile , READONLY , &status) ;
    printError (status , "Error in opening the rasfile" , rasfile) ;
    fits_movabs_hdu (fras , 3 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of  ras file" , rasfile) ;
    fits_get_num_rows (fras , &no_of_records , &status) ;
    printError (status , "Error Reading the number of Rows" , rasfile) ;
//%#Added ON 20July17#% 
if(no_of_records==0){
        LOG(ERROR)<<"***No record found in Drift Series***";
        return(EXIT_FAILURE);
    }
//%#-Till this-20July17#%
    time_drifts = new double[no_of_records];
    fits_read_col (fras , TDOUBLE , 1 , 1 , 1 , no_of_records , NULL , time_drifts , NULL , &status) ;
    printError (status , "***Error reading  centroid x***") ;
    fits_close_file (fras , &status) ;
    printError (status , "Error in closing the file" , rasfile) ;
    //double drift_startTime =time_drifts[0];
    //double drift_stopTime=time_drifts[no_of_records-1];
    //delete[] time_drifts;
    
    
    
    
    
   //creating output directory if it does not exists  
    if (! DirExists ((char *) level1indir.c_str ()))
    {
        LOG(ERROR) << endl << "Input level 1 directory not found " << endl ;
        return (EXIT_FAILURE) ;
    }

    //creating output directory if it does not exists
    if (! DirExists ((char *) level2outdir.c_str ()))
    {
        string cmd = "mkdir -p " + level2outdir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists ((char*) level2outdir.c_str ()) && clobber == YES)
    {
        string cmd = "rm  -r " + level2outdir ;
        system (cmd.c_str ()) ;
        cmd = "mkdir -p " + level2outdir ;
        system (cmd.c_str ()) ;
    }
      else if(DirExists ((char*) level2outdir.c_str ()) && clobber == NO)
    {
       LOG (ERROR) << "***Output Directory " << level2outdir << "    already exists....***" ;
        LOG (ERROR) << "***Use clobber = y for overwriting***" << endl ;
           return(EXIT_FAILURE);
    }
    
    
    
    
    
    
    //    size_rows = new int[avg_Factor];
    //function to get level 1 files
    Directory dirobj ;
    dirobj.setOrbitNo (orbnum);
       int temp_channel_id=-1;
   if(strcmp ((char*)channel.c_str (),"NUV")==0){
       temp_channel_id=NUV;       
               
   }
   else if(strcmp ((char*)channel.c_str (),"FUV")==0){
       temp_channel_id=FUV;
   }
   else if(strcmp ((char*)channel.c_str (),"VIS")==0){
       temp_channel_id=VIS;
   }
   else{
       LOG(INFO)<<"***Invalid channel***";
       return(EXIT_FAILURE);
   }
  
    if (dirobj.setup (level1indir , level2outdir , temp_channel_id))
    { //IM RA  needs VIS data only //setup done only if uvitV directory found
        LOG(ERROR) << endl << "***Error in directory set up***" << endl ;
        LOG(ERROR) << endl << "***This chain runs for NUV/FUV data only...check for uvitN directory in input*** " << endl ;
	LOG(ERROR)<<"CRASH L1 uvtN DIRECTORY ISSUE (uvtLevel2PC.cpp)";
        return (EXIT_FAILURE) ;
    }
    int numofsciencedatafiles = dirobj.sciencedatafile.size () ;//total number of science data file
    if (numofsciencedatafiles <= 00)
    {
        LOG(ERROR) << endl << "No science data file found" << endl ;
        return (EXIT_FAILURE) ;
    }
  //  long timeStarti,timeStopi;
   // float timeStartf,timeStopf;
    string filename_nuvTotal=TOTAL_NUV_FILENAME;
    
    //if(strcmp(channel.c_str(),"NUV")==0){
    if(tar_extracted_flag_PC==TRUE){
       // ofstream of1(filename_nuvTotal.c_str());  
    ofstream of1(filename_nuvTotal.c_str(),ios::app);    
    for(int i=0;i<dirobj.sciencedatafile.size ();i++)
    {
        
        of1<<dirobj.sciencedatafile[i]<<endl;
       
    }
    of1.close();
    }
   // }
//reading TCT file 
    
    
  

    
    //readKeywords((char *) dirobj.sciencedatafile[dataindex].c_str (),1,4,);
    
    
    char obsmode[FLEN_VALUE] ;
    char moduleIndir[FLEN_FILENAME] ; //hold path for input directory for every module, will be updated after every module run
    char outputdir[FLEN_FILENAME] ;
    int division_fact ;

   
    subDivision_size = IMG_DIM_FI ;
    /*==============================================================================================================================*/
     /*+++++====================LOOP FOR NUMBER OF SCIENCE DATA FILE IN LEVEL-1 DIRECTORY=======================================================*/
     /*==============================================================================================================================*/
//method for number of science data file       
    double finalStartTime=0.0f,finalStopTime=0.0f;
  //  LONG  finalStartTime=0,finalStopTime=0;
    fitsfile *sci_data;
    long num_rows_sci_data;
    string clockMaster;
    double time_Overlap=0;
    long nele_tot_l1;
    int vism,nuvm,fuvm;
    float *FinalArray= new float[subDivision_size*subDivision_size];
      vector<string> L1keywords;
      fitsfile *tempfit;
refine_Winsize_track=refine_Winsize;
centroid_Winsize_track=centroid_Winsize;
    for (int dataindex = 0 ; dataindex < numofsciencedatafiles ; dataindex ++)
    {
	flag_Roll_Angleinvalid=FALSE;
	refine_Winsize=refine_Winsize_track;
	centroid_Winsize=centroid_Winsize_track;
        
        Total_exp_time=0.0f;
       
        status=0;
        
         fits_open_file(&tempfit , dirobj.sciencedatafile[dataindex].c_str () , READONLY , &status) ;
         printError (status , "***Error opening  file***" , (char*)dirobj.sciencedatafile[dataindex].c_str ()) ;
         copyUsrkeywrdsTovect(tempfit,L1keywords);
         fits_close_file(tempfit,&status);
         printError (status , "Error in closing  the file" , (char*)dirobj.sciencedatafile[dataindex].c_str ()) ;
        frameIntegrationfrmname.clear();

        LOG(ERROR) << endl << "------------------Data Set " << dataindex + 1 << " : " << dirobj.sciencedatafile[dataindex] << "-----------------------" << endl ;
        /*---finding mode from science data file---*/
        getKeywordVal ((char *) dirobj.sciencedatafile[dataindex].c_str () , "OBS_MODE" , 1 , obsmode) ;
        //readKeywords((char *) dirobj.sciencedatafile[dataindex].c_str (),1,2,TDOUBLE,"TSTART",&finalStartTime,TDOUBLE,"TSTOP",&finalStopTime);
     
        
          status =checkMasterClock((char *) dirobj.sciencedatafile[dataindex].c_str (),clockMaster,vism,fuvm,nuvm,nele_tot_l1);
   if (status)
    {
        LOG (ERROR) << "Error in Finding Master clock for the dataset" ;
        return (EXIT_FAILURE) ;
    }
    fitsfile *ftct ;
    LOG(INFO)<<"Clock Master->"<<clockMaster;

    LOG (INFO) << "Reading TCT file......" ;
    fits_open_file (&ftct , (char*)dirobj.tctfile.c_str() , READONLY , &status) ;
    printError (status , "***Error opening TCT file***" , (char*)dirobj.tctfile.c_str()) ;

    fits_movabs_hdu (ftct , 2 , NULL , &status) ;
    printError (status , "***Error opening TCT file***" , (char*)dirobj.tctfile.c_str()) ;
 int colnum_detector ; //stores column number for columns UVITF, UVITN, UVITV depending on data detector

    
    if (strcasecmp (clockMaster.c_str () , "FUV") == 0) colnum_detector = TCT_FUV_COLNO ;
    else if (strcasecmp (clockMaster.c_str () , "NUV") == 0) colnum_detector = TCT_NUV_COLNO ;
    else if (strcasecmp (clockMaster.c_str () , "VIS") == 0) colnum_detector = TCT_VIS_COLNO ;
    else
    {
        LOG (ERROR) << "***DETECTOR keyword value is  " << datainfo.getDetector () << "\nExpected is FUV, NUV, VIS***" ;
        return (EXIT_FAILURE) ;
    }
   

    long nrow_tct ;
    fits_get_num_rows (ftct , &nrow_tct , &status) ;
    printError (status , "***Error getting rows from TCT file ***" , (char*)dirobj.tctfile.c_str()) ;

    double *uvt_time = new double[nrow_tct] ;
    checkMemoryAvailability (uvt_time , "uvt_time") ;
    double *sps_time = new double[nrow_tct] ;
    checkMemoryAvailability (sps_time , "sps_time") ;

    long firstrow1 = 1 , firstelem1 = 1 ;

    /**reading SPS time column and UVT time column**/
    fits_read_col (ftct , TDOUBLE , 1 , firstrow1 , firstelem1 , nrow_tct , NULL , sps_time , NULL , &status) ;
    printError (status , "***Error  in reading column  value of the sps_time from TCT file ***" , (char*)dirobj.tctfile.c_str()) ;
   fits_read_col (ftct , TDOUBLE , colnum_detector , firstrow1, firstelem1 , nrow_tct , NULL , uvt_time , NULL , &status) ;
    printError (status , "***Error  in reading column  value of the uvt_time from TCT file ***" , (char*)dirobj.tctfile.c_str()) ;
    fits_close_file (ftct , &status) ;
    printError (status , "Error in closing  the TCT File" , (char*)dirobj.tctfile.c_str()) ;
  //tct file reading finished.
    double  drift_startTime,drift_stopTime;
    if(UTC_flag==1){
    LOG(INFO)<<"Drift->"<<setprecision(20)<<time_drifts[0]<<" "<<time_drifts[no_of_records-1];
    for (int i=1;i<nrow_tct;i++)
    {
        if(time_drifts[0]>=sps_time[i-1] && time_drifts[0]<sps_time[i])
        {
            drift_startTime=uvt_time[i-1]*1000;
            break;
        }
        
        
    }
    for (int i=1;i<nrow_tct;i++)
    {
        if(time_drifts[no_of_records-1]>=sps_time[i-1] && time_drifts[no_of_records-1]<sps_time[i])
        {
            drift_stopTime=uvt_time[i-1]*1000;
            break;
        }        
        
    }
    }
    else{
       
        drift_startTime=time_drifts[0]*1000;
        drift_stopTime=time_drifts[no_of_records-1]*1000;
    }
   
    delete[] uvt_time,sps_time;
   
        
        
        
        
        
        
    fits_open_file (&sci_data , (char *) dirobj.sciencedatafile[dataindex].c_str (), READONLY , &status) ;
    printError (status , "Error in opening the rasfile" , rasfile) ;
    fits_movabs_hdu (sci_data , 3 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of  ras file" , rasfile) ;
    fits_get_num_rows (sci_data , &num_rows_sci_data , &status) ;
    printError (status , "Error Reading the number of Rows" , rasfile) ;
    LOG(INFO)<<num_rows_sci_data;
    double *time_scienceData = new  double  [num_rows_sci_data];
     fits_read_col (sci_data , TDOUBLE , 12 , 1 , 1 , num_rows_sci_data , NULL , time_scienceData , NULL , &status) ;
    printError (status , "***Error reading  time from l1 science data file***") ;
    finalStartTime=time_scienceData[0];
    finalStopTime=time_scienceData[num_rows_sci_data-1];
    delete[] time_scienceData;
    fits_close_file(sci_data,&status) ;
    //finalStartTime=timeStarti+timeStartf;
        //finalStopTime=timeStopi+timeStopf;
    //LOG(INFO)<<"Drift_start_time->"<<setprecision(20)<<drift_startTime;
    //LOG(INFO)<<"Drift_End_time->"<<setprecision(20)<<drift_stopTime;
     time_Overlap=abs((drift_stopTime<finalStopTime?drift_stopTime:finalStopTime)-(drift_startTime>finalStartTime?drift_startTime:finalStartTime));
        if(!(((finalStopTime>drift_startTime && finalStopTime <=drift_stopTime ) || 
                (finalStartTime>=drift_startTime && finalStartTime <=drift_stopTime) ||  (drift_startTime>finalStartTime && drift_startTime <=finalStopTime ) ||
               (drift_stopTime>finalStartTime && drift_stopTime <=finalStopTime )) && time_Overlap>20 ))
    {   
        LOG(INFO)<<setprecision(20)<<finalStartTime<<" "<<finalStopTime<<" "<<drift_startTime<<" "<<drift_stopTime;
        LOG(INFO)<<"***Dataset time not matching with RAS file***";

       
       // throw dirobj.sciencedatafile[dataindex];
        
     continue;
    }
    
     
        if (strcasecmp (obsmode , "PC") != 0)
        {
            LOG(ERROR) << endl << "Observation mode is " << obsmode << "  in file  " << dirobj.sciencedatafile[dataindex] ;
            LOG(ERROR) << endl << "Checking next file...." << endl ;
            continue ; //go to next file if obs mode is not IM
        }

        LOG(INFO) << endl << "Data Mode is " << obsmode << endl ;
       
        strcpy (outputdir , dirobj.level2path[dataindex].c_str ()) ;
         //setting path for output directory to be generated 
          sprintf (moduleoutdir_bp , "%s/%s_%s" , outputdir , moduleoutdir_badpix , VERSION) ;
        sprintf (moduleoutdir_uc , "%s/%s_%s" , outputdir , moduleoutdir_unit , VERSION) ;
        sprintf (moduleoutdir_ff , "%s/%s_%s" , outputdir , moduleoutdir_fltfield , VERSION) ;
        sprintf (moduleoutdir_qemcp , "%s/%s_%s" , outputdir , moduleoutdir_qe , VERSION) ;
        sprintf (moduleoutdir_pp , "%s/%s_%s" , outputdir , moduleoutdir_pixpad , VERSION) ;
        sprintf (moduleoutdir_cr , "%s/%s_%s" , outputdir , moduleoutdir_cosmicray , VERSION) ;
        sprintf (moduleoutdir_fi , "%s/%s_%s" , outputdir , moduleoutdir_FrameIntegration , VERSION) ;
        sprintf (moduleoutdir_dd , "%s/%s_%s" , outputdir , moduleoutdir_detectordistortion , VERSION) ;
        sprintf (moduleoutdir_od , "%s/%s_%s" , outputdir , moduleoutdir_opticaldistortion , VERSION) ;
        sprintf (moduleoutdir_centroidCorr , "%s/%s_%s" , outputdir , moduleoutdir_centCorr , VERSION) ;
        sprintf (moduleoutdir_centroidBias , "%s/%s_%s" , outputdir , moduleoutdir_centBias , VERSION) ;
        sprintf (moduleoutdir_snr , "%s/%s_%s" , outputdir , moduleoutdir_shiftNRot , VERSION) ;
        sprintf (moduleoutdir_sd , "%s/%s_%s" , outputdir , moduleoutdir_subDivision , VERSION) ;
        sprintf (moduleoutdir_rav , "%s/%s_%s" , outputdir , moduleoutdir_ravg , VERSION) ;
        sprintf (moduleoutdir_ravFlipped , "%s/%s_%s" , outputdir , moduleoutdir_ravgflipped , VERSION) ;
        sprintf (moduleoutdir_radecimg , "%s/%s_%s" , outputdir , moduleoutdir_RADECIMAGE , VERSION) ;
        sprintf (moduleoutdir_expo , "%s/%s_%s" , outputdir , moduleoutdir_expoframes , VERSION) ;
        strcpy (lbtfile,(char *) dirobj.lbtfile.c_str ());
        //----------DATAINGEST----------//
        LOG(ERROR) << endl << "===================================DATAINGEST===================================================" << endl ;
        DataIngest di_obj ;
        
        di_obj.read ((char *) dirobj.sciencedatafile[dataindex].c_str () , (char*)caldbindir.c_str (),(char *) dirobj.tctfile.c_str () ,(char*)dirobj.mkffile.c_str(),(char *) dirobj.gtifile[dataindex].c_str () ,
                (char *) dirobj.lbtfile.c_str () , (char*)dirobj.attfile.c_str (),(char*) dirobj.darkDirectory.c_str () , att_flag_val,gti_flag , valid_bit , all_or_custom , outputdir , dropframe , parity_flag ,UTC_flag,crc_flag, clobber , history) ;
        di_obj.display () ;
        status = di_obj.DataIngestProcess () ;
        if (status)
        {
           LOG(ERROR)<<  "***Error in Data Ingest Process***" ;
           LOG(ERROR)<<"CRASH DATAINGEST ERROR (uvtLevel2PC.cpp)";
            continue;
        }
      
        strcpy (moduleIndir , di_obj.getModuleOutdir ()) ;
        strcpy (dataIngest_out_dir , di_obj.getModuleOutdir ()) ;
        LOG(INFO) << "Using directory " << moduleIndir << "  as input to uvtUnitConversion" << endl ;

       

//    }
    char caldb_common_dir[PIL_LINESIZE] , caldb_common_dir_qefile[PIL_LINESIZE] ;
    char caldb_temp_dir[PIL_LINESIZE] , caldb_common_dir_ff[FLEN_FILENAME] , caldb_common_dir_cc[FLEN_FILENAME] , caldb_common_dir_optic[FLEN_FILENAME] ;
    strcpy (caldb_temp_dir , (char*) caldbindir.c_str ()) ;
    strcpy (caldb_common_dir , (char*) caldbindir.c_str ()) ;
    strcpy (caldb_common_dir_qefile , (char*) caldbindir.c_str ()) ;
    strcpy (caldb_common_dir_ff , (char*) caldbindir.c_str ()) ;
    strcpy (caldb_common_dir_cc , (char*) caldbindir.c_str ()) ;
    strcpy (caldb_common_dir_optic , (char*) caldbindir.c_str ()) ;
   
     /*Reading information file generated after  DataIngest directory*/
    char infofile_in[PIL_LINESIZE] ;
    string  tempfilepath = searchFile (moduleIndir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG(ERROR) << "***Information file not found in " << moduleIndir << "***"  ;
         continue;
    }
    sprintf (infofile_in , "%s/%s" , moduleIndir , tempfilepath.c_str()) ;
    LOG(INFO) << "Information file :" << infofile_in  ;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (! (FileExists (infofile_in)))
    {
        LOG(ERROR) << endl << "***Input FileList not Found at Specified PATH,Check Input Directory***" ;
        continue;
    }
 vector<string> header_info;
getHistory(header_info);
writeHistory(infofile_in,header_info);

    fitsfile *finfo_in ;
    fits_open_file (&finfo_in , infofile_in , READWRITE , &status) ;
    printError (status , "Error in opening the input information file" , infofile_in) ;
    fits_read_key (finfo_in , TINT , "WIN_X_SZ" , &win_xsize , NULL , &status) ;
    printError (status , "Error in reading the key value of the Window xsize" , infofile_in) ; 
    fits_read_key (finfo_in , TINT , "WIN_Y_SZ" , &win_ysize , NULL , &status) ;
    printError (status , "Error in reading the key value of the Window ysize " , infofile_in) ; 
      fits_read_key (finfo_in , TSTRING , "DATE-OBS" , dateOfObs , NULL , &status) ;
    printError (status , "Error in reading the key value of the DATE-OBS" , infofile_in) ;
        fits_read_key (finfo_in , TSTRING , "TIME-OBS" , timeobs , NULL , &status) ;
    printError (status , "Error in reading the key value of the TIME-OBS " , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU in input information file") ;
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;

    char nameprefix[PIL_LINESIZE] ;
fits_read_key (finfo_in , TDOUBLE , "BZERO_MJD" , &bzero_mjd , NULL , &status) ;
    printError (status , "Error in reading the key value of  BZERO_MJD" , infofile_in) ;
fits_read_key (finfo_in , TDOUBLE , "BSCALE_MJD" , &bscale_mjd , NULL , &status) ;
    printError (status , "Error in reading the key value of  BSCALE_MJD" , infofile_in) ;
fits_read_key (finfo_in , TLONG , "N_TOTAL" , &ntotal, NULL , &status) ;
    printError (status , "Error in reading the key value of  N_TOTAL" , infofile_in) ;
fits_read_key (finfo_in , TLONG , "N_ZERO_CENTROID" , &nzerocentroid , NULL , &status) ;
    printError (status , "Error in reading the key value of  N_ZERO_CENTROID" , infofile_in) ;
    fits_read_key (finfo_in , TDOUBLE , "INTEGRATION_TIME_UTC" , &integrationTime_frmUTC , NULL , &status) ;
    printError (status , "Error in reading the key value of  INTEGRATION_TIME_UTC" , infofile_in) ;
    fits_read_key (finfo_in , TDOUBLE , "INTEGRATION_TIME_CURR_OPTION" , &integrationTime_frm_curr_option , NULL , &status) ;
    printError (status , "Error in reading the key value of the INTEGRATION_TIME_CURR_OPTION " , infofile_in) ;
    fits_read_key (finfo_in , TLONG , "TOTAL INPUT PACKET FROM L1" , &total_InputPacket , NULL , &status) ;
    printError (status , "Error in reading the key value of  TOTAL INPUT PACKET FROM L1" , infofile_in) ;
    fits_read_key (finfo_in , TLONG , "Starting PACKET NUMBER" , &start_pktNo , NULL , &status) ;
    printError (status , "Error in reading the key value of the Starting PACKET NUMBER " , infofile_in) ;
     fits_read_key (finfo_in , TLONG , "CRC_TotalPackets_Failed" , &crcfailedpckts , NULL , &status) ;
    printError (status , "Error in reading the key value of the CRC TotalPackets Failed " , infofile_in) ;
 fits_read_key (finfo_in , TLONG , "PARITY_TotalEvents_Failed" , &parityfailedEvents , NULL , &status) ;
    printError (status , "Error in reading the key value of the PARITY TotalPackets Failed " , infofile_in) ;
fits_read_key (finfo_in , TLONG , "Total remaining packet from L1" , &pcktremained , NULL , &status) ;
    printError (status , "Error in reading the key value of the PARITY TotalPackets Failed " , infofile_in) ;
fits_read_key (finfo_in , TLONG , "Total Events Used" , &eventsremained , NULL , &status) ;
    printError (status , "Error in reading the key value of the PARITY TotalPackets Failed " , infofile_in) ;
fits_read_key (finfo_in , TLONG , "Total Packet from L1" , &pktsFromL1 , NULL , &status) ;
    printError (status , "Error in reading the key value of the PARITY TotalPackets Failed " , infofile_in) ;


    fits_read_key (finfo_in , TLONG , "Total Used packet" , &usedpkts , NULL , &status) ;
    printError (status , "Error in reading the key value of the Total Used packet " , infofile_in) ;
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; 
    fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; 
    fits_read_key (finfo_in , TSTRING , "IMGFILE" , imgfile , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; 
    fits_read_key (finfo_in , TSTRING , "DARKDIR" , darkdir , NULL , &status) ;
    printError (status , "***Error in reading the  key value of the NFILES ***" , infofile_in) ;
    LOG(INFO)<<integrationTime_frmUTC<<" "<<integrationTime_frm_curr_option;
//LOG(INFO)<<crcfailedpckts<<" "<<parityfailedEvents<<" "<<pcktremained<<" "<<eventsremained;exit(1);
       long naxes1[2] ;
    naxes1[0] = naxes1[1] = 0 ;
    int naxis = 2 ;
    int bitpix = FLOAT_IMG ;
    fits_create_img (finfo_in , bitpix , naxis , naxes1 , &status) ;
    printError (status , "Error in Creating the image for Signal Fie" , infofile_in);
    fits_close_file (finfo_in , &status) ;
    printError (status , "Error in reading the column value of the Input Signal List" , infofile_in) ;

    fits_open_file (&finfo_in , infofile_in , READWRITE , &status) ;
    printError (status , "Error in opening the input information file",infofile_in) ;
    fits_movabs_hdu (finfo_in , 3 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU",infofile_in) ;
         for(int i=0;i<L1keywords.size ();i++)
         {   
             if(strstr(L1keywords[i].c_str () ,"NAXIS")==NULL)
                fits_write_record (finfo_in , L1keywords[i].c_str () , &status) ;
         }
        fits_close_file (finfo_in , &status) ;
        printError (status , "Error in closing the information file" , infofile_in) ;
    
  
      
   /*=================================================================================================================================*/
    /*====================================================CALDB READING  STARTED===========================================================*/
    /*=================================================================================================================================*/
    //caldb reading started for QE and MCP file
   string    tempname=caldb_handler.getQEFile (datainfo.getDetector (),datainfo.getObsMode (),caldb_common_dir_qefile);//get filename  of  QE and MCP caldb file 
   if (tempname == " ")
    {
        LOG(ERROR) << "Couldn't find QEMCP file from caldb"  ;
         continue;
    }
    joinStrings (qe_factor_file , 2 , caldb_common_dir_qefile , tempname.c_str()) ;//setting path of QE and MCP caldb file
    
  status= readNumRowsFromFits (qe_factor_file,2,nCalDBTempValues);
   if (status)
            {
                LOG (ERROR) << "Error in reading the number of rows from fits file "<<qe_factor_file ;
                return (EXIT_FAILURE) ;
            }
  
   cal_temp= new float[nCalDBTempValues];   cal_f0 = new float[nCalDBTempValues]; cal_f1=new float[nCalDBTempValues];
   cal_f2=new float[nCalDBTempValues];    cal_f3 = new float[nCalDBTempValues];  cal_f4=new float[nCalDBTempValues];
   cal_f5=new float[nCalDBTempValues];    cal_f6 = new float[nCalDBTempValues],  cal_f7=new float[nCalDBTempValues];
   
   status=readQEMCPFile (qe_factor_file,datainfo.getDetector (),nCalDBTempValues,cal_temp,cal_f0,cal_f1,cal_f2,cal_f3,cal_f4,cal_f5,cal_f6,cal_f7);
   if (status)
            {
                LOG (ERROR) << "Error in reading the QEMCP caldb file ->"<<qe_factor_file ;
                return (EXIT_FAILURE) ;
            }
   status=readNumRowsFromFits ((char *)dirobj.lbtfile.c_str (),2,nrows_lbt);
   if (status)
   {
            LOG (ERROR) << "Error in reading the number of rows from fits file" ;
            return (EXIT_FAILURE) ;
   }
   
   time_lbt= new double[nrows_lbt];
   insideTemp= new float [nrows_lbt],outsideTemp=new float[nrows_lbt];
   //reading inside and outside temperature from LBT file.
    status = getTemp ((char *)dirobj.lbtfile.c_str (),datainfo.getDetector (),time_lbt,insideTemp,outsideTemp,nrows_lbt) ;
    if (status)
    {
        LOG(ERROR) << "***temperature reading from the lbt file unsuccessful***" ;
        continue;
    }

     int filter_coln ;
    int filternumber ;
    sprintf(filter,"%s",datainfo.getFilter ());
   
    qe_mg_factor = new float[nCalDBTempValues] ;
    if (filter == (string) "F0")
    {
        filternumber = 0 ;
        filter_coln = 1 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = cal_f0[q] ;
    }
    else  if (filter == (string) "F1")
    {
        filternumber = 1 ;
        filter_coln = 2 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = cal_f1[q] ;
    }
    else if (filter == (string) "F2")
    {
        filternumber = 2 ;
        filter_coln = 3 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = cal_f2[q] ;
    }
    else if (filter == (string) "F3")
    {
        filternumber = 3 ;
        filter_coln = 4 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = cal_f3[q] ;
    }
    else if (filter == (string) "F4")
    {
        filternumber = 4 ;
        filter_coln = 5 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = cal_f4[q] ;
    }
    else if (filter == (string) "F5")
    {
        filternumber = 5 ;
        filter_coln = 6 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = cal_f5[q] ;
    }
    else if (filter == (string) "F6")
    {
        filternumber = 6 ;
        filter_coln = 7 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = cal_f6[q] ;
    }
    else if (filter == (string) "F7")
    {
        filternumber = 7 ;
        filter_coln = 8 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = cal_f7[q] ;
    }
    
    else{
         LOG(ERROR)<<"***Invalid filter option*** ";
         continue;
    } 
  
    int tfields = 12 ;
    char *ttype[] = {"PacketSequence " , "FrameCount" , "Time" , "Ix" , "Fx" , "Iy" , "Fy" , "Max-Min" , "Min" ,"UVIT_MASTER_TIME" ,"BAD FLAG" , "MULT PHOTON"} ;
    char *tform[] = {"U" , "U" , "D" , "U" , "E" , "U" , "E" , "B" , "B" ,"D", "B" , "B"} ;
    char *tunit[] = {"" , "" , "" , "" , "" , "" , "" , "" , "" , "" , "",""} ;
    char infile[NAMESIZE] ;

    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;

    status = takeDarkinfo () ; //taking dark frame start path and dark end path from  dataIngest 's dark directory.
    if (status)
    {
        LOG (ERROR) << "***Error getting the information of Darks***" << endl ;
        return (EXIT_FAILURE) ;
    }
    LOG(INFO)<<"Reading Dark frames";
    status = readDarkFrame (dstartpath , darkFramestart_data) ;//reading dark start frame

    status = readDarkFrame (dendpath , darkFrameend_data) ;//reading dark end frame


//reading  caldb for flat field correction
    float *finalFlatFielddata ;

    string  tempname1 = caldb_handler.getFlatFieldFile (datainfo.getDetector () , datainfo.getObsMode () , datainfo.getFilter () , caldb_common_dir_ff) ;
    if(tempname1== " ")
    {
     LOG(ERROR)    <<"***Error in finding the flat field file***";
     continue;
    }
    joinStrings (flatfieldfile , 2 , caldb_common_dir_ff , tempname1.c_str()) ;
    flatfielddata = new float[xsize * ysize] ;
    status = readImage (flatfieldfile , 1 , flatfielddata) ;
    finalFlatFielddata = new float[PIX_PADD_SIZE * PIX_PADD_SIZE] ;
    for (int i = 0 ; i < PIX_PADD_SIZE * PIX_PADD_SIZE ; i ++)
    {
        finalFlatFielddata[i] = 0.0f ;
    }

    /*padding is done to  match xsize of event file to that of caldb file*/
    status = Applypadding (flatfielddata , xsize , ysize , finalFlatFielddata , 600 , 600) ;
    if (status)
    {
        LOG (ERROR) << "***Padding Fails For the flatfield CALDB data***" << endl ;
         continue;
    }

    //reading caldb file for masking bad pixels
    tempname = caldb_handler.getBadPixelFile (datainfo.getDetector () , datainfo.getObsMode () , win_xsize+1 , win_ysize+1 , (char*) caldbindir.c_str ()) ;

    if (tempname== " ")
    {
        LOG(ERROR) << "Couldn't find bad pixel file from calDB" ;
        continue;
    }
    char outfile[FLEN_FILENAME] ;
    char file_in[FLEN_FILENAME] ;
    joinStrings (badpixfile , 1 , tempname.c_str()) ;
    badpixdata = new float[xsize * ysize] ;
    status = readImage (badpixfile , 1 , badpixdata) ;

    //reading  caldb file for centroid Correction
    tempname = caldb_handler.getCentroidEAFile (datainfo.getDetector () , caldb_common_dir_cc) ;
    if (tempname==" ")
    {
        LOG (ERROR) <<"Couldn't find CentroidBias file from caldb" ;
        continue;
    }

    joinStrings (centroidEAfile , 2 , caldb_common_dir_cc , tempname.c_str()) ;
    fitsfile *f_ea ;
    long EA , n_ele ;
    fits_open_file (&f_ea , centroidEAfile , READONLY , &status) ;
    printError (status , "Error in closing the info out file" , centroidEAfile) ;
    fits_movabs_hdu (f_ea , 2 , NULL , &status) ;
    printError (status , "Error in closing the info out file" , centroidEAfile) ;
    fits_get_num_rows (f_ea , &n_ele , &status) ;
    printError (status , "Error in closing the info out file" , centroidEAfile) ;
    fits_read_col (f_ea , TINT , 1 , 1 , 1 , n_ele , NULL , (void *) &EA , NULL , &status) ;
    printError (status , "Error in closing the info out file" , centroidEAfile) ;
    fits_close_file (f_ea , &status) ;
    sprintf (file_in , "%s/%s" , moduleIndir , eventfile) ; //taking event file full path



    LOG (INFO) << " Input Event File " << file_in ;
   long nrows = 0 ;
    char centroid_algo[FLEN_FILENAME] ;
    fitsfile *fevt_in , *fout ;
    
    //reading the event file information
    fits_open_file (&fevt_in , file_in , READWRITE , &status) ;
    printError (status , "Error in opening the input event file" , file_in) ;
    fits_read_key (fevt_in , TSTRING , "CENTROID" , centroid_algo , NULL , &status) ;
    printError (status , "Error Reading the Key value for CENTROID" , file_in ) ;
    status=copyAllheaderKeys(file_in);
     if (status)
        {
            LOG (ERROR) << "Error in copying level1 keywords" ;
            continue;
        }       
    if ((strcmp (centroid_algo , "3C") != 0) && (strcmp (centroid_algo , "3S") != 0) && (strcmp (centroid_algo , "5S") != 0))
    {
        LOG (ERROR) << "***Invalid value of Centroid_algo ***" << endl ;
        continue;
    }


    if ((strcmp (centroid_algo , "3C") == 0))
    {
        cent_corr_win_size = 5 ; //3x3 cross  (5 pixels - centre, top, bottom, left & right used)
        cent_bias_win_size = 3 ;
    }
    else if (strcmp (centroid_algo , "5S") == 0)
    {
        cent_corr_win_size = 25 ; //5x5 square
        cent_bias_win_size = 5 ;
    }
    else if ((strcmp (centroid_algo , "3S") == 0))
    {
        cent_corr_win_size = 9 ; //3x3 square
        cent_bias_win_size = 3 ;
    }
    else
    {
        LOG (ERROR) << "***INVALID value for Centroid window.allowed values are (3/5) *** " ;
         continue;
    }


    //reading caldb for centroid Bias correction
    tempname = caldb_handler.getCentroidBiasFile (datainfo.getDetector () , caldb_common_dir , cent_bias_win_size) ;

    if (tempname==" ")
    {
        LOG (ERROR) << endl << "Couldn't find CentroidBias file from caldb" << endl ;
        continue;
    }

    /***For Joining the String ***/

    joinStrings (centroidbiasfile , 2 , caldb_common_dir , tempname.c_str()) ;
    LOG (INFO) << "Centroid bias  file  from calDB :" << centroidbiasfile ;
     /*Reading a bias file From the CALDB dir*/
    status = readcentroidbiasFile () ;
    if (status)
    {
        LOG (ERROR) << "***Error Reading centroid bias File From the Caldb***";
       continue;
    }
    
    
    //reading caldb file for detector Distortion
    tempname = caldb_handler.getDetectorFile (datainfo.getDetector () , caldb_temp_dir) ;
     if (tempname==" ")
    {
        LOG(ERROR) << "Couldn't find detector File From caldb"  ;
        return (EXIT_FAILURE) ;
    }
    joinStrings (detector_distortion_corr_file , 2 , caldb_temp_dir , tempname.c_str()) ;
    LOG(INFO) << "Distortion correction file is " << detector_distortion_corr_file ;
   

    X_detect_distortion = new float[xsize * ysize] ;
    Y_detector_distortion = new float[xsize * ysize] ;
    status = caldb_handler.readCaldbDistFile (X_detect_distortion , Y_detector_distortion , detector_distortion_corr_file) ;
 
    
    //reading caldb for optical distortion correction
    tempname = caldb_handler.getOpticalDetectorFile (datainfo.getDetector () , caldb_common_dir_optic , "F0") ;
    if (tempname==" ")
    {
        LOG(ERROR) << "Couldn't find Optical Dist detector File From caldb" ;
        continue;
    }
    
    
    joinStrings (optical_distortion_corr_file , 2 , caldb_common_dir_optic , tempname.c_str()) ;
    X_optical_distortion = new float[xsize * ysize] ;
    Y_optical_distortion = new float[xsize * ysize] ;
    status = caldb_handler.readCaldbDistFile (X_optical_distortion , Y_optical_distortion , optical_distortion_corr_file) ;
    if(status){
        LOG(ERROR)<<"Error in reading the distortion correction file";
        return(EXIT_FAILURE);
    }
    
    
     tempname = caldb_handler.getTemplateFileForExposure(datainfo.getDetector () , datainfo.getObsMode () , win_xsize+1 , win_ysize+1 , (char*) caldbindir.c_str ()) ;

    if (tempname==" ")
    {
        LOG(ERROR) << "Couldn't find bad pixel file from calDB" ;
        continue;
    }
     xsize=PIX_PADD_SIZE*2;
     ysize =PIX_PADD_SIZE*2;
     joinStrings (ExpTemplateFile , 1 , tempname.c_str()) ;
    float *exp_templateArray = new float[PIX_PADD_SIZE*2*PIX_PADD_SIZE*2] ;
    status = readImage (ExpTemplateFile , 1 , exp_templateArray) ;
    if(status){
        LOG(ERROR)<<"Error in reading the  template Array of the image";
        return(EXIT_FAILURE);
    }
    LOG(INFO)<<"Reading CALDB files completed";
    xsize=datainfo.getXsize();
    ysize =datainfo.getYsize();
   /*=================================================================================================================================*/
    /*====================================================CALDB READING  FINISHED===========================================================*/
    /*=================================================================================================================================*/
    
    fits_movabs_hdu (fevt_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU" , eventfile) ;
    fits_get_num_rows (fevt_in , &nrows , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;

    unsigned short *psc = new unsigned short[nrows] ;
    unsigned short *frame_no = new unsigned short[nrows] ;
    unsigned short *Ix = new unsigned short[nrows] ;
    unsigned short *Iy = new unsigned short[nrows] ;

    double *time_frame = new double[nrows] ;
    float *fx = new float[nrows] ;
    float *fy = new float[nrows] ;
    unsigned char *dmm = new unsigned char[nrows] ;
    unsigned char *Min = new unsigned char[nrows] ;
    unsigned char *multflag = new unsigned char [nrows] ;
    unsigned char *badflag = new unsigned char [nrows] ;
     double *time_frame_uvt = new double[nrows] ;
    fits_read_col (fevt_in , TUSHORT , 1 , 1 , 1 , nrows , NULL , psc , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
    fits_read_col (fevt_in , TUSHORT , 2 , 1 , 1 , nrows , NULL , frame_no , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
    fits_read_col (fevt_in , TDOUBLE , 3 , 1 , 1 , nrows , NULL , time_frame , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
    fits_read_col (fevt_in , TUSHORT , 4 , 1 , 1 , nrows , NULL , Ix , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
    fits_read_col (fevt_in , TFLOAT , 5 , 1 , 1 , nrows , NULL , fx , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
    fits_read_col (fevt_in , TUSHORT , 6 , 1 , 1 , nrows , NULL , Iy , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
    fits_read_col (fevt_in , TFLOAT , 7 , 1 , 1 , nrows , NULL , fy , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
    fits_read_col (fevt_in , TBYTE , 8 , 1 , 1 , nrows , NULL , dmm , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
    fits_read_col (fevt_in , TBYTE , 9 , 1 , 1 , nrows , NULL , Min , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
   fits_read_col (fevt_in , TDOUBLE , 10 , 1 , 1 , nrows , NULL , time_frame_uvt , NULL , &status) ;
    printError (status , "Error in reading the number of rows" , file_in) ;
    fits_close_file (fevt_in , &status) ;


//calculating Total frames came as a input.
cntTotalFrameInput=1;
if(nrows>0){
unsigned short frmcntref= frame_no[0];
for(int i =1;i<nrows;i++)
{
   if(frame_no[i]!=frmcntref){
frmcntref=frame_no[i];
cntTotalFrameInput++;
}
}
}

//till this
    

    //masking bad pixel calculation
     LOG(INFO)<<"Masking Bad Pixel";
    int TotalNum_Effected_Rows = performMaskBadPix (Ix , Iy , badpixdata , dmm , badflag , multflag , nrows , xsize , ysize , thr_multiph) ;
    TotalNum_Effected_Rows=0;
    int current_Frmno=0,nextFrameNo=0,lastFrameNo=0;
    int k=0;
    for (int i =0;i<nrows;i++){
       
        nextFrameNo=i;
        lastFrameNo=i;
        if (multflag[i]==0){
            //TotalNum_Effected_Rows++;
            current_Frmno=frame_no[i];
            if(nextFrameNo+1 <=nrows-1) {
           
            while(current_Frmno==frame_no[++nextFrameNo] && multflag[nextFrameNo]!=0){
                multflag[nextFrameNo]=0;
              //  TotalNum_Effected_Rows++;
                i++;
            }
            //i--;
            }
             if(lastFrameNo-1>=1) {
            while(current_Frmno==frame_no[--lastFrameNo] && multflag[lastFrameNo]!=0){
                multflag[lastFrameNo]=0;
                //TotalNum_Effected_Rows++;
            }
           }
            
        }
        
         //k=i;
    }
    for (int i =0;i<nrows;i++){
        if (multflag[i]==0){
            TotalNum_Effected_Rows++;
        }
    }
    
    
    
    if (wtd_bp == 1)//incase of mask bad pixel  to be written on disk
    {
        status = setDirectoryStructure (moduleoutdir_bp , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return (EXIT_FAILURE) ;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_bp , nameprefix , "bp") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 12 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TUSHORT , 4 , Ix , nrows , TFLOAT , 5 , fx , nrows , TUSHORT , 6 , Iy , nrows ,
                TFLOAT , 7 , fy , nrows , TBYTE , 8 , dmm , nrows , TBYTE , 9 , Min , nrows ,TDOUBLE,10,time_frame_uvt,nrows, TBYTE , 11 , badflag , nrows , TBYTE , 12 , multflag , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  bad pixel FITS  file" << endl ;
             continue;
        }

    }

    LOG(INFO)<<"Performing Pixel padding ";
    status = performPixPadding (Ix , Iy , nrows) ;
    xsize = 600 ;
    ysize = 600 ;
   
    if (wtd_pp == 1)//incase of pix padding to be written on disk
    {

        status = setDirectoryStructure (moduleoutdir_pp , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             continue;
        }

        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_pp , nameprefix , "pp") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 12 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TUSHORT , 4 , Ix , nrows , TFLOAT , 5 , fx , nrows , TUSHORT , 6 , Iy , nrows ,
                TFLOAT , 7 , fy , nrows , TBYTE , 8 , dmm , nrows , TBYTE , 9 , Min , nrows ,TDOUBLE,10,time_frame_uvt,nrows, TBYTE , 11 , badflag , nrows , TBYTE , 12 , multflag , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  Pix padding FITS  file" << endl ;
            continue;
        }

    }
   // LOG(INFO)<<"Performing Sub division";

    float subDiv_fact = subDivision_size / xsize ;

    float temp_x = 0.0f ;
    float temp_y = 0.0f ;

    //performing Subdivision
    for (int j = 0 ; j < nrows ; j ++)
    {
        temp_x = subDiv_fact * (Ix[j] + fx[j]) ;//multiply with  subDiv_fact  to convert 600*600 scheme to 9600*9600 scheme
        temp_y = subDiv_fact * (Iy[j] + fy[j]) ;
        fx[j] = temp_x - (int) (temp_x) ;
        fy[j] = temp_y - (int) temp_y ;
        Ix[j] = (int) temp_x ;
        Iy[j] = (int) temp_y ;
    }

    if (wtd_sd == 1)//incase of subdivision to written on disk 
    {
        status = setDirectoryStructure (moduleoutdir_sd , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return (EXIT_FAILURE) ;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_sd , nameprefix , "sd") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;
        status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TUSHORT , 4 , Ix , nrows , TFLOAT , 5 , fx , nrows , TUSHORT , 6 , Iy , nrows ,
                TFLOAT , 7 , fy , nrows , TBYTE , 8 , dmm , nrows , TBYTE , 9 , Min , nrows ,TDOUBLE,10,time_frame_uvt,nrows, TBYTE , 11 , badflag , nrows , TBYTE , 12 , multflag , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  subDivision  FITS  file" << endl ;
             continue;
        }
    }
    xsize = subDivision_size ;
    ysize = subDivision_size ;
    //subdivision finished.       

    //performing cosmic Ray correction
    long nFrames ;
    LOG(INFO)<<"Performing Cosmic ray correction";
    nFrames = frame_no[nrows - 1] ; //number of frames

    long nFramesCounter = 0 , j = 1 ;
    vector<long> del_CReffctedrows , frameno_crFailed ;
    vector<double> x_crFailed , y_crFailed ;
    long nEventsInFrame[nFrames] , GoodFrames[nFrames] , CRAFFECTED[nFrames] , m = 0 ;

    for (int i=0;i<nFrames;i++)  nEventsInFrame[i]=0;
    
     if(nrows==1){
            nEventsInFrame[nFrames-1]=1;
            GoodFrames[nFrames-1]=frame_no[0];
            j=1;
        }
    
    for (long i = 0 ; i < nrows - 1 ; i ++)
    {
        if ((frame_no[i] == frame_no[i + 1]))//looking for continues frame number 
        {

            nFramesCounter ++ ; //counter for number of events in a frame
            if (i == nrows - 2)
            {
                nEventsInFrame[j-1] = nFramesCounter ;
                GoodFrames[j-1] = frame_no[i] ;
            }

        }
        else //if frame number will change in next iteration
        {
            nEventsInFrame[j-1] = nFramesCounter ;//stores total number of events occure for current frame
            GoodFrames[j-1] = frame_no[i] ;
            j ++ ;
            nFramesCounter = 0;
            if ((frame_no[i + 1] - frame_no[i]) > 1)//incase of missing frame(i.e not continues frames one by one)
            {
                m = (frame_no[i + 1] - frame_no[i]) ;
                for (int index = i ; index < m ; index ++)
                {
                    nEventsInFrame[j-1] = 0 ;
                    //this is added.
                    GoodFrames[j-1] = 0 ;
                    j ++ ;
                }
            }
        }
        
    }
    
    int diff=nFrames-j;
    //nFrames=j;
//    LOG(INFO)<<j;exit(1);
    int cra = 0 ;//cr affected 
    int sum = 0 ;
  //  int cntsum=0;
   // double meanVal,sdval;
   // vector<int> TotEvents;
    vector<bool> flag_Forthrcheck;
    //double previous_frmSD=0;
    double curr_frmthr=0.0;
    sum = 0 ;
    //TotEvents.clear();
    int cnt_ForTotFrames=0;
    for (long i = 0; i < j ; i ++)
    {
              
        if (nEventsInFrame[i] != 0)
        {  
              
                sum = sum + nEventsInFrame[i] ;//adding number of events in total cr_avg_frames_compare frames 
                cnt_ForTotFrames++;
           
        }
    }
    if(cnt_ForTotFrames==0){
            LOG(INFO)<<"***Devide by ZERO,While calculating average value for events***";
            continue;
        }
    double avg_Val_Events=sum/cnt_ForTotFrames;
    LOG(INFO)<<"\033[1;31mCalculated Average events per frame-> "<<avg_Val_Events<<"\033[0m";
    double square_root_of_avg=sqrt(avg_Val_Events);
     curr_frmthr=avg_Val_Events+square_root_of_avg*first_mult_factor_CR+second_mult_Factor_CR/square_root_of_avg; 
     LOG(INFO)<<"\033[1;31mCalculated threshold value for CosmicRay removal -> "<<curr_frmthr<<"\033[0m";
    
      for (long i = 0; i < j ; i ++){
            if ((nEventsInFrame[i]) > curr_frmthr  )
            {
                
                flag_Forthrcheck.push_back(0);
                CRAFFECTED[cra] = GoodFrames[i] ;
                cra ++ ;
            }
            else{
                flag_Forthrcheck.push_back(1);
            }
    
      }
    int flag = 0 ;
//creating list of frames that are CR effected.
    for (long i = 0 ; i < nrows ; i ++)
    {
        flag = 0 ;
        for (long m = 0 ; m < cra ; m ++)
        {
            if (frame_no[i] == CRAFFECTED[m])
            {
                flag = 1 ;
                break ;
            }
        }

        if (flag == 1)
        {
            del_CReffctedrows.push_back (i + 1) ;
            frameno_crFailed.push_back (frame_no[i]) ;
            x_crFailed.push_back (Ix[i] + fx[i]) ;
            y_crFailed.push_back (Iy[i] + fy[i]) ;
        }
    }
    sort(frameno_crFailed.begin(),frameno_crFailed.end());
    frameno_crFailed.erase(unique(frameno_crFailed.begin(),frameno_crFailed.end()),frameno_crFailed.end());
    LOG(INFO)<<"CRAffected row "<<del_CReffctedrows.size();
    // This loop is for common array only.
    vector<unsigned short> frmno , pack_seq , intX , intY ;
    vector<double> tm_frm ,tm_uvit;
    vector<float> fract_X , fract_Y ;
    vector<unsigned char> dmm_temp , min_temp , badflag_temp , multflag_temp ;
   // bool flag_cr = 0 ;
    //int last_index=0;
//removing the records which are CR effected
    
    //for (int j = 0 ; j <= nrows - 1 ; j ++)
   // {
       // LOG(INFO)<<j;
    
    if(del_CReffctedrows.size()==0){
        
        for (int i=0;i<nrows;i++)
        {
            frmno.push_back (frame_no[i]) ;
            pack_seq.push_back (psc[i]) ;
            intX.push_back (Ix[i]) ;
            intY.push_back (Iy[i]) ;
            tm_frm.push_back (time_frame[i]) ;
            fract_X.push_back (fx[i]) ;
            fract_Y.push_back (fy[i]) ;
            dmm_temp.push_back (dmm[i]) ;
            min_temp.push_back (Min[i]) ;
            badflag_temp.push_back (badflag[i]) ;
            multflag_temp.push_back (multflag[i]) ;
tm_uvit.push_back(time_frame_uvt[i]);
        }
    }
    else {
       // flag_cr = 0 ;
        for(int i=0;i<del_CReffctedrows[0];i++) 
        {
            frmno.push_back (frame_no[i]) ;
            pack_seq.push_back (psc[i]) ;
            intX.push_back (Ix[i]) ;
            intY.push_back (Iy[i]) ;
            tm_frm.push_back (time_frame[i]) ;
            fract_X.push_back (fx[i]) ;
            fract_Y.push_back (fy[i]) ;
            dmm_temp.push_back (dmm[i]) ;
            min_temp.push_back (Min[i]) ;
            badflag_temp.push_back (badflag[i]) ;
            multflag_temp.push_back (multflag[i]) ;
tm_uvit.push_back(time_frame_uvt[i]);
        }
        
        for (int i = 0; i < del_CReffctedrows.size ()-1 ; i ++)
        {
            for(int j=del_CReffctedrows[i]+1;j<del_CReffctedrows[i+1];j++)
            {
               // if(del_CReffctedrows[i]+1<del_CReffctedrows[i+1])
               // {
                frmno.push_back (frame_no[j-1]) ;
                pack_seq.push_back (psc[j-1]) ;
                intX.push_back (Ix[j-1]) ;
                intY.push_back (Iy[j-1]) ;
                tm_frm.push_back (time_frame[j-1]) ;
                fract_X.push_back (fx[j-1]) ;
                fract_Y.push_back (fy[j-1]) ;
                dmm_temp.push_back (dmm[j-1]) ;
                min_temp.push_back (Min[j-1]) ;
                badflag_temp.push_back (badflag[j-1]) ;
                multflag_temp.push_back (multflag[j-1]) ;
tm_uvit.push_back(time_frame_uvt[i]);
               // }
               
//            if (del_CReffctedrows[i] == j + 1)
//            {   
//                //LOG(INFO)<<del_CReffctedrows[i]<<" "<<j;
//                flag_cr = 1 ;
//                last_index=i;
//                break ;
//            }
             }
//        if (flag_cr == 0)
//        {
//           
//        }

    }
        for (int i=del_CReffctedrows[del_CReffctedrows.size()-1]+1;i<nrows;i++)
        {
            frmno.push_back (frame_no[i-1]) ;
            pack_seq.push_back (psc[i-1]) ;
            intX.push_back (Ix[i-1]) ;
            intY.push_back (Iy[i-1]) ;
            tm_frm.push_back (time_frame[i-1]) ;
            fract_X.push_back (fx[i-1]) ;
            fract_Y.push_back (fy[i-1]) ;
            dmm_temp.push_back (dmm[i-1]) ;
            min_temp.push_back (Min[i-1]) ;
            badflag_temp.push_back (badflag[i-1]) ;
            multflag_temp.push_back (multflag[i-1]) ;
tm_uvit.push_back(time_frame_uvt[i]);
        }
    
    }
//    for (int j = 0 ; j <= nrows - 1 ; j ++)
//    {
//       // LOG(INFO)<<j;
//        flag_cr = 0 ;
//        for (int i = last_index ; i < del_CReffctedrows.size () ; i ++)
//        {
//            if (del_CReffctedrows[i] == j + 1)
//            {   
//                //LOG(INFO)<<del_CReffctedrows[i]<<" "<<j;
//                flag_cr = 1 ;
//                last_index=i;
//                break ;
//            }
//        }
//        if (flag_cr == 0)
//        {
//            frmno.push_back (frame_no[j]) ;
//            pack_seq.push_back (psc[j]) ;
//            intX.push_back (Ix[j]) ;
//            intY.push_back (Iy[j]) ;
//            tm_frm.push_back (time_frame[j]) ;
//            fract_X.push_back (fx[j]) ;
//            fract_Y.push_back (fy[j]) ;
//            dmm_temp.push_back (dmm[j]) ;
//            min_temp.push_back (Min[j]) ;
//            badflag_temp.push_back (badflag[j]) ;
//            multflag_temp.push_back (multflag[j]) ;
//        }
//
//    }
 //LOG(INFO)<<"HHH "<<cra<<" "<<del_CReffctedrows.size()<<" "<<del_CReffctedrows[del_CReffctedrows.size()-1]<<" "<<frmno.size()<<" "<<nrows<<" "<<frame_no[0]<<" "<<del_CReffctedrows[0]<<" "<<del_CReffctedrows[1]<<" "<<del_CReffctedrows[2];
 
    delete[] frame_no , psc , Ix , Iy , time_frame , fx , fy , dmm , Min , badflag , multflag,time_frame_uvt ;
    psc = new unsigned short[frmno.size ()] ;
    frame_no = new unsigned short[frmno.size ()] ;
    Ix = new unsigned short[frmno.size ()] ;
    Iy = new unsigned short[frmno.size ()] ;
    time_frame = new double[frmno.size ()] ;
    fx = new float[frmno.size ()] ;
    fy = new float[frmno.size ()] ;
    dmm = new unsigned char[frmno.size ()] ;
    Min = new unsigned char[frmno.size ()] ;
    multflag = new unsigned char [frmno.size ()] ;
    badflag = new unsigned char [frmno.size ()] ;
	time_frame_uvt= new double[frmno.size ()] ;;
    //copy to  different Array for the use of  further module
    for (int i = 0 ; i < frmno.size () ; i ++)
    {
        frame_no[i] = frmno[i] ;
        psc[i] = pack_seq[i] ;
        Ix[i] = intX[i] ;
        Iy[i] = intY[i] ;
        fx[i] = fract_X[i] ;
        fy[i] = fract_Y[i] ;
        dmm[i] = dmm_temp[i] ;
        badflag[i] = badflag_temp[i] ;
        multflag[i] = multflag_temp[i] ;
        time_frame[i] = tm_frm[i] ;
        Min[i] = min_temp[i] ;
	time_frame_uvt[i]=tm_uvit[i];
    }

    if (wtd_cr == 1)//incase of cosmic ray output to be written to the disk.
    {

        status = setDirectoryStructure (moduleoutdir_cr , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             continue;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_cr , nameprefix , "cr") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 12 , TUSHORT , 1 , psc , frmno.size () , TUSHORT , 2 , frame_no , frmno.size () , TDOUBLE , 3 , time_frame , frmno.size () , TUSHORT , 4 , Ix , frmno.size () , TFLOAT , 5 , fx , frmno.size () , TUSHORT , 6 , Iy , frmno.size () ,
                TFLOAT , 7 , fy , frmno.size () , TBYTE , 8 , dmm , frmno.size () , TBYTE , 9 , Min , frmno.size () ,TDOUBLE,10,time_frame_uvt,frmno.size (), TBYTE , 11 , badflag , frmno.size () , TBYTE , 12 , multflag , frmno.size ()) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  Pix padding FITS  file" << endl ;
             continue;
        }

    }
//    
    //performing centroid Correction
    float *newXfract = new float[frmno.size ()] ;
    float *newYfract = new float [frmno.size ()] ;
  
 //for (int i=0;i<nrows;i++)
  // {
//	newXfract[i]=Ix[i]+fx[i];
  //      newYfract[i]=Iy[i]+fy[i];
        //new_tempX[i]=Ix[i]+fx[i];
        //new_tempY[i]=Iy[i]+fy[i];
        
    //}
if(centCorrflag==1)
{
    LOG(INFO)<<"Performing Centroid Correction";
    status = performCentroidCorr (time_frame , Ix , Iy , fx , fy , Min , newXfract , newYfract , (float*) darkFramestart_data , (float*) darkFrameend_data , datainfo.getIntegrationTime () , EA , (int) frmno.size ()) ;
    if (status)
    {
        LOG(ERROR) << "Error in Centroid Correction module" << endl ;
        return (EXIT_FAILURE) ;
    }
    char *ttype_1[] = {"PacketSequence " , "FrameCount" , "Time" , "Fx" , "Fy" , "Max-Min" , "Min" ,"UVIT_MASTER_TIME ","BAD FLAG" , "MULT PHOTON"} ;
    char *tform_1[] = {"U" , "U" , "D" , "1D" , "1D" , "B" , "B" ,"D" ,"B" , "B"} ;
char * tunit_1[] = {"" , "" , "" , "" , "" , "" , "" , "" , "",""} ;
    if (wtd_centCorr == 1)//incase of centroid correction to be written on disk
    {
        tfields = 10 ;

        status = setDirectoryStructure (moduleoutdir_centroidCorr , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return (EXIT_FAILURE) ;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_centroidCorr , nameprefix , "cc") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_1 , tform_1 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        // status = writeColumnsToFITS (outfile , 2 , 9 , TUSHORT , 1 , psc , frmno.size () , TUSHORT , 2 , frame_no , frmno.size () , TDOUBLE , 3 , time_frame , frmno.size () , TUSHORT , 4 , Ix , frmno.size () , TFLOAT , 5 , newXfract , frmno.size () , TUSHORT , 6 , Iy , frmno.size () ,
       //         TFLOAT , 7 , newYfract , frmno.size () , TBYTE , 8 , dmm , frmno.size () , TBYTE , 9 , Min , frmno.size () , TBYTE , 10 , badflag , frmno.size () , TBYTE , 11 , multflag , frmno.size ()) ;
 status = writeColumnsToFITS (outfile , 2 , 10 , TUSHORT , 1 , psc , frmno.size () , TUSHORT , 2 , frame_no , frmno.size () , TDOUBLE , 3 , time_frame , frmno.size () , TFLOAT , 4 , newXfract , frmno.size () , TFLOAT , 5 , newYfract ,
                frmno.size () , TBYTE , 6 , dmm , frmno.size () , TBYTE , 7 , Min , frmno.size () , TDOUBLE,8,time_frame_uvt,frmno.size (),TBYTE , 9 , badflag , frmno.size () , TBYTE , 10 , multflag , frmno.size ()) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  Pix padding FITS  file" << endl ;
             continue;
        }

    }
}
else
{
for(int i=0;i<frmno.size();i++)
{
 newXfract[i]=Ix[i]+fx[i];
 newYfract[i]=Iy[i]+fy[i];
}

}
    nrows=frmno.size();
float *effective_NumPhotons = new float[nrows] ;
    for (int i = 0 ; i < nrows ; i ++) effective_NumPhotons[i] = 1.0f ;

    tfields = 11 ;
    char *ttype_2[] = {"PacketSequence " , "FrameCount" , "Time" , "Fx" , "Fy" , "Max-Min" , "Min" ,"UVIT_MASTER_TIME", "BAD FLAG" , "MULT PHOTON" , "EFFECTIVE_NUM_PHOTONS"} ;
    char *tform_2[] = {"U" , "U" , "D" , "1D" , "1D" , "B" , "B" ,"D" ,"B" , "B" , "1D"} ;
    char *tunit_2[] = {"" , "" , "" , "" , "" , "" , "" , "" , "" , ""," "} ;
  //  nrows = frmno.size () ;
    
    float *new_tempX = new float [nrows] ;
    float *new_tempY = new float [nrows] ;
   if(centBiasflag==1)
{
    LOG(INFO)<<"Performing Centroid Bias";
    //performing centroid bias.
    status = performCentroidBias (nrows , newXfract , newYfract , new_tempX , new_tempY) ;
        
    if (wtd_centBias == 1)//incase of centroid bias output to be written on the disk
    {

        status = setDirectoryStructure (moduleoutdir_centroidBias , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             continue;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_centroidBias , nameprefix , "cb") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows ,TDOUBLE,8,time_frame_uvt,nrows, TBYTE , 9 , badflag , nrows , TBYTE , 10 , multflag , nrows , TFLOAT , 11 , effective_NumPhotons , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the Centroid BIAS  FITS  file" << endl ;
             continue;
        }

     }
}
else
{
	for(int i=0;i<nrows;i++)
	{
		new_tempX[i]=newXfract[i];
		new_tempY[i]=newYfract[i];

	}	
}
//performing Flat Field Correction

     //// if(flatfieldFlag){
//    LOG(INFO)<<"Performing Flat field  Correction";
   status = performFlatFieldCorr (effective_NumPhotons , finalFlatFielddata , new_tempX , new_tempY , nrows) ;
    if (wtd_ff == 1)//incase of flat field to be written on disk
    {

        status = setDirectoryStructure (moduleoutdir_ff , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            continue;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_ff , nameprefix , "ff") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows ,TDOUBLE,8,time_frame_uvt,nrows, TBYTE , 9 , badflag , nrows , TBYTE , 10 , multflag , nrows , TFLOAT , 11 , effective_NumPhotons , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  Pix padding FITS  file" << endl ;
             continue;
        }

    }
    ////  }
//    LOG(INFO)<<"Performing QEMCP correction";
     if(UTC_flag==0 && qemcpflag==1){
                LOG(WARNING)<<"***QEMCP correction cant be kept ON if UTC correction not to be done!!,NOW skipping QEMCP correction***";
                qemcpflag=0;
            }
       bool flag_qemcp_success=FALSE;
       double t1 , t2 , x1 , x2 ,factor;
       if(qemcpflag==1){
       
        for (int p = 0 ; p < nrows ; p++)//loop for total number of events 
        {
         temperature=INITIALIZATION_VALUE;        
         for (int i=0;i<nrows_lbt;i++)//loop for number of rows of LBT file for finding the relevant temperature of  current frame based on the frame time   
        {
                    
            if (time_frame[p]>=time_lbt[i] && time_frame[p]<time_lbt[i+1])
            {
                
                temperature=(insideTemp[i]+outsideTemp[i])/2;
                break;
            }
        }
        if(temperature==INITIALIZATION_VALUE)//incase of no relevant temperature found for the current frame 
        {
            LOG(ERROR)<<"No record found in LBT file"<<p<<" "<<time_frame[p];
            break;
        }
	flag_qemcp_success=TRUE;
                for (int j = 0 ; j < nCalDBTempValues - 1 ; j++)//loop for finding the relevant factor value by comparing calculated temperature value and  caldb temperature value 
            {
                if (temperature >= cal_temp[j] && temperature < cal_temp[j + 1])          //temperature - from LBT file ,cal_temp- from caldb file
                {
                    t1 = cal_temp[j] ;
                    t2 = cal_temp[j + 1] ;
                    if ((t2 - t1) == 0)
                    {
                        LOG(INFO) << "***Divide By zero***" << endl ;
                        return (EXIT_FAILURE) ;
                    }
                    x1 = qe_mg_factor[j] ;
                    x2 = qe_mg_factor[j + 1] ;
                    factor = x1 + ((temperature - t1)*((x2 - x1) / (t2 - t1))) ; //calculated factor
                }
            }
//    
         effective_NumPhotons[p]=effective_NumPhotons[p]*factor; //performing QE and MCP correction(by applying factor to effective number of photons )
        }
     if(flag_qemcp_success==FALSE)
            {
	LOG(ERROR)<<"Error in finding the related  factor for QEMCP correction";
         return(EXIT_FAILURE);
            }
    if (wtd_qemcp)//incase of QE and MCp to be written on the disk.
    {
        status = setDirectoryStructure (moduleoutdir_qemcp , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            continue;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_qemcp , nameprefix , "qe") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows , TDOUBLE,8,time_frame_uvt,nrows,TBYTE , 9 , badflag , nrows , TBYTE , 10 , multflag , nrows , TFLOAT , 11 , effective_NumPhotons , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  Pix padding FITS  file" << endl ;
            return(EXIT_SUCCESS);
        }


    }
       }
    LOG(INFO)<<"Performing Unit Conversion";
    status = performUnitConversion (effective_NumPhotons , datainfo.getIntegrationTime () , nrows) ;
    if(status){
        LOG(ERROR)<<"***Integration time is less than or equal to 0***";
	LOG(ERROR)<<"CRASH FRAME INTEGRATION TIME =< 0 (uvtLevel2PC.cpp)";
        continue;
    }
    if (wtd_uc == 1)//incase of unit conversion to be written to tthe disk.
    {

        status = setDirectoryStructure (moduleoutdir_uc , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            continue;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_uc , nameprefix , "uc") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows ,TDOUBLE,8,time_frame_uvt,nrows, TBYTE , 9 , badflag , nrows , TBYTE , 10 , multflag , nrows , TFLOAT , 11 , effective_NumPhotons , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  Pix padding FITS  file" << endl ;
            continue;
        }

    }
//    float *new_badPixarr= new float[PIX_PADD_SIZE*PIX_PADD_SIZE];
//    for (int i=0;i<PIX_PADD_SIZE*PIX_PADD_SIZE;i++)  new_badPixarr[i]=0.0;      
//    status=Applypadding(badpixdata,512,512,new_badPixarr,PIX_PADD_SIZE,PIX_PADD_SIZE);
//    if(status){
//        LOG(ERROR)<<"Error in pix padding on exposure array" ;
//        return(EXIT_FAILURE);
//    }
//    vector<float> badpix_xloc,badpix_yloc,badpix_orig_x,badpix_orig_y;
//    int xlocbadpix=0,ylocbadpix=0;
//    for (int i=0;i<PIX_PADD_SIZE*PIX_PADD_SIZE;i++)
//    {
//        xlocbadpix=(int)i%600;
//        ylocbadpix=(int)i/600;
//        
//        //sqrt((ylocbadpix-300)*(ylocbadpix-300)+(xlocbadpix-300)*(xlocbadpix-300))<255
//        if(new_badPixarr[i]==0.0  &&  sqrt((ylocbadpix-300)*(ylocbadpix-300)+(xlocbadpix-300)*(xlocbadpix-300))<245 )
//        {
//            //LOG(INFO)<<sqrt((ylocbadpix-300)*(ylocbadpix-300)+(xlocbadpix-300)*(xlocbadpix-300));
//            badpix_xloc.push_back(((i%PIX_PADD_SIZE)+1)*1.0);
//            badpix_yloc.push_back(((i/PIX_PADD_SIZE)+1)*1.0);
//            badpix_orig_x.push_back(((i%PIX_PADD_SIZE)+1)*1.0);
//            badpix_orig_y.push_back(((i/PIX_PADD_SIZE)+1)*1.0);
//           // LOG(INFO)<<(i%PIX_PADD_SIZE)+1<<" "<<(i/PIX_PADD_SIZE)+1;
//        }
//    }
//    
if(DetectDistflag==1)
{     
//    status = performDistCorrection (badpix_xloc.size() , badpix_xloc.data() , badpix_yloc.data() , X_detect_distortion , Y_detector_distortion ,IMG_DIM_DI) ;
//    if(status){
//        LOG(ERROR)<<"Error in performing detector distortion correction in Exposure frames";
//        return(EXIT_FAILURE);
//    }
   
    status = performDistCorrection (nrows , new_tempX , new_tempY , X_detect_distortion , Y_detector_distortion ,IMG_DIM_DI) ;
    if (wtd_dd == 1)//incase of detector distortion to be written to the disk.
    {

        status = setDirectoryStructure (moduleoutdir_dd , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            continue;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_dd , nameprefix , "dd") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows ,TDOUBLE,8,time_frame_uvt,nrows, TBYTE , 9 , badflag , nrows , TBYTE , 10 , multflag , nrows , TFLOAT , 11 , effective_NumPhotons , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  Detect Dist  Correction File" << endl ;
             continue;
        }
   //      fits_open_file (&fout , outfile , READWRITE , &status) ;
     //   printError (status , " Error in opening file " , outfile) ;
       // fits_movabs_hdu (fout , 2 , NULL , &status) ;
       // printError (status , "Error in moving to HDU " , outfile) ;
       // fits_delete_rowlist (fout , ListOf_OutsideEvnts.data () , ListOf_OutsideEvnts.size () , &status) ;
       // printError (status , "Error in Deleting the row list " , outfile) ;
        //fits_close_file (fout , &status) ;
        //printError (status , "Error closing the file " , outfile) ;
        

    }
   // LOG(INFO)<<ListOf_OutsideEvnts.size();
  //  ListOf_OutsideEvnts.clear();
//    for (int i=0;i<badpix_orig_x.size();i++)
//    {
//        
//         new_badPixarr[(int)((round(badpix_yloc[i]))*PIX_PADD_SIZE+round(badpix_xloc[i]))]=0.0;  
//        
//        
//        
//    }
//    for (int i=0;i<badpix_orig_x.size();i++)
//    {
//        
//         new_badPixarr[(int)((round(badpix_orig_y[i]))*PIX_PADD_SIZE+round(badpix_orig_x[i]))]=1.0;  
//        
//                
//    }
}
//    badpix_orig_x.clear();
//    badpix_orig_y.clear();
//    for (int i=0;i<badpix_xloc.size();i++)
//    {
//            badpix_orig_x.push_back(badpix_xloc[i]);
//            badpix_orig_y.push_back(badpix_yloc[i]);
//    }
//    status = setDirectoryStructure ("abcd1" , " ") ;
//        if (status)
//        {
//            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
//            continue;
//        }
//    writeOutputImageToDisk ("fi" , "abcd1" , "" , "" , new_badPixarr , "abc" ,0.0 ,1 , 600 , 600) ; //this is for the SignalFram
 LOG(INFO)<<"Performing Optical Assembly distortion  Correction";
if(OpticDistflag==1){
//    status = performDistCorrection (badpix_xloc.size() , badpix_xloc.data() , badpix_yloc.data() , X_optical_distortion , Y_optical_distortion ,IMG_DIM_DI) ;
//    if(status){
//        LOG(ERROR)<<"Error in performing detector distortion correction in Exposure frames";
//        return(EXIT_FAILURE);
//    }
    //LOG(INFO)<<badpix_orig_x[0]<<" "<<badpix_xloc[0];exit(1);
    
    status = performDistCorrection (nrows , new_tempX , new_tempY , X_optical_distortion , Y_optical_distortion , IMG_DIM_DI) ;
    if (wtd_od == 1)//incase of optical distortion to be written to the disk.
    {
        status = setDirectoryStructure (moduleoutdir_od , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             continue;
        }
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_od , nameprefix , "od") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows ,TDOUBLE,8,time_frame_uvt,nrows, TBYTE , 9 , badflag , nrows , TBYTE , 10 , multflag , nrows , TFLOAT , 11 , effective_NumPhotons , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  optical Assmbly Dist  Correction File" << endl ;
             continue;
        }
  //      fits_open_file (&fout , outfile , READWRITE , &status) ;
    //    printError (status , " Error in opening file " , outfile) ;
      //  fits_movabs_hdu (fout , 2 , NULL , &status) ;
       // printError (status , "Error in moving to HDU " , outfile) ;
        //fits_delete_rowlist (fout , ListOf_OutsideEvnts.data () , ListOf_OutsideEvnts.size () , &status) ;
       // printError (status , "Error in Deleting the row list " , outfile) ;
        //fits_close_file (fout , &status) ;
        //printError (status , "Error closing the file " , outfile) ;

    }
     //LOG(INFO)<<ListOf_OutsideEvnts.size();
//    for (int i=0;i<badpix_orig_x.size();i++)
//    {
//       
//         new_badPixarr[(int)((round(badpix_yloc[i]))*PIX_PADD_SIZE+round(badpix_xloc[i]))]=0.0;  
//        
//        
//        
//    }
//    for (int i=0;i<badpix_orig_x.size();i++)
//    {
//       
//         new_badPixarr[(int)(round(badpix_orig_y[i])*PIX_PADD_SIZE+round(badpix_orig_x[i]))]=1.0;   
//        
//        
//        
//    }
    
    
  }
 //badpix_orig_x.clear();
 //badpix_orig_y.clear();
  
    unsigned short *mult_temp = new unsigned short[nrows] ;
    unsigned short *badFlag_temp = new unsigned short[nrows] ;

    for (int i = 0 ; i < nrows ; i ++)
    {
        mult_temp[i] = multflag[i] ;
        badFlag_temp[i] = badflag[i] ;
    }

    float *xi_final = new float[nrows] ;
    float *yi_final = new float[nrows] ;
    //Added-keeping track of exposure frame
    vector<float> xShift_track,yShift_track,thetaShift_track;
    vector<float> xShift_track_final,yShift_track_final,thetaShift_track_final; 
    vector<unsigned short> frm_no_vect,frm_no_vect_final;
    vector<unsigned short> newmultflag,newmultflag_final;
    //till this
    LOG(INFO)<<"Performing Shift and Rotation";
    track_rown.clear ();
    
    status = performShiftNRot (dirobj.sciencedatafile[dataindex],nrows , frame_no,time_frame , new_tempX , new_tempY , xi_final , yi_final,xShift_track,yShift_track,thetaShift_track,frm_no_vect,mult_temp,newmultflag) ;
    if (status)
    {
        LOG (ERROR) << "Error in performing shift and rotate" ;
         continue;
    }
    
    if(xShift_track.size()==0){
        LOG(ERROR)<<"***No shifts found accordingly to RAS file***";
        continue;
        //return(EXIT_FAILURE);
    }
//    for(int i=0;i<nrows;i++)
//    {
//      	
//
//    } 
    
   //LOG(INFO)<<track_rown.size();
   if(fi_flag==FALSE)
        {  
       if(track_rown.size()>0){
       if (track_rown[track_rown.size()-1]<nrows){
          
           nFrameIntegrate_fi=frame_no[nrows-1];
       }
             
       else
       {
           
                            
           long i=nrows;
           long j=track_rown.size()-1;
          // LOG(INFO)<<track_rown[j]<<" "<<i;
//           for (int j=track_rown.size()-1;j>0;j--)
//           {
               while ((track_rown[j]==i))
              { 
                i--;
                j--;
                nFrameIntegrate_fi=frame_no[i-1];
              }   
          // }
        
       }
            
       }
       else
       {
            
       nFrameIntegrate_fi=frame_no[nrows-1];
       }
   }
   
   // LOG(INFO)<<nFrameIntegrate_fi;
   // LOG(INFO)<<xShift_track.size();
    float temp_val_x=xShift_track[0];
    float temp_val_y=yShift_track[0];
    float temp_val_theta=thetaShift_track[0];
    unsigned short temp_frm_no_val=frm_no_vect[0];
    unsigned short temp_val_multFlag = newmultflag[0];
    //logic to remove the duplicates.
    //LOG(INFO)<<frm_no_vect[0];
    for (int i=0;i<xShift_track.size();i++)
    {
       // LOG(INFO)<<temp_val_multFlag;
        if(newmultflag[i]!=0){
        xShift_track_final.push_back(temp_val_x);
        yShift_track_final.push_back(temp_val_y);
        thetaShift_track_final.push_back(temp_val_theta);
        frm_no_vect_final.push_back(temp_frm_no_val);
        newmultflag_final.push_back(temp_val_multFlag);
        }
        else {
            
            continue;
        }
        while((temp_val_x==xShift_track[i] && temp_val_y==yShift_track[i] && temp_val_theta==thetaShift_track[i] && temp_frm_no_val==frm_no_vect[i]  && i <xShift_track.size()) )
        {           
            i++;
        }
        
        temp_val_x=xShift_track[i];
        temp_val_y=yShift_track[i];
        temp_val_theta=thetaShift_track[i];
        temp_frm_no_val=frm_no_vect[i];
        temp_val_multFlag = newmultflag[i];
        if(i==xShift_track.size()-1 && frm_no_vect_final[frm_no_vect_final.size()-1]!=frm_no_vect[i] && newmultflag[i] !=0){//incase of last row frame count is diferent from previous ones.
        
            
        xShift_track_final.push_back(temp_val_x);
        yShift_track_final.push_back(temp_val_y);
        thetaShift_track_final.push_back(temp_val_theta);
        frm_no_vect_final.push_back(temp_frm_no_val);
        newmultflag_final.push_back(temp_val_multFlag);
        }
        //LOG(INFO)<<temp_frm_no_val<<" "<<frm_no_vect_final[frm_no_vect_final.size()-1];
    }

LOG(INFO)<<"Reading attitude data from "<<(char*)dirobj.attfile.c_str ();
 double tstart=datainfo.getTstart ();
 double tstop=datainfo.getTstop ();
    int status=0;
   // LOG(INFO)<<xShift_track_final.size();exit(1);
 
    
    
    float *bad_shifted_final_Arr = new float[2*PIX_PADD_SIZE*2*PIX_PADD_SIZE];
    for(int j=0;j<2*PIX_PADD_SIZE*2*PIX_PADD_SIZE;j++)
    {
       bad_shifted_final_Arr[j] = 0.0;
    }
    
   
//    status= performSubDivisionIM (new_badPixarr,600,600,badPix_Arr_fi,PIX_PADD_SIZE*2,PIX_PADD_SIZE*2);
//   if (status)
//    {
//        LOG (ERROR) << "Error in performing sub division in badpixel data" ;
//        return (EXIT_FAILURE) ;
//    }

    
    unsigned short start_firstFrame=frmno[0];
    int newX=0,newY=0;
    int midd_x=0;
    int midd_y=0;
    int num_frm_integrate=nFrameIntegrate_fi;
    float costheta=0.0,sintheta=0.0;
    vector<FrameIntegration_Arr_strct> frm_exp_track;
    FrameIntegration_Arr_strct obj_exp;
//    float *Expo_Arr= new float[IMG_DIM_FI*IMG_DIM_FI];
//    for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++)
//    {
//    Expo_Arr[i]=0.0f;    
//    }
    float *Expo_Arr= new float[PIX_PADD_SIZE*2*PIX_PADD_SIZE*2];
    for (int i=0;i<PIX_PADD_SIZE*2*PIX_PADD_SIZE*2;i++)
    {
    Expo_Arr[i]=0.0f;    
    }
    //LOG(INFO)<<frm_no_vect_final[xShift_track_final.size()-1];
    if(num_frm_integrate>(int)frm_no_vect_final[xShift_track_final.size()-1])
    {
        LOG(ERROR)<<"***Number of frames to be integrate are larger than total available frames***";
        return(EXIT_FAILURE);
    }
   // bool enough_frm_found=FALSE;
    LOG(INFO)<<"Exposure Frame calculation started....";
    //LOG(INFO)<<frm_no_vect_final[0];
   // int cnty=0;
//    for (int i=1;i<frm_no_vect_final.size();i++){
//        
//        if(frm_no_vect_final[i]!=frm_no_vect_final[i-1]+1)
//      LOG(INFO)<<frm_no_vect_final[i];
//            //cnty++;
//    }
   // LOG(INFO)<<cnty;
    status = setDirectoryStructure (moduleoutdir_expo , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             continue;
        }
    
    int temp_numFrameDiscard=nFrameDiscard_fi;
    int cntx=0;
    float *finalSubDivided_arr_exposure = new float[IMG_DIM_FI*IMG_DIM_FI];
   // LOG(INFO)<<"The flag value is "<<flag_Theta;
    long sizeExpoArr=PIX_PADD_SIZE*2;
    long Half_sizeExpoArr=sizeExpoArr/2;
    float xShifts_New=0.0,yShifts_New=0.0;
    float add_term_X=0.0,add_term_Y=0.0;
    int dim=0;
    int TotFrameAvail=0;
    int cnt_forExp=0;
    vector<int> TotFrameAvai_track;
    int numToAdd = (int)round((1/datainfo.getIntegrationTime()))/2;
    //LOG(INFO)<<numToAdd<<" "<<(1/datainfo.getIntegrationTime());
   // rotAng=0;
    double   ctheta =0.0f;//= cos (1.0 *rotAng* M_PI / 180) ;
   double    stheta =0.0f ;//= sin (1.0 * rotAng * M_PI / 180) ;
	// multFactor_ZeroCentroid=1+((double)nzerocentroid/(double)ntotal);
    multFactor_ZeroCentroid=(double)ntotal/((double)ntotal-(double)nzerocentroid);
    
    for (int i =0;i<xShift_track_final.size();i++)
    {
        if(frm_no_vect_final[i]>=start_firstFrame+nFrameDiscard_fi)
        {
            
            cntx++;
            TotFrameAvail=0;
           // LOG(INFO)<<frm_no_vect_final[i]<<" "<<temp_numFrameDiscard+num_frm_integrate;
           time_t  st =time(NULL);
            while(frm_no_vect_final[i]<temp_numFrameDiscard+num_frm_integrate  && i< xShift_track_final.size())
            { 
                cnt_forExp++;
                TotFrameAvail++;
              //  if((cnt_forExp%numToAdd ==0 && cnt_forExp>0) ||( frm_no_vect_final[i]+numToAdd> temp_numFrameDiscard+num_frm_integrate) )
               // {
//                xShifts_New=xShift_track_final[i]*2;  
//                yShifts_New=yShift_track_final[i]*2;  
                 xShifts_New=xShift_track_final[i]*2;  
                yShifts_New=yShift_track_final[i]*2;  
                costheta=cos(1.0*thetaShift_track_final[i]*M_PI/180);
                sintheta=sin(1.0*thetaShift_track_final[i]*M_PI/180);
                add_term_X=Half_sizeExpoArr-xShifts_New;
                add_term_Y=Half_sizeExpoArr-yShifts_New;
                for (int index_y=0;index_y<sizeExpoArr;index_y++)
                {            
                    dim=index_y*(sizeExpoArr);
                   for (int index_x=0;index_x<sizeExpoArr;index_x++)
                   {
                       
                       
                       if(flag_Theta==1)
                       {
                           midd_x=index_x-Half_sizeExpoArr;
                           midd_y=index_y-Half_sizeExpoArr;
                     //costheta=;
                    // sintheta=sin(-1.0*thetaShift_track_final[i]*M_PI/180);
                    newX =(int) round ((midd_x) * costheta - (midd_y) * sintheta + add_term_X) ; //new index x
                    newY= (int)round ((midd_x) * sintheta + (midd_y) * costheta + add_term_Y) ; //new index y
                       }
                       else
                       {
                    newX = index_x-xShifts_New; //new index x
                    newY=  index_y-yShifts_New; //new index y
                       }
                    //bad_shifted_final_Arr[newY*IMG_DIM_FI+newX]=badPix_Arr_fi[index_y*IMG_DIM_FI+index_x];
//%#Added ON 20July17#% 
                     if(newY>0 && newY<sizeExpoArr && newX>0 && newX<sizeExpoArr)
//%#-Till this-20July17#%
                    bad_shifted_final_Arr[newY*(sizeExpoArr)+newX]=exp_templateArray[dim+index_x];
                   }    
                }
                
                for (int index_expo=0;index_expo<sizeExpoArr*sizeExpoArr;index_expo++)
                {
                    Expo_Arr[index_expo]=Expo_Arr[index_expo]+bad_shifted_final_Arr[index_expo];
                    bad_shifted_final_Arr[index_expo]=0.0;
                }
              //  }
//                else{
//                    for (int index_expo=0;index_expo<sizeExpoArr*sizeExpoArr;index_expo++)
//                     {
//                     Expo_Arr[index_expo]=Expo_Arr[index_expo]+badPix_Arr_fi[index_expo];
//                    }
//                }
               
                i++;
                                
            }
           time_t et=time(NULL);
           LOG(INFO)<<"Total time "<<et-st;
           TotFrameAvai_track.push_back(TotFrameAvail);
           //vector<float> TotFrameCount;
           
//             status = setDirectoryStructure ("abcd" , " ") ;
//        if (status)
//        {
//            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
//            continue;
//        }
//    writeOutputImageToDisk ("fi" , "abcd" , "" , "" , Expo_Arr , "abc" ,0.0 ,1 , 1200 , 1200) ; //this is for the SignalFram
            initArray<float>(finalSubDivided_arr_exposure,IMG_DIM_FI*IMG_DIM_FI,0.0);
            status= performSubDivisionIM (Expo_Arr,PIX_PADD_SIZE*2,PIX_PADD_SIZE*2,finalSubDivided_arr_exposure,IMG_DIM_FI,IMG_DIM_FI);
   if (status)
    {
        LOG (ERROR) << "Error in performing sub division in badpixel data" ;
        return (EXIT_FAILURE) ;
    }
//           
//             float *Rotated_Expo_Array= new float[IMG_DIM_FI*IMG_DIM_FI];
//        initArray<float>(Rotated_Expo_Array,IMG_DIM_FI*IMG_DIM_FI,0.0);
//        int new_Xloc,new_Yloc;
//        int mid_X_new,mid_Y_new;
//      for (int i = 0 ; i < IMG_DIM_FI ; i++)
//    {    
//          mid_X_new=i-IMG_DIM_FI/2;
//          for (int j=0;j<IMG_DIM_FI;j++)
//          {  
//              mid_Y_new=j-IMG_DIM_FI/2;
//             
//              new_Xloc=(-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta)))+IMG_DIM_FI/2;
//              new_Yloc=(-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) +IMG_DIM_FI/2;      
//              if(new_Xloc>0 && new_Xloc<IMG_DIM_FI && new_Yloc>0 && new_Yloc<IMG_DIM_FI)
//              Rotated_Expo_Array[(int)round(new_Yloc*IMG_DIM_FI+new_Xloc)]=finalSubDivided_arr_exposure[j*IMG_DIM_FI+i];
//          }         
//    }
//      
            for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++){

finalSubDivided_arr_exposure[i]=finalSubDivided_arr_exposure[i]*multFactor_ZeroCentroid;
         }
		

            status = writeOutputImageToDisk ("fi" , moduleoutdir_expo , "" , "exp" , finalSubDivided_arr_exposure , nameprefix ,0.0 ,cntx , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output
       // status = writeOutputImageToDisk ("fi" , moduleoutdir_expo , "" , "exp" , Rotated_Expo_Array , nameprefix ,0.0 ,cntx , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output  
        if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                continue;
            }
            //maximum element
            delete[] finalSubDivided_arr_exposure;
            //delete[] Rotated_Expo_Array;
          
               
                temp_numFrameDiscard=0;
            for (int i=0;i<PIX_PADD_SIZE*PIX_PADD_SIZE*2*2;i++)
            {
             //obj_exp.Frame_pixels_exp.push_back(Expo_Arr[i]);
             Expo_Arr[i]=0.0f;    
            }
            //frm_exp_track.push_back(obj_exp);
           // obj_exp.Frame_pixels_exp.clear();
                   
           // }
            num_frm_integrate=num_frm_integrate+nFrameIntegrate_fi+nFrameDiscard_fi;
            
        }       
    }
   delete[] bad_shifted_final_Arr;
   // exit(1);
    
    if((nrows-track_rown.size ())<=0)
    {
        LOG(ERROR)<<"Total number of rows remains is  0 ,NO shifts found accordingly to RAS file time";
        continue;    
    }
  //  LOG(INFO)<<"ENTER"<<" "<<nrows<<" "<<track_rown.size ()<<endl;

delete[] Expo_Arr;
    float *temp_X_newsnr = new float[nrows - track_rown.size ()] ;
    float *temp_Y_newsnr = new float[nrows - track_rown.size ()] ;
    unsigned short *temp_frameno_newsnr = new unsigned short[nrows - track_rown.size ()] ;
    double *temp_ftime_newsnr = new double[nrows - track_rown.size ()] ;
    float *temp_enp__newsnr = new float[nrows - track_rown.size ()] ;
    unsigned short *temp_mult_newsnr = new unsigned short[nrows - track_rown.size ()] ;
    unsigned short *temp_bflag_newsnr = new unsigned short[nrows - track_rown.size ()] ;
 	 double *temp_uvittime_newsnr = new double[nrows - track_rown.size ()] ;
		unsigned char *dmm_newsnr= new unsigned char[nrows - track_rown.size ()]  ;
	unsigned char  *min_newsnr = new unsigned char[nrows - track_rown.size ()] ;

//applyinh CRC factor and parity factor;
//LOG(INFO)<<crcfailedpckts<<" "<<pcktremained<<" "<<parityfailedEvents<<" "<<eventsremained;
if(pcktremained>0 && eventsremained>0){
CrcFactor=(1 / (1 - ((double)crcfailedpckts / (double)pcktremained)));
ParityFactor=(1 / (1 - ((double)parityfailedEvents /(double) eventsremained)));
//LOG(INFO)<<CrcFactor<<" "<<ParityFactor;exit(1);
  for(int i=0;i<nrows;i++){
//effective_NumPhotons[i]=effective_NumPhotons[i]*CrcFactor*ParityFactor;
effective_NumPhotons[i]=effective_NumPhotons[i]*ParityFactor;

}
}

        status = setDirectoryStructure (moduleoutdir_snr , " ") ;
        if (status)
        {
            LOG(ERROR) << "***Error in setting shift and rotate***" << endl ;
             continue;
        }
        //writing
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_snr , nameprefix , "snr") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

        status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , xi_final , nrows , TFLOAT , 5 , yi_final ,
                nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows ,TDOUBLE,8,time_frame_uvt,nrows, TUSHORT , 9 , badFlag_temp , nrows , TUSHORT , 10 , mult_temp , nrows , TFLOAT , 11 , effective_NumPhotons , nrows) ;
        if (status)
        {
            LOG (INFO) << "Error in writing to the  optical Assmbly Dist  Correction File" << endl ;
             continue;
        }
   
        fitsfile *fptr ;
        fits_open_file (&fptr , outfile , READWRITE , &status) ;
        printError (status , " Error in opening file " , outfile) ;
        fits_movabs_hdu (fptr , 2 , NULL , &status) ;
        printError (status , "Error in moving to HDU " , outfile) ;
        fits_delete_rowlist (fptr , track_rown.data () , track_rown.size () , &status) ;
        printError (status , "Error in Deleting the row list " , outfile) ;
        fits_read_col (fptr , TUSHORT , 2 , 1 , 1 , nrows - track_rown.size () , NULL , temp_frameno_newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fptr , TDOUBLE , 3 , 1 , 1 , nrows - track_rown.size () , NULL , temp_ftime_newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fptr , TFLOAT , 5 , 1 , 1 , nrows - track_rown.size () , NULL , temp_Y_newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fptr , TFLOAT , 4 , 1 , 1 , nrows - track_rown.size () , NULL , temp_X_newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
	 fits_read_col (fptr , TBYTE , 6 , 1 , 1 , nrows - track_rown.size () , NULL , dmm_newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
	 fits_read_col (fptr , TBYTE , 7 , 1 , 1 , nrows - track_rown.size () , NULL , min_newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        //           fits_read_col (fptr , TFLOAT , 5 , 1 , 1 , nrows-track_rown.size () , NULL , temp_Y_newsnr, NULL , &status) ;
        //           printError (status , "Error in reading the number of rows",file_in) ;
        fits_read_col (fptr , TDOUBLE , 8 , 1 , 1 , nrows - track_rown.size () , NULL , temp_uvittime_newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
	
	fits_read_col (fptr , TUSHORT , 9 , 1 , 1 , nrows - track_rown.size () , NULL , temp_bflag_newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fptr , TUSHORT , 10 , 1 , 1 , nrows - track_rown.size () , NULL , temp_mult_newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_read_col (fptr , TFLOAT , 11, 1 , 1 , nrows - track_rown.size () , NULL , temp_enp__newsnr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , file_in) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error closing the file " , outfile) ;
     
        filenamesnr=outfile;
        //now creating images
        float *image_snr = new float[xsize*ysize];
        for(int i=0;i<xsize*ysize;i++){
            image_snr[i]=0.0f;
        }
        for (int i=0;i<nrows-track_rown.size();i++)
        {
            if((int)(round(temp_Y_newsnr[i])*xsize+round(temp_X_newsnr[i])) >0 && (int)(round(temp_Y_newsnr[i])*xsize+round(temp_X_newsnr[i]))<xsize*ysize )
         image_snr[(int)(round(temp_Y_newsnr[i])*xsize+round(temp_X_newsnr[i]))]=image_snr[(int)(round(temp_Y_newsnr[i])*xsize+round(temp_X_newsnr[i]))]+temp_enp__newsnr[i]*temp_mult_newsnr[i]*temp_bflag_newsnr[i];    
        
        }
        //image file creation after shiftNRot module
//        status = writeOutputImageToDisk ("fi" , moduleoutdir_snr , "" , "img" , image_snr , nameprefix ,9999 ,1 , xsize , ysize) ; //this is for the SignalFrame output
//            if (status)
//            {
//                LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
//                continue;
//            }
        
        //till this
        
//added
bool flag_rot=FALSE;
float i1=0,j1=0;
double theta1,theta2,y1,y2;
 t1=0.0 , t2=0.0 , theta1=0.0 , theta2=0.0 , x1=0.0 , x2=0.0 , y1=0.0 , y2=0.0 ;
//float  ctheta,stheta;
double curr_frm_time;
double new_delta_theta = 0.0 , new_delta_x = 0.0 , new_delta_y = 0.0 ;

//   
      float *temp_ra_newsnr = new float [nrows-track_rown.size()];
      float * temp_dec_newsnr= new float[nrows-track_rown.size()];
//        float cdelt1,cdelt2;
//      float factor_delta=xsize/600;
//      if(strcmp(datainfo.getDetector (),"VIS")==0)
//    {
//           cdelt1=(3.357/3600)/factor_delta;
//           cdelt2=(3.311/3600)/factor_delta;
//    }
//    else if (strcmp(datainfo.getDetector (),"FUV")==0){
//        //cdelt1=3.357/3600;
//       //  cdelt2=3.311/3600;
//       cdelt1=cdelt2=(3.3373/3600)/factor_delta;
//    }
//    else if(strcmp(datainfo.getDetector (),"NUV")==0){
//       cdelt1=cdelt2=(3.3307/3600)/factor_delta;
//       // cdelt1=cdelt2=(3.3307*factor_delta/3600)/;
//       
//    }
//   // cdelt1=-cdelt1/cos(center_dec*M_PI/180);
//     cdelt1=-cdelt1;///cos(center_dec*M_PI/180);
//       ctheta = cos (1.0 *rotAng* M_PI / 180) ;
//        stheta = sin (1.0 * rotAng * M_PI / 180) ;
//      
//    for (int i = 0 ; i < nrows - track_rown.size() ; i ++)
//    {
//           
//        flag_rot = FALSE ;
//        curr_frm_time = temp_ftime_newsnr[i] ;
//            
//            
//        i1 = temp_X_newsnr[i] - xsize / 2 ;
//        j1 = temp_Y_newsnr[i] - ysize / 2 ;
//        //calculating the new values for Xand Y based on the deltas
//        if ((((i1 * ctheta) - (j1 * stheta)) + xsize / 2 )>- 1 && (((i1 * ctheta) - (j1 * stheta)) + xsize / 2 ) < xsize && (((i1 * stheta) + (j1 * ctheta)) + ysize / 2 )>- 1 &&  (((i1 * stheta) + (j1 * ctheta)) + ysize / 2 ) < xsize)
//        {
//            //for NUV date-01sept16 flip signs of first term(i1) in both the equation
//           
//            temp_X_newsnr[i] = (-(i1 * ctheta) - (j1 * stheta)) + xsize / 2; //-new_delta_x ; //new index x
//            temp_Y_newsnr[i] = (-(i1 * stheta) +(j1 * ctheta)) + ysize / 2; //-new_delta_y ; //new index y
//            temp_ra_newsnr[i]=RA_pnt+(temp_X_newsnr[i]-xsize/2)*(cdelt1/cos(DEC_pnt*M_PI/180));
//            temp_dec_newsnr[i]=DEC_pnt+(temp_Y_newsnr[i]-ysize/2)*cdelt2;
//          
//        }
//    }
//     
         
    delete[] xi_final , yi_final ;
 
    LOG(INFO)<<"Performing Frame integration";
     //Allocating memory for frame integration

      vector<FrameIntegration_Arr_strct> frameIntegration_track ;//vector for storing the output of frame integration to the disk
    float *one_dim_exp ;
    one_dim_exp = new float[IMG_DIM_FI * IMG_DIM_FI] ;
    float *one_dim_img ;
    one_dim_img = new float[IMG_DIM_FI * IMG_DIM_FI] ;
     float *one_dim_exp_final ;
    one_dim_exp_final = new float[IMG_DIM_FI * IMG_DIM_FI] ;
    float *one_dim_img_final ;
    one_dim_img_final = new float[IMG_DIM_FI * IMG_DIM_FI] ;
    float *noise_Array_output = new float[IMG_DIM_FI*IMG_DIM_FI];

  
     
    fpixel[0] = fpixel[1] = 1 ;
    
    fits_open_file (&fptr , (char*)frameIntegrationfrmname[0].c_str() , READWRITE , &status) ;
    printError (status , "Error in reading the fits file" , (char*)frameIntegrationfrmname[0].c_str()) ;
    
    for (int q = 0 ; q < IMG_DIM_FI * IMG_DIM_FI ; q ++) one_dim_exp[q] = 0.0 ;

    fits_read_pix (fptr , TFLOAT , fpixel , IMG_DIM_FI * IMG_DIM_FI , NULL , one_dim_exp , NULL , &status) ;
    printError (status , "Error in reading the pixels from caldb file" , (char*)frameIntegrationfrmname[0].c_str()) ;
   
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the file" , (char*)frameIntegrationfrmname[0].c_str()) ;
    
    
    //for registration and averaging 
    mult_fact = subDiv_fact ;
	
    refine_Winsize = refine_Winsize * mult_fact ;
    centroid_Winsize = centroid_Winsize * mult_fact;
    //registering  and Averaging started
    vector<float> cx_ref , cy_ref , ci_ref ;
    int min_stars_match ;
    float *avg_sigArray = new float[subDivision_size * subDivision_size] ;
    float *avg_expArray = new float[subDivision_size * subDivision_size] ;
    float *noice_mapArray= new float[subDivision_size * subDivision_size];
    frame_sig_Data = new float[subDivision_size * subDivision_size] ;
    frame_exp_Data = new float[subDivision_size * subDivision_size] ;
  
    
    for (int i = 0 ; i < subDivision_size * subDivision_size ; i ++)
    {
        frame_sig_Data[i]=0.0f;
        frame_exp_Data[i]=0.0f;
        avg_sigArray[i] = 0 ;
        avg_expArray[i] = one_dim_exp[i] ;
        noice_mapArray[i]=0.0f;
    }
    bool flag_firstFrame_accomplished=FALSE;
   

   LOG(INFO)<<"Performing Registration and averaging";
   float frm_time_fi_new;
   long frmno_fi_new;
   long num_frm_count=0;
   long tmp_startrow=0;
   long currframe_end_index=0;
   long nextframe_start_index=0;
   vector<float> x_ref_arr , y_ref_arr , x_arr , y_arr , temp_x_arr , temp_y_arr ,int_ref_arr,int_arr;
    double x_dx = 0.0 , y_dy = 0.0 , theta_dt = 0.0 ;
     int cnt = 0 ; 
     int matching_points = 0 ;
       ctheta = 0.0f , stheta = 0.0f ;
       int x_index ,y_index;
         float new_index_x,new_index_y;
//         float t1,t2;
         ofstream fout1("xytheta.txt");
bool flag_onlyoneframe=FALSE;
int cnt_frmExpo=-1;
int Loc_end=nrows-track_rown.size()-1;
float *temp_SigArryWithDiv;
temp_SigArryWithDiv= new float[subDivision_size*subDivision_size];//Array for storing Sig/Exp 
for (int i=0;i<subDivision_size*subDivision_size;i++)
{
 temp_SigArryWithDiv[i]=0.0f;    
}
float sum_OfPixelsOfExp=0.0f;
    for (int i = 0 ; i < nrows-track_rown.size () ; i =i+(currframe_end_index-nextframe_start_index))
    {    cnt_frmExpo++;
	sum_OfPixelsOfExp=0.0f;
        num_frm_count++;
        sd_mul_factor = sd_multi_factor_default ;
       // LOG(INFO)<<nextframe_start_index;
        if (i==0) //incase of reference frame  for registration and averaging
        {
            t1=time_t(NULL);
           status = performFrameIntegration ((long)(nrows - track_rown.size ()) , tmp_startrow,temp_frameno_newsnr , xsize , nFrameDiscard_fi , nFrameIntegrate_fi , temp_X_newsnr , temp_Y_newsnr , temp_mult_newsnr , temp_bflag_newsnr , temp_ftime_newsnr , temp_enp__newsnr , one_dim_img_final , one_dim_exp_final , frameIntegration_track,noise_Array_output,one_dim_img,one_dim_exp,frm_time_fi_new,frmno_fi_new,currframe_end_index) ;
//           for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++)
//            {
//            one_dim_img[i]=one_dim_img[i]/TotFrameAvai_track[cnt_frmExpo]    ;
//            }
          // LOG(INFO)<<TotFrameAvai_track[cnt_frmExpo]<<" "<<moduleoutdir_fi<<" "<<status;
           if (wtd_fi)//incase of frame integration to be written on the disk
           {
        status = setDirectoryStructure (moduleoutdir_fi , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            continue;
        }
        status = setDirectoryStructure (moduleoutdir_fi , "SignalFrames_DividedWithExposure") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             continue;
        }
       status = setDirectoryStructure (moduleoutdir_fi , "NOISE_MAP") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             continue;
        } 
         
        
            status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "SignalFrames" , "sig" , one_dim_img , nameprefix ,frm_time_fi_new ,num_frm_count , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                continue;
            }
          //  LOG(INFO)<<cnt_frmExpo;
         
           
           }      
            t2=time_t(NULL);
           
            if(fi_flag==TRUE){
            status = findStar_algo1 (one_dim_img) ; // for first frame taken as  reference
            if (status)
            {
                 LOG (ERROR) << endl << "***Error in finding star algorithm  for frame  " << infile << "  ***" << endl ;
                 LOG(ERROR)<<"CRASH NO STAR FOUND (POSSIBLY DUE TO SIGMA MULTIPLIER =< 0) (uvtLevel2PC.cpp)";
                return (EXIT_FAILURE) ;
            }
            //LOG(INFO)<<Cx.size ()<<endl;
           
         }
            for (int i = 0 ; i < Cx.size () ; i ++)
            {
                cx_ref.push_back (Cx[i]) ;
                cy_ref.push_back (Cy[i]) ;
                ci_ref.push_back (Ci[i]) ;
             
            }
           
                    
           // fout.close ();
            // exit(1);
            Cx.clear () ;
            Cy.clear () ;
            Ci.clear () ;
           
            for(int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++)
            {
            avg_sigArray[i] =one_dim_img[i];//*avg_expArray[i]; 
//            noice_mapArray[i]=noise_Array_output[i];
            }
            //}
            
              for (int p = 0; p < subDivision_size * subDivision_size; p++) {
                if (avg_expArray[p] != 0) {
                    temp_SigArryWithDiv[p]=0;
                    temp_SigArryWithDiv[p]=one_dim_img[p]/avg_expArray[p];
                    noise_Array_output[p]=sqrt(noise_Array_output[p]*avg_expArray[p]*datainfo.getIntegrationTime())/(avg_expArray[p]*datainfo.getIntegrationTime());
                }
            }
 if (wtd_fi)//incase of frame integration to be written on the disk
            {

                status = writeOutputImageToDisk("Sig_DividedWithExp", moduleoutdir_fi, "SignalFrames_DividedWithExposure", "sig", temp_SigArryWithDiv, nameprefix, frm_time_fi_new, num_frm_count, subDivision_size, subDivision_size); //this is for the SignalFrame output
                if (status) {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl;
                    continue;
                }
                 status = writeOutputImageToDisk("noise_map_sig", moduleoutdir_fi, "NOISE_MAP", "sig", noise_Array_output, nameprefix, frm_time_fi_new, num_frm_count, subDivision_size, subDivision_size); //this is for the SignalFrame output
                if (status) {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl;
                    continue;
                }

            }
            
            //exit(1);
            flag_firstFrame_accomplished=TRUE;
        }
        else if(flag_firstFrame_accomplished==TRUE)//for remaining frame
        {
           
            nextframe_start_index=currframe_end_index;
            nFrameDiscard_fi=0;//bcoz it is already discarded in first iteration
            status = performFrameIntegration ((long)(nrows - track_rown.size ()) , nextframe_start_index,temp_frameno_newsnr , xsize , nFrameDiscard_fi , nFrameIntegrate_fi , temp_X_newsnr , temp_Y_newsnr , temp_mult_newsnr , temp_bflag_newsnr , temp_ftime_newsnr , temp_enp__newsnr , one_dim_img_final , one_dim_exp_final , frameIntegration_track,noise_Array_output,one_dim_img,one_dim_exp,frm_time_fi_new,frmno_fi_new,currframe_end_index) ;
//         for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++)
//            {
//            one_dim_img[i]=one_dim_img[i]/TotFrameAvai_track[cnt_frmExpo]    ;
//            }
          //  LOG(INFO)<<TotFrameAvai_track[cnt_frmExpo];
            if (wtd_fi)//incase of frame integration to be written on the disk
             {
    
            status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "SignalFrames" , "sig" , one_dim_img , nameprefix ,frm_time_fi_new , num_frm_count , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                continue;
            }
          
   
             }
            
            //reading exposure frames
            
            LOG(INFO)<<"\033[1;34mRegistration is applying on frame number "<<num_frm_count<<"\033[0m";
          //  if(num_frm_count==214) goto labelrg;
             
            
             cnt = 0 ;
             matching_points = 0 ;
            status = findStar_algo1 (one_dim_img) ;//for finding stars
            if (status)
            {
                LOG (INFO) << "***Error in finding star algorithm 1 for frame  " << infile << "  ***" ;
		LOG(ERROR)<<"CRASH NO STAR FOUND (POSSIBLY DUE TO SIGMA MULTIPLIER =< 0) (uvtLevel2PC.cpp)";
                return (EXIT_FAILURE) ;
            }
            Rx.clear () ;
            Ry.clear () ;
            Rval.clear () ;
            Fx.clear () ;
            Fy.clear () ;
            Fval.clear () ;
            char name[100];
            sprintf(name,"File_%d",num_frm_count);
            
            if (Cx.size () >= cx_ref.size ())
            {
                matching_points = cx_ref.size () ;
            }
            else
            {
                matching_points = Cx.size () ;
            }

            
            x_ref_arr.clear ();y_ref_arr.clear ();x_arr.clear ();y_arr.clear ();temp_x_arr.clear ();temp_y_arr.clear ();int_ref_arr.clear();int_arr.clear();
            //calculating how many stars to be matched with reference frame.
            min_stars_match = (int) ((100 - PERCENTGE_ERR_ALLOWED) * matching_points / 100) ;
            cnt = 0 ;
            
            //finding stars        
            
            while (cnt < min_stars_match)
            {

                x_ref_arr.clear () ;
                y_ref_arr.clear () ;
                x_arr.clear () ;
                y_arr.clear () ;
                temp_x_arr.clear () ;
                temp_y_arr.clear () ;
                int_ref_arr.clear();
                int_arr.clear();
                cnt = matchStars (cx_ref.size () , Cx.size () , mult_fact , cx_ref.data () , cy_ref.data () , Cx.data () , Cy.data () ,ci_ref.data(),Ci.data(), x_ref_arr , y_ref_arr , x_arr ,
                        y_arr ,int_ref_arr,int_arr, temp_x_arr , temp_y_arr) ;
              //  LOG(INFO)<<cnt;
                diff_Dist = diff_Dist * 2 ;
               
            }
            char name1[100];

            diff_Dist = 1 ;
            //  cx_ref.clear ();cy_ref.clear ();
             x_dx = 0.0 , y_dy = 0.0 , theta_dt = 0.0 ;

            /* option_LeastSquare parameter decide which algorithm to use for finding the Drifts between  points 
           these are  the different techniques for finding shifts between two point of two frames.
             */
            vector<float> New_X_ref,New_Y_ref,New_x_arr,New_y_arr,New_Xdiff,New_Ydiff;
             vector<float> int_new_one,int_new_two;
            status=removeRecords(x_ref_arr,y_ref_arr,x_arr,y_arr,temp_x_arr,temp_y_arr,int_ref_arr.data(),int_arr.data(),New_X_ref,New_Y_ref,New_x_arr,New_y_arr,New_Xdiff,New_Ydiff,int_new_one,int_new_two);
          if (status)
           {
               LOG (ERROR) << "Error in removing the records above mean+sigma" ;
               return (EXIT_FAILURE) ;
           }

           
             sprintf(name1,"File_%d.txt",num_frm_count);
            ofstream fout(name1);
            for (int i=0;i<x_ref_arr.size ();i++)
            {
                fout<<New_X_ref[i]<<" "<<New_Y_ref[i]<<" "<<New_x_arr[i]<<" "<<New_y_arr[i]<<" "<<New_Xdiff[i]<<" "<<New_Ydiff[i]<<endl;
            }
            
            fout.close ();
             
            status = findShiftsNtheta (New_X_ref.size () , New_X_ref , New_Y_ref ,int_new_one,New_x_arr , New_y_arr , int_new_two,New_Xdiff, New_Ydiff , flag_thetaComp,x_dx , y_dy , theta_dt) ;
            if (status)
            {
                LOG (INFO) << "Error in finding shifts n theta " << endl ;
                return (EXIT_FAILURE) ;
            }
           // LOG(INFO)<<x_dx<<" "<<y_dy<<" "<<theta_dt;      
             fout1<<x_dx<<" "<<y_dy<<" "<<theta_dt<<endl;;      
     
                    
            ctheta = 0.0f , stheta = 0.0f ;
            ctheta = cos (-1.0 * theta_dt) ;
            stheta = sin (-1.0 * theta_dt) ;
           
            LOG (INFO) <<"Loop Started for assigning the correction to the frames.."  ;

            for (int i = 0 ; i < subDivision_size * subDivision_size ; i ++)
            {
                frame_sig_Data[i] = 0 ;
                frame_exp_Data[i] = 0 ;
            }

            //loop for applying Shifts to get registered image
             fitsfile *fptr ;
            status = 0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
   // LOG(INFO)<<frameIntegrationfrmname[num_frm_count-1];
    fits_open_file (&fptr , (char*)frameIntegrationfrmname[num_frm_count-1].c_str() , READONLY , &status) ;
    printError (status , "Error in reading the fits file" , (char*)frameIntegrationfrmname[num_frm_count].c_str()) ;

    
    for (int q = 0 ; q < IMG_DIM_FI * IMG_DIM_FI ; q ++) one_dim_exp[q] = 0.0 ;
    fits_read_pix (fptr , TFLOAT , fpixel , IMG_DIM_FI * IMG_DIM_FI , NULL , one_dim_exp , NULL , &status) ;
    printError (status , "Error in reading the pixels from caldb file" , (char*)frameIntegrationfrmname[num_frm_count].c_str()) ;
   // LOG(INFO)<<"File close";
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the file" , (char*)frameIntegrationfrmname[num_frm_count].c_str()) ;
            
    //LOG(INFO)<<"Entered in First stage";
            for (int index = 0 ; index < subDivision_size ; index ++)
            {
                x_index = index - subDivision_size / 2 ;
                for (int jindex = 0 ; jindex < subDivision_size ; jindex ++)
                {
                    if (one_dim_img[jindex * subDivision_size + index] != INVALID_PIX_VALUE && one_dim_exp[jindex * subDivision_size + index] != INVALID_PIX_VALUE)
                    {

                         y_index = jindex - subDivision_size / 2 ;
                         //ROUNDING:Correction needed.previously rounding was only on the subset of full equation.
                         new_index_x = round ((x_index) * ctheta - (y_index) * stheta + subDivision_size / 2 - x_dx* mult_fact) ; //new index x
                         new_index_y = round ((x_index) * stheta + (y_index) * ctheta + subDivision_size / 2 - y_dy* mult_fact) ; //new index y

                        if (round (new_index_x) < subDivision_size && round (new_index_x) > 0 && round (new_index_y) > 0 && round (new_index_y) < subDivision_size)
                        {
                            // cnt_loop ++ ;
                            //Rounding:No correction needed.As rounding is applied on each direction i.e X and Y.
                            frame_sig_Data[(int) (round (new_index_y) * subDivision_size + round (new_index_x))] = one_dim_img[jindex * subDivision_size + index] ;
                            frame_exp_Data[(int) (round (new_index_y) * subDivision_size + round (new_index_x))] = one_dim_exp[jindex * subDivision_size + index];
                        }
                    }
                   
                }
            }
    
           for (int p = 0; p < subDivision_size * subDivision_size; p++) {
                if ( frame_exp_Data[p] != 0) {
                    temp_SigArryWithDiv[p]=0;
                    temp_SigArryWithDiv[p]=frame_sig_Data[p]/frame_exp_Data[p];
                    noise_Array_output[p]=sqrt(noise_Array_output[p]*frame_exp_Data[p]*datainfo.getIntegrationTime())/(frame_exp_Data[p]*datainfo.getIntegrationTime());
                }
            }
 if (wtd_fi)//incase of frame integration to be written on the disk
            {

                status = writeOutputImageToDisk("Sig_DividedWithExp", moduleoutdir_fi, "SignalFrames_DividedWithExposure", "sig", temp_SigArryWithDiv, nameprefix, frm_time_fi_new, num_frm_count, subDivision_size, subDivision_size); //this is for the SignalFrame output
                if (status) {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl;
                    continue;
                }
                  status = writeOutputImageToDisk("noice_map_sig", moduleoutdir_fi, "NOISE_MAP", "sig", noise_Array_output, nameprefix, frm_time_fi_new, num_frm_count, subDivision_size, subDivision_size); //this is for the SignalFrame output
                if (status) {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl;
                    continue;
                }

            }
            Cx.clear () ;
            Cy.clear () ;
            Ci.clear () ;
}

        //Array for storing the Sig/Exp array.
        if(fi_flag==1){
            for (int p = 0 ; p < subDivision_size * subDivision_size ; p ++)
            {               	
                if (frame_sig_Data[p] != INVALID_PIX_VALUE && frame_exp_Data[p] != INVALID_PIX_VALUE)
                {
                    
//                    avg_sigArray[p] = avg_sigArray[p] + frame_sig_Data[p] * frame_exp_Data[p] ;
                   // avg_sigArray[p] = avg_sigArray[p] + frame_sig_Data[p] *frame_exp_Data[p] ;
                     avg_sigArray[p] = avg_sigArray[p] + frame_sig_Data[p];// *frame_exp_Data[p] ;
                    avg_expArray[p] = avg_expArray[p] + frame_exp_Data[p] ;
             //       noice_mapArray[p]=noice_mapArray[p]+noise_Array_output[p];
                   
                }
            }
        }
     //Writing Sig/Exp Array to output disk.
         
        
        
        
        
        
        }
//float Max_value_Exp=0.0f;
//for (int i=0;i<FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG;i++)
//{
  //  if(Max_value_Exp<avg_expArray[i]){
    //    Max_value_Exp=avg_expArray[i];
   // }   
    
//}

float Max_value_Exp=avg_expArray[(FINALFRAMESIZE_REGAVG/2)*FINALFRAMESIZE_REGAVG+(FINALFRAMESIZE_REGAVG/2)];
int cnter_elements=0;
for (int i=(FINALFRAMESIZE_REGAVG/2)-2;i<=(FINALFRAMESIZE_REGAVG/2)+2;i++)
{
	for (int j=(FINALFRAMESIZE_REGAVG/2)-2;j<=(FINALFRAMESIZE_REGAVG/2)+2;j++){

    if(avg_expArray[j*FINALFRAMESIZE_REGAVG+i]==Max_value_Exp){
       cnter_elements++;
    }   
    }
}

//LOG(INFO)<<"Cner type"<<cnter_elements<<" "<<Max_value_Exp;exit(1);

long PerAvgVal=Max_value_Exp*10/100;
float *expArray= new float[FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG];
initArray(expArray,FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG,0.0f);
for (int i=0;i<FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG;i++)
{
    if(avg_expArray[i]>PerAvgVal)
    {
        expArray[i]=avg_expArray[i];
    }
}
vector<float> NonZero_ExpMatrix;
for(int i=0;i<FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG;i++)
{
    if(avg_expArray[i]>0){
        NonZero_ExpMatrix.push_back(avg_expArray[i]);
                
    }
}
//float Mean_ExpTemplate=getmean(NonZero_ExpMatrix.data(),NonZero_ExpMatrix.size());

//status = writeOutputImageToDisk ("rg" , moduleoutdir_rav , "" , "10Per" , expArray , nameprefix , 0 , 1 , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
peakValueExp=INVALID_PIX_VALUE;
if(cnter_elements!=25){
Total_exp_time=INVALID_PIX_VALUE;
LOG(ERROR)<<"CRASH :EXPOSURE TIME CALCULATION FAILED";
}
else{
peakValueExp=Max_value_Exp;
 Total_exp_time =Max_value_Exp*datainfo.getIntegrationTime();
// Total_exp_time =Mean_ExpTemplate*datainfo.getIntegrationTime();
}
    
         fout1.close();
         labelrg:
    delete[] frame_sig_Data , frame_exp_Data ;
    float *Regavg_subSampled_Array_Sig = new float[FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG] ;
    float *Regavg_subSampled_Array_Exp = new float[FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG] ;
    float *Regavg_NoicemapArray= new float[FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG];
    //loop for getting final registration image by dividing by exposure frame
    //if(fi_flag){
    for (int p = 0 ; p < subDivision_size * subDivision_size ; p ++)
    {   //LOG(INFO)<<avg_sigArray[p]<<" "<<avg_expArray[p];
        if (expArray[p] > 0.0 && avg_sigArray[p] != INVALID_PIX_VALUE){
//            if(avg_sigArray[p]>0.0)
//            LOG(INFO)<<avg_sigArray[p]<<" "<<avg_expArray[p]<<" "<<p;
//            avg_sigArray[p] = (avg_sigArray[p] / avg_expArray[p]) ;
             avg_sigArray[p] = (avg_sigArray[p] / expArray[p]) ;
             noice_mapArray[p]=(sqrt(avg_sigArray[p]*expArray[p]*datainfo.getIntegrationTime())/(expArray[p]*datainfo.getIntegrationTime()));
            
        }
        else{
            avg_sigArray[p]=0.0f;
            noice_mapArray[p]=0.0f;
        }
    }
    //}
  
    for (int i = 0 ; i < FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG ; i ++)
    {
        Regavg_subSampled_Array_Sig[i] = INITIALIZATION_VALUE ;
        Regavg_subSampled_Array_Exp[i] = INITIALIZATION_VALUE ;
        Regavg_NoicemapArray[i]=INITIALIZATION_VALUE;
    }
   //  LOG(INFO)<<"Status"<<status;
    fitsfile *fatt;
    fits_open_file(&fatt,dirobj.attfile.c_str (),READONLY,&status);
    printError(status,"Error opening attitude file",(char*)dirobj.attfile.c_str ());
    
    fits_movabs_hdu(fatt,2,NULL,&status);                           //Attitude data is in 2nd HDU
    printError(status,"Error in moving to 2nd HDU ",(char*)dirobj.attfile.c_str ());
    
    long nrows_att;
    fits_get_num_rows(fatt,&nrows_att,&status);
    printError(status,"Error reading number of rows from attitude file  ",(char*)dirobj.attfile.c_str ());
        // LOG(INFO)<<"Status"<<status;
    long firstrow=1;
    long firstelem=1;
    long nelements=nrows_att;                  //for time
    double *time_att = new double[nelements];
    double *RA_arr= new double[nelements];
    double *DEC_arr= new double[nelements];
    int tcoln,qcoln;               //column number for time and quaternion
    double RA_pnt,DEC_pnt;
//    fits_get_colnum(fatt,CASEINSEN,timecol,&tcoln,&status);
//     printError(status,"Error reading time column number from attitude file  ",(char*)dirobj.attfile.c_str ());   
    
     fits_read_col(fatt,TDOUBLE,1, firstrow,firstelem, nelements, NULL,time_att,NULL,&status); 
     printError(status,"Error reading time column from attitude file  ",(char*)dirobj.attfile.c_str ());     
      fits_read_col(fatt,TDOUBLE,3, firstrow,firstelem, nelements, NULL,RA_arr,NULL,&status); 
     printError(status,"Error reading RA column from attitude file  ",(char*)dirobj.attfile.c_str ());   
     fits_read_col(fatt,TDOUBLE,4, firstrow,firstelem, nelements, NULL,DEC_arr,NULL,&status); 
     printError(status,"Error reading DEC column from attitude file  ",(char*)dirobj.attfile.c_str ());     
    long startindex=0, lastindex=nrows_att-1;
  //  LOG(INFO)<<time_att[0]<<" "<<time_att[nrows_att-1]<<" "<<tstart<<" "<<tstop;
    if(tstart>time_att[nrows_att-1] || tstop<time_att[0])
    {
        LOG(ERROR)<<"\033[1;31m***Time information in data file and attitude file do not match***\033[0m";
        return (EXIT_FAILURE);
    }
 //LOG(INFO)<<"Status"<<status;
//    double time_mid =( tstart+tstop)/2;
    double time_mid =tstart;
   long index_attFile=0;
    for (int i=1;i<nelements;i++)
    {
        if(time_mid>=time_att[i-1] && time_mid<time_att[i])
        {
            index_attFile=i-1;
            RA_pnt=RA_arr[i-1];
            DEC_pnt=DEC_arr[i-1];
            break;
        }    
        
    }
   
  
 //  LOG(INFO)<<RA_pnt<<" "<<DEC_pnt;exit(1);
    double rotAng=0.0f;
   //  LOG(INFO)<<"Status"<<status<<" "<<RA_pnt<<" "<<DEC_pnt;
    fits_read_col(fatt,TDOUBLE,5,index_attFile,1, 1, NULL,&rotAng,NULL,&status); 
    printError(status,"Error reading rotation angle  from attitude file  ",(char*)dirobj.attfile.c_str ());
 fits_close_file(fatt,&status);
    printError(status,"Error in closing the file  ",(char*)dirobj.attfile.c_str ());
if(rotAng <-360 || rotAng >360){

flag_Roll_Angleinvalid=TRUE;
fits_open_file(&fatt,dirobj.attfile.c_str (),READONLY,&status);
    printError(status,"Error opening attitude file",(char*)dirobj.attfile.c_str ());
    
    fits_movabs_hdu(fatt,2,NULL,&status);                           //Attitude data is in 2nd HDU
    printError(status,"Error in moving to 2nd HDU ",(char*)dirobj.attfile.c_str ());
    
  //  long nrows_att;
double *RollAngles = new double[nrows_att];
double *RaPntarr= new double[nrows_att];
double *decPntarr= new double[nrows_att];
    fits_get_num_rows(fatt,&nrows_att,&status);
    printError(status,"Error reading number of rows from attitude file  ",(char*)dirobj.attfile.c_str ());
    fits_read_col(fatt,TDOUBLE,5,1,1,nrows_att, NULL,RollAngles,NULL,&status); 
    printError(status,"Error reading rotation angle  from attitude file  ",(char*)dirobj.attfile.c_str ());
    fits_read_col(fatt,TDOUBLE,3,1,1,nrows_att, NULL,RaPntarr,NULL,&status); 
    printError(status,"Error reading RA  from attitude file  ",(char*)dirobj.attfile.c_str ());
    fits_read_col(fatt,TDOUBLE,4,1,1,nrows_att, NULL,decPntarr,NULL,&status); 
    printError(status,"Error reading DEC  from attitude file  ",(char*)dirobj.attfile.c_str ());
    fits_close_file(fatt,&status);
    printError(status,"Error in closing the file  ",(char*)dirobj.attfile.c_str ());

            for (int i = 0; i < nrows_att; i++) {
                if (RollAngles[i] < 360 && RollAngles[i]>-360) {
                    rotAng = RollAngles[i];
                    break;

                }

            }
     for(int i=0;i<nrows_att;i++){
if (RaPntarr[i]<360 && RaPntarr[i]>-360){
RA_pnt=RaPntarr[i];
DEC_pnt=decPntarr[i];
break;

}

}
    
	

}
     RAPNT_FRMATT=RA_pnt;
   DECPNT_FRMATT=DEC_pnt;

    if(strcmp(datainfo.getDetector(),"NUV")==0)
    {
         rotAng=(rotAng+ANGLE_NUV);
    }
    else if(strcmp(datainfo.getDetector(),"FUV")==0)
    {
       // LOG(INFO)<<"INSIDE";
        //)<<rotAng;
        rotAng=180-(rotAng+ANGLE_FUV);
        //)<<rotAng;
        //rotAng=(rotAng+ANGLE_FUV);
    }
    else 
    {
        LOG(ERROR)<<"***Invalid channel***";
        return(EXIT_FAILURE);
    }
   
    //)<<"Roll Angle "<<rotAng; 
    roll_angle_applied=rotAng;
    
    if(strcmp(datainfo.getDetector(),"NUV")==0){
        
        rotAng=-rotAng;
    }
    //)<<"Status"<<status;
    
    //for getting sub-sampled signal  image      
    status = ApplySubSampling (avg_sigArray , subDivision_size , subDivision_size , Regavg_subSampled_Array_Sig , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG (ERROR) << "ERROR in sub sampling "  ;
        return (EXIT_FAILURE) ;
    }
    //for getting sub-sampled exposure image
    status = ApplySubSampling (avg_expArray , subDivision_size , subDivision_size , Regavg_subSampled_Array_Exp , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG (ERROR) << "ERROR in sub sampling "  ;
        return (EXIT_FAILURE) ;
    }
    
    
     status = ApplySubSampling (noice_mapArray , subDivision_size , subDivision_size , Regavg_NoicemapArray , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG (ERROR) << "ERROR in sub sampling "  ;
        return (EXIT_FAILURE) ;
    }
     //)<<"Status"<<status<<" "<<moduleoutdir_ravFlipped;
    //storing the output of Non flipped image to output directoy.
     status = setDirectoryStructure (moduleoutdir_ravFlipped , "") ;//setting output directory structure for registration and averaging 
    if (status)
    {
        LOG(ERROR) << "***Directory Structure has not been successfully set-up***"  ;
        return (EXIT_FAILURE) ;
    };
    // LOG(INFO)<<"Status"<<status;
     status = writeOutputImageToDisk ("regAvg" , moduleoutdir_ravFlipped , "" , "sig" , Regavg_subSampled_Array_Sig , nameprefix , 0 , 1 , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
     status = writeOutputImageToDisk ("regAvg" , moduleoutdir_ravFlipped , "" , "exp" , Regavg_subSampled_Array_Exp , nameprefix , 0 , 1 , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
     status = writeOutputImageToDisk ("regAvg" , moduleoutdir_ravFlipped , "" , "noiseMap" , Regavg_NoicemapArray , nameprefix , 0 , 1 , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
      //)<<"Status"<<status;
     double mid_time = 0.0f ;
    int cnt_frm_wmean = 1 ;
    status = setDirectoryStructure (moduleoutdir_rav , "") ;//setting output directory structure for registration and averaging 
    if (status)
    {
        LOG(ERROR) << "***Directory Structure has not been successfully set-up***"  ;
        return (EXIT_FAILURE) ;
    };
 //)<<"Status"<<status;
    //Bottom part is commented on 01sept16
    status = writeOutputImageToDisk ("rg" , moduleoutdir_rav , "" , "10Per" , expArray , nameprefix , 0 , 1 , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
if(cnter_elements!=25){
Total_exp_time=INVALID_PIX_VALUE;
LOG(ERROR)<<"CRASH :EXPOSURE TIME CALCULATION FAILED";
}
else{

 Total_exp_time =Max_value_Exp*datainfo.getIntegrationTime();
}    
 //)<<"Status"<<status;
if (strcmp(datainfo.getDetector (),"NUV")==0)
    {
        float* temp_sig_regArr= new float[FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
        float* temp_exp_regArr= new float[FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
        float* temp_noicemap_regArr = new float[FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
        for (int i=0;i<FINALFRAMESIZE_REGAVG;i++)
        {
            for (int j=0;j<FINALFRAMESIZE_REGAVG;j++)
            {
            temp_sig_regArr[i*FINALFRAMESIZE_REGAVG+j]=Regavg_subSampled_Array_Sig[i*FINALFRAMESIZE_REGAVG+((FINALFRAMESIZE_REGAVG-1)-j)];
            temp_exp_regArr[i*FINALFRAMESIZE_REGAVG+j]=Regavg_subSampled_Array_Exp[i*FINALFRAMESIZE_REGAVG+((FINALFRAMESIZE_REGAVG-1)-j)];
            temp_noicemap_regArr[i*FINALFRAMESIZE_REGAVG+j]=Regavg_NoicemapArray[i*FINALFRAMESIZE_REGAVG+((FINALFRAMESIZE_REGAVG-1)-j)];
            }
            
        }
        for (int i=0;i<FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG;i++)
        {
            Regavg_subSampled_Array_Sig[i]=temp_sig_regArr[i];
            Regavg_subSampled_Array_Exp[i]=temp_exp_regArr[i];
            Regavg_NoicemapArray[i]=temp_noicemap_regArr[i];
        }
        
    }
    
    //writing sig and exp images to the output file.
    status = writeOutputImageToDisk ("rg" , moduleoutdir_rav , "" , "sigFlipped" , Regavg_subSampled_Array_Sig , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG(ERROR) << "***Writing to Disk Fails***"  ;
        return (EXIT_FAILURE) ;
    }
    status = writeOutputImageToDisk ("rg" , moduleoutdir_rav , "" , "expFlipped" , Regavg_subSampled_Array_Exp , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG(ERROR) << "***Writing to Disk Fails***"  ;
        return (EXIT_FAILURE) ;
    }
     status = writeOutputImageToDisk ("rg" , moduleoutdir_rav , "" , "NoiseMapFlipped" , Regavg_NoicemapArray , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG(ERROR) << "***Writing to Disk Fails***"  ;
        return (EXIT_FAILURE) ;
    }
    
    
    for (int i=0;i<xsize*ysize;i++){
        if (expArray[i] > 0.0 && image_snr[i] != INVALID_PIX_VALUE){        
        image_snr[i]=image_snr[i]/expArray[i];
        }
        else{
            image_snr[i]=0.0f;
        }
    }
    
    
     status = writeOutputImageToDisk ("fi" , moduleoutdir_snr , "" , "img" , image_snr , nameprefix ,9999 ,1 , xsize , ysize) ; //this is for the SignalFrame output
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                continue;
            }
    
    
    LOG(INFO)<<"Now starting Full Frame Astrometry ";
//    if((star_detect_algo_flag==1 || star_detect_algo_flag==3 || star_detect_algo_flag==5) && strcmp (channel.c_str (),"NUV"))
//    {
//        star_detect_algo_flag=5;
//    }
//    else if((star_detect_algo_flag==1 || star_detect_algo_flag==3 || star_detect_algo_flag==5) && strcmp (channel.c_str (),"FUV")){
//        star_detect_algo_flag=3 ;
//    }
//    else if((star_detect_algo_flag==2 || star_detect_algo_flag==4 || star_detect_algo_flag==6 )&& strcmp (channel.c_str (),"NUV"))
//    {
//        star_detect_algo_flag=4;
//    }
//    else if((star_detect_algo_flag==2 || star_detect_algo_flag==4 || star_detect_algo_flag==6 )&& strcmp (channel.c_str (),"FUV"))
//    {
//        star_detect_algo_flag=6;
//    }
    
    //convert the x-y image to RA-DEC image
    
    
//     fitsfile *fatt;
//    fits_open_file(&fatt,dirobj.attfile.c_str (),READONLY,&status);
//    printError(status,"Error opening attitude file",(char*)dirobj.attfile.c_str ());
//    
//    fits_movabs_hdu(fatt,2,NULL,&status);                           //Attitude data is in 2nd HDU
//    printError(status,"Error in moving to 2nd HDU ",(char*)dirobj.attfile.c_str ());
//    
//    long nrows_att;
//    fits_get_num_rows(fatt,&nrows_att,&status);
//    printError(status,"Error reading number of rows from attitude file  ",(char*)dirobj.attfile.c_str ());
//        
//    long firstrow=1;
//    long firstelem=1;
//    long nelements=nrows_att;                  //for time
//    double *time_att = new double[nelements];
//    double *RA_arr= new double[nelements];
//    double *DEC_arr= new double[nelements];
//    int tcoln,qcoln;               //column number for time and quaternion
//    double RA_pnt,DEC_pnt;
////    fits_get_colnum(fatt,CASEINSEN,timecol,&tcoln,&status);
////     printError(status,"Error reading time column number from attitude file  ",(char*)dirobj.attfile.c_str ());   
//    
//     fits_read_col(fatt,TDOUBLE,1, firstrow,firstelem, nelements, NULL,time_att,NULL,&status); 
//     printError(status,"Error reading time column from attitude file  ",(char*)dirobj.attfile.c_str ());     
//      fits_read_col(fatt,TDOUBLE,3, firstrow,firstelem, nelements, NULL,RA_arr,NULL,&status); 
//     printError(status,"Error reading RA column from attitude file  ",(char*)dirobj.attfile.c_str ());   
//     fits_read_col(fatt,TDOUBLE,4, firstrow,firstelem, nelements, NULL,DEC_arr,NULL,&status); 
//     printError(status,"Error reading DEC column from attitude file  ",(char*)dirobj.attfile.c_str ());     
//    long startindex=0, lastindex=nrows_att-1;
//    LOG(INFO)<<time_att[0]<<" "<<time_att[nrows_att-1]<<" "<<tstart<<" "<<tstop;
//    if(tstart>time_att[nrows_att-1] || tstop<time_att[0])
//    {
//        LOG(ERROR)<<"\033[1;31m***Time information in data file and attitude file do not match***\033[0m";
//        return (EXIT_FAILURE);
//    }
////    double time_mid =( tstart+tstop)/2;
//    double time_mid =tstart;
//   long index_attFile=0;
//    for (int i=1;i<nelements;i++)
//    {
//        if(time_mid>=time_att[i-1] && time_mid<time_att[i])
//        {
//            index_attFile=i-1;
//            RA_pnt=RA_arr[i-1];
//            DEC_pnt=DEC_arr[i-1];
//            break;
//        }    
//        
//    }
//   
//    double rotAng=0.0f;
//    
//    fits_read_col(fatt,TDOUBLE,5,index_attFile,1, 1, NULL,&rotAng,NULL,&status); 
//    printError(status,"Error reading rotation angle  from attitude file  ",(char*)dirobj.attfile.c_str ());
//    if(strcmp(datainfo.getDetector(),"NUV")==0)
//    {
//         rotAng=(rotAng+ANGLE_NUV);
//    }
//    else if(strcmp(datainfo.getDetector(),"FUV")==0)
//    {
//        LOG(INFO)<<"INSIDE";
//        LOG(INFO)<<rotAng;
//        rotAng=180-(rotAng+ANGLE_FUV);
//        LOG(INFO)<<rotAng;
//        //rotAng=(rotAng+ANGLE_FUV);
//    }
//    else 
//    {
//        LOG(ERROR)<<"***Invalid channel***";
//        return(EXIT_FAILURE);
//    }
//   
//   
//    fits_close_file(fatt,&status);
//    printError(status,"Error in closing the file  ",(char*)dirobj.attfile.c_str ());
//
//    LOG(INFO)<<"Roll Angle "<<rotAng; 
//    roll_angle_applied=rotAng;
//    
//    if(strcmp(datainfo.getDetector(),"NUV")==0){
//        
//        rotAng=-rotAng;
//    }
    //now update roll angle information  in shift and rotation file.
    long creffected_rows=del_CReffctedrows.size();
        long creffected_frames=frameno_crFailed.size();
    fitsfile *fsnr;
    //LOG(INFO)<<"Entered";
      fits_open_file (&fsnr , snrimageFilename.c_str() , READWRITE , &status) ;
     printError (status , "Error in opening the input information file" , infofile_in) ;
    
     fits_update_key (fsnr , TDOUBLE, "ROLLAPPLIED" , &rotAng, NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ;
      fits_update_key (fsnr , TLONG, "Cosmic-Ray Effected Events " , &creffected_rows , NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ;
     fits_update_key (fsnr , TLONG, "Cosmic-Ray_Effected_Frames " , &creffected_frames , NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ;
     fits_update_key (fsnr , TINT, "bad pixel-Multiple photon Effected Events " , &TotalNum_Effected_Rows , NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ; 
     fits_update_key (fsnr , TDOUBLE, "Average events per frame  " , &avg_Val_Events, NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ; 
     fits_update_key (fsnr , TDOUBLE, "Cosmic Ray calculated threshold value " , &curr_frmthr , NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ; 
     long TotalFrameCount=0;
     for (int i=0;i<TotFrameAvai_track.size();i++)
    {
        TotalFrameCount=TotalFrameCount+TotFrameAvai_track[i];
    }
    fits_update_key (fsnr , TLONG, "Total number of remaining frames" , &TotalFrameCount , NULL , &status) ;
    printError (status , "***Error in updating the nameprefix keyword***") ; 
    
     
     fits_close_file(fsnr,&status);
    
    //till this
    
    
    int mid_X_new,mid_Y_new;
    ctheta=cos((rotAng)*M_PI/180);
    stheta=sin((rotAng)*M_PI/180);
//    if(strcmp(datainfo.getDetector(),"NUV")==0){
//    ctheta=cos((-rotAng)*M_PI/180);
//    stheta=sin((-rotAng)*M_PI/180);
//    }
//    else if (strcmp(datainfo.getDetector(),"FUV")==0){
//        ctheta=cos((rotAng)*M_PI/180);
//    stheta=sin((rotAng)*M_PI/180);
//    }
   
    //LOG(INFO)<<"outside1";
  //  float *Rotated_image = new float[FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
   float *Rotated_image = new float[FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2];
   float *subDividedArray=new float[FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2];
   initArray(subDividedArray,FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2,0.0f);
   
   
   
   float *Rotated_Noicemap= new float[FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2];
   float *subDividedNoice= new float[FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2];
   //SubDividing the image from 4800 to 9600 for bi-linear correction. 
   for (int i=0;i<FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG;i++){
       if(Regavg_subSampled_Array_Sig[i]!=INVALID_PIX_VALUE)
       Regavg_subSampled_Array_Sig[i]=Regavg_subSampled_Array_Sig[i]/4;
   }
    for (int i=0;i<FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG;i++){
        Regavg_NoicemapArray[i]=Regavg_NoicemapArray[i]/4;
    }
   
   performSubDivisionIM(Regavg_NoicemapArray,FINALFRAMESIZE_REGAVG,FINALFRAMESIZE_REGAVG,subDividedNoice,FINALFRAMESIZE_REGAVG*2,FINALFRAMESIZE_REGAVG*2);
    performSubDivisionIM(Regavg_subSampled_Array_Sig,FINALFRAMESIZE_REGAVG,FINALFRAMESIZE_REGAVG,subDividedArray,FINALFRAMESIZE_REGAVG*2,FINALFRAMESIZE_REGAVG*2);
  //  LOG(INFO)<<"outside2";
   //  for(int i=0;i<FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG;i++){
    for(int i=0;i<FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2;i++){
          Rotated_image[i]=0.0f;
      }
    long indexes=0;
     int n1, n2, k1, k2;
     bool cnt_FullFrameAst=0;
    // while(cnt_FullFrameAst<2){
         cnt_FullFrameAst++;
         
    
//      for (int i =0;i<FINALFRAMESIZE_REGAVG*2;i++)
//    {
//        mid_X_new=i-FINALFRAMESIZE_REGAVG;
//          for (int j=0;j<FINALFRAMESIZE_REGAVG*2;j++)
//          {  
//              mid_Y_new=j-FINALFRAMESIZE_REGAVG;
//             
//                x1=(-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta)))+FINALFRAMESIZE_REGAVG;
//              x2=(-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) +FINALFRAMESIZE_REGAVG;
//              
////              x1=(-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta)))+FINALFRAMESIZE_REGAVG;
////              x2=(-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) +FINALFRAMESIZE_REGAVG;  
////              if(x1>0 && x1<FINALFRAMESIZE_REGAVG && x2>0 && x2<FINALFRAMESIZE_REGAVG){
////                  imageLocation_Array_X[(int)(round(x2)*FINALFRAMESIZE_REGAVG+round(x1))]=Regavg_subSampled_Array_Sig[j*FINALFRAMESIZE_REGAVG+i];
////              }
//              
//               k1 = (int)(x1); x1 -= k1;
//	       k2 = (int)(x2); x2 -= k2;
//            
//             if (k1 >=0 && k1 < FINALFRAMESIZE_REGAVG*2-1 && k2 >=0 && k2 < FINALFRAMESIZE_REGAVG*2-1) {
//			//Rotated_image[(FINALFRAMESIZE_REGAVG*2-i)*FINALFRAMESIZE_REGAVG*2+(FINALFRAMESIZE_REGAVG*2-j)] = 
//                 //Rotated_image[x2*(FINALFRAMESIZE_REGAVG*2)+x1] = 
//                  Rotated_image[j*(FINALFRAMESIZE_REGAVG*2)+i] = 
////			    (1.-x1)*(1.-x2)*subDividedArray[k1*(FINALFRAMESIZE_REGAVG*2)+k2]   +
////			    x1     *(1.-x2)*subDividedArray[(k1+1)*(FINALFRAMESIZE_REGAVG*2)+k2] +
////			    (1.-x1)*x2     *subDividedArray[k1*(FINALFRAMESIZE_REGAVG*2)+(k2+1)] +
////			    x1     *x2     *subDividedArray[(k1+1)*(FINALFRAMESIZE_REGAVG*2)+(k2+1)];
//                            (1.-x2)*(1.-x1)*subDividedArray[k2*(FINALFRAMESIZE_REGAVG*2)+k1]   +
//			    x2     *(1.-x1)*subDividedArray[(k2)*(FINALFRAMESIZE_REGAVG*2)+(k1+1)] +
//			    (1.-x2)*x1     *subDividedArray[(k2+1)*(FINALFRAMESIZE_REGAVG*2)+(k1)] +
//			    x2     *x1     *subDividedArray[(k2+1)*(FINALFRAMESIZE_REGAVG*2)+(k1+1)];
//             }
//             else {
//                 Rotated_image[j*(FINALFRAMESIZE_REGAVG*2)+i]=0.0f;
//             }
//          }
//
//        
//        
//     }
         
         for (int i = 0; i < FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG; i++) {
            if (expArray[i] != INVALID_PIX_VALUE)
                expArray[i] = expArray[i] / 4;
            
        }
         if(strcmp(datainfo.getDetector(),"NUV")==0){
         float *flippedExposureArray = new float[FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
          for (int i=0;i<FINALFRAMESIZE_REGAVG;i++)
        {
            for (int j=0;j<FINALFRAMESIZE_REGAVG;j++)
            {
            flippedExposureArray[i*FINALFRAMESIZE_REGAVG+j]=expArray[i*FINALFRAMESIZE_REGAVG+((FINALFRAMESIZE_REGAVG-1)-j)];
           
            }
            
        }

            for (int i = 0; i < FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG; i++) {
               
                expArray[i] = flippedExposureArray[i];

            }

            for (int i = 0; i < nrows - track_rown.size(); i++) {
                temp_X_newsnr[i] = FINALFRAMESIZE_REGAVG - temp_X_newsnr[i];

            }
            
             delete[] flippedExposureArray;
             
         }
         float *Binned_ExpArray = new float[FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
         initArray(Binned_ExpArray,FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG,0.0f);
         float *Rotated_ExpArray= new float[FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2];
         initArray(Rotated_ExpArray,FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2,0.0f);
         
         float *counterArray = new float[FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2];
         initArray(counterArray,FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2,0.0f);
         
         
double locX,locY;
         if(fi_flag==FALSE){
         float *Exparray_Expanded= new float[FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2];
         initArray(Exparray_Expanded,FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2,0.0f);
         
         
         performSubDivisionIM(expArray,FINALFRAMESIZE_REGAVG,FINALFRAMESIZE_REGAVG,Exparray_Expanded,FINALFRAMESIZE_REGAVG*2,FINALFRAMESIZE_REGAVG*2);

            for (int i = 0; i < FINALFRAMESIZE_REGAVG * 2; i++) {
                mid_X_new = i - FINALFRAMESIZE_REGAVG;
                for (int j = 0; j < FINALFRAMESIZE_REGAVG * 2; j++) {
                    mid_Y_new = j - FINALFRAMESIZE_REGAVG;
                    x1 = (-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + FINALFRAMESIZE_REGAVG;
                    x2 = (-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + FINALFRAMESIZE_REGAVG;

                    //                  if((int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))> 0 && (int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1)) <FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2)
                    if (round(x1) > 0 && round(x1) < FINALFRAMESIZE_REGAVG * 2 && round(x2) > 0 && round(x2) < FINALFRAMESIZE_REGAVG * 2) {
                        if (Exparray_Expanded[j * FINALFRAMESIZE_REGAVG * 2 + i] != INVALID_PIX_VALUE && Rotated_ExpArray[(int) (round(x2) * FINALFRAMESIZE_REGAVG * 2 + round(x1))]!=INVALID_PIX_VALUE ) {
                            Rotated_ExpArray[(int) (round(x2) * FINALFRAMESIZE_REGAVG * 2 + round(x1))] = Rotated_ExpArray[(int) (round(x2) * FINALFRAMESIZE_REGAVG * 2 + round(x1))] + Exparray_Expanded[j * FINALFRAMESIZE_REGAVG * 2 + i];
                            counterArray[(int) (round(x2) * FINALFRAMESIZE_REGAVG * 2 + round(x1))]=counterArray[(int) (round(x2) * FINALFRAMESIZE_REGAVG * 2 + round(x1))]+1;
                        } else {
                            Rotated_ExpArray[(int) (round(x2) * FINALFRAMESIZE_REGAVG * 2 + round(x1))] = INVALID_PIX_VALUE;
                        }

                    }
                }
            }
         
         //Applying binning 
         ApplyBinning(Rotated_ExpArray,FINALFRAMESIZE_REGAVG*2,FINALFRAMESIZE_REGAVG*2,Binned_ExpArray,FINALFRAMESIZE_REGAVG,FINALFRAMESIZE_REGAVG,counterArray);
        
       //  initArray(counterArray,FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG,0.0f);
         
//Exposure Array manipulation.


//float temporary_binnedExpArray = new float[FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
for (int i=0;i<FINALFRAMESIZE_REGAVG-8;i=i+8){
for (int j=0;j<FINALFRAMESIZE_REGAVG-8;j=j+8){
 sum_OfPixelsOfExp=0.0f;
for (int k=i;k<i+8;k++){
for (int l=j;l<j+8;l++){
sum_OfPixelsOfExp=sum_OfPixelsOfExp+Binned_ExpArray[l*FINALFRAMESIZE_REGAVG+k];

}
}
//compute average now
sum_OfPixelsOfExp=sum_OfPixelsOfExp/64;
 
for (int k=i;k<i+8;k++){
for (int l=j;l<j+8;l++){
Binned_ExpArray[l*FINALFRAMESIZE_REGAVG+k]=sum_OfPixelsOfExp;

}
}



}
}




//Exposure Array manipulation ends


         for(int i =0;i<nrows-track_rown.size();i++)
         {
             locX=temp_X_newsnr[i]*2-FINALFRAMESIZE_REGAVG;
             locY=temp_Y_newsnr[i]*2-FINALFRAMESIZE_REGAVG;
             
              x1=(-(locX * (ctheta)) - (locY * (stheta)))+FINALFRAMESIZE_REGAVG;
              x2=(-(locX * (stheta)) +(locY * (ctheta))) +FINALFRAMESIZE_REGAVG;
              //if((round(x2)*(FINALFRAMESIZE_REGAVG*2)+round(x1))>0 &&  (round(x2)*(FINALFRAMESIZE_REGAVG*2)+round(x1))<FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2 )
            if(round(x1)>0 && round(x1)<FINALFRAMESIZE_REGAVG*2   && round(x2)>0 && round(x2)<FINALFRAMESIZE_REGAVG*2)
              Rotated_image[(int)(round(x2)*(FINALFRAMESIZE_REGAVG*2)+round(x1))] =Rotated_image[(int)(round(x2)*(FINALFRAMESIZE_REGAVG*2)+round(x1))]+temp_enp__newsnr[i]*temp_mult_newsnr[i]*temp_bflag_newsnr[i];
           //   counterArray[(int)(round(x2)*(FINALFRAMESIZE_REGAVG*2)+round(x1))]=counterArray[(int)(round(x2)*(FINALFRAMESIZE_REGAVG*2)+round(x1))]+1;
         }
         delete[] Exparray_Expanded;//,Rotated_ExpArray;
         }
         else {
            
            for (int i = 0; i < FINALFRAMESIZE_REGAVG * 2; i++) {
                mid_X_new = i - FINALFRAMESIZE_REGAVG;
                for (int j = 0; j < FINALFRAMESIZE_REGAVG * 2; j++) {
                    mid_Y_new = j - FINALFRAMESIZE_REGAVG;

                    x1 = (-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + FINALFRAMESIZE_REGAVG;
                    x2 = (-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + FINALFRAMESIZE_REGAVG;

                   
                    //if ((int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1)) >= 0 && (int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1)) <FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2) {
                      if(round(x1)>0 && round(x1)<FINALFRAMESIZE_REGAVG*2   && round(x2)>0 && round(x2)<FINALFRAMESIZE_REGAVG*2) {
                          if(subDividedArray[j * (FINALFRAMESIZE_REGAVG * 2) + i]!=INVALID_PIX_VALUE && Rotated_image[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))]!=INVALID_PIX_VALUE ){
                        Rotated_image[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))] =Rotated_image[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))]+subDividedArray[j * (FINALFRAMESIZE_REGAVG * 2) + i];
                        counterArray[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))]= counterArray[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))]+1; 
                          }
                          else{
                          Rotated_image[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))] =INVALID_PIX_VALUE;    
                      }
                                
                    } 
                }



            }

        }

//rotating noicemap 
for (int i = 0; i < FINALFRAMESIZE_REGAVG * 2; i++) {
                mid_X_new = i - FINALFRAMESIZE_REGAVG;
                for (int j = 0; j < FINALFRAMESIZE_REGAVG * 2; j++) {
                    mid_Y_new = j - FINALFRAMESIZE_REGAVG;

                    x1 = (-(mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + FINALFRAMESIZE_REGAVG;
                    x2 = (-(mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + FINALFRAMESIZE_REGAVG;

                   
                    //if ((int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1)) >= 0 && (int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1)) <FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2) {
                      if(round(x1)>0 && round(x1)<FINALFRAMESIZE_REGAVG*2   && round(x2)>0 && round(x2)<FINALFRAMESIZE_REGAVG*2) {
                          if(subDividedNoice[j * (FINALFRAMESIZE_REGAVG * 2) + i]!=INVALID_PIX_VALUE && Rotated_Noicemap[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))]!=INVALID_PIX_VALUE ){
                        Rotated_Noicemap[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))] =Rotated_Noicemap[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))]+subDividedNoice[j * (FINALFRAMESIZE_REGAVG * 2) + i];
                        counterArray[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))]= counterArray[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))]+1; 
                          }
                          else{
                          Rotated_Noicemap[(int)(round(x2)*FINALFRAMESIZE_REGAVG*2+round(x1))] =INVALID_PIX_VALUE;    
                      }
                                
                    } 
                }



            }









         
            LOG(INFO)<<"Writing RA_DEC image to the disk";
    int  tfields_3 = 14 ;
    char *ttype_3[] = {"PacketSequence " , "FrameCount" , "Time" , "Fx" , "Fy" , "Max-Min" , "Min" ,"UVIT_MASTER_TIME","BAD FLAG" , "MULT PHOTON" , "EFFECTIVE_NUM_PHOTONS","RA","DEC","MJD_L2"} ;
    char *tform_3[] = {"U" , "U" , "D" , "1D" , "1D" , "B" , "B" ,"D", "B" , "B" , "1D","1D","1D","D"} ;
    char *tunit_3[] = {"" , "" , "" , "" , "" , "" , "" , "" , "" , "","","","",""} ;
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_snr , nameprefix , "radec") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
       
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields_3 , ttype_3 , tform_3 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
          fits_update_key (fout , TSTRING, "FILTER" , datainfo.getFilter(), NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ;
      fits_update_key (fout , TINT, "WIN_X_SZ" , &win_xsize, NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ;
        
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;
        for(int i =0;i<nrows-track_rown.size();i++)
         {
             locX=temp_X_newsnr[i]-FINALFRAMESIZE_REGAVG/2;
             locY=temp_Y_newsnr[i]-FINALFRAMESIZE_REGAVG/2;
             if(((-(locX * (ctheta)) - (locY * (stheta)))+(FINALFRAMESIZE_REGAVG/2)) >0 && ((-(locX * (ctheta)) - (locY * (stheta)))+(FINALFRAMESIZE_REGAVG/2))<4800 &&
                    ((-(locX * (stheta)) +(locY * (ctheta))) +(FINALFRAMESIZE_REGAVG/2)) > 0 && ((-(locX * (stheta)) +(locY * (ctheta))) +(FINALFRAMESIZE_REGAVG/2) )<4800 ){
              temp_ra_newsnr[i]=(-(locX * (ctheta)) - (locY * (stheta)))+(FINALFRAMESIZE_REGAVG/2);
              temp_dec_newsnr[i]=(-(locX * (stheta)) +(locY * (ctheta))) +(FINALFRAMESIZE_REGAVG/2);
             }
             else{
                 temp_ra_newsnr[i]=INVALID_PIX_VALUE;
                 temp_dec_newsnr[i]=INVALID_PIX_VALUE;
             }
              
            
         }
        double * Mjd_l2 = new double[nrows-track_rown.size()];
	for(int i =0;i<nrows-track_rown.size();i++)
		{
		  Mjd_l2[i]=bscale_mjd*temp_uvittime_newsnr[i]+bzero_mjd;
		}
/* .......................................... 27-Oct-2022 (trying FIX)
 + replacing variable-names to post-trimming operation in ShiftNRot
*/    


        status = writeColumnsToFITS (outfile , 2 , 14 , TUSHORT , 1 , psc , nrows-track_rown.size() , TUSHORT , 2 , temp_frameno_newsnr, nrows-track_rown.size() , TDOUBLE , 3 , temp_ftime_newsnr, nrows-track_rown.size() , TFLOAT , 4 , temp_X_newsnr , nrows-track_rown.size() , TFLOAT , 5 , temp_Y_newsnr ,
                nrows-track_rown.size() , TBYTE , 6 , dmm_newsnr , nrows-track_rown.size() , TBYTE , 7 , min_newsnr , nrows-track_rown.size() ,TDOUBLE,8,temp_uvittime_newsnr ,nrows-track_rown.size(), TUSHORT , 9 , temp_bflag_newsnr, nrows-track_rown.size() , TUSHORT , 10 , 
                temp_mult_newsnr, nrows-track_rown.size() , TFLOAT , 11 , temp_enp__newsnr, nrows-track_rown.size(),TFLOAT,12,temp_ra_newsnr,nrows-track_rown.size(),TFLOAT,13,temp_dec_newsnr,nrows-track_rown.size(),TDOUBLE,14,Mjd_l2,nrows-track_rown.size()) ;
//
//        status = writeColumnsToFITS (outfile , 2 , 14 , TUSHORT , 1 , psc , nrows-track_rown.size() , TUSHORT , 2 , frame_no , nrows-track_rown.size() , TDOUBLE , 3 , time_frame , nrows-track_rown.size() , TFLOAT , 4 , temp_X_newsnr , nrows-track_rown.size() , TFLOAT , 5 , temp_Y_newsnr ,
//                nrows-track_rown.size() , TBYTE , 6 , dmm_newsnr , nrows-track_rown.size() , TBYTE , 7 , min_newsnr , nrows-track_rown.size() ,TDOUBLE,8,temp_uvittime_newsnr ,nrows-track_rown.size(), TUSHORT , 9 , badFlag_temp , nrows-track_rown.size() , TUSHORT , 10 , 
//                mult_temp , nrows-track_rown.size() , TFLOAT , 11 , effective_NumPhotons , nrows-track_rown.size(),TFLOAT,12,temp_ra_newsnr,nrows-track_rown.size(),TFLOAT,13,temp_dec_newsnr,nrows-track_rown.size(),TDOUBLE,14,Mjd_l2,nrows-track_rown.size()) ;
// .............................................................................
        if (status)
        {
            LOG (INFO) << "Error in writing to the shift and rotation RADEC eventFile to the disk" << endl ;
             continue;
        }
         
           //free(subDividedArray);
           // free(Rotated_image);
   //  }
//     else {
////
////         ctheta=cos((360-rotAng)*M_PI/180);
////         stheta=sin((360-rotAng)*M_PI/180);
//            ctheta=cos((rotAng)*M_PI/180);
//         stheta=sin((rotAng)*M_PI/180);
//            for (int i = 0; i < nrows - track_rown.size(); i++) 
//            {
//                if(strcmp(datainfo.getDetector(),"NUV")==0){
//              x1=((-(((FINALFRAMESIZE_REGAVG-1)-temp_X_newsnr[i])-FINALFRAMESIZE_REGAVG/2) * (ctheta)) - ((temp_Y_newsnr[i]-FINALFRAMESIZE_REGAVG/2) * (stheta)))+FINALFRAMESIZE_REGAVG/2;
//              x2=((-(((FINALFRAMESIZE_REGAVG-1)-temp_X_newsnr[i])-FINALFRAMESIZE_REGAVG/2) * (stheta)) +((temp_Y_newsnr[i]-FINALFRAMESIZE_REGAVG/2)* (ctheta))) +FINALFRAMESIZE_REGAVG/2;  
//                }
//                else {
//              x1=((-((temp_X_newsnr[i])-FINALFRAMESIZE_REGAVG/2) * (ctheta)) - ((temp_Y_newsnr[i]-FINALFRAMESIZE_REGAVG/2) * (stheta)))+FINALFRAMESIZE_REGAVG/2;
//              x2=((-((temp_X_newsnr[i])-FINALFRAMESIZE_REGAVG/2) * (stheta)) +((temp_Y_newsnr[i]-FINALFRAMESIZE_REGAVG/2)* (ctheta))) +FINALFRAMESIZE_REGAVG/2;  
//                }
//              if((int)((int)(round(x2))*FINALFRAMESIZE_REGAVG+(int)(round(x1)))<FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG){ 
//              //Rotated_image[(int)(FINALFRAMESIZE_REGAVG-round(x2))*FINALFRAMESIZE_REGAVG+(int)(round(x1))]=Rotated_image[(int)(FINALFRAMESIZE_REGAVG-round(x2))*FINALFRAMESIZE_REGAVG+(int)(round(x1))]+effective_NumPhotons[i]*mult_temp[i]*badFlag_temp[i];
//                 //Rotated_image[(int)(round(x1))*FINALFRAMESIZE_REGAVG+(int)(FINALFRAMESIZE_REGAVG-round(x2))]=Rotated_image[(int)(round(x1))*FINALFRAMESIZE_REGAVG+(int)(FINALFRAMESIZE_REGAVG-round(x2))]+effective_NumPhotons[i]*mult_temp[i]*badFlag_temp[i];
//             Rotated_image[(int)(round(x2))*FINALFRAMESIZE_REGAVG+(int)(round(x1))]=Rotated_image[(int)(round(x2))*FINALFRAMESIZE_REGAVG+(int)(round(x1))]+effective_NumPhotons[i]*mult_temp[i]*badFlag_temp[i];
//              }
//            }
//     }
         
         //added 
               
      float *subSampledArray = new float [FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
      initArray(subSampledArray,FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG,0.0f);
      if(fi_flag==FALSE){
          ApplySubSampling_Addition(Rotated_image,FINALFRAMESIZE_REGAVG*2,FINALFRAMESIZE_REGAVG*2,subSampledArray,FINALFRAMESIZE_REGAVG,FINALFRAMESIZE_REGAVG);
      }
      else{
        ApplyBinning(Rotated_image,FINALFRAMESIZE_REGAVG*2,FINALFRAMESIZE_REGAVG*2,subSampledArray,FINALFRAMESIZE_REGAVG,FINALFRAMESIZE_REGAVG,counterArray);
     }
       float *subSampledNoiceArray = new float [FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
      ApplySubSampling_Addition(Rotated_Noicemap,FINALFRAMESIZE_REGAVG*2,FINALFRAMESIZE_REGAVG*2,subSampledNoiceArray,FINALFRAMESIZE_REGAVG,FINALFRAMESIZE_REGAVG);
        if(fi_flag==FALSE){
  // float *temparrforreg = new float[FINALFRAMESIZE_REGAVG*2*FINALFRAMESIZE_REGAVG*2];
        for(int i =0;i<FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG;i++){
	//temparrforreg[i]=subSampledArray[i];
            if(Binned_ExpArray[i]!=0.0 && Binned_ExpArray[i]!=INVALID_PIX_VALUE && subSampledArray[i]!=INVALID_PIX_VALUE){
            subSampledArray[i]=subSampledArray[i]/Binned_ExpArray[i];
            }
            else{
                subSampledArray[i]=0.0f;
            }
        }
        }
       
//          for(int i=0;i<FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG;i++){
//              subSampledArray[i]=Rotated_image[i];
//          }
     
   //  writeHistory((char*)snrimageFilename.c_str(),header_info);
   
     
     // delete[] subSampledArray;
    vector<string> HeaderStrings;
      fits_open_file (&fsnr , snrimageFilename.c_str() , READWRITE , &status) ;
     printError (status , "Error in opening the input information file" , (char*)snrimageFilename.c_str()) ;  
     copyUsrkeywrdsTovect(fsnr,HeaderStrings);
     //fits_movabs_hdu (fsnr , 1 , NULL , &status) ;
           //long naxes1[2] ;
    naxes1[0] = naxes1[1] = 0 ;
     naxis = 2 ;
    bitpix = FLOAT_IMG ;
        fits_create_img (fsnr , bitpix , naxis , naxes1 , &status) ;
    printError (status , "Error in Creating the image for Signal Fie" , infofile_in);
     for(int i=0;i<L1keywords.size ();i++)
         {
            if(strstr(L1keywords[i].c_str () ,"NAXIS")==NULL)
                fits_write_record (fsnr , (char*)L1keywords[i].c_str () , &status) ;
         }
     
     
     
     fits_close_file(fsnr,&status);
    printError (status , "Error in closing the file" , (char*)snrimageFilename.c_str()) ;         
        
    

        
    
    

    fits_open_file (&fsnr , regAvgFilename.c_str() , READWRITE , &status) ;
    printError (status , "Error in opening the input information file" , (char*)regAvgFilename.c_str()) ;  
     //copyUsrkeywrdsTovect(fsnr,HeaderStrings);
     //fits_movabs_hdu (fsnr , 1 , NULL , &status) ;
           //long naxes1[2] ;
    naxes1[0] = naxes1[1] = 0 ;
     naxis = 2 ;
    bitpix = FLOAT_IMG ;
        fits_create_img (fsnr , bitpix , naxis , naxes1 , &status) ;
    printError (status , "Error in Creating the image for Signal Fie" ,(char*) regAvgFilename.c_str());
     for(int i=0;i<L1keywords.size ();i++)
         {
            if(strstr(L1keywords[i].c_str () ,"NAXIS")==NULL)
                fits_write_record (fsnr ,(char*)L1keywords[i].c_str () , &status) ;
         }
     
    fits_close_file(fsnr,&status);
    printError (status , "Error in closing the file" , (char*)regAvgFilename.c_str()) ;        
    
        
    status = setDirectoryStructure (moduleoutdir_radecimg , "") ;//setting output directory structure for registration and averaging 
    if (status)
    {
        LOG(ERROR) << "***Directory Structure has not been successfully set-up***"  ;
        return (EXIT_FAILURE) ;
    };
//if(flag_Roll_Angleinvalid==FALSE){
  status = writeOutputImageToDisk ("ra-dec" , moduleoutdir_radecimg , "" , "sig" , subSampledArray , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
//  status = writeOutputImageToDisk ("ra-dec" , moduleoutdir_radecimg , "" , "expfull" , Rotated_ExpArray , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG*2 , FINALFRAMESIZE_REGAVG*2) ;
  //status = writeOutputImageToDisk ("ra-dec" , moduleoutdir_radecimg , "" , "sigfull" , temparrforreg , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG*2 , FINALFRAMESIZE_REGAVG*2) ;
 //status = writeOutputImageToDisk ("ra-dec" , moduleoutdir_radecimg , "" , "sig" , Rotated_image , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG*2 , FINALFRAMESIZE_REGAVG*2) ;
//  delete[] temparrforreg;
  if (status)
    {
        LOG(ERROR) << "***Writing to Disk Fails***"  ;
        return (EXIT_FAILURE) ;
    }
   status = writeOutputImageToDisk ("ra-decexp" , moduleoutdir_radecimg , "" , "sig" ,Binned_ExpArray  , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG(ERROR) << "***Writing to Disk Fails***"  ;
        return (EXIT_FAILURE) ;
    }
   
   status = writeOutputImageToDisk ("ra-decnoise" , moduleoutdir_radecimg , "" , "noise" ,subSampledNoiceArray  , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG(ERROR) << "***Writing to Disk Fails***"  ;
        return (EXIT_FAILURE) ;
    }
 
//}
//else {
  //status = writeOutputImageToDisk ("NOT-CORRECTsigCHECK_X_Y_IMAGEONLY_ra-dec" , moduleoutdir_radecimg , "" , "sig" , subSampledArray , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;

  //if (status)
    //{
      //  LOG(ERROR) << "***Writing to Disk Fails***"  ;
        //return (EXIT_FAILURE) ;
   // }
   //status = writeOutputImageToDisk ("NOT-CORRECTexp-CHECK_X_Y_IMAGEONLY_ra-decexp" , moduleoutdir_radecimg , "" , "sig" ,Binned_ExpArray  , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
  //  if (status)
    //{
      //  LOG(ERROR) << "***Writing to Disk Fails***"  ;
        //return (EXIT_FAILURE) ;
 //   }




//}

    fits_open_file (&fsnr , radecFilename.c_str() , READWRITE , &status) ;
    printError (status , "Error in opening the input information file" , (char*)radecFilename.c_str()) ;  
     //copyUsrkeywrdsTovect(fsnr,HeaderStrings);
     //fits_movabs_hdu (fsnr , 1 , NULL , &status) ;
           //long naxes1[2] ;
    naxes1[0] = naxes1[1] = 0 ;
     naxis = 2 ;
    bitpix = FLOAT_IMG ;
        fits_create_img (fsnr , bitpix , naxis , naxes1 , &status) ;
    printError (status , "Error in Creating the image for Signal Fie" ,(char*) radecFilename.c_str());
     for(int i=0;i<L1keywords.size ();i++)
         {
            if(strstr(L1keywords[i].c_str () ,"NAXIS")==NULL)
                fits_write_record (fsnr ,(char*)L1keywords[i].c_str () , &status) ;
         }
     
    fits_close_file(fsnr,&status);
    printError (status , "Error in closing the file" ,(char*)radecFilename.c_str()) ;  
//  LOG(INFO)<<"fdvdvdvdv";

     delete[] Binned_ExpArray;
     delete[] Rotated_ExpArray;
     delete[] Rotated_Noicemap;
   
delete[] subSampledArray;
delete[] counterArray;
delete[] Rotated_image;
delete[] subDividedArray;
delete[] subDividedNoice;
delete[]Regavg_subSampled_Array_Sig,Regavg_subSampled_Array_Exp,Regavg_NoicemapArray;
    //writeUsrkeywordsFrmvect((char*)radecFilename.c_str(),HeaderStrings);  
    
      //    LOG(INFO)<<"fdvdvdvdv";exit(1);
        
     fits_open_file (&finfo_in , infofile_in , READWRITE , &status) ;
     printError (status , "Error in opening the input information file" , infofile_in) ;
     fits_update_key (finfo_in , TLONG, "Cosmic-Ray Effected Events " , &creffected_rows , NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ;
     fits_update_key (finfo_in , TLONG, "Cosmic-Ray_Effected_Frames " , &creffected_frames , NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ;
     fits_update_key (finfo_in , TINT, "bad pixel-Multiple photon Effected Events " , &TotalNum_Effected_Rows , NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ; 
     fits_update_key (finfo_in , TDOUBLE, "Average events per frame  " , &avg_Val_Events, NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ; 
     fits_update_key (finfo_in , TDOUBLE, "Cosmic Ray calculated threshold value " , &curr_frmthr , NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ; 
//     long TotalFrameCount=0;
//     forrm co    
//     for (int i = 0; i < size; i++) {
//            Object elem = array[i];
//
//        }
//(int i=0;i<TotFrameAvai_track.size();i++)
//    {
//        TotalFrameCount=TotalFrameCount+TotFrameAvai_track[i];
//    }
     
    fits_update_key (finfo_in , TLONG, "Total number of remaining frames" , &TotalFrameCount , NULL , &status) ;
    printError (status , "***Error in updating the nameprefix keyword***") ; 
   
    vector<string> tempHeaderInfo;
    copyUsrkeywrdsTovect(finfo_in,tempHeaderInfo);
    
    fits_close_file(finfo_in,&status);
    printError (status , "***Error in closing the file***") ; 
    double integration_Time;
    readKeywords (infofile_in, 2 , 1 , TDOUBLE , "INT_TIME" , &integration_Time) ;
    updateKeywords((char*)snrimageFilename.c_str(),1,1,TDOUBLE,"INT_TIME" , &integration_Time);
    writeUsrkeywordsFrmvect((char*)snrimageFilename.c_str(),tempHeaderInfo);
writeUsrkeywordsFrmvect((char*)radecFilename.c_str(),tempHeaderInfo);
    writeHistory((char*)radecFilename.c_str(),header_info);
    //performing full frame astrometry
LOG(INFO)<<"The refineWin->"<<refine_Winsize; 
float crpix1=FINALFRAMESIZE_REGAVG/2;
float crpix2=FINALFRAMESIZE_REGAVG/2;
float crota1=0.0f,cdelt1=0.0f,cdelt2=0.0f;
float crota2=0.0f;
 int factor_delta=FINALFRAMESIZE_REGAVG/600;
 if(strcmp(datainfo.getDetector (),"VIS")==0)
    {
           cdelt1=(3.357/3600)/factor_delta;
           cdelt2=(3.311/3600)/factor_delta;
    }
    else if (strcmp(datainfo.getDetector (),"FUV")==0){
        //cdelt1=3.357/3600;
       //  cdelt2=3.311/3600;
       cdelt1=cdelt2=(3.3373/3600)/factor_delta;
    }
    else if(strcmp(datainfo.getDetector (),"NUV")==0){
       cdelt1=cdelt2=(3.3307/3600)/factor_delta;
       // cdelt1=cdelt2=(3.3307*factor_delta/3600)/;
       
    }
   // cdelt1=-cdelt1/cos(center_dec*M_PI/180);
     cdelt1=-cdelt1;///cos(center_dec*M_PI/180);
     LOG(INFO)<<status;
     fits_open_file (&finfo_in , radecFilename.c_str() , READWRITE , &status) ;
     printError (status , "Error in opening the input information file" , infofile_in) ;
     fits_write_key(finfo_in, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_write_key(finfo_in, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_write_key(finfo_in, TDOUBLE, "CRVAL1", &RA_pnt, "", &status);                                              printError(status,"");  
    fits_write_key(finfo_in, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_write_key(finfo_in, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_write_key(finfo_in, TDOUBLE, "CRVAL2", &DEC_pnt, "", &status);                                    printError(status,"");   
    fits_write_key(finfo_in, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
    
    fits_close_file(finfo_in,&status);
    
      LOG(INFO)<<status;
     fits_open_file (&finfo_in , radecNoiceFilename.c_str() , READWRITE , &status) ;
     printError (status , "Error in opening the input information file" , infofile_in) ;
    fits_write_key(finfo_in, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_write_key(finfo_in, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_write_key(finfo_in, TDOUBLE, "CRVAL1", &RA_pnt, "", &status);                                              printError(status,"");  
    fits_write_key(finfo_in, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_write_key(finfo_in, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_write_key(finfo_in, TDOUBLE, "CRVAL2", &DEC_pnt, "", &status);                                    printError(status,"");   
    fits_write_key(finfo_in, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
    fits_close_file(finfo_in,&status);
    
    
      LOG(INFO)<<status;
    fits_open_file (&finfo_in , radecexpFilename.c_str() , READWRITE , &status) ;
     printError (status , "Error in opening the input information file" , infofile_in) ;
     fits_write_key(finfo_in, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_write_key(finfo_in, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_write_key(finfo_in, TDOUBLE, "CRVAL1", &RA_pnt, "", &status);                                              printError(status,"");  
    fits_write_key(finfo_in, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_write_key(finfo_in, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_write_key(finfo_in, TDOUBLE, "CRVAL2", &DEC_pnt, "", &status);                                    printError(status,"");   
    fits_write_key(finfo_in, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_write_key(finfo_in, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
    fits_close_file(finfo_in,&status);
    
       uvtFullFrameAst ast_obj ;
       ast_obj.read ((char*)moduleoutdir_snr ,(char*)radecFilename.c_str(),(char*)caldbindir.c_str (),(char*)filenamesnr.c_str(),(char*)radecNoiceFilename.c_str(),(char*)radecexpFilename.c_str(),(char*)outputdir ,(char*)dirobj.attfile.c_str (),(char*)att_timecol ,(char*)att_qcol,(char*)caldbindir.c_str (),sd_multi_factor_default,minimum_No_of_Stars,refine_Winsize/mult_fact,centroid_Winsize/mult_fact,databasename,search_algo_ctlg,len_a,len_b,rad_search,clobber , history,0) ;
       ast_obj.display () ;
     
       
       status = ast_obj.uvtFullFrmAstProcess () ;
       if (status)
       {
            LOG(ERROR) << endl << "Error in full frame astrometry  module" ;
LOG(ERROR)<<"CRASH: FULL FRAME ASTROMETRY FAILED (POSSIBLY DUE TO SIGMA MULTIPLIER BEING =< 0) (uvtLevel2PC.cpp)";
            continue;
       }
     
       
     //  float diffAddra=0.0f;
      // float diffAdddec=0.0f;
      float  diffAddra=ast_obj.getDiffRAval();
      float  diffAdddec=ast_obj.getDiffDECval();
       
       double centra=ast_obj.getRAVAL();
       double centdec=ast_obj.getDECVAL();
         LOG(INFO)<<status;
       fits_open_file (&finfo_in , radecFilename.c_str() , READWRITE , &status) ;
     printError (status , "Error in opening the input information file" , infofile_in) ;
fits_update_key(finfo_in, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(finfo_in, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_update_key(finfo_in, TDOUBLE, "CRVAL1", &centra, "", &status);                                              printError(status,"");  
    fits_update_key(finfo_in, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(finfo_in, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(finfo_in, TDOUBLE, "CRVAL2", &centdec, "", &status);                                    printError(status,"");   
    fits_update_key(finfo_in, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
    fits_close_file(finfo_in,&status);
    
      LOG(INFO)<<status;
     fits_open_file (&finfo_in , radecNoiceFilename.c_str() , READWRITE , &status) ;
     printError (status , "Error in opening the input information file" , infofile_in) ;
    fits_update_key(finfo_in, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(finfo_in, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_update_key(finfo_in, TDOUBLE, "CRVAL1", &centra, "", &status);                                              printError(status,"");  
    fits_update_key(finfo_in, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(finfo_in, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(finfo_in, TDOUBLE, "CRVAL2", &centdec, "", &status);                                    printError(status,"");   
    fits_update_key(finfo_in, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
    fits_close_file(finfo_in,&status);
    
    
      LOG(INFO)<<status;
    fits_open_file (&finfo_in , radecexpFilename.c_str() , READWRITE , &status) ;
     printError (status , "Error in opening the input information file" , infofile_in) ;
fits_update_key(finfo_in, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(finfo_in, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
    fits_update_key(finfo_in, TDOUBLE, "CRVAL1", &centra, "", &status);                                              printError(status,"");  
    fits_update_key(finfo_in, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(finfo_in, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(finfo_in, TDOUBLE, "CRVAL2", &centdec, "", &status);                                    printError(status,"");   
    fits_update_key(finfo_in, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(finfo_in, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
    fits_close_file(finfo_in,&status);
    
       
  /*     
    LOG(INFO)<<"Writing RA_DEC image to the disk";
    int  tfields_3 = 12 ;
    char *ttype_3[] = {"PacketSequence " , "FrameCount" , "Time" , "Fx" , "Fy" , "Max-Min" , "Min" , "BAD FLAG" , "MULT PHOTON" , "EFFECTIVE_NUM_PHOTONS","RA","DEC"} ;
    char *tform_3[] = {"U" , "U" , "D" , "1D" , "1D" , "B" , "B" , "B" , "B" , "1D","1D","1D"} ;
    char *tunit_3[] = {"" , "" , "" , "" , "" , "" , "" , "" , "" , "","",""} ;
        sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_snr , nameprefix , "radec") ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;
        fits_create_tbl (fout , BINARY_TBL , 0 , tfields_3 , ttype_3 , tform_3 , NULL , "Events" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error closing the file " , outfile) ;

*/
        LOG(INFO)<<"Diff RA "<<diffAddra<<" "<<diffAdddec;
         float *OldRA = new float[nrows-track_rown.size()];
         float *olddec = new float[nrows-track_rown.size()];
        for(int i=0;i<nrows-track_rown.size();i++)
        {
            OldRA[i]=temp_ra_newsnr[i];
            olddec[i]=temp_dec_newsnr[i];
            temp_ra_newsnr[i]=temp_ra_newsnr[i]+diffAddra;
            temp_dec_newsnr[i]=temp_dec_newsnr[i]+diffAdddec;
                    
        }
             
    char *ttype11 = {"RA-RADEC"} ;
    char *tform11 = {"D"} ;   
    char *ttype12 = {"DEC-RADEC"} ;
    char *tform12 = {"D"} ;   
    fits_open_file(&fout, outfile, READWRITE, &status);
    printError(status, " Error in opening file " , outfile);
    fits_movabs_hdu(fout,2,NULL,&status);
    printError(status, "Error in moving to other hdu ",outfile);
    
    fits_write_key(fout, TDOUBLE, "RA_PNT_RADEC", &RA_pnt, "", &status);                                              printError(status,""); 
    fits_write_key(fout, TDOUBLE, "DEC_PNT_RADEC", &DEC_pnt, "", &status);                                              printError(status,""); 
    fits_write_key(fout, TSTRING, "NAMEPRFX", nameprefix, "", &status);                                              printError(status,""); 
    fits_write_key(fout, TINT, "XSIZE", &xsize, "", &status); 
    fits_write_key(fout, TINT, "YSIZE", &ysize, "", &status); 
    fits_write_key (fout , TDOUBLE, "ROLLAPPLIED" , &rotAng, NULL , &status) ;
    fits_write_key (fout , TSTRING, "DATE-OBS" , dateOfObs, NULL , &status) ;
    fits_write_key (fout , TSTRING, "TIME-OBS" , timeobs, NULL , &status) ;
     printError (status , "***Error in updating the nameprefix keyword***") ;
    fits_insert_col (fout , 15 , ttype11 , tform11 , &status) ; //inserting 3rd column
     fits_insert_col (fout , 16 , ttype12 , tform12 , &status) ; //inserting 3rd column
     
     fits_close_file(fout,&status);
	//writeHistory(outfile,header_info);
//
/* ............................................. 27-Oct-2022 trying fix for
  +  possible error while passing arrays ...
*/
//
        status = writeColumnsToFITS (outfile , 2 , 16 , TUSHORT , 1 , psc , nrows-track_rown.size() , TUSHORT , 2 , temp_frameno_newsnr, nrows-track_rown.size() , TDOUBLE , 3 , temp_ftime_newsnr, nrows-track_rown.size() , TFLOAT , 4 , temp_X_newsnr , nrows-track_rown.size() , TFLOAT , 5 , temp_Y_newsnr ,
                nrows-track_rown.size() , TBYTE , 6 , dmm_newsnr , nrows-track_rown.size() , TBYTE , 7 , min_newsnr , nrows-track_rown.size() ,TDOUBLE,8,temp_uvittime_newsnr, nrows-track_rown.size(), TUSHORT , 9 , temp_bflag_newsnr, nrows-track_rown.size() , TUSHORT , 10 , 
                temp_mult_newsnr, nrows-track_rown.size() , TFLOAT , 11 , temp_enp__newsnr, nrows-track_rown.size(),TFLOAT,12,temp_ra_newsnr,nrows-track_rown.size(),TFLOAT,13,temp_dec_newsnr,nrows-track_rown.size(),TDOUBLE,14,Mjd_l2,nrows-track_rown.size(),TFLOAT,15,OldRA,nrows-track_rown.size(),TFLOAT,16,olddec,nrows-track_rown.size()) ;
//
//        status = writeColumnsToFITS (outfile , 2 , 16 , TUSHORT , 1 , psc , nrows-track_rown.size() , TUSHORT , 2 , frame_no , nrows-track_rown.size() , TDOUBLE , 3 , time_frame , nrows-track_rown.size() , TFLOAT , 4 , temp_X_newsnr , nrows-track_rown.size() , TFLOAT , 5 , temp_Y_newsnr ,
//                nrows-track_rown.size() , TBYTE , 6 , dmm_newsnr , nrows-track_rown.size() , TBYTE , 7 , min_newsnr , nrows-track_rown.size() ,TDOUBLE,8,temp_uvittime_newsnr, nrows-track_rown.size(), TUSHORT , 9 , badFlag_temp , nrows-track_rown.size() , TUSHORT , 10 , 
//                mult_temp , nrows-track_rown.size() , TFLOAT , 11 , effective_NumPhotons , nrows-track_rown.size(),TFLOAT,12,temp_ra_newsnr,nrows-track_rown.size(),TFLOAT,13,temp_dec_newsnr,nrows-track_rown.size(),TDOUBLE,14,Mjd_l2,nrows-track_rown.size(),TFLOAT,15,OldRA,nrows-track_rown.size(),TFLOAT,16,olddec,nrows-track_rown.size()) ;
//
//--------------------------------------------------------------------------------------
        if (status)
        {
            LOG (INFO) << "Error in writing to the shift and rotation RADEC eventFile to the disk" << endl ;
             continue;
        }
  
    fits_open_file(&fout, outfile, READWRITE, &status);
    printError(status, " Error in opening file " , outfile);
    fits_movabs_hdu(fout,2,NULL,&status);
    printError(status, "Error in moving to other hdu ",outfile);
    fits_update_key(fout,TDOUBLE,"BZERO_MJD",&bzero_mjd,NULL,&status);
    printError (status , "***Error in updating the BZERO_MJD  keyword***",outfile) ;
    fits_update_key(fout,TDOUBLE,"BSCALE_MJD",&bscale_mjd,NULL,&status);
    printError (status , "***Error in updating the BSCALE_MJD keyword***",outfile) ;
    fits_close_file(fout,&status);
	
    writeHistory(outfile,header_info);
// updateKeywords(outfile,2,2,TDOUBLE,"INT_TIME" , &integration_Time);

       delete[]  temp_ra_newsnr,temp_dec_newsnr;
      
       if(tar_extracted_flag_PC==FALSE)
       {

      string finalDir=temp_str_l2;
      string cmd="rm -r "+finalDir;
     system (cmd.c_str ());
       cmd ="mkdir -p "+finalDir;
      system (cmd.c_str ());
      LOG(INFO)<<level2outdir;
      cmd ="cp -r  "+(string)level2outdir+"/* "+finalDir;
        system (cmd.c_str ());
       status=compressTars (temp_str,finalDir,outtarpath);
       if (status)
        {
            LOG (ERROR) << "Error in compressing the tar" ;
            return (EXIT_FAILURE) ;
        }
       }
       
        //if(strcmp(datainfo.getDetector(),"NUV")==0){
        scienceDataSuccess.push_back( dirobj.sciencedatafile[dataindex]); 
       // }
          delete[] psc,frame_no,Ix,Iy,time_frame,fx, dmm,Min,multflag,badflag,time_frame_uvt,Mjd_l2 ;
delete[] one_dim_exp,one_dim_img,one_dim_img_final,one_dim_exp_final;
   
         
 }
    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::setDirectoryStructure (char *Dir , const char *subdir)
{
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s/" , Dir , subdir) ;
    LOG(INFO)<<"Directory"<< dir<<endl;
    string cmd ;
   
    if (DirExists (dir) && clobber == YES)
    {
        LOG(ERROR) << "Directory exists and clobber=yes" ;
        cmd = (string) "rm -rf " + (string) dir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists (dir) && clobber == NO)
    {
        LOG(ERROR) << endl << dir << "  already exists " ;
        LOG(ERROR) << endl << "Use clobber=yes for overwriting" ;
        return (EXIT_FAILURE) ;
    }
    cmd = "mkdir -p " + (string) dir ;

   int status= system (cmd.c_str ()) ;
    if (WIFEXITED(status))
        std::cout << "Program returned normally, exit code " << WEXITSTATUS(status) << '\n';
    else
        std::cout << "Program exited abnormaly\n"<<WEXITSTATUS(status) ;
	LOG(INFO)<<cmd;
    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::writeOutputImageToDisk (char *id , char *outDir , char *dir , char *subscript , float *Array , char *namepre , double ftime , unsigned short fno , int sizex , int sizey)
{
    int status = 0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char outfile[NAMESIZE] ;
    long naxes[2] ;
    naxes[0] = naxes[1] = sizex ;
    int naxis = 2 ;
    int bitpix = FLOAT_IMG ;
    fitsfile *fout ;
    sprintf (outfile , "%s/%s/%s_t%.4f_f%d_%s_%s.fits" , outDir , dir , namepre , ftime , fno , subscript , id) ;
    int numhdu = 0 ;
   
    if((strcmp(id,"fi") ==0) && (strcmp(subscript,"exp")==0))
    {
        
        frameIntegrationfrmname.push_back((string)outfile);
    }
    if((strcmp(id,"fi") ==0) && (strcmp(subscript,"img")==0)){
        snrimageFilename=outfile;
    }
    if(((strcmp(id,"ra-dec") ==0) || (strcmp(id,"NOT-CORRECTsigCHECK_X_Y_IMAGEONLY") ==0))&& (strcmp(subscript,"sig")==0)){
        radecFilename=outfile;
    }
    if(((strcmp(id,"ra-decnoise") ==0))){
        radecNoiceFilename=outfile;
    }
    if(((strcmp(id,"ra-decexp") ==0))){
        radecexpFilename=outfile;
    }
    if((strcmp(id,"regAvg") ==0) && (strcmp(subscript,"sig")==0)){
        regAvgFilename=outfile;
    }
    
    
    fits_create_file (&fout , outfile , &status) ;
    printError (status , "Error in creating the output Signal File" , outfile) ;
    fits_create_img (fout , bitpix , naxis , naxes , &status) ;
    printError (status , "Error in Creating the image for Signal Fie" , outfile) ;
    fits_write_pix (fout , TFLOAT , fpixel , sizex*sizey , Array , &status) ;
    printError (status , "***Error in writing the pixels to output***" , outfile) ;
  
    fits_get_num_hdus (fout , &numhdu , &status) ;
    char temp[1000] ;
    writeCommonKeywords (fout , outDir) ;
   if((strcmp(id,"fi") !=0) || (strcmp(subscript,"img")!=0)){
        for (int i = 0 ; i < numhdu ; i ++)
    {
        for (int p = 0 ; p < header_info.size () ; p ++)
        {
	//header_info[p]="HISTORY "+header_info[p];
            sprintf (temp , "%s" , header_info[p].c_str ()) ;
            if(strstr(header_info[p].c_str () ,"NAXIS")==NULL && strstr(header_info[p].c_str () ,"BITPIX")==NULL && strstr(header_info[p].c_str () ,"NAXIS1")==NULL 
                    && strstr(header_info[p].c_str () ,"NAXES2")==NULL && strstr(header_info[p].c_str(),"SIMPLE")==NULL  && strstr(header_info[p].c_str(),"EXTEND")==NULL   && strstr(header_info[p].c_str(),"DATASUM")==NULL   && strstr(header_info[p].c_str(),"CREATOR")==NULL && strstr(header_info[p].c_str(),"XSIZE")==NULL && strstr(header_info[p].c_str(),"YSIZE")==NULL && strstr(header_info[p].c_str(),"EXP_TIME")==NULL && strstr(header_info[p].c_str(),"NAMEPRFX")==NULL && strstr(header_info[p].c_str(),"CHECKSUM")==NULL  &&  strstr(header_info[p].c_str(),"ORIGIN")==NULL ){
         
   fits_write_record (fout , (char*) &temp , &status) ;
             }
        }
    }
   }
 
   
  
       if(strcmp (id,"rg")==0 || strcmp(id,"ra-dec")==0 || strcmp(id,"fi")==0 || strcmp(id,"regAvg")==0 || strcmp(id,"NOT-CORRECTsigCHECK_X_Y_IMAGEONLY")==0){
        fits_update_key (fout , TSTRING , "NAMEPRFX" , namepre , NULL , &status) ;
    printError (status , "***Error in updating the nameprefix keyword***") ;
    fits_update_key (fout , TINT, "XSIZE" , &sizex , NULL , &status) ;
    printError (status , "***Error in updating the XSIZE keyword***") ;
    fits_update_key (fout , TINT , "YSIZE" , &sizey , NULL , &status) ;
    printError (status , "***Error in updating the YSIZE keyword***") ;
    fits_update_key(fout,TFLOAT,"EXP_TIME",&Total_exp_time,NULL,&status);
char flagvalue[FLEN_FILENAME];
if(flag_Roll_Angleinvalid==1)
strcpy(flagvalue,"NOT VALID");
else 
strcpy(flagvalue,"VALID");
    printError (status , "***Error in updating the EXP_TIME keyword***") ;
fits_update_key(fout,TSTRING,"ROLLANGLE_VALID_OR_NOT",flagvalue,NULL,&status);
    printError (status , "***Error in updating the EXP_TIME keyword***") ;
    fits_update_key(fout,TFLOAT,"ROLLAPPLIED",&roll_angle_applied,NULL,&status);
    printError (status , "***Error in updating the EXP_TIME keyword***") ;
fits_update_key(fout,TDOUBLE,"BZERO_MJD",&bzero_mjd,NULL,&status);
    printError (status , "***Error in updating the BZERO_MJD  keyword***") ;
    fits_update_key(fout,TDOUBLE,"BSCALE_MJD",&bscale_mjd,NULL,&status);
    printError (status , "***Error in updating the BSCALE_MJD keyword***") ;
fits_update_key(fout,TLONG,"N_TOTAL",&ntotal,NULL,&status);
    printError (status , "***Error in updating the N_TOTAL  keyword***") ;
    fits_update_key(fout,TLONG,"N_ZERO_CENTROID",&nzerocentroid,NULL,&status);
    printError (status , "***Error in updating the N_ZERO_CENTROID keyword***") ;
   fits_update_key(fout,TDOUBLE,"CORR_FACTOR_ZEROCENTROID",&multFactor_ZeroCentroid,NULL,&status);
    printError (status , "***Error in updating the CORR_FACTOR_ZEROCENTROIDkeyword***") ;
 fits_update_key(fout,TDOUBLE,"INTEGRATION_TIME_UTC",&integrationTime_frmUTC,NULL,&status);
    printError (status , "***Error in updating the INTEGRATION_TIME_UTC keyword***") ;
    fits_update_key(fout,TDOUBLE,"INTEGRATION_TIME_CURR_OPTION",&integrationTime_frm_curr_option,NULL,&status);
    printError (status , "***Error in updating the INTEGRATION_TIME_CURR_OPTION keyword***") ;
 fits_update_key (fout , TINT , "UTC_Switch" , &UTC_flag , NULL , &status) ;
    printError (status , "***Error in updating the YSIZE keyword***") ;
    fits_update_key (fout , TLONG , "TOTAL INPUT PACKET FROM L1" , &total_InputPacket, NULL , &status) ;
    printError (status , "Error in Writing the TOTAL INPUT PACKET FROM L1 " ) ;
    fits_update_key (fout , TLONG , "Starting PACKET NUMBER" , &start_pktNo, NULL , &status) ;
    printError (status , "Error in Writing the Starting PACKET NUMBER ") ;
    fits_update_key (fout , TLONG , "Total Used packet" , &usedpkts, NULL , &status) ;
    printError (status , "Error in Writing the Total Used packet " ) ;
     fits_update_key (fout , TINT , "NFRAMES" , &cntTotalFrameInput, NULL , &status) ;
    printError (status , "Error in Writing the Total Used packet " ) ;
   fits_update_key (fout , TDOUBLE , "CRC_CORRECTION_FACTOR" , &CrcFactor, NULL , &status) ;
    printError (status , "Error in Writing the Total Used packet " ) ;
    fits_update_key (fout , TDOUBLE , "PARITY_CORRECTION_FACTOR" , &ParityFactor, NULL , &status) ;
    printError (status , "Error in Writing the Total Used packet " ) ;
    fits_write_key (fout , TLONG , "Total Packet from L1" , &pktsFromL1 , "" , &status) ;
fits_write_key (fout , TLONG , "CRC_TotalPackets_Failed" , &crcfailedpckts , "" , &status) ;
fits_write_key (fout , TLONG , "Total remaining packet from L1" , &pcktremained , "" , &status) ;
fits_write_key (fout , TLONG , "PARITY_TotalEvents_Failed" , &parityfailedEvents , "" , &status) ;
fits_write_key (fout , TLONG , "Total Events Used" , &eventsremained , "" , &status) ;
fits_write_key (fout , TDOUBLE , "PEAK_OF_EXP" , &peakValueExp , "" , &status) ;
 fits_write_key(fout, TDOUBLE, "RA_PNT_RADEC", &RAPNT_FRMATT, "", &status);                                              printError(status,""); 
    fits_write_key(fout, TDOUBLE, "DEC_PNT_RADEC", &DECPNT_FRMATT, "", &status);
       }
    fits_close_file (fout , &status) ;
    printError (status , "Error in closing the  output Signal fits file" , outfile) ;

    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::copyAllheaderKeys (char* infile)
{

    int status = 0 , keyexist ;
    char record[FLEN_CARD] ;
    header_info.clear () ;
    fitsfile *fin ;
    fits_open_file (&fin , infile , READONLY , &status) ;
    if (status) return (EXIT_FAILURE) ;
    fits_get_hdrspace (fin , &keyexist , NULL , &status) ;
    if (status)
    {
        LOG(INFO) << endl << "***Could not find number of keywords in file " << infile << "***" ;
        fits_report_error (stderr , status) ;
        return (EXIT_FAILURE) ;
    }
    int keyclass ;
    //LOG(INFO)<<"\nNumber of keywords found:"<<keyexist;
    for (int i = 1 ; i <= keyexist ; i ++)
    {
        fits_read_record (fin , i , record , &status) ;
        if (status)
        {
            LOG(INFO) << endl << "***Error in reading record number " << i << " in file " << infile << "***" ;
            fits_report_error (stderr , status) ;
            return (EXIT_FAILURE) ;
        }
        //keyclass = fits_get_keyclass (record) ;
       // if (keyclass == TYP_COMM_KEY)
       //     continue ;
        //else if (keyclass == TYP_USER_KEY || keyclass == TYP_REFSYS_KEY || keyclass == TYP_WCS_KEY)
       // {
            header_info.push_back (record) ;
       // }
    }

    fits_close_file (fin , &status) ;
    // delete[] record;
    return (EXIT_SUCCESS) ;
}





double uvtLevel2PC::readDarkFrame (char * path , float *Array)
{

    fitsfile *fptr ;
    int status = 0 ;
    long fpixel[2] ;
    for (int i = 0 ; i < xsize * ysize ; i ++) Array[i] = (float) 0.0 ;
    fpixel[0] = fpixel[1] = 1 ;
    char errstr[500] ;
    fits_open_file (&fptr , path , READONLY , &status) ;
    printError (status , "Error in opening the Dark frame" , path) ;

    fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , Array , NULL , &status) ;

    printError (status , "Error in reading the pixels from the Dark Frames" , path) ;

    double frame_time = 0 ;
    fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frame_time , NULL , &status) ;
    printError (status , errstr , path) ;

    fits_close_file (fptr , &status) ;
    printError (status , "Error in Closing the  Dark Frame" , path) ;
    return frame_time ;

}





int uvtLevel2PC::readImage (char * caldb_file , int hduno , float *frm_data)
{

    fitsfile *fptr ;
    int status = 0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    fits_open_file (&fptr , caldb_file , READONLY , &status) ;
    printError (status , "Error in readFlatField()" , caldb_file) ;
    fits_movabs_hdu (fptr , hduno , NULL , &status) ;
    printError (status , "Error in  moving to the 2nd HDU of the out information file" , caldb_file) ;
    for (int q = 0 ; q < xsize * ysize ; q ++) frm_data[q] = 0.0 ;
    fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , frm_data , NULL , &status) ;
    printError (status , "Error in reading the pixels from caldb file" , caldb_file) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the file" , caldb_file) ;
    return (EXIT_SUCCESS) ;
}

int uvtLevel2PC::writeOutputTblToDisk (char *id , char *outDir , char *dir , char *subscript , char *namepre , double ftime , unsigned short fno , char **type1 , char**tform1 , int tfields1 , char *extname , vector<float> &X , vector<float> &Y , vector<float> &val)
{
    char out_file[FLEN_FILENAME] ;
    sprintf (out_file , "%s/%s/%s_t%f_f%d_%s_%s.fits" , outDir , dir , namepre , ftime , fno , id , subscript) ;
    if (strcmp (id , "rf") == 0)
    {
        name_track.push_back (basename (out_file)) ;
    }
    fitsfile *fout1 ;
    int status = 0 ;
    fits_create_file (&fout1 , out_file , &status) ;
    printError (status , "Error creating the output file " , out_file) ;
    fits_create_tbl (fout1 , BINARY_TBL , 0 , tfields1 , type1 , tform1 , NULL , extname , &status) ;
    printError (status , "Error in creating the table" , out_file) ;
    fits_update_key (fout1 , TDOUBLE , "FRMTIME" , &ftime , " Average Frame time" , &status) ;
    printError (status , "Error in updating the key value of the FRMTIME" , out_file) ;
    fits_write_col (fout1 , TFLOAT , 1 , 1 , 1 , X.size () , X.data () , &status) ;
    printError (status , "Error in writing  column" , out_file) ;
    fits_write_col (fout1 , TFLOAT , 2 , 1 , 1 , Y.size () , Y.data () , &status) ;
    printError (status , "Error in writing the column" , out_file) ;
    fits_write_col (fout1 , TFLOAT , 3 , 1 , 1 , val.size () , val.data () , &status) ;
    printError (status , "Error in writing the column" , out_file) ;
    fits_close_file (fout1 , &status) ;
    return (EXIT_SUCCESS) ;

}


int uvtLevel2PC::performDistortionCorr (vector<float> &X , vector<float> &Y , float *Xdistr , float *Ydistr , int sizex , long caldbsize)
{
    float multi_factor = sizex / 600 ;
    if (multi_factor == 0)
    {
        LOG(ERROR) << "***Divide by Zero***" << endl ;
        return (EXIT_FAILURE) ;
    }
    float tempx , tempy ;
    for (int p = 0 ; p < X.size () ; p ++)
    {

        tempx = X[p] / multi_factor-44 ;
        tempy = Y[p] / multi_factor-44 ;
//ROUNDING:No correction is needed,As rounding is applied on each direction
        float locate = ((int) round (tempy) - 44) * caldbsize + ((int) round (tempx) - 44) ;
        round (locate) ;
LOG(INFO)<<locate;
	
        tempx = tempx - Xdistr[(int) locate] ;
        tempy = tempy - Ydistr[(int) locate] ;
	

        X[p] = (tempx+44) * multi_factor ;
        Y[p] = (tempy+44) * multi_factor ;


    }
    return (EXIT_SUCCESS) ;
}




int uvtLevel2PC::takeDarkinfo ()
{

    double dark_Firstframetime , dark_Secframetime , start_i , start_f , stop_i , stop_f ;
    char darkpath_First[FLEN_FILENAME] ;
    char darkend_Sec[FLEN_FILENAME] ;
    char dark_temp[FLEN_FILENAME] ;

    //vector for storing the dark frame names
    vector<string> dark_framenames ;

    //setting path for dark frame directory  
    sprintf (dark_temp , "%s/%s" , dataIngest_out_dir , darkdir) ;

    getFileNamesfrmDir (dark_temp , ".fits" , dark_framenames) ; //get dark file names from dark directory

    if (dark_framenames.size () < 2)
    {
        LOG (INFO) << "No enough Dark frames found at input Directory" << endl ;
        return (EXIT_SUCCESS) ;
    }

    //setting first darkframe's path
    //for dark frame 1
    sprintf (darkpath_First , "%s/%s" , dark_temp , dark_framenames[0].c_str ()) ;

    readKeywords (darkpath_First , 1 , 1 , TDOUBLE , "FRMTIME" , &dark_Firstframetime) ;

    //calculate first dark frame's time by adding fraction and integer portion of  time
    //   dark_Firstframetime = start_i + start_f ;

    //for dark frame 2 
    //setting  second darkframe's path
    sprintf (darkend_Sec , "%s/%s" , dark_temp , dark_framenames[1].c_str ()) ;

    readKeywords (darkend_Sec , 1 , 1 , TDOUBLE , "FRMTIME" , &dark_Secframetime) ;

    //calculate second dark frame's time by adding fraction and integer portion of  time
    //dark_Secframetime = start_i + start_f ;

    //comparing darkframe's time  to identify dark begin and dark end.
    if (dark_Firstframetime < dark_Secframetime)
    {
        strcpy (dstartpath , darkpath_First) ;
        strcpy (dendpath , darkend_Sec) ;

    }
    else
    {
        strcpy (dstartpath , darkend_Sec) ;
        strcpy (dendpath , darkpath_First) ;
    }

    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::performQEMCPcorrection (float *sigdata , int size , double fact)
{

    for (int i = 0 ; i < size ; i ++)
    {
        sigdata[i] = sigdata[i] * fact ;
    }


    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::performMaskBadPix (unsigned short* int_x , unsigned short* int_y , float* badpixArr , unsigned char * Max_Min , unsigned char* badflag , unsigned char* multflag , int nrows , int x_size , int y_size , float thr_multphn)
{
    bool WinDow_Mode=FALSE;
    if(x_size !=IMG_DIM_DI || y_size!=IMG_DIM_DI)
    {
        WinDow_Mode=TRUE;
        
    }
    int Total_MultPhoton_effected_frames=0;
    for (int i = 0 ; i < nrows ; i ++)
    {
        int *fg = 0 ;
        badflag[i] = static_cast<unsigned char> (badpixArr[int_y[i] * x_size + int_x[i]]) ;
        if(WinDow_Mode==TRUE){
            if(int_y[i]==datainfo.yoff+1 || int_y[i]==datainfo.yoff+ysize+1 || int_x[i] == datainfo.xoff+1 || int_x[i]==datainfo.xoff+xsize+1)
                badflag[i]=0;
        }
        
        
        if ((float) Max_Min[i] > thr_multphn){
            multflag[i] = 0 ; //if max-min > user_threshold then flag=0 else flag=1
            Total_MultPhoton_effected_frames++;
        } 
        else multflag[i] = 1 ;

    }
    //added
    
    
    
    //till this
    

    return Total_MultPhoton_effected_frames;
}


int uvtLevel2PC::performPixPadding (unsigned short *X_int , unsigned short *Y_int , int nrows)
{
    for (int i = 0 ; i < nrows ; i ++)
    {
        X_int[i] = X_int[i] + 44 ;
        Y_int[i] = Y_int[i] + 44 ;
    }
    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::performCentroidCorr (double *t , unsigned short *xint , unsigned short *yint , float *xf , float *yf , unsigned char *mc , float *newXfract , float *newYfract , float *darkbeginData , float *darkenddata , double integrtn_time , long EA , int nrows)
{
    float dark_value = 0.0 ;
    float sd = 0.0 , sc = 0 ;
    double t1 , t2 ;
    //   double *t = new double[nrows] ;
    t2 = t[nrows - 1] + integrtn_time ;
    int index_x ;
    int index_y ;
    int diff_add = xsize - darkSize ;
    // LOG(INFO)<<"p"<<endl;exit(1);
    for (int i = 0 ; i < nrows ; i ++)
    {
        //            LOG(INFO)<<i<<endl;
        dark_value = 0.0f ;
        t1 = t[i] - integrtn_time ;
        if ((t2 - t1) == 0)
        {
            LOG (ERROR) << "***Divide by Zero***" << endl ;
            return (EXIT_FAILURE) ;
        }
        index_x = (int) xint[i] ;
        index_y = (int) yint[i] ;
        newXfract[i] = 0.0f ;
        newYfract[i] = 0.0f ;
        if ((index_x - diff_add) >= 0 && (index_y - diff_add) >= 0)
        {
            dark_value = darkbeginData[(index_y - diff_add) * darkSize + (index_x - diff_add)] +
                    ((t[nrows - 1] - t1) * ((darkenddata[(index_y - diff_add) * darkSize + (index_x - diff_add)]
                    - darkbeginData[(index_y - diff_add) * darkSize + (index_x - diff_add)]) / (t2 - t1))) ;
        }
        sc = mc[i] ;
        sd = dark_value ;
        newXfract[i] = xf[i]*(1 + ((sd - sc) * cent_corr_win_size) / EA) + index_x ;
        newYfract[i] = yf[i]*(1 + ((sd - sc) * cent_corr_win_size) / EA)+ index_y ;
    }

    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::readcentroidbiasFile ()
{
    fitsfile *fcbias ;
    int status = 0 ;
    long nrows ;
    fits_open_file (&fcbias , centroidbiasfile , READONLY , &status) ;
    printError (status , "Error opening Centroid Bias File" , centroidbiasfile) ;
    fits_movabs_hdu (fcbias , 2 , NULL , &status) ;
    printError (status , "Error Moving to 2nd HDU" , centroidbiasfile) ;
    fits_get_num_rows (fcbias , &nrows , &status) ;
    printError (status , "Error Reading the number of Rows" , centroidbiasfile) ;
    biasRows = nrows ;
    fraction_bias = new double[nrows] , x_corr = new double[nrows] , y_corr = new double[nrows] ;
    fits_read_col (fcbias , TDOUBLE , 1 , 1 , 1 , nrows , NULL , fraction_bias , NULL , &status) ;
    printError (status , "Error reading  centroid x" , centroidbiasfile) ;
    fits_read_col (fcbias , TDOUBLE , 2 , 1 , 1 , nrows , NULL , x_corr , NULL , &status) ;
    printError (status , "Error writing x-correction" , centroidbiasfile) ;
    fits_read_col (fcbias , TDOUBLE , 3 , 1 , 1 , nrows , NULL , y_corr , NULL , &status) ;
    printError (status , "Error writing y-correction" , centroidbiasfile) ;
    fits_close_file (fcbias , &status) ;
    printError (status , "Error closing caldb file" , centroidbiasfile) ;
    LOG (INFO) << "Reading Centroid Bias file from calDB completed" << endl ;
    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::performCentroidBias (long nrows , float *xFrac , float *yFrac , float *new_xFrac , float *new_yFrac)
{
    
   float *xFrac_temp = new float[nrows] ;
float *yFrac_temp = new float[nrows] ;
int *New_X_Val= new int [nrows];
int *New_Y_Val= new int [nrows];
int temp_Xval,temp_Yval;

for (int i=0;i<nrows;i++)
{
    xFrac_temp[i] =( (xFrac[i]-(int) xFrac[i])*32);
    yFrac_temp[i] =( (yFrac[i]-(int) yFrac[i])*32);
    
    New_X_Val[i]=((int) xFrac[i])*32;
    New_Y_Val[i]=((int) yFrac[i])*32;
    //LOG(INFO)<<xFrac_temp[i]<<" "<<yFrac_temp[i];
    if(xFrac_temp[i]<-16){
        xFrac_temp[i]=xFrac_temp[i]+32;
        New_X_Val[i]=New_X_Val[i]-32;
    }
    if(yFrac_temp[i]<-16){
        yFrac_temp[i]=yFrac_temp[i]+32;
         New_Y_Val[i]=New_Y_Val[i]-32;
    }
    if(xFrac_temp[i]>15){
        xFrac_temp[i]=xFrac_temp[i]-32;
        New_X_Val[i]=New_X_Val[i]+32;
    }
    if(yFrac_temp[i]>15){
        yFrac_temp[i]=yFrac_temp[i]-32;
        New_Y_Val[i]=New_Y_Val[i]+32;
    }
    temp_Xval=New_X_Val[i]/32;
    temp_Yval=New_Y_Val[i]/32;
    //For Converting from +-0.5 to 0-1
    xFrac_temp[i]=xFrac_temp[i]+16;
    yFrac_temp[i]=yFrac_temp[i]+16;
  //  LOG(INFO)<<xFrac_temp[i]<<" "<<yFrac_temp[i];
    //rand->nextDouble
    //FPN correction
    for (int j = 0 ; j < biasRows ; j ++)
       {
      
              if (xFrac_temp[i]>= fraction_bias[j] && xFrac_temp[i] < fraction_bias[j + 1])
              {
                       xFrac_temp[i] = xFrac_temp[i]*((x_corr[j]) + ((xFrac_temp[i] - fraction_bias[j]) * (x_corr[j + 1] - x_corr[j])) / (fraction_bias[j + 1] - fraction_bias[j])) ; 
                       break;
              }
              
       }
   // exit(1);
     for (int j = 0 ; j < biasRows ; j ++)
       {
              if (yFrac_temp[i] >= fraction_bias[j] && yFrac_temp[i] < fraction_bias[j + 1])
              {
                       yFrac_temp[i] =yFrac_temp[i]*(( y_corr[j] )+ ((yFrac_temp[i] - fraction_bias[j]) * (y_corr[j + 1] - y_corr[j])) / (fraction_bias[j + 1] - fraction_bias[j])); 
                       break;
              }              
       }

    //For reverse conversion
    xFrac_temp[i]=xFrac_temp[i]-16;
    yFrac_temp[i]=yFrac_temp[i]-16;
//    new_xFrac[i]=xFrac_temp[i]+New_X_Val[i];
//    new_yFrac[i]=yFrac_temp[i]+New_Y_Val[i];
    
    new_xFrac[i]=temp_Xval+xFrac_temp[i]/32;
    new_yFrac[i]=temp_Yval+yFrac_temp[i]/32;
   //LOG(INFO)<<new_xFrac[i]<<" "<<new_yFrac[i];
    //LOG(INFO)<<new_xFrac[i]<<" "<<new_yFrac[i];
}
//exit(1);
return(EXIT_SUCCESS);
}


int uvtLevel2PC::performDistCorrection (long nrows , float *xFrac , float *yFrac , float *x_Distortion , float *y_Distortion , int caldbdim)
{
/*
    float locate ;
    float tempx , tempy , multi_factor = 0.0 ;
    multi_factor = xsize / 600 ;
    LOG (INFO) << "Applying  Distortion Correction..." << nrows<<" "<<caldbdim<<endl ;
    double xFracValue,yFracValue,newXfracVal,newYfracVal;
    double value_X_toAdd,value_Y_toAdd;
    long curr_loc=0;
    long next_loc=0;
    for (int k = 0 ; k < nrows ; k ++)
    {
        if(xFrac[k]!=INVALID_PIX_VALUE && yFrac[k]!=INVALID_PIX_VALUE){
        tempx = (xFrac[k] / multi_factor)-44 ;
        tempy =(yFrac[k] / multi_factor)-44 ;
         
       //ROUNDING:No correction is needed.As round is applied on each  direction.
        locate = ((int) round (tempy)) * caldbdim + ((int) round (tempx) ) ;
     
               curr_loc=(int)locate-1;
               next_loc=curr_loc+1;
                if(curr_loc>0 && next_loc<512*512 && locate > 0 && locate <512*512){  
               value_X_toAdd=x_Distortion[curr_loc]+((x_Distortion[next_loc]-x_Distortion[curr_loc]))*(locate-curr_loc);
               value_Y_toAdd=y_Distortion[curr_loc]+((y_Distortion[next_loc]-y_Distortion[curr_loc]))*(locate-curr_loc);
     
        //ROUNDING:this round is not needed.because it is already a nearest integer value.
       
           tempx = tempx  - (value_X_toAdd) ;
           tempy = tempy -  (value_Y_toAdd) ;
 
        if ((tempx) * multi_factor < xsize && (tempy )* multi_factor < ysize   && (tempx) * multi_factor >0 && (tempy) * multi_factor > 0)
        {
            xFrac[k] = (tempx+44) * multi_factor ;
            yFrac[k] = (tempy+44) * multi_factor ;
            
        }
        else{
            xFrac[k] = INVALID_PIX_VALUE ;
            yFrac[k] = INVALID_PIX_VALUE ;  
        }
        }
        else{
            xFrac[k] = INVALID_PIX_VALUE ;
            yFrac[k] = INVALID_PIX_VALUE ;  
            //ListOf_OutsideEvnts.push_back(k+1);
        }
        }
        else{
             xFrac[k] = INVALID_PIX_VALUE ;
            yFrac[k] = INVALID_PIX_VALUE ;  
        }
     
    }
    return (EXIT_SUCCESS) ;
 */
///*   
     float locate ;
    float tempx , tempy , multi_factor = 0.0 ;
    long k1,k2;
    multi_factor = xsize / 600 ;
    LOG (INFO) << "Applying  Distortion Correction..." << nrows<<" "<<caldbdim<<endl ;
    double xFracValue,yFracValue,newXfracVal,newYfracVal;
    double value_X_toAdd,value_Y_toAdd;
    long curr_loc=0;
    long next_loc=0;
    long nextPixX=0,nextPixY=0;
    double tempAddx_X,tempAddy_X,tempAddx_Y,tempAddy_Y,x1,x2, FinalXinDistAdd_X, FinalXinDistAdd_Y;
    for (int k = 0 ; k < nrows ; k ++)
    {
        if(xFrac[k]!=INVALID_PIX_VALUE && yFrac[k]!=INVALID_PIX_VALUE){
        x1=tempx = (xFrac[k] / multi_factor)-44 ;
        x2=tempy =(yFrac[k] / multi_factor)-44 ;
         
        k1=(int)((xFrac[k] / multi_factor)-44);
        k2=(int)((yFrac[k] / multi_factor)-44);
        
        
        nextPixX=k1+1;
        nextPixY=k2+1;
         
         
         tempAddx_X=(((nextPixX)-(x1))*x_Distortion[(k2)*(512)+(k1)])+(((x1)-(k1))*x_Distortion[(k2)*(512)+(k1+1)]);
         tempAddy_X=(((nextPixX)-(x1))*x_Distortion[(k2+1)*(512)+(k1)])+(((x1)-(k1))*x_Distortion[(k2+1)*(512)+(k1+1)]);    
         
         
         FinalXinDistAdd_X=((k2+1)-(x2))*tempAddx_X+((x2)-(k2))*tempAddy_X;
         
         
         tempAddx_Y=(((nextPixX)-(x1))*y_Distortion[(k2)*(512)+(k1)])+(((x1)-(k1))*y_Distortion[(k2)*(512)+(k1+1)]);
         tempAddy_Y=(((nextPixX)-(x1))*y_Distortion[(k2+1)*(512)+(k1)])+(((x1)-(k1))*y_Distortion[(k2+1)*(512)+(k1+1)]);
         
         
         FinalXinDistAdd_Y=((k2+1)-(x2))*tempAddx_Y+((x2)-(k2))*tempAddy_Y;
         
       //ROUNDING:No correction is needed.As round is applied on each  direction.
        locate = ((int) round (tempy)) * caldbdim + ((int) round (tempx) ) ;
     
               curr_loc=(int)locate-1;
               next_loc=curr_loc+1;
                if(curr_loc>0 && next_loc<512*512 && locate > 0 && locate <512*512){  
               value_X_toAdd=x_Distortion[curr_loc]+((x_Distortion[next_loc]-x_Distortion[curr_loc]))*(locate-curr_loc);
               value_Y_toAdd=y_Distortion[curr_loc]+((y_Distortion[next_loc]-y_Distortion[curr_loc]))*(locate-curr_loc);
     
        //ROUNDING:this round is not needed.because it is already a nearest integer value.
       
           tempx = tempx  - (FinalXinDistAdd_X) ;
           tempy = tempy -  (FinalXinDistAdd_Y) ;
 
        if ((tempx) * multi_factor < xsize && (tempy )* multi_factor < ysize   && (tempx) * multi_factor >0 && (tempy) * multi_factor > 0)
        {
            xFrac[k] = (tempx+44) * multi_factor ;
            yFrac[k] = (tempy+44) * multi_factor ;
            
        }
        else{
            xFrac[k] = INVALID_PIX_VALUE ;
            yFrac[k] = INVALID_PIX_VALUE ;  
        }
        }
        else{
            xFrac[k] = INVALID_PIX_VALUE ;
            yFrac[k] = INVALID_PIX_VALUE ;  
            //ListOf_OutsideEvnts.push_back(k+1);
        }
        }
        else{
             xFrac[k] = INVALID_PIX_VALUE ;
            yFrac[k] = INVALID_PIX_VALUE ;  
        }
      
    }
    return (EXIT_SUCCESS) ;
    
    
 //   */
    
    

}
int uvtLevel2PC::performFrameIntegration (long nrows , long start_row,unsigned short*frame_no , int xsize , int Ndiscard , int Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr_strct> &vect,
        float *noice_map_Sig_Array,float *outputSigArr,float *outputExpArray,float &outputFrmtime,long &outputFrmNo,long &last_index_for_frame)
{
    for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++) 
    {
        outputSigArr[i]=0.0f;
        outputExpArray[i]=0.0f;
    }
    FrameIntegration_Arr_strct obj ;
    int jj = 0 ;
    int x_tem , y_tem ;
    double Start_Time , End_Time , Avg_Time = 0.0 ;

   // long firstpix1[2] ;
   // firstpix1[0] = firstpix1[1] = 1 ;
  //  int naxis = 2 ;
   // long naxes[2] = {IMG_DIM_FI , IMG_DIM_FI} ;
    float Eff_NoPhtn = 0.0f ;
    float multi_factor = IMG_DIM_FI / xsize ;
  //  int file_num = 1 ;
    int frame_check ;
    LOG (INFO) << "MULTIPLICATION FACTOR  for converting image to "<<IMG_DIM_FI<<" 9600 is " << multi_factor << endl ;
   // int frame_init = frame_check = Nacc ;
    double time_frame = 0.0 , time_frame_final = 0.0 ;
    int cnt_frame = 0 ;
  //  vector<string> vhistorystr ;
    
    /****provision for if frames are not starting from 1****/
   // int start_row = 0 ;
    int cmpr_term = 0 ;
 
     if(Nacc>frame_no[nrows-1]){
               LOG(ERROR)<<"***Number of frames to accumulate is greater than total available frames";
                return(EXIT_FAILURE);
    }
    if (frame_no[0] != 1)
    {
        cmpr_term = Ndiscard + (frame_no[0] - 1) ;
    }

    float *badPix_Arr_padded= new float [600*600];
    for(int i=0;i<600*600;i++) badPix_Arr_padded[i]=0.0;
    float *badPix_Arr_fi= new float[IMG_DIM_FI*IMG_DIM_FI];
      for(int i=0;i<600*600;i++) badPix_Arr_fi[i]=0.0;
   int   status = Applypadding (badpixdata,512,512,badPix_Arr_padded,600,600);
    if (status)
    {
        LOG (ERROR) << "Error in padding for bad pix data" ;
        return (EXIT_FAILURE) ;
    }
    status= performSubDivisionIM (badPix_Arr_padded,600,600,badPix_Arr_fi,IMG_DIM_FI,IMG_DIM_FI);
   if (status)
    {
        LOG (ERROR) << "Error in performing sub division in badpixel data" ;
        return (EXIT_FAILURE) ;
    } 
    long int tot_num_frames = 0 ;
    for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++)
    {
    one_dim_img[i]=0.0f;
    //cout<<badpixdata[i]<<endl;
    noice_map_Sig_Array[i]=0.0f;
    one_dim_exp[i]=badPix_Arr_fi[i]*Nacc;
    }
    
    LOG (INFO) << "Performing Frame Integration...." << endl ;
    long p = 0 ;
    
    //calculate final frame count value for the current output frame.
    int nrows_cur_frm=start_row;
    for (int i=start_row;i<nrows;i++)
    {
  if (frame_no[i]<frame_no[start_row]+Ndiscard+Nacc)
    {  
      nrows_cur_frm++;
  //LOG(INFO)<<i;
     }
    }
    last_index_for_frame=nrows_cur_frm;
  // LOG(INFO)<<"FFF "<<nrows<<" "<<nrows_cur_frm<<start_row;
  //LOG(INFO)<<"PPP "<<start_row<<" "<<nrows_cur_frm;
    int NframesCount=0;
    //unsigned short frmno_previous=frame_no[NframesCount];
    for (int k = start_row ; k < nrows_cur_frm ; k ++)
    {
        /**in 'IF condition' accumulated array of  Nacc frames will be  generate   & In 'ELSE condition' array will be written to the frame.
         x and y location  for number of events are taken as a indexes for the array .For each Nacc events one array will be written as a frame. 
         **/
      
        if (frame_no[k] < (Nacc+frame_no[start_row]))
        {
            
            while (jj == 0)
            {
                Start_Time = t[k] ;
                jj ++ ;
            }

            x_tem = (int)round(xFrac[k]) ;
            y_tem = (int)round(yFrac[k]) ;
            x_tem = x_tem * multi_factor ;
            y_tem = y_tem * multi_factor ;


            if (y_tem < IMG_DIM_FI && x_tem < IMG_DIM_FI)
            {

                float mult = (float) mult_phn[k] ;
                float badflg = (float) bad_Flag[k] ;

                 Eff_NoPhtn = (float) ENP[k] ;
//%#Edited ON 20July17#% 
               if(y_tem> 0 && y_tem<IMG_DIM_FI && x_tem>0 && x_tem<IMG_DIM_FI)
//%#-Till this-20July17#%
{
                one_dim_img[y_tem*IMG_DIM_FI+x_tem]=one_dim_img[y_tem*IMG_DIM_FI+x_tem]+Eff_NoPhtn*mult*badflg; 
                noice_map_Sig_Array[y_tem*IMG_DIM_FI+x_tem]=noice_map_Sig_Array[y_tem*IMG_DIM_FI+x_tem]+Eff_NoPhtn*mult*badflg;
               }

            }

            if (frame_no[k - 1] == frame_no[k])
            {
                time_frame = t[k] ;
            }
            else
            {
                time_frame_final = time_frame_final + time_frame ;
            }
                        if (k == nrows_cur_frm - 1)
                        {
                            goto label_else ;
                        }
        } 
        else
        {
             
label_else:

            jj = 0 ;
            End_Time = t[k - 1] ;

            Avg_Time = Start_Time + (End_Time - Start_Time) / 2 ;

            //long temp = 0 ;

            tot_num_frames ++ ;
           LOG(INFO)<< "\rFrame created " << tot_num_frames  ;
            cnt_frame = 0 ;
            time_frame = 0.0 ;
            time_frame_final = 0.0 ;

             for (int i = 0 ; i < IMG_DIM_FI * IMG_DIM_FI ; i ++)
            {
                 outputSigArr[i]=one_dim_img[i];
                 outputExpArray[i]=one_dim_exp[i];                
            }
            outputFrmtime=Avg_Time;
            outputFrmNo=tot_num_frames;

        }

    }
     delete[] badPix_Arr_fi,badPix_Arr_padded;
    return (EXIT_SUCCESS) ;
}

int uvtLevel2PC::performUnitConversion (float *frmdata , double intgrntime , int size)
{
    
    if (intgrntime <= 0)
    {
        LOG(ERROR) << "***Divide by Zero***integration time may be zero or less than zero" ;
        return (EXIT_FAILURE) ;

    }
    LOG(INFO)<<"Integration Time ->"<<intgrntime<<endl;
    for (int i = 0 ; i < size ; i ++)
    {
        frmdata[i] = frmdata[i] / intgrntime ;

    }
    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::performFlatFieldCorr (float* frmsigdata , float* frmflatfield , float *Xfract , float *Yfract , int size)
{
    int div_fact = subDivision_size / 600 ;
    for (int i = 0 ; i < size ; i ++)
    {
        //ROUNDING:No correction  is needed.As rounding is applied on each direction.
        frmsigdata[i] = frmsigdata[i] * frmflatfield[(int) (round ((Yfract[i] / div_fact))*600 + (int) round ((Xfract[i]) / div_fact))] ;
    }
    return (EXIT_SUCCESS) ;

}


int uvtLevel2PC::performShiftNRot (string sciencedataFile,long nrows ,unsigned short *frmno, double *t , float *xf , float *yf , float *xi_fi , float *yi_fi,vector<float> &Shifts_x,vector<float> &Shifts_y,
        vector<float> &Shifts_theta,vector<unsigned short> &frm_no_vect,unsigned short *multphotonFlagIn,vector<unsigned short> &multphotonFlagOut)
{
    
    
//  
    int  status = 0 ;
    //opening attitude file
    
    
    
    float i1 = 0 , j1 = 0 ;
    double time_frame = 0.0 ;
    bool flag = FALSE ;

    //  float *frm_no = new long[nrows];
    for (int i = 0 ; i < nrows ; i ++)
    {
        xi_fi[i] = 0.0f ;
        yi_fi[i] = 0.0f ;
    }
   // LOG(INFO)<<xsize;exit(1);
    float mult_fact = xsize / 600 ;
    double t1 , t2 , theta1 , theta2 , x1 , x2 , y1 , y2 ;
    double new_delta_theta = 0.0 , new_delta_x = 0.0 , new_delta_y = 0.0 ;
    fitsfile *fras ;
    long no_of_records ;
    
   //LOG(INFO)<<"The RAS file used "<<ras_finalFile;
    fits_open_file (&fras ,rasfile , READONLY , &status) ;
    printError (status , "Error in opening the rasfile" , rasfile) ;
    fits_movabs_hdu (fras , 3 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of  ras file" , rasfile) ;
    fits_read_key (fras , TINT , "ThetaFlag" ,&flag_Theta , NULL , &status) ;
    printError (status , "Error reading the key value of the INT_TIME" , rasfile) ;
    fits_movabs_hdu (fras , 2 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of  ras file" , rasfile) ;
    fits_get_num_rows (fras , &no_of_records , &status) ;
    printError (status , "Error Reading the number of Rows" , rasfile) ;


    time_drifts = new double[no_of_records] ;
    roll_ras = new double [no_of_records] ;
    pitch_ras = new double [no_of_records] ;
    yaw_ras = new double[no_of_records] ;

    delta_x = new double [no_of_records] ;
    delta_y = new double [no_of_records] ;
    delta_theta = new double[no_of_records] ;
    fits_read_col (fras , TDOUBLE , 1 , 1 , 1 , no_of_records , NULL , time_drifts , NULL , &status) ;
    printError (status , "***Error reading  centroid x***") ;
    fits_read_col (fras , TDOUBLE , 2 , 1 , 1 , no_of_records , NULL , yaw_ras , NULL , &status) ;
    printError (status , "***Error writing centroid y***") ;
    fits_read_col (fras , TDOUBLE , 3 , 1 , 1 , no_of_records , NULL , pitch_ras , NULL , &status) ;
    printError (status , "***Error writing x-correction***") ;
    fits_read_col (fras , TDOUBLE , 4 , 1 , 1 , no_of_records , NULL , roll_ras , NULL , &status) ;
    printError (status , "***Error writing y-correction***") ;
    fits_close_file (fras , &status) ;
    printError (status , "Error in closing the file" , rasfile) ;



    double temp_pitch , temp_yaw ;
    if (strcasecmp (datainfo.getDetector () , "NUV") == 0)
    {
        for (int i = 0 ; i < no_of_records ; i ++)
        {
            temp_pitch = pitch_ras[i]*3600 ;
            temp_yaw = yaw_ras[i]*3600 ;
            transformRPYtoDXDYDTHETA_NUV (roll_ras[i] , temp_pitch , temp_yaw , delta_x[i] , delta_y[i] , delta_theta[i]) ;      
        }

    }
    else if (strcasecmp (datainfo.getDetector () , "FUV") == 0)
    {
        double tempx,tempy,temptheta;
	
        for (int i = 0 ; i < no_of_records ; i ++)
        {
             temp_pitch = pitch_ras[i]*3600 ;
            temp_yaw = yaw_ras[i]*3600 ;
            transformRPYtoDXDYDTHETA_FUV (roll_ras[i] , temp_pitch , temp_yaw , delta_x[i] , delta_y[i] , delta_theta[i]) ;
 // transformRPYtoDXDYDTHETA_FUV (roll_ras[i] , temp_pitch , temp_yaw , tempx , tempy , temptheta) ;
           
//transformFUVtoNUV(tempx,tempy,temptheta,delta_x[i],delta_y[i],delta_theta[i]);    
//          LOG(INFO)<<delta_x[i]<<" "<<delta_y[i];
        }
    }
    else if (strcasecmp (datainfo.getDetector () , "VIS") == 0)
    {
        for (int i = 0 ; i < no_of_records ; i ++)
        { temp_pitch = pitch_ras[i]*3600 ;
            temp_yaw = yaw_ras[i]*3600 ;            
            transformRPYtoDXDYDTHETA_VIS (roll_ras[i] , temp_pitch , temp_yaw , delta_x[i] , delta_y[i] , delta_theta[i]) ;

        }

    }
    else
    {
        LOG (INFO) << "***Invalid Channel***"  ;
        return (EXIT_FAILURE) ;
    }

    double ctheta , stheta ;
    for (int i = 0 ; i < nrows ; i ++)
    {
           
        flag = FALSE ;
        time_frame = t[i] ;
            
          for (int index = 1 ; index < no_of_records ; index ++)
        {
           //LOG(INFO)<<setprecision (20)<<time[index-1]<<" "<<time[index]<<" "<<time_frame<<" "<<no_of_records; 
            if (time_frame >= time_drifts[index-1] && time_frame < time_drifts[index])
            {
                flag = TRUE ;
                t1 = time_drifts[index-1] ;
                t2 = time_drifts[index ] ;
                theta1 = delta_theta[index-1] ;
                theta2 = delta_theta[index] ;
                x1 = delta_x[index-1] ;
                x2 = delta_x[index ] ;
                y1 = delta_y[index-1] ;
                y2 = delta_y[index ] ;
                new_delta_theta = (theta1 + (time_frame - t1)*(theta2 - theta1) / (t2 - t1)) *mult_fact ;
                new_delta_x = (x1 + (time_frame - t1)*(x2 - x1) / (t2 - t1)) * mult_fact ;
                new_delta_y = (y1 + (time_frame - t1)*(y2 - y1) / (t2 - t1)) * mult_fact ;
                Shifts_x.push_back(x1 + (time_frame - t1)*(x2 - x1) / (t2 - t1));
                Shifts_y.push_back(y1 + (time_frame - t1)*(y2 - y1) / (t2 - t1));
                Shifts_theta.push_back(theta1 + (time_frame - t1)*(theta2 - theta1) / (t2 - t1));
                frm_no_vect.push_back(frmno[i]);
                multphotonFlagOut.push_back(multphotonFlagIn[i]);
                break ;
            }
        }
        //check whether  matching time found or not
        if (flag == FALSE)
        {
            track_rown.push_back (i + 1) ;
            new_delta_theta = 0 ;
            new_delta_x = 0 ;
            new_delta_y = 0 ;
        }
      //  new_delta_x=0.0;new_delta_y=0.0f;new_delta_theta=0.0f;
        ctheta = cos (-1.0 * new_delta_theta * M_PI / 180) ;
        stheta = sin (-1.0 * new_delta_theta * M_PI / 180) ;
        i1 = xf[i] - xsize / 2 ;
        j1 = yf[i] - ysize / 2 ;
        //calculating the new values for Xand Y based on the deltas
        if ((((i1 * ctheta) - (j1 * stheta)) + xsize / 2 - new_delta_x)>=0 && (((i1 * ctheta) - (j1 * stheta)) + xsize / 2 - new_delta_x) < xsize && (((i1 * stheta) + (j1 * ctheta)) + ysize / 2 - new_delta_y)>=0 &&  (((i1 * stheta) + (j1 * ctheta)) + ysize / 2 - new_delta_y) < xsize)
        {
            xi_fi[i] = ((i1 * ctheta) - (j1 * stheta)) + xsize / 2 -new_delta_x ; //new index x
            yi_fi[i] = ((i1 * stheta) + (j1 * ctheta)) + ysize / 2 -new_delta_y ; //new index y
        }
    }
   

    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::findStar_algo1 (float *inputArray) //algorithm for finding the peaks
{

    Fx.clear () ;
    Fy.clear () ;
    Fval.clear () ;
    Rx.clear () ;
    Ry.clear () ;
    Rval.clear () ;
    Cx.clear () ;
    Cy.clear () ;
    Ci.clear () ;
    int r , c ;
    float *temp_array ;
    vector<float> array_temp ;
    if (datainfo.getModeFlag () == PC)
    {


        array_temp.clear () ;

        for (int i = 0 ; i < subDivision_size * subDivision_size ; i ++)
        {

            if (inputArray[i] != 0.0f)
            {
                array_temp.push_back (inputArray[i]) ;

            }
        }
        temp_array = new float[array_temp.size ()] ;
        for (int in = 0 ; in < array_temp.size () ; in ++)
        {
            temp_array[in] = array_temp[in] ;

        }
    }

label:
    Fval.clear () ;
    Fx.clear () ;
    Fy.clear () ;
    Rx.clear () ;
    Ry.clear () ;
    Rval.clear () ;


    // LOG(INFO)<<sd_mul_factor<<endl;exit(1);
    if (sd_mul_factor < 0)
    {
        LOG (ERROR) << "***SD_MULTI_FACTOR is <0***"  ;
	LOG(ERROR)<<"CRASH: PIL CHECK : SIGMA MULTIPLIER =< 0) (uvtLevel2PC.cpp)";
        return (EXIT_FAILURE) ;
    }
    double thr = 0 ;
    float mean_Ofimage , sd_temp ;
    ;
    if (datainfo.getModeFlag () == PC)
    {

        sd_temp = getSD (temp_array , array_temp.size ()) ;
        mean_Ofimage = getmean (temp_array , array_temp.size ()) ;
        //    thr =  sd_temp* sd_mul_factor ;
        thr = mean_Ofimage + sd_temp* sd_mul_factor ;
    }
    else
    {
        sd_temp = getSD (inputArray , subDivision_size * subDivision_size) ;
        mean_Ofimage = getmean (inputArray , subDivision_size * subDivision_size) ;
        //    thr =  sd_temp* sd_mul_factor ;
        thr = mean_Ofimage + sd_temp* sd_mul_factor ;
        ;
    }
   
    LOG (ERROR) << "Threshold for first cut peaks is   " << mean_Ofimage << " + " << sd_temp << " X " << sd_mul_factor << " = " << thr ;
  //  thr=0;//for time being.
    for (int i = 0 ; i < subDivision_size * subDivision_size ; i ++)
    {
        r = (i / subDivision_size) ;
        c = (i % subDivision_size) ;

        if (inputArray[i] > thr)
        {
            Fval.push_back (inputArray[i]) ;
            Fx.push_back (c+1) ; //x is for column
            Fy.push_back (r+1) ; //y is for row
        }

    }

    LOG (INFO) << "SIGMA  Factor::" << sd_mul_factor ;
    LOG (INFO) << " Size of First cut Peaks  " << Fy.size ();

    if (Fy.size () < minimum_No_of_Stars)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor < 0)
        {
            LOG (ERROR) << sd_mul_factor << " less than 0!!!! " ;
	    LOG(ERROR)<<"CRASH NO STAR FOUND DUE TO SIGMA MULTIPLIER BEING =< 0) (uvtLevel2PC.cpp)";
            return (EXIT_FAILURE) ;
        }
        goto label ;

    }

 

   
    //if winsize is even, make it odd
    if (refine_Winsize % 2 == 0)
        refine_Winsize = refine_Winsize - 1 ;

    LOG (INFO) << "Using window size : " << refine_Winsize << " for refining peaks " ;

    //refined peaks
    vector<int> Tx , Ty ;
    vector<float> Tval ;
   
    Tx.reserve (Fx.size ());
    Ty.reserve (Fy.size ());
    Tval.reserve (Fval.size ());
    
   
    Tx = Fx ;
    Ty = Fy ;
    Tval = Fval ;
    Star star1 ;
    star_track.clear () ;
    star_track.reserve (Tx.size ());
    for (int i = 0 ; i < Tx.size () ; i ++)
    {
        star1.x = Tx[i] ;
        star1.y = Ty[i] ;
        star1.intensity = Tval[i] ;
        star_track.push_back (star1) ;
        //LOG(INFO)<<Tx[i]<<" "<<Ty[i]<<" "<<Tval[i]<<endl;
    }
    

    sort (star_track.begin () , star_track.end () , compare1) ;
    //LOG(INFO)<<"sorting  finished";
    /*refining peaks logic
    refined Window size is for the refined  peaks.
   Refined peaks are found by  making window around each of the star(i.e first cut peaks)  and  finding brightest star among that window.*/
    int start_r , end_r , start_c , end_c ;
    //to be removed 
    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;
    //  bool flag_unique=FALSE;
    
    // vector<Star> ::iterator itr =star_track.begin ();
    for (int i = 0 ; i < star_track.size () ; i ++)
    {
        start_r = star_track[i].y - refine_Winsize / 2 ;
        end_r = star_track[i].y + refine_Winsize / 2 ;
        start_c = star_track[i].x - refine_Winsize / 2 ;
        end_c = star_track[i].x + refine_Winsize / 2 ;
        if (start_r < 0) start_r = 0 ;
        if (end_r >= subDivision_size) end_r = subDivision_size - 1 ;
        if (start_c < 0) start_c = 0 ;
        if (end_c >= subDivision_size) end_c = subDivision_size - 1 ;
        for (int fcpeak = i + 1 ; fcpeak < star_track.size () ; fcpeak ++)
        {
            if (star_track[fcpeak].x > start_c && star_track[fcpeak].x < end_c && star_track[fcpeak].y > start_r && star_track[fcpeak].y < end_r)
            {

                star_track.erase (star_track.begin () + fcpeak) ;
                fcpeak -- ;
            }
        }
        Tx.push_back (star_track[i].x) ;
        Ty.push_back (star_track[i].y) ;
        Tval.push_back (star_track[i].intensity) ;
    }
  

    /*--------------Refining peaks completed----------------*/

    float *arr_refine = new float[subDivision_size * subDivision_size] ; //to store refined peaks

    for (int i = 0 ; i < subDivision_size * subDivision_size ; i ++)
        arr_refine[i] = 0 ;

    if (subDivision_size == 0 || subDivision_size == 0)
    {
        LOG (ERROR) << "***Divide by Zero***" << endl ;
        return (EXIT_FAILURE) ;
    }
    for (int i = 0 ; i < Ty.size () ; i ++){
        if((Ty[i] * subDivision_size + Tx[i])>0 && (Ty[i] * subDivision_size + Tx[i])<subDivision_size*subDivision_size)
        arr_refine[Ty[i] * subDivision_size + Tx[i]] = Tval[i] ; //overwriting the same place..
    }
        

    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;

    for (int i = 0 ; i < subDivision_size * subDivision_size ; i ++)
    {
        if (arr_refine[i] != 0)
        {
            Rx.push_back ((i % subDivision_size)) ;
            Ry.push_back ((i / subDivision_size)) ;
            Rval.push_back (arr_refine[i]) ;
        }
    }
    LOG (INFO) << "Number of final peaks is " << Rval.size () ;


    if (Ry.size () < minimum_No_of_Stars)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor < 0)
        {
            LOG (ERROR) << sd_mul_factor << " less than 0!!!! " ;
	    LOG(ERROR)<<"CRASH NO STAR FOUND DUE TO SIGMA MULTIPLIER BEING =< 0) (uvtLevel2PC.cpp)";
            return (EXIT_FAILURE) ;
        }
        delete[]  arr_refine;
        goto label ;
    }
    /**method for  find Centroid within the Stars**/
   doCentroiding (Rx , Ry , centroid_Winsize , inputArray , subDivision_size , subDivision_size) ;
    
    delete[] arr_refine ;
    return (EXIT_SUCCESS) ;
}


void uvtLevel2PC::doCentroiding (vector<int> &X , vector<int> &Y , int centroidwindow , float *arr , int h , int w)
{

    Cx.clear () ;
    Cy.clear () ;
    Ci.clear () ;

    float x , y , val = 0 ;
    //   LOG(INFO) << "THe Centroid Window" << centroidwindow << endl;
    int num = centroidwindow*centroidwindow ;
    double sum_x = 0 , sum_y = 0 , sum = 0 ;


    //  LOG(INFO) << "The X.size is " << X.size() << endl;
    for (int i = 0 ; i < X.size () ; i ++)
    {
        sum_x = 0 ;
        sum_y = 0 , sum = 0 ;
        for (int j = - 1 * (centroidwindow / 2) ; j <= (centroidwindow / 2) ; j ++)
        {
            for (int k = - 1 * (centroidwindow / 2) ; k <= (centroidwindow / 2) ; k ++)
            {
                x = X[i] ;
                y = Y[i] ;
              //  if(j>0 && j<subDivision_size && k>0 && k<subDivision_size )
                val = arr[(Y[i] + j) * w + (X[i] + k)] ;
             //   else
            //        continue;
                if (val == INVALID_PIX_VALUE)
                {
                    continue ;
                }
                sum_x = sum_x + (x + k) * val ;
                sum_y = sum_y + (y + j) * val ;
                sum = sum + val ;
            }
        }
        if (sum <= 0)
        {
            LOG (INFO)<<  "Sum of intensites for (" << X[i] << " , " << Y[i] << ")  is <=0"<<val ;
            LOG (INFO) <<"Divide by zero error" ;
            exit (EXIT_FAILURE) ;
        }
        Cx.push_back ((float) sum_x / (float) sum) ;
        Cy.push_back ((float) sum_y / (float) sum) ;
        Ci.push_back ((float) sum) ;

    }
    float temp , tx , ty ;
    for (int i = 0 ; i < Ci.size () ; i ++)
    {
       
        //LOG(INFO)<<endl<<i;
        for (int j = Ci.size () - 1 ; j > i ; j --)
        {
            if (Ci[j - 1] < Ci[j])
            {
                swap1 (Ci[j] , Ci[j - 1]) ;
                swap1 (Cx[j] , Cx[j - 1]) ;
                swap1 (Cy[j] , Cy[j - 1]) ;
            }
        }
    }
//    Cx.clear();Cy.clear ();Ci.clear ();
//    for(int i=0;i<Rx.size ();i++)
//    {
//    Cx.push_back (Rx[i]);
//    Cy.push_back (Ry[i]);
//    Ci.push_back (Rval[i]);
//    }
}


int uvtLevel2PC::matchStars (int numrowsFirstfile , int numrowsSecfile , float divFact ,
        float *xlocFirst , float *ylocFirst , float *xlocSec , float *ylocSec ,float * int1,float *int2, vector<float> &matchPixelXone ,
        vector<float> &matchPixelYone , vector<float> &matchPixelXtwo , vector<float> &matchPixelYtwo ,vector<float> &int_one,vector<float> &int_two,
        vector<float> &matchPixelDiffX , vector<float> &matchPixelDiffY)
{
    int cnt = 0 ;
    float temp_x1 = 0 , temp_x2 = 0 , temp_y1 = 0 , temp_y2 = 0 ;
    
    //    LOG(INFO) << numrowsFirstfile << " " << numrowsSecfile << endl ;
    for (long int index = 0 ; index < numrowsFirstfile ; index ++) //loop for finding the similer x-cordinates and y-cordinates  from the both frames
    {
        //int pairing_flag = 0 ;
        temp_x1 = xlocFirst[index] / divFact ;
        temp_y1 = ylocFirst[index] / divFact ;
        //  ints1 = Intensity[index] ;

        for (long int i = 0 ; i < numrowsSecfile ; i ++)
        {
            float diff_x = 0.0 , diff_y = 0.0 ;
            temp_x2 = xlocSec[i] / divFact ;
            temp_y2 = ylocSec[i] / divFact ;
            diff_x = temp_x1 - temp_x2 ;
            diff_y = temp_y1 - temp_y2 ;

            if ((diff_x * diff_x) < diff_Dist && (diff_y * diff_y) < diff_Dist) // finding the similar points
            {
                matchPixelXone.push_back (temp_x1) ;
                matchPixelYone.push_back (temp_y1) ;
                matchPixelXtwo.push_back (temp_x2) ;
                matchPixelYtwo.push_back (temp_y2) ;
                matchPixelDiffX.push_back ((- 1.0) * diff_x) ;
                matchPixelDiffY.push_back ((- 1.0) * diff_y) ;
                int_one.push_back(int1[(int)index]);
                int_two.push_back(int2[(int)i]);
              //  ofptr1<<temp_x1*div_fact<<setw(20)<<temp_y1*div_fact<<setw(20)<<ints1<<setw(20)<<temp_x2*div_fact<<setw(20)<<temp_y2*div_fact<<setw(20)<<ints2<<endl;
                cnt ++ ;
                break ;

            }
        }
    }
    return cnt ;
}


int uvtLevel2PC::findShiftsNtheta (int totalelements , vector<float> &Xone , vector<float> &Yone ,vector<float> &int1, vector<float> &Xtwo , vector<float> &Ytwo ,
       vector<float> &int2, vector<float> &DiffOfX , vector<float> &DiffOfY , bool flag_theta_computation,double &Xdx , double &Ydy , double &Theta)
{
    
     if(flag_theta_computation==0)
    {
        LOG(INFO)<<"NO theta value will be calculated";
        vector<double> diff_x_cumm,diff_y_cumm;
          for (int t = 0 ; t < totalelements  ; t++)
          {
             
              diff_x_cumm.push_back (Xtwo[t] - Xone[t]);
              diff_y_cumm.push_back (Ytwo[t] - Yone[t]);
          }
        Xdx=getmean (diff_x_cumm.data (),diff_x_cumm.size ());
        Ydy=getmean (diff_y_cumm.data (),diff_y_cumm.size ());
        Theta=0.0f;
    }
     else   if (shift_N_Rotate_algo == 1)
    {
        spMatrix B ((totalelements) * 2 , 1) ;
        spMatrix A ((totalelements) * 2 , 3) ;
        spMatrix X (3 , 1) ;
        int temp = 0 ;

        for (int t = 0 ; t < totalelements * 2 ; t = t + 2)
        {
            B (t , 0) = (Xtwo[temp] - Xone[temp]) ; //+IMAGE_ARRAYSIZE*0.5;
            B (t + 1 , 0) = (Ytwo[temp] - Yone[temp]) ; //+IMAGE_ARRAYSIZE*0.5;

            A (t , 0) = - 1.0 * (Yone[temp] - (subDivision_size / mult_fact) / 2) ;
            A (t , 1) = 1 ;
            A (t , 2) = 0 ;

            A (t + 1 , 0) = (Xone[temp] - (subDivision_size / mult_fact) / 2) ;
            A (t + 1 , 1) = 0 ;
            A (t + 1 , 2) = 1 ;

            temp ++ ;
        }

        X.ApplyLeastSquare (A , B) ;

        Xdx = X (1 , 0) ;
        Ydy = X (2 , 0) ;
        Theta = X (0 , 0) ;
        //X (0 , 0) = X (0 , 0) ;
        // ofptr << p << setw (15) << X (0 , 0)*180 / M_PI << setw (25) << X (1 , 0) << setw (25) << X (2 , 0) << endl ;
    }
    else if (shift_N_Rotate_algo == 2)
    {
        spMatrix A (totalelements * 2 , 4) ;
        spMatrix B (totalelements * 2 , 1) ;
        spMatrix X (4 , 1) ;
        for (int aindex = 0 ; aindex < totalelements ; aindex ++)
        {
            A (2 * aindex , 0) = Xone[aindex] - (subDivision_size / mult_fact) * 0.5 ;
            A (2 * aindex , 1) = - 1.0 * (Yone[aindex] - (subDivision_size / mult_fact) * 0.5) ;
            A (2 * aindex , 2) = 1.0 ;
            A (2 * aindex , 3) = 0.0 ;
            A (2 * aindex + 1 , 0) = (Yone[aindex] - (subDivision_size / mult_fact) * 0.5) ;
            A (2 * aindex + 1 , 1) = Xone[aindex] - (subDivision_size / mult_fact) * 0.5 ;
            A (2 * aindex + 1 , 2) = 0.0 ;
            A (2 * aindex + 1 , 3) = 1.0 ;
            B (2 * aindex , 0) = Xtwo[aindex] - (subDivision_size / mult_fact) * 0.5 ;
            B (2 * aindex + 1 , 0) = Ytwo[aindex] - (subDivision_size / mult_fact) * 0.5 ;
        }
        X.ApplyLeastSquare (A , B) ;

        double theta = atan2 (X (1 , 0) , X (0 , 0)) ;

        // ofptr << theta * 180 / M_PI << setw (15) << X (0 , 0) << setw (25) << X (1 , 0) << endl ;
        Xdx = X (2 , 0) ;
        Ydy = X (3 , 0) ;
        Theta = theta ;
    }
    else if (shift_N_Rotate_algo == 3)
    {

        double a11 = 0.0 , a12 = 0.0 , a13 = 0.0 , a21 = 0.0 , a22 = 0.0 , a23 = 0.0 , a31 = 0.0 , a32 = 0.0 , a33 = 0.0 , b1 = 0.0 , b2 = 0.0 , b3 = 0.0 ;

        spMatrix A (3 , 3) ;
        spMatrix B (3 , 1) ;
        spMatrix X (3 , 1) ;
         float sum_int_ref=0;
         
        for(int i=0;i<totalelements;i++)
        {
         sum_int_ref=sum_int_ref+int1[i];      
         
        }
        for (int k = 0 ; k < totalelements ; k ++)
        {
            /* weight used for a star*/
            double w =int1[k]/sum_int_ref;
            //     w=pow( (25.0*photonbk_per_pixel+starmag[k]), 2.0 )/
            //        (  (25.0*photonbk_per_pixel)+(0.09*starmag[k])  );


            /* row 1 */
            double y_mod = Ytwo[k]-((subDivision_size / mult_fact) * 0.5) ;
            double x_mod = Xtwo[k]-((subDivision_size / mult_fact) * 0.5) ;

            y_mod = - 1.0 * y_mod ;
            a11 = a11 + 2.0 * w ;
            a12 = 0.0 ;
            a13 = a13 + (2.0 * w * (y_mod)) ;
            b1 = b1 + 2.0 * (w * DiffOfX[k]) ;

            /* row 2*/
            a21 = 0.0 ;
            a22 = a22 + 2.0 * w ;
            a23 = a23 + (2.0 * w * (x_mod)) ;
            b2 = b2 + 2.0 * (w * DiffOfY[k]) ;

            /* row 3 */
            a31 = a31 + (2.0 * w * (y_mod)) ;
            a32 = a32 + (2.0 * w * (x_mod)) ;
            a33 = a33 + (2.0 * w * (pow (x_mod , 2.0) + pow (y_mod , 2.0))) ;
            b3 = b3 + (2.0 * w * ((x_mod) * DiffOfY[k] + (y_mod) * DiffOfX[k])) ;

        }

        A (0 , 0) = a11 ;
        A (0 , 1) = a12 ;
        A (0 , 2) = a13 ;
        A (1 , 0) = a21 ;
        A (1 , 1) = a22 ;
        A (1 , 2) = a23 ;
        A (2 , 0) = a31 ;
        A (2 , 1) = a32 ;
        A (2 , 2) = a33 ;
        B (0 , 0) = b1 ;
        B (1 , 0) = b2 ;
        B (2 , 0) = b3 ;

        X.ApplyLeastSquare (A , B) ;

        // ofptr << X (2 , 0)*180 / M_PI << setw (15) << X (0 , 0) << setw (25) << X (1 , 0) << endl ;
        Xdx = X (0 , 0) ;
        Ydy = X (1 , 0) ;
        Theta = X (2 , 0) ;
    }
    return (EXIT_SUCCESS) ;
}


int uvtLevel2PC::ApplySubSampling (float* inputarray , int in_xsize , int in_ysize , float* outputarray , int out_xsize , int out_ysize)
{

    if (in_xsize % out_xsize != 0)
    {
        LOG (ERROR) << "Can't sub sampling with  this input array ,input Array size is not matching" ;
        return (EXIT_FAILURE) ;
    }
    float devide_term_x = in_xsize / out_xsize ;
    float devide_term_y = in_ysize / out_ysize ;
    int cnt_win = 0 ;
    float sum_win = 0.0f ;
    int index_finalArray = 0 ;
    for (int temp_x = 0 ; temp_x < in_xsize ; temp_x = temp_x + devide_term_x)
    {
        for (int temp_y = 0 ; temp_y < in_ysize ; temp_y = temp_y + devide_term_y)
        {
            if (index_finalArray > out_xsize * out_ysize)
            {
                LOG (ERROR) << "Array is out of bound, EXEED to " << out_xsize << " !!!" ;
                return (EXIT_FAILURE) ;
            }
            cnt_win = 0 ;
            sum_win = 0.0f ;
            for (int i = temp_y ; i < temp_y + devide_term_y ; i ++)
            {
                for (int j = temp_x ; j < temp_x + devide_term_x ; j ++)
                {
                    if (inputarray[j * in_xsize + i] != INVALID_PIX_VALUE)
                    {
                        sum_win = sum_win + inputarray[j * in_xsize + i] ;
                        cnt_win ++ ;
                    }

                }
            }
            if (cnt_win != 0)
            {
                outputarray[index_finalArray ++] = (float) (sum_win/ cnt_win) ;
            }
            else
            {
                outputarray[index_finalArray ++] = INVALID_PIX_VALUE ;
            }


        }

    }

    return (EXIT_SUCCESS) ;

}




int  uvtLevel2PC:: removeRecords(vector<float>  &Xone , vector<float> &Yone ,vector<float> &Xtwo , vector<float> &Ytwo,vector<float> &DiffOfX , vector<float> &DiffOfY ,float *ints1,float*ints2,
   vector<float> &newXone,vector<float> &newYone,vector<float> &newXtwo,vector<float> &newYtwo,vector<float> &newDiffX,vector<float> &newDiffY,vector<float> &new_one_ints,vector<float> &new_two_ints)
 {
     float mean_diff_X=0.0f,mean_diff_Y=0.0f;
     float sd_diff_X=0.0f,sd_diff_Y=0.0f;
     mean_diff_X=getmean (DiffOfX.data (),DiffOfX.size ());
     sd_diff_X=getSD (DiffOfX.data (),DiffOfX.size ());
      mean_diff_Y=getmean (DiffOfY.data (),DiffOfY.size ());
     sd_diff_Y=getSD (DiffOfY.data (),DiffOfY.size ());
     newDiffX.clear (),newDiffY.clear (),newXone.clear ();newXtwo.clear ();newYone.clear ();newYtwo.clear ();new_one_ints.clear (),new_two_ints.clear ();
        vector<float> Distance_vect;
     for (int i=0;i<DiffOfX.size ();i++)
     {
             Distance_vect.push_back (sqrt(DiffOfX[i]*DiffOfX[i]+DiffOfY[i]*DiffOfY[i]));
         
     }
        float mean_dist,sd_dist;
        mean_dist=getmean (Distance_vect.data (),Distance_vect.size ());
        sd_dist=getSD (Distance_vect.data (),Distance_vect.size ());
//     for (int i=0;i<DiffOfX.size ();i++)    
//                {
//                    if (DiffOfX[i]<=mean_diff_X+sd_diff_X && DiffOfY[i]<=mean_diff_Y+sd_diff_Y)
//                    {
//                        newXone.push_back (Xone[i]);
//                        newYone.push_back (Yone[i]);
//                        newXtwo.push_back (Xtwo[i]);
//                        newYtwo.push_back (Ytwo[i]);
//                        newDiffX.push_back (DiffOfX[i]);
//                        newDiffY.push_back (DiffOfY[i]);
//                        new_one_ints.push_back (ints1[i]);
//                        new_two_ints.push_back (ints2[i]);
//                    }
//                    
//                }
         for (int i=0;i<DiffOfX.size ();i++)    
         {          
             //LOG(INFO)<<mean_dist<<" "<<sd_dist<<" "<<Distance_vect[i];
                    if (Distance_vect[i]<=mean_dist+3*sd_dist  && Distance_vect[i]>-1.0*(mean_dist+3*sd_dist))
                    {
                        newXone.push_back (Xone[i]);
                        newYone.push_back (Yone[i]);
                        newXtwo.push_back (Xtwo[i]);
                        newYtwo.push_back (Ytwo[i]);
                        newDiffX.push_back (DiffOfX[i]);
                        newDiffY.push_back (DiffOfY[i]);
                        new_one_ints.push_back (ints1[i]);
                        new_two_ints.push_back (ints2[i]);
                    }
                    
         }
     
       
        return (EXIT_SUCCESS);
     
 }
     
    
 int uvtLevel2PC::getHistory (vector<string> &vhistory)
{
  //  char *user = getlogin () ;
    int cnt = 0 ;
    char validgtiflag_str[FLEN_FILENAME] ;
 //   string str = "Module run by " + (string) user ;
    char temp[PIL_LINESIZE];
   // vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Input Level1 tar  file = " + (string) level1indir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Caldb used= " + (string) caldbindir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Module Output directory = " + (string)level2outdir ) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " channel used= " + (string)channel) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Parity check flag for photon event= " + (string)convertIntToStr(parity_flag)) ;
vhistory.push_back ((string) getSerialNo (cnt) + " CRC flag = " + (string)convertIntToStr(crc_flag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Flag for CRC failure action = " + (string) convertIntToStr (dropframe)) ;
    //vhistory.push_back ((string) getSerialNo (cnt) + " Dark frame subtraction to be done or not = " + (string)convertIntToStr (darkframe_flag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " GTI flag  = " + (string)convertIntToStr(gti_flag)) ;
    
    vhistory.push_back ((string) getSerialNo (cnt) + " Threshold for multi-photon events " + (string)convertIntToStr(thr_multiph)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " CR threshold parameter-1-“_N_” in threshold=AVG+N*sqrt(avg_events)+ST/(sqrt(avg_events)) to identify Cosmic Ray affected frame = " + (string)convertIntToStr( first_mult_factor_CR)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " CR threshold parameter-2-“_ST_” in threshold=AVG+N*sqrt(avg_events)+ST/(sqrt(avg_events)) to identify Cosmic Ray affected frame = " + (string)convertFloatToStr(second_mult_Factor_CR)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " QEMCP Correction Flag= " + (string)convertIntToStr(qemcpflag )) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Centroid Dark Correction Flag= " + (string)convertIntToStr(centCorrflag )) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Centroid Bias Flag = " + (string)convertIntToStr( centBiasflag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Detector Distortion Flag = " + (string)convertIntToStr( DetectDistflag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Optical Distortion Flag = " + (string)convertIntToStr( OpticDistflag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Frame Integration Flag:MULTI/SINGLE = " + (string)convertIntToStr(fi_flag)) ;
    if(fi_flag){
            vhistory.push_back ((string) getSerialNo (cnt) + " Number of initial frames to be discarded for MULTI case = " + (string) convertIntToStr(nFrameDiscard_fi)) ;
            vhistory.push_back ((string) getSerialNo (cnt) + " Number of frames to be combined in MULTI case = " + (string) convertIntToStr(nFrameIntegrate_fi)) ;
    }
    vhistory.push_back ((string) getSerialNo (cnt) + " Frame size after frame integration = " + (string) convertIntToStr(IMG_DIM_FI)) ;
            vhistory.push_back ((string) getSerialNo (cnt) + " Star finding algorithm = " + (string) convertIntToStr(star_detect_algo_flag)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Algorithm to find  offsets of current image w.r.t  reference image = " + (string) convertIntToStr(shift_N_Rotate_algo)) ;
           // if(star_detect_algo_flag==1 || star_detect_algo_flag==3 || star_detect_algo_flag==4)
          //  {
                vhistory.push_back ((string) getSerialNo (cnt) + " Parameter in fullframe astrometry to identify outliers : multiplier to sigma = " + (string) convertFloatToStr (sd_multi_factor_default)) ;
                vhistory.push_back ((string) getSerialNo (cnt) + " Minimum targeted stars requirement = " + (string)convertIntToStr(minimum_No_of_Stars)) ;
          //  }
        //   if(flag_Roll_Angleinvalid==TRUE){
//vhistory.push_back ((string) getSerialNo (cnt) + " Roll Angle valid(1)/invalid(0)  = " + (string)convertIntToStr(flag_Roll_Angleinvalid)) ;
//}
           vhistory.push_back ((string) getSerialNo (cnt) + " Shift and rotation algorithm  = " + (string)convertIntToStr(shift_N_Rotate_algo)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Neighbourhood criteron for identifying stars = " + (string)convertIntToStr( refine_Winsize)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Box Size to compute Centroid for  detected stars = " + (string)convertIntToStr( centroid_Winsize)) ;
             
             vhistory.push_back ((string) getSerialNo (cnt) + " Relative Aspect Series file path  =" + (string) rasfile) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Search algorithm for  star catalogue database =" + (string) convertIntToStr(search_algo_ctlg)) ;
             vhistory.push_back((string) getSerialNo (cnt)+ "Search radius for star catalogue="+(string)rad_search);
//             if(search_algo_ctlg==1 || search_algo_ctlg==3 || search_algo_ctlg==5){
//            vhistory.push_back ((string) getSerialNo (cnt) + " Lengh of rectangle search   = " + (string)(len_a)) ;    
//             vhistory.push_back ((string) getSerialNo (cnt) + " Width of rectangle search  = " + (string) (len_b)) ;  
//             }
             vhistory.push_back ((string) getSerialNo (cnt) + "Registration and averaging image size    = " + (string) convertIntToStr(FINALFRAMESIZE_REGAVG)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + "Flag theta ON/OFF    = " +  (string)convertIntToStr(flag_thetaComp )) ;
             //vhistory.push_back ((string) getSerialNo (cnt) + " frame to be discarded in reference frame calculation = " + (string)convertIntToStr(frames_toDiscard)) ;
            //vhistory.push_back ((string) getSerialNo (cnt) + " Average Factor = " + (string)convertIntToStr( nFrameToAverage)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Minimum search distance for star match in succesive images of FrameIntegration stage. " + (string)convertFloatToStr ( diff_Dist)) ;
             //vhistory.push_back ((string) getSerialNo (cnt) + " Frequency domain filtering flag  " + (string) convertIntToStr(freqDomainFilter_Flag)) ;
             //if(freqDomainFilter_Flag==1)
             //vhistory.push_back ((string) getSerialNo (cnt) + " Type Filtering used   " + (string) convertIntToStr(type_Filtering)) ;
              vhistory.push_back ((string) getSerialNo (cnt) + " Catalogue database path = " + (string)databasename) ;
               vhistory.push_back ((string) getSerialNo (cnt) + "Search radius (in degrees) for star match with star catalogue. = " + (string)rad_search) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Output tar location = " + (string)outtarpath) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for unitConversion= " + (string)convertIntToStr( wtd_uc)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for flat field correction= " + (string)convertIntToStr( wtd_ff)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for mask Bad pixel correction= " + (string)convertIntToStr(wtd_bp)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for Pix padding= " + (string)convertIntToStr( wtd_pp)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for  subDivision= " + (string)convertIntToStr( wtd_sd)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for cosmic ray correction= " + (string)convertIntToStr( wtd_cr)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for QEMCP correction= " + (string)convertIntToStr( wtd_qemcp)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for detector distortion=  " + (string)convertIntToStr( wtd_dd)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for optical distortion= " + (string)convertIntToStr( wtd_od)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for Shift and Rotate= " + (string)convertIntToStr( wtd_snr)) ;
		vhistory.push_back ((string) getSerialNo (cnt) + " Equation for calculation of MJD -> MJD = BZERO_MJD + BSCALE_MJD * LOCAL_TIMEV " ) ;
             
             
             if (clobber == YES)
        vhistory.push_back ((string) getSerialNo (cnt) + " clobber = yes") ;
    else
        vhistory.push_back ((string) getSerialNo (cnt) + "clobber = no") ;
    if (history == YES)
        vhistory.push_back ((string) getSerialNo (cnt) + " history = yes") ;
    else
        vhistory.push_back ((string) getSerialNo (cnt) + "  P9 history = no") ;
    vhistory.push_back ("Parameter List END") ;


    return (EXIT_SUCCESS) ;
}
