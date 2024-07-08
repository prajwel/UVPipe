/* 
 * File:   uvtRelativeAspectPC.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <stdlib.h>

#include "uvtRelativeAspectPC.h"
#include "uvtComputeDrift.h"
#include <uvtQEMCPCorr.h>
#include<sstream>

//#include "uvtUnitConversion.h"
//#include "caldb_Handler.h"
//#include "uvtFilterBadpix.h"
//#include "uvtPixPadding.h"
//#include "uvtSubDivision.h"
//#include<uvtComputeDrift.h>
//#include "uvtDetectStar.h"
#include<Directory.h>
#include<DataInfo.h>
#include<DataIngest.h>
#include <algorithm>
//#include<uvtQEMCPCorr.h>
#include<set>
//#include<uvtDetectDistCorr.h>
#include <vector>
#include<memory.h>
#include<spMatrix.h>
#include<glog/logging.h>
#include <bits/basic_string.h>
#include<uvtDetectStar.h>
#include<uvtComputeJitter.h>
#include<uvtComputeThermal.h>
#include<uvtRelAspCal.h>
//#include<uvtComputeDrift.h>
#define  NBHD_RADIUS 5
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr   //used in doFft() function
#define moduleoutdir_unit "uvtUnitConvertion"
#define moduleoutdir_badpix "uvtMaskBadPix"
#define moduleoutdir_fltfield "uvtFlatFieldCorr"
#define moduleoutdir_qe "uvtQEMCPCorr"
#define moduleoutdir_pixpad "uvtPixPadding"
#define moduleoutdir_cosmicray "uvtCosmicRayCorrection"
#define moduleoutdir_FrameIntegration "uvtFrameIntegration"
#define moduleoutdir_findstarcentroid "uvtDetectStar"
#define moduleoutdir_detectordistortion "uvtDetectDistCorr"
#define moduleoutdir_opticaldistortion "uvtOpticDistCorr"
#define moduleoutdir_driftExercise "uvtComputeDrift"
#define moduleoutdir_refFrameCal "uvtRefFrameCal"
#define moduleoutdir_centCorr "uvtCentroidCorr"
#define moduleoutdir_centBias "uvtCentroidBias"
//#define IMAGE_ARRAYSIZE  4800
#define ter 1

using namespace std ;


uvtRelativeAspectPC::uvtRelativeAspectPC () {
    tar_extracted_flag_PC=FALSE;
    //    sprintf (moduleoutdir_bp , "%s_%s" , moduleoutdir_badpix , VERSION) ;
    //    sprintf (moduleoutdir_uc , "%s_%s" , moduleoutdir_unit , VERSION) ;
    //    sprintf (moduleoutdir_ff , "%s_%s" , moduleoutdir_fltfield , VERSION) ;
    //    sprintf (moduleoutdir_pp , "%s_%s" , moduleoutdir_pixpad , VERSION) ;
    //    sprintf (moduleoutdir_sd , "%s_%s" , moduleoutdir_subdiv , VERSION) ;
    //    sprintf (moduleoutdir_cr , "%s_%s" , moduleoutdir_cosmicray , VERSION) ;
    //    sprintf (moduleoutdir_ac , "%s_%s" , moduleoutdir_AccEverytsec , VERSION) ;
    //    sprintf (moduleoutdir_sc , "%s_%s" , moduleoutdir_findstarcentroid , VERSION) ;
    //    sprintf (moduleoutdir_dd , "%s_%s" , moduleoutdir_detectordistortion , VERSION) ;
    //    sprintf (moduleoutdir_od , "%s_%s" , moduleoutdir_opticaldistortion , VERSION) ;
    //    sprintf (moduleoutdir_de , "%s_%s" , moduleoutdir_driftExercise , VERSION) ;
    //    sprintf (moduleoutdir_rfc , "%s_%s" , moduleoutdir_refFrameCal , VERSION) ;
}


uvtRelativeAspectPC::~ uvtRelativeAspectPC () {

    //    delete[] frame_Data_subdivided ;
    //    delete[] frame_ExpData_subdivided ;
    //    delete[] frame_Data_Padded ;
    //    delete[] frame_ExpoData_padded ;
    //    delete[] frame_fc_data ;
    //delete[] frame_Data;
    // delete[] flatfielddata;
}


int uvtRelativeAspectPC::readPILParameters (int argc , char** argv)
{
    int status = 0 ;
    char temp[PIL_LINESIZE] ;

    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG (INFO) << "***Error Initializing PIL***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetBool("zipFlag" , &zipFlag)))
    {
        LOG (INFO) << endl << "***Error reading caldb directory name***" ;
        return status ;
    }
    //if(paramfile_Varask==FALSE)
    //{ls -tlr
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
    //}
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
//        LOG (INFO) << endl << "***Error Reading att  flag value :" << padding_dim << "***" ;
//        return status ;
//    }

//    if (PIL_OK != (status = PILGetReal4 ("thrVal" , &Threshold_cr)))
//    {
//        LOG (INFO) << endl << status << "***Error reading threshold***" ;
//        return status ;
//    }
    if (PIL_OK != (status = PILGetInt ("cmpFrames" , &nCompareFrames)))
    {
        LOG (INFO) << endl << status << "***Error reading number of frames to be compared***" ;
        return status ;
    }


    if (PIL_OK != (status = PILGetBool ("GTI_FLAG" , &gti_flag)))
    {
        LOG (INFO) << endl << "***Error Reading GTI FLAG value:***" ;
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
    if (PIL_OK != (status = PILGetBool ("flatfieldFlag" , &flatfieldFlag)))
    {
        LOG (INFO) << endl << "***Error Reading flatfieldFlag:" << flatfieldFlag << "***" ;
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
    

    if (PIL_OK != (status = PILGetBool ("qemcpFlag" , &qemcpFlag)))
    {
        LOG (INFO) << endl << "***Error Reading qemcpFlag :" << qemcpFlag << "***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetInt ("algoFlag" , &star_detect_algo_flag)))
    {
        LOG (INFO) << endl << "***Error reading algo_flag ***" ;
        return status ;
    }
    if (star_detect_algo_flag == 1 || star_detect_algo_flag==3 || star_detect_algo_flag==4)
    {
        if (PIL_OK != (status = PILGetReal4 ("threshold" , &rms_mul_factor_default)))
        {
            LOG (INFO) << endl << "***Error reading threshold value for star detection***" ;
            return status ;
        }

        if (PIL_OK != (status = PILGetInt ("minimumTargetedstars" , &minimum_No_of_Stars)))
        {
            LOG (INFO) << endl << "***Error in reading the minimum number of stars ***" ;
            return status ;
        }
    }
    else if (star_detect_algo_flag == 2)
    {

        if (PIL_OK != (status = PILGetInt ("StarDetectionCentroidSqrSize" , &footprint_win_star_detect)))
        {
            LOG (INFO) << endl << "***Error reading  window size for JOE's algorithm ***" ;
            return status ;
        }
        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroidPriThr" , &primary_threshold_Val)))
        {
            LOG (INFO) << endl << "***Error reading primary threshold value for JOE's algorithm***" ;
            return status ;
        }
        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroidSecThr" , &secondary_threshold_Val)))
        {
            LOG (INFO) << endl << "***Error reading secondary threshold value for JOE's algorithm" ;
            return status ;
        }
    }
    else
    {
        LOG (INFO) << "***algo flag must be 1 or 2***" << endl ;
        return (EXIT_FAILURE) ;
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
     if (PIL_OK != (status = PILGetInt ("FrameIntDim" , &IMG_DIM_FI)))
    {
        LOG (INFO) << endl << "***Error reading frame size of frame integration ***" ;
        return status ;
    }    
   
    if(star_detect_algo_flag==1)
    {
     if (PIL_OK != (status = PILGetReal4 ("background_fact" , &backgrd_fact)))
        {
            LOG (INFO) << endl << "***Error reading background_fact ***" ;
            return status ;
        }
    }
    else{
        backgrd_fact=0;//can be ignored.
    }
    if (PIL_OK != (status = PILGetInt ("Window_nb" , &win_size)))
    {
        LOG (INFO) << endl << "***Error reading centroid  window size ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("GenMatchStarsFile_flag" , &match_stars_file_flag)))
    {
        LOG (INFO) << endl << "***Error reading GenMatchStarsFile_flag parameter***" ;
        return status ;
    }


    if (PIL_OK != (status = PILGetInt ("framesToBeDiscard" , &frames_toDiscard)))
    {
        LOG (INFO) << endl << "***Error reading the number of frames to discard in reference frame calculation***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("averageFactor" , &nFrameToAverage)))
    {
        LOG (INFO) << endl << "***Error reading average factor value in reference frame calculation ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetReal4 ("diffDist" , (float*) &pair_nbhd_distance)))
    {
        LOG (INFO) << endl << "***Error reading diffDist parameter***" ;
        return status ;
    }
//    if (PIL_OK != (status = PILGetReal4 ("error_per" , &err_per)))
//    {
//        LOG (INFO) << endl << "***Error reading error percent  ***" ;
//        return status ;
//    }
    if (PIL_OK != (status = PILGetInt ("freqDomainFilterFlag" , &FreqDomainFilterFlag)))
    {
        LOG (INFO) << endl << "***Error reading freqDomainFilter Value***" ;
        return status ;
    }
    if (FreqDomainFilterFlag == 0)
    {
        if (PIL_OK != (status = PILGetInt ("typeFiltering" , &type_Filtering)))
        {
            LOG (INFO) << endl << "***Error reading type Filtering Value***" ;
            return status ;
        }
        if (type_Filtering == 1)
        {

            if (PIL_OK != (status = PILGetBool ("fittingFlag" , &fitting_flag)))
            {
                LOG (INFO) << endl << "***Error reading fitting flag***" ;
                return status ;
            }
            if (fitting_flag)
            {
                if (PIL_OK != (status = PILGetInt ("orderPitch" , &orderpitch)))
                {
                    LOG (INFO) << endl << "***Error reading order of pitch***" ;
                    return status ;
                }
                if (PIL_OK != (status = PILGetInt ("orderYaw" , &orderyaw)))
                {
                    LOG (INFO) << endl << "***Error reading order of yaw***" ;
                    return status ;
                }
                if (PIL_OK != (status = PILGetInt ("orderRoll" , &orderroll)))
                {
                    LOG (INFO) << endl << "***Error reading order of roll***" ;
                    return status ;
                }
                if (PIL_OK != (status = PILGetReal ("deltaTime" , &poly_fit_interval)))
                {
                    LOG (INFO) << endl << "***Error reading the delta_Time***" ;
                    return status ;
                }

            }
        }
        else if (type_Filtering == 0)
        {
            if (PIL_OK != (status = PILGetReal ("freqValue" , (double*) &freqvalue)))
            {
                LOG (INFO) << endl << "***Error reading cut-off frequency value***" ;
                return status ;
            }

        }
    }
    else     if (FreqDomainFilterFlag == 1)
    {
        if (PIL_OK != (status = PILGetReal ("freqValue" , (double*) &freqvalue)))
        {
            LOG (INFO) << endl << "***Error reading freq_value***" ;
            return status ;
        }
    }
    else if (FreqDomainFilterFlag!=2)
    {
        LOG(ERROR)<<"***Value must be 0 0r 1 for FreqDomainFilterFlag ";
        return(EXIT_FAILURE);
    }

    if (PIL_OK != (status = PILGetInt ("shiftRotDetAlgoFlag" , &shift_rotation_algo)))
    {
        LOG (INFO) << endl << "***Error reading shiftRotDetAlgoFlag***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("framesDiscard" , &nFrameDiscard_fi)))
    {
        LOG (INFO) << endl << "***Error reading number of frames to be discarded in frame integration***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("framesCompute" , &nFrameIntegrate_fi)))
    {
        LOG (INFO) << endl << "***Error reading number of frames to integrate ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("flag_thetaComp" , &flag_thetaComp)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobber << "***" ;
        return status ;
    }
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
        LOG (ERROR) << endl << "***Error Reading Writing to Disk Flag For unit Conversion:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskbp" , &wtd_bp)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the bad pixel module:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskff" , &wtd_ff)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the FlatFieldCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskpp" , &wtd_pp)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentCorr" , &wtd_centCorr)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentBias" , &wtd_centBias)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }

    //    if (PIL_OK != (status = PILGetInt ("Write_todisksd" , &wtd_sd)))
    //    {
    //        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the subDivision:" << "***" ;
    //        return status ;
    //    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskcr" , &wtd_cr)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskqe" , &wtd_qemcp)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the QEMCP correction" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskfi" , &wtd_fi)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the frameintegration:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todisksc" ,  &wtd_fsc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the find Star centroid:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskdd" , &wtd_dd)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the Detector  Distortion:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskod" , &wtd_od)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the optical Distortion:" << "***" ;
        return status ;
    }
    PILClose (status) ;
    
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::uvtRelativeAspectPCProcess ()
{
    int status = 0 ; //status flag to check return status of functions
    string temp_str;
    temp_str.assign (level1indir) ;
    if(tar_extracted_flag_PC==FALSE)
    {
    level1indir = "" ;
    if(zipFlag==FALSE){
    status = extractTars (temp_str , level1indir,orbnum) ;
    if (status)
    {
        LOG (INFO) << "Error in extracting tar" ;
        LOG(ERROR)<<"CRASH: TAR EXTRACTION FAILED (uvtRelativeAspectPC.cpp)";
        return (EXIT_FAILURE) ;
    }
    }
     else{
           status = extractZip (temp_str , level1indir , orbnum) ;//extract level-1 data  tar file

            if (status) {
                LOG(INFO) << "Error in extracting tar";
                LOG(ERROR)<<"CRASH ZIP EXTRACTION FAILED (uvtRelativeAspectPC.cpp)";
                return (EXIT_FAILURE);
            } 
         
     }
    }
    if (! DirExists ((char *) level1indir.c_str ()))
    {
        LOG (ERROR) << endl << "Input level 1 directory not found " << endl ;
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
   else
   {
       LOG(INFO)<<"***Invalid channel***";
       return(EXIT_FAILURE);
   }
  
    if (dirobj.setup (level1indir , level2outdir , temp_channel_id))
    { //IM RA  needs VIS data only //setup done only if uvitV directory found
        LOG (ERROR) << endl << "***Error in directory set up***"  ;
        LOG (ERROR) << endl << "***This chain runs for NUV data only...check for uvitN directory in input*** "  ;
        LOG(ERROR)<<"CRASH L1 uvtN DIRECTORY ISSUE  (uvtRelativeAspectPC.cpp)";
        return (EXIT_FAILURE) ;
    }
    int numofsciencedatafiles = dirobj.sciencedatafile.size () ;
    if (numofsciencedatafiles <= 00)
    {
        LOG (ERROR) << endl << "No science data file found" << endl ;
        return (EXIT_FAILURE) ;
    }

    char obsmode[FLEN_VALUE] ;
    char moduleIndir[FLEN_FILENAME] ; //hold path for input directory for every module, will be updated after every module run
    char outputdir[FLEN_FILENAME] ;
    int division_fact ;

    //    uvtDetectStar obj_sc ;
    //    uvtQEMCPCorr qe_obj;
    if (subdivisionFlag == 0)
    {
        subDivision_size = padding_dim ;
    }
    frame_Data = new float[IMG_DIM_FI * IMG_DIM_FI] ;
    frame_fc_data = new float[IMG_DIM_FI * IMG_DIM_FI] ;
    frame_Data_Padded = new float[padding_dim * padding_dim] ;
    frame_ExpoData_padded = new float[padding_dim * padding_dim] ;
   

    char drift_outDir[FLEN_FILENAME] ;
    char jitter_outDir[FLEN_FILENAME] ;
     bool flag_rfc = FALSE ;
          int nFrame_Avg_ref=nFrameToAverage;
          fitsfile *fptr;
          vector<string> L1keywords;
         
    for (int dataindex = 0 ; dataindex < numofsciencedatafiles ; dataindex ++)
    {
        status=0;
         fits_open_file(&fptr , dirobj.sciencedatafile[dataindex].c_str () , READONLY , &status) ;
        copyUsrkeywrdsTovect(fptr,L1keywords);
   
        fits_close_file(fptr,&status);
        ref_frame_module_filename_track.clear();
         nFrameToAverage=nFrame_Avg_ref;
         flag_rfc = FALSE ;
        LOG (ERROR) << endl << "------------------Data Set " << dataindex + 1 << " : " << dirobj.sciencedatafile[dataindex] << "-----------------------" << endl ;
        /*---finding mode from science data file---*/
        getKeywordVal ((char *) dirobj.sciencedatafile[dataindex].c_str () , "OBS_MODE" , 1 , obsmode) ;
        if (strcasecmp (obsmode , "PC") != 0)
        {
            LOG (ERROR) << endl << "Observation mode is " << obsmode << "  in file  " << dirobj.sciencedatafile[dataindex] ;
            LOG (ERROR) << endl << "Checking next file...." << endl ;
            continue ; //go to next file if obs mode is not IM
        }

        cout << endl << "Data Mode is " << obsmode << endl ;

        strcpy (outputdir , dirobj.level2path[dataindex].c_str ()) ;
        strcpy (lbtfile , (char *) dirobj.lbtfile.c_str ()) ;
        //----------DATAINGEST----------//
        LOG (ERROR) << endl << "===================================DATAINGEST==================================================="  ;
        DataIngest di_obj ;
        
        LOG(INFO)<<dirobj.sciencedatafile[dataindex].c_str ();
        LOG(INFO)<<caldbindir.c_str ();
        LOG(INFO)<<(char *) dirobj.tctfile.c_str ();
        LOG(INFO)<<(char *) dirobj.mkffile.c_str ();
        LOG(INFO)<<(char *) dirobj.gtifile[dataindex].c_str ();
        di_obj.read ((char *) dirobj.sciencedatafile[dataindex].c_str () , (char*) caldbindir.c_str () , (char *) dirobj.tctfile.c_str () ,(char *) dirobj.mkffile.c_str ()  ,(char *) dirobj.gtifile[dataindex].c_str () ,
                (char *) dirobj.lbtfile.c_str () , (char*)dirobj.attfile.c_str(), (char*) dirobj.darkDirectory.c_str () , att_flag_val,gti_flag , valid_bit , all_or_custom , outputdir , dropframe , parity_flag ,UTC_flag,crc_flag, clobber , history) ;
        di_obj.display () ;
        status = di_obj.DataIngestProcess () ;
        if (status)
        {
            LOG (ERROR) << endl << "***Error in Data Ingest Process***"  ;
LOG(ERROR)<<" CRASH DATAINGEST FAILED (uvtRelativeAspectPC.cpp)";
            continue ;
        }

        strcpy (moduleIndir , di_obj.getModuleOutdir ()) ;
        strcpy (Indir_dataIngest , di_obj.getModuleOutdir ()) ;
        LOG (INFO) << "Using directory " << moduleIndir << "  as input to uvtUnitConversion" ;

        sprintf (moduleoutdir_bp , "%s/%s_%s" , outputdir , moduleoutdir_badpix , VERSION) ;
        sprintf (moduleoutdir_uc , "%s/%s_%s" , outputdir , moduleoutdir_unit , VERSION) ;
        sprintf (moduleoutdir_ff , "%s/%s_%s" , outputdir , moduleoutdir_fltfield , VERSION ) ;
        sprintf (moduleoutdir_qemcp , "%s/%s_%s" , outputdir , moduleoutdir_qe , VERSION ) ;
        sprintf (moduleoutdir_pp , "%s/%s_%s" , outputdir , moduleoutdir_pixpad , VERSION) ;
        sprintf (moduleoutdir_cr , "%s/%s_%s" , outputdir , moduleoutdir_cosmicray , VERSION) ;
        sprintf (moduleoutdir_fi , "%s/%s_%s" , outputdir , moduleoutdir_FrameIntegration , VERSION) ;
        sprintf (moduleoutdir_sc , "%s/%s_%s" , outputdir , moduleoutdir_findstarcentroid , VERSION) ;
        sprintf (moduleoutdir_dd , "%s/%s_%s" , outputdir , moduleoutdir_detectordistortion , VERSION) ;
        sprintf (moduleoutdir_od , "%s/%s_%s" , outputdir , moduleoutdir_opticaldistortion , VERSION) ;
        sprintf (moduleoutdir_de , "%s/%s_%s" , outputdir , moduleoutdir_driftExercise , VERSION) ;
        sprintf (moduleoutdir_rfc , "%s/%s_%s" , outputdir , moduleoutdir_refFrameCal , VERSION) ;
        sprintf (moduleoutdir_centroidCorr , "%s/%s_%s" , outputdir , moduleoutdir_centCorr , VERSION) ;
        sprintf (moduleoutdir_centroidBias , "%s/%s_%s" , outputdir , moduleoutdir_centBias , VERSION) ;


        char caldb_common_dir[PIL_LINESIZE] , caldb_common_dir_qefile[PIL_LINESIZE],caldb_common_dir_od[PIL_LINESIZE] ;
        char caldb_temp_dir[PIL_LINESIZE] , caldb_common_dir_ff[FLEN_FILENAME] , caldb_common_dir_cc[FLEN_FILENAME] ;
        strcpy (caldb_temp_dir , (char*) caldbindir.c_str ()) ;
        strcpy (caldb_common_dir , (char*) caldbindir.c_str ()) ;
        strcpy (caldb_common_dir_qefile , (char*) caldbindir.c_str ()) ;
        strcpy (caldb_common_dir_ff , (char*) caldbindir.c_str ()) ;
        strcpy (caldb_common_dir_cc , (char*) caldbindir.c_str ()) ;
        strcpy (caldb_common_dir_od , (char*) caldbindir.c_str ()) ;
        double factor ;

        //
        char infofile_in[PIL_LINESIZE] ;
        string tempfilepath = searchFile (moduleIndir , ".info") ;
        if (tempfilepath == " ")
        {
            LOG (ERROR) << "***Information file not found in " << moduleIndir << "***" ;
            continue ;
        }
        sprintf (infofile_in , "%s/%s" , moduleIndir , tempfilepath.c_str()) ;

        /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
        if (! (FileExists (infofile_in)))
        {
            LOG (ERROR) << "***Input FileList not Found at Specified PATH,Check Input Directory***" ;
            continue ;
        }
 vector<string> header_info;
        getHistory (header_info);
        writeHistory (infofile_in,header_info);
        fitsfile *finfo_in , *finfo_out ;
        fits_open_file (&finfo_in , infofile_in , READWRITE , &status) ;
        printError (status , "Error in opening the information file") ;
        fits_read_key (finfo_in , TINT , "WIN_X_SZ" , &win_xsize , NULL , &status) ;
        printError (status , "Error in reading the key value of the Window xsize" , infofile_in) ; //for creating name for output information file
        fits_read_key (finfo_in , TINT , "WIN_Y_SZ" , &win_ysize , NULL , &status) ;
        printError (status , "Error in reading the key value of the Window ysize " , infofile_in) ; //for creating name for output information
        fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
        printError (status , "Error in moving to 2nd HDU") ;
        datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
        xsize = datainfo.getXsize () ;
        ysize = datainfo.getYsize () ;
        int nframes = 0 ;
        ;
        char nameprefix[PIL_LINESIZE] ;
        fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
        printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file
        fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
        printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file
        fits_read_key (finfo_in , TSTRING , "IMGFILE" , imgfile , NULL , &status) ;
        printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file
        //    readKeywords (infofile_in , 2 , 3 , TSTRING , "EVTFILE" , eventfile ,
        //                TSTRING , "IMGFILE" , imgfile) ;
        fits_read_key (finfo_in , TSTRING , "DARKDIR" , darkdir , NULL , &status) ;
        printError (status , "***Error in reading the  key value of the NFILES ***" , infofile_in) ;
        long naxes[2] ;
    naxes[0] = naxes[1] = 0 ;
    int naxis = 2 ;
    int bitpix = FLOAT_IMG ;
        fits_create_img (finfo_in , bitpix , naxis , naxes , &status) ;
    printError (status , "Error in Creating the image for Signal Fie" , infofile_in);
        
           fits_close_file (finfo_in , &status) ;
        printError (status , "Error in closing the information file" , infofile_in) ;
       
        
        
        fits_open_file (&finfo_in , infofile_in , READWRITE , &status) ;
        printError (status , "Error in opening the input information file",infofile_in) ;
        fits_movabs_hdu (finfo_in , 3 , NULL , &status) ;
        printError (status , "Error in moving to 2nd HDU",infofile_in) ;
         for(int i=0;i<L1keywords.size ();i++)
         {
                fits_write_record (finfo_in , L1keywords[i].c_str () , &status) ;
         }
        fits_close_file (finfo_in , &status) ;
        printError (status , "Error in closing the information file" , infofile_in) ;
        
        int tfields = 11 ;
        char *ttype[] = {"PacketSequence " , "FrameCount" , "Time" , "Ix" , "Fx" , "Iy" , "Fy" , "Max-Min" , "Min" , "BAD FLAG" , "MULT PHOTON"} ;
        char *tform[] = {"U" , "U" , "D" , "U" , "E" , "U" , "E" , "B" , "B" , "B" , "B"} ;
        char *tunit[] = {"" , "" , "" , "" , "" , "" , "" , "" , "" , "" , ""} ;
        char infile[NAMESIZE] ;
        char errstr[500] ;
        unsigned short frameno ;
        double frametime , integrationtime ;
        long fpixel[2] ;
        fpixel[0] = fpixel[1] = 1 ;
        //long naxes[2] ;
        naxes[0] = naxes[1] = xsize ;

        status = takeDarkinfo () ;
        if (status)
        {
            LOG (ERROR) << "***Error getting the information of Darks***" ;
            continue ;
        }
        //cout<<dstartpath<<" "<<dendpath<<endl;exit(1);
        status = readDarkFrame (dstartpath , darkFramestart_data) ;

        status = readDarkFrame (dendpath , darkFrameend_data) ;

        float *finalFlatFielddata ;
        
        if(UTC_flag==0 && qemcpFlag==1){
                LOG(WARNING)<<"***QEMCP correction cant be kept ON if UTC correction not to be done!!,NOW skipping QEMCP correction***";
                qemcpFlag=0;
            }
        
        if (qemcpFlag)
        {
            string    tempname = caldb_handler.getQEFile (datainfo.getDetector () , datainfo.getObsMode () , caldb_common_dir_qefile) ;
            if (tempname == " ")
            {
                LOG (ERROR) << endl << "Couldn't find QEMCP file from caldb" << endl ;
                continue ;
            }

            joinStrings (qeFile , 2 , caldb_common_dir_qefile , tempname.c_str()) ;

             status= readNumRowsFromFits (qeFile,2,nCalDBTempValues);
           if (status)
            {
                LOG (ERROR) << "Error in reading the number of rows from fits file "<<qeFile ;
                return (EXIT_FAILURE) ;
            }
  
   temp= new float[nCalDBTempValues];
   f0 = new float[nCalDBTempValues],f1=new float[nCalDBTempValues];f2=new float[nCalDBTempValues];
    f3 = new float[nCalDBTempValues],f4=new float[nCalDBTempValues];f5=new float[nCalDBTempValues];
     f6 = new float[nCalDBTempValues],f7=new float[nCalDBTempValues];
   status=readQEMCPFile (qeFile,datainfo.getDetector (),nCalDBTempValues,temp,f0,f1,f2,f3,f4,f5,f6,f7);
   if (status)
            {
                LOG (ERROR) << "Error in reading the QEMCP caldb file ->"<<qeFile ;
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
    status = getTemp ((char *)dirobj.lbtfile.c_str (),datainfo.getDetector (),time_lbt,insideTemp,outsideTemp,nrows_lbt) ;
    if (status)
    {
        LOG(ERROR) << "***temperature reading from the lbt file unsuccessful***" ;
        continue;
    }
            
//            status = getTemp () ;
//            if (status)
//            {
//                LOG (ERROR) << "***temperature reading from the lbt file unsuccessful***"  ;
//                continue ;
//            }
            int filter_coln ;
            int filternumber ;
            sprintf (filter , "%s" , datainfo.getFilter ()) ;

            qe_mg_factor = new float[nCalDBTempValues] ;
            if (filter == (string) "F0")
            {
                filternumber = 0 ;
                filter_coln = 1 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = f0[q] ;
            }
            else  if (filter == (string) "F1")
            {
                filternumber = 1 ;
                filter_coln = 2 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = f1[q] ;
            }
            else if (filter == (string) "F2")
            {
                filternumber = 2 ;
                filter_coln = 3 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = f2[q] ;
            }
            else if (filter == (string) "F3")
            {
                filternumber = 3 ;
                filter_coln = 4 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = f3[q] ;
            }
            else if (filter == (string) "F4")
            {
                filternumber = 4 ;
                filter_coln = 5 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = f4[q] ;
            }
            else if (filter == (string) "F5")
            {
                filternumber = 5 ;
                filter_coln = 6 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = f5[q] ;
            }
            else if (filter == (string) "F6")
            {
                filternumber = 6 ;
                filter_coln = 7 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = f6[q] ;
            }
            else if (filter == (string) "F7")
            {
                filternumber = 7 ;
                filter_coln = 8 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = f7[q] ;
            }

            else
            {
                LOG (ERROR) << endl << "***Invalid filter option*** " << endl ;
                continue ;
            }

        }
        if (flatfieldFlag)
        {

            string tempname1 = caldb_handler.getFlatFieldFile (datainfo.getDetector () , datainfo.getObsMode () , datainfo.getFilter () , caldb_common_dir_ff) ;
            // tempname = caldb_handler.getFlatFieldFile (datainfo.getDetector () , datainfo.getObsMode () , datainfo.getFilter () , (char*)caldbindir.c_str ()) ;

            joinStrings (flatfieldfile , 2 , caldb_common_dir_ff , tempname1.c_str()) ;
            flatfielddata = new float[xsize * ysize] ;
            status = readImage (flatfieldfile , 1 , flatfielddata) ;
            finalFlatFielddata = new float[600 * 600] ;
            for (int i = 0 ; i < 600 * 600 ; i ++)
            {
                finalFlatFielddata[i] = 0.0f ;
            }

            /*padding is done to  match xsize of event file to caldb file*/
            status = Applypadding (flatfielddata , xsize , ysize , finalFlatFielddata , 600 , 600) ;
            if (status)
            {
                LOG (ERROR) << "***Padding Fails For the flatfield CALDB data***" ;
                continue ;
            }
        }


        //reading caldb file from CALDB directory

        string tempname = caldb_handler.getBadPixelFile (datainfo.getDetector () , datainfo.getObsMode () , win_xsize + 1 , win_ysize + 1 , (char*) caldbindir.c_str ()) ;
        if (tempname == " ")
        {
            LOG (ERROR) << endl << "Couldn't find bad pixel file from calDB" ;
            continue ;
        }
        char outfile[FLEN_FILENAME] ;
        char  file_in[FLEN_FILENAME] ;
        joinStrings (badpixfile , 1 , tempname.c_str()) ;
        badpixdata = new float[xsize * ysize] ;
        status = readImage (badpixfile , 1 , badpixdata) ;

        //reading centroid Corr caldb file
        tempname = caldb_handler.getCentroidEAFile (datainfo.getDetector () , caldb_common_dir_cc) ;
        if (tempname == " ")
        {
            LOG (ERROR) << "Couldn't find CentroidBias file from caldb" << endl ;
            continue ;
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
        // EA= new long[n_ele];
        fits_read_col (f_ea , TINT , 1 , 1 , 1 , n_ele , NULL , (void *) &EA , NULL , &status) ;
        printError (status , "Error in closing the info out file" , centroidEAfile) ;
        fits_close_file (f_ea , &status) ;
        sprintf (file_in , "%s/%s" , moduleIndir , eventfile) ; //taking event file full path



        LOG (INFO) << " Input Event File " << file_in << endl ;
        long  nrows = 0 ;
        char centroid_algo[FLEN_FILENAME] ;
        fitsfile *fevt_in , *fout ;
        //reading the event file information
        fits_open_file (&fevt_in , file_in , READWRITE , &status) ;
        printError (status , "Error in opening the input event file" , file_in) ;
        fits_read_key (fevt_in , TSTRING , "CENTROID" , centroid_algo , NULL , &status) ;
        printError (status , "Error Reading the Key value for CENTROID" , infofile_in) ;
        if ((strcmp (centroid_algo , "3C") != 0) && (strcmp (centroid_algo , "3S") != 0) && (strcmp (centroid_algo , "5S") != 0))
        {
            LOG (ERROR) << "***Invalid value of Centroid_algo ***" ;
            continue ;
        }


        if ((strcmp (centroid_algo , "3C") == 0))
        {
            cent_corr_win_size = 5 ;          //3x3 cross  (5 pixels - centre, top, bottom, left & right used)
            cent_bias_win_size = 3 ;
        }
        else if (strcmp (centroid_algo , "5S") == 0)
        {
            cent_corr_win_size = 25 ;         //5x5 square
            cent_bias_win_size = 5 ;
        }
        else if ( (strcmp (centroid_algo , "3S") == 0))
        {
            cent_corr_win_size = 9 ;             //3x3 square
            cent_bias_win_size = 3 ;
        }
        else
        {
            LOG (ERROR) << "***INVALID value for Centroid window.allowed values are (3/5) *** "  ;
            continue ;
        }


        // cout<<(char*)caldbindir.c_str () <<endl;exit(1);
        tempname = caldb_handler.getCentroidBiasFile (datainfo.getDetector () , caldb_common_dir , cent_bias_win_size) ;

        if (tempname == " ")
        {
            LOG (ERROR) << "Couldn't find CentroidBias file from caldb"  ;
            continue ;
        }

        /***For Joining the String ***/

        // cout<<caldb_common_dir<<endl;exit(1);
        joinStrings (centroidbiasfile , 2 , caldb_common_dir , tempname.c_str()) ;
        LOG (INFO) << "Centroid bias  file  from calDB :" << centroidbiasfile ;
        //
        //    /*Reading a bias file From the CALDB dir*/
        status = readcentroidbiasFile () ;
        if (status)
        {
            LOG (ERROR) << "***Error Reading centroid bias File From the Caldb***" << endl ;
            continue ;
        }
        //detector Distortion File Reading
        tempname = caldb_handler.getDetectorFile (datainfo.getDetector () , caldb_temp_dir) ;
        joinStrings (detector_distortion_corr_file , 2 , caldb_temp_dir , tempname.c_str()) ;
        LOG (INFO) << "Distortion correction file is " << detector_distortion_corr_file ;
        if (tempname == " ")
        {
            LOG (ERROR) << endl << "Couldn't find detector File From caldb" ;
            continue ;
        }

        X_detect_distortion = new float[xsize * ysize] ;
        Y_detector_distortion = new float[xsize * ysize] ;
        status = caldb_handler.readCaldbDistFile (X_detect_distortion , Y_detector_distortion , detector_distortion_corr_file) ;

        tempname = caldb_handler.getOpticalDetectorFile (datainfo.getDetector () , caldb_common_dir_od , "F0") ;

        if (tempname == " ")
        {
            LOG (ERROR) << endl << "Couldn't find Optical Dist detector File From caldb" ;
            continue ;
        }
        joinStrings (optical_distortion_corr_file , 2 , caldb_common_dir_od , tempname.c_str()) ;
        X_optical_distortion = new float[xsize * ysize] ;
        Y_optical_distortion = new float[xsize * ysize] ;
        status = caldb_handler.readCaldbDistFile (X_optical_distortion , Y_optical_distortion , optical_distortion_corr_file) ;

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
        unsigned  char *dmm = new unsigned char[nrows] ;
        unsigned char *Min = new  unsigned char[nrows] ;
        unsigned char *multflag = new unsigned char [nrows] ;
        unsigned char *badflag = new unsigned char [nrows] ;

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
        fits_close_file (fevt_in , &status) ;

        int TotalNum_Effected_Rows = performMaskBadPix (Ix , Iy , badpixdata , dmm , badflag , multflag , nrows , xsize , ysize , thr_multiph) ;
        if (status)
        {
            LOG (ERROR) << "Error in Mask bad pixel module " ;
            continue ;
        }
        if (wtd_bp == 1)
        {

            status = setDirectoryStructure (moduleoutdir_bp , " ") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" ;
                continue ;
            }
            sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_bp , nameprefix , "bp") ;
            fits_create_file (&fout , outfile , &status) ;
            printError (status , "Error creating the output File " , outfile) ;
            fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Events" , &status) ;
            printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
            fits_close_file (fout , &status) ;
            printError (status , "Error closing the file " , outfile) ;

            status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TUSHORT , 4 , Ix , nrows , TFLOAT , 5 , fx , nrows , TUSHORT , 6 , Iy , nrows ,
                    TFLOAT , 7 , fy , nrows , TBYTE , 8 , dmm , nrows , TBYTE , 9 , Min , nrows , TBYTE , 10 , badflag , nrows , TBYTE , 11 , multflag , nrows) ;
            if (status)
            {
                LOG (INFO) << "Error in writing to the  bad pixel FITS  file" ;
                continue ;
            }

        }
            status = performPixPadding (Ix , Iy , nrows) ;
        xsize = 600 ;
        ysize = 600 ;
        if (wtd_pp == 1)
        {

            status = setDirectoryStructure (moduleoutdir_pp , " ") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" ;
                continue ;
            }

            sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_pp , nameprefix , "pp") ;
            fits_create_file (&fout , outfile , &status) ;
            printError (status , "Error creating the output File " , outfile) ;
            fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Events" , &status) ;
            printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
            fits_close_file (fout , &status) ;
            printError (status , "Error closing the file " , outfile) ;

            status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TUSHORT , 4 , Ix , nrows , TFLOAT , 5 , fx , nrows , TUSHORT , 6 , Iy , nrows ,
                    TFLOAT , 7 , fy , nrows , TBYTE , 8 , dmm , nrows , TBYTE , 9 , Min , nrows , TBYTE , 10 , badflag , nrows , TBYTE , 11 , multflag , nrows) ;
            if (status)
            {
                LOG (INFO) << "Error in writing to the  Pix padding FITS  file" ;
                continue ;
            }

        }

        //performing cosmic Ray correction
        long nFrames ;
        nFrames = frame_no[nrows - 1] ; //number of frames

        long nFramesCounter = 0 , j = 1 ;
        vector<long> del_CReffctedrows , frameno_crFailed ;
        vector<double> x_crFailed , y_crFailed ;
        long nEventsInFrame[nFrames] , GoodFrames[nFrames] , CRAFFECTED[nFrames] , m = 0 ;

         for (int i=0;i<nFrames;i++)  nEventsInFrame[i]=0;
        
        //if only one event  is there than
        if(nrows==1){
            nEventsInFrame[nFrames-1]=1;
            GoodFrames[nFrames-1]=frame_no[0];
            j=1;
        }
        
        //
        
        for (long i = 0 ; i < nrows - 1 ; i ++)
        {
            if ((frame_no[i] == frame_no[i + 1]))
            {

                nFramesCounter ++ ; //counter for number of events in a frame
                if (i == nrows - 2)
                {
                    nEventsInFrame[j-1] = nFramesCounter ;
                    GoodFrames[j-1] = frame_no[i] ;
                }

            }
            else
            {
                nEventsInFrame[j-1] = nFramesCounter ;
                GoodFrames[j-1] = frame_no[i] ;
                j ++ ;
                nFramesCounter = 1 ;
                if ((frame_no[i + 1] - frame_no[i]) > 1)
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
        unsigned short Image_Array[xsize][ysize] ;

        for (int ii = 0 ; ii < xsize ; ii ++)
        {
            for (int jj = 0 ; jj < ysize ; jj ++)
            {
                Image_Array[ii][jj] = 0 ;
            }
        }
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
        vector<double> tm_frm ;
        vector<float> fract_X , fract_Y ;
        vector<unsigned char> dmm_temp , min_temp , badflag_temp , multflag_temp ;
        bool flag_cr = 0 ;

        flag_cr = 0 ;
        
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
        }
         }
        
         else{
        
        
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
        }
        
         }
       

        delete[] frame_no , psc , Ix , Iy , time_frame , fx , fy , dmm , Min , badflag , multflag ;
        psc = new unsigned short[frmno.size ()] ;
        frame_no = new unsigned short[frmno.size ()] ;
        Ix = new unsigned short[frmno.size ()] ;
        Iy = new unsigned short[frmno.size ()] ;
        time_frame = new double[frmno.size ()] ;
        fx = new float[frmno.size ()] ;
        fy = new float[frmno.size ()] ;
        dmm = new unsigned char[frmno.size ()] ;
        Min = new  unsigned char[frmno.size ()] ;
        multflag = new unsigned char [frmno.size ()] ;
        badflag = new unsigned char [frmno.size ()] ;

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

        }
        if (wtd_cr == 1)
        {

            status = setDirectoryStructure (moduleoutdir_cr , " ") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_cr , nameprefix , "cr") ;
            fits_create_file (&fout , outfile , &status) ;
            printError (status , "Error creating the output File " , outfile) ;
            fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Events" , &status) ;
            printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
            fits_close_file (fout , &status) ;
            printError (status , "Error closing the file " , outfile) ;

            status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , frmno.size () , TUSHORT , 2 , frame_no , frmno.size () , TDOUBLE , 3 , time_frame , frmno.size () , TUSHORT , 4 , Ix , frmno.size () , TFLOAT , 5 , fx , frmno.size () , TUSHORT , 6 , Iy , frmno.size () ,
                    TFLOAT , 7 , fy , frmno.size () , TBYTE , 8 , dmm , frmno.size () , TBYTE , 9 , Min , frmno.size () , TBYTE , 10 , badflag , frmno.size () , TBYTE , 11 , multflag , frmno.size ()) ;
            if (status)
            {
                LOG (ERROR) << "Error in writing to the  Pix padding FITS  file" ;
                continue ;
            }

        }
         float *newXfract = new float[frmno.size ()] ;
        float *newYfract = new float [frmno.size ()] ;
      //  performing centroid Correction
        //tobe removed
      
//        for (int i=0;i<nFrames;i++){
//                 frmno.push_back (frame_no[i]) ;
//              
//        }
        //till this
       

//        status = performCentroidCorr (time_frame , Ix , Iy , fx , fy , Min , newXfract , newYfract , (float*) darkFramestart_data , (float*) darkFrameend_data , datainfo.getIntegrationTime () , EA , (int) frmno.size ()) ;
//      
//        if (status)
//        { 
//             
//            LOG (INFO) << "Error in Centroid Correction module" << endl ;
//            continue ;
//        }
//      char *ttype_1[] = {"PacketSequence " , "FrameCount" , "Time" ,  "Fx"  , "Fy" , "Max-Min" , "Min" , "BAD FLAG" , "MULT PHOTON"} ;
//      char * tform_1[] = {"U" , "U" , "D" , "1D" , "1D" , "B" , "B" , "B" , "B"} ;
//        char * tunit_1[] = {"" , "" , "" , "" , "" , "" , "" , "" , ""} ;
//        if (wtd_centCorr == 1)
//        {
//            tfields = 11;
//            status = setDirectoryStructure (moduleoutdir_centroidCorr , " ") ;
//            if (status)
//            {
//                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
//                continue ;
//            }
//            sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_centroidCorr , nameprefix , "cc") ;
//            fits_create_file (&fout , outfile , &status) ;
//            printError (status , "Error creating the output File " , outfile) ;
//            fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Events" , &status) ;
//            printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
//            fits_close_file (fout , &status) ;
//            printError (status , "Error closing the file " , outfile) ;
//
//             status = writeColumnsToFITS (outfile , 2 , 11 , TUSHORT , 1 , psc , frmno.size () , TUSHORT , 2 , frame_no , frmno.size () , TDOUBLE , 3 , time_frame , frmno.size () , TUSHORT , 4 , Ix , frmno.size () , TFLOAT , 5 , fx , frmno.size () , TUSHORT , 6 , Iy , frmno.size () ,
//                    TFLOAT , 7 , fy , frmno.size () , TBYTE , 8 , dmm , frmno.size () , TBYTE , 9 , Min , frmno.size () , TBYTE , 10 , badflag , frmno.size () , TBYTE , 11 , multflag , frmno.size ()) ;
//            if (status)
//            {
//                LOG (ERROR) << "Error in writing to the  Pix padding FITS  file" ;
//                continue ;
//            }
//
//        }
        tfields = 10 ;
       char  *ttype_2[] = {"PacketSequence " , "FrameCount" , "Time" ,  "Fx"  , "Fy" , "Max-Min" , "Min" , "BAD FLAG" , "MULT PHOTON" , "EFFECTIVE_NUM_PHOTONS"} ;
        char *tform_2[] = {"U" , "U" , "D" , "1D" , "1D" , "B" , "B" , "B" , "B" , "1D"} ;
       char *tunit_2[] = {"" , "" , "" , "" , "" , "" , "" , "" , "" , ""} ;
       // nrows = frmno.size () ;
        float *new_tempX = new float [frmno.size ()] ;
        float *new_tempY = new float [frmno.size ()] ;
  for (int i=0;i<frmno.size ();i++)
   {
        new_tempX[i]=Ix[i]+fx[i];
        new_tempY[i]=Iy[i]+fy[i];
       //new_tempX[i]=Ix[i]+newXfract[i];
       // new_tempY[i]=Iy[i]+newYfract[i];
        
    }
        //status = performCentroidBias (frmno.size () , newXfract , newYfract , new_tempX , new_tempY) ;

        float *effective_NumPhotons = new float[frmno.size ()] ;
        for (int i = 0 ; i < frmno.size () ; i ++) effective_NumPhotons[i] = 1.0f ;
        if (wtd_centBias == 1)
        {

            status = setDirectoryStructure (moduleoutdir_centroidBias , " ") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***"  ;
                continue ;
            }
            sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_centroidBias , nameprefix , "cb") ;
            fits_create_file (&fout , outfile , &status) ;
            printError (status , "Error creating the output File " , outfile) ;
            fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
            printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
            fits_close_file (fout , &status) ;
            printError (status , "Error closing the file " , outfile) ;

            status = writeColumnsToFITS (outfile , 2 , 10 , TUSHORT , 1 , psc , frmno.size () , TUSHORT , 2 , frame_no , frmno.size () , TDOUBLE , 3 , time_frame , frmno.size () , TFLOAT , 4 , new_tempX , frmno.size () , TFLOAT , 5 , new_tempY ,
                    frmno.size () , TBYTE , 6 , dmm , frmno.size () , TBYTE , 7 , Min , frmno.size () , TBYTE , 8 , badflag , frmno.size () , TBYTE , 9 , multflag , frmno.size () , TFLOAT , 10 , effective_NumPhotons , frmno.size ()) ;
            if (status)
            {
                LOG (ERROR) << "Error in writing to the  Pix padding FITS  file" ;
                continue ;
            }

        }
      
        if (flatfieldFlag)
        {
            status = performFlatFieldCorr (effective_NumPhotons , finalFlatFielddata , new_tempX , new_tempY , frmno.size ()) ;
            if (wtd_ff == 1)
            {

                status = setDirectoryStructure (moduleoutdir_ff , " ") ;
                if (status)
                {
                    LOG (ERROR) << "***Directory Structure has not been successfully set-up***"  ;
                    break ;
                }
                sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_ff , nameprefix , "ff") ;
                fits_create_file (&fout , outfile , &status) ;
                printError (status , "Error creating the output File " , outfile) ;
                fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
                printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
                fits_close_file (fout , &status) ;
                printError (status , "Error closing the file " , outfile) ;

                status = writeColumnsToFITS (outfile , 2 , 10 , TUSHORT , 1 , psc , frmno.size () , TUSHORT , 2 , frame_no , frmno.size () , TDOUBLE , 3 , time_frame , frmno.size () , TFLOAT , 4 , new_tempX , frmno.size () , TFLOAT , 5 , new_tempY ,
                        frmno.size () , TBYTE , 6 , dmm , frmno.size () , TBYTE , 7 , Min , frmno.size () , TBYTE , 8 , badflag , frmno.size () , TBYTE , 9 , multflag , frmno.size () , TFLOAT , 10 , effective_NumPhotons , frmno.size ()) ;
                if (status)
                {
                    LOG (INFO) << "Error in writing to the  Pix padding FITS  file" ;
                    continue ;
                }

            }
        }
        nrows=frmno.size ();
        //float intgrtion_time=datainfo.getIntegrationTime ();cout<<intgrtion_time<<endl;exit (1);
        status = performUnitConversion (effective_NumPhotons , datainfo.getIntegrationTime () , nrows) ;
        if(status){
            LOG(ERROR)<<"***Integration time is coming ZERO***";
            LOG(ERROR)<<" CRASH FRAME INTEGRATION TIME =< ZERO (uvtRelativeAspectPC.cpp)";
            continue;
        }
        if (wtd_uc == 1)
        {

            status = setDirectoryStructure (moduleoutdir_uc , " ") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" ;
                continue ;
            }
            sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_uc , nameprefix , "uc") ;
            fits_create_file (&fout , outfile , &status) ;
            printError (status , "Error creating the output File " , outfile) ;
            fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
            printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
            fits_close_file (fout , &status) ;
            printError (status , "Error closing the file " , outfile) ;

            status = writeColumnsToFITS (outfile , 2 , 10 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                    nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows , TBYTE , 8 , badflag , nrows , TBYTE , 9 , multflag , nrows , TFLOAT , 10 , effective_NumPhotons , nrows) ;
            if (status)
            {
                LOG (INFO) << "Error in writing to the  Pix padding FITS  file" ;
                continue ;
            }

        }
        bool flag_qemcp_success=FALSE;
        if (qemcpFlag)
        {
            double t1 , t2 , x1 , x2 , factor ;
            //ofstream of1("temp2.txt");
            for (int p = 0 ; p < nrows ; p ++)
            {
                temperature = - 9999 ;
                //            double delta_time=time_lbt[1]-time_lbt[0];
                //            int index=(time_frame[p]-time_lbt[0])/delta_time;
                //            temperature=(insideTemp[index+1]+outsideTemp[index+1])/2;
                for (int i = 0 ; i < nrows_lbt ; i ++)
                {
                    //cout<<time_lbt[i]<<" "<<time_lbt[i+1]<<endl;

                    if (time_frame[p] >= time_lbt[i] && time_frame[p] < time_lbt[i + 1])
                    {

                        temperature = (insideTemp[i] + outsideTemp[i]) / 2 ;
                        break ;
                    }
                }
                // of1<<temperature<<endl;;

                if (temperature == - 9999)
                {
                    LOG (ERROR) << "No record found in LBT file" << p << " " << time_frame[p] << endl ;
                    break;
                }
                flag_qemcp_success=TRUE;
                for (int j = 0 ; j < nCalDBTempValues - 1 ; j ++)
                {
                    if (temperature >= temp[j] && temperature < temp[j + 1])          //temperature - from LBT file ,temp- from caldb file
                    {
                        t1 = temp[j] ;
                        t2 = temp[j + 1] ;
                        if ((t2 - t1) == 0)
                        {
                            LOG (INFO) << "***Divide By zero***" << endl ;
                            continue ;
                        }
                        x1 = qe_mg_factor[j] ;
                        x2 = qe_mg_factor[j + 1] ;
                        factor = x1 + ((temperature - t1)*((x2 - x1) / (t2 - t1))) ;
                    }
                }

                effective_NumPhotons[p] = effective_NumPhotons[p] * factor ;


            }
            if(flag_qemcp_success==FALSE)
            {   LOG(ERROR)<<"Error in finding the related  factor for QEMCP correction";
                continue;//next science data file
            }
            
           
            

            // status=performQEMCPcorrection (effective_NumPhotons,nrows,factor);
            if (wtd_qemcp)
            {
                status = setDirectoryStructure (moduleoutdir_qemcp , " ") ;
                if (status)
                {
                    LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                    continue ;
                }
                sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_qemcp , nameprefix , "qe") ;
                fits_create_file (&fout , outfile , &status) ;
                printError (status , "Error creating the output File " , outfile) ;
                fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
                printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
                fits_close_file (fout , &status) ;
                printError (status , "Error closing the file " , outfile) ;

                status = writeColumnsToFITS (outfile , 2 , 10 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                        nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows , TBYTE , 8 , badflag , nrows , TBYTE , 9 , multflag , nrows , TFLOAT , 10 , effective_NumPhotons , nrows) ;
                if (status)
                {
                    LOG (INFO) << "Error in writing to the  Pix padding FITS  file" ;
                    continue ;
                }


            }
        }

        //detector Dist Correction 

        status = performDistCorrection (nrows , new_tempX , new_tempY , X_detect_distortion , Y_detector_distortion , 512) ;
        if (wtd_dd == 1)
        {

            status = setDirectoryStructure (moduleoutdir_dd , " ") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" ;
                continue ;
            }
            sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_dd , nameprefix , "dd") ;
            fits_create_file (&fout , outfile , &status) ;
            printError (status , "Error creating the output File " , outfile) ;
            fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
            printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
            fits_close_file (fout , &status) ;
            printError (status , "Error closing the file " , outfile) ;

            status = writeColumnsToFITS (outfile , 2 , 10 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                    nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows , TBYTE , 8 , badflag , nrows , TBYTE , 9 , multflag , nrows , TFLOAT , 10 , effective_NumPhotons , nrows) ;
            if (status)
            {
                LOG (ERROR) << "Error in writing to the  Detector Distortion  Correction File" ;
                continue ;
            }

        }
        status = performDistCorrection (nrows , new_tempX , new_tempY , X_optical_distortion , Y_optical_distortion , 512) ;
        if (wtd_od == 1)
        {
            status = setDirectoryStructure (moduleoutdir_od , " ") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            sprintf (outfile , "%s/%s_%s.fits" , moduleoutdir_od , nameprefix , "od") ;
            fits_create_file (&fout , outfile , &status) ;
            printError (status , "Error creating the output File " , outfile) ;
            fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype_2 , tform_2 , NULL , "Events" , &status) ;
            printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
            fits_close_file (fout , &status) ;
            printError (status , "Error closing the file " , outfile) ;

            status = writeColumnsToFITS (outfile , 2 , 10 , TUSHORT , 1 , psc , nrows , TUSHORT , 2 , frame_no , nrows , TDOUBLE , 3 , time_frame , nrows , TFLOAT , 4 , new_tempX , nrows , TFLOAT , 5 , new_tempY ,
                    nrows , TBYTE , 6 , dmm , nrows , TBYTE , 7 , Min , nrows , TBYTE , 8 , badflag , nrows , TBYTE , 9 , multflag , nrows , TFLOAT , 10 , effective_NumPhotons , nrows) ;
            if (status)
            {
                LOG (ERROR) << "Error in writing to the  optical Assmbly Distortion  Correction File" ;
                continue ;
            }

        }
        bool flag_fsc = FALSE ;
        //float *one_dim_exp ;
        //one_dim_exp = new float[IMG_DIM_FI * IMG_DIM_FI] ;
       // float *one_dim_img ;
        // one_dim_img = new float[IMG_DIM_FI * IMG_DIM_FI] ;
        //cout<<"KKKK"<<endl;exit(1);
        vector<FrameIntegration_Arr_strct> frameIntegration_track ;
        unsigned short *mult_temp = new unsigned short[nrows] ;
        unsigned short *badFlag_temp = new unsigned short[nrows] ;

        for (int i = 0 ; i < nrows ; i ++)
        {
            mult_temp[i] = multflag[i] ;
            badFlag_temp[i] = badflag[i] ;
        }
        
        //step for calculating RAM needed for run pipelining 
        long tot_frm = frame_no[nrows - 1] - frame_no[0] ;
        //LOG(INFO)<<(tot_frm - nFrameDiscard_fi) / nFrameIntegrate_fi<<" "<<nFrameIntegrate_fi<<endl;
//        try
//        {
//            float *temp_Array_sig = new float[IMG_DIM_FI * IMG_DIM_FI * ((int) (tot_frm - nFrameDiscard_fi) / nFrameIntegrate_fi)] ;
//           // float *temp_Array_exp = new float[IMG_DIM_FI * IMG_DIM_FI * ((int) (tot_frm - nFrameDiscard_fi) / nFrameIntegrate_fi)] ;
////            delete[] temp_Array_sig , temp_Array_exp ;
//            delete[] temp_Array_sig  ;
//        }
//        catch (bad_alloc)
//        {
//            LOG (ERROR) << "***Memory(RAM)  is not sufficient for frame integration .please  run  chain on system with  higher configuration  OR  increase the number of frames to be integrated .***" ;
//            return (EXIT_FAILURE) ;
//        }
   
//        status = performFrameIntegration (nrows , frame_no , xsize , nFrameDiscard_fi , nFrameIntegrate_fi , new_tempX , new_tempY , mult_temp , badFlag_temp , time_frame , effective_NumPhotons , one_dim_img , one_dim_exp , frameIntegration_track) ;
//        if (status)
//        {
//            LOG (ERROR) << "Error in frameIntegration " ;
//            return (EXIT_FAILURE) ;
//        }
     
//        if (wtd_fi)
//        {
//            status = setDirectoryStructure (moduleoutdir_fi , "SignalFrames") ;
//            if (status)
//            {
//                LOG (ERROR) << "***Directory Structure has not been successfully set-up***"  ;
//                continue ;
//            }
//            status = setDirectoryStructure (moduleoutdir_fi , "ExposureFrames") ;
//            if (status)
//            {
//                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" ;
//                continue ;
//            }
//            for (int i = 0 ; i < frameIntegration_track.size () ; i ++)
//            {
//                status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "SignalFrames" , "sig" , frameIntegration_track[i].pixels_sig_fi.data () , nameprefix , frameIntegration_track[i].frame_time_fi , frameIntegration_track[i].frame_number_fi , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output
//                if (status)
//                {
//                    LOG (ERROR) << "***Writing to Disk Fails***" ;
//                    continue ;
//                }
////                status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "ExposureFrames" , "exp" , frameIntegration_track[i].pixels_exp_fi , nameprefix , frameIntegration_track[i].frame_time_fi , frameIntegration_track[i].frame_number_fi , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output
////                if (status)
////                {
////                    LOG (ERROR) << "***Writing to Disk Fails***" ;
////                    continue ;
////                }
//
//            }
//        }
        char *tform2[] = {"E" , "E" , "E"} ;
        tfields = 3 ;
      char *ttype_3[] = {"X" , "Y" , "Intensity"} ;
       // tform = {"I" , "I" , "E"} ;

        if (wtd_fsc == 1)
        {
            status = setDirectoryStructure (moduleoutdir_sc , "Star") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_sc , "Centroid") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***"  ;
                continue ;
            }



        }
        status = setDirectoryStructure (moduleoutdir_rfc , "Centroid") ;
        if (status)
        {
            LOG (ERROR) << "***Directory Structure has not been successfully set-up***"  ;
            continue ;
        }
        int cnt_forRefFrmCal = 0 ;
        int cnt_temp = 0 ;
        int trackcnt_forRefFrame = 0 ;
        double time_RefFrame_Avg = 0.0f ;

        int cnt_finish = 0 ;
        uvtDetectStar obj_sc ;
      long currframe_end_index=0;
        int nextframe_start_index=0;
         long tmp_startrow=0;
          float *one_dim_exp_final ;
    one_dim_exp_final = new float[IMG_DIM_FI * IMG_DIM_FI] ;
    float *one_dim_img_final ;
    one_dim_img_final = new float[IMG_DIM_FI * IMG_DIM_FI] ;
     long frmno_fi_new;
   double frm_time_fi_new;
       long num_frm_count=0;
       frame_Data= new float[IMG_DIM_FI*IMG_DIM_FI];
       frame_ExpData= new float[IMG_DIM_FI*IMG_DIM_FI];
       time_t t1,t2;
        bool flag_firstFrame_accomplished=FALSE;
        lastFrame_flag=FALSE;
//        LOG(INFO)<<" number of frames to integrate "<<nFrameIntegrate_fi;exit(1);
        if (nFrameIntegrate_fi>frame_no[nrows-1]){
            LOG(ERROR)<<"***Number of frames to accumulate is greater than total available frames***";
            continue;
            
        }
         for (int i = 0 ; i <nrows ; i=i+(currframe_end_index-nextframe_start_index))
      // for (int i = 0 ; i < frameIntegration_track.size () ; i ++)
      //  {
         {
            
             num_frm_count++;
            if (i  > 0)    flag_fsc = TRUE ;
             t1=time (NULL) ; //Computes execution start time (in seconds)
            if(i==0)
            {
                 status = performFrameIntegration ((long)nrows , tmp_startrow,frame_no , xsize , nFrameDiscard_fi , nFrameIntegrate_fi , new_tempX , new_tempY , mult_temp , badFlag_temp , time_frame ,effective_NumPhotons , one_dim_img_final , one_dim_exp_final , frameIntegration_track,frame_Data,frame_ExpData,frm_time_fi_new,frmno_fi_new,currframe_end_index) ;
                 if (wtd_fi)//incase of frame integration to be written on the disk
                    {
                        status = setDirectoryStructure (moduleoutdir_fi , "SignalFrames") ;
                        if (status)
                        {
                            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                            continue;
                        }
                        status = setDirectoryStructure (moduleoutdir_fi , "ExposureFrames") ;
                        if (status)
                        {
                            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                            continue;
                        }
                    // for (int i = 0 ; i < frameIntegration_track.size () ; i ++)//loop for total  number of frames generated after frame integration
                    // {
                            status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "SignalFrames" , "sig" , frame_Data , nameprefix ,frm_time_fi_new ,num_frm_count , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output
                            if (status)
                            {
                                LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                                continue;
                            }
                            status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frm_time_fi_new ,num_frm_count , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output
                            if (status)
                            {
                                LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                                continue;
                            }

                    //}
                    }
                 t2 =time (NULL) ; //Computes execution start time (in seconds)
                 //LOG(INFO)<<"Frame time "<<t2-t1;
                 flag_firstFrame_accomplished=TRUE;
            }
            else if (flag_firstFrame_accomplished==TRUE)
            {
//               
             nextframe_start_index=currframe_end_index;
           
             nFrameDiscard_fi=0;//bcoz already deleted in first iteration.
             
            status = performFrameIntegration ((long)nrows , nextframe_start_index,frame_no , xsize , nFrameDiscard_fi , nFrameIntegrate_fi , new_tempX , new_tempY , mult_temp , badFlag_temp , time_frame , effective_NumPhotons , one_dim_img_final , one_dim_exp_final , frameIntegration_track,frame_Data,frame_ExpData,frm_time_fi_new,frmno_fi_new,currframe_end_index) ;
            if(lastFrame_flag==TRUE){
                LOG(INFO)<<"out of method";
                break;
            }
            if(wtd_fi)
            {
             status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "SignalFrames" , "sig" , frame_Data , nameprefix ,frm_time_fi_new ,num_frm_count , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output
             if (status)
            {
                      LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                      continue;
             }
              status = writeOutputImageToDisk ("fi" , moduleoutdir_fi , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frm_time_fi_new ,num_frm_count , IMG_DIM_FI , IMG_DIM_FI) ; //this is for the SignalFrame output
              if (status)
              {
                       LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                       continue;
             }
                             
            }
                
                
                
            }
            
            
            
            rms_mul_factor = rms_mul_factor_default ;

//            for (int index = 0 ; index < IMG_DIM_FI * IMG_DIM_FI ; index ++)
//            {
//                frame_Data[index] = frameIntegration_track[i].pixels_sig_fi[index] ;
//            }
            
            status = obj_sc.getStarCentroid (frame_Data ,nFrameIntegrate_fi,backgrd_fact ,rms_mul_factor , IMG_DIM_FI , IMG_DIM_FI , refine_Winsize * (IMG_DIM_FI / 600) ,
                    centroid_Winsize * (IMG_DIM_FI / 600)  , refined_peaks_x_vect , refined_peaks_y_vect , centroid_x_vect , centroid_y_vect , firstcut_peaks_x_vect , firstcut_peaks_y_vect , R_int ,
                    c_int , F_int , datainfo.getObsMode () , minimum_No_of_Stars , flag_fsc , win_size,star_detect_algo_flag) ;
            if (status)
            {
                LOG (ERROR) << "***Error in star detection algorithm***" ;
                LOG(ERROR)<<"CRASH NO STAR FOUND (POSSIBLY BECAUSE SIGMA MULTIPLIER =< 0) (uvtRelativeAspectPC.cpp)";
                break;
            }
            
                       
            vector<float> fpeak_x,fpeaks_y,fpeaks_int;
         
            int maxele=MAX_FIRSTCUTPIX_CMPR;;
            if(R_int.size ()<MAX_FIRSTCUTPIX_CMPR)
            {
                maxele=R_int.size ();
               
            }
            for(int i=0;i<maxele;i++)
            {
                fpeak_x.push_back (firstcut_peaks_x_vect[i]);
                fpeaks_y.push_back (firstcut_peaks_y_vect[i]);
                fpeaks_int.push_back (F_int[i]);
            }
          //fpeak_x=(float)firstcut_peaks_x_vect;fpeaks_y=(float)firstcut_peaks_y_vect;fpeaks_int=(float)F_int;
            if (wtd_fsc == 1)
            {
                status = writeOutputTblToDisk ("sc" , moduleoutdir_sc , "Star" , "Star" , nameprefix , frm_time_fi_new , i + 1 , ttype_3 , tform2 , tfields , "First-cut peaks" , fpeak_x , fpeaks_y ,fpeaks_int) ;
               //status = writeOutputImageToDisk ("sc" , moduleoutdir_sc , "Star" , "Star" , frame_Data , nameprefix , frameIntegration_track[i].frame_time_fi , i + 1 , IMG_DIM_FI , IMG_DIM_FI) ;
               status = writeOutputImageToDisk ("sc" , moduleoutdir_sc , "Star" , "Star" , frame_Data , nameprefix , frm_time_fi_new, i + 1 , xsize , ysize) ;
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***"  ;
                    continue ;
                }
                status = writeOutputTblToDisk ("sc" , moduleoutdir_sc , "Centroid" , "centroid" , nameprefix , frm_time_fi_new , i + 1 , ttype_3 , tform2 , tfields , "Centroids" , centroid_x_vect , centroid_y_vect , c_int) ;
            }

            //centroid_x_vect.resize (20);
         //   centroid_y_vect.resize (20);
            //c_int.resize (20);
            firstcut_peaks_x_vect.clear () ;
            firstcut_peaks_y_vect.clear () ;
            F_int.clear () ;
            refined_peaks_x_vect.clear () ;
            refined_peaks_y_vect.clear () ;
            R_int.clear () ;

            //bool flag_rfc = FALSE ;
            
          
            if (cnt_forRefFrmCal >= frames_toDiscard)
            {

              //  if ((trackcnt_forRefFrame % nFrameToAverage != 0 || (nFrameToAverage == 1) || trackcnt_forRefFrame >= 0))
              //  {
                    size_rows.push_back (centroid_x_vect.size ()) ;
                    ref_frame_time_data.push_back ((frm_time_fi_new)) ;


                    for (int p = 0 ; p < centroid_x_vect.size () ; p ++)
                    {
                        x_Of_refFrame.push_back (centroid_x_vect[p]) ;
                        y_Of_refFrame.push_back (centroid_y_vect[p]) ;
                        int_Of_refFrame.push_back (c_int[p]) ;
                        //flag = TRUE ;
                    }
                    cnt_temp ++ ;
             //   }
                if (cnt_temp >= nFrameToAverage)
                {
                    if(flag_rfc==TRUE) nFrameToAverage=1;
                    status = calculateRefFrame (size_rows , x_Of_refFrame , y_Of_refFrame , int_Of_refFrame , nFrameToAverage , x_ref_cumm , y_ref_cumm , int_ref_cumm) ;
                    flag_rfc=TRUE;
                    centroid_x_vect = x_ref_cumm ;
                    centroid_y_vect = y_ref_cumm ;
                    c_int = int_ref_cumm ;
                    //clearing all the memory of the VECTORS of Reference Frame.
                    x_Of_refFrame.clear () ;
                    y_Of_refFrame.clear () ;
                    int_Of_refFrame.clear () ;
                    size_rows.clear () ;
                    x_ref_cumm.clear () ;
                    y_ref_cumm.clear () ;
                    int_ref_cumm.clear () ;
                    time_RefFrame_Avg = (ref_frame_time_data[0] + ref_frame_time_data[ref_frame_time_data.size () - 1]) / 2 ;
                    //status = writeOutputTblToDisk ("rf" , moduleoutdir_rfc , "Centroid" , "centroid" , nameprefix , time_RefFrame_Avg , cnt_forRefFrmCal-nFrameToAverage+1 , ttype_3 , tform2 , tfields , "Centroids" , centroid_x_vect , centroid_y_vect , c_int) ;
                    status = writeOutputTblToDisk ("rf" , moduleoutdir_rfc , "Centroid" , "centroid" , nameprefix , time_RefFrame_Avg , trackcnt_forRefFrame+1 , ttype_3 , tform2 , tfields , "Centroids" , centroid_x_vect , centroid_y_vect , c_int) ;
                    if (status)
                    {
                        LOG (ERROR) << "Error in writing the output to the Disk" ;
                        return (EXIT_FAILURE) ;
                    }
                    time_RefFrame_Avg = 0.0 ;
                    ref_frame_time_data.clear () ;
                   // cnt_temp = 0 ;
                    trackcnt_forRefFrame ++ ;
                }
                //trackcnt_forRefFrame ++ ;
            }
            cnt_forRefFrmCal = cnt_forRefFrmCal + 1 ;
            centroid_x_vect.clear () ;
            centroid_y_vect.clear () ;
            c_int.clear () ;
            LOG (INFO) << "\033[1;34mProcess ends for  frame number -> " << cnt_finish << "\033[0m" ;
            cnt_finish ++ ;
            
            

        }

        char infofile_out[NAMESIZE] ;
        sprintf (infofile_out , "%s/%s_rfc.info" , moduleoutdir_rfc , nameprefix) ;
        fits_create_file (&finfo_out , infofile_out , &status) ;
        printError (status , "Error in creating the output information file") ;
        char *ttype1[] = {"centroidFileList"} ;
        char *tform1[] = {"A256"} ;
        fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype1 , tform1 , NULL , "FileList" , &status) ;
        printError (status , "Error in Creating the table ***") ; //for creating name for output information file
        fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
        printError (status , "Error in updating the key value of the NAMEPRFX***") ; //for creating name for output information file
        char dir[100] = "Centroid" ;
        fits_update_key (finfo_out , TSTRING , "CENTDIR" , (char*) dir , NULL , &status) ;
        printError (status , "Error in updating the  key value of the centroidDir") ;
        fits_update_key (finfo_out , TSTRING , "OBS_MODE" , obsmode , NULL , &status) ;
        printError (status , "Error in updating the  key value of the OBS_MODE") ;
        fits_update_key (finfo_out , TDOUBLE , "FRMTIME" , &time_RefFrame_Avg , NULL , &status) ;
        printError (status , "Error in updating the  key value of the OBS_MODE") ;

        int size_vect = ref_frame_module_filename_track.size () ;
        string *name_ofFile = new string[size_vect] ;
        for (int p = 0 ; p < size_vect ; p ++)
        {
            name_ofFile[p] = ref_frame_module_filename_track[p];
            //LOG(INFO)<<name_ofFile[p];
        }
        fits_update_key (finfo_out , TINT , "NFILES" , &size_vect , NULL , &status) ;
        printError (status , "Error in updating the key value of the NFILES") ;
        xsize = ysize = IMG_DIM_FI ;
        datainfo.write (finfo_out) ;
        fits_update_key (finfo_out , TINT , "XSIZE" , &xsize , NULL , &status) ;
        printError (status , "Error in updating the  key value of the XSIZE") ;
        fits_update_key (finfo_out , TINT , "YSIZE" , &ysize , NULL , &status) ;
        printError (status , "Error in updating the  key value of the YSIZE ") ;
for(int i =0;i<size_vect;i++)
{
        fits_write_col (finfo_out , TSTRING , 1 , i+1 , 1 , 1 , &name_ofFile[i] , &status) ;
        printError (status , "Error in writing the column of outout Exposure frame list") ;
}

//        fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , size_vect , name_ofFile , &status) ;
  //      printError (status , "Error in writing the colucdmn of outout Exposure frame list") ;
        fits_close_file (finfo_out , &status) ;
        printError (status , "Error in Closing the output information file") ;
        // delete[] flatfielddata;

        LOG (INFO) << endl << "==================================UVTCOMPUTEDRIFT==================================================" << endl ;
        uvtDriftComputation drift_obj ;
        drift_obj.read (moduleoutdir_rfc ,(char*)dirobj.attfile.c_str (), err_per , pair_nbhd_distance , FreqDomainFilterFlag , type_Filtering , poly_fit_interval , freqvalue , fitting_flag , orderpitch , orderyaw , orderroll , outputdir , shift_rotation_algo , match_stars_file_flag , star_detect_algo_flag,flag_thetaComp,clobber , history) ;
        //  dex_obj.read (moduleIndir ,err_per , diff_Dist , freqDomainFilter_Flag , type_Filtering , delta_time , freqvalue , fitting_flag , orderpitch , orderyaw , orderroll , outputdir , option_LeastSquare , match_stars_file_flag,clobber , history) ;
        drift_obj.display () ;
        status = drift_obj.uvtDriftComputationProcess () ;
        if (status)
        {
            LOG (ERROR) << endl << "***Error in DriftComputation  Module***" << endl ;
	LOG(ERROR)<<"CRASH DRIFT COMPUTATION FAILED (uvtRelativeAspectPC.cpp)";
            continue;
        }
        strcpy (moduleIndir , drift_obj.getModuleOutdir ()) ;
        strcpy (drift_outDir , moduleIndir) ;
//        LOG (INFO) << endl << "==================================UVTCOMPUTEJITTER==================================================" << endl ;
//        uvtComputeJitter jitter_obj ;
//        jitter_obj.read (moduleIndir , (char*) caldbindir.c_str () , (char *) dirobj.gyrofile.c_str () , FreqDomainFilterFlag , freqvalue , fitting_flag , orderpitch , orderyaw , orderroll , outputdir , type_Filtering , clobber , history) ;
//        jitter_obj.display () ;
//        status = jitter_obj.uvtComputeJitterProcess () ;
//        if (status)
//        {
//            LOG (INFO) << endl << "\033[1;31m***Error in ComputeJitter  module***" << endl ;
//              return(EXIT_FAILURE);
//        }
//        strcpy (moduleIndir , jitter_obj.getModuleOutdir ()) ;
//        strcpy (jitter_outDir , moduleIndir) ;
//        LOG (INFO) << "Using directory " << moduleIndir << "  as input to uvtComputeThermalr "  ;
//
//        //----------------------THERMAL SERIES-------------//
//          if (strcasecmp (datainfo.getDetector () , "VIS") != 0){
//        LOG (INFO) << endl << "\033[1;34m==================================UVTCOMPUTETHERMAL==================================================\033[0m"  ;
//
//        uvtComputeThermal  thermal_obj ;
//        thermal_obj.read (moduleIndir , moduleoutdir_rfc , (char*) caldbindir.c_str () , (char*) dirobj.lbtfile.c_str () , outputdir , clobber , history) ;
//        thermal_obj.display () ;
//        status = thermal_obj.uvtThermalCalcProcess () ;
//        if (status)
//        {
//            LOG (INFO) << endl << "\033[1;31m***Error in uvtComputeThermal  module***" ;
//              return(EXIT_FAILURE);
//        }
//        strcpy (moduleIndir , thermal_obj.getModuleOutdir ()) ;
//
//        LOG (INFO) << "Using directory " << moduleIndir << "  as input to uvtComputeThermal "  ;
//          }
//
//        //        
//        //        //----------------------RELATIVE ASPECT CALCULATION--------------//
//        LOG (INFO) << endl << "\033[1;34m==================================UVTCOMPUTERELASPECT==================================================\033[0m"  ;
//        uvtRelAspCal  ras_obj ;
//        ras_obj.read (drift_outDir , jitter_outDir , moduleIndir , outputdir , clobber , history) ;
//        ras_obj.display () ;
//        status = ras_obj.uvtRelAspCalProcess () ;
//        if (status)
//        {
//            LOG (INFO) << endl << "\033[1;31m***Error in uvtRelAspCal  module***" << endl ;
//             return(EXIT_FAILURE);
//        }
//        strcpy (RelAspFilename , ras_obj.getModuleOutdir ()) ;


    }

    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::setDirectoryStructure (char *Dir , const char *subdir)
{
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s/" , Dir , subdir) ;
    //cout<<dir<<endl;
    string cmd ;
    system (cmd.c_str ()) ;
    if (DirExists (dir) && clobber == YES)
    {
        LOG (ERROR) << "Directory exists and clobber=yes" ;
        cmd = (string) "rm -rf " + (string) dir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists (dir) && clobber == NO)
    {
        LOG (ERROR) << endl << dir << "  already exists " ;
        LOG (ERROR) << endl << "Use clobber=yes for overwriting" ;
        return (EXIT_FAILURE) ;
    }
    cmd = "mkdir -p " + (string) dir ;
    system (cmd.c_str ()) ;

    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::writeOutputImageToDisk (char *id , char *outDir , char *dir , char *subscript , float *Array , char *namepre , double ftime , unsigned short fno , int sizex , int sizey)
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
//    LOG(INFO)<<setprecision(20)<<ftime;exit(1);
    sprintf (outfile , "%s/%s/%s_t%.4f_f%d_%s_%s.fits" , outDir , dir , namepre , ftime , fno , id , subscript) ;
    int numhdu = 0 ;
    fits_create_file (&fout , outfile , &status) ;
    printError (status , "Error in creating the output Signal File" , outfile) ;
    fits_create_img (fout , bitpix , naxis , naxes , &status) ;
    printError (status , "Error in Creating the image for Signal Fie" , outfile) ;
    fits_write_pix (fout , TFLOAT , fpixel , sizex*sizey , Array , &status) ;
    printError (status , "***Error in writing the pixels to output***" , outfile) ;
    fits_get_num_hdus (fout , &numhdu , &status) ;
    char temp[1000] ;
    writeCommonKeywords (fout , outDir) ;
    for (int i = 0 ; i < numhdu ; i ++)
    {
        for (int p = 0 ; p < header_info.size () ; p ++)
        {
            sprintf (temp , "%s" , header_info[p].c_str ()) ;
            fits_write_record (fout , (char*) &temp , &status) ;
        }
    }
    fits_close_file (fout , &status) ;
    printError (status , "Error in closing the  output Signal fits file" , outfile) ;

    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::copyAllheaderKeys (char* infile)
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
        cout << endl << "***Could not find number of keywords in file " << infile << "***" ;
        fits_report_error (stderr , status) ;
        return (EXIT_FAILURE) ;
    }
    int keyclass ;
    //cout<<"\nNumber of keywords found:"<<keyexist;
    for (int i = 1 ; i <= keyexist ; i ++)
    {
        fits_read_record (fin , i , record , &status) ;
        if (status)
        {
            cout << endl << "***Error in reading record number " << i << " in file " << infile << "***" ;
            fits_report_error (stderr , status) ;
            return (EXIT_FAILURE) ;
        }
        keyclass = fits_get_keyclass (record) ;
        if (keyclass == TYP_COMM_KEY)
            continue ;
        else if (keyclass == TYP_USER_KEY || keyclass == TYP_REFSYS_KEY || keyclass == TYP_WCS_KEY)
        {
            header_info.push_back (record) ;
        }
    }

    fits_close_file (fin , &status) ;
    // delete[] record;
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::transformToUVITFrame (double *t , double *r , double *p , double *y)
{
    cout << "the 1,1" << endl ;
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::transformToSpacecraftFrame (double *t , double *r , double *p , double *y)
{
    cout << "the 1,1" << endl ;
    return (EXIT_SUCCESS) ;
}


double uvtRelativeAspectPC::readDarkFrame (char * path , float *Array)
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


int uvtRelativeAspectPC::darkFrameComputation (float *outArray)
{
    for (int i = 0 ; i < xsize * ysize ; i ++)
    {
        outArray[i] = 0.0 ;
        float d1 = darkFramestart_data[i] ;
        float d2 = darkFrameend_data[i] ;
        double t1 = t_darkframestart ;
        double t2 = t_darkframeend ;
        //    cout<<"D1::"<<d1<<endl;
        //    cout<<"D2::"<<d2<<endl;
        //    cout<<"t1::"<<t1<<endl;
        //    cout<<"t2::"<<t2<<endl;

        if ((t2 - t1) == 0)
        {
            LOG (ERROR) << "***Divide By zero Error***" << endl ;
            return (EXIT_FAILURE) ;
        }
        float d = d1 + ((d2 - d1) *((t2 - t1) / (t_curr - t1))) ;
        //cout<<d<<endl;

        outArray[i] = (float) d ;
        // outArray[i] = d1;

    }
    return (EXIT_SUCCESS) ;

}


int uvtRelativeAspectPC::darkFrameSubtraction (float *Array , float *frame_data)
{

    for (int p = 0 ; p < xsize * ysize ; p ++)
    {
        frame_data[p] = frame_data[p] - Array[p] ;
        if (frame_data[p] < 0)
        {
            frame_data[p] = 0.0 ;
        }
    }


    return (EXIT_SUCCESS) ;


}


int uvtRelativeAspectPC::readImage (char * caldb_file , int hduno , float *frm_data)
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


int uvtRelativeAspectPC::performUnitConversion (float *frmdata , float *expdata , double  intgrntime , int sizex , int sizey)
{

    if (intgrntime == 0)
    {
        LOG (ERROR) << "***Divide by Zero***" << endl ;
        return (EXIT_FAILURE) ;

    }
    for (int i = 0 ; i < sizex * sizey ; i ++)
    {
        frmdata[i] = frmdata[i] / intgrntime ;

    }
    return (EXIT_SUCCESS) ;
}
//int uvtPcRaCommanArr::performCorrection(float *frmsigdata, float *frmexpdata, float *badpixarry,int size,double  intgrntime){
//    
//     for (int pixno = 0 ; pixno <size  ; pixno++)
//            frmsigdata[pixno] = frmsigdata[pixno] * badpixarry[pixno] ;
//
//        for (int pix = 0 ; pix < size; pix++)
//        {
//            frmexpdata[pix] = 0.0f ;
//            frmexpdata[pix] = badpixarry[pix] * intgrntime ;
//        }
//    return(EXIT_SUCCESS);
//}


int uvtRelativeAspectPC::performCosmicRayCorr (float *frmsigdata , float *frmexpdata , int sizex , int sizey , float threshold_cr)
{
    vector<int> xpix , ypix ;
    int cnt_cosmicAffected = 0 ;
    for (int pixno = 0 ; pixno < sizex * sizey ; pixno ++)
    {
        if (frmsigdata[pixno] > threshold_cr)
        {
            xpix.push_back (pixno % sizey) ;
            ypix.push_back (pixno / sizex) ;
            //frmsigdata[pixno] = frmexpdata[pixno] = 0.0 ;
        }

    }

    for (int i = 0 ; i < xpix.size () ; i ++)
    {
        cnt_cosmicAffected = 0 ;

        for (int j = xpix[i] - 1 ; j <= xpix[i] + 1 ; j ++)
        {
            for (int k = ypix[i] - 1 ; k <= ypix[i] + 1 ; k ++)
            {

                if (j < sizex && j > 0 && k < sizey && k > 0)
                {
                    if (frmsigdata[k * sizey + j] > threshold_cr)
                    {
                        cnt_cosmicAffected ++ ;
                    }

                }

            }

        }
        if (cnt_cosmicAffected == 1)
        {
            frmexpdata[ypix[i] * sizex + xpix[i]] = frmsigdata[ypix[i] * sizey + xpix[i]] = 0.0f ; //cr effected
              }
    }



    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::performAccOFFrame (float *frmsigdata , float *sumdata , int sizex , int sizey , int numoffrmAcc)
{
    if (numoffrmAcc == 0)
    {
        cout << "Divide by Zero" << endl ;
        return (EXIT_FAILURE) ;
    }
    for (int pix = 0 ; pix < sizex * sizey ; pix ++)
    {
        frmsigdata[pix] = 0.0 ;
        frmsigdata[pix] = sumdata[pix] / numoffrmAcc ;
        sumdata[pix] = 0.0 ;

    }

    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::writeOutputTblToDisk (char *id , char *outDir , char *dir , char *subscript  , char *namepre , double ftime , unsigned short fno , char **type1 , char**tform1 , int tfields1 , char *extname , vector<float> &X , vector<float> &Y , vector<float> &val)
{
    char out_file[FLEN_FILENAME] ;
    sprintf (out_file , "%s/%s/%s_t%f_f%d_%s_%s.fits" , outDir , dir , namepre , ftime , fno , id , subscript ) ;
    if (strcmp (id , "rf") == 0)
    {
        ref_frame_module_filename_track.push_back (basename (out_file)) ;
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


int uvtRelativeAspectPC::performDistortionCorr (vector<float> &X , vector<float> &Y , float *Xdistr , float *Ydistr , int sizex , long caldbsize)
{
    float multi_factor = sizex / 600 ;
    if (multi_factor == 0)
    {
        LOG (ERROR) << "***Divide by Zero***" << endl ;
        return (EXIT_FAILURE) ;
    }
    float tempx , tempy ;
    for (int p = 0 ; p < X.size () ; p ++)
    {

        tempx = X[p] / multi_factor ;
        tempy = Y[p] / multi_factor ;

        float  locate = ((int) round (tempy) - 44) * caldbsize + ((int) round (tempx) - 44) ;
        round (locate) ;

        tempx = tempx + Xdistr[(int) locate] ;
        tempy = tempy + Ydistr[(int) locate] ;


        X[p] = tempx * multi_factor ;
        Y[p] = tempy * multi_factor ;


    }
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::calculateRefFrame (vector<int> &sizevect , vector<float> &xref , vector<float> &yref , vector<float> &intref , int avgfact , vector<float> &xrefcumm , vector<float> &yrefcumm , vector<float> &intrefcumm)
{
    double temp_x1 , temp_y1 , temp_int1 , temp_x2 , temp_y2 , temp_int2 ;
    float x_fref[xref.size ()] , y_fref[xref.size ()] , xrefarr[xref.size ()] , yrefarr[xref.size ()] , intrefarr[xref.size ()] ;
    for (int d1 = 0 ; d1 < sizevect[0] ; d1 ++)
    {
        int total = sizevect[0] ;
        int cnt_ref = 0 ;
        temp_x1 = xref[d1] ;
        temp_y1 = yref[d1] ;
        temp_int1 = intref[d1] ;
        // cout<<size_rows[0]<<" "<<size_rows[1]<<" "<<size_rows[2]<<" "<<endl;

        for (int d2 = 1 ; d2 < avgfact ; d2 ++)
        {
            total = total + sizevect[d2] ;
            for (int d3 = total - sizevect[d2] ; d3 < total ; d3 ++)
            {
                double diff_x = 0 , diff_y = 0 ;
                temp_x2 = xref[d3] ;
                temp_y2 = yref[d3] ;
                temp_int2 = intref[d3] ;
                //cout<<d3 << "  " <<temp_x2<<"  "<<temp_y2<<endl;
                diff_x = temp_x1 - temp_x2 ;
                //cout<<"the "<<diff_x<<endl;
                diff_y = temp_y1 - temp_y2 ;

                if (sqrt ((diff_x * diff_x)+(diff_y * diff_y)) < NBHD_RADIUS)
                {
                    x_fref[cnt_ref] = temp_x1 ;
                    y_fref[cnt_ref] = temp_y1 ;
                    xrefarr[cnt_ref] = temp_x2 ;
                    yrefarr[cnt_ref] = temp_y2 ;
                    intrefarr[cnt_ref] = temp_int2 ;
                    cnt_ref ++ ;
                }
            }
        }

        float x = temp_x1 ;
        float y = temp_y1 ;
        double inten = temp_int1 ;
        //  cout<<"outside of the Calcu"<<endl;
        if (cnt_ref >= avgfact - 1)
        {
            for (int d4 = 0 ; d4 < cnt_ref ; d4 ++)
            {
                x = x + xrefarr[d4] ;
                y = y + yrefarr[d4] ;
                inten = inten + intrefarr[d4] ;
            }

            x = x / avgfact ;
            y = y / avgfact ;
            inten = inten / avgfact ;
            xrefcumm.push_back (x) ;
            yrefcumm.push_back (y) ;
            intrefcumm.push_back (inten) ;
        }
    }
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::takeDarkinfo ()
{

    double dark_Firstframetime , dark_Secframetime , start_i , start_f , stop_i , stop_f ;
    char darkpath_First[FLEN_FILENAME] ;
    char darkend_Sec[FLEN_FILENAME] ;
    char dark_temp[FLEN_FILENAME] ;

    //vector for storing the dark frame names
    vector<string> dark_framenames ;

    //setting path for dark frame directory  
    sprintf (dark_temp , "%s/%s" , Indir_dataIngest , darkdir) ;

    getFileNamesfrmDir (dark_temp , ".fits" , dark_framenames) ; //get dark file names from dark directory

    if (dark_framenames.size () < 2)
    {
        LOG (INFO) << "No enough Dark frames found at input Directory" << endl ;
        return (EXIT_SUCCESS) ;
    }

    //setting first darkframe's path
    //for dark frame 1
    sprintf (darkpath_First , "%s/%s" , dark_temp , dark_framenames[0].c_str ()) ;

    readKeywords (darkpath_First , 1 , 1 , TDOUBLE , "FRMTIME" , &dark_Firstframetime ) ;

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


int uvtRelativeAspectPC::performQEMCPcorrection (float *sigdata , int size , double fact)
{

    for (int i = 0 ; i < size ; i ++)
    {
        sigdata[i] = sigdata[i] * fact ;
    }


    return (EXIT_SUCCESS) ;
}
//int uvtPcRaCommanArr::readQEMCPFile (){
//    
//    
//    
//}


int uvtRelativeAspectPC::performMaskBadPix (unsigned short* int_x , unsigned short* int_y , float* badpixArr , unsigned char * Max_Min , unsigned char* badflag , unsigned char* multflag , int nrows , int x_size , int y_size , float thr_multphn)
{
//    for (int i = 0 ; i < nrows ; i ++)
//    {
//        int *fg = 0 ;
//        badflag[i] = static_cast<unsigned char> (badpixArr[int_y[i] * x_size + int_x[i]]) ;
//        if ((float) Max_Min[i] > thr_multphn) multflag[i] = 0 ; //if max-min > user_threshold then flag=0 else flag=1
//        else multflag[i] = 1 ;
//
//    }
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

    return Total_MultPhoton_effected_frames;
    
}


int uvtRelativeAspectPC::performPixPadding (unsigned short *X_int , unsigned short *Y_int , int nrows)
{
    for (int i = 0 ; i < nrows ; i ++)
    {
        X_int[i] = X_int[i] + 44 ;
        Y_int[i] = Y_int[i] + 44 ;
    }
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::performCentroidCorr (double *t , unsigned short *xint , unsigned short *yint , float *xf , float *yf , unsigned char *mc , float *newXfract , float *newYfract , float *darkbeginData , float *darkenddata , double integrtn_time , long EA , int nrows)
{
    float dark_value = 0.0 ;
    float sd = 0.0 , sc = 0 ;
    double t1 , t2 ;
    t2 = t[nrows - 1] + integrtn_time ;
    int index_x ;
    int index_y ;
    int diff_add = 600 - darkSize ;
    
    for (int i = 0 ; i < nrows ; i ++)
    {
        
        dark_value = 0.0f ;
        t1 = t[i] - integrtn_time ;
        if ((t2 - t1) == 0)
        {
            LOG (ERROR) << "***Divide by Zero***"<<integrtn_time ;
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
        newYfract[i] = yf[i]*(1 + ((sd - sc) * cent_corr_win_size) / EA) + index_y ;
    }

    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::readcentroidbiasFile ()
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
    fraction_bias = new double[nrows]  , x_corr = new double[nrows] , y_corr = new double[nrows] ;
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


int uvtRelativeAspectPC::performCentroidBias (long nrows , float *xFrac , float *yFrac , float *new_xFrac , float *new_yFrac)
{
    float *xFrac_temp = new float[nrows] ;
float *yFrac_temp = new float[nrows] ;
int *New_X_Val= new int [nrows];
int *New_Y_Val= new int [nrows];
int temp_Xval,temp_Yval;
float tempXX,tempYY;
bool flag_match_found_X=FALSE;
bool flag_match_found_Y=FALSE;
for (int i=0;i<nrows;i++)
{
    LOG(INFO)<<(xFrac[i]-(int)xFrac[i])*32<<" "<<yFrac_temp[i];
    xFrac_temp[i] =( (xFrac[i]-(int) xFrac[i])*32);
    yFrac_temp[i] =( (yFrac[i]-(int) yFrac[i])*32);
    
    tempXX=xFrac_temp[i]/32;
    tempYY=yFrac_temp[i]/32;
    New_X_Val[i]=((int) xFrac[i])*32;
    New_Y_Val[i]=((int) yFrac[i])*32;
    
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
    
    
    //FPN correction
    for (int j = 0 ; j < biasRows ; j ++)
       {
      //  LOG(INFO)<<(x_corr[j] + ((xFrac_temp[i] - fraction_bias[j]) * (x_corr[j + 1] - x_corr[j])) / (fraction_bias[j + 1] - fraction_bias[j]));
              if (xFrac_temp[i]*4096/32 >= fraction_bias[j] && xFrac_temp[i]*4096/32 < fraction_bias[j + 1])
              {
                       xFrac_temp[i] = xFrac_temp[i]*((x_corr[j]) + ((xFrac_temp[i]*4096/32 - fraction_bias[j]) * (x_corr[j + 1] - x_corr[j])) / (fraction_bias[j + 1] - fraction_bias[j])) ; 
                       flag_match_found_X=TRUE;
                       break;
              }
              
       }
   // exit(1);
     for (int j = 0 ; j < biasRows ; j ++)
       {
              if (yFrac_temp[i]*4096/32 >= fraction_bias[j] && yFrac_temp[i]*4096/32 < fraction_bias[j + 1])
              {
                       yFrac_temp[i] =yFrac_temp[i]*(( y_corr[j] )+ ((yFrac_temp[i]*4096/32 - fraction_bias[j]) * (y_corr[j + 1] - y_corr[j])) / (fraction_bias[j + 1] - fraction_bias[j])); 
                       flag_match_found_Y=TRUE;
                       break;
              }              
       }
    if(flag_match_found_X ==FALSE || flag_match_found_Y== FALSE)
    {
        //LOG(INFO)<<i<<" Not Found"<<" "<<xFrac_temp[i]<<" "<<yFrac_temp[i];
        
    }
        //LOG(INFO)<<xFrac_temp[i]<<" "<<yFrac_temp[i];
    //For reverse conversion
    xFrac_temp[i]=xFrac_temp[i]-16;
    yFrac_temp[i]=yFrac_temp[i]-16;
//    new_xFrac[i]=xFrac_temp[i]+New_X_Val[i];
//    new_yFrac[i]=yFrac_temp[i]+New_Y_Val[i];
    
    new_xFrac[i]=temp_Xval+xFrac_temp[i]/32;
    new_yFrac[i]=temp_Yval+yFrac_temp[i]/32;
    LOG(INFO)<<new_xFrac[i]<<" "<<new_yFrac[i]<<" "<<tempXX<<" "<<tempYY;
}
return(EXIT_SUCCESS);
//
//    float *xFrac_temp = new float[nrows] ;
//    float *yFrac_temp = new float[nrows] ;
//    for (int i = 0 ; i < nrows ; i ++)
//    {
//        xFrac_temp[i] = ( (xFrac[i]-(int) xFrac[i]) + 32) / 2 ;
//        yFrac_temp[i] = ( (yFrac[i]-(int) yFrac[i]) + 32) / 2 ;
//        //cout<<xFrac_temp[i]<<" "<<yFrac_temp[i]<<endl;
//    }
//
//    int x_sign_flag = 0 ;
//    int y_sign_flag = 0 ;
//    for (int i = 0 ; i < nrows ; i ++)
//    {
//        if (xFrac[i] >= 0)
//            x_sign_flag = 1 ;
//        else
//            x_sign_flag = 0 ;
//
//        if (yFrac[i] >= 0)
//            y_sign_flag = 1 ;
//        else
//            y_sign_flag = 0 ;
//
//        //       xFrac[i] = fabs (xFrac[i]) ; //*32.0;
//        //      yFrac[i] = fabs (yFrac[i]) ; //*32.0;
//        xFrac[i] = fabs (xFrac[i]) ;
//        yFrac[i] = fabs (yFrac[i]) ;
//        /*Searching for X Frac value in the Centroid Bias correction file in CalDB**/
//        if (xFrac_temp[i] == fraction_bias[biasRows - 1])
//        {
//            new_xFrac[i] = x_corr[biasRows - 1] ;
//            if (x_sign_flag == 0)
//                new_xFrac[i] = - new_xFrac[i] ;
//
//            if (new_xFrac[i] < 600)
//                new_xFrac[i] = new_xFrac[i] + (int) xFrac[i] ;
//        }
//        else
//        {
//            for (int j = 0 ; j < biasRows ; j ++)
//            {
//                if (xFrac_temp[i] >= fraction_bias[j] && xFrac_temp[i] < fraction_bias[j + 1])
//                {
//                    new_xFrac[i] = x_corr[j] + (((xFrac_temp[i]) - fraction_bias[j]) * (x_corr[j + 1] - x_corr[j])) / (fraction_bias[j + 1] - fraction_bias[j]) ;
//                    if (x_sign_flag == 0)
//                        new_xFrac[i] = - new_xFrac[i] ;
//                    if (new_xFrac[i] < 600)
//                        new_xFrac[i] = new_xFrac[i] + (int) xFrac[i] ;
//
//                    break ;
//                }
//            }
//
//        }
//        /***Searching for Y Frac value in the Centroid Bias correction file in CalDB***/
//        if (yFrac_temp[i] == fraction_bias[biasRows - 1])
//        {
//
//            new_yFrac[i] = y_corr[biasRows - 1] ;
//            if (y_sign_flag == 0)
//                new_yFrac[i] = - new_yFrac[i] ;
//
//            if (new_yFrac[i] < 600)
//                new_yFrac[i] = new_yFrac[i] + (int) yFrac[i] ;
//        }
//        else
//        {
//            for (int j = 0 ; j < biasRows ; j ++)
//            {
//
//                if (yFrac_temp[i] >= fraction_bias[j] && yFrac_temp[i] < fraction_bias[j + 1])
//                {
//                    new_yFrac[i] = y_corr[j] + ((yFrac_temp[i] - fraction_bias[j]) * (y_corr[j + 1] - y_corr[j])) / (fraction_bias[j + 1] - fraction_bias[j]) ;
//                    if (y_sign_flag == 0)
//                        new_yFrac[i] = - new_yFrac[i] ;
//                    if (new_yFrac[i] < 600)
//                        new_yFrac[i] = new_yFrac[i] + (int) yFrac[i] ;
//                    break ;
//                }
//            }
//        }
//        //  new_xFrac[i]=(new_xFrac[i]*2-32)+(int)xFrac[i];
//        //  new_yFrac[i]=new_yFrac[i]*2-32+(int)yFrac[i];
//
//    }
//
//    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::performDistCorrection (long nrows , float *xFrac , float *yFrac , float *x_Distortion , float *y_Distortion , int  caldbdim)
{
 float locate ;
    float tempx , tempy , multi_factor = 0.0 ;
    multi_factor = xsize / 600 ;
    LOG (INFO) << "Applying Correction..." ;
    long curr_loc=0;
    long next_loc=0;
    double value_X_toAdd,value_Y_toAdd;
    for (int k = 0 ; k < nrows ; k ++)
    {
        tempx = (xFrac[k] / multi_factor)-44 ;
        tempy = (yFrac[k] / multi_factor)-44 ;
        if (xFrac[k] != INVALID_PIX_VALUE && yFrac[k] != INVALID_PIX_VALUE)
        {
       // locate = ((int) round (tempy) - 44) * caldbdim + ((int) round (tempx) - 44) ;
        locate = ((int) round (tempy)) * caldbdim + ((int) round (tempx)) ;
        //locate = (tempy) * caldbsize + (tempx) ;
       
               curr_loc=(int)locate-1;
               next_loc=curr_loc+1;
                if(curr_loc>0 && next_loc<512*512 && locate > 0 && locate <512*512){
               value_X_toAdd=x_Distortion[curr_loc]+((x_Distortion[next_loc]-x_Distortion[curr_loc]))*(locate-curr_loc);
               value_Y_toAdd=y_Distortion[curr_loc]+((y_Distortion[next_loc]-y_Distortion[curr_loc]))*(locate-curr_loc);
   
            tempx = tempx - (value_X_toAdd) ;
        tempy = tempy - (value_Y_toAdd) ;
       
        if (tempx * multi_factor < xsize && tempy * multi_factor < ysize && tempx * multi_factor >0 && tempy * multi_factor>0 )
        {
            xFrac[k] = (tempx +44)* multi_factor ;
            yFrac[k] = (tempy +44)* multi_factor ;
        }    
        
        else{
            xFrac[k]=INVALID_PIX_VALUE;
            yFrac[k]=INVALID_PIX_VALUE;
        }
    }
        else{
            xFrac[k]=INVALID_PIX_VALUE;
            yFrac[k]=INVALID_PIX_VALUE;
        }
    }
        else{
             xFrac[k]=INVALID_PIX_VALUE;
            yFrac[k]=INVALID_PIX_VALUE;
        }
    }
    return (EXIT_SUCCESS) ;

}

//int uvtRelativeAspectPC::performFrameIntegration (long nrows , long start_row,unsigned short*frame_no , int xsize , int Ndiscard , int Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr_strct> &vect,
//       float *outputSigArr,float *outputExpArray,double &outputFrmtime,long &outputFrmNo,long &last_index_for_frame)
//{
//    for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++) 
//    {
//        outputSigArr[i]=0.0f;
//        outputExpArray[i]=0.0f;
//    }
//    FrameIntegration_Arr_strct obj ;
//    int jj = 0 ;
//    int x_tem , y_tem ;
//    double Start_Time , End_Time , Avg_Time = 0.0 ;
//
//   // long firstpix1[2] ;
//   // firstpix1[0] = firstpix1[1] = 1 ;
//  //  int naxis = 2 ;
//   // long naxes[2] = {IMG_DIM_FI , IMG_DIM_FI} ;
//    float Eff_NoPhtn = 0.0f ;
//    float multi_factor = IMG_DIM_FI / xsize ;
//  //  int file_num = 1 ;
//    int frame_check ;
//    LOG (INFO) << "MULTIPLICATION FACTOR  for converting image to "<<IMG_DIM_FI<<" 9600 is " << multi_factor << endl ;
//    int frame_init = frame_check = Nacc ;
//    double time_frame = 0.0 , time_frame_final = 0.0 ;
//    int cnt_frame = 0 ;
//  //  vector<string> vhistorystr ;
//    
//    /****provision for if frames are not starting from 1****/
//   // int start_row = 0 ;
//    int cmpr_term = 0 ;
// 
//     if(Nacc>frame_no[nrows-1]){
//               LOG(ERROR)<<"***Number of frames to accumulate is greater than total available frames";
//                return(EXIT_FAILURE);
//    }
//    if (frame_no[0] != 1)
//    {
//        cmpr_term = Ndiscard + (frame_no[0] - 1) ;
//    }
////    for (int j = 0 ; j < nrows ; j ++)
////    {
////        if (frame_no[j] > cmpr_term)
////        {
////            start_row = j ;
////            break ;
////        }
////    }
//    float *badPix_Arr_padded= new float [600*600];
//    for(int i=0;i<600*600;i++) badPix_Arr_padded[i]=0.0;
//    float *badPix_Arr_fi= new float[IMG_DIM_FI*IMG_DIM_FI];
//      for(int i=0;i<600*600;i++) badPix_Arr_fi[i]=0.0;
//   int   status = Applypadding (badpixdata,512,512,badPix_Arr_padded,600,600);
//    if (status)
//    {
//        LOG (ERROR) << "Error in padding for bad pix data" ;
//        return (EXIT_FAILURE) ;
//    }
//    status= performSubDivisionIM (badPix_Arr_padded,600,600,badPix_Arr_fi,IMG_DIM_FI,IMG_DIM_FI);
//   if (status)
//    {
//        LOG (ERROR) << "Error in performing sub division in badpixel data" ;
//        return (EXIT_FAILURE) ;
//    } 
//    long int tot_num_frames = 0 ;
//    for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++)
//    {
//    one_dim_img[i]=0;
//    //cout<<badpixdata[i]<<endl;
//    one_dim_exp[i]=badPix_Arr_fi[i];
//    }
//    
//    LOG (INFO) << "Performing Frame Integration...." << endl ;
//    long p = 0 ;
//    
//    //calculate final frame count value for the current output frame.
//    int nrows_cur_frm=start_row;
//    for (int i=start_row;i<nrows;i++)
//    {
//        if (frame_no[i]<frame_no[start_row]+Ndiscard+Nacc)
//         {  
//         nrows_cur_frm++;
//  //LOG(INFO)<<i;
//          }
//    }
//    last_index_for_frame=nrows_cur_frm;
//  // LOG(INFO)<<"FFF "<<nrows<<" "<<nrows_cur_frm<<start_row;
//    //LOG(INFO)<<"PPP "<<start_row<<" "<<nrows_cur_frm<<" "<<frame_no[0]<<" "<<nrows;exit(1);
//    for (int k = start_row ; k < nrows_cur_frm ; k ++)
//    {
//        /**in 'IF condition' accumulated array of  Nacc frames will be  generate   & In 'ELSE condition' array will be written to the frame.
//         x and y location  for number of events are taken as a indexes for the array .For each Nacc events one array will be written as a frame. 
//         **/
//      
//        if (frame_no[k] < (Nacc+frame_no[start_row]))
//        {
//            
//            while (jj == 0)
//            {
//                Start_Time = t[k] ;
//                jj ++ ;
//            }
//
//            x_tem = (int) xFrac[k] ;
//            y_tem = (int) yFrac[k] ;
//            x_tem = x_tem * multi_factor ;
//            y_tem = y_tem * multi_factor ;
//
//
//            if (y_tem < IMG_DIM_FI && x_tem < IMG_DIM_FI)
//            {
//
//                float mult = (float) mult_phn[k] ;
//                float badflg = (float) bad_Flag[k] ;
//
//                 Eff_NoPhtn = (float) ENP[k] ;
//                //one_dim_exp[y_tem*IMG_DIM_FI+x_tem]=(one_dim_exp[y_tem*IMG_DIM_FI+x_tem]+1.0f)*mult*badflg ;
//                
//                one_dim_img[y_tem*IMG_DIM_FI+x_tem]=(one_dim_img[y_tem*IMG_DIM_FI+x_tem]+Eff_NoPhtn)*mult*badflg; 
//
////                Exposure_Img[y_tem][x_tem] = (Exposure_Img[y_tem][x_tem] + 1.0f) * mult*badflg ;
////
////                Eff_NoPhtn = (float) ENP[k] ;
////                //LOG(INFO)<<Eff_NoPhtn<<endl;exit(1);
////                //                if(Eff_NoPhtn!=28.7053){
////                //                    LOG(INFO)<<k<<" "<<Eff_NoPhtn<<endl;exit(1);
////                //                }
////                //                LOG(INFO)<<Eff_NoPhtn<<" "<<badflg<<" "<<" "<<k<<" "<<mult<<endl;
////                Image_Array[y_tem][x_tem] = (Image_Array[y_tem][x_tem] + Eff_NoPhtn) * mult*badflg ;
//
//            }
//
//            if (frame_no[k - 1] == frame_no[k])
//            {
//                time_frame = t[k] ;
//            }
//            else
//            {
//                time_frame_final = time_frame_final + time_frame ;
//            }
//                        if (k == nrows_cur_frm - 1)
//                        {
//                            goto label_else ;
//                        }
//        } 
//        else
//        {
//            LOG(INFO)<<k;
//label_else:
//
//            jj = 0 ;
//            End_Time = t[k - 1] ;
//
//            Avg_Time = Start_Time + (End_Time - Start_Time) / 2 ;
//
//            long temp = 0 ;
////            for (int i = 0 ; i < IMG_DIM_FI ; i ++)
////            {
////                for (int j = 0 ; j < IMG_DIM_FI ; j ++)
////                {
////                    one_dim_img[temp] = 0.0f ;
////                    one_dim_exp[temp] = 0.0f ;
////                    one_dim_img[temp] = (float) Image_Array[i][j] ;
////                    one_dim_exp[temp] = Exposure_Img[i][j] ;
////                    //                    if(one_dim_img[temp]!=0)
////                    //                    LOG(INFO)<<one_dim_img[temp]<<" "<<one_dim_exp[temp]<<endl;
////                    temp ++ ;
////                }
////            }
//            //            exit(1);
//            //Nacc = Nacc + frame_init ;
//            tot_num_frames ++ ;
//            LOG(INFO)<< "\rFrame created " << tot_num_frames  ;
//            cnt_frame = 0 ;
//            time_frame = 0.0 ;
//            time_frame_final = 0.0 ;
////            if (k != nrows_cur_frm - 2)
////            {
////                k -- ;
////            }
//           // obj.img_pixels_sig= new float[IMG_DIM_FI*IMG_DIM_FI];
//           // obj.Frame_pixels_exp=new float[IMG_DIM_FI*IMG_DIM_FI];
//            
////            for (int i = 0 ; i < IMG_DIM_FI * IMG_DIM_FI ; i ++)
////            {
////                obj.img_pixels_sig[i] = 0.0f ;
////                obj.Frame_pixels_exp[i] = 0.0f ;
////                obj.img_pixels_sig[i] = one_dim_img[i] ;
////                obj.Frame_pixels_exp[i]= badPix_Arr_fi[i]*frame_check;
////                //obj.Frame_pixels_exp[i] = one_dim_exp[i] ;
////                one_dim_img[i]=0;
////              //  one_dim_exp[i]=badPix_Arr_fi[i];
////                
////                
////            }
//             for (int i = 0 ; i < IMG_DIM_FI * IMG_DIM_FI ; i ++)
//            {
//               
//                //obj.img_pixels_sig.push_back (one_dim_img[i]) ;
//              //  obj.Frame_pixels_exp.push_back (badPix_Arr_fi[i]*frame_check);
//                //obj.Frame_pixels_exp[i] = one_dim_exp[i] ;
//                //one_dim_img[i]=0;
//              //  one_dim_exp[i]=badPix_Arr_fi[i];
//                 outputSigArr[i]=one_dim_img[i];
//                 outputExpArray[i]=one_dim_exp[i];
//                
//            }
//            outputFrmtime=Avg_Time;
//            outputFrmNo=tot_num_frames;
////            obj.frameTime_fi = Avg_Time ;
////            obj.frameNo_fi = tot_num_frames ;
//            
//            //vect.push_back (obj) ;
//          //  obj.img_pixels_sig.clear ();
//         //   obj.Frame_pixels_exp.clear ();
//        
//            //  delete[] obj.img_pixels_sig,obj.Frame_pixels_exp;
//
//        }
//
//    }
//     delete[] badPix_Arr_fi,badPix_Arr_padded;
//    return (EXIT_SUCCESS) ;
//    
//}

int uvtRelativeAspectPC::performFrameIntegration (long nrows , long start_row,unsigned short*frame_no , int xsize , int Ndiscard , int Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr_strct> &vect,
        float *outputSigArr,float *outputExpArray,double  &outputFrmtime,long &outputFrmNo,long &last_index_for_frame)
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
    //  int frame_init = frame_check = Nacc ;
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
    //noice_map_Sig_Array[i]=0.0f;
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
   // unsigned short frmno_previous=frame_no[NframesCount];
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
            //%#Edited  ON 20July#%
                 if(y_tem>0 && y_tem <IMG_DIM_FI && x_tem>0 && x_tem <IMG_DIM_FI)
//%#-Till this-20July17#%
                one_dim_img[y_tem*IMG_DIM_FI+x_tem]=one_dim_img[y_tem*IMG_DIM_FI+x_tem]+Eff_NoPhtn*mult*badflg; 
               // noice_map_Sig_Array[y_tem*IMG_DIM_FI+x_tem]=noice_map_Sig_Array[y_tem*IMG_DIM_FI+x_tem]+1.0*mult*badflg;
         

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




//Commneted as on 7-feb-17

//int uvtRelativeAspectPC::performFrameIntegration (long nrows , long start_row,unsigned short*frame_no , int xsize , int Ndiscard , int Nacc , float *xFrac , float *yFrac , unsigned short *mult_phn , unsigned short *bad_Flag , double *t , float *ENP , float *one_dim_img , float *one_dim_exp , vector<FrameIntegration_Arr_strct> &vect,
//        float *outputSigArr,float *outputExpArray,double &outputFrmtime,long &outputFrmNo,long &last_index_for_frame)
//{
// 
//    for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++) 
//    {
//        outputSigArr[i]=0.0f;
//        outputExpArray[i]=0.0f;
//    }
//  
//    FrameIntegration_Arr_strct obj ;
//    int jj = 0 ;
//    int x_tem , y_tem ;
//    double Start_Time , End_Time , Avg_Time = 0.0 ;
//
//  
//    float Eff_NoPhtn = 0.0f ;
//    float multi_factor = IMG_DIM_FI / xsize ;
//  
//    int frame_check ;
//    LOG (INFO) << "MULTIPLICATION FACTOR  for converting image to"<<IMG_DIM_FI<<" is " << multi_factor << endl ;
//    int frame_init = frame_check = Nacc ;
//    double time_frame = 0.0 , time_frame_final = 0.0 ;
//    int cnt_frame = 0 ;
//  //  vector<string> vhistorystr ;
//    
//    /****provision for if frames are not starting from 1****/
//   // int start_row = 0 ;
//    int cmpr_term = 0 ;
// 
//     if(Nacc>frame_no[nrows-1]){
//               LOG(ERROR)<<"***Number of frames to accumulate is greater than total available frames";
//                return(EXIT_FAILURE);
//    }
//    if (frame_no[0] != 1)
//    {
//        cmpr_term = Ndiscard + (frame_no[0] - 1) ;
//    }
////    for (int j = 0 ; j < nrows ; j ++)
////    {
////        if (frame_no[j] > cmpr_term)
////        {
////            start_row = j ;
////            break ;
////        }
////    }
//    float *badPix_Arr_padded= new float [600*600];
//    for(int i=0;i<600*600;i++) badPix_Arr_padded[i]=0.0;
//    float *badPix_Arr_fi= new float[IMG_DIM_FI*IMG_DIM_FI];
//      for(int i=0;i<600*600;i++) badPix_Arr_fi[i]=0.0;
//   int   status = Applypadding (badpixdata,512,512,badPix_Arr_padded,600,600);
//    if (status)
//    {
//        LOG (ERROR) << "Error in padding for bad pix data" ;
//        return (EXIT_FAILURE) ;
//    }
//    status= performSubDivisionIM (badPix_Arr_padded,600,600,badPix_Arr_fi,IMG_DIM_FI,IMG_DIM_FI);
//   if (status)
//    {
//        LOG (ERROR) << "Error in performing sub division in badpixel data" ;
//        return (EXIT_FAILURE) ;
//    } 
//    long int tot_num_frames = 0 ;
//    for (int i=0;i<IMG_DIM_FI*IMG_DIM_FI;i++)
//    {
//    one_dim_img[i]=0;
//    //cout<<badpixdata[i]<<endl;
//    one_dim_exp[i]=badPix_Arr_fi[i];
//    }
//    
//    LOG (INFO) << "Performing Frame Integration...." << endl ;
//    long p = 0 ;
//    
//    //calculate final frame count value for the current output frame.
//    int nrows_cur_frm=start_row;
//    for (int i=start_row;i<nrows;i++)
//    {
//  if (frame_no[i]<frame_no[start_row]+Ndiscard+Nacc)
//    {  
//      nrows_cur_frm++;
//    }
//    }
//    last_index_for_frame=nrows_cur_frm;
//  
//   //  LOG(INFO)<<"PPP "<<start_row<<" "<<nrows_cur_frm;
//    for (int k = start_row ; k < nrows_cur_frm ; k ++)
//    {
//        
//        /**in 'IF condition' accumulated array of  Nacc frames will be  generate   & In 'ELSE condition' array will be written to the frame.
//         x and y location  for number of events are taken as a indexes for the array .For each Nacc events one array will be written as a frame. 
//         **/
//      
//        if (frame_no[k] < (Nacc+frame_no[start_row]))
//        {
//          
//            while (jj == 0)
//            {
//                Start_Time = t[k] ;
//                jj ++ ;
//            }
//
//            x_tem = (int) xFrac[k] ;
//            y_tem = (int) yFrac[k] ;
//            x_tem = x_tem * multi_factor ;
//            y_tem = y_tem * multi_factor ;
//
//
//            if (y_tem < IMG_DIM_FI && x_tem < IMG_DIM_FI)
//            {
//
//                float mult = (float) mult_phn[k] ;
//                float badflg = (float) bad_Flag[k] ;
//
//                 Eff_NoPhtn = (float) ENP[k] ;
//                //one_dim_exp[y_tem*IMG_DIM_FI+x_tem]=(one_dim_exp[y_tem*IMG_DIM_FI+x_tem]+1.0f)*mult*badflg ;
//                
////                one_dim_img[y_tem*IMG_DIM_FI+x_tem]=(one_dim_img[y_tem*IMG_DIM_FI+x_tem]+Eff_NoPhtn)*mult*badflg; 
//                 one_dim_img[y_tem*IMG_DIM_FI+x_tem]=(one_dim_img[y_tem*IMG_DIM_FI+x_tem]+1)*mult*badflg;
//                 
////                Exposure_Img[y_tem][x_tem] = (Exposure_Img[y_tem][x_tem] + 1.0f) * mult*badflg ;
////
////                Eff_NoPhtn = (float) ENP[k] ;
////                //LOG(INFO)<<Eff_NoPhtn<<endl;exit(1);
////                //                if(Eff_NoPhtn!=28.7053){
////                //                    LOG(INFO)<<k<<" "<<Eff_NoPhtn<<endl;exit(1);
////                //                }
////                //                LOG(INFO)<<Eff_NoPhtn<<" "<<badflg<<" "<<" "<<k<<" "<<mult<<endl;
////                Image_Array[y_tem][x_tem] = (Image_Array[y_tem][x_tem] + Eff_NoPhtn) * mult*badflg ;
//
//            }
//
//            if (frame_no[k - 1] == frame_no[k])
//            {
//                time_frame = t[k] ;
//            }
//            else
//            {
//                time_frame_final = time_frame_final + time_frame ;
//            }
//                        if (k == nrows - 1)
//                        {
//                            lastFrame_flag=TRUE;
//                            LOG(INFO)<<"LAST FRAME CAME";
//                            break;
//                            //goto label_else ;
//                        }
//             if (k == nrows_cur_frm - 1)
//                        {
//                           goto label_else ;
//                        }
//        } 
//        else
//        {
//          
//label_else:
//          
//            jj = 0 ;
//            End_Time = t[k - 1] ;
//
//            Avg_Time = Start_Time + (End_Time - Start_Time) / 2 ;
//
//            long temp = 0 ;
//           // LOG(INFO)<<k;
//            tot_num_frames ++ ;
//            //cout<< "\rFrame created " << tot_num_frames  ;
//            
//            cnt_frame = 0 ;
//            time_frame = 0.0 ;
//            time_frame_final = 0.0 ;
//            
//             for (int i = 0 ; i < IMG_DIM_FI * IMG_DIM_FI ; i ++)
//            {
//               
//                //obj.img_pixels_sig.push_back (one_dim_img[i]) ;
//              //  obj.Frame_pixels_exp.push_back (badPix_Arr_fi[i]*frame_check);
//                //obj.Frame_pixels_exp[i] = one_dim_exp[i] ;
//                //one_dim_img[i]=0;
//              //  one_dim_exp[i]=badPix_Arr_fi[i];
//                 outputSigArr[i]=one_dim_img[i];
//                 outputExpArray[i]=one_dim_exp[i];
//                
//            }
//             //LOG(INFO)<<"PPLOPLPL";
//            outputFrmtime=Avg_Time;
//            outputFrmNo=tot_num_frames;
//         
////            obj.frameTime_fi = Avg_Time ;
////            obj.frameNo_fi = tot_num_frames ;
//            
//            //vect.push_back (obj) ;
//          //  obj.img_pixels_sig.clear ();
//         //   obj.Frame_pixels_exp.clear ();
//        
//            //  delete[] obj.img_pixels_sig,obj.Frame_pixels_exp;
//
//        }
//
//    }
//     delete[] badPix_Arr_fi,badPix_Arr_padded;
//    return (EXIT_SUCCESS) ;
//}


int uvtRelativeAspectPC::performUnitConversion (float  *frmdata , double  intgrntime , int size)
{
    LOG(INFO)<<"The integration time "<<intgrntime;
    if (intgrntime <= 0)
    {
        LOG (ERROR) << "***Divide by Zero***" << endl ;
        return (EXIT_FAILURE) ;

    }
    //cout<<intgrntime<<endl;exit(1);
    for (int i = 0 ; i < size ; i ++)
    {
        frmdata[i] = frmdata[i] /  intgrntime ;

    }
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectPC::performFlatFieldCorr (float* frmsigdata , float* frmflatfield , float *Xfract , float *Yfract , int size)
{
    int div_fact = subDivision_size / 600 ;
    for (int i = 0 ; i < size ; i ++)
    {
         frmsigdata[i] = frmsigdata[i] * frmflatfield[(int) (round ((Yfract[i] / div_fact))*600 + (int) round ((Xfract[i]) / div_fact))] ;
        //frmsigdata[i] = frmsigdata[i] * frmflatfield[(int) (round ((Yfract[i]))*600 + (int) round ((Xfract[i])))] ;
    }
    return (EXIT_SUCCESS) ;

}

int uvtRelativeAspectPC::getHistory (vector<string> &vhistory)
{
   // char *user = getlogin () ;
    int cnt = 0 ;
    char validgtiflag_str[FLEN_FILENAME] ;
    //string str = "Module run by " + (string) user ;
    char temp[PIL_LINESIZE];
    //vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Input Level1 tar  file = " + (string) level1indir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Caldb used= " + (string) caldbindir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Module Output directory = " + (string)level2outdir ) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " channel used= " + (string)channel) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Parity check flag for photon event= " + (string)convertIntToStr(parity_flag)) ;
vhistory.push_back ((string) getSerialNo (cnt) + " CRC flag = " + (string) convertIntToStr (crc_flag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Flag for CRC failure action = " + (string) convertIntToStr (dropframe)) ;
    //vhistory.push_back ((string) getSerialNo (cnt) + " Dark frame subtraction to be done or not = " + (string)convertIntToStr (darkframe_flag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " GTI flag  = " + (string)convertIntToStr(gti_flag)) ;
    
    vhistory.push_back ((string) getSerialNo (cnt) + " Threshold for multi-photon events " + (string)convertIntToStr(thr_multiph)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " CR threshold parameter-1-“_N_” in threshold=AVG+N*sqrt(avg_events)+ST/(sqrt(avg_events)) to identify Cosmic Ray affected frame = " + (string)convertIntToStr( first_mult_factor_CR)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " CR threshold parameter-2-“_ST_” in threshold=AVG+N*sqrt(avg_events)+ST/(sqrt(avg_events)) to identify Cosmic Ray affected frame = " + (string)convertFloatToStr(second_mult_Factor_CR)) ;
   
  
 
            vhistory.push_back ((string) getSerialNo (cnt) + " Number of Initial frames to be discarded for MULTI case = " + (string) convertIntToStr(nFrameDiscard_fi)) ;
            vhistory.push_back ((string) getSerialNo (cnt) + " Number of frames to be combined in MULTI case = " + (string) convertIntToStr(nFrameIntegrate_fi)) ;
  
    vhistory.push_back ((string) getSerialNo (cnt) + " Frame size after frame integration = " + (string) convertIntToStr(IMG_DIM_FI)) ;
            vhistory.push_back ((string) getSerialNo (cnt) + " Star finding algorithm = " + (string) convertIntToStr(star_detect_algo_flag)) ;
             
            if(star_detect_algo_flag==1 || star_detect_algo_flag==3 || star_detect_algo_flag==4)
            {
                vhistory.push_back ((string) getSerialNo (cnt) + " First-cut threshold for star detection = " + (string) convertFloatToStr (sd_multi_factor_default)) ;
                vhistory.push_back ((string) getSerialNo (cnt) + " Minimum targeted stars = " + (string)convertIntToStr(minimum_No_of_Stars)) ;
            }
            
             vhistory.push_back ((string) getSerialNo (cnt) + " Neighbourhood criteron for identifying stars = " + (string)convertIntToStr( refine_Winsize)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Box Size to compute Centroid for  detected stars = " + (string)convertIntToStr( centroid_Winsize)) ;
             
            vhistory.push_back ((string) getSerialNo (cnt) + " No. of accumulated frames discarded at beginning = " + (string)convertIntToStr(frames_toDiscard)) ;
            vhistory.push_back ((string) getSerialNo (cnt) + " No. of accumulated frames used to construct reference frame = " + (string)convertIntToStr( nFrameToAverage)) ;
           
           
//             if(search_algo_ctlg==1 || search_algo_ctlg==3 || search_algo_ctlg==5){
//            vhistory.push_back ((string) getSerialNo (cnt) + " Lengh of rectangle search   = " + (string)(len_a)) ;    
//             vhistory.push_back ((string) getSerialNo (cnt) + " Width of rectangle search  = " + (string) (len_b)) ;  
//             }
            
             
             //vhistory.push_back ((string) getSerialNo (cnt) + " frame to be discarded in reference frame calculation = " + (string)convertIntToStr(frames_toDiscard)) ;
            //vhistory.push_back ((string) getSerialNo (cnt) + " Average Factor = " + (string)convertIntToStr( nFrameToAverage)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Minimum search distance for star match in succesive images of FrameIntegration stage. " + (string)convertFloatToStr (pair_nbhd_distance)) ;
             //vhistory.push_back ((string) getSerialNo (cnt) + " Frequency domain filtering flag  " + (string) convertIntToStr(freqDomainFilter_Flag)) ;
             //if(freqDomainFilter_Flag==1)
             //vhistory.push_back ((string) getSerialNo (cnt) + " Type Filtering used   " + (string) convertIntToStr(type_Filtering)) ;
              //vhistory.push_back ((string) getSerialNo (cnt) + " Catalogue database path = " + (string)databasename) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Frequency domain filtering flag  " + (string) convertIntToStr(FreqDomainFilterFlag)) ;
             if(FreqDomainFilterFlag==1)
             vhistory.push_back ((string) getSerialNo (cnt) + " Type Filtering used   " + (string) convertIntToStr(type_Filtering)) ;
             
             if(type_Filtering==1 || type_Filtering==2 ){
                  vhistory.push_back ((string) getSerialNo (cnt) + " order of Pitch   " + (string) convertIntToStr(orderpitch)) ;
                  vhistory.push_back ((string) getSerialNo (cnt) + " order of roll   " + (string) convertIntToStr(orderroll)) ;
                  vhistory.push_back ((string) getSerialNo (cnt) + " order of yaw  " + (string) convertIntToStr(orderyaw)) ;
                   vhistory.push_back ((string) getSerialNo (cnt) + " Delta time    " + (string) convertFloatToStr(poly_fit_interval)) ;
             }
vhistory.push_back ((string) getSerialNo (cnt) + " Frequency value=    " + (string) convertFloatToStr(freqvalue)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Flag for writing matched star pair list" + (string) convertIntToStr(match_stars_file_flag)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Theta switch(ON/OFF)=   " + (string) convertIntToStr(flag_thetaComp)) ;
	     vhistory.push_back ((string) getSerialNo (cnt) + " Shift and rotation algorithm =  " + (string) convertIntToStr(shift_rotation_algo)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for unitConversion= " + (string)convertIntToStr( wtd_uc)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for flat field correction= " + (string)convertIntToStr( wtd_ff)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for mask Bad pixel correction= " + (string)convertIntToStr(wtd_bp)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for Pix padding= " + (string)convertIntToStr( wtd_pp)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for  subDivision= " + (string)convertIntToStr( wtd_sd)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for cosmic ray correction= " + (string)convertIntToStr( wtd_cr)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for QEMCP correction= " + (string)convertIntToStr( wtd_qemcp)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for detector distortion=  " + (string)convertIntToStr( wtd_dd)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for optical distortion= " + (string)convertIntToStr( wtd_od)) ;
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
