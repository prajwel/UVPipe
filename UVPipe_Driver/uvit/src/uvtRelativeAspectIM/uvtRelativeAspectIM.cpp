/* 
 * File:   uvtImRa_commonArray.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <stdlib.h>

#include "uvtRelativeAspectIM.h"


#include<uvtComputeDrift.h>
#include "uvtDetectStar.h"
#include<Directory.h>
#include<DataInfo.h>
#include<DataIngest.h>
#include <algorithm>
#include<uvtQEMCPCorr.h>
#include<uvtComputeJitter.h>
#include<uvtComputeThermal.h>
#include<uvtRelAspCal.h>

#include <vector>

#include<spMatrix.h>
#include<glog/logging.h>
#include<uvtUtils.h>
#include<macro_def.h>

#define  NBHD_RADIUS 5
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr   //used in doFft() function
#define moduleoutdir_unit "uvtUnitConvertion"
#define moduleoutdir_badpix "uvtMaskBadPix"
#define moduleoutdir_fltfield "uvtFlatFieldCorrection"
#define moduleoutdir_qe "uvtQEMCPCorr"
#define moduleoutdir_pixpad "uvtPixPadding"
#define moduleoutdir_subdiv "uvtSubDivision"
#define moduleoutdir_cosmicray "uvtCosmicRayCorrection"
#define moduleoutdir_AccEverytsec "uvtAccEveryTsec"
#define moduleoutdir_findstarcentroid "uvtFindStarCentroid"
#define moduleoutdir_detectordistortion "uvtDetectDistortion"
#define moduleoutdir_opticaldistortion "uvtOpticalDistortion"
#define moduleoutdir_driftExercise "uvtComputeDrift"
#define moduleoutdir_refFrameCal "uvtRefFrameCal"



uvtRelativeAspectIM::uvtRelativeAspectIM () {
    tar_extracted_flag_IM=FALSE;
    
}


uvtRelativeAspectIM::~ uvtRelativeAspectIM () {
    
}


int uvtRelativeAspectIM::readPILParameters (int argc , char** argv)
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
   // if(this->paramfile_Varask==FALSE){
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
   // }
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

    if (PIL_OK != (status = PILGetBool ("darkframeFlag" , &darkframe_flag)))
    {
        LOG (INFO) << endl << "***Error Reading darkframeFlag:" << unitConversionFlag << "***" ;
        return status ;
    }
   if (PIL_OK != (status = PILGetBool ("junkfileflag" , &junkFrameFlag)))
    {
        LOG (INFO) << endl << "***Error Reading GTI flag value***" ;
        return status ;
    }
    if(junkFrameFlag)
    {
       if (PIL_OK != (status = PILGetReal4 ("thresholdforjunkFrame" , &thrJunkFrame)))
        {
            LOG (INFO) << endl << "***Error reading standard deviation multiplication factor value ***" ;
            return status ;
        } 
    }

    if (PIL_OK != (status = PILGetBool ("GTI_FLAG" , &gti_flag)))
    {
        LOG (INFO) << endl << "***Error Reading GTI flag value***" ;
        return status ;
    }
    if (gti_flag)
    {
        all_or_custom = 1 ;
        valid_bit = 1 ;
    }
    if (PIL_OK != (status = PILGetBool ("flatfieldFlag" , &flatfieldFlag)))
    {
        LOG (INFO) << endl << "***Error Reading flatfieldFlag:" << flatfieldFlag << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("paddingDim" , &padding_dim)))
    {
        LOG (INFO) << endl << "***Error Reading padding_dim :" << padding_dim << "***" ;
        return status ;
    }
   // if (PIL_OK != (status = PILGetInt ("att_flagVal" , &att_flag_val)))
    //{
      //  LOG (INFO) << endl << "***Error Reading att  flag value :" << padding_dim << "***" ;
        //return status ;
    //}


    if (PIL_OK != (status = PILGetInt ("Nacc" , &no_ofFramesToAcc)))
    {
        LOG (INFO) << endl << "***Error Reading Nacc :" << no_ofFramesToAcc << "***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetBool ("qemcpFlag" , &qe_mcpFlag)))
    {
        LOG (INFO) << endl << "***Error Reading qemcpFlag :" << qe_mcpFlag << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("subdivisionFlag" , &subdivisionFlag)))
    {
        LOG (INFO) << endl << "***Error Reading subdivisionFlag :" << subdivisionFlag << "***" ;
        return status ;
    }
    if (subdivisionFlag == 1)
    {

        if (PIL_OK != (status = PILGetInt ("subdivisionsize" , &subDivision_size)))
        {
            LOG (INFO) << endl << "***Error Reading subDivision size :" << subDivision_size << "***" ;
            return status ;
        }
    }
    if (PIL_OK != (status = PILGetInt ("algoFlag" , &star_detect_algo_flag)))
    {
        LOG (INFO) << endl << "***Error reading algo_flag ***" ;
        return status ;
    }
    if (star_detect_algo_flag == 1 || star_detect_algo_flag == 3 || star_detect_algo_flag == 4)
    {
        if (PIL_OK != (status = PILGetReal4 ("threshold" , &rms_mul_factor_default)))
        {
            LOG (INFO) << endl << "***Error reading standard deviation multiplication factor value ***" ;
            return status ;
        }

        if (PIL_OK != (status = PILGetReal4 ("minimumTargetedstars" , &min_num_stars)))
        {
            LOG (INFO) << endl << "***Error reading minimum number of stars in star detection ***" ;
            return status ;
        }
    }
    else if (star_detect_algo_flag == 2 )
    {

//        if (PIL_OK != (status = PILGetInt ("StarDetectionCentroidSqrSize" , &algo_Square_Size)))
//        {
//            LOG (INFO) << endl << "***Error reading refine window size ***" ;
//            return status ;
//        }
//        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroidPriThr" , &primary_threshold_Val)))
//        {
//            LOG (INFO) << endl << "***Error reading refine window size ***" ;
//            return status ;
//        }
//        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroidSecThr" , &secondary_threshold_Val)))
//        {
//            LOG (INFO) << endl << "***Error reading refine window size ***" ;
//            return status ;
//        }
    }
    else
    {
        LOG (INFO) << "***algo flag must be 1 or 2***" << endl ;
        return (EXIT_FAILURE) ;
    }

    if (PIL_OK != (status = PILGetInt ("refineWindow" , &refine_Winsize)))
    {
        LOG (INFO) << endl << "***Error reading refine window size ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("centroidWindow" , &centroid_Winsize)))
    {
        LOG (INFO) << endl << "***Error reading centroid  window size ***" ;
        return status ;
    }
    if (star_detect_algo_flag == 1)
    {
        if (PIL_OK != (status = PILGetReal4 ("background_fact" , &backgrd_fact)))
        {
            LOG (INFO) << endl << "***Error reading background_fact ***" ;
            return status ;
        }
    }
    else
    {
        backgrd_fact = 0 ; //can be ignored.
    }
    if (PIL_OK != (status = PILGetInt ("Window_nb" , &search_win_size)))
    {
        LOG (INFO) << endl << "***Error reading centroid  window size ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("GenMatchStarsFile_flag" , &match_stars_file_flag)))
    {
        LOG (INFO) << endl << "***Error reading GenMatchStarsFile_flag parameter***" ;
        return status ;
    }



    if (PIL_OK != (status = PILGetReal ("ThresholdValue" , &threshold)))
    {
        LOG (INFO) << endl << "***Error reading the threshold value for Cosmic Ray correction module***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("framesToBeDiscard" , &frames_toDiscard)))
    {
        LOG (INFO) << endl << "***Error reading the number of frames to discard ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("averageFactor" , &nFrameToAverage)))
    {
        LOG (INFO) << endl << "***Error reading average factor value for reference frame calculation ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetReal4 ("diffDist" , (float*) &nbhd_dist)))
    {
        LOG (INFO) << endl << "***Error reading diffDist parameter***" ;
        return status ;
    }
//    if (PIL_OK != (status = PILGetReal4 ("error_per" , &err_per)))
//    {
//        LOG (INFO) << endl << "***Error reading error percent  ***" ;
//        return status ;
//    }
    if (PIL_OK != (status = PILGetInt ("freqDomainFilterFlag" , &freqDomainFilter_Flag)))
    {
        LOG (INFO) << endl << "***Error reading freqDomainFilter Value***" ;
        return status ;
    }
    if (freqDomainFilter_Flag == 0)
    {
        if (PIL_OK != (status = PILGetInt ("typeFiltering" , &type_Filtering)))
        {
            LOG (INFO) << endl << "***Error reading type Filtering Value***" ;
            return status ;
        }
        if (type_Filtering == 1 || type_Filtering==2)
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
                if (PIL_OK != (status = PILGetReal ("deltaTime" , &delta_time)))
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
    if (freqDomainFilter_Flag == 1)
    {
        if (PIL_OK != (status = PILGetReal ("freqValue" , (double*) &freqvalue)))
        {
            LOG (INFO) << endl << "***Error reading cut-off frequency value***" ;
            return status ;
        }
    }

    if (PIL_OK != (status = PILGetInt ("shiftRotDetAlgoFlag" , &shift_rotation_algo)))
    {
        LOG (INFO) << endl << "***Error reading shiftRotDetAlgoFlag***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetBool ("clobber" , &clobber)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobber << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("flag_thetaComp" , &flag_thetaComp)))
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
    if (flatfieldFlag==1){
    if (PIL_OK != (status = PILGetInt ("Write_todiskff" , &wtd_ff)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the FlatFieldCorrection:" << "***" ;
        return status ;
    }
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskpp" , &wtd_pp)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (subdivisionFlag == 1)
    {
        if (PIL_OK != (status = PILGetInt ("Write_todisksd" , &wtd_sd)))
        {
            LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the subDivision:" << "***" ;
            return status ;
        }
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskcr" , &wtd_cr)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
    if(qe_mcpFlag==1)
    {
    if (PIL_OK != (status = PILGetInt ("Write_todiskqe" , &wtd_qemcp)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the qemcp correction:" << "***" ;
        return status ;
    }
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskac" , &wtd_ac)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the Accumulated module:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todisksc" , &wtd_fsc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the find Star Centroids:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskdd" , &wtd_dd)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the Detect Distortion:" << "***" ;
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


int uvtRelativeAspectIM::uvtRelativeAspectIMProcess ()
{
    int status = 0 ; //status flag to check return status of functions
    //added
   // vector<string> l2output_tarnames;
  //  string l2path_temp="/home/dhruv/Desktop/";
   // ExpandMergeData(level1indir,l2path_temp,l2output_tarnames);
   // int total_tars=0;
   // while(total_tars<l2output_tarnames.size())
   // {
//        level2outdir=level2outdir+"_"+total_tars;
    //till this.ExpandMergeData(level1tar,l2path_temp);
    
    string temp_str ;
  //  temp_str.assign (l2output_tarnames[total_tars]) ;
    temp_str.assign (level1indir) ;
    if(tar_extracted_flag_IM==FALSE)//incase tar is already extracted than no need for extracting again (mainly for driver module)
    {
        
    level1indir = "" ;
    
    if(zipFlag==FALSE)
    {
    status = extractTars (temp_str , level1indir , orbnum) ;//extract level-1 data  tar file
    
    if (status)
    {
        LOG (INFO) << "Error in extracting tar" ;
LOG(ERROR)<<"CRASH: TAR EXTRACTION FAILED (uvtRelativeAspectIM.cpp)";
        return (EXIT_FAILURE) ;
    }
    }
    else{
       status = extractZip (temp_str , level1indir , orbnum) ;//extract level-1 data  tar file
    
    if (status)
    {
        LOG (INFO) << "Error in extracting tar" ;
LOG(ERROR)<<"CRASH: ZIP EXTRACTION FAILED (uvtRelativeAspectIM.cpp)";
        return (EXIT_FAILURE) ;
    } 
    }
    
    }
    
    //history reading
    cout << "level1 " << level1indir << endl ;
    //check Whether  level -1 directory Exist or not.
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



    //function to get level 1 files
    Directory dirobj ;
    dirobj.setOrbitNo (orbnum) ;

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
    {
        LOG (ERROR) << endl << "***Error in directory set up***" << endl ;
        return (EXIT_FAILURE) ;
    }
   
    int numofsciencedatafiles = dirobj.sciencedatafile.size () ; //total number of science data file
    
    if (numofsciencedatafiles <= 0)
    {
        LOG (ERROR) << endl << "No science data file found" << endl ;
        return (EXIT_FAILURE) ;
    }

    char obsmode[FLEN_VALUE] ;
    char moduleIndir[FLEN_FILENAME] ; //hold path for input directory for every module, will be updated after every module run
    char outputdir[FLEN_FILENAME] ;
    int division_fact ;

    uvtDetectStar obj_sc ; //creating the object for finding peaks for the image
    //incase of subdivision  not to be done
    if (subdivisionFlag == 0)       subDivision_size = padding_dim ;
    //Array for storing  pixel padding 
    frame_Data_Padded = new float[padding_dim * padding_dim] ;
    frame_ExpoData_padded = new float[padding_dim * padding_dim] ;
    //Array for storing subdivided pixels incase of subdivision to be done
    if (subdivisionFlag)
    {
        frame_Data_subdivided = new float[subDivision_size * subDivision_size] ;
        frame_ExpData_subdivided = new float[subDivision_size * subDivision_size] ;
        division_fact = subDivision_size / padding_dim ;
        division_fact = division_fact*division_fact ;
    }

    char errstr[500] ;
    unsigned short frameno ;
    double frametime , integrationtime ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    char caldb_common_dir_od[PIL_LINESIZE] , caldb_common_dir_qe[PIL_LINESIZE] ;
    char caldb_common_dir_dd[PIL_LINESIZE] , caldb_common_dir_ff[FLEN_FILENAME] ;
    
    
    //loop for number of science data file
         bool flag_rfc = FALSE ;
         int nFrame_Avg_ref=nFrameToAverage;
         fitsfile *fptr;
         long tot_Num_Rows=0;
         vector<string> L1keywords;
    for (int dataindex = 0 ; dataindex < numofsciencedatafiles ; dataindex ++)
    {
        status =0;
        
        fits_open_file(&fptr , dirobj.sciencedatafile[dataindex].c_str () , READONLY , &status) ;
        copyUsrkeywrdsTovect(fptr,L1keywords);
        fits_movabs_hdu (fptr , 3 , NULL , &status) ;
        fits_get_num_rows (fptr , &tot_Num_Rows, &status) ;
   
    unsigned short *frmCount_Vis = new unsigned short[tot_Num_Rows];
    double *frmTime_Vis = new double[tot_Num_Rows];
   
    fits_read_col (fptr , TUSHORT , 11 , 1 , 1 , tot_Num_Rows , NULL , frmCount_Vis , NULL , &status) ;
    LOG(INFO)<<status;
    fits_read_col (fptr , TDOUBLE , 1 , 1 , 1 , tot_Num_Rows , NULL , frmTime_Vis , NULL , &status) ;
    LOG(INFO)<<status;
    fits_close_file(fptr,&status);
    LOG(INFO)<<frmTime_Vis[tot_Num_Rows-1]-frmTime_Vis[0];
    if(frmTime_Vis[tot_Num_Rows-1]-frmTime_Vis[0]<=20){
        LOG(INFO)<<"NOT CONTAIN ENOUGH DATA ,SKIPPING "<<dirobj.sciencedatafile[dataindex]<<" file";
	LOG(ERROR)<<"CRASH DATA LESS THAN 20 SECONDS ! (uvtRelativeAspectIM.cpp) ";
        delete[] frmCount_Vis,frmTime_Vis;
        continue;
    }   
    delete[] frmCount_Vis,frmTime_Vis;
    LOG(INFO) <<(char *) dirobj.sciencedatafile[dataindex].c_str ();
        nFrameToAverage=nFrame_Avg_ref;
        flag_rfc = FALSE ;
        ref_frame_module_filename_track.clear () ;//clears the vector that stores file name of generated reference frame calculation module
        LOG (ERROR) << endl << "------------------Data Set " << dataindex + 1 << " : " << dirobj.sciencedatafile[dataindex] << "-----------------------" << endl ;
        
        /*---finding mode from science data file---*/
        getKeywordVal ((char *) dirobj.sciencedatafile[dataindex].c_str () , "OBS_MODE" , 1 , obsmode) ;
        //check for mode parameter.
        if (strcasecmp (obsmode , "IM") != 0)
        {
            LOG (ERROR) << endl << "Observation mode is " << obsmode << "  in file  " << dirobj.sciencedatafile[dataindex] ;
	    LOG(ERROR)<<" CRASH DATA NOT IM MODE  ! (uvtRelativeAspectIM.cpp) ";
            LOG (ERROR) << endl << "Checking next file...." << endl ;
            continue ; //go to next file if obs mode is not IM
        }

        LOG (INFO) << endl << "Data Mode is " << obsmode << endl ;
        LOG (INFO) << "Information about level2 path" << dirobj.level2path[dataindex].c_str () ;
        strcpy (outputdir , dirobj.level2path[dataindex].c_str ()) ;
        
        if (wtd_bp == 1)                                                                        sprintf (moduleoutdir_bp , "%s/%s_%s" , outputdir , moduleoutdir_badpix , VERSION) ;
        if (wtd_uc == 1)                                                                         sprintf (moduleoutdir_uc , "%s/%s_%s" , outputdir , moduleoutdir_unit , VERSION) ;
        if (flatfieldFlag == 1 && wtd_ff == 1)                                        sprintf (moduleoutdir_ff , "%s/%s_%s" , outputdir , moduleoutdir_fltfield , VERSION ) ;
        if (qe_mcpFlag == 1 && wtd_qemcp == 1)                               sprintf (moduleoutdir_qemcp , "%s/%s_%s" , outputdir , moduleoutdir_qe , VERSION ) ;
        if (wtd_pp == 1)                                                                        sprintf (moduleoutdir_pp , "%s/%s_%s" , outputdir , moduleoutdir_pixpad , VERSION) ;
        if (subdivisionFlag == 1 && wtd_sd)                                         sprintf (moduleoutdir_sd , "%s/%s_%s" , outputdir , moduleoutdir_subdiv , VERSION) ;
        if (wtd_cr == 1)                                                                         sprintf (moduleoutdir_cr , "%s/%s_%s" , outputdir , moduleoutdir_cosmicray , VERSION) ;
        if (wtd_ac == 1)                                                                         sprintf (moduleoutdir_ac , "%s/%s_%s" , outputdir , moduleoutdir_AccEverytsec , VERSION) ;
        if (wtd_fsc == 1)                                                                        sprintf (moduleoutdir_sc , "%s/%s_%s" , outputdir , moduleoutdir_findstarcentroid , VERSION) ;
        if (wtd_dd == 1)                                                                         sprintf (moduleoutdir_dd , "%s/%s_%s" , outputdir , moduleoutdir_detectordistortion , VERSION) ;
        if (wtd_od == 1)                                                                         sprintf (moduleoutdir_od , "%s/%s_%s" , outputdir , moduleoutdir_opticaldistortion , VERSION) ;
        sprintf (moduleoutdir_de , "%s/%s_%s" , outputdir , moduleoutdir_driftExercise , VERSION) ;
        sprintf (moduleoutdir_rfc , "%s/%s_%s" , outputdir , moduleoutdir_refFrameCal , VERSION) ;

        char *a0,*a1,*a2;
      
        //status=ReadAtd_status ((char*)dirobj.attfile.c_str (),6,a0,a1,a2,listofData);
        
        //----------DATAINGEST----------//
        LOG (ERROR) << "===================================DATAINGEST===================================================" << endl ;
        
        strcpy (lbtfile , (char *) dirobj.lbtfile.c_str ()) ;
     /*process starts for DataIngest */
        DataIngest di_obj ;
        di_obj.read ((char *) dirobj.sciencedatafile[dataindex].c_str () , (char*) caldbindir.c_str () , (char *)dirobj.tctfile.c_str () ,(char *) dirobj.mkffile.c_str () ,(char *) dirobj.gtifile[dataindex].c_str () ,
                (char *) dirobj.lbtfile.c_str () ,(char*)dirobj.attfile.c_str(), (char*) dirobj.darkDirectory.c_str () , att_flag_val,gti_flag , valid_bit , all_or_custom , outputdir , dropframe , parity_flag ,UTC_flag,crc_flag, clobber , history) ;

        di_obj.display () ;
        status = di_obj.DataIngestProcess () ;
        if (status)
        {
            LOG (ERROR) <<  "***Error in Data Ingest Process***"  ;
	    LOG(ERROR)<<" CRASH DATAINGEST FAILED ! (uvtRelativeAspectIM.cpp) ";
            continue;
        }

        strcpy (moduleIndir , di_obj.getModuleOutdir ()) ;
        strcpy (Indir_dataIngest , di_obj.getModuleOutdir ()) ;
        LOG (INFO) << "Using directory " << moduleIndir << "  as input to uvtUnitConversion"  ;
      /*process ends for DataIngest */

        //setting output directory path incase of output to e written
        //setting path for input information file of Dataingest module
      
        char infofile_in[PIL_LINESIZE] ;
        string tempfilepath = searchFile (moduleIndir , ".info") ;
        if (tempfilepath == " ")
        {
            LOG (ERROR) << "***Information file not found in " << moduleIndir << "***" ;
            continue ;
        }
        sprintf (infofile_in , "%s/%s" , moduleIndir , tempfilepath.c_str()) ;
        LOG (INFO) << "Information file :" << infofile_in ;
        /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
        if (! (FileExists (infofile_in)))
        {
            LOG (ERROR) << endl << "***Input FileList not Found at Specified PATH,Check Input Directory***" ;
            continue ;
        }
        /*
     open the .info FITS file  and read the header information from the second HDU.
         */
        vector<string> header_info;
        getHistory (header_info);
        writeHistory (infofile_in,header_info);
        fitsfile *finfo_in , *finfo_out ;
        fits_open_file (&finfo_in , infofile_in , READWRITE , &status) ;
        printError (status , "Error in opening the input information file",infofile_in) ;
        //copyUsrkeywrdsTovect(finfo_in,L1keywords);
        fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
        printError (status , "Error in moving to 2nd HDU",infofile_in) ;
        datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
        xsize = datainfo.getXsize () ;
        ysize = datainfo.getYsize () ;
        int nframes ;
        char nameprefix[PIL_LINESIZE] ;
        fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
        printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file
        fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "Error in reading the  key value of the NFILES " , infofile_in) ;
        fits_read_key (finfo_in , TSTRING , "DARKDIR" , darkdir , NULL , &status) ;
        printError (status , "Error in reading the  key value of the NFILES " , infofile_in) ;
        fits_movabs_hdu (finfo_in , 1 , NULL , &status) ;
        printError (status , "Error in moving to 1st HDU ",infofile_in) ;
        fits_read_key (finfo_in , TINT , "WIN_X_SZ" , &win_xsize , NULL , &status) ;
        printError (status , "Error in reading the key value of the Window xsize" , infofile_in) ; //for creating name for output information file
        fits_read_key (finfo_in , TINT , "WIN_Y_SZ" , &win_ysize , NULL , &status) ;
        printError (status , "Error in reading the key value of the Window ysize " , infofile_in) ; //for creating name for output information file
        fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
        printError (status , "Error in moving to 2nd HDU of information file",infofile_in) ;
     
        
        if (win_xsize != win_ysize)
        {
            LOG (ERROR)    << "***Window is not square***" ;
LOG(ERROR)<<"CRASH WINDOW NOT SQUARE ! (uvtRelativeAspectIM.cpp) ";
            continue;
        }
        char **sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        //reading frame names from information file into vector
        fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) sigframelist , NULL , &status) ;
        printError (status , "Error in reading the signal frame list from the information file" , infofile_in) ;
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
//         status = writeUsrkeywordsFrmvect (infofile_in , L1keywords) ;
//    if (status)
//    {
//        LOG (INFO) << "Information file creation completed" ;
//        return (EXIT_FAILURE) ;
//    }
        char infile[NAMESIZE] ;
        /*===========================================================================================================================================*/
        /*========================================================CALDB reading started==================================================================== */
         /*===========================================================================================================================================*/
        //copying caldb directory path for different module for further use.
        strcpy (caldb_common_dir_dd , (char*) caldbindir.c_str ()) ;
        strcpy (caldb_common_dir_od , (char*) caldbindir.c_str ()) ;
        strcpy (caldb_common_dir_qe , (char*) caldbindir.c_str ()) ;
        strcpy (caldb_common_dir_ff , (char*) caldbindir.c_str ()) ;
        
        //caldb reading started for bad pixel correction 
        string tempname = caldb_handler.getBadPixelFile (datainfo.getDetector () , datainfo.getObsMode () , win_xsize + 1 , win_ysize + 1 , (char*) caldbindir.c_str ()) ;
        
                
        if (tempname ==" ")
        {
            LOG (ERROR) << "Couldn't find bad pixel file from calDB"  ;
            continue ;
        }
        joinStrings (badpixfile , 1 , tempname.c_str()) ; //setting path for badpixel file
        LOG (INFO) << "Caldb BAD pixel file is " << badpixfile << endl ;
        badpixdata = new float[xsize * ysize] ;       
        status = readImage (badpixfile , 1 , badpixdata , xsize , ysize) ; //reading badpixel file.
        if (status)
        {
            LOG (ERROR) << "Error in reading the BAD pixel  image from caldb  " ;
            return (EXIT_FAILURE) ;
        }
        //caldb reading started for flat field correction if flat field to be done 
        if (flatfieldFlag == 1)
        {
            string   tempname1 = caldb_handler.getFlatFieldFile (datainfo.getDetector () , datainfo.getObsMode () , datainfo.getFilter () , caldb_common_dir_ff) ;
             if (tempname1==" "){
                 LOG (ERROR) << endl << "Couldn't find flat field file from caldb" << endl ;
            continue ;
             }
            joinStrings (flatfieldfile , 2 , caldb_common_dir_ff , tempname1.c_str()) ;
            flatfielddata = new float[xsize * ysize] ;
            status = readImage (flatfieldfile , 1 , flatfielddata , xsize , ysize) ; //reading flatfield image from caldb 
            if (status)
            {
                LOG (ERROR) << "Error in reading the flat field correction image from caldb" ;
                return (EXIT_FAILURE) ;
            }
        }

        //caldb reading started for detector distortion
        tempname = caldb_handler.getDetectorFile (datainfo.getDetector () , (char*)caldb_common_dir_dd) ;
         if (tempname == " ")
        {
            LOG (ERROR) << endl << "Couldn't find detector file from caldb" << endl ;
            continue ;
        }
        joinStrings (detector_distortion_corr_file , 2 , (char*)caldb_common_dir_dd ,tempname.c_str()) ;
        LOG (INFO) << "Detector distortion correction file is " << detector_distortion_corr_file<<" ppp"<<tempname ;
        

        //Array for storing the detector distortion file
        X_detect_distortion = new float[xsize * ysize] ;
        Y_detect_distortion = new float[xsize * ysize] ;
        status = caldb_handler.readCaldbDistFile (X_detect_distortion , Y_detect_distortion , detector_distortion_corr_file) ;//reading file for distortion correction file  from CALDB  for detector Distortion
        if (status)
        {
            LOG (ERROR) << "Error in reading distortion correction file from caldb" ;
            return (EXIT_FAILURE) ;
        }
        //caldb reading started for optical distortion correction
        tempname = caldb_handler.getOpticalDetectorFile (datainfo.getDetector () ,(char*) caldb_common_dir_od , datainfo.getFilter ()) ;
        if (tempname==" ")
        {
            LOG (ERROR)  << "Couldn't find Optical Dist detector File From caldb" ;
            continue ;
        }
        joinStrings (optical_distortion_corr_file , 2 , (char*)caldb_common_dir_od , tempname.c_str()) ;
        LOG (INFO) << "Optical distortion correction file is " << optical_distortion_corr_file ;
        //Array for storing CALDB optical distortion correction data
        X_optical_distortion = new float[xsize * ysize] ;
        Y_optical_distortion = new float[xsize * ysize] ;
        status = caldb_handler.readCaldbDistFile (X_optical_distortion , Y_optical_distortion , optical_distortion_corr_file) ;
        if (status)
        {
            LOG (ERROR) << "Error in reading the distortion correction file from CALDB" ;
            return (EXIT_FAILURE) ;
        }
        //caldb reading started for QE and MCP correction if QE and MCP correction to be done
        if (qe_mcpFlag)
        {
            tempname = caldb_handler.getQEFile (datainfo.getDetector () , datainfo.getObsMode () ,(char*) caldb_common_dir_qe) ;
           if (tempname == " ")
            {
                LOG (ERROR) << endl << "Couldn't find QEMCP file from caldb" ;
                continue ;
            }
            joinStrings (qe_factor_file , 2 ,(char*)caldb_common_dir_qe , tempname.c_str()) ;

          
            status = readNumRowsFromFits (qe_factor_file , 2 , nCalDBTempValues) ;
            if (status)
            {
                LOG (ERROR) << "Error in reading the number of rows from fits file " << qe_factor_file ;
                return (EXIT_FAILURE) ;
            }
            //array for storing the caldb values of QE and MCP correction temperature values and filter values.
            cal_temperature = new float[nCalDBTempValues] ;
            cal_f0 = new float[nCalDBTempValues] , cal_f1 = new float[nCalDBTempValues] ;   cal_f2 = new float[nCalDBTempValues] ;
            cal_f3 = new float[nCalDBTempValues] , cal_f4 = new float[nCalDBTempValues] ;   cal_f5 = new float[nCalDBTempValues] ;
            cal_f6 = new float[nCalDBTempValues] , cal_f7 = new float[nCalDBTempValues] ;
            LOG(INFO)<<"QE filename "<<qe_factor_file<<endl;
            status = readQEMCPFile (qe_factor_file , datainfo.getDetector () , nCalDBTempValues , cal_temperature , cal_f0 , cal_f1 , cal_f2 , cal_f3 , cal_f4 , cal_f5 , cal_f6 , cal_f7) ;
                 if (status)
            {
                LOG (ERROR) << "Error in reading the QEMCP caldb file ->" << qe_factor_file ;
                return (EXIT_FAILURE) ;
            }
            status = readNumRowsFromFits ((char *) dirobj.lbtfile.c_str () , 2 , nrows_lbt) ;
            if (status)
            {
                LOG (ERROR) << "Error in reading the number of rows from fits file" << dirobj.lbtfile ;
                return (EXIT_FAILURE) ;
            }

            time_lbt = new double[nrows_lbt] ;//Array for storing time of LBT file.
            insideTemp = new float [nrows_lbt] , outsideTemp = new float[nrows_lbt] ;//Array for storing inside temperatre and outside temperature of LBT file.
            
            
            
             /*===========================================================================================================================================*/
             /*=========================================CALDB READING FINISHED================================================================================*/
             /*===========================================================================================================================================*/
            //reading time,inside temperature and outside temperature from LBT file.
            status = getTemp ((char *) dirobj.lbtfile.c_str () , datainfo.getDetector () , time_lbt , insideTemp , outsideTemp , nrows_lbt) ;
            if (status)
            {
                LOG (ERROR) << "***temperature reading from the lbt file unsuccessful***" ;
                continue ;
            }
         
            int filter_coln ;
            int filternumber ;
            sprintf (filter , "%s" , datainfo.getFilter ()) ;
            LOG (INFO) << "The Filter " << filter ;

            //Array for storing factor values from caldb based on  filter value(i.e column in caldb)
            qe_mg_factor = new float[nCalDBTempValues] ;
           
            if (filter == (string) "F0")
            {
               
                filternumber = 0 ;
                filter_coln = 1 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = cal_f0[q] ;
            }
            else  if (filter == (string) "F1")
            {
                filternumber = 1 ;
                filter_coln = 2 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = cal_f1[q] ;
            }
            else if (filter == (string) "F2")
            {
                filternumber = 2 ;
                filter_coln = 3 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = cal_f2[q] ;
            }
            else if (filter == (string) "F3")
            {
                filternumber = 3 ;
                filter_coln = 4 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = cal_f3[q] ;
            }
            else if (filter == (string) "F4")
            {
                filternumber = 4 ;
                filter_coln = 5 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = cal_f4[q] ;
            }
            else if (filter == (string) "F5")
            {
                filternumber = 5 ;
                filter_coln = 6 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = cal_f5[q] ;
            }
            else if (filter == (string) "F6")
            {
                filternumber = 6 ;
                filter_coln = 7 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = cal_f6[q] ;
            }
            else if (filter == (string) "F7")
            {
                filternumber = 7 ;
                filter_coln = 8 ;
                for (int q = 0 ; q < nCalDBTempValues ; q ++)
                    qe_mg_factor[q] = cal_f7[q] ;
            }

            else
            {
                LOG (ERROR) << endl << "***Invalid filter option*** " << endl ;
                continue ;
            }

        }
        //incase of Darkframe subtraction to be done or not
        if (darkframe_flag)
        {
            status = takeDarkinfo () ;
            t_darkframestart = readDarkFrame (dstartpath , darkFramestart_data) ;
            t_darkframeend = readDarkFrame (dendpath , darkFrameend_data) ;
            darkCompute_array = new float[xsize * ysize] ;
        }

        //setting directory structure for output directory
        if (wtd_uc == 1)//for unit Conversion
        {
            status = setDirectoryStructure (moduleoutdir_uc , "SignalFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
        }
        if (wtd_ac == 1)// for Accumulation
        {
            status = setDirectoryStructure (moduleoutdir_ac , "SignalFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_ac , "ExposureFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
        }
        if (wtd_bp == 1)// for Mask Bad pixel
        {
            status = setDirectoryStructure (moduleoutdir_bp , "SignalFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_bp , "ExposureFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
        }
        if (wtd_cr == 1)// for cosmicRay corr
        {
            status = setDirectoryStructure (moduleoutdir_cr , "SignalFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_cr , "ExposureFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
        }
        if (wtd_ff == 1 && flatfieldFlag==1)// for flatfield corr
        {
            status = setDirectoryStructure (moduleoutdir_ff , "SignalFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_ff , "ExposureFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
        }
        if (wtd_qemcp && qe_mcpFlag==1) // for qemcp corr
        {
            status = setDirectoryStructure (moduleoutdir_qemcp , "SignalFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_qemcp , "ExposureFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
        }
        if (wtd_pp )// for pixpadding 
        {
            status = setDirectoryStructure (moduleoutdir_pp , "SignalFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_pp , "ExposureFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
        }
        if (wtd_sd == 1 && subdivisionFlag == 1)// for sub division
        {
            status = setDirectoryStructure (moduleoutdir_sd , "SignalFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_sd , "ExposureFrames") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
        }
        if (wtd_fsc == 1)// for DetectStar
        {
            status = setDirectoryStructure (moduleoutdir_sc , "Star") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
            status = setDirectoryStructure (moduleoutdir_sc , "Centroid") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }
        }
        if (wtd_dd == 1) // for Detector Distortion correction
        {
            status = setDirectoryStructure (moduleoutdir_dd , "Centroid") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;

            }
        }
        if (wtd_od == 1)// for Optical Distortion
        {
            status = setDirectoryStructure (moduleoutdir_od , "Centroid") ;
            if (status)
            {
                LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
                continue ;
            }

        }
        // for reference frame calculation
        status = setDirectoryStructure (moduleoutdir_rfc , "Centroid") ;
        if (status)
        {
            LOG (ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            continue ;
        }

       
        float *sum_ac = new float[padding_dim * padding_dim] ;
        float *sum_ac_exp = new float[padding_dim * padding_dim] ;
        initArray (sum_ac , padding_dim* padding_dim , (float) 0.0) ;
        initArray (sum_ac_exp , padding_dim* padding_dim , (float) 0.0) ;
        LOG (INFO) << "Total number of frames are " << nframes << endl ;
        if(nframes==0){
            LOG(INFO)<<"***NOT Sufficient frames found for input, Total number of frames Available are "<<nframes<<" ***Cant Proceed further***";
            continue;
            
        }
        double avg_time = 0.0 ;
        int count_frame_ac = 0 ;
        int cnt_forRefFrmCal = 0 ;
        int trackcnt_forRefFrame = 0 ;
        char *tform2[] = {"E" , "E" , "E"} ;
        int   tfields = 3 ;
        char *ttype[] = {"X" , "Y" , "Intensity"} ;
        //char *tform[] = {"I" , "I" , "E"} ;
        int cnt_temp = 0 ;
        double time_RefFrame_Avg = 0.0f ;
        fitsfile *fptr ;//pointer for pointing   fits file
        
        int cnt_finish = 0 ;
        bool flag_fsc = FALSE ;
        
        int Accu_fra_no = no_ofFramesToAcc ;
        int No_discard = frames_toDiscard ;
        double t1 , t2 , x1 , x2 , factor ;
        
        if (nframes <Accu_fra_no)//in case number of frames to accumulate is greater than total available frames.
        {
          LOG(WARNING)<< "Total number of frames to be accumulated are greater  than total available frames,Accumulation will be done on individual frame" << endl ;
            Accu_fra_no = 1 ;
        }
        if (nframes/no_ofFramesToAcc <= No_discard)//incase of number of frames to be discarded is greater than the total available frames.
        {
            LOG(WARNING) << "Total number of frames to be discarded are greater  than total available frames,No  frames will be discarded  " << endl ;
            No_discard = 0 ;
        }
        int tempx = 0 , tempy = 0 ;
        //initialize track_invalidPix to 0 as there is no invalid pixel initially.
        for (int k = 0 ; k < padding_dim ; k ++)
        {
            for (int l = 0 ; l < padding_dim ; l ++) track_invalidPix.push_back (0) ;
        }
       // Rows_tocheck rowscheck;
    int r1=44;
    int r2=108;
    int r3=172;
    int r4=236;
    int r5=300;
    int r6=364;
    int r7=428;
    int r8=492;
   
    int r1h=12;
    int r2h=76;
    int r3h=140;
    int r4h=204;
    int r5h=268;
    int r6h=332;
    int r7h=396;
    int r8h=460;
    bool flag_incompete=FALSE;
    bool flag_junkFrame=FALSE;
   counter_JunkFrame=0;
        //loop for total number of frames in level 1 science data file.     
        for (int i = 0 ; i < nframes ; i ++)
        {
            flag_junkFrame=FALSE;
            sprintf (errstr , "Error at iteration number %d" , i) ;
            sprintf (infile , "%s/%s/%s" , moduleIndir , "SignalFrames" , sigframelist[i]) ;

            frame_Data = (float*) malloc (xsize * ysize * sizeof (float) ) ;
            //frame_Data= new float[xsize*ysize];
           // frame_ExpData= new float[xsize*ysize];
            initArray (frame_Data , xsize*ysize , (float)INITIALIZATION_VALUE) ; //allocating memory to frame_Data to store image pixels of frames
           frame_ExpData = (float*) malloc (xsize * ysize * sizeof (float) ) ;
            initArray (frame_ExpData , xsize*ysize , (float)INITIALIZATION_VALUE) ; //allocating memory to frame_ExpData to store image pixels of  exposure frames

            //read signal frame image and store it in frame_Data.
            status = readImage (infile , 1 , frame_Data , xsize , ysize) ;
                    
            if(junkFrameFlag)
            {
                LOG(INFO)<<"Junk frame pixel correction started...";
              
                  status =editframe(frame_Data,xsize,ysize,thrJunkFrame,flag_junkFrame,r1);

                  if (flag_junkFrame==TRUE){
                      counter_JunkFrame++;
                  }
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,1,r1);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,1,r2);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,1,r3);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,1,r4);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,1,r5);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,1,r6);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,1,r7);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,1,r8);
//                
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,0,r1h);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,0,r2h);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,0,r3h);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,0,r4h);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,0,r5h);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,0,r6h);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,0,r7h);
//                status =editframe(frame_Data,MAX_XLOC_TOCHECKED,ysize,thrJunkFrame,0,r8h);
                
            }
           
             
            //opening DataIngest Signal frame
            fits_open_file (&fptr , infile , READONLY , &status) ;
            printError (status , "Error in opening the input file" , infile) ;
            fits_movabs_hdu (fptr , 1 , NULL , &status) ;
            printError (status , "Error in  moving to the 2nd HDU of the out information file" , infile) ;
            fits_read_key (fptr , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ;
            printError (status , "Error in  reading the FRAMENO keyword" , infile) ;
            fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
            printError (status , "Error in  reading the FRAMETIME keyword" , infile) ;
            fits_read_key (fptr , TDOUBLE , "INT_TIME" , &integrationtime , NULL , &status) ;
            printError (status , "Error in  reading the INT_TIME keyword" , infile) ;
            fits_close_file (fptr , &status) ;
            printError (status , "Error in  closing the file" , infile) ;

            //copy level1 keywords to header_info vector
            status = copyAllheaderKeys (infile) ;

            //for the Dark Subtraction process..      
            if (darkframe_flag)
            {
                t_curr = frametime ; //stores image frame's time
                status = darkFrameComputation (darkCompute_array) ; //intermediate dark value for the  current frame according to t_curr.
                if (status)
                {
                    LOG (ERROR) << "Error in Dark frame computation" ;
                    return (EXIT_FAILURE) ;
                }

                status = darkFrameSubtraction (darkCompute_array , frame_Data , xsize , ysize) ; //performing the dark frame subtraction.
                if (status)
                {
                    LOG (ERROR) << "***Error in Dark Frame Subtraction***" << endl ;
                    return (EXIT_FAILURE) ;
                }
            }

            //performing unit conversion (divide image pixels with the integration time)
            status = performUnitConversionIM (frame_Data , frame_ExpData , integrationtime , xsize , ysize) ;
            if (status)
            {
                LOG (INFO) << "ERROR in unitConversion " << endl ;
LOG(ERROR)<<"CRASH FRAME INTEGRATION TIME =< 0 ! (uvtRelativeAspectIM.cpp) ";
                flag_incompete=TRUE;
                break;
            }
            if (wtd_uc == 1)//in case of unit conversion output to be written to the disk
            {
                status = writeOutputImageToDisk ("uc" , moduleoutdir_uc , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ;
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                    return(EXIT_FAILURE);
                }
            }
           // exit(1);
            //performing badpixel filtering
            status = performCorrectionIM (frame_Data , frame_ExpData , badpixdata , xsize , ysize , integrationtime) ;
            if (status)
            {
                LOG (ERROR) << "***Error in  BAD pixel filtering***" << endl ;
                  flag_incompete=TRUE;
                break;
            }
            if (wtd_bp == 1)//in case of bad pixels output to be written to the disk
            {
                status = writeOutputImageToDisk ("bp" , moduleoutdir_bp , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ; //this is for the SignalFrame output
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                   return(EXIT_FAILURE);
                }
                status = writeOutputImageToDisk ("bp" , moduleoutdir_bp , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , xsize , ysize) ; //this is for the SignalFrame output
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                    return(EXIT_FAILURE);
                }
            }

            //performing cosmic ray correction
            status = performCosmicRayCorrIM (frame_Data , frame_ExpData , xsize , ysize , threshold) ;
            if (status)
            {
                LOG (ERROR) << "Error in Cosmic Ray correctiom" ;
                  flag_incompete=TRUE;
                break;
            }
            if (wtd_cr == 1)//in case of Cosmic ray correction  output to be written to the disk
            {
                status = writeOutputImageToDisk ("cr" , moduleoutdir_cr , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ; //this is for the SignalFrame output
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                  return(EXIT_FAILURE);
                }
                status = writeOutputImageToDisk ("cr" , moduleoutdir_cr , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , xsize , ysize) ; //this is for the SignalFrame output
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                   return(EXIT_FAILURE);
                }
            }

            /*process starts for the FlatField */
            if (flatfieldFlag)//incase of flat field correction to be done or not
            {
                //performing flat filed correction
                status = performFlatFieldCorrIM (frame_Data , flatfielddata , xsize , ysize) ;
                if (status)
                {
                    LOG (ERROR) << "Error in FlatField Correction" ;
                      flag_incompete=TRUE;
                    break;
                }
                if (wtd_ff == 1)//in case of flat field output to be written to the disk
                {
                    status = writeOutputImageToDisk ("ff" , moduleoutdir_ff , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ;
                    if (status)
                    {
                        LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                        return(EXIT_FAILURE);
                    }
                    status = writeOutputImageToDisk ("ff" , moduleoutdir_ff , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , xsize , ysize) ;
                    if (status)
                    {
                        LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                        return(EXIT_FAILURE);
                    }
                }
            }
            /*process Ends For the FlatField */
            if(UTC_flag==0 && qe_mcpFlag==1){
                LOG(WARNING)<<"***QEMCP correction cant be kept ON if UTC correction not to be done!!,NOW skipping QEMCP correction***";
                qe_mcpFlag=0;
            }
            /*process starts for the QE and MCP correction */
            if (qe_mcpFlag)//incase of QE and MCP correction to be done or not
            {

                temperature = INITIALIZATION_VALUE ; //initialization of temperature.
                for (int i = 0 ; i < nrows_lbt ; i ++)//loop for number of rows of LBT file for finding the relevant temperature of  current frame based on the frame time 
                {

                    if (frametime >= time_lbt[i] && frametime < time_lbt[i + 1])
                    {

                        temperature = (insideTemp[i] + outsideTemp[i]) / 2 ;
                        break ;
                    }
                }
                if (temperature == INITIALIZATION_VALUE)//incase of no relevant temperature found for the current frame 
                {
                    LOG (ERROR) << "***No record found in LBT file***" << endl ;
                    continue;
                }

                for (int j = 0 ; j < nCalDBTempValues - 1 ; j ++)//loop for finding the relevant factor value by comparing calculated temperature value and  caldb temperature value  
                {
                    if (temperature >= cal_temperature[j] && temperature < cal_temperature[j + 1])   //temperature - from LBT file ,cal_temperature from caldb file       
                    {
                        t1 = cal_temperature[j] ;
                        t2 = cal_temperature[j + 1] ;
                        if ((t2 - t1) == 0)
                        {
                            LOG (INFO) << "***Divide By zero***" << endl ;
                            return (EXIT_FAILURE) ;
                        }
                        x1 = qe_mg_factor[j] ;
                        x2 = qe_mg_factor[j + 1] ;
                        factor = x1 + ((temperature - t1)*((x2 - x1) / (t2 - t1))) ; //calculated factor
                       
                    }
                }
                
                //performing QE and MCP correction(by applying factor to image pixels)
                status = performQEMCPcorrection (frame_Data , xsize , ysize , factor) ;
                if (status)
                {
                    LOG (ERROR) << "Error in QEMCP module" ;
                     flag_incompete=TRUE;
                    break;
                }
                if (wtd_qemcp == 1)//incase of QE and MCP correction to be written to the disk
                {
                    status = writeOutputImageToDisk ("qe" , moduleoutdir_qemcp , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ;
                    if (status)
                    {
                        LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                        return(EXIT_FAILURE);
                    }
                    status = writeOutputImageToDisk ("qe" , moduleoutdir_qemcp , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , xsize , ysize) ;
                    if (status)
                    {
                        LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                       return(EXIT_FAILURE);
                    }
                }

            }
            /*process ends  for the QE and MCP correction */

            /*process starts for the pixel padding */
            for (int i = 0 ; i < padding_dim * padding_dim ; i ++)//loop for initializing frame_Data_Padded and frame_ExpoData_padded array (i.e 600*600)
            {
                frame_Data_Padded[i] = INITIALIZATION_VALUE ;
                frame_ExpoData_padded[i] = INITIALIZATION_VALUE ;
            }

            //performing pixel padding
            status = Applypadding (frame_Data , xsize , ysize , frame_Data_Padded , padding_dim , padding_dim) ;
            if (status)
            {
                LOG (ERROR) << "Error in applying padding" ;
                  flag_incompete=TRUE;
                break;
            }
            status = Applypadding (frame_ExpData , xsize , ysize , frame_ExpoData_padded , padding_dim , padding_dim) ;
            if (status)
            {
                LOG (ERROR) << "Error in applying padding " ;
                  flag_incompete=TRUE;
                break;
            }
            
           // delete[] frame_Data,frame_ExpData;
            free (frame_Data) ; //releasing the memory
            free (frame_ExpData) ; //releasing the memory

        //    frame_Data= new float[padding_dim*padding_dim];
         //   frame_ExpData= new float[padding_dim]
            frame_Data = (float*) malloc (padding_dim * padding_dim * sizeof (float) ) ;
            initArray (frame_Data , padding_dim*padding_dim , (float) INITIALIZATION_VALUE) ; //allocating memory to frame_Data
            frame_ExpData = (float*) malloc (padding_dim * padding_dim * sizeof (float) ) ;
            initArray (frame_ExpData , padding_dim*padding_dim , (float) INITIALIZATION_VALUE) ; //allocating  memory to frame_ExpData


            for (int pix = 0 ; pix < padding_dim * padding_dim ; pix ++)//loop for assigning the values of arrays to frame_Data and frame_ExpData.
            {
                frame_Data[pix] = frame_Data_Padded[pix] ;
                frame_ExpData[pix] = frame_ExpoData_padded[pix] ;
            }

            if (wtd_pp)//incase of  pixel padding to be written to the disk
            {
                status = writeOutputImageToDisk ("pp" , moduleoutdir_pp , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , padding_dim , padding_dim) ;
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                  return(EXIT_FAILURE);
                }
                status = writeOutputImageToDisk ("pp" , moduleoutdir_pp , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , padding_dim , padding_dim) ;
                if (status)
                {
                    LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                    return(EXIT_FAILURE);
                }
            }
            /*process ends for the pixel padding  */


            tempx = 0 , tempy = 0 ;

            //performing Accumulation on  total NAcc frames 
            for (int k = tempx = (padding_dim - (padding_dim * IMG_DIM_DI) / PIX_PADD_SIZE) / 2 ; k < tempx + IMG_DIM_DI * padding_dim / PIX_PADD_SIZE ; k ++)
            {
                for (int l = tempy = (padding_dim - (padding_dim * IMG_DIM_DI) / PIX_PADD_SIZE) / 2 ; l < tempy + IMG_DIM_DI * padding_dim / PIX_PADD_SIZE ; l ++)
                {
                    if (frame_Data[k * padding_dim + l] != INVALID_PIX_VALUE  && frame_ExpData[k * padding_dim + l] != INVALID_PIX_VALUE)
                    {
                        sum_ac[k * padding_dim + l] = sum_ac[k * padding_dim + l] + frame_Data[k * padding_dim + l] ;
                        sum_ac_exp[k * padding_dim + l] = sum_ac_exp[k * padding_dim + l] + frame_ExpData[k * padding_dim + l] ;
                    }
                    else
                    {
                        track_invalidPix[k * padding_dim + l] = track_invalidPix[k * padding_dim + l] + 1 ;

                    }
                }
            }

            avg_time = avg_time + frametime ;      //Accumulate time of total Nacc frames 

            //check for accumulation.
            if ((i + 1) % Accu_fra_no == 0  || i+1==nframes)
            {
                count_frame_ac ++ ;
                if(i+1==nframes && (i + 1) % Accu_fra_no !=0) Accu_fra_no=(i + 1) % Accu_fra_no;
                LOG(INFO)<<"Accumulated frame count"<<Accu_fra_no<<endl;;
                /*process starts for the Accumulation   */
                //performing average on pixels  of   accumulated Nacc frames.
                status = performAccOFFrame (frame_Data , sum_ac , padding_dim , padding_dim , Accu_fra_no) ; //for Signal frames
                if (status)
                {
                    LOG (ERROR) << "Error in performing Accumulation" ;
                      flag_incompete=TRUE;
                    break;
                }

                status = performAccOFFrame (frame_ExpData , sum_ac_exp , padding_dim , padding_dim , Accu_fra_no) ; //for exposure frames
                if (status)
                {
                    LOG (ERROR) << "Error in performing Accumulation" ;
                      flag_incompete=TRUE;
                    break;
                }
                //Added
                
    for (int k =0 ; k <  xsize ; k ++)
    {
        for (int l = 0 ; l <  ysize; l ++)
        {
//%#Added ON 20July-divide by zero#%  
             if((Accu_fra_no - track_invalidPix[k * xsize + l] )==0){
                frame_Data[k * xsize + l]=INVALID_PIX_VALUE;
            }
//%#-Till this-20July17#%
            else{
            if( frame_Data[k * xsize + l]!=INVALID_PIX_VALUE)
             //   frmsigdata[k * sizex + l] = (float) (sumdata[k * sizex + l] / (numoffrmAcc - track_invalidPix[k * sizex + l])) ;
                 frame_Data[k * xsize + l] = (float) (frame_Data[k * xsize + l] )/(Accu_fra_no - track_invalidPix[k * xsize + l]) ; ;
          //  }
            }
            
        }
    }
                
                //till this
                track_invalidPix.clear () ;

                //no invalid pixel now 
                for (int k = 0 ; k < padding_dim ; k ++)
                {
                    for (int l = 0 ; l < padding_dim ; l ++)
                    {

                        track_invalidPix.push_back (0) ;

                    }
                }
                //vect_cnt_duplicates.clear ();
                frameno = count_frame_ac ;

                frametime = avg_time / Accu_fra_no ; //calculating frametime for Accumulated frame by averaging the avg_time
                avg_time = 0.0 ;

                if (wtd_ac == 1)//incase of Accumulation to be written on disk
                {
                    status = writeOutputImageToDisk ("ac" , moduleoutdir_ac , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , padding_dim , padding_dim) ;
                    if (status)
                    {
                        LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                       return(EXIT_FAILURE);
                    }
                    status = writeOutputImageToDisk ("ac" , moduleoutdir_ac , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , padding_dim , padding_dim) ;
                    if (status)
                    {
                        LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                        return(EXIT_FAILURE);
                    }
                }
                /*process ends for the Accumulation  */

                /*process starts for the Sub Division  */

                if (subdivisionFlag)//flag to decide sub division to be done or not
                {
                    initArray (frame_Data_subdivided , subDivision_size*subDivision_size , (float) INITIALIZATION_VALUE) ; //initialization of frame_Data
                    initArray (frame_ExpData_subdivided , subDivision_size*subDivision_size , (float) INITIALIZATION_VALUE) ; //initialization of frame_ExpData

                    //divide the pixel Array to division factor   because of 1 pixel to be distributed  into 8*8 window(incase of Subdivision to be done)
//%#Added ON 20July-devide by zero#%  .
			  if(division_fact==0){
                        LOG(ERROR)<<"***Divide by ZERO error***";
                        break;
                    }
//%#-Till this-20July17#%
                    for (int p = 0 ; p < padding_dim * padding_dim ; p ++)
                    {
                        if (frame_Data[p] != INVALID_PIX_VALUE)
                        {
                            frame_Data[p] = frame_Data[p] / division_fact ;
                        }
                    }
                    //performing Sub division
                    status = performSubDivisionIM (frame_Data , padding_dim , padding_dim , frame_Data_subdivided , subDivision_size , subDivision_size) ;
                    if (status)
                    {
                        LOG (ERROR) << "Error in subdivision module" << endl ;
                          flag_incompete=TRUE;
                        break;
                    }
                    status = performSubDivisionIM (frame_ExpData , padding_dim , padding_dim , frame_ExpData_subdivided , subDivision_size , subDivision_size) ;
                    if (status)
                    {
                        LOG (ERROR) << "Error in subdivision module" << endl ;
                          flag_incompete=TRUE;
                        break;
                    }

                   free (frame_Data) ; //releasing the memory
                    free (frame_ExpData) ; //releasing the memory
                  //  delete[] frame_Data,frame_ExpData;
                    frame_Data = (float*) malloc (subDivision_size * subDivision_size * sizeof (float) ) ; //Allocating the memory to frame_Data
                    frame_ExpData = (float*) malloc (subDivision_size * subDivision_size * sizeof (float) ) ; //Allocating the memory to frame_ExpData

                    status = initArray (frame_Data , subDivision_size*subDivision_size , (float) INITIALIZATION_VALUE) ;
                    status = initArray (frame_ExpData , subDivision_size*subDivision_size , (float) INITIALIZATION_VALUE) ;

                    for (int pix = 0 ; pix < subDivision_size * subDivision_size ; pix ++)
                    {
                        frame_Data[pix] = frame_Data_subdivided[pix] ;
                        frame_ExpData[pix] = frame_ExpData_subdivided[pix] ;
                    }

                    if (wtd_sd)//incase of Sub division to be written on the Disk.
                    {
                        status = writeOutputImageToDisk ("sd" , moduleoutdir_sd , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                        if (status)
                        {
                            LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                            return(EXIT_FAILURE);
                        }
                        status = writeOutputImageToDisk ("sd" , moduleoutdir_sd , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                        if (status)
                        {
                            LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                           return(EXIT_FAILURE);
                        }
                    }
                }


                rms_mul_factor = rms_mul_factor_default ;

                if (count_frame_ac > 1)            flag_fsc = TRUE ;

                /*process started for Star detection*/
                status = obj_sc.getStarCentroid (frame_Data  , Accu_fra_no , backgrd_fact , rms_mul_factor , subDivision_size , subDivision_size , refine_Winsize * (subDivision_size / padding_dim) ,
                        centroid_Winsize * (subDivision_size / padding_dim)  , Rx_vect , Ry_vect , cx_vect , cy_vect , Fx_vect , Fy_vect , R_int ,
                        c_int , F_int , datainfo.getObsMode () , min_num_stars , flag_fsc , search_win_size , star_detect_algo_flag) ;
                if (status)
                {
                    LOG (ERROR) << "Error in finding stars and centroids for frame" ;
                    LOG(ERROR)<<"CRASH ERROR IN FINDING STARS (POSSIBLY MULTIPLIER TO SIGMA =< 0) (uvtRelativeAspectIM.cpp) ";
                      flag_incompete=TRUE;
                    break ;
                }


                centroid_row = cx_vect.size () ;
            
                if (wtd_fsc == 1)//incase of  star detection to be written on disk or not
                {
                    status = writeOutputImageToDisk ("sc" , moduleoutdir_sc , "Star" , "Star" , frame_Data , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                    if (status)
                    {
                        LOG (ERROR) << "***Writing to Disk Fails***" << endl ;
                       return(EXIT_FAILURE);
                    }
                    status = writeOutputTblToDisk ("sc" , moduleoutdir_sc , "Centroid" , "centroid" , nameprefix , frametime , frameno , ttype , tform2 , tfields , "Centroids" , cx_vect , cy_vect , c_int) ;
                    if (status)
                    {
                        LOG (ERROR) << "Error in writing the table to the disk for DetectStar module" ;
                        return(EXIT_FAILURE);
                    }
                }
              
                /*process ends for Star detection*/
        //      if(i==1)  exit(1);
                //clearing the vectors of refined peaks and First-cut peaks.           
                Fx_vect.clear () ;
                Fy_vect.clear () ;
                F_int.clear () ;
                Rx_vect.clear () ;
                Ry_vect.clear () ;
                R_int.clear () ;

                /*process starts for detector distortion*/
                //applying distortion correction
                status = performDistortionCorr (cx_vect , cy_vect , X_detect_distortion , Y_detect_distortion , subDivision_size , xsize) ;
                if (status)
                {
                    LOG (ERROR) << "***Error in detector distortion module***" ;
                      flag_incompete=TRUE;
                    break ;
                }
                if (wtd_dd == 1)//incase of detector distortion to be written  to disk
                {
                    status = writeOutputTblToDisk ("dd" , moduleoutdir_dd , "Centroid" , "centroid" , nameprefix , frametime , frameno , ttype , tform2 , tfields , "Centroids" , cx_vect , cy_vect , c_int) ;
                    if (status)
                    {
                        LOG (ERROR) << "Error in performing detector distortion correction" ;
                        return(EXIT_FAILURE);
                    }
                }
                 
                /*process ends for detector distortion*/

                /*process starts for soptical Assembly distortion correction*/

                //applying optical  distortion correction
                status = performDistortionCorr (cx_vect , cy_vect , X_optical_distortion , Y_optical_distortion , subDivision_size , xsize) ;
                if (status)
                {
                    LOG (ERROR) << "Error in performing optical distortion correction" ;
                      flag_incompete=TRUE;
                    break;
                }
                if (wtd_od == 1)//incase of optical Distortion correction to be written to the dosk or not
                {
                    status = writeOutputTblToDisk ("od" , moduleoutdir_od , "Centroid" , "centroid" , nameprefix , frametime , frameno , ttype , tform2 , tfields , "Centroids" , cx_vect , cy_vect , c_int) ;
                    if (status)
                    {
                        LOG (ERROR) << "Error in optical distortion correction process" ;
                        return(EXIT_FAILURE);
                    }
                }
                /*process ends for optical Assembly distortion correction*/

                /*process starts for reference frame calculation*/
                if (cnt_forRefFrmCal >= No_discard)
                {

          //         if ((trackcnt_forRefFrame % nFrameToAverage != 0 || (nFrameToAverage == 1) || trackcnt_forRefFrame >= 0))
                //        if (((trackcnt_forRefFrame % nFrameToAverage != 0 || (nFrameToAverage == 1) ) && trackcnt_forRefFrame == 0))
                  //  {
                        size_rows.push_back (cx_vect.size ()) ;
                        ref_frame_time_data.push_back (frametime) ;


                        for (int p = 0 ; p < cx_vect.size () ; p ++)
                        {
                            x_Of_refFrame.push_back (cx_vect[p]) ;
                            y_Of_refFrame.push_back (cy_vect[p]) ;
                            int_Of_refFrame.push_back (c_int[p]) ;

                        }
                        cnt_temp ++ ;
                 //   }
                   // if (cnt_temp == nFrameToAverage  )//after or it is added.
                          if (cnt_temp >= nFrameToAverage  )//after or it is added.
                 //   if (x_Of_refFrame.size ()!=0  )//after or it is added.
                    {
                        if(flag_rfc==TRUE) nFrameToAverage=1;
                        //calculating reference frame
                        status = calculateRefFrame (size_rows , x_Of_refFrame , y_Of_refFrame , int_Of_refFrame , nFrameToAverage , x_ref_cumm , y_ref_cumm , int_ref_cumm) ;
                        if (status)
                        {
                            LOG (ERROR) << "Error in calculating reference frame" ;
                            break ;
                        }
                         flag_rfc=TRUE;
                        cx_vect = x_ref_cumm ;
                        cy_vect = y_ref_cumm ;
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

//                        status = writeOutputTblToDisk ("rf" , moduleoutdir_rfc , "Centroid" , "centroid" , nameprefix , time_RefFrame_Avg , cnt_forRefFrmCal-nFrameToAverage+1 , ttype , tform2 , tfields , "Centroids" , cx_vect , cy_vect , c_int) ;
                         status = writeOutputTblToDisk ("rf" , moduleoutdir_rfc , "Centroid" , "centroid" , nameprefix , time_RefFrame_Avg , trackcnt_forRefFrame+1 , ttype , tform2 , tfields , "Centroids" , cx_vect , cy_vect , c_int) ;
                        if (status)
                        {
                            LOG (ERROR) << "Error in writing out table  to the disk" ;
                            break ;
                        }
                        time_RefFrame_Avg = 0.0 ;
                        ref_frame_time_data.clear () ;
                      //  cnt_temp = 0 ;//commented
                         trackcnt_forRefFrame ++ ;
                    }
                //   trackcnt_forRefFrame ++ ;
                }
               
                cnt_forRefFrmCal = cnt_forRefFrmCal + 1 ;
                cx_vect.clear () ;
                cy_vect.clear () ;
                ci_tempi.clear () ;
                
                LOG (INFO) << "\033[1;34mProcess ends for  frame number -> " << cnt_finish << "\033[0m" ;
                cnt_finish ++ ;
//                if (datainfo.getIntegrationTime () >= 1)
//                {
//                    Accu_fra_no = 1 ;
//                }

            }
            //releasing memory
           free (frame_Data) ;
           free (frame_ExpData) ;
           
            //delete[] frame_Data,frame_ExpData;
        }
        if (flatfieldFlag)
        {
            delete[] flatfielddata ;
        }
    if(  flag_incompete==TRUE){
        LOG(INFO)<<"CHECKING NEXT VIS FILE";
        continue;
    }
        //creating output information file to the reference frame calculation.
        char infofile_out[NAMESIZE] ;
        sprintf (infofile_out , "%s/%s_rfc.info" , moduleoutdir_rfc , nameprefix) ;
        fits_create_file (&finfo_out , infofile_out , &status) ;
        printError (status , "Error in creating the output information file") ;
         char temp[NAMESIZE] ;
//        for (int p = 0 ; p < header_info.size () ; p ++)
//        {
//            
//            sprintf (temp , "%s" , header_info[p].c_str ()) ;
//            fits_write_record (finfo_out , (char*) &temp , &status) ;
//        }
//        
        
       
        char *ttype1[] = {"centroidFileList"} ;
        char *tform1[] = {"A256"} ;
        fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype1 , tform1 , NULL , "FileList" , &status) ;
        printError (status , "Error in Creating the table ***") ; //for creating name for output information file
        fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
        printError (status , "Error in updating the key value of the NAMEPRFX***") ; //for creating name for output information file
        char dir[100] = "Centroid" ;
        fits_update_key (finfo_out , TSTRING , "centroidDir" , (char*) dir , NULL , &status) ;
        printError (status , "Error in updating the  key value of the centroidDir") ;
        fits_update_key (finfo_out , TSTRING , "OBS_MODE" , obsmode , NULL , &status) ;
        printError (status , "Error in updating the  key value of the OBS_MODE") ;

        int size_vect = ref_frame_module_filename_track.size () ; //total number of frames that are generated after reference frame calculation 
        //cout<<"Size of vect"
        string *name_ofFile = new string[size_vect] ;
        for (int p = 0 ; p < size_vect ; p ++)
        {
            name_ofFile[p] = ref_frame_module_filename_track[p] ;
            //LOG(INFO)<<name_ofFile[p]<<endl;
        }
        //updating/writing  the key value of NFILES keyword
        fits_update_key (finfo_out , TINT , "NFILES" , &size_vect , NULL , &status) ;
        printError (status , "Error in updating the key value of the NFILES") ;
        //writing basic information to reference frame calculation information file

        datainfo.write (finfo_out) ;

        //updaiting keyvalue of XSIZE and YSIZE.
        fits_update_key (finfo_out , TINT , "XSIZE" , &subDivision_size , NULL , &status) ;
        printError (status , "Error in updating the  key value of the XSIZE") ;
        fits_update_key (finfo_out , TINT , "YSIZE" , &subDivision_size , NULL , &status) ;
        printError (status , "Error in updating the  key value of the YSIZE ") ;
        
        //Writing  total number of frames to the output information file
        
for(int i =0;i<size_vect;i++)
{
        fits_write_col (finfo_out , TSTRING , 1 , i+1 , 1 , 1 , &name_ofFile[i] , &status) ;
        printError (status , "Error in writing the column of outout Exposure frame list") ;
}

//fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , size_vect , name_ofFile , &status) ;
        //printError (status , "Error in writing the column of outout Exposure frame list") ;
         fits_movabs_hdu (finfo_out , 1 , NULL , &status) ;
        printError (status , "Error in moving to 1st HDU",infofile_out) ;
        fits_update_key (finfo_out , TLONG , "JunkFrames_Counts" , &counter_JunkFrame , NULL , &status) ;
        printError (status , "Error in updating the  key value of the YSIZE ") ;
        
        fits_close_file (finfo_out , &status) ;
        printError (status , "Error in Closing the output information file") ;
// writeHistory (infofile_out,header_info);
        writeUsrkeywordsFrmvect(infofile_out,L1keywords);
        //writeHistory (infofile_out,header_info);
        char drift_outDir[FLEN_FILENAME] ;
        /*process ends for the reference frame calculation*/

        /*process starts for Drift Computation*/
        //Drift Computation
        LOG(INFO)<<level1indir;
        string dir_aux4="AUX4_"+(string)basename(level1indir.c_str())+"_"+convertIntToStr(dataindex+1);
        string cmd_aux4="mkdir -p "+dir_aux4;
        system (cmd_aux4.c_str ());
        LOG (ERROR) << endl << "==================================UVTCOMPUTEDRIFT==================================================" << endl ;
        uvtDriftComputation drift_obj ;
        status = drift_obj.read (moduleoutdir_rfc ,(char*)dirobj.attfile.c_str (),err_per , nbhd_dist , freqDomainFilter_Flag , type_Filtering , delta_time , freqvalue , fitting_flag , orderpitch , orderyaw , orderroll , outputdir , shift_rotation_algo , match_stars_file_flag , star_detect_algo_flag , flag_thetaComp,clobber , history) ;
        if (status)
        {
            LOG (ERROR) << "Error in drift calculation" ;
	    LOG(ERROR)<<"CRASH ERROR IN FINDING DRIFT (uvtRelativeAspectIM.cpp) ";
            break;
        }
        drift_obj.display () ;
        status = drift_obj.uvtDriftComputationProcess () ;
        if (status)
        {
            LOG (ERROR) << endl << "***Error in DriftComputation  Module***" << endl ;
            break;
        }
        strcpy (moduleIndir , drift_obj.getModuleOutdir ()) ;
        strcpy (drift_outDir , moduleIndir) ;
        cmd_aux4=" cp  "+(string)drift_outDir+(string)"/* "+(string)dir_aux4;
        LOG(INFO)<<cmd_aux4;
        system (cmd_aux4.c_str ());
        /*process ends for Drift Computation*/
//       LOG (ERROR) << endl << "==================================UVTCOMPUTEJITTER==================================================" << endl ;
//        /*process starts for Jitter Computation*/
//        char jitter_outDir[FLEN_FILENAME] ;
//        uvtComputeJitter jitter_obj ;
//        jitter_obj.read (moduleIndir , (char*) caldbindir.c_str () , (char *) dirobj.gyrofile.c_str () , freqDomainFilter_Flag , freqvalue , fitting_flag , orderpitch , orderyaw , orderroll , outputdir , type_Filtering , clobber , history) ;
//        jitter_obj.display () ;
//        status = jitter_obj.uvtComputeJitterProcess () ;
//        if (status)
//        {
//            LOG (INFO) << endl << "\033[1;31m***Error in ComputeJitter  module***\033[0m" << endl ;
//            return (EXIT_FAILURE) ;
//        }
//        strcpy (moduleIndir , jitter_obj.getModuleOutdir ()) ;
//        strcpy (jitter_outDir , moduleIndir) ;
//        /*process ends for Drift Computation*/
//        LOG (INFO) << endl << "Using directory " << moduleIndir << "  as input to uvtComputeThermal " << endl ;

        /*process starts for thermal calculation*/
        //----------------------THERMAL SERIES-------------//
        

        //uvtComputeThermal should not be run for VIS mode data.
        
//        if (strcasecmp (datainfo.getDetector () , "VIS") != 0)
//        {
//            LOG (INFO) << endl << "\033[1;34m==================================UVTCOMPUTETHERMAL==================================================\033[0m" << endl ;
//            uvtComputeThermal  thermal_obj ;
//            thermal_obj.read (moduleIndir , Indir_dataIngest , (char*) caldbindir.c_str () , (char*) dirobj.lbtfile.c_str () , outputdir , clobber , history) ;
//            thermal_obj.display () ;
//            status = thermal_obj.uvtThermalCalcProcess () ;
//            if (status)
//            {
//                LOG (INFO) << endl << "\033[1;31m***Error in uvtComputeThermal  module***\033[0m" << endl ;
//                return (EXIT_FAILURE) ;
//            }
//
//            strcpy (moduleIndir , thermal_obj.getModuleOutdir ()) ;
//        }
//
//        /*process ends for thermal Computation*/
//        LOG (INFO) << endl << "Using directory " << moduleIndir << "  as input to uvtRelAspCal " << endl ;

        /*process starts for Relative Aspect calculation*/
//        LOG (INFO) << endl << "\033[1;34m==================================UVTCOMPUTERELASPECT==================================================\033[0m" << endl ;
//        uvtRelAspCal  ras_obj ;
//        ras_obj.read (drift_outDir , jitter_outDir , moduleIndir , outputdir , clobber , history) ;
//        ras_obj.display () ;
//        status = ras_obj.uvtRelAspCalProcess () ;
//        if (status)
//        {
//            LOG (INFO) << endl << "\033[1;31m***Error in uvtRelativeAspCal  module***\033[0m" << endl ;
//            return (EXIT_FAILURE) ;
//        }

//        strcpy (RelAspFilename , ras_obj.getModuleOutdir ()) ;
        /*process ends for RelativeAspect Computation*/
        delete[] sum_ac,sum_ac_exp,flatfielddata,badpixdata,X_detect_distortion,Y_detect_distortion,X_optical_distortion,Y_optical_distortion,cal_f0,cal_f1,cal_f2,cal_f3,cal_f4,cal_f5,
                cal_f6,cal_f7,cal_temperature,time_lbt,insideTemp,outsideTemp,qe_mg_factor;
    }
//         total_tars++;
    
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectIM::setDirectoryStructure (char *Dir , const char *subdir)
{
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s/" , Dir , subdir) ;
    //LOG(INFO)<<dir<<endl;
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


int uvtRelativeAspectIM::writeOutputImageToDisk (char *id , char *outDir , char *dir , char *subscript , float *Array , char *namepre , double ftime , unsigned short fno , int sizex , int sizey)
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
              if(strstr(header_info[p].c_str () ,"NAXIS")==NULL){
            fits_write_record (fout , (char*) &temp , &status) ;
              }
        }
    }
    fits_close_file (fout , &status) ;
    printError (status , "Error in closing the  output Signal fits file" , outfile) ;
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectIM::copyAllheaderKeys (char* infile)
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
        LOG (INFO) << endl << "***Could not find number of keywords in file " << infile << "***" ;
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
            LOG (INFO) << endl << "***Error in reading record number " << i << " in file " << infile << "***" ;
            fits_report_error (stderr , status) ;
            return (EXIT_FAILURE) ;
        }
       // keyclass = fits_get_keyclass (record) ;
       // if (keyclass == TYP_COMM_KEY)
       //     continue ;
       // else if (keyclass == TYP_USER_KEY || keyclass == TYP_REFSYS_KEY || keyclass == TYP_WCS_KEY)
       // {
            header_info.push_back (record) ;
       // }
    }

    fits_close_file (fin , &status) ;
    // delete[] record;
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectIM::transformToUVITFrame (double *t , double *r , double *p , double *y)
{
    LOG (INFO) << "the 1,1" << endl ;
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectIM::transformToSpacecraftFrame (double *t , double *r , double *p , double *y)
{
    LOG (INFO) << "the 1,1" << endl ;
    return (EXIT_SUCCESS) ;
}


double uvtRelativeAspectIM::readDarkFrame (char * path , float *Array)
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


int uvtRelativeAspectIM::darkFrameComputation (float *outArray)
{
    for (int i = 0 ; i < xsize * ysize ; i ++)
    {
        outArray[i] = 0.0 ;
        float d1 = darkFramestart_data[i] ;
        float d2 = darkFrameend_data[i] ;
        double t1 = t_darkframestart ;
        double t2 = t_darkframeend ;
        float d;
        if ((t2 - t1) == 0)
        {
            d=d1;
            //LOG (ERROR) << "***Divide By zero Error in unit conversion ***" << endl ;
            
        }else
        {
         d = d1 + ((d2 - d1) *((t_curr- t1) / (t2 - t1))) ;
        //LOG(INFO)<<d<<endl;
        }
        outArray[i] = (float) d ;
        // outArray[i] = d1;

    }
    return (EXIT_SUCCESS) ;

}







int uvtRelativeAspectIM::performAccOFFrame (float *frmsigdata , float *sumdata , int sizex , int sizey , int numoffrmAcc)
{

    //LOG(INFO)<<track_invalidPix[45*600+45];
    if (numoffrmAcc == 0)
    {
        LOG (INFO) << "Divide by Zero" << endl ;
        return (EXIT_FAILURE) ;
    }
    //    for (int pix = 0 ; pix < sizex*sizey ; pix++)
    //            {
    //                frmsigdata[pix] = -9999;
    //                if(sumdata[pix]!=INVALID_PIX_VALUE){
    //                frmsigdata[pix] = sumdata[pix] / numoffrmAcc ;
    //                }
    //                sumdata[pix] = 0.0 ;
    //                
    //            }
    // LOG(INFO)<<sizex<< ""<<sizey<<" "<<track_invalidPix.size ()<<endl;
    int tempx = 0 , tempy = 0 ;
    int num_times_repeat = 0 ;
    //      for (int k = tempx=PADD_DEFAULT*sizex/PIX_PADD_SIZE ; k < tempx+IMG_DIM_DI*sizex/PIX_PADD_SIZE ; k++)
    //         {
    //         
    //             for(int l=tempy=PADD_DEFAULT*sizey/PIX_PADD_SIZE;l<tempy+IMG_DIM_DI*sizey/PIX_PADD_SIZE;l++)
    //             {
    //                 //LOG(INFO)<<"SIZE::"<<tempy<<" "<<" "<<tempy+IMG_DIM_DI*sizey/PIX_PADD_SIZE;
    ////                 if(sumdata[k*sizex+l]<0){
    ////                 LOG(INFO)<<sumdata[k*sizex+l]<<" "<<k<<" "<<l;
    ////                 exit(1);
    ////                 }
    //                      num_times_repeat = count (track_invalidPix.begin () , track_invalidPix.end () , k*sizex+l) ; //STL method for count
    //                     frmsigdata[k*sizex+l] = (float) (sumdata[k*sizex+l] / (numoffrmAcc- num_times_repeat)) ;
    //                     sumdata[k*sizex+l]=0.0f;
    //             }
    //         }
    //LOG(INFO)<<"EEE"<<endl;exit(1);
    for (int k = tempx = (sizex - (sizex * IMG_DIM_DI) / PIX_PADD_SIZE) / 2 ; k < tempx + IMG_DIM_DI * sizex / PIX_PADD_SIZE ; k ++)
    {
        for (int l = tempy = (sizey - (sizey * IMG_DIM_DI) / PIX_PADD_SIZE) / 2 ; l < tempy + IMG_DIM_DI * sizey / PIX_PADD_SIZE ; l ++)
        {
            //  num_times_repeat = count (track_invalidPix.begin () , track_invalidPix.end () , k*sizex+l) ; //STL method for count
            frmsigdata[k * sizex + l] = INVALID_PIX_VALUE ;
            if (track_invalidPix[k * sizex + l] == numoffrmAcc)
            {
                frmsigdata[k * sizex + l] = INVALID_PIX_VALUE ;
            }
            else
            {             
             //   frmsigdata[k * sizex + l] = (float) (sumdata[k * sizex + l] / (numoffrmAcc - track_invalidPix[k * sizex + l])) ;
                 frmsigdata[k * sizex + l] = (float) (sumdata[k * sizex + l] ) ;
            }
            sumdata[k * sizex + l] = 0.0 ;
        }
    }
   

    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectIM::writeOutputTblToDisk (char *id , char *outDir , char *dir , char *subscript  , char *namepre , double ftime , unsigned short fno , char **type1 , char**tform1 , int tfields1 , char *extname , vector<float> &X , vector<float> &Y , vector<float> &val)
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


int uvtRelativeAspectIM::performDistortionCorr (vector<float> &X , vector<float> &Y , float *Xdistr , float *Ydistr , int sizex , long caldbsize)
{
   float multi_factor = sizex / 600 ;
    if (multi_factor == 0)
    {
        LOG (ERROR) << "***Divide by Zero***" << endl ;
        return (EXIT_FAILURE) ;
    }
    float tempx , tempy ;
    double locate = 0.0f ;
    long curr_loc=0;
    long next_loc=0;
    double value_X_toAdd,value_Y_toAdd;
    double locate1;
    // LOG(INFO)<<caldbsize<<endl;exit(1);
    for (int p = 0 ; p < X.size () ; p ++)
    {
        if (X[p] != INVALID_PIX_VALUE && Y[p] != INVALID_PIX_VALUE)
        {
            tempx = (X[p] / multi_factor)-44 ;
            tempy = (Y[p] / multi_factor)-44;

             locate = ((int) round (tempy)) * caldbsize + ((int) round (tempx)) ;
           
               curr_loc=(int)locate-1;
               next_loc=curr_loc+1;
               if(curr_loc>0 && curr_loc<512*512 && next_loc > 0 && next_loc<512*512){
               value_X_toAdd=Xdistr[curr_loc]+((Xdistr[next_loc]-Xdistr[curr_loc]))*(locate-curr_loc);
               value_Y_toAdd=Ydistr[curr_loc]+((Ydistr[next_loc]-Ydistr[curr_loc]))*(locate-curr_loc);
              
             
           // tempx = tempx  - (Xdistr[(int) locate]) ;
           // tempy = tempy -  (Ydistr[(int) locate]) ;
           tempx = tempx  - (value_X_toAdd) ;
           tempy = tempy -  (value_Y_toAdd) ;
           
             if (tempx*multi_factor>0 && tempx * multi_factor < sizex && tempy * multi_factor>0 && tempy * multi_factor < sizex )
            {
                X[p] = (tempx+44) * multi_factor ;
                Y[p] = (tempy+44) * multi_factor ;
            }
             else{
                  X[p] = INVALID_PIX_VALUE ;
                  Y[p] = INVALID_PIX_VALUE ;
             }
           //LOG(INFO)<<value_X_toAdd<<" "<<Xdistr[(int) locate]<<" "<<value_Y_toAdd<<" "<<Ydistr[(int) locate];
               


            } else {
                X[p] = INVALID_PIX_VALUE;
                Y[p] = INVALID_PIX_VALUE;
            }
        } else {
            X[p] = INVALID_PIX_VALUE;
            Y[p] = INVALID_PIX_VALUE;
        }
    }
    return (EXIT_SUCCESS) ;
}


int uvtRelativeAspectIM::calculateRefFrame (vector<int> &sizevect , vector<float> &xref , vector<float> &yref , vector<float> &intref , int avgfact , vector<float> &xrefcumm , vector<float> &yrefcumm , vector<float> &intrefcumm)
{
    double temp_x1 , temp_y1 , temp_int1 , temp_x2 , temp_y2 , temp_int2 ;
    float x_fref[xref.size ()] , y_fref[xref.size ()] , xrefarr[xref.size ()] , yrefarr[xref.size ()] , intrefarr[xref.size ()] ;
//%#Added ON 20July-divide by zero#%  
if(avgfact==0){
        LOG(ERROR)<<"***Divide by zero***";
        return(EXIT_FAILURE);
    }
//%#-Till this-20July17#%
    for (int d1 = 0 ; d1 < sizevect[0] ; d1 ++)
    {

        int total = sizevect[0] ;
        int cnt_ref = 0 ;
        temp_x1 = xref[d1] ;
        temp_y1 = yref[d1] ;
        temp_int1 = intref[d1] ;
        // LOG(INFO)<<size_rows[0]<<" "<<size_rows[1]<<" "<<size_rows[2]<<" "<<endl;
        for (int d2 = 1 ; d2 < avgfact ; d2 ++)
        {
            total = total + sizevect[d2] ;
            for (int d3 = total - sizevect[d2] ; d3 < total ; d3 ++)
            {
                double diff_x = 0 , diff_y = 0 ;
                temp_x2 = xref[d3] ;
                temp_y2 = yref[d3] ;
                temp_int2 = intref[d3] ;
                //LOG(INFO)<<d3 << "  " <<temp_x2<<"  "<<temp_y2<<endl;
                diff_x = temp_x1 - temp_x2 ;
                //LOG(INFO)<<"the "<<diff_x<<endl;
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
        //  LOG(INFO)<<"outside of the Calcu"<<endl;
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


int uvtRelativeAspectIM::takeDarkinfo ()
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
// int uvtImRa_commonArray:: performQEMCPcorrection(float *sigdata,int xsize,int ysize,double fact){      
//     
//     for(int i=0;i<xsize*ysize;i++){
//         if(sigdata[i]!=INVALID_PIX_VALUE){
//         sigdata[i]=sigdata[i]*fact;
//         }
//     }
//     
//     
//     return(EXIT_SUCCESS);
// }
// 
// int uvtImRa_commonArray::performFlatFieldCorr(float *frmsigdata,float *flatfieldarry,int sizex,int sizey)
//{
//for (int pixno = 0 ; pixno <sizex*sizey  ; pixno++)
//     { 
//         if(frmsigdata[pixno]!=INVALID_PIX_VALUE){
//            frmsigdata[pixno] = frmsigdata[pixno] * flatfieldarry[pixno] ;
//         }
//     }
////        for (int pix = 0 ; pix < sizex * sizey; pix++)
////        {
////            frmexpdata[pix] = 0.0f ;
////            frmexpdata[pix] = badpixarry[pix] * intgrntime ;
////        }    
//    
//    
//    
//    return(EXIT_SUCCESS);
//}

//int uvtImRa_commonArray::readQEMCPFile ()
//{
//    LOG(INFO) << "Reading QE MCP  Temperature vs filter file from calDB........" ;
//    fitsfile *fqemcp ;
//    int status = 0 ;
//    fits_open_file (&fqemcp , qeFile , READONLY , &status) ;
//    printError (status , "Error in opening the qeFile ",qeFile) ;
//    int hdutype ;
//    fits_movabs_hdu (fqemcp , 2 , &hdutype , &status) ;
//    printError (status , "Error in moving to 2nd HDU in qeFile  ",qeFile) ;
//    if (hdutype != BINARY_TBL)
//    {
//        LOG(ERROR) << endl << "***Expected binary table at hdu 2 of temperature vs filter file*** " << endl ;
//        return (EXIT_FAILURE) ;
//    }
//    long nrows ;
//    fits_get_num_rows (fqemcp , &nrows , &status) ;
//    printError (status , "Error in readQEMCP()",qeFile) ;
//    nCalDBTempValues = nrows ;
//    temp = new float[nrows] ;
//    f0 = new float[nrows] ;
//    f1 = new float[nrows] ;
//    f2 = new float[nrows] ;
//    f3 = new float[nrows] ;
//    f5 = new float[nrows] ;
//    f6 = new float[nrows] ;
//    f7 = new float[nrows] ;
//    f4 = new float[nrows] ;
//    
//    fits_read_col (fqemcp , TFLOAT , 1 , 1 , 1 , nrows , NULL , (void*) temp , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb",qeFile) ;
//    fits_read_col (fqemcp , TFLOAT , 2 , 1 , 1 , nrows , NULL , (void*) f0 , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb",qeFile) ;
//   
//    fits_read_col (fqemcp , TFLOAT , 3 , 1 , 1 , nrows , NULL , (void*) f1 , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb",qeFile) ;
//    fits_read_col (fqemcp , TFLOAT , 4 , 1 , 1 , nrows , NULL , (void*) f2 , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb",qeFile) ;
//    fits_read_col (fqemcp , TFLOAT , 5 , 1 , 1 , nrows , NULL , (void*) f3 , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb",qeFile) ;
//    fits_read_col (fqemcp , TFLOAT , 6 , 1 , 1 , nrows , NULL , (void*) f4 , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb",qeFile) ;
//    fits_read_col (fqemcp , TFLOAT , 7 , 1 , 1 , nrows , NULL , (void*) f5 , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb",qeFile) ;
//    fits_read_col (fqemcp , TFLOAT , 8 , 1 , 1 , nrows , NULL , (void*) f6 , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb",qeFile) ;
//    if(strcmp(datainfo.getDetector (),"FUV")==0 || strcmp(datainfo.getDetector (),"VIS")==0)
//    {      fits_read_col (fqemcp , TFLOAT , 9 , 1 , 1 , nrows , NULL , (void*) f7 , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb",qeFile) ;
//    }
//    LOG(INFO) <<"Reading QEMCP Temperature vs filter file from caldb Finished.." ;
//    fits_close_file (fqemcp , &status) ;
//    printError (status , "Error in closing qeFile",qeFile) ;
//    return (EXIT_SUCCESS) ;
//} 


//int uvtImRa_commonArray::getTemp ()
//{
//    int status = 0 ;
//    //reading observation id from data
//   
//    fitsfile *flbt ;
//    fits_open_file (&flbt , lbtfile , READONLY , &status) ;
//    printError (status , "Error opening LBT file",lbtfile) ;
//    
//    char obsid[FLEN_KEYWORD];
//    
//    //reading observation id from header of LBT file
//    fits_read_key(flbt, TSTRING, "OBS_ID",obsid, NULL, &status);
//    printError(status, " Error reading OBS_ID from header " ,lbtfile);
//    
//    fits_movabs_hdu (flbt , 2 , NULL , &status) ;
//    printError (status , "Error moving to HDU 2 of lbtfile",lbtfile) ;
//    //check when  number of Rows in lbt file > 1(TBD)
//    fits_get_num_rows (flbt , &nrows_lbt , &status) ;
//    printError (status , "Error reading the number of rows in lbt file",lbtfile) ;
//     int colinside = 0 , coloutside = 0 ;             //variables to store column numbers to be used from LBT file
//    char *md = datainfo.getDetector();      //read channel for data
//    
//    LOG(INFO)<<"Channel :" <<md<<endl;
//          
//    if (md = (char *) "NUV")
//    {
//        colinside = INSIDE_TEMP_NUV ;
//        coloutside = OUTSIDE_TEMP_NUV ;
//    }
//    else if (md = (char *) "FUV")
//    {
//        colinside =INSIDE_TEMP_FUV ;
//        coloutside = OUTSIDE_TEMP_FUV ;
//    }
//    else if (md = (char *) "VIS")
//    {
//        colinside = INSIDE_TEMP_VIS;
//        coloutside = OUTSIDE_TEMP_VIS;
//    }
//    time_lbt=new double[nrows_lbt];
//   insideTemp = new float [nrows_lbt];
//    outsideTemp= new float[nrows_lbt];
////    int rowno;              
//    fits_read_col (flbt , TDOUBLE , 1 , 1 , 1 ,nrows_lbt , NULL , time_lbt , NULL , &status) ;
//    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
//    fits_read_col (flbt , TFLOAT , colinside , 1 , 1 ,nrows_lbt , NULL , insideTemp , NULL , &status) ;
//    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
//    fits_read_col (flbt , TFLOAT , coloutside , 1 , 1 ,nrows_lbt , NULL , outsideTemp , NULL , &status) ;
//    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
//                   //row number in LBT file from which the temperatures will be read
//     fits_close_file (flbt , &status) ;
//   printError (status , "Error in closing file",lbtfile) ;
//   
//  
//    return (EXIT_SUCCESS) ;
//}



int uvtRelativeAspectIM::getHistory (vector<string> &vhistory)
{
   // char *user = getlogin () ;
    int cnt = 0 ;
    char validgtiflag_str[FLEN_FILENAME] ;
   // string str = "Module run by " + (string) user ;
    char temp[PIL_LINESIZE];
   // vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " input Level1 tar  file = " + (string) level1indir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " caldb used= " + (string) caldbindir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Module Output directory = " + (string)level2outdir ) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " channel used= " + (string)channel) ;
   vhistory.push_back ((string) getSerialNo (cnt) + " CRC flag= " + (string)convertIntToStr(crc_flag)) ;
    //vhistory.push_back ((string) getSerialNo (cnt) + "Total Junk Frames with striping= " + (string)convertIntToStr(counter_JunkFrame)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Flag for CRC failure action = " + (string) convertIntToStr (dropframe)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Dark subtraction Flag = " + (string)convertIntToStr (darkframe_flag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " GTI flag  = " + (string)convertIntToStr(gti_flag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " VIS stripe handler flag  = " + (string)convertIntToStr(junkFrameFlag)) ;
    if(junkFrameFlag==1){
        vhistory.push_back ((string) getSerialNo (cnt) + " Threshold for stripe detection  = " + (string)convertIntToStr(thrJunkFrame)) ;
    }
    vhistory.push_back ((string) getSerialNo (cnt) + " Flat field Flag = " + (string)convertIntToStr( flatfieldFlag)) ;
    //else if (gti_flag == 1)
    //{
     vhistory.push_back ((string) getSerialNo (cnt) + " Final Padded Size  = " + (string) convertIntToStr(padding_dim) );
        //sprintf (validgtiflag_str , "%d" , valid_gtiflag) ;
     vhistory.push_back ((string) getSerialNo (cnt) + " No. of raw frames for accumulation  = " + (string)convertIntToStr( no_ofFramesToAcc)) ;
        //if (all_Or_custom)
        //{
     vhistory.push_back ((string) getSerialNo (cnt) + " QEMCP flag= " + (string) convertIntToStr(qe_mcpFlag)) ;
        //}
       // else if (all_Or_custom == 0)
        //{
     vhistory.push_back ((string) getSerialNo (cnt) + " Sub-division Flag = " + (string) convertIntToStr (subdivisionFlag)) ;
        //}

    //}
            if(subdivisionFlag)
            vhistory.push_back ((string) getSerialNo (cnt) + " Sub-division size = " + (string) convertIntToStr(subDivision_size)) ;
            
            vhistory.push_back ((string) getSerialNo (cnt) + " Star finding algorithm = " + (string) convertIntToStr(star_detect_algo_flag)) ;
            
            if(star_detect_algo_flag==1 || star_detect_algo_flag==3 || star_detect_algo_flag==4){
                vhistory.push_back ((string) getSerialNo (cnt) + " First cut threshold for star detection = " + (string) convertFloatToStr(rms_mul_factor_default)) ;
                  vhistory.push_back ((string) getSerialNo (cnt) + " Minimum targeted stars = " + (string)convertIntToStr(min_num_stars)) ;
            }
            
             vhistory.push_back ((string) getSerialNo (cnt) + " Neighbourhood criteron for identifying stars. = " + (string)convertIntToStr( refine_Winsize)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Box Size to compute Centroid for  detected stars = " + (string)convertIntToStr( centroid_Winsize)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Threshold for CR affected pixels in IM mode = " + (string)convertFloatToStr( threshold)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " No. of accumulated frames to be discarded at beginning = " + (string)convertIntToStr(frames_toDiscard)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " No. of accumulated frames used to construct reference frame = " + (string)convertIntToStr( nFrameToAverage)) ;
             //vhistory.push_back ((string) getSerialNo (cnt) + " Distance from Which matching of stars to be started i.e Diff_dist " + (string)convertFloatToStr( nbhd_dist)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Frequency domain filtering flag  " + (string) convertIntToStr(freqDomainFilter_Flag)) ;
		vhistory.push_back ((string) getSerialNo (cnt) + " frequency value  " + (string) convertFloatToStr(freqvalue)) ;
		//vhistory.push_back ((string) getSerialNo (cnt) + " Fitting Flag  " + ) ;
             if(freqDomainFilter_Flag==1)
             vhistory.push_back ((string) getSerialNo (cnt) + " Type Filtering used   " + (string) convertIntToStr(type_Filtering)) ;
             
             if(type_Filtering==1 || type_Filtering==2 ){
                  vhistory.push_back ((string) getSerialNo (cnt) + " order of Pitch   " + (string) convertIntToStr(orderpitch)) ;
                  vhistory.push_back ((string) getSerialNo (cnt) + " order of roll   " + (string) convertIntToStr(orderroll)) ;
                  vhistory.push_back ((string) getSerialNo (cnt) + " order of yaw  " + (string) convertIntToStr(orderyaw)) ;
                   vhistory.push_back ((string) getSerialNo (cnt) + " Delta time    " + (string) convertFloatToStr(delta_time)) ;

                 
             }
 vhistory.push_back ((string) getSerialNo (cnt) + " Shift And Rotation algorithm used =  " + (string) convertIntToStr(shift_rotation_algo)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Theta switch(ON/OFF)   " + (string) convertIntToStr(flag_thetaComp)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Minimum search distance for star match in successive frames  " + (string) convertIntToStr(diff_dist)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Flag for writing matched star-pair list" + (string) convertIntToStr(match_stars_file_flag)) ;
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
            // vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for Shift and Rotate= " + (string)convertIntToStr( wtd_snr)) ;
             //vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for Weighted mean= " + (string)convertIntToStr( wtd_wm)) ;
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
