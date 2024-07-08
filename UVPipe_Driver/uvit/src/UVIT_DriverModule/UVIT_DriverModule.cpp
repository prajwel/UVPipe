
#include "UVIT_DriverModule.h"
#include<pil.h>
#include<glog/logging.h>
#include<uvtUtils.h>
#include <sstream>
#define EXE_RA_IM "uvtRelativeAspectIM"
#define EXE_L2_PC "uvtLevel2PC"
#define EXE_RA_PC "uvtRelativeAspectPC"

int UVIT_DriverModule::readPILParameters(int argc, char** argv)
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

    if (PIL_OK != (status = PILGetBool ("NUVonNUVflag" , &NUVonNUVflag)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }
    if (NUVonNUVflag==0){
       if (PIL_OK != (status = PILGetBool ("FUVonFUVflag" , &FUVonFUVflag)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }  
        
    }
    
    if (PIL_OK != (status = PILGetBool ("ManualMode" , &manualMode)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }
    if(manualMode==TRUE){
        if (PIL_OK != (status = PILGetFname ("previousOutputL2" , temp)))
    {
        LOG (INFO) << endl << "***Error reading input level1 directory name***" ;
        return status ;
    }
        prev_Output_DirL2.assign(temp);
    
    }
    if (PIL_OK != (status = PILGetString ("channel" , temp)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    channel.assign (temp) ;

    //dataingest param
    if (PIL_OK != (status = PILGetBool ("dropframe" , &dropframe)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetBool ("darkframeFlag" , &darkframe_flag)))
    {
        LOG (INFO) << endl << "***Error Reading darkframeFlag:" << darkframe_flag << "***" ;
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
    if (PIL_OK != (status = PILGetReal4 ("diffDist" ,  &nbhd_dist)))
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
            if (PIL_OK != (status = PILGetReal ("freq_value" , (double*) &freqvalue)))
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
    LOG(INFO)<<"=======================================================";
    LOG(INFO)<<"============NOW ASKING LEVEL2 PC PARAMETERS============";
    LOG(INFO)<<"========================================================";
     if (PIL_OK != (status = PILGetBool ("dropframepc" , &dropframepc)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("parityFlagpc" , &parity_flagpc)))
    {
        LOG (INFO) << "Error reading parity Flag" << endl ;
        return (status) ;
    }
if (PIL_OK != (status = PILGetBool ("crcflagpc" , &crc_flagpc)))
    {
        LOG (INFO) << "Error reading parity Flag" << endl ;
        return (status) ;
    }
    
//    if (PIL_OK != (status = PILGetInt ("att_flagVal" , &att_flag_val)))
//    {
//        LOG (INFO) << endl << "***Error Reading att  flag value :***" ;
//        return status ;
//    }
    
 if (PIL_OK != (status = PILGetReal("startTimepc" , &star_time_ForDatasepc)))
    {
        LOG (INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    
     if (PIL_OK != (status = PILGetReal ("endTimepc" , &end_time_ForDatasepc)))
    {
        LOG (INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thr_One_crpc" , &first_thr_Crpc)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thr_Two_crpc" , &sec_thr_Crpc)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    

    if (PIL_OK != (status = PILGetBool ("GTI_FLAGpc" , &gti_flagpc)))
    {
        LOG (INFO) << endl << "***Error Reading GTI flag :" << gti_flag << "***" ;
        return status ;
    }
    if (gti_flagpc)
    {
        all_or_custom = 1 ;
        valid_bit = 1 ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thresholdMultphpc" , &thr_multiphpc)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    

//    if (PIL_OK != (status = PILGetInt ("algoFlagpc" , &star_detect_algo_flagpc)))
//    {
//        LOG (INFO) << endl << "***Error reading algo_flag ***" ;
//        return status ;
//    }

    if (PIL_OK != (status = PILGetReal4 ("thresholdpc" , &sd_multi_factor_defaultpc)))
    {
        LOG (INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetInt ("minimum_targetedstarspc" , &minimum_No_of_Starspc)))
    {
        LOG (INFO) << endl << "***Error reading algo_flag ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("refinedWinSizepc" , &refine_Winsizepc)))
    {
        LOG (INFO) << endl << "***Error reading refine window size ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("centroidWinSizepc" , &centroid_Winsizepc)))
    {
        LOG (INFO) << endl << "***Error reading centroid  window size ***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetReal4 ("diffDistpc" , (float*) &diff_Distpc)))
    {
        LOG (INFO) << endl << "***Error reading diffDist parameter***" ;
        return status ;
    }
 
   if (PIL_OK != (status = PILGetFname ("catalogpathpc" , catalogpathpc)))
            {
                LOG(ERROR) << endl << "\033[1;31m***Error reading output directory name***" ;
                return status ;
            }
     if (PIL_OK != (status = PILGetString ("att_timecolpc" , att_timecolpc)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("database_namepc" , databasenamepc)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("att_qcolpc" , att_qcolpc)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("shiftRotDetAlgoFlagpc" , &shift_N_Rotate_algopc)))
    {
        LOG (INFO) << endl << "***Error reading shiftRotDetAlgoFlag***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("frameIntFlagpc" , &fi_flagpc)))
    {
        LOG (INFO) << endl << "***Error reading frameIntegration flag***" ;
        return status ;
    }
   
            if (PIL_OK != (status = PILGetInt ("framesDiscardpc" , &nFrameDiscard_fipc)))
            {
                LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
                return status ;
            }
            if (PIL_OK != (status = PILGetInt ("framesComputepc" , &nFrameIntegrate_fipc)))
            {
                LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
                return status ;
            }
    //}
    if (PIL_OK != (status = PILGetInt ("FrameIntDimpc" , &IMG_DIM_FIpc)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
      if (PIL_OK != (status = PILGetInt ("RegAvgfrmsizepc" , &FINALFRAMESIZE_REGAVGpc)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
//Added
     if (PIL_OK != (status = PILGetInt ("QEMCP_tobedonepc" , &qemcpflagpc)))
    {
        LOG (INFO) << endl << "***Error reading QEMCP flag ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("CentCorr_tobedonepc" , &centCorrflagpc)))
    {
        LOG (INFO) << endl << "***Error reading Centroid correction flag ***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetInt ("CentBias_tobedonepc" , &centBiasflagpc)))
    {
        LOG (INFO) << endl << "***Error reading centroid bias flag***" ;
        return status ;
    }
   if (PIL_OK != (status = PILGetInt ("DetectDist_tobedonepc" , &DetectDistflagpc)))
    {
        LOG (INFO) << endl << "***Error reading detector distortion flag ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("OpticDist_tobedonepc" , &OpticDistflagpc)))
    {
        LOG (INFO) << endl << "***Error reading Optical distortion flag ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("pathToOutputTarpc" , outtarpathpc)))
    {
        LOG (INFO) << endl << "***Error reading output tar path***" ;
        return status ;
    }
    if(!DirExists(outtarpathpc)){
        LOG(ERROR)<<"**Tar file path not exists!!!***";
        return(EXIT_FAILURE);
    }

//till this    
if(FINALFRAMESIZE_REGAVGpc>IMG_DIM_FIpc)
   {
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
      if (PIL_OK != (status = PILGetInt ("search_algo_forFullFrameAstpc" , &search_algo_ctlgpc)))
        {
            LOG (INFO) << endl << "***Error Reading search method" << search_algo_ctlgpc << "***" ;
            return status ;
        }
 if (PIL_OK != (status = PILGetString("Radi_searchpc" , (char*)&rad_searchpc)))
    {
        LOG (INFO) << endl << "***Error reading the radius value***" ;
        return status ;
    }
    strcpy(len_apc,"NULL");
    strcpy(len_bpc,"NULL");
    if(search_algo_ctlgpc==1 || search_algo_ctlgpc==3 || search_algo_ctlgpc==5)
    {
      if (PIL_OK != (status = PILGetString ("len_rect_apc" , len_apc)))
        {
            LOG (INFO) << endl << "***Error Reading length of rectangle :" <<len_apc << "***" ;
            return status ;
        }
    
      if (PIL_OK != (status = PILGetString ("len_rect_bpc" , len_bpc)))
        {
            LOG (INFO) << endl << "***Error Reading width of rectangle:" << len_bpc << "***" ;
            return status ;
        }
    
    } 
    
     
     if (PIL_OK != (status = PILGetBool ("flag_thetaComppc" , &flag_thetaComppc)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobberpc << "***" ;
        return status ;
    }
    
    

    if (PIL_OK != (status = PILGetBool ("clobberpc" , &clobberpc)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobberpc << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("historypc" , &historypc)))
    {
        LOG (INFO) << endl << "***Error Reading history parameter:" << history << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("modepc" , modepc)))
    {
        LOG (INFO) << "***Error Reading mode parameter:" << historypc << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskucpc" , &wtd_ucpc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing to Disk Flag For unit Conversion:" << "***" ;
        return status ;
    }
    //    if (PIL_OK != (status = PILGetInt ("Write_todiskuc" , &wtd_sd)))
    //    {
    //        LOG(ERROR) << endl << "***Error Reading Writing to Disk Flag For unit Conversion:" << "***" ;
    //        return status ;
    //    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskbppc" , &wtd_bppc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the bad pixel module:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskffpc" , &wtd_ffpc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the FlatFieldCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskpppc" , &wtd_pppc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentCorrpc" , &wtd_centCorrpc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentBiaspc" , &wtd_centBiaspc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetInt ("Write_todisksdpc" , &wtd_sdpc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the subDivision:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskcrpc" , &wtd_crpc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskqepc" , &wtd_qemcppc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskfipc" , &wtd_fipc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the Accumulated module:" << "***" ;
        return status ;
    }
//    if (PIL_OK != (status = PILGetInt ("Write_todisksc" , &wtd_fsc)))
//    {
//        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the find Star Centroids:" << "***" ;
//        return status ;
//    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskddpc" , &wtd_ddpc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the Detect Distortion:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskodpc" , &wtd_odpc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the optical Distortion:" << "***" ;
        return status ;
}
    
    //for FUV
    
    LOG(INFO)<<"=======================================================";
    LOG(INFO)<<"============NOW ASKING LEVEL2 PC PARAMETERS for FUV channel============";
    LOG(INFO)<<"========================================================";
     if (PIL_OK != (status = PILGetBool ("dropframepcfuv" , &dropframepcfuv)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("parityFlagpcfuv" , &parity_flagpcfuv)))
    {
        LOG (INFO) << "Error reading parity Flag" << endl ;
        return (status) ;
    }
if (PIL_OK != (status = PILGetBool ("crcflagpcfuv" , &crc_flagpcfuv)))
    {
        LOG (INFO) << "Error reading parity Flag" << endl ;
        return (status) ;
    }
    
//    if (PIL_OK != (status = PILGetInt ("att_flagVal" , &att_flag_val)))
//    {
//        LOG (INFO) << endl << "***Error Reading att  flag value :***" ;
//        return status ;
//    }
    
 if (PIL_OK != (status = PILGetReal("startTimepcfuv" , &star_time_ForDatasepcfuv)))
    {
        LOG (INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    
     if (PIL_OK != (status = PILGetReal ("endTimepcfuv" , &end_time_ForDatasepcfuv)))
    {
        LOG (INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thr_One_crpcfuv", &first_thr_For_Crfuv)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thr_Two_crpcfuv" ,&sec_thr_For_Crfuv)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    

    if (PIL_OK != (status = PILGetBool ("GTI_FLAGpcfuv" , &gti_flagpcfuv)))
    {
        LOG (INFO) << endl << "***Error Reading GTI flag :" << gti_flag << "***" ;
        return status ;
    }
    if (gti_flagpc)
    {
        all_or_custom = 1 ;
        valid_bit = 1 ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thresholdMultphpcfuv" , &thr_multiphpcfuv)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    

//    if (PIL_OK != (status = PILGetInt ("algoFlagpcfuv" , &star_detect_algo_flagpcfuv)))
//    {
//        LOG (INFO) << endl << "***Error reading algo_flag ***" ;
//        return status ;
//    }

    if (PIL_OK != (status = PILGetReal4 ("thresholdpcfuv" , &sd_multi_factor_defaultpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetInt ("minimum_targetedstarspcfuv" , &minimum_No_of_Starspcfuv)))
    {
        LOG (INFO) << endl << "***Error reading algo_flag ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("refinedWinSizepcfuv" , &refine_Winsizepcfuv)))
    {
        LOG (INFO) << endl << "***Error reading refine window size ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("centroidWinSizepcfuv" , &centroid_Winsizepcfuv)))
    {
        LOG (INFO) << endl << "***Error reading centroid  window size ***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetReal4 ("diffDistpcfuv" , (float*) &diff_Distpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading diffDist parameter***" ;
        return status ;
    }
 
   if (PIL_OK != (status = PILGetFname ("catalogpathpcfuv" , catalogpathpcfuv)))
            {
                LOG(ERROR) << endl << "\033[1;31m***Error reading output directory name***" ;
                return status ;
            }
     if (PIL_OK != (status = PILGetString ("att_timecolpcfuv" , att_timecolpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("database_namepcfuv" , databasenamepcfuv)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("att_qcolpcfuv" , att_qcolpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("shiftRotDetAlgoFlagpcfuv" , &shift_N_Rotate_algopcfuv)))
    {
        LOG (INFO) << endl << "***Error reading shiftRotDetAlgoFlag***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("frameIntFlagpcfuv" , &fi_flagpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading frameIntegration flag***" ;
        return status ;
    }
   
            if (PIL_OK != (status = PILGetInt ("framesDiscardpcfuv" , &nFrameDiscard_fipcfuv)))
            {
                LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
                return status ;
            }
            if (PIL_OK != (status = PILGetInt ("framesComputepcfuv" , &nFrameIntegrate_fipcfuv)))
            {
                LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
                return status ;
            }
    //}
    if (PIL_OK != (status = PILGetInt ("FrameIntDimpcfuv" , &IMG_DIM_FIpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
      if (PIL_OK != (status = PILGetInt ("RegAvgfrmsizepcfuv" , &FINALFRAMESIZE_REGAVGpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
//Added
     if (PIL_OK != (status = PILGetInt ("QEMCP_tobedonepcfuv" , &qemcpflagpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("CentCorr_tobedonepcfuv" , &centCorrflagpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetInt ("CentBias_tobedonepcfuv" , &centBiasflagpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
   if (PIL_OK != (status = PILGetInt ("DetectDist_tobedonepcfuv" , &DetectDistflagpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("OpticDist_tobedonepcfuv" , &OpticDistflagpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("pathToOutputTarpcfuv" , outtarpathpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading output tar path***" ;
        return status ;
    }
    if(!DirExists(outtarpathpcfuv)){
        LOG(ERROR)<<"**Tar file path not exists!!!***";
        return(EXIT_FAILURE);
    }

//till this    
if(FINALFRAMESIZE_REGAVGpcfuv>IMG_DIM_FIpcfuv){
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
      if (PIL_OK != (status = PILGetInt ("search_algo_forFullFrameAstpcfuv" , &search_algo_ctlgpcfuv)))
        {
            LOG (INFO) << endl << "***Error Reading search method" << search_algo_ctlgpcfuv << "***" ;
            return status ;
        }
 if (PIL_OK != (status = PILGetString("Radi_searchpcfuv" , (char*)&rad_searchpcfuv)))
    {
        LOG (INFO) << endl << "***Error reading the radius value***" ;
        return status ;
    }
    strcpy(len_apc,"NULL");
    strcpy(len_bpc,"NULL");
    if(search_algo_ctlgpcfuv==1 || search_algo_ctlgpcfuv==3 || search_algo_ctlgpcfuv==5)
    {
      if (PIL_OK != (status = PILGetString ("len_rect_apcfuv" , len_apcfuv)))
        {
            LOG (INFO) << endl << "***Error Reading length of rectangle :" <<len_apc << "***" ;
            return status ;
        }
    
      if (PIL_OK != (status = PILGetString ("len_rect_bpcfuv" , len_bpcfuv)))
        {
            LOG (INFO) << endl << "***Error Reading width of rectangle:" << len_bpc << "***" ;
            return status ;
        }
    
    } 
    
     
     if (PIL_OK != (status = PILGetBool ("flag_thetaComppcfuv" , &flag_thetaComppcfuv)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobberpc << "***" ;
        return status ;
    }
    
    

    if (PIL_OK != (status = PILGetBool ("clobberpcfuv" , &clobberpcfuv)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobberpc << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("historypcfuv" , &historypcfuv)))
    {
        LOG (INFO) << endl << "***Error Reading history parameter:" << history << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("modepcfuv" , modepcfuv)))
    {
        LOG (INFO) << "***Error Reading mode parameter:" << historypcfuv << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskucpcfuv" , &wtd_ucpcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing to Disk Flag For unit Conversion:" << "***" ;
        return status ;
    }
    //    if (PIL_OK != (status = PILGetInt ("Write_todiskuc" , &wtd_sd)))
    //    {
    //        LOG(ERROR) << endl << "***Error Reading Writing to Disk Flag For unit Conversion:" << "***" ;
    //        return status ;
    //    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskbppcfuv" , &wtd_bppcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the bad pixel module:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskffpcfuv" , &wtd_ffpcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the FlatFieldCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskpppcfuv" , &wtd_pppcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentCorrpcfuv" , &wtd_centCorrpcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentBiaspcfuv" , &wtd_centBiaspcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetInt ("Write_todisksdpcfuv" , &wtd_sdpcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the subDivision:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskcrpcfuv" , &wtd_crpcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskqepcfuv" , &wtd_qemcppcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskfipcfuv" , &wtd_fipcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the Accumulated module:" << "***" ;
        return status ;
    }
//    if (PIL_OK != (status = PILGetInt ("Write_todisksc" , &wtd_fsc)))
//    {
//        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the find Star Centroids:" << "***" ;
//        return status ;
//    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskddpcfuv" , &wtd_ddpcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the Detect Distortion:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskodpcfuv" , &wtd_odpcfuv)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the optical Distortion:" << "***" ;
        return status ;
}
    
    LOG(INFO)<<"RELATIVE ASPECT PC parameter reading started";
    
//    if (PIL_OK != (status = PILGetFname ("level1indirrapc" , temp)))
//    {
//        LOG (INFO) << endl << "***Error reading input level1 directory name***" ;
//        return status ;
//    }
//    level1indirrapc.assign (temp) ;
//    if (PIL_OK != (status = PILGetFname ("caldbdir" , temp)))
//    {
//        LOG (INFO) << endl << "***Error reading caldb directory name***" ;
//        return status ;
//    }
//    caldbindir.assign (temp) ;
//    if (PIL_OK != (status = PILGetFname ("level2outdirrapc" , temp)))
//    {
//        LOG (INFO) << endl << "***Error reading output level2 directory name***" ;
//        return status ;
//    }
//    level2outdirrapc.assign (temp) ;
    //}
//    if (PIL_OK != (status = PILGetString ("channelrapc" , temp)))
//    {
//        LOG (INFO) << endl << "***Error reading channel***" ;
//        return status ;
//    }
//    channelrapc.assign (temp) ;

    //dataingest param
    if (PIL_OK != (status = PILGetBool ("dropframerapc" , &dropframerapc)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("parityFlagrapc" , &parity_flagrapc)))
    {
        LOG (INFO) << "Error reading parity Flag" << endl ;
        return (status) ;
    }
if (PIL_OK != (status = PILGetBool ("crcflagrapc" , &crc_flagrapc)))
    {
        LOG (INFO) << "Error reading crcFlag" << endl ;
        return (status) ;
    }
//    if (PIL_OK != (status = PILGetInt ("att_flagVal" , &att_flag_val)))
//    {
//        LOG (INFO) << endl << "***Error Reading att  flag value :" << padding_dim << "***" ;
//        return status ;
//    }

//    if (PIL_OK != (status = PILGetReal4 ("thrValrapc" , &Threshold_crrapc)))
//    {
//        LOG (INFO) << endl << status << "***Error reading threshold***" ;
//        return status ;
//    }
    if (PIL_OK != (status = PILGetInt ("cmpFramesrapc" , &nCompareFramesrapc)))
    {
        LOG (INFO) << endl << status << "***Error reading number of frames to be compared***" ;
        return status ;
    }


    if (PIL_OK != (status = PILGetBool ("GTI_FLAGrapc" , &gti_flagrapc)))
    {
        LOG (INFO) << endl << "***Error Reading GTI FLAG value:***" ;
        return status ;
    }
    if (gti_flag)
    {
        all_or_custom = 1 ;
        valid_bit = 1 ;
    }
    if (PIL_OK != (status = PILGetReal4 ("thresholdMultphrapc" , &thr_multiphrapc)))
    {
        LOG (INFO) << endl << status << "***Error reading threshold***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("flatfieldFlagrapc" , &flatfieldFlagrapc)))
    {
        LOG (INFO) << endl << "***Error Reading flatfieldFlag:" << flatfieldFlag << "***" ;
        return status ;
    }


    if (PIL_OK != (status = PILGetBool ("qemcpFlagrapc" , &qemcpFlagrapc)))
    {
        LOG (INFO) << endl << "***Error Reading qemcpFlag :" << qemcpFlagrapc << "***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetInt ("algoFlagrapc" , &star_detect_algo_flagrapc)))
    {
        LOG (INFO) << endl << "***Error reading algo_flag ***" ;
        return status ;
    }
    if (star_detect_algo_flagrapc == 1 || star_detect_algo_flagrapc==3 || star_detect_algo_flagrapc==4)
    {
        if (PIL_OK != (status = PILGetReal4 ("thresholdrapc" , &rms_mul_factor_defaultrapc)))
        {
            LOG (INFO) << endl << "***Error reading threshold value for star detection***" ;
            return status ;
        }

        if (PIL_OK != (status = PILGetInt ("minimumTargetedstarsrapc" , &minimum_No_of_Starsrapc)))
        {
            LOG (INFO) << endl << "***Error in reading the minimum number of stars ***" ;
            return status ;
        }
    }
    else if (star_detect_algo_flag == 2)
    {

        if (PIL_OK != (status = PILGetInt ("StarDetectionCentroidSqrSizerapc" , &footprint_win_star_detectrapc)))
        {
            LOG (INFO) << endl << "***Error reading  window size for JOE's algorithm ***" ;
            return status ;
        }
        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroidPriThrrapc" , &primary_threshold_Valrapc)))
        {
            LOG (INFO) << endl << "***Error reading primary threshold value for JOE's algorithm***" ;
            return status ;
        }
        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroidSecThrrapc" , &secondary_threshold_Valrapc)))
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

    if (PIL_OK != (status = PILGetInt ("refinedWinSizerapc" , &refine_Winsizerapc)))
    {
        LOG (INFO) << endl << "***Error reading refine window size ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("centroidWinSizerapc" , &centroid_Winsizerapc)))
    {
        LOG (INFO) << endl << "***Error reading centroid  window size ***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetInt ("FrameIntDimrapc" , &IMG_DIM_FIrapc)))
    {
        LOG (INFO) << endl << "***Error reading frame size of frame integration ***" ;
        return status ;
    }    
   
    if(star_detect_algo_flag==1)
    {
     if (PIL_OK != (status = PILGetReal4 ("background_factrapc" , &backgrd_factrapc)))
        {
            LOG (INFO) << endl << "***Error reading background_fact ***" ;
            return status ;
        }
    }
    else{
        backgrd_fact=0;//can be ignored.
    }
    if (PIL_OK != (status = PILGetInt ("Window_nbrapc" , &win_sizerapc)))
    {
        LOG (INFO) << endl << "***Error reading centroid  window size ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("GenMatchStarsFile_flagrapc" , &match_stars_file_flagrapc)))
    {
        LOG (INFO) << endl << "***Error reading GenMatchStarsFile_flag parameter***" ;
        return status ;
    }


    if (PIL_OK != (status = PILGetInt ("framesToBeDiscardrapc" , &frames_toDiscardrapc)))
    {
        LOG (INFO) << endl << "***Error reading the number of frames to discard in reference frame calculation***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("averageFactorrapc" , &nFrameToAveragerapc)))
    {
        LOG (INFO) << endl << "***Error reading average factor value in reference frame calculation ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetReal4 ("diffDistrapc" , (float*) &pair_nbhd_distancerapc)))
    {
        LOG (INFO) << endl << "***Error reading diffDist parameter***" ;
        return status ;
    }
//    if (PIL_OK != (status = PILGetReal4 ("error_per" , &err_per)))
//    {
//        LOG (INFO) << endl << "***Error reading error percent  ***" ;
//        return status ;
//    }
    if (PIL_OK != (status = PILGetInt ("freqDomainFilterFlagrapc" , &FreqDomainFilterFlagrapc)))
    {
        LOG (INFO) << endl << "***Error reading freqDomainFilter Value***" ;
        return status ;
    }
    if (FreqDomainFilterFlagrapc == 0)
    {
        if (PIL_OK != (status = PILGetInt ("typeFilteringrapc" , &type_Filteringrapc)))
        {
            LOG (INFO) << endl << "***Error reading type Filtering Value***" ;
            return status ;
        }
        if (type_Filteringrapc == 1)
        {

            if (PIL_OK != (status = PILGetBool ("fittingFlagrapc" , &fitting_flagrapc)))
            {
                LOG (INFO) << endl << "***Error reading fitting flag***" ;
                return status ;
            }
            if (fitting_flagrapc)
            {
                if (PIL_OK != (status = PILGetInt ("orderPitchrapc" , &orderpitchrapc)))
                {
                    LOG (INFO) << endl << "***Error reading order of pitch***" ;
                    return status ;
                }
                if (PIL_OK != (status = PILGetInt ("orderYawrapc" , &orderyawrapc)))
                {
                    LOG (INFO) << endl << "***Error reading order of yaw***" ;
                    return status ;
                }
                if (PIL_OK != (status = PILGetInt ("orderRollrapc" , &orderrollrapc)))
                {
                    LOG (INFO) << endl << "***Error reading order of roll***" ;
                    return status ;
                }
                if (PIL_OK != (status = PILGetReal ("deltaTimerapc" , &poly_fit_intervalrapc)))
                {
                    LOG (INFO) << endl << "***Error reading the delta_Time***" ;
                    return status ;
                }

            }
        }
        else if (type_Filtering == 0)
        {
            if (PIL_OK != (status = PILGetReal ("freqValuerapc" , (double*) &freqvaluerapc)))
            {
                LOG (INFO) << endl << "***Error reading cut-off frequency value***" ;
                return status ;
            }

        }
    }
    else     if (FreqDomainFilterFlagrapc == 1)
    {
        if (PIL_OK != (status = PILGetReal ("freqValuerapc" , (double*) &freqvaluerapc)))
        {
            LOG (INFO) << endl << "***Error reading freq_value***" ;
            return status ;
        }
    }
    else if (FreqDomainFilterFlagrapc!=2)
    {
        LOG(ERROR)<<"***Value must be 0 0r 1 for FreqDomainFilterFlag ";
        return(EXIT_FAILURE);
    }

    if (PIL_OK != (status = PILGetInt ("shiftRotDetAlgoFlagrapc" , &shift_rotation_algorapc)))
    {
        LOG (INFO) << endl << "***Error reading shiftRotDetAlgoFlag***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("framesDiscardrapc" , &nFrameDiscard_firapc)))
    {
        LOG (INFO) << endl << "***Error reading number of frames to be discarded in frame integration***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("framesComputerapc" , &nFrameIntegrate_firapc)))
    {
        LOG (INFO) << endl << "***Error reading number of frames to integrate ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("flag_thetaComprapc" , &flag_thetaComprapc)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobberrapc << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("clobberrapc" , &clobberrapc)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobberrapc << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("historyrapc" , &historyrapc)))
    {
        LOG (INFO) << endl << "***Error Reading history parameter:" << historyrapc << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("moderapc" , moderapc)))
    {
        LOG (INFO) << "***Error Reading mode parameter:" << historyrapc << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskucrapc" , &wtd_ucrapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing to Disk Flag For unit Conversion:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskbprapc" , &wtd_bprapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the bad pixel module:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskffrapc" , &wtd_ffrapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the FlatFieldCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskpprapc" , &wtd_pprapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentCorrrapc" , &wtd_centCorrrapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskCentBiasrapc" , &wtd_centBiasrapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the PixPadding:" << "***" ;
        return status ;
    }

    //    if (PIL_OK != (status = PILGetInt ("Write_todisksd" , &wtd_sd)))
    //    {
    //        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the subDivision:" << "***" ;
    //        return status ;
    //    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskcrrapc" , &wtd_crrapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskqerapc" , &wtd_qemcprapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the QEMCP correction" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskfirapc" , &wtd_firapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the frameintegration:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskscrapc" ,  &wtd_fscrapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the find Star centroid:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskddrapc" , &wtd_ddrapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the Detector  Distortion:" << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskodrapc" , &wtd_odrapc)))
    {
        LOG (ERROR) << endl << "***Error Reading Writing the Disk Flag For the optical Distortion:" << "***" ;
        return status ;
    }
    
   
    PILClose (status) ;
    return (EXIT_SUCCESS) ;   
    
}
int UVIT_DriverModule::RunRelativeAspectChain()
{
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    LOG(INFO)<<"\033[1;31m===============UVIT RELATIVE ASPECT CHAIN STARTED===============\033[0m";
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    stringstream cmd ;
    flag_tarExtracted=FALSE;
    cmd<<EXE_RA_IM<<" zipFlag = "<<ReturnYorN(zipFlag)<<" level1indir ="<<(string)level1indir<<" caldbdir ="<<(string)caldbindir<<" level2outdir ="<<(string)level2outdir<<" channel ="<<channel<<
            " dropframe = "<<ReturnYorN(dropframe)<<" utcFlag = "<<ReturnYorN(UTC_flag)<<" crcflag = "<<ReturnYorN(crc_flag)<<" darkframeFlag= "<<ReturnYorN(darkframe_flag)<<" junkfileflag ="<<ReturnYorN(junkFrameFlag)<<" thresholdforjunkFrame="<<thrJunkFrame<<" GTI_Flag ="<<ReturnYorN(gti_flag)<<
            " flatfieldFlag ="<<ReturnYorN(flatfieldFlag)<<" paddingDim="<<padding_dim<<" Nacc="<<no_ofFramesToAcc<<" qemcpFlag="<<ReturnYorN(qe_mcpFlag)<<" subdivisionFlag ="<<ReturnYorN(subdivisionFlag)<<
            " subdivisionsize ="<<subDivision_size<<" algoFlag ="<<star_detect_algo_flag<<" threshold ="<<rms_mul_factor_default<<" minimumTargetedstars ="<<min_num_stars<<
            " refineWindow ="<<refine_Winsize<<" centroidWindow="<<centroid_Winsize<<" background_fact ="<<backgrd_fact<<" Window_nb ="<<search_win_size<<
            " GenmatchstarsFile_flag ="<<ReturnYorN(match_stars_file_flag)<<" ThresholdValue ="<<threshold<<" framesToBeDiscard="<<frames_toDiscard<<" averageFactor ="<<nFrameToAverage<<
            " diffDist ="<<nbhd_dist<<" freqDomainFilterFlag="<<freqDomainFilter_Flag<<" type_Filtering ="<<type_Filtering<<" fitting_flag ="<<ReturnYorN(fitting_flag)<<" orderPitch="<<orderpitch<<
            " orderYaw ="<<orderyaw<<" orderRoll ="<<orderroll<<" deltaTime ="<<delta_time<<" freqValue="<<freqvalue<<" shiftRotDetAlgoFlag="<<shift_rotation_algo<<
            " clobber ="<<ReturnYorN(clobber)<<" flag_thetaComp ="<<ReturnYorN(flag_thetaComp)<<" history ="<<ReturnYorN(history)<<" mode ="<<(string)mode<<" Write_todiskuc="<<wtd_uc<<
            " Write_todiskbp="<<wtd_bp<<" Write_todiskff="<<wtd_ff<<" Write_todisksd="<<wtd_sd<<" Write_todiskpp="<<wtd_pp<<" Write_todiskcr="<<wtd_cr<<" Write_todiskqe="<<wtd_qemcp<<
            " Write_todiskac="<<wtd_ac<<" Write_todisksc="<<wtd_fsc<<" Write_todiskdd="<<wtd_dd<<" Write_todiskod="<<wtd_od;
   
    
    LOG(INFO)<<cmd.str();
    try{
     system(cmd.str().c_str())   ;
    }
    catch(...)
    {   
        LOG(ERROR)<<"***Error in relativeAspect generation";
        return(EXIT_FAILURE);
    }
    flag_tarExtracted=TRUE;
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    LOG(INFO)<<"\033[1;31m===============UVIT RELATIVE ASPECT CHAIN ENDS===============\033[0m";
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    
    return(EXIT_SUCCESS);
}
 int UVIT_DriverModule::RunL2PCChain(string channel_l2,int cnt_num)
{
     LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    LOG(INFO)<<"\033[1;31m===============UVIT LEVEL 2 PC CHAIN FOR "<<channel_l2<<" STARTED===============\033[0m";
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    stringstream cmd ;
    flag_tarExtracted=FALSE;
    level2outdirc=level2outdir+"_"+channel_l2+"_"+convertIntToStr(cnt_num);
 
    if(strcmp(channel_l2.c_str(),"NUV")==0){
        OutputDir_l2pcNUV.push_back(level2outdirc);
    cmd<<EXE_L2_PC<<" zipFlag = n"<<" level1indir ="<<(string)level1indirpc<<" caldbdir ="<<(string)caldbindir<<" level2outdir ="<<(string)level2outdirc<<" channel ="<<channel_l2<<" startTime= "<<-9999
            <<" endTime="<<-9999<<" unTar_flag= y"
            <<" dropframe = "<<ReturnYorN(dropframepc)<<" utcFlag = "<<ReturnYorN(UTC_flag)<<" crcflag = "<<ReturnYorN(crc_flagpc)<<" parityFlag= "<<parity_flagpc<<" thr_One_cr="<<first_thr_Crpc<<" thr_Two_cr="<<sec_thr_Crpc<<" GTI_FLAG ="<<ReturnYorN(gti_flagpc)<<
            " thresholdMultph ="<<thr_multiphpc<<" threshold= "<<sd_multi_factor_defaultpc<<" minimumTargetedstars="<<minimum_No_of_Starspc<<
            " refinedWinSize="<<refine_Winsizepc<<" centroidWinSize="<<centroid_Winsizepc<<" diffDist="<<diff_Distpc<<" catalogpath="<<catalogpathpc<<" att_timecol="<<att_timecolpc<<
            " database_name="<<databasenamepc<<" att_qcol="<<att_qcolpc<<" shiftRotDetAlgoFlag="<<shift_N_Rotate_algopc<<" frameIntFlag="<<fi_flagpc<<" framesDiscard="<<nFrameDiscard_fipc<<
            " framesCompute="<<nFrameIntegrate_fipc<<" FrameIntDim="<<IMG_DIM_FIpc<<" RegAvgfrmsize="<<FINALFRAMESIZE_REGAVGpc
            <<" QEMCP_tobedone="<<qemcpflagpc<<" CentCorr_tobedone="<<centCorrflagpc<<" CentBias_tobedone="<<centBiasflagpc<<" DetectDist_tobedone="<<DetectDistflagpc<<" OpticDist_tobedone="<<OpticDistflagpc<<
            " pathToOutputTar="<<outtarpathpc<<" search_algo_forFullFrameAst="<<search_algo_ctlgpc<<" Radi_search="<<rad_searchpc<<" len_rect_a="<<(string)len_apc<<" len_rect_b="<<(string)len_bpc<<
            " flag_thetaComp="<<ReturnYorN(flag_thetaComppc)<<" RASfile="<<rasfilepc<<" "<<"clobber ="<<ReturnYorN(clobberpc)<<" history="<<ReturnYorN(historypc)<<" mode="<<modepc<<" Write_todiskuc="<<wtd_ucpc
            <<" Write_todiskbp="<<wtd_bppc<<" Write_todiskff="<<wtd_ffpc<<" Write_todiskpp="<<wtd_pppc<<" Write_todiskCentCorr="<<wtd_centCorrpc<<" Write_todiskCentBias="<<wtd_centBiaspc<<" Write_todisksd="<<wtd_sdpc
            <<" Write_todiskcr="<<wtd_crpc<<" Write_todiskqe="<<wtd_qemcppc<<" Write_todiskfi="<<wtd_fipc<<" Write_todiskdd="<<wtd_ddpc<<" Write_todiskod="<<wtd_odpc<<" LastFileFlag= "<<this->last_FileFlag;
            
     LOG(INFO)<<cmd.str();
    try{
     system(cmd.str().c_str())   ;
    }
    catch(...)
    {   
        LOG(ERROR)<<"***Error in LEVEL2PC generation for NUV";
        return(EXIT_FAILURE);
    }
    }
    else if(strcmp(channel_l2.c_str(),"FUV")==0)
    {
        OutputDir_l2pcFUV.push_back(level2outdirc);
        cmd<<EXE_L2_PC<<" zipFlag = n"<<" level1indir ="<<(string)level1indirpc<<" caldbdir ="<<(string)caldbindir<<" level2outdir ="<<(string)level2outdirc<<" channel ="<<channel_l2<<" startTime= "<<-9999
            <<" endTime="<<-9999<<" unTar_flag= y"
            <<" dropframe = "<<ReturnYorN(dropframepcfuv)<<" utcFlag = "<<ReturnYorN(UTC_flag)<<" crcflag = "<<ReturnYorN(crc_flagpcfuv)<<" parityFlag= "<<parity_flagpcfuv<<" thr_One_cr="<<first_thr_For_Crfuv<<" thr_Two_cr="<<sec_thr_For_Crfuv<<" GTI_FLAG ="<<ReturnYorN(gti_flagpcfuv)<<
            " thresholdMultph ="<<thr_multiphpcfuv<<" algoFlag="<<star_detect_algo_flagpcfuv<<" threshold= "<<sd_multi_factor_defaultpcfuv<<" minimumTargetedstars="<<minimum_No_of_Starspcfuv<<
            " refinedWinSize="<<refine_Winsizepcfuv<<" centroidWinSize="<<centroid_Winsizepcfuv<<" diffDist="<<diff_Distpcfuv<<" catalogpath="<<catalogpathpcfuv<<" att_timecol="<<att_timecolpcfuv<<
            " database_name="<<databasenamepcfuv<<" att_qcol="<<att_qcolpcfuv<<" shiftRotDetAlgoFlag="<<shift_N_Rotate_algopcfuv<<" frameIntFlag="<<fi_flagpcfuv<<" framesDiscard="<<nFrameDiscard_fipcfuv<<
            " framesCompute="<<nFrameIntegrate_fipcfuv<<" FrameIntDim="<<IMG_DIM_FIpcfuv<<" RegAvgfrmsize="<<FINALFRAMESIZE_REGAVGpcfuv
            <<" QEMCP_tobedone="<<qemcpflagpcfuv<<" CentCorr_tobedone="<<centCorrflagpcfuv<<" CentBias_tobedone="<<centBiasflagpcfuv<<" DetectDist_tobedone="<<DetectDistflagpcfuv<<" OpticDist_tobedone="<<OpticDistflagpcfuv<<
            " pathToOutputTar="<<outtarpathpcfuv<<" search_algo_forFullFrameAst="<<search_algo_ctlgpcfuv<<" Radi_search="<<rad_searchpcfuv<<" len_rect_a="<<(string)len_apcfuv<<" len_rect_b="<<(string)len_bpcfuv<<
            " flag_thetaComp="<<ReturnYorN(flag_thetaComppc)<<" RASfile="<<rasfilepc<<" "<<"clobber ="<<ReturnYorN(clobberpc)<<" history="<<ReturnYorN(historypc)<<" mode="<<modepc<<" Write_todiskuc="<<wtd_ucpc
            <<" Write_todiskbp="<<wtd_bppcfuv<<" Write_todiskff="<<wtd_ffpcfuv<<" Write_todiskpp="<<wtd_pppcfuv<<" Write_todiskCentCorr="<<wtd_centCorrpcfuv<<" Write_todiskCentBias="<<wtd_centBiaspcfuv<<" Write_todisksd="<<wtd_sdpcfuv
            <<" Write_todiskcr="<<wtd_crpcfuv<<" Write_todiskqe="<<wtd_qemcppcfuv<<" Write_todiskfi="<<wtd_fipcfuv<<" Write_todiskdd="<<wtd_ddpcfuv<<" Write_todiskod="<<wtd_odpcfuv;
            
     LOG(INFO)<<cmd.str();
    try{
     system(cmd.str().c_str())   ;
     }
    catch(...){
        LOG(ERROR)<<"***Error in LEVEL2PC generation for FUV";
        return(EXIT_FAILURE); 
    }
    }
    else{
        
        LOG(ERROR)<<"***Invalid channel value***"<<channel_l2;
        return(EXIT_FAILURE);
    }
    flag_tarExtracted=TRUE;
    
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    LOG(INFO)<<"\033[1;31m===============UVIT LEVEL2 PC CHAIN FOR "<<channel_l2<<" ENDS===============\033[0m";
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    
return(EXIT_SUCCESS);    
}
 
 int UVIT_DriverModule::RunRAPCChain()
{
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    LOG(INFO)<<"\033[1;31m===============UVIT RELATIVE ASPECT CHAIN  FOR PC STARTED===============\033[0m";
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    stringstream cmd ;
    flag_tarExtracted=FALSE;
    
    
    cmd<<EXE_RA_PC<<" zipFlag ="<<ReturnYorN(zipFlag)<<" level1indir ="<<(string)level1indirrapc<<" caldbdir ="<<(string)caldbindir<<" level2outdir ="<<(string)level2outdirrapc<<" channel ="<<channelRAPC<<" cmpFrames = "<<nCompareFramesrapc<<" thr_one_cr = "<<first_thr_Crpc<<" thr_two_cr ="<<sec_thr_Crpc<<
            " dropframe = "<<ReturnYorN(dropframerapc)<<" utcFlag = "<<ReturnYorN(UTC_flag)<<" crcflag = "<<ReturnYorN(crc_flagrapc)<<" parityFlag = "<<parity_flagrapc<<" GTI_Flag ="<<ReturnYorN(gti_flagrapc)<<" thresholdMultph ="<<thr_multiphrapc<<
            " flatfieldFlag ="<<ReturnYorN(flatfieldFlagrapc)<<" paddingDim="<<padding_dim<<" Nacc="<<no_ofFramesToAcc<<" qemcpFlag="<<ReturnYorN(qemcpFlagrapc)<<" subdivisionFlag ="<<ReturnYorN(subdivisionFlag)<<
            " subdivisionsize ="<<subDivision_size<<" algoFlag ="<<star_detect_algo_flagrapc<<" threshold ="<<rms_mul_factor_defaultrapc<<" minimumTargetedstars ="<<minimum_No_of_Starsrapc<<
            " refinedWinsize ="<<refine_Winsizerapc<<" centroidWinSize="<<centroid_Winsizerapc<<" background_fact ="<<backgrd_factrapc<<" FrameIntDim ="<<IMG_DIM_FIrapc<< " Window_nb ="<<win_sizerapc<<
            " GenmatchstarsFile_flag ="<<ReturnYorN(match_stars_file_flagrapc)<<" ThresholdValue ="<<threshold<<" framesToBeDiscard="<<frames_toDiscardrapc<<" averageFactor ="<<nFrameToAveragerapc<<
            " diffDist ="<<pair_nbhd_distancerapc<<" freqDomainFilterFlag="<<FreqDomainFilterFlagrapc<<" type_Filtering ="<<type_Filteringrapc<<" fitting_flag ="<<ReturnYorN(fitting_flagrapc)<<" orderPitch="<<orderpitchrapc<<
            " orderYaw ="<<orderyawrapc<<" orderRoll ="<<orderrollrapc<<" deltaTime ="<<poly_fit_intervalrapc<<" freqValue="<<freqvaluerapc<<" shiftRotDetAlgoFlag="<<shift_rotation_algorapc<<
            " clobber ="<<ReturnYorN(clobberrapc)<<" flag_thetaComp ="<<ReturnYorN(flag_thetaComprapc)<<" history ="<<ReturnYorN(historyrapc)<<" mode ="<<(string)moderapc<<" Write_todiskuc="<<wtd_ucrapc<<
            " Write_todiskbp="<<wtd_bprapc<<" Write_todiskff="<<wtd_ffrapc<<" Write_todisksd="<<wtd_sd<<" Write_todiskpp="<<wtd_pprapc<<" Write_todiskcr="<<wtd_crrapc<<" Write_todiskqe="<<wtd_qemcprapc<<
            " Write_todiskfi="<<wtd_firapc<<" Write_todisksc="<<wtd_fscrapc<<" Write_todiskdd="<<wtd_ddrapc<<" Write_todiskod="<<wtd_odrapc<<" Write_todiskCentCorr= "<<wtd_centCorrrapc<<
            " Write_todiskCentBias="<<wtd_centBiasrapc<<" framesDiscard ="<<nFrameDiscard_firapc<<" framesCompute= "<<nFrameIntegrate_firapc;
   
    
    LOG(INFO)<<cmd.str();
    try{
     system(cmd.str().c_str())   ;
    }
    catch(...)
    {   
        LOG(ERROR)<<"***Error in relativeAspect generation for PC mode";
        return(EXIT_FAILURE);
    }
    flag_tarExtracted=TRUE;
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    LOG(INFO)<<"\033[1;31m===============UVIT RELATIVE ASPECT CHAIN FOR PC ENDS ===============\033[0m";
    LOG(INFO)<<"\033[1;31m================================================================\033[0m";
    
    return(EXIT_SUCCESS);
}
 int UVIT_DriverModule::copyAllheaderKeys (char* infile)
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
//        if (keyclass == TYP_COMM_KEY)
//            continue ;
       // if (keyclass == TYP_USER_KEY || keyclass == TYP_REFSYS_KEY || keyclass == TYP_WCS_KEY || keyclass == TYP_COMM_KEY)
       // {
            header_info.push_back (record) ;
        //}
    }

    fits_close_file (fin , &status) ;
    // delete[] record;
    return (EXIT_SUCCESS) ;
}






