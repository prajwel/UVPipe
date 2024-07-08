/* 
 * File:   uvtImL2_commonArray.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <stdlib.h>

#include "uvtLevel2IM.h"

#include<Directory.h>
#include<DataInfo.h>
#include<DataIngest.h>
#include <algorithm>

#include<set>

#include <vector>
#include<memory.h>
#include<spMatrix.h>
#include<glog/logging.h>
#include<uvtUtils.h>
#include<macro_def.h>
#include<transform.h>
//#include<uvtComputeDrift.h>
#define hj 1
#define  NBHD_RADIUS 5
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr   //used in doFft() function
#define moduleoutdir_unit "uvtUnitConvertion"
#define moduleoutdir_badpix "uvtMaskBadPix"
#define moduleoutdir_fltfield "uvtFlatFieldCorr"
#define moduleoutdir_qe "uvtQEMCPCorr"
#define moduleoutdir_pixpad "uvtPixPadding"
#define moduleoutdir_subdiv "uvtSubDivision"
#define moduleoutdir_cosmicray "uvtCosmicRayCorrection"
//#define moduleoutdir_AccEverytsec "uvtAccEveryTsec"
//#define moduleoutdir_findstarcentroid "uvtFindStarCentroid"
#define moduleoutdir_detectordistortion "uvtDetectDistCorr"
#define moduleoutdir_opticaldistortion "uvtOpticAssDistCorr"
//#define moduleoutdir_driftExercise "uvtComputeDrift"
//#define moduleoutdir_refFrameCal "uvtRefFrameCal"
#define moduleoutdir_shiftNRot "uvtShiftRot"
#define moduleoutdir_findWtdmean "uvtFindWtdMean"
#define moduleoutdir_RegAvg "uvtRegAvg"
//#define IMAGE_ARRAYSIZE  4800

bool compare1 (struct Star_l2 star1 , struct Star_l2 star2);
bool compare1 (struct Star_l2  star1 , struct Star_l2 star2)
{
    return (star1.intensity > star2.intensity) ;
}
uvtLevel2IM::uvtLevel2IM ()
{
    tar_extracted_flag_IM=FALSE;
}

uvtLevel2IM::~uvtLevel2IM ()
{
   

}

int uvtLevel2IM::readPILParameters (int argc , char** argv)
{
    int status = 0 ;
    char temp[PIL_LINESIZE] ;

    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG (INFO) << "***Error Initializing PIL***" ;
        return status ;
    }
     if(this->paramfile_Varask_iml2==FALSE){
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
     }
    if (PIL_OK != (status = PILGetString ("channel" , temp)))
    {
        LOG (INFO) << endl << "***Error reading channel***" ;
        return status ;
    }
    channel.assign (temp);
 if (PIL_OK != (status = PILGetBool ("utcFlag" , &UTC_flag)))
    {
        LOG (INFO) << endl << "***Error reading drop parameters***" ;
        return status ;
    }
if (PIL_OK != (status = PILGetInt ("crcflag" , &crc_flag)))
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
//     if (PIL_OK != (status = PILGetInt ("att_flagVal" , &att_flag_val)))
//    {
//        LOG (INFO) << endl << "***Error Reading att  flag value :" << padding_dim << "***" ;
//        return status ;
//    }

    if (PIL_OK != (status = PILGetBool ("GTI_FLAG" , &gti_flag)))
    {
        LOG (INFO) << endl << "***Error Reading GTI Flag value***" ;
        return status ;
    }
    if (gti_flag)
    {
        all_or_custom = 1 ;
        valid_bit = 1 ;
    }

    if (PIL_OK != (status = PILGetInt ("paddingDim" , &padding_dim)))
    {
        LOG (INFO) << endl << "***Error Reading padding_dim :" << padding_dim << "***" ;
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
    
    
    if (PIL_OK != (status = PILGetReal ("ThresholdValue" , &cr_threshold)))
    {
        LOG (INFO) << endl << "***Error reading the threshold value for Cosmic Ray correction module***" ;
        return status ;
    }
   
    if (PIL_OK != (status = PILGetReal4 ("diffDist" , (float*) &diff_Dist)))
    {
        LOG (INFO) << endl << "***Error reading diffDist parameter***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetInt ("RegAvgfrmsize" , &FINALFRAMESIZE_REGAVG)))
    {
        LOG (INFO) << endl << "***Error reading size of registered image  ***" ;
        return status ;
    }
    if(subdivisionFlag==1)
    {
        if(FINALFRAMESIZE_REGAVG>subDivision_size)
        {
        LOG(ERROR)    <<"***Subsampling size can not be less than frame size of Accumulated frames..!!!!***";
                return(EXIT_FAILURE);
        }
    }
  if(this->paramfile_Varask_iml2==FALSE)
  {
   if (PIL_OK != (status = PILGetFname ("RASfile" , rasfile)))
            {
                LOG(ERROR) << endl << "\033[1;31m***Error reading Relative Aspect file path***" ;
                return status ;
            }
  }
     if (PIL_OK != (status = PILGetString ("att_timecol" , att_timecol)))
    {
      LOG (INFO) << endl << "***Error reading name of time column in attitude file***" ;
    return status ;
    }
    if (PIL_OK != (status = PILGetString ("att_qcol" , att_qcol)))
    {
        LOG (INFO) << endl << "***Error reading name of quaternion column name in attitude file ***" ;
        return status ;
    }
   
  if (PIL_OK != (status = PILGetReal4 ("threshold" , &sd_multi_factor_default)))
    {
        LOG (INFO) << endl << "***Error reading standard deviation multiplication factor value***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("shiftRotDetAlgoFlag" , &star_detect_algo_flag)))
    {
        LOG (INFO) << endl << "***Error reading shiftRotDetAlgoFlag***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetInt ("no_of_weightedframes" , &no_ofWeigh)))
    {
        LOG(ERROR) << endl << "\033[1;31m***Error reading number of frames to be weighted in weighted mean module ***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetInt ("minimum_targetedstars" , &minimum_No_of_Stars)))
    {
        LOG(ERROR) << endl << "\033[1;31m***Error reading minimum number of stars required ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("database_name" , databasename)))
    {
        LOG (INFO) << endl << "***Error reading  database name of star catalogue***" ;
        return status ;
    }
    
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
    
     if (PIL_OK != (status = PILGetBool ("flag_thetaComp" , &flag_thetaComp)))
    {
        LOG (INFO) << endl << "***Error Reading clobber:" << clobber << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskuc" , &wtd_uc)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing to Disk Flag For unit Conversion:" << "***" ;
        return status ;
    }
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
    if(subdivisionFlag==1){
    if (PIL_OK != (status = PILGetInt ("Write_todisksd" , &wtd_sd)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the subDivision:" << "***" ;
        return status ;
    }
    }
    if (PIL_OK != (status = PILGetInt ("Write_todiskcr" , &wtd_cr)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the cosmicRayCorrection:" << "***" ;
        return status ;
    }
      if (PIL_OK != (status = PILGetInt ("Write_todiskqe" , &wtd_qemcp)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the qemcpcorrection:" << "***" ;
        return status ;
    }
    //if (PIL_OK != (status = PILGetInt ("Write_todiskac" , &wtd_ac)))
    //{
      //  LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For the Accumulated module:" << "***" ;
       // return status ;
    //}
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
    if (PIL_OK != (status = PILGetInt ("Write_todisksnr" , &wtd_snr)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For shift and rotate:" << "***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetInt ("Write_todiskwm" , &wtd_wm)))
    {
        LOG(ERROR) << endl << "***Error Reading Writing the Disk Flag For  weighted mean" << "***" ;
        return status ;
    }


    PILClose (status) ;
    return (EXIT_SUCCESS) ;
}

int uvtLevel2IM::uvtLevel2IMProcess ()
{
  
     int status = 0 ; //status flag to check return status of functions
    string temp_str;
    temp_str.assign (level1indir);
  
    if(tar_extracted_flag_IM==FALSE)
    {
        
    level1indir="";
   status= extractTars (temp_str,level1indir,orbnum);;//extract level-1 data  tar file
    if(status)
    {
       LOG(INFO)<<"Error in extracting tar";
       return(EXIT_FAILURE);
    }
   
    }
    
    if (!DirExists ((char *) level1indir.c_str ()))
    {
        LOG(ERROR) << endl << "Input level 1 directory not found " << endl ;
        return (EXIT_FAILURE) ;
    }
  
    //creating output directory if it does not exists  
    if (!DirExists ((char *) level2outdir.c_str ()))
    {
        string cmd = "mkdir -p " + level2outdir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists ((char*)level2outdir.c_str ()) && clobber==YES){
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
    /*---NOTE : Change 'NUV' to 'VIS'  for valid chain operation, presently using NUV due to non availability of data*/
    if (dirobj.setup (level1indir , level2outdir , temp_channel_id))
    { //IM RA  needs VIS data only //setup done only if uvitV directory found
        LOG(ERROR) << endl << "***Error in directory set up***" << endl ;
        LOG(ERROR) << endl << "***This chain runs for VIS data only...check for uvitV directory in input*** " << endl ;
        return (EXIT_FAILURE) ;
    }
    int numofsciencedatafiles = dirobj.sciencedatafile.size () ;//total number of science data file
    if (numofsciencedatafiles <= 0)
    {
        LOG(ERROR) << endl << "No science data file found" << endl ;
        return (EXIT_FAILURE) ;
    }

    char obsmode[FLEN_VALUE] ;
    char moduleIndir[FLEN_FILENAME] ; //hold path for input directory for every module, will be updated after every module run
    char outputdir[FLEN_FILENAME] ;
     int division_fact;


      if (subdivisionFlag == 0)
    {
        subDivision_size = padding_dim ;
    }
    //memory allocation for pixel padding module
    frame_Data_Padded = new float[padding_dim * padding_dim] ;
    frame_ExpoData_padded = new float[padding_dim * padding_dim] ;
    if(subdivisionFlag)
    {
    frame_Data_subdivided = new float[subDivision_size * subDivision_size] ;
    frame_ExpData_subdivided = new float[subDivision_size * subDivision_size] ;
       division_fact = subDivision_size / padding_dim ;
        division_fact = division_fact*division_fact ;
    }
    //reading RelativeAspect file for input.
   
    
     
     /*==============================================================================================================================*/
     /*====================READING  OF RELATIVE ASPECT FILE GENERATED AFTER  RELATIVE ASPECT CHAIN===================================================*/
     /*==============================================================================================================================*/
    /**Reading Relative Aspect File from RA series**/
     
    
    LOG(INFO) << "\nReading RAS file from RA series..." << endl ;
     fitsfile *fras ;
//      long  no_of_records=0;
//    fits_open_file (&fras , rasfile , READONLY , &status) ;
//    printError (status , "Error in opening the rasfile",rasfile) ;
////    fits_movabs_hdu (fras , 2 , NULL , &status) ;
////    printError (status , "Error in moving the 2nd HDU of  ras file",rasfile) ;
//     fits_movabs_hdu (fras , 2 , NULL , &status) ;
//    printError (status , "Error in moving the 2nd HDU of  ras file",rasfile) ;
//    fits_get_num_rows (fras , &no_of_records , &status) ;
//    printError (status , "Error Reading the number of Rows",rasfile) ;
//
//  double * time_ras = new double[no_of_records] ;
//     double *roll_ras = new double [no_of_records] ;
//    double *pitch_ras = new double [no_of_records] ;
//     double *yaw_ras = new double[no_of_records] ;
//    
//    fits_read_col (fras , TDOUBLE , 1 , 1 , 1 , no_of_records , NULL , time_ras , NULL , &status) ;
//    printError (status , "***Error reading  centroid x***") ;
//    fits_read_col (fras , TDOUBLE , 2 , 1 , 1 , no_of_records , NULL , yaw_ras , NULL , &status) ;
//    printError (status , "***Error writing centroid y***") ;
//    fits_read_col (fras , TDOUBLE , 3 , 1 , 1 , no_of_records , NULL , pitch_ras , NULL , &status) ;
//    printError (status , "***Error writing x-correction***") ;
//    fits_read_col (fras , TDOUBLE , 4 , 1 , 1 , no_of_records , NULL , roll_ras , NULL , &status) ;
//    printError (status , "***Error writing y-correction***") ;
//    fits_close_file (fras , &status) ;
//    printError (status , "Error in closing the file",rasfile) ;
//    
     /*==============================================================================================================================*/
     /*====================READING  OF RELATIVE ASPECT FILE GENERATED AFTER  RELATIVE ASPECT CHAIN FINISHED ===================================================*/
     /*==============================================================================================================================*/
    //variable declared for shift and rotate
    float **subSignalArray , **subExposureArray , **tempSigarr , **tempExparr ;
    //Variables for  Find Weighted mean 
        
    subSignalArray = new float * [subDivision_size] ;
    subExposureArray = new float * [subDivision_size] ;
    tempSigarr = new float * [subDivision_size] ;
    tempExparr = new float * [subDivision_size] ;
    for (int i = 0 ; i < subDivision_size ; i++)
    {
        subExposureArray[i] = new float [subDivision_size] ;
        subSignalArray[i] = new float [subDivision_size] ;
        tempExparr[i] = new float [subDivision_size] ;
        tempSigarr[i] = new float [subDivision_size] ;
    }

    //memory allocation for storing shifts and theta
//    double *delta_x = new double [no_of_records] ;
//    double *delta_y = new double [no_of_records] ;
//    double *delta_theta = new double[no_of_records] ;
    
    double temp_pitch,temp_yaw;
     
     /*==============================================================================================================================*/
     /*+++++====================LOOP FOR NUMBER OF SCIENCE DATA FILE IN LEVEL-1 DIRECTORY=======================================================*/
     /*==============================================================================================================================*/
//method for number of science data file  
    
   for (int dataindex = 0 ; dataindex < numofsciencedatafiles ; dataindex++)
    {

        LOG(ERROR) << endl << "------------------Data Set " << dataindex + 1 << " : " << dirobj.sciencedatafile[dataindex] << "-----------------------" << endl ;
        /*---finding mode from science data file---*/
        getKeywordVal ((char *) dirobj.sciencedatafile[dataindex].c_str () , "OBS_MODE" , 1 , obsmode) ;
        if (strcasecmp (obsmode , "IM") != 0)
        {
            LOG(ERROR) << endl << "Observation mode is " << obsmode << "  in file  " << dirobj.sciencedatafile[dataindex] ;
            LOG(ERROR) << endl << "Checking next file...." << endl ;
            continue ; //go to next file if obs mode is not IM
        }
        string ras_finalFile;
        string rasFile_str=rasfile;
//      status=readpathOfRASfile (dirobj.sciencedatafile[dataindex],(char*)channel.c_str (),rasFile_str,ras_finalFile);
//      if (status)
//        {
//            LOG (ERROR) << "Error in reading the content of RAS file path from Relative Aspect directory" ;
//            return (EXIT_FAILURE) ;
//        }
  //  sprintf(rasfile,"%s",ras_finalFile);  
    LOG(INFO)<<"Relative Aspect file ->"<<ras_finalFile;
    long  no_of_records=0;
    fits_open_file (&fras ,rasfile , READONLY , &status) ;
    printError (status , "Error in opening the rasfile",(char*)ras_finalFile.c_str ()) ;
//    fits_movabs_hdu (fras , 2 , NULL , &status) ;
//    printError (status , "Error in moving the 2nd HDU of  ras file",rasfile) ;
     fits_movabs_hdu (fras , 3 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of  ras file",(char*)ras_finalFile.c_str ()) ;
    fits_get_num_rows (fras , &no_of_records , &status) ;
    printError (status , "Error Reading the number of Rows",(char*)ras_finalFile.c_str ()) ;

  double * time_ras = new double[no_of_records] ;
     double *roll_ras = new double [no_of_records] ;
    double *pitch_ras = new double [no_of_records] ;
     double *yaw_ras = new double[no_of_records] ;
    
    fits_read_col (fras , TDOUBLE , 1 , 1 , 1 , no_of_records , NULL , time_ras , NULL , &status) ;
    printError (status , "***Error reading  centroid x***") ;
    fits_read_col (fras , TDOUBLE , 2 , 1 , 1 , no_of_records , NULL , yaw_ras , NULL , &status) ;
    printError (status , "***Error writing centroid y***") ;
    fits_read_col (fras , TDOUBLE , 3 , 1 , 1 , no_of_records , NULL , pitch_ras , NULL , &status) ;
    printError (status , "***Error writing x-correction***") ;
    fits_read_col (fras , TDOUBLE , 4 , 1 , 1 , no_of_records , NULL , roll_ras , NULL , &status) ;
    printError (status , "***Error writing y-correction***") ;
    fits_close_file (fras , &status) ;
    printError (status , "Error in closing the file",(char*)ras_finalFile.c_str ()) ;
    
    
    double *delta_x = new double [no_of_records] ;
    double *delta_y = new double [no_of_records] ;
    double *delta_theta = new double[no_of_records] ;
    
   LOG(INFO) << endl << "Data Mode is " << obsmode << endl ;
    strcpy (outputdir , dirobj.level2path[dataindex].c_str ()) ;
         
        //setting path for output directory to be generated 
    sprintf (moduleoutdir_bp , "%s/%s_%s", outputdir, moduleoutdir_badpix,VERSION) ;
    sprintf (moduleoutdir_uc , "%s/%s_%s" , outputdir ,moduleoutdir_unit,VERSION) ;
    sprintf (moduleoutdir_ff , "%s/%s_%s" ,outputdir, moduleoutdir_fltfield,VERSION ) ;
    sprintf (moduleoutdir_qemcp , "%s/%s_%s" ,outputdir, moduleoutdir_qe,VERSION ) ;
    sprintf (moduleoutdir_pp , "%s/%s_%s" , outputdir ,moduleoutdir_pixpad,VERSION) ;
    sprintf (moduleoutdir_sd , "%s/%s_%s" , outputdir , moduleoutdir_subdiv,VERSION) ;
    sprintf (moduleoutdir_cr , "%s/%s_%s" , outputdir , moduleoutdir_cosmicray,VERSION) ;

    sprintf (moduleoutdir_dd , "%s/%s_%s" ,outputdir ,moduleoutdir_detectordistortion,VERSION) ;
    sprintf (moduleoutdir_od , "%s/%s_%s" ,outputdir , moduleoutdir_opticaldistortion,VERSION) ;
  
        sprintf (moduleoutdir_snr , "%s/%s_%s" , outputdir , moduleoutdir_shiftNRot,VERSION) ;
        sprintf (moduleoutdir_wm , "%s/%s_%s" , outputdir , moduleoutdir_findWtdmean,VERSION) ;
        sprintf (moduleoutdir_ravg , "%s/%s_%s" , outputdir , moduleoutdir_RegAvg,VERSION) ;
        
        
        //for DataIngest
        //----------DATAINGEST----------//
        LOG(ERROR) << endl << "===================================DATAINGEST===================================================" << endl ;
        DataIngest di_obj ;
         di_obj.read ((char *) dirobj.sciencedatafile[dataindex].c_str () ,(char*)caldbindir.c_str () ,(char *) dirobj.tctfile.c_str () ,(char *) dirobj.mkffile.c_str (), (char *) dirobj.gtifile[dataindex].c_str () ,
                (char *) dirobj.lbtfile.c_str () , (char*)dirobj.attfile.c_str (),(char*) dirobj.darkDirectory.c_str () , att_flag_val,gti_flag , valid_bit , all_or_custom , outputdir , dropframe , parity_flag ,UTC_flag,crc_flag, clobber , history) ;
        di_obj.display () ;
        status = di_obj.DataIngestProcess () ;
        if (status)
        {
            LOG(ERROR) << endl << "***Error in Data Ingest Process***" << endl ;
            continue;
        }

        strcpy (moduleIndir , di_obj.getModuleOutdir ()) ;
         strcpy (dataIngest_out_dir , di_obj.getModuleOutdir ()) ;//copy output path of dataingest directory for taking dark information
        LOG(INFO) << "Using directory " << moduleIndir << "  as input to uvtUnitConversion" ;


   /*Reading information file generated after  DataIngest directory*/
    char infofile_in[PIL_LINESIZE] ;
    string  tempfilepath = searchFile (moduleIndir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG(ERROR) << endl << "***Information file not found in " << moduleIndir << "***" << endl ;
        continue;
    }
    sprintf (infofile_in , "%s/%s" , moduleIndir , tempfilepath.c_str()) ;
    LOG(INFO) << "Information file :" << infofile_in << endl ;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists (infofile_in)))
    {
        LOG(ERROR) << endl << "***Input FileList not Found at Specified PATH,Check Input Direcrory***" << endl ;
       continue;
    }
    /*
 open the .info FITS file  and read the header information from the second HDU.
     */
vector<string> header_info;
getHistory(header_info);
writeHistory(infofile_in,header_info);
    fitsfile *finfo_in  ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input information file",infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU") ;
    datainfo.getInfo(finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    
       int nframes ;
    char nameprefix[PIL_LINESIZE] ;
     char infile[NAMESIZE] ;
    char errstr[FLEN_FILENAME] ;
     unsigned short frameno ;
    double frametime , integrationtime,mid_time ;
     //reading keywords from information file
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file
    fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
    printError (status , "***Error in reading the  key value of the NFILES ***" , infofile_in) ;
    fits_read_key (finfo_in , TSTRING , "DARKDIR" , darkdir , NULL , &status) ;
    printError (status , "***Error in reading the  key value of  DARKDIR ***" , infofile_in) ;
    char **sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
    fits_movabs_hdu (finfo_in , 1 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU") ;
    fits_read_key (finfo_in , TINT , "WIN_X_SZ" , &win_xsize , NULL , &status) ;
    printError (status , "Error in reading the key value of the Window xsize" , infofile_in) ; //for creating name for output information file
    fits_read_key (finfo_in , TINT , "WIN_Y_SZ" , &win_ysize , NULL , &status) ;
    printError (status , "Error in reading the key value of the Window ysize " , infofile_in) ; //for creating name for output information
    //reading frame names from information file into vector
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU") ;
    fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) sigframelist , NULL , &status) ;
    printError (status , "Error in reading  column of signal frame list") ;
     
    fits_close_file (finfo_in , &status) ;
    printError (status , "Error in reading the column value of the Input Signal List" , infofile_in) ;
    
    //added
    double dataset_starttime,dataset_endtime;
     sprintf (infile , "%s/%s/%s" , moduleIndir , "SignalFrames" , sigframelist[0]) ;
     readKeywords (infile,1,1,TDOUBLE,"FRMTIME",&dataset_starttime);
      sprintf (infile , "%s/%s/%s" , moduleIndir , "SignalFrames" , sigframelist[nframes-1]) ;
     readKeywords (infile,1,1,TDOUBLE,"FRMTIME",&dataset_endtime);
     LOG(INFO)<<dataset_endtime<<" "<<dataset_starttime;
     if(time_ras[0]>dataset_endtime || time_ras[no_of_records-1]<dataset_starttime){
         LOG(ERROR)<<"Error in finding the shifts -> No time matching between RAS file and Dataset Time";
         continue;
     }
    //till this
    
    
    //transformations RPY to DXDYDTHETA of RAS file based on Detector
     if (strcasecmp (datainfo.getDetector () , "NUV") == 0)//incase of NUV
        {
            for(int i=0;i<no_of_records;i++)
            {
                temp_pitch=pitch_ras[i]*3600;
                temp_yaw=yaw_ras[i]*3600;
                transformRPYtoDXDYDTHETA_NUV (roll_ras[i],temp_pitch,temp_yaw,delta_x[i],delta_y[i],delta_theta[i]);//for converting RPY to DX,DY and DTHETA
               
            }
           
        }
    
     else if(strcasecmp (datainfo.getDetector () , "FUV") == 0)//incase of FUV
        {
            for(int i=0;i<no_of_records;i++)
            {
                   temp_pitch=pitch_ras[i]*3600;
                temp_yaw=yaw_ras[i]*3600;
                transformRPYtoDXDYDTHETA_FUV (roll_ras[i],temp_pitch,temp_yaw,delta_x[i],delta_y[i],delta_theta[i]);//for converting RPY to DX,DY and DTHETA
            }
        }
     else if(strcasecmp (datainfo.getDetector () , "VIS") == 0)//incase of VIS
        {
            for(int i=0;i<no_of_records;i++)
            {
                   temp_pitch=pitch_ras[i]*3600;
                temp_yaw=yaw_ras[i]*3600;
                transformRPYtoDXDYDTHETA_VIS (roll_ras[i],temp_pitch,temp_yaw,delta_x[i],delta_y[i],delta_theta[i]);//for converting RPY to DX,DY and DTHETA            
            
            }
             
        }
        else
        {
            LOG(INFO)<<"***Invalid Channel***"<<endl;
            continue;
        }
  
    for (int i=0;i<no_of_records;i++)
   {
  delta_theta[i]=roll_ras[i];
 delta_y[i]=pitch_ras[i];
            delta_x[i]=yaw_ras[i];
           LOG(INFO)<<delta_x[i]<<" "<<delta_y[i]<<" "<<delta_theta[i];
    }
//   
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    
    /*==============================================================================================================*/
    /*===================================READING STARTED FOR  CALDB FILES =================================================*/
    /*==============================================================================================================*/
    
    char caldb_common_dir[PIL_LINESIZE] ,caldb_common_dir_qefile[PIL_LINESIZE];
    char caldb_temp_dir[PIL_LINESIZE],caldb_common_dir_ff[FLEN_FILENAME] ;
    //copying caldb directory path for different module for further use.
     strcpy (caldb_temp_dir , (char*)caldbindir.c_str ()) ;
     strcpy (caldb_common_dir ,(char*) caldbindir.c_str ()) ;
     strcpy (caldb_common_dir_qefile ,(char*) caldbindir.c_str ()) ;
     strcpy (caldb_common_dir_ff ,(char*) caldbindir.c_str ()) ;
     
        //caldb reading started for bad pixel correction 
     string tempname = caldb_handler.getBadPixelFile (datainfo.getDetector () , datainfo.getObsMode () ,win_xsize+1 ,win_ysize+1 , (char*)caldbindir.c_str ()) ;
    if (tempname == " ")
    {
        LOG(ERROR) << endl << "Couldn't find bad pixel file from calDB" << endl ;
        continue;
    }
    joinStrings (badpixfile , 1 , tempname.c_str()) ;
    badpixdata=new float[xsize*ysize];
    status =readImage(badpixfile,1,badpixdata,xsize,ysize);
  
     //caldb reading started for flat field correction 
string   tempname1 = caldb_handler.getFlatFieldFile (datainfo.getDetector () , datainfo.getObsMode () , datainfo.getFilter () , caldb_common_dir_ff) ;
    joinStrings (flatfieldfile , 2 ,caldb_common_dir_ff, tempname1.c_str()) ; 
    flatfielddata=new float[xsize*ysize];    
    status =readImage(flatfieldfile,1,flatfielddata,xsize,ysize);
  
    //caldb reading started for detector distortion
    tempname = caldb_handler.getDetectorFile (datainfo.getDetector () , caldb_temp_dir) ;
    joinStrings (detector_distortion_corr_file , 2 , caldb_temp_dir , tempname.c_str()) ;
    LOG(INFO) <<" Detector Distortion correction file is " << detector_distortion_corr_file ;
    if (tempname == " ")
    {
        LOG(ERROR) << endl << "Couldn't find detector  distortion file From caldb" << endl ;
       continue;
    }
    
   X_detect_distortion= new float[xsize*ysize]; 
   Y_detect_distortion= new float[xsize*ysize];   
   status=caldb_handler.readCaldbDistFile(X_detect_distortion,Y_detect_distortion,detector_distortion_corr_file);
   
   //caldb reading started for optical distortion correction
   tempname = caldb_handler.getOpticalDetectorFile (datainfo.getDetector () , caldb_common_dir , "F0") ;
   joinStrings (optical_distortion_corr_file , 2 , caldb_common_dir , tempname.c_str()) ;
   if (tempname == " ")     
    {
        LOG(ERROR) << endl << "Couldn't find Optical Dist detector File From caldb" << endl ;
      continue;
    }
    //Array for storing CALDB optical distortion correction data
   X_optical_distortion= new float[xsize*ysize];    Y_optical_distortion= new float[xsize*ysize];   
   status=caldb_handler.readCaldbDistFile (X_optical_distortion,Y_optical_distortion,optical_distortion_corr_file);
   double factor;
   //caldb reading started for QE and MCP correction if QE and MCP correction to be done
   //if(qe_mcpFlag)
   //{
   tempname=caldb_handler.getQEFile (datainfo.getDetector (),datainfo.getObsMode (),caldb_common_dir_qefile);
   if (tempname == " ")
    {
        LOG(ERROR) << endl << "Couldn't find QEMCP file from caldb" << endl ;
        continue;
    }
    joinStrings (qe_factor_file , 2 , caldb_common_dir_qefile , tempname.c_str()) ;
     status= readNumRowsFromFits (qe_factor_file,2,nCalDBTempValues);
   if (status)
            {
                LOG (ERROR) << "Error in reading the number of rows from fits file "<<qe_factor_file ;
                return (EXIT_FAILURE) ;
            }
  
      cal_temp= new float[nCalDBTempValues];
     cal_f0 = new float[nCalDBTempValues],cal_f1=new float[nCalDBTempValues]; cal_f2=new float[nCalDBTempValues];
    cal_f3 = new float[nCalDBTempValues], cal_f4=new float[nCalDBTempValues]; cal_f5=new float[nCalDBTempValues];
     cal_f6 = new float[nCalDBTempValues],cal_f7=new float[nCalDBTempValues];
     
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
   //reading time,inside temperature and outside temperature from LBT file.
    status = getTemp ((char *)dirobj.lbtfile.c_str (),datainfo.getDetector (),time_lbt,insideTemp,outsideTemp,nrows_lbt) ;
    if (status)
    {
        LOG(ERROR) << "***temperature reading from the lbt file unsuccessful***" ;
        continue;
    }

    /*==============================================================================================================*/
    /*===================================READING ENDED  FOR  CALDB FILES =================================================*/
    /*==============================================================================================================*/
     int filter_coln ;
    int filternumber ;
    sprintf(filter,"%s",datainfo.getFilter ());
    //Array for storing factor values from caldb based on  filter value(i.e column in caldb)
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
        LOG(ERROR)<<endl<<"***Invalid filter option*** "<<endl;
       continue;
    }


    //incase of Darkframe subtraction to be done or not
   if (darkframe_flag)
    { 
       status = takeDarkinfo () ;       
       t_darkframestart = readDarkFrame (dstartpath , darkFramestart_data) ;
        LOG(INFO) << "the start time:: " << t_darkframestart << endl ;
        t_darkframeend = readDarkFrame (dendpath , darkFrameend_data) ;
        LOG(INFO) << "the end  time:: " << t_darkframeend << endl ;
        darkCompute_array = new float[xsize * ysize] ;
    }
     //setting directory structure for output directory
    if (wtd_uc == 1)//for unit Conversion
    {
        status = setDirectoryStructure (moduleoutdir_uc , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
    }
   
    if (wtd_bp == 1)// for Mask Bad pixel
    {
        status = setDirectoryStructure (moduleoutdir_bp , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             return(EXIT_FAILURE);
        }
        status = setDirectoryStructure (moduleoutdir_bp , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
    }
    if (wtd_cr == 1)// for cosmicRay corr
    {
        status = setDirectoryStructure (moduleoutdir_cr , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
        status = setDirectoryStructure (moduleoutdir_cr , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             return(EXIT_FAILURE);
        }
    }
    if (wtd_ff == 1)// for flatfield corr
    {
        status = setDirectoryStructure (moduleoutdir_ff , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             return(EXIT_FAILURE);
        }
        status = setDirectoryStructure (moduleoutdir_ff , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
    }
   if (wtd_qemcp) // for qemcp corr
    {
        status = setDirectoryStructure (moduleoutdir_qemcp , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             return(EXIT_FAILURE);
        }
        status = setDirectoryStructure (moduleoutdir_qemcp , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
    }
    if (wtd_pp )// for pixpadding 
    {
       status = setDirectoryStructure (moduleoutdir_pp , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
        status = setDirectoryStructure (moduleoutdir_pp , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
    }
    if (wtd_sd == 1 && subdivisionFlag==1)// for sub division
    {
        status = setDirectoryStructure (moduleoutdir_sd , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
        status = setDirectoryStructure (moduleoutdir_sd , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             return(EXIT_FAILURE);
        }
    }
  
    if (wtd_dd == 1) // for Detector Distortion correction
    {
        status = setDirectoryStructure (moduleoutdir_dd , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);

        }
        status = setDirectoryStructure (moduleoutdir_dd , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);

        }
    }
    if (wtd_od == 1)// for Optical Distortion
    {
        status = setDirectoryStructure (moduleoutdir_od , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             return(EXIT_FAILURE);
        }
         status = setDirectoryStructure (moduleoutdir_od , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             return(EXIT_FAILURE);
        }

    }  
     if (wtd_snr == 1)//for shift and rotate
    {
        status = setDirectoryStructure (moduleoutdir_snr , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
         status = setDirectoryStructure (moduleoutdir_snr , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }

    }
    if (wtd_wm == 1)//for weighted mean
    {
        status = setDirectoryStructure (moduleoutdir_wm , "SignalFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
             return(EXIT_FAILURE);
        }
         status = setDirectoryStructure (moduleoutdir_wm , "ExposureFrames") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
           return(EXIT_FAILURE);
        }

    }  
   
   //for registration and averaging
        status = setDirectoryStructure (moduleoutdir_ravg , "") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
         status = setDirectoryStructure (moduleoutdir_ravg , "") ;
        if (status)
        {
            LOG(ERROR) << "***Directory Structure has not been successfully set-up***" << endl ;
            return(EXIT_FAILURE);
        }
     
   
 
    float *sum_ac = new float[subDivision_size * subDivision_size] ;
    float *sum_ac_exp = new float[subDivision_size * subDivision_size] ;
     initArray (sum_ac,subDivision_size*subDivision_size,0.0f);
     initArray (sum_ac_exp,subDivision_size*subDivision_size,0.0f);
      float *  caldb_xdist_final=new float[subDivision_size * subDivision_size];
      float *caldb_ydist_final=new float[subDivision_size * subDivision_size];
      float *  caldb_xdist_final_od=new float[subDivision_size * subDivision_size];
      float *caldb_ydist_final_od=new float[subDivision_size * subDivision_size];
      
       float *  caldb_xdist_final_padded=new float[PIX_PADD_SIZE*PIX_PADD_SIZE];
      float *caldb_ydist_final_padded=new float[PIX_PADD_SIZE*PIX_PADD_SIZE];
      float *  caldb_xdist_final_od_padded=new float[PIX_PADD_SIZE*PIX_PADD_SIZE];
      float *caldb_ydist_final_od_padded=new float[PIX_PADD_SIZE*PIX_PADD_SIZE];
      
      for (int p = 0; p < PIX_PADD_SIZE*PIX_PADD_SIZE; p++) {
        caldb_xdist_final_padded[p] = 0.0;
        caldb_ydist_final_padded[p] = 0.0;
        caldb_xdist_final_od_padded[p]=0.0;
        caldb_ydist_final_od_padded[p]=0.0;
    }
      
      
    for (int p = 0; p < subDivision_size * subDivision_size; p++) {
        caldb_xdist_final[p] = 0.0;
        caldb_ydist_final[p] = 0.0;
        caldb_xdist_final_od[p]=0.0;
        caldb_ydist_final_od[p]=0.0;
    }
    /*padding process starts here...padding is applied because 512*512 scheme must be converted into 600*600 scheme***/
      //for detector distortion
      float fact_div=subDivision_size/PIX_PADD_SIZE;
      fact_div=fact_div*fact_div;      
      
      for (int i=0;i<CALDB_DIST_SIZE*CALDB_DIST_SIZE;i++)
      {
          X_detect_distortion[i]= X_detect_distortion[i]/fact_div;
          Y_detect_distortion[i]=Y_detect_distortion[i]/fact_div;
          X_optical_distortion[i]=X_optical_distortion[i]/fact_div;
          Y_optical_distortion[i]=Y_optical_distortion[i]/fact_div;
      }
           status = Applypadding(X_detect_distortion,CALDB_DIST_SIZE,CALDB_DIST_SIZE,caldb_xdist_final_padded,PIX_PADD_SIZE,PIX_PADD_SIZE);
          if (status) {
        LOG(ERROR) << "***Padding unsuccessfull for x-distortion***" ;
           continue;
    }
    status =Applypadding(Y_detect_distortion,CALDB_DIST_SIZE,CALDB_DIST_SIZE, caldb_ydist_final_padded,PIX_PADD_SIZE,PIX_PADD_SIZE);
    if (status) {
        LOG(ERROR) << "***Padding Unsuccessfull for y-distortion***" ;
           continue;
    }
 
 //  for optical distortion correction
    status = Applypadding(X_optical_distortion,CALDB_DIST_SIZE,CALDB_DIST_SIZE,caldb_xdist_final_od_padded,PIX_PADD_SIZE,PIX_PADD_SIZE);
    if (status) {
        LOG(ERROR) << "***Padding unsuccessfull for x-distortion***" ;
           continue;
    }
   
   status =Applypadding(Y_optical_distortion,CALDB_DIST_SIZE,CALDB_DIST_SIZE, caldb_ydist_final_od_padded,PIX_PADD_SIZE,PIX_PADD_SIZE);
    if (status) {
        LOG(ERROR) << "***Padding Unsuccessfull for y-distortion***" ;
        continue;
    }
      
   status = performSubDivisionIM (caldb_xdist_final_padded,PIX_PADD_SIZE,PIX_PADD_SIZE,caldb_xdist_final,subDivision_size,subDivision_size);
     if (status) {
        LOG(ERROR) << "***Padding unsuccessfull for x-distortion***" ;
           continue;
    }
    status =performSubDivisionIM(caldb_ydist_final_padded,PIX_PADD_SIZE,PIX_PADD_SIZE, caldb_ydist_final,subDivision_size,subDivision_size);
    if (status) {
        LOG(ERROR) << "***Padding Unsuccessfull for y-distortion***" ;
           continue;
    }
 
 //  for optical distortion correction
    status = performSubDivisionIM(caldb_xdist_final_od_padded,PIX_PADD_SIZE,PIX_PADD_SIZE,caldb_xdist_final_od,subDivision_size,subDivision_size);
    if (status) {
        LOG(ERROR) << "***Padding unsuccessfull for x-distortion***" ;
           continue;
    }
   
    status =performSubDivisionIM(caldb_ydist_final_od_padded,PIX_PADD_SIZE,PIX_PADD_SIZE, caldb_ydist_final_od,subDivision_size,subDivision_size);
    if (status) {
        LOG(ERROR) << "***Padding Unsuccessfull for y-distortion***" ;
        continue;
    }     
      

     LOG(INFO) << "Total number of frames are " << nframes  ;
  
    
       fitsfile *fptr ;
       int cnt_frm=0,cnt_frm_wmean=0,total_frames_processed=0;
    
       float *sum_weighted_sig,*sum_weighted_exp;
          sum_weighted_sig= new float[subDivision_size*subDivision_size];//Array for storing  Accumulated values of total no_ofWeigh  Signal frames
          sum_weighted_exp= new float[subDivision_size*subDivision_size];//Array for storing Accumulated values of total no_ofweigh Exposure frames.
          for (int i=0;i<subDivision_size*subDivision_size;i++) 
          {
              sum_weighted_sig[i]=0.0f;
              sum_weighted_exp[i]=0.0f;
          }
          refine_Winsize=refine_Winsize*subDivision_size/padding_dim;//setting refined window size for star finding for registration and averaging
          centroid_Winsize=centroid_Winsize*subDivision_size/padding_dim;//setting centroid window size for star finding for registration and averaging 
          
double t1 , t2 , theta1 , theta2 , x1 , x2 , y1 , y2 ;
vector<double> time_track;
vector<float> cx_ref , cy_ref , ci_ref ;
int  min_stars_match ;
float *avg_sigArray = new float[subDivision_size*subDivision_size];//Array for storing final weighted average signal frames
float *avg_expArray = new float[subDivision_size*subDivision_size];//Array for storing final weighted average signal frames
//initialization
for(int i=0;i<subDivision_size*subDivision_size;i++)
{
avg_sigArray[i]=0;
avg_expArray[i]=0;    
}

    float *Regavg_subSampled_Array_Sig = new float[FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG] ;//Array for storing final sub-sampled signal frame computed registering and averaging all the frames.
    float *Regavg_subSampled_Array_Exp = new float[FINALFRAMESIZE_REGAVG * FINALFRAMESIZE_REGAVG] ;//Array for storing final sub-sampled exposure frame computed registering and averaging all the frames.
    
     int Accu_fra_no = no_ofWeigh;
         //incase of total number of frames to be weighted is greater than the total available frames 
     if(nframes/Accu_fra_no<1){
            LOG(INFO)<<"***Total number of frames to be weighted are greater  than total available frames,Accumulation will be done on individual frame***";
            Accu_fra_no=1;
            
        }
     /*==========================================================================================================================*/
       /*==============================LOOP FOR TOTAL AVAILABLE FRAMES=================================================================*/
       /*==========================================================================================================================*/
      int cnt_Shift_applied=0;
//loop for number of frames.
       LOG(INFO)<<"Total number of frames are "<<nframes; 
         int cnt_shifts=0;
    for (int i = 0 ; i <nframes ; i++)
    {
        if (i==0){
            cnt_frm=0;
            cnt_frm_wmean=0;
        }
      
        sprintf (errstr , "Error at iteration number %d" , i) ;    
        sprintf (infile , "%s/%s/%s" , moduleIndir , "SignalFrames" , sigframelist[i]) ;
        
        frame_Data= (float*)malloc (xsize*ysize*sizeof(float));    initArray (frame_Data , xsize*ysize , -9999.0f) ;//allocating memory to frame_Data to store image pixels of frames
        frame_ExpData= (float*)malloc (xsize*ysize*sizeof(float)); initArray (frame_ExpData , xsize*ysize , -9999.0f) ;//allocating memory to frame_ExpData to store image pixels of  exposure frames
         //read signal frame image and store it in frame_Data.
        status = readImage (infile , 1 , frame_Data,xsize,ysize) ;
       
        //opening DataIngest Signal frame
        fits_open_file (&fptr , infile , READONLY , &status) ; printError (status , "Error in opening the input file" , infile) ;
        fits_movabs_hdu (fptr , 1 , NULL , &status) ; printError (status , "Error in  moving to the 2nd HDU of the out information file" , infile) ;
        fits_read_key (fptr , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ; printError (status , "Error in  reading the FRAMENO keyword" , infile) ;
        fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ; printError (status , "Error in  reading the FRAMETIME keyword" , infile) ;
        fits_read_key (fptr , TDOUBLE , "INT_TIME" , &integrationtime , NULL , &status) ; printError (status , "Error in  reading the INT_TIME keyword" , infile) ;
        fits_close_file (fptr , &status) ;printError (status , "Error in  closing the file" , infile) ;
       
        //copy level1 keywords to header_info vector
        status = copyAllheaderKeys (infile) ;     
        //for the Dark Subtraction process..
        
        if (darkframe_flag)
        {
             t_curr = frametime ;//stores image frame's time
             status = darkFrameComputation (darkCompute_array) ; //intermediate dark value for the  current frame according to t_curr.
            if (status)
            {
                LOG(ERROR) << "***Error in DarkFrame Computation***" ;
                return(EXIT_FAILURE);
            }
            status = darkFrameSubtraction(darkCompute_array , frame_Data,xsize,ysize) ;//performing the dark frame subtraction.
            if (status)
            {
                LOG(ERROR) << "***Error in Dark Frame Subtraction***" ;
                return(EXIT_FAILURE);
            }
         }
        
        //performing unit conversion(divide image pixels with the integration time)
        status=performUnitConversionIM (frame_Data,frame_ExpData,integrationtime,xsize,ysize);
        if(status)
        {
            LOG(INFO)<<"ERROR in unitConversion ";
            break;
        }
        if (wtd_uc == 1)//in case of unit conversion output to be written to the disk
        {
            status = writeOutputImageToDisk ("uc" , moduleoutdir_uc , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ;
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***"  ;
                 break;
            }
        }
         
        //performing badpixel filtering
        status=performCorrectionIM(frame_Data,frame_ExpData,badpixdata,xsize,ysize,integrationtime);
        if (status)
        {
            LOG(ERROR) << "***Error in  BAD pixel filtering***" ;
           break;
        }
        if (wtd_bp == 1)//in case of bad pixels output to be written to the disk
        {
            status = writeOutputImageToDisk ("bp" , moduleoutdir_bp , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ; //this is for the SignalFrame output
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***"  ;
                break;
            }
            status = writeOutputImageToDisk ("bp" , moduleoutdir_bp , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , xsize , ysize) ; //this is for the SignalFrame output
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***" ;
              break;
            }
        }
        
        //performing cosmic ray correction
        status=performCosmicRayCorrIM (frame_Data,frame_ExpData,xsize,ysize,cr_threshold);
        if (wtd_cr == 1)//in case of Cosmic ray correction  output to be written to the disk
        {
            status = writeOutputImageToDisk ("cr" , moduleoutdir_cr , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ; //this is for the SignalFrame output
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***"  ;
                break;
            }
            status = writeOutputImageToDisk ("cr" , moduleoutdir_cr , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , xsize , ysize) ; //this is for the SignalFrame output
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***"  ;
                break;
            }
        }
        
      /*process starts for the FlatField */
       
          status=performFlatFieldCorrIM (frame_Data,flatfielddata,xsize,ysize);
             if (wtd_ff == 1)//in case of flat field output to be written to the disk
            {
                status = writeOutputImageToDisk ("ff" , moduleoutdir_ff , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***"  ;
                    break;
                }
                status = writeOutputImageToDisk ("ff" , moduleoutdir_ff , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , xsize , ysize) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***"  ;
                    break;
                }
            }
       
       /*process Ends For the FlatField */
        
        temperature=-9999;//initialization of temperature.
                    for (int i=0;i<nrows_lbt;i++)//loop for number of rows of LBT file for finding the relevant temperature of  current frame based on the frame time 
        {
                             
            if (frametime>=time_lbt[i] && frametime<time_lbt[i+1])
            {
                
                temperature=(insideTemp[i]+outsideTemp[i])/2;
                break;
            }
        }
        if(temperature==-9999)//incase of no relevant temperature found for the current frame 
        {
            LOG(ERROR)<<"No record found in LBT file";
            return(EXIT_FAILURE);
        }
                for (int j = 0 ; j < nCalDBTempValues - 1 ; j++)//loop for finding the relevant factor value by comparing calculated temperature value and  caldb temperature value  
            {
                if (temperature >= cal_temp[j] && temperature < cal_temp[j + 1])   //temperature - from LBT file ,cal_temperature from caldb file    
                {
                    t1 = cal_temp[j] ;
                    t2 = cal_temp[j + 1] ;
                    if ((t2 - t1) == 0)
                    {
                        LOG(INFO) << "***Divide By zero***" ;
                        return (EXIT_FAILURE) ;
                    }
                    x1 = qe_mg_factor[j] ;
                    x2 = qe_mg_factor[j + 1] ;
                    factor = x1 + ((temperature - t1)*((x2 - x1) / (t2 - t1))) ;//calculated factor
                }
            }
            
    //    performing QE and MCP correction(by applying factor to image pixels)
     //   LOG(INFO)<<"The Factor value is "<<factor;
            status=performQEMCPcorrection(frame_Data,xsize,ysize,factor);
              if (wtd_qemcp== 1)//incase of QE and MCP correction to be written to the disk
              {
                status = writeOutputImageToDisk ("qe" , moduleoutdir_qemcp , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , xsize , ysize) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                    break;
                }
                status = writeOutputImageToDisk ("qe" , moduleoutdir_qemcp, "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , xsize , ysize) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                   break;
                }
             }
       
     //loop for initializing frame_Data_Padded and frame_ExpoData_padded array (i.e 600*600)
        for (int i=0;i<padding_dim*padding_dim;i++)
        {
            frame_Data_Padded[i]=-9999;
            frame_ExpoData_padded[i]=-9999;
        }
        
        //performing pixel padding
       status=Applypadding (frame_Data ,xsize,ysize, frame_Data_Padded , padding_dim , padding_dim) ;
        if (status)
            {
                LOG (ERROR) << "Error in applying padding" ;
               return(EXIT_FAILURE);
            }
       status=Applypadding (frame_ExpData ,xsize,ysize, frame_ExpoData_padded , padding_dim , padding_dim) ;
        if (status)
            {
                LOG (ERROR) << "Error in applying padding" ;
               return(EXIT_FAILURE);
            }
       
//       free(frame_Data); //releasing the memory
//       free(frame_ExpData); //releasing the memory
       delete[]  frame_Data,frame_ExpData;
       frame_Data= (float*)malloc (padding_dim*padding_dim*sizeof(float));
       initArray (frame_Data,padding_dim*padding_dim,-9999.0f); //allocating memory to frame_Data
       frame_ExpData= (float*)malloc (padding_dim*padding_dim*sizeof(float));
       initArray (frame_ExpData,padding_dim*padding_dim,-9999.0f); //allocating memory to frame_Data
       
     
        for (int pix = 0 ; pix < padding_dim * padding_dim ; pix++)//loop for assigning the values of arrays to frame_Data and frame_ExpData
        {
            frame_Data[pix] = frame_Data_Padded[pix] ;
            frame_ExpData[pix] = frame_ExpoData_padded[pix] ;
        }
          
        if (wtd_pp)//incase of  pixel padding to be written to the disk
        {
            status = writeOutputImageToDisk ("pp" , moduleoutdir_pp , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , padding_dim , padding_dim) ;
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                break;
            }
            status = writeOutputImageToDisk ("pp" , moduleoutdir_pp , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , padding_dim , padding_dim) ;
            if (status)
            {
                LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                break;
            }
        }

        if (subdivisionFlag)
        {            
           initArray (frame_Data_subdivided,subDivision_size*subDivision_size,-9999.0f);
           initArray (frame_ExpData_subdivided,subDivision_size*subDivision_size,-9999.0f);
          
            for (int p = 0 ; p < padding_dim * padding_dim ; p++) {
                if(frame_Data[p]!=INVALID_PIX_VALUE){
                frame_Data[p] = frame_Data[p] / division_fact ;
                }
            }
         
           status=performSubDivisionIM (frame_Data,padding_dim,padding_dim,frame_Data_subdivided,subDivision_size,subDivision_size);
           if(status)
           {LOG(ERROR)<<"Error in subdivision module"<<endl;
           return(EXIT_FAILURE);
           }
           status=performSubDivisionIM (frame_ExpData,padding_dim,padding_dim,frame_ExpData_subdivided,subDivision_size,subDivision_size);
           if(status) {
               LOG(ERROR)<<"Error in subdivision module"<<endl;
                return(EXIT_FAILURE);
           }
           
//             free(frame_Data);
//             free(frame_ExpData);
             delete[]  frame_Data,frame_ExpData;
             frame_Data= (float*)malloc (subDivision_size*subDivision_size*sizeof(float));
             frame_ExpData= (float*)malloc (subDivision_size*subDivision_size*sizeof(float));

             status=initArray (frame_Data,subDivision_size*subDivision_size,-9999.0f);
             status=initArray (frame_ExpData,subDivision_size*subDivision_size,-9999.0f);
           
             for (int pix = 0 ; pix < subDivision_size * subDivision_size ; pix++)
             {
                frame_Data[pix] = frame_Data_subdivided[pix] ;
                frame_ExpData[pix] = frame_ExpData_subdivided[pix] ;
             }
             
            if (wtd_sd)//incase of sub division  to be written to the disk
            {
                status = writeOutputImageToDisk ("sd" , moduleoutdir_sd , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                    break;
                }
                status = writeOutputImageToDisk ("sd" , moduleoutdir_sd , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                   break;
                }
            }
          }
       float *frm_tempdata_sig= new float[subDivision_size*subDivision_size];
       float *frm_tempdata_exp= new float[subDivision_size*subDivision_size];
       
       for (int i=0;i<subDivision_size*subDivision_size;i++){
           frm_tempdata_sig[i]=INVALID_PIX_VALUE;
           frm_tempdata_exp[i]=INVALID_PIX_VALUE;
       }
 
       /*for detector distortion correction*/
       status=performDistortionCorr (frame_Data, frm_tempdata_sig,caldb_xdist_final,caldb_ydist_final,subDivision_size,subDivision_size);
       if(status)
       {
           LOG(ERROR)<<"Error in performing  distortion correction for signal frame "<<endl;
           break;     
       }
          status=performDistortionCorr (frame_ExpData, frm_tempdata_exp,caldb_xdist_final,caldb_ydist_final,subDivision_size,subDivision_size);
       if(status)
       {
           LOG(ERROR)<<"Error in performing  distortion correction  for exposure frame"<<endl;
           break;        
       }
       status=performDistortionCorr (frame_Data, frm_tempdata_sig,caldb_xdist_final_od,caldb_ydist_final_od,subDivision_size,subDivision_size);
       if(status)
       {
           LOG(ERROR)<<"Error in performing  distortion correction for signal frame "<<endl;
           break;     
       }
          status=performDistortionCorr (frame_ExpData, frm_tempdata_exp,caldb_xdist_final_od,caldb_ydist_final_od,subDivision_size,subDivision_size);
       if(status)
       {
           LOG(ERROR)<<"Error in performing  distortion correction  for exposure frame"<<endl;
           break;        
       }
          for(int i=0;i<subDivision_size*subDivision_size;i++)
          {        
          frame_Data[i]=frm_tempdata_sig[i];
          frame_ExpData[i]=frm_tempdata_exp[i];
          }          
        
         if (wtd_dd)//incase of  detect distortion  to be written to the disk
            {
                status = writeOutputImageToDisk ("dd" , moduleoutdir_dd , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                    break;
                }
                status = writeOutputImageToDisk ("dd" , moduleoutdir_dd , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                    break;
                }
            }  
           for (int i=0;i<subDivision_size*subDivision_size;i++)
           {
           frm_tempdata_sig[i]=INVALID_PIX_VALUE;
           frm_tempdata_exp[i]=INVALID_PIX_VALUE;
           }
          status=performDistortionCorr (frame_Data, frm_tempdata_sig,caldb_xdist_final_od,caldb_ydist_final_od,subDivision_size,subDivision_size);
       if(status)
       {
           LOG(ERROR)<<"Error in performing  distortion correction for signal frame "<<endl;
           break;        
       }
          status=performDistortionCorr (frame_ExpData, frm_tempdata_exp,caldb_xdist_final_od,caldb_ydist_final_od,subDivision_size,subDivision_size);
       if(status)
       {
           LOG(ERROR)<<"Error in performing  distortion correction  for exposure frame"<<endl;
           break;           
       }
          for(int i=0;i<subDivision_size*subDivision_size;i++)
          {
             
          frame_Data[i]=frm_tempdata_sig[i];
          frame_ExpData[i]=frm_tempdata_exp[i];
          }          
         free(frm_tempdata_sig);//releasing memory
         free(frm_tempdata_exp);//releasing memory
         if (wtd_od)//incase of  optical distortion to be written to the disk.
            {
                status = writeOutputImageToDisk ("od" , moduleoutdir_od , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                    break;
                }
                status = writeOutputImageToDisk ("od" , moduleoutdir_od , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                   break;
                }
            }
//         
         
          //Appyig Shift and Rotation
          long  temp1 = 0 ;
         for (int i = 0 ; i < subDivision_size ; i++)
        {
            for (int j = 0 ; j < subDivision_size ; j++)
            {
                tempSigarr[j][i] = -9999 ;//initialize
                tempExparr[j][i] = -9999 ;//initialize
            }
        }
          //copy values of frame_Data and frame_ExpData to the SubSignalArray and subExposureArray
        for (int i = 0 ; i < subDivision_size ; i++)
        {
            for (int j = 0 ; j < subDivision_size ; j++)
           {
                subExposureArray[j][i] =frame_ExpData[temp1] ;
                subSignalArray[j][i] = frame_Data[temp1] ;
                temp1++ ;
           }
        }
          
           mult_fact=subDivision_size/600; 
          
        double ctheta , stheta ;
        double temp_theta_snr,temp_x_snr,temp_y_snr;
        //calculating current value of deltas.
      
//        temp_x_snr=0.0f;
//        temp_y_snr=0.0f;
//        temp_theta_snr=0.0f;
          temp_x_snr=INVALID_PIX_VALUE;
        temp_y_snr=INVALID_PIX_VALUE;
        temp_theta_snr=INVALID_PIX_VALUE;
        //cnt_Shift_applied=0;
      
        //loop for calculating shifts and rotation for current frame based on the frame time from RAS file 
          for (int index2 = 0 ; index2 < no_of_records ; index2++)
        {
//               temp_theta_snr=INVALID_PIX_VALUE;
//                temp_x_snr=INVALID_PIX_VALUE;
//                temp_y_snr=INVALID_PIX_VALUE;
//                 temp_theta_snr=0;
//                temp_x_snr=0;
//                temp_y_snr=0;
            if (frametime>= time_ras[index2] && frametime<=time_ras[index2 + 1])
            {
                cnt_shifts++;
              //  LOG(INFO)<<"Time RAS"<<time_ras[index2];exit(1);
                t1 = time_ras[index2] ;
                t2 = time_ras[index2 + 1] ;
                theta1 = delta_theta[index2] ;
                theta2 = delta_theta[index2 + 1] ;
                x1 = delta_x[index2]*mult_fact ;
                x2 = delta_x[index2 + 1]*mult_fact ;
                y1 = delta_y[index2] *mult_fact;
                y2 = delta_y[index2 + 1]*mult_fact ;
                
                temp_theta_snr= theta1 + (frametime - t1)*(theta2 - theta1) / (t2 - t1) ;
                temp_x_snr= x1 + (frametime - t1)*(x2 - x1) / (t2 - t1) ;
                temp_y_snr = y1 + (frametime - t1)*(y2 - y1) / (t2 - t1) ;
               //LOG(INFO)<<"Shifts ::"<<temp_theta_snr<<" "<<temp_x_snr<<" "<<temp_y_snr;
                cnt_Shift_applied++;  
                break;
           }
        }
       // exit(1);
        
         double index_i = 0.0 , index_j = 0.0 ;
        int i1 , j1;
        int mid=subDivision_size/2;
       //added
        
        if(temp_x_snr==INVALID_PIX_VALUE && temp_y_snr==INVALID_PIX_VALUE && temp_theta_snr==INVALID_PIX_VALUE) 
        {
           
           if (i+1==nframes ) goto labelwm;
            else 
            continue;
        }
//            if(temp_x_snr==0 && temp_y_snr==0 && temp_theta_snr==0) 
//        {
//           
//           if (i+1==nframes ) goto labelwm;
//            else 
//            continue;
//        }
        
        //till this 
        //applying the deltas 
        ctheta = cos ((double)(-temp_theta_snr * M_PI / 180)) ;
        stheta = sin ((double)(-temp_theta_snr * M_PI / 180) );
        

       
        //Applying reverse shifts and rotation on current frame 
        for (int i = 0 ; i < subDivision_size ; i++)
        {
            i1 = i - mid ;
            for (int j = 0 ; j < subDivision_size ; j++)
            {
                
               // if(subSignalArray[i][j]!=INVALID_PIX_VALUE){
                      j1 = j - mid ;
                /*applying new_delta_theta[index] degree rotation to find out new pixel indexes */
                index_i =   ((i1 * ctheta) - (j1 * stheta)) + mid - temp_x_snr ; //new index x
                index_j =   ((i1 * stheta) + (j1 * ctheta)) + mid - temp_y_snr ; //new index y
                //ROUNDING:Need not to be corrected.as it is applied on full number.
                round (index_i) ;
                round (index_j) ;
                if (index_i > -1 && index_i < subDivision_size && index_j > -1 && index_j < subDivision_size)
                {
                    /*appling correction i.e xshift,yshift ,rotation.*/
                   
                    tempSigarr[(int) index_i][(int) index_j] = subSignalArray[i][j] ;
                    tempExparr[(int) index_i][(int) index_j] = subExposureArray[i][j] ;
               
                }
             }         
        }
        //copy values of Array tempSigarr and tempExparr to the frame_Data and frame_ExpData .
        temp1 = 0 ;
        for (int i = 0 ; i < subDivision_size ; i++)
        {
            for (int j = 0 ; j < subDivision_size ; j++)
            {
                frame_Data[temp1] = tempSigarr[j][i] ;
                frame_ExpData[temp1] = tempExparr[j][i] ;
                temp1++ ;
            }
        }
        
        if(wtd_snr==1)//incase of  shift and rotate to be written on the disk
        {
       
            status = writeOutputImageToDisk ("snr" , moduleoutdir_snr , "SignalFrames" , "sig" , frame_Data , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                    break;
                }
                status = writeOutputImageToDisk ("snr" , moduleoutdir_snr , "ExposureFrames" , "exp" , frame_ExpData , nameprefix , frametime , frameno , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                    break;
                }
                
        }
        
        //performing Weighted mean of total  no_ofWeight frames.
   
        cnt_frm++;
      for(int i=0;i<subDivision_size*subDivision_size;i++)
          {
            if (frame_Data[i]!=INVALID_PIX_VALUE && frame_ExpData[i]!=INVALID_PIX_VALUE)
            {
                if(cnt_frm==1){
                     sum_weighted_sig[i]=frame_Data[i];
                     sum_weighted_exp[i]=frame_ExpData[i];
                }
                else
                {                 
                     sum_weighted_sig[i]=(sum_weighted_sig[i]+frame_Data[i]*frame_ExpData[i]);///2;
                     sum_weighted_exp[i]=sum_weighted_exp[i]+frame_ExpData[i];
                }
            }
           }
        
//          for(int i=0;i<subDivision_size*subDivision_size;i++)
//          {
//              if (sum_weighted_exp[i] != 0.0 && sum_weighted_sig[i]!=INVALID_PIX_VALUE )
//              sum_weighted_sig[i]=sum_weighted_sig[i]/sum_weighted_exp[i];
//          }
       
    
      time_track.push_back (frametime);//for track of time 
 
      
    //condition will be satisfiy every no_ofWeigh frames 
      
      
      if(cnt_frm%Accu_fra_no ==0   || ( i+1==nframes) )
      {
           for(int i=0;i<subDivision_size*subDivision_size;i++)
          {
              if (sum_weighted_exp[i] != 0.0 && sum_weighted_sig[i]!=INVALID_PIX_VALUE )
              sum_weighted_sig[i]=sum_weighted_sig[i]/sum_weighted_exp[i];
          }
         // LOG(INFO)<<i;
         // LOG(INFO)<<"The Accumulated Frames are "<<Accu_fra_no<<" "<<cnt_frm;exit(1);
               sd_mul_factor=sd_multi_factor_default;
          cnt_frm_wmean++;
          LOG(INFO)<<"\033[1;34mProcess  started  for  frame "<<cnt_frm_wmean<<"\033[0m"<<endl;
          mid_time=(time_track[0]+time_track[time_track.size ()-1])/2;//calculating  mid time for weighted mean frame.
          cnt_frm=0;//again initialization for next no_of weight
          
          if(wtd_wm==1)//incase of  weighted mean  to be written on disk.
          {
               status = writeOutputImageToDisk ("wm" , moduleoutdir_wm , "SignalFrames" , "sig" , sum_weighted_sig , nameprefix , mid_time, cnt_frm_wmean , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                   break;
                }
                status = writeOutputImageToDisk ("wm" , moduleoutdir_wm , "ExposureFrames" , "exp" , sum_weighted_exp , nameprefix , mid_time , cnt_frm_wmean , subDivision_size , subDivision_size) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                    break;
                }
              
          }          
          
          
          //Registration and averaging started.
          if(cnt_frm_wmean==1)   //incase of reference frame  for registration and averaging
          {
                
                    status = findStar_algo1 (sum_weighted_sig) ; // for first frame taken as  reference
                    if (status)
                    {
                       // LOG (ERROR) << endl << "***Error in finding star algorithm 1 for frame  " << infile << "  ***" << endl ;
                        break;
                    }
                    //LOG(INFO)<<Cx.size ()<<endl;
                    for (int i = 0 ; i < Cx.size () ; i ++)
                    {
                        cx_ref.push_back (Cx[i]) ;
                        cy_ref.push_back (Cy[i]) ;
                        ci_ref.push_back (Ci[i]) ;
                     
                    }
                    ofstream fout("Ref.txt");
                    for (int i=0;i<Cx.size ();i++) fout<<Cx[i]<<" "<<Cy[i]<<" "<<Ci[i]<<endl;
                             
                    fout.close ();
                    Cx.clear () ;
                    Cy.clear ()  ;
                    Ci.clear () ;
             
                    
                            
          }
          else
          { //for remaining frame
        
            int cnt = 0 ;
            int matching_points = 0 ;
            status = findStar_algo1 (sum_weighted_sig) ;//for finding stars
            if (status)
            {
                LOG (INFO) << "***Error in finding star algorithm 1 for frame  " << infile << "  ***" ;
               break;
            }
            Rx.clear () ;
            Ry.clear () ;
            Rval.clear () ;
            Fx.clear () ;
            Fy.clear () ;
            Fval.clear () ;
            char name[100];
            sprintf(name,"File_%d",cnt_frm_wmean);
             ofstream fout(name);
            for (int i=0;i<Cx.size ();i++) fout<<Cx[i]<<" "<<Cy[i]<<" "<<Ci[i]<<endl;;
                             
                    fout.close ();
           
            if (Cx.size () >= cx_ref.size ())
            {
                matching_points = cx_ref.size () ;
            }
            else
            {
                matching_points = Cx.size () ;
            }

            vector<float> x_ref_arr , y_ref_arr , x_arr , y_arr , temp_x_arr , temp_y_arr ;

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

                cnt = matchStars (cx_ref.size () , Cx.size () , mult_fact , cx_ref.data () , cy_ref.data () , Cx.data () , Cy.data () , x_ref_arr , y_ref_arr , x_arr ,
                        y_arr , temp_x_arr , temp_y_arr) ;

                diff_Dist = diff_Dist * 2 ;
            }
            
            
            diff_Dist = 1 ;
       
            double x_dx = 0.0 , y_dy = 0.0 , theta_dt = 0.0 ;

            /* option_LeastSquare parameter decide which algorithm to use for finding the Drifts between  points 
           these are  the different techniques for finding shifts between two point of two frames.
             */
            
            //finding Shift and theta between current frame and reference frame
            status = findShiftsNtheta (cnt , x_ref_arr , y_ref_arr , x_arr , y_arr , temp_x_arr , temp_y_arr ,flag_thetaComp, x_dx , y_dy , theta_dt) ;
            if (status)
            {
                LOG (INFO) << "Error in finding shifts n theta " << endl ;
               break;
            }
            LOG(INFO)<<x_dx<<" "<<y_dy<<" "<<theta_dt;
            x_ref_arr.clear ();y_ref_arr.clear ();x_arr.clear ();y_arr.clear ();temp_x_arr.clear ();temp_y_arr.clear ();
            float ctheta = 0.0f , stheta = 0.0f ;
            ctheta = cos (-1.0 * theta_dt) ;
            stheta = sin (-1.0 * theta_dt) ;
            int cnt_loop = 0 ;
            LOG (INFO) << "Loop Started for assigning the correction to the frames.."  ;

            for (int i = 0 ; i < subDivision_size * subDivision_size ; i ++)
            {
                frame_Data[i] = 0 ;
                frame_ExpData[i] =0 ;
            }

            //loop for applying Shifts to get registered image
            for (int i = 0 ; i < subDivision_size ; i ++)
            {
                int x_index = i - subDivision_size / 2 ;
                for (int j = 0 ; j < subDivision_size ; j ++)
                {
                     if(sum_weighted_sig[j*subDivision_size+i]!=INVALID_PIX_VALUE  && sum_weighted_exp[j*subDivision_size+i]!=INVALID_PIX_VALUE){
                    
                    int y_index = j - subDivision_size / 2 ;
                    //ROUNDING:changed !! before it is applied on only to subpart of equation.
                    float new_index_x = round ((x_index) * ctheta - (y_index) * stheta + subDivision_size / 2 - x_dx* mult_fact) ; //new index x
                    float new_index_y = round ((x_index) * stheta + (y_index) * ctheta + subDivision_size/ 2 - y_dy* mult_fact) ; //new index y
                   
                    if (round (new_index_x) < subDivision_size && round (new_index_x) > 0 && round (new_index_y) > 0 && round (new_index_y) < subDivision_size)
                    {
                        cnt_loop ++ ;
                        
                        //ROUNDING:Correct!! As rounding is applied on both the direction i.e X and Y.
                        frame_Data[(int) (round (new_index_y) * subDivision_size + round (new_index_x))] = sum_weighted_sig[j * subDivision_size + i] ;
                        frame_ExpData[(int) (round (new_index_y) * subDivision_size + round (new_index_x))] = sum_weighted_exp[j * subDivision_size + i] ;
                    }
                     }
                }
            }
            Cx.clear () ;
            Cy.clear () ;
            Ci.clear () ;
            
            for (int p = 0 ; p < subDivision_size*subDivision_size ; p ++)
            {
                        if(frame_Data[p]!=INVALID_PIX_VALUE  && frame_ExpData[p]!=INVALID_PIX_VALUE)
                        {

                            avg_sigArray[p] = avg_sigArray[p] + frame_Data[p] * frame_ExpData[p] ;
                            avg_expArray[p] = avg_expArray[p] + frame_ExpData[p] ;
                            }
            }
          }
          for (int i=0;i<subDivision_size*subDivision_size;i++) 
          {
              sum_weighted_sig[i]=0;
              sum_weighted_exp[i]=0;
          }
          time_track.clear ();
           
      }
       
     free(frame_Data);
     free(frame_ExpData);
     //delete[]  frame_Data,frame_ExpData;
   
     
      } ///loop finished for nframes.
       
   labelwm:     
     //loop for getting final registration image by dividing by exposure frame
       for (int p = 0 ; p < subDivision_size*subDivision_size ; p ++)
         {
                        if (avg_expArray[p] > 0.0 && avg_sigArray[p]!=INVALID_PIX_VALUE )
                          avg_sigArray[p] = (avg_sigArray[p] / avg_expArray[p]) ;
         }
    
    for(int i=0;i<FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG;i++)
    {    
        Regavg_subSampled_Array_Sig[i]=INVALID_PIX_VALUE;
        Regavg_subSampled_Array_Exp[i]=INVALID_PIX_VALUE;        
    }       
   //for getting sub-sampled signal  image                                                                   
    status = ApplySubSampling (avg_sigArray , subDivision_size , subDivision_size , Regavg_subSampled_Array_Sig , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG (ERROR) << "ERROR in sub sampling " << endl ;
       return(EXIT_FAILURE);
    }
    //for getting sub-sampled exposure image
    status = ApplySubSampling (avg_expArray , subDivision_size , subDivision_size, Regavg_subSampled_Array_Exp , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
    if (status)
    {
        LOG (ERROR) << "ERROR in sub sampling " << endl ;
        return(EXIT_FAILURE);
    }   
    
                status = writeOutputImageToDisk ("rg" , moduleoutdir_ravg , "" , "sig" , Regavg_subSampled_Array_Sig, nameprefix , mid_time, cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                   return(EXIT_FAILURE);
                }
                status = writeOutputImageToDisk ("rg" , moduleoutdir_ravg , "" , "exp" , Regavg_subSampled_Array_Exp , nameprefix , mid_time , cnt_frm_wmean , FINALFRAMESIZE_REGAVG , FINALFRAMESIZE_REGAVG) ;
                if (status)
                {
                    LOG(ERROR) << "***Writing to Disk Fails***" << endl ;
                  return(EXIT_FAILURE);
                }
         
  
         delete[] flatfielddata;
         string filenamesnr;
       //performing full frame astrometry
       uvtFullFrameAst ast_obj ;
       ast_obj.read ((char*)moduleoutdir_ravg, (char*)caldbindir.c_str (),(char*)filenamesnr.c_str(),(char*)outputdir ,(char*)dirobj.attfile.c_str (),(char*)att_timecol ,(char*)att_qcol,(char*)caldbindir.c_str (),sd_multi_factor_default,minimum_No_of_Stars,refine_Winsize,centroid_Winsize,databasename,search_algo_ctlg,len_a,len_b,rad_search,clobber , history,1) ;
       ast_obj.display () ;
       status = ast_obj.uvtFullFrmAstProcess () ;
       if (status)
       {
            LOG(ERROR) << endl << "Error in full frame astrometry  module" << endl ;
            return(EXIT_FAILURE);
       }
     
       delete[] delta_theta,delta_x,delta_y;
     
      }
       
   return (EXIT_SUCCESS) ;
}

int uvtLevel2IM::setDirectoryStructure (char *Dir , const char *subdir)
{
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s/" , Dir , subdir) ;
    //LOG(INFO)<<dir<<endl;
    string cmd ;
    system (cmd.c_str ()) ;
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
    system (cmd.c_str ()) ;

    return (EXIT_SUCCESS) ;

}

int uvtLevel2IM::writeOutputImageToDisk (char *id , char *outDir , char *dir , char *subscript , float *Array , char *namepre , double ftime , unsigned short fno , int sizex , int sizey)
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
    fits_create_file (&fout , outfile , &status) ;
    printError(status , "Error in creating the output Signal File" , outfile) ;
    fits_create_img (fout , bitpix , naxis , naxes , &status) ;
    printError(status , "Error in Creating the image for Signal Fie" , outfile) ;
    fits_write_pix (fout , TFLOAT , fpixel , sizex*sizey , Array , &status) ;
    printError(status , "***Error in writing the pixels to output***" , outfile) ;
    if(strcmp (id,"rg")==0){
     fits_write_key (fout , TSTRING , "NAMEPRFX" , namepre , NULL , &status) ;
    printError (status , "***Error in updating the nameprefix keyword***") ;
     fits_write_key (fout , TINT, "XSIZE" , &sizex , NULL , &status) ;
    printError (status , "***Error in writing xsize keyword***") ;
    fits_write_key (fout , TINT , "YSIZE" , &sizey , NULL , &status) ;
    printError (status , "***Error in writing the ysize keyword***") ;
    }
  
    fits_get_num_hdus (fout , &numhdu , &status) ;
    char temp[1000] ;
    writeCommonKeywords (fout , outDir) ;
    for (int i = 0 ; i < numhdu ; i++)
    {
        for (int p = 0 ; p < header_info.size () ; p++)
        {
            sprintf (temp , "%s" , header_info[p].c_str ()) ;
            fits_write_record (fout , (char*) &temp , &status) ;
        }
    }
    fits_close_file (fout , &status) ;
    printError (status , "Error in closing the  output Signal fits file" , outfile) ;
    return (EXIT_SUCCESS) ;
}

int uvtLevel2IM::copyAllheaderKeys (char* infile)
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
    for (int i = 1 ; i <= keyexist ; i++)
    {
        fits_read_record (fin , i , record , &status) ;
        if (status)
        {
            LOG(INFO) << endl << "***Error in reading record number " << i << " in file " << infile << "***" ;
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
    
    return (EXIT_SUCCESS) ;
}

double uvtLevel2IM::readDarkFrame (char * path , float *Array)
{

    fitsfile *fptr ;
    int status = 0 ;
    long fpixel[2] ;
    for (int i = 0 ; i < xsize * ysize ; i++) Array[i] = (float) 0.0 ;
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

int uvtLevel2IM::darkFrameComputation (float *outArray)
{
    for (int i = 0 ; i < xsize * ysize ; i++)
    {
        outArray[i] = 0.0 ;
        float d1 = darkFramestart_data[i] ;
        float d2 = darkFrameend_data[i] ;
        double t1 = t_darkframestart ;
        double t2 = t_darkframeend ;
        //    LOG(INFO)<<"D1::"<<d1<<endl;
        //    LOG(INFO)<<"D2::"<<d2<<endl;
        //    LOG(INFO)<<"t1::"<<t1<<endl;
        //    LOG(INFO)<<"t2::"<<t2<<endl;
           float d ;
        if ((t2- t1) == 0)
        {
            d=d1;
            //LOG(ERROR) << "***Divide By zero Error***" << endl ;
            //return (EXIT_FAILURE) ;
        }
        else
        {       
          d = d1 + ((d2 - d1) *((t_curr - t1) / (t2 - t1))) ;
        //LOG(INFO)<<d<<endl;
        }
        outArray[i] = (float) d ;
        // outArray[i] = d1;

    }
    return (EXIT_SUCCESS) ;

}




//int uvtImRa_commonArray::performCorrection(float *frmsigdata, float *frmexpdata, float *badpixarry,int sizex,int sizey,double  intgrntime)
//{
//    
//     for (int pixno = 0 ; pixno <sizex*sizey  ; pixno++)
//     { 
//         if(frmsigdata[pixno]!=INVALID_PIX_VALUE){
//            frmsigdata[pixno] = frmsigdata[pixno] * badpixarry[pixno] ;
//         }
//     }
//        for (int pix = 0 ; pix < sizex * sizey; pix++)
//        {
//            frmexpdata[pix] = 0.0f ;
//            frmexpdata[pix] = badpixarry[pix] * intgrntime ;
//        }
//    return(EXIT_SUCCESS);
//}
//int uvtImRa_commonArray::performCosmicRayCorr(float *frmsigdata, float *frmexpdata,int sizex,int sizey,float threshold_cr)
//{
//    vector<int> xpix,ypix;
//    int cnt_cosmicAffected = 0 ;
//     for (int pixno = 0 ; pixno <sizex*sizey ; pixno++)
//        {
//         if(frmsigdata[pixno]!=INVALID_PIX_VALUE){
//           
//             if (frmsigdata[pixno] > threshold_cr)
//            {
//                 xpix.push_back (pixno % sizey+1) ;
//                 ypix.push_back (pixno / sizex+1) ;
//                //frmsigdata[pixno] = frmexpdata[pixno] = 0.0 ;
//            }
//         }
//            
//
//        }
//    
//      for (int i = 0 ; i < xpix.size () ; i ++)
//        {
//            cnt_cosmicAffected = 0 ;
//
//            for (int j = xpix[i] - 3/2 ; j <= xpix[i] +3/2 ; j ++)
//            {
//                for (int k = ypix[i] -3/2 ; k <= ypix[i] +3/2 ; k ++)
//                {
//
//                    if (j < sizex && j > 0 && k < sizey && k > 0)
//                    {
//                        if (frmsigdata[k * sizey + j] > threshold_cr)
//                        {
//                            cnt_cosmicAffected ++ ;
//                        }
//
//                    }
//
//                }
//
//            }
//             if (cnt_cosmicAffected == 1)
//            {
//                frmexpdata[ypix[i] * sizex + xpix[i]] = frmsigdata[ypix[i] * sizey + xpix[i]] = -9999 ; //cr effected
//                //cr_effected.push_back (frameno) ;
//               // x_crFailed.push_back (X_pixel[i]) ;
//                //y_crFailed.push_back (Y_pixel[i]) ;
//
//            }
//      }
//        
//   
//    
//    return(EXIT_SUCCESS);
//}
//int uvtImRa_commonArray::performSubDivision(float *frmsigdata,int sizex,int sizey,float *subdivideddata,int size_subdivx,int size_subdivy)
//{
//    
//      if (sizex == 0 || sizey == 0)
//    {
//        LOG(INFO) << "Divide by Zero" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//    int xfactor = size_subdivx / sizex ;
//    int yfactor = size_subdivy / sizey ;
//    for (int i = 0 ; i < size_subdivx ; i++)
//    {
//        for (int j = 0 ; j < size_subdivx ; j++)
//        {
//            subdivideddata[i * size_subdivx + j] = frmsigdata[(i / yfactor) * sizey + (j / xfactor)] ;
//        }
//    }
//    return(EXIT_SUCCESS);
//}



int uvtLevel2IM::writeOutputTblToDisk (char *id , char *outDir , char *dir , char *subscript  , char *namepre , double ftime , unsigned short fno ,char **type1,char**tform1,int tfields1,char *extname,vector<float> &X ,vector<float> &Y,vector<float> &val)
{
    char out_file[FLEN_FILENAME];
    sprintf (out_file , "%s/%s/%s_t%f_f%d_%s_%s.fits" , outDir , dir , namepre , ftime , fno , id,subscript ) ;
    if(strcmp (id,"rf")==0){
        ref_frame_module_filename_track.push_back (basename(out_file));
    }
                fitsfile *fout1 ;
                int status=0;
                fits_create_file (&fout1 , out_file , &status) ;
                printError (status , "Error creating the output file ",out_file) ;
                fits_create_tbl (fout1 , BINARY_TBL , 0 , tfields1 , type1,tform1 ,NULL,extname , &status) ;
                printError (status , "Error in creating the table",out_file) ;
               fits_update_key (fout1 , TDOUBLE , "FRMTIME" , &ftime , " Average Frame time" , &status) ;
               printError (status , "Error in updating the key value of the FRMTIME",out_file) ;
                fits_write_col (fout1 , TFLOAT , 1 , 1 , 1 , X.size () , X.data () , &status) ;
                printError (status , "Error in writing  column",out_file) ;
                fits_write_col (fout1 , TFLOAT , 2 , 1 , 1 , Y.size () , Y.data () , &status) ;
                printError (status , "Error in writing the column",out_file) ;
                fits_write_col (fout1 , TFLOAT , 3 , 1 , 1 , val.size () , val.data () , &status) ;
                printError (status , "Error in writing the column",out_file) ;
                fits_close_file(fout1,&status);                
                return (EXIT_SUCCESS) ;

}
int uvtLevel2IM::performDistortionCorr(float *frmdata, float *frmtempdata,float *Xdistr,float *Ydistr,int sizex,int sizey)
{
    
    for (int i = 0; i < sizex*sizey; i++)
   {
        if(frmdata[i]!=INVALID_PIX_VALUE)
        {
            //ROUNDING:Actually No need for rounding.because it is giving an integer value.e.g 5/512=0.0.if it  casted with  integer  than (int)0.0=0.  
           int rown=(int)round(i/sizex);
           int coln=(int)round(i%sizey);

           rown=rown-Ydistr[rown*sizex+coln];
           coln=coln-Xdistr[rown*sizey+coln];
        
          
           if((int)(round(rown)*sizex+coln)<(sizex*sizey) && (int)(round(rown)*sizex+coln)>0)
           {
               // frmtempdata[(int)(round(rown)*sizex+round(coln))]=INVALID_PIX_VALUE;
               //ROUNDING:No need to correct as rounding is applied on both the direction.
                frmtempdata[(int)round(round(rown)*sizex+round(coln))]=frmdata[i];
           }
        }
     }
//           for (int i = 0; i < sizex*sizey; i++)
//          {
//            
//          int rown=(int)round(i/sizex);
//           int coln=(int)round(i%sizey);
//
//           rown=rown-Ydistr[rown*sizex+coln];
//           coln=coln-Xdistr[rown*sizey+coln];
//        
//          
//           if((int)(round(rown)*sizex+coln)<(sizex*sizey) && (int)(round(rown)*sizex+coln)>0)
//           {
//                frmtempdata[(int)(round(round(rown)*sizex+round(coln)))]=INVALID_PIX_VALUE;
//                frmtempdata[(int)(round(round(rown)*sizex+round(coln)))]=frmdata[i];
//           }
//
//           }
   
         return(EXIT_SUCCESS);
}


int uvtLevel2IM::takeDarkinfo ()
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

    if (dark_framenames.size () <2)
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
//int uvtImRa_commonArray::readQEMCPFile (){
//    
//    
//    
//}
// int uvtImRa_commonArray::readRASfile(char *ras_file,float *time_ras,float *roll_ras,float *pitch_ras,float *yaw_ras)
// {
//     int status=0;
//     int no_of_records=0;
//     fitsfile *fras ;
//    /**Reading Reletive Aspect File from RA series**/
//    LOG(INFO) << "\nReading RAS file from RA series..." << endl ;
//    fits_open_file (&fras , ras_file , READONLY , &status) ;
//    printError (status , "Error in opening the rasfile",ras_file) ;
//    fits_movabs_hdu (fras , 2 , NULL , &status) ;
//    printError (status , "Error in moving the 2nd HDU of  ras file",ras_file) ;
//    fits_get_num_rows (fras , &no_of_records , &status) ;
//    printError (status , "Error Reading the number of Rows",ras_file) ;
//
//    time = new double[no_of_records] ;
//    roll_ras = new double [no_of_records] ;
//    pitch_ras = new double [no_of_records] ;
//    yaw_ras = new double[no_of_records] ;
//    
//    fits_read_col (fras , TDOUBLE , 1 , 1 , 1 , no_of_records , NULL , time_ras , NULL , &status) ;
//    printError (status , "***Error reading  centroid x***") ;
//    fits_read_col (fras , TDOUBLE , 2 , 1 , 1 , no_of_records , NULL , yaw_ras , NULL , &status) ;
//    printError (status , "***Error writing centroid y***") ;
//    fits_read_col (fras , TDOUBLE , 3 , 1 , 1 , no_of_records , NULL , pitch_ras , NULL , &status) ;
//    printError (status , "***Error writing x-correction***") ;
//    fits_read_col (fras , TDOUBLE , 4 , 1 , 1 , no_of_records , NULL , roll_ras , NULL , &status) ;
//    printError (status , "***Error writing y-correction***") ;
//    fits_close_file (fras , &status) ;
//    printError (status , "Error in closing the file",ras_file) ;
// 
//     return(EXIT_SUCCESS);
// }
 
 
 int uvtLevel2IM::findStar_algo1 (float *inputArray) //algorithm for finding the peaks
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
        LOG (ERROR) << "***SD_MULTI_FACTOR is <0***" << endl ;
        return (EXIT_FAILURE) ;
    }
    double thr = 0 ;
    float mean_Ofimage,sd_temp;;
    if (datainfo.getModeFlag () == PC)
    {
        
        sd_temp= getSD (temp_array , array_temp.size ())  ;
         mean_Ofimage=getmean (temp_array,array_temp.size ());
//    thr =  sd_temp* sd_mul_factor ;
        thr = mean_Ofimage+sd_temp* sd_mul_factor ;
    }
    else
    {
        sd_temp= getSD (inputArray , subDivision_size*subDivision_size) ;
        mean_Ofimage=getmean (inputArray,subDivision_size*subDivision_size);
//    thr =  sd_temp* sd_mul_factor ;
   thr = mean_Ofimage+sd_temp* sd_mul_factor ;;
    }
    LOG (ERROR)  << "Threshold for first cut peaks is   " << mean_Ofimage<<" + "<<sd_temp<<" X "<<sd_mul_factor<<" = "<<thr ;
 
    for (int i = 0 ; i <subDivision_size*subDivision_size ; i ++)
    {
        r = (i / subDivision_size) ;
        c = (i % subDivision_size) ;

        if (inputArray[i] > thr)
        {
            Fval.push_back (inputArray[i]) ;
            Fx.push_back (c) ; //x is for column
            Fy.push_back (r) ; //y is for row
        }
        
    }
   
    LOG (INFO) << "SIGMA  Factor::" << sd_mul_factor  ;
    LOG (INFO) << " Size of First cut Peaks  " << Fy.size ()  ;

    if (Fy.size () < minimum_No_of_Stars)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG (ERROR) << sd_mul_factor << " less than 0!!!! " ;
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

    Tx = Fx ;
    Ty = Fy ;
    Tval = Fval ;
Star_l2 star1;
 star_track.clear ();
 for(int i=0;i<Tx.size ();i++)
 {
     star1.x=Tx[i];
     star1.y=Ty[i];
     star1.intensity=Tval[i];
     star_track.push_back (star1);
     //LOG(INFO)<<Tx[i]<<" "<<Ty[i]<<" "<<Tval[i]<<endl;
 }
 
 sort (star_track.begin (),star_track.end (),compare1);
    /*refining peaks logic
    refined Window size is for the refined  peaks.
   Refined peaks are found by  making window around each of the star(i.e first cut peaks)  and  finding brightest star among that window.*/
    int start_r , end_r , start_c , end_c ;
//to be removed 
   Tx.clear ();Ty.clear ();Tval.clear ();
  //  bool flag_unique=FALSE;
   //LOG(INFO)<<star_track.size ()<<endl;
  // vector<Star> ::iterator itr =star_track.begin ();
    for (int i = 0 ; i < star_track.size () ; i ++)
    {
        start_r = star_track[i].y - refine_Winsize / 2 ;
        end_r = star_track[i].y + refine_Winsize / 2 ;
        start_c = star_track[i].x- refine_Winsize / 2 ;
        end_c = star_track[i].x + refine_Winsize / 2 ;
           if (start_r < 0) start_r = 0 ;
        if (end_r >= subDivision_size) end_r = subDivision_size - 1 ;
        if (start_c < 0) start_c = 0 ;
        if (end_c >= subDivision_size) end_c = subDivision_size - 1 ;
         for(int fcpeak=i+1;fcpeak<star_track.size ();fcpeak++)
         {
             if(star_track[fcpeak].x>start_c && star_track[fcpeak].x<end_c && star_track[fcpeak].y>start_r && star_track[fcpeak].y<end_r)
             {
                 
                 star_track.erase (star_track.begin ()+fcpeak);
                 fcpeak--;
             }
          }
        Tx.push_back (star_track[i].x);
        Ty.push_back (star_track[i].y);
        Tval.push_back (star_track[i].intensity);
    }
    
    /*refining peaks logic
    refined Window size is for the refined  peaks.
   Refined peaks are found by  making window around each of the star(i.e first cut peaks)  and  finding brightest star among that window.*/
//    int start_r , end_r , start_c , end_c ;
//
//    for (int i = 0 ; i < Fx.size () ; i ++)
//    {
//        start_r = Ty[i] - refine_Winsize / 2 ;
//        end_r = Ty[i] + refine_Winsize / 2 ;
//        start_c = Tx[i] - refine_Winsize / 2 ;
//        end_c = Tx[i] + refine_Winsize / 2 ;
//        if (start_r < 0) start_r = 0 ;
//        if (end_r >= ysize) end_r = ysize - 1 ;
//        if (start_c < 0) start_c = 0 ;
//        if (end_c >= xsize) end_c = xsize - 1 ;
//        int max = 0 ;
//        for (int k = start_r ; k <= end_r ; k ++)
//        {
//            for (int l = start_c ; l <= end_c ; l ++)
//            {
//
//                if (inputArray[k * xsize + l] > max)
//                {
//                    max = inputArray[k * xsize + l] ;
//                    Tx[i] = l ;
//                    Ty[i] = k ;
//                    Tval[i] = inputArray[k * xsize + l] ;
//                } //  end of if block 
//            } //end of l loop
//        } //end of  k  loop
//    } // end of i loop

    /*--------------Refining peaks completed----------------*/

    float *arr_refine = new float[subDivision_size*subDivision_size] ; //to store refined peaks

    for (int i = 0 ; i < subDivision_size*subDivision_size ; i ++)
        arr_refine[i] = 0 ;

    if (subDivision_size == 0 || subDivision_size == 0)
    {
        LOG (ERROR) << "***Divide by Zero***"  ;
        return (EXIT_FAILURE) ;
    }
    for (int i = 0 ; i < Ty.size () ; i ++)
        arr_refine[Ty[i] * subDivision_size + Tx[i]] = Tval[i] ; //overwriting the same place..

    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;

    for (int i = 0 ; i <subDivision_size*subDivision_size ; i ++)
    {
        if (arr_refine[i] != 0)
        {
            Rx.push_back ((i % subDivision_size)) ;
            Ry.push_back ((i / subDivision_size)) ;
            Rval.push_back (arr_refine[i]) ;
        }
    }
    LOG (INFO) << "Number of final peaks is " << Rval.size ()  ;
    

    if (Ry.size () < minimum_No_of_Stars)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG (ERROR) << sd_mul_factor << " less than 0!!!! "  ;
            return (EXIT_FAILURE) ;
        }
        goto label ;
    }
    /**method for  find Centroid within the Stars**/
    doCentroiding (Rx , Ry , centroid_Winsize , inputArray , subDivision_size , subDivision_size) ;
    delete[] arr_refine ;
    return (EXIT_SUCCESS) ;
}


void uvtLevel2IM::doCentroiding (vector<int> &X , vector<int> &Y , int centroidwindow , float *arr , int h , int w)
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

                val = arr[(Y[i] + j) * w + (X[i] + k)] ;
                if(val==INVALID_PIX_VALUE){
                    continue;
                }
                //  val=arr[(Y[i]+j)*w+(X[i]+k)];
                //   LOG(INFO)<<"The value::"<<x+j<<y+j<<arr[(Y[i]+j)*w+(X[i]+k)]<<endl;
                sum_x = sum_x + (x + k) * val ;
                sum_y = sum_y + (y + j) * val ;
                sum = sum + val ;
            }
        }
        if (sum <= 0)
        {
            LOG (INFO) << "Sum of intensites for (" << X[i] << " , " << Y[i] << ")  is <=0" ;
            LOG (INFO) << "\nDivide by zero error\n" ;
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
  
}
    
    int uvtLevel2IM::matchStars (int numrowsFirstfile , int numrowsSecfile , float divFact ,
        float *xlocFirst , float *ylocFirst , float *xlocSec , float *ylocSec , vector<float> &matchPixelXone ,
        vector<float> &matchPixelYone , vector<float> &matchPixelXtwo , vector<float> &matchPixelYtwo ,
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
                //ofptr1<<temp_x1*div_fact<<setw(20)<<temp_y1*div_fact<<setw(20)<<ints1<<setw(20)<<temp_x2*div_fact<<setw(20)<<temp_y2*div_fact<<setw(20)<<ints2<<endl;
                cnt ++ ;
                break ;

            }
        }
    }
    return cnt ;
}
    
    int uvtLevel2IM::findShiftsNtheta (int totalelements , vector<float> &Xone , vector<float> &Yone , vector<float> &Xtwo , vector<float> &Ytwo ,
        vector<float> &DiffOfX , vector<float> &DiffOfY ,bool flag_theta_computation, double &Xdx , double &Ydy , double &Theta)
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
        else  if (star_detect_algo_flag == 1)
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
    else if (star_detect_algo_flag == 2)
    {
        spMatrix A (totalelements * 2 , 4) ;
        spMatrix B (totalelements * 2 , 1) ;
        spMatrix X (4 , 1) ;
        for (int aindex = 0 ; aindex < totalelements ; aindex ++)
        {
            A (2 * aindex , 0) = Xone[aindex] - (subDivision_size / mult_fact) * 0.5 ;
            A (2 * aindex , 1) = - 1.0 * (Yone[aindex] - (subDivision_size /  mult_fact) * 0.5) ;
            A (2 * aindex , 2) = 1.0 ;
            A (2 * aindex , 3) = 0.0 ;
            A (2 * aindex + 1 , 0) = (Yone[aindex] - (subDivision_size / mult_fact) * 0.5) ;
            A (2 * aindex + 1 , 1) = Xone[aindex] - (subDivision_size /  mult_fact) * 0.5 ;
            A (2 * aindex + 1 , 2) = 0.0 ;
            A (2 * aindex + 1 , 3) = 1.0 ;
            B (2 * aindex , 0) = Xtwo[aindex] - (subDivision_size / mult_fact) * 0.5 ;
            B (2 * aindex + 1 , 0) = Ytwo[aindex] - (subDivision_size /  mult_fact) * 0.5 ;
        }
        X.ApplyLeastSquare (A , B) ;

        double theta = atan2 (X (1 , 0) , X (0 , 0)) ;

        // ofptr << theta * 180 / M_PI << setw (15) << X (0 , 0) << setw (25) << X (1 , 0) << endl ;
        Xdx = X (2 , 0) ;
        Ydy = X (3 , 0) ;
        Theta = theta ;
    }
    else if (star_detect_algo_flag == 3)
    {

        double a11 = 0.0 , a12 = 0.0 , a13 = 0.0 , a21 = 0.0 , a22 = 0.0 , a23 = 0.0 , a31 = 0.0 , a32 = 0.0 , a33 = 0.0 , b1 = 0.0 , b2 = 0.0 , b3 = 0.0 ;

        spMatrix A (3 , 3) ;
        spMatrix B (3 , 1) ;
        spMatrix X (3 , 1) ;

        for (int k = 0 ; k < totalelements ; k ++)
        {
            /* weight used for a star*/
            double w = 1.0 ;
            //     w=pow( (25.0*photonbk_per_pixel+starmag[k]), 2.0 )/
            //        (  (25.0*photonbk_per_pixel)+(0.09*starmag[k])  );


            /* row 1 */
            double y_mod = Ytwo[k]-((subDivision_size /  mult_fact) * 0.5) ;
            double x_mod = Xtwo[k]-((subDivision_size /  mult_fact) * 0.5) ;

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

int uvtLevel2IM::ApplySubSampling (float* inputarray , int in_xsize , int in_ysize , float* outputarray , int out_xsize , int out_ysize)
{
    
    if (in_xsize % out_xsize != 0)
    {
        LOG (ERROR) << "Can't sub sampling with  this input array " << endl ;
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
                LOG (ERROR) << "Array is out of bound, EXEED to " << out_xsize << " !!!" << endl ;
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
                outputarray[index_finalArray ++] = (float) (sum_win / cnt_win) ;
            }
            else
            {
                outputarray[index_finalArray ++] = INVALID_PIX_VALUE ;
            }


        }

    }

    return (EXIT_SUCCESS) ;

}
//int uvtImRa_commonArray::performFlatFieldCorr(float *frmsigdata,float *flatfieldarry,int sizex,int sizey)
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
//    LOG(INFO) <<"Reading QE MCP  Temperature vs filter file from calDB........" ;
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
//    {
//    fits_read_col (fqemcp , TFLOAT , 9 , 1 , 1 , nrows , NULL , (void*) f7 , NULL , &status) ;
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
//   // LOG(INFO)<<"Channel :" <<md<<;
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
//
int uvtLevel2IM::getHistory (vector<string> &vhistory)
{
   // char *user = getlogin () ;
    int cnt = 0 ;
    char validgtiflag_str[FLEN_FILENAME] ;
   // string str = "Module run by " + (string) user ;
    char temp[PIL_LINESIZE];
  //  vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " input Level1 tar  file = " + (string) level1indir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " caldb used= " + (string) caldbindir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Module Output directory = " + (string)level2outdir ) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " channel used= " + (string)channel) ;
    
    vhistory.push_back ((string) getSerialNo (cnt) + " drop frame flag = " + (string) convertIntToStr (dropframe)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Dark frame subtraction to be done or not = " + (string)convertIntToStr (darkframe_flag)) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " GTI flag  = " + (string)convertIntToStr(gti_flag)) ;
    //vhistory.push_back ((string) getSerialNo (cnt) + " Flat field Flag = " + (string)convertIntToStr( flatfieldFlag)) ;
    //else if (gti_flag == 1)
    //{
    vhistory.push_back ((string) getSerialNo (cnt) + " Threshold for Cosmic Ray correction = " + (string)convertIntToStr( cr_threshold)) ;
        vhistory.push_back ((string) getSerialNo (cnt) + " Padding Dimension  = " + (string) convertIntToStr(padding_dim) );
        //sprintf (validgtiflag_str , "%d" , valid_gtiflag) ;
       // vhistory.push_back ((string) getSerialNo (cnt) + " Number of frames to be Accumulated  = " + (string)convertIntToStr( no_ofFramesToAcc)) ;
        //if (all_Or_custom)
        //{
         //  vhistory.push_back ((string) getSerialNo (cnt) + " QEMCP flag= " + (string) convertIntToStr(qe_mcpFlag)) ;
        //}
       // else if (all_Or_custom == 0)
        //{
            vhistory.push_back ((string) getSerialNo (cnt) + " SubDivision Flag = " + (string) convertIntToStr (subdivisionFlag)) ;
        //}

    //}
            if(subdivisionFlag)
            vhistory.push_back ((string) getSerialNo (cnt) + " SubDivision size = " + (string) convertIntToStr(subDivision_size)) ;
            
            vhistory.push_back ((string) getSerialNo (cnt) + " AlgoFlag used for calculating finding stars = " + (string) convertIntToStr(star_detect_algo_flag)) ;
            
            if(star_detect_algo_flag==1 || star_detect_algo_flag==3 || star_detect_algo_flag==4){
                vhistory.push_back ((string) getSerialNo (cnt) + " S.D. multiplication factor value = " + (string) convertFloatToStr (sd_multi_factor_default)) ;
                  vhistory.push_back ((string) getSerialNo (cnt) + " Minimum targeted stars = " + (string)convertIntToStr(minimum_No_of_Stars)) ;
            }
            
             vhistory.push_back ((string) getSerialNo (cnt) + " Refined Window size = " + (string)convertIntToStr( refine_Winsize)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Centroid  Window size = " + (string)convertIntToStr( centroid_Winsize)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " number of frames to be weighted = " + (string) convertIntToStr(no_ofWeigh)) ;
              vhistory.push_back ((string) getSerialNo (cnt) + " Relative Aspect Series file path  =" + (string) rasfile) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Search algorithm for  star catalogue database =" + (string) convertIntToStr(search_algo_ctlg)) ;
             
             if(search_algo_ctlg==1 || search_algo_ctlg==3 || search_algo_ctlg==5){
            vhistory.push_back ((string) getSerialNo (cnt) + "Lengh of rectangle search   = " + (string)(len_a)) ;    
             vhistory.push_back ((string) getSerialNo (cnt) + " Width of rectangle search  = " + (string) (len_b)) ;  
             }
             vhistory.push_back ((string) getSerialNo (cnt) + "Registration and averaging image size    = " + (string) convertIntToStr(FINALFRAMESIZE_REGAVG)) ;
             
             //vhistory.push_back ((string) getSerialNo (cnt) + " frame to be discarded in reference frame calculation = " + (string)convertIntToStr(frames_toDiscard)) ;
            //vhistory.push_back ((string) getSerialNo (cnt) + " Average Factor = " + (string)convertIntToStr( nFrameToAverage)) ;
             vhistory.push_back ((string) getSerialNo (cnt) + " Distance from Which matching of stars to be started i.e Diff_dist = " + (string)convertFloatToStr ( diff_Dist)) ;
             //vhistory.push_back ((string) getSerialNo (cnt) + " Frequency domain filtering flag  " + (string) convertIntToStr(freqDomainFilter_Flag)) ;
             //if(freqDomainFilter_Flag==1)
             //vhistory.push_back ((string) getSerialNo (cnt) + " Type Filtering used   " + (string) convertIntToStr(type_Filtering)) ;
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
             vhistory.push_back ((string) getSerialNo (cnt) + " Write to disk flag for Weighted mean= " + (string)convertIntToStr( wtd_wm)) ;
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
