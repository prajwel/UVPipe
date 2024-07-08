/* 
 * File:   uvtComputeThermal.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#include<iostream>
#include<unistd.h>
#include<stdlib.h>
#include<dirent.h>  //Accessing Directory
#include<string.h> 
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <uvtComputeThermal.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<caldb_Handler.h>
#include<glog/logging.h>
#include<macro_def.h>
#include<transform.h>
//fitsfile *fitin;
// int status =0;
//constructor: called when the object is created


uvtComputeThermal::uvtComputeThermal ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}
//Destructor


uvtComputeThermal::~ uvtComputeThermal () { }

//for parameter file reading


int uvtComputeThermal::read (int argc , char** argv)
{
    int status = 0 ;
    
     
    status = readParams (argc , argv , 7 , FNAME , "caldbpath" , caldb_path ,
            FNAME , "HKpath" , LBTpath ,
            FNAME , "inputdir" , indir ,
            FNAME , "outdir" , outdir ,
            BOOL , "clobber" , &clobber ,
            BOOL , "history" , &history ,
            STRING , "mode" , mode) ;

    if (status) return (EXIT_FAILURE) ;
  string tempfilepath = searchFile (indir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }

    sprintf (infofile_in , "%s/%s" , indir , tempfilepath.c_str()) ;
    char *obs_mode = new char[FLEN_VALUE] ;
    getKeywordVal (infofile_in , "OBS_MODE" , 2 , obs_mode) ;
    if (strcasecmp (obs_mode , "PC") == 0)
    {
        status = readParams (argc , argv , 1 , FNAME , "frameIntegrationDir" ,InframeDir) ;
        if (status) return (EXIT_FAILURE) ;

    }
    else if (strcasecmp (obs_mode , "IM") == 0){
        status = readParams (argc , argv , 1 , FNAME , "DataIngestDir" ,InframeDir) ;
        if (status) return (EXIT_FAILURE) ;
    }
    return (EXIT_SUCCESS) ;
}


int uvtComputeThermal::read (char *indir , char *dataIngestindir,char*caldbdir , char* hkdir , char *outdir , int clobber , int history)
{
    strcpy (this->indir , indir) ;
     strcpy (this->InframeDir , dataIngestindir) ;
    strcpy (this->caldb_path , caldbdir) ;
    strcpy (this->outdir , outdir) ;
    strcpy (this->LBTpath , hkdir) ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_FAILURE) ;
}

//Display the Parameter file content 


void uvtComputeThermal::display ()
{
    LOG (INFO) << endl << "-----------------------------------------------------------------------------------" ;
    LOG (INFO) << endl << "-----------------uvtComputeThermalParameters Display---------" ;
    LOG (INFO) << endl << "Input caldb  Directory                        : " << caldb_path ;
    //LOG(INFO)<<endl<<"Threshold value                                    : "<<Threshold_value;
    LOG (INFO) << endl << "Output Directory                               : " << outdir ;
    LOG (INFO) << endl << "HK  file path                               : " << LBTpath ;
    if (clobber == YES)
        LOG (INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG (INFO) << endl << "Overwrite                                         : NO" ;
    if (history == YES)
        LOG (INFO) << endl << "History                                             : YES" ;
    else
        LOG (INFO) << endl << "History                                              : NO" ;
    LOG (INFO) << endl << "-----------------------------------------------------------------------------------" ;
    LOG (INFO) << endl << "--------uvtComputeThermal  parameters Parameters Display Ends-----------" ;

}

//PixPadding Process started


int uvtComputeThermal::uvtThermalCalcProcess ()
{
    //creating module output directory name and checking whether it already exists

    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG (INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;
    int status = 0 ;
    fitsfile *finfo_in ;
    //check  existence  of output directory if output  exist and 
    //clobber =YES than remove that directory and recreate it,if exist and clobber =NO than exit from module.
    //if directory is not exist than create it .
    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;

    LOG (INFO) << endl << moduleoutdir << "  directory created" ;
    LOG (INFO) << indir << endl ;
    strcpy (caldb_orig_path , caldb_path) ;
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
   
    string tempfilepath = searchFile (indir , ".info") ;
    /**Setting the path for the input information File**/
    sprintf (infofile_in , "%s/%s" , indir , tempfilepath.c_str()) ;

    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input information file" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of the  input information file" , infofile_in) ;
    datainfo.getInfo (finfo_in) ;
    fits_close_file (finfo_in , &status) ;
    printError (status , "Error in closing  the input information file" , infofile_in) ;
    string tempname = caldb_handler.getThermalFile (datainfo.getDetector () , caldb_path) ;
    if (tempname == " ")
    {
        LOG (ERROR) << endl << "Couldn't find thermal  file from caldb" << endl ;
        return (EXIT_FAILURE) ;
    }

    joinStrings (thermalFile , 2 , caldb_path , tempname.c_str()) ;
    LOG (INFO) << endl << "Thermal  file From caldb  is  " << thermalFile ;
   
    status=readThermalFile () ; //reading thermal file from CALDB
    if(status){
        cout<<"Error in reading the thermal CALDB  file"<<endl;
        return(EXIT_SUCCESS);
    }
    if (thermalCalc ()) return (EXIT_FAILURE) ;

    LOG (INFO) << endl << "uvtComputeThermal  process completed successfully" << endl ;

    return (EXIT_SUCCESS) ;
}


int uvtComputeThermal::thermalCalc ()
{
    LOG (INFO) << endl << "Inside Thermal Calculation" ;
    int status = 0 ;
    bool found=FALSE;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    fitsfile *dataIngest_file,*fptr;
    string tempfilepath ;
     // if (datainfo.getModeFlag () == IM){
         tempfilepath  = searchFile (InframeDir , "_di.info") ;
           if (tempfilepath==" "){
               tempfilepath = searchFile(InframeDir, "_rfc.info");
        if (tempfilepath  ==" "){
             LOG(ERROR)<<"Input Directory does not contain required file"<<endl;
             return(EXIT_FAILURE);
         }
         else{
             found=TRUE;
         }
           }
    /**Setting the path for the input information File**/
    sprintf (dataIngestinfofile_in , "%s/%s" , InframeDir , tempfilepath.c_str()) ;
    fits_open_file (&dataIngest_file , dataIngestinfofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input information file" , dataIngestinfofile_in) ;
    fits_movabs_hdu (dataIngest_file , 2 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of the  input information file" , dataIngestinfofile_in) ;
  ///  if(datainfo.getModeFlag ()==IM){
    fits_read_key (dataIngest_file, TINT , "NFILES" , &nframes , NULL , &status) ;
    printError (status , "NFILES not Found" , dataIngestinfofile_in) ;
    if(found==TRUE)
    {
          fits_read_key (dataIngest_file, TSTRING , "CENTDIR" ,sigdir , NULL , &status) ;
    printError (status , "SIGDIR  not Found" , dataIngestinfofile_in) ; 
    }
    else
    {
    fits_read_key (dataIngest_file, TSTRING , "SIGDIR" ,sigdir , NULL , &status) ;
    printError (status , "SIGDIR  not Found" , dataIngestinfofile_in) ; 
    }
    
    char **signalframelist = allocateMemory<char>(nframes , NAMESIZE) ; //for  storing the  input frame list.
    fits_read_col (dataIngest_file , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) signalframelist , NULL , &status) ;
    printError (status , "Error in reading the column " , infofile_in) ;
    fits_close_file(dataIngest_file,&status);
    printError (status , "Error in closing file " , infofile_in) ;
    char infile[FLEN_FILENAME];
     double  frametime;
     vector<double> time_tracker;
     //loop for getting the frame time .
    for(int i=0;i<nframes;i++)
    {
        
    sprintf (infile , "%s/%s/%s" , InframeDir , sigdir, signalframelist[i]) ;
    fits_open_file (&fptr , infile , READONLY , &status) ;
    printError (status , "***Error in opening file***" , infile) ;
    
    if(found==TRUE){
   fits_movabs_hdu (fptr , 2 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of the  input information file" , infile) ;        
    }
    fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
    printError (status , "***FRMTIME keyword not Found***") ;
  
    time_tracker.push_back (frametime);
  
    fits_close_file(fptr,&status);
    printError (status , "Error in closing the file**") ;
    
    }
     
  
    vector<double> yaw_corr_vect , pitch_corr_vect , roll_corr_vect,time_vect ;
    //method for reading  4 temperatures of LBT file
    status = getTempFromLBT () ;
    if (status)
    {
        LOG (ERROR) << "***temperature reading From the HK file fails***" << endl ;
        return (EXIT_FAILURE) ;
    }

    bool flag_matchFound = FALSE ; //flag for checking corresponding  temperature found or not
 
    
    //loop for matching temperature of HK file to temperature of  caldb file 
  
    fitsfile *flbt;
    fits_open_file (&flbt ,LBTpath , READONLY , &status) ;
    printError (status , "Error in opening lbt file" ,LBTpath ) ;
    fits_movabs_hdu (flbt , 2 , NULL , &status) ;
    printError (status , "Error moving to hdu 2 of LBT FILE" ,LBTpath) ;
    fits_get_num_rows (flbt , &num_record_LBT , &status) ; 
    printError (status , "Error in getting number of rows" ,LBTpath) ;
    fits_close_file(flbt,&status);
    printError (status , "Error in closing the file" ,LBTpath) ;
  
    for (int index_lbt = 0 ; index_lbt < num_record_LBT ; index_lbt ++)
    {
        flag_matchFound = FALSE ;
        for (int i = 0 ; i < nCalDBTempValues ; i ++)
        {
          
            if ((MB_top_temp[index_lbt] > temp_One_caldb[i] && MB_top_temp[index_lbt] < temp_One_caldb[i + 1]) || (MB_top_temp[index_lbt] == temp_One_caldb[i]))
            {

                for (int index_secnd_tempHK = i ; index_secnd_tempHK < i + cnt_Fortemp2_repeat ; index_secnd_tempHK ++)
                {

                    if ((TT_top_temp[index_lbt] > temp_Two_caldb[index_secnd_tempHK] && TT_top_temp[index_lbt] < temp_Two_caldb[index_secnd_tempHK + 1]) || (TT_top_temp[index_lbt] == temp_Two_caldb[index_secnd_tempHK]))
                    {

                        for (int index_three_tempHK = index_secnd_tempHK ; index_three_tempHK < index_secnd_tempHK + cnt_Fortemp3_repeat ; index_three_tempHK ++)
                        {


                            if ((TT_bot_temp[index_lbt] > temp_Three_caldb[index_three_tempHK] && TT_bot_temp[index_lbt] < temp_Three_caldb[index_three_tempHK + 1]) || (TT_bot_temp[index_lbt] == temp_Three_caldb[index_three_tempHK]))
                            {

                                for (int index_Four_tempHK = index_three_tempHK ; index_Four_tempHK < index_three_tempHK + cnt_Fortemp4_repeat ; index_Four_tempHK ++)
                                {

                                    if ((BOT_ring_temp[index_lbt] > temp_Four_caldb[index_Four_tempHK] && BOT_ring_temp[index_lbt] < temp_Four_caldb[index_Four_tempHK + 1]) || (BOT_ring_temp[index_lbt] == temp_Four_caldb[index_Four_tempHK]))
                                    {
                                       
                                        time_vect.push_back (time_tracker[index_lbt]);
                                        yaw_corr_vect.push_back (yaw_corr_caldb[index_Four_tempHK]/ARC_PERSEC_TO_DEGREE) ;
                                        pitch_corr_vect.push_back (pitch_corr_caldb[index_Four_tempHK]/ARC_PERSEC_TO_DEGREE) ;
                                        roll_corr_vect.push_back (0.0f) ;
                                        flag_matchFound = TRUE ;
                                        break ;
                                    }

                                }
                                if (flag_matchFound == TRUE)
                                {
                                    break ;
                                }
                            }

                        }
                        if (flag_matchFound == TRUE)
                        {
                            break ;
                        }
                    }
                }
                if (flag_matchFound == TRUE)
                {
                    break ;
                }
            }

        }

    }

    /*Introduce matrix for transformation between star250 to spacecraft.
      * 
      */
    
    writeThermalcorrTotable (time_vect,yaw_corr_vect , pitch_corr_vect , roll_corr_vect) ;

    yaw_corr_vect.clear () ;
    pitch_corr_vect.clear () ;
    return (EXIT_SUCCESS) ;


}


int uvtComputeThermal::getHistory (vector<string> &vhistory)
{
    int cnt = 0 ;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ((string) getSerialNo (cnt) + "Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Caldb Di=" + (string) caldb_path) ;
    //vhistory.push_back("P2 Flat Field file " + (string)flatfieldfile);
    vhistory.push_back ((string) getSerialNo (cnt) + " outdir=" + (string) outdir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Module Output directory=" + (string) moduleoutdir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " HK  file path=" + (string) LBTpath) ;
    if (clobber == YES)
        vhistory.push_back ((string) getSerialNo (cnt) + " clobber=yes") ;
    else
        vhistory.push_back ((string) getSerialNo (cnt) + " clobber=no") ;
    if (history == YES)
        vhistory.push_back ((string) getSerialNo (cnt) + " history=yes") ;
    else
        vhistory.push_back ((string) getSerialNo (cnt) + " history=no") ;
    vhistory.push_back ("Parameter List END") ;
    return (EXIT_SUCCESS) ;
}


int uvtComputeThermal::readThermalFile ()
{

    LOG (INFO) << endl << "Reading Thermal  file from calDB" ;
    int status = 0 , hdutype ;
    long nrows ;

    //opening caldb thermal 
    fitsfile *fthermal ;
    fits_open_file (&fthermal , thermalFile , READONLY , &status) ;
    printError (status , "Error in opening the thermal file" , thermalFile) ;
    fits_movabs_hdu (fthermal , 2 , &hdutype , &status) ;
    printError (status , "Error in readThermalFile" , thermalFile) ;

    if (hdutype != BINARY_TBL)
    {
        LOG (ERROR) << endl << "***Expected binary table at HDU 2 of Thermal File*** " << endl ;
        return (EXIT_FAILURE) ;
    }

    fits_get_num_rows (fthermal , &nrows , &status) ;
    printError (status , "Error in readThermalFile()" , thermalFile) ;

    //allocating memory to arrays 
    nCalDBTempValues = nrows ;
    LOG (INFO) << "Total  number of rows " << nrows << endl ;
    temp_One_caldb = new float[nrows] ;
    temp_Two_caldb = new float[nrows] ;
    temp_Three_caldb = new float[nrows] ;
    temp_Four_caldb = new float[nrows] ;
    yaw_corr_caldb = new float[nrows] ;
    pitch_corr_caldb = new float[nrows] ;
    fits_close_file (fthermal , &status) ;
    printError (status , "***Error in readThermalFile()***") ;

    //reading columns of thermal file
    status = readColumnsFromFITS (thermalFile , 2 , 6 , TFLOAT , 1 , temp_One_caldb , nrows , TFLOAT , 2 , temp_Two_caldb , nrows ,
            TFLOAT , 3 , temp_Three_caldb , nrows , TFLOAT , 4 , temp_Four_caldb , nrows ,
            TFLOAT , 5 , yaw_corr_caldb , nrows , TFLOAT , 6 , pitch_corr_caldb , nrows) ;
    if (status)
    {
        LOG (INFO) << "Error reading column from thermal file" << endl ;
        return (EXIT_FAILURE) ;
    }

    LOG (INFO) << endl << "Reading Thermal   file from caldb Finished........" ;



    cnt_Fortemp1_repeat = 0 ;
    cnt_Fortemp2_repeat = 0 ;
    cnt_Fortemp3_repeat = 0 ;
    cnt_Fortemp4_repeat = 0 ;
    int cnt_caldb = 0 ;

    //loop for calculating number of  counts  each temperature repeats.
    for (int i = 0 ; i < nrows - 1 ; i ++)
    {
        if (temp_One_caldb[i] == temp_One_caldb[i + 1])
        {
            cnt_caldb ++ ;
        }
        else
        {
            cnt_Fortemp1_repeat = cnt_caldb + 1 ;
            break ;
        }
    }
    cnt_caldb = 0 ;
    for (int i = 0 ; i < nrows ; i ++)
    {
        if (temp_Two_caldb[i] == temp_Two_caldb[i + 1])
        {
            cnt_caldb ++ ;
        }
        else
        {
            cnt_Fortemp2_repeat = cnt_caldb + 1 ;
            break ;
        }
    }
    cnt_caldb = 0 ;
    for (int i = 0 ; i < nrows ; i ++)
    {
        if (temp_Three_caldb[i] == temp_Three_caldb[i + 1])
        {
            cnt_caldb ++ ;
        }
        else
        {
            cnt_Fortemp3_repeat = cnt_caldb + 1 ;
            break ;
        }
    }
    cnt_caldb = 0 ;
    for (int i = 0 ; i < nrows ; i ++)
    {
        if (temp_Four_caldb[i] == temp_Four_caldb[i + 1])
        {
            cnt_caldb ++ ;
        }
        else
        {
            cnt_Fortemp4_repeat = cnt_caldb + 1 ;
            break ;
        }
    }

    return (EXIT_SUCCESS) ;
}


int uvtComputeThermal::getTempFromLBT ()
{
    int colnum=0;
     int status = 0 ;
    long total_rows = 0 ;
    char *temp = basename (LBTpath) ;
    string strtemp (temp) ;
    int pos = strtemp.find ("level1.lbt") ;
    if (pos > 0 && pos < strtemp.length ())
    {
        strtemp.replace (pos , string ("level1.lbt").length () , "l2") ;
    }
    strcpy (nameprefix , strtemp.c_str ()) ; //storing filename prefix  
    fitsfile *flbt ;
    fits_open_file (&flbt ,LBTpath , READONLY , &status) ;
    printError (status , "Error in opening LBT file" ,LBTpath ) ;
 
  
 
    fits_open_file (&flbt ,LBTpath , READONLY , &status) ;
    printError (status , "Error in opening LBT file" ,LBTpath ) ;
    fits_movabs_hdu (flbt , 2 , NULL , &status) ;
    printError (status , "Error moving to HDU 2 of HK FILE" ,LBTpath) ;
    fits_get_num_rows (flbt , &total_rows , &status) ; //getting number of rows from HK file
    printError (status , "Error reading the number of rows in HK file" , LBTpath) ;
    
    MB_top_temp= new unsigned short[total_rows];
    TT_top_temp= new unsigned short[total_rows];
    TT_bot_temp= new unsigned short[total_rows];
    BOT_ring_temp= new unsigned short[total_rows];
    
    if(strcmp(datainfo.getDetector (),"NUV")==0)
    {
        
        fits_get_colnum(flbt,CASEINSEN,"NUV_MB_TOP_LOC_AT_MID_Th_91_TS167",&colnum,&status);
        printError (status , "Error reading the column  number of the LBT file" , LBTpath) ;
        fits_read_col (flbt , TUSHORT , colnum , 1 , 1 , total_rows, NULL , MB_top_temp , NULL , &status) ;
        printError (status , "Reading a column Fails in caldb" , LBTpath) ;
        fits_get_colnum(flbt,CASEINSEN,"NUV_TT_TOP_HTR_TOP_M_SN_302_FTS9",&colnum,&status);
        printError (status , "Error reading the column  number of the LBT file" , LBTpath) ;
        fits_read_col (flbt , TUSHORT , colnum , 1 , 1 , total_rows, NULL , TT_top_temp , NULL , &status) ;
        printError (status , "Reading a column Fails in caldb" , LBTpath) ;
        fits_get_colnum(flbt,CASEINSEN,"NUV_TT_BOTTOM_HTR_TOP_M_SN_306_FTS13",&colnum,&status);
        printError (status , "Error reading the column  number of the LBT file" , LBTpath) ;
        fits_read_col (flbt , TUSHORT , colnum , 1 , 1 , total_rows, NULL , TT_bot_temp , NULL , &status) ;
        printError (status , "Reading a column Fails in caldb" , LBTpath) ;
       
        fits_get_colnum(flbt,CASEINSEN,"NUV_BOT_RING_OUT_HTR_M_SN_308_FTS15",&colnum,&status);
        printError (status , "Error reading the column  number of the LBT file" , LBTpath) ;
        fits_read_col (flbt , TUSHORT , colnum , 1 , 1 , total_rows, NULL ,BOT_ring_temp , NULL , &status) ;
        printError (status , "Reading a column Fails in caldb" , LBTpath) ;

    }
    else if(strcmp(datainfo.getDetector (),"FUV")==0)
    {
      
        fits_get_colnum(flbt,CASEINSEN,"FUV_MB_TOP_LOC_AT_MID_Th_34_TS151",&colnum,&status);
        printError (status , "Error reading the column  number of the output file" , LBTpath) ;
        fits_read_col (flbt , TUSHORT , colnum , 1 , 1 , total_rows, NULL , MB_top_temp , NULL , &status) ;
        printError (status , "Reading a column Fails in caldb" , LBTpath) ;
        fits_get_colnum(flbt,CASEINSEN,"FUV_TT_TOP_HTR_TOP_M_SN_294_FTS1",&colnum,&status);
        printError (status , "Error reading the column  number of the output file" , LBTpath) ;
        fits_read_col (flbt , TUSHORT , colnum , 1 , 1 , total_rows, NULL , TT_top_temp , NULL , &status) ;
        printError (status , "Reading a column Fails in caldb" , LBTpath) ;
        fits_get_colnum(flbt,CASEINSEN,"FUV_TT_BOTTOM_HTR_TOP_M_SN_298_FTS5",&colnum,&status);
        printError (status , "Error reading the column  number of the output file" , LBTpath) ;
        fits_read_col (flbt , TUSHORT , colnum , 1 , 1 , total_rows, NULL , TT_bot_temp , NULL , &status) ;
        printError (status , "Reading a column Fails in caldb" , LBTpath) ;
        fits_get_colnum(flbt,CASEINSEN,"FUV_BOT_RING_OUT_HTR_M_SN_300_FTS7",&colnum,&status);
        printError (status , "Error reading the column  number of the output file" , LBTpath) ;
        fits_read_col (flbt , TUSHORT , colnum , 1 , 1 , total_rows, NULL ,BOT_ring_temp , NULL , &status) ;
        printError (status , "Reading a column Fails in caldb" , LBTpath) ;
        
    }
      fits_close_file (flbt , &status) ;
     printError (status , "Error in closing the lbt file" , LBTpath) ;
     
    
    
//    int status = 0 ;
//    long total_rows = 0 ;
//    char *temp = basename (HKpath) ;
//
//    string strtemp (temp) ;
//    int pos = strtemp.find ("level1.hk") ;
//    if (pos > 0 && pos < strtemp.length ())
//    {
//        strtemp.replace (pos , string ("level1.hk").length () , "l2") ;
//    }
//    strcpy (nameprefix , strtemp.c_str ()) ; //storing filename prefix  
//    //opening input HK file 
//    fitsfile *fhk ;
//    fits_open_file (&fhk , HKpath , READONLY , &status) ;
//    printError (status , "Error in opening hk file" , HKpath) ;
//
//    copyUsrkeywrdsTovect (fhk , key_records) ;
//
//    fits_movabs_hdu (fhk , 2 , NULL , &status) ;
//    printError (status , "Error moving to HDU 2 of HK FILE" , HKpath) ;
//    fits_get_num_rows (fhk , &total_rows , &status) ; //getting number of rows from HK file
//    printError (status , "Error reading the number of rows in HK file" , HKpath) ;
//
//    num_record_HK = total_rows ; //setting  total number of records of HK file
//
//    //array for storing monitor 2 information from HK file
//    unsigned char *monitor_Two_data = new unsigned char[num_record_HK * MONITOR2_BYTES] ;
//
//    //reading monitor2 information from HK file to array
//    fits_read_col (fhk , TBYTE , 11 , 1 , 1 , num_record_HK * MONITOR2_BYTES , NULL , monitor_Two_data , NULL , &status) ;
//    printError (status , "Reading a column Fails in caldb" , HKpath) ;
//
//    //memory allocation for the  storing 4  HK temperature 
//    temp_One_HK = new float[num_record_HK] ;
//    temp_Two_HK = new float[num_record_HK] ;
//    temp_Three_HK = new float[num_record_HK] ;
//    temp_Four_HK = new float[num_record_HK] ;
//
//
//    unsigned short data1 , data2 ;
//    float final_word ;
//    ofstream fileout("outHktemp.txt");
//    //calculating temperature values from HK file (MONITOR-2)
//    for (int i = 0 ; i < num_record_HK ; i ++)
//    {
//
//        data1 = monitor_Two_data[MONITOR2_BYTES * i] ;
//        data2 = monitor_Two_data[MONITOR2_BYTES * i + 1] ;
//        final_word = (data1 << 8) | data2 ;
//        temp_One_HK[i] = final_word ;
//        data1 = monitor_Two_data[MONITOR2_BYTES * i + 2] ;
//        data2 = monitor_Two_data[MONITOR2_BYTES * i + 3] ;
//        final_word = (data1 << 8) | data2 ;
//        temp_Two_HK[i] = final_word ;
//        data1 = monitor_Two_data[MONITOR2_BYTES * i + 4] ;
//        data2 = monitor_Two_data[MONITOR2_BYTES * i + 5] ;
//        final_word = (data1 << 8) | data2 ;
//        temp_Three_HK[i] = final_word ;
//        data1 = monitor_Two_data[MONITOR2_BYTES * i + 6] ;
//        data2 = monitor_Two_data[MONITOR2_BYTES * i + 7] ;
//        final_word = (data1 << 8) | data2 ;
//        temp_Four_HK[i] = final_word ;
//        fileout<<temp_One_HK[i]<<setw(20)<<temp_Two_HK[i]<<setw (20)<<temp_Three_HK[i]<<setw(20)<<temp_Four_HK[i]<<endl;
//
//    }
//    fileout.close ();
//
//    fits_close_file (fhk , &status) ;
//    printError (status , "Error in closing the HKfile" , HKpath) ;

    return (EXIT_SUCCESS) ;
}


int uvtComputeThermal::writeThermalcorrTotable (vector<double> &timeData, vector<double> &yawcorr , vector<double> &pitchcorr , vector<double> &rollcorr)
{
    int status = 0 ;
    char outfile[FLEN_FILENAME] ;
    fitsfile *flbt;
    fits_open_file (&flbt ,LBTpath , READONLY , &status) ;
    printError (status , "Error in opening LBT file" ,LBTpath ) ;
    copyUsrkeywrdsTovect (flbt , key_records) ;
    fits_close_file(flbt,&status);
    vector<string> vhistorystr ;
    if (history == YES) getHistory (vhistorystr) ;

    fitsfile *fout ;
    sprintf (outfile , "%s/%s_th.fits" , moduleoutdir , nameprefix) ;
    fits_create_file (&fout , outfile , &status) ;
    printError (status , "Error in creating the output thermal file" , outfile) ;

    int tfields = 4 ;
    char *ttype[] = {"Time","YAW" , "PITCH" , "ROLL"} ;
    char *tform2[] = {"D","E" , "E" , "E"} ;

    fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform2 , NULL , "Corrected Pitch and Yaw" , &status) ;
    printError (status , "Error in creating the table" , outfile) ;

    status=writeColumnsToFITS (outfile , 2 , 4 , TDOUBLE,1,timeData.data (),(int)timeData.size (),TDOUBLE , 2 , yawcorr.data () , (int) yawcorr.size () , TDOUBLE , 3 , pitchcorr.data () , (int) yawcorr.size () , TDOUBLE , 4 , rollcorr.data () , (int) rollcorr.size ()) ;
    if(status){
        LOG(INFO)<<"Error in writing to the  thermal  FITS  file"<<endl;
       return(EXIT_FAILURE);
    }
  
    //cout<<"leba "<<key_records.size ()<<endl;exit(1);
    writeUsrkeywordsFrmvect (outfile , key_records) ;
   // cout<<"leba "<<key_records.size ()<<endl;exit(1);
    if (history == YES) writeHistory (outfile , vhistorystr) ;
  // cout<<"leba "<<key_records.size ()<<endl;exit(1);
    writeCommonKeywords (fout , modulename) ;
      fits_close_file (fout , &status) ;
    printError (status , "Error in closing the output thermal file" , outfile) ;

    
    return (EXIT_SUCCESS) ;

}
