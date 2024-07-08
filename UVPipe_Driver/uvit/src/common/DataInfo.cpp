
/* 
 * File:   DataInfo.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include <vector>
#include"uvtUtils.h"
#include"DataInfo.h"

using namespace std ;

/*
 *Constructor of DataInfo - allocates memory to pointer data types
 */
DataInfo::DataInfo ()
{
    obsmode = new char[FLEN_VALUE] ;
    source = new char[FLEN_VALUE] ;
    filter = new char[FLEN_VALUE] ;
    detector = new char[FLEN_VALUE] ;
}

/*
 *Destructor - deallocates memory 
 */
DataInfo::~DataInfo ()
{
    delete[] obsmode ;
    delete[] source ;
    delete[] filter ;
    delete[] detector ;
}

/**
 * Function to set observation mode of data
 * @param obsmode
 */
void DataInfo::setObsMode (char *obsmode)
{
    strcpy (this->obsmode , obsmode) ;
    if (strcasecmp (obsmode , "PC") == 0){ modeflag = 0 ;
    }
    else  {
        modeflag = 1 ;
    }
    }


void DataInfo::setObsMode (bool modeflag)
{
    if (modeflag == 1) strcpy (obsmode , "IM") ;
    else strcpy (obsmode , "PC") ;
}

/**
 * Function to write OBS_MODE, FILTER, SOURCE, TSTARTI, TSTARTF, TSTOPI, TSTOPF, INT_TIME, 
 * DETECTOR,, XSIZE, YSIZE to information file
 * @param fptr - FITS file pointer pointing to the HDU where keywords are to be written
 */
void DataInfo::write (fitsfile *fptr)
{
    int status = 0 ;
    fits_update_key (fptr , TSTRING , "OBS_MODE" , obsmode , "Observation mode IM or PC" , &status) ;
    printError (status , "Error writing obs mode") ;
    fits_update_key (fptr , TSTRING , "FILTER" , filter , "Filter " , &status) ;
    printError (status , "Error writing filter") ;
    fits_update_key (fptr , TSTRING , "SOURCE" , source , "EC, ET or II" , &status) ;
    printError (status , "Error writing source") ;
    fits_update_key (fptr , TDOUBLE , "TSTARTI" , &tstarti , "Start Time " , &status) ;
    printError (status , "Error writing start time integer part") ;
    fits_update_key (fptr , TDOUBLE , "TSTARTF" , &tstartf , "Start Time " , &status) ;
    printError (status , "Error writing start time fraction part") ;
    fits_update_key (fptr , TDOUBLE , "TSTOPI" , &tstopi , "Stop Time " , &status) ;
    printError (status , "Error writing stop time integer part") ;
    fits_update_key (fptr , TDOUBLE , "TSTOPF" , &tstopf , "Stop Time " , &status) ;
    printError (status , "Error writing stop time fraction part") ;
    fits_update_key (fptr , TDOUBLE , "INT_TIME" , &integrationTime , "Integration Time " , &status) ;
    printError (status , "Error writing integration time") ;
    fits_update_key (fptr , TSTRING , "DETECTOR" , detector , "detector  FUV,NUV or VIS " , &status) ;
    printError (status , "Error writing detector") ;
    fits_update_key (fptr , TINT , "XSIZE" , &xsize , "xsize " , &status) ;
    printError (status , "Error writing xsize") ;
    fits_update_key (fptr , TINT , "YSIZE" , &ysize , "ysize " , &status) ;
    printError (status , "Error writing ysize") ;
}

/**
 * Function to read OBS_MODE, DETECTOR, FILTER, SOURCE, TSTARTI, TSTARTF, TSTOPI, TSTOPF, INT_TIME
 * @param fptr
 */
void DataInfo::getInfo (fitsfile *fptr)
{
    int status = 0 ;
    
    fits_read_key (fptr , TSTRING , "OBS_MODE" , obsmode , NULL , &status) ;
    printError (status , "Error reading obsmode") ;
    setObsMode (obsmode);
    fits_read_key (fptr , TSTRING , "DETECTOR" , detector , NULL , &status) ;
    printError (status , "Error reading detector") ;
    fits_read_key (fptr , TSTRING , "FILTER" , filter , NULL , &status) ;
    printError (status , "Error reading filter") ;
    fits_read_key (fptr , TSTRING , "SOURCE" , source , NULL , &status) ;
    printError (status , "Error reading source") ;
    fits_read_key (fptr , TDOUBLE , "TSTARTI" , &tstarti , NULL , &status) ;
    printError (status , "Error reading tstarti") ;
    fits_read_key (fptr , TDOUBLE , "TSTARTF" , &tstartf , NULL , &status) ;
    printError (status , "Error reading tstartf") ;
    fits_read_key (fptr , TDOUBLE , "TSTOPF" , &tstopf , NULL , &status) ;
    printError (status , "Error reading tstopf") ;
    fits_read_key (fptr , TDOUBLE , "TSTOPI" , &tstopi , NULL , &status) ;
    printError (status , "Error reading tstopi") ;
  
    fits_read_key (fptr , TDOUBLE , "INT_TIME" , &integrationTime , NULL , &status) ;          //integration time won't be their in first file
    if (status == KEY_NO_EXIST)   status=0;                                                                                  //Handled because of DataImgest module where reading 
                                                                                                                                                      //will be done from level-1 science data file                
    printError (status , "Error reading integration time") ;

    fits_read_key (fptr , TINT , "XSIZE" , &xsize , NULL , &status) ;                      //XSIZE will not be there in level-1 science data file
    if (status == KEY_NO_EXIST)
    {
        status = 0 ;
        fits_read_key (fptr , TINT , "WIN_X_SZ" , &x , NULL , &status) ;
        if (status == KEY_NO_EXIST)
        {
            
            status = 0 ;
        }
        printError (status , "Error reading WIN_X_SZ") ;
        fits_read_key (fptr , TINT , "WIN_XOFF" , &xoff , NULL , &status) ;
        if (status == KEY_NO_EXIST)
        {
           
            status = 0 ;
        }
        printError (status , "Error reading WIN_XOFF") ;
     //   LOG(INFO)<<x<<" "<<xoff;;exit(1);
        xsize= x - xoff + 1 ;
     }
    printError (status , "Error reading xsize") ;
    
    fits_read_key (fptr , TINT , "YSIZE" , &ysize , NULL , &status) ;                //YSIZE will not be there in level-1 science data file
    if (status == KEY_NO_EXIST)
    {
        status = 0 ;
        fits_read_key (fptr , TINT , "WIN_Y_SZ" , &y , NULL , &status) ;
        if (status == KEY_NO_EXIST)
        {
             LOG(INFO)<<"INSIDE";
            status = 0 ;
        }
        printError (status , "Error reading WIN_Y_SZ") ;
        fits_read_key (fptr , TINT , "WIN_YOFF" , &yoff , NULL , &status) ;
        if (status == KEY_NO_EXIST)
        {
             LOG(INFO)<<"INSIDE";
            status = 0 ;
        }
        printError (status , "Error reading WIN_YOFF") ;
        ysize = y - yoff + 1 ;
    }
    printError (status , "Error reading ysize") ;
    //LOG(INFO)<<xsize<<" "<<ysize;exit(1);
            
//    if (xsize <= 0 || ysize <= 0)
//    {
//        LOG(ERROR) << endl << "***Unable to set xsize and ysize***" << endl ;
//        exit (EXIT_FAILURE) ;
//    }

}
