/* 
 * File:   main.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on August 2, 2013, 09:06 PM
 */
 

#include<iostream>
#include<unistd.h>
#include<stdlib.h>
#include<dirent.h>  //Accessing Directory
#include<string.h> 
#include <cstdlib>
#include "stdio.h"
#include <fstream>
#include "uvtCentroidBias.h"
//#include "uvtDetectStar.h"
#include<pthread.h>
#include<glog/logging.h>

uvtCentroidBias::uvtCentroidBias ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}

uvtCentroidBias::~uvtCentroidBias ()
{
    //delete[] x_corr , y_corr , fraction_bias  ;
}

int uvtCentroidBias::read (int argc , char** argv)
{
    int status = 0 ;

    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG(ERROR) << "***Error Initializing PIL***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("framelistDir" , inputdatadir)))
    {
        LOG(ERROR) << endl << "***Error reading input data directory name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("caldbDir" , caldbDir)))
    {
        LOG(ERROR) << endl << "***Error reading caldb  directory ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("outdir" , outdir)))
    {
        LOG(ERROR) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("clobber" , &clobber)))
    {
        LOG(ERROR) << "***Error Reading clobber:" << clobber << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("history" , &history)))
    {
        LOG(ERROR) << "***Error Reading history parameter:" << history << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("mode" , mode)))
    {
        LOG(ERROR) << "***Error Reading mode parameter:" << mode << "***" ;
        return status ;
    }
    /*---------Code to find data mode from information file inside input data dir-----------*/

    PILClose (status) ;
    return (EXIT_SUCCESS) ;
}

int uvtCentroidBias::read (char *inputdatadir , char *caldbDir , char *outdir , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->caldbDir , caldbDir) ;
    strcpy (this->outdir , outdir) ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}

void uvtCentroidBias::display ()
{
    LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT CENTROID BIAS PARAMETERS              " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input data directory  : " << inputdatadir ;
    LOG(INFO) << endl << "CALDB directory : " << caldbDir ;
    LOG(INFO) << endl << "Output Directory : " << outdir ;
    if (clobber == YES)
        LOG(INFO) << endl << "Overwrite : YES" ;
    else
        LOG(INFO) << endl << "Overwrite : NO" ;
    if (history == YES)
        LOG(INFO) << endl << "History    : YES" ;
    else
        LOG(INFO) << endl << "History     : NO" ;
    LOG(INFO) << endl << "------------------------------------------------------------------------" << endl ;
}

int uvtCentroidBias::uvtCentroidBiasProcess ()
{
    
    LOG(INFO) << endl << "Centroid Bias process started" << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;

    LOG(INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;
    string cmd ;

    if (DirExists (moduleoutdir) && clobber == YES)
    {
        LOG(INFO) << endl << "Directory exists and clobber=yes" ;
        cmd = (string) "rm -rf " + (string) moduleoutdir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists (moduleoutdir) && clobber == NO)
    {
        LOG(INFO) << endl << moduleoutdir << "  already exists " ;
        LOG(INFO) << endl << "Use clobber=yes for overwriting" ;
        return (EXIT_FAILURE) ;
    }
    nframes = 0 ;
    /**Shell command for creating the output Directory**/
    cmd = "mkdir -p " + (string) moduleoutdir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ; // creating output directory to keep output from unitConversion
    LOG(INFO) << endl << moduleoutdir << "  directory created" ;
    //opening info file in input directory to get data information
    string tempfilepath = searchFile(inputdatadir , "info") ;
    if (tempfilepath == " ")
    {
        LOG(ERROR) << endl << "Error in finding info file" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    LOG(INFO)<<"\nInput information file :"<<infofile_in<<endl;
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input Information File" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to the perticuler HDU " , infofile_in) ;
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ;
    fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
    printError (status , "Error reading the key value of the EVTFILE" , infofile_in) ;
    datainfo.getInfo(finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    if (xsize <= 0 || ysize <= 0)
    {
        LOG(ERROR) << endl << "***Invalid xsize/ysize***\n" ;
        return (EXIT_FAILURE) ;
    }
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error reading the NAMEPREX key in information file" , infofile_in) ; //for creating name for output information file
    fits_close_file(finfo_in,&status);
    
    //creating output information file
    sprintf (infofile_out , "%s/%s_cb.info" , moduleoutdir , nameprefix) ;
    LOG(INFO)<<"\nOutput information file :"<<infofile_out<<endl;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error creating the output information File" , infofile_out) ;
    char *ttype[] = {"FileList"} ;
    char *tform[] = {"A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error creating the table in  out information File" , infofile_out) ;

    datainfo.write (finfo_out) ; //writing basic data information

    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error updating the NAMEPREX key  in out information file" , infofile_out) ; //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    char file_in[NAMESIZE] , file_out[NAMESIZE] ; //for input and output event file
    sprintf (file_in , "%s/%s" , inputdatadir , eventfile) ; //taking event file full path
    sprintf (file_out , "%s/%s_cb.events" , moduleoutdir , nameprefix) ; //setting output event file path
    LOG(INFO)<<"\nInput Event file :"<<file_in<<endl;
    LOG(INFO)<<"\nOutput Event file:"<<file_out<<endl;
    fitsfile *fevt_in , *fevt_out  ; //pending work for IMG File
    fits_open_file (&fevt_in , file_in , READWRITE , &status) ;
    printError (status , "Error in opening input event File" , file_in) ;
    fits_read_key (fevt_in , TSTRING , "CENTROID" , centroid_algo , NULL , &status) ;
    copyUsrkeywrdsTovect (fevt_in,key_record);
    LOG(INFO) << "\nUsing Centroid Algorithm " << centroid_algo << endl ;

    if ((strcmp (centroid_algo , "3C") == 0) || (strcmp (centroid_algo , "3S") == 0))
    {
        windowsize = 3 ;
    }
    else if (strcmp (centroid_algo , "5S") == 0)
    {
        windowsize = 5 ;
    }
    else
    {
        LOG(ERROR) << "***INVALID value for Centroid window.allowed values are (3/5) *** " << endl ;
        return (EXIT_FAILURE) ;
    }
    //printError (status , "***Error Reading the Key value for NAMEPRFX***" , infofile_in) ;
    string tempname = caldb_handler.getCentroidBiasFile (datainfo.getDetector () , caldbDir , windowsize) ;
    if (tempname ==" ")
    {
        LOG(ERROR) << endl << "Couldn't find CentroidBias file from caldb" << endl ;
        return (EXIT_FAILURE) ;
    }

    /***For Joining the String ***/
    joinStrings (centroidbiasfile , 2 , caldbDir , tempname.c_str()) ;
    LOG(INFO) << endl << "\nCentroid bias  file  from calDB :" << centroidbiasfile ;

    /*Reading a bias file From the CALDB dir*/
    status = readcentroidbiasFile () ;
    if (status)
    {
        LOG(ERROR) << "***Error Reading centroid bias File From the Caldb***" << endl ;
        return (EXIT_FAILURE) ;
    }
    /*sorting the X_corr  based on the cent_x */
    double   temp_y_corr ;
    double temp_x , temp_x_corr ;
    LOG(INFO)<<"\nSorting the  data of calDB..."<<endl;
    for (int i = 0 ; i < biasRows - 1 ; i++)
    {
        if (fraction_bias[i] > fraction_bias[i + 1])
        {
            temp_x = fraction_bias[i + 1] ;
            fraction_bias[i + 1] = fraction_bias[i] ;
            fraction_bias[i] = temp_x ;
            temp_x_corr = x_corr[i + 1] ;
            x_corr[i + 1] = x_corr[i] ;
            x_corr[i] = temp_x_corr ;
            temp_y_corr=y_corr[i+1];
             y_corr[i + 1] = y_corr[i] ;
            y_corr[i] = temp_y_corr ;
        }
    }
     /*sorting the y_corr  based on the cent_y */
//    for (int i = 0 ; i < biasRows - 1 ; i++)
//    {
//        if (cent_y[i] > cent_y[i + 1])
//        {
//            temp_y = cent_y[i + 1] ;
//            cent_y[i + 1] = cent_y[i] ;
//            cent_y[i] = temp_y ;
//            temp_y_corr = y_corr[i + 1] ;
//            y_corr[i + 1] = y_corr[i] ;
//            y_corr[i] = temp_y_corr ;
//        }
//    }
    fits_create_file (&fevt_out , file_out , &status) ;
    printError (status , "Error in creating a output event File" , file_out) ;
    fits_copy_file (fevt_in , fevt_out , 1 , 1 , 1 , &status) ;
    printError (status , "***Error in coping the header***" , file_out) ;
    double *t ;
    float *xFrac , *yFrac , *xFrac_temp , *yFrac_temp ;
    long nrows ;
    float *new_xFrac , *new_yFrac,*ff ;
    fits_movabs_hdu (fevt_out , 2 , NULL , &status) ;
    printError (status , "Error in  moving to perticuler HDU" , file_in) ;
    fits_get_num_rows (fevt_out , &nrows , &status) ;
    printError (status , "Error in  getting the nu,mber of rows in the  in input Event File" , file_in) ;
    /***Assigning a memory to Arrays***/
    xFrac = new float[nrows] ;
    yFrac = new float[nrows] ;
    xFrac_temp = new float[nrows] ;
    yFrac_temp = new float[nrows] ;
    new_xFrac = new float[nrows] ;
    new_yFrac = new float[nrows] ;
    ff=  new float[nrows];
    /**Reading XFraction and YFraction values from Event file**/
    fits_read_col (fevt_out , TFLOAT , 4 , 1 , 1 , nrows , NULL , xFrac , NULL , &status) ;
    printError (status , "Error reading xinteger" , file_in) ;
    fits_read_col (fevt_out , TFLOAT , 5 , 1 , 1 , nrows , NULL , yFrac , NULL , &status) ;
    printError (status , "Error reading Yinteger" , file_in) ;
    int x_sign_flag = 0 ;
    int y_sign_flag = 0 ;
    for (int i = 0 ; i < nrows ; i++)
    {
        xFrac_temp[i] = ((xFrac[i]-(int) xFrac[i])+32)/2 ;
        yFrac_temp[i] = ((yFrac[i]-(int) yFrac[i])+32)/2 ;
    }
    /*loop for updating a fractional Part of x and y   */
    LOG(INFO)<<"\nPerforming Centroid Bias..."<<endl;
    for (int i = 0 ; i < nrows ; i++)
    {
        if (xFrac[i] >= 0)
            x_sign_flag = 1 ;
        else
            x_sign_flag = 0 ;

        if (yFrac[i] >= 0)
            y_sign_flag = 1 ;
        else
            y_sign_flag = 0 ;

 //       xFrac[i] = fabs (xFrac[i]) ; //*32.0;
  //      yFrac[i] = fabs (yFrac[i]) ; //*32.0;
        xFrac[i] = fabs (xFrac[i]) ;
        yFrac[i] = fabs (yFrac[i]) ; 
        /*Searching for X Frac value in the Centroid Bias correction file in CalDB**/
        if (xFrac_temp[i] == fraction_bias[biasRows - 1])
        {
            new_xFrac[i] = x_corr[biasRows - 1] ;
            if (x_sign_flag == 0)
                new_xFrac[i] = -new_xFrac[i] ;
            
            if(new_xFrac[i]<xsize)
            new_xFrac[i] = new_xFrac[i] + (int) xFrac[i] ;
        }
        else
        {
            for (int j = 0 ; j < biasRows ; j++)
            {
               if (xFrac_temp[i] >= fraction_bias[j] && xFrac_temp[i] < fraction_bias[j + 1])
                { 
                    new_xFrac[i] = x_corr[j] + ((xFrac_temp[i] - fraction_bias[j]) * (x_corr[j + 1] - x_corr[j])) / (fraction_bias[j + 1] - fraction_bias[j]) ;
                    if (x_sign_flag == 0)
                        new_xFrac[i] = -new_xFrac[i] ;
                    if(new_xFrac[i]<xsize)
                    new_xFrac[i] = new_xFrac[i] + (int) xFrac[i] ;
                                   
                    break ;
                }
            }

        }
         /***Searching for Y Frac value in the Centroid Bias correction file in CalDB***/
        if (yFrac_temp[i] == fraction_bias[biasRows - 1])
        {

            new_yFrac[i] = y_corr[biasRows - 1] ;
            if (y_sign_flag == 0)
                new_yFrac[i] = -new_yFrac[i] ;
            
             if(new_yFrac[i]<ysize)
            new_yFrac[i] = new_yFrac[i] + (int) yFrac[i] ;
        }
        else
        {
            for (int j = 0 ; j < biasRows ; j++)
            {

                if (yFrac_temp[i] >= fraction_bias[j] && yFrac_temp[i] < fraction_bias[j + 1])
                {
                    new_yFrac[i] = y_corr[j] + ((yFrac_temp[i] - fraction_bias[j]) * (y_corr[j + 1] - y_corr[j])) / (fraction_bias[j + 1] - fraction_bias[j]) ;
                    if (y_sign_flag == 0)
                        new_yFrac[i] = -new_yFrac[i] ;
                     if(new_yFrac[i]<ysize)
                    new_yFrac[i] = new_yFrac[i] + (int) yFrac[i] ;
                    break ;
                }
            }
        }
}
    for(int i=0;i<nrows;i++)
        ff[i]=1.0;
    
    LOG(INFO)<<"\nCorrection of fractional part of event (X,Y) is done.writing to output event file.."<<endl;
    fits_write_col (fevt_out , TFLOAT , 4 , 1 , 1 , nrows , new_xFrac , &status) ;
    printError (status , "Error writing new_Xinteger" , file_out) ;
    fits_write_col (fevt_out , TFLOAT , 5 , 1 , 1 , nrows , new_yFrac , &status) ;
    printError (status , "Error writing new_Yinteger " , file_out) ;
    char *ttype1= {FF_COLNAME} ;
    char *tform1 = {"1D"} ;

    /**Inserting the new column of the Effective number of Photons**/
    fits_insert_col (fevt_out , 10, ttype1 , tform1 , &status) ;
    printError (status , "Error inserting to  Column" , file_out) ;
    fits_write_col (fevt_out , TFLOAT, 10 , 1 , 1 , nrows , ff , &status) ;
    printError (status , "Error writing to Column" , file_out) ;

    /*releasing the memory*/
    delete[] new_xFrac ;
    delete[] new_yFrac ;
    delete[] xFrac_temp ;
    delete[] yFrac_temp ;
    delete[] xFrac ;
    delete[] yFrac ;

    /*Closing the input Event File*/
    fits_movabs_hdu (fevt_out , 1 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU" , infofile_out) ;
    /***Writes the  common information to  output information File***/
   
    vector<string> vhistorystr ;
    getHistory (vhistorystr) ;
    if (history == YES) writeHistory (file_out , vhistorystr) ;
      fits_movabs_hdu (fevt_out , 1 , NULL , &status) ;
    printError (status , "Error in  moving to perticular  HDU" , file_in) ;
    writeCommonKeywords (fevt_out , modulename) ;
    fits_close_file (fevt_out , &status) ;
    printError (status , "Error in Closing the output Event File" , file_in) ;
    fits_update_key (finfo_out , TSTRING , "EVTFILE" , basename (file_out) , NULL , &status) ;
    printError (status , "Error in updating the key  value of the EVTFILE" , infofile_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the info out file" , file_out) ;
     writeUsrkeywordsFrmvect(infofile_out,key_record);
    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    fits_close_file (fevt_in , &status) ;
    printError (status , "Error in closing the input Event File" , file_out) ;
    return (EXIT_SUCCESS) ;
}

int uvtCentroidBias::getHistory (vector<string> &vhistory)
{
    int cnt=0; 
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;

    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir) ;

    if (clobber == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+" clobber=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+" clobber=no") ;
    if (history == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+" history=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+" history=no") ;
    vhistory.push_back ("Parameter List END") ;

    return (EXIT_SUCCESS) ;
}

int uvtCentroidBias::readcentroidbiasFile ()
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
    LOG(INFO) << "\nReading Centroid Bias file from calDB completed" << endl ;
    return (EXIT_SUCCESS) ;
}
