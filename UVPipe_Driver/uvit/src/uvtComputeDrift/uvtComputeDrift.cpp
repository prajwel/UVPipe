/* 
 * File:   uvtComputeDrift.cpp
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
#include <uvtComputeDrift.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<fft.h>
#include<spMatrix.h>
#include<spGeneral.h>
#include<glog/logging.h>
#include<transform.h>
#include <map>
#include <algorithm>
#include <iterator>

//#include "uvtComputeJitter.h"


//Constructor -called when object is created


uvtDriftComputation::uvtDriftComputation ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    strcpy (centroidDir , "Centroid") ;
}

//Destructor


uvtDriftComputation::~ uvtDriftComputation ()
{
    //freeMemory (centroidframelist , nframes , NAMESIZE) ;
}

//parameter File reading


int uvtDriftComputation::read (int argc , char** argv)
{
    int status = 0 ;
    status = readParams (argc , argv , 6 , FNAME , "inputdatadir" , inputdatadir ,
            BOOL,"GenMatchStarsFile_flag",&match_starsfile_gen,
            REAL4 , "error_per" , &err_per ,
            REAL4 , "diff_Dist" , &diff_dist ,
            INT , "freqDomainFilter_Flag" , &FreqDomainFilter_Flag ,
            FNAME,"attFile",attFile
           // FNAME,"caldbdir",caldbDir
           // REAL4 , "scaleNUV" , &scaleNUV ,
           // REAL4 , "scaleVIS" , &scaleVIS ,
           // REAL4 , "shiftX_NUVtoVIS" , &Shift_X ,
           // REAL4 , "shiftY_NUVtoVIS" , &Shift_Y ,
            //REAL4 , "angle_NUVtoVIS" , &angle_NUVtoVIS
            ) ;
    
    if (status) return (EXIT_FAILURE) ;

    if (FreqDomainFilter_Flag == 0)
    {
        status = readParams (argc , argv , 1 , INT , "type_Filtering" , &type_Filtering) ;
        if (status) return (EXIT_FAILURE) ;
       
        if (type_Filtering == 1 || type_Filtering==2)
        {
            status = readParams (argc , argv , 1 , INT , "fitting_flag" , &fittingflag) ;
            if (status) return (EXIT_FAILURE) ;
            if (fittingflag)
            {
                status = readParams (argc , argv , 4 , INT , "order_pitch" , &orderpitch ,
                        INT , "order_yaw" , &orderyaw ,
                        INT , "order_roll" , &orderroll ,
                        REAL , "delta_Time" , &delta_time

                        ) ;
                if (status) return (EXIT_FAILURE) ;
            }
        }
        else if (type_Filtering == 0)
        {
            status = readParams (argc , argv , 1 , REAL , "freq_value" , &freqvalue) ;
            if (status) return (EXIT_FAILURE) ;

        }

    }
    else if (FreqDomainFilter_Flag == 1)
    {
        status = readParams (argc , argv , 1 , REAL , "freq_value" , &freqvalue) ;
        if (status) return (EXIT_FAILURE) ;
    }

    status = readParams (argc , argv , 7 , INT , "shiftRotDetAlgoFlag" , &option_LeastSquare ,BOOL,"flagTheta",&flag_theta,
            INT,"algoFlag",&algo_flag_value,
            FNAME , "outdir" , outdir ,
            BOOL , "clobber" , &clobber ,
            BOOL , "history" , &history ,
            STRING , "mode" , mode) ;
    if (status) return (EXIT_FAILURE) ;

    return (EXIT_SUCCESS) ;
}


int uvtDriftComputation::read (char *inputdatadir  ,char * attFile_in,float percent , float diff_Dist , int freqDomainFilter_Flag , int type_filter ,
        double d_time , double freqvalue , int fitting_flag , int orderpitch , int orderyaw , int orderroll ,
        char *outdir , int op_leastsquare , int flag_matchstars, int algoFlag,int  flag_thetacomp,int clobber ,int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->attFile , attFile_in) ;
    strcpy (this->outdir , outdir) ;
    
    this->match_starsfile_gen=flag_matchstars;
    this->err_per = percent ;
    this->clobber = clobber ;
    this->FreqDomainFilter_Flag = freqDomainFilter_Flag ;
    this->delta_time = d_time ;
    this->diff_dist = diff_Dist ;
    this->type_Filtering = type_filter ;
    this->freqvalue = freqvalue ;
    this->algo_flag_value=algoFlag;
    //LOG(INFO)<<"Algo frl" <<algo_flag_value<<" "<<algoFlag;exit(1);
           this->flag_theta=flag_thetacomp; 
    if (fitting_flag)
    {
        this->orderpitch = orderpitch ;
        this->orderyaw = orderyaw ;
        this->orderroll = orderroll ;
    }
    this->fittingflag = fitting_flag ;
    this->option_LeastSquare = op_leastsquare ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}


void uvtDriftComputation::display ()
{

    LOG (INFO) << endl << "----------Display of parameters for Drift Computation  ---------" ;
    LOG (INFO) << endl << "Input Frame List Directory                        : " << inputdatadir ;
    LOG (INFO) << endl << "percentage Error allowed                          : " << err_per ;
    LOG (INFO) << endl << "minimum starting distance for matching stars of two frames                               : " << diff_dist ;
    LOG (INFO) << endl << "freqDomainFilter_flag                               : " << FreqDomainFilter_Flag ;
   
    
    if (FreqDomainFilter_Flag==1)
    {
        LOG (INFO) << endl << "Cut-off frequency                               : " << freqvalue ;
    }
    else if (FreqDomainFilter_Flag == 0)
    {
        LOG (INFO) << endl << "Type Filtering                               : " << type_Filtering ;
        if (type_Filtering == 1 || type_Filtering ==2)
        {
            LOG (INFO) << endl << "Fitting Flag                               : " << fittingflag ;
            LOG (INFO) << endl << "Order pitch                               : " << orderpitch ;
            LOG (INFO) << endl << "Order Yaw                               : " << orderyaw ;
            LOG (INFO) << endl << "Order Roll                               : " << orderroll ;
            LOG (INFO) << endl << "Smoothaning Duration                                : " << delta_time ;
        }
        else
        {
            LOG (INFO) << endl << "Cut-off frequency                               : " << freqvalue ;

        }
    }
    LOG (INFO) << endl << "Shift and Rotate algorithm                               : " << option_LeastSquare ;
    LOG (INFO) << endl << "Output Directory                               : " << outdir ;
    if (clobber == YES)
        LOG (INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG (INFO) << endl << "Overwrite                                         : NO" ;
    if (history == YES)
        LOG (INFO) << endl << "History                                             : YES" ;
    else
        LOG (INFO) << endl << "History                                              : NO" ;

    LOG (INFO) << endl << "----------Display of parameters for Drift Computation  Ends----------\n" ;

}

// Drift Computation


int uvtDriftComputation::uvtDriftComputationProcess ()
{
    LOG (INFO) << endl << "Drift Computation process started" << endl ;
   
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG (INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;
    diff_dist_cp=diff_dist;
    //check  existence  of output directory if output  exist and 
    //clobber =YES than remove that directory and recreate it,if exist and clobber =NO than exit from module.
    //if directory is not exist than create it .
    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;
    
    nframes = 0 ;
    LOG (INFO) << endl << moduleoutdir << "  directory created" ;

    string  tempfilepath = searchFile (inputdatadir , ".info") ;
     if (tempfilepath =="")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;

    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (! (FileExists (infofile_in)))
    {
        LOG (INFO) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY" ;
        return (EXIT_FAILURE) ;        
    }

    LOG (INFO) << endl << "\nInformation File :" << infofile_in ;
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the information file" , infofile_in) ;
     copyUsrkeywrdsTovect(finfo_in,key_records); 
    vector<string> track_keyRecord;
    track_keyRecord= key_records;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in Moving the 2nd HDU" , infofile_in) ;
   
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    if (xsize <= 0 || ysize <= 0)
    {
        LOG (ERROR) << endl << "Invalid xsize/ysize\n" ;
        return (EXIT_FAILURE) ;
    }
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file

    fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
    printError (status , "NFILES not Found" , infofile_in) ;

    centroidframelist = allocateMemory<char>(nframes , NAMESIZE) ; //for  storing the  input frame list.

    fits_read_col (finfo_in , TSTRING ,1 , 1 , 1 , nframes , NULL , (void *) centroidframelist , NULL , &status) ;
    printError (status , "Error in opening the out  information file" , infofile_in) ;
    
    fits_close_file(finfo_in,&status);
    printError (status , "Error in closing the file" , infofile_in) ;

    //setting path for output information file
    sprintf (infofile_out , "%s/%s_dr.info" , moduleoutdir , nameprefix) ;
    LOG (INFO) << endl << "\nOutput information File :" << infofile_out ;

    //creating  output information file
    LOG(INFO)<<"INformation befor status flag" <<status;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information file" , infofile_out) ;
    
    char *ttype[] = {"DUMMY"} ;
    char *tform[] = {"A1"} ;
    
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table" , infofile_out) ;

    datainfo.write (finfo_out) ; //writing basic data information
    
    /*----Info file creating completed, rest of the information will be put by other functions-----------*/
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "NAMEPRFX keyword  not updated/not Found" , infofile_in) ; //for creating name for output information file

    //check for mode 
    if (! (datainfo.getModeFlag () == IM || datainfo.getModeFlag () == PC))
    {
        LOG (ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG (ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
    //LOG(INFO)<<algo_flag_value;exit(1);
    if(algo_flag_value==1  ){
       
     //if (findDrift_WithAlgo_1 ())    return (EXIT_FAILURE) ;//find drift between  frames from given input directory
    }
    else if (algo_flag_value==2 ||algo_flag_value==3  || algo_flag_value==4){
        // LOG(INFO)<<"HHH";exit(1);
        if (findDrift_WithAlgo_3 ())    return (EXIT_FAILURE) ;//find drift between  frames from given input directory
    }
    else{
        LOG(ERROR)<<"***Invalid algo flag value***";
        return(EXIT_FAILURE);
    }
   
    
       writeDrift () ; //writes drift data to output file

    
    //write user keywords from vector to drift  file and  output information file
       fits_open_file (&finfo_in , infile_drift , READWRITE , &status) ;
       printError (status , "Error in opening the information file" , infofile_in) ;
   
      
    writeUsrkeywordsFrmvect (infile_drift , key_records) ;
    writeUsrkeywordsFrmvect (infofile_out , key_records) ;
     writeUsrkeywordsFrmvect (infile_drift , track_keyRecord) ;
    writeUsrkeywordsFrmvect (infofile_out , track_keyRecord) ;

    getHistory (vhistorystr) ; //reading history
    if (history == YES)
    {
        //writing history to  output information file
        writeHistory (infofile_out , vhistorystr) ;
        writeHistory (infile_drift , vhistorystr) ;
    }
    
   
    writeCommonKeywords (finfo_out , modulename) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the File " , infofile_out) ;
    
//   fits_open_file (&finfo_in , infile_drift , READWRITE , &status) ;
//   printError (status , "Error in opening the information file" , infofile_in) ;
 
   writeCommonKeywords (finfo_in, modulename) ;
   
   fits_close_file (finfo_in , &status) ;
   printError (status , "Error in closing the File " , infofile_out) ;
   
   return (EXIT_SUCCESS) ;
}


//int uvtDriftComputation::findDrift ()
//{
//     cntglobal=0;
//    LOG (INFO) << endl << "\nStarted Drift Calculation process..." ;
//    int status = 0 ;
//    long fpixel[2] ;
//    fpixel[0] = 1 , fpixel[1] = 1 ;
//    char errstr[512] ;
//    long naxes[2] ;
//    naxes[0] = naxes[1] = xsize ;
//    double frametime = 0 ;
//    key_records.clear () ;
//    int power = 1 ;
//    while (power < nframes - 1)
//    {
//        power = power * 2 ;
//    }
//    nframes_power2 = power ;    
//    time = new double[nframes_power2] ;
//    Xshift_arr = new double[nframes_power2] ;
//    Yshift_arr = new double[nframes_power2] ;
//    theta_arrfinal = new double[nframes_power2] ;
//    double cumm_x = 0.0 , cumm_y = 0.0 , cumm_theta = 0.0 ;
//    
//    for (int tt = 0 ; tt < nframes_power2 ; tt ++)
//    {
//        time[tt] = 0.0 ;
//        Xshift_arr[tt] = 0.0 ;
//        Yshift_arr[tt] = 0.0 ;
//        theta_arrfinal[tt] = 0.0 ;
//    }
//    
//    long numrows ;
//    LOG (INFO) << "\nTotal number of frames " << nframes << endl ;
//
//    char infile[FLEN_FILENAME] ;
//    char infile1[FLEN_FILENAME] ;
//
//    fitsfile *fptr , *fptr2 ;
//   
//
//    Xshift_arr[0] = 0.0 ;
//    Yshift_arr[0] = 0.0 ;
//    theta_arrfinal[0] = 0.0 ;
//
//    div_fact = xsize / 600 ;
//    int num_final = 0 ;
//    int num_arrSize = 0 ;
//  
//    int cnt = 0 ;
//
//    int numrows_firstFile = 0 ;
//    
//    //opening first file for reading centroid 
//    sprintf (infile , "%s/%s/%s" , inputdatadir , centroidDir , centroidframelist[0]) ;
//    fits_open_file (&fptr , infile , READONLY , &status) ;
//    printError (status , "***Error in opening file***" , infile) ;
//    
//    numrows = 0 ;
//    fits_movabs_hdu (fptr , 2 , NULL , &status) ;
//    printError (status , "***Error in reading HDU 2***" , infile) ;
//    
//    fits_get_num_rows (fptr , &numrows , &status) ;
//    printError (status , "***Error in getting number of rows***" , infile) ;
//
//    fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
//    printError (status , "***FRMTIME keyword not Found***") ;
//
//    //storing time of first frame
//    time[0] = frametime ;
//   
//    float *Xloc = new float[numrows] ;
//    float *Yloc = new float[numrows] ;
//    float *Intensity = new float[numrows] ;
//    numrows_firstFile = numrows ;
//   
//     fits_movabs_hdu (fptr , 1 , NULL , &status) ;
//    printError (status , "***Moving to  1st HDU Fails***") ;
//    
//    //copying user keywords  from first  HDU to vector
//    copyUsrkeywrdsTovect (fptr , key_records) ;
//
//    fits_close_file (fptr , &status) ;
//    printError (status , "Error in closing the input file",infile) ;
//    
//
//    //reading centroids from first file
//      status=readColumnsFromFITS (infile,2,3,TFLOAT,1,Xloc,numrows,TFLOAT,2,Yloc,numrows,
//            TFLOAT,3,Intensity,numrows);
//        if(status)
//         {
//        LOG(INFO)<<"Error reading column from  1st  centroid file"<<endl;
//        return(EXIT_FAILURE);
//        }
//   
//   LOG (INFO) << "\nFinding  similar Pixels from every two consecutive files " << endl ;
//
//    //loop For reading the number of frames in the input Directory
//    for (int p = 0 ; p < nframes - 1 ; p ++) 
//    {
//        LOG(INFO)<<"Pixel matching for "<<p+1<<" and "<<p+2<<" started...."<<endl;
//        sprintf (errstr , "Error at iteration number %d" , p) ;
//        sprintf (infile1 , "%s/%s/%s" , inputdatadir , centroidDir , centroidframelist[p + 1]) ;
//          
//        //opening second file for comparing
//        fits_open_file (&fptr2 , infile1 , READONLY , &status) ;
//        printError (status , "***input File reading Fails***") ;
//        fits_movabs_hdu (fptr2 , 2 , NULL , &status) ;
//        printError (status , "***Moving to  2nd HDU Fails***") ;
//        fits_get_num_rows (fptr2 , &numrows , &status) ;
//        int numrows_secondFile = numrows ;
//        fits_read_key (fptr2 , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
//        printError (status , "***FRMTIME keyword not Found***") ;
//
//         time[p + 1] = frametime ;
//       
//         //for storing the  centroid of the  file
//        float *Xloc1 = new float[numrows] ;
//        float *Yloc1 = new float[numrows] ;
//        float *Intensity1 = new float[numrows] ;
//         fits_close_file (fptr2 , &status) ;
//        
//         //reading the second file 's centroid to Xloc,Yloc and Intensity 
//         status=readColumnsFromFITS (infile1,2,3,TFLOAT,1,Xloc1,numrows,TFLOAT,2,Yloc1,numrows,
//            TFLOAT,3,Intensity1,numrows);
//        if(status)
//        {
//        LOG(INFO)<<"Error reading column from file"<<endl;
//        return(EXIT_FAILURE);
//        }
//
//        //number of match point initialize with zero
//        cnt = 0 ;
//
//        //decide  how many stars to be detected  for matching
//        num_final = numrows_firstFile < numrows_secondFile ? numrows_firstFile : numrows_secondFile ;
//
//        //calculating final array size for star matching
//        num_arrSize = numrows_firstFile > numrows_secondFile ? numrows_firstFile : numrows_secondFile ;
//        
//       vector<float> x_ref_arr, y_ref_arr , x_arr , y_arr ;
//        vector<float> temp_x_arr , temp_y_arr ;
//        
//        /**
//         * err_per-percentage of error allowed (user input)  
//         */
//        minimum_noTargetStars = (int) ((100 - err_per) * num_final / 100) ;
//        //minimum_noTargetStars=100;
//         cntglobal++;   
//       
//        // cout<<minimum_noTargetStars<<endl;
//    //older ones     
//         while(cnt < minimum_noTargetStars)
//        {
//             //cout<<diff_dist<<" "<<xsize<<" "<<div_fact<<endl;
//            x_ref_arr.clear ();y_ref_arr.clear ();x_arr.clear ();
//            y_arr.clear ();temp_x_arr.clear ();temp_y_arr.clear ();
//        
//           LOG(INFO)<<"Process is  repeating with newer value of  distance for  matching stars. New NH  : "<<diff_dist<<endl;
//                cnt = matchStars (numrows_firstFile , numrows_secondFile , div_fact , Xloc , Yloc , Xloc1 , Yloc1 , x_ref_arr , y_ref_arr , x_arr ,
//                y_arr , temp_x_arr , temp_y_arr) ;
//               // cout<<cnt<<endl;
//        diff_dist = diff_dist * 2 ;
//        
//        }
//       
//        // till this
//        
//       
//        diff_dist=diff_dist_cp;
//        
//        
//        double x_dx = 0.0 , y_dy = 0.0 , theta_dt = 0.0 ;
//
//        /* option_LeastSquare parameter decide which algorithm to use for finding the Drifts between  points 
//       these are  the different techniques for finding shifts between two point of two frames.
//         */
//       
//        status = findShiftsNtheta (cnt , x_ref_arr , y_ref_arr , Intensity,x_arr , y_arr ,Intensity1, temp_x_arr , temp_y_arr , x_dx , y_dy , theta_dt) ;
//        if (status)
//        {
//            LOG (INFO) << "Error in finding shifts n theta " << endl ;
//            return (EXIT_FAILURE) ;
//        }
//
//        cumm_x = cumm_x + x_dx ;
//        cumm_y = cumm_y + y_dy ;
//        cumm_theta = cumm_theta + theta_dt ;
//        
//        Xshift_arr[p + 1] = cumm_x; //cumulative X_shifts
//        Yshift_arr[p + 1] = cumm_y; //cumulative Y_shifts
//        theta_arrfinal[p + 1] = cumm_theta ; //cumulative theta
//   
//        delete[] Xloc , Yloc , Intensity ;
//
//        // storing current frame'ss X,Y,theta to  respective array  for the comparison of next frame in  next iteration .
//        Xloc = new float[numrows_secondFile] ;
//        Yloc = new float[numrows_secondFile] ;
//        Intensity = new float[numrows_secondFile] ;
//        for (int index = 0 ; index < numrows_secondFile ; index ++)
//        {
//            Xloc[index] = Xloc1[index] ;
//            Yloc[index] = Yloc1[index] ;
//            Intensity[index] = Intensity1[index] ;
//        }
//        delete[] Xloc1 , Yloc1 , Intensity1 ;
//
//        //setting number of rows of current frame  for next iteration
//        numrows_firstFile = numrows_secondFile ;
//
//    }
//  
//    status = doFiltering (x_arr , y_arr , theta_arr) ; //method for Applying filtering on the cumulative contents
//    if (status)
//    {
//        LOG (ERROR) << "***Filtering failed***" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//
//    return (EXIT_SUCCESS) ;
//}


//int uvtDriftComputation::findDrift_WithAlgo_1 ()
//{//LOG(INFO)<<"JJJJ";exit(1);
//     cntglobal=0;
//    LOG (INFO) << endl << "\nStarted Drift Calculation process..." ;
//    int status = 0 ;
//    long fpixel[2] ;
//    fpixel[0] = 1 , fpixel[1] = 1 ;
//    char errstr[512] ;
//    long naxes[2] ;
//    naxes[0] = naxes[1] = xsize ;
//    double frametime = 0 ;
//   // key_records.clear () ;
//    int power = 1 ;
//    while (power < nframes - 1)
//    {
//        power = power * 2 ;
//    }
//    nframes_power2 = power ;    
//    time = new double[nframes_power2] ;
//    Xshift_arr = new double[nframes_power2] ;
//    Yshift_arr = new double[nframes_power2] ;
//    theta_arrfinal = new double[nframes_power2] ;
//    double cumm_x = 0.0 , cumm_y = 0.0 , cumm_theta = 0.0 ;
//    //added
//    
//    
//     LOG(INFO)<<"Reading attitude data from "<<attFile;
//    
//    status=0;
//    
//    fitsfile *fatt;
//    fits_open_file(&fatt,attFile,READONLY,&status);
//    printError(status,"Error opening attitude file",attFile);
//    
//    fits_movabs_hdu(fatt,2,NULL,&status);                           //Attitude data is in 2nd HDU
//    printError(status,"Error in moving to 2nd HDU ", attFile);
//    
//    long nrows;
//    fits_get_num_rows(fatt,&nrows,&status);
//    printError(status,"Error reading number of rows from attitude file  ",attFile);
//        
//    long firstrow=1;
//    long firstelem=1;
//    long nelements=nrows;                  //for time
//    double *time_att = new double[nelements];
//    int tcoln,qcoln;               //column number for time and quaternion
//    
//     fits_read_col(fatt,TDOUBLE,1, firstrow,firstelem, nelements, NULL,time_att,NULL,&status); 
//     printError(status,"Error reading time column from attitude file  ",attFile);       
//         
//                      //number of elements to be read for quaternions
//   
//    long index_attFile=0;
//    for (int i=1;i<nelements;i++)
//    {
//        if(time[0]>time_att[i-1] && time[0]<time_att[i]){
//            index_attFile=i-1;
//            break;
//        }      
//        
//    }
//     
//   nelements =4;
//   double *qs = new double[nelements];
//   //quaternion reading 
//    fits_read_col(fatt,TDOUBLE,2,1,1, nelements, NULL,(void*)qs,NULL,&status); 
//    printError(status,"Error reading quaternions from attitude file  ",attFile);
//     
//    fits_close_file(fatt,&status);
//    printError(status,"Error in closing the file ",attFile);     
//     
//    double Roll_frmQua,Pitch_frmQua,yaw_frmQua;
//    
//    compute_Rotmatrix (qs[3],qs[0],qs[1],qs[2],Roll_frmQua,Pitch_frmQua,yaw_frmQua);
//    
//     if (strcasecmp (datainfo.getDetector () , "NUV") == 0)
//     {
//         transformRPYtoDXDYDTHETA_NUV (Roll_frmQua,Pitch_frmQua,yaw_frmQua,cumm_x,cumm_y,cumm_theta);
//         
//     }
//    else if (strcasecmp (datainfo.getDetector () , "FUV") == 0)
//     {
//         transformRPYtoDXDYDTHETA_FUV (Roll_frmQua,Pitch_frmQua,yaw_frmQua,cumm_x,cumm_y,cumm_theta);
//    
//     }
//    else if (strcasecmp (datainfo.getDetector () , "VIS") == 0)
//     {
//      transformRPYtoDXDYDTHETA_VIS (Roll_frmQua,Pitch_frmQua,yaw_frmQua,cumm_x,cumm_y,cumm_theta);         
//     }
//    //LOG(INFO)<<cumm_x<<" "<<cumm_y<<" "<<cumm_theta;exit(1);
////     for (int tt = 0 ; tt < nframes_power2 ; tt ++)
////    {
////        time[tt] = 0.0 ;
////        Xshift_arr[tt] = cumm_x ;
////        Yshift_arr[tt] = cumm_y ;
////        theta_arrfinal[tt] = cumm_theta ;
////    }
//  //  x_arr[0]=cumm_x;
//  //  y_arr[0]=cumm_y;
//   // theta_arr[0]=cumm_theta;
//    
//    //till this
//    
//    
//    for (int tt = 0 ; tt < nframes_power2 ; tt ++)
//    {
//        time[tt] = 0.0 ;
//        Xshift_arr[tt] = 0.0 ;
//        Yshift_arr[tt] = 0.0 ;
//        theta_arrfinal[tt] = 0.0 ;
//    }
//    
//    long numrows ;
//    LOG (INFO) << "\nTotal number of frames " << nframes << endl ;
//
//    char infile[FLEN_FILENAME] ;
//    char infile1[FLEN_FILENAME] ;
//
//    fitsfile *fptr , *fptr2 ;
//   
//    LOG(INFO)<<cumm_x<<" "<<cumm_y<<" "<<cumm_theta;
//    Xshift_arr[0] = cumm_x ;
//    Yshift_arr[0] = cumm_y ;
//    theta_arrfinal[0] = cumm_theta ;
//
//    div_fact = xsize / 600 ;
//    int num_final = 0 ;
//    int num_arrSize = 0 ;
//  
//    int cnt = 0 ;
//
//    int numrows_firstFile = 0 ;
//    
//    //opening first file for reading centroid 
//    sprintf (infile , "%s/%s/%s" , inputdatadir , centroidDir , centroidframelist[0]) ;
//    fits_open_file (&fptr , infile , READONLY , &status) ;
//    printError (status , "***Error in opening file***" , infile) ;
//    
//    numrows = 0 ;
//    fits_movabs_hdu (fptr , 2 , NULL , &status) ;
//    printError (status , "***Error in reading HDU 2***" , infile) ;
//    
//    fits_get_num_rows (fptr , &numrows , &status) ;
//    printError (status , "***Error in getting number of rows***" , infile) ;
//
//    fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
//    printError (status , "***FRMTIME keyword not Found***") ;
//
//    //storing time of first frame
//    time[0] = frametime ;
//   
//    float *Xloc = new float[numrows] ;
//    float *Yloc = new float[numrows] ;
//    float *Intensity = new float[numrows] ;
//    numrows_firstFile = numrows ;
//   
//     fits_movabs_hdu (fptr , 1 , NULL , &status) ;
//    printError (status , "***Moving to  1st HDU Fails***") ;
//    
//    //copying user keywords  from first  HDU to vector
//    copyUsrkeywrdsTovect (fptr , key_records) ;
//
//    fits_close_file (fptr , &status) ;
//    printError (status , "Error in closing the input file",infile) ;
//    
//
//    //reading centroids from first file
//      status=readColumnsFromFITS (infile,2,3,TFLOAT,1,Xloc,numrows,TFLOAT,2,Yloc,numrows,
//            TFLOAT,3,Intensity,numrows);
//        if(status)
//         {
//        LOG(INFO)<<"Error reading column from  1st  centroid file"<<endl;
//        return(EXIT_FAILURE);
//        }
//   
//   LOG (INFO) << "\nFinding  similar Pixels from every two consecutive files " << endl ;
//   char temp_path_txt_temp[FLEN_FILENAME];
//   ofstream ofptr1;
//   ofstream of1("Difference.txt");
//   double minvalx,maxvalx,minvaly,maxvaly;
//    //loop For reading the number of frames in the input Directory
//    for (int p = 0 ; p < nframes - 1 ; p ++) 
//    {
//        LOG(INFO)<<"Pixel matching for "<<p+1<<" and "<<p+2<<" started...."<<endl;
//        sprintf (errstr , "Error at iteration number %d" , p) ;
//        sprintf (infile1 , "%s/%s/%s" , inputdatadir , centroidDir , centroidframelist[p + 1]) ;
//          
//        //opening second file for comparing
//        fits_open_file (&fptr2 , infile1 , READONLY , &status) ;
//        printError (status , "***input File reading Fails***") ;
//        fits_movabs_hdu (fptr2 , 2 , NULL , &status) ;
//        printError (status , "***Moving to  2nd HDU Fails***") ;
//        fits_get_num_rows (fptr2 , &numrows , &status) ;
//        int numrows_secondFile = numrows ;
//        fits_read_key (fptr2 , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
//        printError (status , "***FRMTIME keyword not Found***") ;
//
//         time[p + 1] = frametime ;
//       
//         //for storing the  centroid of the  file
//        float *Xloc1 = new float[numrows] ;
//        float *Yloc1 = new float[numrows] ;
//        float *Intensity1 = new float[numrows] ;
//         fits_close_file (fptr2 , &status) ;
//        
//         //reading the second file 's centroid to Xloc,Yloc and Intensity 
//         status=readColumnsFromFITS (infile1,2,3,TFLOAT,1,Xloc1,numrows,TFLOAT,2,Yloc1,numrows,
//            TFLOAT,3,Intensity1,numrows);
//        if(status)
//        {
//        LOG(INFO)<<"Error reading column from file"<<endl;
//        return(EXIT_FAILURE);
//        }
//
//        //number of match point initialize with zero
//        cnt = 0 ;
//
//        //decide  how many stars to be detected  for matching
//        num_final = numrows_firstFile < numrows_secondFile ? numrows_firstFile : numrows_secondFile ;
//
//        //calculating final array size for star matching
//        num_arrSize = numrows_firstFile > numrows_secondFile ? numrows_firstFile : numrows_secondFile ;
//        
//       vector<float> x_ref_arr, y_ref_arr , x_arr , y_arr,int_ref_vect,int_vect ;
//        vector<float> temp_x_arr , temp_y_arr ;
//        
//        /**
//         * err_per-percentage of error allowed (user input)  
//         */
//        //minimum_noTargetStars = (int) ((100 - err_per) * num_final / 100) ;
//        //minimum_noTargetStars=100;
//         cntglobal++;   
//       
//        // cout<<minimum_noTargetStars<<endl;
//    //older ones
////           while(cnt < minimum_noTargetStars)
////        {
////           //  cout<<diff_dist<<" "<<xsize<<" "<<div_fact<<endl;
////            x_ref_arr.clear ();y_ref_arr.clear ();x_arr.clear ();
////            y_arr.clear ();temp_x_arr.clear ();temp_y_arr.clear ();
////        
////           LOG(INFO)<<"Process is  repeating with newer value of  distance for  matching stars. New NH  : "<<diff_dist<<endl;
////           cnt = matchStars (numrows_firstFile , numrows_secondFile , div_fact , Xloc , Yloc , Xloc1 , Yloc1 , x_ref_arr , y_ref_arr , x_arr ,
////                y_arr , temp_x_arr , temp_y_arr) ;
////              //  cout<<cnt<<endl;
////        diff_dist = diff_dist * 2 ;
////        
////        }
//         
//         //remove for new algortihm
////            cnt=numrows;
//            for(int i=0;i<numrows_secondFile;i++)
//            {
//                if((Xloc[i]!=-9999 && Yloc[i]!=-9999) && ( Xloc1[i]!=-9999 &&  Yloc1[i] !=-9999) )
//                {
//                x_arr.push_back (Xloc1[i]);
//                y_arr.push_back (Yloc1[i]);
//                int_vect.push_back (Intensity[i]);
//                }
//               
//            }
//            for(int i=0;i<numrows_firstFile;i++)
//           {
//                 if((Xloc[i]!=-9999 && Yloc[i]!=-9999) && ( Xloc1[i]!=-9999 && Yloc1[i] !=-9999) ){
//             x_ref_arr.push_back (Xloc[i]) ;
//             y_ref_arr.push_back (Yloc[i]);
//              int_ref_vect.push_back (Intensity1[i]);
//                 }
//            }
//             for(int i=0;i<numrows_firstFile;i++)
//           {
//                  if((Xloc[i]!=-9999 && Yloc[i]!=-9999)  && ( Xloc1[i]!=-9999 &&  Yloc1[i] !=-9999) )
//                  {
//             temp_x_arr.push_back (Xloc1[i]-Xloc[i]) ;
//             temp_y_arr.push_back (Yloc1[i]-Yloc[i]);
//                  }
//            }
//            sprintf (temp_path_txt_temp , "%s/matchstars_%d.txt" , moduleoutdir,cntglobal) ;
//            ofptr1.open (temp_path_txt_temp , ios::out) ;
//    for (int i=0;i<x_ref_arr.size ();i++)
//    {
//        ofptr1<<x_ref_arr[i]<<setw(20)<<y_ref_arr[i]<<setw(20)<<setw(20)<<x_arr[i]<<setw(20)<<y_arr[i]<<setw (20)<<temp_x_arr[i]<<setw(20)<<temp_y_arr[i]<<endl;
//        
//    }
////          ofptr1.close ();
////            
//        // till this         
////        diff_dist=diff_dist_cp;
//        double x_dx = 0.0 , y_dy = 0.0 , theta_dt = 0.0 ;
//
//        /* option_LeastSquare parameter decide which algorithm to use for finding the Drifts between  points 
//       these are  the different techniques for finding shifts between two point of two frames.
//         */
//       // cnt=numrows;
//        //cout<<" ppp::"<<x_ref_arr.size ()<<" "<<x_arr.size ()<<" "<<temp_x_arr.size ()<<endl;
//        
//        status = findShiftsNtheta (x_ref_arr.size () , x_ref_arr , y_ref_arr , Intensity,x_arr , y_arr ,Intensity1, temp_x_arr , temp_y_arr ,flag_theta, x_dx , y_dy , theta_dt,minvalx,maxvalx,minvaly,maxvaly) ;
//        if (status)
//        {
//            LOG (INFO) << "Error in finding shifts n theta " << endl ;
//            return (EXIT_FAILURE) ;
//        }
//       
//        //cout<<x_dx<<" "<<y_dy<<" "<<theta_dt<<endl;
//        cumm_x = cumm_x + x_dx ;
//        cumm_y = cumm_y + y_dy ;
//        cumm_theta = cumm_theta + theta_dt ;
//        
//        Xshift_arr[p + 1] = cumm_x; //cumulative X_shifts
//        Yshift_arr[p + 1] = cumm_y; //cumulative Y_shifts
//        theta_arrfinal[p + 1] = cumm_theta ; //cumulative theta
//        LOG(INFO)<<cumm_x<<" "<<cumm_y<<" "<<cumm_theta<<endl;
//        delete[] Xloc , Yloc , Intensity ;
//       
//        // storing current frame'ss X,Y,theta to  respective array  for the comparison of next frame in  next iteration .
//        Xloc = new float[numrows_secondFile] ;
//        Yloc = new float[numrows_secondFile] ;
//        Intensity = new float[numrows_secondFile] ;
//        for (int index = 0 ; index < numrows_secondFile ; index ++)
//        {
//            Xloc[index] = Xloc1[index] ;
//            Yloc[index] = Yloc1[index] ;
//            Intensity[index] = Intensity1[index] ;
//        }
//        delete[] Xloc1 , Yloc1 , Intensity1 ;
//
//        //setting number of rows of current frame  for next iteration
//        numrows_firstFile = numrows_secondFile ;
//        ofptr1.close ();
//    }
//  of1.close();
//    status = doFiltering (x_arr , y_arr , theta_arr) ; //method for Applying filtering on the cumulative contents
//    if (status)
//    {
//        LOG (ERROR) << "***Filtering failed***" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//
//    return (EXIT_SUCCESS) ;
//}

int uvtDriftComputation::findDrift_WithAlgo_3 ()
{
     cntglobal=0;
    LOG (INFO) << endl << "\nStarted Drift Calculation process..." ;
    int status = 0 ;
    long fpixel[2] ;
    fpixel[0] = 1 , fpixel[1] = 1 ;
    char errstr[512] ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    double frametime = 0 ;
    //key_records.clear () ;
    int power = 1 ;
    while (power < nframes - 1)
    {
        power = power * 2 ;
    }
    nframes_power2 = power ;    
    time = new double[nframes_power2] ;
    Xshift_arr = new double[nframes_power2] ;
    Yshift_arr = new double[nframes_power2] ;
    theta_arrfinal = new double[nframes_power2] ;
    double cumm_x = 0.0 , cumm_y = 0.0 , cumm_theta = 0.0 ;
    double temp_cumm_x=0.0,temp_cumm_y=0.0,temp_cumm_theta=0.0f;
    //added 
    
    
    //till
    for (int tt = 0 ; tt < nframes_power2 ; tt ++)
    {
        time[tt] = 0.0 ;
        Xshift_arr[tt] = 0.0 ;
        Yshift_arr[tt] = 0.0 ;
        theta_arrfinal[tt] = 0.0 ;
    }
    
    long numrows ;
    LOG (INFO) << "\nTotal number of frames " << nframes << endl ;

    if(nframes==0){
        LOG(ERROR)<<"***Total number of available frames are 0**";
	LOG(ERROR)<<"CRASH TOO FEW FRAMES TO COMPUTE DRIFT (uvtComputeDrift.cpp)";
        return(EXIT_FAILURE);
    }
    char infile[FLEN_FILENAME] ;
    char infile1[FLEN_FILENAME] ;

  fitsfile *fptr , *fptr2 ;
  
    
    Xshift_arr[0] = 0.0 ;
    Yshift_arr[0] = 0.0 ;
    theta_arrfinal[0] = 0.0 ;

    div_fact = xsize / 600 ;
    int num_final = 0 ;
    int num_arrSize = 0 ;
  
    int cnt = 0 ;

    int numrows_firstFile = 0 ;
    int numrows_refFile=0;
    //opening first file for reading centroid 
    sprintf (infile , "%s/%s/%s" , inputdatadir , centroidDir , centroidframelist[0]) ;
    fits_open_file (&fptr , infile , READONLY , &status) ;
    printError (status , "***Error in opening file***" , infile) ;
    
    numrows = 0 ;
    fits_movabs_hdu (fptr , 2 , NULL , &status) ;
    printError (status , "***Error in reading HDU 2***" , infile) ;
    
    fits_get_num_rows (fptr , &numrows , &status) ;
    printError (status , "***Error in getting number of rows***" , infile) ;

    fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
    printError (status , "***FRMTIME keyword not Found***") ;

    //storing time of first frame
    time[0] = frametime ;
   
    float *Xloc = new float[numrows] ;
    float *Yloc = new float[numrows] ;
    float *Intensity = new float[numrows] ;
    
    float *Xloc_refframe_backup = new float[numrows] ;
    float *Yloc_refframe_backup = new float[numrows] ;
    float *Intensity_refframe_backup = new float[numrows] ;
    
    
    numrows_firstFile = numrows ;
    numrows_refFile=numrows;
    
     fits_movabs_hdu (fptr , 1 , NULL , &status) ;
    printError (status , "***Moving to  1st HDU Fails***") ;
    
    //copying user keywords  from first  HDU to vector
    copyUsrkeywrdsTovect (fptr , key_records) ;

    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the input file",infile) ;
    

    //reading centroids from first file
      status=readColumnsFromFITS (infile,2,3,TFLOAT,1,Xloc,numrows,TFLOAT,2,Yloc,numrows,
            TFLOAT,3,Intensity,numrows);
        if(status)
         {
        LOG(INFO)<<"Error reading column from  1st  centroid file"<<endl;
        return(EXIT_FAILURE);
        }
   //Added new 
            
      for (int i=0;i<numrows;i++)
    {
    Xloc_refframe_backup[i]=Xloc[i];
    Yloc_refframe_backup[i]=Yloc[i];
    Intensity_refframe_backup[i]=Intensity[i];
    }
   //till this
   LOG (INFO) << "\nFinding  similar Pixels from every two consecutive files " << endl ;
 
   ofstream ofptr1;
   
   LOG(INFO)<<"Reading attitude data from "<<attFile;
    
    status=0;
    
//    fitsfile *fatt;
//    fits_open_file(&fatt,attFile,READONLY,&status);
//    printError(status,"Error opening attitude file",attFile);
//    
//    fits_movabs_hdu(fatt,2,NULL,&status);                           //Attitude data is in 2nd HDU
//    printError(status,"Error in moving to 2nd HDU ", attFile);
//    
//    long nrows;
//    fits_get_num_rows(fatt,&nrows,&status);
//    printError(status,"Error reading number of rows from attitude file  ",attFile);
//        
//    long firstrow=1;
//    long firstelem=1;
//    long nelements=nrows;                  //for time
//    double *time_att = new double[nelements];
//    int tcoln,qcoln;               //column number for time and quaternion
//    
//     fits_read_col(fatt,TDOUBLE,1, firstrow,firstelem, nelements, NULL,time_att,NULL,&status); 
//     printError(status,"Error reading time column from attitude file  ",attFile);       
//         
//                      //number of elements to be read for quaternions
//   
//    long index_attFile=0;
// 
//    for (int i=1;i<nelements;i++)
//    {
//        if(time[0]>time_att[i-1] && time[0]<time_att[i]){
//            index_attFile=i-1;
//            break;
//        }      
//        
//    }
//    
//    if(index_attFile==0)
//    {
//        LOG(INFO)<<"***Not related record found for particular time attitude file***";
//        return(EXIT_FAILURE);
//    }
//   nelements =4;
//   double *qs = new double[nelements];
//   //quaternion reading 
//    fits_read_col(fatt,TDOUBLE,2,index_attFile+1,1, nelements, NULL,(void*)qs,NULL,&status); 
//    printError(status,"Error reading quaternions from attitude file  ",attFile);
//     
//    fits_close_file(fatt,&status);
//    printError(status,"Error in closing the file ",attFile);     
//     
//    double Roll_frmQua,Pitch_frmQua,yaw_frmQua;
//    //LOG(INFO)<<qs[3]<<" "<<qs[0]<<" "<<qs[1]<<" "<<qs[2];exit(1);
//    compute_Rotmatrix (qs[3],qs[0],qs[1],qs[2],Roll_frmQua,Pitch_frmQua,yaw_frmQua);
//    
//     if (strcasecmp (datainfo.getDetector () , "NUV") == 0)
//     {
//         transformRPYtoDXDYDTHETA_NUV (Roll_frmQua,Pitch_frmQua,yaw_frmQua,cumm_x,cumm_y,cumm_theta);
//         
//     }
//    else if (strcasecmp (datainfo.getDetector () , "FUV") == 0)
//     {
//         transformRPYtoDXDYDTHETA_FUV (Roll_frmQua,Pitch_frmQua,yaw_frmQua,cumm_x,cumm_y,cumm_theta);
//    
//     }
//    else if (strcasecmp (datainfo.getDetector () , "VIS") == 0)
//     {
//      transformRPYtoDXDYDTHETA_VIS (Roll_frmQua,Pitch_frmQua,yaw_frmQua,cumm_x,cumm_y,cumm_theta);         
//     }
    
 Xshift_arr[0] = cumm_x ;
  Yshift_arr[0] = cumm_y ;
  theta_arrfinal[0] = cumm_theta ;
   
  cumm_x=0.0f;
  cumm_y=0.0f;
  cumm_theta=0.0f;
   Xshift_arr[0] = 0 ;
 Yshift_arr[0] = 0 ;
  theta_arrfinal[0] = 0 ;
  double x_temp_dx,y_temp_dy,theta_temp_theta;
  bool flag_EnoughMatchpointsFound=FALSE;
  int lastFrame_cnt=0;
  int num_matchForExecution=0;
  bool flag_onlyOnePoint_firstTime=FALSE;
  int cnt_ForFirstTime=0;
  double minvalx,maxvalx,minvaly,maxvaly;
  ofstream of1;
  of1.open("Differences.txt",ios::out|ios::trunc);
  of1<<"==========================================================================================================================="<<endl;
  of1<<"FILENUMBER"<<setw(20)<<"MIN_X"<<setw(20)<<"MAX_X"<<setw(20)<<"MIN_Y"<<setw(20)<<"MAX_Y"<<setw(20)<<"RANGE_X"<<setw(20)<<"RANGE_Y"<<endl;
  of1<<"==========================================================================================================================="<<endl;
    //loop For reading the number of frames in the input Directory 
    for (int p = 0 ; p < nframes - 1 ; p ++) 
    {
        LOG(INFO)<<"Pixel matching for "<<p+1<<" and "<<p+2<<" started...."<<endl;
        sprintf (errstr , "Error at iteration number %d" , p) ;
        sprintf (infile1 , "%s/%s/%s" , inputdatadir , centroidDir , centroidframelist[p + 1]) ;
          
        //opening second file for comparing
        fits_open_file (&fptr2 , infile1 , READONLY , &status) ;
        printError (status , "***input File reading Fails***") ;
        fits_movabs_hdu (fptr2 , 2 , NULL , &status) ;
        printError (status , "***Moving to  2nd HDU Fails***") ;
        fits_get_num_rows (fptr2 , &numrows , &status) ;
        int numrows_secondFile = numrows ;
        fits_read_key (fptr2 , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
        printError (status , "***FRMTIME keyword not Found***") ;

         time[p + 1] = frametime ;
       
         //for storing the  centroid of the  file
        float *Xloc1 = new float[numrows] ;
        float *Yloc1 = new float[numrows] ;
        float *Intensity1 = new float[numrows] ;
        
        float *Xloc1_temp = new float[numrows] ;
        float *Yloc1_temp = new float[numrows] ;
        float *Intensity1_temp = new float[numrows] ;
         fits_close_file (fptr2 , &status) ;
        
         //reading the second file 's centroid to Xloc,Yloc and Intensity 
         status=readColumnsFromFITS (infile1,2,3,TFLOAT,1,Xloc1,numrows,TFLOAT,2,Yloc1,numrows,
            TFLOAT,3,Intensity1,numrows);
        if(status)
        {
        LOG(INFO)<<"Error reading column from file"<<endl;
        return(EXIT_FAILURE);
        }

        //number of match point initialize with zero
        cnt = 0 ;

        //decide  how many stars to be detected  for matching
        num_final = numrows_firstFile < numrows_secondFile ? numrows_firstFile : numrows_secondFile ;

        //calculating final array size for star matching
        num_arrSize = numrows_firstFile > numrows_secondFile ? numrows_firstFile : numrows_secondFile ;
        
       vector<float> x_ref_arr, y_ref_arr , x_arr , y_arr,int_ref_vect,int_vect ;
       vector<float> temp_x_arr , temp_y_arr ,New_X_ref,New_Y_ref,New_x_arr,New_y_arr,New_Xdiff,New_Ydiff;
       vector<float> int_new_one,int_new_two;
        /**
         * err_per-percentage of error allowed (user input)  
         */
     //   minimum_noTargetStars = (int) ((100 - err_per) * num_final / 100) ;
     //   minimum_noTargetStars=100;
         cntglobal++;   
       
        //LOG(INFO)<<minimum_noTargetStars<<" "<<num_final<<endl;
        
        //  older ones
         flag_EnoughMatchpointsFound=FALSE;
         if(datainfo.getModeFlag()==IM){
             if(flag_theta==0)
         num_matchForExecution=2;
             else
                 num_matchForExecution=3;
         }
         else{
             if(flag_theta==0)
             num_matchForExecution=2;
             else
                 num_matchForExecution=3;
         }
         flag_onlyOnePoint_firstTime=FALSE;
         cnt_ForFirstTime=0;
        while(cnt <num_matchForExecution )
       {
            
            if(diff_dist>16)
            {   cnt_ForFirstTime++;
                flag_EnoughMatchpointsFound=TRUE;
                if(cnt == 1 && flag_onlyOnePoint_firstTime == FALSE){
                    flag_onlyOnePoint_firstTime=TRUE;
                    cnt=0;
                    num_matchForExecution=1;
                    diff_dist=diff_dist_cp;
                    continue;
                }
                LOG(ERROR)<<"***Stars are not matching ***";
		LOG(ERROR)<<"CRASH NO MATCH OF STARS BETWEEN SUCCESSIVE FRAMES (uvtComputeDrift.cpp)";
                lastFrame_cnt=p;
                break;
            }
           
             //cout<<diff_dist<<" "<<xsize<<" "<<div_fact<<endl;
            x_ref_arr.clear ();y_ref_arr.clear ();x_arr.clear ();
            y_arr.clear ();temp_x_arr.clear ();temp_y_arr.clear ();
            int_ref_vect.clear(),int_vect.clear();
           LOG(INFO)<<"Process is  repeating with newer value of  distance for  matching stars. New NH  : "<<diff_dist<<" "<<cnt<<endl;
           cnt = matchStars (numrows_firstFile , numrows_secondFile , div_fact , Xloc , Yloc , 
                   Xloc1 , Yloc1 ,Intensity,Intensity1, x_ref_arr , y_ref_arr , x_arr ,
                y_arr ,int_ref_vect,int_vect, temp_x_arr , temp_y_arr) ;
                LOG(INFO)<<cnt;
              //  cout<<cnt<<endl;
        diff_dist = diff_dist * 2 ;        
       }        

        diff_dist=diff_dist_cp;
        double x_dx = 0.0 , y_dy = 0.0 , theta_dt = 0.0 ;

        /* option_LeastSquare parameter decide which algorithm to use for finding the Drifts between  points 
       these are  the different techniques for finding shifts between two point of two frames.
         */
        if(flag_onlyOnePoint_firstTime==TRUE && cnt_ForFirstTime==1){
           
           status = findShiftsNtheta (x_ref_arr.size () , x_ref_arr , y_ref_arr , int_ref_vect,x_arr , y_arr ,int_vect, New_Xdiff , New_Ydiff , flag_theta,x_dx , y_dy , theta_dt,minvalx,maxvalx,minvaly,maxvaly) ;
        if (status)
        {
            LOG (INFO) << "Error in finding shifts n theta " << endl ;
            return (EXIT_FAILURE) ;
        } 
        }
        
      else if(flag_EnoughMatchpointsFound==FALSE){
//         status = findShiftsNtheta (x_ref_arr.size () , x_ref_arr , y_ref_arr , Intensity,x_arr , y_arr ,Intensity1, New_Xdiff , New_Ydiff , flag_theta,x_dx , y_dy , theta_dt) ;
//        if (status)
//        {
//            LOG (INFO) << "Error in finding shifts n theta " << endl ;
//            return (EXIT_FAILURE) ;
//        } 
        status=removeRecords(x_ref_arr,y_ref_arr,x_arr,y_arr,temp_x_arr,temp_y_arr,int_ref_vect,int_vect,New_X_ref,New_Y_ref,New_x_arr,New_y_arr,New_Xdiff,New_Ydiff,int_new_one,int_new_two);
        if (status)
        {
            LOG (ERROR) << "Error in removing the records above mean+sigma" ;
            return (EXIT_FAILURE) ;
        }
         status = findShiftsNtheta (New_X_ref.size () , New_X_ref , New_Y_ref , int_new_one,New_x_arr , New_y_arr ,int_new_two, New_Xdiff , New_Ydiff , flag_theta,x_dx , y_dy , theta_dt,minvalx,maxvalx,minvaly,maxvaly) ;
        if (status)
        {
            LOG (INFO) << "Error in finding shifts n theta " << endl ;
            return (EXIT_FAILURE) ;
        }
     }
         else
        {
            LOG(INFO)<<"***Ignoring Last frame because not enough points are found";
            nframes=p+1;
            break;           
        }
//        status = findShiftsNtheta (x_ref_arr.size () , x_ref_arr , y_ref_arr , Intensity,x_arr , y_arr ,Intensity1, temp_x_arr , temp_y_arr , x_dx , y_dy , theta_dt) ;
//        if (status)
//        {
//            LOG (INFO) << "Error in finding shifts n theta " << endl ;
//            return (EXIT_FAILURE) ;
//        }
         
         //New logic introduced here.
        //of1<<"--------------------------------------------------------------------------------------------"<<endl;
        of1<<setw(20)<<p+1<<"-"<<p+2<<setw(20)<<minvalx<<setw(20)<<maxvalx<<setw(20)<<minvaly<<setw(20)<<maxvaly<<setw(20)<<maxvalx-minvalx<<setw(20)<<maxvaly-minvaly<<endl;
         LOG(INFO)<<x_dx<<" "<<y_dy<<" "<<theta_dt;
         x_temp_dx=x_dx;
         y_temp_dy=y_dy;
         theta_temp_theta=theta_dt;
         temp_cumm_x=temp_cumm_x+x_dx;
         temp_cumm_y=temp_cumm_y+y_dy;
         temp_cumm_theta=temp_cumm_theta+theta_dt;
         for(int index=0;index<numrows;index++)
         {
         x_temp_dx=(Xloc1[index]-xsize/2)*cos(-temp_cumm_theta*M_PI/180)-(Yloc1[index]-ysize/2)*sin(-temp_cumm_theta*M_PI/180)-temp_cumm_x+xsize/2;
             y_temp_dy=(Xloc1[index]-xsize/2)*sin(-temp_cumm_theta*M_PI/180)+(Yloc1[index]-ysize/2)*cos(-temp_cumm_theta*M_PI/180)-temp_cumm_y+ysize/2; 
//%#Added ON 20July#%  
             if(round(x_temp_dx)> 0 && round(x_temp_dx)<xsize && round(y_temp_dy)> 0 && round(y_temp_dy)<ysize )
//%#-Till this-20July17#%
{
             Xloc1_temp[index]=x_temp_dx;
             Yloc1_temp[index]=y_temp_dy;
             }
             else{
                 Xloc1_temp[index]=INVALID_PIX_VALUE;
                 Yloc1_temp[index]=INVALID_PIX_VALUE;
             }
             
         }
         cnt=0;
//         for (int i=0;i<numrows;i++){
//             LOG(INFO)<<Xloc_refframe_backup[i]<<" "<<Xloc1_temp[i]<<" "<<Yloc_refframe_backup[i]<<" "<<Yloc1_temp[i];
//         }
         flag_EnoughMatchpointsFound=FALSE;
         cnt_ForFirstTime=0;
          //num_matchForExecution=1;
         while(cnt <num_matchForExecution )
         {             
            
            if(diff_dist>4){
//                if(cnt==1){
//                    cnt=0;
//                    num_matchForExecution=1;
//                    diff_dist=diff_dist_cp;
//                    continue;                    
//                }
                cnt_ForFirstTime++;
                flag_EnoughMatchpointsFound=TRUE;
                LOG(ERROR)<<"***Not matched in second iteration ***"<<" "<<cnt;
                LOG(ERROR)<<"CRASH FAILED TO MATCH STARS WITH REFERENCE FRAME (uvtCompiteDrift.cpp)";
                break;
               // return(EXIT_FAILURE);
            }
             
             //cout<<diff_dist<<" "<<xsize<<" "<<div_fact<<endl;
            x_ref_arr.clear ();y_ref_arr.clear ();x_arr.clear ();
            y_arr.clear ();temp_x_arr.clear ();temp_y_arr.clear ();
            int_ref_vect.clear(),int_vect.clear();
           LOG(INFO)<<"Process is  repeating with newer value of  distance for  matching stars. New NH  : "<<diff_dist<<" "<<cnt;
           cnt = matchStars (numrows_refFile , numrows_secondFile , div_fact , Xloc_refframe_backup , Yloc_refframe_backup , Xloc1_temp , Yloc1_temp ,Intensity_refframe_backup,Intensity1 , x_ref_arr , y_ref_arr , x_arr ,
                y_arr ,int_ref_vect,int_vect, temp_x_arr , temp_y_arr) ;
              //  cout<<cnt<<endl;
            diff_dist = diff_dist * 2 ;        
            
        }        

        diff_dist=diff_dist_cp;
        
//         if((flag_onlyOnePoint_firstTime==TRUE && cnt_ForFirstTime==0) || cnt ==1){
         if( cnt ==1 || cnt ==2){
           status = findShiftsNtheta (x_ref_arr.size () , x_ref_arr , y_ref_arr ,int_ref_vect,x_arr , y_arr ,int_vect, New_Xdiff , New_Ydiff , flag_theta,x_dx , y_dy , theta_dt,minvalx,maxvalx,minvaly,maxvaly) ;
        if (status)
        {
            LOG (INFO) << "Error in finding shifts n theta " << endl ;
            return (EXIT_FAILURE) ;
        } 
        }
        
        
      else  if(flag_EnoughMatchpointsFound==FALSE){
         
        x_dx = 0.0 , y_dy = 0.0 , theta_dt = 0.0 ;
         status=removeRecords(x_ref_arr,y_ref_arr,x_arr,y_arr,temp_x_arr,temp_y_arr,int_ref_vect,int_vect,New_X_ref,New_Y_ref,New_x_arr,New_y_arr,New_Xdiff,New_Ydiff,int_new_one,int_new_two);
        if (status)
        {
            LOG (ERROR) << "Error in removing the records above mean+sigma" ;
            return (EXIT_FAILURE) ;
        }
         
         status = findShiftsNtheta (New_X_ref.size () , New_X_ref , New_Y_ref , int_new_one,New_x_arr , New_y_arr ,int_new_two, New_Xdiff , New_Ydiff , flag_theta,x_dx , y_dy , theta_dt,minvalx,maxvalx,minvaly,maxvaly) ;
        if (status)
        {
            LOG (INFO) << "Error in finding shifts n theta " << endl ;
            return (EXIT_FAILURE) ;
        }
        }
        else
        {
           
            nframes=p+1;
            break;           
        }
        
        temp_cumm_x=temp_cumm_x+x_dx;
        temp_cumm_y=temp_cumm_y+y_dy;
        temp_cumm_theta=temp_cumm_theta+theta_dt;
         //till this.
        LOG(INFO)<<x_dx<<" "<<y_dy<<" "<<theta_dt;
        cumm_x = temp_cumm_x;//+ x_dx-x_temp_dx ;
        cumm_y = temp_cumm_y ;//+ y_dy-y_temp_dy ;
        cumm_theta = temp_cumm_theta;// + theta_dt-theta_temp_theta ;
//        cumm_x = cumm_x + x_dx+x_temp_dx ;
//        cumm_y = cumm_y + y_dy+y_temp_dy ;
//        cumm_theta = cumm_theta + theta_dt ;
        
        Xshift_arr[p + 1] = cumm_x; //cumulative X_shifts
        Yshift_arr[p + 1] = cumm_y; //cumulative Y_shifts
        theta_arrfinal[p + 1] = cumm_theta ; //cumulative theta
        LOG(INFO)<<cumm_x<<" "<<cumm_y<<" "<<cumm_theta<<endl;
        delete[] Xloc , Yloc , Intensity ;
       
        // storing current frame'ss X,Y,theta to  respective array  for the comparison of next frame in  next iteration .
        Xloc = new float[numrows_secondFile] ;
        Yloc = new float[numrows_secondFile] ;
        Intensity = new float[numrows_secondFile] ;
        for (int index = 0 ; index < numrows_secondFile ; index ++)
        {
            Xloc[index] = Xloc1[index] ;
            Yloc[index] = Yloc1[index] ;
            Intensity[index] = Intensity1[index] ;
        }
        delete[] Xloc1 , Yloc1 , Intensity1,Xloc1_temp,Yloc1_temp,Intensity1_temp ;

        //setting number of rows of current frame  for next iteration
        numrows_firstFile = numrows_secondFile ;
        ofptr1.close ();
        
                
    }
  of1.close();
    status = doFiltering (x_arr , y_arr , theta_arr) ; //method for Applying filtering on the cummulative contents
    if (status)
    {
        LOG (ERROR) << "***Filtering failed***" << endl ;
        return (EXIT_FAILURE) ;
    }

    return (EXIT_SUCCESS) ;
}





int uvtDriftComputation::getHistory (vector<string> &vhistory)
{
    vhistory.clear () ;
    int cnt=0;
   // char *user = getlogin () ;
   // string str = "Module run by " + (string) user ;
  //  vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;
    //    vhistory.push_back ("P2 caldbDir=" + (string) caldbDir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir) ;
    char temp[FLEN_FILENAME] ;
    sprintf (temp , "%d" , FreqDomainFilter_Flag) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Frequency Domain Filter Flag=" + (string) temp) ;
    if (FreqDomainFilter_Flag)
    {
        sprintf (temp , "%f" , freqvalue) ;
        vhistory.push_back ((string)getSerialNo (cnt)+" Frequency Value=" + (string) temp) ;

    }
    else if(FreqDomainFilter_Flag==0)
    {
        sprintf (temp , "%d" , type_Filtering) ;
        vhistory.push_back ((string)getSerialNo (cnt)+" Type Filtering used=" + (string) temp) ;
        if (type_Filtering)
        {
            sprintf (temp , "%d" , fittingflag) ;
            vhistory.push_back ((string)getSerialNo (cnt)+" Fitting Flag =" + (string) temp) ;
        }
        if (fittingflag == 1 && type_Filtering == 1 ||type_Filtering==2)
        {
            sprintf (temp , "%d" , orderpitch) ;
            vhistory.push_back ((string)getSerialNo (cnt)+" order pitch =" + (string) temp) ;
            sprintf (temp , "%d" , orderroll) ;
            vhistory.push_back ((string)getSerialNo (cnt)+" order roll =" + (string) temp) ;
            sprintf (temp , "%d" , orderyaw) ;
            vhistory.push_back ((string)getSerialNo (cnt)+" order yaw =" + (string) temp) ;
            sprintf (temp , "%d" , delta_time) ;
            vhistory.push_back ((string)getSerialNo (cnt)+" Delta time=" + (string) temp) ;

        }
    }
    sprintf (temp , "%d" , option_LeastSquare) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Drift Computation Algorithm=" + (string) temp) ;
    sprintf (temp , "%f" , diff_dist) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Minimum search distance for star match in succesive frames=" + (string) temp) ;
    // vhistory.push_back ("Centroid Directory=" + (string) centroidDir) ;
    if (clobber == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+" clobber=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+" clobber=no") ;
    if (history == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+" history=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+" history=no") ;
    vhistory.push_back ((string)getSerialNo (cnt)+"Parameter List END") ;
    return (EXIT_SUCCESS) ;
}


int uvtDriftComputation::doFiltering (vector<double> &xshiftvect , vector<double> &yshiftvect , vector<double> &thetavect)
{
    int status = 0 ;
    LOG (INFO) << "Filtering Started..." << endl ;
    
    xshiftvect.clear () ;
    yshiftvect.clear () ;
    thetavect.clear () ;
 
    //time domain filtering
    if (FreqDomainFilter_Flag == 0)
    {
        
        if (type_Filtering == 0) //spatial domain
        {

            LOG (INFO) << "\nSpatial domain Filtering is performing....." << endl ;
            
            double cutoffFreq = freqvalue ; //cut off frequency for low pass filter           ;
            double RC = 1 / (2 * M_PI * cutoffFreq) ; //Time period
            double delta_t = time[1] - time[0] ;
            double alpha = delta_t / (RC + delta_t) ;
            double *delta_x_out = new double[nframes] ;
            double *delta_y_out = new double[nframes] ;
            double *delta_theta_out = new double[nframes] ;
           
            delta_x_out[0] = Xshift_arr[0] ;
            delta_y_out[0] = Yshift_arr[0] ;
            delta_theta_out[0] = theta_arrfinal[0] ;
           
            char temp_path_txt_temp[FLEN_FILENAME] ;
            sprintf (temp_path_txt_temp , "%s/Delta_x_beforeLowpassFilter.txt" , moduleoutdir) ;
            ofstream ofptr ;
            ofptr.open (temp_path_txt_temp , ios::out) ;

            for (int i = 0 ; i < nframes ; i ++)
            {
                ofptr << time[i] << setw (20) << theta_arrfinal[i] << setw (20) << Xshift_arr[i] << setw (20) << Yshift_arr[i] << endl ;
            }


            ofptr.close () ;

            sprintf (temp_path_txt_temp , "%s/Delta_x_afterLowpassFilter.txt" , moduleoutdir) ; //file For storing the before FFT content.

            ofptr.open (temp_path_txt_temp , ios::out) ;

            ofptr << time[0] << setw (20) << delta_theta_out[0] << setw (20) << delta_x_out[0] << setw (20) << delta_y_out[0] << endl ;
            
            //algorithm of spatial domain filtering 
            for (int i = 1 ; i < nframes ; i ++)
            {
                delta_x_out[i] = alpha * Xshift_arr[i]+(1 - alpha) * delta_x_out[i - 1] ;
                delta_y_out[i] = alpha * Yshift_arr[i]+(1 - alpha) * delta_y_out[i - 1] ;
                delta_theta_out[i] = alpha * theta_arrfinal[i]+(1 - alpha) * delta_theta_out[i - 1] ;
                ofptr << time[i] << setw (20) << delta_theta_out[i] << setw (20) << delta_x_out[i] << setw (20) << delta_y_out[i] << endl ;
                ;
            }
            for (int i = 0 ; i < nframes ; i ++)
            {
                xshiftvect.push_back (delta_x_out[i]) ;
                yshiftvect.push_back (delta_y_out[i]) ;
                thetavect.push_back (delta_theta_out[i]) ;
            }

            ofptr.close () ;


        }
        else if (type_Filtering == 1)
        {       //polynomial Fitting
            
            if (fittingflag) //decide fitting to be done or not
            {
                LOG (INFO) << "\nPolynomial Fitting is performing..." << endl ;
                char temp_path_txt_temp[FLEN_FILENAME] ;
                sprintf (temp_path_txt_temp , "%s/Delta_x_beforelocalSmoothingFilter.txt" , moduleoutdir) ;
                ofstream ofptr ;
                ofptr.open (temp_path_txt_temp , ios::out) ;

                             
                for (int i = 0 ; i < nframes ; i ++)
                {

                    ofptr << time[i] << setw (20) << theta_arrfinal[i] << setw (20) << Xshift_arr[i] << setw (20) << Yshift_arr[i] << endl ;

                }
                ofptr.close () ;
                vector<double> time_temp , X_shift_vect , Y_shift_vect , theta_vect ;
                // cout<<"NFRAMES "<<nframes<<endl;//exit(1);
               
                for (int i = 0 ; i < nframes ; i ++)
                {
                    X_shift_vect.push_back (Xshift_arr[i]) ;
                    Y_shift_vect.push_back (Yshift_arr[i]) ;
                    theta_vect.push_back (theta_arrfinal[i]) ;
                    time_temp.push_back (time[i]) ;
                }
                LOG (INFO) << "\nSmoothening duration : " << delta_time << endl ;
               
                //method for applying smoothening over 'delta_time' duration
                status = local_Smoothening (time_temp , X_shift_vect ,orderpitch, delta_time ) ;
                if (status)
                {
                    LOG (ERROR) << "Error in local smoothening the X_shift" << endl ;
                    return (EXIT_FAILURE) ;
                }
   

                status = local_Smoothening (time_temp , Y_shift_vect , orderyaw , delta_time) ;
                if (status)
                {
                    LOG (ERROR) << "Error in local smoothening the y_shift" << endl ;
                    return (EXIT_FAILURE) ;
                }

                status = local_Smoothening (time_temp , theta_vect , orderroll , delta_time ) ;
                if (status)
                {
                    LOG (ERROR) << "Error in local smoothening the theta_shift" << endl ;
                    return (EXIT_FAILURE) ;
                }

                sprintf (temp_path_txt_temp , "%s/Delta_x_afterlocalSmoothingFilter.txt" , moduleoutdir) ;

                ofptr.open (temp_path_txt_temp , ios::out) ;

               
                delete[] Xshift_arr , Yshift_arr , theta_arrfinal ;
                Xshift_arr = new double[X_shift_vect.size ()] ;
                Yshift_arr = new double[X_shift_vect.size ()] ;
                theta_arrfinal = new double[X_shift_vect.size ()] ;
               // cout<<X_shift_vect.size ()<<endl;
                
                for (int i = 0 ; i < X_shift_vect.size () ; i ++)
                {
                    Xshift_arr[i] = X_shift_vect[i] ;
                    Yshift_arr[i] = Y_shift_vect[i] ;
                    theta_arrfinal[i] = theta_vect[i] ;
                    time[i] = time_temp[i] ;
                }

                for (int i = 0 ; i < X_shift_vect.size () ; i ++)
                {

                    ofptr << time[i] << setw (20) << theta_arrfinal[i] << setw (20) << Xshift_arr[i] << setw (20) << Yshift_arr[i] << endl ;

                }
                ofptr.close () ;
                nframes = X_shift_vect.size () ;                 
                xshiftvect.clear ();   yshiftvect.clear ();    theta_vect.clear ();   time_temp.clear ();
               
            }

          //writing output to final vectors
            for (int i = 0 ; i < nframes ; i ++)
            {
                xshiftvect.push_back (Xshift_arr[i]) ;
                yshiftvect.push_back (Yshift_arr[i]) ;
                thetavect.push_back (theta_arrfinal[i]) ;
            }

        }
         else if (type_Filtering == 2)
        {       //polynomial Fitting for 
            
            if (fittingflag) //decide fitting to be done or not
            {
                LOG (INFO) << "\nPolynomial Fitting is performing..." << endl ;
                char temp_path_txt_temp[FLEN_FILENAME] ;
                sprintf (temp_path_txt_temp , "%s/Delta_x_beforelocalSmoothingFilter.txt" , moduleoutdir) ;
                ofstream ofptr ;
                ofptr.open (temp_path_txt_temp , ios::out) ;

                             
                for (int i = 0 ; i < nframes ; i ++)
                {

                    ofptr << time[i] << setw (20) << theta_arrfinal[i] << setw (20) << Xshift_arr[i] << setw (20) << Yshift_arr[i] << endl ;

                }
                ofptr.close () ;
                vector<double> time_temp , X_shift_vect , Y_shift_vect , theta_vect ;
                // cout<<"NFRAMES "<<nframes<<endl;//exit(1);
               
                for (int i = 0 ; i < nframes ; i ++)
                {
                    X_shift_vect.push_back (Xshift_arr[i]) ;
                    Y_shift_vect.push_back (Yshift_arr[i]) ;
                    theta_vect.push_back (theta_arrfinal[i]) ;
                    time_temp.push_back (time[i]) ;
                }
                LOG (INFO) << "\nSmoothening duration : " << delta_time << endl ;
               
                //method for applying smoothening over 'delta_time' duration
                status = local_Smoothening_new (time_temp , X_shift_vect ,orderpitch, delta_time ) ;
                if (status)
                {
                    LOG (ERROR) << "Error in local smoothening the X_shift" << endl ;
                    return (EXIT_FAILURE) ;
                }
                LOG(INFO)<<"outside";

                status = local_Smoothening_new (time_temp , Y_shift_vect , orderyaw , delta_time) ;
                if (status)
                {
                    LOG (ERROR) << "Error in local smoothening the y_shift" << endl ;
                    return (EXIT_FAILURE) ;
                }
 
                status = local_Smoothening_new (time_temp , theta_vect , orderroll , delta_time ) ;
                if (status)
                {
                    LOG (ERROR) << "Error in local smoothening the theta_shift" << endl ;
                    return (EXIT_FAILURE) ;
                }
 
                sprintf (temp_path_txt_temp , "%s/Delta_x_afterlocalSmoothingFilter.txt" , moduleoutdir) ;

                ofptr.open (temp_path_txt_temp , ios::out) ;

               
                delete[] Xshift_arr , Yshift_arr , theta_arrfinal ;
                Xshift_arr = new double[X_shift_vect.size ()] ;
                Yshift_arr = new double[X_shift_vect.size ()] ;
                theta_arrfinal = new double[X_shift_vect.size ()] ;
               // cout<<X_shift_vect.size ()<<endl;
                
                for (int i = 0 ; i < X_shift_vect.size () ; i ++)
                {
                    Xshift_arr[i] = X_shift_vect[i] ;
                    Yshift_arr[i] = Y_shift_vect[i] ;
                    theta_arrfinal[i] = theta_vect[i] ;
                    time[i] = time_temp[i] ;
                }
 
                for (int i = 0 ; i < X_shift_vect.size () ; i ++)
                {

                    ofptr << time[i] <<X_shift_vect[i]<<setw(20)<<Y_shift_vect[i]<<setw(20)<<theta_vect[i]<< setw (20) << Xshift_arr[i] << setw (20) << Yshift_arr[i] << setw(20)<<theta_arrfinal[i]<<endl ;

                }
  
                ofptr.close () ;
                nframes = X_shift_vect.size () ;                 
                xshiftvect.clear ();   yshiftvect.clear ();    theta_vect.clear ();   time_temp.clear ();
                
            }

          //writing output to final vectors
            for (int i = 0 ; i < nframes ; i ++)
            {
                xshiftvect.push_back (Xshift_arr[i]) ;
                yshiftvect.push_back (Yshift_arr[i]) ;
                thetavect.push_back (theta_arrfinal[i]) ;
            }

        }
        else
        {
            LOG (INFO) << "***Invalid value for type Filtering,it must be 0 or 1***" << endl ;
            return (EXIT_FAILURE) ;

        }
    }
        //frequency domain filtering
    else if (FreqDomainFilter_Flag == 1)
    {
        LOG (INFO) << "\nStarted Frequency domain Filtering..." << endl ;
        float *data , *freq , *delta_t , *fft ;
        long int *sample_no ;
        long int fftN = nframes_power2 * 2 ;
        int No_Of_OATRecords = nframes ;
        delta_t = new float [nframes_power2] ;
        data = new float[nframes_power2] ;
        sample_no = new long int[nframes_power2] ;
        freq = new float[nframes_power2] ;
        fft = new float [nframes_power2] ;

        for (long int i = 0 ; i < nframes_power2 ; i ++)
        {
            freq[i] = 0.0 ;
            sample_no[i] = 0 ;
            data[i] = 0.0 ;
            delta_t[i] = 0.0 ;
        }

        delta_t[0] = time[0] - time[0] ;
        sample_no[0] = 1 ;

        int k = 0 ;

        for (long int i = 1 ; i < nframes_power2 ; i ++)
        {
            if (i >= No_Of_OATRecords)
            {
                delta_t[i] = 0.0 ;
            }
            else
            {
                delta_t[i] = time[i] - time[k] ;
            }
            k ++ ;
            sample_no[i] = i + 1 ;
        }

        for (long int i = 0 ; i < nframes_power2 ; i ++)
        {
            if (i < No_Of_OATRecords - 1)
            {
                freq[i] = sample_no[i + 1] / (nframes_power2 * delta_t[i + 1]) ;
            }
            else
            {
                freq[i] = 0.0 ;
            }
        }

        double *delta_theta_complex = new double[fftN + 1] ;
        double *delta_x_complex = new double[fftN + 1] ;
        double *delta_y_complex = new double[fftN + 1] ;
        delta_theta_complex[0] = 0 ;
        delta_x_complex[0] = 0 ;
        delta_y_complex[0] = 0 ;

        for (long int i = 1 , j = 0 ; i <= fftN ; i ++)
        {
            if (i % 2 != 0)//in case of real part
            {
                delta_theta_complex[i] = theta_arrfinal[j] ;
                delta_x_complex[i] = Xshift_arr[j] ;
                delta_y_complex[i] = Yshift_arr[j] ;
                j ++ ;
            }
            else //incase of the imaginary part
            {

                delta_theta_complex[i] = 0.0 ;
                delta_x_complex[i] = 0.0 ;
                delta_y_complex[i] = 0.0 ;
            }
            //status=TransformNUVtoVIS (scaleNUV,scaleVIS,angle_NUVtoVIS,Shift_X,Shift_Y,delta_x_out,delta_y_out,nframes);
        }

        char temp_path_txt[FLEN_FILENAME] ;
        sprintf (temp_path_txt , "%s/Delta_x_beforeFFT.txt" , moduleoutdir) ; //file For storing the before FFT content.
        ofstream ofptr ;
        ofptr.open (temp_path_txt , ios::out) ;
        for (long int i = 0 ; i <= fftN ; i ++)
        {
            ofptr << delta_x_complex[i] << endl ;
        }
        ofptr.close () ;

        //Applying FFT
        doFft (delta_theta_complex , nframes_power2 , - 1) ;
        doFft (delta_x_complex , nframes_power2 , - 1) ;
        doFft (delta_y_complex , nframes_power2 , - 1) ;

        sprintf (temp_path_txt , "%s/freq.txt" , moduleoutdir) ;

        ofptr.open (temp_path_txt , ios::out) ;
        for (long int i = 0 ; i < nframes_power2 ; i ++)
        {
            ofptr << freq[i] << endl ;
        }
        ofptr.close () ;


        sprintf (temp_path_txt , "%s/Delta_x_afterFFT.txt" , moduleoutdir) ; //file For storing the After  FFT content.

        ofptr.open (temp_path_txt , ios::out) ;
        for (long int i = 0 ; i <= fftN ; i ++)
        {
            ofptr << delta_x_complex[i] << endl ;
        }
        ofptr.close () ;

        for (long int i = 1 , j = 0 ; i < fftN + 1 , j < nframes_power2 ; i = i + 2 , j ++)
        {

            if (freq[j] >= freqvalue || freq[j] == 0.0)
            {
                delta_theta_complex[i] = 0.0 ;
                delta_theta_complex[i + 1] = 0.0 ;
                delta_x_complex[i] = 0.0 ;
                delta_x_complex[i + 1] = 0.0 ;
                delta_y_complex[i] = 0.0 ;
                delta_y_complex[i + 1] = 0.0 ;
            }
        }

        sprintf (temp_path_txt , "%s/Delta_x_afterremovingHIGHcompo.txt" , moduleoutdir) ;

        ofptr.open (temp_path_txt , ios::out) ;
        for (long int i = 0 ; i <= fftN ; i ++)
        {
            ofptr << delta_x_complex[i] << endl ; //file For storing the  FFT contents after removing the higher frequency Component..
        }
        ofptr.close () ;

        //Applying inverse FFT
        doFft (delta_theta_complex , nframes_power2 , 1) ;
        doFft (delta_x_complex , nframes_power2 , 1) ;
        doFft (delta_y_complex , nframes_power2 , 1) ;

        sprintf (temp_path_txt , "%s/Delta_x_afterInverseFFT.txt" , moduleoutdir) ; //file For storing the after inverse content.

        ofptr.open (temp_path_txt , ios::out) ;
        for (long int i = 0 ; i <= fftN ; i ++)
        {
            ofptr << delta_x_complex[i] << endl ;
        }
        ofptr.close () ;

        double temp_roll = 0.0 , temp_pitch = 0.0 , temp_yaw = 0.0 ;

        if (nframes_power2 == 0)
        {
            LOG (ERROR) << "***Divide by zero*** *" << endl ;
            return (EXIT_FAILURE) ;
        }


        for (long int i = 1 , j = 0 ; i <= No_Of_OATRecords * 2 ; i ++)
        {

            if (i % 2 != 0)//taking only real part.
            {
                temp_roll = delta_theta_complex[i] / nframes_power2 ;

                temp_pitch = delta_x_complex[i] / nframes_power2 ;
                temp_yaw = delta_y_complex[i] / nframes_power2 ;
              
                //writing output to final vectors
                xshiftvect.push_back (temp_yaw) ;
                yshiftvect.push_back (temp_pitch) ;
                thetavect.push_back (temp_roll) ;
                j ++ ;
            }
        }
    }
        //no frequency domain, no time domain 
    else if (FreqDomainFilter_Flag != 0 && FreqDomainFilter_Flag != 1)
    {
        LOG (INFO) << "\n No Filtering and Fitting will done." << endl ;
        char temp_path[FLEN_FILENAME] ;
        sprintf (temp_path , "%s/observation.txt" , moduleoutdir) ; //file For storing the before FFT content.
        ofstream ofptr ;
        ofptr.open (temp_path , ios::out) ;
        
        for (long int i = 0 ; i < nframes ; i ++)
        {
            ofptr << time[i] << setw (20) << Xshift_arr[i] << setw (20) << Yshift_arr[i] << setw (20) << theta_arrfinal[i] << endl ;
        }
        ofptr.close () ;

         //writing output to final vectors
        for (int i = 0 ; i < nframes ; i ++)
        {
            xshiftvect.push_back (Xshift_arr[i]) ;
            yshiftvect.push_back (Yshift_arr[i]) ;
            thetavect.push_back (theta_arrfinal[i]) ;
        }
    }
    
    return (EXIT_SUCCESS) ;

}


/*int uvtDriftComputation::local_Smoothening (vector<double> &time_data , vector<double> &delta_term , int order , 
        double delta_time , vector<double> &delta_term_final , vector<double> &time_data_final)
{

    int last_index = 0 ;
    int flag = 0 ;
    double time_temp = time_data[0] ;
    double time_last_temp ;
    int time_last_index = - 1 ;
    
    vector<double> time_data_intermediate , delta_x_intermediate , delta_y_intermediate , delta_theta_intermediate ;
    double *arr_term , *arr_time ;

    delta_term_final.clear () ;

    if (delta_time > (time_data[time_data.size () - 1] - time_data[0]))
    {
        LOG (ERROR) << "***Delta Time is Greater then total time Data ***" << endl ;
        LOG (ERROR) << "***Fitting Will be done on total size of the time Data***" << endl ;
    }

    //loop for fitting  the delta_term  with delta_time duration
    while (time_temp < time_data[time_data.size () - 1])
    {
        time_temp = time_data[time_last_index + 1] ;//first time of fitting period

        delta_x_intermediate.clear () ;
        delta_y_intermediate.clear () ;
        delta_theta_intermediate.clear () ;
        time_data_intermediate.clear () ;
        time_last_temp = time_temp + delta_time ;//last time of  fitting period
      
        //if final time exeed to total time range 
        if (time_last_temp > time_data[time_data.size () - 1])
        {

            for (int i = last_index ; i <= time_data.size () - 1 ; i ++)
            {
                time_data_intermediate.push_back (time_data[i]) ;
                delta_x_intermediate.push_back (delta_term[i]) ;

            }

            goto label ;

        }
        
        //setting last time of fitting period to the  matched time within input time data &also store index of that last time of fitting
        for (int i = last_index ; i < time_data.size () ; i ++)
        {
            if (time_last_temp >= time_data[i] && time_last_temp < time_data[i + 1])
            {
                time_last_temp = time_data[i] ;
                time_last_index = i ;

            }
        }

label2:

         //loop for storing  delta term and time data of  current delta_term range
        for (int i = last_index ; i < time_data.size () ; i ++)
        {

            if (time_data[i] >= time_temp && time_data[i] <= time_last_temp)
            {
                time_data_intermediate.push_back (time_data[i]) ;
                delta_x_intermediate.push_back (delta_term[i]) ;

                last_index = i + 1 ;
            }

        }
        if (flag)
        {
            flag = 0 ;
            goto label3 ;

        }

label:
          //checking  availability of data for next  if only one data left after current delta time duration 
        if (((time_data.size () - 1)-(last_index + 1)) < 2)
        {

            time_last_temp = time_data[time_data.size () - 1] ;
            time_last_index = time_data.size () - 1 ;
            flag = 1 ;
            goto label2 ;

        }
label3:

        arr_term = new double[delta_x_intermediate.size ()] ;

        arr_time = new double[delta_x_intermediate.size ()] ;
        cout<<endl<<"New Set ";
        for (int i = 0 ; i < delta_x_intermediate.size () ; i ++)
        {
            arr_term[i] = 0.0 ;   arr_time[i] = 0.0 ;
            arr_term[i] = delta_x_intermediate[i] ;  arr_time[i] = time_data_intermediate[i] ;
            cout<<endl<<arr_time[i]<<" "<<arr_term[i];
        }
       
        arr_term = DataFitting (arr_time , arr_term , order , delta_x_intermediate.size ()) ;


        for (int i = 0 ; i < delta_x_intermediate.size () ; i ++)
        {
            delta_term_final.push_back (arr_term[i]) ;

            time_data_final.push_back (arr_time[i]) ;
        }
        time_temp = time_last_temp ;

    }
   
    delta_x_intermediate.clear () ;
    delta_y_intermediate.clear () ;
    delta_theta_intermediate.clear () ;
    time_data_intermediate.clear () ;

    return 0 ;


}*/


int uvtDriftComputation::findShiftsNtheta (int totalelements , vector<float>  &Xone , vector<float> &Yone , vector<float> int1, vector<float> &Xtwo , vector<float> &Ytwo ,vector<float> ints2,
                                                                        vector<float> &DiffOfX , vector<float> &DiffOfY ,bool flag_theta_computation, double &Xdx , double &Ydy , double &Theta,double &minEleX,double &maxeEleX,double &minEleY,double &maxEleY)
{
    if(totalelements==0 ){
        LOG(ERROR)<<"***Insufficient number of points found***";
        return(EXIT_SUCCESS);
    }
    if(totalelements<=2){
        flag_theta_computation=0;
    }
   // bool flag_theta_computation;
    if(flag_theta_computation==0)
    { 
        
        
        LOG(INFO)<<"NO theta value will be calculated"<<" "<<totalelements;
        vector<double> diff_x_cumm,diff_y_cumm;
          for (int t = 0 ; t < totalelements  ; t++)
          {
              LOG(INFO)<<" DX "<<Xtwo[t] - Xone[t]<<" DY "<<Ytwo[t] - Yone[t]<<" "<<Xone[t]<<" "<<Xtwo[t]<<" "<<Yone[t]<<" "<<Ytwo[t];
              diff_x_cumm.push_back (Xtwo[t] - Xone[t]);
              diff_y_cumm.push_back (Ytwo[t] - Yone[t]);
             
          }
       minEleX= *min_element(diff_x_cumm.begin(),diff_x_cumm.end());
       maxeEleX=*max_element(diff_x_cumm.begin(),diff_x_cumm.end());
       minEleY=*min_element(diff_y_cumm.begin(),diff_y_cumm.end());
       maxEleY=*max_element(diff_y_cumm.begin(),diff_y_cumm.end());
         //of1<<setw(20)<<min_element(diff_x_cumm.begin(),diff_x_cumm.end())<<setw(20)<<max_element(diff_x_cumm.begin(),diff_x_cumm.end())<<setw(20)<<
         //        min_element(diff_y_cumm.begin(),diff_y_cumm.end())<<setw(20)<<max_element(diff_y_cumm.begin(),diff_y_cumm.end())<<endl;
        Xdx=getmean (diff_x_cumm.data (),diff_x_cumm.size ());
        Ydy=getmean (diff_y_cumm.data (),diff_y_cumm.size ());
        Theta=0.0f;
        LOG(INFO)<<" X "<<Xdx<<" Y "<<Ydy;
        
    }
    
    
    else if(flag_theta_computation==1){
    if (option_LeastSquare == 1)
    {
        
        spMatrix B ((totalelements) * 2 , 1) ;
        spMatrix A ((totalelements) * 2 , 3) ;
        spMatrix X (3 , 1) ;
        int temp = 0 ;

        for (int t = 0 ; t < totalelements * 2 ; t = t + 2)
        {
            B (t , 0) = (Xtwo[temp] - Xone[temp]) ;                //+IMAGE_ARRAYSIZE*0.5;
            B (t + 1 , 0) = (Ytwo[temp]- Yone[temp]) ;            //+IMAGE_ARRAYSIZE*0.5;

            A (t , 0) = - 1.0 * (Yone[temp] - (xsize/div_fact) / 2) ;
            A (t , 1) = 1 ;
            A (t , 2) = 0 ;

            A (t + 1 , 0) = (Xone[temp] - (xsize/div_fact) / 2) ;
            A (t + 1 , 1) = 0 ;
            A (t + 1 , 2) = 1 ;

            temp ++ ;
        }
       
        X.ApplyLeastSquare (A , B) ;

        Xdx = X (1 , 0) ;
        Ydy = X (2 , 0) ;
        Theta = X (0 , 0)*180/M_PI ;
        //X (0 , 0) = X (0 , 0) ;
        // ofptr << p << setw (15) << X (0 , 0)*180 / M_PI << setw (25) << X (1 , 0) << setw (25) << X (2 , 0) << endl ;
    }
    else if (option_LeastSquare == 2)
    {
        spMatrix A (totalelements * 2 , 4) ;
        spMatrix B (totalelements * 2 , 1) ;
        spMatrix X (4 , 1) ;
        for (int aindex = 0 ; aindex < totalelements ; aindex ++)
        {
            A (2 * aindex , 0) = Xone[aindex] - (xsize/div_fact) * 0.5 ;
            A (2 * aindex , 1) = - 1.0 * (Yone[aindex] - (xsize/div_fact) * 0.5) ;
            A (2 * aindex , 2) = 1.0 ;
            A (2 * aindex , 3) = 0.0 ;
            A (2 * aindex + 1 , 0) = (Yone[aindex] - (xsize/div_fact) * 0.5) ;
            A (2 * aindex + 1 , 1) = Xone[aindex] - (xsize/div_fact) * 0.5 ;
            A (2 * aindex + 1 , 2) = 0.0 ;
            A (2 * aindex + 1 , 3) = 1.0 ;
            B (2 * aindex , 0) = Xtwo[aindex] - (xsize/div_fact)* 0.5 ;
            B (2 * aindex + 1 , 0) = Ytwo[aindex] - (xsize/div_fact) * 0.5 ;
        }
        X.ApplyLeastSquare (A , B) ;

        double theta = atan2 (X (1 , 0) , X (0 , 0)) ;

        // ofptr << theta * 180 / M_PI << setw (15) << X (0 , 0) << setw (25) << X (1 , 0) << endl ;
        Xdx = X (2 , 0) ;
        Ydy = X (3 , 0) ;
        Theta = theta*180/M_PI ;
    }
    else if (option_LeastSquare == 3)
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
        if(sum_int_ref==0){
            LOG(ERROR)<<"***Sum of intencity is zero!!!!***";
            return(EXIT_FAILURE);
        }
        for (int k = 0 ; k < totalelements ; k ++)
        {
            /* weight used for a star*/
            
          //  double w = 1.0 ;
            double w =int1[k]/sum_int_ref;
            
            //     w=pow( (25.0*photonbk_per_pixel+starmag[k]), 2.0 )/
            //        (  (25.0*photonbk_per_pixel)+(0.09*starmag[k])  );


            /* row 1 */
            double y_mod = Ytwo[k]-((xsize/div_fact) * 0.5) ;
            double x_mod = Xtwo[k]-((xsize/div_fact) * 0.5) ;
            
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
           // LOG(INFO)<<w<<" "<<sum_int_ref;
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

      //  LOG(INFO)<<a11<<" "<<a12<<" "<<a13<<" "<<a21<<" "<<a22<<" "<<a23<<" "<<a31<<" "<<a32<<" "<<a33<<" "<<b1<<" "<<b2<<" "<<b3; 
        X.ApplyLeastSquare (A , B) ;

        // ofptr << X (2 , 0)*180 / M_PI << setw (15) << X (0 , 0) << setw (25) << X (1 , 0) << endl ;
        Xdx = X (0 , 0) ;
        Ydy = X (1 , 0) ;
        Theta = X (2 , 0)*180/M_PI ;
    }
    }
    return (EXIT_SUCCESS) ;
}


int uvtDriftComputation::matchStars (int numrowsFirstfile , int numrowsSecfile , float divFact , float *xlocFirst , float *ylocFirst ,
        float *xlocSec , float *ylocSec , float * int1,float *int2,vector<float> &matchPixelXone , vector<float>  &matchPixelYone ,
        vector<float> &matchPixelXtwo , vector<float>  &matchPixelYtwo , vector<float> &int_one,vector<float> &int_two, 
        vector<float> &matchPixelDiffX , vector<float> &matchPixelDiffY)
{
     char temp_path_txt_temp[FLEN_FILENAME];
     int cnt = 0 ;
     vector<float> xtrack,ytrack;
     float temp_x1 = 0 , temp_x2 = 0 , temp_y1 = 0 , temp_y2 = 0 ,temp_int_one=0,temp_int_two=0;
//     cout<<"p"<<endl;
//     exit(1);
 //   cout << numrowsFirstfile << " " << numrowsSecfile << endl ;
    ofstream ofptr1 ;
    if(match_starsfile_gen)
    {
    sprintf (temp_path_txt_temp , "%s/matchstars_%d.txt" , moduleoutdir,cntglobal) ;
    ofptr1.open (temp_path_txt_temp , ios::out) ;
    
    }
   
   // cout<<numrowsFirstfile<<" "<<numrowsSecfile<<" "<<divFact<<endl;
    for (long int index = 0 ; index < numrowsFirstfile ; index ++) //loop for finding the similar x-coordinates and y-coordinates  from the both frames
    {        
        //cout<<index<<endl;
        if(xlocFirst[index]!=-9999 && ylocFirst[index]!=-9999)
        {
            temp_x1 = xlocFirst[index] / divFact ;
            temp_y1 = ylocFirst[index] / divFact ;   
          for (long int i = 0 ; i < numrowsSecfile ; i ++)
        {
        if(xlocSec[i]!=-9999  && ylocSec[i]!=-9999){
             float diff_x = 0.0 , diff_y = 0.0 ;
            temp_x2 = xlocSec[i] / divFact ;
            temp_y2 = ylocSec[i] / divFact ;
            
           vector<float>::iterator xdupli=find(xtrack.begin (),xtrack.end (),temp_x2);
           vector<float>::iterator ydupli =find (ytrack.begin (),ytrack.end (),temp_y2);
          
          
           if ( xdupli !=xtrack.end () && ydupli !=ytrack.end ())
           {           
               continue;
           }
            diff_x = temp_x2 - temp_x1 ;
            diff_y = temp_y2 - temp_y1 ;

            if ((diff_x * diff_x) < diff_dist && (diff_y * diff_y) < diff_dist ) // finding the similar points
            {
                 xtrack.push_back (temp_x2);
                 ytrack.push_back (temp_y2);
                matchPixelXone.push_back (temp_x1) ;
                matchPixelYone.push_back (temp_y1) ;
                matchPixelXtwo.push_back (temp_x2) ;
                matchPixelYtwo.push_back (temp_y2) ;
                int_one.push_back(int1[(int)index]);
                int_two.push_back(int2[(int)i]);
                matchPixelDiffX.push_back (diff_x) ;
                matchPixelDiffY.push_back (diff_y) ;
                if(match_starsfile_gen)
                {
                ofptr1<<temp_x1<<setw(20)<<temp_y1<<setw(20)<<setw(20)<<temp_x2<<setw(20)<<temp_y2<<setw(20)<<diff_x<<setw(20)<<diff_y<<endl;
                }
                cnt ++ ;
                break ;

            }
        }
          }
        }
        
//
//        for (long int i = 0 ; i < numrowsSecfile ; i ++)
//        {
//            float diff_x = 0.0 , diff_y = 0.0 ;
//            temp_x2 = xlocSec[i] / divFact ;
//            temp_y2 = ylocSec[i] / divFact ;
//            
//           vector<float>::iterator xdupli=find(xtrack.begin (),xtrack.end (),temp_x2);
//           vector<float>::iterator ydupli =find (ytrack.begin (),ytrack.end (),temp_y2);
//          
//          
//           if ( xdupli !=xtrack.end () && ydupli !=ytrack.end ())
//           {           
//               continue;
//           }
//            diff_x = temp_x2 - temp_x1 ;
//            diff_y = temp_y2 - temp_y1 ;
//
//            if ((diff_x * diff_x) < diff_dist && (diff_y * diff_y) < diff_dist ) // finding the similar points
//            {
//                 xtrack.push_back (temp_x2);
//                 ytrack.push_back (temp_y2);
//                matchPixelXone.push_back (temp_x1) ;
//                matchPixelYone.push_back (temp_y1) ;
//                matchPixelXtwo.push_back (temp_x2) ;
//                matchPixelYtwo.push_back (temp_y2) ;
//                matchPixelDiffX.push_back (diff_x) ;
//                matchPixelDiffY.push_back (diff_y) ;
//                if(match_starsfile_gen)
//                {
//                ofptr1<<temp_x1<<setw(20)<<temp_y1<<setw(20)<<setw(20)<<temp_x2<<setw(20)<<temp_y2<<setw(20)<<diff_x<<setw(20)<<diff_y<<endl;
//                }
//                cnt ++ ;
//                break ;
//
//            }
//        }
    }
     if(match_starsfile_gen)
     {
    ofptr1.close ();
     }
    return cnt ;
}


int uvtDriftComputation::writeDrift ()
{
    
    fitsfile *fdrift ;
    int status = 0 ;
    char caldb_orig_path[FLEN_FILENAME];
    sprintf (infile_drift , "%s/%s_dr.fits" , moduleoutdir , nameprefix) ;
  
    fits_create_file (&fdrift , infile_drift , &status) ;
    printError (status , "Error in creating the file") ;

    char *ttype[] = {"Time" , "X_shift" , "Y_shift" , "Rotation"} ;
     char *ttype1[] = {"Time" , "YAW" , "PITCH" , "ROLL"} ;
    char *tform1[] = {"1D" , "1E" , "1E" , "1E"} ;
    fits_create_tbl (fdrift , BINARY_TBL , 0 , 4 , ttype1 , tform1 , NULL , "BINTABLE" , &status) ;
    fits_create_tbl (fdrift , BINARY_TBL , 0 , 4 , ttype , tform1 , NULL , "BINTABLE" , &status) ;
    fits_write_key (fdrift , TINT , "ThetaFlag" , &flag_theta , "Theta value flag" , &status) ;
    printError (status , "ThetaFlag value not updated" , infile_drift) ; //for creating name for output information file
    fits_close_file (fdrift , &status) ;
    printError (status , "Error in closing the file" , infile_drift) ;
    
    //char finalpath_One[FLEN_FILENAME];
    //char finalpath_Two[FLEN_FILENAME];
    
     double  *pitch_data= new double[nframes];
     double  *yaw_data=new double[nframes];
     double  *roll_data = new double[nframes];
//    if (strcasecmp (datainfo.getDetector () , "VIS") == 0)
//    {
//        //char temp_filepath[FLEN_FILENAME];
//        //strcpy (caldb_orig_path,caldbDir);
//        //const  char *teldefcaldbfile=caldb_handler.getTelDefFile (caldbDir,datainfo.getDetector (),datainfo.getFilter ());
//        //sprintf(temp_filepath,"%s/%s/%s/%s/%s",caldbDir,CALDB_TELDEF_DIR,datainfo.getDetector (),datainfo.getFilter (),basename(teldefcaldbfile));
//        for(int i=0;i<nframes;i++)
//        {
//          //transformtoFUV (x_arr[i],y_arr[i],theta_arr[i],x_vis[i],y_vis[i],theta_vis[i]);
//        }
//        
//        
//                
//        
//       // RawToSat (x_arr.data (),y_arr.data (),nframes,temp_filepath);
//        //cout<<temp_filepath<<" "<<teldefcaldbfile<<endl;
//      
//       
////       TelDef teldef;
////       teldef.read (temp_filepath);
//      // readParamsFrmteldef ((char *)temp_filepath,&Shift_X,&Shift_Y);
//    
////       const char * tempPlatesclfile1=caldb_handler.getPlatScaleFile(caldbDir ,datainfo.getDetector ());
////       joinStrings (finalpath_One, 2 , caldbDir , tempPlatesclfile1) ;
////       readPlateScaleFile((char*)finalpath_One,&scaleXCH1,&scaleYCH1);
////              
////       char * tempPlatesclfile2=caldb_handler.getPlatScaleFile(caldb_orig_path ,"VIS");
////      
////       joinStrings (finalpath_Two, 2 , caldb_orig_path, tempPlatesclfile2) ;
////        
////       readPlateScaleFile((char*)finalpath_Two,&scaleXVIS,&scaleYVIS);
////       angle_NUVtoVIS=0;
////       status = transformCH1toCH2 (scaleXCH1 , scaleXVIS ,scaleYCH1,scaleYVIS, angle_NUVtoVIS , Shift_X , Shift_Y , x_arr.data () , y_arr.data () , nframes) ;
//    }
//    else{
//        
//        for(int i=0;i<nframes;i++)
//        {
//            x_vis[i]=x_arr[i];
//            y_vis[i]=y_arr[i];
//            theta_vis[i]=theta_arr[i];            
//        }
//    }

     //original
     
             if (strcasecmp (datainfo.getDetector () , "NUV") == 0)
        {
            for(int i=0;i<nframes;i++){
               
                transformDXDYDTHETAtoRPY_NUV (x_arr[i],y_arr[i],theta_arr[i],roll_data[i],pitch_data[i],yaw_data[i]);
                yaw_data[i]=yaw_data[i]/3600;
                pitch_data[i]=pitch_data[i]/3600;
            }
            
        }
        else if(strcasecmp (datainfo.getDetector () , "FUV") == 0)
        {
            for(int i=0;i<nframes;i++){
            transformDXDYDTHETAtoRPY_FUV (x_arr[i],y_arr[i],theta_arr[i],roll_data[i],pitch_data[i],yaw_data[i]);
            yaw_data[i]=yaw_data[i]/3600;
                pitch_data[i]=pitch_data[i]/3600;
            }
        }
        else if(strcasecmp (datainfo.getDetector () , "VIS") == 0)
        {
            for(int i=0;i<nframes;i++)
            {
                transformDXDYDTHETAtoRPY_VIS (x_arr[i],y_arr[i],theta_arr[i],roll_data[i],pitch_data[i],yaw_data[i]);
                yaw_data[i]=yaw_data[i]/3600;
                pitch_data[i]=pitch_data[i]/3600;
            }
             
        }
        else
        {
            LOG(INFO)<<"***Invalid Channel***"<<endl;
            return(EXIT_FAILURE);
        }
     /*Introduce matrix for transformation between star250 to spacecraft.
      * 
      */
     
     //Added
     for (int i=0;i<nframes;i++)
     {
      
         //time[i]=time[i]+datainfo.getIntegrationTime ()/2;
         //time[i]=time[i]-datainfo.getIntegrationTime ()/2;
         
         //time[i]=time[i]+datainfo.getIntegrationTime ()*(9.0/16.);
         time[i]=time[i]-datainfo.getIntegrationTime ()*(7.0/16.);
         // chosen 't - 0.5 sec' above ...
     }
     
     
     //till this

   status= writeColumnsToFITS (infile_drift,2,4,TDOUBLE,1,time,nframes,TDOUBLE,2,yaw_data,nframes,TDOUBLE,3,pitch_data,nframes,TDOUBLE,4,roll_data,nframes);
   if(status)
   {
       LOG(INFO)<<"Error in writing to the  Drift FITS  file"<<endl;
       return(EXIT_FAILURE);
   }
   status= writeColumnsToFITS (infile_drift,3,4,TDOUBLE,1,time,x_arr.size (),TDOUBLE,2,x_arr.data (),x_arr.size (),TDOUBLE,3,y_arr.data (),x_arr.size (),TDOUBLE,4,theta_arr.data (),x_arr.size ());
   if(status)
   {
       LOG(INFO)<<"Error in writing to the  Drift FITS  file"<<endl;
       return(EXIT_FAILURE);
   }

     return (EXIT_SUCCESS) ;
}


int uvtDriftComputation::local_Smoothening (vector<double> &time , vector<double> &data, int order,double delta_time)
{
    LOG(INFO)<<"Smoothening data at "<<delta_time<<" sec  intervals";
    vector<double> temp_data = data;
    data.clear();
    
    double starttime = time[0];
    double endtime = time[time.size()-1];
    
    if(endtime <= starttime){
        LOG(ERROR)<<"***Fitting not possible :StartTime : "<<starttime<<"  End time: "<<endtime<<"***";
        return (EXIT_FAILURE);
    }
    
    vector<double>  time_to_fit;       //vector to store time values for delta_time interval  
    vector<double> data_to_fit;       //vector to store data values for delta_time interval
    
    int j=0;                //looping variable
    int index=0;            //index for maintaining index of time in data
    
    cout.precision (10);
    cout.setf (ios::fixed,ios::floatfield);
    while(starttime<endtime)
    {
        
         for( ;time[index]<(starttime+delta_time);index++)
         {
           //  cout<<"Index::"<<index<<" "<<starttime <<" " <<starttime+delta_time<<endl;
            data_to_fit.push_back(temp_data[index]);
            time_to_fit.push_back(time[index]);
            if(index>=time.size()-1)  break;
         }
    //     cout<<"===================="<<endl;
        starttime = time[index];    //update start time after every delta_time 
                     
        double *arr = new double[time_to_fit.size()];
        LOG(INFO)<<"Time FIT SIZE "<<time_to_fit.size ();
        arr = DataFitting (time_to_fit.data() , data_to_fit.data() , order, time_to_fit.size());
        
        for(j=0;j<time_to_fit.size();j++)   
            data.push_back(arr[j]);                    //copy output of fitting to data vector which will hold all the final outputs
        
        time_to_fit.clear();
        data_to_fit.clear();
        delete[] arr;
    
        if(index>=time.size()-1) break;
    }
    temp_data.clear();    
    return (EXIT_SUCCESS);
}

int uvtDriftComputation::local_Smoothening_new (vector<double> &time , vector<double> &data, int order,double delta_time)
{
    LOG(INFO)<<"Smoothening data at "<<delta_time<<" sec  intervals";
    vector<double> temp_data = data;
    data.clear();
    
    double starttime = time[0];
    double endtime = time[time.size()-1];
    
    if(endtime <= starttime){
        LOG(ERROR)<<"***Fitting not possible :StartTime : "<<starttime<<"  End time: "<<endtime<<"***";
        return (EXIT_FAILURE);
    }
    double *data_arr = new double[time.size ()];
    double temp_val;
    for (int i=0;i<time.size ();i++)
    {
        data_arr[i]=temp_data[i];
    }
    vector<double>  time_to_fit;       //vector to store time values for delta_time interval  
    vector<double> data_to_fit;       //vector to store data values for delta_time interval
    
    int j=0;                //looping variable
    int index=0;            //index for maintaining index of time in data
    int next_start_index=0;;
    cout.precision (10);
    cout.setf (ios::fixed,ios::floatfield);
    double mid_point_val;
    double temp_time_store=time[1]-time[0];
    double TimeIntegration=temp_time_store;
    while(temp_time_store<delta_time )
        temp_time_store=temp_time_store+TimeIntegration;
    
    delta_time=temp_time_store;
    LOG(INFO)<<"Delta Time->"<<delta_time<<endl;
//    if(delta_time%datainfo.getIntegrationTime ()==0.0){
//        delta_time=delta_time+1.0;//incase of VIS channel.
//    }
    LOG(INFO)<<"Final Time Duration  used for smoothning is "<<delta_time;
    
            data_to_fit.clear ();time_to_fit.clear ();index=next_start_index;
    while(starttime<endtime-delta_time)
    {
      
         next_start_index=index;
         //LOG(INFO)<<setprecision (20)<<starttime<<" "<<setprecision (20)<<starttime+delta_time<<" "<<delta_time<<" "<<endl;
         for( ;time[index]<(starttime+delta_time);index++)
         {
           //  cout<<"Index::"<<index<<" "<<starttime <<" " <<starttime+delta_time<<endl;
            data_to_fit.push_back(temp_data[index]);
            time_to_fit.push_back(time[index]);
            if(index>=time.size()-1)  break;
         }
         if(time_to_fit.size ()%2 ==0 && index<temp_data.size()-1){
           data_to_fit.push_back(temp_data[index]);
           time_to_fit.push_back(time[index]);            
         }
        // LOG(INFO)<<data_to_fit.size();
    //     cout<<"===================="<<endl;
        index=next_start_index+1;
        if(index<temp_data.size()-1)
        starttime = time[index];    //update start time after every delta_time 
                 
      
        double *arr = new double[time_to_fit.size()];
      // LOG(INFO)<<"Time FIT SIZE "<<time_to_fit.size ();
        arr = DataFitting (time_to_fit.data() , data_to_fit.data() , order, time_to_fit.size());
        mid_point_val=arr[(time_to_fit.size ()/2)];//change it and Add it +1.
      //  temp_val=data_arr[next_start_index+time_to_fit.size ()/2];
       // LOG(INFO)<<temp_val<<" "<<mid_point_val;
        LOG(INFO)<<next_start_index+(time_to_fit.size ()/2)<<" "<<time_to_fit.size()<<" "<<index<<" "<<delta_time<<" "<<datainfo.getIntegrationTime();
//        exit(1);
        data_arr[next_start_index+(time_to_fit.size ()/2)]=mid_point_val;//add +1
       // data.
        //for(j=0;j<time_to_fit.size();j++)   data_arr[]
//            data.push_back(arr[j]);                    //copy output of fitting to data vector which will hold all the final outputs
        
        time_to_fit.clear();
        data_to_fit.clear();
        delete[] arr;
    
        if(index>=time.size()-1) break;
    }
    for(int i=0;i<time.size ();i++)
    {
        data.push_back (data_arr[i]);
    }
    temp_data.clear();    
    return (EXIT_SUCCESS);
}


 int  uvtDriftComputation:: removeRecords(vector<float>  &Xone , vector<float> &Yone ,vector<float> &Xtwo , vector<float> &Ytwo,vector<float> &DiffOfX , vector<float> &DiffOfY ,vector<float> &ints1,vector<float> &ints2,
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
                    if (Distance_vect[i]<=mean_dist+3*sd_dist  && Distance_vect[i]>(mean_dist-3*sd_dist))
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
     ofstream ofptr1;
     char temp_path_txt_temp[FLEN_FILENAME];
       if(match_starsfile_gen)
    {
    sprintf (temp_path_txt_temp , "%s/matchstars_%d.txt" , moduleoutdir,cntglobal) ;
    ofptr1.open (temp_path_txt_temp , ios::out) ;
    
    }
      if(match_starsfile_gen){
               
                //ofptr1<<temp_x1<<setw(20)<<temp_y1<<setw(20)<<setw(20)<<temp_x2<<setw(20)<<temp_y2<<setw(20)<<diff_x<<setw(20)<<diff_y<<endl;
               
     for (int i=0;i<new_two_ints.size ();i++){
        ofptr1<<newXone[i]<<setw(20)<<newYone[i]<<setw(20)<<setw(20)<<newXtwo[i]<<setw(20)<<newYtwo[i]<<setw(20)<<newDiffX[i]<<setw(20)<<newDiffY[i]<<endl;
     }
     
 }
     ofptr1.close ();
     return(EXIT_SUCCESS);
 }