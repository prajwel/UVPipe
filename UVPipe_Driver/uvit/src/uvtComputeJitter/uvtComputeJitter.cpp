/* 
 * File:   uvtComputeJitter.cpp
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
#include <uvtComputeJitter.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<glog/logging.h>
#include "fft.h"
#include "spMatrix.h"
#include "transform.h"



//Constructor - called when object is created


uvtComputeJitter::uvtComputeJitter ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}

//Destructor


uvtComputeJitter::~ uvtComputeJitter ()
{
  //   delete[] time_data , time_gyro , time ;
}

//parameter File reading
int uvtComputeJitter::read (int argc , char** argv)
{
    int status = 0 ;
    
    status = readParams (argc , argv , 3 , FNAME , "inputdatadir" , inputdatadir ,
                                                FNAME, "caldbDir" , caldbDir ,
                                                INT , "freqDomainFilter_Flag" , &FreqDomainFilter_Flag 
                                          // REAL4 , "eulerangleone" , &eulerangle_one ,
                                          // REAL4 , "eulerangletwo" , &eulerangle_two ,
                                          // REAL4 , "euleranglethree" , &eulerangle_three ,
                                          //REAL4 , "scaleVIS" , &scaleChannel
                                        ) ;
    if (status) return (EXIT_FAILURE) ;
    
   
      
     if (FreqDomainFilter_Flag == 0)
    {
        status = readParams (argc , argv , 1 , INT , "type_Filtering" , &type_Filtering) ;
        if (status) return (EXIT_FAILURE) ;
       
        if (type_Filtering == 1)
        {
            status = readParams (argc , argv , 1 , INT , "fitting_flag" , &fittingflag) ;
            if (status) return (EXIT_FAILURE) ;
            if (fittingflag==1)
            {
                status = readParams (argc , argv , 3 , INT , "order_pitch" , &orderpitch ,
                        INT , "order_yaw" , &orderyaw ,
                        INT , "order_roll" , &orderroll 
                       // REAL , "delta_Time" , &delta_time

                        ) ;
                if (status) return (EXIT_FAILURE) ;
            }
        }
   }
      
    status = readParams (argc , argv , 1 , REAL , "freq_value" , &freqvalue) ;
            if (status) return (EXIT_FAILURE) ;
            
     status = readParams (argc , argv , 5 , FNAME , "outdir" , outdir ,
             FNAME,"gyro_file",&gyrofile,
            BOOL , "clobber" , &clobber ,
            BOOL , "history" , &history ,
            STRING , "mode" , mode) ;
    if (status) return (EXIT_FAILURE) ;
     

    return (EXIT_SUCCESS) ;
}


int uvtComputeJitter::read (char *inputdatadir , char *caldbDir,char *gyro_file , int freqDomainFilter_Flag , double freqvalue , int fitting_flag , int orderpitch , int orderyaw , int orderroll , char *outdir , int Type_Filtr,int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;
    strcpy (this->gyrofile , gyro_file) ;
    strcpy (this->caldbDir , caldbDir) ;
    this->FreqDomainFilter_Flag = freqDomainFilter_Flag ;
    this->type_Filtering=Type_Filtr;
    this->freqvalue = freqvalue ;
    if (fitting_flag == 1)
    {
        this->orderpitch = orderpitch ;
        this->orderyaw = orderyaw ;
        this->orderroll = orderroll ;
    }
     this->clobber = clobber ;
    this->fittingflag = fitting_flag ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}
//Parameter file content Display


void uvtComputeJitter::display ()
{
    LOG (INFO) << endl << "----------Display of parameters for Jitter  Computation  ---------" ;
    LOG (INFO) << endl << "Input Frame List Directory                        : " << inputdatadir ;
    LOG (INFO) << endl << "scale of channel                               : " << scaleChannel ;
    LOG (INFO) << endl << "1st of 3 Euler angle between NUV and VIS                               : " <<eulerangle_one ;
    LOG (INFO) << endl << "2nd of 3 Euler angle between NUV and VIS                               : " <<eulerangle_two ;
    LOG (INFO) << endl << "3rd of 3 Euler angle between NUV and VIS                               : " <<eulerangle_three ;
    LOG (INFO) << endl << "freqDomainFilter_flag                               : " << FreqDomainFilter_Flag ;
     if (FreqDomainFilter_Flag == 0)
    {
        LOG (INFO) << endl << "Type Filtering                               : " << type_Filtering ;
        if (type_Filtering == 1)
        {
            LOG (INFO) << endl << "Fitting Flag                               : " << fittingflag ;
            if (fittingflag)
            {
            LOG (INFO) << endl << "Order pitch                               : " << orderpitch ;
            LOG (INFO) << endl << "Order Yaw                               : " << orderyaw ;
            LOG (INFO) << endl << "Order Roll                               : " << orderroll ;
            
            }
        }
    }
    LOG (INFO) << endl << "Cut-off frequency                               : " << freqvalue ;
    LOG (INFO) << endl << "Output Directory                               : " << outdir ;
    if (clobber == YES)
        LOG (INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG (INFO) << endl << "Overwrite                                         : NO" ;
    if (history == YES)
        LOG (INFO) << endl << "History                                             : YES" ;
    else
        LOG (INFO) << endl << "History                                              : NO" ;

    LOG (INFO) << endl << "----------Display of parameters for Jitter Computation  Ends----------\n" ;

}

//Correction for the  Cosmic Ray process
int uvtComputeJitter::uvtComputeJitterProcess ()
{
    
    LOG (INFO) << endl << "Jitter Computation Started" << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG (INFO) << "Module Output Directory : " << moduleoutdir << endl ;
    
    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;

    LOG (INFO) << endl << moduleoutdir << "  directory created" ;

    //searching information file in input directory
    string  tempfilepath = searchFile (inputdatadir , ".info") ;
     if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }
    //setting path of input information file
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;

    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (! (FileExists (infofile_in)))
    {
        LOG (ERROR) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY" ;
        return (EXIT_FAILURE) ;
    }
    LOG (INFO) << endl << "\nInformation File :" << infofile_in ;

   //opening input information file
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the information file") ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in Moving the 2nd HDU") ;
    
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "NAMEPRFX keyword not Found") ; //for creating name for output information file

    fits_close_file(finfo_in,&status);
    //setting path of output information file
    sprintf (infofile_out , "%s/%s_jt.info" , moduleoutdir , nameprefix) ;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information file",infofile_out) ;
    char *ttype[] = {"DUMMY"} ;
    char *tform[] = {"A1"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table",infofile_out) ;
    
    datainfo.write (finfo_out) ; //writing basic data information

    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the key value of NAMEPRFX",infofile_out) ; //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file",infofile_out) ;
 
 if (!(datainfo.getModeFlag () == IM || datainfo.getModeFlag ()==PC))
 {
        LOG (ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG (ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
     
 }      
    //computing jitter  
    if (computeJitter ()) return (EXIT_FAILURE) ;
    
    return (EXIT_SUCCESS) ;

}


int uvtComputeJitter::computeJitter ()
{

    LOG (INFO) << endl << " Inside the Jitter Calculation" << endl ;
    int status = 0 ;

    status = readGyroFile () ; //reading the gyro file
    if (status)
    {
        LOG (ERROR) << "***Gyro file reading fails***" << endl ;
        return (EXIT_FAILURE) ;
    }

    //searching   drift output file from input drift directory
    string  tempfilepath = searchFile (inputdatadir , ".fits") ;
    sprintf (input_driftFile , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    LOG (INFO) << " Drift input  File is " << input_driftFile << endl ;

    G1_P_final.clear () ;
    G1_R_final.clear () ;
    G2_Y_final.clear () ;
    
    fitsfile *fptr ;
    fits_open_file (&fptr , input_driftFile , READONLY , &status) ;
    printError (status , "***input File reading Fails***") ;
    long numrows = 0 ;
    copyUsrkeywrdsTovect (fptr , key_records) ; 
    fits_movabs_hdu (fptr , 2 , NULL , &status) ;
    printError (status , "***Moving to  2nd HDU Fails***") ;
    fits_get_num_rows (fptr , &numrows , &status) ;
    printError (status , "Error in getting the number of rows ",input_driftFile) ;
    
    //check for no data
    if(numrows==0)
    {
    LOG(INFO)<<"Input given Drift file is empty"<<endl;
    return(EXIT_FAILURE);
    }
    
    time_data = new double[numrows] ;
    fits_read_col (fptr , TDOUBLE , 1 , 1 , 1 , numrows , NULL , (void *) time_data , NULL , &status) ;
    printError (status , "time data reading fails",input_driftFile) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the file",input_driftFile) ;
   
    status = extractGyro (numrows) ;//method for gyro extraction as per time in drift file 
    if (status)
    {
        LOG (ERROR) << "Extraction of gyro failed" << endl ;
        return (EXIT_FAILURE) ;
    }
 
    status = doFiltering () ; //Filtering 
    if (status)
    {
        LOG (ERROR) << "Filtering Failed" << endl ;
        return (EXIT_FAILURE) ;
    }
    
    status = doIntegration () ; //integration of the roll pitch and yaw of the gyro file.
    if (status)
    {
        LOG (ERROR) << "Integration Failed" << endl ;
        return (EXIT_FAILURE) ;
    }
   
    LOG (ERROR) << "Integration Process Completed" << endl ;
   //convert it into PIXELS.
    
    
//    for (long int i = 0 ; i < integrated_p.size () ; i ++)
//    {
//        integrated_r[i] = integrated_r[i] ;
//        integrated_p[i] = integrated_p[i] * IMAGE_ARRAYSIZE / 0.5 ;
//        integrated_y[i] = integrated_y[i] * IMAGE_ARRAYSIZE / 0.5 ;
//    }
//       
   //roll pitch yaw to dx dy and theta
    
    
    char caldb_orig_path[FLEN_FILENAME];
     char finalpath_One[FLEN_FILENAME];
     float  pix_persec_X,pix_persec_Y;
  
    //spacecraft to VIS transform
//     if (strcasecmp (datainfo.getDetector () , "VIS") != 0)//if VIS channel than no need to transform
//    {
//      strcpy (caldb_orig_path,caldbDir);
////      const  char *teldefcaldbfile=caldb_handler.getTelDefFile (caldbDir);
////     
////       readKeywords ((char*)teldefcaldbfile, 1, 3 , TFLOAT , "EULERA" , &eulerangle_one ,
////                TFLOAT , "EULERB" , &eulerangle_two ,
////                TFLOAT , "EULERC" , &eulerangle_three) ;    
//      
//     
//      const char * tempPlatesclfile1=caldb_handler.getPlatScaleFile(caldbDir ,datainfo.getDetector ());
//      joinStrings (finalpath_One, 2 , caldbDir , tempPlatesclfile1) ;
//      readPlateScaleFile((char*)finalpath_One,&pix_persec_X,&pix_persec_Y);
//      
//      pix_perdeg_X=pix_persec_X*60*60; //converting to degree
//      pix_perdeg_Y=pix_persec_Y*60*60;
//      
//      status = RPY2dxdy (pix_perdeg_X , pix_perdeg_Y,integrated_y.data () , integrated_p.data () , integrated_p.size ()) ;
//     if (status)
//    {
//        cout << "Error in Roll,pitch yaw to shifts conversion" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//      
//     
//      
//      /*transform dx,dy and dtheta from spacecraft  frame to respective channel.
//      this can be done in two steps.
//       a) spacecraft to UVIT frame
//      b) UVIT frame to respective channel.
//      */
//      
////    status = transformSpacecraft2VIS (pix_perdeg_X , pix_perdeg_Y,scaleChannel , eulerangle_one , eulerangle_two , eulerangle_three , integrated_y.data () , integrated_p.data () , integrated_p.size ()) ;
////      if (status)
////    {
////        cout << "Error in spacecraft to VIS fails" << endl ;
////        return (EXIT_FAILURE) ;
////    }
//    
//     }
        
   writeJitter (); //write x,y,theta  to output jitter file.
   
    return (EXIT_SUCCESS) ;
}


int uvtComputeJitter::getHistory (vector<string> &vhistory)
{
    int cnt=0;
   vhistory.clear () ;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;
    //    vhistory.push_back ("P2 caldbDir=" + (string) caldbDir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir) ;
    char temp[FLEN_FILENAME] ;
    sprintf (temp , "%d" , FreqDomainFilter_Flag) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Frequency Domain Filter Flag=" + (string) temp) ;

    sprintf (temp , "%f" , freqvalue) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Frequency Value=" + (string) temp) ;

    if (FreqDomainFilter_Flag == 0)
    {
        sprintf (temp , "%d" , type_Filtering) ;
        vhistory.push_back ((string)getSerialNo (cnt)+" Type Filtering used=" + (string) temp) ;
        if (type_Filtering)
        {
            sprintf (temp , "%d" , fittingflag) ;
            vhistory.push_back ((string)getSerialNo (cnt)+" Fitting Flag =" + (string) temp) ;
            if (fittingflag)
            {
                sprintf (temp , "%d" , orderpitch) ;
                vhistory.push_back ((string)getSerialNo (cnt)+" order pitch =" + (string) temp) ;
                sprintf (temp , "%d" , orderroll) ;
                vhistory.push_back ((string)getSerialNo (cnt)+" order roll =" + (string) temp) ;
                sprintf (temp , "%d" , orderyaw) ;
                vhistory.push_back ((string)getSerialNo (cnt)+" order yaw =" + (string) temp) ;

            }
        }
    }

    // vhistory.push_back ("Centroid Directory=" + (string) centroidDir) ;
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

int uvtComputeJitter::readGyroFile ()
{
    LOG (INFO) << endl << "Reading Gyro......." << endl ;
    fitsfile *fgyro ;
    int status = 0 ;

    long nrows ;
  
    //opening input GYRO file
   
    fits_open_file (&fgyro , gyrofile , READONLY , &status) ;
    printError (status , "***Error in opening   gyro  File***") ;
    if (status)
    {
        LOG (ERROR) << "Error in opening GYRO file" ;
        return (EXIT_FAILURE) ;
    }
    fits_movabs_hdu (fgyro , 2 , NULL , &status) ;
    printError (status , "***Error Moving to 2nd HDU***") ;
    if (status)
    {
        LOG (ERROR) << "Error in moving to particular HDU  in GYRO file" ;
        return (EXIT_FAILURE) ;
    }
    fits_get_num_rows (fgyro , &nrows , &status) ;
    printError (status , "***Error Reading the number of Rows***") ;
    if (status)
    {
        LOG (ERROR) << "Error in getting number of rows in GYRO file" ;
        return (EXIT_FAILURE) ;
    }

    //allocating memory to store GYRO data
    nrows_gyro = nrows ;
    time = new double[nrows] ;
    G1_R = new double[nrows] ;
    G1_P = new double[nrows] ;
    G2_Y = new double[nrows] ;
    G2_P = new double[nrows] ;
    G3_Y = new double[nrows] ;
    G3_R = new double[nrows] ;
    gyro_rows = nrows ;

    LOG (INFO) << "Total  number of rows In GYRO file : " << nrows << endl ;
    
    fits_close_file (fgyro , &status) ;
    printError (status , "Error in closing the file",gyrofile) ;
    if (status)
    {
        LOG (ERROR) << "Error in closing GYRO file" ;
        return (EXIT_FAILURE) ;
    }    
    status=readColumnsFromFITS (gyrofile,2,4,TDOUBLE,1,time,nrows,TDOUBLE,5,G2_Y,nrows,
            TDOUBLE,6,G2_P,nrows,TDOUBLE,7,G3_R,nrows);
    
    if(status)
    {
        LOG(INFO)<<"Error reading column from gyro file"<<endl;
        return(EXIT_FAILURE);
    }
   
    LOG (INFO) << endl << "Reading GYRO file  completed" << endl ;
    return (EXIT_SUCCESS) ;

}


int uvtComputeJitter::doIntegration ()
{
    double integrated_roll = 0.0 , integrated_yaw = 0.0 , integrated_pitch = 0.0 ;
    LOG (INFO) << "****Integration Process Started****" << endl ;
   
    integrated_p.clear () ;
    integrated_r.clear () ;
    integrated_y.clear () ;
 
    //read first R,P,Y from attitude file by converting QUNS to R,P,Y.
    integrated_p.push_back (0.0) ;
    integrated_r.push_back (0.0) ;
    integrated_y.push_back (0.0) ;
    
    integrated_timeGyro.push_back (time_final[0]) ;
 
    
         
    for (int i = 1 ; i < delta_x_final.size () ; i ++)
    {        
        double time_diff = time_final[i] - time_final[i - 1] ;
        integrated_yaw = integrated_y[i - 1] + (delta_x_final[i - 1] * time_diff) ;
        integrated_pitch = integrated_p[i - 1] + (delta_y_final[i - 1] * time_diff) ;
        integrated_roll = integrated_r[i - 1] + (delta_theta_final[i - 1] * time_diff) ;
        integrated_timeGyro.push_back (time_final[i]) ;
        integrated_p.push_back (integrated_pitch) ;
        integrated_r.push_back (integrated_roll) ;
        integrated_y.push_back (integrated_yaw) ;
        //cout<<i<<" "<<integrated_timeGyro[i-1]<<" "<<integrated_p[i-1]<<" "<<integrated_y[i-1]<<" "<<integrated_r[i-1]<<endl;
        
    }
     
    char temp_path_txt[FLEN_FILENAME] ;
    sprintf (temp_path_txt , "%s/Delta_x_Integrated.txt" , moduleoutdir) ; //file For storing the before FFT content.
    ofstream ofptr ;
    ofptr.precision (8) ;
    ofptr.setf (ios::fixed , ios::floatfield) ;
    ofptr.open (temp_path_txt , ios::out) ;
    for (long int i = 0 ; i < integrated_p.size () ; i ++)
    {
        ofptr << integrated_timeGyro[i] << "   " << integrated_y[i] << "   " << integrated_p[i] << " " << integrated_r[i] << endl ;
    }
    ofptr.close () ;

    return (EXIT_SUCCESS) ;
}

int uvtComputeJitter::doFiltering ()
{

    LOG (INFO) << "Filtering Started" << endl ;
    delta_x_final.clear () ;
    delta_y_final.clear () ;
    delta_theta_final.clear () ;
    time_final.clear () ;
  
    //Data Array size computation
    int power = 1 ;
    while (power < G1_P_final.size ())
    {
        power = power * 2 ;
    }
    gyroData_power2 = power ;
    double *delta_x = new double[gyroData_power2] ;

    double *delta_y = new double[gyroData_power2] ;
    double *delta_theta = new double[gyroData_power2] ;
    double *time_Integrated = new double[gyroData_power2] ;
  
    for (int i = 0 ; i < gyroData_power2 ; i ++)
    {
        delta_x[i] = 0.0 ;
        delta_y[i] = 0.0 ;
        delta_theta[i] = 0.0 ;
        time_Integrated[i] = 0.0 ;
    }

    for (int c = 0 ; c < G1_P_final.size () ; c ++)
    {
        time_Integrated[c] = time_final_gyro[c] ;
        delta_x[c] = G1_P_final[c] ;
        delta_y[c] = G2_Y_final[c] ;
        delta_theta[c] = G1_R_final[c] ; 
    }

    /*----Should we multiply it with FOCUL LENGH?----*/
    //time domain filtering 
    if (FreqDomainFilter_Flag == 0)
    {
        if (type_Filtering == 0)
        {
            double cutof_f = freqvalue ;
            double RC = 1 / (2 * M_PI * cutof_f) ;
            double delta_t = time_final_gyro[1] - time_final_gyro[0] ;
            double alpha = RC / (RC + delta_t) ;
            double *delta_x_out = new double[G1_P_final.size ()] ;
            double *delta_y_out = new double[G1_P_final.size ()] ;
            double *delta_theta_out = new double[G1_P_final.size ()] ;
            delta_x_out[0] = delta_x[0] ;
            delta_y_out[0] = delta_y[0] ;
            delta_theta_out[0] = delta_theta[0] ;
            char temp_path_txt_temp[FLEN_FILENAME] ;
            sprintf (temp_path_txt_temp , "%s/Delta_x_beforeHighpassFilter.txt" , moduleoutdir) ;
            ofstream ofptr ;
            ofptr.open (temp_path_txt_temp , ios::out) ;

            for (int i = 0 ; i < G1_P_final.size () ; i ++)
            {
                ofptr << time_final_gyro[i] << setw (20) << delta_theta[i] << setw (20) << delta_x[i] << setw (20) << delta_y[i] << endl ;
            }
            ofptr.close () ;

            sprintf (temp_path_txt_temp , "%s/Delta_x_afterHighpassFilter.txt" , moduleoutdir) ; //file For storing the before FFT content.

            ofptr.open (temp_path_txt_temp , ios::out) ;

            ofptr << time[0] << setw (20) << delta_theta_out[0] << setw (20) << delta_x_out[0] << setw (20) << delta_y_out[0] << endl ;

            //algorithm for spatial  domain filtering 
            for (int i = 1 ; i < G1_P_final.size () ; i ++)
            {
                delta_x_out[i] = alpha * delta_x_out[i - 1] + alpha * (delta_x[i] - delta_x[i - 1]) ;
                delta_y_out[i] = alpha * delta_y_out[i - 1] + alpha * (delta_y[i] - delta_y[i - 1]) ;
                delta_theta_out[i] = alpha * delta_theta_out[i - 1] + alpha * (delta_theta[i] - delta_theta[i - 1]) ;
                ofptr << time_final_gyro[i] << setw (20) << delta_theta_out[i] << setw (20) << delta_x_out[i] << setw (20) << delta_y_out[i] << endl ;

            }
            ofptr.close () ;
           
            for (int i = 0 ; i < G1_P_final.size () ; i ++)
            {
                delta_x_final.push_back (delta_x_out[i]) ;
                delta_y_final.push_back (delta_y_out[i]) ;
                delta_theta_final.push_back (delta_theta_out[i]) ;
                time_final.push_back (time_Integrated[i]) ;

            }

        }
        else if (type_Filtering == 1)
        {

            if (fittingflag == 1)
            {
                delta_x = DataFitting (time_Integrated , delta_x , orderpitch , G1_P_final.size ()) ;
                delta_y = DataFitting (time_Integrated , delta_y , orderyaw , G1_P_final.size ()) ;
                delta_theta = DataFitting (time_Integrated , delta_theta , orderroll , G1_P_final.size ()) ;
            }
            
            for (int i = 0 ; i < G1_P_final.size () ; i ++)
            {
                delta_x_final.push_back (delta_x[i]) ;
                delta_y_final.push_back (delta_y[i]) ;
                delta_theta_final.push_back (delta_theta[i]) ;
                time_final.push_back (time_Integrated[i]) ;
                //cout<<delta_x_final[i]<<" "<<delta_y_final[i]<<" "<<delta_theta_final[i]<<endl;
            }
        }
        else
        {
            LOG (INFO) << "Invalid value for type Filtering,it must be 0 or 1" << endl ;
            return (EXIT_FAILURE) ;

        }

    }
    //frequency domain filtering 
    else if (FreqDomainFilter_Flag == 1)
    {

        float *data , *freq , *delta_t , *fft ;
        long int *sample_no ;
        long int fftN = gyroData_power2 * 2 ;
        long int No_Of_OATRecords = G1_P_final.size () ;
        delta_t = new float [gyroData_power2] ;
        data = new float[gyroData_power2] ;
        sample_no = new long int[gyroData_power2] ;
        freq = new float[gyroData_power2] ;
        fft = new float [gyroData_power2] ;

        for (long int i = 0 ; i < gyroData_power2 ; i ++)
        {
            freq[i] = 0.0 ;
            sample_no[i] = 0 ;
            data[i] = 0.0 ;
            delta_t[i] = 0.0 ;
        }

        delta_t[0] = time_Integrated[0] - time_Integrated[0] ;
        sample_no[0] = 1 ;

        int k = 0 ;

        for (long int i = 1 ; i < gyroData_power2 ; i ++)
        {
            if (i >= No_Of_OATRecords)
            {
                delta_t[i] = 0.0 ;
            }
            else
            {
                delta_t[i] = time_Integrated[i] - time_Integrated[k] ;
            }
            k ++ ;
            sample_no[i] = i + 1 ;
        }

        for (long int i = 0 ; i < gyroData_power2 ; i ++)
        {
            if (i < No_Of_OATRecords - 1)
            {
                freq[i] = sample_no[i + 1] / (gyroData_power2 * delta_t[i + 1]) ;
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

        LOG (INFO) << "the delta_x value " << delta_x[3] << endl ;
        LOG (INFO) << "the value of FFTN " << fftN << endl ;
      
        for (long int i = 1 , j = 0 ; i <= fftN ; i ++)
        {
            if (i % 2 != 0)//in case of real part
            {
               
                delta_theta_complex[i] = delta_theta[j] ;
                delta_x_complex[i] = delta_x[j] ;
                delta_y_complex[i] = delta_y[j] ;
                j ++ ;
            }
            else//incase of the imaginary part
            {
                delta_theta_complex[i] = 0.0 ;
                delta_x_complex[i] = 0.0 ;
                delta_y_complex[i] = 0.0 ;
            }
        }

        //        LOG(INFO) << "the delta_x value " << delta_x_complex[7] << endl ;
        char temp_path_txt[FLEN_FILENAME] ;
        sprintf (temp_path_txt , "%s/Delta_x_beforeFFT.txt" , moduleoutdir) ; //file For storing the before FFT content.
        ofstream ofptr ;
        ofptr.precision (8) ;
        ofptr.setf (ios::fixed , ios::floatfield) ;
        ofptr.open (temp_path_txt , ios::out) ;
        for (long int i = 1 , j = 0 ; i <= fftN ; i = i + 2 , j ++)
        {
            ofptr << time_Integrated[j] << "   " << delta_x_complex[i] << "   " << delta_x_complex[i + 1] << endl ;
        }
        ofptr.close () ;

        //Applying FFT
        doFft  (delta_theta_complex , gyroData_power2 , - 1) ;
        doFft (delta_x_complex , gyroData_power2 , - 1) ;
        doFft (delta_y_complex , gyroData_power2 , - 1) ;

        sprintf (temp_path_txt , "%s/freq.txt" , moduleoutdir) ;

        ofptr.open (temp_path_txt , ios::out) ;
        for (long int i = 0 ; i < gyroData_power2 ; i ++)
        {
            ofptr << freq[i] << endl ;
        }
        ofptr.close () ;

        sprintf (temp_path_txt , "%s/Delta_x_afterFFT.txt" , moduleoutdir) ; //file For storing the After  FFT content.

        ofptr.open (temp_path_txt , ios::out) ;
        for (long int i = 1 , j = 0 ; i <= fftN ; i = i + 2 , j ++)
        {
            ofptr << time_Integrated[j] << "   " << delta_x_complex[i] << "   " << delta_x_complex[i + 1] << endl ;
        }
        ofptr.close () ;

        //Suppressing the low frequency components
        for (long int i = 1 , j = 0 ; i < fftN + 1 , j < gyroData_power2 ; i = i + 2 , j ++)
        {
            if (freq[j] <= freqvalue || freq[j] == 0)
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
        doFft (delta_theta_complex , gyroData_power2 , + 1) ;
        doFft (delta_x_complex , gyroData_power2 , + 1) ;
        doFft (delta_y_complex , gyroData_power2 , + 1) ;

        double temp_roll = 0.0 , temp_pitch = 0.0 , temp_yaw = 0.0 ;
        double temp_time = 0.0 ;

        if (gyroData_power2 == 0)
        {
            LOG (INFO) << "***Divide by zero*** *" << endl ;
            return (EXIT_FAILURE) ;
        }

        for (long int i = 1 , j = 0 ; i < No_Of_OATRecords * 2 ; i ++)
        {

            if (i % 2 != 0)//taking only real part.
            {
                temp_roll = delta_theta_complex[i] / gyroData_power2 ;

                temp_pitch = delta_x_complex[i] / gyroData_power2 ;
                temp_yaw = delta_y_complex[i] / gyroData_power2 ;
                temp_time = time_Integrated[j] ;
                delta_x_final.push_back (temp_pitch) ;
                delta_y_final.push_back (temp_yaw) ;
                delta_theta_final.push_back (temp_roll) ;
                time_final.push_back (temp_time) ;

                j ++ ;
            }
        }

        sprintf (temp_path_txt , "%s/Delta_x_afterInverse_temp.txt" , moduleoutdir) ; //file For storing the before FFT content.

        ofptr.open (temp_path_txt , ios::out) ;
        for (long int i = 0 ; i < delta_x_final.size () ; i = i + 1)
        {
            ofptr << time_final[i] << "   " << delta_x_final[i] << endl ;
        }
        ofptr.close () ;

    }
    //if nono of above
    else if (FreqDomainFilter_Flag != 0 && FreqDomainFilter_Flag != 1)
    {
       // cout<<G1_P_final.size ()<<endl;exit(1);
        for (int i = 0 ; i < G1_P_final.size () ; i ++)
        {
            delta_x_final.push_back (delta_x[i]) ;
            delta_y_final.push_back (delta_y[i]) ;
            delta_theta_final.push_back (delta_theta[i]) ;
            time_final.push_back (time_Integrated[i]) ;

        }

    }


    return (EXIT_SUCCESS) ;
}


int uvtComputeJitter::extractGyro (long nrows)
{
    int flag = 0 ;
    
    //loop for  finding range of data from gyro according to dataset time range
    for (int i = 0 ; i < nrows_gyro ; i ++)
    {
       
        if ((time_data[0] > time[i] && time_data[0] < time[i + 1]) || ((flag == 1) && time[i] <= time_data[nrows - 1]) || (time_data[0] == time[i]))
        {
            
            G1_P_final.push_back (G2_P[i]) ;
            G2_Y_final.push_back (G3_Y[i]) ;
            G1_R_final.push_back (G3_R[i]) ;
            time_final_gyro.push_back (time[i]) ;
            flag = 1 ;

        }
    }
//    if(G1_P_final.size ()<nrows-1){
//        LOG(ERROR)<<"Not enough record found in GYRO file wrt to Time";
//        return(EXIT_FAILURE);
//    }
//   
    return (EXIT_SUCCESS) ;
}


int uvtComputeJitter::writeJitter ()
{
    
    char infile_jitter[FLEN_FILENAME];
     vector<string> vhistorystr ;
     if(history==YES)    getHistory (vhistorystr) ;
    //setting  the path  for jitter file
    sprintf (infile_jitter , "%s/%s_jt.fits" , moduleoutdir , nameprefix) ;
    int status=0;
    
    //creating the output jitter file
    fitsfile *f_jitter ;
    fits_create_file (&f_jitter , infile_jitter , &status) ;
    printError (status , "Error in creating the jitter file",infile_jitter) ;
    char *ttype[] = {"Time" , "YAW" , "PITCH" , "ROLL"} ;
    char *tform[] = {"1D" , "1E" , "1E" , "1E"} ;

    fits_create_tbl (f_jitter , BINARY_TBL , 0 , 4 , ttype , tform , NULL , "BINTABLE" , &status) ;
    printError (status , "Error in creating table",infile_jitter) ;
        
    status =writeColumnsToFITS (infile_jitter,2,4, TDOUBLE,1,integrated_timeGyro.data (),(int)integrated_r.size (),
                                        TDOUBLE,2,integrated_p.data (),(int)integrated_r.size (),
                                        TDOUBLE,3,integrated_y.data () ,(int)integrated_r.size (),
                                        TDOUBLE,4,integrated_r.data (),(int)integrated_r.size ());
    if(status){
        LOG(INFO)<<"Error in writing to the  jitter  FITS  file"<<endl;
       return(EXIT_FAILURE);
    }
    writeUsrkeywordsFrmvect (infile_jitter , key_records) ;
    
    if (history==YES)  writeHistory (infile_jitter,vhistorystr);
    writeCommonKeywords (f_jitter,modulename);
 
    fits_close_file (f_jitter , &status) ;
    printError (status , "Error in closing the jitter file",infile_jitter) ;
   
    return(EXIT_SUCCESS);
}
