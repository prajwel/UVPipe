/* 
 * File:   uvtQEMCPCorr.cpp
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
#include "stdio.h"
#include <fstream>
#include "uvtQEMCPCorr.h"
#include<pthread.h>
#include<DataInfo.h>
#include<caldb_Handler.h>
#include<glog/logging.h>
//#define IMG_DIM  512

uvtQEMCPCorr::uvtQEMCPCorr ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}
uvtQEMCPCorr::~uvtQEMCPCorr ()
{
    
}

int uvtQEMCPCorr::read (int argc , char** argv)
{
    int status = 0 ;

    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG(ERROR) << "***Error Initializing PIL***" ;
        return status ;
    }
    
    if (PIL_OK != (status = PILGetFname ("inputdatadir" , indir)))
    {
        LOG(ERROR) << endl << "***Error reading input data directory name***" ;
        return status ;
    }
   
    if (PIL_OK != (status = PILGetFname ("caldbDir" , caldbDir)))
    {
        LOG(ERROR) << endl << "***Error reading temperature vs filter file name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("lbtfile" , lbtfile)))
    {
        LOG(ERROR) << endl << "***Error reading lbt file name***" ;
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
        LOG(ERROR) << "***Error Reading history parameter***" << history << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("mode" , mode)))
    {
        LOG(ERROR) << "***Error Reading mode parameter:" << history << "***" ;
        return status ;
    }
    PILClose (status) ;
    return (EXIT_SUCCESS) ;
}

int uvtQEMCPCorr::read (char *indir , char *caldir , char *lbtfile , char*outdir ,
        int clobber , int history)
{
    strcpy (this->indir , indir) ;
    strcpy (this->caldbDir , caldir) ;
    strcpy (this->lbtfile , lbtfile) ;
    strcpy (this->outdir , outdir) ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}

void uvtQEMCPCorr::display ()
{
    LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT QEMCP CORRECTION PARAMETERS      " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input data directory                          :" << this->indir ;
     LOG(INFO) << endl << "caldb Dir  is                                       :" << this->caldbDir ;
    LOG(INFO) << endl << "LBT file                                             :" << this->lbtfile ;
    LOG(INFO) << endl << "Output Directory                               : " << this->outdir ;
    if (clobber == YES)
        LOG(INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG(INFO) << endl << "Overwrite                                         : NO" ;
    if (history == YES)
        LOG(INFO) << endl << "History                                             : YES" ;
    else
        LOG(INFO) << endl << "History                                              : NO" ;
    LOG(INFO) << endl << "------------------------------------------------------------------------" << endl ;
}

int uvtQEMCPCorr::uvtQEMCPCorrectionProcess ()
{

    /**Setting path for the output directory**/
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG(INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;
    string cmd ;
    if (DirExists (moduleoutdir) && clobber == YES)
    {
        LOG(INFO) << "Directory exists and clobber=yes" ;
        cmd = (string) "rm -rf " + (string) moduleoutdir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists (moduleoutdir) && clobber == NO)
    {
        LOG(INFO) << endl << moduleoutdir << "  already exists " ;
        LOG(INFO) << endl << "Use clobber=yes for overwriting" ;
        return (EXIT_FAILURE) ;
    }
/**Shell Command  for the making Directory**/
    cmd = "mkdir -p " + (string) moduleoutdir ;
    /**Executing the shell command**/
    system (cmd.c_str ()) ; // creating output directory to keep output from QEMCPCorrection
    LOG(INFO) << endl << moduleoutdir << "  directory created" ;
   // LOG(INFO) << indir << endl ;
    
    /**Searching the input information File**/
    string tempfilepath = searchFile (indir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG(INFO) << endl << "Error in finding info file" ;
        return (EXIT_FAILURE) ;
    }
       /**Setting input information File path**/
    sprintf (infofile_in , "%s/%s" , indir , tempfilepath.c_str()) ;
    LOG(INFO) << endl << " Input information file :" << infofile_in << endl ;
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    /**opening the input information File**/
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input information File" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU in input information file" , infofile_in) ;
  /**getting the keyword value of the input information File**/
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    LOG(INFO)<<"Filter : "<<datainfo.getFilter ()<<endl;
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    if (xsize <= 0 || ysize <= 0)
    {
        LOG(ERROR) << endl << "Invalid xsize/ysize\n" ;
        return (EXIT_FAILURE) ;
    }

    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the keyword value of the NAMEPRFX" , infofile_in) ; //for creating name for output information file

    //creating output information file
    sprintf (infofile_out , "%s/%s_qe.info" , moduleoutdir , nameprefix) ;
     LOG(INFO) << endl << " Output information file :" << infofile_out << endl ;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information file" , infofile_out) ;
    char *ttype[] = {"SignalFileList" , "ExposureFileList"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table in  the output information file" , infofile_out) ;
   /**writing the keyword value to the output information File**/
    datainfo.write (finfo_out) ; //writing basic data information

    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the keyvalue of the NAMEPRFX" , infofile_out) ; //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
     string tempname = caldb_handler.getQEFile (datainfo.getDetector () , datainfo.getObsMode () , caldbDir) ;
    if (tempname == "")
    {
        LOG(INFO) << endl << "Couldn't find QEMCP  file from caldb" << endl ;
        return (EXIT_FAILURE) ;
    }
    joinStrings (qeFile , 2 , caldbDir , tempname.c_str()) ;
    LOG(INFO) << endl << "QE file is  " << qeFile ;
    /*Reading the Caldb File for QE vs MG relationship*/
    status = readQEMCPFile () ;
    if (status)
    {
        LOG(ERROR) << "***Error reading the CALDB (QE vs MG ) file Data***" << endl ;
        return (EXIT_FAILURE) ;
    }
    status = getTemp () ;
    if (status)
    {
        LOG(ERROR) << "***temperature reading from the lbt file unsuccessful***" << endl ;
        return (EXIT_FAILURE) ;
    }   
    double t1 , t2 , x1 , x2 ;

    int filter_coln ;
    int filternumber ;
    sprintf(filter,"%s",datainfo.getFilter ());
  
    qe_mg_factor = new float[nCalDBTempValues] ;
    if (filter == (string) "F0")
    {
        filternumber = 0 ;
        filter_coln = 1 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = f0[q] ;
    }
    else  if (filter == (string) "F1")
    {
        filternumber = 1 ;
        filter_coln = 2 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = f1[q] ;
    }
    else if (filter == (string) "F2")
    {
        filternumber = 2 ;
        filter_coln = 3 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = f2[q] ;
    }
    else if (filter == (string) "F3")
    {
        filternumber = 3 ;
        filter_coln = 4 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = f3[q] ;
    }
    else if (filter == (string) "F4")
    {
        filternumber = 4 ;
        filter_coln = 5 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = f4[q] ;
    }
    else if (filter == (string) "F5")
    {
        filternumber = 5 ;
        filter_coln = 6 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = f5[q] ;
    }
    else if (filter == (string) "F6")
    {
        filternumber = 6 ;
        filter_coln = 7 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = f6[q] ;
    }
    else if (filter == (string) "F7")
    {
        filternumber = 7 ;
        filter_coln = 8 ;
        for (int q = 0 ; q < nCalDBTempValues ; q++)
            qe_mg_factor[q] = f7[q] ;
    }
    
    else{
        LOG(ERROR)<<endl<<"***Invalid filter option*** "<<endl;
        return (EXIT_FAILURE);
    }
    numfilter = 7 ;
    
//     if (filternumber >= 0 && filternumber < numfilter)
//        {
//            for (int j = 0 ; j < nCalDBTempValues - 1 ; j++)
//            {
//                if (temperature >= temp[j] && temperature < temp[j + 1])          //temperature - from LBT file ,temp- from caldb file
//                {
//                    t1 = temp[j] ;
//                    t2 = temp[j + 1] ;
//                    if ((t2 - t1) == 0)
//                    {
//                        LOG(INFO) << "***Divide By zero***" << endl ;
//                        return (EXIT_FAILURE) ;
//                    }
//                    x1 = qe_mg_factor[j] ;
//                    x2 = qe_mg_factor[j + 1] ;
//                    factor = x1 + ((temperature - t1)*((x2 - x1) / (t2 - t1))) ;
//                }
//            }
//        }
//        else
//        {
//            LOG(ERROR) << endl << "***Incorrect Filter number for Thermal Effects Correction on QE and MG*** " << endl ;
//        }

//        if (factor == -1)
//        {
//            LOG(ERROR) << endl << "***Temperature is out of range***" << endl ;
//            return (EXIT_FAILURE) ;
//        }
   
    

    if (datainfo.getModeFlag () == IM)            //for PC/IM   1 for IM, 0 for PC
    {
        fits_read_key (finfo_in , TSTRING , "SIGDIR" , signalframedir , NULL , &status) ;
        printError (status , "Error in reading the key value of the SIGDIR" , infofile_in) ;
        fits_read_key (finfo_in , TSTRING , "EXPDIR" , expframedir , NULL , &status) ;
        printError (status , "Error in reading the key value of the EXPDIR" , infofile_in) ;
        fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "Error in reading the key value of the NFILES" , infofile_in) ;
        // strcpy(expframedir,"ExposureFrames");        //setting exposure frame dir for output to be created
        LOG(INFO) << endl << "Number of files :" << nframes << endl ;
        //reading frame names from information file into vector
        sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        expframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) sigframelist , NULL , &status) ;
        printError (status , "Error in reading the column of the sigframelist " , infofile_in) ;
        fits_read_col (finfo_in , TSTRING , 2 , 1 , 1 , nframes , NULL , (void *) expframelist , NULL , &status) ;
        printError (status , "Error in reading the column of the expframelist" , infofile_in) ;
        
        if (qemcpCorrectionIM ()) return (EXIT_FAILURE) ;   //calling qemcp correction for IM

        fits_update_key (finfo_out , TSTRING , "SIGDIR" , signalframedir , NULL , &status) ;
        printError (status , "Error in updating the key value of the SIGDIR" , infofile_out) ;
        fits_update_key (finfo_out , TSTRING , "EXPDIR" , expframedir , NULL , &status) ;
        printError (status , "Error in updating the key value of the EXPDIR" , infofile_out) ;
        fits_update_key (finfo_out , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "Error in updating the key value of the NFILES" , infofile_out) ;

        freeMemory (sigframelist , nframes , NAMESIZE) ;
    }
    else if (datainfo.getModeFlag () == PC)
    {
        fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
        printError (status , "Error in reading the key value of the EVTFILE" , infofile_in) ;
//        fits_read_key (finfo_in , TSTRING , "IMGFILE" , imgfile , NULL , &status) ;
//        printError (status , "Error in reading the key value of the IMGFILE" , infofile_in) ;

        if (qemcpCorrectionPC ()) return (EXIT_FAILURE) ;
    }
    else
    {
        LOG(ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG(ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }

    fits_close_file (finfo_in , &status) ;
    printError (status , "Error in  closing the  input information file" , infofile_in) ;
    writeCommonKeywords (finfo_out , modulename) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the out information file" , infofile_out) ;

    LOG(INFO) << endl << "uvtQEMCPCorr  process completed successfully" << endl ;
    return (EXIT_SUCCESS) ;
}

int uvtQEMCPCorr::qemcpCorrectionIM ()
{

    LOG(INFO) << endl << "\nStarted QE MCP correction For IM mode" ;
    char **outsigframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    char **outexpframelist = allocateMemory<char>(nframes , NAMESIZE) ; //to store output exposure frame list  
    //creating signal directory
    char dir[FLEN_FILENAME] ;
    /**setting the signal frame directory path**/
    sprintf (dir , "%s/%s" , moduleoutdir , signalframedir) ;
    /**shell command for creating the Directory**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    //creating exposure directory
    /**setting the exposure frame directory path**/
    sprintf (dir , "%s/%s" , moduleoutdir , expframedir) ;
    cmd = "mkdir -p " + (string) dir ;
     /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    int status = 0 ;
    float framedata[xsize * ysize] ;
    double integrationtime = 0.0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char errstr[512] ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;

    vector<string> vhistorystr ;
    if(history==YES)  getHistory (vhistorystr) ;

    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;
  /*loop for processing of the  number of frames */ 
    LOG(INFO)<<"\nPerforming QEMCP correction..."<<endl;  
    double t1 , t2 , x1 , x2 ;
    for (int i = 0 ; i < nframes ; i++)
    {
        temperature=-9999;
        sprintf (errstr , "Error at iteration number %d" , i) ;
        fitsfile *fptr , *fout ;
        /**setting the input Signal Frame dir path**/
        sprintf (infile , "%s/%s/%s" , indir , signalframedir , sigframelist[i]) ;
        /**Initialization of the PIXEL array**/
        for (int i = 0 ; i < xsize * ysize ; i++)
            framedata[i] = 0.0f ;
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input signal file" , infile) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize * ysize , NULL , framedata , NULL , &status) ;
        printError (status , "Error in reading the pixels from input signal Frame " , infile) ;
        fits_read_key (fptr , TDOUBLE , "INT_TIME" , &integrationtime , NULL , &status) ;
        printError (status , "Error in reading the key value of the INT_TIME" , infile) ;
        fits_read_key (fptr , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ;
        printError (status , "Error in reading the key value of the FRAMENO" , infile) ;
        fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
        printError (status , "Error in  reading the key value of the FRMTIME" , infile) ;
       //factor = -1 ;
       
        for (int i=0;i<nrows_lbt;i++)
        {
          //  cout<<time_lbt[i]<<" "<<time_lbt[i+1]<<endl;
                    
            if (frametime>=time_lbt[i] && frametime<time_lbt[i+1])
            {
                
                temperature=(insideTemp[i]+outsideTemp[i])/2;
                break;
            }
        }
        if(temperature==-9999)
        {
            LOG(ERROR)<<"No record found in LBT file"<<endl;
            return(EXIT_FAILURE);
        }
                for (int j = 0 ; j < nCalDBTempValues - 1 ; j++)
            {
                if (temperature >= temp[j] && temperature < temp[j + 1])          //temperature - from LBT file ,temp- from caldb file
                {
                    t1 = temp[j] ;
                    t2 = temp[j + 1] ;
                    if ((t2 - t1) == 0)
                    {
                        LOG(INFO) << "***Divide By zero***" << endl ;
                        return (EXIT_FAILURE) ;
                    }
                    x1 = qe_mg_factor[j] ;
                    x2 = qe_mg_factor[j + 1] ;
                    factor = x1 + ((temperature - t1)*((x2 - x1) / (t2 - t1))) ;
                }
            }
    
        
        if(i==0){
            copyUsrkeywrdsTovect (fptr,key_record);
        }
        /*calculation of Factor for multiplication*/
       /*QEMCP correction is apply  here*/
        for (int pixno = 0 ; pixno < xsize * ysize ; pixno++){
          if(framedata[pixno]!=INVALID_PIX_VALUE)  
            framedata[pixno] = framedata[pixno] * factor ;
        }
        /*creating output Frame */
        /**setting the path for output Signal Frame**/
        sprintf (outfile , "%s/%s/%s_t%.4f_f%d_sig_qe.fits" , moduleoutdir , signalframedir , nameprefix , frametime , frameno) ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output Signal  file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the image " , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata , &status) ;
        printError (status , "Error in writing the pixels to output Signal  file ") ;
        copyUserKeywords (fptr , fout) ;
           if (history == YES) writeHistory (outfile , vhistorystr) ;
        writeCommonKeywords (fout , modulename) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output  Signal file" , outfile) ;
     
        strcpy (outsigframelist[i] , basename (outfile)) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input Signal file" , outfile) ;

        /*reading the Exposure input Frame */
        /**setting the input Exposure frame file**/
        sprintf (infile , "%s/%s/%s" , indir , expframedir , expframelist[i]) ;
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input Exposure file" , infile) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
        printError (status , "Error in Reading the pixel value of the input  Exposure file" , infile) ;
        
        /*writing  Exposure output Frame */
        
         /**setting the output Exposure frame file**/
        sprintf (outfile , "%s/%s/%s_t%f_f%d_exp_qe.fits" , moduleoutdir , expframedir , nameprefix , frametime , frameno) ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output Exposure file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the output File 1st HDU" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata , &status) ;
        printError (status , "Error in writing the pixels to output Exposure file " , outfile) ;
        copyUserKeywords (fptr , fout) ;
        
        if (history == YES) writeHistory (outfile , vhistorystr) ;
        strcpy (outexpframelist[i] , basename (outfile)) ;
        writeCommonKeywords (fout , modulename) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output Exposure file" , outfile) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input Exposure file" , outfile) ;
      //  cout<< "Total Files written=" << i << "     Remaining files= " << nframes - i <<" \r" ;
    }
   
    /*adding information to info (.info) file for next module Reference*/
    /**opening the output information File which has been generated earlier **/
    LOG(INFO)<<"\nWriting list of output file name to information file"<<endl;
    fitsfile *finfo_out ;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening  the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU" , infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , nframes , (void *) outsigframelist , &status) ;
    printError (status , "Error in  writing the output signal frame list" , infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 2 , 1 , 1 , nframes , (void *) outexpframelist , &status) ;
    printError (status , "Error in writing the output Exposure frame list" , infofile_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file" , infofile_out) ;
    writeUsrkeywordsFrmvect (infofile_out,key_record);
     if (history == YES) writeHistory (infofile_out , vhistorystr) ;
     vhistorystr.clear () ;
    /*releasing the memory*/
    freeMemory (outsigframelist , nframes , NAMESIZE) ;
    freeMemory (outexpframelist , nframes , NAMESIZE) ;

    return (EXIT_SUCCESS) ;
}

int uvtQEMCPCorr::qemcpCorrectionPC ()
{
    
     fitsfile *fptr , *fout , *finfo_out ;
 //   const char * tempfilepath = searchFile (indir , ".event") ;
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    char errstr[512] ;
    int status = 0 ;
    char eventfilename[FLEN_FILENAME] ;
    sprintf (eventfilename , "%s_qemcp.events" , nameprefix) ;
    sprintf (outfile , "%s/%s" , moduleoutdir , eventfilename) ;
    fits_create_file (&fout , outfile , &status) ;
    printError (status , errstr) ;

    long nrows = 0 ;
    /**setting  the input  event file path**/
    sprintf (infile , "%s/%s" , indir , eventfile) ;
     LOG(INFO) << "Input Event File  " << infile << endl ;
     LOG(INFO) << "Output  Event File  " << outfile << endl ;
     fits_open_file (&fptr , infile , READWRITE , &status) ;
     printError (status , "Error in opening the infile" , infile) ;
     copyUsrkeywrdsTovect (fptr,key_record);
     fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
     printError (status , "Error in coping the file" , outfile) ;
      int colnum=0;
     fits_movabs_hdu (fout , 2 , NULL , &status) ;
     printError (status , "***Moving To particular HDU fails***" , outfile) ;
     fits_get_num_rows (fout , &nrows , &status) ;
     printError (status , "Error reading the row number of the output file" , outfile) ;
     fits_get_colnum(fout,CASEINSEN,FF_COLNAME,&colnum,&status);
     printError (status , "Error reading the column  number of the output file" , outfile) ;
     
     float * ENP = new float[nrows];
      double* time_event= new double[nrows];
      fits_read_col (fout , TDOUBLE ,3 , 1 , 1 , nrows , NULL , time_event , NULL , &status) ;
     printError (status , "Error reading Column" , outfile) ;
     fits_read_col (fout , TFLOAT , colnum , 1 , 1 , nrows , NULL , ENP , NULL , &status) ;
     printError (status , "Error reading Column" , outfile) ;
     LOG(INFO)<<"\nPerforming QEMCP correction "<<endl;
    
     double t1 , t2 , x1 , x2 ;
     
     for (int p = 0 ; p < nrows ; p++)
        {
         temperature=-9999;
        
         for (int i=0;i<nrows_lbt;i++)
        {
          //  cout<<time_lbt[i]<<" "<<time_lbt[i+1]<<endl;
            
            if (time_event[p]>=time_lbt[i] && time_event[p]<time_lbt[i+1])
            {
                
                temperature=(insideTemp[i]+outsideTemp[i])/2;
                break;
            }
        }
        if(temperature==-9999)
        {
            LOG(ERROR)<<"No record found in LBT file"<<p<<" "<<time_event[p]<<endl;
            return(EXIT_FAILURE);
        }
                for (int j = 0 ; j < nCalDBTempValues - 1 ; j++)
            {
                if (temperature >= temp[j] && temperature < temp[j + 1])          //temperature - from LBT file ,temp- from caldb file
                {
                    t1 = temp[j] ;
                    t2 = temp[j + 1] ;
                    if ((t2 - t1) == 0)
                    {
                        LOG(INFO) << "***Divide By zero***" << endl ;
                        return (EXIT_FAILURE) ;
                    }
                    x1 = qe_mg_factor[j] ;
                    x2 = qe_mg_factor[j + 1] ;
                    factor = x1 + ((temperature - t1)*((x2 - x1) / (t2 - t1))) ;
                }
            }
    
              
         
            ENP[p] = ENP[p] * factor ;
        }
        fits_write_col (fout , TFLOAT , colnum , 1 , 1 , nrows , ENP , &status) ;
        printError (status , "Error in  writing  the column of multiple photon event" , outfile) ;
        fits_movabs_hdu (fout , 1 , NULL , &status) ;
        printError (status , "Error in moving to 2nd HDU" , infofile_out) ;
    /***Writes the  common information to  output information File***/
  
    vector<string> vhistorystr ;
    
    /*HISTORY Writing*/
    if (history == YES){
        getHistory (vhistorystr) ;
        writeHistory(outfile , vhistorystr) ;
    }
       fits_movabs_hdu (fout , 1 , NULL , &status) ;
     printError (status , "***Moving To particular HDU fails***" , outfile) ;
    writeCommonKeywords (fout , modulename) ;
    fits_close_file (fout , &status) ;
    printError (status , "Error in Closing the output Event File" , outfile) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the input Event File" , outfile) ;
   
    /**update information to the  output info file**/
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU in output information file" , infofile_out) ;
    fits_update_key (finfo_out , TSTRING , "EVTFILE" , basename (outfile) , NULL , &status) ;
    printError (status , "Error in updating the key value of the output information file") ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file  " , infofile_out) ;
      writeUsrkeywordsFrmvect (infofile_out,key_record);
     if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    return (EXIT_SUCCESS) ;
}

int uvtQEMCPCorr::readQEMCPFile ()
{
    LOG(INFO) << endl << "Reading QE MCP  Temperature vs filter file from calDB........" ;
    fitsfile *fqemcp ;
    int status = 0 ;
    fits_open_file (&fqemcp , qeFile , READONLY , &status) ;
    printError (status , "Error in opening the qeFile ",qeFile) ;
    int hdutype ;
    fits_movabs_hdu (fqemcp , 2 , &hdutype , &status) ;
    printError (status , "Error in moving to 2nd HDU in qeFile  ",qeFile) ;
    if (hdutype != BINARY_TBL)
    {
        LOG(ERROR) << endl << "***Expected binary table at hdu 2 of temperature vs filter file*** " << endl ;
        return (EXIT_FAILURE) ;
    }
    long nrows ;
    fits_get_num_rows (fqemcp , &nrows , &status) ;
    printError (status , "Error in readQEMCP()",qeFile) ;
    nCalDBTempValues = nrows ;
    temp = new float[nrows] ;
    f0 = new float[nrows] ;
    f1 = new float[nrows] ;
    f2 = new float[nrows] ;
    f3 = new float[nrows] ;
    f5 = new float[nrows] ;
    f6 = new float[nrows] ;
    f4 = new float[nrows] ;
    fits_read_col (fqemcp , TFLOAT , 1 , 1 , 1 , nrows , NULL , (void*) temp , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 2 , 1 , 1 , nrows , NULL , (void*) f0 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
   
    fits_read_col (fqemcp , TFLOAT , 3 , 1 , 1 , nrows , NULL , (void*) f1 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 4 , 1 , 1 , nrows , NULL , (void*) f2 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 5 , 1 , 1 , nrows , NULL , (void*) f3 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 6 , 1 , 1 , nrows , NULL , (void*) f4 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 7 , 1 , 1 , nrows , NULL , (void*) f5 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    fits_read_col (fqemcp , TFLOAT , 8 , 1 , 1 , nrows , NULL , (void*) f6 , NULL , &status) ;
    printError (status , "Reading a column Fails in caldb",qeFile) ;
    LOG(INFO) << endl << "\nReading QEMCP Temperature vs filter file from caldb Finished.." ;
    fits_close_file (fqemcp , &status) ;
    printError (status , "Error in closing qeFile",qeFile) ;
    return (EXIT_SUCCESS) ;
}

int uvtQEMCPCorr::getTemp ()
{
    int status = 0 ;
    //reading observation id from data
   
    fitsfile *flbt ;
    fits_open_file (&flbt , lbtfile , READONLY , &status) ;
    printError (status , "Error opening LBT file",lbtfile) ;
    
    char obsid[FLEN_KEYWORD];
    
    //reading observation id from header of LBT file
    fits_read_key(flbt, TSTRING, "OBS_ID",obsid, NULL, &status);
    printError(status, " Error reading OBS_ID from header " ,lbtfile);
    
    fits_movabs_hdu (flbt , 2 , NULL , &status) ;
    printError (status , "Error moving to HDU 2 of lbtfile",lbtfile) ;
    //check when  number of Rows in lbt file > 1(TBD)
    fits_get_num_rows (flbt , &nrows_lbt , &status) ;
    printError (status , "Error reading the number of rows in lbt file",lbtfile) ;
     int colinside = 0 , coloutside = 0 ;             //variables to store column numbers to be used from LBT file
    char *md = datainfo.getDetector();      //read channel for data
    
    LOG(INFO)<<endl<<"Channel :" <<md<<endl;
          
    if (md = (char *) "NUV")
    {
        colinside = INSIDE_TEMP_NUV ;
        coloutside = OUTSIDE_TEMP_NUV ;
    }
    else if (md = (char *) "FUV")
    {
        colinside =INSIDE_TEMP_FUV ;
        coloutside = OUTSIDE_TEMP_FUV ;
    }
    else if (md = (char *) "VIS")
    {
        colinside = INSIDE_TEMP_VIS;
        coloutside = OUTSIDE_TEMP_VIS;
    }
    time_lbt=new double[nrows_lbt];
   insideTemp = new float [nrows_lbt];
    outsideTemp= new float[nrows_lbt];
//    int rowno;              
    fits_read_col (flbt , TDOUBLE , 1 , 1 , 1 ,nrows_lbt , NULL , time_lbt , NULL , &status) ;
    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
    fits_read_col (flbt , TFLOAT , colinside , 1 , 1 ,nrows_lbt , NULL , insideTemp , NULL , &status) ;
    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
    fits_read_col (flbt , TFLOAT , coloutside , 1 , 1 ,nrows_lbt , NULL , outsideTemp , NULL , &status) ;
    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
                   //row number in LBT file from which the temperatures will be read
     fits_close_file (flbt , &status) ;
   printError (status , "Error in closing file",lbtfile) ;
   // obsid=69//this is hardcoded.
    //Extracting last 10 digits from observation id
  //  string obsid_str = obsid;
    //int pos_last_underscore = obsid_str.find_last_of("_");
    //string obsid_no_str = obsid_str.substr(pos_last_underscore+1,10);
    //int obsid_no = atol(obsid_no_str.c_str());
        
   // unsigned int  observation_id[total_rows];      
    /**Reading the inside TEMP**/
    //int colnum_obsid = 2;                                      //column number for observation id in LBT file
//    fits_read_col (flbt , TUINT , colnum_obsid , 1 , 1 ,total_rows , NULL ,observation_id  , NULL , &status) ;  
//    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
//    bool flag_lbt=FALSE;                            //flag to check whether required observation id is found in LBT or not
//    char * temp_obsid;
//     for (int i=0;i<total_rows;i++)
//    {
//         
//        if(observation_id[i]==obsid_no)
//        {
//            rowno=i;
//            flag_lbt=TRUE;
//            break;
//            
//        }
//    }
    
 //   if(flag_lbt==TRUE)
   // {
    
   //     fits_read_col (flbt , TFLOAT , colinside , rowno+1 , 1 ,1 , NULL , &insideTemp , NULL , &status) ;
    //    printError (status , "Error in reading the  column value of the Inside temp",lbtfile) ;
        /**Reading the outside TEMP**/
     //  fits_read_col (flbt , TFLOAT , coloutside , rowno+1 , 1 , 1 , NULL , &outsideTemp , NULL , &status) ;
      //  printError (status , "Error in reading the column value of the Outside temp",lbtfile) ;
       // temperature = (insideTemp + outsideTemp) / 2 ;
       
  //  }
//    else
//    {
//        LOG(INFO)<<"***Temperature not found***"<<endl;
//        return(EXIT_FAILURE);
//    }
        
  
    return (EXIT_SUCCESS) ;
}

int uvtQEMCPCorr::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    char  name_factor[FLEN_FILENAME];
    sprintf(name_factor,"%f",factor);
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) indir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Temperature vs Filter file=  " + (string)qeFile) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
     vhistory.push_back ((string)getSerialNo (cnt)+" lbtfile=" + (string) lbtfile) ;
     vhistory.push_back ((string)getSerialNo (cnt)+" Factor=" + (string)name_factor) ;
     vhistory.push_back ((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir) ;
     if (clobber == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+" clobber=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+" clobber=no") ;
    if (history == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+" history=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+" history=no") ;
    vhistory.push_back (" Parameter List END") ;
    return (EXIT_SUCCESS) ;
}
int uvtQEMCPCorr::readQEMCPFileforCmnArray (char* filename,char *lbtfilename,int &caldb_counts,double &temperature,char *detector, vector<float> &temp,vector<float> &f0,vector<float> &f1,vector<float> &f2,vector<float> &f3,vector<float> &f4,vector<float> &f5,vector<float> &f6)
{
    //char name;
    //name =filename;
    int status =0;
    fitsfile *fqemcp;
//    sprintf(name,"%s",filename);
    strcpy(this->qeFile,filename);
     strcpy(this->lbtfile,lbtfilename);
   this->datainfo.setDetector (detector);
    readQEMCPFile ();
    caldb_counts=this->nCalDBTempValues;
     fits_open_file (&fqemcp , qeFile , READONLY , &status) ;
    printError (status , "Error in readQEMCP() ",qeFile) ;
    int hdutype ;
    fits_movabs_hdu (fqemcp , 2 , &hdutype , &status) ;
    printError (status , "Error in readQEMCP() ",qeFile) ;
    if (hdutype != BINARY_TBL)
    {
        LOG(ERROR) << endl << "***Expected binary table at hdu 2 of temperature vs filter file*** " << endl ;
        return (EXIT_FAILURE) ;
    }
    long nrows ;
    fits_get_num_rows (fqemcp , &nrows , &status) ;
    printError (status , "Error in readQEMCP()",qeFile) ;
   
    fits_close_file(fqemcp,&status);
    printError (status , "Error in closing the qemcp file",qeFile) ;
   //f7=new float [nrows];
    for(int i=0;i<nrows;i++)
    {
        f0.push_back (this->f0[i]);f1.push_back (this->f1[i]);f2.push_back (this->f2[i]);f3.push_back (this->f3[i]);
        f4.push_back (this->f4[i]);f5.push_back (this->f5[i]);f6.push_back (this->f6[i]);temp.push_back (this->temp[i]);
    }
    getTemp ();
    temperature=this->temperature;
    //f0=this->f0; //f1=this->f1; f2=this->f2; f3=this->f3; f4=this->f4; f5=this->f5; f6=this->f6;temp=this->temp;
 //  f0[0]=this->f0[0];
//  cout<<"FF"<<endl;
//cout<<f0[1]<<endl;
  return(EXIT_SUCCESS);
    
}
