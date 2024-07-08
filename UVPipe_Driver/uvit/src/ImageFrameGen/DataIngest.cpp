/* 
 * File:   DataIngest.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#include<cstdlib>
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string.h>
#include<ctime>
#include<vector>
#include<map>
//#include<libgen.h>
//#include <fitsio.h>
#include "pil.h"
#include "pil_error.h"
#include "spMatrix.h"
#include "lib_crc.c"
#include "DataIngest.h"
#include "uvtUtils.h"
#include "Directory.h"
#include<glog/logging.h>
#include<algorithm>
#include<macro_def.h>
//#define MAX_STRING_SIZE	2048            // For CRC_Check Only
//#define P_CCITT     0x1021
//#define FRAME_NO_OF_EVENTS 336         //number of events in a frame of PC mode
//#define IMG_DIM 512                                     //Image dimensions
//#define START_WORD 30
//#define BYTES 2269

#define MODULENAME  "DataIngest"


using namespace std ;
bool compare_endtime (struct darkframe_info info1 , struct darkframe_info info2) ;
bool compare_starttime (struct darkframe_info info1 , struct darkframe_info info2) ;


bool compare_bitPosition (GTIParams P1 , GTIParams P2)
{
    return (P1.bitpos < P2.bitpos) ;
}


bool validate (char flag)
{
    if (flag == 'y' || flag == 'n') return true ;
}


typedef struct darkframe_info
{

    string darkfrm_name ;
    double darkframe_starttime , darkframe_endtime ;
} ;


DataIngest::DataIngest ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    sprintf (dark_outputDir , "%s" , "Dark") ;
}


DataIngest::~ DataIngest ()
{
    frm_count.clear () ;
    frm_time_int.clear () ;
} //destructor


/***Function for reading the parameters   when DataIngest  module is executed***/
int DataIngest::read (int argc , char** argv)
{

    int status ;

    status = readParams (argc , argv , 7 ,
            FNAME , "DataIngest_inFile" , inFile ,
            FNAME , "DataIngest_inFile_GTI" , inFile_GTI ,
            FNAME , "DarkFrameDir" , darkFramesPath ,
            FNAME , "DataIngest_inFile_LBT" , inFile_LBT ,
            FNAME , "DataIngest_inFile_TCT" , inFile_TCT ,
            FNAME , "caldbdir" , caldbindir ,
            FNAME , "DataIngest_output_dirpath" , outdir) ;

    if (status) return (EXIT_FAILURE) ;

    /*check whether the Data is of PC mode or IM mode,In case of PC mode take Extra parameter of parity check */

    char *obs_mode = new char[FLEN_VALUE] ;
    getKeywordVal (inFile , "OBS_MODE" , 1 , obs_mode) ;
    if (strcasecmp (obs_mode , "PC") == 0)
    {
        status = readParams (argc , argv , 1 , INT , "parityFlag" , &parity_flag) ;


    }
    status = readParams (argc , argv , 1 ,
            BOOL , "GTI_FLAG" , &gti_flag) ;
    if (status) return (EXIT_FAILURE) ;

    if (gti_flag == 1)
    {

        status = readParams (argc , argv , 2 , STRING , "All_or_custom" , param_all_or_cust ,
                //  FNAME ,"gtiparam_file",inGtifile,
                INT , "validGTI_flag" , &valid_gtiflag
                ) ;
        if (status) return (EXIT_FAILURE) ;

        if (strcasecmp (param_all_or_cust , "ALL") == 0)
        {
            all_Or_custom = 1 ;
        }
        else if ((strcasecmp (param_all_or_cust , "CUSTOM")) == 0)
        {
            all_Or_custom = 0 ;
        }

    }
    status = readParams (argc , argv , 4 ,
            BOOL , "dropFrame" , &dropFrame ,
            BOOL , "clobber" , &clobber ,
            BOOL , "history" , &history ,
            STRING , "mode" , &mode
            ) ;
    if (status) return (EXIT_FAILURE) ;

    /*dropFrame-if set to 1 then frame will be discarded  if CRC fails for particular frame  */

    return (EXIT_SUCCESS) ;
}


/***Function for reading the parameters   when DataIngest  module is executed in chain***/
int DataIngest::read (char* scidatafile , char *caldb_dir , char* tctfile , char* gtifile , char* lbtfile , char *darkDir , int gtiflag , int validbit , int allOrcust , char* outdir , int dropframe , int parityflag , int clobber , int history)
{
    strcpy (this->inFile , scidatafile) ;
    strcpy (this->inFile_GTI , gtifile) ;
    strcpy (this->inFile_LBT , lbtfile) ;
    strcpy (this->inFile_TCT , tctfile) ;
    strcpy (this->outdir , outdir) ;
    strcpy (this->darkFramesPath , darkDir) ;
    strcpy (this->caldbindir , caldb_dir) ;
    this->gti_flag = gtiflag ;
    this->valid_gtiflag = validbit ;
    this->all_Or_custom = allOrcust ;
    this->parity_flag = parityflag ;
    this->dropFrame = dropframe ;
    this->clobber = clobber ;
    this->history = history ;
}


/***Function to display DataIngest parameter ***/
void DataIngest::display ()
{
    LOG (INFO) << "-----------------------------------------------------------" ;
    LOG (INFO) << "             DATA INGEST PARAMETERS          " ;
    LOG (INFO) << "-----------------------------------------------------------" ;
    LOG (INFO) << " Input File : " << inFile ;
    LOG (INFO) << " Input GTI file : " << inFile_GTI ;
    LOG (INFO) << " Input TCT file : " << inFile_TCT ;
    LOG (INFO) << " Input LBT file : " << inFile_LBT ;
    LOG (INFO) << " Data Ingest Output Path : " << outdir ;
    LOG (INFO) << "GTI Flag: " << gti_flag ;
    if (gti_flag)
    {
        LOG (INFO) << "All or custom: " << all_Or_custom ;
        LOG (INFO) << "valid GTI flag  " << valid_gtiflag ;
    }

    if (dropFrame)
        LOG (INFO) << " Full Frame will be dropped for failed CRC" ;
    else
        LOG (INFO) << " Packet will be dropped for failed CRC" ;
    if (clobber)
        LOG (INFO) << " Clobber : File OverWrite " ;
    else
        LOG (INFO) << " Clobber : No Overwrite " ;
    if (history == YES)
        LOG (INFO) << " History : YES " ;
    else
        LOG (INFO) << " History : NO " ;
    LOG (INFO) << "----------------------------------------------------------" << endl ;

}

//readLevel1Darkfiles()
//findMedianofDarkAtBegandEnd();
//createDarkDirectoryInDataIngest
//float findMedianValue(unsigned short *input, int noOfValues);

//void medianDarksIM(char *filenameprefix,char *dirpath,int m,int n);
//void getMedianFrame(unsigned short *outputImage,unsigned short *allImagePixels,int N,int imageHeight,int imageWidth);


int DataIngest::DataIngestProcess ()
{
    LOG (INFO) << "Data Ingest Process started" << endl ;

    int status = 0 ;

    status = setFilePaths () ; //function to set all output file paths and create output directory

    if (status)
    {
        LOG (ERROR) << status << "***Error in setting filepaths***" << endl ;
        return (EXIT_FAILURE) ;
    }

    fitsfile *time_correction_file_ptr ; //for TCT level-1 file
    fitsfile *di_input_file_ptr , *inputfile_primary_ptr ; //for input level-1 file

    //inputfile_primary_ptr - First HDU in Science data file 
    //di_input_file_ptr - Third HDU in Science data file 

    long felement = 1 , nelement = 1 ;

    if (history == YES) getHistory (vhistorystr) ; //function to create history for output event file

    /*-----------------------------------Opening  and Reading Level-1 Input Science Data File---------------------------------*/
    fits_open_file (&inputfile_primary_ptr , inFile , READONLY , &status) ;
    printError (status , "Error in Opening the input file" , inFile) ; //required for keyword copying

    copyUsrkeywrdsTovect (inputfile_primary_ptr , key_record) ;

    datainfo.getInfo (inputfile_primary_ptr) ; //gets basic information about data like mode, source, filter, detector

    fits_open_file (&di_input_file_ptr , inFile , READONLY , &status) ;
    printError (status , "Error in opening the data ingest input file" , inFile) ;

    LOG (INFO) << "Observation Mode of data  is    " << datainfo.getObsMode () ;

    /***Move to HDU containing science data **/
    fits_movnam_hdu (di_input_file_ptr , BINARY_TBL , SCIENCEDATA_HDUNAME , 0 , &status) ;
    printError (status , "Error in moving to 2nd HDU in input File" , inFile) ;

    fits_get_num_rows (di_input_file_ptr , &nrow_l1 , &status) ; //Get number of rows in level 1 science data file in 'nrow_l1'
    printError (status , "Error in getting the number of rows of input file" , inFile) ;

    LOG (INFO) << "Number of rows in input file : " << nrow_l1 ;

    //check whether GTI filtering to be done or not
    if (gti_flag)
    {
        char * fgtiname ;
        fgtiname = getenv ("PFILES") ;
        sprintf (inGtifile , "%s/%s" , fgtiname , GTI_FILE_NAME) ;
        LOG (INFO) << "GTI  parameter File : " << inGtifile << endl ;
        if (! (FileExists (inGtifile)))
        {
            LOG (INFO) << endl << "GTI parameter file  not found at PFILES area." ;
            return (EXIT_FAILURE) ;
        }

    }

    /*-----------------------------------Creating Data Ingest Output fits file----------------------------------*/
    int frow = 1 ;
    int start_rowNo=1;
    /**creating data ingest output File**/
    fitsfile *di_output_file_ptr ;

    fits_create_file (&di_output_file_ptr , outFile , &status) ;
    printError (status , "***Error in creating the output file***" , outFile) ;

    fits_create_img (di_output_file_ptr , BYTE_IMG , 0 , NULL , &status) ;
    printError (status , "***Error in creating the output file***" , outFile) ;
    copyUserKeywords (inputfile_primary_ptr , di_output_file_ptr) ; //Copy keywords from input level-1 file to data ingest output file
    //writes date, origin, checksum and creator to file
    LOG (INFO) << outFile << "    file created" ;

    if (datainfo.getModeFlag () == PC)
    { //for Photon Counting Mode
        char *ttype1[] = {"PacketSequence " , "FrameCount" , "Time" , "2016 Centroids"} ;
        char *tform1[] = {"1U" , "1U" , "1D" , "2016B"} ;
        char *tunit1[] = {"" , "" , "" , ""} ;
        fits_create_tbl (di_output_file_ptr , BINARY_TBL , 0 , 4 , ttype1 , tform1 , tunit1 , SCIENCEDATA_HDUNAME , &status) ;
        printError (status , "***Error in creating the table For PC mode(dataIngestOut)***" , outFile) ;
    }
    else
    { //for Integrating Mode
        char *ttype1[] = {"PacketSequence " , "FrameCount" , "Time" , "1008 Pixels"} ;
        char *tform1[] = {"1U" , "1U" , "1D" , "1008U"} ;
        char *tunit1[] = {"" , "" , "" , ""} ;
        fits_create_tbl (di_output_file_ptr , BINARY_TBL , 0 , 4 , ttype1 , tform1 , tunit1 , SCIENCEDATA_HDUNAME , &status) ;
        printError (status , "***Error in creating the table for the IM mode(dataIngestOut)" , outFile) ;
    }
    fits_movnam_hdu (di_output_file_ptr , ANY_HDU , SCIENCEDATA_HDUNAME , 0 , &status) ;
    printError (status , "***Error in moving to 2nd HDU inn  output file(.dataIngestOut) ***" , outFile) ;

    /*-----------------------------------Creating Time Correction Fits File---------------------------------*/
    fits_create_file (&time_correction_file_ptr , timeCorrectionFile , &status) ;
    fits_create_img (time_correction_file_ptr , BYTE_IMG , 0 , NULL , &status) ;
    printError (status , "Error in creating the Image in the time_correction File" , timeCorrectionFile) ;

    copyUserKeywords (inputfile_primary_ptr , time_correction_file_ptr) ;

    char *ttype2[] = {"FrameCount" , "Time" , "Corrected time"} ;
    char *tform2[] = {"1U" , "1D" , "1D"} ;
    char *tunit2[] = {"" , "" , ""} ;
    fits_create_tbl (time_correction_file_ptr , BINARY_TBL , 0 , 3 , ttype2 , tform2 , tunit2 , "CORRECTED_TIME" , &status) ;
    printError (status , "Error in creating the table im timeCorrecton File" , timeCorrectionFile) ;
    fits_movnam_hdu (time_correction_file_ptr , ANY_HDU , "CORRECTED_TIME" , 0 , &status) ;
    printError (status , "Error in creating the Image in the time_correction File" , timeCorrectionFile) ;

    LOG (INFO) << timeCorrectionFile << "  file created" ;

    /*-----------------declaring arrays to read from input science data file------------------------------*/
    unsigned short *pktSeqCtrl = new unsigned short[BUFFERSIZE] ;
    checkMemoryAvailability (pktSeqCtrl , "pktSeqCtrl") ;

    //Packet Sequence control used for CRC and also to be copied to output file
    unsigned short *frameCount = new unsigned short[BUFFERSIZE] ;
    checkMemoryAvailability (frameCount , "frameCount") ; //for CRC and output file

    //for CRC and output file
    unsigned int *frameTime = new unsigned int[BUFFERSIZE] ;
    checkMemoryAvailability (frameTime , "frameTime") ;
     double *frameTime_temp = new double[BUFFERSIZE] ;
    checkMemoryAvailability (frameTime , "frameTime_temp") ;

    unsigned char *centroids = new unsigned char[BUFFERSIZE * CENTROID_BYTES] ;
    checkMemoryAvailability (centroids , "centroids") ;

    double *correctedTime = new double[BUFFERSIZE] ;
    checkMemoryAvailability (correctedTime , "correctedTime") ;

    unsigned short *pixel = new unsigned short[BUFFERSIZE * PIXELNO] ;
    checkMemoryAvailability (pixel , "pixel") ;

    unsigned char *GTI_arr = new unsigned char[BUFFERSIZE * GTI_BYTES] ;
    checkMemoryAvailability (GTI_arr , "GTI_arr") ;


    unsigned short *pktSeqCtrl_temp = new unsigned short[nrow_l1] ;
    checkMemoryAvailability (pktSeqCtrl , "pktSeqCtrl") ;
    /*--------------------Matrix declaration for time correction computations-------------------------*/

    /*Loop for-
     1. Reading input file by BUFFERSIZE at a time
     2. Populating map for time and frame
     3. Writing data to output file, time column is written after fitting
     4. Finding missing frames with frame counts
     */
    vector<unsigned short> missingFrameCounts ;
    vector<double> missingFrameTimes ;

    LOG (INFO)<< "Reading data from science data file" ;
    bool before_start_flag = FALSE ;
    int count_missingFrame_beforeStart = 0 ;
    bool first_Continues_frm = FALSE ; //checking and discarding first n continues frames;
    int segflag = 0 ;
    fits_read_col (di_input_file_ptr , TUSHORT , 8 , frow , felement , nrow_l1 , NULL , pktSeqCtrl_temp , NULL , &status) ;
    printError (status , "Error in reading the column " , inFile) ;

    //discarding first starting  continues packets;
  //  int start_row_loc=0;
    // int start_row_loc=0;
  // if(datainfo.getModeFlag ()==IM){
    int  start_row_loc = -1 ;
    int index = -1 ;
    do
    {

        start_row_loc ++ ;
        index ++ ;
        if (index > nrow_l1)
        {
            LOG (ERROR) << "No starting found for data  ,packet id  ->" << pktSeqCtrl_temp[index] << endl ;
            ;
            break ;
        }
    }
    while (pktSeqCtrl_temp[index] != 1 && datainfo.getModeFlag ()==IM ) ;

    //}
    delete[] pktSeqCtrl_temp ;
    
    LOG(INFO) << "Starting from row number " << start_row_loc+1 << "  in level1 science data file "  ;
    for (long i = start_row_loc + 1 ; i <= nrow_l1 ; i = i + BUFFERSIZE)
    {
       start_rowNo = i ;
    
        felement = 1 ;
        if ((nrow_l1 + 1 - i) >= BUFFERSIZE) nelement = BUFFERSIZE ;
        else nelement = nrow_l1 + 1 - i ;
        //  cout<<"File "<<frow;

      
        /*-----------------reading data---------------------*/
        fits_read_col (di_input_file_ptr , TUSHORT , 8 , start_rowNo , felement , nelement , NULL , pktSeqCtrl , NULL , &status) ;
        printError (status , "Error in reading the column " , inFile) ;
        fits_read_col (di_input_file_ptr , TUSHORT , 11 , start_rowNo , felement , nelement , NULL , frameCount , NULL , &status) ;
        printError (status , "Error in reading the column value of the framecount" , inFile) ;
        fits_read_col (di_input_file_ptr , TUINT , 12 , start_rowNo , felement , nelement , NULL , frameTime , NULL , &status) ;
        printError (status , "Error in reading the column of the frametime" , inFile) ;
        fits_read_col (di_input_file_ptr , TBYTE , 24 , start_rowNo , felement , nelement*GTI_BYTES , NULL , GTI_arr , NULL , &status) ;
        printError (status , "Error in reading the column of the frametime" , inFile) ;

            /**loop for finding missing frame**/
        long frameindex ;
        double tmptime1 , tmptime2 , missingtime ;
        if (i == start_row_loc + 1)
        {
            frameindex = frameCount[0] ; //store first framecount in science data 'frameindex'
        }

        if (frameCount[0] != 1 && i == start_row_loc + 1) //to check missing framecounts before the first frame in science data
        {
            for (int i = 1 ; i < frameCount[0] ; i ++) //to store missing framecounts before the first frame in science data
            {
                missingFrameCounts.push_back (i) ;
                missingFrameTimes.push_back (0.0) ;
                count_missingFrame_beforeStart ++ ; //counter for missing frame numbers before first frame in science data
            }
            before_start_flag = TRUE ; //record that there are missing frames before first frame
        }
        LOG(INFO)<<missingFrameCounts.size ()<<endl;
        for (long j = 1 ; j <= nelement ; j ++) //to find missing frame counts between the data
        {
            while (frameCount[j - 1] == frameCount[j]) j ++ ; //compare each framecount with its previous one, increase the loop variable till it is same
            if (j >= nelement) break ; //break when loop variable exceeds number of rows read from data
            if (frameindex != frameCount[j] - 1) //
            {
                tmptime1 = frameTime[j - 1] ;
                tmptime2 = frameTime[j] ;

                missingtime = (tmptime1 + tmptime2) / 2.0 ;
                for (int m = frameindex + 1 ; m < frameCount[j] ; m ++)
                    missingFrameCounts.push_back (m) ;
                missingFrameTimes.push_back (missingtime) ;
            }

            frameindex = frameCount[j] ;
        }
 LOG(INFO)<<missingFrameCounts.size ()<<endl;

        /*--------------writing data to output file----------------------*/
        fits_write_col (di_output_file_ptr , TUSHORT , 1 ,frow , felement , nelement , pktSeqCtrl , &status) ;
        printError (status , "Error in writing the column of the pktSeqCtrl" , outFile) ;
        fits_write_col (di_output_file_ptr , TUSHORT , 2 , frow , felement , nelement , frameCount , &status) ;
        printError (status , "Error in writing the column of the frameCount" , outFile) ;

        if (datainfo.getModeFlag () == PC) //for PC Mode
        {
            fits_read_col (di_input_file_ptr , TBYTE , 15 , start_rowNo , felement , nelement*CENTROID_BYTES , NULL , centroids , NULL , &status) ;
            printError (status , "Error in reading the column of the centroids in inputFile" , inFile) ;
            fits_write_col (di_output_file_ptr , TBYTE , 4 ,frow , felement , nelement*CENTROID_BYTES , centroids , &status) ;
            printError (status , "Error in writing the column of the centroids in output File" , outFile) ;
        }
        else //for Integration Mode
        {
            fits_read_col (di_input_file_ptr , TUSHORT , 15 , start_rowNo , felement , nelement*PIXELNO , NULL , pixel , NULL , &status) ;
            printError (status , "Error in reading the column of the pixels from  input File" , inFile) ;
            fits_write_col (di_output_file_ptr , TUSHORT , 4 , frow , felement , nelement*PIXELNO , pixel , &status) ;
            printError (status , "Error in writing the column of the pixels in outFile" , outFile) ;
        }
        /*-----------------writing data-----------------------*/
        for (int j=0;j<BUFFERSIZE;j++){
            frameTime_temp[j]=frameTime[j]/1000.0;
        }
        
         fits_write_col (di_output_file_ptr , TDOUBLE , 3 , frow , felement , nelement , frameTime_temp , &status) ;
          printError (status , "Error in writing the column of the pixels in outFile" , outFile) ;
        for (long k = 0 ; k < nelement * GTI_BYTES ; k ++)
        {

            GTI_vect.push_back (GTI_arr[k]) ;

        }
        for (long k = 0 ; k < nelement ; k ++)
        {
            frm_count.push_back (frameCount[k]) ;
            frm_time_int.push_back (frameTime[k]) ;
            FrameTimeMap.insert (make_pair (frameCount[k] , frameTime[k])) ;
            //            GTI_vect.push_back (GTI_arr[k]);
            
        }
        //   cout<<"frame time is "<<frm_time_int[0];exit(1);
        //FrameTimeMap.
            frow=frow+BUFFERSIZE;
    }  //end of i loop
    LOG (INFO) << "Added columns PacketSequenceControl, FrameCount and  Pixel/Centroid to " << outFile << " file" ;

    //Write HISTORY to output frame file

    writeCommonKeywords (di_output_file_ptr , modulename) ;
    fits_close_file (di_output_file_ptr , &status) ; //closing data ingest output file
    printError (status , "Error in closing the file" , outFile) ;


    //discarding  for  first n Continues frames

    fits_close_file (di_input_file_ptr , &status) ;
    printError (status , "Error in closing the file" , inFile) ;
    /*---------------Creating Missing Frame List file-----------------*/
    fitsfile *missingFrameList_ptr ;
    fits_create_file (&missingFrameList_ptr , missingFrameList , &status) ;
    printError (status , "Error in creating the file " , missingFrameList) ;
    fits_create_img (missingFrameList_ptr , BYTE_IMG , 0 , NULL , &status) ;
    printError (status , "Error in creating the image " , missingFrameList) ;

    copyUserKeywords (inputfile_primary_ptr , missingFrameList_ptr) ;

    char *ttype_m[] = {"FrameCount" , "Time"} ;
    char *tform_m[] = {"1U" , "D"} ;
    char *tunit_m[] = {"" , ""} ;

    fits_create_tbl (missingFrameList_ptr , BINARY_TBL , 0 , 2 , ttype_m , tform_m , tunit_m , "MISSING_FRAME_NO" , &status) ;
    printError (status , "Error in creating the table in missing frame list" , missingFrameList) ;
    fits_movnam_hdu (missingFrameList_ptr , ANY_HDU , "MISSING_FRAME_NO" , 0 , &status) ;
    printError (status , "Error in moving to the MISSING _FRAME_NO HDU" , missingFrameList) ;
    fits_write_col (missingFrameList_ptr , TUSHORT , 1 , 1 , 1 , missingFrameCounts.size () , missingFrameCounts.data () , &status) ;
    printError (status , "Error in writing the column of the missingFrameCount" , missingFrameList) ;
    LOG (INFO) << missingFrameList << " file created and data written to it" ;
    fits_close_file (inputfile_primary_ptr , &status) ;

    /*------------fitting for time correction------------------*/
    nframes = frm_count.size () ;

    spMatrix A (nframes , 2) ; //Matrix A will have 1 in first column and frameCount in column 2 , values will be

    //populated from the level-1 science data file from which 4 columns has already
    //been copied to dataIngest output file
    spMatrix B (nframes , 1) ; //will contain frame time which is read from level-1 science data file, it is single column matrix
    spMatrix X (2 , 1) ; //Resulting coefficients , first value will offset(factor to be added) and second will be slope (multiplying factor)
    LOG(INFO)<<"Number of frames "<<nframes<<endl;
    if((frm_count[nframes-1]-frm_count[0])>1)
    {
    
    for (long i = 0 ; i < nframes ; i ++)
    {
        A (i , 0) = 1 ;
        A (i , 1) = frm_count[i] ; //assigning framecount to A
        B (i , 0) = frm_time_int[i] ; //assigning time to B
        B (i , 0) = B (i , 0) / 1000.0 ; //for conversion to seconds
    }
    LOG (INFO) << "Applying least square for time correction......." ;
    X.ApplyLeastSquare (A , B) ;
    }
    else{
        X(0,0)=0;
        X(1,0)=1;
    }
    unsigned short temp_time_fitting = 0 ;
    if (before_start_flag)
    {
        for (int i = 0 ; i < count_missingFrame_beforeStart ; i ++)
        {
            temp_time_fitting = missingFrameCounts[i] ;
            missingFrameTimes[i] = X (0 , 0)+ (X (1 , 0) * temp_time_fitting) ;
        }

    }
   // LOG(INFO)<<"NNN "<<missingFrameCounts.size (); exit(1);
    delete[] pktSeqCtrl , centroids , frameTime , frameCount ;
    fits_write_col (missingFrameList_ptr , TDOUBLE , 2 , 1 , 1 , missingFrameCounts.size () , missingFrameTimes.data () , &status) ;
    printError (status , "Error in Writing the column of missingFrameTimes" , missingFrameList) ;
    if (history == YES) writeHistory (missingFrameList , vhistorystr) ;
    writeCommonKeywords (missingFrameList_ptr , modulename) ; //writes date, origin, checksum and creator to file
    fits_close_file (missingFrameList_ptr , &status) ;
    printError (status , "Error in closing the missing frame list file" , missingFrameList) ;

    /*----------------computing corrected time and writing to file-----------------*/

    double fittedtime , temptime ;
    unsigned short tempfcount ;
    felement = 1 ;
    nelement = 1 ;
  
    for (long i = 1 ; i <= nframes ; i ++)
    {
        tempfcount = frm_count[i - 1] ;
        temptime = frm_time_int[i - 1] ;

        temptime = temptime / 1000 ; //for  converting time to seconds.

        fits_write_col (time_correction_file_ptr , TUSHORT , 1 , i , felement , nelement , &tempfcount , &status) ;
        printError (status , "Error in Writing the column of FrameCount" , timeCorrectionFile) ; // 1st column containing frame count
        fits_write_col (time_correction_file_ptr , TDOUBLE , 2 , i , felement , nelement , &temptime , &status) ;
        printError (status , "Error in Writing the column of FrameTime(converting in to  milli seconds" , timeCorrectionFile) ; //2nd column containing time read from file

       fittedtime = X (0 , 0)+ (X (1 , 0) * tempfcount) ;

        fits_write_col (time_correction_file_ptr , TDOUBLE , 3 , i , felement , nelement , &fittedtime , &status) ;
        printError (status , "Error in writing the column of the fitted Time" , timeCorrectionFile) ; //3rd column containing time after fitting
        frm_time_double.push_back (fittedtime) ; //updating map with fitted time
        
    }
  
    //writes date, origin, checksum and creator to file
    fits_close_file (time_correction_file_ptr , &status) ; //closing time correction fits file
    printError (status , "Error in closing the File" , timeCorrectionFile) ;
    LOG (INFO) << "Applying UTC correction......." ;
   
//    /**applying UTC correction**/
//    status = doUTCCorrection () ; //
//    if (status)
//    {
//        LOG (ERROR) << "Error in UTC correction" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//    LOG (INFO) << "Time Correction completed" << endl ;
double  t1 = frm_time_double[0] ;
    double f1 = frm_count[0] ;
    double t2 = frm_time_double[nframes - 1] ;
    double f2 = frm_count[nframes - 1] ;
  double integrationTime = (double) (t2 - t1) / (double) (f2 - f1 + 1) ;
updateKeywords (outFile , 2 ,1 , TDOUBLE , "INT_TIME" , &integrationTime) ;

    /*--------------CRC Correction-----------------*/
    /**Method For the CRC check**/
    doCRCCorrection () ;

    //Generating LEVEL2 GTI file(LEVEL1 code to be integrated)
    LOG (INFO) << "Generating Level-2 GTI file using LBT and GTI file" ;
    vector<long> rowno ;
    //gti_flag=1;
    if (gti_flag)
    {
        if (doGTIfiltering (rowno))
        {
            LOG (ERROR) << "Error in GTI filtering" << endl ;
            return (EXIT_FAILURE) ;
        }
    }
    else
    {
        LOG (ERROR) << "GTI Filtering has not been done" ;
    }
    LOG (INFO) << "GTI Failed for   " << rowno.size () << " rows : " ;


    //deleting rows where GTI failed.
    fits_open_file (&di_output_file_ptr , outFile , READWRITE , &status) ;
    printError (status , "Error in opening the  output file " , outFile) ;

    fits_movnam_hdu (di_output_file_ptr , BINARY_TBL , SCIENCEDATA_HDUNAME , 0 , &status) ;
    printError (status , "Error in moving to particular HDU " , outFile) ;

    fits_delete_rowlist (di_output_file_ptr , rowno.data () , rowno.size () , &status) ;
    printError (status , "Error in Deleting the row list " , outFile) ;

    fits_close_file (di_output_file_ptr , &status) ;
    printError (status , "Error in closing the output file " , outFile) ;
     LOG(INFO)<<"The File "<<tempinfo<<endl;
        LOG(INFO)<<"The PCImage File "<<PCImageFile<<endl;
    //Checking DATAMODE for frames/event list creation    
    if (datainfo.getModeFlag () == IM)
    {
        LOG (INFO) << "Creating separate frames for Integration Mode......." ;

        if (createFrames ()) //Function to extract individual frames in case of IM mode
        {
            return (EXIT_FAILURE) ;
        }

    }
    else //for PC mode
    {
        LOG (INFO) << "Creating event data for Photon Counting Mode......." ;
        if (getEvents ()) //Decode events in case of PC mode
        {
            LOG(ERROR)<<"Error in decoding events ";
            return (EXIT_FAILURE) ;
        } //function creates event file and image file for PC mode
    }
    datainfo.setXsize (IMG_DIM_DI) ;
    datainfo.setYsize (IMG_DIM_DI) ;
   LOG (INFO) << "Reading Dark frames........" ;

//    if (readDarkFrames ()) //Function  to read Dark frame time and filename
//    {
//        return (EXIT_FAILURE) ;
//    }

    /**Method for INFO file generation. Information File contains  the listing of all the frames  generated in this module and keyword necessary for next module**/
  
 //  LOG(INFO)<<"The File "<<infoFilename<<endl;
      //LOG(INFO)<<"The PCImage File "<<PCImageFile<<endl;exit(1);
    if (createInfoFile ())
    {
        LOG (ERROR) << "***Error creating info file***" << endl ;
        return (EXIT_FAILURE) ;
    }

    if (history == YES) writeHistory (infoFilename , vhistorystr) ;

    LOG (INFO) << "Data Ingest process is completed successfully." << endl ;

    return 0 ;
}

//function does CRC correction and produces out CRC failed list file


void DataIngest::doCRCCorrection ()
{
    int status = 0 ;
    LOG (INFO) << "Applying CRC filtering......." ;
    fitsfile *di_input_file_ptr , *di_output_file_ptr ;
    fitsfile *crc_log_file_ptr ;

    //opening LEVEL1 input File

    fits_open_file (&di_input_file_ptr , inFile , READONLY , &status) ;
    printError (status , "Error in the opening the input File" , inFile) ;
    fits_movnam_hdu (di_input_file_ptr , BINARY_TBL , SCIENCEDATA_HDUNAME , 0 , &status) ;
    printError (status , "Error in moving to the SCIENCE DATA HDU in the input File " , inFile) ;

    //Creating output CRC log file
    fits_create_file (&crc_log_file_ptr , crcFailedList , &status) ;
    printError (status , "Error in creating the File CRCFailed List" , crcFailedList) ;
    fits_create_img (crc_log_file_ptr , BYTE_IMG , 0 , NULL , &status) ;
    printError (status , "Error in creating the image for  the CRC file" , crcFailedList) ;

    writeUsrkeywordsFrmvect (crcFailedList , key_record) ;

    char *ttype12[] = {"Frame Count" , "Actual CRC" , "Computed CRC"} ;
    char *tform12[] = {"U" , "U" , "U"} ;
    char *tunit12[] = {"" , "" , ""} ;
    fits_create_tbl (crc_log_file_ptr , BINARY_TBL , 0 , 3 , ttype12 , tform12 , tunit12 , "FAILED CRC" , &status) ;
    printError (status , "Error in creating the Tablein crc File" , crcFailedList) ;
    fits_movnam_hdu (crc_log_file_ptr , ANY_HDU , "FAILED CRC" , 0 , &status) ;
    printError (status , "Error in moving to the FAILED CRC  HDU" , crcFailedList) ;

    unsigned char Readarr[PACKET_SIZE] , Bytes2int[PACKET_SIZE] , a1 , b1 ;
    int Bytes2int_INT[PACKET_SIZE] ;
    unsigned int pid , psc , seg_len , sec_hdr ;
    unsigned short pid1 , psc1 , seg_len1 , sec_hdr1 , frmcount1 ;
    unsigned short pid2 , psc2 , seg_len2 , sec_hdr2 , frmcount2 ;
    long frmtime1 , frmtime2 , frmtime3 , frmtime4 , frmtime ;
    long ticktime1 , ticktime2 , ticktime3 , ticktime4 , ticktime ;
    long tickcount1 , tickcount2 , tickcount3 , tickcount4 , tickcount ;

    unsigned short *frameCount = new unsigned short[nrow_l1] ;
    checkMemoryAvailability (frameCount , "frameCount") ;
    //short *actualCRC = new short[nrow_l1] ;
    //checkMemoryAvailability (actualCRC , "actualCRC") ;
    // unsigned   short *computedCRC = new short[nrow_l1] ;
    // checkMemoryAvailability (computedCRC , "computedCRC") ;
    unsigned short *actualCRC = new unsigned short[nrow_l1] ;
    checkMemoryAvailability (actualCRC , "actualCRC") ;
    unsigned short *computedCRC = new unsigned short[nrow_l1] ;
    checkMemoryAvailability (computedCRC , "computedCRC") ;
    unsigned char read_2016[CENTROID_BYTES] ;
    unsigned short read_1008[PIXELNO] ;
    //short read_1008[PIXELNO] ;
    //list of rows to be deleted if CRC fails
    vector<long> list1 , list2 ;
    unsigned short crc_ccitt_ffff , crc_ccitt_0000 ;
    //  short crc_ccitt_ffff , crc_ccitt_0000 ;
    crc_ccitt_ffff = 0xffff ;


    //  short crc_ccitt_ffff , crc_ccitt_0000 ;
    //  crc_ccitt_ffff = 0xffff ;
    //reading all at a time as array is declared for them
    fits_read_col (di_input_file_ptr , TUSHORT , 11 , 1 , 1 , nrow_l1 , NULL , frameCount , NULL , &status) ;
    printError (status , "Error in reading the frame Count column" , inFile) ;
    fits_read_col (di_input_file_ptr , TUSHORT , 16 , 1 , 1 , nrow_l1 , NULL , actualCRC , NULL , &status) ;
    printError (status , "Error in reading the actual CRC column " , inFile) ;
    int temp_ans ;
    int temp2 ;
    int temp_1 ;
    long nelement = 1 , felement = 1 ;
    list.clear () ;
    //ofstream fout ("crcFile2" , ios::out) ;
    // fout.setf (ios::hex,ios::basefield);
    for (long frow = 1 ; frow <= nrow_l1 ; frow ++)
    {
        fits_read_col (di_input_file_ptr , TUSHORT , 7 , frow , felement , nelement , NULL , &pid , NULL , &status) ;
        printError (status , "Error in reading the column of the pid For CRC check" , inFile) ;
        fits_read_col (di_input_file_ptr , TUSHORT , 8 , frow , felement , nelement , NULL , &psc , NULL , &status) ;
        printError (status , "Error in reading the column of the psc For CRC check" , inFile) ;
        fits_read_col (di_input_file_ptr , TUSHORT , 9 , frow , felement , nelement , NULL , &seg_len , NULL , &status) ;
        printError (status , "Error in reading the column of the seg_len For CRC check" , inFile) ;
        fits_read_col (di_input_file_ptr , TUSHORT , 10 , frow , felement , nelement , NULL , &sec_hdr , NULL , &status) ;
        printError (status , "Error in reading the column of the sec_hdr For CRC check" , inFile) ;
        fits_read_col (di_input_file_ptr , TUINT , 12 , frow , felement , nelement , NULL , &frmtime , NULL , &status) ;
        printError (status , "Error in reading the column of the frmtime For CRC check" , inFile) ;
        fits_read_col (di_input_file_ptr , TUINT , 13 , frow , felement , nelement , NULL , &ticktime , NULL , &status) ;
        printError (status , "Error in reading the column of the ticktime For CRC check" , inFile) ;
        fits_read_col (di_input_file_ptr , TUINT , 14 , frow , felement , nelement , NULL , &tickcount , NULL , &status) ;
        printError (status , "Error in reading the column of the tick_count For CRC check" , inFile) ;

        if (datainfo.getModeFlag () == PC)
        {
            fits_read_col (di_input_file_ptr , TBYTE , 15 , frow , felement , nelement*CENTROID_BYTES , NULL , read_2016 , NULL , &status) ;
            printError (status , "Error in reading the column of the centroid bytes  For CRC check" , inFile) ;
        }
        else
        {
            fits_read_col (di_input_file_ptr , TUSHORT , 15 , frow , felement , nelement*PIXELNO , NULL , read_1008 , NULL , &status) ;
            printError (status , "Error in reading the column of the pixels For CRC check" , inFile) ;
        }
        //    cout<<read_1008[0]<<" "<<read_1008[1]<<" "<<read_1008[2];exit(1);
        pid1 = (pid & 0xFF00) >> 8 ;
        pid2 = (pid & 0x00FF) ;
        psc1 = (psc & 0xFF00) >> 8 ;
        psc2 = (psc & 0x00FF) ;
        seg_len1 = (seg_len & 0xFF00) >> 8 ;
        seg_len2 = (seg_len & 0x00FF) ;
        sec_hdr1 = (sec_hdr & 0xFF00) >> 8 ;
        sec_hdr2 = (sec_hdr & 0x00FF) ;
        frmcount1 = (frameCount[frow - 1] & 0xFF00) >> 8 ;
        frmcount2 = (frameCount[frow - 1] & 0x00FF) ;
        frmtime1 = (frmtime & 0xFF000000) >> 24 ;
        frmtime2 = (frmtime & 0x00FF0000) >> 16 ;
        frmtime3 = (frmtime & 0x0000FF00) >> 8 ;
        frmtime4 = (frmtime & 0x000000FF) ;
        ticktime1 = (ticktime & 0xFF000000) >> 24 ;
        ticktime2 = (ticktime & 0x00FF0000) >> 16 ;
        ticktime3 = (ticktime & 0x0000FF00) >> 8 ;
        ticktime4 = (ticktime & 0x000000FF) ;
        tickcount1 = (tickcount & 0xFF000000) >> 24 ;
        tickcount2 = (tickcount & 0x00FF0000) >> 16 ;
        tickcount3 = (tickcount & 0x0000FF00) >> 8 ;
        tickcount4 = (tickcount & 0x000000FF) ;

        for (int i = 0 ; i < PACKET_SIZE ; i ++)
        {
            Readarr[i] = 0 ;
            Bytes2int[i] = 0 ;
        }

        Readarr[0] = pid1 ;
        Readarr[1] = pid2 ;
        Readarr[2] = psc1 ;
        Readarr[3] = psc2 ;
        Readarr[4] = seg_len1 ;
        Readarr[5] = seg_len2 ;
        Readarr[6] = sec_hdr1 ;
        Readarr[7] = sec_hdr2 ;
        Readarr[8] = frmcount1 ;
        Readarr[9] = frmcount2 ;
        Readarr[10] = frmtime1 ;
        Readarr[11] = frmtime2 ;
        Readarr[12] = frmtime3 ;
        Readarr[13] = frmtime4 ;
        Readarr[14] = ticktime1 ;
        Readarr[15] = ticktime2 ;
        Readarr[16] = ticktime3 ;
        Readarr[17] = ticktime4 ;
        Readarr[18] = tickcount1 ;
        Readarr[19] = tickcount2 ;
        Readarr[20] = tickcount3 ;
        Readarr[21] = tickcount4 ;

        if (datainfo.getModeFlag () == PC)
        {
            for (int mn = 0 ; mn < CENTROID_BYTES ; mn ++)
                Readarr[mn + 26] = read_2016[mn] ;
        }
        else
        {
            for (int i = 0 , mn = 26 ; i < PIXELNO ; i ++ , mn = mn + 2)
            {
                unsigned short msb = (read_1008[i] >> 8) ;
                unsigned short lsb = read_1008[i] & 0xff ;
                Readarr[mn] = (unsigned char) msb ;
                Readarr[mn + 1] = (unsigned char) lsb ;

            }
        }

        for (int bytes = 0 ; bytes < 2042 ; bytes ++)
            Bytes2int[bytes] = Readarr[bytes] ;

        //        for (int bytes = 0 ; bytes < 2042 ; bytes=bytes +1)
        //        {
        //             Bytes2int_INT[bytes]=Bytes2int[bytes];
        //        }
        //remove
//        char p ;
//        int ans_int = 0 ;
  //      char temp_str[FLEN_FILENAME] ;
        //       vector<int> Eight_BITs;
        //        for (int bytes = 0 ; bytes < 2042 ; bytes=bytes +1){
        //           ans_int=(int)Bytes2int[bytes];
        //           //cout<<ans_int;
        //           //exit(1);
        //           // ans_int=100;
        //           if(ans_int<16)
        //           {
        //               if(ans_int==10) {p='A';  sprintf (temp_str,"%c",p);    }
        //               else if  (ans_int==11){ p='B';  sprintf (temp_str,"%c",p);    }
        //              else if  (ans_int==12) {p='C';  sprintf (temp_str,"%c",p);    }
        //              else if  (ans_int==13) {p='D';  sprintf (temp_str,"%c",p);    }
        //              else if  (ans_int==14){ p='E';  sprintf (temp_str,"%c",p);    }
        //              else if  (ans_int==15){ p='F';  sprintf (temp_str,"%c",p);    }
        //              
        //              else{
        //              sprintf (temp_str,"%d",ans_int);    
        //              }
        //               fout<<temp_str;
        //           }
        //           else{
        //                fout<<ans_int;
        //           }
        //           
        //            
        ////            for (int i=1;i<=8;i++){
        ////                p=((int)Bytes2int[bytes]>>(8-i)) && 0x0001 ;
        ////                for(int j=1;j<=4;j++){    
        ////                ans_int=(ans_int|p)<<(4-j);
        ////                }
        //////                Eight_BITs.push_back (p);              
        ////            }
        ////            fout<<Bytes2int[bytes];
        ////             fout<<p;
        //        
        //        }

        // fout;
        // fout.close ();exit(1);
        // fout<<temp_ans;
        crc_ccitt_ffff = 0xffff ;

        //Compute CRC for 2042 bytes of a datapacket and store in  'computedCRC[frow - 1]' for row number 'frow'
        for (int a = 0 ; a < 2042 ; a ++)
        {
            crc_ccitt_ffff = update_crc_ccitt (crc_ccitt_ffff , (char) Bytes2int[a]) ;
        }

        computedCRC[frow - 1] = crc_ccitt_ffff ;

    }
    // cout<<computedCRC[0]<<" "<<computedCRC[1]<<" "<<computedCRC[2];exit(1);
    //comparing Computed CRC with Actual CRC and storing failed row indices  in list
    for (long frow = 1 ; frow <= nrow_l1 ; frow ++)
    {

        if (computedCRC[frow - 1] != actualCRC[frow - 1])
        {

            if (! dropFrame)
            {
                //list.push_back (frow -1) ;
                list.push_back (frow) ;
            }

            if (dropFrame)
            {
                int temprowno = frow - 1 ;
                while (frameCount[temprowno] == frameCount[temprowno - 1])
                {
                    list.push_back (temprowno - 1 + 1) ;
                    temprowno -- ;
                }
                temprowno = frow - 1 ;
                list.push_back (temprowno + 1) ;
                while (frameCount[temprowno] == frameCount[temprowno + 1])
                {
                    list.push_back (temprowno + 1 + 1) ;
                    temprowno ++ ;
                }
                frow = temprowno + 1 ;
            }

            //sorting the list
            if (dropFrame)
            {

                sort (list.begin () , list.end ()) ;

            }

            LOG (INFO) << "CRC failed at " << list.size () << " rows" ;
            //            if (nrow_l1 == 0)
            //            {
            //                LOG (ERROR) << "***Divide by Zero Error in in CRC check***" ;
            //            }
            //errorRate = list.size ()*100.0 / (double) (nrow_l1) ; //will give percentage of CRC fail
        }



    }

    LOG (INFO) << "Writing to CRC log file......" ;

    for (long crcfailedrow = 0 ; crcfailedrow < list.size () ; crcfailedrow ++)
    {
        fits_write_col (crc_log_file_ptr , TUSHORT , 1 , crcfailedrow + 1 , felement , 1 , &frameCount[list[crcfailedrow] - 1] , &status) ;
        printError (status , "Error in writing the column of the frameCount" , crcFailedList) ;

        fits_write_col (crc_log_file_ptr , TUSHORT , 2 , crcfailedrow + 1 , felement , 1 , &actualCRC[list[crcfailedrow] - 1] , &status) ;
        printError (status , "Error in writing the column of the actualCRC" , crcFailedList) ;

        fits_write_col (crc_log_file_ptr , TUSHORT , 3 , crcfailedrow + 1 , felement , 1 , &computedCRC[list[crcfailedrow] - 1] , &status) ;
        printError (status , "Error in writing the column of the computedCRC" , crcFailedList) ;
    }
    if (history == YES) writeHistory (crcFailedList , vhistorystr) ;

    writeCommonKeywords (crc_log_file_ptr , modulename) ; //writes date, origin, checksum and creator to file

    fits_close_file (crc_log_file_ptr , &status) ;
    printError (status , "Error in closing the file" , crcFailedList) ;

    delete[] actualCRC , computedCRC , frameCount ;

    fits_close_file (di_input_file_ptr , &status) ;
    printError (status , "Error in closing the File" , crcFailedList) ;

    LOG (INFO) << "Deleting rows where CRC failed from dataingest output file......" ;

    char tempstr[FLEN_FILENAME] ;

    sprintf (tempstr , "%s[DETECTOR_DATA]" , outFile) ;

    fits_open_file (&di_output_file_ptr , tempstr , READWRITE , &status) ;
    printError (status , "Error in opening the  output file in CRC check" , outFile) ;

    fits_delete_rowlist (di_output_file_ptr , list.data () , list.size () , &status) ;
    printError (status , "Error in Deleting the row list " , outFile) ;

    writeCommonKeywords (di_output_file_ptr , modulename) ;

    fits_close_file (di_output_file_ptr , &status) ; //closing data ingest output file
    printError (status , "Error in Closing the output File " , outFile) ;

    if (dropFrame)
    {
        list.clear () ;
    }
    cerr << endl ;

}

//function for performing UTC correction on the data using TCT file


int DataIngest::doUTCCorrection ()
{

    LOG (INFO) << "Time correction using TCT......" ;

    fitsfile *ftct ;
    int status = 0 ;

    LOG (INFO) << "Reading TCT file......" ;
    fits_open_file (&ftct , inFile_TCT , READONLY , &status) ;
    printError (status , "***Error opening TCT file***" , inFile_TCT) ;

    fits_movabs_hdu (ftct , 2 , NULL , &status) ;
    printError (status , "***Error opening TCT file***" , inFile_TCT) ;

    int colnum_detector ; //stores column number for columns UVITF, UVITN, UVITV depending on data detector

    if (strcasecmp (datainfo.getDetector () , "FUV") == 0) colnum_detector = TCT_FUV_COLNO ;
    else if (strcasecmp (datainfo.getDetector () , "NUV") == 0) colnum_detector = TCT_NUV_COLNO ;
    else if (strcasecmp (datainfo.getDetector () , "VIS") == 0) colnum_detector = TCT_VIS_COLNO ;
    else
    {
        LOG (ERROR) << "***DETECTOR keyword value is  " << datainfo.getDetector () << "\nExpected is FUV, NUV, VIS***" ;
        return (EXIT_FAILURE) ;
    }
    long nrow_tct ;
    fits_get_num_rows (ftct , &nrow_tct , &status) ;
    printError (status , "***Error getting rows from TCT file ***" , inFile_TCT) ;

    double *uvt_time = new double[nrow_tct] ;
    checkMemoryAvailability (uvt_time , "uvt_time") ;

    double *sps_time = new double[nrow_tct] ;
    checkMemoryAvailability (sps_time , "sps_time") ;

    long firstrow = 1 , firstelem = 1 ;

    /**reading SPS time column and UVT time column**/
    fits_read_col (ftct , TDOUBLE , 1 , firstrow , firstelem , nrow_tct , NULL , sps_time , NULL , &status) ;
    printError (status , "***Error  in reading column  value of the sps_time from TCT file ***" , inFile_TCT) ;

    fits_read_col (ftct , TDOUBLE , colnum_detector , firstrow , firstelem , nrow_tct , NULL , uvt_time , NULL , &status) ;
    printError (status , "***Error  in reading column  value of the uvt_time from TCT file ***" , inFile_TCT) ;

    fits_close_file (ftct , &status) ;
    printError (status , "Error in closing  the TCT File" , inFile_TCT) ;

    double temp_time_uvt , temp_time_sps ;

    //sorting UVT time and SPS time wrt to UVT time
    for (int i = 0 ; i < nrow_tct ; i ++)
    {
        for (int j = nrow_tct - 1 ; j > i ; j --)
        {
            if (uvt_time[i] > uvt_time[j])
            {
                temp_time_uvt = uvt_time[j] ;
                uvt_time[j] = uvt_time[i] ;
                uvt_time[i] = temp_time_uvt ;
                temp_time_sps = sps_time[j] ;
                sps_time[j] = sps_time[i] ;
                sps_time[i] = temp_time_sps ;
            }
        }

    }

    unsigned int t1 , t2 ; // to store start and end time
    unsigned short f1 , f2 ; //to store start and end framecount

    t1 = frm_time_double[0] ;
    f1 = frm_count[0] ;
    t2 = frm_time_double[nframes - 1] ;
    f2 = frm_count[nframes - 1] ;
    
    bool tct_flag = FALSE ;

    //computing UTC time for frames using UVT time and SPS time from TCT
    // cout<<frm_time_double[0];exit(1);
    LOG(INFO)<<"Printing Frame time after time correction for missing frames" << endl;
 
    
    
    for (long i = 0 ; i < nframes ; i ++)
    {
       //cout<<"INDDDD "<<uvt_time[1]<<" "<<frm_time_double[i];exit(1);
        if(frm_time_double[i]<uvt_time[0])
        {            
             tct_flag = TRUE ;
            frm_time_double[i] = sps_time[0]+((sps_time[1] - sps_time[0]) / (uvt_time[1] - uvt_time[0]))*(frm_time_double[i] - uvt_time[0]) ;     
        }
        else if(frm_time_double[i]>uvt_time[nrow_tct-1])
        {            
            tct_flag = TRUE ;
            frm_time_double[i] = sps_time[nrow_tct-2]+((sps_time[nrow_tct-1] - sps_time[nrow_tct-2]) / (uvt_time[nrow_tct-1] - uvt_time[nrow_tct-2]))*(frm_time_double[i] - uvt_time[nrow_tct-2]) ;            
        }
        else
        {
        for (int index = 1 ; index < nrow_tct ; index ++)
        {
          
            if (frm_time_double[i] < uvt_time[index] && frm_time_double[i] >= uvt_time[index - 1])
            {
               // LOG(INFO)<<index<< " "<<i<<"  " <<setprecision (12)<<sps_time[index-1]<<" "<<setprecision (12)<<sps_time[index]<<" "<<uvt_time[index-1]<<" " <<uvt_time[index]<<" "<<frm_time_double[i];
                tct_flag = TRUE ;
                frm_time_double[i] = sps_time[index - 1]+((sps_time[index] - sps_time[index - 1]) / (uvt_time[index] - uvt_time[index - 1]))*(frm_time_double[i] - uvt_time[index - 1]) ;
             //   new_time[i] = sps_time[index - 1]+((sps_time[index] - sps_time[index - 1]) / (uvt_time[index] - uvt_time[index - 1]))*(frm_time_double[i] - uvt_time[index - 1]) ;
              
                break ;
            }

        }
        }        
       
    }
    if (! tct_flag)
    {
        LOG (INFO) << "***Data time and TCT time does not match***" ;
        return (EXIT_FAILURE) ;
    }
    
    double tstarti , tstartf , tstopi , tstopf ;
    
   /// long temp_One;
  //  temp_One=(long ) frm_time_double[0] ;
    tstarti = (long ) frm_time_double[0] ;
    tstartf = frm_time_double[0]-(long) frm_time_double[0] ;
    tstopi = (long) frm_time_double[nframes - 1] ;
    tstopf = frm_time_double[nframes - 1]-(long) frm_time_double[nframes - 1] ;
//    tstarti = (long ) frm_time_double[0] ;
//    tstartf = frm_time_double[0]-(long) frm_time_double[0] ;
//    tstopi = (long) frm_time_double[nframes - 1] ;
//    tstopf = frm_time_double[nframes - 1]-(long) frm_time_double[nframes - 1] ;
  //  setprecision (20);
    //LOG(INFO).
            

    LOG (INFO) << "\033[1;34mtstarti ->" <<setprecision (20)<< tstarti<<"\033[0m" ;
           LOG (INFO) << "\033[1;34mtstartf -> " <<setprecision (20)<< tstartf<<"\033[0m";
            LOG (INFO)<< "\033[1;34mtstopi -> " <<setprecision (20)<< tstopi<<"\033[0m" ;
            LOG (INFO)<< "\033[1;34mtstopf -> " <<setprecision (20)<< tstopf<<"\033[0m" ;

    datainfo.setTime (tstarti , tstartf , tstopi , tstopf) ; //updating datainfo

    LOG (INFO) << "Updating time correction file....." ;

    fitsfile *ftime ;

    fits_open_file (&ftime , timeCorrectionFile , READWRITE , &status) ;
    printError (status , "Error in opening the timeCorrection File" , timeCorrectionFile) ;

    fits_movabs_hdu (ftime , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU " , timeCorrectionFile) ;

    char *ttype = {"UTC"} ;
    char *tform = {"D"} ;
    fits_insert_col (ftime , 4 , ttype , tform , &status) ; //inserting 3rd column
    printError (status , "Error in inserting the 4th column in the time Correction File" , timeCorrectionFile) ;
    fits_write_col (ftime , TDOUBLE , 4 , firstrow , firstelem , nframes , frm_time_double.data () , &status) ;
    printError (status , "Error in Writing the column new Time" , timeCorrectionFile) ;


    updateKeywords (timeCorrectionFile , 2 , 4 , TDOUBLE , "TSTARTI" , &tstarti ,
            TDOUBLE , "TSTARTF" , &tstartf ,
            TDOUBLE , "TSTOPI" , &tstopi ,
            TDOUBLE , "TSTOPF" , &tstopf) ;

    if (history == YES) writeHistory (timeCorrectionFile , vhistorystr) ;

    writeCommonKeywords (ftime , modulename) ;

    fits_close_file (ftime , &status) ;
    printError (status , "Error in closing the File" , timeCorrectionFile) ;

    /*---------updating time correction fits file completed----------*/

    /*----------------------finding integration time----------------------*/
    // divide by Zero checking
    if ((f2 - f1 + 1) == 0)
    {

        LOG (ERROR) << "***Divide by Zero Error in Finding the Integration Time***" ;
        return (EXIT_FAILURE) ;
    }
    double integrationTime = (double) (t2 - t1) / (double) (f2 - f1 + 1) ;
    //integrationTime = integrationTime / 1000 ; //to convert into seconds
    LOG(INFO) << "Integration time is " << integrationTime << " seconds" ;
    datainfo.setIntegrationTime (integrationTime) ; //updating datainfo

    LOG(INFO) << "Updating output file...." ;

    /*-----updating data ingest output file------------*/
    fitsfile *fdataingestout ;

    fits_open_file (&fdataingestout , outFile , READWRITE , &status) ;
    printError (status , "Error in opening the output File" , outFile) ;

    updateKeywords (outFile , 1 , 5 , TDOUBLE , "TSTARTI" , &tstarti ,
            TDOUBLE , "TSTARTF" , &tstartf ,
            TDOUBLE , "TSTOPI" , &tstopi ,
            TDOUBLE , "TSTOPF" , &tstopf ,
            TDOUBLE , "INT_TIME" , &integrationTime) ;

    fits_movabs_hdu (fdataingestout , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU" , outFile) ;

    long nrow_dout ;

    fits_get_num_rows (fdataingestout , &nrow_dout , &status) ;
    printError (status , "Error in getting the number of rows in output File" , outFile) ;

    // reading time column, cloumn number is 3 for time column
//    fits_write_col (fdataingestout , TDOUBLE , 3 , 1 , firstelem , nframes , frm_time_double.data () , &status) ;
//    printError (status , "Error in writing the column of the Frametime (in Double)" , outFile) ;

     fits_write_col (fdataingestout , TDOUBLE , 3 , 1 , firstelem , nrow_dout , frm_time_double.data () , &status) ;
  printError (status , "Error in writing the column of the Frametime (in Double)" , outFile) ;
    
    // for second hdu
    updateKeywords (outFile , 2 , 5 , TDOUBLE , "TSTARTI" , &tstarti ,
            TDOUBLE , "TSTARTF" , &tstartf ,
            TDOUBLE , "TSTOPI" , &tstopi ,
            TDOUBLE , "TSTOPF" , &tstopf ,
            TDOUBLE , "INT_TIME" , &integrationTime) ;

    fits_close_file (fdataingestout , &status) ;
    printError (status , "Error in closing the File" , outFile) ;


    /*---------updating time in data ingest output file complete-----------*/
    return (EXIT_SUCCESS) ;
}

//function to extract events in case of PC mode data


int DataIngest::getEvents ()
{
    LOG (INFO) << "Decoding Events....." ;
    ;
    fitsfile *fptr , *fout , *fptr_primary , *fout_h3 , *fimg ;
    int status = 0 ;
    strcpy(tempinfo,infoFilename);
    vector<string> vhistorystr ; //Vector for storing HISTORY
    getHistory (vhistorystr) ; //function to create history for output event file

    //opening dataIngestout file
    fits_open_file (&fptr_primary , outFile , READWRITE , &status) ; //for pointing to primary header
    printError (status , "Error in opening the output File" , outFile) ;

    fits_open_file (&fptr , outFile , READONLY , &status) ; // for pointing to second header
    printError (status , "Error in opening the output File" , outFile) ;
    fits_movnam_hdu (fptr , BINARY_TBL , SCIENCEDATA_HDUNAME , 0 , &status) ; //moving to HDU containing science data
    printError (status , "Error moving to the  SCIENCEDATA HDU in getEvent() Function" , outFile) ;

    long nrows = 0 ;
    fits_get_num_rows (fptr , &nrows , &status) ; //getting number of rows from dataIngest out file
    printError (status , "Error getting the number of Rows in output File" , outFile) ;

    //centroiddata - buffer to store centroid data
    unsigned char *centroiddata = new unsigned char[nrows * CENTROID_BYTES] ;
    checkMemoryAvailability (centroiddata , "centroiddata") ;

    //psc - buffer to store packet sequence control
    unsigned short *psc = new unsigned short[nrows] ;
    checkMemoryAvailability (psc , "psc") ;

    //framecount -  buffer to store frame count
    unsigned short *framecount = new unsigned short[nrows] ;
    checkMemoryAvailability (framecount , "framecount") ;

    // time  - buffer to store frame time
    double *time = new double[nrows] ;
    checkMemoryAvailability (time , "time") ;

    //Reading packet Seq Ctrl ,framecount ,time, centroid data from dataIngest out file
    fits_read_col (fptr , TUSHORT , 1 , 1 , 1 , nrows , NULL , psc , NULL , &status) ;
    printError (status , "Error reading the column value of the psc in getEvents()" , outFile) ;
    fits_read_col (fptr , TUSHORT , 2 , 1 , 1 , nrows , NULL , framecount , NULL , &status) ;
    printError (status , "Error reading the column value of the framecount" , outFile) ;
    fits_read_col (fptr , TDOUBLE , 3 , 1 , 1 , nrows , NULL , time , NULL , &status) ;
    printError (status , "Error reading the column value of the time" , outFile) ;
    fits_read_col (fptr , TBYTE , 4 , 1 , 1 , nrows*CENTROID_BYTES , NULL , centroiddata , NULL , &status) ;
    printError (status , "Error reading the column value of the centroiddata" , outFile) ;

    //Creating output event file
    LOG (ERROR) << "Creating event file......." ;
    fits_create_file (&fout , eventFile , &status) ;
    printError (status , "Error creating the output Event file" , eventFile) ;
    fits_create_img (fout , BYTE_IMG , 0 , NULL , &status) ;
    printError (status , "Error creating the image " , eventFile) ;

    //Copying keywords from dataIngestout file to output EVENT file
    writeUsrkeywordsFrmvect (eventFile , key_record) ;
    //copyUserKeywords (fptr_primary , fout) ;

    //Declaring columns for output EVENT file
    char *ttype[] = {"PacketSequence " , "FrameCount" , "Time" , "Ix" , "Fx" , "Iy" , "Fy" , "Max-Min" , "Min" ,} ;
    char *tform[] = {"U" , "U" , "D" , "U" , "E" , "U" , "E" , "B" , "B"} ;
    char *tunit[] = {"" , "" , "" , "" , "" , "" , "" , "" , ""} ;
    int tfields = 9 ;

    //creating event table in output event file
    fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , tunit , EVENTDATA , &status) ;
    printError (status , "Error in creating the table in output event File" , eventFile) ; //table for valid events

    //creating table for storing parity failed events in output event file
    fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , tunit , "ParityFailedEvents" , &status) ;
    printError (status , "Error in creating the table  parityFailed Events in the output event File" , eventFile) ; //table for invalid events

    // moving to 2nd HDU (EVENT table) in output EVENT file
    fits_movabs_hdu (fout , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU" , eventFile) ;

    //for pointing to 3rd HDU (parity failed events) in output EVENT file
    fits_open_file (&fout_h3 , eventFile , READWRITE , &status) ;
    printError (status , "Error in opening to the eventFile" , eventFile) ; //for writing parity failed events
    fits_movabs_hdu (fout_h3 , 3 , NULL , &status) ;
    printError (status , "Error in moving to the 3rd HDU in the output Event File" , eventFile) ; //for writing parity failed events

    //Creating output IMAGE file
    
    fits_create_file (&fimg , PCImageFile , &status) ;
    printError (status , "Error in creating the PCimage File" , PCImageFile) ;

    //Declaring image parameters for FITS
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = IMG_DIM_DI;
    naxes[1] =IMG_DIM_DI ;

    
    //LOG(INFO)<<naxes[0]<<" "<<naxes[1];exit(1);
    //Creating image HDU in output image file
    fits_create_img (fimg , bitpix , naxis , naxes , &status) ;
    printError (status , "Error in creating the image " , PCImageFile) ;

    //Copying keywords from dataIngestout file to output IMAGE file
    writeUsrkeywordsFrmvect (PCImageFile , key_record) ;
    //copyUserKeywords (fptr_primary , fimg) ;

    float *imageArray = new float[naxes[0] * naxes[1]] ; //buffer to store image
    initArray<float>(imageArray , naxes[0] * naxes[1] , 0.0f) ; //initializing array with zero

    long nevent = MAX_EVENTS_PER_PACKET*nrows ; //336 is the max number of events in one row/packet data
    unsigned char Dmm ; //Max-Min
    unsigned char Min ;

    long eventindex = 0 ;
    int Ix , Fx , Iy , Fy ;
    float fx , fy ;
    unsigned char signbit = 0 ;
    unsigned short word1 , word2 , word3 ;
    long invalideventcount = 0 ;
    bool removeInvalidRecords = FALSE ;
    
    
    //Loop to extract event data from 'centroiddata' buffer and also check for parity of events
    for (long i = 0 ; i < nrows * CENTROID_BYTES ; i = i + 6)
    {
        removeInvalidRecords = FALSE ; //Flag to check for INVALID EVENTS (events which are 0)
        word1 = (((int) centroiddata[i]) << 8) + ((int) centroiddata[i + 1]) ;
        word2 = (((int) centroiddata[i + 2]) << 8) + ((int) centroiddata[i + 3]) ;
        word3 = (((int) centroiddata[i + 4]) << 8) + ((int) centroiddata[i + 5]) ;

        unsigned char p1 = getParity (word1) ;
        unsigned char p2 = getParity (word2) ;
        unsigned char p3 = getParity (word3) ;

        Ix = ((int) word1 >> 7) ;
        Fx = (((int) centroiddata[i + 1] >> 1) & 0x3f) ;
        Iy = ((int) word2 >> 7) ;
        Fy = (((int) centroiddata[i + 3] >> 1) & 0x3f) ;


        if (Ix == 0 && Iy == 0 && Fx == 0.0f && Fy == 0.0f)
        {
            removeInvalidRecords = TRUE ;
        }

        signbit = (Fx >> 5) & 0x1 ;
        Fx = Fx & 0x1f ; //taking 5 bits
        if (signbit == 0)
            fx = (float) Fx / 32.0 ;
        else
            fx = (float) (Fx - 32) / 32.0 ;

        signbit = (Fy >> 5) & 0x1 ;
        Fy = Fy & 0x1f ; //taking 5 bits
        if (signbit == 0)
            fy = (float) Fy / 32.0 ;
        else
            fy = (float) (Fy - 32) / 32.0 ;

        if (parity_flag == 1)
        {

            if (p1 != 0 || p2 != 0 || p3 != 0)
            { //Even parity


                if (! removeInvalidRecords)
                {
                    invalideventcount ++ ;
                    //writing invalid events to 3rd extension of the output event file
                    fits_write_col (fout_h3 , TUSHORT , 1 , invalideventcount , 1 , 1 , &psc[i / CENTROID_BYTES] , &status) ;
                    printError (status , "Error in writing the  column psc to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TUSHORT , 2 , invalideventcount , 1 , 1 , &framecount[i / CENTROID_BYTES] , &status) ;
                    printError (status , "Error in writing the  column framecount to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TDOUBLE , 3 , invalideventcount , 1 , 1 , &time[i / CENTROID_BYTES] , &status) ;
                    printError (status , "Error in writing the  column time to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TUSHORT , 4 , invalideventcount , 1 , 1 , &Ix , &status) ;
                    printError (status , "Error in writing the  column IX to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TFLOAT , 5 , invalideventcount , 1 , 1 , &fx , &status) ;
                    printError (status , "Error in writing the  column fx to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TUSHORT , 6 , invalideventcount , 1 , 1 , &Iy , &status) ;
                    printError (status , "Error in writing the  column Iy to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TFLOAT , 7 , invalideventcount , 1 , 1 , &fy , &status) ;
                    printError (status , "Error in writing the  column fy to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TBYTE , 8 , invalideventcount , 1 , 1 , &Dmm , &status) ;
                    printError (status , "Error in writing the  column Dmm to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TBYTE , 9 , invalideventcount , 1 , 1 , &Min , &status) ;
                    printError (status , "Error in writing the  column Min to the output Event File" , eventFile) ;
                }
                continue ;
            }
        }
        else if (parity_flag == 2)
        {
            if (p1 != 0 || p2 != 0)
            {
                //Even parity
                if (! removeInvalidRecords)
                {
                    invalideventcount ++ ;
                    //writing invalid events to 3rd extension of the output event file
                    fits_write_col (fout_h3 , TUSHORT , 1 , invalideventcount , 1 , 1 , &psc[i / CENTROID_BYTES] , &status) ;
                    printError (status , "Error in writing the  column psc to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TUSHORT , 2 , invalideventcount , 1 , 1 , &framecount[i / CENTROID_BYTES] , &status) ;
                    printError (status , "Error in writing the  column framecount to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TDOUBLE , 3 , invalideventcount , 1 , 1 , &time[i / CENTROID_BYTES] , &status) ;
                    printError (status , "Error in writing the  column time to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TUSHORT , 4 , invalideventcount , 1 , 1 , &Ix , &status) ;
                    printError (status , "Error in writing the  column IX to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TFLOAT , 5 , invalideventcount , 1 , 1 , &fx , &status) ;
                    printError (status , "Error in writing the  column fx to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TUSHORT , 6 , invalideventcount , 1 , 1 , &Iy , &status) ;
                    printError (status , "Error in writing the  column Iy to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TFLOAT , 7 , invalideventcount , 1 , 1 , &fy , &status) ;
                    printError (status , "Error in writing the  column fy to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TBYTE , 8 , invalideventcount , 1 , 1 , &Dmm , &status) ;
                    printError (status , "Error in writing the  column Dmm to the output Event File" , eventFile) ;
                    fits_write_col (fout_h3 , TBYTE , 9 , invalideventcount , 1 , 1 , &Min , &status) ;
                    printError (status , "Error in writing the  column Min to the output Event File" , eventFile) ;

                }
                continue ;
            }
        }

        Dmm = centroiddata[i + 4] >> 1 ;
        Min = ((word3 >> 1) & 0xff) ;

        if (! removeInvalidRecords)
        {
            eventindex ++ ; //Counter for number of VALID events
            fits_write_col (fout , TUSHORT , 1 , eventindex , 1 , 1 , &psc[i / CENTROID_BYTES] , &status) ;
            printError (status , "***Error Writing the column of packet sequence***") ;
            fits_write_col (fout , TUSHORT , 2 , eventindex , 1 , 1 , &framecount[i / CENTROID_BYTES] , &status) ;
            printError (status , "Error in writing the  column framecount to the output Event File" , eventFile) ;
            fits_write_col (fout , TDOUBLE , 3 , eventindex , 1 , 1 , &time[i / CENTROID_BYTES] , &status) ;
            printError (status , "Error in writing the  column time to the output Event File" , eventFile) ;
            fits_write_col (fout , TUSHORT , 4 , eventindex , 1 , 1 , &Ix , &status) ;
            printError (status , "Error in writing the  column IX to the output Event File" , eventFile) ;
            fits_write_col (fout , TFLOAT , 5 , eventindex , 1 , 1 , &fx , &status) ;
            printError (status , "Error in writing the  column fx to the output Event File" , eventFile) ;
            fits_write_col (fout , TUSHORT , 6 , eventindex , 1 , 1 , &Iy , &status) ;
            printError (status , "Error in writing the  column Iy to the output Event File" , eventFile) ;
            fits_write_col (fout , TFLOAT , 7 , eventindex , 1 , 1 , &fy , &status) ;
            printError (status , "Error in writing the  column fy to the output Event File" , eventFile) ;
            fits_write_col (fout , TBYTE , 8 , eventindex , 1 , 1 , &Dmm , &status) ;
            printError (status , "Error in writing the  column Dmm to the output Event File" , eventFile) ;
            fits_write_col (fout , TBYTE , 9 , eventindex , 1 , 1 , &Min , &status) ;
            printError (status , "Error in writing the  column Min to the output Event File" , eventFile) ;

            //setting image values
            imageArray[Iy * naxes[0] + Ix] ++ ;
        }
     //LOG(INFO)<<"The PC image File " <<infoFilename<<" "<<i;
    }
  
    //  LOG(INFO)<<"The PC image File " <<infoFilename;exit(1);
    LOG (INFO) << "Total events   : " << nevent ;
    LOG (INFO) << "Invalid Events : " << invalideventcount ;
    LOG (INFO) << "Valid Events   :  " << nevent - invalideventcount ;

  
    //Reading integration time from dataINgest out file
    double integrationTime ;
    fits_read_key (fptr , TDOUBLE , "INT_TIME" , &integrationTime , NULL , &status) ;
    printError (status , "Error reading the key value of the INT_TIME" , outFile) ;

    fits_movabs_hdu (fout , 1 , NULL , &status) ; //for writing INT_TIME to primary header
    printError (status , "Error in Moving to 1st HDU in outFile" , eventFile) ;

    //Writing integration time to output EVENT file
    fits_write_key (fout , TDOUBLE , "INT_TIME" , &integrationTime , "Integration time (time of one frame)" , &status) ;
    printError (status , "Error writing The key value of the INT_TIME" , eventFile) ;

    fits_close_file (fptr , &status) ;
    printError (status , "Error closing the File" , outFile) ;

    //Writing HISTORY to output EVENT file
    // if (history == YES) writeHistory (eventFile , vhistorystr) ;
    fits_close_file (fout_h3 , &status) ;
    printError (status , "Error closing the File" , eventFile) ;

    /*--------------creating event file completed -----------------*/

   

 
 if (history == YES)  writeHistory ((char*)eventFile , vhistorystr) ;
    
    //writing image to output image file
      writeCommonKeywords (fout , modulename) ;

    fits_close_file (fout , &status) ;
    printError (status , "Error closing the File" , eventFile) ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;

    fits_write_pix (fimg , TFLOAT , fpixel , naxes[0] * naxes[1] , imageArray , &status) ;
    printError (status , "Error writing the PIXELS to the output image file" , PCImageFile) ;
    //copyUserKeywords (fptr_primary , fimg) ;
   
    
    if (history == YES) writeHistory (PCImageFile , vhistorystr) ;

    writeCommonKeywords (fimg , modulename) ; //writes date, origin, checksum and creator to file
     fits_close_file (fimg , &status) ;
    printError (status , "Error Closing the image File" , PCImageFile) ;
    

    if (history == YES) writeHistory (outFile , vhistorystr) ;

    fits_close_file (fptr_primary , &status) ;
    printError (status , "Error closing the File" , outFile) ;

    delete[] psc , centroiddata , time , framecount ;

    LOG (INFO) << "Event data decoding completed" ;

    return (EXIT_SUCCESS) ;
}


int DataIngest::createFrames ()
{
    /**Creating  Signal Frame Directory **/
    string command = (string) "mkdir -p " + (string) IMframespath ;
    system (command.c_str ()) ;

    LOG (INFO) << IMframespath << " directory created" ;

    /**reading the File for Creating Frames from the row wise data**/
    fitsfile *IM_mode_file_ptr , *IM_mode_file_primary_ptr , *IM_mode_file_priPtr ;
    int status = 0 ;
    long nrows_pri_HDU = 0 ;


    fits_open_file (&IM_mode_file_primary_ptr , inFile , READONLY , &status) ;
    printError (status , "Error in opening the output File For IM mode" , outFile) ;
    fits_movnam_hdu (IM_mode_file_primary_ptr , BINARY_TBL , SCIENCEDATA_HDUNAME_FIRST , 0 , &status) ;
    printError (status , "Error moving to DETECTOR SETTING   hdu" , outFile) ;
    fits_get_num_rows (IM_mode_file_primary_ptr , &nrows_pri_HDU , &status) ;
    unsigned short telemetry_info[PIXEL_DARK] ; //for using telemetry for dark frame identification

    fits_open_file (&IM_mode_file_ptr , outFile , READWRITE , &status) ;
    printError (status , "Error opening data ingest output file" , outFile) ;
    fits_movnam_hdu (IM_mode_file_ptr , BINARY_TBL , SCIENCEDATA_HDUNAME , 0 , &status) ;
    printError (status , "Error moving to science  hdu" , outFile) ;

    long nrows = 0 ;
    fits_get_num_rows (IM_mode_file_ptr , &nrows , &status) ;
    printError (status , "Error in getting the number of rows " , outFile) ;

    double integrationTime = 0 ;
    fits_read_key (IM_mode_file_ptr , TDOUBLE , "INT_TIME" , &integrationTime , NULL , &status) ;
    printError (status , "Error reading the key value of the INT_TIME" , outFile) ;

    /**Getting the Size of the Frame**/
//    int xsize = datainfo.getXsize () ;
//    int ysize = datainfo.getYsize () ;

    
    
    int xsize = datainfo.x+1 ;
    int ysize = datainfo.y+1 ;
    LOG (INFO) << "Window size : " << xsize << "X" << ysize ;

    float imageframe[ IMG_DIM_DI * IMG_DIM_DI] ;
    unsigned short framecount=0 , psc ;
    unsigned short pixel[PIXELNO] ;

    int segflag = 1 ;
    char filename[NAMESIZE] ; //file name for the signal frame  file

    fitsfile *fframe ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = IMG_DIM_DI ;
    naxes[1] = IMG_DIM_DI ;
    ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;

    double temptime = 0 ;
    long pixelindex = 0 ;
    long firstelem = 1 , nelement = 1 ;

    int p = datainfo.xoff ;
    int q = datainfo.yoff ;
    int i , j ;

    int pointer = 0 ;
    sort (list.begin () , list.end ()) ;
    unsigned short darkframe_id ;
    int pri_hdu_rowno = 0 ;
   
    /**LOOP for reading the Row wise data and creating the pixel array for each frame based on the frame number and pixels(1008 pixels)**/
    for (long firstrow = 1 ; firstrow <= nrows ; firstrow ++)
    {
        fits_read_col (IM_mode_file_ptr , TDOUBLE , 3 , firstrow , firstelem , nelement , NULL , &temptime , NULL , &status) ;
        printError (status , "Error in reading the time from the  file " , outFile) ;

        for (int j = 0 ; j < IMG_DIM_DI * IMG_DIM_DI ; j ++) imageframe[j] = -9999 ; //loop to initialize each image frame
        darkframe_id = 1 ;
        pixelindex = 0 ;

        //long startlocation = (datainfo.yoff) * xsize + datainfo.xoff ;
        long startlocation = (datainfo.yoff) * IMG_DIM_DI + datainfo.xoff ;
       long First_img_pix_loc=startlocation;
        long temp_location = startlocation ;
      //  LOG(INFO)<<temp_location;exit(1);
        int increase_location = 0 ;
        pri_hdu_rowno ++ ;
        
        do
        {
            //INVALID _PIX_VALUE to be put for packets for which CRC fails .
            for (int i = pointer ; i < list.size () ; i ++)
            {
                if (firstrow == list[i])
                {
                    i = i - 1 ;

                    do
                    {
                        i ++ ;
                        pointer ++ ;
                        for (int k = 0 ; k < PIXELNO ; k ++ , pixelindex ++)
                        {
                            if (pixelindex >= (IMG_DIM_DI*IMG_DIM_DI)) break ; //full pixels are not to be read from last packet of frame
                            imageframe[pixelindex] = INVALID_PIX_VALUE ;
                        }

                    }
                    while (list[i] == list[i + 1] - 1) ;
                }
            }

            fits_read_col (IM_mode_file_ptr , TUSHORT , 1 , firstrow , firstelem , nelement , NULL , &psc , NULL , &status) ;
            printError (status , "Error Reading the column of the psc" , outFile) ;
            fits_read_col (IM_mode_file_ptr , TUSHORT , 2 , firstrow , firstelem , nelement , NULL , &framecount , NULL , &status) ;
            printError (status , "Error Reading the Column of the framecount" , outFile) ;
            fits_read_col (IM_mode_file_ptr , TUSHORT , 4 , firstrow , firstelem , (nelement * PIXELNO) , NULL , (void *) pixel , NULL , &status) ;
            printError (status , "Error Reading the Column of Pixels" , outFile) ;

//            if (framecount < nrows_pri_HDU)
//            {
//                fits_read_col (IM_mode_file_primary_ptr , TUSHORT , 15 , pri_hdu_rowno , firstelem , (nelement * PIXEL_DARK) , NULL , (void *) telemetry_info , NULL , &status) ;
//                printError (status , "Error Reading the Column of Pixels" , outFile) ;
//
//                darkframe_id = (telemetry_info[32] << 8) | telemetry_info[33] ; //check for Dark frame
//            }
//            else
//            {
//                darkframe_id = 9999.0 ;
//            }
           // LOG(INFO)<<pixel[0]<<endl;exit(1);
           
            
            for (int k = 0 ; k < PIXELNO ; k++ , startlocation ++)
            {
                   
                
//                if (((startlocation - temp_location) % xsize) == 0 && startlocation != 0)
                   //  if (((startlocation - temp_location) % xsize) == 0 &&startlocation!=First_img_pix_loc)
                   if (((startlocation - temp_location) %xsize) == 0 && startlocation!=First_img_pix_loc)
                {
                  
                   increase_location ++ ; //for the next row to be filled;
                    startlocation = temp_location = (datainfo.yoff + increase_location) * IMG_DIM_DI + datainfo.xoff ;
                    //increase_location ++ ; //for the next row to be filled;
                }
                if (startlocation >= (IMG_DIM_DI * IMG_DIM_DI)) break ;
                if (increase_location < ysize)
                {
                   // LOG(INFO)<<startlocation<<" "<<increase_location<<" "<<temp_location<<" "<<pixel[k];
                    imageframe[startlocation] = pixel[k] ;
                }
                else
                {
                    break ;
                }

            }
           
            segflag = (int) (psc) >> 14 ;
            firstrow ++ ;
            if (firstrow > nrows) break ;
            if (segflag == 2) break ; //come out of loop after reading last segment in frame
        }
        while (segflag == 0 || segflag == 1) ; //segflag = 0 for continuation , 1 for first, 2 for last segment

        firstrow -- ;

        char *basefilename ;
        /*creating Frames and writing pixel Array(512*512) to it*/
        sprintf (filename , "%s/%s_t%.4f_f%d_di.fits" , IMframespath , filePrefix , (float) temptime , framecount) ;
        fits_create_file (&fframe , filename , &status) ;
        printError (status , "Cant Create the new File" , filename) ;
        fits_create_img (fframe , bitpix , naxis , naxes , &status) ;
        printError (status , "Error creating the Image" , filename) ;
        fits_write_pix (fframe , TFLOAT , fpixel , IMG_DIM_DI*IMG_DIM_DI , imageframe , &status) ;
        printError (status , "Error in writing the Pixels" , filename) ;
        writeUsrkeywordsFrmvect (filename , key_record) ;
        //copyUserKeywords (IM_mode_file_priPtr , fframe) ;

        basefilename = basename (filename) ;
        framelistvector.push_back ((string) basefilename) ;

        updateKeywords (filename , 1 , 3 , TDOUBLE , "FRMTIME" , &temptime ,
                TDOUBLE , "INT_TIME" , &integrationTime ,
                TUSHORT , "FRAMENO" , &framecount) ;


        writeCommonKeywords (fframe , modulename) ; //writes date, origin, checksum and creator to file
        fits_close_file (fframe , &status) ;
        printError (status , "Error in closing the File" , filename) ;
        if (history == YES) writeHistory (filename , vhistorystr) ; //Write HISTORY to output frame file

    }
    list.clear () ;
    if (history == YES) writeHistory (outFile , vhistorystr) ;

    // writeCommonKeywords (IM_mode_file_priPtr,modulename);
    fits_close_file (IM_mode_file_ptr , &status) ;
    printError (status , "Error in closing the File" , filename) ;


    LOG (INFO) << "Frame creation completed" ;
    fits_close_file (IM_mode_file_primary_ptr , &status) ;
    printError (status , "Error in closing the File" , inFile) ;
    return (EXIT_SUCCESS) ;
}

//Function to create Info file which is to be used by all other modules for reading information about input data


int DataIngest::createInfoFile ()
{

    fitsfile *finfo ; //FITS pointer for creating output information file
    int status = 0 ;
  
    fits_create_file (&finfo , infoFilename , &status) ; //Creating output Information file
    printError (status , "***Error creating information file***") ;
   
    char *ttype[] = {"SignalFileList"} ;
    char *tform[] = {"A256"} ;
    int tfields = 1 ;
    //LOG(INFO)<<framelistvector.size ();;exit(1);
    long nrows = framelistvector.size () ; //size of number of frames generated as a output
    fits_create_tbl (finfo , ASCII_TBL , nrows , tfields , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table" , infoFilename) ;
    fits_write_key (finfo , TSTRING , "DARKDIR" , (void*) dark_outputDir , "Directory containing dark frames" , &status) ;
    printError (status , "DARK directory keyword not Found" , infoFilename) ;
    if (datainfo.getModeFlag () == IM) //Writing  keywords of NFILES nd SIGDIR to information file  in case of IM mode
    {
        fits_write_key (finfo , TLONG , "NFILES" , &nrows , "Number of files" , &status) ;
        printError (status , "Error in writing the key value of the NFILES" , infoFilename) ;
        fits_write_key (finfo , TSTRING , "SIGDIR" , basename (IMframespath) , "Directory containing signal frames" , &status) ;
        printError (status , "SIGDIR keyword not Found" , infoFilename) ;

    }
    datainfo.write (finfo) ; //Writing basic keywords  to the output information file
    fits_write_key (finfo , TSTRING , "NAMEPRFX" , filePrefix , "Filename Prefix" , &status) ;
    printError (status , "Error in writing the key value of the NAMEPRFX") ;
    fits_write_key (finfo , TSTRING , "MISSINGF" , basename (missingFrameList) , "Missing frame list file" , &status) ;
    printError (status , "Error  in writing the key value of the MISSINGF") ;
    fits_write_key (finfo , TSTRING , "CRCLOGF" , basename (crcFailedList) , "CRC failed list file" , &status) ;
    printError (status , "Error in Writing the CRCLOGF " , infoFilename) ;
    fits_write_key (finfo , TSTRING , "FRMTIMEF" , basename (timeCorrectionFile) , "frame time and file" , &status) ;
    printError (status , "Error in Writing the FRMTIMEF " , infoFilename) ;

    if (datainfo.getModeFlag () == PC) //Writing  keywords of EVTFILE  nd IMGFILE to information file  in case of PC mode
    {
        fits_write_key (finfo , TSTRING , "EVTFILE" , basename (eventFile) , "event file generated" , &status) ;
        printError (status , "Error in writing the key value of the EVTFILE" , infoFilename) ;
        fits_write_key (finfo , TSTRING , "IMGFILE" , basename (PCImageFile) , "PC mode Image file" , &status) ;
        printError (status , "Error in writing the key value of the infoFile" , infoFilename) ;
    }

    char *names[nrows] ;
    //Writing output frame's name to output information file
    if (datainfo.getModeFlag () == IM)
    {
        long firstrow = 1 , firstelem = 1 ;
        int colnum = 1 ;
        for (long i = 0 ; i < nrows ; i ++)
        {
            const char *temp1 = framelistvector[i].c_str () ;
            names[i] = const_cast<char *> (temp1) ;
        }
        fits_write_col (finfo , TSTRING , colnum , firstrow , firstelem , nrows , (void *) names , &status) ;
        printError (status , "***Error writing filename***" , infoFilename) ;
    }
    writeCommonKeywords (finfo , modulename) ; //writes date, origin, checksum and creator to file

    fits_close_file (finfo , &status) ;
    printError (status , "Error in closing the File" , infoFilename) ;
    //LOG (INFO) << "Successfully closed" << endl ;
   // LOG (INFO) << "HHH" << endl ;
    status = writeUsrkeywordsFrmvect (infoFilename , key_record) ;
    if (status)
    {
        LOG (INFO) << "Information file creation completed" ;
        return(EXIT_FAILURE);
    }
    //LOG (INFO) << "HHH1" << endl ;
    //LOG (INFO) << "Information file creation completed" ;

    return (EXIT_SUCCESS) ;
}

//function to generate level2 gti using input gti file and lbt file


//int DataIngest::genL2gti ()
//{
//    fitsfile *flbt , *fgti ;
//    int status = 0 ;
//
//    fits_open_file (&fgti , inFile_GTI , READONLY , &status) ;
//    printError (status , "") ;
//
//    fits_movabs_hdu (fgti , 2 , NULL , &status) ;
//    printError (status , "") ;
//    LOG (WARNING) << endl << "$$$$   Code for generating level 2 gti is yet to be implemented  $$$$$" << endl ;
//
//
//    fits_close_file (fgti , &status) ;
//
//    return (EXIT_SUCCESS) ;
//}


int DataIngest::setFilePaths ()
{
    LOG (INFO) << "Setting file paths" ;

    //getting file prefix for output files based on input file
    char *infilepath = strdup (inFile) ;
    char *temp = basename (infilepath) ; //retrieve onlyfilename
    char *ext = strchr (temp , '.') ; //extract name before extension
    if (ext != NULL) strtok (temp , ".") ;

    //changing level1 in fileprefix to level2
    string strtemp (temp) ;
    int pos = strtemp.find ("level1") ;
    if (pos > 0 && pos < strtemp.length ())
    {
        strtemp.replace (pos , string ("level1").length () , "l2") ;
    }

    strcpy (filePrefix , strtemp.c_str ()) ; //storing filename prefix
    sprintf (outpath , "%s/%s/" , outdir , modulename) ; //creating name for output directory path
    sprintf (IMframespath , "%s/SignalFrames" , outpath) ; //creating name for output path for IM frames
    sprintf (Darkframes_path , "%s/Darkframes/" , "/tmp/tempo") ;
    //sprintf (darkFramesoutPath , "%s/%s" , outpath , basename (darkFramesPath)) ; //creating name for output paths for dark frames
    sprintf (darkFramesoutPath , "%s/%s/" , outpath , dark_outputDir) ; //creating name for output paths for dark frames
    /* check for existence of output directory ,if output directory already exist and clobber is yes than delete the output directory, & if output directory already exist and
       clobber is No than  No need to overwrite and exit.  */
    eventFile = new char[FLEN_FILENAME];
    infoFilename= new char[FLEN_FILENAME];
    PCImageFile= new char[PIL_LINESIZE];
    if (createOutputDirectory (clobber , outpath)) return (EXIT_FAILURE) ;

    //cerr << endl << "Output Files to be generated are-------" ;
    LOG (INFO) << "Output Files to be generated are-------" ;

    //creating output filenames
    sprintf (outFile , "%s%s_di.dataIngestOut" , outpath , filePrefix) ;
    LOG (INFO) << outFile ;
    sprintf (crcFailedList , "%s%s_di.crcFailedList" , outpath , filePrefix) ;
    LOG (INFO) << crcFailedList ;
    sprintf (missingFrameList , "%s%s_di.missingFrameList" , outpath , filePrefix) ;
    LOG (INFO) << missingFrameList ;
    sprintf (timeCorrectionFile , "%s%s_di.timeCorrection" , outpath , filePrefix) ;
    LOG (INFO) << timeCorrectionFile ;
    sprintf (eventFile , "%s%s_di.events" , outpath , filePrefix) ;
    LOG (INFO) << eventFile ;
    sprintf (PCImageFile , "%s/%s_di.image" , outpath , filePrefix) ;
    LOG (INFO) << PCImageFile ;
    sprintf (infoFilename , "%s/%s_di.info" , outpath , filePrefix) ;
    LOG (INFO) << infoFilename;

//    sprintf (dark_infoFile , "%s%s_di.darkinfo" , outpath , filePrefix) ;
//    LOG (INFO) << dark_infoFile ;

    return (EXIT_SUCCESS) ;
}


int DataIngest::getHistory (vector<string> &vhistory)
{
    char *user = getlogin () ;
    int cnt = 0 ;
    char validgtiflag_str[FLEN_FILENAME] ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " input Level1 file = " + (string) inFile) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Module Output directory = " + (string) outdir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Input TCT file = " + (string) inFile_TCT) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Input GTI file = " + (string) inFile_GTI) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Input LBT file = " + (string) inFile_LBT) ;
    if (datainfo.getModeFlag () == IM)
    {
        vhistory.push_back ((string) getSerialNo (cnt) + " Signal Frame Directory = " + (string) "SignalFrames") ;
    }
    if (gti_flag == 0)
    {
        vhistory.push_back ((string) getSerialNo (cnt) + " GTI filtering  = " + (string) "NO") ;
    }
    else if (gti_flag == 1)
    {
        vhistory.push_back ((string) getSerialNo (cnt) + " GTI filtering  = " + (string) "YES") ;
        sprintf (validgtiflag_str , "%d" , valid_gtiflag) ;
        vhistory.push_back ((string) getSerialNo (cnt) + " Valid GTI flag  = " + (string) validgtiflag_str) ;
        if (all_Or_custom)
        {
            vhistory.push_back ((string) getSerialNo (cnt) + " Parameters  checked(all/custom)= " + (string) "All") ;
        }
        else if (all_Or_custom == 0)
        {
            vhistory.push_back ((string) getSerialNo (cnt) + " Parameters checked(all/custom)  = " + (string) "Custom") ;
        }

    }
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


float DataIngest::findMedianValue (vector<float> &input , int noOfValues)
{
    //cout<<"INSIDE";
    float temp ;
    int index , j ;
    sort (input.begin () , input.end ()) ;
    // Sort the input data
    //        for(index=0;index<noOfValues;index++)
    //                for(j=index+1;j<noOfValues;j++)
    //                {
    //                        if(input[index]>input[j])
    //                        {
    //                                temp=input[j];
    //                                input[j]=input[index];
    //                                input[index]=temp;
    //                        }
    //                }

    //    If the no of values are even number, take the average of middle and middle-1 value.
    if (noOfValues % 2 == 0)
        return (input[noOfValues / 2] + input[(noOfValues / 2) - 1]) / 2.0 ;
    else
        return input[noOfValues / 2] ;
}


int DataIngest::readDarkFrames ()
{
    fitsfile *fout ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    int bitpix = FLOAT_IMG ;
    long naxes[2] ;
    bool darkstart_found = FALSE ;
    bool darkend_found = FALSE ;
    naxes[0] = naxes[1] = IMG_DIM_DI ;
    int naxis = 2 ;
    char filename[FLEN_FILENAME] ; //dark frame input file
    char outfilename[FLEN_FILENAME] ; //output dark frame file
    double prod_StartTime , prod_sfractTime ; //start and end  fractional time of the product
    double prod_EndTime , prod_efractTime ; //start and end integer time of the product
    double dark_framestarttime , dark_frameendtime ; //dark frame time
    int status = 0 ;

    vector<string> darkfilenames ; //to store filenames of dark frame files in level-1 data

    status = getFileNamesfrmDir (darkFramesPath , ".fits" , darkfilenames) ; //getting files  from darkframedPath directory those are ends with ".fits" and store their names to darkfilenames
    if (status)
    {
        LOG (ERROR) << "Error in getting  dark files from directory" << endl ;
        return (EXIT_FAILURE) ;
    }
    LOG (INFO) << "Total  Dark Frames - " << darkfilenames.size () ;

    //if only one dark at a input than considering it as a median of the dark.
    if (darkfilenames.size () == 0)
    {
        char full_caldbpath_forDarkfiles[FLEN_FILENAME] ; //path of dark frames from caldb directory.
        vector<string> darkfile_frm_caldb ;
        char temp_Dark_dir[FLEN_FILENAME] ;
        LOG (INFO) << "No dark frames found at input " ;
        LOG (INFO) << "Taking Dark frames from CALDB " << endl ;
        LOG (INFO) << "CALDB directory -> " << caldbindir << endl ;
        sprintf (temp_Dark_dir , "%s/%s" , caldbindir , "DARK") ;

        status = getFileNamesfrmDir (temp_Dark_dir , ".fits" , darkfile_frm_caldb) ;
        if (darkfile_frm_caldb.size () == 0)
        {
            LOG (ERROR) << "Darks are also not available in CALDB directory......!!!!! " << endl ;
            return (EXIT_FAILURE) ;
        }
        string command = (string) "mkdir -p " + (string) darkFramesoutPath ; //
        system (command.c_str ()) ;

        //copyo darkbegin from caldb to output DataIngest directory.
        sprintf (full_caldbpath_forDarkfiles , "%s/%s" , temp_Dark_dir , (char*) darkfile_frm_caldb[0].c_str ()) ;
        sprintf (outfilename , "%s/%s_darkbegin.fits" , darkFramesoutPath , filePrefix) ;
        command = "cp " + (string) full_caldbpath_forDarkfiles + " " + outfilename ;
        LOG (INFO) << "Executing command " << command << endl ;
        system (command.c_str ()) ;

        //copy darkend from caldb  to output DataIngest directory.
        sprintf (full_caldbpath_forDarkfiles , "%s/%s" , temp_Dark_dir , (char*) darkfile_frm_caldb[1].c_str ()) ;
        sprintf (outfilename , "%s/%s_darkend.fits" , darkFramesoutPath , filePrefix) ;
        command = "cp " + (string) full_caldbpath_forDarkfiles + " " + outfilename ;
        LOG (INFO) << "Executing command " << command << endl ;
        system (command.c_str ()) ;


        return (EXIT_SUCCESS) ;
    }
//        if(darkfilenames.size ()==1)//in case of only one dark frame 
//        {
//            LOG(ERROR)<<"Total number of dark frame is "<<darkfilenames.size ()<<" ,  "<<"No extra procees will be done on dark";
//            LOG(ERROR)<<"Now taking median of  dark";
//             sprintf (outfilename , "%s/%s_dark.fits" , darkFramesoutPath , filePrefix) ;
//             string cmd="cp "+darkfilenames[0]+" "+outfilename;
//             LOG(INFO)<<"Executing "<<cmd<<" command";
//             system (cmd.c_str ());
//            return(EXIT_SUCCESS);               
//        }


    //  vector of structure dark frame_info
    vector<darkframe_info> darkframe_track , darkframe_track_temp ;
    darkframe_info darkframe_infoobj ;

    double start_i , start_f , stop_i , stop_f ;
char * filtr_dark= new char[10];
        char *obs_mode_dark = new char[10];
char *high_vol_dark_status= new char[10];
    //  loop for reading name and time of darkframes
 
    for (int i = 0 ; i < darkfilenames.size () ; i ++)
    {
        sprintf (filename , "%s/%s" , darkFramesPath , darkfilenames[i].c_str ()) ;
         LOG(INFO)<<filename;
        readKeywords (filename , 1 , 7, TDOUBLE , "TSTARTI" , &start_i ,
                TDOUBLE , "TSTARTF" , &start_f ,
                TDOUBLE , "TSTOPI" , &stop_i ,
                TDOUBLE , "TSTOPF" , &stop_f,
                TSTRING,"FILTER",filtr_dark,
                TSTRING ,"OBS_MODE",obs_mode_dark,
                TSTRING,"HVSTATUS"  ,high_vol_dark_status) ;
                  //LOG(INFO)<<filename;
                  LOG(INFO)<<filtr_dark<<" "<<obs_mode_dark<<" "<<" "<<high_vol_dark_status;
        if(strcmp (filtr_dark,"F0")!=0 ||strcmp (obs_mode_dark,"IM")!=0|| strcmp(high_vol_dark_status,"OFF")!=0)
        {
            LOG(WARNING)<<"***"<<filename<<" is not  a dark file ,FILE   INFO: FILTER -> "<<filtr_dark<<" Mode -> "<<obs_mode_dark<<" HVSTATUS ->"<<high_vol_dark_status<<" ***";            
        }

        //      calculating the time for Dark frame
        dark_framestarttime = start_i + start_f ;
        dark_frameendtime = stop_i + stop_f ;
        darkframe_infoobj.darkframe_starttime = dark_framestarttime / 1000 ; //storing darkframe start  time (in seconds)
        darkframe_infoobj.darkframe_endtime = dark_frameendtime / 1000 ; //storing darkframe end time (in converting to seconds)
        darkframe_infoobj.darkfrm_name = filename ; //storing darkframe name
        darkframe_track.push_back (darkframe_infoobj) ;

    }
    darkframe_track_temp = darkframe_track ;
   
    //   reading  product start time and  product end  time  from the level1 file
    datainfo.getTime (&prod_StartTime , &prod_sfractTime , &prod_EndTime , &prod_efractTime) ;

    prod_StartTime = prod_StartTime + prod_sfractTime ;
    prod_EndTime = prod_EndTime + prod_efractTime ;

    //    sort dark frame time
    sort (darkframe_track.begin () , darkframe_track.end () , compare_endtime) ;

    string darkframe_startPath = "" , darkframe_endPath = "" ;
    double darkframe_startTime , darkframe_endTime ;

    bool darkFlag = FALSE ;

    //  compare product start time with dark frame's time  to identify dark begin

    for (int i = 0 ; i < darkframe_track.size () ; i ++)
    {

        if (prod_StartTime > darkframe_track[i].darkframe_endtime)
        {
            darkframe_startPath = (darkframe_track[i].darkfrm_name.c_str ()) ;
          //  cout << darkframe_track[i].darkfrm_name.c_str () << endl ;
            darkframe_startTime = darkframe_track[i].darkframe_starttime ;
            darkFlag = TRUE ;
            break ;
        }

    }

    //    if (! darkFlag)
    //    {
    //        LOG (INFO) << "No Start dark frame found at input"  ;
    //        return(EXIT_FAILURE);
    //
    //    }
    //darkFlag = FALSE ;


    sort (darkframe_track.begin () , darkframe_track.end () , compare_starttime) ;
    //  compare product end time with dark frame's time  to identify dark end
    // cout<<prod_EndTime<<" "<<darkframe_track[0].darkframe_starttime;exit(1);

    for (int i = 0 ; i < darkframe_track.size () ; i ++)
    {
        if (prod_EndTime < darkframe_track[i].darkframe_starttime)
        {
            darkframe_endPath = (darkframe_track[i].darkfrm_name.c_str ()) ;
            darkframe_endTime = (darkframe_track[i + 1].darkframe_starttime) ;
            darkFlag = TRUE ;
            break ;
        }
    }
    if (! darkFlag)
    {
        LOG (INFO) << "NO  start or end dark frame found" ;
        return (EXIT_FAILURE) ;
    }

    //LOG (INFO) << "Dark Frame startPath  -" << darkframe_startPath ;
    //LOG (INFO) << "Dark Frame endPath    -" << darkframe_endPath  ;

    float *median_pixelArray_darkstart = new float [IMG_DIM_DI * IMG_DIM_DI] ;
    float *median_pixelArray_darkend = new float [IMG_DIM_DI * IMG_DIM_DI] ;

    //float *temp_start1= new float [IMG_DIM_DI*IMG_DIM_DI];
    //float *temp_start0= new float [IMG_DIM_DI*IMG_DIM_DI];

    vector<float> imageframe_begin , imageframe_end ;
    float *Diff_First2_image_begin = new float [IMG_DIM_DI * IMG_DIM_DI] ;
    float *Diff_First2_image_end = new float [IMG_DIM_DI * IMG_DIM_DI] ;
    vector <float> Diff_First2_image_begin_vect ;
    vector<float>Diff_First2_image_end_vect ;
    double darkBegin_time , darkEnd_time ;
    int total_frames_darkbegin , total_frames_darkend ;
    int full_frames_darkstart ;
   
    if (strcmp (darkframe_startPath.c_str () , "") != 0)
    {
       
        total_frames_darkbegin = full_frames_darkstart = decodeDarkFrames ((char*) darkframe_startPath.c_str () , imageframe_begin , darkBegin_time) ; //decoding dark frame begin.
        darkstart_found = TRUE ;
    }
    vector<float> pixel_array_darkbegin , pixel_array_darkend ;

    if (strcmp (darkframe_endPath.c_str () , "") != 0)
    {
        
        total_frames_darkend = decodeDarkFrames ((char*) darkframe_endPath.c_str () , imageframe_end , darkEnd_time) ; //decoding dark frame end.
        darkend_found = TRUE ;
    }

    float *sum_begin = new float[IMG_DIM_DI * IMG_DIM_DI] ;
    float *sum_end = new float[IMG_DIM_DI * IMG_DIM_DI] ;
   
    double thr_sd_begin ;


    if (strcmp (darkframe_startPath.c_str () , "") != 0)
    {
        LOG (INFO) << "Dark begin Computation started" << endl ;
        for (int filenum = 1 ; filenum < total_frames_darkbegin ; filenum = filenum + 10)
        {

            for (int i = 0 ; i < IMG_DIM_DI * IMG_DIM_DI ; i ++)
            {
                Diff_First2_image_begin[i] = imageframe_begin[(filenum) * IMG_DIM_DI * IMG_DIM_DI + i] - imageframe_begin[(filenum - 1) * IMG_DIM_DI * IMG_DIM_DI + i] ;
            }
            //taking sd of diff image
            double sd_begin = getSD (Diff_First2_image_begin , IMG_DIM_DI * IMG_DIM_DI) ;


            thr_sd_begin = SIGMA_FAC_DARK*sd_begin ;

            Diff_First2_image_begin_vect.clear () ;

            //removing noice i.e above the n*SD
            for (int i = 0 ; i < IMG_DIM_DI * IMG_DIM_DI ; i ++)
            {
                if (Diff_First2_image_begin[i] > thr_sd_begin && Diff_First2_image_begin[i]<- thr_sd_begin)
                {
                    Diff_First2_image_begin[i] = 0.0f ;
                }
                else
                {
                    Diff_First2_image_begin_vect.push_back (Diff_First2_image_begin[i]) ;
                }
            }

            //again calculating SD
            sd_begin = getSD (Diff_First2_image_begin_vect.data () , Diff_First2_image_begin_vect.size ()) ; //newer values

            thr_sd_begin = SIGMA_FAC_DARK*sd_begin ;

            for (int i = 0 ; i < IMG_DIM_DI * IMG_DIM_DI ; i ++) sum_begin[i] = 0.0f ;

            for (int pix = 0 ; pix < IMG_DIM_DI * IMG_DIM_DI ; pix ++)
            {
                for (int j = filenum ; j < filenum + 10 ; j ++)
                {
                    sum_begin[pix] = sum_begin[pix] + imageframe_begin[j * IMG_DIM_DI * IMG_DIM_DI + pix] ;
                }
                sum_begin[pix] = sum_begin[pix] / 10 ;
                if (sum_begin[pix] > thr_sd_begin && sum_begin[pix]< - thr_sd_begin)
                {
                    sum_begin[pix] = 0.0f ;
                }
                pixel_array_darkbegin.push_back (sum_begin[pix]) ;

            }

        } //finish of for loop
        LOG (INFO) << "Dark begin Computation completed" << endl ;

    }

    LOG (INFO) << "Dark end Computation started...." << endl ;
    //for dark end
    if (strcmp (darkframe_endPath.c_str () , "") != 0)
    {
        for (int filenum = 1 ; filenum < total_frames_darkend ; filenum = filenum + 10)
        {

            for (int i = 0 ; i < IMG_DIM_DI * IMG_DIM_DI ; i ++)
            {
                Diff_First2_image_end[i] = imageframe_end[(filenum) * IMG_DIM_DI * IMG_DIM_DI + i] - imageframe_end[(filenum - 1) * IMG_DIM_DI * IMG_DIM_DI + i] ;
            }
            //taking sd of diff image
            double sd_begin = getSD (Diff_First2_image_end , IMG_DIM_DI * IMG_DIM_DI) ;
            //   double sd_end = getSD (Diff_First2_image_end,IMG_DIM_DI*IMG_DIM_DI);

            thr_sd_begin = SIGMA_FAC_DARK*sd_begin ;
            //  double thr_sd_end=3*sd_end; 
            //          fits_create_file (&fout , "Diffend.fits" , &status) ;
            //    printError (status , "Cant Create the new File" , outfilename) ;
            //    fits_create_img (fout , bitpix , naxis , naxes , &status) ;
            //    printError (status , "Error creating the image for the output File" , outfilename) ;
            //    fits_write_pix (fout , TFLOAT , fpixel , IMG_DIM_DI*IMG_DIM_DI , (void*)Diff_First2_image_end , &status) ;
            //    printError (status , "Error in writing the pixels to output Accumulated file" , outfilename) ;
            //    fits_close_file (fout , &status) ;
            //    printError (status , "Error in closing the output File" , outfilename) ;
            Diff_First2_image_end_vect.clear () ;

            //removing noice i.e above the n*SD
            for (int i = 0 ; i < IMG_DIM_DI * IMG_DIM_DI ; i ++)
            {
                if (Diff_First2_image_end[i] > thr_sd_begin && Diff_First2_image_end[i]<- 1 * thr_sd_begin)
                {

                    Diff_First2_image_end[i] = 0.0f ;
                }
                else
                {
                    Diff_First2_image_end_vect.push_back (Diff_First2_image_end[i]) ;
                }
            }
            //         fits_create_file (&fout , "Darkend_DiffImage_byremoving_SD.fits" , &status) ;
            //    printError (status , "Cant Create the new File" , outfilename) ;
            //    fits_create_img (fout , bitpix , naxis , naxes , &status) ;
            //    printError (status , "Error creating the image for the output File" , outfilename) ;
            //    fits_write_pix (fout , TFLOAT , fpixel , IMG_DIM_DI*IMG_DIM_DI , (void*)Diff_First2_image_end_vect.data () , &status) ;
            //    printError (status , "Error in writing the pixels to output Accumulated file" , outfilename) ;
            //    fits_close_file (fout , &status) ;
            //    printError (status , "Error in closing the output File" , outfilename) ;
            //again calculating SD
            sd_begin = getSD (Diff_First2_image_end_vect.data () , Diff_First2_image_end_vect.size ()) ; //newer values

            thr_sd_begin = SIGMA_FAC_DARK*sd_begin ;
            for (int i = 0 ; i < IMG_DIM_DI * IMG_DIM_DI ; i ++) sum_end[i] = 0.0f ;

            for (int pix = 0 ; pix < IMG_DIM_DI * IMG_DIM_DI ; pix ++)
            {
                for (int j = filenum ; j < filenum + 10 ; j ++)
                {
                    sum_end[pix] = sum_end[pix] + imageframe_end[j * IMG_DIM_DI * IMG_DIM_DI + pix] ;
                }
                sum_end[pix] = sum_end[pix] / 10 ;
                if (sum_end[pix] > thr_sd_begin && sum_end[pix]< - thr_sd_begin)
                {
                    sum_end[pix] = 0.0f ;
                }
                pixel_array_darkend.push_back (sum_end[pix]) ;

            }


        } //end of for loop
        LOG (INFO) << "Dark end Computation completed...." << endl ;
    }

    int cnt_pix = 0 ;
    float temp = 0 ;
    vector<float> median_begin ;
    if (strcmp (darkframe_startPath.c_str () , "") != 0)
    {
        do
        {
            median_begin.clear () ;
            for (int i = 0 ; i < total_frames_darkbegin / 10 ; i ++)
            {
                //pixel_array_darkend[i]=imageframe_end[i*nele+cnt_pix];               
                median_begin.push_back (pixel_array_darkbegin[i * IMG_DIM_DI * IMG_DIM_DI + cnt_pix]) ;
            }
            temp = findMedianValue (median_begin , median_begin.size ()) ;
            median_pixelArray_darkstart[cnt_pix] = temp ;
            cnt_pix ++ ;
        }
        while (cnt_pix < IMG_DIM_DI * IMG_DIM_DI) ;
    }


    cnt_pix = 0 ;
    temp = 0 ;
    vector<float> median_end ;
    if (strcmp (darkframe_endPath.c_str () , "") != 0)
    {
        do
        {
            median_end.clear () ;
            for (int i = 0 ; i < total_frames_darkend / 10 ; i ++)
            {
                //pixel_array_darkend[i]=imageframe_end[i*nele+cnt_pix];               
                median_end.push_back (pixel_array_darkend[i * IMG_DIM_DI * IMG_DIM_DI + cnt_pix]) ;
            }
            temp = findMedianValue (median_end , median_end.size ()) ;
            median_pixelArray_darkend[cnt_pix] = temp ;
            cnt_pix ++ ;
        }
        while (cnt_pix < IMG_DIM_DI * IMG_DIM_DI) ;

    }



    // making dark frame's directory at output path

    string command = (string) "mkdir -p " + (string) darkFramesoutPath ; //
    system (command.c_str ()) ;


    //writing darkbegin and darkend frame to output directory.
    if (strcmp (darkframe_startPath.c_str () , "") != 0)
    {
        sprintf (outfilename , "%s/%s_darkbegin.fits" , darkFramesoutPath , filePrefix) ;
        fits_create_file (&fout , outfilename , &status) ;
        printError (status , "Cant Create the new File" , outfilename) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error creating the image for the output File" , outfilename) ;
        fits_write_pix (fout , TFLOAT , fpixel , IMG_DIM_DI*IMG_DIM_DI , (void*) median_pixelArray_darkstart , &status) ;
        printError (status , "Error in writing the pixels to output Accumulated file" , outfilename) ;

        fits_update_key (fout , TDOUBLE , "FRMTIME" , &darkBegin_time , NULL , &status) ; //median time writing
        printError (status , "Error in updating the key value of the FRMTIME" , outfilename) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output File" , outfilename) ;
    }
    if (strcmp (darkframe_endPath.c_str () , "") != 0)
    {
        sprintf (outfilename , "%s/%s_darkend.fits" , darkFramesoutPath , filePrefix) ;
        fits_create_file (&fout , outfilename , &status) ;
        printError (status , "Cant Create the new File" , outfilename) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error creating the image for the output File" , outfilename) ;
        fits_write_pix (fout , TFLOAT , fpixel , IMG_DIM_DI*IMG_DIM_DI , (void*) median_pixelArray_darkend , &status) ;
        printError (status , "Error in writing the pixels to output Accumulated file" , outfilename) ;
        fits_update_key (fout , TDOUBLE , "FRMTIME" , &darkEnd_time , NULL , &status) ; //median time writing
        printError (status , "Error in updating the key value of the FRMTIME" , outfilename) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output File" , outfilename) ;
    }

    if (darkstart_found == FALSE)
    {
        char outdark_begin[FLEN_FILENAME] ;
        sprintf (outdark_begin , "%s/%s_darkbegin.fits" , darkFramesoutPath , filePrefix) ;
        command = "cp " + (string) outfilename + "  " + (string) outdark_begin ;
        system (command.c_str ()) ;
        cout << "inside " << endl ;
        double temp_time = darkEnd_time - 10 ;
        fitsfile *fptr ;
        fits_open_file (&fptr , outdark_begin , READWRITE , &status) ;
        printError (status , "Error opening data ingest output file" , filename) ;
        fits_update_key (fptr , TDOUBLE , "FRMTIME" , &temp_time , NULL , &status) ; //median time writing
        printError (status , "Error in updating the key value of the FRMTIME" , outfilename) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing file" , filename) ;
    }
    else if (darkend_found == FALSE)
    {
        char outdark_end[FLEN_FILENAME] ;
        sprintf (outdark_end , "%s/%s_darkend.fits" , darkFramesoutPath , filePrefix) ;
        command = "cp " + (string) outfilename + "  " + (string) outdark_end ;
        system (command.c_str ()) ;
         double temp_time = darkBegin_time + 10 ;
        fitsfile *fptr ;
        fits_open_file (&fptr , outdark_end , READWRITE , &status) ;
        printError (status , "Error opening data ingest output file" , filename) ;
        fits_update_key (fptr , TDOUBLE , "FRMTIME" , &temp_time , NULL , &status) ; //median time writing
        printError (status , "Error in updating the key value of the FRMTIME" , outfilename) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing file" , filename) ;
    }
    //LOG(INFO)<<"Dark Execution Time : "<<et-st<<" seconds";  
    return (EXIT_SUCCESS) ;
}

//int DataIngest::readDarkFrames ()
//{
//    char filename[FLEN_FILENAME] ; //dark frame input file
//    char outfilename[FLEN_FILENAME] ; //output dark frame file
//    double prod_StartTime , prod_sfractTime ; //start and end  fractional time of the product
//    double prod_EndTime , prod_efractTime ; //start and end integer time of the product
//    double dark_framestarttime ,dark_frameendtime; //dark frame time
//    int status = 0 ;
//    vector<string> darkfilenames ; //to store filenames of dark frame files in level-1 data
//
//
//    getFileNamesfrmDir (darkFramesPath , ".fits" , darkfilenames) ; //getting files  from darkframedPath directory those are ends with ".fits" and store their names to darkfilenames
//
//    LOG (INFO) << "Total  Dark Frames - " << darkfilenames.size () ;
//
//    //vector of structure dark frame_info
//    vector<darkframe_info> darkframe_track ;
//    darkframe_info darkframe_infoobj ;
//
//    double start_i , start_f , stop_i , stop_f ;
//
//    //loop for reading name and time of darkframes
//    for (int i = 0 ; i < darkfilenames.size () ; i ++)
//    {
//        sprintf (filename , "%s/%s" , darkFramesPath , darkfilenames[i].c_str ()) ;
//        readKeywords (filename , 1 , 4 , TDOUBLE , "TSTARTI" , &start_i ,
//                TDOUBLE , "TSTARTF" , &start_f ,
//                TDOUBLE , "TSTOPI" , &stop_i ,
//                TDOUBLE , "TSTOPF" , &stop_f) ;
//
//        //calculating the time for Dark frame
//        dark_framestarttime = start_i + start_f ;
//        dark_frameendtime=stop_i + stop_f;
//        darkframe_infoobj.darkframe_starttime = dark_framestarttime ; //storing darkframe start  time
//        darkframe_infoobj.darkframe_endtime=dark_frameendtime;
//        darkframe_infoobj.darkfrm_name = filename ; //storing darkframe name
//        darkframe_track.push_back (darkframe_infoobj) ;
//
//
//    }
//  
//
//
//    //reading  product start time and  product end  time  from the level1 file
//    datainfo.getTime (&prod_StartTime , &prod_sfractTime , &prod_EndTime , &prod_efractTime) ;
//
//    prod_StartTime = prod_StartTime + prod_sfractTime ;
//    prod_EndTime = prod_EndTime + prod_efractTime ;
//
//    //sort dark frame time
//    sort (darkframe_track.begin () , darkframe_track.end () , compare_starttime) ;
//
//    string darkframe_startPath , darkframe_endPath ;
//    double darkframe_startTime , darkframe_endTime ;
//
//    bool darkFlag = FALSE ;
//    //compare product start time with dark frame's time  to identify dark begin
//    for (int i = 0 ; i < darkframe_track.size () ; i ++)
//    {
//        if (prod_StartTime > darkframe_track[i].darkframe_starttime && prod_StartTime < darkframe_track[i + 1].darkframe_starttime)
//        {
//            darkframe_startPath = (darkframe_track[i].darkfrm_name.c_str ()) ;
//            darkframe_startTime = darkframe_track[i].darkframe_starttime ;
//            darkFlag = TRUE ;
//            break ;
//        }
//
//    }
//    if (! darkFlag)
//    {
//        LOG (INFO) << "NO Start dark frame found at input"  ;
//
//    }
//    darkFlag = FALSE ;
//    //compare product end time with dark frame's time  to identify dark end
//    for (int i = 0 ; i < darkframe_track.size () ; i ++)
//    {
//        if (prod_EndTime > darkframe_track[i].darkframe_starttime && prod_EndTime < darkframe_track[i + 1].darkframe_starttime)
//        {
//            darkframe_endPath = (darkframe_track[i + 1].darkfrm_name.c_str ()) ;
//            darkframe_endTime = (darkframe_track[i + 1].darkframe_starttime) ;
//            darkFlag = TRUE ;
//            break ;
//        }
//    }
//    if (! darkFlag)
//    {
//        LOG (INFO) << "NO End dark  frame found at input"  ;
//    }
//
//    LOG (INFO) << "Dark Frame startPath  -" << darkframe_startPath ;
//    LOG (INFO) << "Dark Frame endPath    -" << darkframe_endPath  ;
//
//
//    //making dark frame's directory at output path
//    string command = (string) "mkdir -p " + (string) darkFramesoutPath ; //
//    system (command.c_str ()) ;
//
//    fitsfile *fptr , *fout ;
//
//    //open dark begin from input directory  and copy to output directory
//    sprintf (outfilename , "%s/%s_darkbegin.fits" , darkFramesoutPath , filePrefix) ;
//    fits_create_file (&fout , outfilename , &status) ;
//    printError (status , "Cant Create the new File" , outfilename) ;
//    fits_open_file (&fptr , darkframe_startPath.c_str () , READONLY , &status) ;
//    printError (status , "Error in opening the output File For IM mode" , outFile) ;
//    fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
//    printError (status , "Error in coping the input event file to the output event file" , outfilename) ;
//
//    if (history == YES) writeHistory (outfilename , vhistorystr) ;
//    writeCommonKeywords (fout , modulename) ;
//    fits_close_file (fptr , &status) ;
//    printError (status , "Error in closing the input File" , (char*)darkframe_startPath.c_str ()) ;
//    fits_close_file (fout , &status) ;
//    printError (status , "Error in closing the input File" , outfilename) ;
//    //open dark end  from input directory  and copy to output directory
//    sprintf (outfilename , "%s/%s_darkend.fits" , darkFramesoutPath , filePrefix) ;
//    fits_create_file (&fout , outfilename , &status) ;
//    printError (status , "Cant Create the new File" , outfilename) ;
//    fits_open_file (&fptr , darkframe_endPath.c_str () , READONLY , &status) ;
//    printError (status , "Error in opening the output File For IM mode" , outFile) ;
//    fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
//    printError (status , "Error in coping the input event file to the output event file" , outfilename) ;
//
//    if (history == YES) writeHistory (outfilename , vhistorystr) ;
//
//    writeCommonKeywords (fout , modulename) ;
//    fits_close_file (fptr , &status) ;
//    printError (status , "Error in closing the input File" , (char*)darkframe_endPath.c_str ()) ;
//
//    fits_close_file (fout , &status) ;
//    printError (status , "Error in closing the output File" , outfilename) ;
//
//    return (EXIT_SUCCESS) ;
//}


bool compare_endtime (struct darkframe_info info1 , struct darkframe_info info2)
{
    return (info1.darkframe_endtime > info2.darkframe_endtime) ;
    //   return (info1.darkframe_starttime <info2.darkframe_starttime) ;

}


bool compare_starttime (struct darkframe_info info1 , struct darkframe_info info2)
{
    return (info1.darkframe_starttime < info2.darkframe_starttime) ;
    //   return (info1.darkframe_starttime <info2.darkframe_starttime) ;

}


int DataIngest::doGTIfiltering (vector<long> &rowno)
{
    
  int status =readGTIparameters();
  if(status)
  {
      LOG(ERROR)<<"Error in reading GTI parameters";
      return(EXIT_FAILURE);
  }
    const int N = 72 ;
    const int NBYTES = 9 ;
    unsigned char final_ByteArray[NBYTES] ;
    int p = 0 ;
    GTIParams P ;

    vector<GTIParams> paramvector ;

    rowno.clear () ;

    ifstream fin (inGtifile) ;

    cout << "File opened " << inGtifile << endl ;
    int i ;
    string temp_line ;

    fin >> temp_line >> temp_line >> temp_line >> temp_line>>temp_line>>temp_line ;
    //cout<<temp_line;exit(1);
    //reading the parameter name,bit positions and allowed values.
    while (! fin.eof ())
    {
        //fin >> num ;
        fin >> P.name >> P.bitpos >> P.min>>P.max>>P.incl_Or_Excl ;
        // LOG(INFO)<< "outputs " << P.name << P.bitpos << P.acceptCriteria ;
        
        if (fin.eof ()) break ;
        paramvector.push_back (P) ;
    }
    fin.close () ;
    
    //Sorting according to bit position
    sort (paramvector.begin () , paramvector.end () , compare_bitPosition) ;


    char flag[N] ;
    for (int i = 0 ; i < paramvector.size () ; i ++) flag[i] = 'y' ;

    for (int i = paramvector.size () ; i < N ; i ++) flag[i] = 'n' ;


    bool flagcheck = true ;
    
    //loop for asking  user for  GTI parameter to be checked or not
    //if param_all_or_cust is 1 than all parameters to be checked.if 0 than ask user which parameter to be check
    if (all_Or_custom==0)
    { // custom

        for (int i = 0 ; i < paramvector.size () ; i ++)
        {
            do
            {
                cout << endl << "Check " << paramvector[i].name << " (y/n) ?" ;
                cin >> flag[i] ;
                cout.flush () ;
                flagcheck = validate (flag[i]) ; //validate the user input whether it is right or not
            }
            while (! flagcheck) ;
        }

    }

    //set 1 for y value of flag
    if (! valid_gtiflag && all_Or_custom)
    {
        for (int i = 0 ; i < paramvector.size () ; i ++)
        {
            flag[i] == 'y' ? 'n' : 'y' ;
            flag[i] == 'n' ? 'y' : 'n' ;
        }
    }
    else if (! valid_gtiflag && ! all_Or_custom)
    {
        for (int i = 0 ; i < paramvector.size () ; i ++)
        {
            flag[i] == 'y' ? 'n' : 'y' ;
            flag[i] == 'n' ? 'y' : 'n' ;
        }

    }

    //converting to BYTE
    for (i = 0 ; i < NBYTES ; i ++) final_ByteArray[i] = 0 ;

    cerr << "GTI flag to be matched  -  " ;
    for (i = 0 ; i < N ; i ++)
    {
        if (flag[i] == 'y')
        {
            cout << '1' ;
        }
        else
        {
            cout << '0' ;
        }

    }
    cerr << endl ;
   
    for (int i=0;i<cpu_VIS_temp.size ();i++)
    {
        
        if(flag[0]==1)
        {
        if(!is_in_range (cpu_VIS_temp[i],paramvector[0].min,paramvector[0].max))
        {
            
        rowno.push_back (i +1);
        continue;
        }
        
        }
        if(flag[1]==1)
        {
            if(!is_in_range (cpu_NUV_temp[i],paramvector[1].min,paramvector[1].max))
            {
                 
             rowno.push_back (i +1);
        continue;
           }
            
        }
        if(flag[2]==1)
        {
            if(!is_in_range (cpu_FUV_temp[i],paramvector[2].min,paramvector[2].max)){
                 
            rowno.push_back (i +1);
        continue;
            }
        }
        
    }
    
    
    
    
    
    
    
    //combine 8 bit and create a byte
//    for (int i = 0 ; i < NBYTES ; i ++)
//    {
//        p = 1 ;
//        for (int j = NO_OF_BITS_PER_BYTE * i ; j < NO_OF_BITS_PER_BYTE * i + NO_OF_BITS_PER_BYTE ; j ++ , p ++)
//        {
//
//            final_ByteArray[i] = (flag[j] << (NO_OF_BITS_PER_BYTE - p)) | final_ByteArray[i] ;
//        }
//
//    }
//    
//
//
//
//    int numBits = 0 ;


    //loop for anding 9 byte  GTI flag from level 1 file with 9 byte of user gti flag
//    for (i = 0 ; i < GTI_vect.size () / NBYTES ; i ++)//loop for number of rows of level1 science data file
//    {
//        numBits = 0 ;
//        for (int j = 0 ; j < NBYTES ; j ++)
//        {
//            if ((GTI_vect[i * NBYTES + j] & final_ByteArray[j]) == final_ByteArray[j]) numBits ++ ;
//        }
//        //  cout<<"the numbits are"<<numBits;
//        if (numBits != NBYTES) rowno.push_back (i + 1) ; //pushing row number where GTI flag is not satisfied
//
//    }


    return (EXIT_SUCCESS) ;
}


int DataIngest::decodeDarkFrames (char * filename , vector<float> &Arr , double &timeMedian)
{

    //    string path_f=(string)Darkframes_path+"/"+basename(filename);
    //     if(DirExists((char*)path_f.c_str ())){
    //                string command1 = "rm -rf  " + (string) path_f;
    //                LOG(INFO)<<"Executing......... "<<command1;
    //                system (command1.c_str ()) ;    //Executing the shell command for removing the output directory 
    //                LOG(INFO)<<path_f<<" directory deleted";
    //        }
    //    
    //    string command = (string) "mkdir -p " + (string) Darkframes_path+"/"+basename(filename); //
    //    system (command.c_str ()) ;

    /**reading the File for Creating Frames from the row wise data**/
    fitsfile *IM_mode_file_ptr ;
    int status = 0 ;
    fits_open_file (&IM_mode_file_ptr , filename , READONLY , &status) ;
    printError (status , "Error opening data ingest output file" , filename) ;
    fits_movnam_hdu (IM_mode_file_ptr , BINARY_TBL , SCIENCEDATA_HDUNAME , 0 , &status) ;
    printError (status , "Error moving to science  hdu" , filename) ;

    long nrows = 0 ;
    fits_get_num_rows (IM_mode_file_ptr , &nrows , &status) ;
    printError (status , "Error in getting the number of rows " , filename) ;
    unsigned short *frmcnt = new unsigned short[nrows];
    fits_read_col (IM_mode_file_ptr , TUSHORT , 11 , 1 ,1 , nrows , NULL, frmcnt , NULL , &status) ;
   printError (status , "Error in reading the column value of the framecount" , inFile) ;
    
   long totfrm= frmcnt[nrows-1]-frmcnt[0];
   try
   {
   float *temp_memcheck = new float[totfrm*IMG_DIM_DI*IMG_DIM_DI];
   delete[] temp_memcheck;
   }
   catch(bad_alloc)
   {
       LOG(ERROR)<<"***RAM+SWP memory is not sufficient  for Dark frame reading*** ";
       return(EXIT_FAILURE);
   }
    unsigned short psc ;
    unsigned short pixel[PIXELNO] ;

    int segflag = 1 ;

    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;

    unsigned int temptime ;

    long firstelem = 1 , nelement = 1 ;

    Arr.clear () ;
    int counter_pixel = 0 ;
    int num_frames = 0 ;
    vector<long> time_median ;
    ;
    /**LOOP for reading the Row wise data and creating the pixel array for each frame based on the frame number and pixels(1008 pixels)**/
    for (long firstrow = 1 ; firstrow <= nrows ; firstrow ++)
    {
        counter_pixel = 0 ;
        num_frames ++ ;

        fits_read_col (IM_mode_file_ptr , TUINT , 12 , firstrow , firstelem , nelement , NULL , &temptime , NULL , &status) ;
        printError (status , "Error in reading the time from the  file " , outFile) ;
        time_median.push_back (temptime) ;
        do
        {

            fits_read_col (IM_mode_file_ptr , TUSHORT , 8 , firstrow , firstelem , nelement , NULL , &psc , NULL , &status) ;
            printError (status , "Error Reading the column of the psc" , outFile) ;
            fits_read_col (IM_mode_file_ptr , TUSHORT , 15 , firstrow , firstelem , (nelement * PIXELNO) , NULL , (void *) pixel , NULL , &status) ;
            printError (status , "Error Reading the Column of Pixels" , outFile) ;

            for (int k = 0 ; k < PIXELNO ; k ++)
            {

                Arr.push_back (pixel[k]) ; //Arr contains the vector of  all the pixels(512*512) from all the dark frames.
                counter_pixel ++ ;
                if (counter_pixel == (IMG_DIM_DI * IMG_DIM_DI)) break ;

            }

            segflag = (int) (psc) >> 14 ;
            firstrow ++ ;
            if (firstrow > nrows) break ;
            if (segflag == 2) break ; //come out of loop after reading last segment in frame
        }
        while (segflag == 0 || segflag == 1) ; //segflag = 0 for continuation , 1 for first, 2 for last segment

        firstrow -- ;
    }
    if (time_median.size () % 2 == 0)
        timeMedian = ((time_median[time_median.size () / 2] + time_median[(time_median.size () / 2) - 1]) / 2.0) / 1000 ;
    else
        timeMedian = (time_median[time_median.size () / 2]) / 1000 ;

    fits_close_file (IM_mode_file_ptr , &status) ;
    printError (status , "Error in closing  file " , filename) ;
    return num_frames ;

}


int DataIngest:: readGTIparameters()
{
    int status=0;
    
    fitsfile *fptr;
     fits_open_file (&fptr , inFile_LBT , READONLY , &status) ;
    printError (status , "Error opening data ingest output file" ,  inFile_LBT) ;
    fits_movabs_hdu (fptr , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU in input information File" ,  inFile_LBT) ;
     fits_get_num_rows (fptr , &nrow_lbt , &status) ;
     printError (status , "Error in getting number of rows" ,  inFile_LBT) ;
   double *time_lbt=new double[nrow_lbt];      
     fits_read_col (fptr , TDOUBLE , 1 , 1 , 1 , nrow_lbt , NULL , time_lbt, NULL , &status) ;
    printError (status , "Error in reading the column of Exposure",inFile_LBT) ;
    float cpu_temp_vis,cpu_temp_nuv,cpu_temp_fuv;
    int count=0;
    for(int i=1;i<nrow_lbt;i++)
    {
        if(time_lbt[i-1]<=frm_time_double[count] && time_lbt[i]>frm_time_double[count])
        {           
    fits_read_col (fptr , TFLOAT , 25 , i , 1 , 1 , NULL , &cpu_temp_vis, NULL , &status) ;
    printError (status , "Error in reading the column of Exposure",inFile_LBT) ;
    fits_read_col (fptr , TFLOAT , 17 , i , 1 , 1 , NULL , &cpu_temp_nuv, NULL , &status) ;
    printError (status , "Error in reading the column of Exposure",inFile_LBT) ;
    fits_read_col (fptr , TFLOAT , 9 , i , 1 , 1 , NULL , &cpu_temp_fuv, NULL , &status) ;
    printError (status , "Error in reading the column of Exposure",inFile_LBT) ;         
    cpu_FUV_temp.push_back (cpu_temp_fuv);
    cpu_NUV_temp.push_back (cpu_temp_nuv);
    cpu_VIS_temp.push_back (cpu_temp_vis);
    count++;
        }
        //cout<<time_lbt[i-1]<<" "<<cpu_temp_vis<<" "<<cpu_temp_nuv<<" "<<cpu_temp_fuv;
    }
   
   
   
    fits_close_file(fptr,&status);
    printError (status , "Error in closing the lbt file" ,  inFile_LBT) ;
    return(EXIT_SUCCESS);
}
