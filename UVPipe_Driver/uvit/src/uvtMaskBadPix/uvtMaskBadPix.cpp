/* 
 * File:   uvtMaskBadPix.cpp
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
#include "uvtMaskBadPix.h"
#include<glog/logging.h>
#include<pthread.h>


uvtMaskBadPix::uvtMaskBadPix ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    
}


uvtMaskBadPix::~ uvtMaskBadPix ()
{
    //delete[] badpixdata ;
}


int uvtMaskBadPix::read (int argc , char** argv)
{
    int status = 0 ;
    readParams (argc , argv , 1 , FNAME , "inputdatadir" , inputdatadir) ;

    /*---------Code to find data mode from information file inside input data dir-----------*/
    string tempfilepath = searchFile (inputdatadir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG (INFO) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }

    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    char *obs_mode = new char[FLEN_VALUE] ;
    getKeywordVal (infofile_in , "OBS_MODE" , 2 , obs_mode) ;
    if (strcasecmp (obs_mode , "PC") == 0)
    {
        //LOG(INFO)<<endl<<"data mode is "<<obs_mode<<endl;     
        readParams (argc , argv , 1 , REAL4 , "threshold_multph" , &threshold_multph) ;
        if (status) return (EXIT_FAILURE) ;
    }
    status = readParams (argc , argv , 6 ,
            FNAME , "caldbDir" , caldbDir ,
            FNAME , "outdir" , outdir ,
            BOOL , "clobber" , &clobber ,
            BOOL , "history" , &history ,
            STRING , "mode" , &mode
            ) ;

    if (status) return (EXIT_FAILURE) ;

    return (EXIT_SUCCESS) ;
}


int uvtMaskBadPix::read (char *inputdatadir , char *caldbDir , char *outdir , float thr , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->caldbDir , caldbDir) ;
    strcpy (this->outdir , outdir) ;
    this->threshold_multph = thr ;
    this->clobber = clobber ;
    this->history = history ;
    //    this->exposureFlag = exposure_Flag ;
    return (EXIT_SUCCESS) ;
}


void uvtMaskBadPix::display ()
{
    LOG (INFO) << endl ;
    LOG (INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG (INFO) << "             UVT FILTER BADPIX PARAMETERS              " << endl ;
    LOG (INFO) << "------------------------------------------------------------------------" ;
    LOG (INFO) << endl << "Input data directory  : " << inputdatadir ;
    LOG (INFO) << endl << "CALDB directory : " << caldbDir ;
    LOG (INFO) << endl << "Output Directory : " << outdir ;
    string tempfilepath = searchFile (inputdatadir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG (INFO) << endl << "***Error in finding info file***" ;
        exit(1);
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    char *obs_mode = new char[FLEN_VALUE] ;
    getKeywordVal (infofile_in , "OBS_MODE" , 2 , obs_mode) ;
    if (strcasecmp (obs_mode , "PC") == 0)
    {
        LOG (INFO) << endl << "Multiple photon event : " << threshold_multph ;
    }
    if (clobber == YES)
        LOG (INFO) << endl << "Overwrite : YES" ;
    else
        LOG (INFO) << endl << "Overwrite : NO" ;
    if (history == YES)
        LOG (INFO) << endl << "History    : YES" ;
    else
        LOG (INFO) << endl << "History     : NO" ;
    LOG (INFO) << endl << "------------------------------------------------------------------------" << endl ;
}


int uvtMaskBadPix::uvtMaskBadPixProcess ()
{
    LOG (INFO) << endl << "Filter Bad Pixel process started" << endl ;

    /**Setting the path for module output directory**/
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
   
     //check  existence  of output directory if output  exist and 
    //clobber =YES than remove that directory and recreate it,if exist and clobber =NO than exit from module.
    //if directory is not exist than create it .
    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;
    LOG (INFO) << endl << moduleoutdir << "  directory created" ;

    /*Searching the Input Directory for the input information file*/
    string tempfilepath = searchFile (inputdatadir , "info") ;
    if (tempfilepath == " ")
    {
        LOG (INFO) << endl << "Error in finding info file" ;
        return (EXIT_FAILURE) ;
    }

    /**setting path for input  information File**/
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    LOG (INFO) << "\nInput Information File " << infofile_in << endl ;


    /*----Opening the Input information file and read needed information from it ----*/
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY, &status) ;
    printError (status , "Error in opening the input info file" , infofile_in) ;
    fits_read_key (finfo_in , TINT , "WIN_X_SZ" , &win_xsize , NULL , &status) ;
    printError (status , "Error in reading the key value of the WIN_X_SZ" , infofile_in) ; //for creating name for output information file
    fits_read_key (finfo_in , TINT , "WIN_X_SZ" , &win_ysize , NULL , &status) ;
    printError (status , "Error in reading the key value of the  WIN_X_SZ" , infofile_in) ; //for creating name for output information file
    LOG(INFO)<<"Successfully opened"<<endl;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU in input information file" , infofile_in) ;

    /**Reading the information from Keywords **/
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ; //Width of frame
    ysize = datainfo.getYsize () ; //Height of frame
    if (xsize <= 0 || ysize <= 0)
    {
        LOG (INFO) << endl << "***Invalid xsize/ysize : xsize/ysize <= 0 ***\n" ;
        return (EXIT_FAILURE) ;
    }

    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX" , infofile_in) ; //for creating name for output information file
     
    /*----creating the output information path----*/
    sprintf (infofile_out , "%s/%s_bp.info" , moduleoutdir , nameprefix) ;
    LOG (INFO) << "\nOutput Information File :" << infofile_out << endl ;

    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the file" , infofile_out) ; //for creating name for output information file

    char *ttype[] = {"SignalFileList"} ;
    char *tform[] = {"A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table" , infofile_out) ;

    
    /*----Reading caldb badpixel file path from caldb Directory----*/
    LOG (INFO) << "Getting Path of badpixel correction  file of CalDB directory" << endl ;
    
    string tempname = caldb_handler.getBadPixelFile (datainfo.getDetector () , datainfo.getObsMode () , xsize , ysize , caldbDir) ;
    if (tempname == " ")
    {
        LOG (INFO) << endl << "Couldn't find bad pixel file from caldb" << endl ;
        return (EXIT_FAILURE) ;
    }
    /**Getting the absolute path for the badpixfile **/
    joinStrings (badpixfile , 1 , tempname.c_str()) ;
    LOG (INFO) << endl << "\nBad pixel  Correction file :" << badpixfile ;
    
    
    /*----reading the bad pixel file from caldb----*/
    status = readBadpixFile () ;
    if (status)
    {
        LOG (INFO) << endl << "Reading of the bad pixel correction file from the calDB fails" << endl ;
        return (EXIT_FAILURE) ;
    }
 fits_close_file (finfo_in , &status) ;
 printError (status , "Error in closing the input information File" , infofile_in) ;
 LOG(INFO)<<"Successfully closed"<<endl;
    /*----check whether the MODE is IM or PC----*/
    if (datainfo.getModeFlag () == IM)
    {
        readKeywords (infofile_in , 2 , 2 , TSTRING , "SIGDIR" , sigframedir ,
                TINT , "NFILES" , &nframes) ;

        strcpy (expframedir , "ExposureFrames") ; //setting exposure frame dir for output to be created
           fits_open_file (&finfo_in , infofile_in , READWRITE, &status) ;
           printError (status , "Error in opening the input info file" , infofile_in) ;
           fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
                printError (status , "Error in moving to 2nd HDU in input information file" , infofile_in) ;

        //reading signal  frame names from information file into vector 'sigframelist'
        sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) sigframelist , NULL , &status) ;
        printError (status , "Error in reading the column of input SignalFrame list" , infofile_in) ;
        
        fits_close_file(finfo_in,&status);
        printError (status , "Error in closing file" , infofile_in) ;
        if (maskBadpixIM ()) //Function to filter bad pixels for IM
            return (EXIT_FAILURE) ;

        //Updating keywords 'SIGDIR' , 'EXPDIR', 'NFILES' in output information file
        updateKeywords (infofile_out , 2 , 3 , TSTRING , "SIGDIR" , sigframedir ,
                TSTRING , "EXPDIR" , expframedir ,
                TINT , "NFILES" , &nframes) ;

        freeMemory (sigframelist , nframes , NAMESIZE) ;
    }
    else if (datainfo.getModeFlag () == PC) //incase of the PC mode
    {
        //Reading event filename and image filename from input information file
        
        readKeywords (infofile_in , 2 , 2, TSTRING , "EVTFILE" , eventfile ,
                TSTRING , "IMGFILE" , imgfile) ;
        //LOG(INFO)<<"Succ"<<endl;
        if (maskBadPixPC ()) //Function to filter bad pixels for PC mode
            return (EXIT_FAILURE) ;
    }
    else
    { //invalid mode
        LOG (ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG (ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
    
   

    /*writing the keyword information to the output information File*/
    datainfo.write (finfo_out) ; //writing basic data information
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the key value of the NAMEPRFX" , infofile_out) ; //for creating name for output information file

    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information File" , infofile_out) ;

    LOG (INFO) << endl << "uvtMaskBadPix process completed successfully" << endl ;
    return (EXIT_SUCCESS) ;
}


int uvtMaskBadPix::maskBadpixIM ()
{
    LOG (INFO) << endl << "Started Filtering bad pixels for IM mode" << endl ;
    char **outsigframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    char **outexpframelist = allocateMemory<char>(nframes , NAMESIZE) ; //to store output exposure frame list  

    //creating signal directory
    char dir[FLEN_FILENAME] ;

    /**Setting the output directory path for signal frames**/
    sprintf (dir , "%s/%s" , moduleoutdir , sigframedir) ;

    /**Shell command for creating the output Directory**/
    string cmd = "mkdir -p " + (string) dir ;

    /**Executing the shell command **/
    system (cmd.c_str ()) ;
    LOG (INFO) << endl << dir << " directory created" << endl ;

    /*----Creating the Exposure Directory----*/
    /**Setting the exposure frame path**/
    sprintf (dir , "%s/%s" , moduleoutdir , expframedir) ;
    cmd = "mkdir -p " + (string) dir ;
    system (cmd.c_str ()) ;
    LOG (INFO) << endl << dir << " directory created" << endl ;

    LOG (INFO) << endl << "\nTotal number of frames :" << nframes << endl ;
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
    if (history == YES) getHistory (vhistorystr) ;

    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;

    badExparray = new float[xsize * ysize] ;
    //for exposure frames
    for (int p = 0 ; p < xsize * ysize ; p ++)
    {
        badExparray[p] = 0.0 ;
        badExparray[p] = badpixdata[p] * datainfo.getIntegrationTime () ;
    }

    LOG (INFO) << "Performing Bad Pixel Filtering.." << endl ;

    /**LOOP for the processing the  total number of frames **/
    for (int i = 0 ; i < nframes ; i ++)
    {
        sprintf (errstr , "Error at iteration number %d" , i) ;
        fitsfile *fptr , *fout ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[i]) ;
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input  file" , infile) ;
        copyUsrkeywrdsTovect (fptr , key_record) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
        printError (status , "Error in reading the pixels from input file" , infile) ;

        readKeywords (infile , 1 , 3 , TDOUBLE , "INT_TIME" , &integrationtime ,
                TUSHORT , "FRAMENO" , &frameno ,
                TDOUBLE , "FRMTIME" , &frametime) ;

        /*----Applying the Filterbadpix correction----*/
        for (int pixno = 0 ; pixno < xsize * ysize ; pixno ++)
        {
            if (framedata[pixno] != INVALID_PIX_VALUE){
                framedata[pixno] = framedata[pixno] * badpixdata[pixno] ;
                if(badpixdata[pixno]==0) framedata[pixno]=INVALID_PIX_VALUE;
            }
        }
        /*----Setting the output signal  frame path----*/
        sprintf (outfile , "%s/%s/%s_t%.4f_f%d_sig_bp.fits" , moduleoutdir , sigframedir , nameprefix , frametime , frameno) ;

        /**Creating the output File**/
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the img" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata , &status) ;
        printError (status , "Error in writing the pixel  value to the output file" , outfile) ;

        if (history == YES) writeHistory (outfile , vhistorystr) ;
        copyUserKeywords (fptr , fout) ;
        writeCommonKeywords (fout , modulename) ;

        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the file output file" , outfile) ;
        /**writing the HISTORY **/


        strcpy (outsigframelist[i] , basename (outfile)) ; //storing output signal frame name to output signal frame list

        /*----Setting the output exposure   frame path----*/
        sprintf (outfile , "%s/%s/%s_t%f_f%d_exp_bp.fits" , moduleoutdir , expframedir , nameprefix , frametime , frameno) ;

        /**Creating  the output Exposure Frame**/
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the out file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the IMG file" , outfile) ;

        /**writing the new Exposure array to the output Exposure  Frame**/
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , badExparray , &status) ;
        printError (status , "Error in writing the pixels to the output file" , outfile) ;

        /**Coping the USER kEYWORD to the output Frame**/
        if (history == YES) writeHistory (outfile , vhistorystr) ;
        copyUserKeywords (fptr , fout) ;
        writeCommonKeywords (fout , modulename) ;

        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output file" , outfile) ;
        /**Writing HISTORY to the Exposure Frame**/

        strcpy (outexpframelist[i] , basename (outfile)) ; //storing output exposure frame name to output exposure frame list

        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input file" , outfile) ;
        cout<< "Total Files written=" << i << "     Remaining files= " << nframes - i <<" \r" ;
    } //end of loop for number of frames

    vhistorystr.clear () ;

    /*----Updating the Information to the information file----*/
    /**Open the output information File which has been created before and update some of the information to the information File**/
    LOG (INFO) << "\nUpdating the output information file... " << endl ;
    fitsfile *finfo_out ;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "***Error in opening the out information file***" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "***Error in moving to the 2nd HDU***" , infofile_out) ;
   
    char *ttype = {"ExposureFileList"} ;
    char *tform = {"A256"} ;
    fits_insert_col (finfo_out , 2 , ttype , tform , &status) ;
    printError (status , "***Error in inserting the column in out information file***" , infofile_out) ; //inserting column to add exposure frame list
    
     status=writeColumnsToFITS (infofile_out,2,2,TSTRING,1,outsigframelist,nframes,TSTRING ,2,outexpframelist,nframes);
     if(status){
         LOG(INFO)<<"Error writing  the columns to  the  fits file"<<endl;
        return(EXIT_FAILURE);
     }
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in output information file" , infofile_out) ;
    writeUsrkeywordsFrmvect (infofile_out , key_record) ;
    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    

    /**Releasing the Memory**/
    freeMemory (outsigframelist , nframes , NAMESIZE) ;
    freeMemory (outexpframelist , nframes , NAMESIZE) ;

    delete[] badExparray ;

    return (EXIT_SUCCESS) ;
}

//this function identifies bad pixels in the data using bad pixel information from caldb and adds a column for bad flag to event file


int uvtMaskBadPix::maskBadPixPC ()
{
    LOG (INFO) << endl << "Started Filtering bad pixel for PC mode" << endl ;
    char file_in[NAMESIZE] , file_out[NAMESIZE] , file_out1[NAMESIZE] ; //for input and output event file

    //opening eventfile
    /**Setting the input Event File path**/
    sprintf (file_in , "%s/%s" , inputdatadir , eventfile) ; //taking event file full path
    LOG (INFO) << "Input Event File " << file_in << endl ;

    /**setting the output event file path**/
    sprintf (file_out , "%s/%s_bp.events" , moduleoutdir , nameprefix) ;
    fitsfile *fevt_in , *fevt_out , *fimg , *fimg_out , *finfo_out ;
    int status = 0 ;

    fits_open_file (&fevt_in , file_in , READWRITE , &status) ;
    printError (status , "Error in opening the input Event file" , file_in) ;
    
    //coping keywords from input event file to the vector
    copyUsrkeywrdsTovect (fevt_in , key_record) ;

    //creating the output event file
    fits_create_file (&fevt_out , file_out , &status) ;
    printError (status , "Error in creating the  output event file" , file_out) ;
    //copy input event filt to output event file
    fits_copy_file (fevt_in , fevt_out , 1 , 1 , 1 , &status) ;
    printError (status , "Error in coping the input event file to the output event file" , file_out) ;
    fits_close_file (fevt_in , &status) ;
    printError (status , "Error in closing the input event file" , file_in) ;



    int ncols = 0 ;
    long nrows = 0 ;
    int x = 0 , y = 0 ;
    unsigned char MaxMin , multflag ; //MaxMin is Max minus Min column value
    int colnum_Ix , colnum_Iy , colnum_MaxMin ;
    unsigned char badval = 0 ;
    vector<string> vhistorystr ;
    LOG (INFO) << "Performing BadPixel Filtering..." << endl ;
  
    fits_movabs_hdu (fevt_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU" , file_out) ;

    char *ttype[] = {"BAD_FLAG" , "MULT_PH_EVT"} ;
    char *tform[] = {"B" , "B"} ;

    fits_get_num_cols (fevt_out , &ncols , &status) ;
    printError (status , "Error in getting the number of  columns" , file_out) ;
    fits_get_num_rows (fevt_out , &nrows , &status) ;
    printError (status , " Error in getting the number of rows in output event file" , file_out) ;

    //adding the column to output event fjle
    fits_insert_cols (fevt_out , ncols + 1 , 2 , ttype , tform , &status) ;
    printError (status , "Error in inserting the column" , file_out) ;

    //get column numbers for Ix , Iy and MaxMin
    fits_get_colnum (fevt_out , CASEINSEN , "Ix" , &colnum_Ix , &status) ;
    printError (status , "Error in getting the column number of the Ix column" , file_out) ;
    fits_get_colnum (fevt_out , CASEINSEN , "Iy" , &colnum_Iy , &status) ;
    printError (status , "Error in getting thr column number of the Iy column" , file_out) ;
    fits_get_colnum (fevt_out , CASEINSEN , "Max-Min" , &colnum_MaxMin , &status) ;
    printError (status , "Error in getting the column number of the Max-Min" , file_out) ;
    LOG (INFO) << endl << "Columns are " << colnum_Ix << "  , " << colnum_Iy << "    No of Rows " << nrows << " No of Columns :" << ncols ;

    //for all rows
    for (long rowno = 0 ; rowno < nrows ; rowno ++)
    {
        //reading columns of x,y               
        fits_read_col (fevt_out , TUSHORT , colnum_Ix , rowno + 1 , 1 , 1 , NULL , &x , NULL , &status) ;
        printError (status , "Error in reading the column x" , file_out) ;
        fits_read_col (fevt_out , TUSHORT , colnum_Iy , rowno + 1 , 1 , 1 , NULL , &y , NULL , &status) ;
        printError (status , "Error in reading the column y" , file_out) ;
       
        if (((y * xsize + x) > (xsize * ysize)) || (y * xsize + x) < 0)
        {
            LOG (ERROR) << endl << "(" << x << " , " << y << ")" ;
            LOG (ERROR) << endl << "***Index out of range exception***" << endl ;
            return (EXIT_FAILURE) ;
        }
        badval = static_cast<unsigned char> (badpixdata[y * xsize + x]) ; //float value of badpixel will be casted to unsigned char 

        //writing badval column to output event file
        fits_write_col (fevt_out , TBYTE , ncols + 1 , rowno + 1 , 1 , 1 , &badval , &status) ;
        printError (status , "Error in writing the column" , file_out) ;

        //reading column of MaxMin(for multiple photon event) from output event file
        fits_read_col (fevt_out , TBYTE , colnum_MaxMin , rowno + 1 , 1 , 1 , NULL , &MaxMin , NULL , &status) ;
        printError (status , "Error in reading the column of MAX-MIN" , file_out) ;
        
       
        if ((float)MaxMin > threshold_multph) multflag = 0 ; //if max-min > user_threshold then flag=0 else flag=1
        else multflag = 1 ;
        fits_write_col (fevt_out , TBYTE , ncols + 2 , rowno + 1 , 1 , 1 , &multflag , &status) ;
        printError (status , "Error in writing the multiple photon event column" , file_out) ;
    } //end of loop for row number
    // } // end of loop for hdu number


    LOG (INFO) << "Writing History" << endl ;
    if (history == YES)
    {
        getHistory (vhistorystr) ; //reading history and storing it to the vhistorystr vector
        writeHistory (file_out , vhistorystr) ; //write history to the output event file
    }
    fits_movabs_hdu (fevt_out , 1 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU" , file_out) ;
    //write origin,creator ,date checksum to output event file
    writeCommonKeywords (fevt_out , modulename) ;
    fits_close_file (fevt_out , &status) ;
    printError (status , "Error in closing the output  fits file" , file_out) ;

    //updating "EVTFILE" in information file        
    LOG (INFO) << "Updating the output Information file" << endl ;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU of the output information file" , infofile_out) ;

    //----------------------------------------------------------------------------------------
    //Updating Image file with bad pixel information
    sprintf (file_in , "%s/%s" , inputdatadir , imgfile) ; //taking event file full path
    sprintf (file_out1 , "%s/%s_bp.image" , moduleoutdir , nameprefix) ;

    fits_open_file (&fimg , file_in , READWRITE , &status) ;
    printError (status , "Error in opening the img file " , file_in) ;
    fits_create_file (&fimg_out , file_out1 , &status) ;
    printError (status , "Error in creating the output img file" , file_out1) ;
    fits_copy_header (fimg , fimg_out , &status) ;
    printError (status , "Error in coping the header" , file_out1) ;

    float imagedata[xsize * ysize] ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;

    fits_read_pix (fimg , TFLOAT , fpixel , xsize*ysize , NULL , imagedata , NULL , &status) ;
    printError (status , "Error in reading the  pixels from input image file" , file_in) ;

    fits_close_file (fimg , &status) ;
    printError (status , "Error in closing the img file" , file_in) ;

    for (int pixno = 0 ; pixno < xsize * ysize ; pixno ++)
        imagedata[pixno] = imagedata[pixno] * badpixdata[pixno] ;

    LOG (INFO) << endl << "Writing Image..." ;
    fits_write_pix (fimg_out , TFLOAT , fpixel , xsize*ysize , imagedata , &status) ;
    printError (status , "Error in writing the pixels to the output IMG file" , file_out1) ;
    //writing history and common keywords to the ourput image file.

    if (history == YES) writeHistory (file_out1 , vhistorystr) ;
    writeCommonKeywords (fimg_out , modulename) ;

    fits_close_file (fimg_out , &status) ;
    printError (status , "Error in closing the img file" , file_out1) ;
    LOG (INFO) << endl << "Image writing completed" ;

    //updating "IMGFILE"  and"EVTFILE" in output information file
    updateKeywords (infofile_out , 2 , 2 , TSTRING , "EVTFILE" , basename (file_out) ,
            TSTRING , "IMGFILE" , basename (file_out1)) ;

    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file" , infofile_out) ;
    
    //writing level 1 keywords and history to the output information file
    writeUsrkeywordsFrmvect (infofile_out , key_record) ;
    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    LOG (INFO) << endl << "Event filtering for bad pixel completed" << endl ;
    return (EXIT_SUCCESS) ;
}


int uvtMaskBadPix::readBadpixFile ()
{
    LOG (INFO) << endl << "Reading bad pixel file from caldb........" ;
    fitsfile *fbadpix ;
    int status = 0 ;
    badpixdata = new float[xsize * ysize] ;
    //openig badpixel file from caldb.
    fits_open_file (&fbadpix , badpixfile , READONLY , &status) ;
    printError (status , "Error in opening the badPixel file" , badpixfile) ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;

    /**Reading the BADPIXEL DATA from the badpixel File**/
    fits_read_pix (fbadpix , TFLOAT , fpixel , xsize*ysize , NULL , badpixdata , NULL , &status) ;
    printError (status , "Error in reading the pixels from caldb file" , badpixfile) ;
    fits_close_file (fbadpix , &status) ;
    printError (status , "Error in closing the file" , badpixfile) ;
    return (EXIT_SUCCESS) ;
}


int uvtMaskBadPix::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" badpix file " + (string) badpixfile) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Signal Frame Directory=" + (string) sigframedir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Exposure frame directory=" + (string) expframedir) ;
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


