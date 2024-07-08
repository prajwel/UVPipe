/* 
 * File:   uvtCosmicRayCorr.cpp
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
#include <uvtCosmicRayCorr.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<glog/logging.h>
#define WINSIZE 3
//Constructor -called when object is created


uvtCosmicRayCorr::uvtCosmicRayCorr ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}

//Destructor


uvtCosmicRayCorr::~ uvtCosmicRayCorr () 
{

    
 }

//parameter File reading 


int uvtCosmicRayCorr::read (int argc , char** argv)
{
    int status = 0 ;
    status = readParams (argc , argv , 1 , FNAME , "inputdatadir" , inputdatadir) ;
    if (status) return (EXIT_FAILURE) ;

    string tempfilepath = searchFile (inputdatadir , ".info") ;
    if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }

    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    char *obs_mode = new char[FLEN_VALUE] ;
    getKeywordVal (infofile_in , "OBS_MODE" , 2 , obs_mode) ;
    if (strcasecmp (obs_mode , "PC") == 0)
    {
        status = readParams (argc , argv , 1 , INT , "numFrames" , &nCompareFrames) ;
        if (status) return (EXIT_FAILURE) ;
    }

    status = readParams (argc , argv , 5 , REAL , "threshold_cr" , &Threshold ,
            FNAME , "outdir" , outdir ,
            BOOL , "clobber" , &clobber ,
            BOOL , "history" , &history ,
            STRING , "mode" , mode) ;
    if (status) return (EXIT_FAILURE) ;
    return (EXIT_SUCCESS) ;
}


int uvtCosmicRayCorr::read (char *input_datadir , float threshold , int compare_frames , char *out_dir , int clobber , int history)
{
    strcpy (inputdatadir , input_datadir) ;
    // this->Threshold = threshold ;
    this->Threshold = threshold ;
    this->nCompareFrames = compare_frames ;
    strcpy (outdir , out_dir) ;
    this->clobber = clobber ;
    this->history = history ;
}
//Parameter file content Display


void uvtCosmicRayCorr::display ()
{
    LOG (INFO) << endl ;
    LOG (INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG (INFO) << "             UVT COSMIC RAY  CORRECTION PARAMETERS      " << endl ;
    LOG (INFO) << "------------------------------------------------------------------------" ;
    ;
    LOG (INFO) << endl << "Input Frame List Directory                        : " << inputdatadir ;
    LOG (INFO) << endl << "Threshold value                                    : " << Threshold ;
    LOG (INFO) << endl << "Output Directory                               : " << outdir ;
    if (clobber == YES)
        LOG (INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG (INFO) << endl << "Overwrite                                         : NO" ;
    if (history == YES)
        LOG (INFO) << endl << "History                                             : YES" ;
    else
        LOG (INFO) << endl << "History                                              : NO" ;

    LOG (INFO) << endl << "------------------------------------------------------------------------" << endl ;

}

//Correction for the  Cosmic Ray process


int uvtCosmicRayCorr::uvtCosmicRayCorrProcess ()
{

    LOG (INFO) << "\nStarted Cosmic Ray Correction process" << endl ;
    ;

    /*---setting the output Directory path ----*/
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG (INFO) << endl << "\nModule Output Directory : " << moduleoutdir << endl ;

    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;

    /*----Searches the information file from the input Directory*/
    string  tempfilepath = searchFile (inputdatadir , ".info") ;
    
    /**Setting the path for the input information File**/
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (! (FileExists (infofile_in)))
    {
        LOG (ERROR) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY" ;
        return (EXIT_FAILURE) ;
    }
    LOG (INFO) << endl << "\nInformation File :" << infofile_in ;
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    /*----Reading the needed information from the information file----*/
    /**Opening the input information File**/
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input information file" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of the  input information file" , infofile_in) ;
    /**Reading the keyword  from the input information file**/
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    if (xsize <= 0 || ysize <= 0)
    {
        LOG (ERROR) << endl << "***Invalid xsize/ysize : xsize/ysize <= 0 ***\n" ;
        return (EXIT_FAILURE) ;
    }
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX" , infofile_in) ; //for creating name for output information file
    /**Setting thr output information  file**/
    sprintf (infofile_out , "%s/%s_cr.info" , moduleoutdir , nameprefix) ;
    LOG (INFO) << "\nOutput Information File " << infofile_out << endl ;
    /*----creating the output information file----*/
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information file" , infofile_out) ;
    char *ttype[] = {"SignalFrames" , "ExposureFrames"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "***Error in creating the table in out info file***" , infofile_out) ;
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "***Error in updating the key value of the NAMEPRFX***" , infofile_out) ;
    //writing the basic information to the output information file read from  input information file.
    datainfo.write (finfo_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "***Error in closing the out information file***" , infofile_out) ;
    //in case of IM mode 
    if (datainfo.getModeFlag () == IM)
    {
        readKeywords(infofile_in , 2 , 3 , TSTRING , "SIGDIR" , sigframedir ,
                TSTRING , "EXPDIR" , expframedir ,
                TINT , "NFILES" , &nframes 
                ) ;
 
        sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        expoframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        
        status=readColumnsFromFITS (infofile_in,2,2,TSTRING,1,sigframelist,nframes,TSTRING ,2,expoframelist,nframes);
        if(status)
        {
        LOG(INFO)<<"Error reading  the columns from the file"<<endl;
        return(EXIT_FAILURE);
        }

        /*in case of IM mode*/
        if (cosmicRayCorrIM ()) return (EXIT_FAILURE) ;
        //open  & update keywords of SIGDIR,EXPDIR,NFILES to output information file     
        //SIGDIR-path of output signal directory name
        //EXPDIR-path of output exposure directory name
        //NFILES-number of frames at output
          status=updateKeywords(infofile_out , 2 , 3 , TSTRING , "SIGDIR" , sigframedir ,
                TSTRING , "EXPDIR" , expframedir ,
                TINT , "NFILES" , &nframes 
                ) ;
          if(status) printError (status,"keyword not found in header",infofile_out);
        
           freeMemory (sigframelist , nframes , NAMESIZE) ;
    } //in case of the PC 
    else if (datainfo.getModeFlag () == PC)
    {//incase of PC mode
        fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
        printError (status , "***Error in reading the key value of the EVTFILES ***" , infofile_in) ;

        if (cosmicRayCorrPC ()) return (EXIT_FAILURE) ; //incase of PC mode
    }
    else
    {//invalid mode
        LOG (ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG (ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
    fits_close_file(finfo_in,&status);
    printError (status , "Error in closing the input  info file" , infofile_in) ;
    return (EXIT_SUCCESS) ;

}


int uvtCosmicRayCorr::cosmicRayCorrIM ()
{

    LOG (INFO) << endl << "Started Cosmic ray Correction  for IM mode \n" ;

    char **outsigframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    char **outexpframelist = allocateMemory<char>(nframes , NAMESIZE) ; //to store output exposure frame list  
    //creating signal directory
    char dir[FLEN_FILENAME] ;
    /**Setting the path for the output Signal Directory**/
    sprintf (dir , "%s%s" , moduleoutdir , sigframedir) ;
    /**Shell command for  creating the Directory**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell Command**/
    system (cmd.c_str ()) ;
    LOG (INFO) << endl << dir << " directory created" << endl ;

    //creating exposure directory

    /**Setting the output Exposure Frame directory **/
    sprintf (dir , "%s%s" , moduleoutdir , expframedir) ;
    /**Shell command for creating the Directory**/
    cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG (INFO) << endl << dir << "directory created" << endl ;
    LOG (INFO) << endl << "\nTotal number of frames - " << nframes ;
    int status = 0 ;
    /**
     * framedata-array for storing the input signal frame image pixels
     * framedata_cp-backup array of framedata
     * tempexpdata-array for storing the exposure frame data
     * tempexpdata_cp-backup array of framedata
     * @return 
     */
    float framedata[xsize * ysize] ;
    float framedata_cp[xsize * ysize] ;
    float tempexpdata[xsize * ysize] ;
    float tempexpdata_cp[xsize * ysize] ;


    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char errstr[512] ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;


    vector<int> X_pixel , Y_pixel ;
    if (history == YES)
    {
        getHistory (vhistorystr) ; //put history to vhistorystr vector
    }
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    char outfile1[FLEN_FILENAME] ;
    unsigned short frameno = 0 , frameno1 = 0 ;
    double frametime = 0 , frametime1 = 0 ;
    vector<long> cr_effected ;
    vector<double> x_crFailed , y_crFailed ;

    bool crflag = FALSE ;
    int cnt_cosmicAffected = 0 ;
    //processe for the  number of files in the Input Directory
    /**LOOP for processing for the number of frames **/
    LOG (INFO) << "\nPerforming Cosmic Ray Correction..." << endl ;
    //loop for number of frames
    for (int i = 0 ; i < nframes ; i ++)
    {
        sprintf (errstr , "Error at iteration number %d" , i) ;
        fitsfile *fptr , *fout , *fptr2 , *fout2 ;

        /*Signal frame reading*/
        crflag = FALSE ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[i]) ; //storing name of signal frame in infile
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input file" , infile) ;
        if (i == 0)
        {
            copyUsrkeywrdsTovect (fptr , key_records) ; //copying level1 keywords to the vector
        }
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ; //store signal frame data in 'framedata'
        printError (status , "Error in reading the pixel value of the input file" , infile) ;

        readKeywords(infile , 1 , 2 , TUSHORT , "FRAMENO" , &frameno ,TDOUBLE , "FRMTIME" , &frametime) ;
        
        /*Exposure frame reading*/
        sprintf (infile , "%s/%s/%s" , inputdatadir , expframedir , expoframelist[i]) ; //storing name of exposure frame in infile
        fits_open_file (&fptr2 , infile , READONLY , &status) ;
        printError (status , "Error in opening the input file" , infile) ;
        fits_read_pix (fptr2 , TFLOAT , fpixel , xsize*ysize , NULL , tempexpdata , NULL , &status) ; //store exposure frame data in 'tempexpdata'
        printError (status , "Error in reading the pixel value of the exposure input file" , infile) ;

        readKeywords(infile , 1 , 2 , TUSHORT , "FRAMENO" , &frameno1 , TDOUBLE , "FRMTIME" , &frametime1) ;

        for (int i = 0 ; i < xsize * ysize ; i ++)
        {
            framedata_cp[i] = framedata[i] ;
            tempexpdata_cp[i] = tempexpdata[i] ;
        }

        X_pixel.clear () ;
        Y_pixel.clear () ;
        /**
         loop for  identifying pixels which are above the threshold value
         * X_pixel-it contains x locations of identified pixels(above the threshold)
         * Y_pixel- it contains y locations of  identified pixels(above the threshold)
         * 
         */
        for (int j = 0 ; j < xsize * ysize ; j ++)
        {
            if (framedata[j] != INVALID_PIX_VALUE) //check whether pixel is declared invalid by previous module or not,.if yes then do nothing to it
            {

                if (framedata[j] > Threshold)
                {
                    X_pixel.push_back (j % ysize+1) ;
                    Y_pixel.push_back (j / ysize+1) ;
                    // tempexpdata[j] = framedata[j] = 0.0f ;
                    //crflag=TRUE;
                }
            }
        }

        /**
         * loop for checking the window around the identified pixel(above threshold),if one of the window pixel is also above the threshold value than identified pixel is
         * not CR Affected.
         * cnt_cosmicAffected-number of  above threshold pixels within window.
         * WINSIZE-window size  
         * @return 
         */
        for (int i = 0 ; i < X_pixel.size () ; i ++)
        {
            cnt_cosmicAffected = 0 ;
            for (int j = X_pixel[i] - WINSIZE / 2 ; j <= X_pixel[i] + WINSIZE / 2 ; j ++)
            {
                for (int k = Y_pixel[i] - WINSIZE / 2 ; k <= Y_pixel[i] + WINSIZE / 2 ; k ++)
                {

                    if (j < xsize && j > 0 && k < ysize && k > 0)
                    {
                        if (framedata[k * ysize + j] > Threshold)
                        {
                            cnt_cosmicAffected ++ ;
                        }
                    }
                }
            }
            //if only one pixel is above the threshold than it is CR effected.
            if (cnt_cosmicAffected == 1)
            {
                tempexpdata_cp[Y_pixel[i] * ysize + X_pixel[i]] = framedata_cp[Y_pixel[i] * ysize + X_pixel[i]] = -9999 ; //cr effected
                cr_effected.push_back (frameno) ;
                x_crFailed.push_back (X_pixel[i]) ;
                y_crFailed.push_back (Y_pixel[i]) ;

            }
        }

        /**Setting the path for the output Signal Frame **/
        sprintf (outfile , "%s/%s/%s_t%f_f%d_sig_cr.fits" , moduleoutdir , sigframedir , nameprefix , frametime , frameno) ;
        /**creating the output Sig frame**/
        fits_create_file (&fout , outfile , &status) ; //outfile - output signal frame
        printError (status , "Error in creating the output signal  frame" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the IMG for  Signal  frame" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata_cp , &status) ;
        printError (status , "Error in writing the pixel value of the output file" , outfile) ;

        //write history ,level-1 keywords ,origin,creator ,checksum and date to output signal  frame
        if (history == YES)
            writeHistory (outfile , vhistorystr) ;
       // writeUsrkeywordsFrmvect (outfile , key_records) ;
        copyUserKeywords (fptr,fout);
        writeCommonKeywords (fout , modulename) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output signal file" , outfile) ;
        /**Setting the path for the output Exposure Frame**/
        sprintf (outfile1 , "%s/%s/%s_t%f_f%d_exp_cr.fits" , moduleoutdir , expframedir , nameprefix , frametime1 , frameno1) ;

        //exposure frames processing 
        fits_create_file (&fout2 , outfile1 , &status) ; //outfile1 - output exposure frame
        printError (status , "Error in creating the output Exposure frame" , outfile1) ;
        fits_create_img (fout2 , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the IMG for  Exposure frame" , outfile1) ;
        fits_write_pix (fout2 , TFLOAT , fpixel , xsize*ysize , tempexpdata_cp , &status) ;
        printError (status , "Error in Writing  pixels the output Exposure  file" , outfile1) ;

        //write history ,level-1 keywords ,origin,creator ,checksum and date to output exposure frame
        if (history == YES)
            writeHistory (outfile1 , vhistorystr) ;

        // writeUsrkeywordsFrmvect (outfile1 , key_records) ;
        copyUserKeywords (fptr2,fout2);
        writeCommonKeywords (fout2 , modulename) ;
        fits_close_file (fout2 , &status) ;
        printError (status , "Error in closing the output Exposure  file " , outfile) ;

        //coping the output signal frame  to outputsigframelist  
        strcpy (outsigframelist[i] , basename (outfile)) ;

        //coping the output exposure frame  to outputexpframelist 
        strcpy (outexpframelist[i] , basename (outfile1)) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input Signal  file" , outfile) ;
        fits_close_file (fptr2 , &status) ;
        printError (status , "Error in closing the input  Exposure file" , outfile) ;
        cout<< "Total Files written=" << i << "     Remaining files= " << nframes - i <<" \r" ;

    } //end of loop for number of frames

    //LOG(INFO)<<"The Status value is "<<status<<endl;
    /*Writing  the information in output information File*/
    fitsfile *finfo_out ;
    /**Opening the output informaton File which has been generated earlier**/
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening  output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in  moving to the 2nd HDU of the out information file" , infofile_out) ;
    writeCommonKeywords (finfo_out , modulename) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the out information file" , infofile_out) ;
    status=writeColumnsToFITS (infofile_out,2,2,TSTRING,1,outsigframelist,nframes,TSTRING ,2,outexpframelist,nframes);
        if(status)
        {
        LOG(INFO)<<"Error writing  the columns to the fits file"<<endl;
        return(EXIT_FAILURE);
        }
   
    
    //write level-1 keywords to vector
    writeUsrkeywordsFrmvect (infofile_out , key_records) ;

    //write history to output information file
    if (history == YES)
        writeHistory (infofile_out , vhistorystr) ;

 
    /**Release Memory**/
    freeMemory (outsigframelist , nframes , NAMESIZE) ;
    freeMemory (outexpframelist , nframes , NAMESIZE) ;

    fitsfile *fit_crfailed ;
    char crFailedFrames[FLEN_FILENAME] ;
    /**Setting the path for the CRFailed file
     CRfailed.fits- Fiile containing frame number & pixel's  X and Y location where CR Failed. 
    
     **/

    sprintf (crFailedFrames , "%s/%s_CRFailed.fits" , moduleoutdir , nameprefix) ;
    LOG (INFO) << "Creating Cosmic Ray affected pixel's  list" << endl ;
    LOG (INFO) << "CRFailed file  " << crFailedFrames << endl ;

    /**Creating the CR file**/
    fits_create_file (&fit_crfailed , crFailedFrames , &status) ;
    printError (status , "Error in creating the CRFailed file" , crFailedFrames) ;
    char *exte ;
    char *ttypee[3] ;
    char *tforme[3] ;
    char *tunite[3] ;
    exte = "CRAffected" ;
    tunite[3] ;
    ttypee[0] = "FrameNo" ;
    tforme[0] = "1J" ;
    tunite[0] = "" ;
    ttypee[1] = "X" ;
    tforme[1] = "D" ;
    tunite[1] = "" ;
    ttypee[2] = "Y" ;
    tforme[2] = "D" ;
    tunite[2] = "" ;

    //writing Cr failed pixels to the table
    fits_create_tbl (fit_crfailed , BINARY_TBL , 0 , 3 , ttypee , tforme , tunite , exte , &status) ;
    printError (status , "Error in  creating the table of the Exposure file" , crFailedFrames) ;
    status=writeColumnsToFITS (crFailedFrames,2,3,TLONG,1,cr_effected.data (),cr_effected.size (),TDOUBLE ,2,x_crFailed.data (),x_crFailed.size (),TDOUBLE,3,y_crFailed.data (),y_crFailed.size ());
    if(status)
   {
    LOG(INFO)<<"Error writing   the columns to  the fits  file"<<endl;
    return(EXIT_FAILURE);
    }

    writeUsrkeywordsFrmvect (crFailedFrames , key_records) ;

    if (history == YES)
        writeHistory (crFailedFrames , vhistorystr) ;

    writeCommonKeywords (fit_crfailed , modulename) ;
    fits_close_file (fit_crfailed , &status) ;
    printError (status , "Error in  closing the file" , crFailedFrames) ;

    return (EXIT_SUCCESS) ;
}


int uvtCosmicRayCorr::cosmicRayCorrPC ()
{
    LOG (INFO) << endl << "Started Cosmic ray Correction  for PC mode\n" << endl ;
    fitsfile *fptr , *fout ;

    if (history == YES)
    {
        getHistory (vhistorystr) ; //putting history to vector vhistorystr
    }

    int status = 0 ;
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    char errstr[512] ;
    char eventfilename[FLEN_FILENAME] ;
    /**Setting the output EVENT file path**/
    sprintf (eventfilename , "%s_cr.events" , nameprefix) ;
    sprintf (outfile , "%s/%s" , moduleoutdir , eventfilename) ;
    /**Creating the output Eventfile **/
    fits_create_file (&fout , outfile , &status) ;
    printError (status , errstr) ;

    long nrows = 0 ;
    /**Setting the path for the input EVENT file**/

    sprintf (infile , "%s/%s" , inputdatadir , eventfile) ;
    LOG (INFO) << "Input Event File is  " << infile << endl ;
    LOG (INFO) << "output Event File is  " << outfile << endl ;

    //opening the input event file
    fits_open_file (&fptr , infile , READONLY , &status) ;
    printError (status , "Error in opening the input Event File" , infile) ;
    //copying the input event file to output event file
    fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
    printError (status , "Error in coping the input event file to the output event file" , outfile) ;
    //copying the level-1 keywords to the vector 
    copyUsrkeywrdsTovect (fptr , key_records) ;

    vector<long> del_CReffctedrows , frameno_crFailed ;
    vector<double> x_crFailed , y_crFailed ;

    long *frame_no ;
    unsigned short *xi , *yi ;
    double *t ;
    float *xf , *yf ;
    int cra , crg ;


    fits_movabs_hdu (fout , 2 , NULL , &status) ;
    printError (status , "Error in  moving to the 2nd HDU of the input  file" , infile) ;
    fits_get_num_rows (fout , &nrows , &status) ;
    printError (status , "Error in  reading the number of rows" , infile) ;

    /**
     * declaration of array for storing the column value of inout event file
     * @return 
     */
    frame_no = new long[nrows] ;
    xi = new unsigned short[nrows] , yi = new unsigned short[nrows] ;


    xf = new float[nrows] , yf = new float[nrows] ;
    /**Reading columns from input Event file**/
//    status=readColumnsFromFITS (outfile,2,5,TLONG,1,frame_no,nrows,TUSHORT ,4,xi,nrows,TUSHORT,6,yi,nrows,TFLOAT,5,xf,nrows,TFLOAT,7,yf,nrows);
//        if(status)
//        {
//        LOG(INFO)<<"Error reading  the columns from the file"<<endl;
//        return(EXIT_FAILURE);
//        }
       
        fits_read_col (fout , TLONG , 2 , 1 , 1 , nrows , NULL , frame_no, NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
         fits_read_col (fout , TUSHORT, 4 , 1 , 1 , nrows , NULL , xi, NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
         fits_read_col (fout , TUSHORT , 6 , 1 , 1 , nrows , NULL , yi, NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
         fits_read_col (fout , TFLOAT , 5 , 1 , 1 , nrows , NULL , xf, NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
         fits_read_col (fout , TFLOAT , 7 , 1 , 1 , nrows , NULL , yf, NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
               
    /*
     nEventsInFrame-> Array of number of Events in each frame;
     * GoodFrames->Array of frame number of frame which has minimum 1 event
     *CRAFFECTED-> Array of frame number whose number of Events are higher than the threshold
     *          
     */

    LOG (INFO) << "\nPerforming the CosmicRay Correction...." << endl ;
    long nFrames ;
    nFrames = frame_no[nrows - 1] ; //number of frames

    long nFramesCounter = 0 , j = 0 ;

    long nEventsInFrame[nFrames] , GoodFrames[nFrames] , CRAFFECTED[nFrames] , m = 0 ;
 
    for (long i = 0 ; i < nrows - 1 ; i ++)
    {
        if ((frame_no[i] == frame_no[i + 1]))
        {

            nFramesCounter ++ ; //counter for number of events in a frame
            if (i == nrows - 2)
            {
                nEventsInFrame[j] = nFramesCounter ;
                GoodFrames[j] = frame_no[i] ;
            }

        }
        else
        {
            nEventsInFrame[j] = nFramesCounter ;
            GoodFrames[j] = frame_no[i] ;
            j ++ ;
            nFramesCounter = 1 ;
            if ((frame_no[i + 1] - frame_no[i]) > 1)
            {
                m = (frame_no[i + 1] - frame_no[i]) ;
                for (int index = i ; index < m ; index ++)
                {
                    nEventsInFrame[j] = 0 ;
                    //this is added.
                    GoodFrames[j]=0;
                    j ++ ;
                }
            }
        }
    }
   
    cra = 0 ;
    crg = 0 ;
    int sum = 0 ;
    
    for (long i = nCompareFrames ; i < nFrames ; i ++)
    {
        sum = 0 ;
        if (nEventsInFrame[i] != 0)
        {
            for (int j = i - nCompareFrames ; j < i ; j ++)
            {
                sum = sum + nEventsInFrame[j] ;
            }
            sum = sum / nCompareFrames ;

            if ((nEventsInFrame[i] - sum) > Threshold)
            {
                CRAFFECTED[cra] = GoodFrames[i] ;
                cra ++ ;
            }
        }
    }
 

    int flag = 0 ;
    unsigned short Image_Array[xsize][ysize] ;

    for (int ii = 0 ; ii < xsize ; ii ++)
    {
        for (int jj = 0 ; jj < ysize ; jj ++)
        {
            Image_Array[ii][jj] = 0 ;
        }
    }
    for (long i = 0 ; i < nrows ; i ++)
    {
        flag = 0 ;
        for (long m = 0 ; m < cra ; m ++)
        {
            if (frame_no[i] == CRAFFECTED[m])
            {
                flag = 1 ;
                break ;
            }
        }

        if (flag == 1)
        {
            del_CReffctedrows.push_back (i + 1) ;
            frameno_crFailed.push_back (frame_no[i]) ;
            x_crFailed.push_back (xi[i] + xf[i]) ;
            y_crFailed.push_back (yi[i] + yf[i]) ;
        }
    }
    
    //deleting the rows from output event file where frame is found to be CR Affected  (if one of the event of the frame is CR affected)    
    fits_delete_rowlist (fout , del_CReffctedrows.data () , del_CReffctedrows.size () , &status) ;
    printError (status , "Error in Deleting the row list " , outfile) ;

    char expFilename[FLEN_FILENAME] ;
    fitsfile *fit_exp ;

    LOG (INFO) << "Creating Exposure File" << endl ;

    /**Setting the path for the exposure file**/
//    sprintf (expFilename , "%s/%s_Exposure.fits" , moduleoutdir , nameprefix) ;
//    /* Creation  of exposure File  
//     *Exposure.fits-file containing event(X,Y) and total number of exposures for that event
//     */
//
//    fits_create_file (&fit_exp , expFilename , &status) ;
//    printError (status , "Error in creating the Exposure file" , expFilename) ;
//    LOG (INFO) << "Exposure File name " << expFilename << endl ;
    char *exte = "Exposure" ;
    char *ttypee[3] ;
    char *tforme[3] ;
    char *tunite[3] ;
//    ttypee[0] = "X_Val" ;
//    tforme[0] = "1U" ;
//    tunite[0] = "" ;
//    ttypee[1] = "Y_Val" ;
//    tforme[1] = "1U" ;
//    tunite[1] = "" ;
//    ttypee[2] = "Counts" ;
//    tforme[2] = "1U" ;
//    tunite[2] = "" ;
//    fits_create_tbl (fit_exp , BINARY_TBL , 1 , 3 , ttypee , tforme , tunite , exte , &status) ;
//    printError (status , "Error in  creating the table of the Exposure file" , expFilename) ;
//
//    fits_movnam_hdu (fit_exp , ANY_HDU , exte , 0 , &status) ;
//    printError (status , "Error in  moving to 2nd HDU" , expFilename) ;
//    int row_num = 1 ;
//    
//    for (int i = 0 ; i < xsize ; i ++)
//    {
//        for (int j = 0 ; j < ysize ; j ++)
//        {
//            fits_write_col (fit_exp , TINT , 1 , row_num , 1 , 1 , &i , &status) ;
//            printError (status , "Error in  writing  the column of x(Exposure file)" , expFilename) ;
//            fits_write_col (fit_exp , TINT , 2 , row_num , 1 , 1 , &j , &status) ;
//            printError (status , "Error in  writing  the column of y(Exposure file)" , expFilename) ;
//            fits_write_col (fit_exp , TUSHORT , 3 , row_num , 1 , 1 , &Image_Array[j][i] , &status) ;
//            printError (status , "Error in  writing  the column of Image_Array(Exposure file)" , expFilename) ;
//            row_num ++ ;
//        }
//    }
//    
//    if (history == YES)
//        writeHistory (expFilename , vhistorystr) ; //write history to exposure file
//    writeUsrkeywordsFrmvect (expFilename , key_records) ; //write level-1 keywords  from vector to exposure file
//    fits_movabs_hdu (fit_exp , 1 , NULL , &status) ;
//    printError (status , "Error in  moving to the 2nd HDU of the input  file" , infile) ;
//    writeCommonKeywords (fit_exp , modulename) ; //writing origin,creator,checksum,date to output exposure file
//    fits_close_file (fit_exp , &status) ;
//    printError (status , "Error in  closing the  Exposure file" , expFilename) ;

    char crFailedFrames[FLEN_FILENAME] ;
    /**Setting the path for the CRFailed file**/
    sprintf (crFailedFrames , "%s/%s_CRFailed.fits" , moduleoutdir , nameprefix) ;
    LOG (INFO) << "Creating CR fail event list" << endl ;
    LOG (INFO) << "CRFailed file  " << crFailedFrames << endl ;
   
    fitsfile *fit_crfailed ;
    /**Setting the path for the CRFailed file
   CRfailed.fits- File containing frame number & all  event(x,y) corresponds to the frame for which  CR Failed. 
    
     **/
    fits_create_file (&fit_crfailed , crFailedFrames , &status) ;
    printError (status , "Error in creating the CRFailed file" , crFailedFrames) ;
    exte = "CRFailed" ;
    tunite[3] ;
    ttypee[0] = "frameNo" ;
    tforme[0] = "1J" ;
    tunite[0] = "" ;
    ttypee[1] = "X" ;
    tforme[1] = "D" ;
    tunite[1] = "" ;
    ttypee[2] = "Y" ;
    tforme[2] = "D" ;
    tunite[2] = "" ;
  
    fits_create_tbl (fit_crfailed , BINARY_TBL , 0 , 3 , ttypee , tforme , tunite , exte , &status) ;
    printError (status , "Error in  creating the table of the Exposure file" , crFailedFrames) ;
   
    status=writeColumnsToFITS (crFailedFrames,2,3,TLONG,1,frameno_crFailed.data (),frameno_crFailed.size (),
            TDOUBLE ,2,x_crFailed.data (),x_crFailed.size (),TDOUBLE,3,y_crFailed.data (),y_crFailed.size ());
        if(status)
        {
        LOG(INFO)<<"Error writing  the columns to  the file"<<endl;
        return(EXIT_FAILURE);
        }

    if (history == YES)
        writeHistory (crFailedFrames , vhistorystr) ;
    writeUsrkeywordsFrmvect (crFailedFrames , key_records) ;
    fits_movabs_hdu (fit_crfailed , 1 , NULL , &status) ;
    printError (status , "Error in  moving to the 1st  HDU of the input  file" , infile) ;
    
    writeCommonKeywords (fit_crfailed , modulename) ;
   
    fits_close_file (fit_crfailed , &status) ;
    printError (status , "Error in  closing the  Exposure file" , crFailedFrames) ;
   
    delete[] frame_no ;
    
    if (history == YES)  writeHistory (outfile , vhistorystr) ;
        
    //write origin,creator,checksum,date to output event file
    writeCommonKeywords (fout , modulename) ;
    fits_close_file (fout , &status) ;
    printError (status , "Error in closing the output event file" , infile) ;


    /**open output information File and update EVTFILE  which has been generated earlier.**/
    status=updateKeywords(infofile_out , 2 , 1 , TSTRING , "EVTFILE" , basename(outfile)) ;
          if(status) printError (status,"keyword not found in header",infofile_out);

    writeUsrkeywordsFrmvect (infofile_out , key_records) ;
   
    if (history == YES)
        writeHistory (infofile_out , vhistorystr) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the input file" , infile) ;
    
    return (EXIT_SUCCESS) ;
}


int uvtCosmicRayCorr::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char * Threshold = new char[20] ;
    sprintf (Threshold , "%f" , Threshold) ;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"inputdatadir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Threshold value " + (string) Threshold) ;
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
