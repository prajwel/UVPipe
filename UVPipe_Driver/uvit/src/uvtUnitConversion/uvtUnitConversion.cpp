/* 
 * File:   uvtUnitConversion.cpp
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
#include "fitsio.h"  /* required by every program that uses CFITSIO  */
#include "pil.h"
#include "pil_error.h"
#include "uvtUnitConversion.h"
#include<pthread.h>
#include<glog/logging.h>


/*Constructor*/
uvtUnitConversion::uvtUnitConversion ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}


int uvtUnitConversion::read (int argc , char** argv)
{
    int status = 0 ;

    status = readParams (argc , argv , 1 , FNAME , "inputdatadir" , inputdatadir) ;
    if (status) return (EXIT_FAILURE) ;

    //Finf observation mode of data from input data directory 
    string tempfilepath = searchFile (inputdatadir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    char *obs_mode = new char[FLEN_VALUE] ;
    getKeywordVal (infofile_in , "OBS_MODE" , 2 , obs_mode) ;

    if (strcasecmp (obs_mode , "IM") == 0)
    {
        status = readParams (argc , argv , 1 , BOOL , "darkframe_flag" , &darkFrame_Flag) ;
        if (status) return (EXIT_FAILURE) ;

    }

    status = readParams (argc , argv , 4 , FNAME , "outdir" , outdir ,
            BOOL , "clobber" , &clobber ,
            BOOL , "history" , &history ,
            STRING , "mode" , mode) ;
    
    return (EXIT_SUCCESS) ;
}


int uvtUnitConversion::read (char *inputdatadir , char *outdir , int clobber , int history , int dark_frame_flag , char *dstart_path , char *dend_path)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;
    if (dstart_path != NULL)
    {
        strcpy (this->dstartpath , dstart_path) ;
        strcpy (this->dendpath , dend_path) ;
    }
    this->clobber = clobber ;
    this->history = history ;
    this->darkFrame_Flag = dark_frame_flag ;
    return (EXIT_SUCCESS) ;
}


void uvtUnitConversion::display ()
{
    LOG (INFO) << endl ;
    LOG (INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG (INFO) << "             UVT UNIT CONVERSION PARAMETERS       " << endl ;
    LOG (INFO) << "------------------------------------------------------------------------" ;
    LOG (INFO) << endl << "Input Data Directory    : " << inputdatadir ;
    LOG (INFO) << endl << "Output Directory   : " << outdir ;
    if (clobber == YES)
        LOG (INFO) << endl << "Overwrite                              : YES" ;
    else
        LOG (INFO) << endl << "Overwrite                              : NO" ;
    if (history == YES)
        LOG (INFO) << endl << "History                                  : YES" ;
    else
        LOG (INFO) << endl << "History                                   : NO" ;
    LOG (INFO) << endl << "------------------------------------------------------------------------" << endl ;
}


int uvtUnitConversion::uvtUnitConversionProcess ()
{

    /**Creating Output Directory**/
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG (INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;

    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;

    LOG (INFO) << endl << moduleoutdir << "  directory created" ;

    /**Serching for the information file from the input directory**/
    string tempfilepath = searchFile (inputdatadir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG (INFO) << endl << "Error in finding info file" ;
        return (EXIT_FAILURE) ;
    }

    /*Setting the input information file path*/
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    LOG (INFO) << endl << "Input Information File :" << infofile_in ;
    int status = 0 ;

    /**Opening the input Information  file**/
    fitsfile *finfo_in ;
    fitsfile *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input information file" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU" , infofile_in) ;
    /*Reading the keyword information from the input File*/
    datainfo.getInfo (finfo_in) ;
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    if (xsize <= 0 || ysize <= 0)
    {
        LOG (INFO) << endl << "***Invalid xsize/ysize***\n" ;
        return (EXIT_FAILURE) ;
    }

    readKeywords (infofile_in , 2 , 1 , TSTRING , "NAMEPRFX" , nameprefix) ; //Read NAMEPREFIX from input info file

    /**Setting the output information file path**/
    sprintf (infofile_out , "%s/%s_uc.info" , moduleoutdir , nameprefix) ;
    LOG (INFO) << "\nOutput information File :" << infofile_out << endl ;

    /**Creating the output information file**/
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the out information file" , infofile_out) ;
    char *ttype[] = {"SignalFileList"} ;
    char *tform[] = {"A100"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table in output information file" , infofile_out) ;
    datainfo.write (finfo_out) ;
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the key value of the NAMEPRFX" , infofile_out) ; //for creating name for output information file

    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    /*----check whether th MODE is IM or PC----*/
    if (datainfo.getModeFlag () == IM)
    {
        readKeywords (infofile_in , 2 , 3 , TSTRING , "SIGDIR" , sigframedir ,
                TINT , "NFILES" , &nframes ,
                TSTRING , "DARKDIR" , darkdir) ;

        //allocating the memory  to two dimension Array 
        sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;

        //reading frame names  columns from the information file.
        fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) sigframelist , NULL , &status) ;
        printError (status , "Error in reading the column value of the Input Signal List" , infofile_in) ;
        LOG (INFO) << endl << "Number of files :" << nframes ;

        updateKeywords (infofile_out , 2 , 2 , TSTRING , "SIGDIR" , sigframedir ,
                TINT , "NFILES" , &nframes) ;

        outframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list

        /*----Calling the Filterbad pix  correction  for the IM mode----*/
        if (unitConvertIM ())
            return (EXIT_FAILURE) ;

        fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , nframes , (void *) outframelist , &status) ;
        printError (status , "Error in Writing the outframe list in output info file" , infofile_out) ;
        freeMemory (sigframelist , nframes , NAMESIZE) ;
        freeMemory (outframelist , nframes , NAMESIZE) ;
    }
    else if (datainfo.getModeFlag () == PC)
    { //incase of the PC mode
        /**Calling the Filterbadpix Correction for the PC mode**/
        if (unitConvertPC ())
            return (EXIT_FAILURE) ;
    }
    else
    { //Invalid mode 
        LOG (ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG (ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }

    /**Closing the both input and output information File**/
    fits_close_file (finfo_in , &status) ;
    printError (status , "Error in closing the input information file" , infofile_in) ;
    writeCommonKeywords (finfo_out , modulename) ; //write DATE, ORIGIN, CREATOR, CHECKSUM
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the out information file" , infofile_in) ;
    return (EXIT_SUCCESS) ;
}


int uvtUnitConversion::unitConvertIM ()
{
    LOG (INFO) << endl << "\nStarted  unit conversion For IM mode.... " ;

    /**Setting the path for the   Signal output Directory**/
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s" , moduleoutdir , sigframedir) ;
    /**Shell command for  creating the output signal frame directory**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the  shell command **/
    system (cmd.c_str ()) ;
    LOG (INFO) << endl << dir << " directory created" << endl ;

    /*------------Dividing all frames by integration time--------------*/
    LOG (INFO) << endl << "\nTotal number of frames - " << nframes ;
    int status = 0 ;
    float *framedata ;
    framedata = new float[xsize * ysize] ;
    double integrationtime = 0.0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char errstr[512] ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    vector<string> vhistorystr ;
    getHistory (vhistorystr) ;
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;

    fitsfile *fptr , *fout ;

    //process for dark frames
    if (darkFrame_Flag)
    { //if Darkframe Flag is set
        LOG (INFO) << "Reading Dark frames........" << endl ;

        status = takeDarkinfo () ;
       
        if (status)
        {
            LOG (INFO) << "***Reading Dark fails***" << endl ;
            return(EXIT_FAILURE);
        }

        t_darkframestart = readDarkFrame (dstartpath , darkFramestart_data) ;
        LOG (INFO) << "The Start  Dark frametime  is - " << t_darkframestart << endl ;

        t_darkframeend = readDarkFrame (dendpath , darkFrameend_data) ;
        LOG (INFO) << "The End  Dark frametime is - " << t_darkframeend << endl ;

        darkCompute_array = new float[xsize * ysize] ;
    }

    /**LOOP for  processing the all the Frames of INPUT Directory**/
    LOG (INFO) << "\nPerforming Unit Conversion..." << endl ;

    for (int i = 0 ; i < nframes ; i ++)
    {
        sprintf (errstr , "Error at iteration number %d" , i) ;
        t_curr = 0 ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[i]) ;
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in Opening the File" , infile) ;

        if (i == 0)
            copyUsrkeywrdsTovect (fptr , key_record) ;

        for (int q = 0 ; q < xsize * ysize ; q ++) framedata[q] = 0.0f ;

        /**Reading the PIXEL array from the input Frame**/
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
        printError (status , "Error in Reading the pixels " , infile) ;

        readKeywords (infile , 1 , 3 , TDOUBLE , "INT_TIME" , &integrationtime ,
                TUSHORT , "FRAMENO" , &frameno ,
                TDOUBLE , "FRMTIME" , &frametime) ;

        t_curr = frametime ;
        if (integrationtime == 0)
        {
            LOG (INFO) << "***Divide by Zero(integration time =0)!!!***" << endl ;
            return (EXIT_FAILURE) ;
        }
        //Dark Computation and Subtraction if the flag is set
        if (darkFrame_Flag)
        {
            status = darkFrameComputation (darkCompute_array) ;
            if (status)
            {
                LOG (INFO) << "***Error in DarkFrame Computation***" << endl ;
                return (EXIT_FAILURE) ;
            }
             
            darkFrameSubtraction (darkCompute_array , framedata) ;
            
        }

        /**Performing unit conversion**/
        for (int pixno = 0 ; pixno < xsize * ysize ; pixno ++)
        {
            if (framedata[pixno] != INVALID_PIX_VALUE)
            {
               framedata[pixno] = framedata[pixno] / integrationtime ;
                
            }
        }

        /**Setting the output path  **/
        sprintf (outfile , "%s/%s/%s_t%.4f_f%d_uc.fits" , moduleoutdir , sigframedir , nameprefix , frametime , frameno) ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output Signal File" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in Creating the image for Signal File" , outfile) ;
        /**Writing the new PIXEL array to the output**/
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata , &status) ;
        printError (status , "***Error in writing the pixels to output***" , outfile) ;


        updateKeywords (outfile , 1 , 1 , TDOUBLE , "INT_TIME" , &integrationtime );
                //TUSHORT , "FRAMENO" , &frameno ,
                //TDOUBLE , "FRMTIME" , &frametime) ;


        /* Copy keyword of  LEVEL-1 file  from vector to output Frame*/
        //writeUsrkeywordsFrmvect (outfile , key_record) ;
        copyUserKeywords (fptr,fout);
        writeCommonKeywords (fout , modulename) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the  output Signal file" , outfile) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input Signal  file" , infile) ;
        /**writing the history**/
        if (history == YES) writeHistory (outfile , vhistorystr) ;
       
        /**Maintaining the output frame list**/
        strcpy (outframelist[i] , basename (outfile)) ; //taking only filename
        cout<< "Total Files written=" << i << "     Remaining files= " << nframes - i <<" \r" ;
    }

    writeUsrkeywordsFrmvect (infofile_out , key_record) ;

    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    vhistorystr.clear () ;

    LOG (INFO) << endl << "Completed unit conversion of frames successfully" << endl ;
    return (EXIT_SUCCESS) ;
}


int uvtUnitConversion::unitConvertPC ()
{
    LOG (INFO) << endl << "\nStated unit Conversion process for PC mode" << endl ;
    char file_in[NAMESIZE] , file_out[NAMESIZE] ; //for input and output event file

    readKeywords (infofile_in , 2 , 1 , TSTRING , "EVTFILE" , eventfile) ;

    /**Setting the path for the input Event File**/
    sprintf (file_in , "%s/%s" , inputdatadir , eventfile) ; //taking event file full path

    /**Setting the path for the output Event File**/
    sprintf (file_out , "%s/%s_uc.events" , moduleoutdir , nameprefix) ;

    fitsfile *fevt_in , *fevt_out , *finfo_out ;
    int status = 0 ;
    double intgtn_Time = datainfo.getIntegrationTime () ;
    /**Opening the input event file**/
    LOG (INFO) << "\nInput Event File  " << file_in << endl ;
    LOG (INFO) << "\noutput Event File " << file_out << endl ;
    int colnum=0;
    long nrows=0;;
    fits_open_file (&fevt_in , file_in , READWRITE , &status) ;
    printError (status , "Error in opening the input file" , file_out) ;

    //copy keywords of level 1 fits file to vector
    copyUsrkeywrdsTovect (fevt_in , key_record) ;

    //creating the output event file 
    fits_create_file (&fevt_out , file_out , &status) ;
    printError (status , "Error in creating the file" , file_out) ;

    //coping input event file to output event file
    fits_copy_file (fevt_in , fevt_out , 1 , 1 , 1 , &status) ;
    printError (status , "Error in coping the file " , file_out) ;
    fits_close_file (fevt_in , &status) ;
    printError (status , "Error in closing the file" , file_out) ;

    fits_movabs_hdu (fevt_out , 2 , NULL , &status) ;
    printError (status , "***Moving To particular HDU fails***" , file_out) ;
     fits_get_num_rows (fevt_out , &nrows , &status) ;
    printError (status , "Error reading the row number of the output file" , file_out) ;
   // char integrationtime_str[50] ; //variable to store integration time
    float *ENP= new float[nrows];
    //sprintf (integrationtime_str , "%.7f" , intgtn_Time) ; //storing integration time in a string
 //   cout<<integrationtime_str<<endl;
 //   string expr = (string) ENPCOLNAME + (string) "/" + (string) integrationtime_str ; //creating expression for calculation  (ENP/integrationtime)
//    LOG (INFO) << "Calculator Expression  " << expr ;
    fits_get_colnum(fevt_out,CASEINSEN,ENPCOLNAME,&colnum,&status);
    printError (status , "Error reading the column  number of the output file" , file_out) ;
    fits_read_col (fevt_out , TFLOAT , colnum , 1 , 1 , nrows , NULL , ENP , NULL , &status) ;
    printError (status , "Error reading Column" , file_out) ;
    for (int i=0;i<nrows;i++)
    {
        ENP[i]=ENP[i]/datainfo.getIntegrationTime ();
    }
      fits_write_col (fevt_out , TFLOAT, colnum , 1 , 1 , nrows , ENP , &status) ;
      printError (status , "Error writing to Column" , file_out) ;
 //   fits_calculator (fevt_out , (char *) expr.c_str () , fevt_out , ENPCOLNAME , "1D" , &status) ; //calculate and update ENPCOLNAME


    vector<string> vhistorystr ;
    if (history == YES)
    {
        LOG (INFO) << "Writing HISTORY to output Event file.." << endl ;
        getHistory (vhistorystr) ; //write history  to the vhistorystr vector
        writeHistory (file_out , vhistorystr) ; //write history if history set to yes
    }

    fits_movabs_hdu (fevt_out , 1 , NULL , &status) ;
    printError (status , "***Moving To particular HDU fails***" , file_out) ;


    writeCommonKeywords (fevt_out , modulename) ; //coping origin,creator,checksum and date to output event file
    fits_close_file (fevt_out , &status) ;
    printError (status , "Error in closing the out event file" , file_in) ;

    //updating  some keywords of output information file.
    updateKeywords (infofile_out , 2 , 2 , TSTRING , "EVTFILE" , basename (file_out) ,
            TSTRING , "IMGFILE" , "") ;

    //writing userkeywords from vector to the output information file
    writeUsrkeywordsFrmvect (infofile_out , key_record) ;

    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    vhistorystr.clear () ;

    return (EXIT_SUCCESS) ;
}


int uvtUnitConversion::getHistory (vector<string> &vhistory)
{
    char *user = getlogin () ;
    int cnt=0;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir =" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir) ;
    if(datainfo.getModeFlag ()==IM)
    {
        char nameDarkFrame[FLEN_FILENAME];
        sprintf(nameDarkFrame,"%d",darkFrame_Flag);
        vhistory.push_back ((string)getSerialNo (cnt)+" Dark frame flag =" + (string)nameDarkFrame) ;   
    }
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


double uvtUnitConversion::readDarkFrame (char * path , float *Array)
{

    fitsfile *fptr ;

    /***Opening the Dark Frame***/
    int status = 0 ;
    long fpixel[2] ;

    /**INITIALIZE with ZERO**/
    for (int i = 0 ; i < xsize * ysize ; i ++) Array[i] = (float) 0.0 ;

    fpixel[0] = fpixel[1] = 1 ;

    /**Opening dark frame */
    fits_open_file (&fptr , path , READONLY , &status) ;
    printError (status , "Error in opening the Dark frame" , path) ;

    /**Reading the  image from the frame**/
    fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , Array , NULL , &status) ;
    printError (status , "Error in reading the pixels from the Dark Frames" , path) ;

    double start_i , start_f , stop_i , stop_f ;
    double frame_time = 0 ;

    /**Reading the  frame time From the frame**/
     readKeywords (path , 1 , 1 , TDOUBLE , "FRMTIME" , &frame_time ) ;


    //frame_time = start_i + start_f ;

    fits_close_file (fptr , &status) ;
    printError (status , "Error in Closing the  Dark Frame" , path) ;

    return frame_time ;

}


int uvtUnitConversion::darkFrameComputation (float *outArray)
{

    float t1 = t_darkframestart ; //time for beginning dark frame
    float t2 = t_darkframeend ; //time for ending dark frame
    float t_curr_temp = t_curr ; //time for current frame

    for (int i = 0 ; i < xsize * ysize ; i ++)
    {
        outArray[i] = 0.0f ;
        float d1 = darkFramestart_data[i] ;
        float d2 = darkFrameend_data[i] ;

        /**Check Divide by ZERO**/
        if ((t_curr_temp - t1) == 0)
        {
            LOG (ERROR) << "***Divide By zero Error***" << endl ;
            return (EXIT_FAILURE) ;
        }
        float d = d1 + ((d2 - d1) *((t_curr_temp - t1) / (t2 - t1))) ; // interpolated value for dark frame pixel

        outArray[i] = (float) d ; //value for dark frame pixel corresponding to current frame pixel

    }

    return (EXIT_SUCCESS) ;
}


void  uvtUnitConversion::darkFrameSubtraction (float *Array , float *frame_data)
{
   
    for (int p = 0 ; p < xsize * ysize ; p ++)
    {

        frame_data[p] = frame_data[p] - Array[p] ;
//        if (frame_data[p] < 0)
//        {    
//           frame_data[p] = 0.0 ;
//        }
    }
}


int uvtUnitConversion::takeDarkinfo ()
{
    double dark_Firstframetime , dark_Secframetime , start_i , start_f , stop_i , stop_f ;
    char darkpath_First[FLEN_FILENAME] ;
    char darkend_Sec[FLEN_FILENAME] ;
    char dark_temp[FLEN_FILENAME] ;

    //vector for storing the dark frame names
    vector<string> dark_framenames ;

    //setting path for dark frame directory  
    sprintf (dark_temp , "%s/%s" , inputdatadir , darkdir) ;

    getFileNamesfrmDir (dark_temp , ".fits" , dark_framenames) ; //get dark file names from dark directory

    if (dark_framenames.size () <2)
    {
        LOG (INFO) << "***No enough Dark frames found at input Directory***" << endl ;
        return (EXIT_FAILURE) ;
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
