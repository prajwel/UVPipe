/* 
 * File:   uvtPixPadding.cpp
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
#include <uvtPixPadding.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<glog/logging.h>
#include<macro_def.h>

uvtPixPadding::uvtPixPadding ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}
//Destructor

uvtPixPadding::~uvtPixPadding () {
 }

//for parameter file reading

int uvtPixPadding::read (int argc , char** argv)
{
    int status = 0 ;
    status = readParams (argc , argv , 5 , FNAME , "inputdatadir" ,  inputdatadir , 
            FNAME , "outdir" , outdir , BOOL , "clobber" , &clobber ,BOOL , "history" , &history , STRING , "mode" ,  mode) ;
    
    return (EXIT_SUCCESS) ;
}

int uvtPixPadding::read (char* inputdatadir , int outputdimension , char* outdir , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;
    this-> padding_dimension = outputdimension ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_FAILURE) ;
}

//Display the Parameter file content 

void uvtPixPadding::display ()
{
    LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT PIX PADDING PARAMETERS      " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input Frame List Directory                        : " << inputdatadir ;
    LOG(INFO) << endl << "Output Directory                               : " << outdir ;
    if (clobber == YES)
        LOG(INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG(INFO) << endl << "Overwrite                                         : NO" ;
    if (history == YES)
        LOG(INFO) << endl << "History                                             : YES" ;
    else
        LOG(INFO) << endl << "History                                              : NO" ;
    LOG(INFO) << endl << "-----------------------------------------------------------------------------------" ;
   }

//PixPadding Process started
int uvtPixPadding::uvtPixPaddingProcess ()
{
    padding_dimension=PIX_PADD_SIZE;

    LOG(INFO) << "\nInside the PixPadding  Module" << endl ;
    /*if the Input Directory of the Filelist is not available then exit from module
       else read Filelist 
    
     
     */
    
    LOG(INFO) << " Input Directory is " << inputdatadir << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG(INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;
    
    //check  existence  of output directory if output  exist and 
    //clobber =YES than remove that directory and recreate it,if exist and clobber =NO than exit from module.
    //if directory is not exist than create it .
    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;

    /*Searching For the .info  file in the Framelist Directory for the Filelist, '.info' file  contains the listing 
    of the files and some usefull info in the header.
     */
    string  tempfilepath = searchFile (inputdatadir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG(INFO) << endl << "***Information file not found in " << inputdatadir << "***" << endl ;
        return (EXIT_FAILURE) ;
    }
    /*setting the path for input information file **/
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    LOG(INFO) << endl << " Input information file :" << infofile_in << endl ;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists (infofile_in)))
    {
        LOG(INFO) << endl << "***Input FileList not Found at Specified PATH,Check Input Directory***" << endl ;
        return (EXIT_FAILURE) ;
    }

    int status = 0 ;
    /* 
   open the .info FITS file  and read the header information from the second HDU. 
     */

    fitsfile *finfo_in , *finfo_out ;
    /*Openig the input information file**/
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input information File" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to particular HDU   " , infofile_in) ;

    /**Reading the some keyword information from the input File**/
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file

    xsize = datainfo.getXsize () ;//xsize-x dimension of the input image
    ysize = datainfo.getYsize () ;//ysize-y dimension of the input image
    
    if(padding_dimension<xsize){
        LOG(ERROR) <<"padding size must be greater than input frame size" ;
        return(EXIT_FAILURE);
    }
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file

    /**setting the path for the output information path**/
    sprintf (infofile_out , "%s/%s_pp.info" , moduleoutdir , nameprefix) ;

    /**Creating the output information file**/
    LOG(INFO)<< endl << " Output information file :" << infofile_out<< endl ;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information File " , infofile_out) ;
    char *ttype[] = {"SignalFrames" , "ExposureFrames"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table" , infofile_out) ;
    
    //opening the output information file.
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU in output information file" , infofile_out) ;

    datainfo.setXsize (padding_dimension) ;//setting xsize to the padding_dimension
    datainfo.setYsize (padding_dimension) ;//setting ysize to  the padding dimention
    /**writing the  keywords to the output information file**/
    datainfo.write (finfo_out) ;
    //fits_close_file (finfo_out , &status) ;
   // in case of IM
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the key value of the NAMEPRFX" , infofile_out) ; //for creating name for output information file
    fits_close_file (finfo_out , &status) ;
   printError (status , "Error in closing the out information  file" , infofile_out) ;
    if (datainfo.getModeFlag () == IM)
    {
               
         readKeywords(infofile_in , 2 , 3 , TSTRING , "SIGDIR" , sigframedir ,
                TSTRING , "EXPDIR" , expframedir ,
                TINT , "NFILES" , &nframes 
                ) ;
         
        //reading frame names from information file into vector
        //sigframlist-array for storing the input signal frame names from input information file.
         //expframlist-array for storing the input exposure frame names from input information file.
        sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        expoframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        /**reading signal frame names and exposure frame names  from input information file to the sigframelist and exposure frame names */
        status=readColumnsFromFITS (infofile_in,2,2,TSTRING,1,sigframelist,nframes,TSTRING ,2,expoframelist,nframes);
        if(status){
        LOG(INFO)<<"Error reading  the columns from the file"<<endl;
        return(EXIT_FAILURE);
        }
        
        //Performing Pixel Padding for IM
        if (pixPaddingIM ()) return (EXIT_FAILURE) ;

        //updating the key value of SIGDIR ,EXPDIR,nframes to putput information file.
        /**
         * SIGDIR-output signal frame directory name
         * EXPDIR-output Exposure frame directory name
         * NFILES- total number of frames at output 
         * @return 
         */
          status=updateKeywords(infofile_out , 2 , 3 , TSTRING , "SIGDIR" , sigframedir ,
                TSTRING , "EXPDIR" , expframedir ,
                TINT , "NFILES" , &nframes 
                ) ;
          if(status) printError (status,"keyword not found in header",infofile_out);
        
        freeMemory (sigframelist , nframes , NAMESIZE) ;
        freeMemory (expoframelist , nframes , NAMESIZE) ;
    } //in case of PC
    else if (datainfo.getModeFlag () == PC)
    {
        fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;//updating the key value of EVENT file for future module reference
        printError (status , "Error in reading the key value of the EVTFILE" , infofile_in) ;

        //Performing pixel padding
        if (pixPaddingPC ()) return (EXIT_FAILURE) ;

    } //in case of neither PC or IM
    else
    {
        LOG(INFO) << endl << "Invalid input for operating mode parameter" ;
        LOG(INFO) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
   fits_close_file(finfo_in,&status);
   printError (status , "Error in closing  the  information file" , infofile_in) ;
    LOG(INFO) << endl << "Pixel Padding Process Completed Successfully" << endl ;

    return (EXIT_SUCCESS) ;
}

int uvtPixPadding::pixPaddingIM ()
{

    LOG(INFO) << endl << "Started Pixel Padding  for IM mode"<<endl ;

    /**allocating the memory to the output Frame list**/
    /**
     * outputsigframelist-array for storing the output signal frame names
     * outexpframelist-array for storing the output exposure frame names
     * @return 
     */
    
    char **outsigframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    char **outexpframelist = allocateMemory<char>(nframes , NAMESIZE) ; //to store output exposure frame list  
    //creating signal directory
    char dir[FLEN_FILENAME] ;
    /**setting the path for the output Signal Frame Directory**/
    sprintf (dir , "%s%s" , moduleoutdir , sigframedir) ;
    /**Shell command for setting the output Signal Frame path**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    //creating exposure directory

    /**setting the path for the output Exposure Frame Directory**/
    sprintf (dir , "%s%s" , moduleoutdir , expframedir) ;
    /**Shell command for setting the output Exposure Frame path**/
    cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;

    LOG(INFO) << endl << "Total number of frames - " << nframes ;
    int status = 0 ;
    //frame_data-storing the input image pixels
    float frame_data[xsize * ysize] ;
    //frameoutsigdata-array for storing the output signal frame data
    //frameoutexpdate-array for storing the output exposure frame data
    float *frameoutsigdata = new float[padding_dimension * padding_dimension] ;
    float *frameoutexpdata = new float[padding_dimension * padding_dimension] ;
    double integrationtime = 0.0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char errstr[512] ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = padding_dimension ;
    vector<string> vhistorystr ;
    if (history == YES)
    {
       getHistory (vhistorystr) ;//storing history to vhistorystr vector 
    }
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;
    /**LOOP for processing the number of frames **/
    LOG(INFO)<<"\nApplying Pixel Padding..."<<endl;
    //loop for number of frames 
    for (int i = 0 ; i < nframes ; i++)
    {
        sprintf (errstr , "Error at iteration number %d" , i) ;
        fitsfile *fptr , *fout ;
        //setting the path for input signal  frame 
        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[i]) ;
        //opening the input signal  frame
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input Signal Files" , infile) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , frame_data , NULL , &status) ;
        printError (status , "Error in reading pixels from input frame" , infile) ;
        
       status=readKeywords (infile , 1 , 3 , TDOUBLE , "INT_TIME" , &integrationtime ,
                TUSHORT , "FRAMENO" , &frameno ,
                TDOUBLE , "FRMTIME" , &frametime ) ;
       if(status)   printError (status , "keyword not found in header" , infile) ;
     
       if(i==0)
       {
            copyUsrkeywrdsTovect (fptr,key_record);//copying the  level-1 keywords to the   vector
       }
        for(int p=0;p<padding_dimension*padding_dimension;p++)
            frameoutsigdata[p]=-9999;
          /*
         padding is applied for converting the 512*512 frame to 600*600  frame
         */
        status = Applypadding (frame_data , xsize , ysize , frameoutsigdata , padding_dimension , padding_dimension) ;
        if (status)
        {
            LOG(ERROR) << "Padding failed for SignalFrame " << endl ;
            return (EXIT_FAILURE) ;
        }
        //signal frames
        /**Setting the path for the output Signal Frame**/
        sprintf (outfile , "%s/%s/%s_t%f_f%d_sig_pp.fits" , moduleoutdir , sigframedir , nameprefix , frametime , frameno) ;
        /**creating the output Signal Frame**/
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output signal file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the img" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , padding_dimension*padding_dimension , frameoutsigdata , &status) ;
        printError (status , "Error in writing the pixels to output signal file" , outfile) ;
       
        copyUserKeywords (fptr , fout) ;//copying the level-1 keywords to the output signal frame from input signal  frame
        if (history == YES)writeHistory (outfile , vhistorystr) ;
        writeCommonKeywords (fout , modulename) ;
       
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output signal file" , outfile) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input signal file" , outfile) ;
      
        strcpy (outsigframelist[i] , basename (outfile)) ; //taking only filename

        //for exposure frames
        /**setting the path for the input Exposure Frame **/
        sprintf (infile , "%s/%s/%s" , inputdatadir , expframedir , expoframelist[i]) ;
        /**Reading the input Exposure Frame**/
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input Exposure file" , infile) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , frame_data , NULL , &status) ;
        printError (status , "Error in reading the pixels from the input exposure file" , infile) ;
        
        status=readKeywords (infile , 1 , 3 , TDOUBLE , "INT_TIME" , &integrationtime ,
                TUSHORT , "FRAMENO" , &frameno ,
                TDOUBLE , "FRMTIME" , &frametime ) ;
        if(status)   printError (status , "keyword not found in header" , infile) ;

        /*padding is applied on Exposure frame */
         for(int p=0;p<padding_dimension*padding_dimension;p++)
            frameoutexpdata[p]=-9999;
       
        /*
         padding is applied for converting the 512*512 frame to 600*600  frame
         */
        status = Applypadding (frame_data , xsize , ysize , frameoutexpdata , padding_dimension , padding_dimension) ;
        if (status)
        {
            LOG(ERROR) << "Padding failed For ExposureFrame " << endl ;
            return (EXIT_FAILURE) ;
        }
        /**setting the path for the output Exposure frame **/
        sprintf (outfile , "%s/%s/%s_t%f_f%d_exp_pp.fits" , moduleoutdir , expframedir , nameprefix , frametime , frameno) ;
        /**Creating the output Exposure Frame**/
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output  Exposure file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the image" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , padding_dimension*padding_dimension , frameoutexpdata , &status) ;
        printError (status , "Error in writing  the pixels from output Exposure file" , outfile) ;
        
        copyUserKeywords (fptr , fout) ;//copying the level-1 keywords to the output exposure  frame from input  exposure frame
        if(history==YES)    writeHistory (outfile , vhistorystr) ;
        writeCommonKeywords (fout , modulename) ;
        
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output Exposure File" , outfile) ;
        strcpy (outexpframelist[i] , basename (outfile)) ; //taking only filename
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing  the input Exposure file" , outfile) ;
        
       
           // system (cmd.c_str ());
        cout<< "Total Files written=" << i << "     Remaining files= " << nframes - i <<" \r" ;
         
    } //end of loop for number of frames
    
    /*----Updating the information to the output Information File----*/
    /**opening the output information file which has been generated earlier and update information**/
    fitsfile *finfo_out ;
    LOG(INFO)<<"\nUpdating the  output information file .."<<endl;
    status=writeColumnsToFITS (infofile_out,2, 2,TSTRING,1,outsigframelist,nframes,TSTRING ,2,outexpframelist,nframes);
        if(status)
        {
        LOG(INFO)<<"Error writing  the columns to the fits file"<<endl;
        return(EXIT_FAILURE);
        }

    /**Write user keywords and  history to the vector**/
    writeUsrkeywordsFrmvect (infofile_out,key_record);
    
      if(history==YES)
        writeHistory (infofile_out , vhistorystr) ;
    /**memory release **/
    freeMemory (outsigframelist , nframes , NAMESIZE) ;
    freeMemory (outexpframelist , nframes , NAMESIZE) ;
    vhistorystr.clear () ;
    return (EXIT_SUCCESS) ;
}

int uvtPixPadding::pixPaddingPC ()
{
    LOG(INFO) << endl << "Started Pixel padding for PC mode\n" ;
    fitsfile *fptr , *fout , *finfo_out ;
    int status = 0 ;
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    char eventfilename[FLEN_FILENAME] ;
    unsigned short *xi , *yi ;
     vector<string> vhistorystr ;//vector  for storing the history
   
    /**Setting the path for the output EVENT file**/
    sprintf (eventfilename , "%s_pp.events" , nameprefix) ;
    sprintf (outfile , "%s/%s" , moduleoutdir , eventfilename) ;
   
    /**creating the output File**/
    fits_create_file (&fout , outfile , &status) ;
    printError (status , "Error in creating  the output event file" , outfile) ;
    long nrows = 0 ;
  
    /**Setting the path of input EVENT file**/
    sprintf (infile , "%s/%s" , inputdatadir , eventfile) ;
    LOG(INFO) << "\nInput Event file is  " << infile << endl ;
   
    /**Open the  input EVENT File **/
    LOG(INFO)<<"\nInput Event file :"<<infile<<endl;
    LOG(INFO)<<"\nOutput Event file :"<<outfile<<endl;
    fits_open_file (&fptr , infile , READONLY , &status) ;
    printError (status , "Error in opening the input event file " , infile) ;
  
    copyUsrkeywrdsTovect (fptr,key_record);
    /*copy the input event file to the output event file*/
    fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
    printError (status , "Error in coping the event file " , outfile) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the input event file" , outfile) ;
   
    fits_movabs_hdu (fout , 2 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU in input file" , infile) ;
    fits_get_num_rows (fout , &nrows , &status) ;
    printError (status , "Error in getting the number of rows in input Event file" , infile) ;
    LOG(INFO) << " Total Number of Events are " << nrows << endl ;
 
    //xi-array for storing the xi column from the copied output event file
    //yi array for storing the yi column from the copied output event file
    xi = new unsigned short[nrows] ;
    yi = new unsigned short[nrows] ;
    
    
      fits_read_col (fout , TUSHORT , 4 , 1 , 1 , nrows , NULL , xi , NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
         fits_read_col (fout , TUSHORT , 6 , 1 , 1 , nrows , NULL , yi , NULL , &status) ;
        printError (status , "Error in reading the column y" , outfile) ;
        
//    status=readColumnsFromFITS (outfile,2,2,TUSHORT,4,(unsigned short*)xi,(int)nrows,TUSHORT ,6,(unsigned short*)yi,(int)nrows);
//        if(status){
//        LOG(INFO)<<"Error reading  the columns from the file"<<endl;
//        return(EXIT_FAILURE);
//        }
         
    /**PIXEL padding apply here**/
    //xpadding-x padding to be done
    //ypadding-y padding to be done  
    int xpadding = (padding_dimension - xsize) / 2 ;
    int ypadding = (padding_dimension - ysize) / 2 ;
    LOG(INFO)<<"\nApplying  Pixel  Padding...."<<endl;
    for (int j = 0 ; j < nrows ; j++)
    {
        xi[j] = xi[j] + xpadding ;
        yi[j] = yi[j] + ypadding ;
    }
    
    /**writing the columns  to the output EVENT file**/
   status=writeColumnsToFITS (outfile,2,2,TUSHORT,4,xi,nrows,TUSHORT ,6,yi,nrows);
        if(status)
        {
        LOG(INFO)<<"Error writing  the columns to  the fits file"<<endl;
        return(EXIT_FAILURE);
        }
    
    //updating output information file
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU" , infofile_out) ;
    fits_update_key (finfo_out , TSTRING , "EVTFILE" , basename (outfile) , NULL , &status) ;
    printError (status , "Error in updating the key  value of the EVTFILE" , infofile_out) ;
   
    //writing history to the output event file
    if (history == YES)
    {
        getHistory (vhistorystr) ;
        LOG(INFO)<<"\nWriting  History.."<<endl;
        writeHistory (outfile , vhistorystr) ;
    }
   
    //write origin,creator ,checksum and date to the output event file
    writeCommonKeywords (fout , modulename) ;
    fits_close_file (fout , &status) ;
    printError (status , "Error in closing the output  Event file" , outfile) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file" , outfile) ;
    //writing level1  keywords  and history from vector to the output information file.
    writeUsrkeywordsFrmvect (infofile_out,key_record);
      if(history==YES)
        writeHistory (infofile_out , vhistorystr) ;
    return (EXIT_SUCCESS) ;
}

int uvtPixPadding::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    char str_dim[10] ;
    sprintf (str_dim , "%d" , padding_dimension) ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"inputdatadir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"output Dimension " + (string) str_dim) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"Module Output directory=" + (string) moduleoutdir) ;
    if (clobber == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+"clobber=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+"clobber=no") ;
    if (history == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+"history=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+"history=no") ;
    vhistory.push_back ("Parameter List END") ;
    return (EXIT_SUCCESS) ;
}
