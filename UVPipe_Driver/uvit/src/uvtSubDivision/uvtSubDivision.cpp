/* 
 * File:   uvtSubDivision.cpp
 * Authors:: Dhruv, Sanjay K Singh, Arvind K Singh
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
#include <uvtSubDivision.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<glog/logging.h>
fitsfile *fitin ;
int status = 0 ;
//Constructor:called When the object is  created

uvtSubDivision::uvtSubDivision ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}

uvtSubDivision::~uvtSubDivision () { }

int
uvtSubDivision::read (int argc , char** argv)
{
   int status = 0 ;
   status=readParams (argc , argv , 6 , FNAME , "inputdatadir" , inputdatadir,INT,"subdivisionsize",&subdivision_dimension,FNAME,"outdir",outdir,BOOL,"clobber",&clobber,
            BOOL,"history",&history,STRING ,"mode",&mode) ;
   
   if (status) return (EXIT_FAILURE) ;

    return (EXIT_SUCCESS) ;
}

int uvtSubDivision::read (char* inputdatadir , int outputdimension , char* outdir , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;
    this->subdivision_dimension = outputdimension ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_FAILURE) ;
}

void uvtSubDivision::display ()
{
       LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT SUB DIVISION PARAMETERS      " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input Frame List Directory                        : " << inputdatadir ;

    LOG(INFO) << endl << "Output Directory                               : " << outdir ;
    LOG(INFO) << endl << "Subdivision dimension                               : " << subdivision_dimension ;
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

int uvtSubDivision::uvtSubDivisionProcess ()
{
    /*if the Input Directory of the Filelist is not available then exit from module
     else read Filelist 
    
     */
    LOG(INFO) << "SubDivision Process Started\n" << endl ;
    /*----setting path for the output Directory Structure----*/
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG(INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;
    
    //check  existence  of output directory if output  exist and 
    //clobber =YES than remove that directory and recreate it,if exist and clobber =NO than exit from module.
    //if directory is not exist than create it .
    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;


    /**Searching the information  file  from the input Directory**/
    string  tempfilepath = searchFile (inputdatadir , "info") ; //searching the .info  file in the perticuler directory
    if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }
    /**Setting the  path for  input information File**/
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;

    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists (infofile_in)))
    {
        LOG(INFO) << endl << "***Input FileList not Found at Specified PATH,Check INPUT DIRECTORY***" ;
        return (EXIT_FAILURE) ;
    }
    LOG(INFO) << endl << " Input information File :" << infofile_in ;
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    /**Open input information File**/
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the inputInformation File" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU" , infofile_in) ;

    /**Reading the  keyword information **/
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;//xsize - x dimension of the image
    ysize = datainfo.getYsize () ;//ysize - y dimension of the image
    if(subdivision_dimension<xsize){
        LOG(ERROR)<<"subdivision dimensions must be greater than input image dimensions";
        return(EXIT_FAILURE);
    }
    if (xsize <= 0 || ysize <= 0)
    {
        LOG(INFO) << endl << "Invalid xsize or ysize" << endl ;
        return (EXIT_FAILURE) ;
    }
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "NAMEPRFX  keyword not Found" , infofile_in) ; //for creating name for output information file
    /**Creating the output information PATH**/
    sprintf (infofile_out , "%s/%s_sd.info" , moduleoutdir , nameprefix) ;
      LOG(INFO) << "\nOutput information File " << infofile_out<< endl ;
    /**creating the output information File**/
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information File" , infofile_out) ;
    char *ttype[] = {"SignalFrames" , "ExposureFrames"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table in output information file" , infofile_out) ;

    datainfo.setXsize (subdivision_dimension) ;//setting the output image x-dimension  to 4800
    datainfo.setYsize (subdivision_dimension) ;//setting the output image y-dimension  to 4800
    /**Writing the Keyword information to the information File**/
    datainfo.write (finfo_out) ;
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "NAMEPRFX keyword not updated/Found" , infofile_out) ; //for creating name for output information file
   
    //incase of IM
    
    if (datainfo.getModeFlag () == IM)
    {
        /*reading key value SIGDIR,EXPDIR and NFILES from the input information file
         *SIGDIR-directory path for reading the signal frames
         * EXPDIR-directory path for reading the exposure frames 
         * NFILES-number of frames  to be read
         */
   
          readKeywords (infofile_in , 2, 3 , TSTRING , "SIGDIR" , sigframedir ,
                TSTRING , "EXPDIR" , expframedir ,
                TINT , "NFILES" , &nframes) ;        
        
        LOG(INFO) << endl << "Number of files :" << nframes ;
       
        //reading frame names from information file into array
        /*
         *sigframelist=array for storing the input signal frame's names
         * expoframelist= array for storing the input exposure frame's  names
         */
        
        sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        expoframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        
        // reading the signal frame and exposure frame's name to arrays.
        status=readColumnsFromFITS (infofile_in,2,2,TSTRING,1,sigframelist,nframes,TSTRING ,2,expoframelist,nframes);
        if(status){
        LOG(INFO)<<"Error reading  the columns from the file"<<endl;
        return(EXIT_FAILURE);
        }
         
         //method for IM mode execution
        if (subDivisionIM ()) return (EXIT_FAILURE) ;
      
        updateKeywords (infofile_out , 2 , 3 , TSTRING , "SIGDIR" , sigframedir ,TSTRING,"EXPDIR",expframedir,TINT,"NFILES",&nframes);//updating the keywords
        writeCommonKeywords (finfo_out , modulename) ;
        freeMemory (sigframelist , nframes , NAMESIZE) ;
    }//incase of PC
    else if (datainfo.getModeFlag () == PC)
    {
        fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
        printError (status , "Error reading the key value of the EVTFILE" , infofile_in) ;
        //method for IM mode execution
        if (subDivisionPC ()) return (EXIT_FAILURE) ;
   }
    else//incase of not IM or PC mode
    {
        LOG(ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG(ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
     fits_close_file (finfo_out , &status) ;
    return (EXIT_SUCCESS) ;
}

int uvtSubDivision::subDivisionIM ()
{
   LOG(INFO) << endl << "\nStarted Subdivision process  for IM mode.....\n" ;
    /**memory allocation For output frame listing**/
    char **outsigframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    char **outexpframelist = allocateMemory<char>(nframes , NAMESIZE) ; //to store output exposure frame list  
    //creating signal directory
    char dir[FLEN_FILENAME] ;
    /**Setting the output signal directory path**/
    sprintf (dir , "%s%s" , moduleoutdir , sigframedir) ;
    /**Shell command for the directory creation**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    //creating exposure directory
    /**setting the  path for the Exposure frame Directory**/
    sprintf (dir , "%s%s" , moduleoutdir , expframedir) ;
    /**Shell command for the  creating the Directory**/
    cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    
    LOG(INFO) << endl << "Total number of frames - " << nframes ;
    int status = 0 ;
    float *framedata = new float[xsize * ysize] ;
    float *frameoutsigdata = new float[subdivision_dimension * subdivision_dimension] ;
    float *frameoutexpdata = new float[subdivision_dimension * subdivision_dimension] ;
    double integrationtime = 0.0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char errstr[512] ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = subdivision_dimension ;

    vector<string> vhistorystr ;//vector  for storing the history
    if(history==YES)  getHistory (vhistorystr) ;//reading history and store it in the vector
        
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;
    //factor for subdivision.
    int sub_division_factor = subdivision_dimension / xsize ;
    int division_factor = sub_division_factor*sub_division_factor ;
    if (division_factor == 0)
    {
        LOG(ERROR) << "***Divide by Zero***" << endl ;
        return (EXIT_FAILURE) ;
    }
    LOG(INFO)<<"Applying Subdivision..."<<endl;
    LOG(INFO) << endl << "Processing Started \nPlease wait....\n" ;
    //cmd = "clear" ;
    /**Loop for processing  NUMBER of FRAMES **/
    for (int i = 0 ; i < nframes ; i++)
    {
        if (i % 50 == 0 && i > 0)
        {
           // system (cmd.c_str ());
          cout<< "Total Files written=" << i << "     Remaining files= " << nframes - i <<" \r" ;
         }
        sprintf (errstr , "Error at iteration number %d" , i) ;
        fitsfile *fptr , *fout ;

        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[i]) ;
        //opening the input  signal frame and reading image
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input  Signal  file" , infile) ;
        
        if(i==0) copyUsrkeywrdsTovect (fptr,key_records); //coping the keywords from fits file to vector
        
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
        printError (status , "Error in reading the Pixels from input Signal file" , infile) ;
        
        readKeywords (infile , 1, 3 , TDOUBLE , "INT_TIME" , &integrationtime ,
                TUSHORT , "FRAMENO" ,& frameno ,
                TDOUBLE , "FRMTIME" , &frametime) ;   

        for (int p = 0 ; p < xsize * ysize ; p++)
        {
            if (framedata[p] != INVALID_PIX_VALUE)//identify whether pixel is declared invalid by previous module or not.if yes then don't do anything to it
                framedata[p] = framedata[p] / division_factor ;//subdivide the pixels
        }

        /*method For Applying SubDivision,
         * method takes the input Array as a First Argument ,Apply  subDivision 
         * on it,and update output array passed as a Second Array 
         *  */

        int check = applySubDivision (framedata , frameoutsigdata) ;
        if (check)
        {
            LOG(ERROR) << "SubDivision Failed For SignalFrames" ;
            return (EXIT_FAILURE) ;
        }

        /** setting the path for output  signal frame **/
        sprintf (outfile , "%s/%s/%s_t%f_f%d_sig_sd.fits" , moduleoutdir , sigframedir , nameprefix , frametime , frameno) ;
        /**Creating the output Signal Frame  **/
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output Signal file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
         printError (status , "Error in creating the image " , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , subdivision_dimension*subdivision_dimension , frameoutsigdata , &status) ;
        printError (status , "Error Writing the pixels in output signal frame list " , outfile) ;
        fits_update_key (fout , TUSHORT , "FRAMENO" , &frameno , "File name prefix" , &status) ;
        printError (status , "NAMEPRFX keyword not updated/Found" , outfile) ; //for creating name for output information file
     
        /**Copy user Keyword from input frame to the output Frame **/
        copyUserKeywords (fptr , fout) ;
         if(history==YES) writeHistory (outfile , vhistorystr) ;
        /**write common keywords**/
       
        writeCommonKeywords (fout , modulename) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output Signal File" , outfile) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input signal File" , outfile) ;
      
        strcpy (outsigframelist[i] , basename (outfile)) ; //taking only filename
        /**Setting the path for the  input Exposure frame **/
        sprintf (infile , "%s/%s/%s" , inputdatadir , expframedir , expoframelist[i]) ;
        /**Opening the input Exposure Frame and reading image**/
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input  Exposure file" , outfile) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
        printError (status , "Error in reading the pixels from the input  Exposure file" , outfile) ;
      
        readKeywords (infile , 1, 3 , TDOUBLE , "INT_TIME" , &integrationtime ,
                TUSHORT , "FRAMENO" ,& frameno ,
                TDOUBLE , "FRMTIME" , &frametime) ;   


        check = applySubDivision (framedata , frameoutexpdata) ;
        if (check)
        {
            LOG(ERROR) << "SubDivision Failed for Exposure Frames" ;
            return (EXIT_FAILURE) ;
        }

        /**setting the path for the output Exposure frame **/
        sprintf (outfile , "%s/%s/%s_t%f_f%d_exp_sd.fits" , moduleoutdir , expframedir , nameprefix , frametime , frameno) ;
      /**creating the output exposure frame**/
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output Exposure File" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        fits_write_pix (fout , TFLOAT , fpixel , subdivision_dimension*subdivision_dimension , frameoutexpdata , &status) ;
        printError (status , "Error in creating the output Exposure File" , outfile) ;

        /**Copy user Keyword from input frame to the output Frame **/
        copyUserKeywords (fptr , fout) ;
       
        if(history==YES)writeHistory (outfile , vhistorystr) ;
        /**write common keywords**/
        writeCommonKeywords (fout , modulename) ;
        /**Closing the output EVENT file**/
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output Exposure file" , outfile) ;
       
        strcpy (outexpframelist[i] , basename (outfile)) ; //taking only filename
        /**Closing the input EVENT file**/
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input Exposure file" , infile) ;
    } //end of loop for number of frames
    
    //writing history and level1 keywords to output information file 
     writeUsrkeywordsFrmvect (infofile_out,key_records);
     if(history==YES)
                writeHistory (infofile_out , vhistorystr) ;
    
    vhistorystr.clear () ;
    /**open output information File which has been generated earlier  and update the information in it**/
    LOG(INFO)<<"\nUpdating output information file..."<<endl;
  
    status=writeColumnsToFITS (infofile_out,2,2,TSTRING,1,outsigframelist,nframes,TSTRING ,2,outexpframelist,nframes);
    if(status){
    LOG(INFO)<<"Error writing  the columns to  the fits  file"<<endl;
        return(EXIT_FAILURE);
    }
    /**memory release of the  output Signal frame list**/
    freeMemory (outsigframelist , nframes , NAMESIZE) ;
    freeMemory (outexpframelist , nframes , NAMESIZE) ;

    return (EXIT_SUCCESS) ;

}

int uvtSubDivision::subDivisionPC ()
{
    LOG(INFO) << endl << "\nStarted Subdivision process  for PC mode\n" ;
    fitsfile *fptr , *fout , *finfo_out ;
    int status = 0 ;
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    char eventfilename[FLEN_FILENAME] ;
    /**setting the path for the output EVENT file**/
    sprintf (eventfilename , "%s_sd.events" , nameprefix) ;
    sprintf (outfile , "%s/%s" , moduleoutdir , eventfilename) ;
    /**creating the output EVENT file**/
    fits_create_file (&fout , outfile , &status) ;
    printError (status , "Error in creating the output event file" , outfile) ;
    long nrows = 0 ;
    /**setting the path for the  input EVENT file**/
    sprintf (infile , "%s/%s" , inputdatadir , eventfile) ;
    LOG(INFO) << "\nInput Event file is " << infile << endl ;
      LOG(INFO) << "\nOutput Event file is " << outfile << endl ;
    fits_open_file (&fptr , infile , READONLY , &status) ;
    printError (status , "Error in opening the input event file " , infile) ;
    unsigned short *xi , *yi ;
    float *xf , *yf ;
    //copy input event file's keywords to the  vector
    copyUsrkeywrdsTovect (fptr,key_records);
    /*Copy the input Event file to output Event File*/
    fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
    printError (status , "Error in Coping the input Event file to the output Event file" , infile) ;
    fits_movabs_hdu (fout , 2 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU in input file" , infile) ;
    fits_get_num_rows (fout , &nrows , &status) ;
    printError (status , "Error in getting the number of rows in input Event file" , infile) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the  input  Event file") ;
    LOG(INFO) << "\nTotal number of EVENTS are " << nrows << endl ;
    
    //arrays for storing values of the input event file column values. 
    xi = new unsigned short[nrows] , yi = new unsigned short[nrows] ;
    xf = new float[nrows] , yf = new float[nrows] ;
    
        fits_read_col (fout , TUSHORT, 4 , 1 , 1 , nrows , NULL , xi, NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
         fits_read_col (fout , TUSHORT , 6 , 1 , 1 , nrows , NULL , yi, NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
         fits_read_col (fout , TFLOAT , 5 , 1 , 1 , nrows , NULL , xf, NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
         fits_read_col (fout , TFLOAT , 7 , 1 , 1 , nrows , NULL , yf, NULL , &status) ;
        printError (status , "Error in reading the column x" , outfile) ;
    
//    status=readColumnsFromFITS (outfile,2,2,TUSHORT,4,xi,nrows,TUSHORT ,6,yi,nrows,TFLOAT,5,xf,nrows,TFLOAT,7,yf,nrows);
//    if(status)
//    {
//    LOG(INFO)<<"Error reading  the columns from the file"<<endl;
//    return(EXIT_FAILURE);
//    }
    
    float temp_x = 0.0f ;
    float temp_y = 0.0f ;
    float multi_factor = subdivision_dimension / xsize ;
    
    LOG(INFO)<<"Performing SubDivision ... "<<endl;
    /**multiplying  x and y columns with multiplication factor to convert it into subdivision dimention**/
    LOG(INFO)<< "\nMultiplying each event (x,y) with "<< multi_factor << " to convert " << xsize << "*"<< ysize << "  data to " << subdivision_dimension << "*" << subdivision_dimension << " data" << endl ;
  //combining integer and fraction part of events and  multiply it with multi_factor & again divide it into integer and fraction part
    for (int j = 0 ; j < nrows ; j++)
    {
        temp_x = multi_factor * (xi[j] + xf[j]) ;
        temp_y = multi_factor * (yi[j] + yf[j]) ;
        xf[j] = temp_x - (int) (temp_x) ;
        yf[j] = temp_y - (int) temp_y ;
        xi[j] = (int) temp_x ;
        yi[j] = (int) temp_y ;
    }
    /**updating the columns of X and Y  to the output EVENT file**/    
    status=writeColumnsToFITS (outfile,2,4,TUSHORT,4,xi,nrows,TUSHORT ,6,yi,nrows,TFLOAT,5,xf,nrows,TFLOAT,7,yf,nrows);
     if(status){
    LOG(INFO)<<"Error writing  the columns to  the fits  file"<<endl;
        return(EXIT_FAILURE);
    }
    
    delete[] xi , yi , xf , yf ;    
    
    LOG(INFO)<<"Update output information file "<<endl;
    /**open output information File  which has been generated earlier and update some of the information to it**/
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in Moving to 2nd HDU in Event file" , infofile_out) ;
    fits_update_key (finfo_out , TSTRING , "EVTFILE" , basename (outfile) , NULL , &status) ;
    printError (status , "Error in updating the key  value of the EVTFILE" , infofile_out) ;
  
    vector<string> vhistorystr ;
    if (history == YES)
    {        
        getHistory (vhistorystr) ;//read history and store it in vhistorystr vector.
        LOG(INFO)<<"Writing History.."<<endl;
        writeHistory (outfile , vhistorystr) ;//write history to output event file
    }
    //write origin ,creator ,checksum and date to output event file
    writeCommonKeywords (fout , modulename) ;
    fits_close_file (fout , &status) ;
    printError (status , "Error in closing the output  Event file" , outfile) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file" , outfile) ;
     
   //write user keywords to the output information file 
     writeUsrkeywordsFrmvect (infofile_out,key_records);
     if(history==YES)
                writeHistory (infofile_out , vhistorystr) ;
    return (EXIT_SUCCESS) ;
}

//method for Apply SubDivision
int uvtSubDivision::applySubDivision (float *inputArray , float *outputArray) //size of inputArray is (xsize*ysize)                                                                                                                      
{ //size of outputArray is (subdivision_dimension*subdivision_dimension)                                         
    if (xsize == 0 || ysize == 0)
    {
        LOG(ERROR) << "Divide by Zero" << endl ;
        return (EXIT_FAILURE) ;
    }
    int xfactor = subdivision_dimension / xsize ;
    int yfactor = subdivision_dimension / ysize ;
    for (int i = 0 ; i <subdivision_dimension ; i++)
    {
        for (int j = 0 ; j <subdivision_dimension ; j++)
        {
            outputArray[i * subdivision_dimension + j] = inputArray[(i / yfactor) * ysize + (j / xfactor)] ;
        }
    }
    return (EXIT_SUCCESS) ;
}

int uvtSubDivision::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char str_subdivsize[50] ;
    sprintf (str_subdivsize , "%d" , subdivision_dimension) ;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Subdivision dimension =" + (string) str_subdivsize) ;
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

