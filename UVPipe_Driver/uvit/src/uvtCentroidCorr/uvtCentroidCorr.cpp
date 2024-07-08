/* 
 * File:   uvtCentroidCorr.cpp
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
//#include "uvtCentroid.h"
#include "uvtCentroidCorr.h"
#include<pthread.h>
#include<glog/logging.h>

#define DEFAULT_FRAMESIZE 512
#define DARKFRAME_SIZE 512 

uvtCentroidCorr::uvtCentroidCorr ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}

uvtCentroidCorr::~uvtCentroidCorr ()
{
   // delete[] darkFrameend_data , darkFramestart_data ;
}
int uvtCentroidCorr::read (int argc , char** argv)
{
    int status = 0 ;
      
    status = readParams (argc , argv , 7 , FNAME , "inputdatadir" , inputdatadir,FNAME,"caldbDir",caldbDir,FNAME,"outdir",outdir,FNAME,"dataIngestdir",dataIngestDir,BOOL , "clobber" , &clobber,BOOL,"history",&history,STRING,"mode",mode) ;
    if (status) return (EXIT_FAILURE) ;

    return (EXIT_SUCCESS) ;
}

int uvtCentroidCorr::read (char *inputdatadir , char *caldbDir , char *inputdataIngestDir,char *outdir , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->caldbDir , caldbDir) ;
    strcpy (this->outdir , outdir) ;
     strcpy (this->dataIngestDir , inputdataIngestDir) ;
    //   this->threshold_multph=thr;
    this->clobber = clobber ;
    this->history = history ;
    // this->centroid_algo = cent_algo;
    return (EXIT_SUCCESS) ;
}

void uvtCentroidCorr::display ()
{
    LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT CENTROID CORRECTION PARAMETERS              " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input data directory  : " << inputdatadir ;
    LOG(INFO) << endl << "CALDB directory : " << caldbDir ;
    LOG(INFO) << endl << "Output Directory : " << outdir ;
    LOG(INFO) << endl << "Centroid Algorithm used : " << centroid_algo ;
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

int uvtCentroidCorr::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Centroid Algorithm=" + (string) centroid_algo) ;
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

int uvtCentroidCorr::uvtCentroidCorrProcess (){
    LOG(INFO) << endl << "Started Centroid Correction  process.." << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG(INFO) << endl << "\nModule Output Directory : " << moduleoutdir << endl ;

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
    cmd = "mkdir -p " + (string) moduleoutdir ;
    system (cmd.c_str ()) ; // creating output directory to keep output from unitConversion
    LOG(INFO) << endl << moduleoutdir << "  directory created" ;

    
    char dataIngestInfo[PIL_LINESIZE];   
    string tempfilepath = searchFile (dataIngestDir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG(ERROR) << endl << "Error in finding info file" ;
        return (EXIT_FAILURE) ;
    }
    int status = 0 ;
    fitsfile *fdi;
    sprintf (dataIngestInfo , "%s/%s" , dataIngestDir , tempfilepath.c_str()) ;
    fits_open_file (&fdi , dataIngestInfo , READONLY , &status) ;
    printError (status , "Error Opening the information File" , dataIngestInfo) ;
    LOG(INFO)<<"Successfully opened"<<endl;
    fits_movabs_hdu (fdi , 2 , NULL , &status) ;
    printError (status , "Error Moving to 2nd HDU in  information File" , dataIngestInfo) ;
    fits_read_key (fdi , TSTRING , "DARKDIR" , darkdir , NULL , &status) ;
    printError (status , "Error Reading the Key value for DARKDIR" , dataIngestInfo) ;
    fits_close_file(fdi,&status);
    printError (status , "Error in closing the file" , dataIngestInfo) ;
    LOG(INFO)<<"Successfully closed"<<endl;
    
    //opening info file in input directory to get data information
    tempfilepath = searchFile (inputdatadir , ".info") ;
    if (tempfilepath == "")
    {
        LOG(ERROR) << endl << "\033[1;31m Error in finding info file" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    LOG(INFO)<<"\nInput information file :"<<infofile_in<<endl;
  
    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error Opening the information File" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error Moving to 2nd HDU in  information File" , infofile_in) ;
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error Reading the Key value for NAMEPRFX" , infofile_in) ;
    fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
    printError (status , "Error Reading the Key value for EVTFILE" , infofile_in) ;

    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file

    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;

    if (xsize <= 0 || ysize <= 0)
    {
        LOG(ERROR) << endl << "***Invalid xsize/ysize***\n" ;
        return (EXIT_FAILURE) ;
    }

    //creating output information file
    sprintf (infofile_out , "%s/%s_cc.info" , moduleoutdir , nameprefix) ;
    LOG(INFO)<<"\nOutput information file:"<<infofile_out<<endl;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error Creating the output  information File" , infofile_out) ;
    char *ttype[] = {"LIST"} ;
    char *tform[] = {"A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error creating the table in output information file" , infofile_out) ;
    datainfo.write (finfo_out) ; //writing basic data information
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the NAMEPRFX key " , infofile_out) ; //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    fits_close_file (finfo_in , &status) ;
    printError (status , "Error in closing the  input  information File" , infofile_in) ;
   
    string tempname = caldb_handler.getCentroidEAFile (datainfo.getDetector () , caldbDir) ;
    if (tempname == "")
    {
        LOG(ERROR) << endl << "Couldn't find Centroid correction  file from caldb" << endl ;
        return (EXIT_FAILURE) ;
    }
    joinStrings (centroidEAfile , 2 , caldbDir , tempname.c_str()) ;
    LOG(INFO) << endl << "centroid EA  file is " << centroidEAfile ;

    /***Reading EA(Average Energy )from caldb file***/
    status = readcentroidEAFile () ;
    if (status)
    {
        LOG(ERROR) << "***Error Reading a EA value From the Caldb***" << endl ;
        return (EXIT_FAILURE) ;
    }
    if (EA == 0)
        {
            LOG(ERROR) << "*** Divide by zero***" << endl ;
            return (EXIT_FAILURE) ;
        }
   
    /***Method  For reading the Dark Frame generated from the pre-Processing Block
     for further processing ***/
    darkFramestart_data = new float[DARKFRAME_SIZE * DARKFRAME_SIZE] ;
    darkFrameend_data = new float[DARKFRAME_SIZE * DARKFRAME_SIZE] ;
    status=takeDarkinfo ();
     if (status)
    {
        LOG(ERROR) << "***Error getting the information of Darks***" << endl ;
        return (EXIT_FAILURE) ;
    }
    status = readDarkFrame (dstartpath , darkFramestart_data) ;
    if (status)
    {
        LOG(ERROR) << "***Error Reading a start DarkFrame***" << endl ;
        return (EXIT_FAILURE) ;
    }
    status = readDarkFrame (dendpath , darkFrameend_data) ;
    if (status)
    {
        LOG(ERROR) << "***Error Reading a  End DarkFrame***" << endl ;
        return (EXIT_FAILURE) ;
    }
   
    char file_in[NAMESIZE] , file_out[NAMESIZE] ; //for input and output event file
    sprintf (file_in , "%s/%s" , inputdatadir , eventfile) ;
    sprintf (file_out , "%s/%s_cc.events" , moduleoutdir , nameprefix) ;
    LOG(INFO)<<"\nInput Event File :"<<file_in<<endl;
    LOG(INFO)<<"\nOutput Event File :"<<file_out<<endl;
    fits_update_key (finfo_out , TSTRING , "EVTFILE" , basename (file_out) , "File name prefix" , &status) ;
    printError (status , "Error  in updating the Key value for  EVTFILE " , infofile_out) ;
    fitsfile *fevt_in , *fevt_out;
    fits_open_file (&fevt_in , file_in , READWRITE , &status) ;
    printError (status , "Error in opening the  input event file" , file_in) ;
    copyUsrkeywrdsTovect (fevt_in,key_record);
    fits_read_key (fevt_in , TSTRING , "CENTROID" , centroid_algo , NULL , &status) ;
    printError (status , "Error Reading the Key value for CENTROID" , file_in) ;
    LOG(INFO)<<"CentAlgo value is "<<centroid_algo<<endl;
    if ((strcmp (centroid_algo , "3C") != 0) && (strcmp (centroid_algo , "3S") != 0) && (strcmp(centroid_algo,"5S")!=0)){
           LOG(ERROR)<<"***Invalid value of Centroid_algo ***"<<endl;
           return(EXIT_FAILURE);           
    }
    if ((strcmp (centroid_algo , "3C") == 0))         windowsize = 5 ;          //3x3 cross  (5 pixels - centre, top, bottom, left & right used)
    else if (strcmp (centroid_algo , "5S") == 0)      windowsize = 25 ;         //5x5 square
    else if( (strcmp (centroid_algo , "3S") == 0))    windowsize=9;             //3x3 square
    else
    {
        LOG(ERROR) << "***INVALID value for Centroid window.allowed values are (3/5) *** " << endl ;
        return (EXIT_FAILURE) ;
    }
    fits_create_file (&fevt_out , file_out , &status) ;
    printError (status , "Error in creating the output information File" , file_out) ;
    fits_copy_file(fevt_in ,fevt_out , 1,1,1,&status) ;
    printError (status , "Error in coping the header" , file_out) ;
    /*declaration of variables used for getting column value for  input event file*/
    float *xFrac , *yFrac , *new_xFrac , *new_yFrac ;
    double *t ;
    unsigned short *xi , *yi,*mc;
    long nrows ;
    double time_Interval ;
    fits_movabs_hdu (fevt_out , 2 , NULL , &status) ;
        printError (status , "Error in Moving to 2nd HDU",file_out) ;
        fits_get_num_rows (fevt_out , &nrows , &status) ;
        printError (status , "Error in getting the number of rows from input event File",file_out) ;
        LOG(INFO)<<"\nTotal number of Events :"<<nrows<<endl;
        /***Assigning memory to Arrays***/
        xi = new unsigned short[nrows] ;
        yi = new unsigned short[nrows] ;
        xFrac = new float[nrows] ;
        yFrac = new float[nrows] ;

        new_xFrac = new float[nrows] ;
        new_yFrac = new float[nrows] ;
        t = new double[nrows] ;
        mc = new unsigned short[nrows] ;
  
        for (int i = 0 ; i < nrows ; i++)
        {
            xi[i] = (unsigned short) 0 ;
            yi[i] = (unsigned short) 0 ;
            new_xFrac[i] = 0.0 ;
            new_yFrac[i] = 0.0 ;
        }
        //added..
        fits_read_col (fevt_out , TDOUBLE , 3 , 1 , 1 , nrows , NULL , t , NULL , &status) ;
        printError (status , "Error reading xinteger" , file_out) ;
        //
        fits_read_col (fevt_out , TUSHORT , 4 , 1 , 1 , nrows , NULL , xi , NULL , &status) ;
        printError (status , "Error reading xinteger" , file_out) ;
        fits_read_col (fevt_out , TUSHORT , 6 , 1 , 1 , nrows , NULL , yi , NULL , &status) ;
        printError (status , "Error reading yinteger" , file_out) ;
        fits_read_col (fevt_out , TFLOAT , 5 , 1 , 1 , nrows , NULL , xFrac , NULL , &status) ;
        printError (status , "Error reading the xfractional" , file_out) ;
        fits_read_col (fevt_out , TFLOAT , 7 , 1 , 1 , nrows , NULL , yFrac , NULL , &status) ;
        printError (status , "Error reading the yfractional" , file_out) ;
        fits_read_col (fevt_out , TUSHORT , 9 , 1 , 1 , nrows , NULL , mc , NULL , &status) ;
        printError (status , "Error reading the minimum corner value" , file_out) ;
        time_Interval = datainfo.getIntegrationTime () ;
        LOG(INFO) << "Time interval(Integration time)" << time_Interval << endl ;
        double t1 , t2 ;
        float dark_value = 0.0 ;
        float sd = 0.0 , sc = 0 ;
         t2 = t[nrows - 1] + time_Interval ;
          int index_x;
          int index_y;
          LOG(INFO)<<"\nPerforming Centroid Correction......"<<endl;
        //  long sqr_winsize =windowsize*windowsize;
        int diff_add=xsize-DEFAULT_FRAMESIZE;        
          for (int i = 0 ; i < nrows ; i++)
        {
           dark_value = 0.0f ;
            t1 = t[i] - time_Interval ; 
            if ((t2 - t1) == 0)
            {
                LOG(ERROR) << "***Divide by Zero***" << endl ;
                return (EXIT_FAILURE) ;
            }
            index_x = (int) xi[i] ;
            index_y = (int) yi[i] ;
            new_xFrac[i] = 0.0f ;
            new_yFrac[i] = 0.0f ;
            if ((index_x - diff_add) >= 0 && (index_y - diff_add) >= 0)
            {
               dark_value = darkFramestart_data[(index_y - diff_add)*DEFAULT_FRAMESIZE + (index_x - diff_add)] +
                                      ((t[nrows - 1] - t1) * ((darkFrameend_data[(index_y - diff_add)*DEFAULT_FRAMESIZE + (index_x - diff_add)] 
                                       - darkFramestart_data[(index_y - diff_add)*DEFAULT_FRAMESIZE + (index_x - diff_add)]) / (t2 - t1))) ;
            }
            sc = mc[i] ;
            sd = dark_value ;
            new_xFrac[i] = xFrac[i]*(1 + ((sd - sc)*windowsize) / EA) + index_x ;
            new_yFrac[i] = yFrac[i]*(1 + ((sd - sc)*windowsize) / EA) + index_y ;
         }
          
        char *ttype1="X" ,*ttype2="Y";
        char *tform1=" 1D",*tform2="1D" ;
        /**deleting the cloumns of xi, yi, xfract, and yfrac**/
       for(int drown=1;drown<5;drown++)
        fits_delete_col(fevt_out,4,&status);
        
        fits_insert_col(fevt_out,4,(char*)ttype1,(char*)tform1,&status);
        fits_insert_col(fevt_out,5,(char*)ttype2,(char*)tform2,&status);
        fits_write_col (fevt_out , TFLOAT , 4 , 1 , 1 , nrows , new_xFrac , &status) ;
        printError (status , "Error writing the Xinteger" , file_out) ;
        fits_write_col (fevt_out , TFLOAT , 5 , 1 , 1 , nrows , new_yFrac , &status) ;
        printError (status , "Error writing the Yinteger" , file_out) ;
        delete[] new_xFrac ;
        delete[] new_yFrac ;
        delete[] xFrac , yFrac , xi , yi , mc ;
    fits_movabs_hdu (fevt_out , 1 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU" , file_out) ;
    /***Writes the  common information to  output information File***/
    
    vector<string> vhistorystr ;
    getHistory (vhistorystr) ;
    if (history == YES) writeHistory (file_out , vhistorystr) ;
      fits_movabs_hdu (fevt_out , 1 , NULL , &status) ;
    printError (status , "Error in  moving to perticuler HDU" , file_in) ;
    writeCommonKeywords (fevt_out , modulename) ;
    fits_close_file (fevt_out , &status) ;
    printError (status , "***Error in Closing the output Event File***" , file_in) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "***Error in closing the info out file***" , file_out) ;
    writeUsrkeywordsFrmvect(infofile_out,key_record);
    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    fits_close_file (fevt_in , &status) ;
    printError (status , "***Error in closing the input Event File***" , file_out) ;
    return (EXIT_SUCCESS) ;
}
int uvtCentroidCorr::readcentroidEAFile ()
{
    LOG(INFO) << endl << "\nReading Centroid EA file from calDB." ;
    fitsfile *f_ea ;
    int status = 0 ;
    long n_ele ;
    fits_open_file (&f_ea , centroidEAfile , READONLY , &status) ;
    printError (status , "Error in opening the centroid caldb file" , centroidEAfile) ;
    fits_movabs_hdu (f_ea , 2 , NULL , &status) ;
    printError (status , "Error in moving to the HDU" , centroidEAfile) ;
    fits_get_num_rows (f_ea , &n_ele , &status) ;
    printError (status , "Error in getting the number of rows" , centroidEAfile) ;
    fits_read_col (f_ea , TINT , 1 , 1 , 1 , n_ele , NULL , (void *) &EA , NULL , &status) ;
    printError (status , "Error in closing the centroid corr caldb  file" , centroidEAfile) ;
    fits_close_file (f_ea , &status) ;
    return (EXIT_SUCCESS) ;
}
int uvtCentroidCorr::readDarkFrame (char * path , float *Array)
{
    LOG(INFO)<<"\nReading Dark frame : "<<path<<endl;
    fitsfile *fptr ;
    int status = 0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    for (int i = 0 ; i < DARKFRAME_SIZE * DARKFRAME_SIZE ; i++)
    {
        Array[i] = (float) 0.0 ;
    }
    fits_open_file (&fptr , path , READONLY , &status) ;
    printError (status , "Error in opening the Dark frame" , path) ;
    fits_read_pix (fptr , TFLOAT , fpixel , DARKFRAME_SIZE*DARKFRAME_SIZE , NULL , Array , NULL , &status) ;
    printError (status , "Error in reading the pixels from the Dark Frames" , path) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in Closing the  Dark Frame" , path) ;
    return(EXIT_SUCCESS);
}
int uvtCentroidCorr::takeDarkinfo ()
{
    
    double dark_Firstframetime , dark_Secframetime , start_i , start_f , stop_i , stop_f ;
    char darkpath_First[FLEN_FILENAME] ;
    char darkend_Sec[FLEN_FILENAME] ;
    char dark_temp[FLEN_FILENAME] ;

    //vector for storing the dark frame names
    vector<string> dark_framenames ;

    //setting path for dark frame directory  
    sprintf (dark_temp , "%s/%s" , dataIngestDir , darkdir) ;

    getFileNamesfrmDir (dark_temp , ".fits" , dark_framenames) ; //get dark file names from dark directory

    if (dark_framenames.size () <2)
    {
        LOG (INFO) << "No enough Dark frames found at input Directory" << endl ;
        return (EXIT_SUCCESS) ;
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
//    
//    double dark_framestarttime,dark_frameendtime,start_i,start_f,stop_i,stop_f;
//    double timeDark;
//    vector<string> dark_framenames;
//    char dark_temp[FLEN_FILENAME];
//    sprintf(dark_temp,"%s/%s",dataIngestDir,darkdir);
//    getFileNamesfrmDir (dark_temp,".fits",dark_framenames);
//    if (dark_framenames.size ()<2){
//        LOG(INFO)<<"No enough Dark frames found at input Directory"<<endl;
//        return(EXIT_SUCCESS);
//    }
//    char darkstart_path[FLEN_FILENAME];
//    char darkend_path[FLEN_FILENAME];
//    int status=0;
//    fitsfile *darkfile; 
//    sprintf(darkstart_path,"%s/%s",dark_temp, dark_framenames[0].c_str ());
//    cout<<darkstart_path<<endl;
//    
//    fits_open_file (&darkfile , darkstart_path , READONLY , &status) ;
//    printError (status , "Error in opening the output File For IM mode" , (char*)dark_framenames[0].c_str ()) ;
//    fits_read_key (darkfile , TDOUBLE , "FRMTIME" , &dark_framestarttime , NULL , &status) ;
//    printError (status , "Error reading the key value of the FRMTIME" , (char*)dark_framenames[0].c_str ()) ;
//    fits_close_file(darkfile,&status);
//    printError (status , "Error in closing the dark file" , (char*)dark_framenames[0].c_str ()) ;
//    
//    sprintf(darkend_path,"%s/%s",dark_temp, dark_framenames[1].c_str ());
//    fits_open_file (&darkfile , darkend_path , READONLY , &status) ;
//    printError (status , "Error in opening the output File For IM mode" , (char*)dark_framenames[1].c_str ()) ;
//    fits_close_file(darkfile,&status);
//    printError (status , "Error in closing the dark file" , (char*)dark_framenames[0].c_str ()) ;
//    
//
//    
//    if(dark_framestarttime<dark_frameendtime){
//        strcpy(dstartpath,darkstart_path);
//        strcpy(dendpath,darkend_path);
//        
//    }
//    else{
//        
//        strcpy(dstartpath,darkend_path);
//        strcpy(dendpath,darkstart_path);
//    }
//    
    return(EXIT_SUCCESS);
}
