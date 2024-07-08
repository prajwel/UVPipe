/* 
 * File:   uvtOpticAssDistCorr.cpp
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
#include <uvtOpticAssDistCorr.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<vector>
#include<algorithm>
#include<iterator>
#include<glog/logging.h>


//Constructor -called when object is created

uvtOpticAssDistCorr::uvtOpticAssDistCorr ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    //strcpy(centroidDir, "Centroid");
}

//Destructor

uvtOpticAssDistCorr::~uvtOpticAssDistCorr () {
 }

//parameter File reading 

int uvtOpticAssDistCorr::read (int argc , char** argv)
{
    int status = 0 ;
    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG(INFO) << "***Error Initializing PIL***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("inputdatadir" , inputdatadir)))
    {
        LOG(INFO) << endl << "***Error reading input Directory***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("caldbDir" , caldbDir)))
    {
        LOG(INFO) << endl << "***Error reading Caldb  Directory Path ***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("outdir" , outdir)))
    {
        LOG(INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("clobber" , &clobber)))
    {
        LOG(INFO) << "***Error Reading clobber:" << clobber << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("history" , &history)))
    {
        LOG(INFO) << "***Error Reading history parameter:" << history << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("mode" , mode)))
    {
        LOG(INFO) << "***Error Reading mode parameter:" << history << "***" ;
        return status ;
    }
    PILClose (status) ;
    return (EXIT_SUCCESS) ;
}

int uvtOpticAssDistCorr::read (char *inputdatadir , char * caldbDir , char *outdir , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;
    strcpy (this->caldbDir , caldbDir) ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}
//Parameter file content Display

void uvtOpticAssDistCorr::display ()
{

   LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT OPTICAL DISTORTION PARAMETERS      " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input Frame List Directory                        : " << inputdatadir ;
    LOG(INFO) << endl << "Caldb Directory                                           : " << caldbDir ;
    LOG(INFO) << endl << "Output Directory                               : " << outdir ;
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

//Correction for the  Cosmic Ray process

int uvtOpticAssDistCorr::uvtOpticalDistcorrProcess ()
{
    /*if the Input Directory of the Filelist is not available then exit from module
     else read Filelist 
     */
    LOG(INFO) << endl << "Optical Distrotion process started" << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG(INFO)<<endl<<"Module Output Directory : "<<moduleoutdir<<endl;
    string cmd ;
    //check For The Directory Existence
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
    /**Shell command for creating output directory**/
    cmd = "mkdir -p " + (string) moduleoutdir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ; // creating output directory to keep output from unitConversion
    LOG(INFO) << endl << moduleoutdir << "  directory created" ;
    /*For Searching The Information File  (.info) From the Directory ,It Returns the filename of the .info file which
     *  contain framelist to be processed */
    string  tempfilepath = searchFile (inputdatadir , ".info") ;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (tempfilepath == " ")
    {
        LOG(INFO) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY" << endl ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    LOG(INFO) << endl << "\nInformation File :" << infofile_in ;
    /*Information File Processing Started*/
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    /**opening  the input information file**/
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Information File not Opened" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Moving to 2nd Hdu Fails For Information File Fails" , infofile_in) ;
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file

    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    if (xsize <= 0 || ysize <= 0)
    {
        LOG(INFO) << endl << "***Invalid xsize/ysize***\n" ;
        return (EXIT_FAILURE) ;
    }
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "NAMEPRFX keyword not Found" , infofile_in) ; //for creating name for output information file
    fits_read_key (finfo_in , TSTRING , "FILTER" , filter , NULL , &status) ;
    printError (status , "Filter  keyword not Found" , infofile_in) ;
    sprintf (infofile_out , "%s/%s_od.info" , moduleoutdir , nameprefix) ;
    /*creating the output information file*/
    LOG(INFO) << endl << "\n Output Information File :" << infofile_out ;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    char *ttype[] = {"CentroidFrames"} ;
    char *tform[] = {"A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in Creating the table") ;
    datainfo.setXsize (xsize) ;
    datainfo.setYsize (ysize) ;
    datainfo.write (finfo_out) ; //writing basic data information
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the NAMEPRFX keyword" , infofile_out) ; //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    /*For  Further Proccessing Take input From Caldb Directory (here take Optical Detector File From Caldb) correction on the Frames Will be applied based
     on this File(caldb) */
    string tempname = caldb_handler.getOpticalDetectorFile (datainfo.getDetector () , caldbDir , filter) ;
    if (tempname == "")
    {
        LOG(INFO) << endl << "Couldn't find Correction file from caldb" << endl ;
        return (EXIT_FAILURE) ;
    }
    /**Setting the path of file for distortion correction from calDB directory **/
    joinStrings (distortionCorrfile , 2 , caldbDir , tempname.c_str()) ;
    LOG(INFO) << endl << "\nDistortion correction file is " << distortionCorrfile<<endl ;
    status = caldb_handler.readCaldbOpticDistFile (x_Distortion , y_Distortion , distortionCorrfile) ; //Distortion File Reading From caldb.
    if (status)
    {
        LOG(ERROR) << "Error in reading the Distortion file" << endl ;
        return (EXIT_FAILURE) ;
    }
    if (datainfo.getModeFlag () == IM)
    {//incase of IM mode
        fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "Error in reading the  key value of NFILES " , infofile_in) ;
        fits_read_key (finfo_in , TSTRING , "CENTDIR" , centroidDir , NULL , &status) ;
        printError (status , "Error  reading the key value of the CENTDIR " , infofile_in) ;
        centroidframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) centroidframelist , NULL , &status) ;
        printError (status , "Error in reading the  column of centroidframe list in input information file" , infofile_in) ;
        /*For IM mode */
        if (opticalDistortionIM ()) return (EXIT_FAILURE) ;
        /* updating of some Information open information file for this module*/
        fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
        printError (status , "Error in opening the out  information file" , infofile_in) ;
        fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
        printError (status , "Error in moving the 2nd HDU of  information file" , infofile_in) ;
        fits_update_key (finfo_out , TSTRING , "CENTDIR" , centroidDir , NULL , &status) ;
        printError (status , "Error in updating the key value of the CENTDIR" , infofile_in) ;
        fits_update_key (finfo_out , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "Error in updating the key valu of the NFILES" , infofile_in) ;
        fits_close_file (finfo_out , &status) ;
        printError (status , "Error in closing the file" , infofile_in) ;
        freeMemory (centroidframelist , nframes , NAMESIZE) ; //for releasing the memory
    } //in case of the PC 
    else if (datainfo.getModeFlag () == PC)
    { //incase of PC mode
        fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
        printError (status , "Error in reading the key value of the EVTFILE" , infofile_in) ;
        if (opticalDistortionPC ()) return (EXIT_FAILURE) ;
    }
    else//else neither PC or IM(i.e invalid mode )
    {
        LOG(ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG(ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
    fits_close_file(finfo_in,&status);
    printError (status , "Error in closing  file", infofile_in) ;
    
    return (EXIT_SUCCESS) ;
}

int uvtOpticAssDistCorr::opticalDistortionPC ()
{
    LOG(INFO) << "\nStarted Correction for Distortion in optical assembly for PC mode " << endl ;
    char file_in[NAMESIZE] , file_out[NAMESIZE] ; //for input and output event file
    //opening eventfile
    sprintf (file_in , "%s/%s" , inputdatadir , eventfile) ; //taking event file full path
    sprintf (file_out , "%s/%s_od.events" , moduleoutdir , nameprefix) ;
    fitsfile *fevt_in , *fevt_out , *finfo_out ;
    int status = 0 ;
    /**opening the input EVENT file**/
    fits_open_file (&fevt_in , file_in , READWRITE , &status) ;
    printError (status , "Error in opening the input event file" , file_in) ;
    copyUsrkeywrdsTovect (fevt_in,key_records);    
    fits_create_file (&fevt_out , file_out , &status) ;
    printError (status , "Error in creating the out event file" , file_in) ;
    /**Copying the input event file to  output event file**/
    fits_copy_file (fevt_in , fevt_out , 1 , 1 , 1 , &status) ;
    printError (status , "Error in coping the file" , file_out) ;

    /***Difining vars for reading a event file columns***/
    float *xFrac , *yFrac ;
    long nrows ;
    fits_movabs_hdu (fevt_out , 2 , NULL , &status) ;
    printError (status , "Error in  moving to perticuler HDU" , file_out) ;
    fits_get_num_rows (fevt_out , &nrows , &status) ;
    printError (status , "Error in  getting the number of rows" , file_out) ;
    /***Assigning a memory to Arrays***/
    xFrac = new float[nrows] ;
    yFrac = new float[nrows] ;

    fits_read_col (fevt_out , TFLOAT , 4 , 1 , 1 , nrows , NULL , xFrac , NULL , &status) ;
    printError (status , "Error reading xfractional" , file_out) ;
    fits_read_col (fevt_out , TFLOAT , 5 , 1 , 1 , nrows , NULL , yFrac , NULL , &status) ;
    printError (status , "Error reading yfractional" , file_out) ;
    float locate , tempx , tempy , multi_factor = 0.0 ;
    multi_factor = xsize / 600 ;
    /**Appling correction on XFract and yFrac**/
    LOG(INFO)<<"Applying Correction..."<<endl;
    for (int k = 0 ; k < nrows ; k++)
    {
        tempx = xFrac[k] / multi_factor ;
        tempy = yFrac[k] / multi_factor ;

        locate = ((int) round (tempy) - 44) * 512 + ((int) round (tempx) - 44) ;
        round (locate) ;

        tempx = tempx + x_Distortion[(int) locate] ;
        tempy = tempy + y_Distortion[(int) locate] ;
  if(tempx*multi_factor <xsize && tempy*multi_factor<ysize) { 
        xFrac[k] = tempx * multi_factor ;
        yFrac[k] = tempy * multi_factor ;
  }

    }
    /***Writing to output event file***/
    fits_write_col (fevt_out , TFLOAT , 4 , 1 , 1 , nrows , xFrac , &status) ;
    printError (status , "Error writing Xfractional " , file_out) ;
    fits_write_col (fevt_out , TFLOAT , 5 , 1 , 1 , nrows , yFrac , &status) ;
    printError (status , "Error writing Yfractional " , file_out) ;
    /*releasing the memory*/
    delete[] xFrac ;
    delete[] yFrac ;
    fits_movabs_hdu (fevt_out , 1 , NULL , &status) ;
    printError (status , "Moving to 1st HDU fails for output event file" , file_out) ;
    /***Writes the  common information to  output information File***/
   
    vector<string> vhistorystr ;
    if (history == YES)
    {
        getHistory (vhistorystr) ;
        writeHistory (file_out , vhistorystr) ;
    }
    writeCommonKeywords (fevt_out , modulename) ;
    fits_close_file (fevt_out , &status) ;
    printError (status , "Error in Closing the output Event File" , file_out) ;
    fits_close_file (fevt_in , &status) ;
    printError (status , "Error in closing the input Event File" , file_out) ;
    LOG(INFO)<<"\nUpdating  output information file.."<<endl;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in openig the output information file" , file_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU of output information file" , file_out) ;
    fits_update_key (finfo_out , TSTRING , "EVTFILE" , basename (file_out) , NULL , &status) ;
    printError (status , "Error in updating the key value of the EVTFILE" , file_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the file" , file_out) ;
    writeUsrkeywordsFrmvect (infofile_out,key_records);
    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    return (EXIT_SUCCESS) ;
}

int uvtOpticAssDistCorr::opticalDistortionIM ()
{
   LOG(INFO) << "\nStarted Correction for Distortion in optical assembly for IM mode " << endl ;
    char **outcentroidframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    int status = 0 ;
    fitsfile *finfo_out ;
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s%s" , moduleoutdir , centroidDir) ;
    /**Shell command for creating the directory**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" ;
    long fpixel[2] ;
    fpixel[0] = 1 , fpixel[1] = 1 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    vector<string> vhistorystr ;
    if (history == YES)
    {
        getHistory (vhistorystr) ;
    }
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;
   
    float *Xloc , *Yloc ;
     LOG(INFO) << "\nApplying Correction..." << endl ;
    for (int i = 0 ; i < nframes ; i++)
    {
        fitsfile *fptr , *fout ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , centroidDir , centroidframelist[i]) ;
        /**opening the input Centroid file**/
        
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Reading input file failed" , infile) ;
        long numrows = 0 ;
       

        fits_movabs_hdu (fptr , 2 , NULL , &status) ;
        printError (status , "Moving To 2nd HDU for the Frame Fails" , infile) ;
        fits_get_num_rows (fptr , &numrows , &status) ;
        printError (status , "Error in getting the number of rows" , infile) ;
        
        fits_read_key (fptr , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ;
        printError (status , "FRAMENO not Found" , infile) ;
        fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
        printError (status , "FRMTIME not found" , infile) ;
        /**Creating the output file **/
        sprintf (outfile , "%s/%s/%s_t%.4f_f%d_centroid_od.fits" , moduleoutdir , centroidDir , nameprefix , frametime , frameno) ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output Event File ",outfile) ;
        fits_movabs_hdu (fptr , 1 , NULL , &status) ;
        printError (status , "Moving to  2nd HDU Fails",outfile) ;
        fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
        printError (status , "Error in coping the file" , outfile) ;
        Xloc = new float[numrows] ;
        Yloc = new float[numrows] ;
        fits_read_col (fout , TFLOAT , 1 , 1 , 1 , numrows , NULL , (void *) Xloc , NULL , &status) ;
        printError (status , "Error in reading the column of Xloc" , outfile) ;
        fits_read_col (fout , TFLOAT , 2 , 1 , 1 , numrows , NULL , (void *) Yloc , NULL , &status) ;
        printError (status , "Error in reading the column of  Yloc" , outfile) ;
        float tempx , tempy ;
        double locate ;
        float multi_factor = xsize / 600 ;
        if (multi_factor == 0)
        {
            LOG(INFO) << "Divide by zero" << endl ;
            return (EXIT_FAILURE) ;
        }
        for (int i = 0 ; i < numrows ; i++)
        {
            if(Xloc[i] !=INVALID_PIX_VALUE && Yloc[i] !=INVALID_PIX_VALUE){
            tempx = Xloc[i] / multi_factor ;
            tempy = Yloc[i] / multi_factor ;

            locate = ((int) round (tempy) - 44) * 512 + ((int) round (tempx) - 44) ;
            round (locate) ;

            tempx = tempx + x_Distortion[(int) locate] ;
            tempy = tempy + y_Distortion[(int) locate] ;
            
         if(tempx*multi_factor <xsize && tempy*multi_factor<ysize) { 
            Xloc[i] = tempx * multi_factor ;
            Yloc[i] = tempy * multi_factor ;
         }
            }
        }
        fits_write_col (fout , TFLOAT , 1 , 1 , 1 , numrows , (void *) Xloc , &status) ;
        printError (status , "Error in writing the column of Xloc " , outfile) ;
        fits_write_col (fout , TFLOAT , 2 , 1 , 1 , numrows , (void *) Yloc , &status) ;
        printError (status , "Error in writing the column of Yloc" , outfile) ;
       // copyUserKeywords (fptr , fout) ;
       // writeCommonKeywords (fout , modulename) ;
        // writeUsrkeywordsFrmvect (outfile,key_records);
        if (history == YES) writeHistory(outfile , vhistorystr) ;
         writeCommonKeywords (fout,modulename);
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output file" , outfile) ;
       
         
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input file" , outfile) ;
        strcpy (outcentroidframelist[i] , basename (outfile)) ;
         delete[] Xloc , Yloc ;
         cout<< "Total Files written=" << i << "     Remaining files= " << nframes - i <<" \r" ;
    }
     LOG(INFO)<<"\nwriting list of output frame names to the output information file"<<endl;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening  the output information path" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU of output information file" , infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , nframes , (void *) outcentroidframelist , &status) ;
    printError (status , "Writing To column of centroid listing Fails" , infofile_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file" , infofile_out) ;
    writeUsrkeywordsFrmvect (infofile_out,key_records);
    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    return (EXIT_SUCCESS) ;
}
int uvtOpticAssDistCorr::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Centroid Directory=" + (string) centroidDir) ;
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
