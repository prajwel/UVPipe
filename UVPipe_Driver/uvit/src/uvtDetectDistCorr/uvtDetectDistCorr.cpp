/* 
 * File:   uvtDetectDistCorr.cpp
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
#include <uvtDetectDistCorr.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<caldb_Handler.h>
#include<glog/logging.h>


//Constructor -called when object is created

uvtDetectDistCorr::uvtDetectDistCorr ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    strcpy (centroidDir , "Centroid") ;
}

//Destructor

uvtDetectDistCorr::~uvtDetectDistCorr () {
    
 }
//parameter File reading
int uvtDetectDistCorr::read (int argc , char** argv)
{
    int status = 0 ;
    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG(ERROR) << "***Error Initializing PIL***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("inputdatadir" , inputdatadir)))
    {
        LOG(ERROR) << endl << "***Error reading input directory***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("caldbDir" , caldbDir)))
    {
        LOG(ERROR) << endl << "***Error reading CalDB  Directory Path ***" ;
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
        LOG(ERROR) << "***Error Reading history parameter:" << history << "***" ;
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

int uvtDetectDistCorr::read (char *inputdatadir , char * caldbDir , char *outdir , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;
    strcpy (this->caldbDir , caldbDir) ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}
//Parameter file content Display

void uvtDetectDistCorr::display ()
{ LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT DETECTOR DISTORTION PARAMETERS      " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input Frame List Directory                        : " << inputdatadir ;
    LOG(INFO) << endl << "CaldDB Directory                                           : " << caldbDir ;
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

int uvtDetectDistCorr::uvtDetectDistCorrProcess ()
{
    /*if the Input Directory of the Filelist is not available then exit from module
     else read Filelist
     */
    LOG(INFO) << endl << "Detect Distortion process started" << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
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
    /**Shell Command For creating the Directory**/
    cmd = "mkdir -p " + (string) moduleoutdir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ; // creating output directory to keep output from unitConversion
    LOG(INFO) << endl << moduleoutdir << "  directory created" ;
    string  tempfilepath = searchFile (inputdatadir , ".info") ;
     if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }
    /**Setting the path for the input information File**/
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
   
    if (!(FileExists (infofile_in)))
    {
        LOG(ERROR) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY" ;
        return (EXIT_FAILURE) ;
    }
  
    LOG(INFO) << endl << "\nInformation File :" << infofile_in ;
    
    int status = 0 ;
    /**open a input Information File**/
    fitsfile *finfo_in , *finfo_out ;
   
    /**Reading the information File from information file**/
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the information file",infofile_in) ;
      
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in Moving the 2nd HDU",infofile_in) ;
    /**Reading the keyword information from the input information File**/
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
   
   
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    if (xsize <= 0 || ysize <= 0)
    {
        LOG(ERROR) << endl << "***Invalid xsize/ysize***\n" ;
        return (EXIT_FAILURE) ;
    }
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "NAMEPRFX keyword not Found",infofile_in) ; //for creating name for output information file
    /**Setting a centroid Directory**/
    strcpy (centroidDir , "Centroid") ;
    sprintf (infofile_out , "%s/%s_dd.info" , moduleoutdir , nameprefix) ;
    LOG(INFO)<<"Output Information File :"<<infofile_out<<endl;
    /**Creating the output information File**/
    
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information file",infofile_out) ;
    char *ttype[] = {"CentroidFrames"} ;
    char *tform[] = {"A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table",infofile_out) ;
    /**Write keyword information to the output information File**/
    datainfo.write (finfo_out) ; //writing basic data information
    
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the key value of NAMEPRFX",infofile_out) ; //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    fits_close_file (finfo_out , &status) ;
    
    
    string tempname = caldb_handler.getDetectorFile (datainfo.getDetector () , caldbDir) ;
    if (tempname == " ")
    {
        LOG(ERROR) << endl << "Couldn't find Detector Distortion  file from calDB"<<endl ;
        return (EXIT_FAILURE) ;
    }

    joinStrings (distortionCorrfile , 2 , caldbDir , tempname.c_str()) ;
    LOG(INFO) << endl << "\nDistortion correction file :" << distortionCorrfile ;
 
    status = caldb_handler.readCaldbDistFile (x_Distortion , y_Distortion , distortionCorrfile) ; //Distortion File Reading From caldb.
    if (status)
    {
        LOG(ERROR) << "***Error in reading the Distortion file***" << endl ;
        return (EXIT_FAILURE) ;
    }
    if (datainfo.getModeFlag () == IM) //For IM mode
    {
        fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "NFILES keyword not found ",infofile_in) ;
        centroidframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        fits_read_col (finfo_in , TSTRING , 2 , 1 , 1 , nframes , NULL , (void *) centroidframelist , NULL , &status) ;
       //method for the Detector Distortion  For IM
        if (detectDistortionIM ()) return (EXIT_FAILURE) ;
      
        fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
        printError (status , "Error in opening the out  information file",infofile_out) ;
        fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
        printError (status , "Error in moving the 2nd HDU of  information file",infofile_out) ;
        fits_update_key (finfo_out , TSTRING , "CENTDIR" , &centroidDir , NULL , &status) ;
        printError (status , "Error in updating the key value of the CENTDIR",infofile_out) ;
        fits_update_key (finfo_out , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "NFILES keyword  not updated/not Found",infofile_out) ;
        fits_close_file (finfo_out , &status) ;
        printError (status , "Error in closing the file",infofile_out) ;
        freeMemory (centroidframelist , nframes , NAMESIZE) ; //for releasing the memory
    }//in case of the PC
    else if (datainfo.getModeFlag () == PC)
    {
        fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
        printError (status , " EVTFILE keyword not found ",infofile_in) ;
        if (detectDistortionPC ()) return (EXIT_FAILURE) ;
    }
    else //else neither PC or IM(i.e invalid mode )
    {
        LOG(ERROR) << endl << "Invalid operating mode found" ;
        LOG(ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
    fits_close_file(finfo_in,&status);
    printError (status , "Error in closing the file ",infofile_in) ;
    
    return (EXIT_SUCCESS) ;
}

int uvtDetectDistCorr::detectDistortionIM ()
{
     LOG(INFO)<<"Started  correction for  detector Distortion process for IM mode "<<endl;
    char **outcentroidframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    int status = 0 ;
    fitsfile *finfo_out ;
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s%s" , moduleoutdir , centroidDir) ;
    string cmd = "mkdir -p " + (string) dir ;
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" ;
    long fpixel[2] ;
    fpixel[0] = 1 , fpixel[1] = 1 ;
     long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    int numHDUS=0 ;
    vector<string> vhistorystr ;
    if (history==YES)
    {
            getHistory (vhistorystr) ;
    }
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;
    LOG(INFO)<<"Applying  correction..."<<endl;
     //nframes To be Processed
    for (int i = 0 ; i < nframes ; i++)
    {
         fitsfile *fptr , *fout ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , centroidDir , centroidframelist[i]) ;
        /**opening the input Centroid frame **/
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input Centroid file" , infile) ;
        fits_get_num_hdus (fptr , &numHDUS , &status) ;
        printError (status , "Error in getting the num of HDUs in file" , infile) ;
         if (i==0){
            copyUsrkeywrdsTovect (fptr,key_records);                        
        }
         long numrows = 0 ;
        fits_movabs_hdu (fptr , 2 , NULL , &status) ;
        printError (status , "Error in Moving to  2nd HDU" , infile) ;
        fits_get_num_rows (fptr , &numrows , &status) ;
        printError (status , "Error in getting the number of HDUs" , infile) ;
        fits_read_key (fptr , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ;
        printError (status , "FRAMENO keyword not Found") ;
        fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
        printError (status , "FRMTIME keyword not Found",infile) ;
        fits_movabs_hdu (fptr , 1 , NULL , &status) ;
        printError (status , "Error in moving to 2nd HDU" , infile) ;
      
        /**creating the output event file**/
        sprintf (outfile , "%s/%s/%s_t%.4f_f%d_centroid_dd.fits" , moduleoutdir , centroidDir , nameprefix , frametime , frameno) ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating  the output Centroid frame " , outfile) ;
       
         
        /**coping the  input event file to the output event file**/
        fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
        printError (status , "Error in coping the file" , infile) ;

        float *Xloc = new float[numrows] ;
        float *Yloc = new float[numrows] ;

        fits_movabs_hdu (fout , 2 , NULL , &status) ;
        printError (status , "Error in Moving to  2nd HDU" , outfile) ;
        fits_read_col (fout , TFLOAT , 1 , 1 , 1 , numrows , NULL , (void *) Xloc , NULL , &status) ;
        printError (status , "Error in Reading the Xloc " , outfile) ;
        fits_read_col (fout , TFLOAT , 2 , 1 , 1 , numrows , NULL , (void *) Yloc , NULL , &status) ;
        printError (status , "Error in reading the Yloc" , outfile) ;
        /*Correcton on Centroid based on caldb*/
        float tempx , tempy ;
        double locate =0.0f;
        float multi_factor = xsize / PADDED_FRAMESIZE ; //if the frame is subdivided, then distortion file also needs to be subdivided
        if (multi_factor == 0)
        {
            LOG(ERROR) << "***Divide by Zero***" << endl ;
            return (EXIT_FAILURE) ;
        }
        /**Applying Correction  on xloc and yloc**/
        for (int i = 0 ; i < numrows ; i++)
        {
            if(Xloc[i]!=INVALID_PIX_VALUE && Yloc[i]!=INVALID_PIX_VALUE){
            tempx = Xloc[i] / multi_factor ; //taking dimension in unpadded frame
            tempy = Yloc[i] / multi_factor ; //taking dimension in unpadded frame
            locate = ((int) round (tempy) - 44) *512 + ((int) round (tempx) - 44) ; //44 is subtracted as frame size was padded to 600x600 from initial 512x512
            round (locate) ;
          
            tempx = tempx + x_Distortion[(int) locate] ;
            tempy = tempy + y_Distortion[(int) locate] ;
           if(tempx*multi_factor <xsize && tempy*multi_factor<ysize) 
           { 
            Xloc[i] = tempx * multi_factor ;
            Yloc[i] = tempy * multi_factor ;
           }
         }
        }
        /*corrected Centroid written to the output File.*/
        fits_write_col (fout , TFLOAT , 1 , 1 , 1 , numrows , (void *) Xloc , &status) ;
        printError (status , "Error in reading the Xloc " , outfile) ;
        fits_write_col (fout , TFLOAT , 2 , 1 , 1 , numrows , (void *) Yloc , &status) ;
        printError (status , "Error in reading the Yloc " , outfile) ;
  
        //writeUsrkeywordsFrmvect (outfile,key_records);
        if (history == YES) writeHistory (outfile , vhistorystr) ;
         writeCommonKeywords (fout,modulename);
        
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output Centroid file" , outfile) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the input Centroid file" , outfile) ;
    
        strcpy (outcentroidframelist[i] , basename (outfile)) ;
        delete[] Xloc , Yloc ;
        cout<<"Total Files written=" << i << "     Remaining files= " << nframes - i <<" \r" ;
    }
    /*centroid File list written to the INFO file*/
    LOG(INFO)<<"\nUpdating the output information file..."<<endl;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the information file",infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Moving to 2nd HDU in out Info File Fails",infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , nframes , (void *) outcentroidframelist , &status) ;
    printError (status , "Writing  centroid List for output INFO File Fails",infofile_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the file",infofile_out) ;
    writeUsrkeywordsFrmvect(infofile_out,key_records);
     if(history==YES)writeHistory (infofile_out , vhistorystr) ;
    return (EXIT_SUCCESS) ;
}

int uvtDetectDistCorr::detectDistortionPC ()
{
    LOG(INFO)<<"\nStarted  correction for  detector Distortion process for PC mode "<<endl;
    char file_in[NAMESIZE] , file_out[NAMESIZE] ; //for input and output event file
     /**Setting path for input event file**/
    sprintf (file_in , "%s/%s" , inputdatadir , eventfile) ; //taking event file full path
    /**setting path for output event **/
    sprintf (file_out , "%s/%s_dd.events" , moduleoutdir , nameprefix) ;
    fitsfile *fevt_in , *fevt_out , *finfo_out ;
    int status = 0 ;
    /**opening the input Event file**/
    fits_open_file (&fevt_in , file_in , READONLY , &status) ;
    printError (status , "Error in opening the input event file" , file_in) ;
    copyUsrkeywrdsTovect (fevt_in,key_records);
    fits_create_file (&fevt_out , file_out , &status) ;
    printError (status , "Error in creating the output information file") ;
   /**Coping the input Event file to output Event file**/
    fits_copy_file (fevt_in , fevt_out , 1 , 1 , 1 , &status) ;
    printError (status , "Error in coping the input event file to output event file" , file_out) ;
    /***Defining pointers for reading a event file columns***/
    float *xFrac , *yFrac ;
    long nrows=0 ;
    fits_movabs_hdu (fevt_out , 2 , NULL , &status) ;
    printError (status , "Error in  moving to 2nd HDU in output event file" , file_out) ;
    fits_get_num_rows (fevt_out , &nrows , &status) ;
    printError (status , "Error in  getting the number of rows in output event file" , file_out) ;
    /***Assigning a memory to Arrays***/
    xFrac = new float[nrows] ;
    yFrac = new float[nrows] ;
    fits_read_col (fevt_out , TFLOAT , 4 , 1 , 1 , nrows , NULL , xFrac , NULL , &status) ;
    printError (status , "Error in reading the xFrac " , file_out) ;
    fits_read_col (fevt_out , TFLOAT , 5 , 1 , 1 , nrows , NULL , yFrac , NULL , &status) ;
    printError (status , "Error in  reading the yFrac" , file_out) ;
    float locate ;
    float tempx , tempy , multi_factor = 0.0 ;
    multi_factor = xsize / 600 ;
     LOG(INFO)<<"Applying Correction..."<<endl;
//     cout<<"The caldbDim"<<caldbdim<<endl;exit(1);
      for (int k = 0 ; k < nrows ; k++)
    {
        tempx = xFrac[k] / multi_factor ;
        tempy = yFrac[k] / multi_factor ;

        locate = ((int) round (tempy) - 44) * 512 + ((int) round (tempx) - 44) ;
        round (locate) ;

        tempx = tempx + x_Distortion[(int) locate] ;
        tempy = tempy + y_Distortion[(int) locate] ;

     if(tempx*multi_factor <xsize && tempy*multi_factor<ysize) 
     { 
        xFrac[k] = tempx * multi_factor ;
        yFrac[k] = tempy * multi_factor ;
      }
    }
    /***Writing to output event file***/
    fits_write_col (fevt_out , TFLOAT , 4 , 1 , 1 , nrows , xFrac , &status) ;
    printError (status , "Error writing Xfractional",file_out) ;
    fits_write_col (fevt_out , TFLOAT , 5 , 1 , 1 , nrows , yFrac , &status) ;
    printError (status , "Error writing Yfractional",file_out) ;
    /*releasing the memory*/
    delete[] xFrac ;
    delete[] yFrac ;
    fits_movabs_hdu (fevt_out , 1 , NULL , &status) ;
    printError (status , "Moving to 1st HDU fails for output event file ",file_out) ;
    /***Writes the  common information to  output information File***/
   
    vector<string> vhistorystr ;
    if (history == YES)
    {
        getHistory (vhistorystr) ;
        writeHistory (file_out , vhistorystr) ;
    }
     writeCommonKeywords (fevt_out , modulename) ;
    fits_close_file (fevt_out , &status) ;
    printError (status , "Error in Closing the output Event File",file_out) ;
    fits_close_file (fevt_in , &status) ;
    printError (status , "Error in closing the input Event File",file_out) ;
    LOG(INFO)<<"\nUpdating information to output information File"<<endl;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the output information File",file_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in opening the output information File",file_out) ;
    fits_update_key (finfo_out , TSTRING , "EVTFILE" , basename (file_out) , NULL , &status) ;
    printError (status , "Error in updating the key value for the EVTFILE",file_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error  closing the output information file",file_out) ;
    writeUsrkeywordsFrmvect(infofile_out,key_records);
    if(history==YES)writeHistory (infofile_out , vhistorystr) ;
    return (EXIT_SUCCESS) ;
}

int uvtDetectDistCorr::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" caldbDir=" + (string) caldbDir) ;
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
//int uvtDetectDistCorr::readCaldbDistoFile(float *Xdist, float *Ydist, char *distrotionCorrfile, long &dim) {
//
//    this->x_Distortion = Xdist;
//    this->y_Distortion = Ydist;
//    strcpy(this->distortionCorrfile, distrotionCorrfile);
//    this->caldbdim = dim;
//    readDistortionFile();
//    strcpy(distrotionCorrfile, this->distortionCorrfile);
//    dim = caldbdim;
//    // Xdist=this->x_Distortion;
//    //Ydist=this->y_Distortion;
//    LOG(INFO) << "The dim is " << caldbdim;
//    //  Xdist= new float[caldbdim*caldbdim];
//    // Ydist= new float[caldbdim*caldbdim];
//    for (int i = 0; i < caldbdim * caldbdim; i++) {
//        Xdist[i] = 0.0;
//        Ydist[i] = 0.0;
//        Xdist[i] = x_Distortion[i];
//        Ydist[i] = y_Distortion[i];
//    }
//    //  Xdist=this->x_Distortion;
//    // Ydist=this->y_Distortion;
//    //dim=caldbdim;
//    LOG(INFO) << "The dim" << dim << endl;
//   
//    
//    return (EXIT_SUCCESS);
//}
