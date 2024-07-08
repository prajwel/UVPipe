/* 
 * File:   uvtFindWtdMean.cpp
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
#include <uvtFindWtdMean.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<glog/logging.h>

//Constructor -called when object is created

uvtFindWtdMean::uvtFindWtdMean ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    strcpy (centroidDir , "Centroid") ;
}

//Destructor

uvtFindWtdMean::~uvtFindWtdMean ()
{
    
}

//parameter File reading

int uvtFindWtdMean::read (int argc , char** argv)
{
    int status = 0 ;

    if (PIL_OK != (status = PILInit (argc , argv)))
    {

        LOG(ERROR) << "***Error Initializing PIL***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("framelistDir" , inputdatadir)))
    {
        LOG(ERROR) << endl << "***Error reading framelist file name***" ;
        return status ;
    }
     if (PIL_OK != (status = PILGetInt ("no_of_weightedframes" , &no_ofWeigh)))
    {
        LOG(ERROR) << endl << "***Error reading Caldb  Directory Path ***" ;
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

int uvtFindWtdMean::read (char *inputdatadir , char *outdir , int num_ofweight , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;
    //    strcpy (this->caldbDir , caldbDir) ;
    this->no_ofWeigh = num_ofweight ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}
//Parameter file content Display

void uvtFindWtdMean::display ()
{
    LOG(INFO) << endl << "----------Finding Weighted Mean Parameter Display---------" ;
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

    LOG(INFO) << endl << "----------Finding Weighted Mean parameter  Display Ends----------\n" ;
}

//Correction for the  Cosmic Ray process

int uvtFindWtdMean::uvtFindWeightedMeanProcess ()
{
    /*if the Input Directory of the Filelist is not available then exit from module
     else read Filelist
     */
    LOG(INFO) << endl << "Started Find Weighted Mean process.." << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    //LOG(INFO)<<endl<<"Module Output Directory : "<<moduleoutdir<<endl;
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
    /**Shell command for creating the output Directory**/
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
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists (infofile_in)))
    {
        LOG(ERROR) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY" ;
        return (EXIT_FAILURE) ;
    }
    LOG(INFO) << endl << "\nInformation File :" << infofile_in ;
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "***Error in opening the information file***") ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "***Error in Moving the 2nd HDU***") ;
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "***NAMEPRFX keyword not Found***") ; //for creating name for output information file
    sprintf (infofile_out , "%s/%s_wm.info" , moduleoutdir , nameprefix) ;
    LOG(INFO) << "\nOutput information file : " << infofile_out << endl ;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "***Error in creating the output information file***") ;
    char *ttype[] = {"SignalFrames" , "ExposureFrames"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table",infofile_out) ;
    datainfo.write (finfo_out) ; //writing basic data information
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the key value of NAMEPRFX",infofile_out) ; //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    fits_close_file (finfo_out , &status) ;
    if (datainfo.getModeFlag () == IM)
    {//For IM mode
        fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "NFILES not Found" , infofile_in) ;
        fits_read_key (finfo_in , TSTRING , "SIGDIR" , sigframedir , NULL , &status) ;
        printError (status , "SIGDIR  not Found" , infofile_in) ;
        fits_read_key (finfo_in , TSTRING , "EXPDIR" , expoframedir , NULL , &status) ;
        printError (status , "EXPDIR  not Found" , infofile_in) ;

        signalframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        exposureframelist = allocateMemory<char>(nframes , NAMESIZE) ;

        fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) signalframelist , NULL , &status) ;
        fits_read_col (finfo_in , TSTRING , 2 , 1 , 1 , nframes , NULL , (void *) exposureframelist , NULL , &status) ;

        //method for the Detector Distortion  For IM
        if (findWtdMeanIM ()) return (EXIT_FAILURE) ;
        fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
        printError (status , "Error in opening the out  information file" , infofile_out) ;
        fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
        printError (status , "Error in moving the 2nd HDU of  information file" , infofile_out) ;
        fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
        printError (status , "NAMEPRFX keyword  not updated/not Found" , infofile_out) ; //for creating name for output information file
        fits_update_key (finfo_out , TINT , "NFILES" , &index1 , NULL , &status) ;
        printError (status , "NFILES keyword  not updated/not Found" , infofile_out) ;
        fits_close_file (finfo_out , &status) ;
        printError (status , "Error in closing the file" , infofile_out) ;
        freeMemory (signalframelist , nframes , NAMESIZE) ; //for releasing the memory
        freeMemory (exposureframelist , nframes , NAMESIZE) ; //for releasing the memory
    }
    else
    {
        LOG(ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG(ERROR) << endl << "Allowed value is  Only IM" ;
        return (EXIT_FAILURE) ;
    }
    fits_close_file(finfo_in,&status);
    printError (status , "Error in closing the file" , infofile_in) ;
    return (EXIT_SUCCESS) ;
}

int uvtFindWtdMean::findWtdMeanIM ()
{
    LOG(INFO) << "\nStarted weighted mean process  for IM mode.." << endl ;
    char **outsigframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    char **outexpframelist = allocateMemory<char>(nframes , NAMESIZE) ; //to store output exposure frame list  
    char infile[NAMESIZE] , outfile[NAMESIZE] ;
    float framedata[ysize * xsize] ;
    int status = 0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    fitsfile *fptr , *fout ;
    float sum_sig ;
    float sum_exp ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s" , moduleoutdir , sigframedir) ;
    /**Shell command for creating the directory**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    //creating exposure directory
    sprintf (dir , "%s/%s" , moduleoutdir , expoframedir) ;
    cmd = "mkdir -p " + (string) dir ;
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    vector<string> vhistorystr ;
    //reading history
     if (history == YES)
        getHistory (vhistorystr) ;
    int image_no ;
    int cnt_frames = 0 ;
    int no_of_images = nframes ;
    index1 = 0 ;
    double frametime;
    float *outputArray_sig , *outputArray_exp ;
    vector<double> time_track;
    outputArray_sig = new float [xsize * ysize] ;
    outputArray_exp = new float [xsize * ysize] ;
    if(no_ofWeigh>nframes)
    {
         LOG(INFO) << "\nnumber of frames to be weighted is greater than total number of frames - " <<  endl ;
         return(EXIT_FAILURE);
    }
    LOG(INFO) << "\nTotal number of frames - " << nframes << endl ;
    LOG(INFO) << "\nPerforming  Weighted Mean of frames...." << endl ;
    
    //for total number of frames to be weighted
    for (int x = 0 ; x < nframes ; x = x + no_ofWeigh)
    {
        cout<<x<<endl;
        for(int i=0;i<xsize*ysize;i++)
        {
            outputArray_sig[i]=-9999;
            outputArray_exp[i]=-9999;
         }
         image_no = 0 ;
        cnt_frames++ ;
        float ***tempSignalArray , ***tempExposureArray ;
        tempSignalArray = new float**[no_of_images] ;
        tempExposureArray = new float**[no_of_images] ;
        if (nframes - x < no_ofWeigh)
        {
            break;
            //no_ofWeigh = nframes - x ;
        }
//        
        //for allocating the memory
        for (int i = 0 ; i < no_ofWeigh ; i++)
        {
            tempSignalArray[i] = new float*[xsize] ;
            tempExposureArray[i] = new float*[xsize] ;

            for (int j = 0 ; j < xsize ; j++)
            {
                tempSignalArray[i][j] = new float[ysize] ;
                tempExposureArray[i][j] = new float[ysize] ;
            }
        }
        
        //allocating the memory to output array.
        float **opSigArray ;
        opSigArray = new float *[xsize] ;
        for (int j = 0 ; j < xsize ; j++)
        {
            opSigArray[j] = new float [ysize] ;
        }
        float **opExpArray ;
        opExpArray = new float *[xsize] ;
        for (int j = 0 ; j < xsize ; j++)
        {
            opExpArray[j] = new float [ysize] ;
        }
        time_track.clear ();
        //loop for storing the pixel  of total number of weighted frame to Array i.e tempSignalArray and  tempExposureArray  .
        for (int i = x ; i < x + no_ofWeigh ; i++)
        {
            for (int q = 0 ; q < xsize * ysize ; q++)
            {
                framedata[q] = 0.0 ;
            }
            sprintf (infile , "%s/%s/%s" , inputdatadir , "SignalFrames" , signalframelist[i]) ;
            float *oneDimSignalArray , *oneDimExposureArray ;
            oneDimSignalArray = new float [xsize * ysize] ;
            oneDimExposureArray = new float [xsize * ysize] ;
            for (int q = 0 ; q < xsize * ysize ; q++)
            {
                oneDimSignalArray[q] = 0.0 ;
                oneDimExposureArray[q] = 0.0 ;
            }
            fits_open_file (&fptr , infile , READONLY , &status) ;
            printError (status , "Error in opening the input File" , infile) ;
            fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , oneDimSignalArray , NULL , &status) ;
            printError (status , "Error in reading the pixels from the input File" , infile) ;
            fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
            printError (status , "Error reading the key value of the FRAMRTIME" , infile) ;
            time_track.push_back (frametime);
            copyUsrkeywrdsTovect (fptr,key_records);  
         
            fits_close_file (fptr , &status) ;
            printError (status , "Error in closing the file" , infile) ;
            
            sprintf (infile , "%s/%s/%s" , inputdatadir , "ExposureFrames" , exposureframelist[i]) ;
            fits_open_file (&fptr , infile , READONLY , &status) ;
            printError (status , "Error in opening the input File" , infile) ;
            fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , oneDimExposureArray , NULL , &status) ;
            printError (status , "Error in reading the pixels from the input File" , infile) ;
            fits_close_file (fptr , &status) ;
           
            int temp_index_in = 0 ;
            for (int d1 = 0 ; d1 < xsize ; d1++)
            {
                for (int d2 = 0 ; d2 < ysize ; d2++)
                {
                    tempSignalArray[image_no][d1][d2] = oneDimSignalArray[temp_index_in] ;
                    tempExposureArray[image_no][d1][d2] = oneDimExposureArray[temp_index_in] ;
                    temp_index_in++ ;
                }
            }
            image_no++ ;
            free (oneDimSignalArray) ;
            free (oneDimExposureArray) ;
        }
        double time_middle=(time_track[0]+time_track[time_track.size ()-1])/2;
        long temp_index_out = 0 ;
        for (int ii = 0 ; ii < xsize ; ii++)
        {
            for (int jj = 0 ; jj < ysize ; jj++)
            {
                sum_sig = 0;
                sum_exp =0 ;
                for (int kk = 0 ; kk < no_ofWeigh ; kk++)
                {
                    if(tempSignalArray[kk][ii][jj]!=INVALID_PIX_VALUE  && tempExposureArray[kk][ii][jj] !=INVALID_PIX_VALUE)
                    {
                          if (kk == 0)
                    {
                        sum_sig = tempSignalArray[kk][ii][jj] ;
                        sum_exp = tempExposureArray[kk][ii][jj] ;
                    }
                    else
                    {
                        sum_sig = (sum_sig + (tempSignalArray[kk][ii][jj]) * tempExposureArray[kk][ii][jj]) / 2.0 ;
                        sum_exp = sum_exp + tempExposureArray[kk][ii][jj] ;
                    }
                    }
                  
                }
                opSigArray[ii][jj] =-9999;
                opExpArray[ii][jj]=-9999;
                opSigArray[ii][jj] = sum_sig ;
                opExpArray[ii][jj] = sum_exp ;
                outputArray_sig[temp_index_out]=-9999;
                outputArray_exp[temp_index_out]=-9999;
                outputArray_sig[temp_index_out] = opSigArray[ii][jj] ;
                outputArray_exp[temp_index_out] = opExpArray[ii][jj] ;
                temp_index_out++ ;
            }
        }
        for (int j = 0 ; j < ysize ; j++)
        {
            free (opSigArray[j]) ;
            free (opExpArray[j]) ;
        }
        
        free (opSigArray) ;
        free (opExpArray) ;
        for (int i = 0 ; i < no_ofWeigh ; i++)
        {
            for (int j = 0 ; j < xsize ; j++)
            {
                free (tempSignalArray[i][j]) ;
                free (tempExposureArray[i][j]) ;
            }
        }
        free (tempSignalArray) ;
        free (tempExposureArray) ;
        
        //setting path of output  signal frame.
        sprintf (outfile , "%s/%s/%s_t%f_f%d_wm_sig.fits" , moduleoutdir , sigframedir , nameprefix , time_middle,cnt_frames) ;
        strcpy (outsigframelist[index1] , basename (outfile)) ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the out file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the IMG file" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , outputArray_sig , &status) ;
        printError (status , "Error in writing the pixels to the output file" , outfile) ;           
        fits_close_file (fout , &status) ;
      
          fits_open_file (&fout , outfile , READWRITE , &status) ;
          printError (status , "Error in opening the input File" , outfile) ;
          fits_update_key (fout , TDOUBLE , "FRMTIME" , &time_middle , NULL , &status) ;//updating frame time 
         printError (status , "Error in writing the key value of the FRMTIME" , outfile) ;
          fits_update_key (fout , TUSHORT , "FRAMENO" , &cnt_frames , NULL , &status) ;//updating frame number.
         printError (status , "Error in writing the key value of the FRMENO" , outfile) ;
         fits_close_file (fout , &status) ;
         printError (status , "Error in closing the output signal file" , outfile) ;  
        
        if(history==YES){
               writeHistory(outfile , vhistorystr) ; //write history to each file
        }
         writeUsrkeywordsFrmvect (outfile,key_records);
        updateKeywords (outfile,modulename);
        sprintf (outfile , "%s/%s/%s_t%f_f%d_wm_exp.fits" , moduleoutdir , expoframedir , nameprefix , time_middle,cnt_frames) ;
        strcpy (outexpframelist[index1] , basename (outfile)) ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the out file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the IMG file" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , outputArray_exp , &status) ;
        printError (status , "Error in writing the pixels to the output file" , outfile) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the File" , infile) ;
        
        
        fits_open_file (&fout , outfile , READWRITE , &status) ;
        printError (status , "Error in opening the input File" , outfile) ;       
        fits_update_key (fout , TDOUBLE , "FRMTIME" , &time_middle , NULL , &status) ;
        printError (status , "Error in writing the key value of the FRMTIME" , outfile) ;
        fits_update_key (fout , TUSHORT , "FRAMENO" , &cnt_frames , NULL , &status) ;
        printError (status , "Error in writing the key value of the FRMENO" , outfile) ;
        fits_close_file (fout , &status) ;
         printError (status , "Error in closing the output file" , outfile) ;
        
         if(history==YES)
               writeHistory(outfile , vhistorystr) ; //write history to each file
        
         writeUsrkeywordsFrmvect (outfile,key_records);
        updateKeywords (outfile,modulename);
        index1++ ;
      }
    vhistorystr.clear () ;
    //adding framelist to info file
    LOG(INFO) << "\nWriting list of output file namelist to output information file " << endl ;
    fitsfile *finfo_out ;
    
    //updating the output information file.
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the out information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU" , infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , cnt_frames , (void *) outsigframelist , &status) ;
    printError (status , "Error in writing the  output Signal framelist to the  out information file" , infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 2 , 1 , 1 , cnt_frames , (void *) outexpframelist , &status) ;
    printError (status , "Error in writing the  output Exposure framelist to the  out information file " , infofile_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in output information File" , infofile_out) ;
    cout<<"OURCOFLMF"<<endl;
    freeMemory (outsigframelist , nframes , NAMESIZE) ;
    freeMemory (outexpframelist , nframes , NAMESIZE) ;
    return (EXIT_SUCCESS) ;
}

int uvtFindWtdMean::getHistory (vector<string> &vhistory)
{
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ("P1 framelistDir=" + (string) inputdatadir) ;
    vhistory.push_back ("P3 outdir=" + (string) outdir) ;
    vhistory.push_back ("Module Output directory=" + (string) moduleoutdir) ;
    vhistory.push_back ("Centroid Directory=" + (string) centroidDir) ;
    if (clobber == YES)
        vhistory.push_back ("P4 clobber=yes") ;
    else
        vhistory.push_back ("P4 clobber=no") ;
    if (history == YES)
        vhistory.push_back ("P5 history=yes") ;
    else
        vhistory.push_back ("P5 history=no") ;
    vhistory.push_back ("Parameter List END") ;
    return (EXIT_SUCCESS) ;
}

