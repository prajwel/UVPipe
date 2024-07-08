/* 
 * File:   uvtAccEveryTsec.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include<uvtAccEveryTsec.h>
#include <algorithm>
#include<glog/logging.h>
#include<fitsio.h>
#include<pil.h>
#include<map>
#include<macro_def.h>

//#include "uvtDetectDistCorr.h"


uvtAccEveryTsec::uvtAccEveryTsec ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}


uvtAccEveryTsec::~ uvtAccEveryTsec () {
 }


int uvtAccEveryTsec::read (int argc , char** argv)
{

    int status = 0 ;
    status = readParams (argc , argv , 6 , FNAME , "inputdatadir" , inputdatadir , INT , "Nacc" , &numOfFramesToAcc , FNAME , "outdir" , outdir , BOOL , "clobber" , &clobber ,
            BOOL , "history" , &history , STRING , "mode" , &mode) ;

    if (status) return (EXIT_FAILURE) ;

    return (EXIT_SUCCESS) ;
}


int uvtAccEveryTsec::read (char* inputdatadir , char* outdir , int Nacc , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;
    this->numOfFramesToAcc = Nacc ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}


void uvtAccEveryTsec::display ()
{
    LOG (INFO) << endl ;
    LOG (INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG (INFO) << "             UVT ACC EVERY T SEC PARAMETERS       " << endl ;
    LOG (INFO) << "------------------------------------------------------------------------" ;
    LOG (INFO) << endl << "Input Data Directory    : " << inputdatadir ;
    LOG (INFO) << endl << "Output Directory   : " << outdir ;
    LOG (INFO) << endl << "Number of frames to be Accumulated   : " << numOfFramesToAcc ;
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


int uvtAccEveryTsec::uvtAccEveryTsecProcess ()
{
    sprintf (moduleoutdir , "%s/%s" , outdir , modulename) ;
    LOG (INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;

    //check  existence  of output directory if output  exist and 
    //clobber =YES than remove that directory and recreate it,if exist and clobber =NO than exit from module.
    //if directory is not exist than create it .
    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;

    LOG (INFO) << endl << moduleoutdir << "  directory created" ;

    //opening info file in input directory to get data information
    string tempfilepath = searchFile (inputdatadir , ".info") ; //check for the .info file in the perticuler Directory
    if (tempfilepath == " ")
    {
        LOG (INFO) << endl << "Error in finding info file" ;
        return (EXIT_FAILURE) ;
    }
    /**setting the path of input information file**/
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    LOG (INFO) << endl << "\nInput Information File :" << infofile_in ;
    int status = 0 ;

    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the  input information file " , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU in input information file" , infofile_in) ;
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "NAMEPRFX  keyword not Found" , infofile_in) ; //for creating name for output information file
    fits_close_file (finfo_in , &status) ;
    printError (status , "Error in closing te input information file" , infofile_in) ;
    
    if (strcmp(datainfo.getObsMode (),"PC")==0)
    {
        LOG (INFO) << endl << "This operation is not valid for PC mode data\nExiting" << endl ;
        return (EXIT_FAILURE) ;
    }
    xsize = datainfo.getXsize () ; //xsize-x dimention of the image
    ysize = datainfo.getYsize () ; //ysize-y dimention of the image

    LOG (INFO) << endl << "Xsize :" << xsize << "  Ysize:" << ysize ;
  

    readKeywords (infofile_in , 2 , 3 , TSTRING , "SIGDIR" , sigframedir ,
            TSTRING , "EXPDIR" , expframedir ,
            TINT , "NFILES" , &nframes) ;

    printError (status , "Error in reading the key value of the NFILES " , infofile_in) ;
    LOG (INFO) << endl << "Number of frames :" << nframes ;
  
    if(numOfFramesToAcc>nframes){
        LOG(ERROR)<<"***Number of frames to Accumulate is greater than  available frames***"<<endl;
        return(EXIT_FAILURE);
    }
    //computing number of output frames
    //numOfFramesToAcc=int(timebin/datainfo.getIntegrationTime());
    if (numOfFramesToAcc == 0)
    {
        LOG (INFO) << "***Divide by zero***" << endl ;
        return (EXIT_FAILURE) ;
    }
    /**
     * calculating total  number of frames to be generated after Accumulation
     * nframes_out-total number of frames to be generated at output after accumulation
     * @return 
     */
    int temp_numOFframes =numOfFramesToAcc;
   // if (nframes % numOfFramesToAcc == 0)
        nframes_out = nframes / numOfFramesToAcc ;
        int temp_numFramesout=nframes_out;
//    else
//        nframes_out = (nframes / numOfFramesToAcc) + 1 ;
    LOG (INFO) << endl << "Number of frame to accumulate :" << numOfFramesToAcc ;
    LOG (INFO) << endl << "Number of output frames : " << nframes_out ;
    LOG (INFO) << endl << "File Time duration: " << datainfo.getTimeDuration () << "  sec" ;

    /*
     *sigframelist-array for storing input signal frame names
     * expframelist-array for storing input exposure frame 
     *     
     *  */
    sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
    expframelist = allocateMemory<char>(nframes , NAMESIZE) ;

    status = readColumnsFromFITS (infofile_in , 2 , 2 , TSTRING , 1 , sigframelist , nframes , TSTRING , 2 , expframelist , nframes) ;
    if (status)
    {
        LOG (INFO) << "Error reading  the columns from the file" << endl ;
        return (EXIT_FAILURE) ;
    }

    //LOG(INFO)<<"NFRAMES "<<nframes<<endl;
    //creating output information file
    sprintf (infofile_out , "%s/%s_ac.info" , moduleoutdir , nameprefix) ;
    LOG (INFO) << "\nOutput Informaton File " << infofile_out << endl ;
    LOG (INFO) << "\nCreating file  " << infofile_out ;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information file") ;
    char *ttype[] = {"SignalFileList" , "ExposureFileList"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the output information file") ;
    //cloumn values will be written in accumulateFrames functoin
    datainfo.write (finfo_out) ; //writing basic data information
    updateKeywords (infofile_out , 2 , 4 , TSTRING , "SIGDIR" , sigframedir , TSTRING , "EXPDIR" , expframedir , TINT , "NFILES" , &nframes_out,TSTRING,"NAMEPRFX",nameprefix) ; //updating the keywords

    //LOG(INFO)<<endl<<"Closing information file created";
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file") ;

    double int_time=datainfo.getIntegrationTime ();
    
    flag_Acc=FALSE;
    if(int_time>=1){
        flag_Acc=TRUE;        
    }
    LOG (INFO) << "\nAccumulating Signal frames........." ;
    /**
     * method for accumulate  frames
     * sigframedir-input  signal  frame directory
     * expframedir-input exposure frame directory
     * @return 
     */
    accumulateFrames (sigframedir , sigframelist , "sig") ;

    numOfFramesToAcc=temp_numOFframes;
    nframes_out=temp_numFramesout;
    LOG (INFO) << "\nAccumulating Exposure frames.........." ;
  
    accumulateFrames (expframedir , expframelist , "exp") ;

    LOG (INFO) << endl << "Creating frame time file....." ;
   
    freeMemory (sigframelist , nframes , NAMESIZE) ;
    freeMemory (expframelist , nframes , NAMESIZE) ;

    return (EXIT_SUCCESS) ;
}


int uvtAccEveryTsec::getHistory (vector<string>& vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    char s_numOffrmToAcc[10] ;
    sprintf (s_numOffrmToAcc , "%d" , numOfFramesToAcc) ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Number of Frames to Accumulate=" + (string) s_numOffrmToAcc) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
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


/*function to accumulate frames 
dir - name of directory for  containing  output frames
namelist - List of all the input frames
frameidentifier - it  is   "sig"  for Signal Frame and "exp"  for Exposure Frame*/

int uvtAccEveryTsec::accumulateFrames (char *dir , char **namelist , char *frameIdentifier)
{
    cout.flush ();
    string cmd = "mkdir -p " + (string) moduleoutdir + (string) "/" + (string) dir ;
    system (cmd.c_str ()) ;
    char **outlist;
    if (flag_Acc==TRUE){
        outlist = allocateMemory<char>(nframes-numOfFramesToAcc+1 , NAMESIZE) ; // to hold ouput filenames
    }
    else {
        outlist = allocateMemory<char>(nframes_out , NAMESIZE) ; // to hold ouput filenames
    }
    
    char outfile[FLEN_FILENAME] ;
    char infile[FLEN_FILENAME] ;
    vector<string> vhistorystr ;
    if (history == YES) getHistory (vhistorystr) ;

    fitsfile *fptr , *fout ;
    int status = 0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    float framedata[xsize * ysize] ;
    //double sum[xsize * ysize] ;
    float sum[xsize * ysize] ;
    float average[xsize * ysize] ;

    int inputframeindex = 0 ;
    int count_track = 0 ;
    double frametime = 0 , avgftime = 0 ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    int frameno ;
    LOG (INFO) << endl << "Loop start to accumulate frames " << endl ;
    vector<int> track_invalidPix ; //vector for storing  position of  invalid pixel 
    int tempx=0;int tempy=0;
    //loop for number of frames to be generated at output
  
        
    int temp_Acc=numOfFramesToAcc;
           // if (inputframeindex >= nframes_out) break ;
    for (int i = 0 ; i < nframes_out ; i ++)
    {
       
       // initArray<double>(sum , xsize*ysize , 0.0f) ;
        initArray<float>(sum , xsize*ysize , 0.0f) ;
        initArray<float>(average , xsize*ysize , -9999) ;

        frametime = 0 ;
        avgftime = 0 ;
        count_track = 0 ;
        //loop for  accumulating  number of frames 
         for (int k =0 ; k <xsize ; k++)
            {
             for(int l=0;l<ysize;l++)
             { 
                // int ghj=l*padding_dim+k;
                 track_invalidPix.push_back (0);
                 //invalid_pix_duplicates_indices.insert (pair(ghj, 0));
             }
            }
       
        for (int j = 0 ; j < numOfFramesToAcc ; j ++)
        {
            count_track ++ ; //counter for number of frames averaged, will be different in last loop
            sprintf (infile , "%s/%s/%s" , inputdatadir , dir , namelist[inputframeindex]) ;
            fits_open_file (&fptr , infile , READONLY , &status) ;
            printError (status , "Error in opening the  fits file" , infile) ;
            if (j == 0)
                copyUsrkeywrdsTovect (fptr , key_records) ;
            fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
            printError (status , "Error  in reading the pixels from the input fits file" , infile) ;
            fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
            printError (status , "Error reading the key value of the FRMTIME " , infile) ;
            fits_close_file (fptr , &status) ;
            printError (status , "Error in closing the file" , infile) ;
            avgftime = avgftime + frametime ;
            /**loop for getting track of invalid pixels from numOfFramesToacc frames **/
            for (int k =tempx= (xsize-(xsize*IMG_DIM_DI)/PIX_PADD_SIZE)/2 ; k < tempx+IMG_DIM_DI*xsize/PIX_PADD_SIZE ; k++)
            {
             for(int l=tempy=(ysize-(ysize*IMG_DIM_DI)/PIX_PADD_SIZE)/2;l<tempy+IMG_DIM_DI*xsize/PIX_PADD_SIZE;l++){
                 if(framedata[k*xsize+l]!=INVALID_PIX_VALUE)
                   sum[k*xsize+l] = sum[k*xsize+l] + framedata[k*xsize+l] ;
                 else track_invalidPix[k*xsize+l]= track_invalidPix[k*xsize+l]+1;
             }
            }
//            for (int k = 0 ; k < (xsize * ysize) ; k ++)
//            {
//                if (framedata[k] != INVALID_PIX_VALUE)//check whether pixel is identified invalid from previous modules or not,if yes than do nothing to it
//                    sum[k] = sum[k] + framedata[k] ;
//                else
//                {
//                    track_invalidPix.push_back (k) ;
//                }
//            }
            inputframeindex ++ ; //index for total frames input
            if (inputframeindex >= nframes) break ;
        } //end of j loop for numOfFrameToAcc
      
        
        if(flag_Acc==TRUE)
        {
            nframes_out=nframes-temp_Acc+1;
            numOfFramesToAcc=1;
            
           // cout<<"sgf "<<nframes_out<<endl;exit(1);
        }
        
        if (count_track == 0)
        {
            LOG (ERROR) << "***Divide by zero ***" << endl ;
            return (EXIT_FAILURE) ;
        }
        int num_times_repeat = 0 ;
         for (int k = tempx=(xsize-(xsize*IMG_DIM_DI)/PIX_PADD_SIZE)/2 ; k < tempx+IMG_DIM_DI*xsize/PIX_PADD_SIZE ; k++)
         {
             for(int l=tempy=(ysize-(ysize*IMG_DIM_DI)/PIX_PADD_SIZE)/2;l<tempy+IMG_DIM_DI*xsize/PIX_PADD_SIZE;l++)
             {
//                     num_times_repeat = count (track_invalidPix.begin () , track_invalidPix.end () , k*xsize+l) ; //STL method for count
//                     average[k*xsize+l] = (float) (sum[k*xsize+l] / (count_track - num_times_repeat)) ;
                  average[k*xsize +l]=INVALID_PIX_VALUE;
                    if(track_invalidPix[k*xsize+l] == numOfFramesToAcc) 
                    {
                        average[k*xsize+l]=INVALID_PIX_VALUE;
                    }else{
                 //   frmsigdata[k*sizex+l] = (float) (sumdata[k*sizex+l] / (numoffrmAcc- num_times_repeat)) ;
                  average[k*xsize+l] = (float) (sum[k*xsize+l] / (numOfFramesToAcc- track_invalidPix[k*xsize+l])) ;
                     
                    }
                    sum[k*xsize+l]=0.0;
             }
          }
//        {
//           
//            if((count_track - num_times_repeat)==0){
//                   average[k]=INVALID_PIX_VALUE;
//            }
//            else{
//                   average[k] = (float) (sum[k] / (count_track - num_times_repeat)) ; 
//            }//average is used do identify the average of the pixel value of the  numOfFramesToAcc
//        }
//        for (int k = 0 ; k < xsize * ysize ; k ++)
//        {
//            num_times_repeat = count (track_invalidPix.begin () , track_invalidPix.end () , k) ; //STL method for count
//            if((count_track - num_times_repeat)==0){
//                   average[k]=INVALID_PIX_VALUE;
//            }
//            else{
//                   average[k] = (float) (sum[k] / (count_track - num_times_repeat)) ; 
//            }//average is used do identify the average of the pixel value of the  numOfFramesToAcc
//        }
        ///cout<<"PPPP111111"<<endl;exit(1);
        avgftime = avgftime / count_track ; //calculating average time for accumulated frames
        track_invalidPix.clear () ;
        
        frameno = i + 1 ;
        //setting path for  output frame 
        sprintf (outfile , "%s/%s/%s_t%.4f_f%d_%s_ac.fits" , moduleoutdir , dir , nameprefix , avgftime , frameno , frameIdentifier) ;
        //LOG(INFO)<<endl<<"Outfile "<<i<<"  :  "<<outfile<<"    Status:"<<status;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error creating the image for the output File" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , average , &status) ;
        printError (status , "Error in writing the pixels to output Accumulated file" , outfile) ;      
        writeUsrkeywordsFrmvect (outfile,key_records);
       
        //copying level-1 keywords to output file
        updateKeywords (outfile , 1 , 2 , TINT , "FRAMENO" , &frameno , TDOUBLE , "FRMTIME" , &avgftime) ; //updating the keywords

        FrameTimeMap.insert (make_pair ((unsigned short) frameno , avgftime)) ;

        //write history to the outfile.
        if (history == YES) writeHistory (outfile , vhistorystr) ;

        //write creator,origin,checksum and date to output frame.
        writeCommonKeywords (fout , modulename) ;

        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output Accumulated file" , outfile) ;

        strcpy (outlist[i] , basename (outfile)) ; //taking only filename
        cout<<"Accumulated frame "<<frameno<<" created "<<"\r";
//        cout<<"i "<<i<<"j "<<j<<endl;
    }// end of i loop for number of output frames
    
   //  }
   // while(flag_Acc==FALSE);
    
    int colnum = 0 ;
    if (strcmp (frameIdentifier , "sig") == 0) colnum = 1 ;
    else colnum = 2 ;
    
    // open  & update information file with framelist
    status = writeColumnsToFITS (infofile_out , 2 , 1 , TSTRING , colnum , (void*) outlist , frameno) ;
    if (status)
    {
        LOG (INFO) << "Error writing the columns to  the fits  file" << endl ;
        return (EXIT_FAILURE) ;
    }
  
    // writing level-1 keywords  and  history from vector to output information file
   
    if (colnum == 1)
    {
        writeUsrkeywordsFrmvect (infofile_out , key_records) ;
        if (history == YES)
            writeHistory (infofile_out , vhistorystr) ;
    }
 
    freeMemory (outlist ,frameno , NAMESIZE) ;
    LOG (INFO) << endl << "\nAccumulation complete.." ;
    return (EXIT_SUCCESS) ;
}



