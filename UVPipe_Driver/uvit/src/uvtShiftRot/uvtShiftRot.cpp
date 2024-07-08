/* 
 * File:   uvtShiftRot.cpp
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
#include <uvtShiftRot.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<vector>
#include<glog/logging.h>
#include<transform.h>

//Constructor -called when object is created

uvtShiftRot::uvtShiftRot ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}

//Destructor
//int uvtCosmicRayCorr;

uvtShiftRot::~uvtShiftRot () {
    //delete[] flatfielddata;
}

//parameter File reading 

int uvtShiftRot::read (int argc , char** argv)
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

    if (PIL_OK != (status = PILGetFname ("RASfile" , rasfile)))
    {
        LOG(ERROR) << endl << "***Error reading output directory name***" ;
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

int uvtShiftRot::read (char* inputdatadir , char* rasFile , char* outdir , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->rasfile , rasFile) ;
    strcpy (this->outdir , outdir) ;

    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}
//Parameter file content Display

void uvtShiftRot::display ()
{
    LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT SHIFT ROTATE  PARAMETERS      " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input Frame List Directory                        : " << inputdatadir ;
    LOG(INFO) << endl << "Output Directory                               : " << outdir ;
     LOG(INFO) << endl << "RAS filename                               : " << rasfile;
    if (clobber == YES)
        LOG(INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG(INFO) << endl << "Overwrite                                         : NO" ;
    if (history == YES)
        LOG(INFO) << endl << "History                                             : YES" ;
    else
        LOG(INFO) << endl << "History                                              : NO" ;
    LOG(INFO) << endl << "---------------------------------------------------------------------" ;
}

//Correction for the  Cosmic Ray process

int uvtShiftRot::uvtShiftRotProcess ()
{
    /*if the Input Directory of the Filelist is not available then exit from module
     else read Filelist 
     */
    LOG(INFO) << "\nInside the Shift and Rotate method" << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG(INFO) << endl << "\nModule Output Directory : " << moduleoutdir << endl ;
    string cmd ;
     //check  existence  of output directory if output  exist and 
    //clobber =YES than remove that directory and recreate it,if exist and clobber =NO than exit from module.
    //if directory is not exist than create it .
    if (DirExists (moduleoutdir) && clobber == YES)
    {
        LOG(INFO) << endl << "Directory exists and clobber=yes" ;
        cmd = (string) "rm -rf " + (string) moduleoutdir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists (moduleoutdir) && clobber == NO)
    {
        LOG(INFO) << endl << moduleoutdir << "  already exists " ;
        LOG(INFO) << endl << "\nUse clobber=yes for overwriting" ;
        return (EXIT_FAILURE) ;
    }
    cmd = "mkdir -p " + (string) moduleoutdir ;
    system (cmd.c_str ()) ;
    LOG(INFO) << "\nInput Directory is " << inputdatadir << endl ;
    
    //search for the information file.
    string tempfilepath = searchFile (inputdatadir , ".info") ;
     if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    //LOG(INFO)<<"the fitname "<<Fitsname;

    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists (infofile_in)))
    {
        LOG(ERROR) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY" ;
        return (EXIT_FAILURE) ;
    }
    LOG(INFO) << endl << "\nInformation file :" << infofile_in ;
   
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "***Error in opening the input information file***") ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "***Error in moving the 2nd HDU of the  input information file***") ;
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "***Error in reading the key value of the NAMEPRFX*** ") ; //for creating name for output information file

  //setting path for output information file
    sprintf (infofile_out , "%s/%s_sh.info" , moduleoutdir , nameprefix) ;
    LOG(INFO) << "\nOutput information file :" << infofile_out << endl ;
   
     //creating  output information file
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "***Error in creating the output information file*** ") ;
    char *ttype[] = {"SignalFrames" , "ExposureFrames"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "***Error in creating the table in out info file***") ;
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "***Error in updating the key value of the NAMEPRFX***") ;
    datainfo.write (finfo_out) ; //writing basic information 
    fits_close_file (finfo_out , &status) ;
    printError (status , "***Error in closing the out information file***") ;
  
    //fits pointer for opening the Releative aspect file from user directory
    fitsfile *fras ;
    /**Reading Reletive Aspect File from RA series**/
    LOG(INFO) << "\nReading RAS file from RA series..." << endl ;
    fits_open_file (&fras , rasfile , READONLY , &status) ;
    printError (status , "Error in opening the rasfile",rasfile) ;
    fits_movabs_hdu (fras , 2 , NULL , &status) ;
    printError (status , "Error in moving the 2nd HDU of  ras file",rasfile) ;
    fits_get_num_rows (fras , &no_of_records , &status) ;
    printError (status , "Error Reading the number of Rows",rasfile) ;
    

    time = new double[no_of_records] ;
    roll_ras = new double [no_of_records] ;
    pitch_ras = new double [no_of_records] ;
    yaw_ras = new double[no_of_records] ;

    

    delta_x = new double [no_of_records] ;
    delta_y = new double [no_of_records] ;
    delta_theta = new double[no_of_records] ;
    fits_read_col (fras , TDOUBLE , 1 , 1 , 1 , no_of_records , NULL , time , NULL , &status) ;
    printError (status , "***Error reading  centroid x***") ;
    fits_read_col (fras , TDOUBLE , 2 , 1 , 1 , no_of_records , NULL , yaw_ras , NULL , &status) ;
    printError (status , "***Error writing centroid y***") ;
    fits_read_col (fras , TDOUBLE , 3 , 1 , 1 , no_of_records , NULL , pitch_ras , NULL , &status) ;
    printError (status , "***Error writing x-correction***") ;
    fits_read_col (fras , TDOUBLE , 4 , 1 , 1 , no_of_records , NULL , roll_ras , NULL , &status) ;
    printError (status , "***Error writing y-correction***") ;
    fits_close_file (fras , &status) ;
    printError (status , "Error in closing the file",rasfile) ;
    //logic for integraton will be here.
    LOG(INFO) << "\n Integrating Roll ,pitch,yaw to rotation,xshift,yshift....." << endl ;
    
    double temp_pitch,temp_yaw;
     if (strcasecmp (datainfo.getDetector () , "NUV") == 0)
        {
            for(int i=0;i<no_of_records;i++){
                temp_pitch=pitch_ras[i]*3600;
                temp_yaw=yaw_ras[i]*3600;
                transformRPYtoDXDYDTHETA_NUV (roll_ras[i],temp_pitch,temp_yaw,delta_x[i],delta_y[i],delta_theta[i]);
                cout<<" p"<<delta_x[i]<<endl;
            }
            
        }
        else if(strcasecmp (datainfo.getDetector () , "FUV") == 0)
        {
            for(int i=0;i<no_of_records;i++)
            {
                transformRPYtoDXDYDTHETA_NUV (roll_ras[i],temp_pitch,temp_yaw,delta_x[i],delta_y[i],delta_theta[i]);
            }
        }
        else if(strcasecmp (datainfo.getDetector () , "VIS") == 0)
        {
            for(int i=0;i<no_of_records;i++)
            {
              transformRPYtoDXDYDTHETA_NUV (roll_ras[i],temp_pitch,temp_yaw,delta_x[i],delta_y[i],delta_theta[i]);
            }
             
        }
        else
        {
            LOG(INFO)<<"***Invalid Channel***"<<endl;
            return(EXIT_FAILURE);
        }
    //cout<<"Final value is  "<<delta_y[no_of_records-1]<<endl;exit(1);
    
//    for (int i = 0 ; i < no_of_records ; i++)
//    {
//        delta_x[i] = delta_Xras[i] ;
//        delta_y[i] = delta_Yras[i] ;
//        
//        delta_theta[i] = theta_ras[i] ; //*M_PI/180
//     
//       
//    }
  
    //in case of IM mode 
    if (datainfo.getModeFlag () == IM)
    {
        fits_read_key (finfo_in , TSTRING , "SIGDIR" , sigframedir , NULL , &status) ;
        printError (status , "***Error in finding the key value of SIGDIR***") ;
        fits_read_key (finfo_in , TSTRING , "EXPDIR" , expframedir , NULL , &status) ;
        printError (status , "***Error in  reading the key value of the EXPDIR ***") ;
        fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "***Error in reading the key value of the NFILES***") ;
        sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        expoframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) sigframelist , NULL , &status) ;
        printError (status , "***Error in reading the column  for input Signal frame list***") ;
        fits_read_col (finfo_in , TSTRING , 2 , 1 , 1 , nframes , NULL , (void *) expoframelist , NULL , &status) ;
        printError (status , "***Error in reading the column  for input Exposure  frame list***") ;

        //method for the cosmic ray correction  
        if (uvtshiftrotateIM ()) return (EXIT_FAILURE) ;
        fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
        printError (status , "***Error in opening the out information file***") ;
        fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
        printError (status , "***Error in moving  the out information file***") ;
        datainfo.setXsize (xsize) ;
        datainfo.setYsize (ysize) ;
        datainfo.write (finfo_out) ;
        fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
        printError (status , "***Error in key value  of the NAMEPRFX ***") ; //for creating name for output information file
        fits_update_key (finfo_out , TSTRING , "SIGDIR" , sigframedir , NULL , &status) ;
        printError (status , "***Error in updating the  key value of the SIGDIR ***") ;
        fits_update_key (finfo_out , TSTRING , "EXPDIR" , expframedir , NULL , &status) ;
        printError (status , "***Error in updating the key value of the EXPDIR***") ;
        fits_update_key (finfo_out , TINT , "NFILES" , &nframes , NULL , &status) ;
        printError (status , "***Error in updating the NFILES keyword  not updated/not Found") ;
        fits_close_file (finfo_out , &status) ;
        fits_close_file (finfo_in , &status) ;
        printError (status , "***Error in closing the ot info file***") ;
        freeMemory (sigframelist , nframes , NAMESIZE) ;
        freeMemory (expoframelist , nframes , NAMESIZE) ;
    } //in case of the PC 
    else if (datainfo.getModeFlag () == PC)
    {
        fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
        printError (status , "***Error in reading the key value of the EVTFILES ***") ;
        if (uvtshiftrotatePC ()) return (EXIT_FAILURE) ;
        fits_close_file (finfo_in , &status) ;
    } //else neither PC or IM(i.e invalid mode )
    else
    {
        LOG(ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG(ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
    return (EXIT_SUCCESS) ;

}

int uvtShiftRot::uvtshiftrotateIM ()
{
    LOG(INFO) << endl << "\nStarted Shift and Rotate process..." << endl ;
    //creating signal directory
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s%s" , moduleoutdir , sigframedir) ;
    /**Shell command for creating the directory**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    //creating exposure directory
    sprintf (dir , "%s%s" , moduleoutdir , expframedir) ;
    cmd = "mkdir -p " + (string) dir ;
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    LOG(INFO) << endl << "\nTotal number of frames - " << nframes ;
    int status = 0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    int bitpix = FLOAT_IMG ;
    vector<string> vhistorystr ;
    if (history == YES)
        getHistory (vhistorystr) ;
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    vector<string> sigfilename , expfilename ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;
    int naxis = 2 ;

   
    vector<unsigned short> frame_no ;
    vector<double> frame_time ;
   
    fitsfile *fptr,*fptr1 ;
    //loop for getting the information of frame times related to frame no.
    for (int i = 0 ; i < nframes ; i++)
    {
        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[i]) ;
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input File") ;
        fits_read_key (fptr , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ;
        printError (status , "***Error in reading the key value of the FRAMENO***") ;
        fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
        printError (status , "****Error reading the key value of the FRMTIME") ;
        frame_no.push_back (frameno) ;
        frame_time.push_back (frametime) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the File",infile) ;
    }
    
    double t1 , t2 , theta1 , theta2 , x1 , x2 , y1 , y2 ;
    double *new_delta_theta , *new_delta_x , *new_delta_y ;
    new_delta_theta = new double[nframes] ;
    new_delta_x = new double[nframes] ;
    new_delta_y = new double[nframes] ;
    for (int p = 0 ; p < nframes ; p++)
    {
        new_delta_x[p] = 0.0f ;
        new_delta_y[p] = 0.0f ;
        new_delta_theta[p] = 0.0f ;
    }
   // int point =0;
    
    //processed done for the  number of files in the Input Directory
    //getting median deltas values for all these parameters depending on the time matching
  
    float mult_fact=xsize/600;
    for (int index = 0 ; index < nframes ; index++)// 
    {
        for (int index2 = 0 ; index2 < no_of_records ; index2++)
        {
            if (frame_time[index] >= time[index2] && frame_time[index] < time[index2 + 1])
            {
           //     point++;
               
                t1 = time[index2] ;
                t2 = time[index2 + 1] ;
                theta1 = delta_theta[index2] ;
                theta2 = delta_theta[index2 + 1] ;
                x1 = delta_x[index2]*mult_fact ;
                x2 = delta_x[index2 + 1]*mult_fact ;
                y1 = delta_y[index2] *mult_fact;
                y2 = delta_y[index2 + 1]*mult_fact ;
                new_delta_theta[index] = theta1 + (frame_time[index] - t1)*(theta2 - theta1) / (t2 - t1) ;
                new_delta_x[index] = x1 + (frame_time[index] - t1)*(x2 - x1) / (t2 - t1) ;
                new_delta_y[index] = y1 + (frame_time[index] - t1)*(y2 - y1) / (t2 - t1) ;
               
                break;
           }
        }
    }
  
    
    float **subSignalArray , **subExposureArray , **tempSigarr , **tempExparr ;
    subSignalArray = new float * [xsize] ;
    subExposureArray = new float * [xsize] ;
    tempSigarr = new float * [xsize] ;
    tempExparr = new float * [xsize] ;
    for (int i = 0 ; i < xsize ; i++)
    {
        subExposureArray[i] = new float [xsize] ;
        subSignalArray[i] = new float [xsize] ;
        tempExparr[i] = new float [xsize] ;
        tempSigarr[i] = new float [xsize] ;
    }
    float *one_dim_img ;
    one_dim_img = new float[xsize * ysize] ;
    float *one_dim_exp ;
    one_dim_exp = new float[xsize * ysize] ;
    long naxes[2] = {xsize , ysize} ;
    int SIZE = xsize / 2 ;
    long temp1 = 0 ;
    fitsfile *outfptr1 ;
    LOG(INFO) << "\nApplying Shift and Rotation..." << endl ;
    for (int index = 0 ; index < nframes ; index++)
    {
        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[index]) ;
        for (int p = 0 ; p < xsize * ysize ; p++)
        {
            one_dim_exp[p] = 0 ;
            one_dim_img[p] = 0 ;
        }
       
        /***Reading the signal frame array***/
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input Frame") ;
        fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
        printError (status , "Error reading the key value of the FRMTIME",infile) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , one_dim_img , NULL , &status) ;
        printError (status , "Error in reading the pixel value of the input file" , infile) ;
         if (index==0)
        {
            copyUsrkeywrdsTovect (fptr,key_records);
        }
        temp1 = 0 ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , expframedir , expoframelist[index]) ;
        /***Reading the Exposure frame array***/
        fits_open_file (&fptr1 , infile , READONLY , &status) ;
        printError (status , "Error in opening the input File",infile) ;
        fits_read_pix (fptr1 , TFLOAT , fpixel , xsize*ysize , NULL , one_dim_exp , NULL , &status) ;
        printError (status , "Error in reading the pixel value of the input file" , infile) ;

        for (int i = 0 ; i < xsize ; i++)
        {
            for (int j = 0 ; j < ysize ; j++)
            {
                tempSigarr[j][i] = -9999 ;
                tempExparr[j][i] = -9999 ;
            }
        }
        for (int i = 0 ; i < xsize ; i++)
        {
            for (int j = 0 ; j < ysize ; j++)
            {
                subExposureArray[j][i] = one_dim_exp[temp1] ;
                subSignalArray[j][i] = one_dim_img[temp1] ;
                temp1++ ;
            }
        }
        double ctheta , stheta ;
        //applying the deltas 
        
        ctheta = cos ((double)(new_delta_theta[index] * M_PI / 180)) ;
        stheta = sin ((double)(new_delta_theta[index] * M_PI / 180) );
        

        double index_i = 0.0 , index_j = 0.0 ;
        int i1 , j1;
        for (int i = 0 ; i < xsize ; i++)
        {
            i1 = i - SIZE ;
            for (int j = 0 ; j < ysize ; j++)
            {
                
               // if(subSignalArray[i][j]!=INVALID_PIX_VALUE){
                      j1 = j - SIZE ;
                /*applying new_delta_theta[index] degree rotation to find out new pixel indexes */
                index_i =   ((i1 * ctheta) - (j1 * stheta)) + SIZE-new_delta_x[index] ; //new index x
                index_j =   ((i1 * stheta) + (j1 * ctheta)) + SIZE-new_delta_y[index] ; //new index y
                round (index_i) ;
                round (index_j) ;
                if (index_i > -1 && index_i < xsize && index_j > -1 && index_j < ysize)
                {
                    /*appling correction i.e xshift,yshift ,rotation.*/
                   
                    tempSigarr[(int) index_i][(int) index_j] = subSignalArray[i][j] ;
                    tempExparr[(int) index_i][(int) index_j] = subExposureArray[i][j] ;
               
                }
             }
              
            //}
        }
        temp1 = 0 ;
        for (int i = 0 ; i < ysize ; i++)
        {
            for (int j = 0 ; j < xsize ; j++)
            {
                one_dim_img[temp1] = tempSigarr[j][i] ;
                one_dim_exp[temp1] = tempExparr[j][i] ;
                temp1++ ;
            }
        }
       
        /*generating the output Signal Directory*/
        sprintf (outfile , "%s/%s/%s_t%f_f%d_sig_sh.fits" , moduleoutdir , sigframedir , nameprefix , frametime , frame_no[index]) ;
        fits_create_file (&outfptr1 , outfile , &status) ;
        printError (status , "Error in creating the output Signal Frame ") ;
        fits_create_img (outfptr1 , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the image") ;
        fits_write_pix (outfptr1 , TFLOAT , fpixel , xsize*ysize , one_dim_img , &status) ;
        printError (status , "Error in writing the pixels to the output Signal frame") ;
     
        //copy level-1 keywords to output file
        copyUserKeywords (fptr , outfptr1) ;
        if(history==YES)
         {
               writeHistory(outfile , vhistorystr) ; //write history to each file
         }
        writeCommonKeywords (outfptr1 , modulename) ;
       
        fits_close_file (outfptr1 , &status) ;
        printError (status , "Error in closing the input Signal Frame") ;    
                
        sigfilename.push_back (basename (outfile)) ;
        /*generating the output Exposure  Directory*/
        sprintf (outfile , "%s/%s/%s_t%f_f%d_exp_sh.fits" , moduleoutdir , expframedir , nameprefix , frametime , frame_no[index]) ;
        fits_create_file (&outfptr1 , outfile , &status);
        printError (status , "Error in creating the file",outfile) ;
        fits_create_img (outfptr1 , FLOAT_IMG , naxis , naxes , &status);
        printError (status , "Error in creating the img ",outfile) ;
        fits_write_pix (outfptr1 , TFLOAT , fpixel , xsize*ysize , one_dim_exp , &status) ;
        printError (status , "Error in writing the pixels to Exposure Frame") ;
        //copy level-1 keywords to output file
        copyUserKeywords (fptr1 , outfptr1) ;
        expfilename.push_back (basename (outfile)) ;
        if(history==YES)
        {
               writeHistory(outfile , vhistorystr) ; //write history to each file
        }
         writeCommonKeywords (outfptr1 , modulename) ;
         fits_close_file (outfptr1 , &status);
         printError (status , "Error in closing the input Exposure File") ;
         fits_close_file (fptr1 , &status) ;
         printError (status , "Error in Closing the File",infile) ;
        fits_close_file (fptr, &status) ;
        printError (status , "Error in Closing the File",infile) ;
        
    }
    //updating the information file
    LOG(INFO)<<"\nWriting to output  information file... "<<endl;
    fitsfile *finfo_out ;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "***Error in opening  the output information file***") ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "***Error in  moving to the 2nd HDU of the out information file***") ;
    fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , nframes , (void *) sigfilename.data () , &status) ;
    printError (status , "***Error in writing the column of outout signal frame list***") ;
    fits_write_col (finfo_out , TSTRING , 2 , 1 , 1 , nframes , (void *) expfilename.data () , &status) ;
    printError (status , "***Error in writing the column of outout Exposure frame list***") ;
      if(history==YES)
      {
               writeHistory(infofile_out , vhistorystr) ; //write history to each file
       }
    writeUsrkeywordsFrmvect (infofile_out,key_records);//writing user keywords to the information file
    writeCommonKeywords (finfo_out, modulename) ;//writing common keywords 
    fits_close_file (finfo_out , &status) ;
    printError (status , "***Error in closing the out information file***") ;
    
    sigfilename.clear () ;
    expfilename.clear () ;
    for (int i = 0 ; i < xsize ; i++)
    {
        free (subSignalArray[i]) ;
        free (subExposureArray[i]) ;
        free (tempSigarr[i]) ;
        free (tempExparr[i]) ;
    }
    //Free The memory
    free (subSignalArray) ;
    free (subExposureArray) ;
    free (one_dim_img) ;
    free (one_dim_exp) ;
    free (tempSigarr) ;
    free (tempExparr) ;
    return (EXIT_SUCCESS) ;
}

int uvtShiftRot::uvtshiftrotatePC ()
{
    LOG(INFO) << "\nStarted Shift and Rotation process for PC mode...." << endl ;
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    int status = 0 ;
    vector<string> vhistorystr ;
    if(history==YES)
        getHistory(vhistorystr);
    
    char eventfilename[FLEN_FILENAME] ;
    //setting path for output event file.
    sprintf (eventfilename , "%s_sh.events" , nameprefix) ;
    sprintf (outfile , "%s/%s" , moduleoutdir , eventfilename) ;
    long nrows = 0 ;
    sprintf (infile , "%s/%s" , inputdatadir , eventfile) ;
    LOG(INFO) << "\nInput Event file :" << infile << endl ;
    LOG(INFO) << "\nOutput Event file :" << outfile << endl ;
    //opening the input event file.
    fitsfile *fptr ;
    fits_open_file (&fptr , infile , READONLY , &status) ;
    printError (status , "***Error in opening the input Event File***") ;
  
    fitsfile *fout ;
    fits_create_file (&fout , outfile , &status) ;
    printError (status , "Error in creating the file" , outfile) ;
    fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
    printError (status , "Error in copying the  file" , infile) ;
  
    //copy user keywords to vector
    copyUsrkeywrdsTovect (fptr,key_records);
    double *t ;
    long *frm_no;
    float *xf , *yf , *xi_final , *yi_final ;

    double t1 , t2 , theta1 , theta2 , x1 , x2 , y1 , y2 ;
    double new_delta_theta=0.0 , new_delta_x = 0.0 , new_delta_y = 0.0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
   
    fits_movabs_hdu (fout , 2 , NULL , &status) ;
    printError (status , "Error in  moving to the 2nd HDU of the input  file" , infile) ;
    fits_get_num_rows (fout , &nrows , &status) ;
    printError (status , "Error in  reading the number of rows" , infile) ;

    t = new double[nrows] , xf = new float[nrows] , yf = new float[nrows] ;
    xi_final = new float[nrows] , yi_final = new float[nrows] ;
    frm_no = new long[nrows];
    for (int i = 0 ; i < nrows ; i++)
    {
        xi_final[i] = 0.0f ;
        yi_final[i] = 0.0f ;
    }
    fits_read_col (fout , TLONG,2 , 1 , 1 , nrows , NULL , frm_no , NULL , &status) ;
    printError (status , "Error in  reading the column of t " , infile) ;
    fits_read_col (fout , TDOUBLE , 3 , 1 , 1 , nrows , NULL , t , NULL , &status) ;
    printError (status , "Error in  reading the column of t " , infile) ;
    fits_read_col (fout , TFLOAT , 4 , 1 , 1 , nrows , NULL , xf , NULL , &status) ;
    printError (status , "Error in  reading the column of xi" , infile) ;
    fits_read_col (fout , TFLOAT , 5 , 1 , 1 , nrows , NULL , yf , NULL , &status) ;
    printError (status , "Error in  reading the column of yi" , infile) ;
    
    float i1 = 0 , j1 = 0 ;
   float time_frame = 0.0 ;
    bool flag =FALSE;
  float mult_fact=xsize/600;
  vector<long> track_rown;
    for (int i = 0 ; i < nrows ; i++)
    { 
        flag =FALSE;
        time_frame = t[i] ;
        
        //calculating the median of deltas from Relative aspect file.
        for (int index = 0 ; index < no_of_records ; index++)
        {
            if (time_frame >= time[index] && time_frame < time[index + 1])
            {
                flag=TRUE;
                t1 = time[index] ;
                t2 = time[index + 1] ;
                theta1 = delta_theta[index] ;
                theta2 = delta_theta[index + 1] ;
                x1 = delta_x[index] ;
                x2 = delta_x[index + 1] ;
                y1 = delta_y[index] ;
                y2 = delta_y[index + 1] ;
                new_delta_theta = (theta1 + (time_frame - t1)*(theta2 - theta1) / (t2 - t1));//*mult_fact ;
                new_delta_x = (x1 + (time_frame - t1)*(x2 - x1) / (t2 - t1))*mult_fact ;
                new_delta_y = (y1 + (time_frame - t1)*(y2 - y1) / (t2 - t1))*mult_fact ;
                break ;
            }
        }
        //check whether  matching time found or not
        if(flag==FALSE){
            track_rown.push_back (i+1);
            new_delta_theta=0;
             new_delta_x=0;
             new_delta_y=0;
        }
      
        double ctheta , stheta ;
        ctheta = cos (1.0*new_delta_theta * M_PI / 180) ;
        stheta = sin (1.0*new_delta_theta * M_PI / 180) ;
        i1 = xf[i] - xsize / 2 ;
        j1 = yf[i] - ysize / 2 ;
        //calculating the new values for Xand Y based on the deltas
        if ((((i1 * ctheta) - (j1 * stheta)) + xsize / 2-new_delta_x)>-1 && (xi_final[i] = ((i1 * ctheta) - (j1 * stheta)) + xsize / 2-new_delta_x)< xsize && (yi_final[i] = ((i1 * stheta) + (j1 * ctheta)) + ysize / 2-new_delta_y)>-1 && (yi_final[i] = ((i1 * stheta) + (j1 * ctheta)) + ysize / 2-new_delta_y) < xsize)
        {
            xi_final[i] = ((i1 * ctheta) - (j1 * stheta)) + xsize / 2-new_delta_x ; //new index x
            yi_final[i] = ((i1 * stheta) + (j1 * ctheta)) + ysize / 2-new_delta_y ; //new index y
        }
    }
  //writing the updated X and Y values to the output event file.
    fits_write_col (fout , TDOUBLE , 3 , 1 , 1 , nrows , t , &status) ;
    printError (status , "Error in  writing  the column of t" , outfile) ;
    fits_write_col (fout , TFLOAT , 4 , 1 , 1 , nrows , xi_final , &status) ;
    printError (status , "Error in  writing  the column of X-integer" , outfile) ;
    fits_write_col (fout , TFLOAT , 5 , 1 , 1 , nrows , yi_final , &status) ;
    printError (status , "Error in  writing  the column of  Y-integer" , outfile) ;
     fits_delete_rowlist (fout , track_rown.data () , track_rown.size () , &status) ;
    printError (status , "Error in Deleting the row list " , outfile) ;
    if(history==YES)
     {
               writeHistory(outfile , vhistorystr) ; //write history to each file
     }
    //writeUsrkeywordsFrmvect (outfile,key_records);
    writeCommonKeywords (fout , modulename) ;
    
    
    fits_close_file (fout , &status) ;
    printError (status , "Error in closing the file" , infile) ;
    LOG(INFO) << "\nUpdating output information file......" << endl ;
    fitsfile *finfo_out ;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU in output information file" , infofile_out) ;
    fits_update_key (finfo_out , TSTRING , "EVTFILE" , basename (outfile) , NULL , &status) ;
    printError (status , "Error in updating the key value of the output information file") ;
    
        if(history==YES){
               writeHistory(infofile_out , vhistorystr) ; //write history to each file
       }
    writeUsrkeywordsFrmvect (infofile_out,key_records);//writing user keywords from vector to the output information file.
    writeCommonKeywords (finfo_out , modulename) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file  " , infofile_out) ;
    return (EXIT_SUCCESS) ;
}

int uvtShiftRot::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputDataDir=" + (string)inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" RAS file name  " + (string)rasfile) ;
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

