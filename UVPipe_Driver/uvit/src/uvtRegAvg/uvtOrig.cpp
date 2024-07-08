#include<iostream>
#include<unistd.h>
#include<stdlib.h>
#include<dirent.h>  //Accessing Directory
#include<string.h>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include "uvtRegAvg.h"
#include<pthread.h>
#include<uvtUtils.h>
#include<spMatrix.h>
#include <math.h>
#include<fftw3.h>
#include<algorithm>
#include<glog/logging.h>
#include "macro_def.h"


void  roundoff(float &val){
    val=val-0.5;
    val=ceil(val);
}
//Constructor -called when object is created
bool compare (struct Star star1 , struct Star star2) ;
typedef struct Star
{
    float intensity ;
    int x , y ;
    void print ()
    {
        LOG(INFO) << endl << "(" << x << "," << y << ") : " << intensity ;
    }
}; 
uvtRegAvg::uvtRegAvg ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    strcpy (centroidDir , "Centroid") ;
}
//Destructor
uvtRegAvg::~uvtRegAvg ()
{
    //freeMemory (sigframelist , nframes , NAMESIZE) ;
    //freeMemory (expoframelist , nframes , NAMESIZE) ;
}
//parameter File reading
int uvtRegAvg::read (int argc , char** argv)
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
    if (PIL_OK != (status = PILGetInt ("algo_flag" , &algo_flag)))
    {
        LOG(ERROR) << endl << "***Error reading algo_flag ***" ;
        return status ;
    }
    if (algo_flag == 0)
    {
        if (PIL_OK != (status = PILGetReal4 ("threshold" , &sd_multi_factor_default)))
        {
            LOG(INFO) << endl << "***Error reading output directory name***" ;
            return status ;
        }
        if (PIL_OK != (status = PILGetInt ("centroidlimit" , &centroidlimit)))
        {
            LOG(INFO) << endl << "***Error reading output directory name***" ;
            return status ;
        }
        if (PIL_OK != (status = PILGetInt ("minimum_targetedstars" , &minimum_No_of_Stars)))
        {
            LOG(INFO) << endl << "***Error reading algo_flag ***" ;
            return status ;
        }

        if (PIL_OK != (status = PILGetInt ("refine_Window" , &refine_Winsize)))
        {
            LOG(INFO) << endl << "***Error reading refine window size ***" ;
            return status ;
        }
        if (PIL_OK != (status = PILGetInt ("centroid_Window" , &centroid_Winsize)))
        {
            LOG(INFO) << endl << "***Error reading refine window size ***" ;
            return status ;
        }
    }
    else if (algo_flag == 1)
    {

        if (PIL_OK != (status = PILGetInt ("StarDetectionCentroid_square_size" , &algo_Square_Size)))
        {
            LOG(INFO) << endl << "Error reading refine window size" ;
            return status ;
        }
        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroid_threshold" , &primary_threshold_Val)))
        {
            LOG(INFO) << endl << "Error reading refine window size" ;
            return status ;
        }
        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroid_secondary_threshold" , &secondary_threshold_Val)))
        {
            LOG(INFO) << endl << "Error reading refine window size" ;
            return status ;
        }
    }
    else
    {
        LOG(ERROR) << "algo_flag must be 1 or 2" << endl ;
        return (EXIT_FAILURE) ;
    }

    if (PIL_OK != (status = PILGetInt ("option_LeastSquare" , &option_LeastSquare)))
    {
        LOG(ERROR) << endl << "Error reading the option_Leastsquare" ;
        return status ;
    }
      if (PIL_OK != (status = PILGetReal4  ("diff_Dist" , (float*) &diff_dist)))
   {
          LOG(ERROR) << endl << "***Error reading Distance value for matching***" ;
          return status ;
   }
    
   PILClose (status) ;
    return (EXIT_SUCCESS) ;
}

int uvtRegAvg::read (char *inputdatadir  , char *outdir , int clobber , int history ,int algoflag, float sdMutltifactor ,int centlimit,int min_nostars,int refinewinsize,int centroidwinsize,int algosqrsize,float prithr,float secthr,float opLeastsqr,float diffDist)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;
    this->algo_flag=algoflag;
     this->sd_multi_factor_default=sdMutltifactor;
      this->centroidlimit=centlimit;
       this->minimum_No_of_Stars=min_nostars;
           this->refine_Winsize=refinewinsize;
       this->centroid_Winsize=centroidwinsize;
       this->algo_Square_Size=algosqrsize;
       this->primary_threshold_Val=prithr;
       this->secondary_threshold_Val=secthr;
       this->option_LeastSquare=opLeastsqr;
       this->diff_dist=diffDist;
     //strcpy (this->caldbDir , caldbDir) ;
    //this->numFrames2Gen = no_frames_gen ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}
//Parameter file content Display

void uvtRegAvg::display ()
{
    LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT REG AVG  PARAMETERS      " << endl ;
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
    LOG(INFO) << endl << "----------------------------------------------------------------" ;
}
//Correction for the  Cosmic Ray process

int uvtRegAvg::uvtRegAvgProcess ()
{
    LOG(INFO) << endl << "Registration and Average process started" << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    //LOG(INFO)<<endl<<"Module Output Directory : "<<moduleoutdir<<endl;
    string cmd ;
  //  diff_dist=4096;
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
    /**Shell commad for creating the output directory**/
    cmd = "mkdir -p " + (string) moduleoutdir ;
    /**Executing  the Shell command**/
    system (cmd.c_str ()) ; // creating output directory to keep output from unitConversion
    LOG(INFO) << endl << moduleoutdir << "  directory created" ;
    const char * tempfilepath = searchFile (inputdatadir, ".info") ;
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath) ;
    LOG(INFO) << endl <<" \nInformation File : " << infofile_in ;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists (infofile_in)))
    {
        LOG(ERROR) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY" ;
        return (EXIT_FAILURE) ;
    }
    LOG(INFO) << endl <<" \nInformation File : " << infofile_in ;
    int status = 0 ;
    fitsfile *finfo_in  ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the information file",infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in Moving the 2nd HDU",infofile_in) ;
    datainfo.getInfo(finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "NAMEPRFX keyword not Found",infofile_in) ; //for creating name for output information file
    fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
    printError (status , "NFILES  keyword not Found",infofile_in) ; //for creating name for output information file

    sprintf (infofile_out , "%s/%s_rav.info" , moduleoutdir , nameprefix) ;
//    fits_create_file (&finfo_out , infofile_out , &status) ;
//    printError (status , "Error in creating the output information file",infofile_in) ;
//    char *ttype[] = {"SignalFrames" , "ExposureFrames"} ;
//    char *tform[] = {"A256" , "A256"} ;
//    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
//    printError (status , "Error in creating the table",infofile_out) ;
//    datainfo.write (finfo_out) ; //writing basic data information
//    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
//    printError (status , "Error in updating the key value of NAMEPRFX ",infofile_out) ; //for creating name for output information file
//    /*----info file creating completed, rest of the information will be put by other functions-----------*/
//    fits_close_file (finfo_out , &status) ;
//    printError (status , "Error in closing the output information file ",infofile_out) ;
     sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
    expoframelist = allocateMemory<char>(nframes , NAMESIZE) ;
    fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) sigframelist , NULL , &status) ;
    printError (status , "Error in reading the reading the signal frame list",infofile_in) ;
    fits_read_col (finfo_in , TSTRING , 2 , 1 , 1 , nframes , NULL , (void *) expoframelist , NULL , &status) ;
    printError (status , "Error in reading the reading the Exposure frame list ",infofile_in) ;
    sprintf (sigframedir , "%s" , "SignalFrames") ;
   sprintf (expframedir , "%s" , "ExposureFrames") ;
    if (datainfo.getModeFlag () == IM || datainfo.getModeFlag ()==PC)
    {//For IM mode
        //Performing RegAvg
        if (reg_Avg ()) return (EXIT_FAILURE) ;
   }//in case of the PC and IM
    else//else neither PC or IM(i.e invalid mode )
    {
        LOG(ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG(ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
    /**Updating the output information file**/
//    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
//    printError (status , "Error in opening the output infofile",infofile_out) ;
//    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
//    printError (status , "Error in moving the 2nd HDU in out information file",infofile_out) ;
//    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
//    printError (status , "Error in updating the  key value of the NAMEPRFX",infofile_out) ; //for creating name for output information file
//    fits_update_key (finfo_out , TSTRING , "SIGDIR" , sigframedir , NULL , &status) ;
//    printError (status , "Error in updating the key value of the SIGDIR",infofile_out) ;
//    fits_update_key (finfo_out , TSTRING , "EXPDIR" , expframedir , NULL , &status) ;
//    printError (status , "Error updating the key value of the EXPDIR",infofile_out) ;
//    fits_close_file (finfo_out , &status) ;
//    printError (status , "Error in closing the out informaton  file",infofile_out) ;
   
    return (EXIT_SUCCESS) ;
}

int uvtRegAvg::reg_Avg()
{
    LOG(INFO) << endl << "\nStarted Registration and Averaging of frames..."<<endl;
    long total_size=xsize*ysize ;
    int status = 0 ;
    fitsfile *infptr2,*infptr1 ;
    fitsfile *outfptr2 ;  
     vector<string> vhistorystr ;
    if (history==YES){
        getHistory (vhistorystr) ;
    }
    int naxis = 2  ;
  long naxes[2] = {FINALFRAMESIZE_REGAVG, 600} ;
 //    long naxes[2] = {xsize , ysize} ;
    long nelements = total_size ;
    char infile[NAMESIZE] , outfile[NAMESIZE] ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    LOG(INFO) << "\nTotal number of frames - "<< nframes << endl ;
    // int addition =nframes / no_of_Frames_toGen;
    sd_mul_factor=sd_multi_factor_default;
    float avg_sigArray[xsize*ysize];
    float avg_expArray[xsize*ysize];
    for(int i=0;i<nelements;i++){
        avg_sigArray[i]=0.0f;
        avg_expArray[i]=0.0f;            
    }
   if (numFrames2Gen <= 1)
    {
        numFrames2Gen = nframes ;
    }
    LOG(INFO)<<"\nPerforming Registration and Averaging..."<<endl;
    float *sig_array_ref = new float[nelements] ;
    sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[0]) ;
    fits_open_file (&infptr2 , infile, READONLY , &status) ;
    printError (status , "Error in opening the input File" , infile) ;
    fits_read_pix (infptr2 , TFLOAT , fpixel , nelements , NULL , sig_array_ref , NULL , &status) ;
    printError (status , "Error in reading the pixels  values of the input File") ;
    copyUsrkeywrdsTovect (infptr2,key_records);    
    fits_close_file (infptr2 , &status) ;
    printError (status , "Error in closing the  input File" , infile) ;
    vector<float> cx_ref , cy_ref , ci_ref ;
    if (algo_flag == 0)
    {
        status = findStar_algo1 (sig_array_ref) ;              // for first frame taken as  reference
        if (status)
        {
            LOG(ERROR) << endl << "***Error in finding star algorithm 1 for frame  " << infile << "  ***" << endl ;
            return (EXIT_FAILURE) ;
        }

        for (int i = 0 ; i < Cx.size () ; i++)
        {
            cx_ref.push_back (Cx[i]) ;
            cy_ref.push_back (Cy[i]) ;
            ci_ref.push_back (Ci[i]) ;
        }
        Cx.clear () ;
        Cy.clear () ;
        Ci.clear () ;
  }

    float *sig_array = new float[xsize*ysize] ;
    float *exp_array = new float[xsize*ysize] ;
    float *sig_array_temp = new float[xsize*ysize] ;
    float *exp_array_temp = new float[xsize*ysize] ;
    vector<float> corr_vect_first_x , corr_vect_second_x , corr_vect_first_y , corr_vect_second_y ;

    double x_dx = 0.0 , y_dy = 0.0 , theta_dt = 0.0 ;
    int win_range = TEMPLATE_SIZE ;
    if (win_range % 2 == 0)
    {
        win_range = win_range + 1 ;
    }
    float *arr_final = new float[win_range * win_range] ;
    vector<Star> framevector ;
    //Star array[10];
    Star starobj ;
    int min_size=0;
    vector<float> x_corr_mf , y_corr_mf ;
    vector<string> name_sigfile , name_expfile ;
    char temp_path_txt[FLEN_FILENAME];
    sprintf (temp_path_txt , "%s/shifts.txt" , moduleoutdir) ; //file For storing the After  FFT content.
    //ofstream ofptr;
      //  ofptr.open (temp_path_txt , ios::out) ;
    for (int i = 1 ; i < nframes ; i++)
    {
       LOG(INFO)<<i<<endl;
        for (int i = 0 ; i < nelements ; i++)
        {
            sig_array_temp[i] = 0.0f ;
            exp_array_temp[i] = 0.0f ;
            sig_array[i] = 0.0f ;
            exp_array[i] = 0.0f ;
        }
        corr_vect_first_x.clear () ;
        corr_vect_first_y.clear () ;
        corr_vect_second_x.clear () ;
        corr_vect_second_y.clear () ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[i]) ;
        fits_open_file (&infptr1 , infile , READONLY , &status) ;
        printError (status , "Error in opening the input  Signal File" , infile) ;
        fits_read_pix (infptr1 , TFLOAT , fpixel , nelements , NULL , sig_array , NULL , &status) ;
        printError (status , "Error in reading the pixels  values of the input File") ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , expframedir , expoframelist[i]) ;
        fits_open_file (&infptr2 , infile , READONLY , &status) ;
        printError (status , "Error in opening the input File" , infile) ;
        fits_read_pix (infptr2 , TFLOAT , fpixel , nelements , NULL , exp_array , NULL , &status) ;
        printError (status , "Error in reading the pixels  values of the input File") ;
        fits_close_file (infptr1 , &status) ;
        printError (status , "Error in closing the  input File" , infile) ;
        fits_close_file (infptr2 , &status) ;
        printError (status , "Error in closing the  input File" , infile) ;
        int cnt_input = 0 ;
        for (int i = 0 ; i < nelements ; i++)
        {
            if (sig_array[i] != 0.0)
            {
                cnt_input++ ;
            }
        }
        x_corr = 0.0 , y_corr = 0.0 ;
        if (algo_flag == 1)
        {
            framevector.clear () ;
            for (int i = 0 ; i < nelements ; i++)
            {
                if (sig_array[i] >= 1)
                {
                    starobj.intensity = sig_array[i] ;
                    starobj.x = (int) round (i % xsize) ;
                    starobj.y = (int) round (i / ysize) ;
                    framevector.push_back (starobj) ;
             }
            }
            LOG(INFO) << endl << "Insertion into vector complete" ;
            LOG(INFO) << endl << "Vector size: " << framevector.size () << endl ;
            sort (framevector.begin () , framevector.end () , compare) ;
            LOG(INFO) << endl << "Selected Peaks for template window" ;
            for (int i = 0 ; i < 5 ; i++)
                framevector[i].print () ;
            LOG(INFO) << endl ;
            LOG(INFO) << "Sorting completed" << endl ;
            int x , y ;
            double val = 0.0 ;
            for (int i = 0 ; i < NO_OF_ELEMENTS_TO_CMPR ; i++)
            {
                x = framevector[i].x ;
                y = framevector[i].y ;
                LOG(INFO) << x << " " << y << endl ;
                int cnt = 0 ;
                for (int j = -1 * (win_range / 2) ; j <= (win_range / 2) ; j++)
                {
                    for (int k = -1 * (win_range / 2) ; k <= (win_range / 2) ; k++)
                    {

                        if ((y + win_range / 2) < xsize && (x + win_range / 2) < ysize)
                        {

                            val = sig_array[(y + j) * ysize + (x + k)] ;
                        }
                        else
                        {

                            if (y + win_range / 2 >= ysize && x + win_range / 2 < xsize)
                            {

                                while (y + win_range / 2 >= ysize)
                                {
                                    y = y - 1 ;
                                }

                            }
                            else if (x + win_range / 2 >= xsize && y + win_range / 2 < ysize)
                            {
                                while (x + win_range / 2 >= xsize)
                                {
                                    x = x - 1 ;
                                }

                            }
                            else if (x + win_range / 2 >= xsize && y + win_range / 2 >= ysize)
                            {
                                while (x + win_range / 2 >= xsize)
                                {
                                    x = x - 1 ;
                                }
                                while (y + win_range / 2 >= ysize)
                                {
                                    y = y - 1 ;
                                }
                            }
                            else
                            {
                                LOG(INFO) << "ERROR***" << endl ;
                                return (EXIT_FAILURE) ;
                            }                          
                            val = sig_array[(y + j) * ysize + (x + k)] ;
                         }
                        arr_final[cnt++] = val ;

                    }
                }
            int status = match<float>(sig_array_ref , xsize , ysize , arr_final , win_range , win_range , &x_corr , &y_corr) ;
            if (status)
                {
                    LOG(ERROR) << "***Correlation Failed.***" << endl ;
                    return (EXIT_FAILURE) ;
                }
               LOG(INFO) << x_corr<< " " << y_corr << " " << (x - win_range / 2) << " " << (y - win_range / 2) << endl ;
               LOG(INFO) << "The xshift is " << (x_corr - (x - win_range / 2)) << "the yshift is " << (y_corr - (y - win_range / 2)) << endl ;
                x_corr_mf.push_back (x_corr - (x - win_range / 2)) ;
                y_corr_mf.push_back (y_corr - (y - win_range / 2)) ;
                LOG(INFO) << "Window Calculation is completed...." ;
            }
            x_corr = x_corr_mf[NO_OF_ELEMENTS_TO_CMPR / 2] ;
            y_corr = y_corr_mf[NO_OF_ELEMENTS_TO_CMPR / 2] ;
            x_corr_mf.clear () ;
            y_corr_mf.clear () ;
            LOG(INFO) << "The x_corr is " << x_corr << endl ;
            for (int i = 0 ; i < xsize ; i++)
            {
                for (int j = 0 ; j < ysize ; j++)
                {
                    int temp_x = (i + x_corr) ;
                    int temp_y = (j + y_corr) ;
                    if ((temp_x) < xsize && (temp_x) > -1 && (temp_y) > -1 && (temp_y) < ysize)
                    {
                        sig_array_temp[(int) (temp_y * ysize + temp_x)] = sig_array[j * ysize + i] ;
                        exp_array_temp[(int) (temp_y * ysize + temp_x)] = exp_array[j * ysize + i] ;

                    }
                }
            }
        }
       else if (algo_flag == 0)
        {
            LOG(INFO) <<i<< " Algorithm - find peaks" <<endl ;
          //  diff_dist = 4096.0 ;
            float temp_x1 = 0 , temp_x2 = 0 , temp_y1 = 0 , temp_y2 = 0 ;
            int cnt= 0 ;
            int size_arr = 0 ;
           // float *sig_array_peak = new float[xsize * ysize] ;
            status = findStar_algo1 (sig_array) ;
            if (status)
            {
                LOG(INFO) << endl << "***Error in finding star algorithm 1 for frame  " << infile << "  ***" << endl ;
                return (EXIT_FAILURE) ;
            }
           if (Cx.size () >= cx_ref.size ())
            {
                size_arr = Cx.size () ;
            }            
            else
            {
                size_arr = cx_ref.size () ;
            }
          //  int temp_refsize=cx_ref.size ();
         //   int temp_finsize=Cx.size ();
         //   min_size = temp_refsize<temp_finsize ? temp_refsize : temp_finsize ;
            //minimum_No_of_Stars = (int) ((100 - 20) * min_size / 100) ;
            float x_ref_arr[size_arr] , y_ref_arr[size_arr] , x_arr[size_arr] , y_arr[size_arr] ;
            float temp_x_arr[size_arr] , temp_y_arr[size_arr] ;
            for (long int index = 0 ; index < cx_ref.size () ; index++)//loop for finding the similer x-cordinates and ycordinates  from the both frames
            {
                int pairing_flag = 0 ;
                temp_x1 = cx_ref[index] ;
                temp_y1 = cy_ref[index] ;
                for (long int i = 0 ; i < Cx.size () ; i++)
                {
                    float diff_x = 0.0 , diff_y = 0.0 ;
                    temp_x2 = Cx[i] ;
                    temp_y2 = Cy[i] ;
                    diff_x = temp_x1 - temp_x2 ;
                    diff_y = temp_y1 - temp_y2 ;

                    if ((diff_x * diff_x) < diff_dist && (diff_y * diff_y) < diff_dist)// finding the similer points
                    {
                        pairing_flag = 1 ;
                        x_ref_arr[cnt] = temp_x1 ;
                        y_ref_arr[cnt] = temp_y1 ;
                        x_arr[cnt] = temp_x2 ;
                        y_arr[cnt] = temp_y2 ;
                        temp_x_arr[cnt] = (-1.0) * diff_x ;
                        temp_y_arr[cnt] = (-1.0) * diff_y ;
                        cnt++ ;
                        break ;
                       }
                }
            }
        
               if(cnt<=minimum_No_of_Stars &&cnt<10000)
               {
                
                     diff_dist =diff_dist*2;
                     cx_ref.clear ();
                     cy_ref.clear ();
                     ci_ref.clear ();
                     Cx.clear ();
                     Cy.clear ();
                     Ci.clear ();
                     restart ();
                     return(EXIT_SUCCESS);
              }
         
            if (option_LeastSquare == 1)
            {  
                 spMatrix B1 ((cnt) * 2 , 1) ;
                 spMatrix A1 ((cnt) * 2 , 3) ;
                 spMatrix X1 (3 , 1) ;
                 int temp = 0 ;
                for (int t = 0 ; t < cnt * 2 ; t = t + 2)
                {
                    B1(t , 0) = (x_arr[temp] - x_ref_arr[temp]) ; //+IMAGE_ARRAYSIZE*0.5;
                    B1 (t + 1 , 0) = (y_arr[temp] - y_ref_arr[temp]) ; //+IMAGE_ARRAYSIZE*0.5;
                    A1 (t , 0) = -1.0 * (y_ref_arr[temp] - xsize / 2) ;
                    A1 (t , 1) = 1 ;
                    A1 (t , 2) = 0 ;
                    A1 (t + 1 , 0) = (x_ref_arr[temp] - xsize / 2) ;
                    A1 (t + 1 , 1) = 0 ;
                    A1 (t + 1 , 2) = 1 ;
                    temp++ ;
                }
                X1.ApplyLeastSquare (A1 , B1) ;
                x_dx = X1 (1 , 0) ;
                y_dy = X1 (2 , 0) ;
                theta_dt = X1(0 , 0) ;
             }
            else if (option_LeastSquare == 2)
            {
                spMatrix A (cnt * 2 , 4) ;
                spMatrix B (cnt * 2 , 1) ;
                spMatrix X (4 , 1) ;
                for (int aindex = 0 ; aindex < cnt ; aindex++)
                {
                    A (2 * aindex , 0) = x_ref_arr[aindex] - xsize * 0.5 ;
                    A (2 * aindex , 1) = -1.0 * (y_ref_arr[aindex] - xsize * 0.5) ;
                    A (2 * aindex , 2) = 1.0 ;
                    A (2 * aindex , 3) = 0.0 ;
                    A (2 * aindex + 1 , 0) = (y_ref_arr[aindex] - xsize * 0.5) ;
                    A (2 * aindex + 1 , 1) = x_ref_arr[aindex] - xsize * 0.5 ;
                    A (2 * aindex + 1 , 2) = 0.0 ;
                    A (2 * aindex + 1 , 3) = 1.0 ;
                    B (2 * aindex , 0) = x_arr[aindex] - xsize * 0.5 ;
                    B (2 * aindex + 1 , 0) = y_arr[aindex] - xsize * 0.5 ;
                }
                X.ApplyLeastSquare (A , B) ;
               double theta = atan2 (X (1 , 0) , X (0 , 0)) ;
                x_dx = X (2 , 0) ;
                y_dy = X (3 , 0) ;
                theta_dt = theta ;
            }
            else if (option_LeastSquare == 3)
            {
                double a11 = 0.0 , a12 = 0.0 , a13 = 0.0 , a21 = 0.0 , a22 = 0.0 , a23 = 0.0 , a31 = 0.0 , a32 = 0.0 , a33 = 0.0 , b1 = 0.0 , b2 = 0.0 , b3 = 0.0 ;
                spMatrix A (3 , 3) ;
                spMatrix B (3 , 1) ;
                spMatrix X (3 , 1) ;
                for (int k = 0 ; k < cnt ; k++)
                {
                    /* weight used for a star*/
                    double w = 1.0 ;
                    //     w=pow( (25.0*photonbk_per_pixel+starmag[k]), 2.0 )/
                    //        (  (25.0*photonbk_per_pixel)+(0.09*starmag[k])  );
                       /* row 1 */
                    double y_mod = y_arr[k]-(xsize * 0.5) ;
                    double x_mod = x_arr[k]-(xsize * 0.5) ;
       
                    y_mod = -1.0 * y_mod ;
                    a11 = a11 + 2.0 * w ;
                    a12 = 0.0 ;
                    a13 = a13 + (2.0 * w * (y_mod)) ;
                    b1 = b1 + 2.0 * (w * temp_x_arr[k]) ;

                    /* row 2*/
                    a21 = 0.0 ;
                    a22 = a22 + 2.0 * w ;
                    a23 = a23 + (2.0 * w * (x_mod)) ;
                    b2 = b2 + 2.0 * (w * temp_y_arr[k]) ;

                    /* row 3 */
                    a31 = a31 + (2.0 * w * (y_mod)) ;
                    a32 = a32 + (2.0 * w * (x_mod)) ;
                    a33 = a33 + (2.0 * w * (pow (x_mod , 2.0) + pow (y_mod , 2.0))) ;
                    b3 = b3 + (2.0 * w * ((x_mod) * temp_y_arr[k] + (y_mod) * temp_x_arr[k])) ;

                }
                A (0 , 0) = a11 ;
                A (0 , 1) = a12 ;
                A (0 , 2) = a13 ;
                A (1 , 0) = a21 ;
                A (1 , 1) = a22 ;
                A (1 , 2) = a23 ;
                A (2 , 0) = a31 ;
                A (2 , 1) = a32 ;
                A (2 , 2) = a33 ;
                B (0 , 0) = b1 ;
                B (1 , 0) = b2 ;
                B (2 , 0) = b3 ;
                X.ApplyLeastSquare (A , B) ;
                x_dx = X (0 , 0) ;
                y_dy = X (1 , 0) ;
                theta_dt = X (2 , 0) ;
            }
          //  ofptr<<x_dx<<setw(20)<<y_dy<<setw(20)<<theta_dt<<endl;            
            float ctheta = 0.0f , stheta = 0.0f ;
            ctheta = cos (-1.0*theta_dt) ;
            stheta = sin (-1.0*theta_dt) ;
            int cnt_loop = 0 ;
            LOG(INFO) << "\nLoop Started for assigning the correction to the frames..\n" <<endl ;
            for(int i=0;i<xsize*ysize;i++){
                sig_array_temp[i]=0.0f;
                exp_array_temp[i]=0.0f;
            }

            for (int i = 0 ; i < xsize ; i++)
            {
                int x_index = i - xsize / 2 ;
                for (int j = 0 ; j < ysize ; j++)
                {
                  int y_index = j - ysize / 2 ;
                  float new_index_x = round( (x_index) * ctheta - y_index * stheta + xsize / 2)-x_dx ; //new index x
                  float new_index_y = round( (x_index) * stheta + (y_index) * ctheta + ysize / 2)-y_dy ; //new index y
               
                   if ((new_index_x) < xsize && (new_index_x) > 0 && (new_index_y) > 0 && (new_index_y) < ysize)
                    {
                        cnt_loop++ ;
                        sig_array_temp[(int) (round (new_index_y) * ysize + round (new_index_x))] =  sig_array[j * ysize + i] ;
                        exp_array_temp[(int) (round (new_index_y) * ysize + round (new_index_x))] = exp_array[j * ysize + i] ;
                    }
                }
            }
            Cx.clear () ;
            Cy.clear () ;
            Ci.clear () ;
         }
        else
        {
            LOG(ERROR) << "***Invalid choice for the algorithm***" << endl ;
            return (EXIT_FAILURE) ;
        }                                 

         for(int p=0;p<nelements;p++)
         {
         avg_sigArray[p]=avg_sigArray[p]+sig_array_temp[p]*exp_array_temp[p];
         avg_expArray[p]=avg_expArray[p]+exp_array_temp[p];         
         }                     
    }
      // ofptr.close();
   // cout<<"ppp"<<endl;
    //    exit(1);   
       for(int p=0;p<nelements;p++)
      {
        if(avg_expArray[p]>0.0)
         avg_sigArray[p]=(avg_sigArray[p]/avg_expArray[p]);
      }     
        //for converting into 600*600 frame
         float *new_sigArr= new float[FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
         float *new_expArr=new float[FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG];
      
status= ApplySubSampling (avg_sigArray,xsize,ysize,new_sigArr,FINALFRAMESIZE_REGAVG,FINALFRAMESIZE_REGAVG);
if(status)
{
LOG(ERROR)<<"ERROR in sub sampling "<<endl;
return(EXIT_FAILURE);
}
status=ApplySubSampling (avg_expArray,xsize,ysize,new_expArr,FINALFRAMESIZE_REGAVG,FINALFRAMESIZE_REGAVG);
if(status)
{
LOG(ERROR)<<"ERROR in sub sampling "<<endl;
return(EXIT_FAILURE);
}
          
        
        sprintf (outfile , "%s/%s_sig_rav.fits" , moduleoutdir , nameprefix ) ;
        fits_create_file (&outfptr2 , outfile , &status) ;
        printError (status , "Error in creating the out file" , outfile) ;
        fits_create_img (outfptr2 , FLOAT_IMG , naxis , naxes , &status) ;
        printError (status , "Error in creating the image file" , outfile) ;
        fits_write_pix (outfptr2 , TFLOAT , fpixel ,FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG,new_sigArr , &status) ;
        printError (status , "Error in writing the pixels to the output Exposure File" , outfile) ;
        fits_update_key (outfptr2 , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
        printError (status , "Error in updating the key value of the TIMEFILE" , outfile) ;
        fits_close_file (outfptr2, &status) ;
        writeUsrkeywordsFrmvect (outfile,key_records);
        if(history==YES)writeHistory (outfile,vhistorystr);
        updateKeywords (outfile,modulename);
       
        	
        fits_open_file (&outfptr2 , outfile , READWRITE , &status) ;
        printError (status , "Error in opening the input  Signal File" , outfile) ;
        fits_delete_key(outfptr2,"FRAMENO",&status);
        printError (status , "Error in DELETING  the FRAMENO" , outfile) ;
        fits_delete_key(outfptr2,"FRMTIME",&status);
        printError (status , "Error in DELETING  the FRM_TIME" , outfile) ;
        fits_close_file(outfptr2,&status);
        printError (status , "Error in closing output file" , outfile) ;
        
        sprintf (outfile , "%s/%s_exp_rav.fits" , moduleoutdir  , nameprefix ) ;
        fits_create_file (&outfptr2 , outfile , &status) ;
        printError (status , "Error in creating the out file" , outfile) ;
        fits_create_img (outfptr2 , FLOAT_IMG , naxis , naxes , &status) ;
        printError (status , "Error in creating the image file" , outfile) ;

        fits_write_pix (outfptr2 , TFLOAT , fpixel , FINALFRAMESIZE_REGAVG*FINALFRAMESIZE_REGAVG, new_expArr, &status) ;
        printError (status , "Error in writing the pixels to the output Exposure File" , outfile) ;
        fits_update_key (outfptr2 , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
        printError (status , "Error in updating the key value of the TIMEFILE" , infofile_out) ;
        fits_close_file (outfptr2, &status) ;
        printError (status , "Error in closing the file" , infofile_out) ;
       
        writeUsrkeywordsFrmvect (outfile,key_records);
        if(history==YES)writeHistory (outfile,vhistorystr);
        updateKeywords (outfile,modulename);
         
        fits_open_file (&outfptr2 , outfile , READWRITE , &status) ;
        printError (status , "Error in opening the input  Signal File" , outfile) ;
        fits_delete_key(outfptr2,"FRAMENO",&status);
        printError (status , "Error in DELETING  the FRAMENO" , outfile) ;
        fits_delete_key(outfptr2,"FRMTIME",&status);
        printError (status , "Error in DELETING  the FRM_TIME" , outfile) ;
        fits_close_file(outfptr2,&status);
        printError (status , "Error in closing output file" , outfile) ;
    return (EXIT_SUCCESS) ;
}
int uvtRegAvg::getHistory (vector<string> &vhistory)
{

    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ("P1 framelistDir=" + (string) inputdatadir) ;
    vhistory.push_back ("P3 outdir=" + (string) outdir) ;
    vhistory.push_back ("Module Output directory=" + (string) moduleoutdir) ;
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

void uvtRegAvg::StarDetectionAndCentroids (char infile[] , char outfile[] , int sqrSize , double primaryThreshold , double secondaryThreshold)
{
    LOG(INFO) << "\nInside the Star Detection method" << endl ;
    int i , j , m , n ;
    int algo_Square_Size = sqrSize ;
    int s1 , s2 ;
    s1 = sqrSize ;
    s2 = sqrSize ;

#if 1   

    /* Variable declaration for stardetection and finding centroid */
    double center_pix_val = 0.0 ; //center pixel value of algo.square
    double temp , temp1 ;
    double sum_Of_SqrPixVal ; //sum of all pixel values of algo.square
    double fcmin ; //minimum count value from four corner value of user defined square
    int x , y ;
    double a , b , c , d ; //variables for corner values of 5*5 square
    int flag ;
    double check1 ; //center pixel value - fcmin
    double check2 ; /*sum of all pixel values of 5*5 matrix - (N*fcmin)
                            where N is number of pixels(e.g if 5*5 pixels are 25,if 3*3 pixels are 9 )*/
    int c1 , c2 ;
    int sqr_x[s1][s2] , sqr_y[s1][s2] ;
    int sum_x , sum_y ;
    int temp_m , temp_n , minus_val ;
    double **sqr_pix , **image_Data , *data , *temp_image_Data ;
    //char allocate_Memory_double;
    sqr_pix = allocateMemory<double>(s1 , s2) ;

    /*Fits file variable declaration*/
    int status = 0 , hdutype ;
    long fpixel = 1 , nelements ;
    int nulval = 0 ;
    int anynull = 0 ;
    int naxis ;
    long *naxes ;
    int maxdim = 0 , bitpix ;
    int pixel = 0 , rows = 0 , height = 0 , width = 0 ;

    if (algo_Square_Size == 24)
    {
        minus_val = 12 ;
    }

    if (algo_Square_Size == 3)
    {
        minus_val = 2 ;
    }

    if (algo_Square_Size == 5)
    {
        minus_val = 3 ;
    }
#endif

#if 1   
    fitsfile* infptr ;

    /**********************************************************************/
    /*************** INPUT FILE READING  ********************/
    /**********************************************************************/
    // LOG(INFO)<<"The infile is "<<infile<<endl;
    if (fits_open_file (&infptr , infile , READONLY , &status))
        // printError(status);
        printError (status , "***Error in moving to 3rd HDU in input Event File1***") ;
    if (fits_open_image (&infptr , infile , READONLY , &status))
        printError (status , "***Error in moving to 3rd HDU in input Event File2***") ;
    //printError(status);
    //  LOG(INFO)<<"The outside of first fits conversion"<<endl;
    //*get image dimension (naxis)
    if (fits_get_img_dim (infptr , &naxis , &status))
        printError (status , "***Error in moving to 3rd HDU in input Event File3***") ;
    //  printError(status);
    //LOG(INFO)<<"The outside of first fits conversion"<<endl;
    naxes = new long[naxis] ; //**allocate memory to naxes
    naxes[0] = 0 ;
    naxes[1] = 0 ;
    maxdim = naxis ;
    //  LOG(INFO)<<"The outside of first fits conversion"<<endl;
    //** get image parameters (used for geting naxes[0] & naxes[1])
    if (fits_get_img_param (infptr , maxdim , &bitpix , &naxis , naxes , &status))
        printError (status , "***Error in moving to 3rd HDU in input Event File4***") ;
    //printError(status);
    //    LOG(INFO)<<"The outside of first fits conversion3"<<endl;
    pixel = naxes[0] ;
    rows = naxes[1] ;
    LOG(INFO) << "The rows are " << naxes[1] << endl ;
    height = naxes[0] ;
    width = naxes[1] ;
    nelements = pixel * rows ; //* number of pixels to write 

    image_Data = allocateMemory<double>(pixel , rows) ;
    temp_image_Data = new double[nelements] ;
    //    LOG(INFO)<<"The outside of first fits conversion0"<<endl;
    if (fits_read_img (infptr , TDOUBLE , fpixel , nelements , &nulval , temp_image_Data , &anynull , &status))
        printError (status , "***Error in moving to 3rd HDU in input Event File5***") ;
    //  printError(status);
    // LOG(INFO)<<"The outside of first fits conversion1"<<endl;
    if (fits_close_file (infptr , &status))
        printError (status , "***Error in moving to 3rd HDU in input Event File6***") ;
    //  printError(status);

    /* Assign one dimension array values to two dimension array*/

    /*    for(int d1=0;d1<pixel;d1++)
        {
            for(int d2=0;d2<rows;d2++)
            {
                image_Data[d1][d2]=temp_image_Data[index];
                index++;
            }
        }
     */
    //     int index1=0;
    data = new double[nelements] ;
    int temp_index = 0 , index1 = 0 ;
    for (int d1 = 0 ; d1 < pixel ; d1++)
    {
        for (int d2 = 0 ; d2 < rows ; d2++)
        {
            image_Data[d1][d2] = temp_image_Data[index1] ;
            index1++ ;
            if (image_Data[d1][d2] > 0)
            {
                data[temp_index] = image_Data[d1][d2] ;
                temp_index++ ;
            }
        }
    }

    long int p ;
    double max_pix = data[0] ;
    double min_pix = data[0] ;
    /*      for(p=1;p<temp_index;p++)
            {  
                            if(data[p]>max_pix)
                            {
                                    max_pix = data[p];
                            }
            }
     */
    for (p = 1 ; p < temp_index ; p++)
    {
        if (data[p] > max_pix)
        {
            max_pix = data[p] ;
        }
        if (data[p] < min_pix)
        {
            min_pix = data[p] ;
        }
    }

    double sum_xh , sum_yh ;
    double **Xc , **Yc ;
    int **X , **Y ;
    double **Center_PixelValue ;
    Center_PixelValue = allocateMemory<double>(pixel , rows) ;
    Xc = allocateMemory<double>(pixel , rows) ;
    Yc = allocateMemory<double>(pixel , rows) ;
    X = allocateMemory<int>(pixel , rows) ;
    Y = allocateMemory<int>(pixel , rows) ;

    for (int a = 0 ; a < pixel ; a++)
    {
        for (int b = 0 ; b < rows ; b++)
        {
            Xc[a][b] = 0.0 ;
            Yc[a][b] = 0.0 ;
            X[a][b] = 0 ;
            Y[a][b] = 0 ;
            Center_PixelValue[a][b] = 0.0 ;
        }

    }
    FILE *optr1 = fopen (outfile , "w") ;
    int count ;
    int xy_index = 0 ;
    //  LOG(INFO)<<"The number of pixels are "<<pixel<<endl;
    for (i = 0 ; i < pixel ; i++)
    {
        if (i == 5000)
            LOG(INFO) << i ;
        //   LOG(INFO)<<"The outside of first fits conversion21"<<endl;
        count = 0 ;

        for (j = 0 ; j < rows ; j++)
        {
            // LOG(INFO)<<"The outside of first fits conversion22"<<endl;
            //              exit(1);
            sum_Of_SqrPixVal = 0 ;
            if (i + s1 < height && j + s2 < width)
            {
                //  LOG(INFO)<<"The outside of first fits conversion23"<<endl;

                for (m = i , c1 = 0 ; m < i + s1 ; m++ , c1++)
                {
                    // LOG(INFO)<<"The outside of first fits conversion24"<<endl;
                    for (n = j , c2 = 0 ; n < j + s2 ; n++ , c2++)
                    {
                        //**read required pixel values and x,y locations
                        sqr_pix[c1][c2] = image_Data[m][n] ;
                        sqr_x[c1][c2] = m ;
                        sqr_y[c1][c2] = n ;
                        sum_Of_SqrPixVal = sum_Of_SqrPixVal + sqr_pix[c1][c2] ;

                    }
                }
                //read square corner value & define minimum of them as "fcmin"
                x = m - 1 ;
                y = n - 1 ;

                temp_m = x - minus_val ;
                temp_n = y - minus_val ;

                a = image_Data[i][j] ;
                b = image_Data[x][j] ;
                c = image_Data[i][y] ;
                d = image_Data[x][y] ;

                //find fcmin
                temp = (a <= b) ? a : b ;
                temp1 = (temp <= c) ? temp : c ;
                fcmin = (temp1 <= d) ? temp1 : d ;

                center_pix_val = image_Data[temp_m][temp_n] ; //*finding center pix value
                check1 = (center_pix_val - fcmin) ;
                check2 = sum_Of_SqrPixVal - (fcmin * (s1 * s2)) ; //*for step7 

                if (fcmin != 0)
                {
                    flag = 0 ;
                    if (check1 > primaryThreshold)
                    {

                        flag = check_MaxCenterPixVal (center_pix_val , sqr_pix , algo_Square_Size) ; //check center pixel value is local maximum

                        if (flag == 1)
                        {
                            if (check2 > secondaryThreshold)
                            {

                                //centroid of square
                                sum_xh = 0.0 ;
                                sum_x = 0 ;
                                sum_yh = 0.0 ;
                                sum_y = 0 ;
                                for (int ii = 0 ; ii < s1 ; ii++)
                                {
                                    for (int jj = 0 ; jj < s2 ; jj++)
                                    {
                                        sum_xh = sum_xh + (sqr_x[ii][jj] * sqr_pix[ii][jj]) ;
                                        sum_yh = sum_yh + (sqr_y[ii][jj] * sqr_pix[ii][jj]) ;
                                    }
                                }

                                Xc[i][j] = sum_xh / sum_Of_SqrPixVal ; //divide by intensity
                                Yc[i][j] = sum_yh / sum_Of_SqrPixVal ;
                                X[i][j] = temp_m + 1 ;
                                Y[i][j] = temp_n + 1 ;
                                Center_PixelValue[i][j] = center_pix_val ;
                                xy_index++ ;
                            }
                        }
                    }
                }
            }//end of if(i+s1 && j+s2 <512)
        }
    }
    double **temp_Xc , **temp_Yc , **temp_centroid ;
    int **temp_X , **temp_Y ;
    temp_Xc = allocateMemory<double>(16 , 16) ;
    temp_Yc = allocateMemory<double>(16 , 16) ;
    temp_X = allocateMemory<int>(16 , 16) ;
    temp_Y = allocateMemory<int>(16 , 16) ;
    temp_centroid = allocateMemory<double>(16 , 16) ;
    int r1 , tt1 , ss1 , counter ;
   for (i = 0 ; i < pixel ; i = i + 16)
    {        
        for (j = 0 ; j < rows ; j = j + 16)
        { // if(j==0) LOG(INFO)<<"The outside of first fits conversion2q"<<endl;

            counter = 0 ;
            for (tt1 = i , r1 = 0 ; tt1 < i + 16 ; tt1++ , r1++)
            {

                for (ss1 = j , c1 = 0 ; ss1 < j + 16 ; ss1++ , c1++)
                {
                    // LOG(INFO)<<i<<" "<<j<<endl;	 

                    temp_Xc[r1][c1] = Xc[tt1][ss1] ;

                    temp_Yc[r1][c1] = Yc[tt1][ss1] ;
                    temp_X[r1][c1] = X[tt1][ss1] ;
                    temp_Y[r1][c1] = Y[tt1][ss1] ;
                    temp_centroid[r1][c1] = Center_PixelValue[tt1][ss1] ;
                }
           }

            //#if 0v

            for (r1 = 0 ; r1 < 16 ; r1++)
            {
                //   LOG(INFO)<<"Inside the 16loop "<<endl;
                for (c1 = 0 ; c1 < 16 ; c1++)
                {
                    if (r1 == 16 / 2)
                    {
                        if (c1 == 16 / 2)
                        {
                            if (temp_centroid[r1][c1] > 0.0)
                            {
                                //counter++;
                                ////if(counter ==8)
                                fprintf (optr1 , "%f %d %d %f %f\n" , temp_centroid[r1][c1] , temp_X[r1][c1] , temp_Y[r1][c1] , temp_Xc[r1][c1] , temp_Yc[r1][c1]) ;

                            }
                        }
                    }

                }
            }
            //#endif

        }
    }
    //exit(1);
    free (temp_image_Data) ;
    free (data) ;

    for (int i = 0 ; i < 16 ; i++)
    {

        free (temp_centroid[i]) ;
        free (temp_X[i]) ;
        free (temp_Y[i]) ;
        free (temp_Xc[i]) ;
        free (temp_Yc[i]) ;
    }
    free (temp_Xc) ;
    free (temp_Yc) ;
    free (temp_X) ;
    free (temp_Y) ;
    free (temp_centroid) ;

    for (int d2 = 0 ; d2 < rows ; d2++)
    {

        free (image_Data[d2]) ;
        free (Xc[d2]) ;
        free (X[d2]) ;
        free (Yc[d2]) ;
        free (Y[d2]) ;
    }

    free (image_Data) ;
    free (Xc) ;
    free (Yc) ;
    free (X) ;
    free (Y) ;
    LOG(INFO) << "Output Centroid File  " << outfile << "  is created." << endl ;
    fclose (optr1) ;


#endif   
    LOG(INFO) << "outside of  the Star Ditection method" << endl ;
}

int uvtRegAvg::check_MaxCenterPixVal (double center , double **sqr , int sqr_size)
{
    int flag_count = 0 ;
    int sqr_elem = sqr_size*sqr_size ;
    //int flag_set=sqr_elem-1;
    int center_x = (sqr_size / 2) - 1 ;
    int center_y = (sqr_size / 2) - 1 ;
    for (int d1 = 0 ; d1 < sqr_size ; d1++)
    {
        for (int d2 = 0 ; d2 < sqr_size ; d2++)
        {
            //if(d1!= center_x && d2!= center_y && center>sqr[d1][d2])
            if (center >= sqr[d1][d2])
            {
                flag_count = flag_count + 1 ;
            }
            else
            {
                flag_count = flag_count + 0 ;
            }
        }
    }
    if (flag_count == sqr_elem)
        return 1 ;
    else return 0 ;

}

int uvtRegAvg::correlate (float *array1 , int h1 , int w1 , float *array2 , int h2 , int w2 , float xshift , float yshift , int bp_p1 , int bp_p2)
{

    int h = h1 + h2 ;
    int w = w1 + w2 ;

    LOG(INFO) << endl << "Height :" << h << "   Width:" << w ;

    float **image1 = new float*[h] ;
    float **image2 = new float*[h] ;
    //	float **corr = new float*[h];
    for (int i = 0 ; i < h ; i++)
    {
        image1[i] = new float[w] ;
        image2[i] = new float[w] ;
        //		corr[i]= new float[w];
    }

    LOG(INFO) << endl << "Memory allocated for arrays\n" ;

    for (int i = 0 ; i < h ; i++)
    {
        for (int j = 0 ; j < w ; j++)
        {
            image1[i][j] = 0.0 ;
            image2[i][j] = 0.0 ;
            //			corr[i][j]=0.0;
        }
    }

    LOG(INFO) << endl << "Reading images....\n" ;
    //reading image1
    if (bp_p1 <= 8)
    {
    }
    else if (bp_p1 <= 16)
    {
    }
    else if (bp_p1 <= 32)
    {
        int index = 0 ;
        for (int i = (h - h1) / 2 ; i < (h + h1) / 2 ; i++)
            for (int j = (w - w1) / 2 ; j < (w + w1) / 2 ; j++)
                image1[i][j] = array1[index++] ;

    }
    else
    {
        LOG(INFO) << endl << "Bits per pixel " << bp_p1 << " not supported" << endl ;
    }

    //reading image 2
    if (bp_p2 <= 8)
    {
    }
    else if (bp_p2 <= 16)
    {
    }
    else if (bp_p2 <= 32)
    {
        int index = 0 ;
        for (int i = (h - h2) / 2 ; i < (h + h2) / 2 ; i++)
            for (int j = (w - w2) / 2 ; j < (w + w2) / 2 ; j++)
                image2[i][j] = array2[index++] ;

    }
    else
    {
        LOG(INFO) << endl << "Bits per pixel " << bp_p2 << " not supported" << endl ;
    }

    LOG(INFO) << endl << "Completed reading images\n" ;

    LOG(INFO) << endl << "Correlating images...\n" ;

    fftw_complex *im1_fft , *im2_fft ;
    fftw_plan p1 , p2 ;

    fftw_complex *corr_fft ;

    double *im1 = new double[h * w] ;
    double *im2 = new double[h * w] ;
    double *corr = new double[h * w] ;

    im1_fft = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * h * (w / 2 + 1)) ;
    im2_fft = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * h * (w / 2 + 1)) ;
    corr_fft = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * h * (w / 2 + 1)) ;

    for (int i = 0 ; i < h ; i++)
    {
        for (int j = 0 ; j < w ; j++)
        {
            im1[i * h + j] = image1[i][j] * pow (-1 , i + j) ;
            im2[i * h + j] = image2[i][j] * pow (-1 , i + j) ;
        }
    }

    //  for(int i=0;i<h;i++) delete[] image1[i],image2[i];
    //       delete[] image1,image2;
    free (image1) ;
    free (image2) ;

    ofstream fim1 ("/tmp/im1.bin" , ios::binary) ;
    ofstream fim2 ("/tmp/im2.bin" , ios::binary) ;
    for (int i = 0 ; i < h * w ; i++)
    {
        float val1 = (float) im1[i] ;
        fim1.write ((char *) &val1 , sizeof (float)) ;
        float val2 = (float) im2[i] ;
        fim2.write ((char *) &val2 , sizeof (float)) ;
    }
    fim1.close () ;
    fim2.close () ;

    p1 = fftw_plan_dft_r2c_2d (h , w , im1 , im1_fft , FFTW_ESTIMATE) ;
    p2 = fftw_plan_dft_r2c_2d (h , w , im2 , im2_fft , FFTW_ESTIMATE) ;

    fftw_execute (p1) ;
    fftw_execute (p2) ;

    LOG(INFO) << endl << "FFT computation completed" << endl ;

    //	//writing fft1
    //	ofstream fft1("/tmp/fft1.bin",ios::binary);
    //
    //        for(int i=0;i<h*w;i++){
    //                float val = (float)sqrt(im1_fft[i][0]*im1_fft[i][0]+im1_fft[i][1]*im1_fft[i][1]);
    //                fft1.write((char *)&val,sizeof(float));
    //        }
    //        fft1.close();
    //
    //	//writing fft2
    //	ofstream fft2("/tmp/fft2.bin",ios::binary);
    //
    //        for(int i=0;i<h*w;i++){
    //                float val = (float)sqrt(im2_fft[i][0]*im2_fft[i][0]+im2_fft[i][1]*im2_fft[i][1]);
    //                fft2.write((char *)&val,sizeof(float));
    //        }
    //        fft2.close();


    //taking conjugate of second image
    for (int i = 0 ; i < h * (w / 2 + 1) ; i++)
        im2_fft[i][1] = im2_fft[i][1]*(-1) ;

    //element by element multiplication for correlation
    for (int i = 0 ; i < h * (w / 2 + 1) ; i++)
    {
        int x = i / h ;
        int y = i % h ;
        corr_fft[i][0] = (im1_fft[i][0] * im2_fft[i][0] - im1_fft[i][1] * im2_fft[i][1]) ;
        corr_fft[i][1] = (im1_fft[i][0] * im2_fft[i][1] + im1_fft[i][1] * im2_fft[i][0]) ;
    }

    fftw_plan p3 ;
    p3 = fftw_plan_dft_c2r_2d (h , w , corr_fft , corr , FFTW_ESTIMATE) ;

    double *fftimage = new double[h * w] ;
    for (int i = 0 ; i < h * w ; i++) fftimage[i] = 0 ;
    int index = 0 ;
    for (int i = 0 ; i < h ; i++)
    {
        for (int j = 0 , k = w - 1 ; j < (w / 2 + 1) ; j++ , k--)
        {
            fftimage[i * h + j] = sqrt (corr_fft[index][0] * corr_fft[index][0] + corr_fft[index][1] * corr_fft[index][1]) ;
            fftimage[i * h + k] = fftimage[i * h + j] ;
            index++ ;
        }
    }

    ofstream fftcorr ("/tmp/corrfft.bin" , ios::binary) ;

    for (int i = 0 ; i < h * w ; i++)
    {
        float val = (float) fftimage[i] ;
        fftcorr.write ((char *) &val , sizeof (float)) ;
    }
    fftcorr.close () ;

    delete[] fftimage ;

    fftw_execute (p3) ;

    fftw_destroy_plan (p1) ;
    fftw_destroy_plan (p2) ;
    fftw_destroy_plan (p3) ;
    fftw_free (im1_fft) ;
    fftw_free (im2_fft) ;
    fftw_free (corr_fft) ;

    delete[] im1 , im2 ;

    double max = 0 ;
    float x , y ;
    //ofstream fout(argv[9],ios::binary);
    for (int i = 0 ; i < h * w ; i++)
    {
        float val = (float) corr[i] ;
        if (corr[i] > max)
        {
            max = corr[i] ;
            y = i / h ;
            x = i % h ;
        }
        val = val * pow (-1 , (i / h + i % h)) ;
        val = val / (h * w) ;
        //fout.write((char *)&val,sizeof(float));
    }

    //fout.close();

    delete[] corr ;
    if (x > h / 2 && y > w / 2)
    {
        LOG(INFO) << "INSIDE " ;
        x = x - h / 2 ;
        y = y - w / 2 ;
    }
    else if (x > h / 2 && y < w / 2)
    {
        x = x - h / 2 ;
        y = y + w / 2 ;
    }
    else if (x < h / 2 && y < w / 2)
    {
        x = x + h / 2 ;
        y = y + w / 2 ;
    }
    else if (x < h / 2 && y > w / 2)
    {
        x = x + h / 2 ;
        y = y - w / 2 ;
    }

    //testing	
    LOG(INFO) << endl << "Peak observed at (" << x << "," << y << ")" ;
    x = x - h / 2 ;
    y = y - w / 2 ;
    if (y < 0 && x < 0)
    {

        y = y*-1 ;
        x = x*-1 ;
    }
    LOG(INFO) << "The X_shift " << x << " " << "The Y_shift " << y << endl ;
    xshift = x ;
    yshift = y ;
    LOG(INFO) << endl ;
    return 0 ;
}

template<class T>
int uvtRegAvg::match (T *searchwindow , int h1 , int w1 , T *templatewindow , int h2 , int w2 , float *x , float *y)
{
    int h = h1 + h2 ;
    int w = w1 + w2 ;

    LOG(INFO) << endl << "Height :" << h << "   Width:" << w ;

    fftw_complex *im1_fft , *im2_fft ;
    fftw_plan p1 , p2 ;

    fftw_complex *corr_fft ;

    double *im1 = new double[h * w] ;
    if (im1 == NULL)
    {
        LOG(INFO) << endl << "Out of Memory Error" << endl ;
        return (EXIT_FAILURE) ;
    }

    double *im2 = new double[h * w] ;
    if (im2 == NULL)
    {
        LOG(INFO) << endl << "Out of Memory Error" << endl ;
        return (EXIT_FAILURE) ;
    }

    double *corr = new double[h * w] ;
    if (corr == NULL)
    {
        LOG(INFO) << endl << "Out of Memory Error" << endl ;
        return (EXIT_FAILURE) ;
    }

    im1_fft = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * h * (w / 2 + 1)) ;
    if (im1_fft == NULL)
    {
        LOG(INFO) << endl << "Out of Memory Error" << endl ;
        return (EXIT_FAILURE) ;
    }

    im2_fft = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * h * (w / 2 + 1)) ;
    if (im2_fft == NULL)
    {
        LOG(INFO) << endl << "Out of Memory Error" << endl ;
        return (EXIT_FAILURE) ;
    }

    corr_fft = (fftw_complex*) fftw_malloc (sizeof (fftw_complex) * h * (w / 2 + 1)) ;
    if (corr_fft == NULL)
    {
        LOG(INFO) << endl << "Out of Memory Error" << endl ;
        return (EXIT_FAILURE) ;
    }

    for (int i = 0 ; i < h ; i++)
    {
        for (int j = 0 ; j < w ; j++)
        {
            //copying serach window to padded double array
            if (i < h1 && j < w1)
                im1[i * w + j] = searchwindow[i * w1 + j] * pow (-1 , i + j) ;
            else
                im1[i * w + j] = 0 ;
            //copying template window to padded double array
            if (i < h2 && j < w2)
                im2[i * w + j] = templatewindow[i * w2 + j] * pow (-1 , i + j) ;
            else
                im2[i * w + j] = 0 ;
        }
    }

    ofstream fim1 ("/tmp/im1.bin" , ios::binary) ;
    ofstream fim2 ("/tmp/im2.bin" , ios::binary) ;
    for (int i = 0 ; i < h * w ; i++)
    {
        float val1 = (float) im1[i] ;
        fim1.write ((char *) &val1 , sizeof (float)) ;
        float val2 = (float) im2[i] ;
        fim2.write ((char *) &val2 , sizeof (float)) ;
    }
    fim1.close () ;
    fim2.close () ;

    p1 = fftw_plan_dft_r2c_2d (h , w , im1 , im1_fft , FFTW_ESTIMATE) ;
    p2 = fftw_plan_dft_r2c_2d (h , w , im2 , im2_fft , FFTW_ESTIMATE) ;

    fftw_execute (p1) ;
    fftw_execute (p2) ;

    fftw_destroy_plan (p1) ;
    fftw_destroy_plan (p2) ;
    //fftw_free(p1);
    //   fftw_free(p2);
    delete[] im1 , im2 ;
    LOG(INFO) << endl << "FFT computation completed" << endl ;

    //	//writing fft1
    //	ofstream fft1("/tmp/fft1.bin",ios::binary);
    //
    //        for(int i=0;i<h*w;i++){
    //                float val = (float)sqrt(im1_fft[i][0]*im1_fft[i][0]+im1_fft[i][1]*im1_fft[i][1]);
    //                fft1.write((char *)&val,sizeof(float));
    //        }
    //        fft1.close();
    //
    //	//writing fft2
    //	ofstream fft2("/tmp/fft2.bin",ios::binary);
    //
    //        for(int i=0;i<h*w;i++){
    //                float val = (float)sqrt(im2_fft[i][0]*im2_fft[i][0]+im2_fft[i][1]*im2_fft[i][1]);
    //                fft2.write((char *)&val,sizeof(float));
    //        }
    //        fft2.close();


    //taking conjugate of second image
    for (int i = 0 ; i < h * (w / 2 + 1) ; i++)
        im2_fft[i][1] = im2_fft[i][1]*(-1) ;

    //element by element multiplication for correlation
    for (int i = 0 ; i < h * (w / 2 + 1) ; i++)
    {
        int x = i / h ;
        int y = i % h ;
        corr_fft[i][0] = (im1_fft[i][0] * im2_fft[i][0] - im1_fft[i][1] * im2_fft[i][1]) ;
        corr_fft[i][1] = (im1_fft[i][0] * im2_fft[i][1] + im1_fft[i][1] * im2_fft[i][0]) ;
    }

    fftw_free (im1_fft) ;
    fftw_free (im2_fft) ;

    double *fftimage = new double[h * w] ;
    for (int i = 0 ; i < h * w ; i++) fftimage[i] = 0 ;
    int index = 0 ;
    for (int i = 0 ; i < h ; i++)
    {
        for (int j = 0 , k = w - 1 ; j < (w / 2 + 1) ; j++ , k--)
        {
            fftimage[i * h + j] = sqrt (corr_fft[index][0] * corr_fft[index][0] + corr_fft[index][1] * corr_fft[index][1]) ;
            fftimage[i * h + k] = fftimage[i * h + j] ;
            index++ ;
        }
    }

    ofstream fftcorr ("/tmp/corrfft.bin" , ios::binary) ;

    for (int i = 0 ; i < h * w ; i++)
    {
        float val = (float) fftimage[i] ;
        fftcorr.write ((char *) &val , sizeof (float)) ;
    }
    fftcorr.close () ;

    delete[] fftimage ;


    fftw_plan p3 ;
    p3 = fftw_plan_dft_c2r_2d (h , w , corr_fft , corr , FFTW_ESTIMATE) ;

    fftw_execute (p3) ;
    fftw_destroy_plan (p3) ;
    // fftw_free(p3);
    fftw_free (corr_fft) ;

    double max = 0 ;

    ofstream fout ("/tmp/corr.bin" , ios::binary) ;
    for (int i = 0 ; i < h1 ; i++)
    {
        for (int j = 0 ; j < w1 ; j++)
        {
            float val = (float) corr[i * w + j] ;
            if (corr[i * w + j] > max)
            {
                max = corr[i * w + j] ;
                *x = j ;
                *y = i ;
            }
            val = val * pow (-1 , (i / h + i % h)) ;
            val = val / (h * w) ;
            fout.write ((char *) &val , sizeof (float)) ;
        }
    }

    fout.close () ;

    delete[] corr ;
    LOG(INFO) << endl << "Peak observed at (" << *x << "," << *y << ")" ;
    LOG(INFO) << endl ;

    return (EXIT_SUCCESS) ;
} //end of match function

bool compare (struct Star star1 , struct Star star2)
{
    return (star1.intensity > star2.intensity) ;
}

int uvtRegAvg::findStar_algo1 (float *inputArray) //algorithm for finding the peaks
{

    Fx.clear () ;
    Fy.clear () ;
    Fval.clear () ;
    Rx.clear () ;
    Ry.clear () ;
    Rval.clear () ;
    Cx.clear () ;
    Cy.clear () ;
    Ci.clear () ;
    int r , c ;
     float *temp_array;
     vector<float> array_temp;     
    if(datainfo.getModeFlag ()==PC){
        
           array_temp.clear ();
           
      for (int i=0;i<xsize*ysize;i++){
        
        if(inputArray[i]!=0.0f){
            array_temp.push_back (inputArray[i]);
            
        }
    }
         temp_array  = new float[array_temp.size ()];
           for (int in=0;in<array_temp.size ();in++)
           {
               temp_array[in]=array_temp[in];     
             
           }
    }
          
label:
    Fval.clear () ;
    Fx.clear () ;
    Fy.clear () ;
    Rx.clear () ;
    Ry.clear () ;
    Rval.clear () ;
   
    

    if (sd_mul_factor < 0)
    {
        LOG(ERROR) << "***SD_MULTI_FACTOR is <0***" << endl ;
        return (EXIT_FAILURE) ;
    }
     double thr=0;
    if(datainfo.getModeFlag ()==PC){
//        cout<<"111 "<<array_temp.size ()<<endl;
//        exit(1);
         thr = getSD (temp_array , array_temp.size()) * sd_mul_factor ;
        // cout<<getSD (temp_array , array_temp.size())<<endl;
        
    }
    else{
          thr = getSD (inputArray , xsize * ysize) * sd_mul_factor ;
    }
    LOG(ERROR) << endl << "\nThreshold for first cut peaks is   " << thr ;

    for (int i = 0 ; i < xsize * ysize ; i++)
    {
        r = (i / xsize) ;
        c = (i % xsize) ;

        if (inputArray[i] > thr)
        {
            Fval.push_back (inputArray[i]) ;
            Fx.push_back (c) ; //x is for column
            Fy.push_back (r) ; //y is for row
        }
    }

    LOG(INFO) << "SIGMA  Factor::" << sd_mul_factor << endl ;
    LOG(INFO) << " Size of First cut Peaks  " << Fy.size () << endl ;

    if (Fy.size () < minimum_No_of_Stars)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG(ERROR) << sd_mul_factor << " less than 0!!!! " << endl ;
            return (EXIT_FAILURE) ;
        }
        goto label ;
        //        LOG(INFO) << endl << "First cut peaks detected : " << Fy.size () << endl ;
        //        LOG(INFO) << endl << "***No peaks found ***" << endl ;
    }

//    for (int i = 0 ; i < xsize * ysize ; i++)
//        peakImage[i] = 0 ;
//
//    for (int i = 0 ; i < Fy.size () ; i++)
//        peakImage[Fy[i] * xsize + Fx[i]] = Fval[i] ;

    
     //if winsize is even, make it odd
    if (refine_Winsize % 2 == 0)
        refine_Winsize = refine_Winsize - 1 ;

    LOG(INFO) << endl << "Using window size : " <<refine_Winsize << " for refining peaks " ;

    //refined peaks
    vector<int> Tx , Ty ;
    vector<float> Tval ;

    Tx = Fx ;
    Ty = Fy ;
    Tval = Fval ;

    /*refining peaks logic
    refined Window size is for the refined  peaks.
   Refined peaks are found by  making window around each of the star(i.e first cut peaks)  and  finding brightest star among that window.*/
    int start_r , end_r , start_c , end_c ;

    for (int i = 0 ; i < Fx.size () ; i++)
    {
        start_r = Ty[i] - refine_Winsize / 2 ;
        end_r = Ty[i] + refine_Winsize / 2 ;
        start_c = Tx[i] - refine_Winsize / 2 ;
        end_c = Tx[i] + refine_Winsize / 2 ;
        if (start_r < 0) start_r = 0 ;
        if (end_r >= ysize) end_r = ysize - 1 ;
        if (start_c < 0) start_c = 0 ;
        if (end_c >= xsize) end_c = xsize - 1 ;
        int max = 0 ;
        for (int k = start_r ; k <= end_r ; k++)
        {
            for (int l = start_c ; l <= end_c ; l++)
            {

                if (inputArray[k * xsize + l] > max)
                {
                    max = inputArray[k * xsize + l] ;
                    Tx[i] = l ;
                    Ty[i] = k ;
                    Tval[i] = inputArray[k * xsize + l] ;
                } //  end of if block 
            } //end of l loop
        } //end of  k  loop
    } // end of i loop

    /*--------------Refining peaks completed----------------*/

    float *arr_refine = new float[xsize * ysize] ; //to store refined peaks

    for (int i = 0 ; i < xsize * ysize ; i++)
        arr_refine[i] = 0 ;

    if (xsize == 0 || ysize == 0)
    {
        LOG(ERROR) << "***Divide by Zero***" << endl ;
        return (EXIT_FAILURE) ;
    }
    for (int i = 0 ; i < Ty.size () ; i++)
        arr_refine[Ty[i] * xsize + Tx[i]] = Tval[i] ; //overwriting the same place..

    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;

    for (int i = 0 ; i < xsize * ysize ; i++)
    {
        if (arr_refine[i] != 0)
        {
            Rx.push_back ((i % xsize)) ;
            Ry.push_back ((i / xsize)) ;
            Rval.push_back (arr_refine[i]) ;
        }
    }
    LOG(INFO) << "Number of final peaks is " << Rval.size () << endl ;
//    if (Rval.size () > centroidlimit)
//    {
//        LOG(ERROR) << endl << "Number of refined peaks found:" << Rval.size () << endl ;
//        LOG(ERROR) << endl << "***Number of peaks found is greater than centroidlimit*** \n***Run with greater threshold value***" << endl ;
//        return (EXIT_FAILURE) ;
//    }

    if (Ry.size () < minimum_No_of_Stars)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG(ERROR) << sd_mul_factor << " less than 0!!!! " << endl ;
            return (EXIT_FAILURE) ;
        }
        goto label ;
    }
    /**method for  find Centroid within the Stars**/
    doCentroiding (Rx , Ry , centroid_Winsize , inputArray , ysize , xsize) ;
    delete[] arr_refine ;
    return (EXIT_SUCCESS) ;
}

void uvtRegAvg::doCentroiding (vector<int> &X , vector<int> &Y , int centroidwindow , float *arr , int h , int w)
{

    Cx.clear () ;
    Cy.clear () ;
    Ci.clear () ;

    float x , y , val = 0 ;
    //   LOG(INFO) << "THe Centroid Window" << centroidwindow << endl;
    int num = centroidwindow*centroidwindow ;
    double sum_x = 0 , sum_y = 0 , sum = 0 ;
    //  LOG(INFO) << "The X.size is " << X.size() << endl;
    for (int i = 0 ; i < X.size () ; i++)
    {
        sum_x = 0 ;
        sum_y = 0 , sum = 0 ;
        for (int j = -1 * (centroidwindow / 2) ; j <= (centroidwindow / 2) ; j++)
        {
            for (int k = -1 * (centroidwindow / 2) ; k <= (centroidwindow / 2) ; k++)
            {
                x = X[i] ;
                y = Y[i] ;

                val = arr[(Y[i] + j) * w + (X[i] + k)] ;
                //  val=arr[(Y[i]+j)*w+(X[i]+k)];
                //   LOG(INFO)<<"The value::"<<x+j<<y+j<<arr[(Y[i]+j)*w+(X[i]+k)]<<endl;
                sum_x = sum_x + (x + k) * val ;
                sum_y = sum_y + (y + j) * val ;
                sum = sum + val ;
            }
        }
        if (sum <= 0)
        {
            LOG(INFO) << endl << "Sum of intensites for (" << X[i] << " , " << Y[i] << ")  is <=0" << endl ;
            LOG(INFO) << endl << "\nDivide by zero error\n" ;
            exit (EXIT_FAILURE) ;
        }
        Cx.push_back ((float) sum_x / (float) sum) ;
        Cy.push_back ((float) sum_y / (float) sum) ;
        Ci.push_back ((float) sum) ;

    }

    //sorting the list
    float temp , tx , ty ;
    for (int i = 0 ; i < Ci.size () ; i++)
    {
        //LOG(INFO)<<endl<<i;
        for (int j = Ci.size () - 1 ; j > i ; j--)
        {
            if (Ci[j - 1] < Ci[j])
            {
                swap1 (Ci[j] , Ci[j - 1]) ;
                swap1 (Cx[j] , Cx[j - 1]) ;
                swap1 (Cy[j] , Cy[j - 1]) ;
            }
        }
    }
}


int uvtRegAvg::restart (){
    LOG(INFO)<<"\nuvtComputeDrift will be restarted with newer "<<diff_dist<<" distance as a NB "<<endl;
    uvtRegAvgProcess ();
    return(EXIT_SUCCESS);
}

int uvtRegAvg :: ApplySubSampling (float* inputarray, int in_xsize, int in_ysize, float* outputarray, int out_xsize, int out_ysize)
{
    if(in_xsize%out_xsize!=0){
        LOG(ERROR)<<"Can't sub sampling with  this input array "<<endl;
        return(EXIT_FAILURE);
    }
    float devide_term_x=in_xsize/out_xsize;   
     float devide_term_y=in_ysize/out_ysize;   
    int cnt_win=0;
    float sum_win=0.0f;
    int index_finalArray=0;
    for(int temp_x=0;temp_x<in_xsize;temp_x=temp_x+devide_term_x)
    {
        for(int temp_y=0;temp_y<in_ysize;temp_y=temp_y+devide_term_y)
        {
            if(index_finalArray>out_xsize*out_ysize){
                LOG(ERROR)<<"Array is out of bound, EXEED to "<<out_xsize<<" !!!"<<endl;
                return(EXIT_FAILURE);
            }
            cnt_win=0;
            sum_win=0.0f;
                   for(int i=temp_y;i<temp_y+devide_term_y;i++)
                   {        
                        for(int  j=temp_x;j<temp_x+devide_term_x;j++)
                        {
                            if(inputarray[j*in_xsize+i]!=0)
                            {
                                sum_win=sum_win+inputarray[j*in_xsize+i];
                                cnt_win++;                                
                            }
                                                        
                        }
                   }
            if(cnt_win!=0)
            {   
               outputarray[index_finalArray++]=(float)(sum_win/cnt_win);
            }
            else{
                 outputarray[index_finalArray++]=0.0f;
            }
            
            
        }
       
    }        
    
    return(EXIT_SUCCESS);
    
}
