/* 
 * File:   uvtDetectStar.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */

#include<pil.h>
#include<fitsio.h>
#include<iostream>
#include<unistd.h>

#include<dirent.h>  //Accessing Directory
#include<string.h> 
#include<cstdlib>
#include<stdio.h>
#include<fstream>
#include<uvtDetectStar.h>
#include<uvtUtils.h>
#include<algorithm>
#include<glog/logging.h>
#include<macro_def.h>
#include<spMatrix.h>

#define MODULENAME "uvtDetectStar"

bool compare (struct Star vect1 , struct Star vect2) ;

//typedef struct Star
//{
//    float intensity ;
//   float  x , y ;
//};
bool compare (struct Star vect1 , struct Star vect2)
{
return (vect1.intensity > vect2.intensity ) ;
    //return (vect1.intensity > vect2.intensity) ;
}

float findMedianValue (vector<float> &input , int noOfValues);
uvtDetectStar::uvtDetectStar ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    strcpy (starDir , "Star") ;
    strcpy (centroidDir , "Centroid") ;
}


uvtDetectStar::~ uvtDetectStar () {
 }


int uvtDetectStar::read (int argc , char** argv)
{
    int status = 0 ;


    status = readParams (argc , argv , 3 , FNAME , "inputdatadir" , inputdatadir , FNAME , "outdir" , outdir , INT , "algo_flag" , &algo_flag) ;
    if (status) return (EXIT_FAILURE) ;
    if (algo_flag == 1 || algo_flag == 3 || algo_flag == 4)
    {
        status = readParams (argc , argv , 5 , INT , "Nacc" , &Nacc , REAL4 , "background_fact" , &backgrnd_fact , REAL4 , "threshold" , &sd_multi_factor_default , REAL4 , "minimum_targetedstars" , &thr_intensity_refinement , INT , "refine_Window" , &refine_Winsize) ;
        if (status) return (EXIT_FAILURE) ;
    }
    else  if (algo_flag == 2)
    {
        status = readParams (argc , argv , 3 , INT , "StarDetectionCentroid_square_size" , &algo_Square_Size , REAL4 , "StarDetectionCentroid_threshold" , &primary_threshold_Val , REAL4 , "StarDetectionCentroid_secondary_threshold" , &secondary_threshold_Val) ;
        if (status) return (EXIT_FAILURE) ;
    }
    else
    {
        LOG (INFO) << "***algo_flag must be 1 or 2***" << endl ;
        return (EXIT_FAILURE) ;
    }
    status = readParams (argc , argv , 5 , INT , "centroid_Window" , &centroid_Winsize , INT , "Window_nb" , &win_search , BOOL , "clobber" , &clobber , BOOL , "history" , &history , STRING , "mode" , mode) ;
    if (status) return (EXIT_FAILURE) ;
    //    if (PIL_OK != (status = PILInit (argc , argv)))
    //    {
    //        LOG(INFO) << "***Error Initializing PIL***" ;
    //        return status ;
    //    }
    //    if (PIL_OK != (status = PILGetFname ("inputdatadir" , inputdatadir)))
    //    { //change to inputdir
    //        LOG(INFO) << endl << "***Error reading input directory***" ;
    //        return status ;
    //    }
    //    if (PIL_OK != (status = PILGetFname ("outdir" , outdir)))
    //    {
    //        LOG(INFO) << endl << "***Error reading output directory name***" ;
    //        return status ;
    //    }
    //    if (PIL_OK != (status = PILGetInt ("algo_flag" , &algo_flag)))
    //    {
    //        LOG(INFO) << endl << "***Error reading algo_flag ***" ;
    //        return status ;
    //    }
    //    if (algo_flag == 1)
    //    {
    //        if (PIL_OK != (status = PILGetReal4 ("threshold" , &sd_multi_factor_default)))
    //        {
    //            LOG(INFO) << endl << "***Error reading output directory name***" ;
    //            return status ;
    //        }
    //
    //        if (PIL_OK != (status = PILGetInt ("minimum_targetedstars" , &minimum_No_of_Stars)))
    //        {
    //            LOG(INFO) << endl << "***Error reading algo_flag ***" ;
    //            return status ;
    //        }
    //        if (PIL_OK != (status = PILGetInt ("refine_Window" , &refine_Winsize)))
    //        {
    //            LOG(INFO) << endl << "***Error reading refine window size ***" ;
    //            return status ;
    //        }
    //    }
    //    else if (algo_flag == 2)
    //    {
    //
    //        if (PIL_OK != (status = PILGetInt ("StarDetectionCentroid_square_size" , &algo_Square_Size)))
    //        {
    //            LOG(INFO) << endl << "***Error reading refine window size ***" ;
    //            return status ;
    //        }
    //        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroid_threshold" , &primary_threshold_Val)))
    //        {
    //            LOG(INFO) << endl << "***Error reading refine window size ***" ;
    //            return status ;
    //        }
    //        if (PIL_OK != (status = PILGetReal4 ("StarDetectionCentroid_secondary_threshold" , &secondary_threshold_Val)))
    //        {
    //            LOG(INFO) << endl << "***Error reading refine window size ***" ;
    //            return status ;
    //        }
    //    }
    //    else
    //    {
    //        LOG(INFO) << "***algo_flag must be 1 or 2***" << endl ;
    //        return (EXIT_FAILURE) ;
    //    }
    //
    //    if (PIL_OK != (status = PILGetInt ("centroid_Window" , &centroid_Winsize)))
    //    {
    //        LOG(INFO) << endl << "***Error reading centroid window size ***" ;
    //        return status ;
    //    }
    //
    //    if (PIL_OK != (status = PILGetBool ("clobber" , &clobber)))
    //    {
    //        LOG(INFO) << "***Error Reading clobber:" << clobber << "***" ;
    //        return status ;
    //    }
    //    if (PIL_OK != (status = PILGetBool ("history" , &history)))
    //    {
    //        LOG(INFO) << "***Error Reading history parameter:" << history << "***" ;
    //        return status ;
    //    }
    //    if (PIL_OK != (status = PILGetString ("mode" , mode)))
    //    {
    //        LOG(INFO) << "***Error Reading mode parameter:" << history << "***" ;
    //        return status ;
    //    }
    //
    //    PILClose (status) ;

    return (EXIT_SUCCESS) ;
}


int uvtDetectStar::read (char* inputdatadir , char* outdir , int algo_flag , float threshold , int refine_window , int centroid_window , float num_min_stars , float prm_thr , float sec_thr , int algo_square , int windw_search , int clobber , int history)
{

    strcpy (this->inputdatadir , inputdatadir) ;
    this->centroid_Winsize = centroid_window ;

    this->win_search = windw_search ;
    this->algo_flag = algo_flag ;
    if (algo_flag == 1)
    {
        this->refine_Winsize = refine_window ;
        this->sd_multi_factor_default = threshold ;
        this->thr_intensity_refinement = num_min_stars ;
    }
    else if (algo_flag == 2)
    {
        this->primary_threshold_Val = prm_thr ;
        this->secondary_threshold_Val = sec_thr ;
        this->algo_Square_Size = algo_square ;
    }
    //  this->centroidlimit = centroidlimit ;
    strcpy (this->outdir , outdir) ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}


void uvtDetectStar::display ()
{
    LOG (INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG (INFO) << "             UVT DETECT STAR PARAMETERS              " << endl ;
    LOG (INFO) << "------------------------------------------------------------------------" ;
    LOG (INFO) << endl << "Input Frame List Directory              : " << inputdatadir ;
    LOG (INFO) << endl << "Output Directory                               : " << outdir ;
    LOG (INFO) << endl << "Threshold Value                               : " << sd_multi_factor_default ;
    LOG (INFO) << endl << "Algo Flag                                          : " << algo_flag ;
    LOG (INFO) << endl << "Centroid Window size                               : " << centroid_Winsize ;
    if (algo_flag == 1)
    {
        LOG (INFO) << endl << "Refined Window size                               : " << refine_Winsize ;
        LOG (INFO) << endl << "SD multiplication Factor              : " << sd_multi_factor_default ;
        //LOG(INFO) << endl << "maximum Centroid Limit                               : " << centroidlimit ;
        LOG (INFO) << endl << "Minimum number of Stars                              : " << thr_intensity_refinement ;

    }
    else if (algo_flag == 2)
    {
        LOG (INFO) << endl << "primary threshold value                               : " << primary_threshold_Val ;
        LOG (INFO) << endl << "secondary threshold value                               : " << secondary_threshold_Val ;
        LOG (INFO) << endl << "Window Square size                                : " << algo_Square_Size ;
    }
    if (clobber == 1)
        LOG (INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG (INFO) << endl << "Overwrite                                         : NO" ;
    if (history == 1)
        LOG (INFO) << endl << "History                                             : YES" ;
    else
        LOG (INFO) << endl << "History                                              : NO" ;

    LOG (INFO) << endl << "--------------------------------------------------------------------\n" ;

}


int uvtDetectStar::uvtDetectStarProcess ()
{
    string cmd ;
    nframes = 0 ;
    int status = 0 ;

    if (! DirExists (inputdatadir))
    {
        LOG (INFO) << "***Input directory " << inputdatadir << "  not found***" << endl ;
        return (EXIT_FAILURE) ;
    }
    sprintf (moduleoutdir , "%s/%s" , outdir , modulename) ;

    if (DirExists (moduleoutdir) && clobber == YES)
    {
        LOG (INFO) << endl << "Directory exists and clobber=yes" << endl ;
        cmd = (string) "rm -rf " + (string) moduleoutdir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists (moduleoutdir) && clobber == NO)
    {
        LOG (INFO) << endl << moduleoutdir << "  already exists " << endl ;
        LOG (INFO) << endl << "Use clobber=yes for overwriting" << endl ;
        return (EXIT_FAILURE) ;
    }
    /**Shell command for creating the output directory**/
    cmd = "mkdir -p " + (string) moduleoutdir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;

    string  infofilename = searchFile (inputdatadir , ".info") ;

    if (infofilename == " ")
    {
        LOG (ERROR) << "***Error in finding the info file " << inputdatadir << " ***" << endl ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , infofilename.c_str()) ;

    //open Information file
    fitsfile *finfo_in , *finfo_out ;
    LOG (INFO) << "Information File ::" << infofile_in << endl ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input information file" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU of input info file" , infofile_in) ;
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file

    xsize = datainfo.getXsize () ; //xsize- x dimention of the image
    ysize = datainfo.getYsize () ; //ysize-y dimention of the image
    if (xsize <= 0 || ysize <= 0)
    {
        LOG (INFO) << endl << "***Invalid xsize/ysize***\n" ;
        return (EXIT_FAILURE) ;
    }

    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX" , infofile_in) ; //for creating name for output information file

    //creating output information file
    sprintf (infofile_out , "%s/%s_sc.info" , moduleoutdir , nameprefix) ;
    LOG (INFO) << " output Information File " << infofile_out << endl ;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    char *ttype[] = {"StarFileList" , "CentroidFileList"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table in output information file" , infofile_out) ;
    datainfo.write (finfo_out) ; //writing basic data information

    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the key value of the NAMEPRFX" , infofile_out) ; //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file" , infofile_out) ;
    fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
    printError (status , "Error in reading the key value of the NFILES" , infofile_in) ;
    LOG (INFO) << endl << "Number of files :" << nframes ;
    //reading frame names from information file into vector
    /*
     *sigframelist- array for storing the signal frame names
     * expframelist-array for storing the exposure frame names
     */
    sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
    expoframelist = allocateMemory<char>(nframes , NAMESIZE) ;
    fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) sigframelist , NULL , &status) ;
    printError (status , "Error in reading the list of signal frame" , infofile_in) ;
    fits_read_col (finfo_in , TSTRING , 2 , 1 , 1 , nframes , NULL , (void *) expoframelist , NULL , &status) ;
    printError (status , "Error in reading the list of exposure frame" , infofile_in) ;
    fits_read_key (finfo_in , TSTRING , "SIGDIR" , sigframedir , NULL , &status) ;
    printError (status , "Error in reding the key value of the SIGDIR" , infofile_in) ;
    fits_read_key (finfo_in , TSTRING , "EXPDIR" , expoframedir , NULL , &status) ;
    printError (status , "Error in reading the key value of the EXPDIR" , infofile_in) ;
    fits_close_file (finfo_in , &status) ;
    printError (status , "Error in closing the input information file" , infofile_in) ;
    //method for finding the stars and centroids
    if (detectStars ()) return (EXIT_FAILURE) ;

    LOG (INFO) << "Updating the output information file" << endl ;

    /**updating the keywords NFILES,STARDIR,CENTDIR
     *STARDIR-directory path for output  star directory 
     * EXPDIR-directory path for ourpur centroid directory
     
     */
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the  input information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU of input information file" , infofile_out) ;
    fits_update_key (finfo_out , TINT , "NFILES" , &nframes , "Number of frames" , &status) ;
    printError (status , "Error in updating the key value of the NFILES " , infofile_out) ;
    fits_update_key (finfo_out , TSTRING , "STARDIR" , starDir , "Star directory" , &status) ;
    printError (status , "Error in updating the key value of the STARDIR" , infofile_out) ;
    fits_update_key (finfo_out , TSTRING , "CENTDIR" , centroidDir , "Centroid Directory" , &status) ;
    printError (status , "Error in updating the key value of the CENTDIR" , infofile_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file" , infofile_out) ;


    freeMemory (sigframelist , nframes , NAMESIZE) ;
    freeMemory (expoframelist , nframes , NAMESIZE) ;
    return (EXIT_SUCCESS) ;
}


int uvtDetectStar::detectStars ()
{

    LOG (INFO) << endl << "\nStarted  Finding Stars ..." << endl ;
    fitsfile *finfo_out ;
    //creating star directory
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s" , moduleoutdir , starDir) ;
    /**Shell command for creating the Star directory**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;
    LOG (INFO) << endl << dir << " directory created" << endl ;
    sprintf (dir , "%s/%s" , moduleoutdir , centroidDir) ;
    /**Shell command for creating the Centroid Directory**/
    cmd = "mkdir -p " + (string) dir ;
    /**Executing  the shell command**/
    system (cmd.c_str ()) ;
    LOG (INFO) << endl << dir << " directory created" << endl ;
    LOG (INFO) << endl << "\nTotal number of frames - " << nframes << endl ;
    int status = 0 ;
    float *framedata = new float[xsize * ysize] ; //for signal frame
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char errstr[512] ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    vector<string> vhistorystr ;

    vector<string> starfilelist , centroidfilelist ;
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;
    char temp_path_txt[FLEN_FILENAME] ;
    sprintf (temp_path_txt , "%s/peaks_withoutMeanSubtrcted.txt" , moduleoutdir) ;
    //    ofstream ofptr ;
    //    ofptr.open (temp_path_txt , ios::out) ;
    //    ofptr << "Frame_no" << setw (20) << "First_cut peaks" << setw (20) << "Refined_peaks" << endl ;
    //    ofptr << "===================================================" << endl ;

    //float *peakImage = new float[xsize * ysize] ;
    fitsfile *fptr , *fout ;
    int tfields = 3 ;
    cnt_det = 1 ;
    char *ttype[] = {"X" , "Y" , "Intensity"} ;
    char *tform[] = {"I" , "I" , "E"} ;
    LOG (INFO) << "\nFinding Stars and Centroid  process started for " << nframes << endl ;
    //loop for number of frames
    /**
     * if algo_flag is 1 than use SD based algorithm(findStar_algo1)
     * if algo flag is 2 than use JOE's algorithm(findStar_algo2)
     * @return 
     */
    //    if(xsize%PIX_PADD_SIZE==0)
    //    {
    int mult_term = xsize / PIX_PADD_SIZE ;
    refine_Winsize = refine_Winsize*mult_term ;
    centroid_Winsize = centroid_Winsize*mult_term ;
    //  }
    //  else
    // {
    // cout<<"Invalid Xsize !!!  xsize = "<<xsize<<", not a factor of 600 "<<endl;
    //  return(EXIT_SUCCESS);
    // }
    if (history == YES) getHistory (vhistorystr) ;
    first_frameFlag = FALSE ;
    vector<float> cx_ref , cy_ref , ci_ref ;
    for (int i = 0 ; i < nframes ; i ++)
    {

        sprintf (errstr , "Error at iteration number %d" , i) ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[i]) ;
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input File" , infile) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
        printError (status , "Error in reading the pixels from the input File" , infile) ;
        fits_read_key (fptr , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ;
        printError (status , "Error reading the key value of the FRAMRNO " , infile) ;
        fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
        printError (status , "Error reading the key value of the FRAMRTIME" , infile) ;
        //copy level-1 keywords to the vector
        if (i == 0)
        {
            copyUsrkeywrdsTovect (fptr , key_records) ;
        }
        if (algo_flag == 1)
        {
            sd_mul_factor = sd_multi_factor_default ;
            LOG (INFO) << "\nSD multiplication Factor :" << sd_multi_factor_default ;
            // method for finding firstcut peaks,refined peaks and centroids
            /*enable  following 2 lines  for new algorithm */
            if (first_frameFlag == FALSE)
            {
                status = findStar_algo1 (framedata) ;
                if (status)
                {
                    LOG (ERROR) << endl << "***Error in finding star algorithm 1(RMS based) for frame  " << infile << "  ***" << endl ;
                    return (EXIT_FAILURE) ;
                }
                cx_ref = Cx ;
                cy_ref = Cy ;
                ci_ref = Ci ;
                first_frameFlag = TRUE ;
                /*enable following lines  for new algorithm*/
            }

            else
            {
                status = matchStars_TorefFrame (cx_ref , cy_ref , ci_ref , framedata) ;
                if (status)
                {
                    LOG (ERROR) << endl << "***Error in matching stars wrt first(Reference frame)   " << infile << "  ***" << endl ;
                    return (EXIT_FAILURE) ;
                }

                cx_ref.clear () ;
                cy_ref.clear () ;
                ci_ref.clear () ;
                cx_ref = Cx ;
                cy_ref = Cy , ci_ref = Ci ;
            }

        }
        else if (algo_flag == 2)
        {
            sd_mul_factor = sd_multi_factor_default ;
            //method for finding refined peaks and centroids(this algorithm directly detect refined peaks)
            status = findStar_algo2 (framedata) ;
            if (status)
            {
                LOG (ERROR) << endl << "***Error in finding star algorithm 2(JOE's ) for frame  " << infile << "  ***" << endl ;
                return (EXIT_FAILURE) ;
            }
        }
        else if (algo_flag == 3)
        {
            sd_mul_factor = sd_multi_factor_default ;
            status = findStar_algo3 (framedata ) ;
            if (status)
            {
                LOG (ERROR) << "Error in find Star algorithm" << endl ;
                LOG(ERROR)<<"CRASH FAILED TO FIND STAR (uvtDetectStar.cpp)";
                return (EXIT_FAILURE) ;
            }
        }
        else if (algo_flag == 4)
        {

            vector<float> cint_background , diff_intensity , final_cx , final_cy , final_cint , diff_int_x , diff_int_y ;
            vector <float> cx_background , cy_background ;
            spMatrix A (512 * 512 , 6) ;
            spMatrix B (512 * 512 , 1) ;
            spMatrix X (6 , 1) ;
            //   LOG(INFO)<<"ERR";
            for (int i = 44 ; i < xsize - 44 ; i ++)
            {
                for (int j = 44 ; j < ysize - 44 ; j ++)
                {
                    if (framedata[i * xsize + j] != INVALID_PIX_VALUE)
                    {
                        A (i , 0) = 1 ;
                        A (i , 1) = i ;
                        A (i , 2) = j ;
                        A (i , 3) = i*i ;
                        A (i , 4) = j*j ;
                        A (i , 5) = i*j ;
                        B (i , 0) = framedata[i * xsize + j] ;
                    }
                }
            }
            X.ApplyLeastSquare (A , B) ;
            for (int i = 44 ; i < xsize - 44 ; i ++)
            {
                for (int j = 44 ; j < ysize - 44 ; j ++)
                {
                    if (framedata[i * xsize + j] != INVALID_PIX_VALUE)
                    {
                        cint_background.push_back (X (0 , 0) + X (1 , 0) * i + X (2 , 0) * j + X (3 , 0) * i * i + X (4 , 0) * j * j + X (5 , 0) * i * j) ;
                        cx_background.push_back (i) ;
                        cy_background.push_back (j) ;
                    }
                    else
                    {
                        cint_background.push_back (INVALID_PIX_VALUE) ;
                        cx_background.push_back (INVALID_PIX_VALUE) ;
                        cy_background.push_back (INVALID_PIX_VALUE) ;
                    }
                }
            }
            for (int i = 44 ; i < xsize - 44 ; i ++)
            {
                for (int j = 44 ; j < ysize - 44 ; j ++)
                {
                    if (framedata[i * xsize + j] != INVALID_PIX_VALUE  && cint_background[i * xsize + j] != INVALID_PIX_VALUE  )
                    {
                        diff_intensity.push_back (framedata[i * IMG_DIM_DI + j] - cint_background[i * IMG_DIM_DI + j]) ;
                        diff_int_x.push_back (i) ;
                        diff_int_y.push_back (j) ;
                    }
                    else
                    {
                        diff_intensity.push_back (INVALID_PIX_VALUE) ;
                        diff_int_x.push_back (INVALID_PIX_VALUE) ;
                        diff_int_y.push_back (INVALID_PIX_VALUE) ;
                    }
                }
            }
            //    for(int i=0;i<cx.size ();i++)
            //    {    
            //    diff_intensity.push_back (cint[i]-cint_background[i]);    
            //    }
            double sd_diff_centroid_pix = getSD (diff_intensity.data () , diff_intensity.size ()) ;
            //  cx.clear ();cy.clear ();cint.clear ();
            for (int i = 44 ; i < xsize - 44 ; i ++)
            {
                for (int j = 44 ; j < ysize - 44 ; j ++)
                {
                    if (framedata[i * xsize + j] > cint_background[i * xsize + j] + 2 * sd_diff_centroid_pix)
                    {
                        final_cx.push_back (i) ;
                        final_cy.push_back (j) ;
                        final_cint.push_back (framedata[i * xsize + j]) ;
                    }
                }
            }



            //    diff_intensity.clear ();cint_background.clear ();
            //    cx=final_cx;
            //    cy=final_cy;
            //    cint=final_cint;

            Cx = cx_background ;
            Cy = cy_background ;
            Ci = cint_background ;

            final_cx.clear () ;
            final_cy.clear () ;
            final_cint.clear () ;
        }
        else
        {
            LOG (ERROR) << "***Invalid input for the algo flag value***" << endl ;
        }
        //creating output star File
        //   ofptr << i + 1 << setw (20) << Fx.size () << setw (20) << Rx.size () << endl ;
        //setting path for output star frame
        sprintf (outfile , "%s/%s/%s_t.%f_f%d_sc_star.fits" , moduleoutdir , starDir , nameprefix , frametime , frameno) ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output File " , outfile) ;

        if (algo_flag == 1)
        {
            /*to enable peaks  image file  writing.
            fits_create_img (fout , bitpix , naxis , naxes , &status) ;
            printError (status , "Error in creating the image for the output File " , outfile) ;
            //write first cut peaks to image
            fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , peakImage , &status) ;
            printError (status , "Error Writing the pixels to the output File " , outfile) ;
             */

            /**creating table for first cut peaks
             * write first cut peak's X,Y and intensity values to the output star file
             **/
            fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Firstcut Peaks" , &status) ;
            printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
            fits_write_col (fout , TINT , 1 , 1 , 1 , Fx.size () , (void *) Fx.data () , &status) ;
            printError (status , "Error Writing the firstcut pixel's X-cordinates" , outfile) ;
            fits_write_col (fout , TINT , 2 , 1 , 1 , Fy.size () , (void *) Fy.data () , &status) ;
            printError (status , "Error Writing the firstcut pixel's Y-cordinates" , outfile) ;
            fits_write_col (fout , TFLOAT , 3 , 1 , 1 , Fval.size () , (void *) Fval.data () , &status) ;
            printError (status , "Error Writing the firstcut pixel's Intensity" , outfile) ;
        }

        /**Creating the table for refined peaks
         * write refined  peak's X,Y and intensity values to the output star file
         **/

        fits_create_tbl (fout , BINARY_TBL , 0 , tfields , ttype , tform , NULL , "Refined Peaks" , &status) ;
        printError (status , "Error creating the table for the FirstCut pixels " , outfile) ;
        fits_write_col (fout , TINT , 1 , 1 , 1 , Rx.size () , (void *) Rx.data () , &status) ;
        printError (status , "Error Writing the Refined pixels's X-cordinates" , outfile) ;
        fits_write_col (fout , TINT , 2 , 1 , 1 , Ry.size () , (void *) Ry.data () , &status) ;
        printError (status , "Error Writing the Refined pixels Y-cordinates" , outfile) ;
        fits_write_col (fout , TFLOAT , 3 , 1 , 1 , Rval.size () , (void *) Rval.data () , &status) ;
        printError (status , "Error Writing the Refined  pixel's Intensity" , outfile) ;

        //writing level-1 keywords  and history to the output star file
        //copyUserKeywords(fptr,fout);
        writeUsrkeywordsFrmvect (outfile , key_records) ;
        //       fits_movabs_hdu (fout , 1 , NULL , &status) ;
        //       printError (status , "Error moving to 2nd HDU in the output information file" , infofile_out) ;
        fits_update_key (fout , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ;
        printError (status , "***Error writing Creator***") ;
        fits_update_key (fout , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
        printError (status , "***Error writing Creator***") ;
        //       fits_movabs_hdu (fout , 2 , NULL , &status) ;
        //       printError (status , "Error moving to 2nd HDU in the output information file" , infofile_out) ;
        if (history == YES) writeHistory (outfile , vhistorystr) ;

        fits_movabs_hdu (fout , 1 , NULL , &status) ;
        printError (status , "Error in moving to 2nd HDU of input information file" , infofile_out) ;
        writeCommonKeywords (fout , modulename) ; //write origin ,checksum ,date and creator to the output star file
        fits_close_file (fout , &status) ;
        printError (status , "Error in Closing the output  Star  File " , outfile) ;


        const char *temp = basename (outfile) ;

        starfilelist.push_back ((string) temp) ;
        //setting path for output centroid file
        sprintf (outfile , "%s/%s/%s_t%f_f%d_sc_centroid.fits" , moduleoutdir , centroidDir , nameprefix , frametime , frameno) ;
        temp = basename (outfile) ;
        centroidfilelist.push_back ((string) temp) ;

        //creating output centroid file 
        fitsfile *fout1 ;
        fits_create_file (&fout1 , outfile , &status) ;
        printError (status , "Error creating the output File for Centroids" , outfile) ;

        char *tform2[] = {"E" , "E" , "E"} ;
        fits_create_tbl (fout1 , BINARY_TBL , 0 , tfields , ttype , tform2 , NULL , "Centroids" , &status) ;
        printError (status , errstr) ;
        //write Centroid's X,Y and intensity value to the output centroid file
        fits_write_col (fout1 , TFLOAT , 1 , 1 , 1 , Cx.size () , (void *) Cx.data () , &status) ;
        printError (status , "Error Writing the Centroid's X-cordinates" , outfile) ;

        fits_write_col (fout1 , TFLOAT , 2 , 1 , 1 , Cy.size () , (void *) Cy.data () , &status) ;
        printError (status , "Error Writing the Centroid's Y-cordinates" , outfile) ;
        fits_write_col (fout1 , TFLOAT , 3 , 1 , 1 , Ci.size () , (void *) Ci.data () , &status) ;
        printError (status , "Error Writing the Centroid's Intencity" , outfile) ;

        //write level-1 keywords & history from vector to the output star  file.
        //copyUserKeywords (fptr,fout1);
        writeUsrkeywordsFrmvect (outfile , key_records) ;
        //        fits_movabs_hdu (fout1 , 1 , NULL , &status) ;
        //       printError (status , "Error moving to 2nd HDU in the output information file" , infofile_out) ;
        fits_update_key (fout1 , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ;
        printError (status , "***Error writing Creator***") ;
        fits_update_key (fout1 , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
        printError (status , "***Error writing Creator***") ;
        //       fits_movabs_hdu (fout1 , 2 , NULL , &status) ;
        //       printError (status , "Error moving to 2nd HDU in the output information file" , infofile_out) ;
        if (history == YES) writeHistory (outfile , vhistorystr) ;
        writeCommonKeywords (fout1 , modulename) ; //write origin,creator,checksum and date to the output centroid file

        fits_close_file (fout1 , &status) ;
        printError (status , "Error in Closing the output  Centroid File " , outfile) ;
        fits_close_file (fptr , &status) ;
        cout << endl << "Star detected for written=" << i << "  Remaining files= " << nframes - i << " \r" ;

    }
    //    delete[] peakImage ;
    // ofptr.close () ;
    LOG (INFO) << "\nWrite listing of  output file names to the output  information file" << endl ;
    /**write the listing of frames to the output information file**/
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error opening the  output information File " , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error moving to 2nd HDU in the output information file" , infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , starfilelist.size () , (void *) starfilelist.data () , &status) ;
    printError (status , "Error Writing the column of the starfilelist " , infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 2 , 1 , 1 , centroidfilelist.size () , (void *) centroidfilelist.data () , &status) ;
    printError (status , "Error Writing the column of the centroidfilelist " , infofile_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the outout information file " , infofile_out) ;
    //write level-1 keywords and history  to the output information file.
    writeUsrkeywordsFrmvect (infofile_out , key_records) ;
    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    return (EXIT_SUCCESS) ;
}


int uvtDetectStar::getHistory (vector<string> &vhistory)
{
    int cnt = 0 ;
   // char *user = getlogin () ;
    char s_algo[10] ;
    char s_cent_win[10] ;
  //  string str = "Module run by " + (string) user ;
   // vhistory.push_back (str) ;
    sprintf (s_algo , "%d" , algo_flag) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " P1 inputdatadir =" + (string) inputdatadir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Algorithm Used   =" + (string) s_algo) ;
    if (algo_flag == 1)
    {
        char sthreshold[10] ;
        char s_refine_win[10] ;
        char min_stars[10] ;
        sprintf (sthreshold , "%f" , sd_multi_factor_default) ;
        sprintf (s_refine_win , "%d" , refine_Winsize) ;
        sprintf (min_stars , "%d" , thr_intensity_refinement) ;
        vhistory.push_back ((string) getSerialNo (cnt) + " SD multiplication Factor =" + (string) sthreshold) ;
        vhistory.push_back ((string) getSerialNo (cnt) + " Refine window Size  =" + (string) s_refine_win) ;
        vhistory.push_back ((string) getSerialNo (cnt) + " minimum number of Stars Detected  =" + (string) min_stars) ;
    }
    else if (algo_flag == 2)
    {
        char algo_sqr_size[50] ;
        char pri_Thr[50] ;
        char sec_Thr[50] ;
        sprintf (algo_sqr_size , "%d" , algo_Square_Size) ;
        sprintf (pri_Thr , "%f" , primary_threshold_Val) ;
        sprintf (sec_Thr , "%f" , secondary_threshold_Val) ;
        vhistory.push_back ((string) getSerialNo (cnt) + " Algorithm Square Size =" + (string) algo_sqr_size) ;
        vhistory.push_back ((string) getSerialNo (cnt) + " Primary Threshold value  =" + (string) pri_Thr) ;
        vhistory.push_back ((string) getSerialNo (cnt) + " Secondary Threshold value  =" + (string) sec_Thr) ;
    }
    sprintf (s_cent_win , "%d" , centroid_Winsize) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Centroid Window Size  = " + (string) s_cent_win) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " outdir=" + (string) outdir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Module Output directory=" + (string) moduleoutdir) ;
    if (clobber == YES)
        vhistory.push_back ((string) getSerialNo (cnt) + " clobber=yes") ;
    else
        vhistory.push_back ((string) getSerialNo (cnt) + " clobber=no") ;
    if (history == YES)
        vhistory.push_back ((string) getSerialNo (cnt) + " history=yes") ;
    else
        vhistory.push_back ((string) getSerialNo (cnt) + " history=no") ;
    vhistory.push_back ("Parameter List END") ;
    return (EXIT_SUCCESS) ;
}


int uvtDetectStar::check_MaxCenterPixVal (double center , double **sqr , int sqr_size)
{
    int flag_count = 0 ;
    int sqr_elem = sqr_size*sqr_size ;
    //int flag_set=sqr_elem-1;
    int center_x = (sqr_size / 2) - 1 ;
    int center_y = (sqr_size / 2) - 1 ;
    for (int d1 = 0 ; d1 < sqr_size ; d1 ++)
    {
        for (int d2 = 0 ; d2 < sqr_size ; d2 ++)
        {
            //if(d1!= center_x && d2!= center_y && center>sqr[d1][d2])
            if (center >= sqr[d1][d2])
                flag_count = flag_count + 1 ;
        }
    }
    if (flag_count == sqr_elem)
        return 1 ;
    else return 0 ;

}


int uvtDetectStar::matchStars_TorefFrame (vector<float> &X , vector<float> &Y , vector<float> &intensity , float *inputArray)
{

    star_track.clear () ;
    Star star1 ;
    Cx.clear () ;
    Cy.clear () ;
    Ci.clear () ;
    vector<int>  tx , ty ;
    vector<float> tint ;
    int centroidwindow = win_search * xsize / 600 ;
    //cout<<centroidwindow<<endl;exit(1);
    float x , y ;
    float temp_high_x , temp_high_y ;
    double mean_img = getmean (inputArray , xsize * ysize) ;
    double sd_img = getSD (inputArray , xsize * ysize) ;
    float max_intensity_pixel = 0.0 ;

    for (int i = 0 ; i < X.size () ; i ++)
    {

        if (X[i] == - 9999 || Y[i] == - 9999)
        {
            temp_high_x = - 9999 , temp_high_y = - 9999 , max_intensity_pixel = 0.0 ;
        }

        else
        {
            x = X[i] ;
            y = Y[i] ;


            // cout<<"X:: "<<x<<" "<<"Y::"<<" "<<y<<endl;
            max_intensity_pixel = 0.0 ;
            temp_high_x = - 9999 , temp_high_y = - 9999 ;
            for (int j = x - (centroidwindow / 2) ; j <= x + (centroidwindow / 2) ; j ++)
            {
                for (int k = y - (centroidwindow / 2) ; k <= y + (centroidwindow / 2) ; k ++)
                {
                    if (j < 0 || k < 0 || j > xsize - 1 || k > ysize - 1)
                    {
                        cout << "Index out of range " << i << " " << x << " " << y << " " << j << " " << k << endl ;
                        exit (1) ;
                    }
                    if (max_intensity_pixel < inputArray[k * xsize + j])
                    {
                        max_intensity_pixel = inputArray[k * xsize + j] ;
                        temp_high_x = j ;
                        temp_high_y = k ;
                    }
                }
            }

        }
                     star1.x=temp_high_x;
                     star1.y=temp_high_y;
                     star1.intensity=max_intensity_pixel;
//        if (max_intensity_pixel > (mean_img + sd_mul_factor * sd_img))
//        {
//            star1.x = temp_high_x ;
//            star1.y = temp_high_y ;
//            star1.intensity = max_intensity_pixel ;
//        }
//        else
//        {
//            max_intensity_pixel = 0.0 ;
//            temp_high_x = - 9999 , temp_high_y = - 9999 ;
//            star1.x = temp_high_x ;
//            star1.y = temp_high_y ;
//            star1.intensity = max_intensity_pixel ;
//        }

        star_track.push_back (star1) ;
    }

    //sort (star_track.begin (),star_track.end (),compare);

    tx.clear () ;
    ty.clear () ;
    tint.clear () ;

    for (int i = 0 ; i < star_track.size () ; i ++)
    {
        tx.push_back (star_track[i].x) ;
        ty.push_back (star_track[i].y) ;
        tint.push_back (star_track[i].intensity) ;
        //LOG(INFO)<<star_track[i].x<<" "<<star_track[i].y<<" "<<star_track[i].intensity;
    }  
    doCentroiding (tx , ty , centroid_Winsize ,  inputArray , ysize , xsize) ;
    return (EXIT_SUCCESS) ;
}
//{
//    float mean_Ofimage=0.0;
//     if (xsize == 0 || ysize == 0)
//    {
//        LOG(ERROR) << "***Divide by Zero***" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//    Fx.clear () ;
//    Fy.clear () ;
//    Fval.clear () ;
//    Rx.clear () ;
//    Ry.clear () ;
//    Rval.clear () ;
//    Cx.clear () ;
//    Cy.clear () ;
//    Ci.clear () ;
//    int r , c ;
//     float *temp_array;
//     vector<float> array_temp;     
//     
//      if(datainfo.getModeFlag ()==PC)
//    {
//       array_temp.clear ();
//      for (int i=0;i<xsize*ysize;i++)
//      {
//        
//        if(inputArray[i]!=0.0f){
//            array_temp.push_back (inputArray[i]);
//            
//        }
//      }
//         temp_array  = new float[array_temp.size ()];
//           for (int in=0;in<array_temp.size ();in++)
//           {
//               temp_array[in]=array_temp[in];     
//             
//           }
//    }
//    label:
//    Fval.clear () ;
//    Fx.clear () ;
//    Fy.clear () ;
//    Rx.clear () ;
//    Ry.clear () ;
//    Rval.clear () ;
//   
//    
//
//    if (sd_mul_factor < 0)
//    {
//        LOG(ERROR) << "***SD_MULTI_FACTOR is <0***" << endl ;
//        return (EXIT_FAILURE) ;
//    }
//     double thr=0;
//     double sd_temp=0.0f;
// 
//   if(datainfo.getModeFlag ()==PC)
//    {
//        sd_temp=getSD (temp_array , array_temp.size ());
//  mean_Ofimage=getmean (temp_array,array_temp.size ());
//   thr =  sd_temp* sd_mul_factor ;
//     
//    }
//    else
//    {
//       
//         sd_temp=getSD (inputArray , xsize * ysize);
//  mean_Ofimage=getmean (inputArray,xsize*ysize);
//  thr =  sd_temp* sd_mul_factor ;
//    }
//    return(EXIT_SUCCESS);
//}

//


int uvtDetectStar::findStar_algo1 (float *inputArray ) //algorithm for finding the peaks
{
    float mean_Ofimage = 0.0 ;
    if (xsize == 0 || ysize == 0)
    {
        LOG (ERROR) << "***Divide by Zero***"  ;
        return (EXIT_FAILURE) ;
    }
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
    float *temp_array ;
    vector<float> array_temp ;

    if (datainfo.getModeFlag () == PC)
    {
        array_temp.clear () ;
        for (int i = 0 ; i < xsize * ysize ; i ++)
        {

            if (inputArray[i] != 0.0f)
            {
                array_temp.push_back (inputArray[i]) ;

            }
        }
        temp_array  = new float[array_temp.size ()] ;
        for (int in = 0 ; in < array_temp.size () ; in ++)
        {
            temp_array[in] = array_temp[in] ;

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
        LOG (ERROR) << "***SD_MULTI_FACTOR is <0***"  ;
        return (EXIT_FAILURE) ;
    }
    double thr = 0 ;
    double sd_temp = 0.0f ;

    if (datainfo.getModeFlag () == PC)
    {
        sd_temp = getSD (temp_array , array_temp.size ()) ;
        mean_Ofimage = getmean (temp_array , array_temp.size ()) ;
        //   thr =  sd_temp* sd_mul_factor ;

    }
    else
    {

        sd_temp = getSD (inputArray , xsize * ysize) ;
        mean_Ofimage = getmean (inputArray , xsize * ysize) ;
        //thr =  sd_temp* sd_mul_factor ;
    }


    // sd_temp=getSD (inputArray , xsize * ysize);
    // mean_Ofimage=getmean (inputArray,xsize*ysize);
    //    thr =  sd_temp* sd_mul_factor ;
    //added
    
   // thr = Nacc * backgrnd_fact + sd_temp* sd_mul_factor ;//Nacc 
    thr = backgrnd_fact + sd_temp* sd_mul_factor ;//Nacc 
    //thr = mean_Ofimage+sd_temp* sd_mul_factor ;
    LOG (INFO) << "Threshold for first cut peaks is   " << Nacc << "*" << backgrnd_fact << " + " << sd_temp << " X " << sd_mul_factor << " = " << thr ;

    //Stores those  pixels whose va;ues are higher than 'thr'.
    for (int i = 0 ; i < xsize * ysize ; i ++)
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

    LOG (INFO) << "SIGMA  Factor::" << sd_mul_factor ;
    LOG (INFO) << " Size of First cut Peaks  " << Fy.size () ;

    if (Fy.size () < thr_intensity_refinement)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG (ERROR) << sd_mul_factor << " less than 0!!!! "  ;
            return (EXIT_FAILURE) ;
        }
        goto label ;
        //        LOG(INFO) << endl << "First cut peaks detected : " << Fy.size () << endl ;
        //        LOG(INFO) << endl << "***No peaks found ***" << endl ;
    }


    /*    for (int i = 0 ; i < xsize * ysize ; i++)
           peakImage[i] = 0 ;
        for (int i = 0 ; i < Fy.size () ; i++)
            peakImage[Fy[i] * xsize + Fx[i]] = Fval[i] ;
     */



    //if winsize is even, make it odd
    if (refine_Winsize % 2 == 0)
        refine_Winsize = refine_Winsize - 1 ;

    LOG (INFO) <<  "Using window size : " << refine_Winsize << " for refining peaks " ;

    //refined peaks
    vector<int> Tx , Ty ;
    vector<float> Tval ;

    Tx = Fx ;
    Ty = Fy ;
    Tval = Fval ;

    Star star1 ;
    star_track.clear () ;
    for (int i = 0 ; i < Tx.size () ; i ++)
    {
        star1.x = Tx[i] ;
        star1.y = Ty[i] ;
        star1.intensity = Tval[i] ;
        star_track.push_back (star1) ;

    }

    sort (star_track.begin () , star_track.end () , compare) ;

    /*refining peaks logic
       refined Window size is for the refined  peaks.
      Refined peaks are found by  making window around each of the star(i.e first cut peaks)  and  finding brightest star among that window.*/
    int start_r , end_r , start_c , end_c ;
    //to be removed 
    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;
    //  bool flag_unique=FALSE;
    //int maxEle;
    //added accordingly to GHOSH suggested.
    int maxEle = MAX_FIRSTCUTPIX_CMPR ;


    if (star_track.size () < MAX_FIRSTCUTPIX_CMPR)
    {
        maxEle = star_track.size () ;
    }
        
        
    //     for(int i=0;i<maxEle;i++)
    //   {
    //       LOG(INFO)<<"FirstCut  Peaks->"<<star_track[i].x<<" "<<star_track[i].y;
    //       
    //   }

    //till this
    //   for (int i = 0 ; i < star_track.size () ; i ++)
    //  LOG(INFO)<<maxEle<<endl;
    for (int i = 0 ; i < maxEle ; i ++)
    {
        start_r = star_track[i].y - refine_Winsize / 2 ;
        end_r = star_track[i].y + refine_Winsize / 2 ;
        start_c = star_track[i].x - refine_Winsize / 2 ;
        end_c = star_track[i].x + refine_Winsize / 2 ;
        if (start_r < 0) start_r = 0 ;
        if (end_r >= ysize) end_r = ysize - 1 ;
        if (start_c < 0) start_c = 0 ;
        if (end_c >= xsize) end_c = xsize - 1 ;
        for (int fcpeak = i + 1 ; fcpeak < star_track.size () ; fcpeak ++)
        {
            if (star_track[fcpeak].x > start_c && star_track[fcpeak].x < end_c && star_track[fcpeak].y > start_r && star_track[fcpeak].y < end_r)
            {

                star_track.erase (star_track.begin () + fcpeak) ;
                fcpeak -- ;
            }
        }
        Tx.push_back (star_track[i].x) ;
        Ty.push_back (star_track[i].y) ;
        Tval.push_back (star_track[i].intensity) ;
        maxEle = star_track.size () ;
    }
    //LOG(INFO)<<Tx.size ()<<" "<<maxEle<<endl;exit(1);
    //   for(int i=0;i<Tx.size ();i++)
    //   {
    //       LOG(INFO)<<"Refined Peaks"<<Tx[i]<<" "<<Ty[i];     
    //   }

    //till this

    //    for (int i = 0 ; i < Fx.size () ; i++)
    //    {
    //        start_r = Ty[i] - refine_Winsize / 2 ;
    //        end_r = Ty[i] + refine_Winsize / 2 ;
    //        start_c = Tx[i] - refine_Winsize / 2 ;
    //        end_c = Tx[i] + refine_Winsize / 2 ;
    //        if (start_r < 0) start_r = 0 ;
    //        if (end_r >= ysize) end_r = ysize - 1 ;
    //        if (start_c < 0) start_c = 0 ;
    //        if (end_c >= xsize) end_c = xsize - 1 ;
    //        int max = 0 ;
    //        for (int k = start_r ; k <= end_r ; k++)
    //        {
    //            for (int l = start_c ; l <= end_c ; l++)
    //            {
    //
    //                if (inputArray[k * xsize + l] > max)
    //                {
    //                    max = inputArray[k * xsize + l] ;
    //                    Tx[i] = l ;
    //                    Ty[i] = k ;
    //                    Tval[i] = inputArray[k * xsize + l] ;
    //                } //  end of if block 
    //            } //end of l loop
    //        } //end of  k  loop
    //    } // end of i loop

    /*--------------Refining peaks completed----------------*/

    float *arr_refine = new float[xsize * ysize] ; //to store refined peaks
    for (int i = 0 ; i < xsize * ysize ; i ++)
        arr_refine[i] = 0.0f ;


    for (int i = 0 ; i < Ty.size () ; i ++)
        arr_refine[Ty[i] * xsize + Tx[i]] = Tval[i] ; //overwriting the same place..

    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;

    for (int i = 0 ; i < xsize * ysize ; i ++)
    {
        if (arr_refine[i] != 0)
        {
            Rx.push_back ((i % xsize)) ;
            Ry.push_back ((i / xsize)) ;
            Rval.push_back (arr_refine[i]) ;
        }
    }

    //Added this 
    //int temp_x,temp_y,temp_int;
    float t = 0 ;
    //temp_int=0;

    //    remove for newer algorithm 
    for (int i = 0 ; i < Rval.size () ; i ++)
    {
        for (int j = Rval.size () - 1 ; j > i ; j --)
        {
            if (Rval[j - 1] < Rval[j])
            {
                swap1 (Rval[j] , Rval[j - 1]) ;
                swap1 (Rval[j] , Rval[j - 1]) ;
                swap1 (Rval[j] , Rval[j - 1]) ;
            }
        }
    }
    
    int maxEle_refined = MAX_REFINED_PIX_SIZE ;


    if (Rx.size () < MAX_REFINED_PIX_SIZE)
    {
        maxEle_refined = Rx.size () ;
    }
    LOG(INFO)<<MAX_REFINED_PIX_SIZE<<" "<<Rx.size ();
    //for(int i=0;i<Rx.size ();i++)
    for (int i = 0 ; i < maxEle_refined ; i ++)
    {
        Tx.push_back (Rx[i]) ;
        Ty.push_back (Ry[i]) ;
        Tval.push_back (Rval[i]) ;
    }

    Rx.clear () ;
    Ry.clear () ;
    Rval.clear () ;
    Rx = Tx ;
    Ry = Ty ;
    Rval = Tval ;


    //till this

    LOG (INFO) << "Number of final peaks is " << Rval.size ()  ;


    if (Ry.size () < maxEle_refined)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG (ERROR) << sd_mul_factor << " less than 0!!!! "  ;
            return (EXIT_FAILURE) ;
        }
        goto label ;
    }
    /**method for  find Centroid within the Stars**/
    
    doCentroiding (Rx , Ry , centroid_Winsize , inputArray , ysize , xsize) ;
    
    
    for (int i = 0 ; i < Ci.size () ; i ++)
    {
        for (int j = Ci.size () - 1 ; j > i ; j --)
        {
            if (Ci[j - 1] < Ci[j])
            {
                swap1 (Ci[j] , Ci[j - 1]) ;
                swap1 (Cx[j] , Cx[j - 1]) ;
                swap1 (Cy[j] , Cy[j - 1]) ;
            }
        }
    }
    
    delete[] arr_refine ;
    return (EXIT_SUCCESS) ;
}       

//method For doing Centroiding


void uvtDetectStar::doCentroiding (vector<int> &X , vector<int> &Y , int centroidwindow , float *arr , int h , int w)
{
    Cx.clear () ;
    Cy.clear () ;
    Ci.clear () ;
    float x_temp , y_temp , int_temp ;
    Star star1 ;
    star_track.clear () ;
    float x , y , val = 0 ;
    float val_Envelope;
    //if centroidwindowsize  is even, make it odd
    if (centroidwindow % 2 == 0)
        centroidwindow = centroidwindow - 1 ;

    LOG (INFO) << "Using  window size : " << centroidwindow << " for finding Centroids " ;
    double sum_x = 0 , sum_y = 0 , sum = 0 ;
    int cnt_envelope=0;
    vector<float> ints_envelopePix;
    float median_envelope=0;
    /**Centroidwindow is  window size for the centroiding.
    Centroid is  done by creating the window around each point of X and Y  vector and finding the centroid by adding the pixel intensity of 
     *  each point of window. **/
    
  
    for (int i = 0 ; i < X.size () ; i ++)
    {
      
        //LOG(INFO)<<X[i]<<" "<<Y[i]<<" "<<""<<arr[(Y[i]) * w + (X[i])] <<endl;
        sum=0;
        if (X[i] != INVALID_PIX_VALUE && Y[i] != INVALID_PIX_VALUE)
        {
            sum_x = 0 ;
            sum_y = 0 , sum = 0 ;
            cnt_envelope=0;
            ints_envelopePix.clear ();
            //median calculation for envelope
              for (int j = - 1 *( (centroidwindow+1) / 2) ; j <=( (centroidwindow+1) / 2) ; j ++)
            {
                for (int k = - 1 * ( (centroidwindow+1) / 2) ; k <= ( (centroidwindow+1) / 2) ; k ++)
                {
                   // x = X[i] ;
                   // y = Y[i] ;
                    if((j<- 1 * (centroidwindow/ 2) || j>1 * (centroidwindow / 2)) || ((k<- 1 * (centroidwindow / 2) || k>1 * (centroidwindow / 2))))
                    {
                       //%#Added ON 20July#% 
                       if((Y[i] + j) >0 && (Y[i] + j)<h &&  (X[i] + k)>0 && (X[i] + k)<w)
//%#-Till this-20July17#%
{
                        if (arr[(Y[i] + j) * w + (X[i] + k)] == INVALID_PIX_VALUE) continue ;
                        cnt_envelope++;
                        ints_envelopePix.push_back ( arr[(Y[i] + j) * w + (X[i] + k) ]);
                        //LOG(INFO)<<X[i]<<" "<<Y[i]<<" "<<arr[(Y[i]) * w + (X[i]) ]<<" "<<X[i] + k<<" "<<Y[i] + j<<" "<<arr[(Y[i] + j) * w + (X[i] + k) ]<<" "<<k<<" "<<j;
}                    
}
                    //if (arr[(Y[i] + j) * w + (X[i] + k)] == INVALID_PIX_VALUE) continue ;
                    //val = arr[(Y[i] + j) * w + (X[i] + k)] ;
                    //val_Envelope=val_Envelope+val;
                    //                else val=0;
                    // LOG(INFO)<<val;
                    
                }
            }
            sort(ints_envelopePix.begin (),ints_envelopePix.end ());
            median_envelope=findMedianValue (ints_envelopePix,ints_envelopePix.size ());
            
            //LOG(INFO)<<cnt_envelope<<" "<<median_envelope;exit(1);
//            LOG(INFO)<<median_envelope;
            //till this
            for (int j = - 1 * (centroidwindow / 2) ; j <= (centroidwindow / 2) ; j ++)
            {
                for (int k = - 1 * (centroidwindow / 2) ; k <= (centroidwindow / 2) ; k ++)
                {
                    x = X[i] ;
                    y = Y[i] ;
//%#Added ON 20July#%
			if((Y[i] + j) >0 && (Y[i] + j)<h &&  (X[i] + k)>0 && (X[i] + k)<w)
//%#-Till this-20July17#%	
	{               
		    if (arr[(Y[i] + j) * w + (X[i] + k)] == INVALID_PIX_VALUE) continue ;
		
                    val = arr[(Y[i] + j) * w + (X[i] + k)] ;
                    //                else val=0;
                    // LOG(INFO)<<val;

                    sum_x = sum_x + (x + k) * (val-median_envelope) ;
			
                    sum_y = sum_y + (y + j) * (val-median_envelope) ;
                    sum = sum + val-median_envelope ;
			//LOG(INFO)<<val<<" "<<X[i]<<" "<<Y[i]<<" "<<X[i] + k<<" "<<Y[i] + j<<" "<<sum<<" "<<median_envelope;

		}
                }
            }

//exit(1);
            if (sum <= 0)
           // if (sum == 0)
            {
                LOG (ERROR) <<  "***Sum of intensities for (" << X[i] << " , " << Y[i] << ")  is <=0 ***" ;
                LOG (ERROR) <<  "***Divide by zero error***" ;
                continue;
            }
            /**Average value of x and y of the  Window is find out and put it in the CX,CY and Ci**/
  //%#Editeded ON 20July-Float casting removal #%         
            x_temp =  sum_x /  sum ;
            y_temp =  sum_y /  sum ;
//%#-Till this-20July17#%
            int_temp =  sum ;
            Cx.push_back (x_temp) ;
            Cy.push_back (y_temp) ;
            Ci.push_back (int_temp) ;
            star1.x = x_temp ;
            star1.y = y_temp ;
            star1.intensity = int_temp ;
            star_track.push_back (star1) ;
        }
        else
        {
            star1.x = INVALID_PIX_VALUE ;
            star1.y = INVALID_PIX_VALUE ;
            star1.intensity = 0 ;
            star_track.push_back (star1) ;
        }

    }
    /**Sorting the list on the basis of the intensity**/
    // sort (star_track.begin (),star_track.end (),compare);
    // LOG(INFO)<<"=============================================="<<endl;
    Cx.clear () , Cy.clear () ;
    Ci.clear () ;

    for (int i = 0 ; i < star_track.size () ; i ++)
    {
        Cx.push_back ((float) star_track[i].x) ;
        Cy.push_back ((float) star_track[i].y) ;
        Ci.push_back (star_track[i].intensity) ;
        //LOG(INFO)<<star_track[i].x<<" "<<star_track[i].y<<" "<<star_track[i].intensity;
    }
    //LOG(INFO)<<"=============================================="<<endl;

    //for (int i = 0 ; i < Ci.size () ; i++)
    //{
    //    
    //        for (int j = Ci.size () - 1 ; j > i ; j--)
    //        {
    //            if (Ci[j - 1] < Ci[j])
    //            {
    //                swap1 (Ci[j] , Ci[j - 1]) ;
    //                swap1 (Cx[j] , Cx[j - 1]) ;
    //                swap1 (Cy[j] , Cy[j - 1]) ;
    //            }
    //        }
    //}
//exit(1);
}


int uvtDetectStar::getStarCentroid (float* inputArray  , int nacc , float bckgrd_fct , float rms_factor , int size_x ,
        int size_y , float winsize , float cent_winsize  , vector<int> &rx , vector<int> &ry ,
        vector<float> &cx , vector<float> &cy , vector<int> &fx , vector<int> &fy , vector<float> &rint ,
        vector<float> &cint , vector<float> &fint , char *mode , float min_No_of_stars , bool flag_check , int Window_nh , int algoflg)
{
    //  LOG(INFO)<<"The flag value "<<flag_check<<endl;
    this->xsize = size_x ;
    this->ysize = size_y ;
    this->refine_Winsize = winsize ;
    this->backgrnd_fact = bckgrd_fct ;
    this->Nacc = nacc ;
    this->centroid_Winsize = cent_winsize ;
    this->sd_mul_factor = rms_factor ;
    this->datainfo.setObsMode (mode) ;
    this->thr_intensity_refinement = min_No_of_stars ;
    this->win_search = Window_nh ;
    //LOG(INFO)<< this->refine_Winsize << " " << this->sd_mul_factor << " " << this->centroid_Winsize << endl ;
    //  this->centroidlimit = cent_limit ;
    this->Cx = cx ;
    this->Cy = cy ;
    this->Rx = rx ;
    this->Ry = ry ;
    this->Ci = cint ;
    this->Rval = rint ;
    this->Fx = fx ;
    this->Fy = fy ;
    this->Fval = fint ;
    int status ;
    if (algoflg == 1)
    {
        if (flag_check == FALSE)//for checking first frame
        {
            status = findStar_algo1 (inputArray ) ;
            if (status)
            {
                LOG (ERROR) << "Error in find Star algorithm" << endl ;
                return (EXIT_FAILURE) ;
            }

            this->cx_ref_vect = this->Cx ;
            this->cy_ref_vect = this->Cy ;
            this->ci_ref_vect = this->Ci ;
            //cout<<"Inside the findpeaks "<<cx_ref_vect.size ()<<endl;

        }
        else
        {
            LOG (INFO) << "Window search criteria for matching pixel ->" << win_search ;
            status = matchStars_TorefFrame (this->cx_ref_vect , this->cy_ref_vect , this->ci_ref_vect , inputArray) ;
            if (status)
            {
                LOG (ERROR) << endl << "***Error in matching stars wrt first(Reference frame)   " << endl ;
                return (EXIT_FAILURE) ;
            }

            this->cx_ref_vect.clear () ;
            this->cy_ref_vect.clear () ;
            this->ci_ref_vect.clear () ;
            this->cx_ref_vect = this->Cx ;
            this->cy_ref_vect = this->Cy , this->ci_ref_vect = this->Ci ;

        }

        cx = this->Cx ;
        cy = this->Cy ;
        rx = this->Rx ;
        ry = this->Ry ;
        cint = this->Ci ;
        rint = this->Rval ;
        fint = this->Fval ;
        fx = this->Fx ;
        fy = this->Fy ;
        //    obj_dataInfo=this->datainfo;
    }
    else if(algoflg==2){
         if (flag_check == FALSE)//for checking first frame
        {
          status = findStar_algo_temp (inputArray ) ;
        if (status)
        {
            LOG (ERROR) << "Error in find Star algorithm 2" << endl ;
            return (EXIT_FAILURE) ;
        }
           this->cx_ref_vect = this->Cx ;
            this->cy_ref_vect = this->Cy ;
            this->ci_ref_vect = this->Ci ;
         }
         else{
             LOG (INFO) << "Window search criteria for matching pixel ->" << win_search ;
            status = matchStars_TorefFrame (this->cx_ref_vect , this->cy_ref_vect , this->ci_ref_vect , inputArray) ;
            if (status)
            {
                LOG (ERROR) << endl << "***Error in matching stars wrt first(Reference frame)   " << endl ;
                return (EXIT_FAILURE) ;
            }

            this->cx_ref_vect.clear () ;
            this->cy_ref_vect.clear () ;
            this->ci_ref_vect.clear () ;
            this->cx_ref_vect = this->Cx ;
            this->cy_ref_vect = this->Cy , this->ci_ref_vect = this->Ci ;
         }
        cx = this->Cx ;
        cy = this->Cy ;
        rx = this->Rx ;
        ry = this->Ry ;
        cint = this->Ci ;
        rint = this->Rval ;
        fint = this->Fval ;
        fx = this->Fx ;
        fy = this->Fy ;
        
    }
    else if (algoflg == 3)
    {
        status = findStar_algo3 (inputArray ) ;
        if (status)
        {
            LOG (ERROR) << "Error in find Star algorithm" << endl ;
            LOG(ERROR)<<"CRASH FAILED FIND STAR (uvtDetectStar.cpp)";
            return (EXIT_FAILURE) ;
        }
        cx = this->Cx ;
        cy = this->Cy ;
        rx = this->Rx ;
        ry = this->Ry ;
        cint = this->Ci ;
        rint = this->Rval ;
        fint = this->Fval ;
        fx = this->Fx ;
        fy = this->Fy ;
    }
    else if (algoflg == 4)
    {
        //      status = findStar_algo3 (inputArray ) ;
        //    if(status)
        //    {
        //        LOG(ERROR)<<"Error in find Star algorithm"<<endl;
        //        return(EXIT_FAILURE);
        //    }
        //    cx = this->Cx ;
        //    cy = this->Cy ;
        //    rx = this->Rx ;
        //    ry = this->Ry ;
        //    cint = this->Ci ;
        //    rint = this->Rval ;
        //    fint = this->Fval ;
        //    fx = this->Fx ;
        //    fy = this->Fy ;
        vector<float> cint_background , diff_intensity , final_cx , final_cy , final_cint , diff_int_x , diff_int_y ;
        vector <float> cx_background , cy_background ;
        spMatrix A (512 * 512 , 6) ;
        spMatrix B (512 * 512 , 1) ;
        spMatrix X (6 , 1) ;
        //   LOG(INFO)<<"ERR";
        for (int i = 44 ; i < xsize - 44 ; i ++)
        {
            for (int j = 44 ; j < ysize - 44 ; j ++)
            {
                if (inputArray[i * xsize + j] != INVALID_PIX_VALUE)
                {
                    A (i , 0) = 1 ;
                    A (i , 1) = i ;
                    A (i , 2) = j ;
                    A (i , 3) = i*i ;
                    A (i , 4) = j*j ;
                    A (i , 5) = i*j ;
                    B (i , 0) = inputArray[i * xsize + j] ;
                }
            }
        }
        X.ApplyLeastSquare (A , B) ;
        for (int i = 44 ; i < xsize - 44 ; i ++)
        {
            for (int j = 44 ; j < ysize - 44 ; j ++)
            {
                if (inputArray[i * xsize + j] != INVALID_PIX_VALUE)
                {
                    cint_background.push_back (X (0 , 0) + X (1 , 0) * i + X (2 , 0) * j + X (3 , 0) * i * i + X (4 , 0) * j * j + X (5 , 0) * i * j) ;
                    cx_background.push_back (i) ;
                    cy_background.push_back (j) ;
                }
                else
                {
                    cint_background.push_back (INVALID_PIX_VALUE) ;
                    cx_background.push_back (INVALID_PIX_VALUE) ;
                    cy_background.push_back (INVALID_PIX_VALUE) ;
                }
            }
        }
        for (int i = 44 ; i < xsize - 44 ; i ++)
        {
            for (int j = 44 ; j < ysize - 44 ; j ++)
            {
                if (inputArray[i * xsize + j] != INVALID_PIX_VALUE  && cint_background[i * xsize + j] != INVALID_PIX_VALUE  )
                {
                    diff_intensity.push_back (inputArray[i * IMG_DIM_DI + j] - cint_background[i * IMG_DIM_DI + j]) ;
                    diff_int_x.push_back (i) ;
                    diff_int_y.push_back (j) ;
                }
                else
                {
                    diff_intensity.push_back (INVALID_PIX_VALUE) ;
                    diff_int_x.push_back (INVALID_PIX_VALUE) ;
                    diff_int_y.push_back (INVALID_PIX_VALUE) ;
                }
            }
        }
        //    for(int i=0;i<cx.size ();i++)
        //    {    
        //    diff_intensity.push_back (cint[i]-cint_background[i]);    
        //    }
        double sd_diff_centroid_pix = getSD (diff_intensity.data () , diff_intensity.size ()) ;
        //  cx.clear ();cy.clear ();cint.clear ();
        for (int i = 44 ; i < xsize - 44 ; i ++)
        {
            for (int j = 44 ; j < ysize - 44 ; j ++)
            {
                if (inputArray[i * xsize + j] > cint_background[i * xsize + j] + 2 * sd_diff_centroid_pix)
                {
                    final_cx.push_back (i) ;
                    final_cy.push_back (j) ;
                    final_cint.push_back (inputArray[i * xsize + j]) ;
                }
            }
        }



        //    diff_intensity.clear ();cint_background.clear ();
        //    cx=final_cx;
        //    cy=final_cy;
        //    cint=final_cint;

        cx = cx_background ;
        cy = cy_background ;
        cint = cint_background ;

        final_cx.clear () ;
        final_cy.clear () ;
        final_cint.clear () ;
    }


    return (EXIT_SUCCESS) ;
}


int uvtDetectStar::findStar_algo2 (float *inputArray)
{
    /**findstar_algo2-this algorithm create a window around a each pixel of  first cut peaks based on the square size.
     *  minimum of 4 corner values and center pixel value of each window is identified .
     *calculate summation of each pixel value of window and also calculate  (minimum corner value*window size*window size)[considering the value if we put minimum corner value in each pixel of window]
     * calculation of subtraction of (centerpix value - minimum corner value )
     *perform (sum of all the pixel values of window)  - (minimum corner value*window size*window size).
     * this parameters are compared with primary  and secondary threshold and stars and centroids are identified..
     */
    Rx.clear () ;
    Ry.clear () ;
    Rval.clear () ;
    int sr_no = 0 ;
    double center_pix_val ; //center pixel value of algo.square
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
    int s1 , s2 ;
    s1 = s2 = algo_Square_Size ;
    int sqr_x[s1][s2] , sqr_y[s1][s2] ;
    int sum_x , sum_y ;
    int temp_m , temp_n , minus_val ;
    double **sqr_pix ;
    float *data ;
    int count ;
    int xy_index = 0 ;
    int m , n ;
    sqr_pix = allocateMemory<double> (s1 , s2) ;
    double sum_xh , sum_yh ;
    double **Xc , **Yc ;
    int **X , **Y ;
    double **Center_PixelValue ;
    Center_PixelValue = allocateMemory<double>(xsize , ysize) ;
    Xc = allocateMemory<double>(xsize , ysize) ;
    Yc = allocateMemory<double>(xsize , ysize) ;
    X = allocateMemory<int>(xsize , ysize) ;
    Y = allocateMemory<int>(xsize , ysize) ;
    minus_val = (algo_Square_Size % 2 == 0) ? (algo_Square_Size / 2) : ((algo_Square_Size / 2) + 1) ;
    data = new float[xsize * ysize] ;
    int temp_index = 0 ;
    for (int i = 0 ; i < xsize * ysize ; i ++)
    {
        if (inputArray[i] > 0.0)
        {
            data[i] = inputArray[i] ;
            temp_index ++ ;
        }
    }
    long int p ;
    double max_pix = data[0] ;
    double min_pix = data[0] ;
    for (p = 1 ; p < temp_index ; p ++)
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
    LOG (INFO) << "Using Primary threshold " << primary_threshold_Val << " & secondary  threshold  " << secondary_threshold_Val << endl ;

    for (int i = 0 ; i < xsize ; i ++)
    {
        count = 0 ;
        for (int j = 0 ; j < ysize ; j ++)
        {
            sum_Of_SqrPixVal = 0 ;
            if (i + s1 < xsize && j + s2 < ysize)
            {
                for (m = i , c1 = 0 ; m < i + s1 ; m ++ , c1 ++)
                {
                    for (n = j , c2 = 0 ; n < j + s2 ; n ++ , c2 ++)
                    {
                        //**read required pixel values and x,y locations
                        sqr_pix[c1][c2] = inputArray[m * xsize + n] ;
                        sqr_x[c1][c2] = n ;
                        sqr_y[c1][c2] = m ;
                        //**make total of all pixel within algorithm shape (square of 3*3 or 5*5 or 24*24) 
                        sum_Of_SqrPixVal = sum_Of_SqrPixVal + sqr_pix[c1][c2] ;
                    }
                }
                //read square corner value & define minimum of them as "fcmin"
                x = n - 1 ;
                y = m - 1 ;

                temp_m = x - minus_val ;
                temp_n = y - minus_val ;

                a = inputArray[j * xsize + i] ;
                b = inputArray[j * xsize + x] ;
                c = inputArray[y * xsize + i] ;
                d = inputArray[y * xsize + x] ;

                temp = (a <= b) ? a : b ;
                temp1 = (temp <= c) ? temp : c ;
                fcmin = (temp1 <= d) ? temp1 : d ;

                center_pix_val = inputArray[temp_n * xsize + temp_m] ; //*finding center pix value//

                check1 = (center_pix_val - fcmin) ;
                check2 = sum_Of_SqrPixVal - (fcmin * (s1 * s2)) ; //*for step7 

                if (center_pix_val != 0)
                {
                    if (fcmin != 0)
                    {
                        flag = 0 ;
                        if (check1 > primary_threshold_Val)
                        {

                            flag = check_MaxCenterPixVal (center_pix_val , sqr_pix , algo_Square_Size) ; //check center pixel value is local maximum

                            if (flag == 1)
                            {
                                if (check2 > secondary_threshold_Val)
                                {
                                    //centroid of square
                                    sum_xh = 0.0 ;
                                    sum_x = 0 ;
                                    sum_yh = 0.0 ;
                                    sum_y = 0 ;
                                    for (int ii = 0 ; ii < s1 ; ii ++)
                                    {
                                        for (int jj = 0 ; jj < s2 ; jj ++)
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
                                    Rx.push_back (temp_m + 1) ;
                                    Ry.push_back (temp_n + 1) ;
                                    Rval.push_back (center_pix_val) ;
                                    sr_no ++ ;
                                    xy_index ++ ;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 0 ; i < xsize ; i ++)
        delete[] Xc[i] , Yc[i] , X[i] , Y[i] ;
    doCentroiding (Rx , Ry , centroid_Winsize , inputArray , ysize , xsize) ;
    return (EXIT_SUCCESS) ;
}


//algorithm 3 


int uvtDetectStar::findStar_algo3 (float *inputArray ) //algorithm for finding the peaks
{
    float mean_Ofimage = 0.0 ;
    float mean_Ofimage_Final=0.0f;
    
    if (xsize == 0 || ysize == 0)
    {
        LOG (ERROR) << "***Divide by Zero***"  ;
        return (EXIT_FAILURE) ;
    }
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
    float *temp_array ;
    vector<float> array_temp ;

    if (datainfo.getModeFlag () == PC)
    {
        array_temp.clear () ;
        for (int i = 0 ; i < xsize * ysize ; i ++)
        {

            if (inputArray[i] != 0.0f)
            {
                array_temp.push_back (inputArray[i]) ;

            }
        }
        temp_array  = new float[array_temp.size ()] ;
        for (int in = 0 ; in < array_temp.size () ; in ++)
        {
            temp_array[in] = array_temp[in] ;

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
        LOG (ERROR) << "***SD_MULTI_FACTOR is <0***"  ;
        return (EXIT_FAILURE) ;
    }
    double thr = 0 ;
    double thr_second=0;
    double sd_temp = 0.0f ;
    double sd_temp_Final=0.0f;
    if (datainfo.getModeFlag () == PC)
    {
        sd_temp = getSD (temp_array , array_temp.size ()) ;
        mean_Ofimage = getmean (temp_array , array_temp.size ()) ;
        //   thr =  sd_temp* sd_mul_factor ;

    }
    else
    {
        sd_temp = getSD (inputArray , xsize * ysize) ;
        mean_Ofimage = getmean (inputArray , xsize * ysize) ;
        //thr =  sd_temp* sd_mul_factor ;
    }

if (datainfo.getModeFlag () == PC)
    {
    delete[] temp_array;
    }
   
    thr = mean_Ofimage + sd_temp* 3.5;
    thr_second = mean_Ofimage - sd_temp* 3.5 ;
    
     // thr = mean_Ofimage + sd_temp* sd_mul_factor;
    //thr_second = mean_Ofimage - sd_temp* sd_mul_factor ;
  
    LOG (INFO) << "Threshold for first cut peaks is   " << mean_Ofimage << " + " << sd_temp << " X 3.5" <<  " = " << thr ;

    //Stores those  pixels whose va;ues are higher than 'thr'.
    
//    for (int i = 0 ; i < xsize * ysize ; i ++)
//    {
//        r = (i / xsize) ;
//        c = (i % xsize) ;
//
//        if (inputArray[i] >thr )
//        {
//            Fval.push_back (inputArray[i]) ;
//            Fx.push_back (c) ; //x is for column
//            Fy.push_back (r) ; //y is for row
//        }
//    }
//    
    
    for (int i = 0 ; i < xsize * ysize ; i ++)
    {
        r = (i / xsize) ;
        c = (i % xsize) ;

        if (inputArray[i] <thr &&  inputArray[i]>thr_second  && inputArray[i]!=INVALID_PIX_VALUE )
        {
            Fval.push_back (inputArray[i]) ;
            Fx.push_back (c) ; //x is for column
            Fy.push_back (r) ; //y is for row
        }
    }
    if(Fx.size()==0){
        LOG(WARNING)<<"***0 star found***";
	LOG(ERROR)<<"CRASH:NO PIXELS OUTSIDE +-3.5 SIGMA (uvtDetectStar.cpp)"; 
        return(EXIT_SUCCESS);
        
    }
    
    
    
    //till this
//    Fx=FtempX;
//    Fy=FtempY;
//    Fval=FtempInt;
    

    LOG (INFO) << "SIGMA  Factor::" << sd_mul_factor ;
    LOG (INFO) << " Size of First cut Peaks  " << Fy.size () ;

    
    //Now removing further
    sd_temp_Final = getSD (Fval.data() , Fval.size()) ;
    mean_Ofimage_Final = getmean (Fval.data() , Fval.size()) ;
    
    vector<int> Fx_final,Fy_final;
    vector<float> Fval_final;
    double thr_final= mean_Ofimage_Final+sd_mul_factor*sd_temp_Final;
     LOG (INFO) << "Final threshold for first cut peaks is   " << mean_Ofimage_Final << " + " << sd_temp_Final << " X " <<sd_mul_factor << " = " << thr_final ;
    for (int i=0;i< xsize*ysize ;i++)
    {
        r = (i / xsize) ;
        c = (i % xsize) ;
        if(inputArray[i]>thr_final)
        {
            Fval_final.push_back (inputArray[i]) ;
            Fx_final.push_back (c) ; //x is for column
            Fy_final.push_back (r) ; //y is for row
        }
        
    }
    Fx.clear();Fy.clear();Fval.clear();
    Fx=Fx_final;
    Fy=Fy_final;
    Fval=Fval_final;
  //  LOG(INFO)<<Fx.size();
    //remove those stars which is a stripe .
    vector<int> FtempX,FtempY;
    vector<float> FtempInt;
   // FtempX=Fx;
   // FtempY=Fy;
   // FtempInt=Fval;
    int count_NeighbourPix=0;
    int count_ValidPix=0;
   
    float val_plusY2,val_plusY1,val_plusX2,val_plusX1,val_minY2,val_minY1,val_minX2,val_minX1;
    //LOG(INFO)<<"RRR";
    float curr_Pix_value=0.0f;
    vector<int> index_ToBedelete;
    int count_lesserIntensity_Pix;
    float value_ToCompare=NEIGHBOUR_PIX_MULTFACTOR*thr_final+mean_Ofimage_Final;
    if(datainfo.getModeFlag()==IM){
    for (int i=0;i<Fval.size();i++){
//        LOG(INFO)<<i<<" "<<Fx[i]<<" "<<Fy[i];
        count_lesserIntensity_Pix=0;
        curr_Pix_value=inputArray[Fy[i]*xsize+Fx[i]];
        
        count_NeighbourPix=0;
        count_ValidPix=0;
       
       for (int j=Fx[i]-2;j<=Fx[i]+2;j++) 
       {
           for (int k=Fy[i]-2;k<=Fy[i]+2;k++)
           {
               if((k*xsize+j)>0 && (k*xsize+j)<xsize*ysize )
               {
                   count_NeighbourPix++;
                   if(inputArray[k*xsize+j]!=INVALID_PIX_VALUE){
                       count_ValidPix++;
                   }
                  
                   if(curr_Pix_value>inputArray[k*xsize+j]){
                       count_lesserIntensity_Pix++;
                   }
                   
               }               
           }
           
       }
       // LOG(INFO)<<i<<" "<<count_lesserIntensity_Pix<<" "<<count_NeighbourPix;
        if(count_lesserIntensity_Pix!=(count_NeighbourPix-1)){//not maximum within 5X5 window
             Fx[i]=INVALID_PIX_VALUE;
             Fy[i]=INVALID_PIX_VALUE;
             Fval[i]=INVALID_PIX_VALUE;
     //        break;
             continue;
        }
       
       if(count_NeighbourPix!=count_ValidPix)//Found one of the pixel as a INVALID PIXEL in neighborhood
       {
          // inputArray[Fy[i]*xsize+Fx[i]]=INVALID_PIX_VALUE;
           Fx[i]=INVALID_PIX_VALUE;
           Fy[i]=INVALID_PIX_VALUE;
           Fval[i]=INVALID_PIX_VALUE;
           //break;
          continue;
          
       }
      // else//checking for star like thing
      // {
       
        if(((Fy[i]-2)*xsize+(Fx[i]-2))<0 || ((Fy[i]+2)*xsize+(Fx[i]+2))>xsize*ysize){
            Fx[i]=INVALID_PIX_VALUE;
           Fy[i]=INVALID_PIX_VALUE;
           Fval[i]=INVALID_PIX_VALUE;
           //break;
          continue;
        }
              
           val_minY1=inputArray[(Fy[i]-1)*xsize+Fx[i]];
           val_minY2=inputArray[(Fy[i]-2)*xsize+Fx[i]]    ;    
           val_plusY1=inputArray[(Fy[i]+1)*xsize+Fx[i]];
           val_plusY2=inputArray[(Fy[i]+2)*xsize+Fx[i]];
           val_minX1=inputArray[(Fy[i])*xsize+(Fx[i]-1)];
           val_minX2=inputArray[(Fy[i])*xsize+(Fx[i]-2)];
           val_plusX1=inputArray[(Fy[i])*xsize+(Fx[i]+1)];
           val_plusX2=inputArray[(Fy[i])*xsize+(Fx[i]+2)];
         //  LOG(INFO)<<"PPP";
          
           
           if(!((val_minY1 >= val_minY2 && val_plusY1>=val_plusY2 && val_minX1 >= val_minX2 && val_plusX1>=val_plusX2) && (val_minY1 >value_ToCompare &&
               val_plusY1>value_ToCompare  && val_minX1>value_ToCompare && val_plusX1>value_ToCompare)))
         
            // if(((val_minY1 > val_minY2 && val_plusY1>val_plusY2 && val_minX1 > val_minX2 && val_plusX1>val_plusX2) ))
           {
           //inputArray[Fy[i]*xsize+Fx[i]]=INVALID_PIX_VALUE;
            Fx[i]=INVALID_PIX_VALUE;         
           Fy[i]=INVALID_PIX_VALUE;          
           Fval[i]=INVALID_PIX_VALUE;
          
           }
          
              
          
     //  }      
        
        //check for brightest one withing 5X5
        
        
    }    
    
   
    for (int i = 0; i < Fx.size(); i++) {
        if (Fx[i] != INVALID_PIX_VALUE) {
            FtempX.push_back(Fx[i]);
            FtempY.push_back(Fy[i]);
            FtempInt.push_back(Fval[i]);
           
        }

    }
   // LOG(INFO)<<"TEP->"<<FtempX.size();
    Fx=FtempX;
    Fy=FtempY;
    Fval=FtempInt;
        
    Rx=Fx;
    Ry=Fy;
    Rval=Fval;
    LOG(INFO)<<"Total final star found->"<<Rval.size();
    }
//    exit(1);
    Fval_final.clear();Fx_final.clear();Fy_final.clear();
    
//    if (Fy.size () < 2)
//    {
//        sd_mul_factor = sd_mul_factor - 0.25 ;
//        if (sd_mul_factor <= 0)
//        {
//            LOG (ERROR) << sd_mul_factor << " less than 0!!!! "  ;
//            return (EXIT_FAILURE) ;
//        }
//        goto label ;
//        
//    }


    /*    for (int i = 0 ; i < xsize * ysize ; i++)
           peakImage[i] = 0 ;
        for (int i = 0 ; i < Fy.size () ; i++)
            peakImage[Fy[i] * xsize + Fx[i]] = Fval[i] ;
     */

    //commnted as per need

    //if winsize is even, make it odd
    if(datainfo.getModeFlag()==PC)
    {
    if (refine_Winsize % 2 == 0)
        refine_Winsize = refine_Winsize - 1 ;

    LOG (INFO) <<  "Using window size : " << refine_Winsize << " for refining peaks " ;

    //refined peaks
    vector<int> Tx , Ty ;
    vector<float> Tval ;

    Tx = Fx ;
    Ty = Fy ;
    Tval = Fval ;

    //LOG(INFO)<<"SIZE->"<<Tval.size();
    Star star1 ;
    star_track.clear () ;
    for (int i = 0 ; i < Tx.size () ; i ++)
    {
        star1.x = Tx[i] ;
        star1.y = Ty[i] ;
        star1.intensity = Tval[i] ;
        star_track.push_back (star1) ;

    }

    sort (star_track.begin () , star_track.end () , compare) ;

for (int i=0;i<star_track.size();i++){
for (int j=star_track.size()-1;j>i;j--){

if(star_track[i].intensity ==star_track[j].intensity)
{
if(star_track[i].x>star_track[j].x){
   swap1(star_track[i].intensity,star_track[j].intensity);
swap1(star_track[i].x,star_track[j].x);
swap1(star_track[i].y,star_track[j].y);


}
}
}
}

    /*refining peaks logic
       refined Window size is for the refined  peaks.
      Refined peaks are found by  making window around each of the star(i.e first cut peaks)  and  finding brightest star among that window.*/
    int start_r , end_r , start_c , end_c ;
    //to be removed 
    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;
 

    for (int i = 0 ; i < star_track.size () ; i ++)       
    {
        start_r = star_track[i].y - refine_Winsize / 2 ;
        end_r = star_track[i].y + refine_Winsize / 2 ;
        start_c = star_track[i].x - refine_Winsize / 2 ;
        end_c = star_track[i].x + refine_Winsize / 2 ;
        if (start_r < 0) start_r = 0 ;
        if (end_r >= ysize) end_r = ysize - 1 ;
        if (start_c < 0) start_c = 0 ;
        if (end_c >= xsize) end_c = xsize - 1 ;
        for (int fcpeak = i + 1 ; fcpeak < star_track.size () ; fcpeak ++)
        {
            if (star_track[fcpeak].x > start_c && star_track[fcpeak].x < end_c && star_track[fcpeak].y > start_r && star_track[fcpeak].y < end_r)
            {

                star_track.erase (star_track.begin () + fcpeak) ;
                fcpeak -- ;
            }
        }
        Tx.push_back (star_track[i].x) ;
        Ty.push_back (star_track[i].y) ;
        Tval.push_back (star_track[i].intensity) ;
        // maxEle=star_track.size ();
    }
    

    /*--------------Refining peaks completed----------------*/

    float *arr_refine = new float[xsize * ysize] ; //to store refined peaks
    for (int i = 0 ; i < xsize * ysize ; i ++)
        arr_refine[i] = 0.0f ;


    for (int i = 0 ; i < Ty.size () ; i ++)
        arr_refine[Ty[i] * xsize + Tx[i]] = Tval[i] ; //overwriting the same place..

    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;

    for (int i = 0 ; i < xsize * ysize ; i ++)
    {
        if (arr_refine[i] != 0)
        {
            Rx.push_back ((i % xsize)) ;
            Ry.push_back ((i / xsize)) ;
            Rval.push_back (arr_refine[i]) ;
        }
    }


    LOG (INFO) << "Number of final peaks is " << Rval.size ()  ;

     //Added logic
    double sum_refinePeaks=0.0f;
    vector<float> track_SumIntensity;
    vector<int> Rx_final,Ry_final;
    vector<float> Rval_final;
  // LOG(INFO)<<Rval.size();
    
    for (int i=0;i<Rval.size();i++)
    { 
//          if(Rx[i]<0){
//                LOG(INFO)<<"ERROR "<<Rx[i];exit(1);
//            }
        start_r = Ry[i] - 1 ;
        end_r = Ry[i] + 1 ;
        start_c = Rx[i] - 1 ;
        end_c = Rx[i] + 1 ;
        sum_refinePeaks=0.0f;
        for (int indexx=start_c;indexx<end_c;indexx++)
        {
            for (int indexy=start_r;indexy<end_r;indexy++)
            {
                if(inputArray[indexy*xsize+indexx]!=INVALID_PIX_VALUE)
                sum_refinePeaks=sum_refinePeaks+inputArray[indexy*xsize+indexx];
                
            }
            
        }
        track_SumIntensity.push_back(sum_refinePeaks);     
   }
    
  
    for(int i=0;i<track_SumIntensity.size();i++)
    {
        if(track_SumIntensity[i]> thr_intensity_refinement*sd_temp_Final)
        {
//            if(Rx[i]<0){
//                LOG(INFO)<<"ERROR "<<Rx[i];exit(1);
//            }
            Rx_final.push_back(Rx[i]);
            Ry_final.push_back(Ry[i]);
            Rval_final.push_back(Rval[i]);           
        }
        
    }
   
    Rx.clear();Ry.clear();Rval.clear();
    Rx=Rx_final;
    Ry=Ry_final;
    Rval=Rval_final;
    Rx_final.clear(); Ry_final.clear(); Rval_final.clear();
}
//   
    
    
    
//    if (Ry.size () < 2)
//    {
//        sd_mul_factor = sd_mul_factor - 0.25 ;
//        if (sd_mul_factor <= 0)
//        {
//            LOG (ERROR) << sd_mul_factor << " less than 0!!!! "  ;
//            return (EXIT_FAILURE) ;
//        }
//        delete[] arr_refine;
//        goto label ;
//    }
    /**method for  find Centroid within the Stars**/
   

   doCentroiding (Rx , Ry , centroid_Winsize , inputArray , ysize , xsize) ;
//  Cx.clear ();Cy.clear ();Ci.clear ();
//    for (int i=0;i<Rx.size ();i++)
//    {
//       Cx.push_back (Rx[i]);
//       Cy.push_back (Ry[i]);
//       Ci.push_back (Rval[i]);
//    }
   // LOG(INFO)<<"DDD111";exit(1);
    for (int i = 0 ; i < Ci.size () ; i ++)
    {
      
        for (int j = Ci.size () - 1 ; j > i ; j --)
        {
            if (Ci[j - 1] < Ci[j])
            {
                swap1 (Ci[j] , Ci[j - 1]) ;
                swap1 (Cx[j] , Cx[j - 1]) ;
                swap1 (Cy[j] , Cy[j - 1]) ;
            }
        }
    }
    //delete[] arr_refine ;
    return (EXIT_SUCCESS) ;
}


int uvtDetectStar::findStar_algo4 (float *inputArray ) //algorithm for finding the peaks
{
    float mean_Ofimage = 0.0 ;
    if (xsize == 0 || ysize == 0)
    {
        LOG (ERROR) << "***Divide by Zero***"  ;
        return (EXIT_FAILURE) ;
    }
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
    float *temp_array ;
    vector<float> array_temp ;

    if (datainfo.getModeFlag () == PC)
    {
        array_temp.clear () ;
        for (int i = 0 ; i < xsize * ysize ; i ++)
        {

            if (inputArray[i] != 0.0f)
            {
                array_temp.push_back (inputArray[i]) ;

            }
        }
        temp_array  = new float[array_temp.size ()] ;
        for (int in = 0 ; in < array_temp.size () ; in ++)
        {
            temp_array[in] = array_temp[in] ;

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
        LOG (ERROR) << "***SD_MULTI_FACTOR is <0***"  ;
        return (EXIT_FAILURE) ;
    }
    double thr = 0 ;
    double sd_temp = 0.0f ;

    if (datainfo.getModeFlag () == PC)
    {
        sd_temp = getSD (temp_array , array_temp.size ()) ;
        mean_Ofimage = getmean (temp_array , array_temp.size ()) ;
        //   thr =  sd_temp* sd_mul_factor ;

    }
    else
    {

        sd_temp = getSD (inputArray , xsize * ysize) ;
        mean_Ofimage = getmean (inputArray , xsize * ysize) ;
        //thr =  sd_temp* sd_mul_factor ;
    }


    //  sd_temp=getSD (inputArray , xsize * ysize);
    //  mean_Ofimage=getmean (inputArray,xsize*ysize);
    //    thr =  sd_temp* sd_mul_factor ;
    thr = mean_Ofimage + sd_temp* sd_mul_factor ;
    LOG (INFO) << "Threshold for first cut peaks is   " << mean_Ofimage << " + " << sd_temp << " X " << sd_mul_factor << " = " << thr ;

    //Stores those  pixels whose va;ues are higher than 'thr'.
    for (int i = 0 ; i < xsize * ysize ; i ++)
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

    LOG (INFO) << "SIGMA  Factor::" << sd_mul_factor ;
    LOG (INFO) << " Size of First cut Peaks  " << Fy.size () ;

    if (Fy.size () < thr_intensity_refinement)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG (ERROR) << sd_mul_factor << " less than 0!!!! "  ;
            return (EXIT_FAILURE) ;
        }
        goto label ;
        //        LOG(INFO) << endl << "First cut peaks detected : " << Fy.size () << endl ;
        //        LOG(INFO) << endl << "***No peaks found ***" << endl ;
    }


    /*    for (int i = 0 ; i < xsize * ysize ; i++)
           peakImage[i] = 0 ;
        for (int i = 0 ; i < Fy.size () ; i++)
            peakImage[Fy[i] * xsize + Fx[i]] = Fval[i] ;
     */



    //if winsize is even, make it odd
    if (refine_Winsize % 2 == 0)
        refine_Winsize = refine_Winsize - 1 ;

    LOG (INFO) <<  "Using window size : " << refine_Winsize << " for refining peaks " ;

    //refined peaks
    vector<int> Tx , Ty ;
    vector<float> Tval ;

    Tx = Fx ;
    Ty = Fy ;
    Tval = Fval ;

    Star star1 ;
    star_track.clear () ;
    for (int i = 0 ; i < Tx.size () ; i ++)
    {
        star1.x = Tx[i] ;
        star1.y = Ty[i] ;
        star1.intensity = Tval[i] ;
        star_track.push_back (star1) ;

    }

    sort (star_track.begin () , star_track.end () , compare) ;

    /*refining peaks logic
       refined Window size is for the refined  peaks.
      Refined peaks are found by  making window around each of the star(i.e first cut peaks)  and  finding brightest star among that window.*/
    int start_r , end_r , start_c , end_c ;
    //to be removed 
    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;
   
    for (int i = 0 ; i < star_track.size () ; i ++)
        //  LOG(INFO)<<maxEle<<endl;
        //for (int i = 0 ; i <maxEle ; i ++)
    {
        start_r = star_track[i].y - refine_Winsize / 2 ;
        end_r = star_track[i].y + refine_Winsize / 2 ;
        start_c = star_track[i].x - refine_Winsize / 2 ;
        end_c = star_track[i].x + refine_Winsize / 2 ;
        if (start_r < 0) start_r = 0 ;
        if (end_r >= ysize) end_r = ysize - 1 ;
        if (start_c < 0) start_c = 0 ;
        if (end_c >= xsize) end_c = xsize - 1 ;
        for (int fcpeak = i + 1 ; fcpeak < star_track.size () ; fcpeak ++)
        {
            if (star_track[fcpeak].x > start_c && star_track[fcpeak].x < end_c && star_track[fcpeak].y > start_r && star_track[fcpeak].y < end_r)
            {

                star_track.erase (star_track.begin () + fcpeak) ;
                fcpeak -- ;
            }
        }
        Tx.push_back (star_track[i].x) ;
        Ty.push_back (star_track[i].y) ;
        Tval.push_back (star_track[i].intensity) ;
        // maxEle=star_track.size ();
    }
  

    /*--------------Refining peaks completed----------------*/

    float *arr_refine = new float[xsize * ysize] ; //to store refined peaks
    for (int i = 0 ; i < xsize * ysize ; i ++)
        arr_refine[i] = 0.0f ;


    for (int i = 0 ; i < Ty.size () ; i ++)
        arr_refine[Ty[i] * xsize + Tx[i]] = Tval[i] ; //overwriting the same place..

    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;

    for (int i = 0 ; i < xsize * ysize ; i ++)
    {
        if (arr_refine[i] != 0)
        {
            Rx.push_back ((i % xsize)) ;
            Ry.push_back ((i / xsize)) ;
            Rval.push_back (arr_refine[i]) ;
        }
    }


    LOG (INFO) << "Number of final peaks is " << Rval.size ()  ;
    if (Ry.size () < thr_intensity_refinement)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG (ERROR) << sd_mul_factor << " less than 0!!!! "  ;
            return (EXIT_FAILURE) ;
        }
        goto label ;
    }
    /**method for  find Centroid within the Stars**/
    doCentroiding (Rx , Ry , centroid_Winsize , inputArray , ysize , xsize) ;
    for (int i = 0 ; i < Ci.size () ; i ++)
    {
        for (int j = Ci.size () - 1 ; j > i ; j --)
        {
            if (Ci[j - 1] < Ci[j])
            {
                swap1 (Ci[j] , Ci[j - 1]) ;
                swap1 (Cx[j] , Cx[j - 1]) ;
                swap1 (Cy[j] , Cy[j - 1]) ;
            }
        }
    }
    delete[] arr_refine ;
    return (EXIT_SUCCESS) ;
}


int uvtDetectStar::findStar_algo_temp (float *inputArray ) //algorithm for finding the peaks
{
//    float mean_Ofimage = 0.0 ;
//    if (xsize == 0 || ysize == 0)
//    {
//        LOG (ERROR) << "***Divide by Zero***"  ;
//        return (EXIT_FAILURE) ;
//    }
//    Fx.clear () ;
//    Fy.clear () ;
//    Fval.clear () ;
//    Rx.clear () ;
//    Ry.clear () ;
//    Rval.clear () ;
//    Cx.clear () ;
//    Cy.clear () ;
//    Ci.clear () ;
//    int r , c ;
//    float *temp_array ;
//    vector<float> array_temp ;
//
//    if (datainfo.getModeFlag () == PC)
//    {
//        array_temp.clear () ;
//        for (int i = 0 ; i < xsize * ysize ; i ++)
//        {
//
//            if (inputArray[i] != 0.0f)
//            {
//                array_temp.push_back (inputArray[i]) ;
//
//            }
//        }
//        temp_array  = new float[array_temp.size ()] ;
//        for (int in = 0 ; in < array_temp.size () ; in ++)
//        {
//            temp_array[in] = array_temp[in] ;
//
//        }
//    }
//
//
//
//label:
//    Fval.clear () ;
//    Fx.clear () ;
//    Fy.clear () ;
//    Rx.clear () ;
//    Ry.clear () ;
//    Rval.clear () ;
//
//    if (sd_mul_factor < 0)
//    {
//        LOG (ERROR) << "***SD_MULTI_FACTOR is <0***"  ;
//        return (EXIT_FAILURE) ;
//    }
//    double thr = 0 ;
//    double sd_temp = 0.0f ;
//
//    if (datainfo.getModeFlag () == PC)
//    {
//        sd_temp = getSD (temp_array , array_temp.size ()) ;
//        mean_Ofimage = getmean (temp_array , array_temp.size ()) ;
//        //   thr =  sd_temp* sd_mul_factor ;
//
//    }
//    else
//    {
//
//        sd_temp = getSD (inputArray , xsize * ysize) ;
//        mean_Ofimage = getmean (inputArray , xsize * ysize) ;
//        //thr =  sd_temp* sd_mul_factor ;
//    }

    vector<int > X_pix_OFwindow,Y_pix_OFwindow;
    vector<int> X_pix_OFwindow_final,Y_pix_OFwindow_final;
    vector<float> Int_pix_OFwindow_final,Int_pix_OFwindow;
    double  sd_window,mean_window;
    float X_win,Y_win,int_Win;
    int cnt_invalid=0;
    for (int i = 0 ; i < xsize ; i = i + 15)
    {     
        for (int j = 0 ; j < ysize ; j = j + 15)
        {
            cnt_invalid=0;
            for (int p = i ; p < i + 15 ; p ++)
            {
                for (int q = j ; q < j + 15 ; q ++)
                {
                    
                     if(inputArray[q*ysize+p]==INVALID_PIX_VALUE)
                    {
                         cnt_invalid++;
                    }
                    Int_pix_OFwindow.push_back (inputArray[q * ysize + p]) ;
                    X_pix_OFwindow.push_back (p);
                    Y_pix_OFwindow.push_back (q);

                }
           
            }           
            if (cnt_invalid==15*15)
            {
                mean_window=0;
                sd_window=0;
            }
            else
            {
          mean_window=getmean (Int_pix_OFwindow.data (),Int_pix_OFwindow.size ());
          sd_window=getSD (Int_pix_OFwindow.data (),Int_pix_OFwindow.size ());
            }
      //    LOG(INFO)<<mean_window<<" "<<sd_window<<" "<<Int_pix_OFwindow.size ();
          X_win=X_pix_OFwindow[0];
          Y_win=Y_pix_OFwindow[0];
          int_Win=Int_pix_OFwindow[0];
         int cnt=0;
         for (int p=0;p<Int_pix_OFwindow.size ();p++)
         {
          //   if(Int_pix_OFwindow[p]>int_Win)
             if(Int_pix_OFwindow[p]>mean_window+sd_window)
             {
                  int_Win=Int_pix_OFwindow[cnt];
                       X_pix_OFwindow_final.push_back (X_pix_OFwindow[p]);
                      Y_pix_OFwindow_final.push_back (Y_pix_OFwindow[p]);
                      Int_pix_OFwindow_final.push_back (Int_pix_OFwindow[p]);
                 
             }             
         }
//          for (int p = i ; p < i + 15 ; p ++)
//            {
//              for (int q = j ; q < j + 15 ; q ++)
//               {
//                  if(Int_pix_OFwindow[cnt]>mean_window+sd_window)
//                 //  if(Int_pix_OFwindow[cnt]>int_Win)
//                   {
//                       int_Win=Int_pix_OFwindow[cnt];
//                      X_pix_OFwindow_final.push_back (p);
//                      Y_pix_OFwindow_final.push_back (q);
//                      Int_pix_OFwindow_final.push_back (Int_pix_OFwindow[cnt]);
//                   }   
//                   cnt++;
//               }               
//           }
          Int_pix_OFwindow.clear ();
          X_pix_OFwindow.clear ();
          Y_pix_OFwindow.clear ();

        } 
    
    }
    double mean_final=getmean (Int_pix_OFwindow_final.data (),Int_pix_OFwindow_final.size ());
    double sd_final=getSD (Int_pix_OFwindow_final.data (),Int_pix_OFwindow_final.size ());
     //double mean_final=getmean (inputArray,xsize*ysize);
  //  double sd_final=getSD (inputArray,xsize*ysize);
    Fx.clear (),Fy.clear (),Fval.clear ();Rx.clear (),Ry.clear (),Rval.clear ();
    Fx=X_pix_OFwindow_final;
    Fy=Y_pix_OFwindow_final;
    Fval=Int_pix_OFwindow_final;
    Cx.clear ();Cy.clear ();Ci.clear ();
    for (int i=0;i<Int_pix_OFwindow_final.size ();i++)
    {
       if(Int_pix_OFwindow_final[i]>mean_final+sd_final)
       {
        Rx.push_back (X_pix_OFwindow_final[i]);
        Ry.push_back (Y_pix_OFwindow_final[i]);
        Rval.push_back (Int_pix_OFwindow_final[i]);
        }    
    }
     for (int i = 0 ; i < Rval.size () ; i ++)
    {
        for (int j = Rval.size () - 1 ; j > i ; j --)
        {
            if (Rval[j - 1] < Rval[j])
            {
                swap1 (Rval[j] , Rval[j - 1]) ;
                swap1 (Rval[j] , Rval[j - 1]) ;
                swap1 (Rval[j] , Rval[j - 1]) ;
            }
        }
    }
    
    int maxEle_refined = MAX_REFINED_PIX_SIZE ;
   if (Rx.size () < MAX_REFINED_PIX_SIZE)
    {
        maxEle_refined = Rx.size () ;
    }
    Rx.resize (maxEle_refined);
    Ry.resize (maxEle_refined);
    Rval.resize (maxEle_refined);
    doCentroiding (Rx,Ry,centroid_Winsize,inputArray,xsize,ysize);
    for (int i = 0 ; i < Ci.size () ; i ++)
    {
        for (int j = Ci.size () - 1 ; j > i ; j --)
        {
            if (Ci[j - 1] < Ci[j])
            {
                swap1 (Ci[j] , Ci[j - 1]) ;
                swap1 (Cx[j] , Cx[j - 1]) ;
                swap1 (Cy[j] , Cy[j - 1]) ;
            }
        }
    }
  
    LOG(INFO)<<Rx.size ()<<" "<<Cx.size ();
            
    
    return (EXIT_SUCCESS) ;
}       
float findMedianValue (vector<float> &input , int noOfValues)
{
    //cout<<"INSIDE";
    float temp ;
    int index , j ;
    sort (input.begin () , input.end ()) ;
    // Sort the input data
    //        for(index=0;index<noOfValues;index++)
    //                for(j=index+1;j<noOfValues;j++)
    //                {
    //                        if(input[index]>input[j])
    //                        {
    //                                temp=input[j];
    //                                input[j]=input[index];
    //                                input[index]=temp;
    //                        }
    //                }

    //    If the no of values are even number, take the average of middle and middle-1 value.
    if (noOfValues % 2 == 0)
        return (input[noOfValues / 2] + input[(noOfValues / 2) - 1]) / 2.0 ;
    else
        return input[noOfValues / 2] ;
}