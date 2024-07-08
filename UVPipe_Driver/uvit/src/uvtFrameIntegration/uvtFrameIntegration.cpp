/* 
 * File:   uvtFrameIntegration.cpp
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
#include "uvtFrameIntegration.h"
#include<uvtUtils.h>
#include<vector>
#include<algorithm>
#include<iterator>
#include<glog/logging.h>
//Constructor -called when object is created
//bool compare (struct Star vect1 , struct Star vect2) ;
//bool compare (struct Star vect1 , struct Star vect2)
//{
//    return (vect1.intensity > vect2.intensity) ;
//}

uvtFrameIntegration::uvtFrameIntegration ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}

//Destructor
uvtFrameIntegration::~uvtFrameIntegration () {
 }

int uvtFrameIntegration::read (int argc , char** argv)
{
    int status = 0 ;

    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG(ERROR) << "***Error Initializing PIL***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("inputdatadir" , inputdatadir)))
    {
        LOG(ERROR) << endl << "***Error reading input data directory***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("frames2discard" , &Ndiscard)))
    {
        LOG(ERROR) << endl << "***Error reading number of frames to discard at beginning***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("frames2integrate" , &Nacc)))
    {
        LOG(ERROR) << endl << "***Error reading number of frames to integrate ***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetFname ("outdir" , outdir)))
    {
        LOG(ERROR) << endl << "***Error reading output directory path***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("clobber" , &clobber)))
    {
        LOG(ERROR) << "***Error reading clobber:" << clobber << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetBool ("history" , &history)))
    {
        LOG(ERROR) << "***Error reading history parameter:" << history << "***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetString ("mode" , mode)))
    {
        LOG(ERROR) << "***Error reading mode parameter:" << history << "***" ;
        return status ;
    }
    PILClose (status) ;
    return (EXIT_SUCCESS) ;
}

int uvtFrameIntegration::read (char *inputdata_dir , int Ndiscard , int Nacc , char* out_dir , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdata_dir) ;
    strcpy (this->outdir , out_dir) ;
    this->Ndiscard = Ndiscard ;
    this->Nacc = Nacc ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_FAILURE) ;
}

void uvtFrameIntegration::display ()
{
     LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT FRAME INTEGRATION PARAMETERS              " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input Directory                                      : " << inputdatadir ;
    LOG(INFO) << endl << "Output Directory                               : " << outdir ;
    LOG(INFO) << endl << "Number of frames to discard      : " << Ndiscard ;
    LOG(INFO) << endl << "Number of frames to accumulate      : " << Nacc ;
    if (clobber == YES)
        LOG(INFO) << endl << "Overwrite                                        : YES" ;
    else
        LOG(INFO) << endl << "Overwrite                                        : NO" ;
    if (history == YES)
        LOG(INFO) << endl << "History                                             : YES" ;
    else
        LOG(INFO) << endl << "History                                             : NO" ;
    LOG(INFO) << endl << "----------------------------------------------------------------------------------------------------------\n" ;

}

int uvtFrameIntegration::uvtFrameIntProcess ()
{
    LOG(INFO) << endl << "Starting FrameIntegration process................." << endl ;
    sprintf (moduleoutdir , "%s/%s" , outdir , modulename) ;

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

    nframes = 0 ; //Number of frames
    /**Shell command for creating the output Directory**/
    cmd = "mkdir -p " + (string) moduleoutdir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ; // creating output directory 

    LOG(INFO) << endl << "Created output directory  " << moduleoutdir ;

    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s" , moduleoutdir , "SignalFrames") ;
    cmd = "mkdir -p " + (string) dir ;
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << "Created signal frame directory " << dir ;

    sprintf (dir , "%s/%s" , moduleoutdir , "ExposureFrames") ;
    cmd = "mkdir -p " + (string) dir ;
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << "Created exposure frame directory " << dir ;

    //opening info file in input directory to get data information
    string tempfilepath = searchFile (inputdatadir , "info") ;
    if (tempfilepath == " ")
    {
        LOG(ERROR) << endl << "Error in finding info file" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;

    LOG(INFO) << endl << "\nInput information file :" << infofile_in ;

    long int tot_num_frames = 0 ;
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "***Error in opening the input information file***") ;

    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "***Error in moving to 2nd HDU in input information file***") ;

    datainfo.getInfo(finfo_in) ; //reading basic information for data from information file
    if (datainfo.getModeFlag () != PC)
    {
        LOG(ERROR) << "This process is not applicable to modes other than Photon Counting" ;
        return (EXIT_FAILURE) ;
    }

    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    if (xsize <= 0 || ysize <= 0)
    {
        LOG(ERROR) << endl << "***Invalid xsize/ysize***\n" ;
        return (EXIT_FAILURE) ;
    }
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "***Error in reading the key value of the NAMEPRFX***") ;
    fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
    printError (status , "***Error reading the key value of the EVTFILE***") ;
    fits_close_file(finfo_in,&status);
    sprintf (infofile_out , "%s/%s_fi.info" , moduleoutdir , nameprefix) ;
    LOG(INFO) << endl << "\nOutput information file :"<< infofile_out ;

    //creating output information file
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the file" , infofile_out) ; //for creating name for output information file
    char *ttype[] = {"ExposureFileList" , "ImageFilelist"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating tht table" , infofile_out) ;
 
    datainfo.setXsize (IMG_DIM_FI) ;
    datainfo.setYsize (IMG_DIM_FI) ;
    datainfo.write (finfo_out) ;

    //writeCommonKeywords (finfo_out , modulename) ;
    fits_close_file (finfo_out , &status) ;

    /*----info file creating completed, rest of the information will be put by other functions-----------*/

    vector<string> imglist , explist ;

    char exp_filename[NAMESIZE] , img_filename[NAMESIZE] ;
    char file_Name[NAMESIZE] ;
    char file_in[NAMESIZE]; //for input and output event file
    sprintf (file_in , "%s/%s" , inputdatadir , eventfile) ; //taking event file full path
    char outfile3[FLEN_FILENAME] ;
    LOG(INFO)<<"\nInput Event file: "<<file_in<<endl;
    fitsfile *fevt_in ;
   fits_open_file (&fevt_in , file_in , READWRITE , &status) ;
    printError (status , "Error in opening the input event file",file_in) ;

     long *frame_no ;
    double *t , *ENP ;
    float *xFrac , *yFrac ;
    unsigned short *bad_Flag,*mult_phn;
    long nrows ;
    int tfields3 = 4 ;
    char *ttype3[] = {"Average time" , "X" , "Y" , "Intensity"} ;
    char *tform3[] = {"D" , "I" , "I" , "D"} ;
    sprintf (outfile3 , "%s/%s_fi.time" , moduleoutdir , nameprefix) ;
//    fitsfile *finptr ;
//    fits_create_file (&finptr , outfile3 , &status) ;
//    printError (status , "Error creating the output File " , outfile3) ;
//    fits_create_tbl (finptr , BINARY_TBL , 0 , tfields3 , ttype3 , tform3 , NULL , "" , &status) ;
//    printError (status , "Error creating the table for the FirstCut pixels " , outfile3) ;
//    int rown = 1 ;
    copyUsrkeywrdsTovect (fevt_in,key_records);      
    fits_movabs_hdu (fevt_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU",file_in) ;
    fits_get_num_rows (fevt_in , &nrows , &status) ;
     printError (status , "Error in reading the number of rows",file_in) ;
    
     frame_no = new long[nrows] ;
    ENP = new double[nrows] ;
    bad_Flag = new unsigned short[nrows] ;
    mult_phn = new unsigned short[nrows] ;
    t = new double[nrows] , xFrac = new float[nrows] , yFrac = new float[nrows] ;
    fits_read_col (fevt_in , TLONG , 2 , 1 , 1 , nrows , NULL , frame_no , NULL , &status) ;
    printError (status , "Error in reading the column of frame_no",file_in) ;
    fits_read_col (fevt_in , TDOUBLE , 3 , 1 , 1 , nrows , NULL , t , NULL , &status) ;
    printError (status , "Error in reading the column of t",file_in) ;
    fits_read_col (fevt_in , TFLOAT , 4 , 1 , 1 , nrows , NULL , xFrac , NULL , &status) ;
    printError (status , "Error in reading the column of Xfrac",file_in) ;
    fits_read_col (fevt_in , TFLOAT , 5 , 1 , 1 , nrows , NULL , yFrac , NULL , &status) ;
    printError (status , "Error in reading the  column of  yFrac",file_in) ;
    fits_read_col (fevt_in , TDOUBLE , 10 , 1 , 1 , nrows , NULL , ENP , NULL , &status) ;
    printError (status , "Error in reading the column of Exposure",file_in) ;
    fits_read_col (fevt_in , TUSHORT , 8 , 1 , 1 , nrows , NULL , bad_Flag , NULL , &status) ;
    printError (status , "Error in reading the column of Exposure",file_in) ;
    fits_read_col (fevt_in , TUSHORT , 9 , 1 , 1 , nrows , NULL , mult_phn, NULL , &status) ;
    printError (status , "Error in reading the column of Exposure",file_in) ;

    int jj = 0 ;
    int x_tem , y_tem ;
    LOG(INFO) << "\nTotal number of Events :" << nrows << endl ;
    double Start_Time , End_Time , Avg_Time = 0.0 ;
    float **Image_Array ;
    float **Exposure_Img ;
    Image_Array = new float*[IMG_DIM_FI] ;
    Exposure_Img = new float*[IMG_DIM_FI] ;
    float *one_dim_exp ;
    one_dim_exp = new float[IMG_DIM_FI * IMG_DIM_FI] ;
   // unsigned short *one_dim_exp_temp ;
    //one_dim_exp_temp = new unsigned short[IMG_DIM_FI * IMG_DIM_FI] ;
    float *one_dim_img ;
    fitsfile *fptr ;
    one_dim_img = new float[IMG_DIM_FI * IMG_DIM_FI] ;
    for (int i = 0 ; i < IMG_DIM_FI ; i++)
    {
        Image_Array[i] = new float[IMG_DIM_FI] ;
        Exposure_Img[i] = new float[IMG_DIM_FI] ;
        for (int j = 0 ; j < IMG_DIM_FI ; j++)
        {
            Image_Array[i][j] = 0.0f ;
            Exposure_Img[i][j] =  0.0f ;
        }
    }
    long firstpix1[2] ;
    firstpix1[0] = firstpix1[1] = 1 ;
    int naxis = 2 ;
    long naxes[2] = {IMG_DIM_FI , IMG_DIM_FI} ;
    float Eff_NoPhtn = 0.0f ;
    float multi_factor = IMG_DIM_FI / xsize ;
    int file_num = 1 ;
    int frame_check ;
    LOG(INFO) << "\nMULTIPLICATION FACTOR  for converting image to 9600 is " << multi_factor << endl ;
    int frame_init = frame_check = Nacc ;
    double time_frame = 0.0 , time_frame_final = 0.0 ;
    int cnt_frame = 0 ;
    vector<string> vhistorystr ;
    if(history==YES)
    getHistory (vhistorystr) ;
    /****provision for if frames are not starting from 1****/
    int start_row = 0 ;
    int cmpr_term = 0 ;
    if (frame_no[0] != 1)
    {
        cmpr_term = Ndiscard + (frame_no[0] - 1) ;
    }
    for (int j = 0 ; j < nrows ; j++)
    {
        if (frame_no[j] > cmpr_term)
        {
            start_row = j ;
            break ;
        }
    }
    
  
    LOG(INFO)<<"\nPerforming Frame Integration...."<<endl;
    //int cnt_juzcheck=0;
   // cout<<start_row<<endl;exit(1);
//    bool flag_greater=FALSE;
//    if(frame_no[0]>Nacc)
//    {
//        flag_greater=TRUE;
//    }
    if(Nacc>frame_no[nrows-1]){
               LOG(ERROR)<<"***Number of frames to accumulate is greater than total available frames";
                return(EXIT_FAILURE);
    }
 
    //in case of integrate frame number is less than the starting frame number.
    while(frame_no[0]>Nacc)
    {
            Nacc=Nacc+frame_init;            
    }
    vector<float> x_track_exp,y_track_exp,Eff_noPhtn_track_exp;
      vector<float> x_track_sig,y_track_sig,Eff_noPhtn_track_sig;
    int tfields = 3 ;
    char *ttype1[] = {"X" , "Y" , "Intensity"} ;
    char *tform1[] = {"E" , "E" , "E"} ;
   //  int cntr=0;
//    Star star1;   
 //star_track_exp.clear ();
     if (frame_no[0] != 1)
    {
        cmpr_term = Ndiscard + (frame_no[0] - 1) ;
    }
     for (int j = 0 ; j < nrows ; j ++)
    {
        if (frame_no[j] > cmpr_term)
        {
            start_row = j ;
            break ;
        }
    }
      Nacc= Nacc+frame_no[0]-1;
    for (int k = start_row ; k < nrows ; k++)
    {
        /**in 'IF condition' accumulated array of  Nacc frames will be  generate   & In 'ELSE condition' array will be written to the frame.
         x and y location  for number of events are taken as a indexes for the array .For each Nacc events one array will be written as a frame. 
         **/
//        if(frame_no[k]>Nacc)
//        {
//            flag_greater=FALSE;
//        }
        
       
        if (frame_no[k] <= (Nacc + Ndiscard))
        {
//            if (Nacc + Ndiscard > frame_no[nrows - 1])
//            {
//                frame_init = frame_no[nrows - 1 ] % frame_check ;
//            }
            while (jj == 0)
            {
                Start_Time = t[k] ;
                jj++ ;
            }        
          
            x_tem = (int) xFrac[k] ;
            y_tem = (int) yFrac[k] ;
           
            x_tem = x_tem * multi_factor ;
            y_tem = y_tem * multi_factor ;
           
           
            if (round(y_tem) < IMG_DIM_FI && round(x_tem) < IMG_DIM_FI)
            {
              
                float mult=(float)mult_phn[k];              
                float badflg=(float)bad_Flag[k];
             
                Exposure_Img[y_tem][x_tem] = (Exposure_Img[y_tem][x_tem]+1.0f)*mult*badflg ;
               
                Eff_NoPhtn = (float) ENP[k] ;
                
                Image_Array[y_tem][x_tem] = (Image_Array[y_tem][x_tem] + Eff_NoPhtn)*mult*badflg;            
              
            }
         
            if (frame_no[k - 1] == frame_no[k])
            {
                time_frame = t[k] ;
            }
            else
            {
                time_frame_final = time_frame_final + time_frame ;
            }
//            if (k == nrows - 1)
//            {
//                goto label_else ;
//top
           
                    
        }
        else //frame will be written here.
        { 
          //  x_track_exp.clear ();y_track_exp.clear ();Eff_noPhtn_track_exp.clear ();
          //  x_track_sig.clear ();y_track_sig.clear ();Eff_noPhtn_track_sig.clear ();
          //  star_track_exp.clear ();star_track_sig.clear ();
             jj = 0 ;
            End_Time = t[k - 1] ;

            Avg_Time = Start_Time + (End_Time - Start_Time) / 2 ;
            //cout<<Avg_Time<<endl;
            //sprintf (file_Name , "%s/Frame_Pixels_%d.bin" , moduleoutdir , file_num) ;

            //to be removed
            //cout<<Image_Array[pix_y][pix_x]<<endl
//            cntr=0;
//            for(int pix_x=0;pix_x<IMG_DIM_FI;pix_x++)
//            {
//                for(int pix_y=0;pix_y<IMG_DIM_FI;pix_y++)
//                {
//                    one_dim_exp[cntr]=Exposure_Img[pix_x][pix_y];
//                    one_dim_img[cntr]=Image_Array[pix_x][pix_y];
//                    cntr++;
//                    if(Exposure_Img[pix_y][pix_x]!=0.0f)
//                    {
//                        x_track_sig.push_back (pix_x);
//                        y_track_sig.push_back (pix_y);
//                        Eff_noPhtn_track_sig.push_back (Exposure_Img[pix_y][pix_x]);
////                        star1.x=pix_x;
////                        star1.y=pix_y;
////                        star1.intensity=Exposure_Img[pix_y][pix_x];
////                        star_track_exp.push_back (star1);
//
//                    }
//                     if(Image_Array[pix_y][pix_x]!=0.0f)
//                    {
//                         x_track_exp.push_back (pix_x);
//                        y_track_exp.push_back (pix_y);
//                        Eff_noPhtn_track_exp.push_back (Image_Array[pix_y][pix_x]);
//                         
////                        star1.x=pix_x;
////                        star1.y=pix_y;
////                        star1.intensity=Image_Array[pix_y][pix_x];
////                        star_track_sig.push_back (star1);
//                    }
//                    
//                    
//                }
//                
//            }
//            doCentroiding (x_track_sig,y_track_sig,Eff_noPhtn_track_sig,48,one_dim_img,IMG_DIM_FI,IMG_DIM_FI);
         //   doCentroiding (x_track_exp, y_track_exp,Eff_noPhtn_track_exp,48,one_dim_exp,IMG_DIM_FI,IMG_DIM_FI);
            
//            sort (star_track_sig.begin (),star_track_sig.end (),compare);
//            sort (star_track_exp.begin (),star_track_exp.end (),compare);
//            //till this
//            
//            for(int i=0;i<star_track_sig.size ();i++)
//            {
//                x_track_sig.push_back (star_track_sig[i].x);
//                y_track_sig.push_back (star_track_sig[i].y);
//                Eff_noPhtn_track_sig.push_back (star_track_sig[i].intensity);
//                
//            }
//             for(int i=0;i<star_track_exp.size ();i++)
//            {
//                x_track_exp.push_back (star_track_exp[i].x);
//                y_track_exp.push_back (star_track_exp[i].y);
//                Eff_noPhtn_track_exp.push_back (star_track_exp[i].intensity);
//                
//            }
            long temp = 0 ;
            for (int i = 0 ; i < IMG_DIM_FI ; i++)
            {
                for (int j = 0 ; j < IMG_DIM_FI ; j++)
                {
                    one_dim_img[temp] = 0.0f ;
                    one_dim_exp[temp] = 0.0f ;
                    one_dim_img[temp] = (float) Image_Array[i][j] ;
                    one_dim_exp[temp] =  Exposure_Img[i][j] ;
                    temp++ ;
                }
            }
            // exit(1);
           // ofstream fout (file_Name , ios::binary) ;
           // fout.write ((char*) one_dim_img , (IMG_DIM_FI * IMG_DIM_FI) * sizeof (float)) ;
           // fout.close () ;
        sprintf (exp_filename , "%s/%s/%s_t%f_f%d_%s" , moduleoutdir , "ExposureFrames" , nameprefix , file_num ,Avg_Time, "Centroid_exp.fits") ;
//        fits_create_file (&fptr , exp_filename , &status) ;
//        printError (status , "Error in creating the Exposure File" , exp_filename) ;
//        fits_create_tbl (fptr , BINARY_TBL , 0 , tfields , ttype1 , tform1 , NULL , "Exposures" , &status) ;
//        printError (status , "Error creating the table for the FirstCut pixels " , exp_filename) ;
//        fits_write_col (fptr , TFLOAT , 1 , 1 , 1 , x_track_exp.size () , (void *) x_track_exp.data () , &status) ;
//        printError (status , "Error Writing the Refined pixels's X-cordinates" , exp_filename) ;
//        fits_write_col (fptr , TFLOAT , 2 , 1 , 1 , y_track_exp.size () , (void *) y_track_exp.data () , &status) ;
//        printError (status , "Error Writing the Refined pixels Y-cordinates" , exp_filename) ;
//        fits_write_col (fptr , TFLOAT , 3 , 1 , 1 , Eff_noPhtn_track_exp.size () , (void *)Eff_noPhtn_track_exp.data () , &status) ;
//        printError (status , "Error Writing the Refined  pixel's Intensity" , exp_filename) ;
//        fits_write_key (fptr , TINT , "FRAMENO" , &file_num , "Frame number" , &status) ;
//        printError (status , "Error writing the key value  of the FRAMENO" , exp_filename) ;
//        fits_write_key (fptr , TDOUBLE , "FRMTIME" , &Avg_Time , "Frame time" , &status) ;
//        printError (status , "Error writing the key value of the FRMTIME" , exp_filename) ;
//        fits_close_file (fptr , &status) ;
//        printError (status , "Error closing the File" , exp_filename) ;
            fits_create_file (&fptr , exp_filename , &status) ;
            printError (status , "Error in creating the Exposure File" , exp_filename) ;
            fits_create_img (fptr ,FLOAT_IMG , naxis , naxes , &status) ;
            printError (status , "Error creating the Image for Exposure File" , exp_filename) ;
            fits_write_pix (fptr , TFLOAT , firstpix1 , IMG_DIM_FI*IMG_DIM_FI , one_dim_exp , &status) ;
            printError (status , "Error Writing the pixels " , exp_filename) ;
            fits_write_key (fptr , TINT , "FRAMENO" , &file_num , "Frame number" , &status) ;
            printError (status , "Error writing the key value  of the FRAMENO" , exp_filename) ;
            fits_write_key (fptr , TDOUBLE , "FRMTIME" , &Avg_Time , "Frame time" , &status) ;
            printError (status , "Error writing the key value of the FRMTIME" , exp_filename) ;
            writeUsrkeywordsFrmvect (exp_filename,key_records);
            if (history == YES) writeHistory (exp_filename , vhistorystr) ;
            writeCommonKeywords (fptr,modulename);
                
            fits_close_file (fptr , &status) ;
            printError (status , "Error closing the File" , exp_filename) ;
   
            updateKeywords (exp_filename,modulename);
            explist.push_back (basename (exp_filename)) ;
            sprintf (img_filename , "%s/%s/%s_t%f_f%d_%s" , moduleoutdir , "SignalFrames" , nameprefix , file_num ,Avg_Time, "Centroid_img.fits") ;
//           fits_create_file (&fptr , img_filename , &status) ;
//           printError (status , "Error in creating the Exposure File" , img_filename) ;
//           fits_create_tbl (fptr , BINARY_TBL , 0 , tfields , ttype1 , tform1 , NULL , "Centroid" , &status) ;
//           printError (status , "Error creating the table for the FirstCut pixels " , img_filename) ;
//           fits_write_col (fptr , TFLOAT , 1 , 1 , 1 , x_track_sig.size () , (void *) x_track_sig.data () , &status) ;
//            printError (status , "Error Writing the Refined pixels's X-cordinates" , img_filename) ;
//        fits_write_col (fptr , TFLOAT , 2 , 1 , 1 , y_track_sig.size () , (void *) y_track_sig.data () , &status) ;
//        printError (status , "Error Writing the Refined pixels Y-cordinates" , img_filename) ;
//        fits_write_col (fptr , TFLOAT , 3 , 1 , 1 , Eff_noPhtn_track_sig.size () , (void *)Eff_noPhtn_track_sig.data () , &status) ;
//        printError (status , "Error Writing the Refined  pixel's Intensity" , img_filename) ;
//         fits_write_key (fptr , TINT , "FRAMENO" , &file_num , "Frame number" , &status) ;
//            printError (status , "Error writing the key value of the FRAMENO" , img_filename) ;
//            fits_write_key (fptr , TDOUBLE , "FRMTIME" , &Avg_Time , "Frame time" , &status) ;
//            printError (status , "Error writing the key value of the FRMTIME" , img_filename) ;
//         fits_close_file (fptr , &status) ;
//           printError (status , "Error closing the File" , exp_filename) ;
            fits_create_file (&fptr , img_filename , &status) ;
            printError (status , "Error creating the image File" , outfile3) ;
            fits_create_img (fptr , FLOAT_IMG , naxis , naxes , &status) ;
            printError (status , "Error creating the image in IMAGE  File" , outfile3) ;
            fits_write_pix (fptr , TFLOAT , firstpix1 , IMG_DIM_FI*IMG_DIM_FI , one_dim_img , &status) ;
            printError (status , "Error writing the pixels to the output Signal file" , outfile3) ;
            fits_write_key (fptr , TINT , "FRAMENO" , &file_num , "Frame number" , &status) ;
            printError (status , "Error writing the key value of the FRAMENO" , img_filename) ;
            fits_write_key (fptr , TDOUBLE , "FRMTIME" , &Avg_Time , "Frame time" , &status) ;
            printError (status , "Error writing the key value of the FRMTIME" , img_filename) ;
              writeUsrkeywordsFrmvect (img_filename,key_records);
            if (history == YES) writeHistory (img_filename , vhistorystr) ;
              writeCommonKeywords (fptr,modulename);
            fits_close_file (fptr , &status) ;
            printError (status , "Error closing the Signal File" , outfile3) ;
          
            //updateKeywords (img_filename,modulename);
            imglist.push_back (basename (img_filename)) ;
            file_num++ ;
//           
            for (int i = 0 ; i < IMG_DIM_FI ; i++)
            {
                for (int j = 0 ; j < IMG_DIM_FI ; j++)
                {
                    Image_Array[i][j] = (float)0.0f ;
                    Exposure_Img[i][j] = (float)0.0f ;
                }
            }
            Nacc = Nacc + frame_init ;
            tot_num_frames++ ;
            //LOG(INFO) << "Centroid written for  " << tot_num_frames << endl ;
            cnt_frame = 0 ;
            time_frame = 0.0 ;
            time_frame_final = 0.0 ;
            if (k != nrows - 1)
            {
                k-- ;
            }
            cout<< "Total Files written=" << file_num-1<< "     Remaining files= " <<((int)frame_no[nrows-1]/frame_check)-(file_num-1)<<" \r" ;
        }
        
    
    }
    delete[] one_dim_exp , one_dim_img ;
    delete[] Image_Array ;
    delete[] Exposure_Img ;
    delete[] xFrac ;
    delete[] yFrac ;
    delete[] t ;
    delete[] frame_no ;

   /***Writes the  common information to  output information File***/
   
     LOG(INFO)<<"\nUpdating and writing the list of output frame names to output information file.. "<<endl;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the out information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU" , infofile_out) ;
    fits_update_key (finfo_out , TINT , "NFILES" , &tot_num_frames , "File name prefix" , &status) ;
     printError (status , "***Error in updating the key value of the NFILES***") ;
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "***Error in updating the key value of the NAMEPRFX***") ;
    fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , imglist.size () , imglist.data () , &status) ;
    printError (status , "Error in writing the  output Signal framelist to the  out information file" , infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 2 , 1 , 1 , explist.size () , explist.data () , &status) ;
    printError (status , "Error in writing the  output Exposure framelist to the out information file " , infofile_out) ;
    fits_write_key (finfo_out , TSTRING , "SIGDIR" , (void*) "SignalFrames" , "Directory Name" , &status) ;
    fits_write_key (finfo_out , TSTRING , "EXPDIR" , (void*) "ExposureFrames" , "Directory Name" , &status) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in output information File" , infofile_out) ;
     
    writeUsrkeywordsFrmvect (infofile_out,key_records);          
    if (history == YES) writeHistory (infofile_out , vhistorystr) ;
    return (EXIT_SUCCESS) ;
}

int uvtFrameIntegration::getHistory (vector<string> &vhistory)
{
    char *user = getlogin () ;
    // string cent_alg_str;
    int cnt=0;
    char frm_compute[100] ;
    char frm_discard[100] ;
    sprintf (frm_compute , "%d" , Nacc) ;
    sprintf (frm_discard , "%d" , Ndiscard) ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Number of frames to be discarded" + (string) frm_discard) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Number of frames to compute" + (string) frm_compute) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir) ;
    if (clobber == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+" clobber=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+" clobber=no") ;
    if (history == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+"  history=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+" history=no") ;
    vhistory.push_back ("Parameter List END") ;

    return (EXIT_SUCCESS) ;

}
void uvtFrameIntegration::doCentroiding (vector<float> &X , vector<float> &Y , vector<float> Intensity,int centroidwindow , float *arr , int h , int w)
{
    //Cx.clear () ;
   // Cy.clear () ;
   // Ci.clear () ;
    vector<float> Cx,Cy,Ci;
    float x , y , val = 0 ;
     //if centroidwindowsize  is even, make it odd
    if (centroidwindow % 2 == 0)
        centroidwindow = centroidwindow - 1 ;

    // LOG(INFO) << endl << "Using  window size : " <<centroidwindow<< " for finding Centroids " ;
    double sum_x = 0 , sum_y = 0 , sum = 0 ;
    /**Centroidwindow is  window size for the centroiding.
    Centroid is  done by creating the window around each point of X and Y  vector and finding the centroid by adding the pixel intensity of 
     *  each point of window. **/
    
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
                val = arr[(int)(Y[i] + j) * w + (int)(X[i] + k)] ;
                sum_x = sum_x + (x + k) * val ;
                sum_y = sum_y + (y + j) * val ;
                sum = sum + val ;
            }
        }
        //if (sum <= 0)
        if(sum==0)
        {
            LOG(ERROR) << endl << "Sum of intensities for (" << X[i] << " , " << Y[i] << ")  is <=0" << endl ;
            LOG(ERROR) << endl << "\nDivide by zero error\n" ;
            exit (EXIT_FAILURE) ;
        }
        /**Average value of x and y of the  Window is find out and put it in the CX,CY and Ci**/
        Cx.push_back ((float) sum_x / (float) sum) ;
        Cy.push_back ((float) sum_y / (float) sum) ;
        Ci.push_back ((float) sum) ;
    }
/**Sorting the list on the basis of the intensity**/
    for (int i = 0 ; i < Ci.size () ; i++)
    {

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
    X=Cx;
    Y=Cy;
    Intensity=Ci;
    Cx.clear (),Cy.clear ();Ci.clear ();
}
