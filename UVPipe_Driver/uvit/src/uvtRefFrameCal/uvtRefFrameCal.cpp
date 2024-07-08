/* 
 * File: uvtRefFrameCal.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#include<stdlib.h>
//#include<pil.h>
//#include <fitsio.h>
#include<iostream>
#include<unistd.h>

#include<dirent.h>  //Accessing Directory
#include<string.h> 
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include<glog/logging.h>
#include <uvtRefFrameCal.h>

#define  MODULENAME "uvtRefFrameCal"
#define  NBHD_RADIUS 5

uvtRefFrameCal::uvtRefFrameCal ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    strcpy (centroidDir , "Centroid") ;
}

uvtRefFrameCal::~uvtRefFrameCal () { 

}

int uvtRefFrameCal::read (int argc , char** argv)
{
    int status = 0 ;

    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG(INFO) << "***Error Initializing PIL***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("inputdatadir" , inputdatadir)))
    { //change to inputdir
        LOG(INFO) << endl << "***Error reading framelist file name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("outdir" , outdir)))
    {
        LOG(INFO) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetInt ("framesToBeDiscard" , &framesToDiscard)))
    {
        LOG(INFO) << endl << "***Error reading Number of frames to discard***" ;
        return status ;
    }

    if (PIL_OK != (status = PILGetInt ("averageFactor" , &averageFactor)))
    {
        LOG(INFO) << endl << "***Error reading the average Factor***" ;
        return status ;
    }
  if (averageFactor == 0)
    {
        LOG(ERROR) << "***invalid Input for average Factor,divide by zero error***" << endl ;
        return (EXIT_FAILURE) ;
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

int uvtRefFrameCal::read (char* inputdatadir , char* outdir , int frames_toDiscard , int avg_Factor , int clobber , int history)
{

    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->outdir , outdir) ;

    this->framesToDiscard = frames_toDiscard ;
    this->averageFactor = avg_Factor ;

    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}

void uvtRefFrameCal::display ()
{
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT REF FRAME CAL PARAMETERS              " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input Frame List Directory              : " << inputdatadir ;
    LOG(INFO) << endl << "Output Directory                               : " << outdir ;
    LOG(INFO) << endl << "Average Factor                               : " << averageFactor ;
    LOG(INFO) << endl << "No Of Frames to be Discarded                               : " << framesToDiscard ;
    if (clobber == 1)
        LOG(INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG(INFO) << endl << "Overwrite                                         : NO" ;
    if (history == 1)
        LOG(INFO) << endl << "History                                             : YES" ;
    else
        LOG(INFO) << endl << "History                                              : NO" ;

    LOG(INFO) << endl << "--------------------------------------------------------------------\n" ;

}

int uvtRefFrameCal::uvtRefFrameCalProcess ()
{

    LOG(INFO) << "\nInside the uvtRefFrameCal  Module" << endl ;
    /*if the Input Directory of the Filelist is not available then exit from module
       else read Filelist     
     */
    LOG(INFO) << "Input Directory : " << inputdatadir << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG(INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;
    string cmd ;
    /*check for the output Directory ,clobber is used for the overwriting  
   if clobber is yes then  overwriting is done.else (i.e if outdirectory is already available & dont want to overwrite)
    exit
     */
    if (DirExists (moduleoutdir) && clobber == YES)
    {
        LOG(INFO) << "Directory exists and clobber=yes" ;
        cmd = (string) "rm -rf " + (string) moduleoutdir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists (moduleoutdir) && clobber == NO)
    {
        LOG(INFO) << endl << moduleoutdir << "  already exists " ;
        LOG(INFO) << endl << "Use clobber=yes for overwriting" ;
        return (EXIT_FAILURE) ;
    }
    /**shell command for Creating the output Directory**/
    cmd = "mkdir -p " + (string) moduleoutdir ;
    /**Executing the Shell command**/
    system (cmd.c_str ()) ;

    /*Searching For the .info  file in the Framelist Directory for the Filelist, '.info' file  contains the listing 
    of the files and some usefull info in the header.
     */
    string tempfilepath = searchFile (inputdatadir , ".info") ;
    if (tempfilepath == " ")
    {
        LOG(ERROR) << endl << "***Information file not found in " << inputdatadir << "***" << endl ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;
    LOG(INFO) << endl << "Information file : " << infofile_in << endl ;
    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists (infofile_in)))
    {
        LOG(INFO) << endl << "***Input FileList not Found at Specified PATH,Check Input Direcrory***" << endl ;
        return (EXIT_FAILURE) ;
    }
    int status = 0 ;
    fitsfile *finfo_in , *finfo_out ;
    /**Opening the input information file**/
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the input information file" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving the  2nd HDU in  input information file" , infofile_in) ;
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
    if (xsize <= 0 || ysize <= 0)
    {
        LOG(ERROR) << endl << "***Invalid xsize/ysize***\n" ;
        return (EXIT_FAILURE) ;
    }
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file
    sprintf (infofile_out , "%s/%s_rfc.info" , moduleoutdir , nameprefix) ;
    LOG(INFO) << "\nOutput information File :" << infofile_out << endl ;
    /**Creating  output information file**/
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information file" , infofile_out) ;
    char *ttype[] = {"centroidFileList"} ;
    char *tform[] = {"A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table in output information file" , infofile_out) ;
    datainfo.write (finfo_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the out information file" , infofile_out) ;
    fits_read_key (finfo_in , TINT , "NFILES" , &nframes , NULL , &status) ;
    printError (status , "Error reading NFILES " , infofile_out) ;
   // if(datainfo.getObsMode ()==PC)
    //reading frame names from information file into vector
    centroidframelist = allocateMemory<char>(nframes , NAMESIZE) ;
  
   if(datainfo.getModeFlag () == IM)
   {
    fits_read_col (finfo_in , TSTRING , 1 , 1 , 1 , nframes , NULL , (void *) centroidframelist , NULL , &status) ;
    printError (status , "Error in reading the centroid  file listing " , infofile_in) ;
      
   }
    else if(datainfo.getModeFlag () == PC)
    {
        fits_read_col (finfo_in , TSTRING , 2 , 1 , 1 , nframes , NULL , (void *) centroidframelist , NULL , &status) ;
    printError (status , "Error in reading the centroid  file listing " , infofile_in) ;
  }
    number = (nframes - framesToDiscard) / averageFactor ; //number of output frames             
    if (datainfo.getModeFlag () == IM || datainfo.getModeFlag () == PC)
    {
        if (refFrameCal ()) return (EXIT_FAILURE) ; //referenceFrameCalculation
    }
    else //in case of neither PC or IM
    {
        LOG(ERROR) << endl << "***Invalid input for operating mode parameter***" ;
        LOG(ERROR) << endl << "***Allowed values are pc/PC/im/IM***" ;
        return (EXIT_FAILURE) ;
    }
    LOG(INFO) << "\nUpdating  output information file" << endl ;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to second HDU in out information file" , infofile_out) ;
    //writing the basic information
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the key value of the NAMEPRFX" , infofile_out) ; //for creating name for output information file
    
    fits_update_key (finfo_out , TSTRING , "CENTDIR" , centroidDir , NULL , &status) ;
    printError (status , "Error in updating the  key value of the centroidDir " , infofile_out) ;
       fits_update_key (finfo_out , TINT , "NFILES" , &number , NULL , &status) ;
    printError (status , "Error in updating the key value of the NFILES " , infofile_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file" , infofile_out) ;
    freeMemory (centroidframelist , nframes , NAMESIZE) ;
    LOG(INFO) << endl << "Reference Frame Calculation  Process Completed Successfully" << endl ;
    return (EXIT_SUCCESS) ;

}

int uvtRefFrameCal::refFrameCal ()
{
    LOG(INFO) << endl << "\nStarted Reference Frame calculation process" << endl ;
    char **outframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    fitsfile *finfo_out ;
    //creating output  directory
    char dir[FLEN_FILENAME] ;
    sprintf (dir , "%s/%s" , moduleoutdir , centroidDir) ;
    /**Shell command for Creating the output Directory**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Execute the Shell command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    LOG(INFO) << endl << "\nTotal number of frames - " << nframes << endl ;
    int status = 0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char errstr[512] ;
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

    float *Xloc ;
    float *Yloc ;
    float *Intensity ;
   
    int size[averageFactor] ;
    for (int j = 0 ; j < averageFactor ; j++)
    {
        size[j] = 0 ;
    }
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
//    long naxes[2] ;
//    naxes[0] = naxes[1] = xsize ;
    vector<float> vect_x , vect_y , vect_int ;
    vector<double> time_Data ;
    int cnt_t = 0 ;
    LOG(INFO) << "\nNumber of Frames to be Averaged :" << averageFactor << endl ;
    LOG(INFO) << "Reference frame calculation started  " << endl ;
    for (int j = framesToDiscard ; j < nframes ; j = j + averageFactor)
    {
        int count = 0 ;
        if ((j + averageFactor) > nframes - 1)
        {
            averageFactor = nframes - j ;
            LOG(INFO) << " Taking group of remaining  " << averageFactor << " frames" << endl ;
        }
        time_Data.clear () ;
        /**Loop for  number of Average frames**/
        for (int i = j ; i < j + averageFactor ; i++)
        {
            sprintf (errstr , "Error at iteration number %d" , i) ;
            fitsfile *fptr ;
            sprintf (infile , "%s/%s/%s" , inputdatadir , centroidDir , centroidframelist[i]) ;
            long numrows = 0 ;
            fits_open_file (&fptr , infile , READONLY , &status) ;
            printError (status , "***Error in opening the input file***") ;

            copyUsrkeywrdsTovect (fptr,key_records);                
            fits_movabs_hdu (fptr , 2 , NULL , &status) ;
            printError (status , "***Error in  moving to the 2nd HDU ***") ;
            fits_read_key (fptr , TUSHORT , "FRAMENO" , &frameno , NULL , &status) ;
            printError (status , "Error reading the key value of the FRAMENO",infile) ;
            fits_read_key (fptr , TDOUBLE , "FRMTIME" , &frametime , NULL , &status) ;
            printError (status , "Error reading the key value of the FRMTIME",infile) ;
            fits_get_num_rows (fptr , &numrows , &status) ;
            printError (status , "***Error in  reading the number of rows***") ;

            size[count] = numrows ;
            Xloc = new float[numrows] ;
            Yloc = new float[numrows] ;
            Intensity = new float[numrows] ;
            
            fits_read_col (fptr , TFLOAT , 1 , 1 , 1 , numrows , NULL , (void *) Xloc , NULL , &status) ;
            printError (status , "***Error in reading the column Xloc of  input frame***") ;
            fits_read_col (fptr , TFLOAT , 2 , 1 , 1 , numrows , NULL , (void *) Yloc , NULL , &status) ;
            printError (status , "Error in reading the column Yloc of input frame***") ;
            fits_read_col (fptr , TFLOAT , 3 , 1 , 1 , numrows , NULL , (void *) Intensity , NULL , &status) ;
            printError (status , "Error in reading the Column Intensity of input frame") ;
           
            fits_close_file (fptr , &status) ;
            time_Data.push_back (frametime) ;
            for (int ii = 0 ; ii < numrows ; ii++)
            {
                vect_x.push_back (Xloc[ii]) ;
                vect_y.push_back (Yloc[ii]) ;
                vect_int.push_back (Intensity[ii]) ;
            }
            count++ ;
        }
        double time_middle ;
        cnt_t++ ;
        time_middle = (time_Data[0] + time_Data[time_Data.size () - 1]) / 2 ;
        float x_ref_arr[vect_x.size ()] , y_ref_arr[vect_y.size ()] , x_arr[vect_x.size ()] , y_arr[vect_y.size ()] ;
        double int_array[vect_x.size ()] ;
        for (int index = 0 ; index < vect_x.size () ; index++)
        {
            x_arr[index] = y_arr[index] = int_array[index] = 0.0 ;
        }
        float temp_x1 , temp_x2 , temp_y1 , temp_y2 , temp_int1 , temp_int2 ;
        vector<float> Xref , Yref , Intref ;
        //loop for finding similar(within neighbourhood ) pixels from all consicutive number of  averagefactor file.
        /**
         *vect_x-contains x locations of centroids of number of  average factor files
         * vect_y-contains y locations of centroids of number of  average factor files
         * vect_int-contains intensity of centroids of number of average factor files
         * size-each location contains number of stars in each file to be processed for find ref frame.
         */
        for (int d1 = 0 ; d1 < size[0] ; d1++)
        {
            int total = size[0] ;
            int cnt = 0 ;
            temp_x1 = vect_x[d1] ;
            temp_y1 = vect_y[d1] ;
            temp_int1 = vect_int[d1] ;
            for (int d2 = 1 ; d2 < averageFactor ; d2++)
            {
                total = total + size[d2] ;
                for (int d3 = total - size[d2] ; d3 < total ; d3++)
                {
                    double diff_x = 0 , diff_y = 0 ;
                    temp_x2 = vect_x[d3] ;
                    temp_y2 = vect_y[d3] ;
                    temp_int2 = vect_int[d3] ;

                    diff_x = temp_x1 - temp_x2 ;
                    diff_y = temp_y1 - temp_y2 ;
                    
                    //finding similar pixels
                    if (sqrt ((diff_x * diff_x)+(diff_y * diff_y)) < NBHD_RADIUS)//NBHD_RADIUS-distance value for finding similar pixel.
                    {
                        x_ref_arr[cnt] = temp_x1 ;
                        y_ref_arr[cnt] = temp_y1 ;
                        x_arr[cnt] = temp_x2 ;
                        y_arr[cnt] = temp_y2 ;
                        int_array[cnt] = temp_int2 ;
                        cnt++ ;
                    }
                }
            }

            float x = temp_x1 ;
            float y = temp_y1 ;
            double inten = temp_int1 ;
            if (cnt >= averageFactor - 1)
            {
                for (int d4 = 0 ; d4 < cnt ; d4++)
                {
                    x = x + x_arr[d4] ;
                    y = y + y_arr[d4] ;
                    inten = inten + int_array[d4] ;
                }
                x = x / averageFactor ;
                y = y / averageFactor ;
                inten = inten / averageFactor ;
                Xref.push_back (x) ;
                Yref.push_back (y) ;
                Intref.push_back (inten) ;
            }
        }
        /**Creating the output Reference frame**/
        fitsfile *fout1 ;
        sprintf (outfile , "%s/%s/%s_%f_f%d_rfc.fits" , moduleoutdir , centroidDir , nameprefix , time_middle , cnt_t) ;
        fits_create_file (&fout1 , outfile , &status) ;
        printError (status , "Error in creating the output  Reference Frame" , outfile) ;
        //cout<<"LLL "<<key_records.size ()<<endl;
//        float *array_final= new float [xsize*ysize];
//        cout<<Xref.size ()<<" "<<xsize<<ysize<<endl;
//        
//        int p=0;
//        for(int i=0;i<Yref.size ();i++)
//        {
//          
//                //if((round(Yref[i]*xsize+Xref[j]))<xsize*ysize){
//            array_final[(int)(round(Yref[i]*xsize+Xref[i]))]=Intref[i];
//           // p++;
//            //    }
//                 
//        }
////        cout<<"LL"<<endl;
////                exit(1);
//        fits_create_img (fout1 , bitpix , naxis , naxes , &status) ;
//        printError (status , "Error in creating the image for the output File " , outfile) ;
//            //write first cut peaks to image
//        fits_write_pix (fout1 , TFLOAT , fpixel , xsize*ysize , array_final , &status) ;
//        printError (status , "Error Writing the pixels to the output File " , outfile) ;
//        delete[] array_final;
        char *tform2[] = {"E" , "E" , "E"} ;
        int tfields = 3 ;
        double t1 = time_Data[0] ;
        double t2 = time_Data[averageFactor - 1] ;
        double interpolated_time = 0.0 ;
        interpolated_time = t1 + ((t2 - t1) / 2) ;
        char *ttype[] = {"X" , "Y" , "Intensity"} ; //char *tform[]={"I","I","E"};
        fits_create_tbl (fout1 , BINARY_TBL , 0 , tfields , ttype , tform2 , NULL , "Reference Frame Peaks" , &status) ;
        printError (status , "Error in creating the table " , outfile) ;
      //write calculated Xref,Yref and Intref to outout file
      
        fits_write_col (fout1 , TFLOAT , 1 , 1 , 1 , Xref.size () , (void *) Xref.data () , &status) ;
        printError (status , "Error in writing the column of Xref" , outfile) ;
        fits_write_col (fout1 , TFLOAT , 2 , 1 , 1 , Yref.size () , (void *) Yref.data () , &status) ;
        printError (status , "Error in writing the column of Yref" , outfile) ;
        fits_write_col (fout1 , TFLOAT , 3 , 1 , 1 , Intref.size () , (void *) Intref.data () , &status) ;
        printError (status , "Error in writing the column of Intensity" , outfile) ;
        fits_write_key (fout1 , TDOUBLE , "TIME" , &interpolated_time , NULL , &status) ;
        printError (status , "Error in writing the TIME in output Centroid file" , outfile) ;
        
       
        
       fits_close_file (fout1 , &status) ;
       printError (status , "Error in closing the file" , outfile) ;
       
       fits_open_file (&fout1 , outfile , READWRITE , &status) ;
       printError (status , "Error in opening the input File" , outfile) ;
        //update frametime and frameno to output file       
         for(int i=1;i<=2;i++)
         {
             
           fits_movabs_hdu (fout1 , i , NULL , &status) ;
           printError (status , "Error in  moving to the  HDU") ;  
           fits_update_key (fout1 , TDOUBLE , "FRMTIME" , &time_middle , NULL , &status) ;
           printError (status , "Error in writing the key value of the FRMTIME" , outfile) ;
           fits_update_key (fout1 , TUSHORT , "FRAMENO" , &cnt_t , NULL , &status) ;
           printError (status , "Error in writing the key value of the FRMENO" , outfile) ;

         }
       //write level-1 keywords from vector,history,origin,checksum,creator and date to output file.
          writeUsrkeywordsFrmvect (outfile,key_records);
          if (history == YES) writeHistory (outfile , vhistorystr) ;
          fits_movabs_hdu (fout1, 1 , NULL , &status) ;
          printError (status , "***Error in  moving to the 2nd HDU ***") ;
          writeCommonKeywords (fout1,modulename);  
          fits_close_file (fout1 , &status) ;
           printError (status , "***Error in  closing the file ***") ;
        vect_x.clear () ;
        vect_y.clear () ;
        vect_int.clear () ;
        strcpy (outframelist[cnt_t - 1] , basename (outfile)) ;
        
    }
  
    LOG(INFO) << "\nReference Frame calculation completed" << endl ;
    /**write columns of the output framelist**/
    /**opening output information file**/
    LOG(INFO) << "\nWriting list of output  frame names to output  information file  " << endl ;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening  the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in  moving to the 2nd HDU of the out information file" , infofile_out) ;
    fits_write_col (finfo_out , TSTRING , 1 , 1 , 1 , number , (void*) outframelist , &status) ;
    printError (status , "Error in writing the column of outout Exposure frame list" , infofile_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in Closing the  output information file" , infofile_out) ;
    //write level1 keywords from vector,history  to putput information file.
    writeUsrkeywordsFrmvect (infofile_out,key_records);
     if (history == YES) writeHistory (infofile_out , vhistorystr) ;
      vhistorystr.clear () ;
    return (EXIT_SUCCESS) ;
}

int uvtRefFrameCal::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    char average_f[3] ;
    char s_numOfframe_Disc[3] ;

    sprintf (average_f , "%d" , averageFactor) ;
    sprintf (s_numOfframe_Disc , "%d" , framesToDiscard) ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" inputdatadir =" + (string) inputdatadir) ;
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

    vhistory.push_back ((string)getSerialNo (cnt)+" Average-Factor=" + (string) average_f) ;
    vhistory.push_back ((string)getSerialNo (cnt)+" Frame Discarded=" + (string) s_numOfframe_Disc) ;
    vhistory.push_back ("Parameter List END") ;
    return (EXIT_SUCCESS) ;
}

