/* 
 * File:   uvtRelAspCal.cpp
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
#include<cstdlib>
#include<stdio.h>
#include<fstream>
#include<uvtRelAspCal.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<spMatrix.h>
#include<glog/logging.h>


#define DRIFT_PATH "uvtDriftExersize_1.0/"
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr   //used in doFft() function

//Constructor -called when object is created


uvtRelAspCal::uvtRelAspCal ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
    // strcpy(centroidDir, "Centroid");
}

//Destructor


uvtRelAspCal::~ uvtRelAspCal () {
    // delete[] x_Distortion;
    // delete[] y_Distortion;
 }

//parameter File reading


int uvtRelAspCal::read (int argc , char** argv)
{
    int status = 0 ;


    if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG(ERROR) << "***Error Initializing PIL***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("drift_dir" , driftdir)))
    {
        LOG(ERROR) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("jitter_dir" , jitterdir)))
    {
        LOG(ERROR) << endl << "***Error reading output directory name***" ;
        return status ;
    }
    if (PIL_OK != (status = PILGetFname ("thermal_dir" , thermaldir)))
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


int uvtRelAspCal::read (char *drift_dir , char *jitter_dir , char *thermal_dir , char *outdir , int clobber , int history)
{
    strcpy (this->driftdir , drift_dir) ;
    strcpy (this->jitterdir , jitter_dir) ;
    strcpy (this->thermaldir , thermal_dir) ;
    //strcpy(this->caldbDir, caldbDir);
    strcpy (this->outdir , outdir) ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}
//Parameter file content Display


void uvtRelAspCal::display ()
{

    LOG(ERROR) << endl << "---------- Relative Aspect Series Parameter  Display---------" ;
    LOG(ERROR) << endl << "Input Frame List Directory                        : " << inputdatadir ;
    LOG(ERROR) << endl << "Caldb Directory                                           : " << caldbDir ;
    LOG(ERROR) << endl << "Output Directory                               : " << outdir ;
    if (clobber == YES)
        LOG(ERROR) << endl << "Overwrite                                         : YES" ;
    else
        LOG(ERROR) << endl << "Overwrite                                         : NO" ;
    if (history == YES)
        LOG(ERROR) << endl << "History                                             : YES" ;
    else
        LOG(ERROR) << endl << "History                                              : NO" ;

    LOG(ERROR) << endl << "----------Relative Aspect Series Display Ends----------\n" ;

}


int uvtRelAspCal::uvtRelAspCalProcess ()
{
    /*if the Input Directory of the Filelist is not available then exit from module
     else read Filelist

     */
    LOG(ERROR) << endl << "Relative Aspect Series Calculation process started" << endl ;
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    //LOG(INFO)<<endl<<"Module Output Directory : "<<moduleoutdir<<endl;
    string cmd ;
    if (DirExists (moduleoutdir) && clobber == YES)
    {
        LOG(ERROR) << endl << "Directory exists and clobber=yes" ;
        cmd = (string) "rm -rf " + (string) moduleoutdir ;
        system (cmd.c_str ()) ;
    }
    else if (DirExists (moduleoutdir) && clobber == NO)
    {
        LOG(ERROR) << endl << moduleoutdir << "  already exists " ;
        LOG(ERROR) << endl << "Use clobber=yes for overwriting" ;
        return (EXIT_FAILURE) ;
    }

    nframes = 0 ;
    cmd = "mkdir -p " + (string) moduleoutdir ;
    system (cmd.c_str ()) ; // creating output directory to keep output from unitConversion
    LOG(ERROR) << endl << moduleoutdir << "  directory created" ;

    string tempfilepath = searchFile (driftdir , ".info") ;
     if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , driftdir , tempfilepath.c_str()) ;
    int status = 0 ;

    fitsfile *finfo_in , *finfo_out ;

    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in opening the information file" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in Moving the 2nd HDU" , infofile_in) ;
    datainfo.getInfo (finfo_in) ;
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "Error in reading the key value of the NAMEPRFX " , infofile_in) ; //for creating name for output information file

    tempfilepath = searchFile (driftdir , ".fits") ;
     if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding drift file***" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (driftfile , "%s/%s" , driftdir , tempfilepath.c_str()) ;

    tempfilepath = searchFile (jitterdir , ".fits") ;
     if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding jitter file***" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (jitterfile , "%s/%s" , jitterdir , tempfilepath.c_str()) ;

    if(strcmp (datainfo.getDetector (), "VIS")!=0)
    {
    tempfilepath = searchFile (thermaldir , ".fits") ;
     if (tempfilepath =="")
    {
        LOG (ERROR) << endl << "***Error in finding thermal  file***" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (thermalfile , "%s/%s" , thermaldir , tempfilepath.c_str()) ;
    }
    sprintf (infofile_out , "%s/%s_ras.info" , moduleoutdir , nameprefix) ;
    LOG(INFO) << "The information file  out is " << infofile_out << endl ;
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "***Error in creating the output information file***") ;
    char *ttype[] = {"CentroidFrames"} ;
    char *tform[] = {"A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 1 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "***Error in creating the table***") ;
    datainfo.write (finfo_out) ; //writing basic data information

    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "***Error in updating the key value of NAMEPRFX*** ") ; //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file " , infofile_out) ; //for creating name for output information file
    int mode = 1 ;
    if (datainfo.getModeFlag () == IM || datainfo.getModeFlag () == PC)
    {//For IM mode
        if (computeRelAsp ()) return (EXIT_FAILURE) ;
        fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
        printError (status , "***Error in opening the out  information file***") ;
        fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
        printError (status , "***Error in moving the 2nd HDU of  information file***") ;

        fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
        printError (status , "***NAMEPRFX keyword  not updated/not Found***") ; //for creating name for output information file

        fits_close_file (finfo_out , &status) ;
        printError (status , "***Error in closing the file***") ;

        //freeMemory(centroidframelist, nframes, NAMESIZE);//for releasing the memory
    }//in case of the PC
    else
    {
        LOG(ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG(ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }


    return (EXIT_SUCCESS) ;

}


int uvtRelAspCal::computeRelAsp ()
{

    LOG(ERROR) << endl << " Inside the Relative aspect calculation" << endl ;
    int status = 0 ;

    fitsfile *frelasp ;
    vector<string> vhistorystr ;
    char infile_ras[FLEN_FILENAME] ;
    sprintf (infile_ras , "%s/%s_ras.fits" , moduleoutdir , nameprefix) ;
    // " << infile_ras << endl ;
    bool vis_flag =FALSE;

    fitsfile *fptr ;
    fits_open_file (&fptr , driftfile , READONLY , &status) ;
    printError (status , "input file reading fails",driftfile) ;
    long thermal_noRecords =0;
    copyUsrkeywrdsTovect (fptr , key_record) ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the file" , driftfile) ;
   //  this->setRelAspFilename (infile_ras);
    //cout<<"After"<<endl;exit(1);
   //strcpy (RASfile,infile_ras);
    //cout<<RASfile<<endl;exit(1);
    fits_create_file (&frelasp , infile_ras , &status) ;
    printError (status , "Error in creating the file" , infile_ras) ;
    char *ttype[] = {"Time" , "Yaw" , "Pitch" , "Roll"} ;
    char *tform[] = {"1D" , "1D" , "1D" , "1D"} ;
    fits_create_tbl (frelasp , BINARY_TBL , 0 , 4 , ttype , tform , NULL , "BINTABLE" , &status) ;
    printError (status , "Error in creating the table" , infile_ras) ;

    double *time_data , *Xshiftdata , *Yshiftdata , *thetadata ;
    double *time_data1 , *Xshiftdata1 , *Yshiftdata1 , *thetadata1 ;
    double *time_data2 , *Xshiftdata2 , *Yshiftdata2 , *thetadata2 ;

    vector<double> timedataVect , timedataVect1 , timedataVect2,xshiftVect , yshiftVect , thetaVect , xshiftVect1 , yshiftVect1 , thetaVect1,xshiftVect2 , yshiftVect2 , thetaVect2  ;
    //   double *time_data2,*Xshiftdata2,*Yshiftdata2,*thetadata2;
    long drift_noRecords = readFitsFile (driftfile , time_data1 , Xshiftdata1 , Yshiftdata1 , thetadata1 , timedataVect1 , xshiftVect1 , yshiftVect1 , thetaVect1) ;
    long jitter_noRecords = readFitsFile (jitterfile , time_data , Xshiftdata , Yshiftdata , thetadata , timedataVect , xshiftVect , yshiftVect , thetaVect) ;
    
    if(strcmp (datainfo.getDetector (), "VIS")!=0)
    {
       thermal_noRecords = readFitsFile (thermalfile , time_data2 , Xshiftdata2 , Yshiftdata2 , thetadata2 , timedataVect2 , xshiftVect2 , yshiftVect2 , thetaVect2) ;
       vis_flag=TRUE;
    }
   
    
 
    //added
    vector<double> jitter_xshift , jitter_yshift , jitter_theta ;
    vector<double> thermal_xshift , thermal_yshift , thermal_theta ;
    double temp_x , temp_y , temp_theta ;
    bool flag_jitterFound=0,flag_thermalFound=0;
   
    for (int i = 0 ; i < drift_noRecords ; i ++)
    {
        flag_jitterFound=0;
        for (int j = 1 ; j < jitter_noRecords - 1 ; j ++)
        {
            if (timedataVect[i] < timedataVect1[j] && timedataVect[i] >= timedataVect1[j - 1])
            {
                flag_jitterFound=1;
                temp_x = xshiftVect[j - 1]+((xshiftVect[j] - xshiftVect[j - 1]) / (timedataVect[j] - timedataVect[j - 1]))*(timedataVect1[i] - timedataVect[j - 1]) ;
                temp_y = yshiftVect[j - 1]+((yshiftVect[j] - yshiftVect[j - 1]) / (timedataVect[j] - timedataVect[j - 1]))*(timedataVect1[i] - timedataVect[j - 1]) ;
                temp_theta = thetaVect[j - 1]+((thetaVect[j] - thetaVect[j - 1]) / (timedataVect[j] - timedataVect[j - 1]))*(timedataVect1[i] - timedataVect[j - 1]) ;
                jitter_xshift.push_back (temp_x) ;
                jitter_yshift.push_back (temp_y) ;
                jitter_theta.push_back (temp_theta) ;
//                drift_xshift.push_back (Xshiftdata1[i]);
//                drift_yshift.push_back (Yshiftdata1[i]);
//                drift_theta.push_back (thetaVect1[i]);
                
                break ;
            }
        }
        if(flag_jitterFound==0){
            jitter_xshift.push_back (0.0f);
            jitter_yshift.push_back (0.0f);
            jitter_theta.push_back (0.0f);
        }
    }
    
    //for matching thermal content as per time series of drift ;
    
    if(vis_flag==FALSE)
    {
     for (int i = 0 ; i < drift_noRecords ; i ++)
    {
        flag_thermalFound=0;
        
        for (int j = 1 ; j < thermal_noRecords - 1 ; j ++)
        {
            if (timedataVect2[i] < timedataVect1[j] && timedataVect2[i] >= timedataVect1[j - 1])
            {
                flag_thermalFound=1;
                temp_x = xshiftVect2[j - 1]+((xshiftVect2[j] - xshiftVect2[j - 1]) / (timedataVect2[j] - timedataVect2[j - 1]))*(timedataVect1[i] - timedataVect2[j - 1]) ;
                temp_y = yshiftVect2[j - 1]+((yshiftVect2[j] - yshiftVect2[j - 1]) / (timedataVect2[j] - timedataVect2[j - 1]))*(timedataVect1[i] - timedataVect2[j - 1]) ;
                temp_theta = thetaVect2[j - 1]+((thetaVect2[j] - thetaVect2[j - 1]) / (timedataVect2[j] - timedataVect2[j - 1]))*(timedataVect1[i] - timedataVect2[j - 1]) ;
                thermal_xshift.push_back (temp_x) ;
                thermal_yshift.push_back (temp_y) ;
                thermal_theta.push_back (temp_theta) ;
               
                break ;
            }
        }
        if(flag_thermalFound==0)
        {
            thermal_xshift.push_back (0.0f);
            thermal_yshift.push_back (0.0f);
            thermal_theta.push_back (0.0f);
        }
    }
    }
    
    

    double *sum_roll_data1 = new double[drift_noRecords] ;
    double *sum_pitch_data1 = new double[drift_noRecords] ;
    double *sum_yaw_data1 = new double[drift_noRecords] ;

    for (int i = 0 ; i < drift_noRecords ; i ++)
    {
        sum_yaw_data1[i] = 0.0 ;
        sum_pitch_data1[i] = 0.0 ;
        sum_roll_data1[i] = 0.0 ;
        //cout<<jitter_yshift[i]<<" "<<thermal_yshift[i]<<" "<<yshiftVect1[i]<<endl;
        sum_yaw_data1[i] = jitter_xshift[i] +xshiftVect1[i];
        sum_pitch_data1[i] = jitter_yshift[i] +yshiftVect1[i];
        sum_roll_data1[i] = jitter_theta[i] +thetaVect1[i];
        
        if(vis_flag==FALSE)
        {
        sum_yaw_data1[i] = sum_yaw_data1[i]+ thermal_xshift[i] ;
        sum_pitch_data1[i] =sum_pitch_data1[i]+ thermal_yshift[i] ;
        sum_roll_data1[i] = sum_roll_data1[i]+thermal_theta[i] ;
        }

    }

    //convert it into dx,dy and dtheta..

    fits_write_col (frelasp , TDOUBLE , 1 , 1 , 1 , drift_noRecords , timedataVect1.data () , &status) ;
    printError (status , "Error in writing the  time data to the relative aspect file" , infile_ras) ;
    fits_write_col (frelasp , TDOUBLE , 2 , 1 , 1 , drift_noRecords , sum_yaw_data1 , &status) ;
    printError (status , "Error in writing the yaw  data to the relative aspect file" , infile_ras) ;
    fits_write_col (frelasp , TDOUBLE , 3 , 1 , 1 , drift_noRecords , sum_pitch_data1 , &status) ;
    printError (status , "Error in writing the pitch data to the relative aspect file" , infile_ras) ;
    fits_write_col (frelasp , TDOUBLE , 4 , 1 , 1 , drift_noRecords , sum_roll_data1 , &status) ;
    printError (status , "Error in writing the roll  data to the relative aspect file" , infile_ras) ;

   
    writeUsrkeywordsFrmvect (infile_ras , key_record) ;
    writeCommonKeywords (frelasp , modulename) ;
    fits_close_file (frelasp , &status) ;
    printError (status , "***Error in closing the File***") ;
 if (history == YES)
    {
        getHistory (vhistorystr) ;
        writeHistory (infile_ras , vhistorystr) ;
    }
    return (EXIT_SUCCESS) ;

}


int uvtRelAspCal::getHistory (vector<string> &vhistory)
{
    int cnt = 0 ;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Input Drift directory=" + (string) driftdir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Input Jitter directory=" + (string) jitterdir) ;
    vhistory.push_back ((string) getSerialNo (cnt) + " Input thermal  directory=" + (string) jitterdir) ;
   ;
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


int uvtRelAspCal::readFitsFile (char* filename , double *timedata , double* Xdata , double* Ydata ,
        double* thetadata , vector<double> &timeterm , vector<double> &xterm ,
        vector<double> &yterm , vector<double> &thetaterm)
{
    int status = 0 ;
    fitsfile *fptr ;
    fits_open_file (&fptr , filename , READONLY , &status) ;
    printError (status , "***input file reading fail***",filename) ;
    long numrows = 0 ;
    fits_movabs_hdu (fptr , 2 , NULL , &status) ;
    printError (status , "***Moving to  2nd HDU Fails***",filename) ;
    fits_get_num_rows (fptr , &numrows , &status) ;
    timedata = new double[numrows] ;
    Xdata = new double[numrows] ;
    Ydata = new double[numrows] ;
    thetadata = new double[numrows] ;
    fits_close_file (fptr , &status) ;
    printError (status , "Error in closing the file " , filename) ;
    // exit(1);

    status = readColumnsFromFITS (filename , 2 , 4 , TDOUBLE , 1 , timedata , (int) numrows , TDOUBLE , 2 ,
            Xdata , (int) numrows , TDOUBLE , 3 , Ydata , (int) numrows , TDOUBLE , 4 , thetadata , (int) numrows) ;
    if (status)
    {
        LOG (INFO) << "Error reading  the columns from the file" << endl ;
        return (EXIT_FAILURE) ;
    }
    xterm.clear () ;
    yterm.clear () ;
    thetaterm.clear () ;

    for (int i = 0 ; i < numrows ; i ++)
    {
        timeterm.push_back (timedata[i]) ;
        xterm.push_back (Xdata[i]) ;
        yterm.push_back (Ydata[i]) ;
        thetaterm.push_back (thetadata[i]) ;
    }
     return numrows ;
}
