/* 
 * File:   uvtFlatFieldCorr.cpp
 * Authors:: Dhruv, Sanjay K Singh, Arvind K Singh
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
#include <uvtFlatFieldCorr.h>
#include<pthread.h>
#include<uvtUtils.h>
#include<glog/logging.h>

//int status=0;
//fitsfile *fitin;

//constructor-called when the object is created 

uvtFlatFieldCorr::uvtFlatFieldCorr ()
{
    sprintf (modulename , "%s_%s" , MODULENAME , VERSION) ;
}
//Destructor-called when the obj is no more necessary

uvtFlatFieldCorr::~uvtFlatFieldCorr () {
   
}

/*For Reading the parameter file content*/
int uvtFlatFieldCorr::read (int argc , char** argv)
{
    int status = 0 ;
   status=readParams (argc , argv , 6 , FNAME , "inputdatadir" , inputdatadir,FNAME,"caldbDir",calDir,FNAME,"outdir",outdir,BOOL,"clobber",&clobber,
            BOOL,"history",&history,STRING ,"mode",&mode) ;
   
   if (status) return (EXIT_FAILURE) ;
  }

int uvtFlatFieldCorr::read (char *inputdatadir , char *caldb_dir , char *out_dir , int clobber , int history)
{
    strcpy (this->inputdatadir , inputdatadir) ;
    strcpy (this->calDir , caldb_dir) ;
    strcpy (this->outdir , out_dir) ;
    this->clobber = clobber ;
    this->history = history ;
    return (EXIT_SUCCESS) ;
}

/*For Displaying  the content of parameter file*/
void uvtFlatFieldCorr::display ()
{
     LOG(INFO) << endl ;
    LOG(INFO) << "-----------------------------------------------------------------------" << endl ;
    LOG(INFO) << "             UVT FLAT FIELD CORRECTION PARAMETERS      " << endl ;
    LOG(INFO) << "------------------------------------------------------------------------" ;
    LOG(INFO) << endl << "Input Frame List Directory                        : " << inputdatadir ;
//    LOG(INFO) << endl << "Flat Field  File                                     : " << flatfieldfile ;
    LOG(INFO) << endl << "Output Directory                               : " << outdir ;
    if (clobber == YES)
        LOG(INFO) << endl << "Overwrite                                         : YES" ;
    else
        LOG(INFO) << endl << "Overwrite                                         : NO" ;
    if (history == YES)
        LOG(INFO) << endl << "History                                             : YES" ;
    else
        LOG(INFO) << endl << "History                                              : NO" ;
    //  LOG(INFO)<<endl<<"-----------------------------------------------------------------------------------";
     LOG(INFO) << endl << "------------------------------------------------------------------------" << endl ;

}

int uvtFlatFieldCorr::uvtFlatFieldCorrProcess ()
{
    /*if the Input Directory of the Filelist is not available then exit from module
     else read Filelist 
     */
    LOG(INFO) << endl << "\nInside  FlatFieldCorrection Module" << endl ;
    /*--------*/
    sprintf (moduleoutdir , "%s/%s/" , outdir , modulename) ;
    LOG(INFO) << endl << "Module Output Directory : " << moduleoutdir << endl ;
    
     //check  existence  of output directory if output  exist and 
    //clobber =YES than remove that directory and recreate it,if exist and clobber =NO than exit from module.
    //if directory is not exist than create it .
    if (createOutputDirectory (clobber , moduleoutdir))
        return (EXIT_FAILURE) ;
   
    /*Searching For the .info  file in the Framelist Directory for the Filelist, '.info' file  contains the listing 
     of the files and some usefull info in the header.
     */
   /*----Searching the input information File from the input Directory----*/
    string  tempfilepath = searchFile (inputdatadir , "info") ;
     if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }
    sprintf (infofile_in , "%s/%s" , inputdatadir , tempfilepath.c_str()) ;

    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists (infofile_in)))
    {
        LOG(INFO) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY" ;
        return (EXIT_FAILURE) ;
    }

    LOG(INFO) << endl << "\nInformation File :" << infofile_in << endl ;
    int status = 0 ;
    /*----Reading needed information from the input information file----*/
    fitsfile *finfo_in , *finfo_out ;
    /*  open the .info FITS file  and read the header information from the second HDU.  */
    fits_open_file (&finfo_in , infofile_in , READONLY , &status) ;
    printError (status , "Error in Opening the input information File" , infofile_in) ;
    fits_movabs_hdu (finfo_in , 2 , NULL , &status) ;
    printError (status , "Error in moving to the 2nd HDU in input information File" , infofile_in) ;

    /**Reading the keyword values from the  input information File**/
    datainfo.getInfo (finfo_in) ; //reading basic information for data from information file
    xsize = datainfo.getXsize () ;
    ysize = datainfo.getYsize () ;
   
    fits_read_key (finfo_in , TSTRING , "NAMEPRFX" , nameprefix , NULL , &status) ;
    printError (status , "NAMEPRFX keyword not Found" , infofile_in) ; //for creating name for output information file
      /*   Creating a output list  FITS File which will be used in the next Module   */
    /**Setting the output information File**/
    sprintf (infofile_out , "%s/%s_ff.info" , moduleoutdir , nameprefix) ;
    LOG(INFO)<<"\nOutout information File "<<infofile_out<<endl;
    /*Creating the output information File*/
    fits_create_file (&finfo_out , infofile_out , &status) ;
    printError (status , "Error in creating the output information File" , infofile_out) ;
    char *ttype[] = {"SignalFrames" , "ExposureFrames"} ;
    char *tform[] = {"A256" , "A256"} ;
    fits_create_tbl (finfo_out , ASCII_TBL , 0 , 2 , ttype , tform , NULL , "FileList" , &status) ;
    printError (status , "Error in creating the table in output information File" , infofile_out) ;
    fits_update_key (finfo_out , TSTRING , "NAMEPRFX" , nameprefix , "File name prefix" , &status) ;
    printError (status , "Error in updating the key value of the NAMEPRFX " , infofile_out) ;
     /**writing the keyword values to the output information File**/
    datainfo.write (finfo_out) ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in Closing the output information file" , infofile_out) ;
   /*
     getFlatFieldFile-method finds  the path for the flatfield FIT file from the arguments ,
      this FIT file will be used for the correction in  the input File.
     */
    LOG(INFO)<<"Setting  the path of input flat field Correction file from calDB directory"<<endl;
    string  namepath = caldb_handler.getFlatFieldFile (datainfo.getDetector () , datainfo.getObsMode () , datainfo.getFilter () , calDir) ;
    if (namepath == " ")
    {
        LOG(ERROR) << endl << "Couldn't find flat field correction  file from caldb" << endl ;
        return (EXIT_FAILURE) ;
    }
       joinStrings (flatfieldfile , 2 , calDir , namepath.c_str()) ;
    // cout<<"OUTdie "<<endl;exit(1);
    /* reading FlatField File into 'flatfielddata' buffer */
    LOG(INFO)<<"\nReading FlatField file from calDB  directory"<<endl;
    readFlatFieldFile () ;
    
   //allocating memory to array  for padded caldb data 
    
           
          if (datainfo.getModeFlag () == IM)
    {
          finalFlatFielddata = new float[xsize * ysize] ;
    for(int i=0;i<xsize*ysize;i++){
        finalFlatFielddata[i]=0.0f;
    }
   
    /*padding is done to  match xsize of event file to caldb file*/
    //status = Applypadding (flatfielddata ,flatfield_dim1,flatfield_dim2, finalFlatFielddata,xsize,ysize) ;
    status = Applypadding (flatfielddata ,flatfield_dim1,flatfield_dim2, finalFlatFielddata,xsize,ysize) ;
    if (status)
    {
        LOG(ERROR) << "***Padding Fails For the flatfield CALDB data***" << endl ;
        return (EXIT_FAILURE) ;
    }
       /*SIGDIR-directory path for reading the signal frames
         * EXPDIR-directory path for reading the exposure frames 
         * NFILES-number of frames  to be read
         */
             readKeywords (infofile_in , 2, 3 , TSTRING , "SIGDIR" , sigframedir ,
                TSTRING , "EXPDIR" , expframedir ,
                TINT , "NFILES" , &nframes) ;        
        LOG(INFO) << endl << "Number of files :" << nframes ;
        //reading frame names from information file into vector
        sigframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        expoframelist = allocateMemory<char>(nframes , NAMESIZE) ;
        
        status=readColumnsFromFITS (infofile_in,2,2,TSTRING,1,sigframelist,nframes,TSTRING ,2,expoframelist,nframes);
        if(status){
        LOG(INFO)<<"Error reading  the columns from the file"<<endl;
        return(EXIT_FAILURE);
        }
        /*
         steps to be taken for IM mode
         */
        if (flatFieldIM ()) return (EXIT_FAILURE) ;
        fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
        printError (status , "Error in opening the output information file" , infofile_out) ;
        fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
        printError (status , "Error in moving to the 2nd HDU in output information file" , infofile_out) ;
         datainfo.write (finfo_out) ;
         fits_close_file (finfo_out , &status) ;
        printError (status , "Error in closing the output information file" , infofile_out) ;
        
        updateKeywords (infofile_out , 2 , 4 , TSTRING,"NAMEPRFX",nameprefix,TSTRING , "SIGDIR" , sigframedir ,TSTRING,"EXPDIR",expframedir,TINT,"NFILES",&nframes);//updating the keywords
        
          fits_close_file (finfo_in , &status) ;
         printError (status , "Error in closing the input information File" , infofile_out) ;
        freeMemory (sigframelist , nframes , NAMESIZE) ;
    }
    else if (datainfo.getModeFlag () == PC)
    { //if PC mode 
        finalFlatFielddata = new float[xsize * ysize] ;
    for(int i=0;i<xsize*ysize;i++){
        finalFlatFielddata[i]=0.0f;
    }
   
    /*padding is done to  match xsize of event file to caldb file*/
    //status = Applypadding (flatfielddata ,flatfield_dim1,flatfield_dim2, finalFlatFielddata,xsize,ysize) ;
    status = Applypadding (flatfielddata ,flatfield_dim1,flatfield_dim2, finalFlatFielddata,600,600) ;
    if (status)
    {
        LOG(ERROR) << "***Padding Fails For the flatfield CALDB data***" << endl ;
        return (EXIT_FAILURE) ;
    }
        fits_read_key (finfo_in , TSTRING , "EVTFILE" , eventfile , NULL , &status) ;
        printError (status , "Error in reading the key value of the EVTFILE" , infofile_in) ;
        /*steps to be taken for pc mode */
       if (flatFieldPC ()) return (EXIT_FAILURE) ;
         fits_close_file (finfo_in , &status) ;
         printError (status , "Error in closing the input information File" , infofile_out) ;
     }
    else
    {
        LOG(ERROR) << endl << "Invalid input for operating mode parameter" ;
        LOG(ERROR) << endl << "Allowed values are pc/PC/im/IM" ;
        return (EXIT_FAILURE) ;
    }
    return (EXIT_SUCCESS) ;
}

int uvtFlatFieldCorr::flatFieldPC ()
{
    LOG(INFO) << endl << "\nStarted FlatField correction for PC mode" << endl ;    
    fitsfile *fptr , *fout , *finfo_out ;
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    char errstr[512] ;
    int colnum=0;
    int status = 0 ;
    char eventfilename[FLEN_FILENAME] ;
    sprintf (eventfilename , "%s_ff.events" , nameprefix) ;
    sprintf (outfile , "%s/%s" , moduleoutdir , eventfilename) ;
    fits_create_file (&fout , outfile , &status) ;
    printError (status , errstr) ;

    long nrows = 0 ;
    /**Setting the input EVENT file path**/
    sprintf (infile , "%s/%s" , inputdatadir , eventfile) ;
    LOG(INFO) <<"\nInput Event File " << infile << endl ;
    LOG(INFO)<<"\nOutput Event File "<<outfile<<endl;
    /**opening the inpur EVENT File**/
    fits_open_file (&fptr , infile , READWRITE , &status) ;
    printError (status , "Error in opening the infile" , infile) ;
    copyUsrkeywrdsTovect (fptr,key_records);
    /**Copy the input File to the output File**/
    fits_copy_file (fptr , fout , 1 , 1 , 1 , &status) ;
    printError (status , "Error in coping the file" , outfile) ;
    float *xFrac, *yFrac, *ENP;
    fits_movabs_hdu (fout , 2 , NULL , &status) ;
    printError (status , "***Moving To particular HDU fails***" , outfile) ;
    fits_get_num_rows (fout , &nrows , &status) ;
    printError (status , "Error reading the row number of the output file" , outfile) ;
    LOG(INFO)<<"Total number of Events are "<<nrows<<endl;
    xFrac = new float[nrows];
    yFrac = new float[nrows];
    ENP = new float[nrows];
    for(int i=0;i<nrows;i++){
        ENP[i]=1.0f;
        }
      /**Reading the input event File from  the input Directory**/
//     fits_read_col (fout , TFLOAT , 3 , 1 , 1 , nrows , NULL , xFrac , NULL , &status) ;
//    printError (status , "Error reading Column" , outfile) ;
//    fits_read_col (fout , TFLOAT , 4 , 1 , 1 , nrows , NULL , yFrac , NULL , &status) ;
//    printError (status , "Error reading Column" , outfile) ;
    fits_read_col (fout , TFLOAT , 4 , 1 , 1 , nrows , NULL , xFrac , NULL , &status) ;
    printError (status , "Error reading Column" , outfile) ;
    fits_read_col (fout , TFLOAT , 5 , 1 , 1 , nrows , NULL , yFrac , NULL , &status) ;
    printError (status , "Error reading Column" , outfile) ;
    fits_get_colnum(fout,CASEINSEN,FF_COLNAME,&colnum,&status);
    printError (status , "Error reading the column  number of the output file" , outfile) ;
    fits_read_col (fout , TFLOAT , colnum , 1 , 1 , nrows , NULL , ENP , NULL , &status) ;
    printError (status , "Error reading Column" , outfile) ;
    LOG(INFO)<<"\nApplying flat field Correction on Effective number of Photons for each row... "<<endl;
/**Calculating The Effective number of photons for each Centroid**/
   
    int  divi_fact=xsize/600;
    for (int i = 0 ; i < nrows ; i++)
    {
                  // ENP[i] =ENP[i]*finalFlatFielddata[(int)(round((yFrac[i]))*xsize + (int)round((xFrac[i])))] ;
        ENP[i] =ENP[i]*finalFlatFielddata[(int)(round((yFrac[i])/divi_fact)*600 + (int)round((xFrac[i])/divi_fact))] ;
    }
   /**Write updated value of the  ENP to output event file**/
    fits_write_col (fout , TFLOAT, colnum , 1 , 1 , nrows , ENP , &status) ;
    printError (status , "Error writing to Column" , outfile) ;
    delete[] xFrac,yFrac,ENP,finalFlatFielddata;
    fits_movabs_hdu(fout, 1, NULL, &status);
    printError(status, "Error in moving to 1st HDU",infofile_out);
    /***Writes the  common information to  output information File***/
   
     vector<string> vhistorystr;
    if (history == YES){       
        getHistory(vhistorystr);
        LOG(INFO)<<"\nWriting History"<<endl;
        writeHistory(outfile, vhistorystr);
    }
    writeCommonKeywords(fout, modulename);
    fits_close_file(fout, &status);
    printError(status, "Error in Closing the output Event File",outfile);
    fits_close_file(fptr, &status);
    printError(status, "Error in closing the input Event File",outfile);
    /**open output information File  which has been generated earlier**/
    LOG(INFO)<<"\nUpdating output information file"<<endl;
    fits_open_file (&finfo_out , infofile_out , READWRITE , &status) ;
    printError (status , "Error in opening the output information file" , infofile_out) ;
    fits_movabs_hdu (finfo_out , 2 , NULL , &status) ;
    printError (status , "Error in moving to 2nd HDU in output information file" , infofile_out) ;
    fits_update_key (finfo_out , TSTRING , "EVTFILE" , basename (outfile) , NULL , &status) ;
    printError (status , "Error in updating the key value of the output information file") ;
    fits_close_file (finfo_out , &status) ;
    printError (status , "Error in closing the output information file  ",infofile_out) ;
     writeUsrkeywordsFrmvect (infofile_out,key_records);
     if(history==YES)
     writeHistory (infofile_out , vhistorystr) ;
     vhistorystr.clear ();
    return (EXIT_SUCCESS) ;
}

/* For IM mode */
int uvtFlatFieldCorr::flatFieldIM ()
{
   LOG(INFO) << endl << "\nStarted FlatField correction for IM mode" << endl ;    
    /**allocating the memory to the output signal frame list and Exposure frame list**/
    char **outsigframelist = allocateMemory<char >(nframes , NAMESIZE) ; // to store output frame list
    char **outexpframelist = allocateMemory<char>(nframes , NAMESIZE) ; //to store output exposure frame list  
    //creating signal directory
    char dir[FLEN_FILENAME] ;
    /**setting the path for the output signal Directory**/
    sprintf (dir , "%s%s" , moduleoutdir , sigframedir) ;
    /**Shell command for creating the Directory**/
    string cmd = "mkdir -p " + (string) dir ;
    /**Executing the Shell Command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;
    //creating exposure directory
    /**setting the path for the output Exposure  Directory**/
    sprintf (dir , "%s%s" , moduleoutdir , expframedir) ;
    /**Shell command for creating the Directory**/
    cmd = "mkdir -p " + (string) dir ;
     /**Executing the Shell Command**/
    system (cmd.c_str ()) ;
    LOG(INFO) << endl << dir << " directory created" << endl ;

    LOG(INFO) << endl << "Total number of frames - " << nframes << endl ;
    int status = 0 ;
    float framedata[xsize * ysize] ;
    double integrationtime = 0.0 ;
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    char errstr[512] ;
    int bitpix = FLOAT_IMG ;
    int naxis = 2 ;
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;

    vector<string> vhistorystr ;
    if (history==YES){
        getHistory (vhistorystr) ;
    }
    char infile[FLEN_FILENAME] , outfile[FLEN_FILENAME] ;
    unsigned short frameno = 0 ;
    double frametime = 0 ;
    LOG(INFO)<<"Applying FlatField Correction..."<<endl;
    for (int i = 0 ; i < nframes ; i++)
    {
        sprintf (errstr , "Error at iteration number %d" , i) ;
        fitsfile *fptr , *fout ;
        sprintf (infile , "%s/%s/%s" , inputdatadir , sigframedir , sigframelist[i]) ;
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening the input File" , infile) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
        printError (status , "Error in reading the pixels  values of the input File") ;
        readKeywords (infile , 1, 3 , TDOUBLE , "INT_TIME" , &integrationtime ,
                TUSHORT , "FRAMENO" ,& frameno ,
                TDOUBLE , "FRMTIME" , &frametime) ;   

        /**Applying FlatField Correction **/
        for (int pixno = 0 ; pixno < xsize * ysize ; pixno++){
            if(framedata[pixno]!=INVALID_PIX_VALUE)
            framedata[pixno] = framedata[pixno] * finalFlatFielddata[pixno] ;
        }
        //signal frames
        /**Setting the path for the output Signal  frame */
        sprintf (outfile , "%s/%s/%s_t%f_f%d_sig_ff.fits" , moduleoutdir , sigframedir , nameprefix , frametime , frameno) ;
        /**Creating the output Signal Frame**/
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error creating the output Signal file" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error creating the image Signal file" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata , &status) ;
        printError (status , "Error in writing the pixels to output file" , outfile) ;
        copyUserKeywords (fptr , fout) ;
        if(history==YES)  writeHistory (outfile , vhistorystr) ;
        writeCommonKeywords (fout , modulename) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output Signal File" , outfile) ;
        fits_close_file (fptr , &status) ;
        printError (status , "Error in closing the  input f Signal file" , outfile) ;
        
        strcpy (outsigframelist[i] , basename (outfile)) ; //taking only filename

        /**Setting the path for the exposure Frame directory**/
        sprintf (infile , "%s/%s/%s" , inputdatadir , expframedir , expoframelist[i]) ;
        /**opening the input Exposure File**/
        fits_open_file (&fptr , infile , READONLY , &status) ;
        printError (status , "Error in opening  the input Exposure File" , infile) ;
        fits_read_pix (fptr , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
        printError (status , "Error in reading the pixels for input Exposure file" , infile) ;
        readKeywords (infile , 1, 3 , TDOUBLE , "INT_TIME" , &integrationtime ,
                TUSHORT , "FRAMENO" ,& frameno ,
                TDOUBLE , "FRMTIME" , &frametime) ;   
        sprintf (outfile , "%s/%s/%s_t%f_f%d_exp_ff.fits" , moduleoutdir , expframedir , nameprefix , frametime , frameno) ;
        fits_create_file (&fout , outfile , &status) ;
        printError (status , "Error in creating the output Exposure File" , outfile) ;
        fits_create_img (fout , bitpix , naxis , naxes , &status) ;
        printError (status , "Error in creating the image in output Exposure File" , outfile) ;
        fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata , &status) ;
        printError (status , "Error in writing the pixels to the output Exposure File" , outfile) ;
        copyUserKeywords (fptr , fout) ;
         if(history==YES){
               writeHistory (outfile , vhistorystr) ; //write history to each file
        }
        writeCommonKeywords (fout , modulename) ;
        fits_close_file (fout , &status) ;
        printError (status , "Error in closing the output File" , outfile) ;

       
        strcpy (outexpframelist[i] , basename (outfile)) ; //taking only filename

        fits_close_file (fptr , &status) ;
        printError (status , "Error in Closing the  input File" , infile) ;
      cout<< "Total Files written=" << i << "     Remaining files= " << nframes - i <<" \r" ;
    } //end of loop for number of frames
      writeUsrkeywordsFrmvect (infofile_out,key_records);
     if(history==YES)
                writeHistory (infofile_out , vhistorystr) ;
      vhistorystr.clear () ;
    /*----writing generated frame list to output information file----*/
/**Open output information file  which has been generated **/
    LOG(INFO)<<"Updating the output information file..."<<endl;
   status=writeColumnsToFITS (infofile_out,2,2,TSTRING,1,outsigframelist,nframes,TSTRING ,2,outexpframelist,nframes);
   if(status)
   {
       LOG(INFO)<<"Error writing  the columns to  the fits  file"<<endl;
        return(EXIT_FAILURE);
   }
    LOG(INFO) << " \nInformation File of FlatField module  Created"<< endl ;
    /**Releasing the memory of list of output file**/
    freeMemory (outsigframelist , nframes , NAMESIZE) ;
    freeMemory (outexpframelist , nframes , NAMESIZE) ;
    return (EXIT_SUCCESS) ;
}

int uvtFlatFieldCorr::readFlatFieldFile ()
{
   fitsfile *flatfield ;
    long naxes[2];
    int status = 0 ;
    fits_open_file (&flatfield , flatfieldfile , READONLY , &status) ;
    printError (status , "Error in opening flatfield file : " , flatfieldfile) ;
    //   int hdutype;
    int dim =2;
    fits_get_img_size(flatfield,dim,naxes,&status);
   
    flatfield_dim1=naxes[0];
    flatfield_dim2=naxes[1];
    
    flatfielddata = new float[flatfield_dim1 * flatfield_dim2] ;                        //in case of PC mode;
   
    for (int s = 0 ; s < flatfield_dim1*flatfield_dim2 ; s++)
        flatfielddata[s] = 0.0 ;
    
    long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    fits_read_pix (flatfield , TFLOAT , fpixel ,flatfield_dim1*flatfield_dim2 , NULL , flatfielddata , NULL , &status) ;
    printError (status , "Error in reading the pixels from the caldb FlatfieldFile" , flatfieldfile) ;
  
    fits_close_file (flatfield , &status) ;
    printError (status , "Error in  closing the caldb FlatField File" , flatfieldfile) ;
    LOG(INFO) << "\n----------Reading of caldb  File is Completed----------" << endl ;
    return (EXIT_SUCCESS) ;
}

int uvtFlatFieldCorr::getHistory (vector<string> &vhistory)
{
    int cnt=0;
    char *user = getlogin () ;
    string str = "Module run by " + (string) user ;
    vhistory.push_back (str) ;
    vhistory.push_back ("Parameter List START for " + (string) modulename) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"P1 framelistDir=" + (string) inputdatadir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"P2 Flat Field file " + (string) flatfieldfile) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"P3 outdir=" + (string) outdir) ;
    vhistory.push_back ((string)getSerialNo (cnt)+"Module Output directory=" + (string) moduleoutdir) ;
    if (clobber == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+"P4 clobber=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+"P4 clobber=no") ;
    if (history == YES)
        vhistory.push_back ((string)getSerialNo (cnt)+"P5 history=yes") ;
    else
        vhistory.push_back ((string)getSerialNo (cnt)+"P5 history=no") ;
    vhistory.push_back ("Parameter List END") ;
    return (EXIT_SUCCESS) ;
}

