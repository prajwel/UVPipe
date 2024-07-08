
/* 
 * File:   uvtDetectDistCorrL2.cpp
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
#include "uvtDetectDistCorrL2.h"
#include<pthread.h>
#include<uvtUtils.h>
#include<glog/logging.h>

//#include "uvtFindWtdMean.h"



//Constructor -called when object is createdr


uvtDetectDistCorrl2::uvtDetectDistCorrl2() {
    sprintf(modulename, "%s_%s", MODULENAME, VERSION);
    strcpy(centroidDir, "Centroid");
}

//Destructor

uvtDetectDistCorrl2::~uvtDetectDistCorrl2() {
//    delete[] x_Distortion;
//    delete[] y_Distortion;

}

//parameter File reading

int uvtDetectDistCorrl2::read(int argc, char** argv) {
    int status = 0;


    if (PIL_OK != (status = PILInit(argc, argv))) {

        LOG(ERROR) << "***Error Initializing PIL***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("framelistDir", inputdatadir))) {
        LOG(ERROR) << endl << "***Error reading framelist file name***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("caldbDir", caldbDir))) {
        LOG(ERROR) << endl << "***Error reading Caldb  Directory Path ***";
        return status;
    }


    if (PIL_OK != (status = PILGetFname("outdir", outdir))) {
        LOG(ERROR) << endl << "***Error reading output directory name***";
        return status;
    }
    if (PIL_OK != (status = PILGetBool("clobber", &clobber))) {
        LOG(ERROR) << "***Error Reading clobber:" << clobber << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetBool("history", &history))) {
        LOG(ERROR) << "***Error Reading history parameter:" << history << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetString("mode", mode))) {
        LOG(ERROR) << "***Error Reading mode parameter:" << history << "***";
        return status;
    }
    PILClose(status);
    return (EXIT_SUCCESS);
}

int uvtDetectDistCorrl2::read(char *inputdatadir, char * caldbDir, char *outdir, int clobber, int history) {
    strcpy(this->inputdatadir, inputdatadir);
    strcpy(this->outdir, outdir);
    strcpy(this->caldbDir, caldbDir);
    this->clobber = clobber;
    this->history = history;
    return (EXIT_SUCCESS);
}
//Parameter file content Display

void uvtDetectDistCorrl2::display() {

    LOG(INFO) << endl << "----------Detect Distortion Parameter  Display---------";
    LOG(INFO) << endl << "Input Frame List Directory                        : " << inputdatadir;
    LOG(INFO) << endl << "Caldb Directory                                           : " << caldbDir;
    LOG(INFO) << endl << "Output Directory                               : " << outdir;
    if (clobber == YES)
        LOG(INFO) << endl << "Overwrite                                         : YES";
    else
        LOG(INFO) << endl << "Overwrite                                         : NO";
    if (history == YES)
        LOG(INFO) << endl << "History                                             : YES";
    else
        LOG(INFO) << endl << "History                                              : NO";

    LOG(INFO) << endl << "-------Detect Distortion Parameters Display Ends-------\n";

}

//Correction for the  Cosmic Ray process

int uvtDetectDistCorrl2::uvtDetectDestroProcess() {
    /*if the Input Directory of the Filelist is not available then exit from module
     else read Filelist

     */
    LOG(INFO) << endl << "Detect Destortion process started" << endl;
    sprintf(moduleoutdir, "%s/%s/", outdir, modulename);
    //LOG(INFO)<<endl<<"Module Output Directory : "<<moduleoutdir<<endl;
    string cmd;
    if (DirExists(moduleoutdir) && clobber == YES) {
        LOG(INFO) << endl << "Directory exists and clobber=yes";
        cmd = (string) "rm -rf " + (string) moduleoutdir;
        system(cmd.c_str());
    } else if (DirExists(moduleoutdir) && clobber == NO) {
        LOG(INFO) << endl << moduleoutdir << "  already exists ";
        LOG(INFO) << endl << "Use clobber=yes for overwriting";
        return (EXIT_FAILURE);
    }

    nframes = 0;
    cmd = "mkdir -p " + (string) moduleoutdir;
    system(cmd.c_str()); // creating output directory to keep output from unitConversion
    LOG(INFO) << endl << moduleoutdir << "  directory created";

    string  tempfilepath = searchFile(inputdatadir, ".info");
     if (tempfilepath ==" ")
    {
        LOG (ERROR) << endl << "***Error in finding info file***" ;
        return (EXIT_FAILURE) ;
    }
    sprintf(infofile_in, "%s/%s", inputdatadir, tempfilepath.c_str());

    /*check whether the given Input Directory contains the the list file or not if not then exit from module*/
    if (!(FileExists(infofile_in))) {
        LOG(ERROR) << endl << "Input FileList not Found at Specified PATH,Check INPUT DIRECTORY";
        return (EXIT_FAILURE);
    }
    LOG(INFO) << endl << "\nInformation File :" << infofile_in;
    int status = 0;
    fitsfile *finfo_in, *finfo_out;
    fits_open_file(&finfo_in, infofile_in, READONLY, &status);
    printError(status, "***Error in opening the information file***");
    fits_movabs_hdu(finfo_in, 2, NULL, &status);
    printError(status, "***Error in Moving the 2nd HDU***");
    datainfo.getInfo (finfo_in); //reading basic information for data from information file
    xsize = datainfo.getXsize();
    ysize = datainfo.getYsize();
    fits_read_key(finfo_in, TSTRING, "NAMEPRFX", nameprefix, NULL, &status);
    printError(status, "***NAMEPRFX keyword not Found***"); //for creating name for output information file


    sprintf(infofile_out, "%s/%s_dd.info", moduleoutdir, nameprefix);
    fits_create_file(&finfo_out, infofile_out, &status);
    printError(status, "***Error in creating the output information file***");
    char *ttype[] = {"SignalFrames","ExposureFrames"};
    char *tform[] = {"A256","A256"};
    fits_create_tbl(finfo_out, ASCII_TBL, 0, 2, ttype, tform, NULL, "FileList", &status);
    printError(status, "***Error in creating the table***");
    datainfo.write(finfo_out); //writing basic data information

    fits_update_key(finfo_out, TSTRING, "NAMEPRFX", nameprefix, "File name prefix", &status);
    printError(status, "***Error in updating the key value of NAMEPRFX*** "); //for creating name for output information file
    /*----info file creating completed, rest of the information will be put by other functions-----------*/
    fits_close_file(finfo_out, &status);
    string tempname = caldb_handler.getDetectorFile(datainfo.getDetector(), caldbDir);
    if (tempname == " ") {
        LOG(INFO) << endl << "Couldn't find bad pixel file from caldb" << endl;
        return (EXIT_FAILURE);
    }

    joinStrings(distortionCorrfile, 2, caldbDir, tempname.c_str());
    LOG(INFO) << endl << "Destortion correction file is " << distortionCorrfile;

    x_Distortion= new float[CALDB_DIST_SIZE*CALDB_DIST_SIZE];
    y_Distortion= new float[CALDB_DIST_SIZE*CALDB_DIST_SIZE];
   status = caldb_handler.readCaldbDistFile(x_Distortion, y_Distortion, distortionCorrfile) ; //Distortion File Reading From caldb.
    if (status)
    {
        LOG(ERROR) << "***Error in reading the Distortion file***" << endl ;
        return (EXIT_FAILURE) ;
    }
    
    caldb_xdist_final=new float[xsize*ysize];
    caldb_ydist_final=new float[xsize*ysize];
    for (int p = 0; p < caldbfinalsize * caldbfinalsize; p++) {
        caldb_xdist_final[p] = 0.0;
        caldb_ydist_final[p] = 0.0;
    }
    /*padding process starts here...padding is applied because 512*512 scheme must be converted into 600*600 scheme***/
//LOG(INFO)<<"The Data info flag is "<<datainfo.getModeFlag()<<endl;
  status = Applypadding(x_Distortion,CALDB_DIST_SIZE,CALDB_DIST_SIZE,caldb_xdist_final,xsize,ysize);
    if (status) {
        LOG(ERROR) << "***Padding unsuccessfull for x-distortion***" << endl;
        return (EXIT_FAILURE);
    }
  
  
  //LOG(INFO)<<"y_dist::"<<y_Distortion[0]<<" "<<y_Distortion[1]<<" "<<y_Distortion[2]<<" "<<y_Distortion[512]<<" "<<y_Distortion[1024]<<" "<<endl;
 status =Applypadding(y_Distortion,CALDB_DIST_SIZE,CALDB_DIST_SIZE, caldb_ydist_final,xsize,ysize);
    if (status) {
        LOG(ERROR) << "***Padding Unsuccessfull for y-distortion***" << endl;
        return (EXIT_FAILURE);
    }
 

    if (datainfo.getModeFlag() == 1) {//For IM mode

        fits_read_key(finfo_in, TSTRING, "SIGDIR", sigframedir, NULL, &status);
        printError(status, "***Error in reading the key value of the SIGDIR***");
        fits_read_key(finfo_in, TSTRING, "EXPDIR", expframedir, NULL, &status);
        printError(status, "***Error in reading the key value of the EXPDIR***");
        fits_read_key(finfo_in, TINT, "NFILES", &nframes, NULL, &status);
        printError(status, "***Error in reading the  key value of the NFILES ***");

        LOG(INFO) << endl << "Number of files :" << nframes;
        //reading frame names from information file into vector

        sigframelist = allocateMemory<char>(nframes, NAMESIZE);
        expoframelist = allocateMemory<char>(nframes, NAMESIZE);
        fits_read_col(finfo_in, TSTRING, 1, 1, 1, nframes, NULL, (void *) sigframelist, NULL, &status);
        printError(status, "***Error in reading the reading the signal frame list ***");
        fits_read_col(finfo_in, TSTRING, 2, 1, 1, nframes, NULL, (void *) expoframelist, NULL, &status);
        printError(status, "***Error in reading the reading the Exposure frame list ***");

        //Performing Pixel Padding
        if (detectDestrotionIM()) return (EXIT_FAILURE);

        fits_open_file(&finfo_out, infofile_out, READWRITE, &status);
        printError(status, "***Error in opening the output infofile***");
        fits_movabs_hdu(finfo_out, 2, NULL, &status);
        printError(status, "***Error in moving the 2nd HDU in out information file***");
        fits_update_key(finfo_out, TSTRING, "NAMEPRFX", nameprefix, "File name prefix", &status);
        printError(status, "4"); //for creating name for output information file
        fits_update_key(finfo_out, TSTRING, "SIGDIR", sigframedir, NULL, &status);
        printError(status, "***Error in updating the key value of the SIGDIR***");
        fits_update_key(finfo_out, TSTRING, "EXPDIR", expframedir, NULL, &status);
        printError(status, "***Error updating the key value of the EXPDIR***");
        fits_update_key(finfo_out, TINT, "NFILES", &nframes, NULL, &status);
        printError(status, "***Error in updating the key value of the NFILES ***");
        //        LOG(INFO)<<"after the PixPAdding"<<endl;
        //       exit(1);
        fits_close_file(finfo_out, &status);
        printError(status, "***Error in closing the out information  file***");
        freeMemory(sigframelist, nframes, NAMESIZE);
        freeMemory(expoframelist, nframes, NAMESIZE);



    }//in case of the PC
    else if (datainfo.getModeFlag() == 0) {
        fits_read_key(finfo_in, TSTRING, "EVTFILE", eventfile, NULL, &status);
        printError(status, "***Error in  reading the input information file***");
        //fits_read_key(finfo_in,TSTRING,"IMGFILE",imgfile,NULL,&status);           printError(status,"");
       if (detectDestrotionPC()) return (EXIT_FAILURE);

    }//else neither PC or IM(i.e invalid mode )
    else {
        LOG(ERROR) << endl << "Invalid input for operating mode parameter";
        LOG(ERROR) << endl << "Allowed values are pc/PC/im/IM";
        return (EXIT_FAILURE);
    }
        fits_close_file(finfo_in,&status);
        printError(status, "***Error in  reading the input information file***",infofile_in);
    return (EXIT_SUCCESS);

}

int uvtDetectDistCorrl2::detectDestrotionIM() {

    LOG(INFO) << endl << "\nCorrection of Detector  Distortion  for IM mode";
    char **outsignalframelist = allocateMemory<char >(nframes, NAMESIZE); // to store output frame list
      char **outexposureframelist = allocateMemory<char >(nframes, NAMESIZE); // to store output frame list
    int status = 0;
    fitsfile *finfo_out;

    char dir[FLEN_FILENAME];
    sprintf(dir, "%s%s", moduleoutdir, sigframedir);

    string cmd = "mkdir -p " + (string) dir;

    system(cmd.c_str());
   
    LOG(INFO) << endl << dir << " directory created";
    sprintf(dir, "%s/%s", moduleoutdir, expframedir);
 
    cmd = "mkdir -p " + (string) dir;

    system(cmd.c_str());
 
   
    LOG(INFO) << endl << dir << " directory created ";
  
    long fpixel[2];
    fpixel[0] = 1, fpixel[1] = 1;

    char errstr[512];

    long naxes[2];
    naxes[0] = naxes[1] = xsize;
    vector<string> vhistorystr;
    getHistory(vhistorystr);
    char infile[FLEN_FILENAME], outfile[FLEN_FILENAME];

    unsigned short frameno = 0;
    double frametime = 0;
    float frame_data[xsize*ysize];
    float frame_data_temp[xsize*ysize];
    float frame_data_exp[xsize*ysize];
    float frame_data_temp_exp[xsize*ysize];
    int bitpix = FLOAT_IMG;
    int naxis = 2;
    //ofstream of1("one.txt");
    //nframes To be Processed
    for (int i = 0; i < nframes; i++) 
    {
        sprintf(errstr, "Error at iteration number %d", i);
        fitsfile *fptr,*fout;
        sprintf(infile, "%s/%s/%s", inputdatadir,sigframedir , sigframelist[i]);
         fits_open_file(&fptr, infile, READONLY, &status);
         printError(status, "***input File reading Fails***");
         fits_read_pix(fptr, TFLOAT, fpixel, xsize*ysize, NULL, frame_data, NULL, &status);
         printError(status, "***Error in reading pixels from input frame***");
        fits_read_key(fptr, TUSHORT, "FRAMENO", &frameno, NULL, &status);
        printError(status, "***FRAMENO keyword not Found***");
        fits_read_key(fptr, TDOUBLE, "FRMTIME", &frametime, NULL, &status);
        printError(status, "***FRMTIME keyword not Found***");
        //LOG(INFO)<<xsize*ysize<<endl;exit(1);
        
        for (int i=0;i<xsize*ysize;i++){
            frame_data_temp[i]=INVALID_PIX_VALUE;
            frame_data_temp_exp[i]=INVALID_PIX_VALUE;
        }
         for (int i = 0; i < xsize*ysize; i++)
         {   
            
          int rown=i/xsize;
           int coln=i%xsize;

           rown=rown+caldb_ydist_final[rown*xsize+coln];
           coln=coln+caldb_xdist_final[rown*xsize+coln];
         //  of1<<rown<<" "<<coln<<endl;
           frame_data_temp[(int)(round((rown)*xsize+coln))]=INVALID_PIX_VALUE;
           
           if((int)(round(rown)*xsize+coln)<(xsize*ysize))
           {
           frame_data_temp[(int)(round(rown)*xsize+coln)]=frame_data[i];
           }

        }
         /*corrected Centroid written to the output File.*/
        sprintf(outfile, "%s/%s/%s_t%.4f_f%d_sig_dd.fits", moduleoutdir, sigframedir, nameprefix, frametime, frameno);
        fits_create_file(&fout, outfile, &status);
        printError(status, "Error in creating the level2 dd file");
        fits_create_img(fout, bitpix, naxis, naxes, &status);
        printError(status, "***Error in creating the img***");
        fits_write_pix(fout, TFLOAT, fpixel,xsize*ysize, frame_data_temp, &status);
        printError(status, "***Error in writing the pixels to output signal file***");
        strcpy(outsignalframelist[i], basename(outfile));
        copyUserKeywords(fptr, fout);
        writeCommonKeywords(fout, modulename);
        fits_close_file(fout, &status);
        fits_close_file(fptr, &status);
        sprintf(infile, "%s/%s/%s", inputdatadir,expframedir , expoframelist[i]);
         
        
        
        fits_open_file(&fptr, infile, READONLY, &status);
        printError(status, "***input File reading Fails***");
        fits_read_pix(fptr, TFLOAT, fpixel, xsize*ysize, NULL, frame_data_exp, NULL, &status);
        printError(status, "***Error in reading pixels from input frame***");
        fits_read_key(fptr, TUSHORT, "FRAMENO", &frameno, NULL, &status);
        printError(status, "***FRAMENO keyword not Found***");
        fits_read_key(fptr, TDOUBLE, "FRMTIME", &frametime, NULL, &status);
        sprintf(outfile, "%s/%s/%s_t%.4f_f%d_exp_dd.fits", moduleoutdir, expframedir, nameprefix, frametime, frameno);
      
        for (int i = 0; i < xsize*ysize; i++) {
           
              
           int rown=i/xsize;
           int coln=i%xsize;

           rown=rown+caldb_ydist_final[(rown)*xsize+(coln)];
           coln=coln+caldb_xdist_final[(rown)*xsize+(coln)];
 frame_data_temp_exp[(int)(round(rown)*xsize+coln)]=INVALID_PIX_VALUE;
           if((int)(round(rown)*xsize+coln)<(xsize*ysize)){
               
           frame_data_temp_exp[(int)(round(rown)*xsize+coln)]=frame_data_exp[i];
                   
           
           }

        }
        fits_create_file(&fout, outfile, &status);
        printError(status, "Error creating the output file ",outfile);
        fits_create_img(fout, bitpix, naxis, naxes, &status);
        printError(status, "Error in creating the img ",outfile);
        fits_write_pix(fout, TFLOAT, fpixel,xsize*ysize, frame_data_temp_exp, &status);
        printError(status, "Error in writing the pixels to output signal file ",outfile);
        copyUserKeywords(fptr, fout);
        writeCommonKeywords(fout, modulename);
        fits_close_file(fout, &status);
        fits_close_file(fptr, &status);
        
        if (history == YES) writeHistory(outfile, vhistorystr);
        strcpy(outexposureframelist[i], basename(outfile));
       
    }
   // of1.close();
    /*centroid File list written to the INFO file*/
    fits_open_file(&finfo_out, infofile_out, READWRITE, &status);
    fits_movabs_hdu(finfo_out, 2, NULL, &status);
    printError(status, "***Moving to 2nd HDU in out Info File Fails***");
    fits_write_col(finfo_out, TSTRING, 1, 1, 1, nframes, (void *) outsignalframelist, &status);
    printError(status, "***Writing  centroid List for output INFO File Fails***");
      fits_write_col(finfo_out, TSTRING, 2, 1, 1, nframes, (void *) outexposureframelist, &status);
    printError(status, "***Writing  centroid List for output INFO File Fails***");
    fits_close_file(finfo_out, &status);
    return (EXIT_SUCCESS);
}

int uvtDetectDistCorrl2::detectDestrotionPC() 
{
    char file_in[NAMESIZE], file_out[NAMESIZE]; //for input and output event file
    //opening eventfile

    sprintf(file_in, "%s/%s", inputdatadir, eventfile); //taking event file full path
    LOG(INFO) << eventfile << endl;
    sprintf(file_out, "%s/%s_dd.events", moduleoutdir, nameprefix);
    fitsfile *fevt_in, *fevt_out, *fimg, *fimg_out, *finfo_out;
    int status = 0;

    fits_open_file(&fevt_in, file_in, READWRITE, &status);
    printError(status, "");

    fits_create_file(&fevt_out, file_out, &status);
    printError(status, "");


    /***Difining vars for reading a event file columns***/
    long *frame_no;
    unsigned short *xi, *yi, *mc, *b, *mult, *exp, *max;
    double *t;
    float *xFrac, *yFrac;
    long nrows;

    for (int hdu_no = 2; hdu_no < 3; hdu_no++) {
        fits_movabs_hdu(fevt_in, 2, NULL, &status);
        printError(status, "***Error in  moving to perticuler HDU");
        fits_get_num_rows(fevt_in, &nrows, &status);

        /***Assigning a memory to Arrays***/

        frame_no = new long[nrows];
        xi = new unsigned short[nrows];
        yi = new unsigned short[nrows];
        max = new unsigned short[nrows];
        mc = new unsigned short[nrows];
        b = new unsigned short[nrows];
        mult = new unsigned short[nrows];
        exp = new unsigned short[nrows];
        t = new double[nrows];
        xFrac = new float[nrows];
        yFrac = new float[nrows];




        fits_read_col(fevt_in, TLONG, 1, 1, 1, nrows, NULL, frame_no, NULL, &status);
        printError(status, "***Error reading Frame number ***");
        fits_read_col(fevt_in, TDOUBLE, 2, 1, 1, nrows, NULL, t, NULL, &status);
        printError(status, "***Error reading time***");
        fits_read_col(fevt_in, TUSHORT, 3, 1, 1, nrows, NULL, xi, NULL, &status);
        printError(status, "***Error reading xinteger***");
        fits_read_col(fevt_in, TUSHORT, 4, 1, 1, nrows, NULL, yi, NULL, &status);
        printError(status, "Error reading Yinteger");
        fits_read_col(fevt_in, TFLOAT, 5, 1, 1, nrows, NULL, xFrac, NULL, &status);
        printError(status, "Error reading xfractional");
        fits_read_col(fevt_in, TFLOAT, 6, 1, 1, nrows, NULL, yFrac, NULL, &status);
        printError(status, "Error reading yfractional");
        fits_read_col(fevt_in, TUSHORT, 7, 1, 1, nrows, NULL, max, NULL, &status);
        printError(status, "Error reading maximum ");
        fits_read_col(fevt_in, TUSHORT, 8, 1, 1, nrows, NULL, mc, NULL, &status);
        printError(status, "Error reading minimum corner value");
        fits_read_col(fevt_in, TUSHORT, 9, 1, 1, nrows, NULL, exp, NULL, &status);
        printError(status, "Error reading exposure");
        fits_read_col(fevt_in, TUSHORT, 10, 1, 1, nrows, NULL, b, NULL, &status);
        printError(status, "Error reading Bad flag");
        fits_read_col(fevt_in, TUSHORT, 11, 1, 1, nrows, NULL, mult, NULL, &status);
        printError(status, "***Error reading multiple photon events***");


        float locate;
        for (int k = 0; k < nrows; k++) {
            locate = xi[k]*600 + yi[k];
            ceil(locate);

            xi[k] = xi[k] + caldb_xdist_final[(int) locate];
            yi[k] = yi[k] + caldb_ydist_final[(int) locate];

        }
        char *ttype2[11];
        char *tform[11];
        char *tunit[11];
        char *ext = "EVENT_DATA";
        ttype2[0] = "FrameNo";
        tform[0] = "1J";
        tunit[0] = "";
        ttype2[1] = "Time";
        tform[1] = "1D";
        tunit[1] = "";
        ttype2[2] = "X_INT";
        tform[2] = "1U";
        tunit[2] = "";
        ttype2[3] = "Y_INT";
        tform[3] = "1U";
        tunit[3] = "";
        ttype2[4] = "X_FRAC";
        tform[4] = "1E";
        tunit[4] = "";
        ttype2[5] = "Y_FRAC";
        tform[5] = "1E";
        tunit[5] = "";
        ttype2[6] = "Max-Min";
        tform[6] = "1U";
        tunit[6] = "";
        ttype2[7] = "MinC";
        tform[7] = "1U";
        tunit[7] = "";
        ttype2[8] = "BAD_FLAG";
        tform[8] = "1U";
        tunit[8] = "";
        ttype2[9] = "MULTIPLE_PHOTON_FLAG";
        tform[9] = "1U";
        tunit[9] = "";
        ttype2[10] = "Exposure";
        tform[10] = "1U";
        tunit[10] = "";

        /***Writing to output event file***/

        fits_create_tbl(fevt_out, BINARY_TBL, 1, 11, ttype2, tform, tunit, ext, &status);
        fits_movnam_hdu(fevt_out, ANY_HDU, ext, 0, &status);
        fits_write_col(fevt_out, TINT, 1, 1, 1, nrows, frame_no, &status);
        printError(status, "***Error writing Frame number ***");
        fits_write_col(fevt_out, TDOUBLE, 2, 1, 1, nrows, t, &status);
        printError(status, "***Error writing Frame number ***");
        fits_write_col(fevt_out, TUSHORT, 3, 1, 1, nrows, xi, &status);
        printError(status, "***Error writing Xinteger ***");
        fits_write_col(fevt_out, TUSHORT, 4, 1, 1, nrows, yi, &status);
        printError(status, "***Error writing Yinteger ***");
        fits_write_col(fevt_out, TFLOAT, 5, 1, 1, nrows, xFrac, &status);
        printError(status, "***Error writing Xfractional ***");
        fits_write_col(fevt_out, TFLOAT, 6, 1, 1, nrows, yFrac, &status);
        printError(status, "***Error writing Yfractional ***");
        fits_write_col(fevt_out, TUSHORT, 7, 1, 1, nrows, max, &status);
        printError(status, "***Error writing maximum ***");
        fits_write_col(fevt_out, TUSHORT, 8, 1, 1, nrows, mc, &status);
        printError(status, "***Error writing minimum corner value ***");
        fits_write_col(fevt_out, TUSHORT, 9, 1, 1, nrows, b, &status);
        printError(status, "***Error writing bad Flag ***");
        fits_write_col(fevt_out, TUSHORT, 10, 1, 1, nrows, mult, &status);
        printError(status, "***Error writing multiple Photon events ***");
        fits_write_col(fevt_out, TUSHORT, 11, 1, 1, nrows, exp, &status);
        printError(status, "***Error writing Exposure***");





        /*releasing the memory*/

        delete[] xi;
        delete[] yi;
        delete[] xFrac;
        delete[] yFrac;
        delete[] b;
        delete[] exp;
        delete[] max;
        delete[] mc;
        delete[] mult;
        delete[] t;
        delete[] frame_no;
    }

    fits_movabs_hdu(fevt_out, 1, NULL, &status);
    printError(status, "***Moving to 1st HDU fails for output event file*** ");

    /***Writes the  common information to  output information File***/
    writeCommonKeywords(fevt_out, modulename);
    vector<string> vhistorystr;
    getHistory(vhistorystr);

    if (history == YES) writeHistory(file_out, vhistorystr);
    fits_movabs_hdu(fevt_in, 3, NULL, &status);
    printError(status, "***Error in moving to 3rd HDU in input Event File***");
    fits_copy_hdu(fevt_in, fevt_out, NULL, &status);
    printError(status, "***Error in Coping the 3rd HDU to output Event File***");
    fits_close_file(fevt_out, &status);
    printError(status, "***Error in Closing the output Event File***");
    fits_close_file(fevt_in, &status);
    printError(status, "***Error in closing the input Event File***");
    fits_open_file(&finfo_out, infofile_out, READWRITE, &status);
    printError(status, "***Error in opening the output information File***");
    fits_movabs_hdu(finfo_out, 2, NULL, &status);
    printError(status, "***Error in opening the output information File***");
    fits_update_key(finfo_out, TSTRING, "EVTFILE", basename(file_out), NULL, &status);
    printError(status, "***Error in updating the key value for the EVTFILE***");
    // fits_update_key(finfo_out, TSTRING, "EVTFILE", basename(file_out), "File name prefix", &status);

    fits_close_file(finfo_out, &status);
    printError(status, "***Error  closing the output information file***");
    return (EXIT_SUCCESS);
}

int uvtDetectDistCorrl2::getHistory(vector<string> &vhistory) {

    int cnt=0;
    char *user = getlogin();
    string str = "Module run by " + (string) user;
    vhistory.push_back(str);
    vhistory.push_back("Parameter List START for " + (string) modulename);
    vhistory.push_back((string)getSerialNo (cnt)+" framelistDir=" + (string) inputdatadir);
    vhistory.push_back((string)getSerialNo (cnt)+" caldbDir=" + (string) caldbDir);
    vhistory.push_back((string)getSerialNo (cnt)+" outdir=" + (string) outdir);
    vhistory.push_back((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir);

    vhistory.push_back((string)getSerialNo (cnt)+" Centroid Directory=" + (string) centroidDir);
    if (clobber == YES)
        vhistory.push_back((string)getSerialNo (cnt)+" clobber=yes");
    else
        vhistory.push_back((string)getSerialNo (cnt)+" clobber=no");
    if (history == YES)
        vhistory.push_back((string)getSerialNo (cnt)+" history=yes");
    else
        vhistory.push_back((string)getSerialNo (cnt)+" history=no");
    vhistory.push_back("Parameter List END");
    return (EXIT_SUCCESS);
}

//int uvtDetectDistCorrl2::readDestortionFile() {
//     LOG(INFO) << endl << "Reading Distortion Correction file from caldb........" ;
//    fitsfile *fdist ;
//    int status = 0 ;
//    char errstr[500] ;
//
//    fits_open_file (&fdist , distortionCorrfile , READONLY , &status) ;
//    printError (status , "***Error in readDestortionFile()***") ;
//    fits_movabs_hdu (fdist , 2 , NULL , &status) ;
//    printError (status , "***Moving to 2nd HDU in caldb file Fails***") ;
//    fits_get_num_rows (fdist , &caldbdim , &status) ;
//    LOG(INFO) << "The caldbdim is inside the caldb is " << caldbdim << endl ;
//    printError (status , "***Number of Rows counting Fails***") ;
//
//    //check to compare number of rows and columns to frame xsize & ysize
//    
//    x_Distortion = new float[caldbdim * caldbdim] ; //Array For Storing  x_distortion(512*512) From caldb Directory
//    LOG(INFO)<<"The caldb Dimension is "<<caldbdim<<endl;
////    float err[caldbdim] ;
////     float err1[caldbdim] ;
////     for(int i=0;i<caldbdim;i++){
////         err[i]=0.0f;
////         err1[i]=0.0f;
////     }
//
//     status=0;
////     fits_read_col (fdist , TFLOAT , 2 , 1 , 1 , caldbdim , NULL , err1 , NULL , &status) ;
////   printError (status , "***Reading a column Fails in caldb***") ;
////    fits_read_col (fdist , TFLOAT , 1 , 1 , 1 , caldbdim , NULL , err , NULL , &status) ;
////    printError (status , "***Reading a column Fails in caldb***") ;
//    
//    //LOG(INFO)<<"value for the x-dist is "<<err1[1]<<" "<<err[1]<<endl;
//
//    //LOG(INFO)<<"\nValues from caldb file\n";
//    
//  //  ofstream fout("/tmp/distortionfile");
//    
//  float value=0;
//  int index=0;
//    for(int i=0;i<caldbdim;i++){
//        for(int j=0;j<caldbdim;j++){
//            status=0;
//            fits_read_col (fdist , TFLOAT , j+1 , i+1 , 1 ,1 , NULL , &value , NULL , &status) ;
//            printError (status , "***Reading a column Fails in caldb***") ;
//            x_Distortion[index]=value;
//            index++;
//            //  fout<<"\t"<<setprecision(6)<<value;
//         }
//        
//        //fout<<endl;
//    }
////    fout.close();
////    exit(1);
//    
//  index=0;
//    
//    
////    for (int i = 0 ; i < caldbdim ; i++)
////    {
////        for (int j = 0 ; j < caldbdim ; j++)
////        {
////           
////            x_Distortion[i * caldbdim + j] = err[i] ;
////
////        }
////    }
//
//    //LOG(INFO) << "x distr in class " << x_Distortion[539] << endl ;
//
//   
//    fits_movabs_hdu (fdist , 3 , NULL , &status) ;
//    printError (status , "***Moving to 3rd HDU in caldb file Fails***") ; //for distortion in y
//    fits_get_num_rows (fdist , &caldbdim , &status) ;
//    LOG(INFO) << "The caldbdim is inside the caldb is " << caldbdim << endl ;
//    printError (status , "***Number of Rows counting Fails***") ; //getting no of rows in the perticuler HDU
//    y_Distortion = new float[caldbdim * caldbdim] ; //Array For Storing  y_distortion(512*512) From caldb Directory
//
//    
//     for(int i=0;i<caldbdim;i++){
//        for(int j=0;j<caldbdim;j++){
//            status=0;
//            fits_read_col (fdist , TFLOAT , j+1 , i+1 , 1 ,1 , NULL , &value , NULL , &status) ;
//            printError (status , "***Reading a column Fails in caldb***") ;
//            y_Distortion[index]=value;
//            index++;
//            //  fout<<"\t"<<setprecision(6)<<value;
//         }
//        
//        //fout<<endl;
//    }
//    
//    
//    LOG(INFO)<<"x_distr:::"<<x_Distortion[0]<<" "<<x_Distortion[1]<<endl;
//      LOG(INFO)<<"y_distr:::"<<y_Distortion[0]<<" "<<y_Distortion[1]<<endl;
//      //exit(1);
////    fits_read_col (fdist , TFLOAT , 1 , 1 , 1 , caldbdim , NULL , err , NULL , &status) ;
////    printError (status , "***Reading a column Fails in caldb***") ;
////
////    for (int i = 0 ; i < caldbdim ; i++)
////    {
////        for (int j = 0 ; j < caldbdim ; j++)
////        {
////            y_Distortion[i * caldbdim + j] = err[i] ;
////        }
////    }
//
//
//    fits_close_file (fdist , &status) ;
//
//    LOG(INFO) << endl << "Reading Distortion Correction file from caldb completed" ;
//    return (EXIT_SUCCESS) ;
//}

//int uvtDetectDistCorrl2::Applypadding(float *inputArray, float *outputArray) {
//
//    /*input Array is Of (512*512) and outputArray is of (600*600).*/
//
//    LOG(INFO)<<"INSIDE PADDDING "<<endl;
//    int count = 0;
//    int final = caldbfinalsize - caldbdim;
//    LOG(INFO) << final << endl;
//    int eachside = final / 2;
//    int p;
//    LOG(INFO)<<"TTTT"<<endl;
//    /*This Loop substitute ' inputArray' to the middle of  the 'output Array' i.e same Distance from each side of the inputArray to outputArray*/
//
//    for (long r = 0; r <caldbdim; r++) {
//        for (long j = p = (eachside + r) * caldbfinalsize + eachside; j < p + caldbdim; j++) {
//            //LOG(INFO)<<count<<endl;
////            if(j>)
//            outputArray[j] = inputArray[count];
//            count++;
//
//        }
//    }
//    LOG(INFO)<<count<<endl;
//    LOG(INFO)<<"outside of the paddingng"<<endl;
//    return (EXIT_SUCCESS);
//
//}
