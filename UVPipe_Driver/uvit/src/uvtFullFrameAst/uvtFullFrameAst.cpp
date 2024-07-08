/* 
 * File:   uvtFullFrameAst.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#include <uvtFullFrameAst.h>
#include<uvtUtils.h>
//#include <Attitude.h>
#include <glog/logging.h>
#include <caldb_Handler.h>
#include<algorithm>
#include<sqlite3.h>
#include<spMatrix.h>
#include<fstream>
#include<caldb_Handler.h>
#include<uvtUtils.h>
#include<transform.h>
#include <iomanip>
using namespace std;
int writeUsrkeywordsFrmvectNew(char *filename,vector<string> &strvect);
bool compare (struct Star1 vect1 , struct Star1 vect2) ;
bool compare (struct Star1 vect1 , struct Star1 vect2)
{
    return (vect1.intensity > vect2.intensity) ;
}
//Constructor -called when object is created
uvtFullFrameAst:: uvtFullFrameAst() {
    sprintf(modulename, "%s_%s", MODULENAME, VERSION);
   
}

//Destructor
uvtFullFrameAst::~ uvtFullFrameAst() {
}

//parameter reading
int uvtFullFrameAst::read(int argc, char** argv) {
    int status = EXIT_SUCCESS;

    status = readParams(argc,argv,20,
            FNAME,"indir",inputdatadir,
            FNAME,"RADECfilename",RADECfileinput,
            FNAME,"caldbDir",caldbDir,                                //CALDB dir is required for teldef file
            FNAME,"outdir",outdir,
            FNAME,"attitudefile",attitudefile,
            FNAME,"EXPFILENAME",expfile,
            FNAME,"NOICEMAPFILENAME",noiceMapFile,
            FNAME,"ShiftNRotfile",shiftRotfile,
            STRING,"att_timecol",&att_timecol,
            STRING,"att_qcol",&att_qcol,
            FNAME,"catalogpath",catalogpath,
            BOOL,"clobber",&clobber,
            REAL4,"threshold",&sd_multi_factor_default,
            STRING,"Radi_search",rad_search,
            INT,"minimum_targetedstars",&minimum_No_of_Stars,
            INT,"refine_Window",&refine_Winsize,
            INT,"centroid_Window",&centroid_Winsize,
            STRING,"database_name",databasename,
            BOOL,"history",&history,
            STRING,"mode",mode);
    if(status){
        LOG(ERROR)<<"***Error reading parameters***";
        return(EXIT_FAILURE);
    }
  if (PIL_OK != (status = PILInit (argc , argv)))
    {
        LOG (INFO) << "***Error Initializing PIL***" ;
        return status ;
    }
      if (PIL_OK != (status = PILGetInt ("search_algo_forFullFrameAst" , &search_algo_ctlg)))
        {
            LOG (INFO) << endl << "***Error Reading search method" << search_algo_ctlg << "***" ;
            return status ;
        }
    if(search_algo_ctlg==1 || search_algo_ctlg==3 || search_algo_ctlg==5)
    {
      if (PIL_OK != (status = PILGetString ("len_rect_a" , (char*)&len_a)))
        {
            LOG (INFO) << endl << "***Error Reading length of rectangle :" <<len_a << "***" ;
            return status ;
        }
    
      if (PIL_OK != (status = PILGetString ("len_rect_b" , (char*)&len_b)))
        {
            LOG (INFO) << endl << "***Error Reading width of rectangle:" << len_b << "***" ;
            return status ;
        }
    
    }
    PILClose (status) ;
    return status;
}

int uvtFullFrameAst::read(char *inputdatadir, char* RADECfilename,char * caldbDir, char *shiftNRotfile, char *NoiceMapfile,char *exposureFile,char *outdir, char *attitudefile,
                                                char *att_timecol, char *att_qcol,char *catalogpath,float sd_mult_fact,int min_stars,int refine_win,int centroid_win,char *database_name,int algoval,char *len,char *wid, char *rad,int clobber, int history,int flag_Imageinput) {
   
    strcpy(this->inputdatadir, inputdatadir);
    strcpy(this->RADECfileinput,RADECfilename);
    strcpy(this->outdir, outdir);
    strcpy(this->shiftRotfile,shiftNRotfile);
    strcpy(this->expfile,exposureFile);
    strcpy(this->noiceMapFile,NoiceMapfile);
    strcpy(this->caldbDir, caldbDir);
    strcpy(this->attitudefile,attitudefile);
    strcpy(this->att_qcol,att_qcol);
    strcpy(this->att_timecol,att_timecol);
    strcpy(this->catalogpath,catalogpath);
    strcpy(this->databasename,database_name);
    strcpy(this->rad_search, rad);
  
    this->sd_multi_factor_default=sd_mult_fact;
    this->refine_Winsize=refine_win;
    this->centroid_Winsize=centroid_win;
    this->minimum_No_of_Stars=min_stars;
    this->search_algo_ctlg=algoval;
    strcpy(this->len_a,len);
    strcpy(this->len_b,wid);
    this->clobber = clobber;
    this->history = history;
    this->flag_inputImage=flag_Imageinput;
    return (EXIT_SUCCESS);
}

int uvtFullFrameAst::read(char *inputdatadir, char * caldbDir, char *shiftNRotfile, char *outdir, char *attitudefile,
                                                char *att_timecol, char *att_qcol,char *catalogpath,float sd_mult_fact,int min_stars,int refine_win,int centroid_win,char *database_name,int algoval,char *len,char *wid, char *rad,int clobber, int history,int flag_Imageinput) {
   
    strcpy(this->inputdatadir, inputdatadir);
    strcpy(this->outdir, outdir);
    strcpy(this->shiftRotfile,shiftNRotfile);
    //strcpy(this->expfile,exposureFile);
    //strcpy(this->noiceMapFile,NoiceMapfile);
    strcpy(this->caldbDir, caldbDir);
    strcpy(this->attitudefile,attitudefile);
    strcpy(this->att_qcol,att_qcol);
    strcpy(this->att_timecol,att_timecol);
    strcpy(this->catalogpath,catalogpath);
    strcpy(this->databasename,database_name);
    strcpy(this->rad_search, rad);
  
    this->sd_multi_factor_default=sd_mult_fact;
    this->refine_Winsize=refine_win;
    this->centroid_Winsize=centroid_win;
    this->minimum_No_of_Stars=min_stars;
    this->search_algo_ctlg=algoval;
    strcpy(this->len_a,len);
    strcpy(this->len_b,wid);
    this->clobber = clobber;
    this->history = history;
    this->flag_inputImage=flag_Imageinput;
    return (EXIT_SUCCESS);
}



//Parameter file content Display

void uvtFullFrameAst::display() {

    LOG(INFO)<< "----------FULL FRAME ASTROMETRY PARAMETERS---------";
    LOG(INFO)<< "Input Directory   : "<<inputdatadir;
    LOG(INFO)<< "Caldb Directory  : "<<caldbDir;
    LOG(INFO) <<"Attitude File  :  "<<attitudefile;
    LOG(INFO)<< "Time column name   : "<<att_timecol;
    LOG(INFO)<< "Quaternion column name : "<<att_qcol;
    LOG(INFO)<< "Catalog Path  :  "<<catalogpath;
    LOG(INFO)<< "Output Directory : " << outdir;
     LOG(INFO)<< "Minimum Number of star requirement : " <<minimum_No_of_Stars;
    if (clobber == YES)
       LOG(INFO)<< "Overwrite  : YES";
    else
        LOG(INFO)<< "Overwrite : NO";
    if (history == YES)
        LOG(INFO)<<"History  : YES";
    else
        LOG(INFO)<<"History  : NO";

    LOG(INFO) <<"----------------------------------------------------------------------------------\n";

}
int uvtFullFrameAst::uvtFullFrmAstProcess() {
    
    LOG(INFO)<< "Full Frame Astrometry process started"; 
    sprintf(moduleoutdir, "%s/%s/", outdir, modulename);
    
    if(createOutputDirectory(clobber, moduleoutdir))  
        return (EXIT_FAILURE);
 
    LOG(INFO)<<"\033[1;34m"<<inputdatadir<<"\033[0m";
    string tempimagefile=" ";
    tempimagefile = searchFile(inputdatadir, "_radec.fits");
    
     if (tempimagefile==" "){
tempimagefile = searchFile(inputdatadir, "_sig_rg.fits");
    if (tempimagefile==" ")
    {
      tempimagefile = searchFile(inputdatadir, "_sig_snl.fits");
         if (tempimagefile==" "){
             tempimagefile = searchFile(inputdatadir, "_ra-dec.fits");
		if (tempimagefile==" "){
		  tempimagefile = searchFile(inputdatadir, "NOT-CORRECTsigCHECK_X_Y_IMAGEONLY.fits");
             if (tempimagefile==" "){
			tempimagefile = searchFile(inputdatadir, "_Sig.fits");
           

	 if (tempimagefile==" "){
                tempimagefile = searchFile(inputdatadir, ".fits");
                 if (tempimagefile==" "){
             LOG(ERROR)<<"Input Directory does not contain required file"<<endl;
             return(EXIT_FAILURE);
                 }
             }
             
         }
}
    }
}
     }
 
   
    // printError (status , "Error in closing the input File" , noiceMapFile) ;
     LOG(INFO)<<"INIIIII";
    sprintf (imagefile_in, "%s/%s", inputdatadir, tempimagefile.c_str());              //input image file path
    double  RA_pnt,DEC_pnt;
    char date_obs[FLEN_FILENAME],time_obs[FLEN_FILENAME];
//    readKeywords(imagefile_in,2,5,TSTRING,"NAMEPRFX",nameprefix,TINT,"XSIZE",&xsize,TINT,"YSIZE",&ysize,TSTRING,"DATE-OBS",date_obs,TSTRING,"TIME-OBS",time_obs,TFLOAT,"ROLLAPPLIED",&roll_angle);
//    readKeywords(imagefile_in,2,2,TDOUBLE ,"RA_PNT_RADEC",&RA_pnt,TDOUBLE,"DEC_PNT_RADEC",&DEC_pnt);
    readKeywords(RADECfileinput,1,5,TSTRING,"NAMEPRFX",nameprefix,TINT,"XSIZE",&xsize,TINT,"YSIZE",&ysize,TSTRING,"DATE-OBS",date_obs,TSTRING,"TIME-OBS",time_obs,TFLOAT,"ROLLAPPLIED",&roll_angle);
    readKeywords(RADECfileinput,1,2,TDOUBLE ,"RA_PNT_RADEC",&RA_pnt,TDOUBLE,"DEC_PNT_RADEC",&DEC_pnt);
    
    LOG(INFO)<<RA_pnt<<" "<<DEC_pnt;
    rapnt_att=RA_pnt;
    decpnt_att=DEC_pnt;
    center_ra_prev=(float)rapnt_att;
    center_dec_prev=(float)decpnt_att;
    //LOG(INFO)<<date_obs;
    
    sprintf(imagefile_out,"%s/%s_as_Sig.fits",moduleoutdir,nameprefix);                 //ouptut image file path
    sprintf(imagefile_out_expName,"%s/%s_as_Exp.fits",moduleoutdir,nameprefix); 
    sprintf(imagefile_out_noicemapName,"%s/%s_as_NoiseMap.fits",moduleoutdir,nameprefix); 
    //Read Attitude file 
    
    //till this. 
   int status = 0;
  long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
     fitsfile *expPointer;
     float *Expdata = new float[xsize*ysize];
     
    fitsfile *fin_ra;
    fits_open_file(&fin_ra,RADECfileinput,READONLY,&status);
    printError(status, "File cannot be opened",imagefile_in);
     vector<string> tempL1header;
         copyUsrkeywrdsTovect(fin_ra,tempL1header);
    //datainfo.getInfo(fin);
    fits_close_file(fin_ra,&status);
     
     
     //opening exposure file.
     fits_open_file(&expPointer,expfile,READONLY,&status);
     printError(status, "File cannot be opened",expfile);
     fits_read_pix (expPointer , TFLOAT , fpixel , xsize*ysize , NULL , Expdata , NULL , &status) ;
     printError (status , "Error in reading the pixels from the input File" , expfile) ;
      datainfo.getInfo(expPointer);
     //fits_close_file(expPointer,&status);
     //printError (status , "Error in closing the input File" , expfile) ;
     
     //opening the noicemap file.
     fitsfile *noicemapPointer;
     float *noicemapData= new float[xsize*ysize];
     fits_open_file(&noicemapPointer,noiceMapFile,READONLY,&status);
     printError(status, "File cannot be opened",noiceMapFile);
     fits_read_pix (noicemapPointer , TFLOAT , fpixel , xsize*ysize , NULL , noicemapData , NULL , &status) ;
      printError (status , "Error in reading the pixels from the input File" , noiceMapFile) ;
     //fits_close_file(noicemapPointer,&status);
    
     fitsfile  *fin, *fout,*foutexp,*foutnoicemap;    //, *fteldef;
     long nrows_snr;
     float *framedata= new float[xsize*ysize];
      double  x1,x2;
       float *X_loc_snr,*Y_loc_snr,*mult_phtn_snr,*enp_snr,*bad_flag_snr;
     if(flag_inputImage==0){
     fits_open_file (&fin , imagefile_in , READWRITE , &status) ;
        printError (status , " Error in opening file " , imagefile_in) ;
        fits_movabs_hdu (fin , 2 , NULL , &status) ;
        printError (status , "Error in moving to HDU " ,imagefile_in) ;
        fits_get_num_rows (fin , &nrows_snr , &status) ;
        printError (status , "Error in reading the number of rows" ,imagefile_in) ;
      X_loc_snr= new  float[nrows_snr];
        Y_loc_snr= new float[nrows_snr];
        bad_flag_snr = new float[nrows_snr];
        mult_phtn_snr= new float[nrows_snr];
        enp_snr= new float[nrows_snr];
        fits_read_col (fin , TFLOAT , 13 , 1 , 1 , nrows_snr  , NULL , Y_loc_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" ,imagefile_in) ;
        fits_read_col (fin , TFLOAT , 12 , 1 , 1 , nrows_snr  , NULL , X_loc_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" ,imagefile_in) ;
        fits_read_col (fin , TFLOAT , 9 , 1 , 1 , nrows_snr  , NULL , bad_flag_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" ,imagefile_in) ;
        fits_read_col (fin , TFLOAT , 10 , 1 , 1 , nrows_snr  , NULL , mult_phtn_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , imagefile_in) ;
         fits_read_col (fin , TFLOAT , 11 , 1 , 1 , nrows_snr  , NULL , enp_snr , NULL , &status) ;
        printError (status , "Error in reading the number of rows" , imagefile_in) ;
        fits_close_file (fin , &status) ;
        printError (status , "Error closing the file " , imagefile_in) ;
     
       
     
    for (int i = 0; i < nrows_snr; i++) {
        x1 = X_loc_snr[i];
        x2 = Y_loc_snr[i];


        if (round(x1) > 0 && round(x1) < xsize && round(x2) > 0 && round(x2) < ysize)
            framedata[(int) (round(x2) * xsize + round(x1))] = framedata[(int) (round(x2) * xsize + round(x1))] + bad_flag_snr[i] * mult_phtn_snr[i] * enp_snr[i];
        

    }
          for (int i = 0; i < xsize * ysize; i++) {
        if (Expdata[i] != 0.0f)
            framedata[i] = framedata[i] / Expdata[i];
        else
            framedata[i] = 0.0f;
    }
        
     }
     else{

    fits_open_file(&fin,imagefile_in,READONLY,&status);
    printError(status, "File cannot be opened",imagefile_in);
   
    fits_read_pix (fin , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
    printError (status , "Error in reading the pixels from the input File" , imagefile_in) ;
    
    fits_close_file (fin , &status) ;
    printError (status , "Error closing the file " , imagefile_in) ;
         
     }
        
        
        
     
     
     
     
//    fits_open_file(&fin,imagefile_in,READONLY,&status);
//    printError(status, "File cannot be opened",imagefile_in);
//   
//    fits_read_pix (fin , TFLOAT , fpixel , xsize*ysize , NULL , framedata , NULL , &status) ;
//     printError (status , "Error in reading the pixels from the input File" , imagefile_in) ;
    // xsize=ysize=600;
    sd_mul_factor = sd_multi_factor_default ;
    refine_Winsize=refine_Winsize*xsize/PIX_PADD_SIZE;
    centroid_Winsize=centroid_Winsize*xsize/PIX_PADD_SIZE;

    
  float Max_value_Exp=Expdata[(xsize/2)*xsize+(xsize/2)];
 
   status=findStar_algo1 (framedata,Expdata,Max_value_Exp);
   if (status)
    {
        LOG (ERROR) << "Error in star detection algorithm" ;
	LOG(ERROR)<<"CRASH: NO STARS FOUND IN UV IMAGE (uvtFullFraneAst.cpp)";
        return (EXIT_FAILURE) ;
    }
   
   //added
   LOG(INFO)<<Cx.size ();
  
   

   
   
   
   //opening shift and Rotation file
   vector<double> cx_ref,cy_ref,ci_ref;
   int indexx,indexy,cntpixels;
float Sum_ofpixels,valuetoCmpr;
   for (int i = 0; i < Cx.size(); i++) {
//LOG(INFO)<<uvitObj.Cx[i]<<" "<<uvitObj.Cy[i]<<" "<<Frame_Exposure_data[(int)(round(uvitObj.Cy[i])*uvitObj.IMG_DIM_FIpc+round(uvitObj.Cx[i]))];
cntpixels=0;Sum_ofpixels=0.0f;
indexx=round(Cx[i]);
indexy=round(Cy[i]);
	for (int j=indexx-12;j<=indexx+12;j++){
for (int k =indexy-12;k<=indexy+12;k++){
if((int)(round(k)*xsize+round(j)) >=0 && (int)(round(k)*xsize+round(j))<xsize*ysize && Expdata[(int)(round(k)*xsize+round(j))]!= INVALID_PIX_VALUE ){
Sum_ofpixels=Sum_ofpixels+Expdata[(int)(round(k)*xsize+round(j))];
cntpixels++;
}
}
}
valuetoCmpr=Sum_ofpixels/cntpixels;
//LOG(INFO)<<"TT:"<<Max_value_Exp<<" "<<valuetoCmpr<<" "<<cntpixels;
if(valuetoCmpr>(20*Max_value_Exp/100)){
cx_ref.push_back(Cx[i]);
                                    cy_ref.push_back(Cy[i]);
                                    ci_ref.push_back(Ci[i]);
}

}

Cx.clear();
Cy.clear();
Ci.clear();

Cx=cx_ref;
Cy=cy_ref;
Ci=ci_ref;
LOG(INFO)<<Cx.size();
   
 
   if(RA_pnt == -9999 || DEC_pnt== -9999)
   {
       LOG(ERROR)<<"***NOT found RA_PNT and DEC_PNT direction in header***";
       LOG(ERROR)<<"CRASH NO TIME OVERLAP OF UV IMAGE WITH ATTUTUDE FILE (uvtFullFrameAst.cpp)";
       return(EXIT_FAILURE);
   }
   
   LOG(INFO)<<"TWIST angle "<<setprecision(20)<<" "<<RA_pnt<<" "<<DEC_pnt;
   
      double cdelt1,cdelt2;
      float factor_delta=4800/600;
      if(strcmp(datainfo.getDetector (),"VIS")==0)
    {
        cdelt1=(XSCALE_ARCSEC_VIS/3600)/factor_delta;
        cdelt2=(YSCALE_ARCSEC_VIS/3600)/factor_delta;
    }
    else if (strcmp(datainfo.getDetector (),"FUV")==0)
    {
       cdelt1=(XSCALE_ARCSEC_FUV/3600)/factor_delta;
       cdelt2=(YSCALE_ARCSEC_FUV/3600)/factor_delta;
    }
    else if(strcmp(datainfo.getDetector (),"NUV")==0)
    {
       cdelt1=(XSCALE_ARCSEC_NUV/3600)/factor_delta;
       cdelt2=(YSCALE_ARCSEC_NUV/3600)/factor_delta;       
    }
  
     cdelt1=-cdelt1;///cos(center_dec*M_PI/180);
  
     
      fits_create_file(&foutexp, imagefile_out_expName, &status);
      printError(status," ",imagefile_out);
     
      fits_create_file(&foutnoicemap, imagefile_out_noicemapName, &status);
       printError(status," ",imagefile_out);
      
      fits_create_file(&fout, imagefile_out, &status);
      printError(status," ",imagefile_out);
     
   
     
        long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    int naxis = 2 ;
    int bitpix = FLOAT_IMG ;
        // fits_create_img (fout , bitpix , naxis , naxes , &status) ;
         //fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , Rotated_frmData , &status) ;
        // fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata , &status) ;
        // fits_close_file(fout,&status);
        // exit(1);
    // LOG(INFO)<<date_obs;exit(1);
     string Date_Of_Obs=date_obs;
     string Time_Of_Obs=time_obs;
    string year_ofOBS=Date_Of_Obs.substr(0,4);
    string  month_ofOBS =Date_Of_Obs.substr(5,2);
    string date_ofOBS = Date_Of_Obs.substr(8,2);
    string time_ofOBS_hh=Time_Of_Obs.substr(0,2);
    string time_ofOBS_mm=Time_Of_Obs.substr(3,2);
    string time_ofOBS_ss=Time_Of_Obs.substr(6,2);
    LOG(INFO)<<year_ofOBS<<" "<<month_ofOBS<<" "<<date_ofOBS<<" "<<time_ofOBS_hh<<" "<<time_ofOBS_mm<<" "<<time_ofOBS_ss;
    float time_hh=atof(time_ofOBS_hh.c_str())/24;
    float time_mm=atof(time_ofOBS_mm.c_str())/(24*60);
    float time_ss=atof(time_ofOBS_ss.c_str())/(24*3600);
    
    float final_Date_fract=atof(date_ofOBS.c_str())-1+time_hh+time_mm+time_ss;
    
    LOG(INFO)<<year_ofOBS<<" "<<month_ofOBS<<" "<<date_ofOBS<<" "<<final_Date_fract;
    int monthDays[12]={31, 28, 31, 30, 31, 30,
                           31, 31, 30, 31, 30, 31};
    bool flagLeapyear=FALSE;
    int yearint=atoi(year_ofOBS.c_str());
    
  if (yearint%400 == 0)
         flagLeapyear=TRUE;
  else if ( yearint%100 == 0)
   flagLeapyear=FALSE;
  else if ( yearint%4 == 0 )
     flagLeapyear=TRUE;
  else
    flagLeapyear=FALSE;
    
    
    float daysCountPeryear;
    if(flagLeapyear==TRUE){
        daysCountPeryear=366;
        monthDays[1]=monthDays[1]+1;
        }
    else{
        daysCountPeryear=365;
         
    }
    float TotalDays_count=0.0f;
    int numMonth= atoi(month_ofOBS.c_str());
    for (int i=0;i<numMonth-1;i++){
        TotalDays_count=TotalDays_count+monthDays[i];        
    }
    TotalDays_count=TotalDays_count+final_Date_fract;
    
    float numyearFract=yearint+(TotalDays_count/daysCountPeryear);
    LOG(INFO)<<"OBS_EPOCH->"<<setprecision(10)<<numyearFract<<" "<<TotalDays_count<<" "<<daysCountPeryear;
    //exit(1);
  float MJD_curr_Date=  ConvertDateToMJD(final_Date_fract,atof(month_ofOBS.c_str()),atof(year_ofOBS.c_str()));
  LOG(INFO)<<MJD_curr_Date;
  
  
  dayRefMJD=MJD_curr_Date-51543.5;
  
  string starRADECfile=(string)moduleoutdir+"/"+"star_radec.txt";
  
  ofstream temp(starRADECfile.c_str());
   float Ra_diff_toAdd,dec_diff_toAdd;
   vector<float> RA_img_Stars,DEC_img_Stars;
   vector<float> RA_img_Stars_OBSEpoch,DEC_img_Stars_OBSEpoch;
   vector<float> RA_img_Stars_J2000Epoch,DEC_img_Stars_J2000Epoch;
    
      LOG(INFO)<<setprecision(10)<<RA_pnt<<" "<<DEC_pnt;
    temp<<"=============================================================="<<endl;
    float rapnt_j2000 = RA_pnt;
    float rapnt_j2000_starting=RA_pnt;
    float decpnt_j2000_starting=DEC_pnt;
    float  decpnt_j2000 = DEC_pnt;
     //converting RA_PNT and DEC_PNT to jOBS epoch.
    LOG(INFO)<<RA_pnt<<" "<<DEC_pnt;
    Convert_J2000_To_ObervationRADECVal(dayRefMJD,RA_pnt,DEC_pnt,Ra_diff_toAdd,dec_diff_toAdd);
    
  //Add difference to get JOBS coordinates for the  pointing direction
   RA_pnt=RA_pnt+Ra_diff_toAdd*((dayRefMJD));
   DEC_pnt=DEC_pnt+dec_diff_toAdd*(dayRefMJD);
//     LOG(INFO)<<RA_pnt<<" "<<DEC_pnt;exit(1);
   rapnt_obs=RA_pnt;
   decpnt_obs=DEC_pnt;
  
   LOG(INFO)<<rapnt_obs<<" "<<decpnt_obs;
   //calculate Star RADEC values in JOBS because image is in JOBS. 
   for(int i=0;i<Cx.size();i++)
    {
         star_dec=DEC_pnt+(Cy[i]-2400)*cdelt2;
         star_ra=RA_pnt+(Cx[i]-2400)*(cdelt1/cos(DEC_pnt*M_PI/180));
         LOG(INFO)<<i<<" "<<star_ra<<" "<<star_dec;
         temp<<Cx[i]<<" "<<Cy[i]<<" "<<setprecision (20)<<star_ra<<setprecision (20)<<"  "<<star_dec<<endl;
         RA_img_Stars_OBSEpoch.push_back (star_ra);
         DEC_img_Stars_OBSEpoch.push_back (star_dec);  
    }
   
   //convert back star RADEC value from JOBS to J2000
    for(int i=0;i<Cx.size();i++){
       Convert_J2000_To_ObervationRADECVal(-dayRefMJD,RA_img_Stars_OBSEpoch[i],DEC_img_Stars_OBSEpoch[i],Ra_diff_toAdd,dec_diff_toAdd); 
       RA_img_Stars_J2000Epoch.push_back(RA_img_Stars_OBSEpoch[i]+Ra_diff_toAdd*((-dayRefMJD)));
       DEC_img_Stars_J2000Epoch.push_back(DEC_img_Stars_OBSEpoch[i]+dec_diff_toAdd*(-dayRefMJD));
       //temp<<Cx[i]<<" "<<Cy[i]<<" "<<setprecision (20)<<RA_img_Stars_OBSEpoch[i]+Ra_diff_toAdd*((-dayRefMJD))<<setprecision (20)<<"  "<<DEC_img_Stars_OBSEpoch[i]+dec_diff_toAdd*(-dayRefMJD)<<endl;
     }
   
   vector<float> RA_img_Stars_J2000Epoch_backup,DEC_img_Stars_J2000Epoch_backup;
   RA_img_Stars_J2000Epoch_backup=RA_img_Stars_J2000Epoch;
   DEC_img_Stars_J2000Epoch_backup=DEC_img_Stars_J2000Epoch;
          
   LOG(INFO)<<" Now Opening Database...."<<databasename;
   if (!FileExists (databasename)) 
   {
       LOG(ERROR)<<"***DataBase not found***";
       return(EXIT_FAILURE);
   }
        
        int numStars=Cx.size();
        string newRad;
        double mean_Of_diffra,mean_Of_diffdec,mean_Of_diffra_uv,mean_Of_diffdec_uv;
         string file1=(string)moduleoutdir+"/"+"star_raDec_frmOptics_catalogueWith_5Stars.txt";
         string file2=(string)moduleoutdir+"/"+"star_raDec_frmUV_catalogue.txt";
        string  file3=(string)moduleoutdir+"/"+"star_raDec_frmOptics_catalogueWith_10Stars.txt";
        //LOG(INFO)<<numStars;exit(1);
        diff_dec_add_opt=0.0f;
        diff_ra_add_opt=0.0f;
   if(numStars/2>1){
        LOG(INFO)<<"Searching OPTIC catalogue";
        flag_NOT_FOUND_CATA_MATCH=0;
        newRad=getRaDECmatch(RA_img_Stars_J2000Epoch,DEC_img_Stars_J2000Epoch,2,len_a,len_b,rad_search,numStars/2,mean_Of_diffra,mean_Of_diffdec,file1,&numStars_frmCatamatch_optic,cos((DEC_pnt)*M_PI/180),0);
      
        if(flag_NOT_FOUND_CATA_MATCH==1)
       {       
       flag_Optic_Catalogue=0;
       
//       center_ra=rapnt_obs;
//       center_dec=decpnt_obs;     
//       
       center_ra=rapnt_j2000_starting;
       center_dec=decpnt_j2000_starting;     
       
//       LOG(INFO)<<center_ra<<" "<<center_dec;exit(1);
       fits_create_img (fout , bitpix , naxis , naxes , &status) ;
       fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata, &status) ;
      // fits_close_file(fout,&status);
       writeUsrkeywordsFrmvectNew(imagefile_out,tempL1header);
    fits_copy_file(expPointer,foutexp,1,1,1,&status);
   printError(status,"File cannot be copied",imagefile_in);
     fits_copy_file(noicemapPointer,foutnoicemap,1,1,1,&status);
   printError(status,"File cannot be copied",imagefile_in);
       }
       else{
           flag_Optic_Catalogue=1;
            LOG(INFO)<<mean_Of_diffra<<" "<<mean_Of_diffdec;
            diff_ra_add=mean_Of_diffra/cos((decpnt_j2000)*M_PI/180);
            diff_dec_add=mean_Of_diffdec;
          
            rapnt_j2000=rapnt_j2000+diff_ra_add;
            decpnt_j2000=decpnt_j2000+diff_dec_add;
          
            float rapnt_j2000_afterFirstComparison = rapnt_j2000;
            float decpnt_j2000_afterFirstComparison=decpnt_j2000;
            //again converting RADEC value of center from J2000 to JOBS;
             Convert_J2000_To_ObervationRADECVal(dayRefMJD,rapnt_j2000,decpnt_j2000,Ra_diff_toAdd,dec_diff_toAdd);
    
  //Add difference to get JOBS coordinates for the  pointing direction
   rapnt_j2000=rapnt_j2000+Ra_diff_toAdd*((dayRefMJD));
   decpnt_j2000=decpnt_j2000+dec_diff_toAdd*(dayRefMJD);
            
                        
            //Re-Calculating the centroid value.
            vector<float> ratemp,dectemp;
            ratemp=RA_img_Stars_OBSEpoch;
            dectemp=DEC_img_Stars_OBSEpoch;
                        
            RA_img_Stars_OBSEpoch.clear();
            DEC_img_Stars_OBSEpoch.clear();
            
          
         //temp<<"=========================================================="<<endl;
            for (int i = 0; i < Cx.size(); i++) {
              
                star_dec = decpnt_j2000 + (Cy[i] - 2400) * cdelt2;
                star_ra = rapnt_j2000 + (Cx[i] - 2400)*(cdelt1 / cos(decpnt_j2000 * M_PI / 180));

               // temp << Cx[i] << " " << Cy[i] << " " << setprecision(20) << star_ra << setprecision(20) << "  " << star_dec << endl;
                RA_img_Stars_OBSEpoch.push_back(star_ra);
                DEC_img_Stars_OBSEpoch.push_back(star_dec);
            }
         
            temp<<"================================================================"<<endl;
          RA_img_Stars_J2000Epoch.clear();
          DEC_img_Stars_J2000Epoch.clear(); 
            
         for(int i=0;i<Cx.size();i++){
       Convert_J2000_To_ObervationRADECVal(-dayRefMJD,RA_img_Stars_OBSEpoch[i],DEC_img_Stars_OBSEpoch[i],Ra_diff_toAdd,dec_diff_toAdd); 
       RA_img_Stars_J2000Epoch.push_back(RA_img_Stars_OBSEpoch[i]+Ra_diff_toAdd*((-dayRefMJD)));
       DEC_img_Stars_J2000Epoch.push_back(DEC_img_Stars_OBSEpoch[i]+dec_diff_toAdd*(-dayRefMJD));
         temp<<Cx[i]<<" "<<Cy[i]<<" "<<setprecision (20)<<RA_img_Stars_OBSEpoch[i]+Ra_diff_toAdd*((-dayRefMJD))<<setprecision (20)<<"  "<<DEC_img_Stars_OBSEpoch[i]+dec_diff_toAdd*(-dayRefMJD)<<endl;
         }
        temp.close(); 
        
          flag_NOT_FOUND_CATA_MATCH=0;
          newRad=getRaDECmatch(RA_img_Stars_J2000Epoch,DEC_img_Stars_J2000Epoch,2,len_a,len_b,rad_search,numStars,mean_Of_diffra,mean_Of_diffdec,file3,&numStars_frmCatamatch_optic,cos((DEC_pnt)*M_PI/180),0);  
        
            diff_ra_add=diff_ra_add+mean_Of_diffra/cos((DEC_pnt)*M_PI/180);
            diff_dec_add=diff_dec_add+mean_Of_diffdec;
            diff_ra_add_opt=diff_ra_add;
            diff_dec_add_opt=diff_dec_add;
            rapnt_j2000_afterFirstComparison=rapnt_j2000_afterFirstComparison+diff_ra_add;
            decpnt_j2000_afterFirstComparison=decpnt_j2000_afterFirstComparison+diff_dec_add;
           
            
            //now Comparing List 1 of UVIT-JOBS-radec with J2000 catalogue RADEC ;
            vector<float> JOBS_UVIT_CENT_X,JOBS_UVIT_CENT_Y,J2000_USNO_CENT_X,J2000_USNO_CENT_Y;
            vector<float>  trackRA_frmCata,trackDEC_frmCata;
         
            for (int i = 0; i < track_decback.size(); i++) {
                    if (track_decback[i] != -9999 && track_raback[i] != -9999) {
                        for (int j = 0; j < validStar_index_filtered.size(); j++) {
                            if(i==validStar_index_filtered[j]){
                               // JOBS_UVIT_CENT_X.push_back(Cx[i]);
                               // JOBS_UVIT_CENT_Y.push_back(Cy[i]);
                                JOBS_UVIT_CENT_X.push_back(Cx[i]);
                                JOBS_UVIT_CENT_Y.push_back(Cy[i]);
                                J2000_USNO_CENT_X.push_back(((track_raback[i]-rapnt_j2000_starting)/(cdelt1/cos(decpnt_j2000_starting * M_PI / 180)))+2400);
                                J2000_USNO_CENT_Y.push_back(((track_decback[i]-decpnt_j2000_starting)/(cdelt2))+2400);
                                //J2000_USNO_CENT_X.push_back(((track_raback[i]-rapnt_obs)/(cdelt1/cos(decpnt_obs * M_PI / 180)))+2400);
                               // J2000_USNO_CENT_Y.push_back(((track_decback[i]-decpnt_obs)/(cdelt2))+2400);
                                
                                trackRA_frmCata.push_back(track_raback[i]);
                                trackDEC_frmCata.push_back(track_decback[i]);
                            }
                        }
                    }
                 }
            
            for (int i=0;i<J2000_USNO_CENT_Y.size();i++){
                LOG(INFO)<<"OPT-> "<<trackRA_frmCata[i]<<" "<<trackDEC_frmCata[i]<<" "<<JOBS_UVIT_CENT_X[i]<<" "<<JOBS_UVIT_CENT_Y[i]<<" "<<J2000_USNO_CENT_X[i]<<" "<<J2000_USNO_CENT_Y[i];
                
            }
            
            
            long totalelements=JOBS_UVIT_CENT_X.size();
             double    Xdx,Ydy,Theta;
            if(totalelements>2){
             spMatrix B((totalelements) * 2, 1);
        spMatrix A((totalelements) * 2, 3);
        spMatrix X(3, 1);
        int temp = 0;

        for (int t = 0; t < totalelements * 2; t = t + 2) {
            B(t, 0) = (J2000_USNO_CENT_X[temp]/8 - JOBS_UVIT_CENT_X[temp]/8); //+IMAGE_ARRAYSIZE*0.5;
            B(t + 1, 0) = (J2000_USNO_CENT_Y[temp]/8 - JOBS_UVIT_CENT_Y[temp]/8); //+IMAGE_ARRAYSIZE*0.5;

            A(t, 0) = -1.0 * (JOBS_UVIT_CENT_Y[temp]/8 - ((4800/8) / 2));
            A(t, 1) = 1;
            A(t, 2) = 0;

            A(t + 1, 0) = (JOBS_UVIT_CENT_X[temp]/8 - ((4800/8) / 2));
            A(t + 1, 1) = 0;
            A(t + 1, 2) = 1;

            temp++;
        }

        X.ApplyLeastSquare(A, B);

       Xdx = X(1, 0)*8;//
        Ydy = X(2, 0)*8;
        Theta = X(0, 0);
      LOG(INFO)<<Xdx<<" "<<Ydy<<" "<<Theta;
            }
            else if (totalelements<=2 && totalelements>0){
                 vector<double> diff_x_cumm,diff_y_cumm;
          for (int t = 0 ; t < totalelements  ; t++)
          {
             
              diff_x_cumm.push_back (J2000_USNO_CENT_X[t] - JOBS_UVIT_CENT_X[t]);
              diff_y_cumm.push_back (J2000_USNO_CENT_Y[t] - JOBS_UVIT_CENT_Y[t]);
          }
        Xdx=getmean (diff_x_cumm.data (),diff_x_cumm.size ());
        Ydy=getmean (diff_y_cumm.data (),diff_y_cumm.size ());
        Theta=0.0f;
                
            }
            
      //convert to RADEC shift.
      double RASHIFT =Xdx*(cdelt1 / cos(decpnt_obs * M_PI / 180));
      double DECSHIFT=Ydy*cdelt2;
              
      //LOG(INFO)<<RASHIFT<<" "<<DECSHIFT;exit(1);
      
      //
      
      //Now Apply Rotation and shift
      float *ExpandedframeData = new float[xsize*2*ysize*2];
      float *ExpandedExpData = new float[xsize*2*ysize*2];
      float *ExpandedNoiseData= new float[xsize*2*ysize*2];      
      
      float *counterArrayFrm = new float[xsize*2*ysize*2];//initArray(counterArrayFrm ,xsize*2*ysize*2,1.0f);
      float *counterArrayExp = new float[xsize*2*ysize*2];
      float *counterArraynoise= new float[xsize*2*ysize*2];
      
      for (int i =0;i<xsize*ysize;i++){
          framedata[i]=framedata[i]/4;
          Expdata[i]=Expdata[i]/4;
          noicemapData[i]=noicemapData[i]/4;         
      }
      
      
      performSubDivisionIM(framedata,xsize,ysize,ExpandedframeData,xsize*2,ysize*2);
      performSubDivisionIM(Expdata,xsize,ysize,ExpandedExpData,xsize*2,ysize*2);
      performSubDivisionIM(noicemapData,xsize,ysize,ExpandedNoiseData,xsize*2,ysize*2);
      
      float *Rotated_FrmArray = new float[xsize*2*ysize*2];initArray(Rotated_FrmArray,xsize*2*ysize*2,0.0f);
      float *Rotated_ExpArray=new float[xsize*2*ysize*2];initArray(Rotated_ExpArray,xsize*2*ysize*2,0.0f);
      float *Rotated_NoiseArray = new float[xsize*2*ysize*2];initArray(Rotated_NoiseArray,xsize*2*ysize*2,0.0f);
      
      float *BinnedFrmData =new float[xsize*ysize];
      float *BinnedExpData=new float[xsize*ysize];
      float *BinnedNoiseData=new float[xsize*ysize];
      
      double mid_X_new,mid_Y_new;
      double ctheta =cos(Theta);
      double stheta =sin(Theta);
      LOG(INFO)<<flag_inputImage;
   
      if(flag_inputImage==1){
        for (int i = 0; i < xsize * 2; i++) {
                mid_X_new = i - xsize;
                for (int j = 0; j < xsize * 2; j++) {
                    mid_Y_new = j - ysize;
                    x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + xsize;//+Xdx;
                    x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + ysize;//+Ydy;

                    //                  if((int)(round(x2)*xsize*2+round(x1))> 0 && (int)(round(x2)*xsize*2+round(x1)) <xsize*2*xsize*2)
                    if (round(x1) > 0 && round(x1) < xsize * 2 && round(x2) > 0 && round(x2) < ysize * 2) {
                        if (ExpandedframeData[j * xsize * 2 + i] != INVALID_PIX_VALUE && Rotated_FrmArray[(int) (round(x2) * xsize * 2 + round(x1))]!=INVALID_PIX_VALUE ) {
                            Rotated_FrmArray[(int) (round(x2) * xsize * 2 + round(x1))] = Rotated_FrmArray[(int) (round(x2) * xsize * 2 + round(x1))] + ExpandedframeData[j * xsize * 2 + i];
                           
                        } else {
                            Rotated_FrmArray[(int) (round(x2) * xsize * 2 + round(x1))] = INVALID_PIX_VALUE;
                        }

                    }
                }
            }
      }
      else {
             for (int i = 0; i < nrows_snr; i++) {
        mid_X_new = X_loc_snr[i]-xsize/2;
        mid_Y_new = Y_loc_snr[i]-ysize/2;

        x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + xsize/2;//+Xdx;
        x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + ysize/2;//+Ydy;

        if (round(x1) > 0 && round(x1) < xsize && round(x2) > 0 && round(x2) < ysize) {
            BinnedFrmData[(int) (round(x2) * xsize  + round(x1))]=BinnedFrmData[(int) (round(x2) * xsize  + round(x1))]+bad_flag_snr[i]*mult_phtn_snr[i]*enp_snr[i];
        }
       
       
        

    }
         
          
      }
      
      for (int i = 0; i < xsize * 2; i++) {
                mid_X_new = i - xsize;
                for (int j = 0; j < xsize * 2; j++) {
                    mid_Y_new = j - ysize;
                    x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + xsize;//+Xdx;
                    x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + ysize;//+Ydy;

                    //                  if((int)(round(x2)*xsize*2+round(x1))> 0 && (int)(round(x2)*xsize*2+round(x1)) <xsize*2*xsize*2)
                    if (round(x1) > 0 && round(x1) < xsize * 2 && round(x2) > 0 && round(x2) < ysize * 2) {
                        if (ExpandedExpData[j * xsize * 2 + i] != INVALID_PIX_VALUE && Rotated_ExpArray[(int) (round(x2) * xsize * 2 + round(x1))]!=INVALID_PIX_VALUE ) {
                            Rotated_ExpArray[(int) (round(x2) * xsize * 2 + round(x1))] = Rotated_ExpArray[(int) (round(x2) * xsize * 2 + round(x1))] + ExpandedExpData[j * xsize * 2 + i];
                           
                        } else {
                            Rotated_ExpArray[(int) (round(x2) * xsize * 2 + round(x1))] = INVALID_PIX_VALUE;
                        }

                    }
                }
            }
      
      
       for (int i = 0; i < xsize * 2; i++) {
                mid_X_new = i - xsize;
                for (int j = 0; j < xsize * 2; j++) {
                    mid_Y_new = j - ysize;
                    x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + xsize;//+Xdx;
                    x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + ysize;//+Ydy;

                    //                  if((int)(round(x2)*xsize*2+round(x1))> 0 && (int)(round(x2)*xsize*2+round(x1)) <xsize*2*xsize*2)
                    if (round(x1) > 0 && round(x1) < xsize * 2 && round(x2) > 0 && round(x2) < ysize * 2) {
                        if (ExpandedNoiseData[j * xsize * 2 + i] != INVALID_PIX_VALUE && Rotated_NoiseArray[(int) (round(x2) * xsize * 2 + round(x1))]!=INVALID_PIX_VALUE ) {
                            Rotated_NoiseArray[(int) (round(x2) * xsize * 2 + round(x1))] = Rotated_NoiseArray[(int) (round(x2) * xsize * 2 + round(x1))] + ExpandedNoiseData[j * xsize * 2 + i];
                           
                        } else {
                            Rotated_NoiseArray[(int) (round(x2) * xsize * 2 + round(x1))] = INVALID_PIX_VALUE;
                        }

                    }
                }
            }
      if(flag_inputImage==1){
      ApplyBinning(Rotated_FrmArray,xsize*2,ysize*2,BinnedFrmData,xsize,ysize,counterArrayFrm);
      }
    
          
             
          
          
    
      ApplyBinning(Rotated_ExpArray,xsize*2,ysize*2,BinnedExpData,xsize,ysize,counterArrayExp);
      ApplyBinning(Rotated_NoiseArray,xsize*2,ysize*2,BinnedNoiseData,xsize,ysize,counterArraynoise);
      if(flag_inputImage==0){
       for (int i =0;i<xsize*ysize;i++)  {
           if(BinnedExpData[i]!=0)
              BinnedFrmData[i]=BinnedFrmData[i]/BinnedExpData[i];              
           else 
               BinnedFrmData[i]=0.0f;
              
          }
      }
      
    long naxes[2] ;
    naxes[0] = naxes[1] = xsize ;
    int naxis = 2 ;
    int bitpix = FLOAT_IMG ;
         
          fits_create_img (fout , bitpix , naxis , naxes , &status) ;
//         fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , Rotated_frmData , &status) ;
          fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , BinnedFrmData, &status) ;
        
            fits_create_img (foutexp , bitpix , naxis , naxes , &status) ;

           fits_write_pix (foutexp , TFLOAT , fpixel , xsize*ysize , BinnedExpData , &status) ;
         
          fits_create_img (foutnoicemap , bitpix , naxis , naxes , &status) ;

          fits_write_pix (foutnoicemap , TFLOAT , fpixel , xsize*ysize , BinnedNoiseData , &status) ;
         // LOG(INFO)<<RASHIFT<<" "<<DECSHIFT;exit(1);
//          center_ra=rapnt_obs+RASHIFT;
//          center_dec=decpnt_obs+DECSHIFT;
          
          center_ra=rapnt_j2000_starting+RASHIFT;
          center_dec=decpnt_j2000_starting+DECSHIFT;
          
          LOG(INFO)<<center_ra<<" "<<center_dec;
            
          
            
       }
        
        
       }
   else{
       flag_Optic_Catalogue=0;
   LOG(ERROR)<<"***Total number of stars found is less than 3 ,using RA_PNT and DEC_PNT***";
   flag_NOT_FOUND_CATA_MATCH=1;   
//   center_ra=rapnt_obs;
  // center_dec=decpnt_obs;
   center_ra=rapnt_j2000_starting;
   center_dec=decpnt_j2000_starting;
   
   
   fits_create_img (fout , bitpix , naxis , naxes , &status) ;
//         fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , Rotated_frmData , &status) ;
          fits_write_pix (fout , TFLOAT , fpixel , xsize*ysize , framedata, &status) ;
          //fits_close_file(fout,&status);
              writeUsrkeywordsFrmvectNew(imagefile_out,tempL1header);
     fits_copy_file(expPointer,foutexp,1,1,1,&status);
   printError(status,"File cannot be copied",imagefile_in);
     fits_copy_file(noicemapPointer,foutnoicemap,1,1,1,&status);
   printError(status,"File cannot be copied",imagefile_in);
   
//   fits_copy_file(fin,fout,1,1,1,&status);
//    printError(status,"File cannot be copied",imagefile_in);
   }
        
    //NOw searching UV CATALOGUE    
        LOG(INFO)<<"NOW UV CATALOGUE....";
        RA_img_Stars_J2000Epoch.clear();
        DEC_img_Stars_J2000Epoch.clear();
        RA_img_Stars_J2000Epoch=RA_img_Stars_J2000Epoch_backup;
        DEC_img_Stars_J2000Epoch=DEC_img_Stars_J2000Epoch_backup;
        rapnt_j2000=rapnt_j2000_starting;//resetting the value;
        decpnt_j2000=decpnt_j2000_starting;
        if(numStars/2>3){
      
        newRad=getRaDECmatch(RA_img_Stars_J2000Epoch,DEC_img_Stars_J2000Epoch,6,len_a,len_b,rad_search,numStars/2,mean_Of_diffra,mean_Of_diffdec,file1,&numStars_frmCatamatch_UV,cos((DEC_pnt)*M_PI/180),1);
    if(flag_NOT_FOUND_CATA_MATCH==1)
       { 
        flag_UV_Catalogue=0;        
    }
    else{
        flag_UV_Catalogue=1;
        
        diff_ra_add=mean_Of_diffra/cos((decpnt_j2000)*M_PI/180);
            diff_dec_add=mean_Of_diffdec;
          
            rapnt_j2000=rapnt_j2000+diff_ra_add;
            decpnt_j2000=decpnt_j2000+diff_dec_add;
            
         float rapnt_j2000_afterFirstComparison = rapnt_j2000;
            float decpnt_j2000_afterFirstComparison=decpnt_j2000;
            
             //again converting RADEC value of center from J2000 to JOBS;
                 Convert_J2000_To_ObervationRADECVal(dayRefMJD,rapnt_j2000,decpnt_j2000,Ra_diff_toAdd,dec_diff_toAdd);
    
  //Add difference to get JOBS coordinates for the  pointing direction
   rapnt_j2000=rapnt_j2000+Ra_diff_toAdd*((dayRefMJD));
   decpnt_j2000=decpnt_j2000+dec_diff_toAdd*(dayRefMJD);
          RA_img_Stars_OBSEpoch.clear();
          DEC_img_Stars_OBSEpoch.clear();
                
   
    for (int i = 0; i < Cx.size(); i++) {
              
                star_dec = decpnt_j2000 + (Cy[i] - 2400) * cdelt2;
                star_ra = rapnt_j2000 + (Cx[i] - 2400)*(cdelt1 / cos(decpnt_j2000 * M_PI / 180));

             //   temp << Cx[i] << " " << Cy[i] << " " << setprecision(20) << star_ra << setprecision(20) << "  " << star_dec << endl;
                RA_img_Stars_OBSEpoch.push_back(star_ra);
                DEC_img_Stars_OBSEpoch.push_back(star_dec);
            }
            
            RA_img_Stars_J2000Epoch.clear();
          DEC_img_Stars_J2000Epoch.clear(); 
          
             for(int i=0;i<Cx.size();i++){
       Convert_J2000_To_ObervationRADECVal(-dayRefMJD,RA_img_Stars_OBSEpoch[i],DEC_img_Stars_OBSEpoch[i],Ra_diff_toAdd,dec_diff_toAdd); 
       RA_img_Stars_J2000Epoch.push_back(RA_img_Stars_OBSEpoch[i]+Ra_diff_toAdd*((-dayRefMJD)));
       DEC_img_Stars_J2000Epoch.push_back(DEC_img_Stars_OBSEpoch[i]+dec_diff_toAdd*(-dayRefMJD));
     }
          
        newRad=getRaDECmatch(RA_img_Stars_J2000Epoch,DEC_img_Stars_J2000Epoch,6,len_a,len_b,rad_search,numStars,mean_Of_diffra,mean_Of_diffdec,file1,&numStars_frmCatamatch_UV,cos((DEC_pnt)*M_PI/180),1);
   
        diff_ra_add=mean_Of_diffra/cos((DEC_pnt)*M_PI/180);
        diff_dec_add=mean_Of_diffdec;
          
        rapnt_j2000_afterFirstComparison=rapnt_j2000_afterFirstComparison+diff_ra_add;
        decpnt_j2000_afterFirstComparison=decpnt_j2000_afterFirstComparison+diff_dec_add;
           
         vector<float> JOBS_UVIT_CENT_X,JOBS_UVIT_CENT_Y,J2000_USNO_CENT_X,J2000_USNO_CENT_Y;
           vector<float>  trackRA_frmCata,trackDEC_frmCata;
            for (int i = 0; i < track_decback.size(); i++) {
                    if (track_decback[i] != -9999 && track_raback[i] != -9999) {
                        for (int j = 0; j < validStar_index_filtered.size(); j++) {
                            if(i==validStar_index_filtered[j]){
//                                JOBS_UVIT_CENT_X.push_back(Cx[i]);
//                                JOBS_UVIT_CENT_Y.push_back(Cy[i]);
                                JOBS_UVIT_CENT_X.push_back(Cx[i]);
                                JOBS_UVIT_CENT_Y.push_back(Cy[i]);
//                                J2000_USNO_CENT_X.push_back(((track_raback[i]-rapnt_j2000_starting)/(cdelt1/cos(decpnt_j2000_starting * M_PI / 180)))+2400);
//                                J2000_USNO_CENT_Y.push_back(((track_decback[i]-decpnt_j2000_starting)/(cdelt2))+2400);
                                J2000_USNO_CENT_X.push_back(((track_raback[i]-rapnt_j2000_starting)/(cdelt1/cos(decpnt_j2000_starting * M_PI / 180)))+2400);
                                J2000_USNO_CENT_Y.push_back(((track_decback[i]-decpnt_j2000_starting)/(cdelt2))+2400);
                                
                                trackRA_frmCata.push_back(track_raback[i]);
                                trackDEC_frmCata.push_back(track_decback[i]);
                            
                            }
                        }
                    }
            }
            
            for (int i=0;i<J2000_USNO_CENT_Y.size();i++){
                LOG(INFO)<<"UV->"<<trackRA_frmCata[i]<<" "<<trackDEC_frmCata[i]<<" "<<JOBS_UVIT_CENT_X[i]<<" "<<JOBS_UVIT_CENT_Y[i]<<" "<<J2000_USNO_CENT_X[i]<<" "<<J2000_USNO_CENT_Y[i];
                
            }
            
            
            long totalelements=JOBS_UVIT_CENT_X.size();
  double    Xdx,Ydy,Theta;
if(totalelements >2){
             spMatrix B((totalelements) * 2, 1);
        spMatrix A((totalelements) * 2, 3);
        spMatrix X(3, 1);
        int temp = 0;

        for (int t = 0; t < totalelements * 2; t = t + 2) {
            B(t, 0) = (J2000_USNO_CENT_X[temp]/8 - JOBS_UVIT_CENT_X[temp]/8); //+IMAGE_ARRAYSIZE*0.5;
            B(t + 1, 0) = (J2000_USNO_CENT_Y[temp]/8 - JOBS_UVIT_CENT_Y[temp]/8); //+IMAGE_ARRAYSIZE*0.5;

            A(t, 0) = -1.0 * (JOBS_UVIT_CENT_Y[temp]/8 - ((4800/8) / 2));
            A(t, 1) = 1;
            A(t, 2) = 0;

            A(t + 1, 0) = (JOBS_UVIT_CENT_X[temp]/8 - ((4800/8) / 2));
            A(t + 1, 1) = 0;
            A(t + 1, 2) = 1;

            temp++;
        }

        X.ApplyLeastSquare(A, B);

       Xdx = X(1, 0)*8;//
         Ydy = X(2, 0)*8;
         Theta = X(0, 0);

     LOG(INFO)<<Xdx<<" "<<Ydy<<" "<<Theta;
}  else if (totalelements<=2 && totalelements>0){
                 vector<double> diff_x_cumm,diff_y_cumm;
          for (int t = 0 ; t < totalelements  ; t++)
          {
             
              diff_x_cumm.push_back (J2000_USNO_CENT_X[t] - JOBS_UVIT_CENT_X[t]);
              diff_y_cumm.push_back (J2000_USNO_CENT_Y[t] - JOBS_UVIT_CENT_Y[t]);
          }
        Xdx=getmean (diff_x_cumm.data (),diff_x_cumm.size ());
        Ydy=getmean (diff_y_cumm.data (),diff_y_cumm.size ());
        Theta=0.0f;
                
            }





 
      
      //convert to RADEC shift.
      double RASHIFT =Xdx*(cdelt1 / cos(DEC_pnt * M_PI / 180));
      double DECSHIFT=Ydy*cdelt2;
      double THETADEG=Theta*180/M_PI;
      
      center_RA_UV=rapnt_obs+RASHIFT;
      center_DEC_UV=decpnt_obs+DECSHIFT;
      center_ROLL_UV=THETADEG;
      
     // x1 = ((mid_X_new * (ctheta)) - (mid_Y_new * (stheta))) + xsize;//+Xdx;
                    //x2 = ((mid_X_new * (stheta)) +(mid_Y_new * (ctheta))) + ysize;//+Ydy;
      
        
    }
        
        
        
        }
        else{
            flag_UV_Catalogue=0;            
        }
        
        
     fits_close_file(expPointer,&status);
     fits_close_file(noicemapPointer,&status);
       
            writeUsrkeywordsFrmvectNew(imagefile_out,tempL1header);
             writeUsrkeywordsFrmvectNew(imagefile_out_expName,tempL1header);
            writeUsrkeywordsFrmvectNew(imagefile_out_noicemapName,tempL1header);
        addWCS(fout);
         
         addWCS(foutexp);
         addWCS(foutnoicemap);

         fits_close_file(fout,&status);
          fits_close_file(foutexp,&status);
           fits_close_file(foutnoicemap,&status);
     
    return (EXIT_SUCCESS);
        
}


int uvtFullFrameAst::addWCS(fitsfile *fout)
{
    
    int status=0;
    //const char *teldeffile = caldb_handler.getTelDefFile(caldbDir);
    int xsize, ysize;
    double xscale, yscale;    
//    readKeywords((char *)teldeffile,1,5, TINT, (char *)"DET_XSIZ", &xsize,
//                                                TINT,(char *)"DET_YSIZ",&ysize,
//                                                TDOUBLE,(char *)"DET_XSCL",&xscale,
//                                                TDOUBLE,(char *)"DET_YSCL", &yscale,
//                                                TDOUBLE, (char *)"FOCALLEN",&focallength);
    
    xsize =4800; //teldef.det_xsiz;
    ysize = 4800;//teldef.det_ysiz;
     float cdelt1,cdelt2;
    int factor_delta=xsize/600;
    if(strcmp(datainfo.getDetector (),"VIS")==0)
    {
           cdelt1=(3.357/3600)/factor_delta;
           cdelt2=(3.311/3600)/factor_delta;
    }
    else if (strcmp(datainfo.getDetector (),"FUV")==0){
        //cdelt1=3.357/3600;
       //  cdelt2=3.311/3600;
       cdelt1=cdelt2=(3.3373/3600)/factor_delta;
    }
    else if(strcmp(datainfo.getDetector (),"NUV")==0){
       cdelt1=cdelt2=(3.3307/3600)/factor_delta;
       // cdelt1=cdelt2=(3.3307*factor_delta/3600)/;
       
    }
   // cdelt1=-cdelt1/cos(center_dec*M_PI/180);
     cdelt1=-cdelt1;///cos(center_dec*M_PI/180);
//    xscale =3.32/3600; //teldef.det_xscl;
//    yscale= 3.32/3600;//teldef.det_yscl;
    
    float crpix1 =  xsize/2.0;
    float crpix2 =  ysize/2.0; 
     //cdelt1 = (float) (atan(xscale / teldef.focal_len));               //in radians
    // cdelt2 = (float) (atan(yscale / teldef.focal_len));               //in radians
   
   
    //cdelt1 = cdelt1*180/M_PI;
   // cdelt2 = cdelt2 *180/M_PI;
    
   // cdelt1=cdelt2=(float)atan((5.0/1000)/(512.0*FOCAL_LENGTH))*180.0/M_PI;
    
    //cdelt1=cdelt2=3.32/3600;
    //twist=0;
   // cdelt1=cdelt1*8;
    //cdelt2=cdelt2*8;
    
    LOG(INFO)<<"cdelt1 :"<<cdelt1<<"  cdelt2:"<<cdelt2;
    LOG(INFO)<<"Adding WCS keywords";
    LOG(INFO)<<"CRPIX1 :"<<crpix1;
    LOG(INFO)<<"CDELT1 :"<<cdelt1;
    LOG(INFO)<<"CRVAL1 :"<<center_ra;
    LOG(INFO)<<"CRPIX2 :"<<crpix2;
    LOG(INFO)<<"CDELT2 :"<<cdelt2;
    LOG(INFO)<<"CRVAL2 :"<<center_dec;   
    LOG(INFO)<<"CROTA2 :"<<twist;
    
   fits_update_key(fout, TFLOAT, "CENTER_RA_fromQuartions", &center_ra_prev, "Right Ascension", &status);     printError(status,"Error in writing the key value of CENTER_RA_fromQuartions");
   fits_update_key(fout, TFLOAT, "CENTER_DEC_fromQuartions", &center_dec_prev, "Declination", &status);     printError(status,"Error in writing the key value of CENTER_DEC_fromQuartions");
   
   if(flag_Optic_Catalogue==1){
        fits_update_key(fout, TSTRING, "Optic_Catalogue_flag ", (char*)"SUCCESS", "flag for Optic catalogue", &status);     printError(status,"Error in writing the key value of Optic_Catalogue_flag");
         fits_update_key(fout, TINT, "TOTAL MATCH FOUND FROM CATALOGUE OPTIC ", &numStars_frmCatamatch_optic, "total matched stars", &status);   printError(status,"Error in writing the key value of UV_Catalogue_flag");
     }
     else if(flag_Optic_Catalogue==0){
         fits_update_key(fout, TSTRING, "Optic_Catalogue_flag ", (char*)"FAILURE", "flag for Optic catalogue", &status);     printError(status,"Error in writing the key value of Optic_Catalogue_flag");
    }
   if(flag_UV_Catalogue==1)
   {
       fits_update_key(fout, TSTRING, "UV_Catalogue_flag ", (char*)"SUCCESS", "flag for UV catalgue", &status);   printError(status,"Error in writing the key value of UV_Catalogue_flag");
       fits_update_key(fout, TFLOAT, "CENTER_RA_UV ", &center_RA_UV, "", &status);   printError(status,"Error in writing the key value of UV_Catalogue_flag");
       fits_update_key(fout, TFLOAT, "CENTER_DEC_UV ", &center_DEC_UV, "", &status);   printError(status,"Error in writing the key value of UV_Catalogue_flag");
       fits_update_key(fout, TFLOAT, "CENTER_ROLL_UV ", &center_ROLL_UV, "", &status);   printError(status,"Error in writing the key value of UV_Catalogue_flag");
       fits_update_key(fout, TINT, "TOTAL MATCH FOUND FROM CATALOGUE UV ", &numStars_frmCatamatch_UV, "total matched stars", &status);   printError(status,"Error in writing the key value of UV_Catalogue_flag");
   }
   else if(flag_UV_Catalogue==0){
       fits_update_key(fout, TSTRING, "UV_Catalogue_flag ", (char*)"FAILURE", "flag for UV catalgue", &status);   printError(status,"Error in writing the key value of UV_Catalogue_flag");
   }
   float  crota1=0;
   float  crota2=0;
  
   fits_update_key(fout, TSTRING, "CTYPE1", (char *) "RA---TAN", "Right Ascension", &status);     printError(status,"Error in writing the key value of RA-TAN");
    fits_update_key(fout, TSTRING, "CUNIT1",(char *) "deg", "Unit", &status);                                   printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX1", &crpix1, "Reference Pixel", &status);                            printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT1", &cdelt1, "", &status);                                                   printError(status,"");
   // if((center_ra<360 && center_ra>-360)|| (center_dec<360 && center_dec>-360)){
   //     center_ra=0;
   //     center_dec=0;
   // }
    fits_update_key(fout, TFLOAT, "CRVAL1", &center_ra, "", &status);                                              printError(status,"");  
    fits_update_key(fout, TSTRING, "CTYPE2", (char *) "DEC--TAN", "", &status);                     printError(status,"");
    fits_update_key(fout, TSTRING, "CUNIT2", (char *)"deg", "Unit", &status);                           printError(status,"");
    fits_update_key(fout, TFLOAT, "CRPIX2", &crpix2, "Reference Pixel", &status);                    printError(status,"");
    fits_update_key(fout, TFLOAT, "CDELT2", &cdelt2, "", &status);                                           printError(status,"");
    fits_update_key(fout, TFLOAT, "CRVAL2", &center_dec, "", &status);                                    printError(status,"");   
    fits_update_key(fout, TFLOAT, "CROTA2", &crota2, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(fout, TFLOAT, "CROTA1", &crota1, "Twist Angle", &status);                          printError(status,"");
    fits_update_key(fout, TFLOAT, "ROLLangleapplied  ", &twist, "", &status);                          printError(status,"");
    return(EXIT_SUCCESS);
}

//int uvtFullFrameAst::getRaDecTwist()
//{
//        
//    int status=0;
//    
//     //Reading UVIT alignment from teldef file
//    uvitAlign.q1=(sqrt(1+teldef.m11+teldef.m22+teldef.m33))/2.0;
//    uvitAlign.q2=(teldef.m32-teldef.m23)/(4*uvitAlign.q1);
//    uvitAlign.q3=(teldef.m13-teldef.m31)/(4*uvitAlign.q1);
//    uvitAlign.q4=(teldef.m21-teldef.m12)/(4*uvitAlign.q1);
//        
//    LOG(INFO)<<"UVIT alignment quaternion is ("<<uvitAlign.q1<<","<<
//                        uvitAlign.q2<<","<<uvitAlign.q3<<","<<uvitAlign.q4<<")";
//
//    //Reading UVIT attitude from attitude file
//    vector<Attitude> attvect;
//    status = readAttitude(attitudefile,att_timecol,att_qcol,datainfo.getTstart(),datainfo.getTstop(),attvect);
//    if(status) {
//        LOG(ERROR)<<"\033[1;31m***Error reading attitude data ***\033[0m";
//        return(EXIT_FAILURE);
//    }
//       
//    Axis uvitnormal(0,1,0);                              //considering UVIT center normal axis be Y axis
//    uvitnormal.normalize();
//    
//    Axis uvitX(1,0,0);                      //x axis of UVIT payload, used for twist angle computation
//    uvitX.normalize();
//    
//    Axis uvitnormal_inertial;                    //UVIT normal in inertial coordinates
//    Axis uvitX_inertial;
//    
//    double sum_x=0,sum_y=0,sum_z=0,x,y,z;
//    double sum_xt=0, sum_yt=0, sum_zt=0;
//    
//    for(long i=0;i<attvect.size();i++){
//        
//         Q qbi;                   //body to inertial quaternion at time t;
//         Q qatt(attvect[i].q1,attvect[i].q2,attvect[i].q3,attvect[i].q4);
//         
//         //LOG(INFO)<<"Attitude Quaternion  ";       qatt.display();
//                
//         quaternion_product(qatt,uvitAlign,qbi);
//         qbi.normalize();
//         //LOG(INFO)<<"Qbi ";                     qbi.display();
//         
//         rotate(uvitnormal,qbi,uvitnormal_inertial);                        //for UVIT center vector
//         
//         sum_x+=uvitnormal_inertial.x;
//         sum_y+=uvitnormal_inertial.y;
//         sum_z+=uvitnormal_inertial.z;
//         
//         rotate(uvitX,qbi,uvitX_inertial);
//         
//         sum_xt+=uvitX_inertial.x;
//         sum_yt+=uvitX_inertial.y;
//         sum_zt+=uvitX_inertial.z;
//     }
//
//    x = sum_x/attvect.size();        y=sum_y/attvect.size();       z= sum_z/attvect.size();
//    uvitnormal_inertial.update(x,y,z);
//    uvitnormal_inertial.normalize();
//    LOG(INFO)<<"UVIT Inertial ";           uvitnormal_inertial.display();
//  
//    center_ra = acos( uvitnormal_inertial.y / sqrt(uvitnormal_inertial.x*uvitnormal_inertial.x+ uvitnormal_inertial.y*uvitnormal_inertial.y)) ;                         
//    if(uvitnormal_inertial.y<0)                 center_ra = 2*M_PI-center_ra;
//    center_dec = asin(uvitnormal_inertial.z) ;
//    
//    center_ra = center_ra*180/M_PI;
//    center_dec = center_dec * 180/M_PI;
//        
//    //for twist angle
//    x=sum_xt/attvect.size();          y=sum_yt/attvect.size();         z=sum_zt/attvect.size();
//    uvitX_inertial.update(x,y,z);
//    uvitX_inertial.normalize ();
//    
//    double cos_ra = cos(center_ra * M_PI/180);
//    double sin_dec = sin(center_dec * M_PI/180);
//    
//   Axis ref(cos_ra,sin_dec,0);                                 //for twist angle computation
//   // Axis ref(1,0,0);                                 //for twist angle computation
//   LOG(INFO)<<"Reference  axes"<<endl;
//   ref.display ();
//    twist=(acos((ref.x*uvitX_inertial.x+ref.y*uvitX_inertial.y+ref.z*uvitX_inertial.z)/(uvitX_inertial.getMod()*ref.getMod()))) * 180/M_PI ;             //in degrees 
//    if(uvitX_inertial.z<0) twist=360-twist;
//    if(twist==NAN) {
//             twist=0;
//    }
//    
//    LOG(INFO)<<"Center RA : "<<center_ra;
//    LOG(INFO)<<"Center DEC : "<<center_dec;
//    LOG(INFO)<<"Twist: "<<twist<<"  deg";
//    
//    return (EXIT_SUCCESS);
//}

int uvtFullFrameAst::getHistory(vector<string> &vhistory) {
    int cnt=0;
   // char *user = getlogin();
   // string str = "Module run by " + (string) user;
    char *cdir = getcwd(NULL,PATH_MAX);
    string rundir = "Run directory : " + (string) cdir;
   // vhistory.push_back(str);
    vhistory.push_back(rundir);
    vhistory.push_back("Parameter List START for " + (string) modulename);
    vhistory.push_back((string)getSerialNo (cnt)+" inputdatadir=" + (string) inputdatadir);
    vhistory.push_back((string)getSerialNo (cnt)+" caldbDir= " + (string) caldbDir);
    vhistory.push_back((string)getSerialNo (cnt)+" catalogpath="+(string)catalogpath);
    vhistory.push_back((string)getSerialNo (cnt)+" attitudefile="+(string)attitudefile);
    vhistory.push_back((string)getSerialNo (cnt)+" outdir= " + (string) outdir);
    vhistory.push_back((string)getSerialNo (cnt)+" Module Output directory=" + (string) moduleoutdir);
    vhistory.push_back((string)getSerialNo (cnt)+" Time col in attitude file = " + (string) att_timecol);
    vhistory.push_back((string)getSerialNo (cnt)+" Quaternion col in attitude file = " + (string) att_qcol);
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

int uvtFullFrameAst::readcatalogueFile() {
    
    return (EXIT_SUCCESS);
}


int uvtFullFrameAst::findStar_algo1(float *inputArray, float *expdata ,float peakOFexp) //algorithm for finding the peaks
{
    float mean_Ofimage=0.0;
     if (xsize == 0 || ysize == 0)
    {
        LOG(ERROR) << "***Divide by Zero***" << endl ;
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
     float *temp_array;
     vector<float> array_temp;     
     
      if(datainfo.getModeFlag ()==PC)
    {
       array_temp.clear ();
      for (int i=0;i<xsize*ysize;i++)
      {
         if(inputArray[i]!=0.0f && expdata[i]>(80*peakOFexp/100)){
          // if(inputArray[i]!=0.0f ){
            array_temp.push_back (inputArray[i]);
         }
         
            
      }
          temp_array  = new float[array_temp.size ()];
           for (int in=0;in<array_temp.size ();in++)
           {
               temp_array[in]=array_temp[in];     
           }
     }
      //}
     LOG(INFO)<<array_temp.size();
          
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
     double sd_temp=0.0f;
     
   if(datainfo.getModeFlag ()==PC)
    {
       LOG(INFO)<<"INSIDE PC mode";
        sd_temp=getSD (temp_array , array_temp.size ());
  mean_Ofimage=getmean (temp_array,array_temp.size ());
   thr =  sd_temp* sd_mul_factor ;
     
    }
    else
    {
       
         sd_temp=getSD (inputArray , xsize * ysize);
  mean_Ofimage=getmean (inputArray,xsize*ysize);
  thr =  sd_temp* sd_mul_factor ;
    }
    
     
  // sd_temp=getSD (inputArray , xsize * ysize);
  // mean_Ofimage=getmean (inputArray,xsize*ysize);
//    thr =  sd_temp* sd_mul_factor ;
//   thr = mean_Ofimage+sd_temp* sd_mul_factor ;
     thr=sd_temp* sd_mul_factor;
   LOG(ERROR) << endl << "\nThreshold for first cut peaks is   " << mean_Ofimage<<" + "<<sd_temp<<" X "<<sd_mul_factor<<" = "<<thr ;
   
    //Stores those  pixels whose va;ues are higher than 'thr'.
    for (int i = 0 ; i < xsize * ysize ; i++)
    {
        r = (i / xsize) ;
        c = (i % xsize) ;
       
             if (inputArray[i] > thr  )
                 {
            Fval.push_back (inputArray[i]) ;
            Fx.push_back (c) ; //x is for column
            Fy.push_back (r) ; //y is for row
                }
      
    }
   
   vector<int> FX1temp,FY1temp;
   vector<float> Fval1temp;
   bool flag_found=FALSE;
   LOG(INFO)<<Fx.size();
   for (int i=0;i<Fx.size();i++){
       flag_found=FALSE;
       for (int j=0;j<array_temp.size();j++){
           if(Fval[i]==array_temp[j]){
               flag_found=TRUE;
               break;
           }
       }
       if(flag_found==TRUE){
           FX1temp.push_back(Fx[i]);
           FY1temp.push_back(Fy[i]);
           Fval1temp.push_back(Fval[i]);
       }
       
   }
   
   Fx.clear();Fy.clear();Fval.clear();
   Fx=FX1temp;
   Fy=FY1temp;
   Fval=Fval1temp;
   LOG(INFO)<<Fval.size();
//   exit(1);
    LOG(INFO) << "SIGMA  Factor::" << sd_mul_factor << endl ;
    LOG(INFO) << " Size of First cut Peaks  " << Fy.size () << endl ;

    if (Fy.size () < minimum_No_of_Stars)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG(ERROR) << sd_mul_factor << " less than 0!!!! " << endl ;
            LOG(ERROR)<<"CRASH NO STAR FOUND (FIRST CUT LEVEL) TILL SIGMA MULTIPLIER BECAME ZERO OR LESS (uvtFullFrameAst.cpp)";
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
   LOG(INFO)<<refine_Winsize<<endl;
    if (refine_Winsize % 2 == 0)
        refine_Winsize = refine_Winsize - 1 ;

    LOG(INFO) << endl << "Using window size : " <<refine_Winsize << " for refining peaks " ;

    //refined peaks
    vector<int> Tx , Ty ;
    vector<float> Tval ;

    Tx = Fx ;
    Ty = Fy ;
    Tval = Fval ;

 Star1 star1;
 star_track.clear ();
 for(int i=0;i<Tx.size ();i++)
 {
     star1.x=Tx[i];
     star1.y=Ty[i];
     star1.intensity=Tval[i];
     star_track.push_back (star1);
     
 }
 sort (star_track.begin (),star_track.end (),compare);
    /*refining peaks logic
    refined Window size is for the refined  peaks.
   Refined peaks are found by  making window around each of the star(i.e first cut peaks)  and  finding brightest star among that window.*/
    int start_r , end_r , start_c , end_c ;
//to be removed 
   Tx.clear ();Ty.clear ();Tval.clear ();
  //  bool flag_unique=FALSE;
   LOG(INFO)<<star_track.size ()<<endl;
  // vector<Star> ::iterator itr =star_track.begin ();
    for (int i = 0 ; i < star_track.size () ; i ++)
    {
        start_r = star_track[i].y - refine_Winsize / 2 ;
        end_r = star_track[i].y + refine_Winsize / 2 ;
        start_c = star_track[i].x- refine_Winsize / 2 ;
        end_c = star_track[i].x + refine_Winsize / 2 ;
           if (start_r < 0) start_r = 0 ;
        if (end_r >= ysize) end_r = ysize - 1 ;
        if (start_c < 0) start_c = 0 ;
        if (end_c >= xsize) end_c = xsize - 1 ;
         for(int fcpeak=i+1;fcpeak<star_track.size ();fcpeak++)
         {
             if(star_track[fcpeak].x>start_c && star_track[fcpeak].x<end_c && star_track[fcpeak].y>start_r && star_track[fcpeak].y<end_r)
             {
                 
                 star_track.erase (star_track.begin ()+fcpeak);
                 fcpeak--;
             }
          }
        Tx.push_back (star_track[i].x);
        Ty.push_back (star_track[i].y);
        Tval.push_back (star_track[i].intensity);
    }
    
   

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
    for (int i = 0 ; i < xsize * ysize ; i++)
        arr_refine[i] = 0.0f ;

  
    for (int i = 0 ; i < Ty.size () ; i++)
        arr_refine[Ty[i] * xsize + Tx[i]] = Tval[i] ; //overwriting the same place..
    
    Tx.clear () ;
    Ty.clear () ;
    Tval.clear () ;

    for (int i = 0 ; i < xsize * ysize ; i++)
    {
        if (arr_refine[i] != 0)
        {
            Rx.push_back ((i % xsize)+1) ;
            Ry.push_back ((i / xsize)+1) ;
            Rval.push_back (arr_refine[i]) ;
        }
    }
   
    LOG(INFO) << "Number of final peaks is " << Rval.size () << endl ;

// for (int i=0;i<Rx.size();i++){
//       LOG(INFO)<<Rx[i]<<" "<<Ry[i];       
//   }
    if (Ry.size () < minimum_No_of_Stars)
    {
        sd_mul_factor = sd_mul_factor - 0.25 ;
        if (sd_mul_factor <= 0)
        {
            LOG(ERROR) << sd_mul_factor << " less than 0!!!! " << endl ;
	    LOG(ERROR)<<"CRASH NO STAR FOUND (REFINED PEAK) TILL SUGMA MULTIPLIER LESS THAN ZERO (uvtFullFrameAst.cpp)";
            return (EXIT_FAILURE) ;
        }
         delete[] arr_refine ;
        goto label ;
    }
    /**method for  find Centroid within the Stars**/
    doCentroiding (Rx , Ry , centroid_Winsize , inputArray , ysize , xsize) ;
    delete[] arr_refine ;
    return (EXIT_SUCCESS) ;
}

//method For doing Centroiding

void uvtFullFrameAst::doCentroiding (vector<int> &X , vector<int> &Y , int centroidwindow , float *arr , int h , int w)
{
    Cx.clear () ;
    Cy.clear () ;
    Ci.clear () ;
    double x_temp,y_temp,int_temp;
    Star1 star1;
    star_track.clear ();
     
    float x , y , val = 0 ;
     //if centroidwindowsize  is even, make it odd
    if (centroidwindow % 2 == 0)
        centroidwindow = centroidwindow - 1 ;

     LOG(INFO) << endl << "Using  window size : " <<centroidwindow<< " for finding Centroids " ;
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
//%#Added ON 20July#%
        if((Y[i] + j)>0 && (Y[i] + j)<h  &&  (X[i] + k)>0 && (X[i] + k)<w)
//%#-Till this-20July17#%
{       
		val = arr[(Y[i] + j) * w + (X[i] + k)] ;
                sum_x = sum_x + (x + k) * val ;
                sum_y = sum_y + (y + j) * val ;
                sum = sum + val ;
            }
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
//%#Editeded ON 20July#%
        x_temp= sum_x /  sum;
        y_temp= sum_y /  sum;
//%#-Till this-20July17#%
        int_temp=(float) sum;
        Cx.push_back (x_temp) ;
        Cy.push_back (y_temp) ;
        Ci.push_back (int_temp) ;
        star1.x=x_temp;
        star1.y=y_temp;
        star1.intensity=int_temp;
        star_track.push_back (star1);
    }
/**Sorting the list on the basis of the intensity**/
   sort (star_track.begin (),star_track.end (),compare);
    Cx.clear (),Cy.clear ();Ci.clear ();
   
    for (int i = 0 ; i < star_track.size () ; i++)
    {
        Cx.push_back (star_track[i].x);
        Cy.push_back (star_track[i].y);
        Ci.push_back (star_track[i].intensity);
    }

}
int uvtFullFrameAst::convertToUVITAxes (float &x ,float  &y,double &new_x,double &new_y,double &new_z)
{
    //cout<<tan ((NUV_COEFF/IMG_DIM_DI*HOUR_MACRO)*M_PI/180)<<endl;exit(1);
//    new_x=FOCAL_LENGTH*tan ((x-IMG_DIM_DI/2)*(PIX_PER_DEG_NUV*DEG_TO_RAD));
//    new_y =FOCAL_LENGTH*tan ((y-IMG_DIM_DI/2)*(PIX_PER_DEG_NUV*DEG_TO_RAD));
 double focal_length,detector_size;
 if(strcmp(datainfo.getDetector (),"VIS")==0)
 {
 focal_length=4.4877;
 detector_size=0.000025;
     
 }else if (strcmp(datainfo.getDetector (),"FUV")==0)
 {
 focal_length=4.6378;
 detector_size=0.000025;
 }
 else if (strcmp(datainfo.getDetector (),"NUV")==0)
 {
  focal_length=4.647;
 detector_size=0.000025;
 }
//  new_x=-1.0*(x-IMG_DIM_DI/2)*(sqrt(0.025)/(IMG_DIM_DI*100));
//  new_y =-1.0*(y-IMG_DIM_DI/2)*(sqrt(0.025)/(IMG_DIM_DI*100));
//  new_z=FOCAL_LENGTH;
  new_x=-1.0*(x-2400)*detector_size/8;
 new_y =-1.0*(y-2400)*detector_size/8;
 new_z=focal_length;
 //new_x=-1.0*(x-IMG_DIM_DI/2)*detector_size;
 //new_y =-1.0*(y-IMG_DIM_DI/2)*detector_size;
 //new_z=focal_length;
  //cout<<x-IMG_DIM_DI/2<<" "<<new_x<<" "<<" "<<y-IMG_DIM_DI/2<<" "<<new_y<<" "<<new_z<<endl;
return(EXIT_SUCCESS);
}
int uvtFullFrameAst::toNormalize (double &x, double &y, double &z,double &nor_x,double &nor_y,double &nor_z)
{
double rad=sqrt(x*x+y*y+z*z)    ;
nor_x=x/rad;
nor_y=y/rad;
nor_z=z/rad;
  //  cout<<"normalize::"<<nor_x<<" "<<nor_y<<" "<<nor_z<<endl;
    return(EXIT_SUCCESS);
}
string uvtFullFrameAst::getRaDECmatch(vector<float> &RA_img_Stars,vector<float> &DEC_img_Stars,int search_algo_ctlg,string len_a,string len_b,string rad_search,int no_of_stars,
        double &Abs_Val_DiffRA, double &Abs_Val_DiffDEC ,string nameFile,int *numStarsmatch,double decangle,bool flag_channel){
   
    
    double ra_temp=0.0f,dec_temp=0.0f;
    double ra_temp_currYear=0.0f,dec_temp_currYear=0.0f;
    LOG(INFO)<<rad_search;
    bool enough_points_found=FALSE;
    vector<string> usno_CAT2;
    vector<string> radeg_CAT2 ;
    vector<string>dedeg_CAT2;
    vector<string>bmag_CAT2;
    vector<string> rmag_CAT2;
    vector<string>epoch_CAT2;
   vector<string>nuvmag ;
    vector<string>errnuv;
    vector<string> fuvmag;
    vector<string> errfuv;
    vector<string> nuvcounts;
    vector<string> errcounts;
    vector<double> Catalogue_stars_RA;
    vector<double>Catalogue_stars_DEC;
        vector<double> diff_RA;
        vector<double> diff_DEC;
    
    char temp_ra[FLEN_FILENAME];
   char temp_dec[FLEN_FILENAME];
   double max_bmag=16;
   double rad_temp=0.0f;
   Database db;
  db.openDatabase((string)databasename);
   
 ofstream file;
 float diff_add_cata_ra,diff_add_cata_dec;
    while(enough_points_found==FALSE)
   {   
       file.open(nameFile.c_str(),ios::out| ios::trunc);
       file<<"RA_image"<<setw(20)<<"RA_catalogue"<<setw(20)<<"DEC_image"<<setw(20)<<"DEC_Catalogue"<<setw(20)<<"DIFF_RA"<<setw(20)<<"DIFF_DEC"<<endl;;
       file<<"====================================================================================================================================="<<endl;  
       
        usno_CAT2.clear(), radeg_CAT2.clear(), dedeg_CAT2.clear(), bmag_CAT2.clear(), rmag_CAT2.clear(), epoch_CAT2.clear(), nuvmag.clear(), errnuv.clear(), fuvmag.clear(), errfuv.clear(), nuvcounts.clear(), errcounts.clear();
        Catalogue_stars_RA.clear(), Catalogue_stars_DEC.clear();
        diff_RA.clear(), diff_DEC.clear();
        track_raback.clear();
        track_decback.clear();
        validStar_index_filtered.clear();
        valid_Star_index.clear();
     
        for (int i = 0; i < no_of_stars; i++) 
        {
     
                    usno_CAT2.clear(), radeg_CAT2.clear(), dedeg_CAT2.clear(), bmag_CAT2.clear(),
                    rmag_CAT2.clear(), epoch_CAT2.clear(), nuvmag.clear(), errnuv.clear(), fuvmag.clear(),
                    errfuv.clear(), nuvcounts.clear(), errcounts.clear();
            max_bmag = 16;
            ra_temp = -9999;
            dec_temp = -9999;
            sprintf(temp_ra, " %f", RA_img_Stars[i]);
            sprintf(temp_dec, " %f", DEC_img_Stars[i]);
          //  LOG(INFO)<<temp_ra<<" "<<temp_dec;
            db.select(usno_CAT2, radeg_CAT2, dedeg_CAT2, bmag_CAT2, rmag_CAT2, epoch_CAT2, (string) temp_ra, (string) temp_dec, nuvmag, errnuv, fuvmag, errfuv, nuvcounts, errcounts, search_algo_ctlg, len_a, len_b, rad_search,decangle);
           
            if(flag_channel==1){
                bmag_CAT2=nuvmag;
            }
           
            for (int j = 0; j < bmag_CAT2.size(); j++) 
            {
                
                if (atof(bmag_CAT2[j].c_str()) < max_bmag) 
             
                {
                    diff_add_cata_dec=0.0f;
                    diff_add_cata_ra=0.0f;
                    ra_temp = atof(radeg_CAT2[j].c_str());
                    dec_temp = atof(dedeg_CAT2[j].c_str());
            
                    for(int i=0;i<Catalogue_stars_RA.size();i++)
                    {
                    if(ra_temp==Catalogue_stars_RA[i] && dec_temp==Catalogue_stars_DEC[i])
                    {
                        ra_temp=-9999;dec_temp=-9999;                              
                        break;    
                    }
                    }
                    max_bmag = atof(bmag_CAT2[j].c_str());
                 
                }
            }
            

             if(ra_temp!=-9999 && dec_temp!=-9999)
             {
                 
                file<<temp_ra<<setw(20)<<setprecision(10)<<ra_temp <<setw(20)<<setprecision(10)<< temp_dec<<setw(20)<<dec_temp<<setw(20)<<(ra_temp-RA_img_Stars[i])*decangle<<setw(20)<<dec_temp-DEC_img_Stars[i]<<endl;
                 
             }
            track_raback.push_back(ra_temp);
            track_decback.push_back(dec_temp);
            Catalogue_stars_RA.push_back(ra_temp);
            Catalogue_stars_DEC.push_back(dec_temp);

        }
        
        
   
       for (int k=0;k<Catalogue_stars_RA.size ();k++)      
        {
           if(Catalogue_stars_RA[k]!=-9999  && Catalogue_stars_DEC[k]!=-9999)
           {
               valid_Star_index.push_back(k);              
               diff_RA.push_back ((Catalogue_stars_RA[k]-RA_img_Stars[k])*decangle);//multily by cos DEC_PNT
               diff_DEC.push_back (Catalogue_stars_DEC[k]-DEC_img_Stars[k]);     

       }

       
       }
      
        double *temp_Array_RA = new double[diff_RA.size()];
        double *temp_Array_DEC = new double[diff_RA.size()];
        vector<double> diffRAbackup =diff_RA;
    
        if(diff_RA.size()>0){
       int status=getNearByValues(diff_RA.data(),temp_Array_RA,diff_RA.size(),0.00138888);
        if(status)
        {
        LOG(ERROR)<<"Error in finding the near by values for RA";
        exit(1);
        }
     
      status=getNearByValues(diff_DEC.data(),temp_Array_DEC,diff_DEC.size(),0.00138888);
        if(status)
        {
        LOG(ERROR)<<"Error in finding the near by values for DEC";
        exit(1);
        }

        
        vector<double > diffRA_matched,diffDEC_matched;
        for (int i=0;i<diff_RA.size();i++)
        {
            if(temp_Array_RA[i]==1 && temp_Array_DEC[i]==1){
               
                diffRA_matched.push_back(diff_RA[i]);
                diffDEC_matched.push_back(diff_DEC[i]);
                validStar_index_filtered.push_back(valid_Star_index[i]);
            }

        }
      
        diff_RA.clear();
        diff_DEC.clear();
        diff_RA=diffRA_matched;
        diff_DEC=diffDEC_matched;
//         LOG(INFO)<<diff_RA.size();
        
        *numStarsmatch=diff_RA.size();
        
        
        }
        if(diff_RA.size()>=2)            
        {
           
           enough_points_found=TRUE;
        }
        else 
        {
       rad_temp=atof(rad_search.c_str ());
       rad_temp=rad_temp*1.2;
       sprintf((char*)rad_search.c_str (),"%f",rad_temp);
       if(rad_temp>1.2){
           LOG(ERROR)<<"***Search radius become greater than 3 arc-min***";
	   LOG(ERROR)<<"CRASH SEARCH RADIUS BECAME LARGER THAN 3 ARC-MIN (uvtFullFrameAst.cpp)";
           break;
       }
       file.close();
        
       }
 
        
   }
   //added
   
   if(enough_points_found==FALSE)
   {
       LOG(INFO)<<"Error in getting points for RA and DEC matching";
       flag_NOT_FOUND_CATA_MATCH=1;
   }
   else{
   //till this
      
       float ra_max=0;
       float dec_max=0;
       for (int i=0;i<diff_DEC.size();i++)
       {
        if(ra_max<abs(diff_RA[i])){
            ra_max=abs(diff_RA[i]);
        }
        if(dec_max<abs(diff_DEC[i]))
        {
            dec_max=abs(diff_DEC[i]);
        }
       
       }
 
       Abs_Val_DiffRA=ra_max;//change name of variable that implies diffrence
      Abs_Val_DiffDEC=dec_max;// ||
     
  
   double new_searchRadius=sqrt(Abs_Val_DiffRA*Abs_Val_DiffRA+Abs_Val_DiffDEC*Abs_Val_DiffDEC);//*1.2;
 
   Abs_Val_DiffRA=getmean(diff_RA.data(),diff_RA.size());
   Abs_Val_DiffDEC=getmean(diff_DEC.data(),diff_DEC.size());
 
       file<<"====================================================================="<<endl;
       file<<"Abs max DIFF_RA value "<<Abs_Val_DiffRA<<endl;
       file<<"Abs max DIFF_DEC value "<<Abs_Val_DiffDEC<<endl;
  
   rad_search=convertFloatToStr(new_searchRadius);

   }
  file.close();
   return rad_search;
    
}

int uvtFullFrameAst::calculateShiftsAndRoll( int totalelements,vector<double> &Xone,vector<double> &Yone,vector<double> &Xtwo,vector<double> &Ytwo,int xsize,int ysize,double &Xdx,double  &Ydy,double &Theta)
{
        spMatrix B ((totalelements) * 2 , 1) ;
        spMatrix A ((totalelements) * 2 , 3) ;
        spMatrix X (3 , 1) ;
        int temp = 0 ;

        for (int t = 0 ; t < totalelements * 2 ; t = t + 2)
        {
            B (t , 0) = (Xtwo[temp]/8 - Xone[temp]/8) ; //+IMAGE_ARRAYSIZE*0.5;
            B (t + 1 , 0) = (Ytwo[temp]/8 - Yone[temp]/8) ; //+IMAGE_ARRAYSIZE*0.5;

            A (t , 0) = - 1.0 * (Yone[temp]/8 - 300) ;
            A (t , 1) = 1 ;
            A (t , 2) = 0 ;

            A (t + 1 , 0) = (Xone[temp]/8 - 300) ;
            A (t + 1 , 1) = 0 ;
            A (t + 1 , 2) = 1 ;

            temp ++ ;
        }

        X.ApplyLeastSquare (A , B) ;

        Xdx = X (1 , 0) ;
        Ydy = X (2 , 0) ;
        Theta = X (0 , 0)*180/M_PI ;
        return (EXIT_SUCCESS);
}

int writeUsrkeywordsFrmvectNew(char *filename,vector<string> &strvect){
            
    fitsfile *fout;
            int status=0;
            int numhdu=0;
             fits_open_file (&fout , filename , READWRITE , &status) ;
              printError (status , "Error in opening the input File" , filename) ;
           // LOG(INFO)<<"Successfully opened"<<endl;
               fits_get_num_hdus (fout , &numhdu , &status) ;
                               
               for(int hno=1;hno<=1;hno++){
                    fits_movabs_hdu (fout , hno , NULL , &status) ;

        for (int i = 0; i < strvect.size(); i++) {
            if (strstr(strvect[i].c_str(), "NAXIS") == NULL && strstr(strvect[i].c_str () ,"BITPIX")==NULL && strstr(strvect[i].c_str () ,"NAXIS1")==NULL 
                     && strstr(strvect[i].c_str () ,"NAXES2")==NULL   && strstr(strvect[i].c_str(),"SIMPLE")==NULL  && strstr(strvect[i].c_str(),"EXTEND")==NULL  ) {
                fits_write_record(fout, strvect[i].c_str(), &status);


                //LOG(INFO)<<strvect[i].c_str ()<<endl;
                if (status) {
                    LOG(ERROR) << endl << "***Error in writing keyword " << strvect[i].c_str() << "***";
                    fits_report_error(stderr, status);
                    fits_close_file(fout, &status);
                    printError(status, "Error in closing the input File", filename);
                    //     LOG(INFO)<<"Successfully closed with error"<<endl;
                    return (EXIT_FAILURE);
                }
            }
        }
                           
               }
               
            fits_close_file(fout,&status);
            printError (status , "Error in closing the input File" , filename) ;
           // LOG(INFO)<<"Successfully closed"<<endl;
//            fits_open_file (&fout , filename , READWRITE , &status) ;
//              printError (status , "Error in opening the input File" , filename) ;
//                fits_close_file(fout,&status);
//            printError (status , "Error in closing the input File" , filename) ;
            
            return(EXIT_SUCCESS);
        }
