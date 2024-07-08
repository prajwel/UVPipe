/* 
 * File:   caldb_Handler.cpp
 * Authors:: Preeti Tahlani ,Dhruv Shelat,Sanjay K Singh, Arvind K Singh
 * Scientific Contributers: Swarna k Ghosh, Shyam.N Tandon
 *
 * Created on September 2, 2013, 3:26 PM
 */
#include <vector>
#include <uvtUtils.h>
#include <caldb_Handler.h>
#include<glog/logging.h>
using namespace std;


 string  caldb_Handler:: getFlatFieldFile(char *detector,char *mode,char *filter, char *path)
{
     if(!DirExists(path)){
        LOG(ERROR)<<endl<<"CALDB path  "<<path<<"  does not exist"<<endl;
        return NULL;
    }
     //char* tpath= new char[400];
  char tpath[400];
    char tpath1[400];
    char tfitpath[500];
   
    sprintf(tpath,"%s/%s/%s/%s/%s/",path,CALDB_FLATFIELD_DIR,detector,mode,filter);
//    strcpy(tpath1,tpath);
//       const char *fitpath= searchFile(tpath1,".fits");
//  
//   sprintf(tfitpath,"%s",fitpath);
//    g
       joinStrings (path,1,tpath);
      string fitpath= searchFile(path,".fits");
      if(strcmp(fitpath.c_str(),"")==0)
          return NULL;
      joinStrings (tfitpath,1,fitpath.c_str());
//    delete [] tpath;
    return  tfitpath;
}
// string caldb_Handler:: getCatalogue_LookupFile(char *path)
//{
//     if(!DirExists(path)){
//        LOG(ERROR)<<endl<<"CALDB path  "<<path<<"  does not exist"<<endl;
//        return NULL;
//    }
//     //char* tpath= new char[400];
//  char tpath[400];
//    char tpath1[400];
//    char *tfitpath = new char[500];
//   
//    sprintf(tpath,"%s/%s/",path,CALDB_CATALOGUEFILE_PATH);
////    strcpy(tpath1,tpath);
////       const char *fitpath= searchFile(tpath1,".fits");
////  
////   sprintf(tfitpath,"%s",fitpath);
////    g
//       joinStrings (path,1,tpath);
//      const char *fitpath= searchFile(path,".fits");
//     joinStrings (tfitpath,1,(char*)fitpath);
////    delete [] tpath;
//     LOG(INFO)<<tfitpath;
//    return  (string)tfitpath;
//}

string  caldb_Handler::getBadPixelFile(char *detector,char *mode,int xsize, int ysize, char *path)
{
    if(!DirExists(path))
    {
        LOG(ERROR)<<endl<<"***CALDB path  "<<path<<"  does not exist***"<<endl;
        exit(EXIT_FAILURE);
    }
    char tpath[FLEN_FILENAME];
    char tfitpath[FLEN_FILENAME]=" ";
    char dimensiondir[FLEN_FILENAME];
   
    if (xsize == 100 && ysize == 100) strcpy (dimensiondir , "100X100") ;
//    else if (xsize == 100 && ysize == 150) strcpy (dimensiondir , "100X150") ;
//    else if (xsize == 100 && ysize == 200) strcpy (dimensiondir , "100X200") ;
//    else if (xsize == 100 && ysize == 250) strcpy (dimensiondir , "100X250") ;
//    else if (xsize == 100 && ysize == 300) strcpy (dimensiondir , "100X300") ;
//    else if (xsize == 100 && ysize == 350) strcpy (dimensiondir , "100X350") ;
//    else if (xsize == 100 && ysize == 512) strcpy (dimensiondir , "100X512") ;

   // else if (xsize == 150 && ysize == 100) strcpy (dimensiondir , "150X100") ;
    else if (xsize == 150 && ysize == 150) strcpy (dimensiondir , "150X150") ;
//    else if (xsize == 150 && ysize == 200) strcpy (dimensiondir , "150X200") ;
//    else if (xsize == 150 && ysize == 250) strcpy (dimensiondir , "150X250") ;
//    else if (xsize == 150 && ysize == 300) strcpy (dimensiondir , "150X300") ;
//    else if (xsize == 150 && ysize == 350) strcpy (dimensiondir , "150X350") ;
//    else if (xsize == 150 && ysize == 512) strcpy (dimensiondir , "150X512") ;

//    else if (xsize == 512 && ysize == 100) strcpy (dimensiondir , "512X100") ;
//    else if (xsize == 512 && ysize == 150) strcpy (dimensiondir , "512X150") ;
//    else if (xsize == 512 && ysize == 200) strcpy (dimensiondir , "512X200") ;
//    else if (xsize == 512 && ysize == 250) strcpy (dimensiondir , "512X250") ;
//    else if (xsize == 512 && ysize == 300) strcpy (dimensiondir , "512X300") ;
//    else if (xsize == 512 && ysize == 350) strcpy (dimensiondir , "512X350") ;
    else if (xsize == 512 && ysize == 512) strcpy (dimensiondir , "512X512") ;
    else if(xsize==200 && ysize==200 )  strcpy (dimensiondir,"200X200");
    else if(xsize==250 && ysize==250 )  strcpy (dimensiondir,"250X250");
    else if(xsize==300 && ysize==300 )  strcpy (dimensiondir,"300X300");
    else if(xsize==350 && ysize==350 )  strcpy (dimensiondir,"350X350");
else{
        LOG(ERROR)<<"***No directory found for  xsize :"<<xsize<<"  ysize: "<<ysize<<"***"<<endl;
        return tfitpath;
    }
 
    sprintf(tpath,"%s/%s/%s/%s/%s/",path,CALDB_BADPIXEL_DIR,detector,mode,dimensiondir);
    string fitpath= searchFile(tpath,".fits");
    
    if(strcmp(fitpath.c_str(),"")==0){
        LOG(ERROR)<<endl<<"***No directory found for  xsize :"<<xsize<<"  ysize: "<<ysize<<"***"<<endl;
        return tfitpath;
    }
    sprintf(tfitpath,"%s/%s",tpath,fitpath.c_str());
    return  tfitpath;
       
}

string caldb_Handler:: getDetectorFile(char *detector, char *path)
{
     if(!DirExists(path)){
        LOG(ERROR)<<"CALDB path  "<<path<<"  does not exist";
        return NULL;
    }
    char tpath[100];
    char tfitpath[200]=" ";
    sprintf(tpath,"%s/%s/%s/",path,CALDB_DETECTORDISTO_DIR,detector);
    strcpy(path,tpath);
 
    string fitpath= searchFile(path,".fits");
     if(strcmp(fitpath.c_str(),"")==0){
        LOG(ERROR)<<endl<<"***caldb searching failed for detector distortion***"<<endl;
        return NULL;
    }
   // LOG(INFO)<<fitpath<<endl;
    sprintf(tfitpath,"%s",fitpath.c_str());
    return  fitpath;
    //return  tfitpath;
}
string  caldb_Handler:: getOpticalDetectorFile(char *detector, char *path,char *filter)
{
     if(!DirExists(path)){
        LOG(ERROR)<<"CALDB path  "<<path<<"  does not exist";

    }
    char tpath[100];
    char tfitpath[200]="empty";
    sprintf(tpath,"%s/%s/%s/%s/",path,CALDB_OPTICALDISTO_DIR,detector,filter);
    strcpy(path,tpath);
 
 string fitpath= searchFile(path,".fits");
    if(strcmp(fitpath.c_str(),"")==0){
     LOG(ERROR)<<endl<<"***caldb searching failed for optical distortion***"<<endl;
        return NULL;
    }
    sprintf(tfitpath,"%s",fitpath.c_str());
    
    //return  tfitpath;
    return  fitpath;
}

string caldb_Handler:: getCentroidBiasFile(char* detector, char* path, int winsize)
{
    char windir[10];
    if(winsize==3)
        strcpy(windir,"3X3");
    else if(winsize==5)
        strcpy(windir,"5X5");
    else{
        LOG(ERROR)<<"***Invalid Value For Window Size***";
        return NULL;
    }        
        
    if(!DirExists(path)){
        LOG(ERROR)<<"CALDB path  "<<path<<"  does not exist";
        return NULL;
    }
    char tpath[400];
    char tfitpath[400]=" ";
   
    sprintf(tpath,"%s/%s/%s/%s/",path,CALDB_CENTROIDBIAS_DIR,detector,windir);
    strcpy(path,tpath);
 
   string fitpath= searchFile(path,".fits");
   
   if(strcmp(fitpath.c_str(),"")==0){
     LOG(ERROR)<<endl<<"***caldb searching failed for centroid bias***"<<endl;
        return NULL;
    }
    
    sprintf(tfitpath,"%s",fitpath.c_str());
   // LOG(INFO)<<endl<<"The caldb centroid bias file is  "<<tfitpath;
    return  tfitpath;
}
string  caldb_Handler:: getCentroidEAFile(char *detector, char *path)
{
     if(!DirExists(path))
    {
        LOG(ERROR)<<"CALDB path  "<<path<<"  does not exist";
        return NULL;
    }
    char tpath[100];
    char tfitpath[200]=" ";
    sprintf(tpath,"%s/%s/%s/",path,CALDB_CENTROIDCORR_DIR,detector);
    strcpy(path,tpath);
 
    string fitpath= searchFile(path,".fits");
    if(strcmp(fitpath.c_str(),"")==0){
     LOG(ERROR)<<endl<<"***caldb searching failed for EA file ***"<<endl;
        return NULL;
    }
    sprintf(tfitpath,"%s",fitpath.c_str());
     return  fitpath;
}


string caldb_Handler:: getQEFile(char* detector, char* mod, char* path)
{
     if(!DirExists(path)){
        LOG(ERROR)<<"CALDB path  "<<path<<"  does not exist";
        return NULL;
    }
    char tpath[100];
    char tfitpath[200]="";
    sprintf(tpath,"%s/%s/%s/%s/",path,CALDB_QE_DIR,detector,mod);
    strcpy(path,tpath);
 
    string fitpath= searchFile(path,".fits");
    if(strcmp(fitpath.c_str(),"")==0){
     LOG(ERROR)<<endl<<"***caldb searching failed for QE file***"<<endl;
        return NULL;
    }
    sprintf(tfitpath,"%s",fitpath.c_str());
   
    //return  tfitpath;
    return  fitpath;
}
string  caldb_Handler:: getThermalFile(char *detector, char *path)
{
     if(!DirExists(path)){
        LOG(ERROR)<<"CALDB path  "<<path<<"  does not exist";
        return NULL;
    }
    char tpath[100];
    char tfitpath[200]="";
    sprintf(tpath,"%s/%s/%s/",path,CALDB_THERMAL_DIR,detector);

    strcpy(path,tpath);
 
string  fitpath= searchFile(path,".fits");
    if(strcmp(fitpath.c_str(),"")==0){
     LOG(ERROR)<<endl<<"***caldb searching failed for Thermal file***"<<endl;
        return NULL;
    }
    //LOG(INFO)<<fitpath;
    sprintf(tfitpath,"%s",fitpath.c_str());
  
    return  tfitpath;
}

int caldb_Handler::readCaldbDistFile(float *Xdist, float *Ydist, char *distortionCorrfile) 
{
     long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    fitsfile *fdist;
    int status = 0;
    long dim;
     fits_open_file(&fdist, distortionCorrfile, READONLY, &status);
    printError(status, "Error in opening the distortion correction file",distortionCorrfile);
    fits_movabs_hdu(fdist, 2, NULL, &status);
    printError(status, "Moving to 2nd HDU in calDB file Fails",distortionCorrfile);
//    fits_get_num_rows(fdist, &dim, &status);
//    printError(status, "Number of Rows counting Fails",distortionCorrfile);
    
//    if(dim!=CALDB_DIST_SIZE){
//        LOG(ERROR)<<"\nCALDB distortion size is not equal to  "<<CALDB_DIST_SIZE;
//        return (EXIT_FAILURE);
//    }
//    
   status = 0;
   fits_read_pix (fdist , TFLOAT , fpixel , CALDB_DIST_SIZE*CALDB_DIST_SIZE , NULL , Xdist , NULL , &status) ;
   printError (status , "Error in reading the pixels from input signal Frame " , distortionCorrfile) ;
  //  int index = 0;
    
//    for (int i = 0; i < dim; i++) {
//        for (int j = 0; j < dim; j++) {
//            fits_read_col(fdist, TFLOAT, j + 1, i + 1, 1, 1, NULL, &Xdist[index], NULL, &status);
//            printError(status, "Reading a column Fails in calDB",distortionCorrfile);
//            Xdist[index]=0.0f;
//            index++;
//         }   
//    }
    //  index = 0;
    fits_movabs_hdu(fdist, 3, NULL, &status);
    printError(status, "Moving to 3rd HDU in calDB file Fails"); //for distortion in y
      fits_read_pix (fdist , TFLOAT , fpixel , CALDB_DIST_SIZE*CALDB_DIST_SIZE , NULL , Ydist , NULL , &status) ;
   printError (status , "Error in reading the pixels from input signal Frame " , distortionCorrfile) ;
   
//    fits_get_num_rows(fdist, &dim, &status);
//    printError(status, "Number of Rows counting Fails",distortionCorrfile); //getting no of rows in the perticuler HDU
//
//    for (int i = 0; i < dim; i++) {
//        for (int j = 0; j < dim; j++) {
//            fits_read_col(fdist, TFLOAT, j + 1, i + 1, 1, 1, NULL, & Ydist[index], NULL, &status);
//            printError(status, "Reading a column Fails in calDB ",distortionCorrfile);
//            Ydist[index]=0.0f;
//            index++;
//            }
//         }
       fits_close_file(fdist, &status);
       printError(status, "Error in closing the file  ",distortionCorrfile);
      LOG(INFO) << "Reading of  Distortion Correction file from CALDB completed"<<endl;
      return (EXIT_SUCCESS);
}
int caldb_Handler::readCaldbOpticDistFile(float *Xdist, float *Ydist, char *distortionCorrfile) 
{
     long fpixel[2] ;
    fpixel[0] = fpixel[1] = 1 ;
    fitsfile *fdist;
    int status = 0;
    long dim;
     fits_open_file(&fdist, distortionCorrfile, READONLY, &status);
    printError(status, "Error in opening the distortion correction file",distortionCorrfile);
    fits_movabs_hdu(fdist, 2, NULL, &status);
    printError(status, "Moving to 2nd HDU in calDB file Fails",distortionCorrfile);
    fits_get_num_rows(fdist, &dim, &status);
    printError(status, "Number of Rows counting Fails",distortionCorrfile);
    
    if(dim!=CALDB_DIST_SIZE){
        LOG(ERROR)<<"\nCALDB distortion size is not equal to  "<<CALDB_DIST_SIZE;
        return (EXIT_FAILURE);
    }
    
  // status = 0;
  // fits_read_pix (fdist , TFLOAT , fpixel , CALDB_DIST_SIZE*CALDB_DIST_SIZE , NULL , Xdist , NULL , &status) ;
  // printError (status , "Error in reading the pixels from input signal Frame " , distortionCorrfile) ;
   int index = 0;
    
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            fits_read_col(fdist, TFLOAT, j + 1, i + 1, 1, 1, NULL, &Xdist[index], NULL, &status);
            printError(status, "Reading a column Fails in calDB",distortionCorrfile);
         //   Xdist[index]=0.0f;
            index++;
         }   
    }
     index = 0;
    fits_movabs_hdu(fdist, 3, NULL, &status);
    printError(status, "Moving to 3rd HDU in calDB file Fails"); //for distortion in y
//      fits_read_pix (fdist , TFLOAT , fpixel , CALDB_DIST_SIZE*CALDB_DIST_SIZE , NULL , Ydist , NULL , &status) ;
//   printError (status , "Error in reading the pixels from input signal Frame " , distortionCorrfile) ;
    fits_get_num_rows(fdist, &dim, &status);
    printError(status, "Number of Rows counting Fails",distortionCorrfile); //getting no of rows in the perticuler HDU

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            fits_read_col(fdist, TFLOAT , j + 1, i + 1, 1, 1, NULL, & Ydist[index], NULL, &status);
            printError(status, "Reading a column Fails in calDB ",distortionCorrfile);
           // Ydist[index]=0.0f;
            index++;
            }
         }
       fits_close_file(fdist, &status);
       printError(status, "Error in closing the file  ",distortionCorrfile);
      LOG(INFO) << "Reading of  Distortion Correction file from CALDB completed"<<endl;
      return (EXIT_SUCCESS);
}

char* caldb_Handler::getTelDefFile(char *path ,char *detector,char *filter)
{ 
     if(!DirExists(path)){
        LOG(ERROR)<<endl<<"CALDB path  "<<path<<"  does not exist"<<endl;
        return NULL;
    }
     
     char filepath[FLEN_FILENAME]="";
     sprintf(filepath,"%s/%s/%s/%s",path,CALDB_TELDEF_DIR,detector,filter);
   //  LOG(INFO)<<"File path "<<endl;
     string teldeffile= searchFile(filepath,".fits");
     if(strcmp(teldeffile.c_str(),"")==0){
     LOG(ERROR)<<endl<<"***caldb searching failed for teldef file ***"<<endl;
        return NULL;
    }
     char teldefpath[FLEN_FILENAME]="";
     //sprintf(teldefpath,"%s/%s",filepath,teldeffile);
     
     char finalpath[FLEN_FILENAME]="";
     strcpy(finalpath,teldefpath);
   
     return finalpath;
}

 char* caldb_Handler :: getPlatScaleFile(char *path ,char *detector)
{
     if(!DirExists(path)){
        LOG(ERROR)<<endl<<"CALDB path  "<<path<<"  does not exist"<<endl;
        return NULL;
    }
     
    char tpath[100];
    char tfitpath[200]="";
    sprintf(tpath,"%s/%s/%s/",path,CALDB_PLATESCALE_DIR,detector);

    strcpy(path,tpath);
 
    string fitpath= searchFile(path,".fits");
    sprintf(tfitpath,"%s", fitpath.c_str());
    return  tfitpath;
}

 int TelDef::read(char* teldeffile)
 {
       
   int flag = readKeywords(teldeffile,1,29,TFLOAT,"COE_X0_A",&coef_x0_a,
                                                TFLOAT,"COE_X0_B",&coef_x0_b,
                                                TFLOAT,"COE_X0_C",&coef_x0_c,
                                                TFLOAT,"COE_Y0_A",&coef_y0_a,
                                                TFLOAT,"COE_Y0_B",&coef_y0_b,
                                                TFLOAT,"COE_Y0_C",&coef_y0_c,
                                                TFLOAT,"DETXPIX1",&det_xpix1,
                                                TFLOAT,"DETYPIX1",&det_ypix1,
                                                TFLOAT,"DET_XSIZ",&det_xsiz,
                                                TFLOAT,"DET_YSIZ",&det_ysiz,
                                                TFLOAT,"DET_XSCL",&det_xscl,
                                                TFLOAT,"DET_YSCL",&det_yscl,
                                                TFLOAT,"INT_XCEN",&int_xcen,
                                                TFLOAT,"INT_YCEN",&int_ycen,
                                                TINT,"DETXFLIP",&detxflip,
                                                TINT,"DETYFLIP",&detyflip,
                                                TFLOAT,"DET_XOFF",&det_xoff,
                                                TFLOAT,"DET_YOFF",&det_yoff,
                                                TFLOAT,"DET_SCAL",&det_scal,
                                                TFLOAT,"FOCALLEN",&focal_len,
                                                TFLOAT,"ALIGNM11",&m11,
                                                TFLOAT,"ALIGNM12",&m12,
                                                TFLOAT,"ALIGNM13",&m13,
                                                TFLOAT,"ALIGNM21",&m21,
                                                TFLOAT,"ALIGNM22",&m22,
                                                TFLOAT,"ALIGNM23",&m23,
                                                TFLOAT,"ALIGNM31",&m31,
                                                TFLOAT,"ALIGNM32",&m32,
                                                TFLOAT,"ALIGNM33",&m33);
   
   if(flag)  return (EXIT_FAILURE);
   
   det_xcen = det_xpix1+ (det_xsiz-1)/2.0;
   det_ycen = det_ypix1 + (det_ysiz-1)/2.0;
   return (EXIT_SUCCESS);
   
 }
 string  caldb_Handler::getTemplateFileForExposure(char *detector,char *mode,int xsize, int ysize, char *path)
{
    if(!DirExists(path))
    {
        LOG(ERROR)<<endl<<"***CALDB path  "<<path<<"  does not exist***"<<endl;
        exit(EXIT_FAILURE);
    }
    char tpath[FLEN_FILENAME];
    char tfitpath[FLEN_FILENAME]=" ";
    char dimensiondir[FLEN_FILENAME];
   
//    if (xsize == 100 && ysize == 100) strcpy (dimensiondir , "100X100") ;
//    else if (xsize == 100 && ysize == 150) strcpy (dimensiondir , "100X150") ;
//    else if (xsize == 100 && ysize == 200) strcpy (dimensiondir , "100X200") ;
//    else if (xsize == 100 && ysize == 250) strcpy (dimensiondir , "100X250") ;
//    else if (xsize == 100 && ysize == 300) strcpy (dimensiondir , "100X300") ;
//    else if (xsize == 100 && ysize == 350) strcpy (dimensiondir , "100X350") ;
//    else if (xsize == 100 && ysize == 512) strcpy (dimensiondir , "100X512") ;
//
//    else if (xsize == 150 && ysize == 100) strcpy (dimensiondir , "150X100") ;
//    else if (xsize == 150 && ysize == 150) strcpy (dimensiondir , "150X150") ;
//    else if (xsize == 150 && ysize == 200) strcpy (dimensiondir , "150X200") ;
//    else if (xsize == 150 && ysize == 250) strcpy (dimensiondir , "150X250") ;
//    else if (xsize == 150 && ysize == 300) strcpy (dimensiondir , "150X300") ;
//    else if (xsize == 150 && ysize == 350) strcpy (dimensiondir , "150X350") ;
//    else if (xsize == 150 && ysize == 512) strcpy (dimensiondir , "150X512") ;
//
//    else if (xsize == 512 && ysize == 100) strcpy (dimensiondir , "512X100") ;
//    else if (xsize == 512 && ysize == 150) strcpy (dimensiondir , "512X150") ;
//    else if (xsize == 512 && ysize == 200) strcpy (dimensiondir , "512X200") ;
//    else if (xsize == 512 && ysize == 250) strcpy (dimensiondir , "512X250") ;
//    else if (xsize == 512 && ysize == 300) strcpy (dimensiondir , "512X300") ;
//    else if (xsize == 512 && ysize == 350) strcpy (dimensiondir , "512X350") ;
//    else if (xsize == 512 && ysize == 512) strcpy (dimensiondir , "512X512") ;
//    else if(xsize==200 && ysize==200 )  strcpy (dimensiondir,"200X200");
    
    if (xsize == 100 && ysize == 100) strcpy (dimensiondir , "100X100") ;
//    else if (xsize == 100 && ysize == 150) strcpy (dimensiondir , "100X150") ;
//    else if (xsize == 100 && ysize == 200) strcpy (dimensiondir , "100X200") ;
//    else if (xsize == 100 && ysize == 250) strcpy (dimensiondir , "100X250") ;
//    else if (xsize == 100 && ysize == 300) strcpy (dimensiondir , "100X300") ;
//    else if (xsize == 100 && ysize == 350) strcpy (dimensiondir , "100X350") ;
//    else if (xsize == 100 && ysize == 512) strcpy (dimensiondir , "100X512") ;

   // else if (xsize == 150 && ysize == 100) strcpy (dimensiondir , "150X100") ;
    else if (xsize == 150 && ysize == 150) strcpy (dimensiondir , "150X150") ;
//    else if (xsize == 150 && ysize == 200) strcpy (dimensiondir , "150X200") ;
//    else if (xsize == 150 && ysize == 250) strcpy (dimensiondir , "150X250") ;
//    else if (xsize == 150 && ysize == 300) strcpy (dimensiondir , "150X300") ;
//    else if (xsize == 150 && ysize == 350) strcpy (dimensiondir , "150X350") ;
//    else if (xsize == 150 && ysize == 512) strcpy (dimensiondir , "150X512") ;

//    else if (xsize == 512 && ysize == 100) strcpy (dimensiondir , "512X100") ;
//    else if (xsize == 512 && ysize == 150) strcpy (dimensiondir , "512X150") ;
//    else if (xsize == 512 && ysize == 200) strcpy (dimensiondir , "512X200") ;
//    else if (xsize == 512 && ysize == 250) strcpy (dimensiondir , "512X250") ;
//    else if (xsize == 512 && ysize == 300) strcpy (dimensiondir , "512X300") ;
//    else if (xsize == 512 && ysize == 350) strcpy (dimensiondir , "512X350") ;
    else if (xsize == 512 && ysize == 512) strcpy (dimensiondir , "512X512") ;
    else if(xsize==200 && ysize==200 )  strcpy (dimensiondir,"200X200");
    else if(xsize==250 && ysize==250 )  strcpy (dimensiondir,"250X250");
    else if(xsize==300 && ysize==300 )  strcpy (dimensiondir,"300X300");
    else if(xsize==350 && ysize==350 )  strcpy (dimensiondir,"350X350");
    
else{
        LOG(ERROR)<<"***No directory found for  xsize :"<<xsize<<"  ysize: "<<ysize<<"***"<<endl;
        exit(EXIT_FAILURE);
    }
 
    sprintf(tpath,"%s/%s/%s/%s/%s/",path,CALDB_TEMPLATE_EXP_DIR,detector,mode,dimensiondir);
    string fitpath= searchFile(tpath,".fits");
    
    if(strcmp(fitpath.c_str(),"")==0){
        LOG(ERROR)<<endl<<"***No directory found for  xsize :"<<xsize<<"  ysize: "<<ysize<<"***"<<endl;
        return NULL;
    }
    sprintf(tfitpath,"%s/%s",tpath,fitpath.c_str());
    return  tfitpath;
       
}
