
#include<fitsio.h>
#include<gti.h>
#include<macro_def.h>
#include<iostream>
//#include<uvtUtils.h>
#include<fstream>
#include <string.h>

using namespace std;

struct gti_paraminfo
{
    string  param_name ;
    string   value_paramStr ;
    int  datatype_param;
} ;
int checkValidRange(vector<gti_paraminfo> &gtiInfo ,vector<string> allowedval);
int genL2gti(char *level1Sciencedatafile ,char *gtiparamfiles, char *mkffile ,char *hkfile,vector<unsigned char> &new_gtiFlag)
{
    int status =0;
    long nrow_l1=0;
    vector<float> time_final;
    int numhdu_mkf,numhdu_hk,keyexist;
    int nelement=1,felement=1;
    int  term_datatype;//[FLEN_FILENAME];
    //string  record;//[FLEN_CARD] ;
    char record[FLEN_FILENAME];
    int temp_num;
    string record_final;
     ifstream fin(gtiparamfiles);

    cout<<"File opened "<<gtiparamfiles<<endl;
    int i;
    unsigned short bitPos;
    string temp_line,paramnames,acceptCriteria;
    vector<string> param_vect,param_allowedValues,finalparam_allowedValues;
  
    fin>>temp_line>>temp_line>>temp_line>>temp_line;
    
    //reading the parameter name,bit positions and allowed values.
    while (! fin.eof ())
    {
        //fin >> num ;
        fin >>paramnames >> bitPos >> acceptCriteria ;
       // LOG(INFO)<< "outputs " << P.name << P.bitpos << P.acceptCriteria ;
        if (fin.eof ()) break ;
        param_vect.push_back (paramnames) ;
        param_allowedValues.push_back (acceptCriteria);
    }
    fin.close();
    
    
    fitsfile *fptr,*fmkf,*fhk;
    fits_open_file (&fptr , level1Sciencedatafile , READONLY , &status) ;
    //printError (status , "Error in the opening the input File" , level1Sciencedatafile) ;
    fits_movnam_hdu (fptr , BINARY_TBL , SCIENCEDATA_HDUNAME , 0 , &status) ;
   // printError (status , "Error in moving to the SCIENCE DATA HDU in the input File " , level1Sciencedatafile) ;
    fits_get_num_rows (fptr , &nrow_l1 , &status) ; //Get number of rows in level 1 science data file in 'nrow_l1'
   // printError (status , "Error in getting the number of rows of input file" , level1Sciencedatafile) ;
    int *time_l1= new int[nrow_l1];
    fits_read_col (fptr , TINT , 12 , 1 , 1 , nrow_l1 , NULL ,time_l1  , NULL , &status) ;
    
    cout<<"Total no of rows are "<<nrow_l1<<endl;
    
    for(int i=0;i<nrow_l1;i++)
    {
     time_final.push_back (time_l1[i]);    
    }   
    
    
    param_vect.clear ();
    param_vect.push_back ("SUN_ANGLE");
    param_vect.push_back ("MOON_ANGLE");
    
   //openig mkf file
    fits_open_file (&fmkf , mkffile , READONLY , &status) ;
    //printError (status , "Error in the opening the input File" , mkffile) ;
   
  // printError (status , "Error in moving to the SCIENCE DATA HDU in the input File " , mkffile) ;
     
    //opening hkfile
    cout<<status<<endl;
//    fits_open_file (&fhk , hkfile , READONLY , &status) ;
//  
//    //printError (status , "Error in the opening the input File" , mkffile) ;
//    fits_movnam_hdu (fhk , ASCII_TBL ,HK_HDU_NAME, 0 , &status) ;
//      cout<<status<<endl;
  //printError (status , "Error in moving to particular HDU" , mkffile) ;
     fits_get_num_hdus (fmkf , &numhdu_mkf , &status) ;
     cout<<numhdu_mkf<<endl;
     
     //fits_get_num_hdus (fhk , &numhdu_hk , &status) ;
     char *unit;
     bool flag =FALSE;
     char readunit_key[FLEN_FILENAME];
     int col_of_term;
     char param_unit[FLEN_FILENAME];
     string fkeyname;
     gti_paraminfo gtiparams_info;
     fits_movabs_hdu (fmkf , 2 , NULL , &status) ;
     fits_get_hdrspace (fmkf , &keyexist , NULL , &status) ;
     string term_str;
     vector<gti_paraminfo> gtiparams_vect;
  // for (int i=0; i<time_final.size ();i++)
    for (int i=0; i<5;i++)
   {
        gtiparams_vect.clear (); 
    
        for(int paramno=0;paramno<param_vect.size ();paramno++)
        {
                
           for(int numrecords=1;numrecords<=keyexist;numrecords++)
             {
                
                fits_read_record (fmkf , numrecords, record, &status) ;
               
                  if(strstr (record,param_vect[paramno].c_str())!=NULL)
                  {
                      const   char *temp=strstr(record,"=");
                      //cout<<temp<<endl;
                      int substr1=strlen(temp);
                      int main_str_len=strlen (record);
                      int diff_index=main_str_len-substr1-1-TTYPE_KEYWORDNAME_SKIP;
                      record_final=record;
                      fkeyname=record_final.substr (TTYPE_KEYWORDNAME_SKIP,diff_index);
                      sprintf(readunit_key,"TFORM%s",fkeyname.c_str());
                      
                      fits_read_keyword (fmkf , readunit_key ,param_unit , NULL , &status) ;
                            
                      col_of_term=atoi(fkeyname.c_str());
                      temp_num=0;
                      
                      findDataType(param_unit,term_datatype,temp_num);
             
                      if(temp_num==1) //incase of float datatype
                      {
                     
                      float term=0;                     
                      fits_read_col (fmkf , term_datatype, col_of_term, i+1 , felement , nelement , NULL , &term , NULL , &status) ;
                      cout<<term<<" "<<status<<" "<<col_of_term<<endl;
                      sprintf((char*)term_str.c_str (),"%f",term);
                      gtiparams_info.param_name=param_vect[paramno];
                      gtiparams_info.value_paramStr=term_str;
                      gtiparams_info.datatype_param=1;
                      gtiparams_vect.push_back (gtiparams_info);
                     // finalparam_allowedValues.push_back (param_allowedValues[paramno]);
                      
                      }
                      else   
                          int term;//to be changed
                      
                      break;
    
                  }
                       
             }
       
            
            
        }
          cout<<gtiparams_vect.size ()<<endl;
          checkValidRange(gtiparams_vect,param_allowedValues);
      
    }  
    
     //cout<<gtiparams_vect.size ()<<endl;
             
     
     
    fits_close_file(fptr,&status);
 //   printError (status , "Error in closing the level1 science data file" , mkffile) ;
    fits_close_file(fmkf,&status);
  // printError (status , "Error in closing the mkf file " , mkffile) ;
}

int  findDataType(char *tformterm ,int  &dataType,int &num)
{
    num=0;
if(strstr (tformterm,"E")!=NULL)
{

dataType=TFLOAT;
num=1;

}
   
return (EXIT_SUCCESS);
}

int checkValidRange (vector<gti_paraminfo>& gtiInfo, vector<string> allowedval){
    cout<<"Inside the Range checking"<<endl;
    
    for(int i=0;i<gtiInfo.size ();i++)
    {
       cout<<"parameter name"<<gtiInfo[i].param_name<<endl;  
       cout<<"Allowed Values  "<<allowedval[i]<<endl;
       if(strstr (allowedval[i].c_str (),"-")!=NULL)
       {
           cout<<"MIN-MAX found"<<endl;
           
       }
       else if(strstr(allowedval[i].c_str (),"<")!=NULL)
       {
           
       }
       else if (strstr(allowedval[i].c_str (),">")!=NULL){
           
           
       }
       else if (strstr(allowedval[i].c_str (),">") && strstr(allowedval[i].c_str (),">") !=NULL)
       {
                      
       }
    }
    
    
    return(EXIT_SUCCESS);
}