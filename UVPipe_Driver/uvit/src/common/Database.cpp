/* 
 * File:   Database.cpp
 * Author: priyanka
 * 
 * Created on July 30, 2015, 10:24 AM
 */

#include <cstdlib>
#include <iosfwd>
#include<fstream>
#include<vector>
#include<iostream>
#include<sqlite3.h>
#include<stdio.h>
#include <string.h>
#include "Database.h"
#include<math.h>

using namespace std;

Database::Database() {
}

Database::Database(const Database& orig) {
}

Database::~Database() {
}
//Function to open database
int Database::openDatabase(string databaseName){
           
    rc=sqlite3_open((char *) databaseName.c_str(),&db);
    
    if(rc){
    
        fprintf(stderr,"Cant open database: %s \n",sqlite3_errcode(db));
        
        exit(0);
    }
    
    
    else{
    
        fprintf(stderr,"Opened Database successfully");
    }
    
    
    return 0;
}
//Function to close database
int Database::closeDatabase(){
    
   rc=sqlite3_close(db);
        if(rc){
        
        fprintf(stderr,"Cant close database: %s \n",sqlite3_errcode(db));
        
        exit(0);
    }
    
    else{
     
        fprintf(stderr,"Closed Database successfully");
    }
    
    return 0;
}

//Function for rectangular search
int Database::RectSelect(string RAdeg, string DEdeg,string a,string b,vector<string> &usno,vector<string> &radeg,vector<string> &dedeg,vector<string> &bmag,vector<string> &rmag,vector<string> &epoch)
{

    string sql="";
    
    float ra=0,a1=0,dec=0,b1=0;
    
    string i="",j="",k="",l="";
    
    char temp[200];
    
    ra=atof(RAdeg.c_str());
    
    a1=atof(a.c_str());
    
    dec=atof(DEdeg.c_str());
    
    b1=atof(b.c_str());
    
    sprintf(temp,"%f",(ra+a1));
    
    i = (string) temp;
    
    sprintf(temp,"%f",(ra-a1));
    
    j = (string)temp;
    
    sprintf(temp,"%f",(dec+b1));
    
    k=(string)temp;
    
    sprintf(temp,"%f",(dec-b1));
    
    l=(string)temp;
    
    sql="select * from Visible INDEXED BY uvit_index where RAdeg<"+i+" AND RAdeg>"+j+" AND DEdeg<"+k+" AND DEdeg>"+l+";";
    
    vector<vector<string> > result=query((char*)sql.c_str());  
    
    cout.setf(ios::fixed,ios::floatfield);
    
    for(vector<vector<string> >:: iterator it=result.begin();it<result.end();++it)
    {
    
        vector<string> row=*it;
        usno.push_back(row.at(0).c_str());
            
            radeg.push_back(row.at(1).c_str());
            
            dedeg.push_back(row.at(2).c_str());
            
            bmag.push_back(row.at(3).c_str());
            
            rmag.push_back(row.at(4).c_str());
            
            epoch.push_back(row.at(5).c_str());
       // cout<<endl<<row.at(0)<<"\t"<<row.at(1)<<"\t"<<row.at(2)<<"\t"<<row.at(3)<<"\t"<<row.at(4)<<"\t"<<row.at(5)<<endl;
      }
    
    
    return 0;
}

//Function for rectangular search for FUV
int Database::RectSelectFUV(string RA, string DEC,string a,string b,vector<string> &usno,vector<string> &rafuv,vector<string> &decfuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts, vector<string> &errcounts){

    string sql="";
    
    float ra=0,a1=0,dec=0,b1=0;
    
    string i="",j="",k="",l="";
    
    char temp[200];
    
    ra=atof(RA.c_str());
    
    a1=atof(a.c_str());
    
    dec=atof(DEC.c_str());
    
    b1=atof(b.c_str());
    
    sprintf(temp,"%f",(ra+a1));
    
    i = (string) temp;
    
    sprintf(temp,"%f",(ra-a1));
    
    j = (string)temp;
    
    sprintf(temp,"%f",(dec+b1));
    
    k=(string)temp;
    
    sprintf(temp,"%f",(dec-b1));
    
    l=(string)temp;
    
    sql="select * from FUV INDEXED BY fuv_index where RA<"+i+" AND RA>"+j+" AND DEC<"+k+" AND DEC>"+l+";";
        
    vector<vector<string> > result=query((char*)sql.c_str());  
    
    cout.setf(ios::fixed,ios::floatfield);
    
    for(vector<vector<string> >:: iterator it=result.begin();it<result.end();++it){
    
        vector<string> row=*it;
            rafuv.push_back(row.at(0).c_str());
            
            decfuv.push_back(row.at(1).c_str());
            
            nuvMag.push_back(row.at(2).c_str());
            
            errnuv.push_back(row.at(3).c_str());
            
            fuvmag.push_back(row.at(4).c_str());
            
            errfuv.push_back(row.at(5).c_str());
            
            nuvcounts.push_back(row.at(6).c_str());
            
            errcounts.push_back(row.at(7).c_str());
       // cout<<endl<<row.at(0)<<"\t"<<row.at(1)<<"\t"<<row.at(2)<<"\t"<<row.at(3)<<"\t"<<row.at(4)<<"\t"<<row.at(5)<<endl;
      }
    return 0;
}

//Function for rectangular search for NUV
int Database::RectSelectNUV(string RA, string DEC,string a,string b,vector<string> &usno,vector<string> &rafuv,vector<string> &decfuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts, vector<string> &errcounts){

    string sql="";
    
    float ra=0,a1=0,dec=0,b1=0;
    
    string i="",j="",k="",l="";
    
    char temp[200];
    
    ra=atof(RA.c_str());
    
    a1=atof(a.c_str());
    
    dec=atof(DEC.c_str());
    
    b1=atof(b.c_str());
    
    sprintf(temp,"%f",(ra+a1));
    
    i = (string) temp;
    
    sprintf(temp,"%f",(ra-a1));
    
    j = (string)temp;
    
    sprintf(temp,"%f",(dec+b1));
    
    k=(string)temp;
    
    sprintf(temp,"%f",(dec-b1));
    
    l=(string)temp;
    
    sql="select * from NUV INDEXED BY nuv_index where RA<"+i+" AND RA>"+j+" AND DEC<"+k+" AND DEC>"+l+";";   
    
    vector<vector<string> > result=query((char*)sql.c_str());  
    
    cout.setf(ios::fixed,ios::floatfield);
    
    for(vector<vector<string> >:: iterator it=result.begin();it<result.end();++it){
    
        vector<string> row=*it;
            rafuv.push_back(row.at(0).c_str());
            
            decfuv.push_back(row.at(1).c_str());
            
            nuvMag.push_back(row.at(2).c_str());
            
            errnuv.push_back(row.at(3).c_str());
            
            fuvmag.push_back(row.at(4).c_str());
            
            errfuv.push_back(row.at(5).c_str());
            
            nuvcounts.push_back(row.at(6).c_str());
            
            errcounts.push_back(row.at(7).c_str());
            
        //cout<<endl<<row.at(0)<<"\t"<<row.at(1)<<"\t"<<row.at(2)<<"\t"<<row.at(3)<<"\t"<<row.at(4)<<"\t"<<row.at(5)<<endl;
      }
    
    
    return 0;
}
//Function for conic search
int Database::CirSelect(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &radeg,vector<string> &dedeg,vector<string> &bmag,vector<string> &rmag,
        vector<string> &epoch,double dec_pnt_value)
{
   // dec_pnt_value=1;
    float ra=0,dec=0,rad=0,r=0,s=0,t=0,u=0;
    
    char temp[200];
    
    string sql="";
    
    ra=atof(RAdeg.c_str());
    
    dec=atof(Dedeg.c_str());
    
    rad=atof(radius.c_str());
   
    //r = (ra+rad)*dec_pnt_value;
    
   // s=(ra-rad)*dec_pnt_value;
    //cout<<dec_pnt_value;exit(1);
    // r = ra+(rad);
    
    //s=ra-(rad);
    r = ra+(rad/(float)dec_pnt_value);
//    
    s=ra-(rad/(float)dec_pnt_value);
    // r = ra+(rad);
    
    //s=ra-(rad);
    t=dec+rad;
    
    u=dec-rad;
    
    string x="",y="",z="",w="";
    
    sprintf(temp,"%f",r);
    
    x=temp;
    
    sprintf(temp,"%f",s);
    
    y=temp;
    
    sprintf(temp,"%f",t);
    
    z=temp;
    
    sprintf(temp,"%f",u);
    
    w=temp;
    
    rad=rad*rad;
   // sql=".schema Visible";
    
//    sql="SELECT max() FROM Visible;"
    sql="select * from Visible INDEXED BY uvit_index where RAdeg<"+x+" AND RAdeg>"+y+" AND DEdeg<"+z+" AND DEdeg>"+w+";";
   
    //cout<<sql;
    vector<vector<string> > result=query((char*)sql.c_str());
    
    vector<string> row;
    
   // vector<string> usno,radeg,dedeg,bmag,rmag,epoch;
    //cout<<result.size()<<endl;
    for(vector<vector<string> >:: iterator it=result.begin();it<result.end();++it){
       
        row=*it;
       
        float x=0,y=0,a1=0;
       
        x=atof(row.at(1).c_str());
       
        y=atof(row.at(2).c_str());
       
        a1=((ra-x)*dec_pnt_value*(ra-x)*dec_pnt_value)+((dec-y)*(dec-y));
        //cout<<"FOR RA::: "<<ra<<" "<<dec<<" "<<row.at(1)<<" "<<row.at(2)<<" "<<row.at(3)<<" "<<row.at(4)<<" "<<a1<<" "<<rad<<endl;
        if(a1<rad)
        {
            //cout<<"FOR RA::: "<<ra<<" "<<dec<<row.at(1)<<" "<<row.at(2)<<" "<<row.at(3)<<" "<<row.at(4)<<endl;
            usno.push_back(row.at(0).c_str());
            
            radeg.push_back(row.at(1).c_str());
            
            dedeg.push_back(row.at(2).c_str());
            
            bmag.push_back(row.at(3).c_str());
            
            rmag.push_back(row.at(4).c_str());
            
            epoch.push_back(row.at(5).c_str());
        }
        
    }
    //cout<<"============================================================";
   // cout<<usno.size();
    cout.setf(ios::fixed,ios::floatfield);
    
//    for(int i=0;i<usno.size();i++){
//    
//        cout<<usno[i]<<"\t"<<radeg[i]<<"\t"<<dedeg[i]<<"\t"<<bmag[i]<<"\t"<<rmag[i]<<"\t"<<epoch[i]<<endl;
//    }
    
    
//    usno.clear();
//    
//    radeg.clear();
//    
//    dedeg.clear();
//    
//    bmag.clear();
//    
//    rmag.clear();
//    
//    epoch.clear();
    
   // cout<<"Success"<<endl;
    
    return 0;
}

//Function for conic search FUV
int Database::CirSelectFUV(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &rafuv,vector<string> &decfuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts, vector<string> &errcounts,double dec_pnt_value)
{

    float ra=0,dec=0,rad=0,r=0,s=0,t=0,u=0;
    
    char temp[200];
    
    string sql="";
    
    ra=atof(RAdeg.c_str());
    
    dec=atof(Dedeg.c_str());
    
    rad=atof(radius.c_str());
    
    r = ra+rad/dec_pnt_value;
    
    s=ra-rad/dec_pnt_value;
    
    t=dec+rad;
    
    u=dec-rad;
    
    string x="",y="",z="",w="";
    
    sprintf(temp,"%f",r);
    
    x=temp;
    
    sprintf(temp,"%f",s);
    
    y=temp;
    
    sprintf(temp,"%f",t);
    
    z=temp;
    
    sprintf(temp,"%f",u);
    
    w=temp;
    
    rad=rad*rad;
    
    sql="select * from FUV INDEXED BY fuv_index where RA<"+x+" AND RA>"+y+" AND DEC<"+z+" AND DEC>"+w+";";
    
    vector<vector<string> > result=query((char*)sql.c_str());
    
    vector<string> row;
    
  //  vector<string> RAFUV,DECFUV,NUVmag,errNUVmag,FUVmag,errFUVmag,NUVcounts,errcounts;
   cout<<"The result is "<<result.size();
    for(vector<vector<string> >:: iterator it=result.begin();it<result.end();++it){
    
        row=*it;
       
        float x=0,y=0,a1=0;
       
        x=atof(row.at(0).c_str());
       
        y=atof(row.at(1).c_str());
       
        a1=((ra-x)*dec_pnt_value*(ra-x)*dec_pnt_value)+((dec-y)*(dec-y));
        
        if(a1<rad){
        
//            usno.push_back(row.at(0).c_str());
//            
//            radeg.push_back(row.at(1).c_str());
//            
//            dedeg.push_back(row.at(2).c_str());
//            
//            bmag.push_back(row.at(3).c_str());
//            
//            rmag.push_back(row.at(4).c_str());
//            
//            epoch.push_back(row.at(5).c_str());
            rafuv.push_back(row.at(0).c_str());
            
            decfuv.push_back(row.at(1).c_str());
            
            nuvMag.push_back(row.at(2).c_str());
            
            errnuv.push_back(row.at(3).c_str());
            
            fuvmag.push_back(row.at(4).c_str());
            
            errfuv.push_back(row.at(5).c_str());
            
            nuvcounts.push_back(row.at(6).c_str());
            
            errcounts.push_back(row.at(7).c_str());
        }
        
    }
    
    cout.setf(ios::fixed,ios::floatfield);
    
//    for(int i=0;i<rafuv.size();i++)
//    {
//        cout<<rafuv[i]<<"\t"<<decfuv[i]<<"\t"<<nuvMag[i]<<"\t"<<errnuv[i]<<"\t"<<fuvmag[i]<<"\t"<<errfuv[i]<<"\t"<<nuvcounts[i]<<"\t"<<errcounts[i]<<endl;
//    }
    
//    rafuv.clear ();
//    decfuv.clear ();
//    nuvMag.clear ();
//    errnuv.clear ();
//    fuvmag.clear ();
//    errfuv.clear ();
//    fuvmag.clear ();
//    errfuv.clear ();
//    nuvcounts.clear ();
//    errcounts.clear();
//    RAFUV.clear();
//    
//    DECFUV.clear();
//    
//    NUVmag.clear();
//    
//    errNUVmag.clear();
//    
//    FUVmag.clear();
//    
//    errFUVmag.clear();
//    
//    NUVcounts.clear();
    
    
    
    //cout<<"Success"<<endl;
    
    return 0;
}

//Function for conic search FUV
int Database::CirSelectNUV(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &ranuv,vector<string> &decnuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts,vector<string> &errcounts,double dec_pnt_value)
{

    float ra=0,dec=0,rad=0,r=0,s=0,t=0,u=0;
    
    char temp[200];
    
    string sql="";
    
    ra=atof(RAdeg.c_str());
    
    dec=atof(Dedeg.c_str());
    
    rad=atof(radius.c_str());
    
    r = ra+rad/dec_pnt_value;
    
    s=ra-rad/dec_pnt_value;
    
    t=dec+rad;
    
    u=dec-rad;
    
    string x="",y="",z="",w="";
    
    sprintf(temp,"%f",r);
    
    x=temp;
    
    sprintf(temp,"%f",s);
    
    y=temp;
    
    sprintf(temp,"%f",t);
    
    z=temp;
    
    sprintf(temp,"%f",u);
    
    w=temp;
    
    rad=rad*rad;
    
    sql="select * from NUV INDEXED BY nuv_index where RA<"+x+" AND RA>"+y+" AND DEC<"+z+" AND DEC>"+w+";";
    
    vector<vector<string> > result=query((char*)sql.c_str());
    
    vector<string> row;
    //cout<<"The size is "<<result.size()<<endl;
    //vector<string> RANUV,DECNUV,NUVmag,errNUVmag,FUVmag,errFUVmag,NUVcounts,errcounts;
    
    for(vector<vector<string> >:: iterator it=result.begin();it<result.end();++it)
    {
    
        row=*it;
       
        float x=0,y=0,a1=0;
       
        x=atof(row.at(0).c_str());
       
        y=atof(row.at(1).c_str());
       
        a1=((ra-x)*dec_pnt_value*(ra-x)*dec_pnt_value)+((dec-y)*(dec-y));
        
        if(a1<rad)
        {
        
            ranuv.push_back(row.at(0).c_str());
            
            decnuv.push_back(row.at(1).c_str());
            
            nuvMag.push_back(row.at(2).c_str());
            
            errnuv.push_back(row.at(3).c_str());
            
            fuvmag.push_back(row.at(4).c_str());
            
            errfuv.push_back(row.at(5).c_str());
            
            nuvcounts.push_back(row.at(6).c_str());
            
            errcounts.push_back(row.at(7).c_str());
        }
        
    }
    
//    for(int i=0;i<ranuv.size();i++){
//    
//        cout<<ranuv[i]<<"\t"<<decnuv[i]<<"\t"<<nuvMag[i]<<"\t"<<errnuv[i]<<"\t"<<fuvmag[i]<<"\t"<<errfuv[i]<<"\t"<<nuvcounts[i]<<"\t"<<errcounts[i]<<endl;
//    }
    
    
//    ranuv.clear();
//    
//    decnuv.clear();
//    
//    nuvMag.clear();
//    
//    errnuv.clear();
//    
//    fuvmag.clear();
//    
//    errfuv.clear();
//    
//    nuvcounts.clear();
//    
//    errcounts.clear();
//    
    //cout<<"Success"<<endl;
    
    return 0;
}
//General function to execute query
vector<vector<string> > Database::query(char* query){

    vector<vector<string> > results;
    
    sqlite3_stmt *statement;
    
    if(sqlite3_prepare_v2(db,query,-1,&statement,0)==SQLITE_OK){
    
        int cols=sqlite3_column_count(statement);
        
        int result=0;
        
        while(true){
        
            result=sqlite3_step(statement);
            
            if(result==SQLITE_ROW){
            
                vector<string> values;
                
                for(int col=0;col<cols;col++){
                
                    values.push_back((char*)sqlite3_column_text(statement,col));
                }
                
                
                results.push_back(values);
            }
            
            
            else{
            
                break;
            }
            
        }
        
        
        sqlite3_finalize(statement);
    }
    
    
    string error=sqlite3_errmsg(db);
    
    if(error!="not an error")
        cout<<query<<" "<<error<<endl;
    
    return results;
}

int  Database::select (vector<string> &usno , vector<string> &radeg , vector<string> &dedeg , vector<string> &bmag , vector<string> &rmag , vector<string> &epoch , string ra_cent , string dec_cent  ,vector<string> &nuvMag,vector<string> &errNuv,vector<string> &fuvMag,vector<string> &errFuv,vector<string> &nuvCnts,vector<string> &errCnts
, int  search_algo,string lengh_a,string width_b,string rad_in ,double decangle ) 
{
    
    int i=0;
    
    string a="",b="",r="",ra="",dec="";
    
    //cout<<endl<<"Select the value"<<endl<<"1.Rectangular Select"<<endl<<"2.Circular Select"<<endl<<"3.Rectangular Select for FUV"<<endl
      //         <<"4.Circular Select for FUV"<<endl<<"5.Rectangular Select for NUV"<<endl<<"6.Circular Select for NUV"<<endl;
    
   // cout<<"Which type of select is required :";
    
    //cin>>i;
    
    switch(search_algo)
    {
    
        case 1:
            ra =ra_cent;dec=dec_cent;
//        
//            cout<<endl<<"RAdeg:";
//            
//            cin>>ra;
//            
//            cout<<endl<<"DEdeg:";
//            
//            cin>>dec;
//            
           //cout<<endl<<"a:";
            
           // cin>>a;
            
           // cout<<endl<<"b:";
            
          //  cin>>b;
            a=lengh_a;
            b=width_b;
            RectSelect (ra , dec , a , b , usno , radeg , dedeg , bmag , rmag , epoch) ;
            
            break;
        
        case 2:            
            ra =ra_cent;
            dec=dec_cent;
            r=rad_in;
            
//            cout<<endl<<"RAdeg:";
//            
//            cin>>ra;
//            
//            cout<<endl<<"DEdeg:";
//            
//            cin>>dec;
            
//            cout<<endl<<"Radius:";
//            
//            cin>>r;
            
            CirSelect(ra,dec,r,usno , radeg , dedeg , bmag , rmag , epoch,decangle);
            
            break; 
        
        case 3:
        
//            cout<<endl<<"RA:";
//            
//            cin>>ra;
            ra=ra_cent;
            dec=dec_cent;
//            cout<<endl<<"DEC:";
//            
//            cin>>dec;
//            
           // cout<<endl<<"a:";
            
          //  cin>>a;
            
         //   cout<<endl<<"b:";
            
         //   cin>>b;
            a=lengh_a;
            b=width_b;
            RectSelectFUV(ra,dec,a,b,usno,radeg,dedeg,nuvMag,errNuv,fuvMag,errFuv,nuvCnts,errCnts);
            
            break;
        case 4:
              ra =ra_cent;
            dec=dec_cent;
            r=rad_in;
//            cout<<endl<<"RA:";
//            
//            cin>>ra;
//            
//            cout<<endl<<"DEC:";
//            
//            cin>>dec;
//            
//            cout<<endl<<"Radius:";
//            
//            cin>>r;
            
            CirSelectFUV(ra,dec,r,usno,radeg,dedeg,nuvMag,errNuv,fuvMag,errFuv,nuvCnts,errCnts,decangle);
            
            break; 
        case 5:
            ra=ra_cent;
            dec=dec_cent;
//            cout<<endl<<"RA:";
//            
//            cin>>ra;
//            
//            cout<<endl<<"DEC:";
//            
//            cin>>dec;
            
           // cout<<endl<<"a:";
            
          //  cin>>a;
            
          //  cout<<endl<<"b:";
            
          //  cin>>b;
            a=lengh_a;
            b=width_b;
            RectSelectNUV(ra,dec,a,b,usno,radeg,dedeg,nuvMag,errNuv,fuvMag,errFuv,nuvCnts,errCnts);
            
            break;
        case 6:
              ra =ra_cent;
            dec=dec_cent;
            r=rad_in;
//            cout<<endl<<"RA:";
//            
//            cin>>ra;
//            
//            cout<<endl<<"DEC:";
//            
//            cin>>dec;
//            
//            cout<<endl<<"Radius:";
//            
//            cin>>r;
            
            CirSelectNUV(ra,dec,r,usno,radeg,dedeg,nuvMag,errNuv,fuvMag,errFuv,nuvCnts,errCnts,decangle);
            
            break; 
        default:
            
            if(i>6)
                cout<<"No such choice exist.."<<endl;
            
            break;
    }
    
    
    return 0;
}
static int callback(void *NotUsed, int argc, char **argv, char **azColName)
{
    
    //cout<<"EXIT"<<endl;exit(1);
    
    int i;
    
    for(i=0;i<argc;i++){
    
        printf("%s = %s\n",azColName[i],argv[i]?argv[i]:"Null");
        
        printf("\n");
        
        return 0;
    }
    
    
}
