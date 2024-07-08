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
    
    sqlite3_close(db);
    
    return 0;
}
//Function to create table
int Database::createTable(){
    
    string sql="";
    
    sql ="Create table ABC(USNO_A2 text primary key not null,"
            "RAdeg real not null,"
            "DEdeg real not null,"
            "BMag real not null,"
            "RMag real not null,"
            "Epoch real not null,"
            "ACT text,"
            "Wrong_Value text);";
    
    rc=sqlite3_exec(db,sql.c_str(),callback,0,&zErrMsg);
    
    if(rc!=SQLITE_OK){
    
        fprintf(stderr,"\n SQL error: %s\n",zErrMsg);
        
        sqlite3_free(zErrMsg);
    }
    
    else{
    
        fprintf(stdout,"\nTable created successfully\n");
    }
    return 0;
}
//Function to create index
int Database::createIndex(){
    
    string sql="";
    
    sql="create index uvit_index on ABC(RAdeg,DEdeg);";
    
    rc=sqlite3_exec(db,sql.c_str(),callback,0,&zErrMsg);
    
    if(rc!=SQLITE_OK){
        
        fprintf(stderr,"\n SQL Error:%s\n",zErrMsg);
    
        sqlite3_free(zErrMsg);
    }
    else{
        
        fprintf(stdout,"\n Index created successfully \n");
    }
    return 0;
}
//Function to insert data to database
int Database::Insert(string sql){
       
    sqlite3_exec(db,"BEGIN TRANSACTION",NULL,NULL,&zErrMsg);
    
    rc = sqlite3_exec(db, sql.c_str(), callback, 0, &zErrMsg);
    
    if( rc != SQLITE_OK ){
    
        fprintf(stderr, "SQL error Insert: %s\n", zErrMsg);
        
        sqlite3_free(zErrMsg);
        
    }
    
    else{
    
        fprintf(stdout, "Records created successfully\n");
        
    }  
    
    sqlite3_exec(db,"COMMIT TRANSACTION",NULL,NULL,&zErrMsg);
    
    closeDatabase();
}

//Function for rectangular search
int Database::RectSelect(string RAdeg, string DEdeg,string a,string b,vector<string> &usno,vector<string> &radeg,vector<string> &dedeg,vector<string> &bmag,vector<string> &rmag,vector<string> &epoch){

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
    
    for(vector<vector<string> >:: iterator it=result.begin();it<result.end();++it){
    
        vector<string> row=*it;
        
        cout<<endl<<row.at(0)<<"\t"<<row.at(1)<<"\t"<<row.at(2)<<"\t"<<row.at(3)<<"\t"<<row.at(4)<<"\t"<<row.at(5)<<endl;
      }
    
    
    return 0;
}


//Function for conic search
int Database::CirSelect(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &radeg,vector<string> &dedeg,vector<string> &bmag,vector<string> &rmag,vector<string> &epoch){

    float ra=0,dec=0,rad=0,r=0,s=0,t=0,u=0;
    
    char temp[200];
    
    string sql="";
    
    ra=atof(RAdeg.c_str());
    
    dec=atof(Dedeg.c_str());
    
    rad=atof(radius.c_str());
    
    r = ra+rad;
    
    s=ra-rad;
    
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
    
    sql="select * from Visible INDEXED BY uvit_index where RAdeg<"+x+" AND RAdeg>"+y+" AND DEdeg<"+z+" AND DEdeg>"+w+";";
    
    vector<vector<string> > result=query((char*)sql.c_str());
    
    vector<string> row;
    
   // vector<string> usno,radeg,dedeg,bmag,rmag,epoch;
    
    for(vector<vector<string> >:: iterator it=result.begin();it<result.end();++it){
    
        row=*it;
       
        float x=0,y=0,a1=0;
       
        x=atof(row.at(1).c_str());
       
        y=atof(row.at(2).c_str());
       
        a1=((ra-x)*(ra-x))+((dec-y)*(dec-y));
        
        if(a1<rad){
        
            usno.push_back(row.at(0).c_str());
            
            radeg.push_back(row.at(1).c_str());
            
            dedeg.push_back(row.at(2).c_str());
            
            bmag.push_back(row.at(3).c_str());
            
            rmag.push_back(row.at(4).c_str());
            
            epoch.push_back(row.at(5).c_str());
        }
        
    }
    
    
    for(int i=0;i<usno.size();i++){
    
        cout<<usno[i]<<"\t"<<radeg[i]<<"\t"<<dedeg[i]<<"\t"<<bmag[i]<<"\t"<<rmag[i]<<"\t"<<epoch[i]<<endl;
    }
    
    
    usno.clear();
    
    radeg.clear();
    
    dedeg.clear();
    
    bmag.clear();
    
    rmag.clear();
    
    
    epoch.clear();
      
    cout<<"Success"<<endl;
    
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


int Database::select (vector<string> &usno , vector<string> &radeg , vector<string> &dedeg , vector<string> &bmag , vector<string> &rmag , vector<string> &epoch , string ra_cent , string dec_cent)
{

    int i = 0 ;

    string a = "" , b = "" , r = "" , ra = "" , dec = "" ;

    cout << endl << "Select the value" << endl << "1.Rectangular Select" << endl << "2.Circular Select" << endl ;

    cout << "Which type of select is required :" ;

    cin >> i ;

    switch (i)
    {

        case 1:

            cout << endl << "RAdeg:" ;

            cin >> ra ;

            cout << endl << "DEdeg:" ;

            cin >> dec ;

            cout << endl << "a:" ;

            cin >> a ;

            cout << endl << "b:" ;

            cin >> b ;

            RectSelect (ra , dec , a , b , usno , radeg , dedeg , bmag , rmag , epoch) ;

            break ;
        case 2:


            ra = ra_cent , dec = dec_cent ;
         //   cout << endl << "Radius:" ;
            r=0.3;
          //  cin >> r ;
            //ra =;dec=;
            CirSelect (ra , dec , r , usno , radeg , dedeg , bmag , rmag , epoch) ;

            break ;
    }


    return 0 ;
}

static int callback(void *NotUsed, int argc, char **argv, char **azColName){

    cout<<"EXIT"<<endl;exit(1);
    
    int i;
    
    for(i=0;i<argc;i++){
    
        printf("%s = %s\n",azColName[i],argv[i]?argv[i]:"Null");
        
        printf("\n");
        
        return 0;
    }
    
}
