/* 
 * File:   Database.h
 * Author: priyanka
 *
 * Created on July 30, 2015, 10:24 AM
 */

#ifndef DATABASE_H
#define	DATABASE_H
#include<vector>
#include<sqlite3.h>
using namespace std;
class Database {
public:
    Database();
    Database(const Database& orig);
    virtual ~Database();
    
    char *zErrMsg;
    
    int  rc;
    
    sqlite3 *db;
    
    int openDatabase(string databaseName);
    
    int createTable();
    
    int createIndex();
    
    int closeDatabase();
    
    int Insert(string sql);
    
  int RectSelect(string RAdeg, string DEdeg,string a,string b,vector<string> &usno,vector<string> &radeg,vector<string> &dedeg,vector<string> &bmag,vector<string> &rmag,vector<string> &epoch);
  
    //int RectSelect(string RAdeg,string DEdeg,string a,string b);
    
    int RectSelectFUV(string RA, string DEC,string a,string b,vector<string> &usno,vector<string> &rafuv,vector<string> &decfuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts, vector<string> &errcounts);
    
    int RectSelectNUV(string RA, string DEC,string a,string b,vector<string> &usno,vector<string> &rafuv,vector<string> &decfuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts, vector<string> &errcounts);
    
  //  int CirSelect(string RAdeg,string Dedeg,string radius);
    
   // int CirSelect(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &radeg,vector<string> &dedeg,vector<string> &bmag,vector<string> &rmag,vector<string> &epoch);
     int CirSelect(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &radeg,vector<string> &dedeg,vector<string> &bmag,vector<string> &rmag,
        vector<string> &epoch,double dec_pnt_value); 
    //int CirSelectFUV(string RA,string DEC,string radius);
int  CirSelectFUV(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &rafuv,vector<string> &decfuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts, vector<string> &errcounts,double dec_pnt_value);
     //int CirSelectFUV(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &rafuv,vector<string> &decfuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts,vector<string> &errcounts);
  //  int CirSelectNUV(string RA,string DEC,string radius);
int CirSelectNUV(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &ranuv,vector<string> &decnuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts,vector<string> &errcounts,double dec_pnt_value);   
    //  int CirSelectNUV(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &rafuv,vector<string> &decfuv,vector<string> &nuvMag,vector<string> &errnuv,vector<string> &fuvmag,vector<string> &errfuv,vector<string> &nuvcounts,vector<string> &errcounts);
    vector<vector<string> > query(char *query);
    
    int count();
  // int  select (vector<string> &usno , vector<string> &radeg , vector<string> &dedeg , vector<string> &bmag , vector<string> &rmag , vector<string> &epoch , string ra_cent , string dec_cent  ,vector<string> nuvMag,vector<string> errNuv,vector<string> fuvMag,vector<string> errFuv,vector<string> nuvCnts,vector<string> errCnts
//, int  search_algo,string lengh_a,string width_b);
  
//int   select (vector<string> &usno , vector<string> &radeg , vector<string> &dedeg , vector<string> &bmag , vector<string> &rmag , vector<string> &epoch , string ra_cent , string dec_cent  ,vector<string> nuvMag,vector<string> errNuv,vector<string> fuvMag,vector<string> errFuv,vector<string> nuvCnts,vector<string> errCnts
//, int  search_algo,string lengh_a,string width_b,string rad_in);
int select (vector<string> &usno , vector<string> &radeg , vector<string> &dedeg , vector<string> &bmag , vector<string> &rmag , vector<string> &epoch , string ra_cent , string dec_cent  ,vector<string> &nuvMag,vector<string> &errNuv,vector<string> &fuvMag,vector<string> &errFuv,vector<string> &nuvCnts,vector<string> &errCnts
, int  search_algo,string lengh_a,string width_b,string rad_in ,double decangle );   

 
    // int select (vector<string> &usno , vector<string> &radeg , vector<string> &dedeg , vector<string> &bmag , vector<string> &rmag , vector<string> &epoch , string ra_cent , string dec_cent  ,vector<string> nuvMag,vector<string> errNuv,vector<string> fuvMag,vector<string> errFuv,vector<string> nuvCnts,vector<string> errCnts);
private:

};

//INDEPENDENT FUNCTIONS
static int callback(void *NotUsed, int argc, char **argv, char **azColName);
#endif	/* DATABASE_H */

