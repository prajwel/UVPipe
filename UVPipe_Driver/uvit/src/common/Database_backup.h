/* 
 * File:   Database.h
 * Author: priyanka
 *
 * Created on July 30, 2015, 10:24 AM
 */

#ifndef DATABASE_H
#define	DATABASE_H
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
    
    int CirSelect(string RAdeg, string Dedeg, string radius,vector<string> &usno,vector<string> &radeg,vector<string> &dedeg,vector<string> &bmag,vector<string> &rmag,vector<string> &epoch);
    
    vector<vector<string> > query(char *query);
    
    int select(vector<string> &usno,vector<string> &radeg,vector<string> &dedeg,vector<string> &bmag,vector<string> &rmag,vector<string> &epoch,string ra_cent,string dec_cent);
private:

};

//INDEPENDENT FUNCTIONS
static int callback(void *NotUsed, int argc, char **argv, char **azColName);
#endif	/* DATABASE_H */

