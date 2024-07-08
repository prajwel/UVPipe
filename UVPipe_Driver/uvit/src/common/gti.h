/* 
 * File:   gti.h
 * Author: uvit
 *
 * Created on November 19, 2014, 11:13 AM
 */

#ifndef GTI_H
#define	GTI_H



#include <vector>
#include <algorithm>
#include<string.h>

using namespace std;

int genL2gti(char *level1Sciencedatafile ,char *gtiparamfiles, char *mkffile ,char *hkfile,vector<unsigned char> &new_gtiFlag);

int   findDataType(char *tformterm ,int  &dataType,int &num);


//void readCPUTemp(); //NUV VIS FUV
//void readHUVTemp(); //NUV VIS FUV
//void readTTTempNUV(); //NUV VIS FUV
//void readTTTempFUV();
//void readMCPvolt();
//void readAnodeVolt();
//void readcathodVolt();
//void readSunAngle();
//void readMoonAngle();
//void brightLimbEarth();
//void readDriftRateAsp();
//void readMagnitudeJitter();
//void readSAAFluxCalculated();
//void readCPMcountRate();
//void readTimeSinceSAOexit();
//void readTimeToSAO();
//void readBITErrorRate();
//void readSSMmotion();
//void readLAXPCputrification();
//void readSolarPanalmotion();

#endif	/* GTI_H */

