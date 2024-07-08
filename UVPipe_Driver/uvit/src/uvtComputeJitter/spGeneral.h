// File  Name           :       spGeneral.h
//
// ANALYST/CONCEPTUAL
// DESIGNER             :       Amit Gupta
//
// TECHNICAL DESIGNER   :       Amit Gupta
//
// IMPLEMENTER          :       Amit Gupta
//
// Version              :       1.0

#ifndef SP_GENERAL_H
#define SP_GENERAL_H

#define SGI_IRIX        1
//#define WINDOWS_NT    1

#ifndef FALSE
#define FALSE   0
#define TRUE    1
#endif


#ifdef WINDOWS_NT
#define DLLExport _declspec(dllexport)
#endif

#ifdef SGI_IRIX
#define DLLExport
typedef int  boolean;
#endif

#define ERROR_NO_BIAS           150
#define MODEL_ERROR_CODE_BIAS   200
#endif // SP_GENERAL_H
