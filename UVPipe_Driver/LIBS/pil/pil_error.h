/**************************************************************************************

	I N T E G R A L   S C I E N C E   D A T A   C E N T E R

	P A R A M E T E R   I N T E R F A C E   L I B R A R Y

Copyright:	(c) 1998 ISDC http://obswww.unige.ch/isdc/
File:		pil_error.h
Description:	Parameter Interface Library - PIL errors header file
Authors:	Jerzy.Borkowski@obs.unige.ch (JB)
History:	13-May-98 (JB) : version 0.7 released. Code cleanup, unnecesary
			dependencies removed. C++ safe
		17-Aug-98 (JB) : version 1.0 released. Minor changes, error
			range adapted do -3000 - -3999

**************************************************************************************/

#ifndef	PIL_ERROR_H_INCLUDED
#define	PIL_ERROR_H_INCLUDED

#ifndef EXPSYM
#ifdef WIN32
#define EXPSYM __declspec(dllexport)
#else
#define EXPSYM
#endif
#endif

#if defined(__cplusplus) && ! defined(__CINT__)
extern "C" {
#endif

#define	PIL_OK				(0)

#define	PIL_ERR_BASE			(-3000)
#define	PIL_ERR_MAX_IDX			(PIL_ERR_BASE)

#define	PIL_NUL_PTR			(PIL_ERR_BASE - 0)
#define	PIL_BAD_ARG			(PIL_ERR_BASE - 1)
#define	PIL_NO_MEM			(PIL_ERR_BASE - 2)
#define	PIL_NO_FILE			(PIL_ERR_BASE - 3)
#define	PIL_ERR_FREAD			(PIL_ERR_BASE - 4)
#define	PIL_ERR_FWRITE			(PIL_ERR_BASE - 5)
#define	PIL_EOS				(PIL_ERR_BASE - 6)
#define	PIL_BAD_NAME			(PIL_ERR_BASE - 7)
#define	PIL_BAD_TYPE			(PIL_ERR_BASE - 8)
#define	PIL_BAD_MODE			(PIL_ERR_BASE - 9)
#define	PIL_BAD_LINE			(PIL_ERR_BASE - 10)
#define	PIL_NOT_IMPLEMENTED		(PIL_ERR_BASE - 11)
#define	PIL_FILE_NOT_EXIST		(PIL_ERR_BASE - 12)
#define	PIL_FILE_EXIST			(PIL_ERR_BASE - 13)
#define	PIL_FILE_NO_RD			(PIL_ERR_BASE - 14)
#define	PIL_FILE_NO_WR			(PIL_ERR_BASE - 15)
#define	PIL_LINE_BLANK			(PIL_ERR_BASE - 16)
#define	PIL_LINE_COMMENT		(PIL_ERR_BASE - 17)
#define	PIL_LINE_ERROR			(PIL_ERR_BASE - 18)
#define	PIL_NOT_FOUND			(PIL_ERR_BASE - 19)
#define	PIL_PFILES_TOO_LONG		(PIL_ERR_BASE - 20)
#define	PIL_PFILES_FORMAT		(PIL_ERR_BASE - 21)
#define	PIL_LOCK_FAILED			(PIL_ERR_BASE - 22)
#define	PIL_BOGUS_CMDLINE		(PIL_ERR_BASE - 23)
#define	PIL_NO_LOGGER			(PIL_ERR_BASE - 24)
#define	PIL_LINE_TOO_MANY		(PIL_ERR_BASE - 25)
#define	PIL_LINE_TOO_FEW		(PIL_ERR_BASE - 26)
#define	PIL_LINE_UNMATCHED_QUOTE	(PIL_ERR_BASE - 27)
#define	PIL_LINE_NO_LF			(PIL_ERR_BASE - 28)
#define	PIL_LINE_EXTRA_SPACES		(PIL_ERR_BASE - 29)
#define	PIL_BAD_VAL_BOOL		(PIL_ERR_BASE - 30)
#define	PIL_BAD_VAL_INT			(PIL_ERR_BASE - 31)
#define	PIL_BAD_VAL_REAL		(PIL_ERR_BASE - 32)
#define	PIL_BAD_VAL_INT_VAR_VECTOR	(PIL_ERR_BASE - 33)
#define	PIL_BAD_VAL_INT_VECTOR		(PIL_ERR_BASE - 34)
#define	PIL_BAD_VAL_REAL_VAR_VECTOR	(PIL_ERR_BASE - 35)
#define	PIL_BAD_VAL_REAL_VECTOR		(PIL_ERR_BASE - 36)
#define	PIL_OFF_RANGE			(PIL_ERR_BASE - 37)
#define	PIL_BAD_ENUM_VALUE		(PIL_ERR_BASE - 38)
#define	PIL_BAD_FILE_ACCESS		(PIL_ERR_BASE - 39)

#define	PIL_ERR_MIN_IDX			(PIL_ERR_MAX_IDX - 39)

			/* variables */
    
extern	const char	*PIL_err_table[PIL_ERR_MAX_IDX - PIL_ERR_MIN_IDX + 1];

			/* function prototypes */
#ifndef __CINT__
EXPSYM
#endif
const char	*PIL_err_handler(int r);

#if defined(__cplusplus) && ! defined(__CINT__)
}
#endif

#endif
