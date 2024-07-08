/**************************************************************************************

	I N T E G R A L   S C I E N C E   D A T A   C E N T E R

	P A R A M E T E R   I N T E R F A C E   L I B R A R Y

Copyright:	(c) 1998 ISDC http://obswww.unige.ch/isdc/
File:		pil_error.c
Description:	Parameter Interface Library - PIL errors definitions
Authors:	Jerzy.Borkowski@obs.unige.ch (JB)
History:	13-May-98 (JB) : version 0.7 released. Code cleanup, unnecesary
			dependencies removed. C++ safe
		17-Aug-98 (JB) : version 1.0 released. Minor changes, error
			range adapted do -3000 - -3999

**************************************************************************************/


#include <stdio.h>
#include "pil.h"


const char	*PIL_err_table[PIL_ERR_MAX_IDX - PIL_ERR_MIN_IDX + 1] =
      {	"NULL pointer passed",
	"bad argument(s) passed",
	"no enough memory",
	"file does not exist",
	"error reading file",
	"error writing file",
	"end of string encountered",
	"bad name",
	"bad type",
	"bad mode",
	"bad line",
	"not implemented",
	"file does not exist",
	"file exists",
	"file not readable",
	"file not writeable",
	"blank line",
	"comment line",
	"error in line format",
	"not found",
	"PFILES env. variable too long",
#ifdef WIN32
	"PFILES env. missing '|'",
#else
	"PFILES env. missing ';'",
#endif
	"cannot lock/unlock parameter for exclusive access",
	"bogus parameters in command line",
	"no logger function",
	"too many tokens found in the line",
	"not enough tokens found in the line",
	"unmatched quote/doublequote found in the line",
	"line without terminating character",
	"extra space after trailing quote",
	"cannot convert string to boolean",
	"cannot convert string to integer",
	"cannot convert string to real",
	"cannot convert string to a var. vector of integers",
	"cannot convert string to a vector of integers",
	"cannot convert string to a var. vector of reals",
	"cannot convert string to a vector of reals",
	"value is out of range set by min/max fields",
	"value does not match any in enumerated list (min field)",
	/* SPR 3410: Error message for file access problem should be more sensible. */
	"file not found (or has wrong access type)"
      };


const char	*PIL_err_handler(int r)
 { if ((r >= PIL_ERR_MIN_IDX) && (r <= PIL_ERR_MAX_IDX))
     return(PIL_err_table[PIL_ERR_MAX_IDX - r]);
   return(NULL);
 }
