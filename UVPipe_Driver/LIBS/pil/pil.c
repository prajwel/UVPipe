/**************************************************************************************

	I N T E G R A L   S C I E N C E   D A T A   C E N T E R

	P A R A M E T E R   I N T E R F A C E   L I B R A R Y

Copyright:	(c) 1998-2002 ISDC http://isdc.unige.ch
File:		pil.c
Description:	Parameter Interface Library - core & API routines
Authors:	Jerzy.Borkowski@obs.unige.ch (JB)
History:	06-Oct-97 (JB) : support for command line arguments added
		10-May-98 (JB) : PILGetReal4 function added 
		13-May-98 (JB) : version 0.7 released. Code cleanup, unnecesary
			dependencies removed. Uses/links only PIL library. Changes
			to Makefile.
		13-Aug-98 (JB) : version 1.0 released. More compatible with Common
			library some new functions: PILReloadParameters,
			PILFlushParameters, locking of parameter file during
			read/write implemented.
		17-Aug-98 (JB) : reformatted, more comments added
		25-Aug-98 (JB) : undefined flock_t bug fixed (linux, osf)
		15-Mar-99 (JB) : small changes in header file reflecting changes
			in F90 bindings (for better portability)
		25-Jun-99 (JB) : added "" -> empty string substitution.
			added vector parameters support (3 new API functions)
		08-Dec-99 (JB) : fixed problem with modification time comparison
			between user and system parameter files (if both are present)
		27-Jan-00 (JB) : code const-ized in order to support ANSI C++
		04-Feb-00 (JB) : PILGetNumParameters/PILGetParameter added to API.
			changed logic in build_names(). Now if PFILES specifies
			both system and user part, and there is no system file
			and there is user file function succeeds.
		16-Feb-00 (JB) : bugfix in build_names().
		23-Feb-00 (JB) : bugfix in build_names(). (write to read-only data)
		17-Jul-00 (JB) : bugfix in PIL_verify_access(). _CREATE was creating
			0 length parameter file.
		22-Aug-00 (JB) : added PILVerifyCmdLine().
		27-Sep-00 (JB) : NUL PTR bug fixes in xxxget/xxxput, xxxVarVector added
		31-Oct-00 (JB) : added expansion of ${...}xxcx${...}
		09-Feb-01 (JB) : version 1.6.3 : enumerated values implemented. min 
			field must contain '|' separated list of allowable values,
			max field must be empty. String compare is case insensitive,
			returned value is uppercase.
		20-Apr-01 (JB) : version 1.7.0 : readline support added
		03-Jul-01 (JB) : version 1.7.1 : access renamed to pil_access (was 
			clashing with C++ (SPR 479), PILGetBool - case insensitive 
			now (SCREW 289), PILGetFname now trims spaces (SCREW 255)
		14-Aug-01 (JB) : version 1.7.2 : even more access variables renamed 
			to pil_access (SPR 0479, 0559)
		08-Mar-02 : version 1.8.0, fixes in trim_spaces, added root_name
			support, fixes in flush_params, reload_params and handling
			of the default mode (bugs were reported by James Peachey
			(peachey@bigband.gsfc.nasa.gov) and Bill Pence
			(pence@tetra.gsfc.nasa.gov)
		03-Apr-02 (JB): version 1.8.1-beta, SPR 1153 implemented, Also 2 extra 
			modes for better compatibility with other software :

			1) special READLINE mode : edit buffer is initially set to
			the default value. Pressing ENTER accepts default value (as
			before). However, pressing Ctrl-W (or DEL/BACKSPACE several times, 
			or ctrl+a then ctrl+K) allows empty string to be returned without
			setenv PIL_EMPTY_STRING or similar tricks. This mode is 
			set by calling PILSetReadlinePromptMode(PIL_RL_PROMPT_IEB)
			or by compiling library with STARTUP_PIL_RL_PROMPT_IEB
			defined. Note that default value is _NOT_ printed inside
			square brackets before ':' character anymore.

			2) special cmndline args mode : When this mode is in effect
			the command line arguments must be given using the following
			syntax :

				/path/program_name  name1 = value1 name2 = value2 [...]

			This mode is set by calling PILSetCmndArgMode(PIL_CMNDARG_NSESA)
			or by compiling library with STARTUP_PIL_CMNDARG_NSESA defined.
		19-Apr-02 (JB): version 1.8.1 : SPR 1153, SCREW 00507
		22-Apr-02 (JB): version 1.8.2-beta : renamed some symbols, PIL_CMNDARG_NSESA
			code removed, new parameter syntax implemented :
				when testing if argv[i] contains '=', the following symbols
				are ignored: ==, !=, <=, >=, =<, =>
			Also "named" parameter's syntax is relaxed. "named" parameters can
			span multiple cmdline args. Recognized formats :
				name=value
				name= vale
				name =value
				name = value
		06-May-02 (MT/BP/JP): version 1.8.2-rc1, incorporated 3 fixes by MT/BP/JP.
		08-May-02 (JB): implemented SPR 00687, SCREW 00370, 00505, this means that
			pget/pset/plist/pil_lint/pil_gen_c_code moved to new package :
			pil_tools, furthermode pget -> pil_get, pset -> pil_set, plist
			-> pil_list (but also temporarily there's a link with old name)
			PIL_LINT_MODE implemented, so that pil_lint is much more paranoid
			about parameter file syntax.
		20-May-02 (JB) : reverted some changes made in version 1.8.2 (SPR 1367)
		30-May-02 (JB) : restored N=V syntax as default (SPR 1442)
		23-Aug-02 (JB) : SPR 1618 : uninit read 

**************************************************************************************/


#include <ctype.h>
#include <stdio.h>

/* SCREW 1145: malloc.h not needed because stdlib.h is included.
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#else
#ifndef __APPLE__
#include <malloc.h>
#endif
#endif
*/

#include <string.h>
/* SCREW 1145: system-dependent functions are now in pil_sys.c
#include <unistd.h>
*/
#include <stdlib.h>
/* SCREW 1145: system-dependent functions are now in pil_sys.c
#include <sys/types.h>
#include <unistd.h>

#include <fcntl.h>
*/

/* SCREW 1145: time_t is used in this file; it must have been
   included by accident by some other header file. */
#include <time.h>

/* SPR 2592: Need SCHAR_MAX for handling escape sequences. */
#include <limits.h>

#ifndef HAVE_FLOCK_T				/* some systems (linux, osf) have flock and not flock_t */
typedef struct flock flock_t;			/* anyway, this should be somehow better handled, by */
#define	HAVE_FLOCK_T				/* autoconf stuff. Maybe global isdcautoconf.h file ?? */
#endif						/* for me, better place for ifndef ... would be pil.h, but ... */


#include "pil.h"				/* include our stuff */


PIL_TYPE_TABLE PIL_ttab[] = 			/* conversion table between constants */
      {						/* and names of parameter types */
	{ PIL_TYPE_BOOL,	"boolean" },
	{ PIL_TYPE_INT4,	"integer" },
	{ PIL_TYPE_REAL8,	"real8  " },
	{ PIL_TYPE_STRING,	"string " },
	{ PIL_TYPE_FNAME,	"filenam" },
	{ PIL_TYPE_STRUCT,	"structr" },
	{ PIL_TYPE_IMCUR,	"imgcrsr" },
	{ PIL_TYPE_GCUR,	"gphcrsr" },
	{ PIL_TYPE_D,		"D" },
	{ PIL_TYPE_G,		"G" },
	{ PIL_TYPE_UNKNOWN,	" -\?\?\?- " }	/* this has to be the last item in table */
      };

/* SPR 3169: Make sure PILDefaultPFile is constructed with reasonable values. */
EXPSYM PIL_PFILE	PILDefaultPFile = { 0, NULL, NULL, 0, NULL, 0, 0, NULL }; /* handle of default PFILE for applic */
char		PILSysPFile[PIL_LINESIZE];	/* name of system parameter file */
char		PILUsrPFile[PIL_LINESIZE];	/* name of user parameter file */
int		PILRunMode = ISDC_SINGLE_MODE;	/* mode of operation */
int		PILF90argc = 0;			/* argc from F90 main() */
char		**PILF90argv = NULL;		/* argv from F90 main() */
char		PILModuleName[PIL_LINESIZE];	/* name of module */
char		PILModuleVersion[PIL_LINESIZE];	/* version of module */
int		PILServerState = ISDC_TERMINATE;/* state of server mode */
int		PILSpecialMode = 0;		/* extras, like PSET mode */
int		PILQueryMode = PIL_QUERY_DEFAULT; /* do we want to disable prompting user ? */
char		*PILEmptyStringTemplate = NULL; /* value for empty string substitution */
int		PILVectorType = PIL_TYPE_STRING; /* used to verify vector values */
int		PILVectorLen = 0;		/* used to verify vector values */
int		PILVectorVarMode = 0;		/* nonzero - variable length vector */
char		PILLastNotFoundName[PIL_LINESIZE]; /* name of last not found parameter */
int		(*PILRootFnameFunction)(char *s) = PIL_cfitsio_root_name;
#ifdef	STARTUP_PIL_RL_PROMPT_IEB
int		PILReadlinePromptMode = PIL_RL_PROMPT_IEB;	/* Prompting mode when using READLINE */
#else
int		PILReadlinePromptMode = PIL_RL_PROMPT_PIL;	/* Prompting mode when using READLINE */
#endif
#ifdef	STARTUP_PIL_CMNDARG_ALLOW_SPACES
int		PILCmndArgMode = PIL_CMNDARG_ALLOW_SPACES;  /* Name [ space] = [ space ] Value */
#else
int		PILCmndArgMode = PIL_CMNDARG_PIL;  /* 'Name=Value' or positional */
#endif
int		(*PILLoggerFunction)(char *msg) = NULL;

#ifdef WIN32
static const char sPathDelim = '\\';
static const char sPathDelimString[] = "\\";
static const char sPathFieldSep[] = ";";
static const char sPfilesDelim = '|';
#else
static const char sPathDelim = '/';
static const char sPathDelimString[] = "/";
static const char sPathFieldSep[] = ":";
static const char sPfilesDelim = ';';
#endif

		/* support and conversion routines */


#ifdef WIN32
/* SCREW 1145: Provide strcasecmp implementation for Windows. */
/**************************************************************************************
Function:	strcasecmp
Description:	Windows replacement for UNIX standard strcasecmp utility.
		Provides identical behavior to the standard utility
		as implemented on Linux and Solaris. It does not check
		the pointers passed!
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
s1 (In)		first string to compare
s2 (Out)	second string to compare

Return code:	negative if s1 < s2, 0 if s1 == s2, positive if s1 > s2
**************************************************************************************/
static int strcasecmp(const char *s1, const char *s2)
 { const char *p1 = s1;
   const char *p2 = s2;
   int diff = 0;
   while (!diff) {
     diff = toupper(*p1) - toupper(*p2);
     if (*p1 && *p2) { ++p1; ++p2; }
     else { diff = *p1 - *p2; break; }
   }
   return diff;
 }
#endif

/**************************************************************************************

Function:	PIL_log_error
Description:	prints error via PILLoggerFunction
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
s (In)		String to print before error code message
errcode (In)	error code

Return code:	Function returns PIL_OK, or error code : PIL_NO_LOGGER if there is
		no function defined

**************************************************************************************/

int	PIL_log_error(char *s, int errcode)
 { char	buf[PIL_LINESIZE * 3 + 1000];

   if (NULL == PILLoggerFunction)  return(PIL_NO_LOGGER);
   if (PIL_OK == errcode)  return(PIL_OK);

   if (NULL != s)  { sprintf(buf, "%s : %s", s, PIL_err_handler(errcode)); }
   else  { sprintf(buf, "%s", PIL_err_handler(errcode)); }

   return((*PILLoggerFunction)(buf));
 }


/**************************************************************************************

Function:	PIL_root_name
Description:	finds name of file from input URL. Input URL is replaced by
		its file name component (if evaluates to file)
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
s (In/Out)		name of URL to analyse

Return code:	PIL_ROOTNAME_FILE, PIL_ROOTNAME_STDIN, PIL_ROOTNAME_STDOUT or
		PIL_ROOTNAME_STDINOUT if URL evaluates to file, stdin, stdout or
		stdin/stdout ('-' case). PIL_ROOTNAME_NOTFILE otherwise.
		If verification function is NULL then  it returns PIL_ROOTNAME_FILE.

**************************************************************************************/

int	PIL_root_name(char *s)
 {
   if (NULL == s)  return(PIL_NUL_PTR);
   if (NULL == PILRootFnameFunction)  return(PIL_ROOTNAME_FILE);
   return((*PILRootFnameFunction)(s));
 }


/**************************************************************************************

Function:	PIL_cfitsio_root_name
Description:	Parse the input URL or filename, which may use CFITSIO's extended
		file name syntax, to determine if it refers to a file on magnetic
		disk (as opposed to a remote ftp:// or http:// file specifier).  If
		so, then this routines returns a string containing the root disk
		file name, minus any filtering qualifiers that may have been
		appended. Input URL is replaced by its file name component (if
		evaluates to file).
		Based on original code from: pence@tetra.gsfc.nasa.gov
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
s (In/Out)		name of URL to analyse

Return code:	PIL_ROOTNAME_FILE, PIL_ROOTNAME_STDIN, PIL_ROOTNAME_STDOUT or
		PIL_ROOTNAME_STDINOUT if URL evaluates to file, stdin, stdout or
		stdin/stdout ('-' case). PIL_ROOTNAME_NOTFILE otherwise.

**************************************************************************************/

int	PIL_cfitsio_root_name(char *s)
 { size_t	slen;
   char		*cptr;
   int		ii;


   if (NULL == s)  return(PIL_NUL_PTR);
   for (cptr = s; ' ' == *cptr; cptr++);	/* ignore leading blanks in the name */

   if (*cptr == '\0')  return(PIL_ROOTNAME_NOTFILE); /* blank file name is illegal */

						/* '-' means either stdin or stdout */
   if (!strncmp(cptr,  "-",     1))  return(PIL_ROOTNAME_STDINOUT);

   if (!strncmp(cptr, "stdin",  5) ||		/* stdin case only */
       !strncmp(cptr, "STDIN",  5))  return(PIL_ROOTNAME_STDIN);

   if (!strncmp(cptr, "ftp:",   4) ||		/* various URLs */
       !strncmp(cptr, "http:",  5) ||
       !strncmp(cptr, "root:",  5) ||
       !strncmp(cptr, "mem:",   4) ||
       !strncmp(cptr, "shmem:", 6)) return(PIL_ROOTNAME_FILE);

	/*  skip over optional 'disk:' or 'disk://' or 'file:' or 'file://' prefix */

   if ((!strncmp(cptr, "disk:", 5)) || (!strncmp(cptr, "file:", 5)))
     { cptr += 5;
       if (!strncmp(cptr, "//", 2))  cptr += 2;
     }

	/*  find length of root name, up to the first '(', '[', or space char */
   slen = strcspn(cptr, "([ ");
   if (0 == slen)  return(PIL_ROOTNAME_NOTFILE);  /* degenerate case */

   memmove(s, cptr, slen);			/* copy the root name to the output string */
   s[slen] = 0;

	/* Check if the last characters in the file name are '+n' where n is a
	   1 to 5 digit integer.  If so, then strip off this qualifier from
	   the root name. (This is used to specify an HDU within the FITS file). */

   cptr = s + slen - 1;
   for (ii = 0; ii < 5 && ii < slen - 2; ii++)
    {
      if (((*cptr) < '0') || ((*cptr) > '9'))
        break;					/* last characters are not digits */

      cptr--;
      if ('+' == *cptr)				/* look for a plus sign   */
        { *cptr = 0;				/* terminate the string here */
	  break;
        }
    }

   return(PIL_ROOTNAME_FILE);
 }


/**************************************************************************************

Function:	PIL_file_copy
Description:	copies src to dest.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
dest (In)		full path to destination file
src (In)		full path to source file

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	Binary copy is performed. No locking during copy. Function closes all files
	and deallocates all temporary buffers before returning.

**************************************************************************************/

int	PIL_file_copy(const char *dest, const char *src)
 { FILE *inf = NULL, *outf = NULL;
   char buf[PIL_COPYCHUNK];
   size_t fs, i, delta;

/* SCREW 1145: Use of stat is encapsulated in PIL_get_file_size now.
   struct stat fst;
*/
   int r;

   r = PIL_OK;
   if (NULL == src) return(PIL_NUL_PTR);	/* check for valid args */
   if (NULL == dest) return(PIL_NUL_PTR);

/* SCREW 1145: Get file size before file is officially opened. */
	/* get file size */
   r = PIL_get_file_size(src, &fs);
   if (PIL_OK != r) goto EEX1;

	/* open src & dst files */

   if (NULL == (inf = fopen(src, "r")))  { r = PIL_NO_FILE; goto EEX1; }
   if (NULL == (outf = fopen(dest, "w")))  { r = PIL_NO_FILE; goto EEX1; }

/* SCREW 1145: Use of stat is encapsulated in PIL_get_file_size now.
   if (fstat(fileno(inf), &fst))  { r = PIL_BAD_ARG; goto EEX4; }
   fs = (size_t)fst.st_size;
*/
	/* copy data */

   for (i=0; i<fs; i+= PIL_COPYCHUNK)
    { if ((i + PIL_COPYCHUNK) > fs) delta = fs  - i;
      else delta = PIL_COPYCHUNK;
      if (delta != fread(buf, 1, delta, inf))  { r = PIL_ERR_FREAD; break; }
      if (delta != fwrite(buf, 1, delta, outf))  { r = PIL_ERR_FWRITE; break; }
    }

EEX1:   if (NULL != outf) fclose(outf);		/* cleanup routines */
        if (NULL != inf) fclose(inf);
   return(r);
 }


/**************************************************************************************

Function:	PIL_cmndarg2value
Description:	return value of given parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
numpar (In)		position of parameter in question (used when no parameter with
			'name' can be found)
*result (Out)		returned value of parameter (this is always string)
argc			number of entries in argument list (usually from main)
argv			argument list (usually from main)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	PIL_cmndarg2value scans argument list (given by argc and argv parameters) and
	returns in *result value of parameter with name 'name' (1st parameter) or if
	there is no parameter with such name value of numpar-th parameter. If numpar
	is out of range PIL_NOT_FOUND error is returned. Parameters in parameter
	list should be in the following format :

		named      : Name=Value
		positional : Value

	It is possible to mix named and positional arg. If 2 or more named args have
	the same name, value from first occurence is returned.

Warning:  *result returns pointer to an argv[] element (possibly with some offset).

**************************************************************************************/

int	PIL_cmndarg2value(const char *name, int numpar, char **result, int argc, char **argv)
 { char	buf[PIL_LINESIZE], *p;
   int	i;

   if (NULL == name) return(PIL_NUL_PTR);
   if (NULL == result) return(PIL_NUL_PTR);
   if (argc <= 0) return(PIL_NOT_FOUND);
   if (NULL == argv) return(PIL_NUL_PTR);
   
   for (i = 1; i < argc; i++)
    { if (NULL == argv[i]) continue;
      strncpy(buf, argv[i], PIL_LINESIZE - 1);
      buf[PIL_LINESIZE - 1] = 0;
      if (NULL == (p = strchr(buf, '='))) continue;
      *p = 0;
      if (strcmp(buf, name)) continue;			/* named arg: Name=Value */
      *result = argv[i] + (strlen(buf) + 1);
      return(PIL_OK);
    }
   if (PIL_NO_POSITIONAL & PILSpecialMode) return(PIL_NOT_FOUND);
   if ((numpar >= 1) && (numpar < argc))		/* positional arg: Value */
     { if (NULL == argv[numpar]) return(PIL_NOT_FOUND);
       if (NULL != strchr(argv[numpar], '=')) return(PIL_NOT_FOUND);
       *result = argv[numpar];
       return(PIL_OK);       
     }
   return(PIL_NOT_FOUND);
 }


/**************************************************************************************

Function:	PIL_allow_spaces_cmndarg2value
Description:	return value of given parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
numpar (In)		position of parameter in question (used when no parameter with
			'name' can be found)
*result (Out)		returned value of parameter (this is always string)
argc			number of entries in argument list (usually from main)
argv			argument list (usually from main)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	flexible SPACE/TAB syntax

Warning:  *result returns pointer to an argv[] element (possibly with some offset).

**************************************************************************************/


int	PIL_find_equal(char *buf)	/* check if buf contains '=' (excepting some cases) */
 { int	i;				/* returns: -1: not found, >=0: offset at which found */
   char	c0, c1;

   if (NULL == buf)  return(-1);

   for (i = 0;; i++)
    {
      c0 = buf[i];
      if (0 == c0)  break;		/* EOS ?, yes so terminate, and return not found */

      c1 = buf[i+1];
      if ('=' == c0)			/* ignore cases:  ==  =<  =>  */
        {
          if (('=' != c1) && ('<' != c1) && ('>' != c1))  return(i);
          i++;
        }
      else if ('=' == c1)		/* ignore cases:  ==  <=  >=  !=  */
        {
          if (('=' != c0) && ('<' != c0) && ('>' != c0) && ('!' != c0))  continue;
          i++;
        }
    }

   return(-1);
 }
    
/* Codes to classify arguments. */
enum { UNKNOWN = 0, PAR = 1, EQUALS = 2, VALUE = 4, DONE = 8 };

static void classify_arg(PIL_PFILE * fp, int arg_index, unsigned char * type, int * par_index, char ** value) {
  int end_par = fp->pcnt;
  int index;
  char * arg_p = fp->argv[arg_index];

  /* Skip leading white space in argument. */
  while (0 != isspace(*arg_p)) ++arg_p;

  /* Initially classify as unknown, with parameter index one-past-last par, and value is the first significant character. */
  *type = UNKNOWN;
  *par_index = end_par;
  *value = arg_p;

  /* Check argument against array of parameters to see if argument is of one of the forms PAR, PAR =, or PAR = VALUE. */
  for (index = 0; index < end_par; ++index) {
    char * par_p = 0;

    /* Skip any parameter which is not in the standard named parameter format. */
    if (PIL_FORMAT_OK != fp->pp[index].format) continue;

    /* Get pointer to the name of this parameter. */
    par_p = fp->pp[index].strname;

    /* Always start looking for parameter at the first significant character. */
    arg_p = *value;

    /* Compare parameter to argument. Break when a difference is deteected or the end of either string. */
    for (; '\0' != *arg_p && '\0' != *par_p && toupper(*arg_p) == toupper(*par_p); ++arg_p, ++par_p) {}

    /* If the parameter name was consumed, this argument *may* match this parameter. */
    if ('\0' == *par_p) {
      /* To match for sure the next character in the argument must be either white space or equals sign. */
      while (0 != isspace(*arg_p)) ++arg_p;
      if ('\0' == *arg_p || '=' == *arg_p) {
        /* This argument holds the name of a parameter. */
        *type |= PAR;
        *par_index = index;
        /* The value should be reset to point to the first significant character after the parameter name. */
        *value = arg_p;
        break;
      }
    }
  }

  /* Look for an equals sign. */
  if ('=' == *arg_p) {
    *type |= EQUALS;

    /* Look past equals sign as well as following white space. */
    ++arg_p;
    while (0 != isspace(*arg_p)) ++arg_p;
  }

  /* If there are any significant characters in the argument, the argument contains a value. */
  if ('\0' != *arg_p) *type |= VALUE;
}

/* Return flag indicating whether this character type contains a non-trivial string
   which is *not* a parameter name. */
static char is_significant_value(unsigned char type) {
  return EQUALS == type || (EQUALS | VALUE) == type || VALUE == type;
}

int PIL_allow_spaces_cmndarg2value(PIL_PFILE * fp) {
  int status = PIL_OK;
  if (0 == fp) {
    status = PIL_NUL_PTR;
  } else if (fp->argc > 1 && fp->pcnt > 0) { /* Only parse if command line has > 1 arg and if there is > 0 pars). */
    int end_arg = fp->argc;
    int end_par = fp->pcnt;
    /* Classifications of each argument. Include room for two extra arguments, because explicit parameters
       can involve an argument plus up to two arguments after it. This is initialized to 0 == UNKNOWN type. */
    unsigned char * arg_type = (unsigned char *) calloc(end_arg + 2, sizeof(unsigned char));
    /* Argument number indexed by parameter number. These must be initialized! */
    int * arg_par = (int *) calloc(end_par, sizeof(int));
    /* Parameter number indexed by argument number. These must be initialized! */
    int * par_arg = (int *) calloc(end_arg, sizeof(int));
    /* Pointers to the values. These are initialized to 0. */
    char ** value = (char **) calloc(end_arg, sizeof(char *));

    if (0 == arg_type || 0 == arg_par || 0 == par_arg || 0 == value) {
      free(value);
      free(par_arg);
      free(arg_par);
      free(arg_type);
      status = PIL_NO_MEM;
    } else {
      int arg_index;
      int par_index;

      /* Initialize arguments indexed by parameter number and parameters indexed by argument number to be
         one-past-last of each list. (Default is no matching has occurred. */
      for (arg_index = 0; arg_index < end_arg; ++arg_index) par_arg[arg_index] = end_par;
      for (par_index = 0; par_index < end_par; ++par_index) arg_par[par_index] = end_arg;

      /* Where possible, match up sequences of arguments to form expressions of the form
         PAR EQUALS VALUE. Any argument which fails to form part of such an expression is re-classified
         as a pure value. Go backwards, because that way a re-classified value may still be assigned to
         a valid expression from earlier in the argument list. */
      for (arg_index = end_arg - 1; arg_index > 0 && PIL_OK == status; --arg_index) {
        /* Classify argument into a form like PAR=VALUE, PAR=, =VALUE, etc. based only on the contents
           of the argument and the set of parameters in this file. */
        classify_arg(fp, arg_index, arg_type + arg_index, par_arg + arg_index, value + arg_index);

        /* If an argument appears to involve a parameter name, further checking is needed to make sure
           that makes sense in context. */
        if (0 != (PAR & arg_type[arg_index])) {
          if ((PAR | EQUALS | VALUE) == arg_type[arg_index] || ((PAR | EQUALS) == arg_type[arg_index] &&
            0 == is_significant_value(arg_type[arg_index + 1]))) {
            /* PAR=VALUE. Adjust value to point to the first significant character after EQUALS. */
            ++value[arg_index];
            while (0 != isspace(*value[arg_index])) ++value[arg_index];
          } else if ((PAR | EQUALS) == arg_type[arg_index] && 0 != is_significant_value(arg_type[arg_index + 1])) {
            /* PAR= VALUE or PAR= =VALUE or PAR= =. */
            arg_type[arg_index + 1] = DONE;
            value[arg_index] = value[arg_index + 1];
          } else if (PAR == arg_type[arg_index] && (EQUALS | VALUE) == arg_type[arg_index + 1]) {
            /* PAR =VALUE. */
            ++value[arg_index + 1];
            while (0 != isspace(*value[arg_index + 1])) ++value[arg_index + 1];
            arg_type[arg_index + 1] = DONE;
            value[arg_index] = value[arg_index + 1];
          } else if (PAR == arg_type[arg_index] && EQUALS == arg_type[arg_index + 1] &&
            0 != is_significant_value(arg_type[arg_index + 2])) {
            /* PAR = VALUE or PAR = =VALUE or PAR = =. */
            arg_type[arg_index + 1] = DONE;
            arg_type[arg_index + 2] = DONE;
            value[arg_index] = value[arg_index + 2];
          } else if ((PAR | EQUALS) == arg_type[arg_index]) {
            /* PAR= (blank). */
            ++value[arg_index];
            while (0 != isspace(*value[arg_index])) ++value[arg_index];
          } else if (PAR == arg_type[arg_index] && EQUALS == arg_type[arg_index + 1]) {
            /* PAR = (blank). */
            ++value[arg_index + 1];
            while (0 != isspace(*value[arg_index + 1])) ++value[arg_index + 1];
            arg_type[arg_index + 1] = DONE;
            value[arg_index] = value[arg_index + 1];
          } else if (0 != (PAR & arg_type[arg_index])) {
            /* None of the possible explicit cases were matched, so re-classify this parameter as a value. */
            arg_type[arg_index] = VALUE;
            arg_par[par_arg[arg_index]] = end_arg;
            par_arg[arg_index] = end_par;
            value[arg_index] = fp->argv[arg_index];
            continue;
          } else {
            continue;
          }
          arg_type[arg_index] = DONE;
          arg_par[par_arg[arg_index]] = arg_index;
        }
      }

      if (PIL_OK == status && 0 == (PIL_NO_POSITIONAL & PILSpecialMode)) {
        /* Match up positional arguments with left-over (unassigned) parameters. */
        /* Find first eligible parameter by iterating over all parameters and skipping those which were not
           already assigned, or which are hidden. */
        for (par_index = 0; par_index != end_par &&
          (end_arg != arg_par[par_index] || 0 != (PIL_MODE_HIDDEN & fp->pp[par_index].mode)); ++par_index) {}

        /* Iterate over arguments, and assign each one which was not already assigned to a positional parameter. */
        for (arg_index = 1; arg_index < end_arg && par_index < end_par; ++arg_index) {
          if (DONE != arg_type[arg_index]) {
            arg_type[arg_index] = DONE;
            arg_par[par_index] = arg_index;
            par_arg[arg_index] = par_index;
            value[arg_index] = fp->argv[arg_index];
            /* Find next unassigned parameter. */
            for (; par_index != end_par &&
              (end_arg != arg_par[par_index] || 0 != (PIL_MODE_HIDDEN & fp->pp[par_index].mode)); ++par_index) {}
          }
        }
      }

      /* Interpretation is done. Perform assignments of arguments to parameters. */
      for (arg_index = 1; arg_index < end_arg; ++arg_index) {
        int ii;

        /* Check if arguments were left over. */
        if (DONE != arg_type[arg_index]) {
          status = PIL_BOGUS_CMDLINE;
          continue;
        }

        /* Check if parameter was already assigned; if so, one of the other par_arg's which came before
           will equal the current par_arg. */
        for (ii = 1; ii < arg_index; ++ii) {
          if (par_arg[arg_index] == par_arg[ii] && end_par != par_arg[arg_index]) {
            status = PIL_BOGUS_CMDLINE;
            break;
          }
        }

        /* If this argument should be assigned to a parameter, but was not already assigned, assign it. */
        if (end_par != par_arg[arg_index] && ii == arg_index) {
          int local_status = PIL_OK;
          PIL_free(fp->pp[par_arg[arg_index]].strvalue);
          local_status = PIL_dup(&fp->pp[par_arg[arg_index]].strvalue, value[arg_index]);
          if (PIL_OK != local_status) {
            status = PIL_OK == status ? local_status : status;
          }
          /* Flag this parameter as changed on the command line. */
          if (PIL_OK == status) {
            fp->pp[par_arg[arg_index]].attrib |= PIL_VALUE_FROM_CMDLINE;
            fp->pp[par_arg[arg_index]].modified |= PIL_MODIFIED_VALUE;
          }
        }
      }

      /* Flag overall parameter file as modified. */
      if (PIL_OK == status && 0 == (PIL_FLUSH_PERFORMED & fp->mode))
        fp->mode |= PIL_PARFILE_MODIFIED;

      /* Clean up. */
      free(value);
      free(par_arg);
      free(arg_par);
      free(arg_type);
    }
  }
  return status;
}


/**************************************************************************************

Function:	PIL_verify_access
Description:	verify access mode of the file specified by path/name
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of file in question
path (In)		path of file in question
mode (In)		access mode to be checked for.
modtime (Out)		file's last modification time

Return code:	If given access mode is granted returns PIL_OK, PIL_NO_FILE,
		PIL_FILE_NO_RD, PIL_FILE_NO_WR if access is denied.
		Other errors mean bad arguments, out of memory, etc

Notes:	Possible value for mode are :

		PIL_FILECHECK_EXIST	check whether file exists
		PIL_FILECHECK_WRITE	check whether file is writable (plus mod time)
		PIL_FILECHECK_READ	check whether file readable (returns also
					last modification time)
		PIL_FILECHECK_CREATE	check whether file can be created (does
					_NOT_ return last modification time !!!)

	This function is called by PIL_build_names.

Warning:PIL_FILECHECK_XXX are also used to describe test mode for file-type
	variables (see PIL_parse_line)

**************************************************************************************/

int	PIL_verify_access(const char *name, const char *path, int mode, time_t *modtime)
 { char *p;
   int l, r;
   FILE *fp;
/* SCREW 1145: Use of stat is encapsulated in PIL_get_file_size now.
   struct stat statbuf;
*/

   if (NULL == name) return(PIL_NUL_PTR);
   if (NULL == path) return(PIL_NUL_PTR);
   if (NULL == modtime) return(PIL_NUL_PTR);
   if (0 == name[0]) return(PIL_NO_FILE);	/* name must be specified */
   if (0 == path[0]) return(PIL_NO_FILE);	/* path must be specified */

   r = 0;
   l = strlen(name) + strlen(path) + 2;		/* +2 is for '/' and EOS */
   if (NULL == (p = (char *)PIL_malloc(l))) return(PIL_NO_MEM);
   strcpy(p, path);
   strcat(p, sPathDelimString);
   strcat(p, name);
/* SCREW 1145: Use of access is encapsulated in PIL_file_exists now.
   if (mode & PIL_FILECHECK_EXIST)  if (access(p, F_OK)) r = PIL_NO_FILE;
*/
   if (mode & PIL_FILECHECK_EXIST) if (!PIL_file_exists(p)) r = PIL_NO_FILE;
   if (mode & PIL_FILECHECK_READ)
     { fp = fopen(p, "r");
       if (NULL == fp) r = PIL_FILE_NO_RD;
       else
        { fclose(fp);
/* SCREW 1145: Use of stat is encapsulated in PIL_get_mod_time.
          if (0 == stat(p, &statbuf)) { *modtime = statbuf.st_mtime; } 
          else { r = PIL_NO_FILE; }
*/
          r = PIL_get_mod_time(p, modtime);
        }
     }
   if (mode & PIL_FILECHECK_WRITE)
     { fp = fopen(p, "r+");
       if (NULL == fp) r = PIL_FILE_NO_WR;
       else
        { fclose(fp);
/* SCREW 1145: Use of stat is encapsulated in PIL_get_mod_time.
          if (0 == stat(p, &statbuf)) { *modtime = statbuf.st_mtime; } 
          else { r = PIL_NO_FILE; } 
*/
          r = PIL_get_mod_time(p, modtime);
        }
     }
   if (mode & PIL_FILECHECK_CREATE)
     { fp = fopen(p, "w");
       if (NULL == fp) r = PIL_FILE_NO_WR;
       else
        { fclose(fp);
/* SCREW 1145: Use ISO-standard remove in lieu of unlink. */
          remove(p);		/* bugfix, JB, 17-Jul-2000 */
        }
     }
   PIL_free(p);
   return(r);   
 }


/**************************************************************************************

Function:	PIL_traverse_token
Description:	find a file in a path
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of file to be found
path (In/Out)		comma separated list of directories (modified)
*result (Out)		directory in which file was found
mode (In)		required access mode of file

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	'path' is scanned from left to right for a file 'name' with access mode 'mode'.
First matching occurence is returned in *result (directory name only).
If no match is found function returns PIL_NOT_FOUND error code.

Warning:returned string in *result is actually part of 'path'

**************************************************************************************/

int	PIL_traverse_token(const char *name, char *path, char **result, int mode)
 { char		*p;
   time_t	modtime;

   if (NULL == name) return(PIL_NUL_PTR);
   if (NULL == path) return(PIL_NUL_PTR);
   if (NULL == result) return(PIL_NUL_PTR);

   for (p = strtok(path, sPathFieldSep); NULL != p; p = strtok(NULL, sPathFieldSep))
    { if (PIL_OK == PIL_verify_access(name, p, mode, &modtime))
        { *result = p;
          return(PIL_OK);
        }
    }
   return(PIL_NOT_FOUND);
 }


/**************************************************************************************

Function:	PIL_build_names
Description:	find parameter file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter file
argc (In)		number of items in argument list (usually from main)
argv			argument list (usually from main)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	Name of found (and to be used) parameter file is stored in static variable
	PILUsrPFile. If no parameter file exists it is copied from system R/O parameter
	file. PFILES environment variable stores search path for system and user
	parameter files. Paths are semicolon separated. Each path is colon separated
	list of directories. First is user path, second is system path. If only 1
	search path is specified it is assumed that user path = system path. Both
	relative and absolute directories are allowed. If no match is found error
	code is returned. Usually

		parameter filename = argv[0] + ".par"
	
	but which can be overriden when 'name' is non-NULL. Also application can
	disable scanning of PFILES by setting PILSpecialMode variable with PIL_BYBASS_NAMING
	flag set on.

Warning:results are stored in global static variables

**************************************************************************************/

int	PIL_build_names(char *name, int argc, char **argv)
 { char 	psys[2 * PIL_LINESIZE];
   char 	pusr[2 * PIL_LINESIZE];
   time_t	psystime, pusrtime;
   char		buf2[PIL_LINESIZE];
   char		buf[PIL_LINESIZE];
   char		pfilesbuf[PIL_LINESIZE];
   char		*p, *p2, *p3;
   int		l, i, r, copyflg, sysparexists;

   if (argc <= 0) return(PIL_BAD_ARG);	/* some simple arg validation */
   if (NULL == argv) return(PIL_NUL_PTR);

   if (PIL_BYPASS_NAMING & PILSpecialMode)
     {
       strncpy(PILUsrPFile, PILModuleName, PIL_LINESIZE);
       PILUsrPFile[PIL_LINESIZE - 1] = 0;
       return(PIL_OK);
     }

   if (NULL == name) p = argv[0];	/* typically we take par file name from exec name */
   else p = name;

   p2 = strrchr((char *)p, sPathDelim);   
   if (NULL != p2) p = p2;		/* filter out path, leave only filename */

   strncpy(buf, p, PIL_LINESIZE - 10);	/* 10 less - to accomodate future .par suffix */
   buf[PIL_LINESIZE - 10] = 0;		/* terminate string */

   i = 1;				/* signal append */
   l = strlen(buf);
   if (l > 4)				/* check if .par already appended */
     if (0 == memcmp(".par", buf + (l - 4), 4))
       i = 0;				/* signal do not append .par */

   if (i) strcat(buf, ".par");
   
   	/* ok, now we have valid filename, lets continue and calc path */

   copyflg = 0;				/* assume no copy needed */
   sysparexists = 0;			/* assume no system R/O parameter file exists */

   p = getenv("PFILES");		/* this is read only argument */
   if (NULL == p)  p = ".";		/* default to current dir if PFILES not set (again read-only) */
   strncpy(pfilesbuf, p, PIL_LINESIZE - 1); /* copy read-only arg to read-write buffer - strtok overwrites it !!! */
   pfilesbuf[PIL_LINESIZE - 1] = 0;	/* force end of string */
   p = pfilesbuf;			/* set pointer to the read-write copy */
   p2 = strchr(p, sPfilesDelim);	/* locate sys/usr separator */
   if (NULL != p2)			/* found separator (sys != usr) ? */
     { *p2 = 0;				/* yes, so process system part first */
       p2++;				/* p - usr part (left part), p2 - sys part (right part) */
     
       if (PIL_OK == (r = PIL_traverse_token(buf, p2, &p3, PIL_FILECHECK_EXIST)))
         { strncpy(psys, p3, PIL_LINESIZE - 1);
	   psys[PIL_LINESIZE - 1] = 0;
	   if (PIL_OK == (r = PIL_verify_access(buf, psys, PIL_FILECHECK_READ, &psystime)))
	     sysparexists = 1;		/* signal that sys parameter file can be read */
         }
     }

   strcpy(buf2, p);			/* process usr part */
   if (PIL_OK != (r = PIL_traverse_token(buf, p, &p3, PIL_FILECHECK_EXIST)))
     { if (PIL_OK != (r = PIL_traverse_token(buf, buf2, &p3, PIL_FILECHECK_CREATE))) return(r);
       if (0 == sysparexists) return(PIL_NO_FILE);
       copyflg = 1;			/* we need to copy sys -> usr, no usr file */
     }
   else					/* found parameter file, check if readable */
     { if (PIL_OK != (r = PIL_verify_access(buf, p3, PIL_FILECHECK_WRITE, &pusrtime))) return(r);
       if (sysparexists)		/* we need to copy sys -> usr, sys exists and modified later */
         if (psystime > pusrtime) copyflg = 1;
     }
   strncpy(pusr, p3, PIL_LINESIZE - 1);	/* write final path for parameter file */
   pusr[PIL_LINESIZE - 1] = 0;
/* jp needed temporarily because args to PIL_file_copy changed. */
#ifdef FOO
   if (copyflg)				/* psys/buf guaranteed to exist at this point */
     { if (PIL_OK != (r = PIL_file_copy(pusr, psys, buf))) return(r); /* cannot copy */
     }
#endif

   strcat(pusr, sPathDelimString);	/* append file name to path */
   strcat(pusr, buf);

   strncpy(PILUsrPFile, pusr, PIL_LINESIZE);  /* write output data */
   PILUsrPFile[PIL_LINESIZE - 1] = 0;

   return(PIL_OK);
 }

/**************************************************************************************

Function:	PIL_find_pfiles
Description:	find system and user parameter files in PFILES
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter file
userpfile (Out)		full path to user parameter file
syspfile (Out)		full path to system parameter file

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	PFILES environment variable stores search path for system and user
	parameter files. Paths are semicolon separated. Each path is colon separated
	list of directories. First is user path, second is system path. If only 1
	search path is specified it is assumed that user path = system path. Both
	relative and absolute directories are allowed. If no match is found error
	code is returned. Usually

		parameter filename = argv[0] + ".par"
	
	but which can be overriden when 'name' is non-NULL. Also application can
	disable scanning of PFILES by setting PILSpecialMode variable with PIL_BYBASS_NAMING
	flag set on.
	This function is based on PIL_build_names, which will be phased out
	as its functionality is repacked in PIL_find_pfiles and other functions.

**************************************************************************************/

int	PIL_find_pfiles(char *name, char *userpfile, char *syspfile)
 { char 	psys[2 * PIL_LINESIZE];
   char 	pusr[2 * PIL_LINESIZE];
   char		buf2[PIL_LINESIZE];
   char		buf[PIL_LINESIZE];
   char		pfilesbuf[PIL_LINESIZE];
   char		*p, *p2, *p3;
   int		l, i, r, sysparexists;

   if (NULL == name) return(PIL_NUL_PTR);
   if (NULL == userpfile) return(PIL_NUL_PTR);
   if (NULL == syspfile) return(PIL_NUL_PTR);

   syspfile[0] = 0;

   if (PIL_BYPASS_NAMING & PILSpecialMode)
     {
       strncpy(userpfile, PILModuleName, PIL_LINESIZE);
       userpfile[PIL_LINESIZE - 1] = 0;
       return(PIL_OK);
     }

   /* See if name contains a path delimiter, which would indicate it is a fully qualified
      parameter file name. */
   *pfilesbuf = '\0';
   *userpfile = '\0';
   if (0 != (p = strrchr(name, sPathDelim)))
     { size_t len = strlen(name);
       time_t modtime;
       char * dir_buf = (char *) malloc((1 + len) * sizeof(char));
       if (0 == dir_buf) return(PIL_NO_MEM);
       strncpy(dir_buf, name, p - name);

       /* To be a parameter file name, it must end in .par. */
       if (len >= 4 && ((name + len - 4) == strstr(name, ".par")))
         /* To be used conclusively as the user file at this point, the file must exist or be creatable. */
         { if (PIL_OK == PIL_verify_access(p, dir_buf, PIL_FILECHECK_EXIST, &modtime) ||
               PIL_OK == PIL_verify_access(p, dir_buf, PIL_FILECHECK_CREATE, &modtime))
             { len = len < PIL_LINESIZE - 1 ? len : PIL_LINESIZE - 1;
               strncpy(userpfile, name, len);
               userpfile[len] = '\0';
             }
         }
       free(dir_buf);
       /* Prepend directory to the path which will be scanned for parameter files; effectively, prepend it to PFILES. */
       len = p - name < PIL_LINESIZE - strlen(sPathFieldSep) - 1 ? p - name : PIL_LINESIZE - strlen(sPathFieldSep) - 1;
       strncat(pfilesbuf, name, len);
       strcat(pfilesbuf + len, sPathFieldSep); /* Explicitly add path delimiter. */
     }

   p = name;

   p2 = strrchr((char *)p, sPathDelim);   
   if (NULL != p2) p = p2;		/* filter out path, leave only filename */

   strncpy(buf, p, PIL_LINESIZE - 10);	/* 10 less - to accomodate future .par suffix */
   buf[PIL_LINESIZE - 10] = 0;		/* terminate string */

   i = 1;				/* signal append */
   l = strlen(buf);
   if (l > 4)				/* check if .par already appended */
     if (0 == memcmp(".par", buf + (l - 4), 4))
       i = 0;				/* signal do not append .par */

   if (i) strcat(buf, ".par");
   
   	/* ok, now we have valid filename, lets continue and calc path */

   sysparexists = 0;			/* assume no system R/O parameter file exists */

   p = getenv("PFILES");		/* this is read only argument */
   if (NULL == p)  p = ".";		/* default to current dir if PFILES not set (again read-only) */
   strncat(pfilesbuf, p, PIL_LINESIZE - 1 - strlen(pfilesbuf)); /* copy read-only arg to read-write buffer - strtok overwrites it !!! */
   pfilesbuf[PIL_LINESIZE - 1] = 0;	/* force end of string */
   p = pfilesbuf;			/* set pointer to the read-write copy */
   p2 = strchr(p, sPfilesDelim);	/* locate sys/usr separator */
   if (NULL != p2)			/* found separator (sys != usr) ? */
     { *p2 = 0;				/* yes, so process system part first */
       p2++;				/* p - usr part (left part), p2 - sys part (right part) */
     
       if (PIL_OK == (r = PIL_traverse_token(buf, p2, &p3, PIL_FILECHECK_EXIST)))
         { strncpy(psys, p3, PIL_LINESIZE - 1);
	   psys[PIL_LINESIZE - 1] = 0;
	   sysparexists = 1;		/* signal that sys parameter file exists */
         }
     }

   strcpy(buf2, p);			/* process usr part */
   if (PIL_OK != (r = PIL_traverse_token(buf, p, &p3, PIL_FILECHECK_EXIST)))
     { if (PIL_OK != (r = PIL_traverse_token(buf, buf2, &p3, PIL_FILECHECK_CREATE))) return(r);
       if (0 == sysparexists) return(PIL_NO_FILE);
     }

   strncpy(pusr, p3, PIL_LINESIZE - 1);	/* write final path for parameter file */
   pusr[PIL_LINESIZE - 1] = 0;

   strcat(pusr, sPathDelimString);	/* append file name to user par file path */
   strcat(pusr, buf);

   if ('\0' == *userpfile)
     { strncpy(userpfile, pusr, PIL_LINESIZE);  /* write output data */
       userpfile[PIL_LINESIZE - 1] = 0;
     }

   if (sysparexists)
     { strcat(psys, sPathDelimString);	/* append file name to system par file path */
       strcat(psys, buf);

       strncpy(syspfile, psys, PIL_LINESIZE); /* write more output data */
       syspfile[PIL_LINESIZE - 1] = 0;
     }

   return(PIL_OK);
 }

/**************************************************************************************

Function:	PIL_trim_spaces
Description:	remove leading/trailing spaces
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
*s (In)			string to check (pointer to a string)

Return code:	PIL_OK on success, err code otherwise

**************************************************************************************/

int	PIL_trim_spaces(char *s)
 { int	i, l;

   if (NULL == s) return(PIL_NUL_PTR);

   for (i=0; ' ' == s[i]; i++);
   l = strlen(s + i);
   if (i > 0)  memmove(s, s + i, l + 1);
   for (i = l - 1; i >= 0; i--)
    { 
      if (' ' != s[i])  break;			/* BUGFIX : 06-Mar-2002 pence@tetra.gsfc.nasa.gov */
      s[i] = 0;
    }

   return(PIL_OK);
 }


/**************************************************************************************

Function:	PIL_stricmp
Description:	compare string ignoring case
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
*s1 (In)		1st string (pointer to a string)
*s2 (In)		2nd string

Return code:	zero if strings are equal, nonzero otherwise

Notes:	C locale is assumed and 7bit ASCII charset

**************************************************************************************/

int	PIL_stricmp(const char *s1, const char *s2)
 { int	c1, c2;

   if (NULL == s1) return(2);
   if (NULL == s2) return(2);

   for (;;)
    { c1 = *s1;
      if ((c1 >= 'A') && (c1 <= 'Z'))  c1 += 'a' - 'A';
      c2 = *s2;
      if ((c2 >= 'A') && (c2 <= 'Z'))  c2 += 'a' - 'A';
      
      if (c1 < c2) return(-1);
      if (c1 > c2) return(1);
      if (0 == c1) break;
      s1++;
      s2++;
    }
   return(0);
 }


/**************************************************************************************

Function:	PIL_dup
Description:	safe string copy
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
*dst (Out)		destination string (pointer to a string)
src (In)		source string

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	Checks for various error conditions. Result returned in *dst. Use PIL_free
	to deallocate *dst.

**************************************************************************************/

int	PIL_dup(char **dst, const char *src)
 { int l;

   if (NULL == dst) return(PIL_NUL_PTR);
   if (NULL == src)  { *dst = NULL; return(PIL_OK); }
   l = strlen(src);
   *dst = (char *)PIL_malloc(l + 1);
   if (NULL == *dst) return(PIL_NO_MEM);
   memcpy(*dst, src, l + 1);	/* this is string copy, but we know length already */
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PIL_type2string
Description:	return ascii string name for given parameter data type
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
type (In)		data type (integer)

Return code:	ascii name for given data type or " -???- " string

Notes:	uses global conversion table. " -???- " should be the last entry in this table

**************************************************************************************/

const char *PIL_type2string(int type)
 { int i;
   
   for (i=0; PIL_ttab[i].type != PIL_TYPE_UNKNOWN; i++)
     if (type == PIL_ttab[i].type) break;
   return(PIL_ttab[i].name);
 }


/**************************************************************************************

Function:	PIL_mode2string
Description:	return ascii string for given parameter mode
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
mode (In)		parameter mode (integer)
but (Out)		parameter mode in ASCII representation


Return code:	PIL_OK or error code

Notes:	Output string is exactly 4 bytes long (plus terminating zero). Letters denote :

		a - auto
		h - hidden
		l - learn
		q - query

	Dash (minus sign) in place of letter means given mode is not set.

**************************************************************************************/

int	PIL_mode2string(int mode, char *buf)
 {
   if (NULL == buf)  return(PIL_NUL_PTR);

   buf[0] = ((mode & PIL_MODE_AUTO) ? 'a' : '-');
   buf[1] = ((mode & PIL_MODE_HIDDEN) ? 'h' : '-');
   buf[2] = ((mode & PIL_MODE_LEARN) ? 'l' : '-');
   buf[3] = ((mode & PIL_MODE_QUERY) ? 'q' : '-');
   buf[4] = 0;

   return(PIL_OK);
 }


/**************************************************************************************

Function:	PIL_indir2string
Description:	dereference indirections
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			parameter file
what (In)		type of operation
*dst (Out)		result after all indirections are resolved
src (In)		input string to be resolved

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	given input string (usually value of parameter as read from parameter file)
	dereference any indirections. Currently the following syntax is supported :

		$xxx	is substituted with value of environment variable xxx
		)xxx	is substituted with value/minvalue/maxvalue (given by 'what')
			of parameter xxx. Substitutions may be disabled, when
			pil_ext_indir parameter is not present in parameter file.

Warning:if error occurs during attempt to resolve indirection (say $xxx points to
	nonexisting envirionment variable) original string is returned.

**************************************************************************************/


int     pil_chridx(char *p, char c)
 { int  l;

   for (l = 0; p[l]; l++)  if (c == p[l]) return(l);
   return(-1);
 }


int	pil_curly_append(char **p, char *s, int nchars)
 { int	l, nl;
   char	*np;

   if (NULL == p)  return(PIL_NUL_PTR);
   if (nchars < 0) return(PIL_BAD_ARG);
   if (0 == nchars)
     { if (NULL != *p) return(PIL_OK); }
   else
     { if (NULL == s)  return(PIL_NUL_PTR); }

   if (NULL == *p)  { l = 0; }
   else  { l = strlen(*p); }
   
   nl = nchars + l + 1;						/* + 1 - for EOS */
   
   if (NULL == *p)  { np = (char *)PIL_malloc(nl); }
   else  { np = (char *)PIL_realloc((void *)*p, nl); }
   if (NULL == np)
     { if (NULL != *p)  { PIL_free(*p); *p = NULL; }		/* on error deallocate string */
       return(PIL_NO_MEM);
     }

   if (nchars > 0) memcpy(np + l, s, nchars);			/* append new string */
   np[nl - 1] = 0;						/* signal EOS */
   
   *p = np;
   return(PIL_OK);
 }


int     pil_curly_expand(char *src, char **dst)
 { char		ev[PIL_CURLYSIZE];
   int		sl, l, l2, l3, tl, r;
   char		*p;


   if (NULL == src) return(PIL_NUL_PTR);
   if (NULL == dst) return(PIL_NUL_PTR);
   *dst = NULL;

   sl = tl = 0;
   for (;;)
    { if (-1 == (l = pil_chridx(src + sl + tl, '$'))) break;	/* no '$' from current pos, so terminate */
      if ('{' != src[sl + tl + l + 1])
        { tl += l + 1;
          continue;						/* '$' not followed by '{' */
        }

      l += tl;							/* total number of bytes _BEFORE_ '${...}' */
      if (-1 == (l2 = pil_chridx(src + sl + l + 2, '}'))) break;  /* no closing '}' */
      if (PIL_OK != (r = pil_curly_append(dst, src + sl, l))) return(r);  /* copy part _BEFORE_ '$' */
      sl += l + 2;						/* position on 1st char after '${' */

      if (l2 > 0)						/* if non-empty env.var name */
        { memcpy(ev, src + sl, l2);				/* then copy its name to the temporary buffer */
          ev[l2] = 0;						/* signal end of string */
	  p = getenv(ev);					/* get environment variable */
	  l3 = 0;
	  if (NULL != p)  l3 = strlen(p);			/* ... and compute length of that variable (0 - if notdef) */
	  if (PIL_OK != (r = pil_curly_append(dst, p, l3))) return(r);  /* copy value of environment variable */
	}
      sl += l2 + 1;						/* env.var.name + '}' */
      tl = 0;
    }

   return(pil_curly_append(dst, src + sl, strlen(src + sl)));	/* copy remaining part */
 }


int	PIL_indir2string(PIL_PFILE *fp, int what, char **dst, char *src)
 { char		*p, *p2;
   int		i, idx;


   if (NULL == dst) return(PIL_NUL_PTR);
   if (NULL == src) return(PIL_NUL_PTR);
   
   if (0 == (fp->mode & PIL_SUBST_CURLY))  return(pil_curly_expand(src, dst));  /* default is expand curly ! */

   p = src;
   
   for (i=0; i<PIL_MAX_INDIR; i++)
    { if ('$' == p[0])
        { if (0 == (fp->mode & PIL_SUBST_ENV)) break;
          if (NULL == (p2 = getenv(p+1))) break; /* if no env var, return previous indir */
          p = p2;
        }
      else if (')' == p[0])
        { if (0 == (fp->mode & PIL_SUBST_INDIR)) break;
          if (PIL_OK != PIL_find_name(fp, p+1, &idx)) break;
          switch (what)
           { case PIL_WHAT_MIN: p = fp->pp[idx].strmin; break;
             case PIL_WHAT_MAX: p = fp->pp[idx].strmax; break;
             default: p = fp->pp[idx].strvalue; break;
           }
        }
      else break;
    }   

   return(PIL_dup(dst, p));
 }


/**************************************************************************************

Function:	PIL_string2value
Description:	convert ascii string to typed value (like sscanf(...))
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
buf (In)		string to be converted
type (In)		type of output value
v (Out)			structure holding typed value

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	type can be one of the following values :

		PIL_TYPE_BOOL	accepts yes, y, no, n
		PIL_TYPE_INT4	any integer value
		PIL_TYPE_REAL8	any real value
		PIL_TYPE_STRING	any string
		PIL_TYPE_FNAME	any string (does not check whether file exists)

	See also PIL_value2string

Warning:other XPI/IRAF data types like GCursor, IMCur are not supported

**************************************************************************************/

int	PIL_string2value(const char *buf, int type, PIL_VALUE *v)
 { char c;
   int	my_len, r;

   if (NULL == buf) return(PIL_NUL_PTR);
   if (NULL == v) return(PIL_NUL_PTR);
   r = PIL_OK;
   switch (type)
    { case PIL_TYPE_BOOL:
		if ((!PIL_stricmp(buf, "yes")) || (!PIL_stricmp(buf, "y")))  { v->b = 1; }
		else if ((!PIL_stricmp(buf, "no")) || (!PIL_stricmp(buf, "n")))  { v->b = 0; }
		else r = PIL_BAD_VAL_BOOL;
		break;
      case PIL_TYPE_INT4:
		if (1 != sscanf(buf, "%d %c", &(v->i), &c))  r = PIL_BAD_VAL_INT;
		break;
      case PIL_TYPE_REAL8:
		if (1 != sscanf(buf, "%lg %c", &(v->d), &c))  r = PIL_BAD_VAL_REAL;
		break;
      case PIL_TYPE_STRING:
		switch (PILVectorType)
		 { case PIL_TYPE_INT4:
				if (PILVectorVarMode)
				  { my_len = PILVectorLen;
				    if (PIL_OK != PILCheckIntVarVector(buf, &my_len, NULL)) r = PIL_BAD_VAL_INT_VAR_VECTOR;
				  }
				else
				  { if (PIL_OK != PILCheckIntVector(buf, PILVectorLen, NULL)) r = PIL_BAD_VAL_INT_VECTOR; }
				break;
		   case PIL_TYPE_REAL4:
				if (PILVectorVarMode)
				  { my_len = PILVectorLen;
				    if (PIL_OK != PILCheckReal4VarVector(buf, &my_len, NULL)) r = PIL_BAD_VAL_REAL_VAR_VECTOR;
				  }
				else
				  { if (PIL_OK != PILCheckReal4Vector(buf, PILVectorLen, NULL)) r = PIL_BAD_VAL_REAL_VECTOR; }
				break;
		   case PIL_TYPE_REAL8:
				if (PILVectorVarMode)
				  { my_len = PILVectorLen;
				    if (PIL_OK != PILCheckRealVarVector(buf, &my_len, NULL)) r = PIL_BAD_VAL_REAL_VAR_VECTOR;
				  }
				else
				  { if (PIL_OK != PILCheckRealVector(buf, PILVectorLen, NULL)) r = PIL_BAD_VAL_REAL_VECTOR; }
				break;
		 }
		strncpy(v->s, buf, PIL_LINESIZE);
		v->s[PIL_LINESIZE - 1] = 0;
		break;
      case PIL_TYPE_FNAME:
      case PIL_TYPE_DOL:
		strncpy(v->s, buf, PIL_LINESIZE);
		v->s[PIL_LINESIZE - 1] = 0;
		PIL_trim_spaces(v->s);
		break;
      case PIL_TYPE_STRUCT:
      case PIL_TYPE_GCUR:
      case PIL_TYPE_IMCUR:
      case PIL_TYPE_D:
      case PIL_TYPE_G:
      default:	r = PIL_NOT_IMPLEMENTED;	/* we silently ignore unrecognized types ... */
    }
   return(r);
 }


/**************************************************************************************

Function:	PIL_value2string
Description:	convert typed value to ascii string (like sprintf(...))
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
v (In)			structure holding typed value to be converted
type (In)		type of output value
buf (Out)		output buffer
maxlen (In)		size of output buffer in bytes

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	type can be one of the following values :

		PIL_TYPE_BOOL	accepts yes, y, no, n
		PIL_TYPE_INT4	any integer value
		PIL_TYPE_REAL8	any real value
		PIL_TYPE_STRING	any string
		PIL_TYPE_FNAME	any string (does not check whether file exists)

	See also PIL_string2value

Warning:other XPI/IRAF data types like GCursor, IMCur are not supported

**************************************************************************************/

int	PIL_value2string(PIL_VALUE *v, int type, char *buf, int maxlen)
 { int	r;

   if (NULL == buf) return(PIL_NUL_PTR);
   if (NULL == v) return(PIL_NUL_PTR);
   if (maxlen <= 0) return(PIL_BAD_ARG);
   r = PIL_OK;
   switch (type)
    { case PIL_TYPE_BOOL:
		strcpy(buf, (v->b ? "yes" : "no"));  break;
      case PIL_TYPE_INT4:
		sprintf(buf, "%d", v->i);  break;
      case PIL_TYPE_REAL8:
		sprintf(buf, "%.15g", v->d);  break;
      case PIL_TYPE_STRING:
      case PIL_TYPE_FNAME:
		strncpy(buf, v->s, maxlen);
		buf[maxlen - 1] = 0;
		break;
      case PIL_TYPE_STRUCT:
      case PIL_TYPE_GCUR:
      case PIL_TYPE_IMCUR:
      case PIL_TYPE_D:
      case PIL_TYPE_G:
      default:	r = PIL_NOT_IMPLEMENTED;	/* we silently ignore unrecognized types ... */
    }
   return(r);
 }


/**************************************************************************************

Function:	PIL_check_mode
Description:	convert string representing mode for parameter to binary representation
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
ptab (In)		parameter mode (ascii string)
*mode (Out)		parameter mode (integer)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	Input string can be in one of 2 allowable formats:

	Short format : only following characters are allowed :

		'a'		auto
		'h'		hidden
		'l'		learn
		'q'		query
		' '		separator
		'\t'		separator

	Any other characters cause function to terminate with error PIL_BAD_MODE.
	Allowable characters can be in any order, and can be repeated many times.

	Long format : ptab string should contain the following words separated by
	whitespaces :
	
		'auto'
		'hidden'
		'learn'
		'query'
	
	Words can be in any order and can be repeated many times.

**************************************************************************************/

int	PIL_check_mode(char *ptab, int *mode)
 { int		i, r, nbyt;
   char 	mytok[PIL_LINESIZE];
   char		*p;
	
		/* check short mode */
   r = PIL_OK;
   for (p = ptab;;p++)
    { if (0 == *p) break;
      switch (*p)
       { case 'a' : *mode |= PIL_MODE_AUTO; break;
         case 'h' : *mode |= PIL_MODE_HIDDEN; break;
         case 'l' : *mode |= PIL_MODE_LEARN; break;
         case 'q' : *mode |= PIL_MODE_QUERY; break;
	 case ' ' : break;
	 case '\t': break;
         default  : r = PIL_BAD_MODE; 
       }
    }
		/* check long mode */

   if (PIL_OK != r)	/* not in short format, try long format */
     { r = PIL_OK;
       *mode = 0;
       p = ptab;	/* start from beginning */
       for (i=0;;i++)	/* translate '+' into ' ' */
        { if (0 == p[i]) break;
          if ('+' == p[i]) p[i] = ' ';
        }
       for (;;)
        { if (1 != sscanf(p, "%s%n", mytok, &nbyt)) break;
	  if (!strcmp(mytok, "auto")) { *mode |= PIL_MODE_AUTO; }
	  else if (!strcmp(mytok, "hidden")) { *mode |= PIL_MODE_HIDDEN; }
	  else if (!strcmp(mytok, "learn")) { *mode |= PIL_MODE_LEARN; }
	  else if (!strcmp(mytok, "query")) { *mode |= PIL_MODE_QUERY; }
	  else return(PIL_BAD_MODE);
	  p += nbyt;
        }
     }
   return(r);
 }


/**************************************************************************************

Function:	PIL_check_pfile
Description:	check parameter file structure
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This functions performs generic test on structure pointed to by fp. It tests
	whether fields within structure have "reasonable" values. Many other PIL
	functions call this function. Successful call to this functions means that
	fp pointer is more or less (probably more more than less) valid.

Warning:test is by no means exhaustive.

**************************************************************************************/

int	PIL_check_pfile(PIL_PFILE *fp)
 { if (NULL == fp) return(PIL_NUL_PTR);
/*   if (NULL == fp->fp) return(PIL_BAD_ARG); */
   if (fp->pcnt < 0) return(PIL_BAD_ARG);
   if (fp->argc < 0) return(PIL_BAD_ARG);
   if (fp->pcnt) if (NULL == fp->pp) return(PIL_BAD_ARG);
   if (fp->argc) if (NULL == fp->argv) return(PIL_BAD_ARG);
   /* SCREW 1145: Check new file name field for validity. */
   if (NULL == fp->fname || 0 == *fp->fname) return(PIL_BAD_ARG);
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PIL_check_par
Description:	generic check of parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
idx (In)		index of parameter in question

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	function tests whether all fields for given parameter (in table of parameters
	managed by PIL in memory) are non-NULL.

**************************************************************************************/

int	PIL_check_par(PIL_PFILE *fp, int idx)
 { int r;

   if (PIL_OK != (r = PIL_check_pfile(fp))) return(r);
   if ((idx < 0) || (idx >= fp->pcnt)) return(PIL_BAD_ARG);
   if (PIL_FORMAT_OK != fp->pp[idx].format) return(PIL_BAD_ARG);
   if (NULL == fp->pp[idx].strname) return(PIL_BAD_ARG);
   if (NULL == fp->pp[idx].strmode) return(PIL_BAD_ARG);
   if (NULL == fp->pp[idx].strtype) return(PIL_BAD_ARG);
   if (NULL == fp->pp[idx].strvalue) return(PIL_BAD_ARG);
   if (NULL == fp->pp[idx].strmin) return(PIL_BAD_ARG);
   if (NULL == fp->pp[idx].strmax) return(PIL_BAD_ARG);
   if (NULL == fp->pp[idx].strprompt) return(PIL_BAD_ARG);
   return(r);
 }

/**************************************************************************************

Function:	PIL_merge_pfiles
Description:	merge (if possible) the system and user parameter files
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
spfile (In)		pointer to first parameter file structure
upfile (In/Out)		pointer to second parameter file structure

Return code:	standard PIL error code or PIL_OK (zero)

**************************************************************************************/

int	PIL_merge_pfiles(PIL_PFILE *srcfile, PIL_PFILE *destfile)
  { int r = PIL_OK;
    int modified = 0;
    PIL_PARAM **src = NULL;
    PIL_PARAM **dest = NULL;
    PIL_PARAM **srcend;
    PIL_PARAM **destend;
    PIL_PARAM **sp;
    PIL_PARAM **dp;
    int comparison;

    if (NULL == srcfile) return(PIL_NUL_PTR);
    if (NULL == destfile) return(PIL_NUL_PTR);

    do {
      /* Create sorted arrays for the system and user parameter files. */
      r = PIL_sort_parameters(srcfile, &src);
      if (PIL_OK != r) continue;

      r = PIL_sort_parameters(destfile, &dest);
      if (PIL_OK != r) continue;

      /* Determine start and end points. */
      sp = src; dp = dest;
      srcend = src + srcfile->pcnt; destend = dest + destfile->pcnt;
      
      /* Loop through parameters until the end, or until the last
         one with a valid format is reached. */
      while (sp < srcend && PIL_FORMAT_OK == (*sp)->format &&
        dp < destend && PIL_FORMAT_OK == (*dp)->format)
        { comparison = strcmp((*sp)->strname, (*dp)->strname);
          if (comparison > 0)
            { modified = 1;
              /* System parameter > user parameter, means user
                 parameter not valid in the system file, so remove it. */
              ++dp;
              continue;
            }
          else if (comparison < 0)
            { modified = 1;
              /* System parameter < user parameter, means user
                 file is missing a parameter which is in the system file,
                 so add it. */
              ++sp;
              continue;
            }
          /* Name of the two parameters is the same, so compare
             values. If the values are different, make sure everything
             else about the parameters is the same before copying the
             source value to the destination. */
          if (strcmp((*sp)->strvalue, (*dp)->strvalue))
            { if ((*sp)->type == (*dp)->type && (*sp)->mode == (*dp)->mode &&
                !strcmp((*sp)->strmin, (*dp)->strmin) &&
                !strcmp((*sp)->strmax, (*dp)->strmax))
                { PIL_free((*dp)->strvalue);
                  (*dp)->strvalue = (char *) malloc((strlen((*sp)->strvalue)+1)*sizeof(char));
                  if (NULL == (*dp)->strvalue) { r = PIL_NO_MEM; break; }
                  strcpy((*dp)->strvalue, (*sp)->strvalue);
                  /* SPR 3548: fix a bug which caused learned parameters to
                     be unlearned if they were not modified but others were.
                     The cause of this is that the strline field is used verbatim
                     to write such parameters if and when the output parameter
                     file is written. (For parameters which actually are modified
                     the line in the parameter file for that parameter is
                     reconstructed from scratch.) Before the following addition,
                     the strline field was not copied when merging parameter files.
                     This caused the whole line to revert to the value in the
                     system parameter file. */
                  PIL_free((*dp)->strline);
                  (*dp)->strline = (char *) malloc((strlen((*sp)->strline)+1)*sizeof(char));
                  if (NULL == (*dp)->strline) { r = PIL_NO_MEM; break; }
                  strcpy((*dp)->strline, (*sp)->strline);
                }
            }
          ++sp; ++dp;
        }
    } while (0);

    PIL_free(dest);
    PIL_free(src);

    return(r);
  }

		/* functions dealing with single parameter */


/**************************************************************************************

Function:	PIL_init_param
Description:	initialize parameter entry in parameter table.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
pp (In/Out)		pointer to parameter entry structure

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This function should be called when new parameter is to be added to parameter
	table, and table has been already reallocated to hold new entry. Subsequently
	fields should be filled in by application (see PIL_append_param)

**************************************************************************************/

int	PIL_init_param(PIL_PARAM *pp)
 { if (NULL == pp) return(PIL_NUL_PTR);
   pp->strline = NULL;
   pp->strname = NULL;
   pp->strtype = NULL;
   pp->strmode = NULL;
   pp->strvalue = NULL;
   pp->strmin = NULL;
   pp->strmax = NULL;
   pp->strprompt = NULL;
   pp->type = 0;
   pp->mode = 0;
   pp->attrib = 0;
   pp->modified = 0;
   pp->format = PIL_FORMAT_ERR;
   pp->reprompt = 0; /* SCREW 1496: new reprompt flag. */
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PIL_free_entry
Description:	make parameter entry ready to be deleted
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
pp (In/Out)		pointer to parameter entry structure

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	all non-NULL fields within parameter entry structure are deallocated.
	Application may subsequently remove this entry from parameter table (see
	also PIL_delete_param)

**************************************************************************************/

int	PIL_free_entry(PIL_PARAM *pp)
 { if (NULL == pp) return(PIL_NUL_PTR);
   if (NULL != pp->strline)  PIL_free(pp->strline);
   if (NULL != pp->strname)  PIL_free(pp->strname);
   if (NULL != pp->strtype)  PIL_free(pp->strtype);
   if (NULL != pp->strmode)  PIL_free(pp->strmode);
   if (NULL != pp->strvalue)  PIL_free(pp->strvalue);
   if (NULL != pp->strmin)  PIL_free(pp->strmin);
   if (NULL != pp->strmax)  PIL_free(pp->strmax);
   if (NULL != pp->strprompt)  PIL_free(pp->strprompt);
   return(PIL_init_param(pp));
 }


/**************************************************************************************

Function:	PIL_delete_param
Description:	delete given parameter from parameter table
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In/Out)		pointer to parameter file structure
idx (In)		index of parameter to be deleted

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	for performance reasons it is advisable to delete parameter from last to first

**************************************************************************************/

int	PIL_delete_param(PIL_PFILE *fp, int idx)
 { PIL_PARAM *p;
   int r;

   if (NULL == fp) return(PIL_NUL_PTR);
   if ((idx < 0) || (idx >= fp->pcnt)) return(PIL_BAD_ARG);
   if (NULL == fp->pp) return(PIL_BAD_ARG);	/* should never happen */
   
   if (PIL_OK != (r = PIL_free_entry(&(fp->pp[idx])))) return(r);
   if (idx < (fp->pcnt - 1))
     { memmove(fp->pp + (idx), fp->pp + (idx + 1), sizeof(PIL_PARAM) * (fp->pcnt - 1 - idx));
     }
   if (fp->pcnt > 1)			/* is this the last parameter ? */
     {  p = (PIL_PARAM *)PIL_realloc((void *)fp->pp, (fp->pcnt - 1) * sizeof(PIL_PARAM));
        if (NULL == p) r = PIL_NO_MEM;
        else fp->pp = p;
     }
   else					/* yes, so delete whole list ... */
     { PIL_free((void *)fp->pp);
       fp->pp = NULL;
     }
   (fp->pcnt)--;			/* signal one parameter less */
   return(r);
 }


/**************************************************************************************

Function:	PIL_append_param
Description:	append new parameter to parameter table
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In/Out)		pointer to parameter file structure
pp (In)			pointer to parameter structure

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	All fields within pp parameter structure are copied into fp, so application
	is free to call PIL_free_entry. Function is atomic, it is either succeeds
	or nothing is done.

**************************************************************************************/

int	PIL_append_param(PIL_PFILE *fp, PIL_PARAM *pp)
 { PIL_PARAM *p;
   int r;

   if ((NULL == fp) || (NULL == pp)) return(PIL_NUL_PTR);

   if (0 != fp->pcnt)
     { if (NULL == fp->pp) return(PIL_BAD_ARG);	/* should never happen */
       p = (PIL_PARAM *)PIL_realloc((void *)fp->pp, (fp->pcnt + 1) * sizeof(PIL_PARAM));
     }
   else
     { p = (PIL_PARAM *)PIL_malloc(sizeof(PIL_PARAM));
     }
   if (NULL == p) return(PIL_NO_MEM);
   if (PIL_OK != (r = PIL_init_param(p+(fp->pcnt)))) return(r);
   if (PIL_OK != (r = PIL_dup(&(p[fp->pcnt].strline), pp->strline))) goto EEX1;
   if (PIL_OK != (r = PIL_dup(&(p[fp->pcnt].strname), pp->strname))) goto EEX1;
   if (PIL_OK != (r = PIL_dup(&(p[fp->pcnt].strmode), pp->strmode))) goto EEX1;
   if (PIL_OK != (r = PIL_dup(&(p[fp->pcnt].strtype), pp->strtype))) goto EEX1;
   if (PIL_OK != (r = PIL_dup(&(p[fp->pcnt].strvalue), pp->strvalue))) goto EEX1;
   if (PIL_OK != (r = PIL_dup(&(p[fp->pcnt].strmin), pp->strmin))) goto EEX1;
   if (PIL_OK != (r = PIL_dup(&(p[fp->pcnt].strmax), pp->strmax))) goto EEX1;
   if (PIL_OK != (r = PIL_dup(&(p[fp->pcnt].strprompt), pp->strprompt))) goto EEX1;
   p[fp->pcnt].type = pp->type;
   p[fp->pcnt].mode = pp->mode;
   p[fp->pcnt].attrib = pp->attrib;
   p[fp->pcnt].modified = pp->modified;
   p[fp->pcnt].format = pp->format;
   p[fp->pcnt].reprompt = pp->reprompt; /* SCREW 1496: new reprompt flag. */

EEX1:
   if (PIL_OK != r) PIL_free_entry(p + (fp->pcnt));
   else
     { fp->pp =  p;
       (fp->pcnt)++;
     }
   return(r);
 }



	/* functions dealing with parameter file */


/**************************************************************************************

Function:	PIL_lock_pfile
Description:	lock parameter file for exclusive access
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This function is necessary when application runs in server mode and wants to
	flush/reread parameter file. It should be called just before application
	wants to read/write parameter file.

Warning:this function may hang application. In this case any signal catched will help
	This function assures mutual exclusion as long as other processes accessing
	parameter file use the same technique.

**************************************************************************************/

/* SCREW 1145: This function moved to pil_sys.c
int	PIL_lock_pfile(PIL_PFILE *fp, int mode)
 { flock_t flk;

   if (NULL == fp) return(PIL_NUL_PTR);

   if (PIL_RDWRITE & mode)  { flk.l_type = F_WRLCK; }
   else  { flk.l_type = F_RDLCK; } 
   flk.l_whence = 0;
   flk.l_start = 0;
   flk.l_len = 0;
   if (-1 == fcntl(fileno(fp->fp), F_SETLKW, &flk))  return(PIL_LOCK_FAILED);

   return(PIL_OK);
 }
*/


/**************************************************************************************

Function:	PIL_unlock_file
Description:	unlock parameter file for exclusive access
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)                 pointer to parameter file structure

Return code:    standard PIL error code or PIL_OK (zero)

Notes:  This function is necessary when application runs in server mode and wants to
        flush/reread parameter file. It should be called when application is done
        reading/writing parameter file
        
Warning:this function may hang application. In this case any signal catched will help
	This function assures mutual exclusion as long as other processes accessing
	parameter file use the same technique.

**************************************************************************************/

/* SCREW 1145: This function moved to pil_sys.c
int	PIL_unlock_pfile(PIL_PFILE *fp)
 { flock_t flk;

   if (NULL == fp) return(PIL_NUL_PTR);

   flk.l_type = F_UNLCK;
   flk.l_whence = 0;
   flk.l_start = 0;
   flk.l_len = 0;
   if (-1 == fcntl(fileno(fp->fp), F_SETLKW, &flk))  return(PIL_LOCK_FAILED);

   return(PIL_OK);
 }
*/

/* Prototype for new open utility. */
static int	PIL_pfile_open_only(char *fname, int argc, char **argv, int pfile_acc_mode, int numtries, PIL_PFILE *fp);

/**************************************************************************************

Function:	PIL_reload_parameters
Description:	reload parameters from parameter file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	this function is called during initial read of parameter file. Additionally
	it is called whenever application running in server mode wants to reread
	parameter file.

Warning:Function unconditionally deletes current parameter list in memory, without
	flushing it to disk, so it is necessary for server applications to call
	PIL_flush_parameters before calling this function

**************************************************************************************/

int	PIL_reload_parameters(PIL_PFILE *fp)
 { char 	buf[PIL_LINESIZE], buf2[PIL_LINESIZE];
   int		r, r2, i, idx, numpar, savemode, j;
   int		numtries = 5;
   PIL_PARAM	pp;
   
   if (PIL_OK != (r = PIL_check_pfile(fp))) return(r);

   for (i = fp->pcnt; i > 0; i--)  PIL_delete_param(fp, i - 1); /* delete parameters (if any) */

   numpar = 1;						/* 1st par is at idx 1, since argv[0] is exec name */
   if (PIL_PSET_MODE & PILSpecialMode)  numpar = 2;	/* but during pset syntax is: execname filename par1 par2 .. */

   savemode = fp->mode & PIL_FLUSH_PERFORMED;		/* save flush flag */

   r = PIL_lock_pfile(fp, PIL_RDONLY);			/* to be on the safe side, lock file for excl. access */
   if (NULL == fp->fp)
     { char * fname = fp->fname;
       r = PIL_pfile_open_only(fp->fname, fp->argc, fp->argv, PIL_RDONLY, numtries, fp);
       if (PIL_OK != r)
         { fp->fname = fname;				/* problem occurred, so restore state of fp */
           return(r);
         }
       free(fname);					/* fp->fname was changed, so prevent mem leak */
     }
   if (PIL_OK != r) { return(r); }

   fseek(fp->fp, 0L, SEEK_SET);				/* rewind to the beginning of file */
   for (j = 0;; j++)	
    { if (PIL_OK != PIL_read_line(fp, buf, j)) break;	/* break on EOF or error */
      switch (PIL_parse_line(fp, buf, numpar, &pp, fp->argc, fp->argv, buf2))
       { case PIL_OK : ;
		PIL_append_param(fp, &pp);		/* this may fail, but we ignore errcode */
		numpar++;				/* this counts number of args (do not count comment lines) */
		break;
         case PIL_LINE_BLANK: ;
         case PIL_LINE_COMMENT:
         case PIL_LINE_ERROR:
         case PIL_LINE_TOO_MANY:
         case PIL_LINE_TOO_FEW:
         case PIL_LINE_UNMATCHED_QUOTE:
         case PIL_LINE_NO_LF:
         case PIL_LINE_EXTRA_SPACES:
         case PIL_BAD_NAME:
         case PIL_BAD_TYPE:
         case PIL_BAD_MODE:
		PIL_append_param(fp, &pp);		/* this may fail, but we ignore errcode */
		break;
         default:
		break;					/* we silently ignore error input lines */
       }
    }
   r2 = PIL_unlock_pfile(fp);				/* file i/o done, unlock file */
   if (PIL_OK == r)  r = r2;				/* propagate error code */

   r = PIL_allow_spaces_cmndarg2value(fp);
   if (PIL_OK != r) { return(r); }

   savemode |= (fp->mode & (PIL_PARFILE_MODIFIED));	/* save modify/flush flags (flush_flag from previous call) */
   fp->mode = 0;
   if (PIL_OK == PIL_find_name(fp, "mode", &idx))	/* now get global mode for whole file */
     { if (PIL_TYPE_STRING == fp->pp[idx].type)
         { if (PIL_OK != PIL_check_mode(fp->pp[idx].strvalue, &(fp->mode)))
		fp->mode = PIL_MODE_HIDDEN;
         }
       else fp->mode = PIL_MODE_HIDDEN;			/* errors are ignored ... */
     }
   else fp->mode = PIL_MODE_HIDDEN;
   if (PIL_MODE_AUTO == fp->mode)  fp->mode = PIL_MODE_HIDDEN;  /* if it is auto convert to hidden */
   fp->mode |= savemode;				/* add modify/flush flags (from PIL_parse_line) */

   if (PIL_OK == PIL_find_name(fp, "pil_ext_indir", &idx))  /* now get indir mode for whole file */
     { fp->mode &= ~(PIL_SUBST_ENV | PIL_SUBST_INDIR | PIL_SUBST_EMPTYSTR | PIL_SUBST_CURLY); /* clear indir modes */
       if (PIL_TYPE_STRING == fp->pp[idx].type)
         { if (NULL != strchr(fp->pp[idx].strvalue, '$')) fp->mode |= PIL_SUBST_ENV;
           if (NULL != strchr(fp->pp[idx].strvalue, ')')) fp->mode |= PIL_SUBST_INDIR;
           if (NULL != strchr(fp->pp[idx].strvalue, 'E')) fp->mode |= PIL_SUBST_EMPTYSTR;
           if (NULL != strchr(fp->pp[idx].strvalue, '}')) fp->mode |= PIL_SUBST_CURLY;
         }
     }

   if (0 == (fp->mode & PIL_SUBST_EMPTYSTR))		/* load empty string translation value */
     { PILEmptyStringTemplate = getenv("PIL_EMPTY_STRING");
     }

   if (PIL_PSET_MODE & PILSpecialMode) fp->mode |= PIL_PSET_MODE;	/* specials for pset */
   if (PIL_NO_POSITIONAL & PILSpecialMode) fp->mode |= PIL_NO_POSITIONAL; /* specials for pget */

   if (NULL != fp->fp)
     { fclose(fp->fp);
       fp->fp = NULL;
     }
   return(r);
 }

/**************************************************************************************

Function:	PIL_parcmp
Description:	Compare two parameters based on their name strings.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
par1 (In)		pointer to first parameter being compared.
par2 (In)		pointer to second parameter being compared.

Return code:	-1 if first parameter's name is lexically less than second,
		+1 if first parameter's name is lexically greater, and
		 0 if they are equal. Null pointers and parameters with a
		format different from PIL_FORMAT_OK (e.g. comments)
		are always considered greater than parameters with
		format PIL_FORMAT_OK. All such non-OK parameters are
		furthermore considered equal to each other.
		Note: NOT standard PIL error code!

**************************************************************************************/

static int	PIL_parcmp(const void *par1, const void *par2)
  { const PIL_PARAM * p1 = *(const PIL_PARAM * const *) par1;
    const PIL_PARAM * p2 = *(const PIL_PARAM * const *) par2;

    if (NULL == p1 && NULL == p2) return 0;
    else if (NULL == p1) return 1;
    else if (NULL == p2) return -1;

    if (PIL_FORMAT_OK != p1->format && PIL_FORMAT_OK != p2->format) return 0;
    else if (PIL_FORMAT_OK != p1->format) return 1;
    else if (PIL_FORMAT_OK != p2->format) return -1;

    return strcmp(p1->strname, p2->strname);
  }


/**************************************************************************************

Function:	PIL_sort_parameters
Description:	Sort parameters in a PIL_PFILE structure into an
		array of pointers to PIL_PARAM. Uses qsort and PIL_parcmp.
ParameterName (I/O)	ParameterDescription:

Notes:		The input parameter array contained by the PIL_PFILE structure
		is not affected. The output array must be freed by the caller.
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file.
parray (Out)		pointer to array of pointers to PIL_PARAMs.

Return code:	standard PIL error code or PIL_OK (zero)

**************************************************************************************/

int	PIL_sort_parameters(PIL_PFILE *fp, PIL_PARAM ***parray)
  { int ii;

    if (NULL == fp) return(PIL_NUL_PTR);
    if (NULL == parray) return(PIL_NUL_PTR);
    if (NULL == fp->pp) return(PIL_OK);
    if (0 >= fp->pcnt) return(PIL_OK);

    *parray = (PIL_PARAM **) calloc(fp->pcnt, sizeof(PIL_PARAM *));
    if (NULL == *parray) return(PIL_NO_MEM);

    for (ii = 0; ii < fp->pcnt; ++ii) (*parray)[ii] = fp->pp + ii;

    qsort(*parray, fp->pcnt, sizeof(PIL_PARAM *), &PIL_parcmp);
    
    return(PIL_OK);
  }


/**************************************************************************************

Function:	PIL_pfile_open_only
Description:	open parameter file and do nothing else
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fname (In)		name of parameter file
argc (In)		number of items in parameter list (usually from main)
argv (In)		parameter list (usally from main)
access (In)		required access mode (PIL_RDWRITE or PIL_RDONLY)
numtries (In)		number of attempts to make to open the file
fp (Out)		parameter file structure to be created

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	if fname is NULL then filename is taken from argv[0]

**************************************************************************************/

static int	PIL_pfile_open_only(char *fname, int argc, char **argv, int pfile_acc_mode, int numtries, PIL_PFILE *fp)
 { static const char	*suffix = ".par";
   static const int	sleeptime = 1;
   char			*p;
   int			ii;

   if (NULL == fp) return(PIL_NUL_PTR);
   fp->argv = argv;
   fp->argc = argc;
   if (fp->argc <= 0) fp->argv = NULL;
   fp->pil_access = ((PIL_RDWRITE == pfile_acc_mode) ? PIL_RDWRITE : PIL_RDONLY);

   if (NULL == fname)					/* if no fname specified, compute one from argv[0] */
     { if (NULL == argv) return(PIL_BAD_ARG);
       p = (char *)PIL_malloc(strlen(argv[0]) + strlen(suffix) + 1);
       if (NULL == p) return(PIL_NO_MEM);
       strcpy(p, argv[0]);
       strcat(p, suffix);
     }
   /* SCREW 1145: Make a copy of the file name no matter what. */
   /* else { p = fname; } */					/* ... else you supplied one */
   else
     { p = (char *)PIL_malloc(strlen(fname) + 1);
       if (NULL == p) return(PIL_NO_MEM);
       strcpy(p, fname);
     }

   for (ii = 0; ii < numtries; ++ii)
     { fp->fp = fopen(p, ((0 == (PIL_OPEN_RDONLY & PILSpecialMode)) && (PIL_RDWRITE == pfile_acc_mode)) ? "r+" : "r");
       if (NULL != fp->fp) break;
       PIL_sleep(sleeptime);
     }

   /* SCREW 1145: Free copy of file name only if file fails to open. */
   /* if (NULL == fname) if (NULL != p) PIL_free(p); */	/* we do not need fname anymore, so do some cleanup */

   if (NULL == fp->fp)
     { if (NULL != p) { PIL_free(p); p = NULL; }
       fp->fname = NULL;
       if ((0 == (PIL_OPEN_RDONLY & PILSpecialMode)) && (PIL_RDWRITE == pfile_acc_mode))
         { return(PIL_FILE_NO_WR); }
       else
         { return(PIL_FILE_NO_RD); }
     }
   else { fp->fname = p; }

   fp->pp = NULL;					/* no parameters at the moment */
   fp->pcnt = 0;

   /* SPR 3169: Return PIL_OK, not r, which hadn't been set anywhere. */
   return(PIL_OK);
 }

/**************************************************************************************

Function:	PIL_pfile_open
Description:	open parameter file and load its parameters
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fname (In)		name of parameter file
argc (In)		number of items in parameter list (usually from main)
argv (In)		parameter list (usally from main)
access (In)		required access mode (PIL_RDWRITE or PIL_RDONLY)
numtries (In)		number of attempts to make to open the file
fp (Out)		parameter file structure to be created

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	if fname is NULL then filename is taken from argv[0]

**************************************************************************************/

int	PIL_pfile_open(char *fname, int argc, char **argv, int pfile_acc_mode, int numtries, PIL_PFILE *fp)
 { int r = PIL_pfile_open_only(fname, argc, argv, pfile_acc_mode, numtries, fp);
   if (PIL_OK != r) return(r);
   r = PIL_reload_parameters(fp);
   return(r);
 }

/**************************************************************************************

Function:	PIL_expand_cat
Description:	Concatenates source string onto destination string, preceding quotes
		and backslashes by a backslash. Used when writing strings to output
		parameter file. A maximum total  of PIL_LINESIZE characters will ever
		be written to the destination.

ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
s (In)		Source string to be appended.
d (In/Out)	Destination string to which to append s. Must hold at least PIL_LINESIZE
		characters.
index (In/Out)	Position in destination string where the next character from s should be
		written.

Return code:	None.

**************************************************************************************/
static void	PIL_expand_cat(const char *s, char *d, int *index)
 { int ii;
   for (ii = 0; s[ii] && *index < PIL_LINESIZE; ++ii, ++*index)
    { if ('"' == s[ii] || '\\' == s[ii])
	{ d[*index] = '\\';
	  ++*index;
	}
      d[*index] = s[ii];
    }
   if (*index < PIL_LINESIZE) { d[*index] = '\0'; }
 }

/**************************************************************************************

Function:	PIL_noexpand_cat
Description:	Concatenates source string onto destination string. Used when writing
		strings to output parameter file. A maximum total of PIL_LINESIZE
		characters will ever be written to the destination.

ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
s (In)		Source string to be appended.
d (In/Out)	Destination string to which to append s. Must hold at least PIL_LINESIZE
		characters.
index (In/Out)	Position in destination string where the next character from s should be
		written.

Return code:	None.

**************************************************************************************/
static void	PIL_noexpand_cat(const char *s, char *d, int *index)
 { int ii;
   for (ii = 0; s[ii] && *index < PIL_LINESIZE; ++ii, ++*index)
    { d[*index] = s[ii]; }
   if (*index < PIL_LINESIZE) { d[*index] = '\0'; }
 }


/**************************************************************************************

Function:	PIL_flush_parameters
Description:	flush parameters to disk file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	All lines are written back to parameter file, including comment lines, and in
	the order they were read from it. Parameter file is locked during write
	operation. Function flushes buffers to disk.

**************************************************************************************/

int	PIL_flush_parameters(PIL_PFILE *fp)
 { int		i, r, rewriteflag;
   PIL_PARAM	*pp;

   if (PIL_OK != (r = PIL_check_pfile(fp))) return(r);

		/* check whether anything was modified, if not do nothing */

   rewriteflag = 0;
   for (i=0; i<fp->pcnt; i++)	/* check if there are any modifications that require rewriting of the parameter file */
    { pp = &(fp->pp[i]);
      if (PIL_FORMAT_OK == fp->pp[i].format)
	if (   (PIL_PSET_MODE & fp->mode)	/* in pset mode always update all parameters */
	    || (   (PIL_MODIFIED_ANY & pp->modified)  /* or if this parameter has been modified */
	        && (   (PIL_MODE_LEARN & pp->mode)  /* and its (effective) mode has learn attrib */
	            || ((PIL_MODE_AUTO & pp->mode) && (PIL_MODE_LEARN & fp->mode))
	           )
	       )
	   )
	  { rewriteflag = 1;
	    break;
	  }
    }
   
   if ((0 == (PIL_OPEN_RDONLY & PILSpecialMode))		/* no update if applic requested RDONLY mode */
      && (rewriteflag)						/* no updates if not really necessary */
      && (PIL_RDWRITE == fp->pil_access)			/* no update if no write permission */
      && ((PIL_PARFILE_MODIFIED | PIL_PSET_MODE) & fp->mode))	/* no update if no modifications made to pfile */
     { char *tmpfile = PIL_uniq_fname(fp->fname);		/* temporary file used instead of directly writing */
       if (NULL == tmpfile) tmpfile = fp->fname;		/* use the regular file name if tmpfile invalid */
       if (NULL != (fp->fp = fopen(tmpfile, "w")))
        {
          for (i=0; i<fp->pcnt; i++)	/* and write the parameters one by one ... */
           { pp = &(fp->pp[i]);
             r = 0;
             switch (fp->pp[i].format)
	      { case PIL_FORMAT_OK:
			if (   (PIL_PSET_MODE & fp->mode)	/* in pset mode always update all parameters */
			    || (   (PIL_MODIFIED_ANY & pp->modified)  /* or if this parameter has been modified */
			        && (   (PIL_MODE_LEARN & pp->mode)  /* and its (effective) mode has learn attrib */
			            || ((PIL_MODE_AUTO & pp->mode) && (PIL_MODE_LEARN & fp->mode))
			           )
			       )
			   )
			  { /* SPR 2592: Escape " and \ with a backslash. */
			    char line[PIL_LINESIZE] = "";
			    int index = 0;

			    PIL_expand_cat(pp->strname, line, &index);
			    PIL_noexpand_cat(",", line, &index);
			    PIL_expand_cat(pp->strtype, line, &index);
			    PIL_noexpand_cat(",", line, &index);
			    PIL_expand_cat(pp->strmode, line, &index);
			    PIL_noexpand_cat(",", line, &index);

			    switch (pp->type)
			     { case PIL_TYPE_STRING:
			       case PIL_TYPE_FNAME:
					if (*pp->strvalue)
					  { PIL_noexpand_cat("\"", line, &index);
					    PIL_expand_cat(pp->strvalue, line, &index);
					    PIL_noexpand_cat("\"", line, &index);
					  }
			       		break;
			       default:
					PIL_expand_cat(pp->strvalue, line, &index);
			       		break;
			     }
			    PIL_noexpand_cat(",", line, &index);

			    PIL_expand_cat(pp->strmin, line, &index);
			    PIL_noexpand_cat(",", line, &index);

			    PIL_expand_cat(pp->strmax, line, &index);
			    PIL_noexpand_cat(",", line, &index);

			    PIL_noexpand_cat("\"", line, &index);
			    PIL_expand_cat(pp->strprompt, line, &index);
			    PIL_noexpand_cat("\"", line, &index);

			    r = fprintf(fp->fp, "%s\n", line);
			    break;
			  }
	        case PIL_FORMAT_BLANK:
	        case PIL_FORMAT_COMMENT:
			r = fprintf(fp->fp, "%s\n", pp->strline);
			break;
	      }
	     if (r < 0) { r = PIL_ERR_FWRITE; break; }
	     else { r = PIL_OK; }
	   }
       if (fflush(fp->fp))		/* flush dirty buffers to disk */
         { if (PIL_OK == r) r = PIL_ERR_FWRITE;
         }

       if (PIL_OK == r)
         { fp->mode = ((fp->mode | PIL_FLUSH_PERFORMED) & (~PIL_PARFILE_MODIFIED));
         }
       fclose(fp->fp);
       fp->fp = NULL;

       /* Move temporary unique file onto parameter file, and cleanup tmp file name. */
       if (strcmp(tmpfile, fp->fname))
         { remove(fp->fname); /* Visual Studio compiler is non-ANSI compliant; need to rm destination explicitly first. */
	       rename(tmpfile, fp->fname);
           free(tmpfile);
         }

       PIL_unlock_pfile(fp);
        }

     }

   return(r);			/* may return err code if ftruncate/write failed */
 }
 

/**************************************************************************************

Function:	PIL_pfile_close
Description:	close parameter file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
status (In)		if nonzero (status != PIL_OK) function does not call
			PIL_flush_parameters, it is all changes made so far
			are discarded

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	function clears/deallocates all fields within parameter file structure.

**************************************************************************************/

int	PIL_pfile_close(PIL_PFILE *fp, int status)
 { int i, r;

   if (PIL_OK != (r = PIL_check_pfile(fp))) return(r);

   if (PIL_OK == status)  r = PIL_flush_parameters(fp);

   for (i = fp->pcnt; i > 0; i--)  PIL_delete_param(fp, i - 1);
   if (NULL != fp->fp)
     { fclose(fp->fp);
       fp->fp = NULL;
     }
   fp->argc = 0;
   fp->argv = NULL;
   fp->pil_access = 0;
   /* SCREW 1145: Clean up file name field. */
   if (NULL != fp->fname)
     { PIL_free(fp->fname);
       fp->fname = NULL;
     }
   return(r);			/* may return err code if ftruncate/write failed */
 }

	/* reading / writing / parsing stuff */


/**************************************************************************************

Function:	PIL_read_line
Description:	read one line from opened parameter file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
buf (Out)		buffer to store line

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	buf has to have at least PIL_LINESIZE bytes. Any \n or \r characters are
	filtered out and converted to zero (end of string).
	The sleep() is to work around NFS argument caching. The flush_parameters()
	function, first truncates file to 1 (writing 0 at offset 0), and immediately
	it rewrites files contents with the new data. NFS clients sometimes sees
	the param file of length 1 for several tens of seconds.

Warning:function does not behave cleverly when it encounters line longer then
	PIL_LINESIZE. In this case it simply splits this line into pieces.
	A backslash followed by any character is replaced by a single character
	whose (unsigned) integer value is SCHAR_MAX more than the original character.
	Thus, the string returned by this function is not literally what was in the file!

**************************************************************************************/

int	PIL_read_line(PIL_PFILE *fp, char *buf, int lineno)
 { int r, i;
   unsigned char * bptr;
   char raw[PIL_LINESIZE];

   if (NULL == buf) return(PIL_NUL_PTR);
   if (PIL_OK != (r = PIL_check_pfile(fp))) return(r);
   if (NULL == fp->fp) return(PIL_BAD_ARG);

   for (i = 0; i<100; i++)
    { if (NULL == fgets(raw, PIL_LINESIZE, fp->fp)) return(PIL_ERR_FREAD);
      if (lineno > 0)  break;
      if (1L != ftell(fp->fp))  break;
      if (0 != raw[0])  break;
      fseek(fp->fp, 0L, SEEK_SET);
/* SCREW 1145: Use of sleep is encapsulated now in PIL_sleep.
      sleep(1);
*/
      PIL_sleep(1);
    }

   /* SPR 2592: Convert escape sequences: \X becomes X+128. */
   bptr = (unsigned char *) buf;
   for (i = 0; raw[i] && i < PIL_LINESIZE - 1; ++i)
    { if ('\\' == raw[i] && 0 != raw[i + 1])
        { *bptr++ = (unsigned char) raw[i + 1] + (unsigned char) SCHAR_MAX;
          ++i;
        }
      else
        { *bptr++ = raw[i];
        }
    }

   /* SPR 3580: make sure buffer is null terminated. */
   *bptr = '\0';

   return(PIL_OK);
 }


/**************************************************************************************

Function:	PIL_mark_EOS
Description:	find and mark end of string
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
*p (In/Out)		string
term (In)		end of string marker

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	function scans string for a 1st [term\n\r\0]. If found, strips off trailing
	spaces/TABS (if any) and writes end of string (zero) character

**************************************************************************************/

int	PIL_mark_eos(char **p, char term)	/* scan up to 1st [term\n\r\0], mark EOS with \0 */
 { char *p1, *p2;

   if (NULL == p) return(PIL_NUL_PTR);
   if (NULL == *p) return(PIL_BAD_ARG);
   for (p1 = *p;;(*p)++)
    { if (term == **p)
        { p2 = *p;
	  if (p1 != p2)				/* check for empty string, and if no */
	    for (p2--;;)			/* wipe out trailing spaces, if any */
	     { if ((' ' != *p2) && ('\t' != *p2)) break;
	       *p2 = 0;
	       if (p1 == p2) break;		/* was it the first character in string ? */
	       p2--;
	     }
          **p = 0;				/* normal exit, mark EOS */
	  (*p)++;				/* and position at the next character */
          return(PIL_OK);
        }
      if (('\r' == **p) || ('\n' == **p) || (0 == **p))
        { p2 = *p;
	  if (p1 != p2)				/* check for empty string, and if no */
	    for (p2--;;)			/* wipe out trailing spaces, if any */
	     { if ((' ' != *p2) && ('\t' != *p2)) break;
	       *p2 = 0;
	       if (p1 == p2) break;		/* was it the first character in the string ? */
	       p2--;
	     }
	  **p = 0;				/* we have reached EOS, so mark it */
	  return(PIL_EOS);			/* and return errcode */
        }
    }
 }


/**************************************************************************************

Function:	PIL_parse_line
Description:	parse line accoring to parameter file format
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
buf (In/Out)		raw line from parameter file
numpar (In)		position of this parameter in parameter file
pp (Out)		parameter structure (all results go there)
argc (In)		number of items in parameter list (usually from main)
argv (In)		parameter list (usually from main)

Return code:	standard PIL error code or the following code :
			PIL_OK			parameter line parsed successfully
			PIL_LINE_BLANK		blank line
			PIL_LINE_COMMENT	comment line
			PIL_LINE_ERROR		format error

Notes:	indirections are not handled here. Line (buf) is split into 7 fields and
	pointers in pp structure are set to point to those fields. Then rudimentary
	checks are performed : parameter name, type, mode are checked. If value was
	specified in command line it overrides original value in parameter file.
	Characters whose value exceeds SCHAR_MAX, which are not normally printed, will
	have SCHAR_MAX subtracted from their integer values. This is because PIL_read_line
	converts escape sequences beginning with a backslash into a single character
	which is not printable to allow PIL_parse_line to handle quotes as a delimiter.

Warning:buf is modified. fields in pp point to it. Be careful.

**************************************************************************************/

int	PIL_parse_line(PIL_PFILE *fp, char *buf, int numpar, PIL_PARAM *pp, int argc, char **argv, char *buf2)
 { char c, *p, *tmptr, *ptab[7];
   unsigned char *ubuf;
   int	r, i, extra_spaces;

   if (NULL == fp) return(PIL_NUL_PTR);		/* foolproof checks first */
   if (NULL == buf) return(PIL_NUL_PTR);	/* first check args and init pp (mostly with zeroes) */
   if (NULL == buf2) return(PIL_NUL_PTR);	/* this temporary buffer is to support reentrancy */
   if (PIL_OK != (r = PIL_init_param(pp))) return(r);

   extra_spaces = 0;
   i = 0;					/* assume no CR/LF character */
   if (NULL != (p = strchr(buf, '\n')))  { *p = 0; i = 1; }	/* filter out EOL characters */
   if (NULL != (p = strchr(buf, '\r')))  { *p = 0; i = 1; }	/* filter out EOL characters */

   strncpy(buf2, buf, PIL_LINESIZE);
   buf2[PIL_LINESIZE - 1] = 0;			/* make a copy of input line, we will need it later */

   pp->strline = buf2;

   if (PIL_LINT_MODE & PILSpecialMode)  if (0 == i)	/* if in lint mode, check if the last line has terminating LF */
     { pp->format = PIL_FORMAT_NO_LF;		/* signal line without terminating LF character */
       return(PIL_LINE_NO_LF);
     }
   
   if (0 == buf[0])
     { pp->format = PIL_FORMAT_BLANK;		/* signal blank line */
       return(PIL_LINE_BLANK);
     }
   else if ('#' == buf[0])
     { pp->format = PIL_FORMAT_COMMENT;		/* signal comment line */
       return(PIL_LINE_COMMENT);
     }

   pp->format = PIL_FORMAT_ERR;

   for (i=0; i<7; i++) ptab[i] = (char *)"";	/* discard of const attribute */
   r = PIL_OK;

   ubuf = (unsigned char *) buf; /* SPR 2592: save beginning of buffer. */

   for (i=0; i<7; i++)
    { for (;; buf++) if ((' ' != *buf) && ('\t' != *buf)) break;
      switch (*buf)
       { case '\'': ;
         case '"' : ptab[i] = ++buf;	/* found string in "xxxx" or 'xxxx' */
		    r = PIL_mark_eos(&buf, buf[-1]);
		    if (PIL_EOS == r)
		      { r = PIL_LINE_UNMATCHED_QUOTE;
		        break;
		      }
		    if (PIL_OK == r)
		      { for (;; buf++)
		         {
		           if ((' ' != *buf) && ('\t' != *buf)) break;
		           if (PIL_LINT_MODE & PILSpecialMode)
		             if (i < 6)
		               extra_spaces = 1;
		         }
			switch (*buf)
			 { case ',' : buf++; break;
			   case '\n': ;
			   case '\r': ;
			   case 0   : *buf = 0;
			              r = PIL_EOS;
			              break;
			 }
		      }
		    break;	
         case '\n': ;
         case '\r': ;
         case 0   : *buf = 0;		/* found EOS */
		    r = PIL_EOS;
		    break;
         default  : ptab[i] = buf;	/* found regular string (from ',' to ',') */
		    r = PIL_mark_eos(&buf, ',');
		    break;
       }
      if (6 == i)
        { switch (r)
           { case PIL_EOS:
			r = PIL_OK;
			break;
	     case PIL_LINE_UNMATCHED_QUOTE:
	     case PIL_LINE_EXTRA_SPACES:
			break;
	     default:
			r = PIL_LINE_TOO_MANY;
			break;
	   }
          break;
        }
      if (PIL_OK != r) break;
    }

   if (0 != extra_spaces)
     { pp->format = PIL_FORMAT_EXTRA_SPACES;
       return(r);
     }
   if (i < 6)
     { pp->format = PIL_FORMAT_TOO_FEW;
       return(PIL_LINE_TOO_FEW);
     }
   if (i > 6)
     { if (PIL_LINT_MODE & PILSpecialMode)
         { pp->format = PIL_FORMAT_TOO_MANY;
           return(PIL_LINE_TOO_MANY);
         }
     }
   if (PIL_LINE_UNMATCHED_QUOTE == r)
     { if (PIL_LINT_MODE & PILSpecialMode)
         { pp->format = PIL_FORMAT_QUOTE;
           return(r);
         }
       r = PIL_OK;
     }
   if (PIL_LINE_TOO_MANY == r)
     { if (PIL_LINT_MODE & PILSpecialMode)
         { pp->format = PIL_FORMAT_TOO_MANY;
           return(r);
         }
       r = PIL_OK;
     }
   if (PIL_OK != r)
     { pp->format = PIL_FORMAT_ERR;
       return(PIL_LINE_ERROR);
     }
		/* now the line is parsed, tokens stored, analyse further */

   /* SPR 2592: Map characters which were escaped with \ back to their
      original values. */
   for (i = 0; i < 6; ++i)
    { ubuf = (unsigned char *) ptab[i];
      while (*ubuf)
       { if (*ubuf >= (unsigned char) SCHAR_MAX)
           { *ubuf -= (unsigned char) SCHAR_MAX; }
         ++ubuf;
       }
    }
   ubuf = (unsigned char *) pp->strline;
   while (*ubuf)
    { if (*ubuf >= (unsigned char) SCHAR_MAX)
        { *ubuf -= (unsigned char) SCHAR_MAX; }
      ++ubuf;
    }

   pp->strname = ptab[0];
   pp->strtype = ptab[1];
   pp->strmode = ptab[2];
   pp->strvalue = ptab[3];
   pp->strmin = ptab[4];
   pp->strmax = ptab[5];
   pp->strprompt = ptab[6];

   		/* check parameter name */
   
   c = ptab[0][0];
   if ((!((c >= 'a') && (c <= 'z'))) && (!((c >= 'A') && (c <= 'Z'))) && (c != '$')) return(PIL_BAD_NAME);
   for (i=1;;i++)
    { c = ptab[0][i];
      if (0 == c) break;
      if (!((c >= ' ') && (c <= '~'))) return(PIL_BAD_NAME);	/* we are very liberal here .... */
    } 
   		/* check type */
   
   for (;; buf++) if ((' ' != *buf) && ('\t' != *buf)) break;	/* skip blanks */
   pp->type = PIL_TYPE_UNKNOWN;
   pp->attrib = 0;
   p = ptab[1];
   if ('*' == *p)
     { pp->attrib |= PIL_ATTRIB_REDIRECTION;
       p++;
     }
   if (!strcmp(p, "b")) pp->type = PIL_TYPE_BOOL;
   else if (!strcmp(p, "i")) pp->type = PIL_TYPE_INT4;
   else if (!strcmp(p, "r")) pp->type = PIL_TYPE_REAL8;
   else if (!strcmp(p, "s")) pp->type = PIL_TYPE_STRING;
   else if (!strcmp(p, "struct")) pp->type = PIL_TYPE_STRUCT;
   else if (!strcmp(p, "gcur")) pp->type = PIL_TYPE_GCUR;
   else if (!strcmp(p, "imcur")) pp->type = PIL_TYPE_IMCUR;
   else if (!strcmp(p, "d")) pp->type = PIL_TYPE_D;
   else if (!strcmp(p, "g")) pp->type = PIL_TYPE_G;
   if ('f' == *p)
     { for (i=1;;i++)
        { if (0 == p[i]) break;
          switch (p[i])
           { case 'r' : pp->attrib |= PIL_FILECHECK_READ; break;
             case 'w' : pp->attrib |= PIL_FILECHECK_WRITE; break;
             case 'e' : pp->attrib |= PIL_FILECHECK_EXIST; break;
             case 'n' : pp->attrib |= PIL_FILECHECK_NOEXIST; break;
	     default  : return(PIL_BAD_TYPE);
           }
          if (PIL_OK != r) break;
        }
       if (PIL_OK == r) pp->type = PIL_TYPE_FNAME;
     }

   if (PIL_TYPE_UNKNOWN == pp->type) return(PIL_BAD_TYPE);

		/* if cmd line arg was specified override value, and signal it (only on initial reading) */

   switch (PILCmndArgMode)
    { 
      case PIL_CMNDARG_ALLOW_SPACES:
/*		r = PIL_allow_spaces_cmndarg2value(fp); */
		r = PIL_NOT_FOUND;
		break;
      default:
		r = PIL_cmndarg2value(pp->strname, numpar, &tmptr, argc, argv);
		break;
    }
   if (PIL_OK == r)
     {
       pp->attrib |= PIL_VALUE_FROM_CMDLINE;
       pp->modified |= PIL_MODIFIED_VALUE;

       if (0 == (PIL_FLUSH_PERFORMED & fp->mode))
         { pp->strvalue = tmptr;
           fp->mode |= PIL_PARFILE_MODIFIED;
         }
     }   
		/* check mode - we are a bit paranoid here, and fail if bad mode */

   pp->mode = 0;
   if (PIL_OK != (r = PIL_check_mode(ptab[2], &(pp->mode)))) return(r);

   pp->format = PIL_FORMAT_OK;		/* still, further processing may find errors, but format is ok */
   return(r);
 }

	/* get / put / find functions */

/**************************************************************************************

Function:	PIL_find_name
Description:	find parameter named 'name'
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter to be found
idx (Out)		position of parameter (if found)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	functions scans list of parameters in parameter file structure. Index of
	found parameter is returned in idx. Scan is not case sensitive.

**************************************************************************************/

int	PIL_find_name(PIL_PFILE *fp, const char *name, int *idx)
 { int i, r;

   if (PIL_OK != (r = PIL_check_pfile(fp)))  return(r);
   for (i=0; i<fp->pcnt; i++)
    { if (PIL_FORMAT_OK != fp->pp[i].format) continue;	/* ignore invalid lines */
      if (NULL == fp->pp[i].strname) continue;		/* ignore invalid lines */
      if (!strcasecmp(name, fp->pp[i].strname))
        { *idx = i;					/* line returned may have errors */
          return(PIL_OK);
        }
    }
   return(PIL_NOT_FOUND);
 }

/**************************************************************************************

Function:	PIL_get_range
Description:	get allowable range of parameter value (if any)
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
idx (In)		index of parameter
vmin (Out)		minimum value
vmax (Out)		maximum value

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	functions returns min and max value limits for given parameter. Both min and
	max fields must be present in parameter file. Function resolves any
	indirections found in min and max fields. Min value must be less than max
	value. Limits can be specified for: PIL_TYPE_INT4, PIL_TYPE_REAL8,
	PIL_TYPE_STRING and PIL_TYPE_FNAME parameter types. If the above conditions 
	are not met functions returns with error PIL_BAD_ARG and PIL library assumes 
	that no range limit is specified for given parameter.

**************************************************************************************/

int	PIL_get_range(PIL_PFILE *fp, int idx, PIL_VALUE *vmin, PIL_VALUE *vmax)
 { char *smin, *smax;
   PIL_PARAM *pp;
   int r;

   if (PIL_OK != (r = PIL_check_par(fp, idx))) return(r);
   if ((NULL == vmin) || (NULL == vmax)) return(PIL_NUL_PTR);
   pp = &(fp->pp[idx]);

   smin = NULL;
   if (PIL_OK != (r = PIL_indir2string(fp, PIL_WHAT_MIN, &smin, pp->strmin))) return(r); /* handle min indirections, if any */
   r = PIL_string2value(smin, pp->type, vmin);
   PIL_free(smin);
   smin = NULL;
   if (PIL_OK != r) return(r);

   smax = NULL;
   if (PIL_OK != (r = PIL_indir2string(fp, PIL_WHAT_MAX, &smax, pp->strmax))) return(r); /* handle max indirections, if any */
   r = PIL_string2value(smax, pp->type, vmax);
   PIL_free(smax);
   smax = NULL;
   if (PIL_OK != r) return(r);

   switch (pp->type)
    { case PIL_TYPE_INT4:
		if (vmax->i < vmin->i)  r = PIL_BAD_ARG;
		break;
      case PIL_TYPE_REAL8:
		if (vmax->d < vmin->d)  r = PIL_BAD_ARG;
		break;
      case PIL_TYPE_STRING:
      case PIL_TYPE_FNAME:
		if ((0 == vmin->s[0]) || (0 == vmax->s[0])) { r = PIL_BAD_ARG; break; }
		if (strcmp(vmax->s, vmin->s) < 0)  r = PIL_BAD_ARG;
		break;
      default:	break;
    }
   return(r);
 }

/**************************************************************************************

Function:	PIL_get_enum_range
Description:	get allowable range or enumerated list of values for parameter (if any)
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
idx (In)		index of parameter
vmin (Out)		minimum value
vmax (Out)		maximum value
enumlist (Out)		list of enumerated values (filled only in enum mode)
unumcount (Out)		number of items in enumlist (-1 if not enum mode)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	functions returns min and max value limits for given parameter. Both min and
	max fields must be present in parameter file. Function resolves any
	indirections found in min and max fields. Min value must be less than max
	value. Limits can be specified for: PIL_TYPE_INT4, PIL_TYPE_REAL8,
	PIL_TYPE_STRING and PIL_TYPE_FNAME parameter types. If above conditions are not
	met function returns with error PIL_BAD_ARG and PIL library assumes that no
	range limit is specified for given parameter.

        enum mode is effective when : max field is empty, min field contains
        '|' separated list of items (strings/ints/reals/filenames)

	this function is more general then PIL_get_range() and should be called instead.
**************************************************************************************/

int	PIL_get_enum_range(PIL_PFILE *fp, int idx, PIL_VALUE *vmin, PIL_VALUE *vmax, PIL_VALUE **enumlist, int *enumcount)
 { char		*p, *smin, *smax, *stmp;
   PIL_VALUE	tmpvmin, tmpvmax, tmpenum, *vp, *nvp;
   PIL_PARAM	*pp;
   int		j, r;

   if (NULL == enumlist) return(PIL_NUL_PTR);
   if (NULL == enumcount) return(PIL_NUL_PTR);
   if (PIL_OK != (r = PIL_check_par(fp, idx))) return(r);
   if ((NULL == vmin) || (NULL == vmax)) return(PIL_NUL_PTR);
   pp = &(fp->pp[idx]);

   if (PIL_OK == (r = PIL_get_range(fp, idx, vmin, vmax)))
     { *enumlist = NULL;
       *enumcount = -1;
       return(r);
     }

   smin = NULL;
   if (PIL_OK != (r = PIL_indir2string(fp, PIL_WHAT_MIN, &smin, pp->strmin))) return(r); /* handle min indirections, if any */
   r = PIL_string2value(smin, PIL_TYPE_STRING, &tmpvmin);
   PIL_free(smin);
   smin = NULL;
   if (PIL_OK != r) return(r);

   if (NULL == strchr(tmpvmin.s, '|'))  return(PIL_BAD_ARG);	/* at least 1 '|' is needed for enum list */
		/* the line above means that "1.23|" specifies valid enumlist with 1 element only !!! */

   smax = NULL;
   if (PIL_OK != (r = PIL_indir2string(fp, PIL_WHAT_MAX, &smax, pp->strmax))) return(r); /* handle max indirections, if any */
   r = PIL_string2value(smax, PIL_TYPE_STRING, &tmpvmax);
   PIL_free(smax);
   smax = NULL;
   if (PIL_OK != r) return(r);

   if (0 != tmpvmax.s[0])  return(PIL_BAD_ARG);			/* max field must be empty in enum mode */

   j = 0;
   nvp = vp = NULL;						/* initialize enum value vector (empty) */

   for (p = strtok(&(tmpvmin.s[0]), "|"); NULL != p; p = strtok(NULL, "|"))	/* break string into pieces */
    { stmp = NULL;
      if (PIL_OK != (r = PIL_indir2string(fp, PIL_WHAT_MIN, &stmp, p))) break;	/* handle min indirections, if any */
      r = PIL_string2value(stmp, pp->type, &tmpenum);		/* convert to typed value */
      PIL_free(stmp);
      stmp = NULL;
      if (PIL_OK != r) break;

      if (NULL == vp)  { nvp = (PIL_VALUE *)PIL_malloc(sizeof(PIL_VALUE)); }	/* resize buffer */
      else  { nvp = (PIL_VALUE *)PIL_realloc((void *)vp, (j + 1) * sizeof(PIL_VALUE)); }
      if (NULL == nvp)
        { r = PIL_NO_MEM;
          break;
        }
      vp = nvp;							/* update pointers and store new value */
      vp[j] = tmpenum;
      j++;
    }

   if (PIL_OK == r)
     { if (0 == j)  r = PIL_BAD_ARG;				/* signal both min/max empty */
     }

   if (PIL_OK != r)						/* clean up in case of error */
     {
       if (NULL != vp)
         { PIL_free((void *)vp);				/* in error condition, we do not leave any bufs allocated */
           vp = NULL;
         }
       j = -1;							/* signal no enum mode also means enums are ignored on error */
     }

   *enumlist = vp;						/* store results - WARNING: enumlist must be dealloc later */
   *enumcount = j;

   return(r);
 }


/**************************************************************************************

Function:	PIL_put_value
Description:	set new value of parameter in parameter table
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
idx (In)		index of parameter in question
type (In)		data type of new value
vs (In)			typed new value of parameter
strv (In)		ascii representation of new value
mode (In)		mode of operation

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	function first checks whether 'type' equals type of parameter. If this is not
	the case function returns error. Then if mode has PIL_CHECK_RANGE_YES bit
	set on function performs range checking. For data type PIL_TYPE_FNAME
	file access mode is checked. If range check fails function returns with
	PIL_BAD_ARG error code.
	Finally function writes typed value and ascii string value to parameter
	table and sets PIL_MODIFIED_VALUE and PIL_PARFILE_MODIFIED flags in
	the parameter file structure.
	This function is a low-level function, called internally by PIL_put_by_value

Warning:Parameter file is not updated until PILFlushParameters is called. All changes
	done by PIL_get_by_value are in memory (parameter table in parameter file
	structure)

**************************************************************************************/

int	PIL_put_value(PIL_PFILE *fp, int idx, int type, PIL_VALUE *vs, const char *strv, int mode)
 { char		*p, rootfname[PIL_LINESIZE];
   int		i, j, r, r2, venumcount;
   PIL_VALUE	vmin, vmax, *venum;
   PIL_PARAM	*pp;

   if (NULL == vs) return(PIL_NUL_PTR);
   if (NULL == strv) return(PIL_NUL_PTR);
   if (PIL_OK != (r = PIL_check_par(fp, idx))) return(r);
   pp = &(fp->pp[idx]);
   if (type != pp->type) return(PIL_BAD_ARG);
   
   if (mode & PIL_CHECK_RANGE_YES) switch (pp->type)
    { case PIL_TYPE_BOOL:
		break;
      case PIL_TYPE_STRING:
		if (NULL == vs->s) { r = PIL_BAD_ARG; break; }
		if (PIL_OK == PIL_get_enum_range(fp, idx, &vmin, &vmax, &venum, &venumcount))
		  { if (-1 == venumcount)
		      { if ((strcmp(vs->s, vmin.s) < 0) || (strcmp(vs->s, vmax.s) > 0))  r = PIL_OFF_RANGE;
		      }
		    else
		      { r = PIL_BAD_ENUM_VALUE;
		        for (i=0; i<venumcount; i++)  if (0 == strcasecmp(venum[i].s, vs->s))
		           { for (j=0; j<PIL_LINESIZE; j++)	/* convert string to uppercase */
		              { if (0 == vs->s[j]) break;
		                vs->s[j] = toupper(vs->s[j] & 0xFF);
		              }
		             r = PIL_OK; 
		             break;
		           }
		        if (NULL != venum)  { PIL_free((void *)venum); venum = NULL; }
		      }
		  }
		break;
      case PIL_TYPE_INT4:
		if (PIL_OK == PIL_get_enum_range(fp, idx, &vmin, &vmax, &venum, &venumcount))
		  { if (-1 == venumcount)
		      { if ((vs->i < vmin.i) || (vs->i > vmax.i)) r = PIL_OFF_RANGE;
		      }
		    else
		      { r = PIL_BAD_ENUM_VALUE;
		        for (i=0; i<venumcount; i++)  if (venum[i].i == vs->i)  { r = PIL_OK; break; }
		        if (NULL != venum)  { PIL_free((void *)venum); venum = NULL; }
		      }
		  }
		break;
      case PIL_TYPE_REAL8:
		if (PIL_OK == PIL_get_enum_range(fp, idx, &vmin, &vmax, &venum, &venumcount))
		  { if (-1 == venumcount)
		      { if ((vs->d < vmin.d) || (vs->d > vmax.d)) r = PIL_OFF_RANGE;
		      }
		    else
		      { r = PIL_BAD_ENUM_VALUE;
		        for (i=0; i<venumcount; i++)  if (venum[i].d == vs->d)  { r = PIL_OK; break; }
		        if (NULL != venum)  { PIL_free((void *)venum); venum = NULL; }
		      }
		  }
		break;
      case PIL_TYPE_FNAME:
		if (NULL == vs->s) { r = PIL_BAD_ARG; break; }
		if (PIL_OK == PIL_get_enum_range(fp, idx, &vmin, &vmax, &venum, &venumcount))
		  { if (-1 == venumcount)
		      { if ((strcmp(vs->s, vmin.s) < 0) || (strcmp(vs->s, vmax.s) > 0))  r = PIL_OFF_RANGE;
		      }
		    else
		      { r = PIL_BAD_ENUM_VALUE;
		        for (i=0; i<venumcount; i++)  if (0 == strcmp(venum[i].s, vs->s))  { r = PIL_OK; break; }
		        if (NULL != venum)  { PIL_free((void *)venum); venum = NULL; }
		      }
		  }
		if (PIL_OK != r) break;

		if (PIL_PSET_MODE & PILSpecialMode) break; /* In pset mode, skip all checks for file existence. */

		strncpy(rootfname, vs->s, PIL_LINESIZE);	/* convert URL to file name */
		rootfname[PIL_LINESIZE - 1] = 0;
		r2 = PIL_root_name(rootfname);	/* remove file://, [2], and the stuff like that */
		switch (r2)
		 { case PIL_ROOTNAME_FILE:
/* SCREW 1145: Use of access is encapsulated now in PIL_file_* functions.
			if (pp->attrib & PIL_FILECHECK_EXIST)  if (access(rootfname, F_OK)) { r = PIL_BAD_FILE_ACCESS; break; }
			if (pp->attrib & PIL_FILECHECK_NOEXIST)  if (0 == access(rootfname, F_OK)) { r = PIL_BAD_FILE_ACCESS; break; }
			if (pp->attrib & PIL_FILECHECK_READ)  if (access(rootfname, R_OK)) { r = PIL_BAD_FILE_ACCESS; break; }
			if (pp->attrib & PIL_FILECHECK_WRITE)  if (access(rootfname, W_OK)) { r = PIL_BAD_FILE_ACCESS; break; }
*/
			if (pp->attrib & PIL_FILECHECK_EXIST)  if (!PIL_file_exists(rootfname)) { r = PIL_BAD_FILE_ACCESS; break; }
			if (pp->attrib & PIL_FILECHECK_NOEXIST)  if (PIL_file_exists(rootfname)) { r = PIL_BAD_FILE_ACCESS; break; }
			if (pp->attrib & PIL_FILECHECK_READ)  if (!PIL_file_readable(rootfname)) { r = PIL_BAD_FILE_ACCESS; break; }
			if (pp->attrib & PIL_FILECHECK_WRITE)  if (!PIL_file_writable(rootfname)) { r = PIL_BAD_FILE_ACCESS; break; }
			break;
		   case PIL_ROOTNAME_STDIN:
			if (pp->attrib & PIL_FILECHECK_NOEXIST)  { r = PIL_BAD_FILE_ACCESS; break; }
			if (pp->attrib & PIL_FILECHECK_WRITE)  { r = PIL_BAD_FILE_ACCESS; break; }
			break;
		   case PIL_ROOTNAME_STDOUT:
			if (pp->attrib & PIL_FILECHECK_NOEXIST)  { r = PIL_BAD_FILE_ACCESS; break; }
			if (pp->attrib & PIL_FILECHECK_READ)  { r = PIL_BAD_FILE_ACCESS; break; }
			break;
		   case PIL_ROOTNAME_STDINOUT:
			if (pp->attrib & PIL_FILECHECK_NOEXIST)  { r = PIL_BAD_FILE_ACCESS; break; }
			break;
		   default:
			r = PIL_BAD_FILE_ACCESS;		/* http://, ftp:// shmem:// and similar cases ... */
			break;
		 }
		break;
      default:	r = PIL_NOT_IMPLEMENTED;	/* we silently ignore unrecognized types ... */
    }
   if (PIL_OK == r)
     { if (PIL_OK == (r = PIL_dup(&p, strv)))
	 { if (NULL != pp->strvalue) PIL_free(pp->strvalue);	/* free old value string if any */
	   pp->strvalue = p;			/* store new one */
	   if (PIL_MODIFIED_YES & mode)
	     { pp->modified |= PIL_MODIFIED_VALUE;
	       fp->mode |= PIL_PARFILE_MODIFIED;
	     }
	 }
     }
   return(r);
 }


/**************************************************************************************

Function:	PIL_put_by_value
Description:	set new value of named parameter in parameter table (typed value)
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
type (In)		data type of new value
p (In)			typed new value of parameter (pointer to it)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	function first calls PIL_find_name to get index of parameter, then it checks
	status of this parameter. If everything seems to be OK, functions fills
	parameter value structure with data pointed to by 'p', then calculates its
	ascii string representation (PIL_value2string) and finally calls
	PIL_put_value
	This function is an intermediate level function.
	See companion function PIL_put_by_string

Warning:Parameter file is not updated until PILFlushParameters is called. All changes
	done by PIL_get_by_value are in memory (parameter table in parameter file
	structure)

**************************************************************************************/

int     PIL_put_by_value(PIL_PFILE *fp, const char *name, int type, void *p)
 { char buf[PIL_LINESIZE];
   PIL_PARAM *pp;
   PIL_VALUE v;
   int r, idx;
       
   if (NULL == name) return(PIL_NUL_PTR);
   if (0 == name[0]) return(PIL_BAD_ARG);
   if (PIL_OK != (r = PIL_find_name(fp, name, &idx))) return(r);
   if (PIL_OK != (r = PIL_check_par(fp, idx))) return(r);
   if (NULL == p) return(PIL_BAD_ARG);
   pp = &(fp->pp[idx]);
   if (type != pp->type) return(PIL_BAD_ARG);
   switch (type)
    { case PIL_TYPE_BOOL:
		v.b = *(int *)p;
		break;
      case PIL_TYPE_INT4:
		v.i = *(int *)p;
		break;
      case PIL_TYPE_REAL8:
		v.d = *(double *)p;
		break;
      case PIL_TYPE_STRING:
      case PIL_TYPE_FNAME:
		strncpy(v.s, (char *)p, PIL_LINESIZE);
		v.s[PIL_LINESIZE - 1] = 0;
		break;
      default:	return(PIL_BAD_ARG);
    }
   if (PIL_OK != (r = PIL_value2string(&v, pp->type, buf, PIL_LINESIZE))) return(r);
   return(PIL_put_value(fp, idx, pp->type, &v, buf, PIL_MODIFIED_YES | PIL_CHECK_RANGE_YES));
 }


/**************************************************************************************

Function:	PIL_put_by_string
Description:	set new value of named parameter in parameter table (string value)
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
src (In)		new value of parameter (ascii string)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	function first calls PIL_find_name to get index of parameter, then it checks
	status of this parameter. If everything seems to be OK, functions converts
	ascii value to typed value (PIL_string2value) and finally calls
	PIL_put_value.
	This function is an intermediate level function.
	See companion function PIL_put_by_value

Warning:Parameter file is not updated until PILFlushParameters is called. All changes
	done by PIL_get_by_value are in memory (parameter table in parameter file
	structure)

**************************************************************************************/

int     PIL_put_by_string(PIL_PFILE *fp, const char *name, const char *src)
 { PIL_PARAM *pp;
   PIL_VALUE v;
   int r, idx;
       
   if (PIL_OK != (r = PIL_find_name(fp, name, &idx))) return(r);
   if (PIL_OK != (r = PIL_check_par(fp, idx))) return(r);
   if (NULL == src) return(PIL_BAD_ARG);
   pp = &(fp->pp[idx]);
   if (PIL_OK != (r = PIL_string2value(src, pp->type, &v))) return(r);
   return(PIL_put_value(fp, idx, pp->type, &v, src, PIL_MODIFIED_YES | PIL_CHECK_RANGE_YES));
 }


/**************************************************************************************

Function:	PIL_put_bool
Description:	set new value of named boolean parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
b (In)			new value of parameter (integer)

Return code:	standard PIL error code or PIL_OK (zero)

**************************************************************************************/

int     PIL_put_bool(PIL_PFILE *fp, const char *name, int b)
 { return(PIL_put_by_value(fp, name, PIL_TYPE_BOOL, (void *)(&b)));  }


/**************************************************************************************

Function:	PIL_put_int
Description:	set new value of named integer parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
i (In)			new value of parameter (integer)

Return code:	standard PIL error code or PIL_OK (zero)

**************************************************************************************/

int     PIL_put_int(PIL_PFILE *fp, const char *name, int i)
 { return(PIL_put_by_value(fp, name, PIL_TYPE_INT4, (void *)(&i)));  }


/**************************************************************************************

Function:	PIL_put_real
Description:	set new value of named double (real*8) parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
d (In)			new value of parameter (double)

Return code:	standard PIL error code or PIL_OK (zero)

**************************************************************************************/

int     PIL_put_real(PIL_PFILE *fp, const char *name, double d)
 { return(PIL_put_by_value(fp, name, PIL_TYPE_REAL8, (void *)(&d)));  }


/**************************************************************************************

Function:	PIL_put_string
Description:	set new value of named string parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
s (In)			new value of parameter (ascii string)

Return code:	standard PIL error code or PIL_OK (zero)

**************************************************************************************/

int     PIL_put_string(PIL_PFILE *fp, const char *name, const char *s)
 { return(PIL_put_by_value(fp, name, PIL_TYPE_STRING, (void *)s));  }


/**************************************************************************************

Function:	PIL_put_fname
Description:	set new value of named string parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
s (In)			new value of parameter (ascii string)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	file access mode is verified by PIL_put_by_value

**************************************************************************************/

int     PIL_put_fname(PIL_PFILE *fp, const char *name, const char *s)
 { return(PIL_put_by_value(fp, name, PIL_TYPE_FNAME, (void *)s));  }


/**************************************************************************************

Function:	PIL_get_by_value
Description:	get (ask) for a value of parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
v (Out)			value of parameter
type (In)		data type

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	function first calls PIL_find_name to get index of parameter, then it checks
        status of this parameter. If everything seems to be OK, functions resolves
	indirections (if any). Now, when Query mode is off, functions either
	immediately returns parameter value, or if query mode is on or value of
	parameter id out of range prompts user to enter correct value. It keeps
	prompting user for a new value until correct value is entered.
	Once correct value is entered it is converted to proper type, then new value
	is written to parameter file structure (in memory) and function returns
	Function updates PIL_MODIFIED_XXX flags when necessary.
	When PIL_QUERY_OVERRIDE mode is in effect stdin/stdout i/o is disabled and
	function returns error when current parameter value is invalid or out of
	range.

        This function is an intermediate level function. It is called by
        family of PIL_get_xxx routines.

Warning:Parameter file is not updated until PILFlushParameters is called. All changes
	done by PIL_get_by_value are in memory (parameter table in parameter file
	structure)

**************************************************************************************/

#ifdef HAVE_READLINE_LIBS
static	char *rl_def_value;

static	int	pil_set_rl_def_value(void)
 { return(rl_insert_text(rl_def_value));
 }
#endif


int	PIL_get_by_value(PIL_PFILE *fp, const char *name, PIL_VALUE *v, int type)
 { char		buf[PIL_LINESIZE], *p, *eval, pr2[PIL_LINESIZE];
#ifdef HAVE_READLINE_LIBS
   char		prompt[PIL_LINESIZE], *line_read;
#endif
   int		emode, etype, changed, ft, r, idx, venumcount;
   PIL_VALUE	vmin, vmax, *venum;
   PIL_PARAM	*pp;

   if (NULL == v) return(PIL_NUL_PTR);
   if (NULL == name) return(PIL_NUL_PTR);
   if (0 == name[0]) return(PIL_BAD_ARG);
   if (PIL_OK != (r = PIL_find_name(fp, name, &idx)))
     { if (PIL_NOT_FOUND == r)
         { strncpy(PILLastNotFoundName, name, PIL_LINESIZE);
           PILLastNotFoundName[PIL_LINESIZE - 1] = 0;
         }
       return(r);
     }
   if (PIL_OK != (r = PIL_check_par(fp, idx))) return(r);
   pp = &(fp->pp[idx]);

   etype = type;
   if (PIL_TYPE_DOL == type)  type = PIL_TYPE_STRING;
   if (type != pp->type) return(PIL_BAD_ARG);

   emode = pp->mode;
   if (PIL_MODE_AUTO == emode) emode = fp->mode;
   if (PIL_MODE_AUTO == emode) emode = PIL_MODE_HIDDEN;

   /* for pset we want every parameter to be queried */
   if (PIL_PSET_MODE & PILSpecialMode)  emode = PIL_MODE_QUERY;

   eval = NULL;

   for (ft = 1;; ft=0)
    {
      if (NULL != eval)  { PIL_free(eval); eval = NULL; }

      changed = 0;
      if (PIL_QUERY_OVERRIDE != PILQueryMode)
	/* SCREW 1496: Check reprompt flag to see whether to force reprompting even for hidden parameters, and
	   for parameters supplied on the command line. */
        if ((0 == ft) || (0 != pp->reprompt) || ((PIL_MODE_QUERY & emode) && (0 == (PIL_VALUE_FROM_CMDLINE & pp->attrib))))
	  { pp->reprompt = 0;
            sprintf(pr2, "%s ", (pp->strprompt[0] ? pp->strprompt : pp->strname) );
	    if (PIL_OK == PIL_get_enum_range(fp, idx, &vmin, &vmax, &venum, &venumcount))
	      { if (-1 == venumcount)
		  { sprintf(pr2 + strlen(pr2), "<%s - %s> ", pp->strmin, pp->strmax);  /* print min/max range */
		  }
		else
		  { sprintf(pr2 + strlen(pr2), "<%s> ", pp->strmin);		/* print enum's list */
		    if (NULL != venum)  { PIL_free((void *)venum); venum = NULL; }
		  }
	      }
#ifndef HAVE_READLINE_LIBS					/* simple printf/fgets case */
            sprintf(pr2 + strlen(pr2), "[%s] : ", pp->strvalue);
            printf(pr2);					/* print prompt string */
            if (NULL == fgets(buf, PIL_LINESIZE, stdin))	/* and read user input */
              { r = PIL_ERR_FREAD; break; }
	    if (NULL != (p = strchr(buf, '\n'))) *p = 0;	/* filter out EOL characters */
	    if (NULL != (p = strchr(buf, '\r'))) *p = 0;	/* filter out EOL characters */
	    if (buf[0]) changed = 1;
	    if (0 == (fp->mode & PIL_SUBST_EMPTYSTR))		/* if enabled, translate empty string */
	      { 
	        if (0 == strcmp(buf, (PILEmptyStringTemplate ? PILEmptyStringTemplate : "\"\""))) buf[0] = 0;
	      }
#else								/* readline library case */
	    if (PIL_RL_PROMPT_IEB == PILReadlinePromptMode)
	      {
                sprintf(prompt, "%s: ", pr2);			/* construct prompt string for readlinelib */
	        rl_def_value = pp->strvalue;
	        rl_startup_hook = pil_set_rl_def_value;		/* force readline to have nonempty buffer value */
	        line_read = readline(prompt);			/* get a line from the user. */
	        rl_def_value = NULL;				/* not really necessary, but just in case */
	        if (NULL == line_read)  { r = PIL_ERR_FREAD; break; }

	        strncpy(buf, line_read, PIL_LINESIZE - 1);	/* copy to output variable */
	        buf[PIL_LINESIZE - 1] = 0;

	        if (0 != line_read[0])  add_history(line_read);	/* add nonempty line to the history list */ 
	        free(line_read);
	        line_read = (char *)NULL;

	        if (NULL != (p = strchr(buf, '\n'))) *p = 0;	/* filter out EOL characters */
	        if (NULL != (p = strchr(buf, '\r'))) *p = 0;	/* filter out EOL characters */
	        changed = 1;
	      }
	    else
	      {
                sprintf(prompt, "%s[%s] : ", pr2, pp->strvalue); /* construct prompt string for readlinelib */
	        line_read = readline(prompt);			/* get next line from the user. */
	        if (NULL == line_read)  { r = PIL_ERR_FREAD; break; }

	        strncpy(buf, line_read, PIL_LINESIZE - 1);	/* copy new line to the output variable */
	        buf[PIL_LINESIZE - 1] = 0;

	        if (0 != line_read[0])  add_history(line_read);	/* add nonempty line to the history list */ 
	        free(line_read);				/* free readline allocated resources */
	        line_read = (char *)NULL;

	        if (NULL != (p = strchr(buf, '\n'))) *p = 0;	/* filter out EOL characters */
	        if (NULL != (p = strchr(buf, '\r'))) *p = 0;	/* filter out EOL characters */
	        if (buf[0]) changed = 1;
	        if (0 == (fp->mode & PIL_SUBST_EMPTYSTR))	/* if enabled, translate empty string */
	          { 
	            if (0 == strcmp(buf, (PILEmptyStringTemplate ? PILEmptyStringTemplate : "\"\""))) buf[0] = 0;
	          }
	      }
#endif
          }

      if (0 == changed)
        { strncpy(buf, pp->strvalue, PIL_LINESIZE);		/* so the buf contains always value to be written (nonexpanded) */
	  buf[PIL_LINESIZE - 1] = 0;
	}

      if (PIL_OK != (r = PIL_indir2string(fp, PIL_WHAT_VALUE, &eval, buf))) break;  /* indirections */
      if (PIL_OK == (r = PIL_string2value(eval, etype, v)))	/* JB 03-Jul-01: was pp->type, but now there is PIL_TYPE_DOL */
        { 
          if (PIL_OK == (r = PIL_put_value(fp, idx, pp->type, v, buf,
               (changed ? (PIL_MODIFIED_YES | PIL_CHECK_RANGE_YES) : PIL_CHECK_RANGE_YES))))
            break;	/* value assigned properly, but may loop infinitely here .... */
        }
      if (PIL_QUERY_OVERRIDE == PILQueryMode) break;
      PIL_log_error("Invalid parameter value", r);
    }      

   if (NULL != eval)  { PIL_free(eval); eval = NULL; }
   return(r);
 }


/**************************************************************************************

Function:	PIL_get_bool
Description:	get (ask) for a value of boolean parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'

**************************************************************************************/

int	PIL_get_bool(PIL_PFILE *fp, const char *name, int *result)
 { char		buf[PIL_LINESIZE + 1000];
   PIL_VALUE	v;
   int		r;

   if (NULL == result) return(PIL_NUL_PTR);
   if (PIL_OK == (r = PIL_get_by_value(fp, name, &v, PIL_TYPE_BOOL)))  *result = v.b;
   sprintf(buf, "PIL_get_bool(%s)", name);
   PIL_log_error(buf, r);
   return(r);
 }


/**************************************************************************************

Function:	PIL_get_int
Description:	get (ask) for a value of integer parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'

**************************************************************************************/

int	PIL_get_int(PIL_PFILE *fp, const char *name, int *result)
 { char		buf[PIL_LINESIZE + 1000];
   PIL_VALUE	v;
   int		r;

   if (NULL == result) return(PIL_NUL_PTR);
   if (PIL_OK == (r = PIL_get_by_value(fp, name, &v, PIL_TYPE_INT4)))  *result = v.i;
   sprintf(buf, "PIL_get_int(%s)", name);
   PIL_log_error(buf, r);
   return(r);
 }


/**************************************************************************************

Function:	PIL_get_real
Description:	get (ask) for a value of double (real*8) parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'

**************************************************************************************/

int	PIL_get_real(PIL_PFILE *fp, const char *name, double *result)
 { char		buf[PIL_LINESIZE + 1000];
   PIL_VALUE	v;
   int		r;

   if (NULL == result) return(PIL_NUL_PTR);
   if (PIL_OK == (r = PIL_get_by_value(fp, name, &v, PIL_TYPE_REAL8)))  *result = v.d;
   sprintf(buf, "PIL_get_real(%s)", name);
   PIL_log_error(buf, r);
   return(r);
 }


/**************************************************************************************

Function:	PIL_get_string
Description:	get (ask) for a value of string parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	result has to have at least PIL_LINESIZE bytes long.

**************************************************************************************/

int	PIL_get_string(PIL_PFILE *fp, const char *name, char *result)
 { char		buf[PIL_LINESIZE + 1000];
   PIL_VALUE	v;
   int		r;

   if (NULL == result) return(PIL_NUL_PTR);
   if (PIL_OK == (r = PIL_get_by_value(fp, name, &v, PIL_TYPE_STRING)))
     { strcpy(result, v.s); }
   sprintf(buf, "PIL_get_string(%s)", name);
   PIL_log_error(buf, r);
   return(r);
 }


/**************************************************************************************

Function:	PIL_get_fname
Description:	get (ask) for a value of file name parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	result has to have at least PIL_LINESIZE bytes long.

**************************************************************************************/

int	PIL_get_fname(PIL_PFILE *fp, const char *name, char *result)
 { char		buf[PIL_LINESIZE + 1000];
   PIL_VALUE	v;
   int		r;

   if (NULL == result) return(PIL_NUL_PTR);
   if (PIL_OK == (r = PIL_get_by_value(fp, name, &v, PIL_TYPE_FNAME)))
     { strcpy(result, v.s); }
   sprintf(buf, "PIL_get_fname(%s)", name);
   PIL_log_error(buf, r);
   return(r);
 }


/**************************************************************************************

Function:	PIL_get_dol
Description:	get (ask) for a value of DOL parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	result has to have at least PIL_LINESIZE bytes long.

**************************************************************************************/

int	PIL_get_dol(PIL_PFILE *fp, const char *name, char *result)
 { char		buf[PIL_LINESIZE + 1000];
   PIL_VALUE	v;
   int		r;

   if (NULL == result) return(PIL_NUL_PTR);
   if (PIL_OK == (r = PIL_get_by_value(fp, name, &v, PIL_TYPE_DOL)))
     { strcpy(result, v.s); }
   sprintf(buf, "PIL_get_dol(%s)", name);
   PIL_log_error(buf, r);
   return(r);
 }

/**************************************************************************************

Function:	PIL_get_as_string
Description:	get (ask) for a value of any parameter type, return it as a string.
		This function was added to address SCREW 1497.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
fp (In)			pointer to parameter file structure
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	result has to have at least PIL_LINESIZE bytes long.

**************************************************************************************/

int	PIL_get_as_string(PIL_PFILE *fp, const char *name, char *result)
 { char		buf[PIL_LINESIZE + 1000];
   PIL_VALUE	v;
   int		r = PIL_OK;
   int		idx = 0;

   if (NULL == result) return(PIL_NUL_PTR);
   /* Determine the index of the parameter. */
   if (PIL_OK != (r = PIL_find_name(fp, name, &idx)))
     { if (PIL_NOT_FOUND == r)
	{ strncpy(PILLastNotFoundName, name, PIL_LINESIZE);
	  PILLastNotFoundName[PIL_LINESIZE - 1] = 0;
	}
     }
   else
     /* Use PIL_get_by_value to prompt for "q" parameters, using the type which is
        known to be correct. The returned value v is ignored, because what is needed
        is the string representation of the parameter, regardless of its type. */
     { int par_type = fp->pp[idx].type;
       fp->pp[idx].type = PIL_TYPE_STRING;
       if (PIL_OK == (r = PIL_get_by_value(fp, name, &v, fp->pp[idx].type)))
        /* The parameter struct now has the correct value, so copy it to the output. */
	{ strncpy(result, fp->pp[idx].strvalue, PIL_LINESIZE);
	  result[PIL_LINESIZE-1] = '\0';
        }
       fp->pp[idx].type = par_type;
     }

   sprintf(buf, "PIL_get_as_string(%s)", name);
   PIL_log_error(buf, r);
   return(r);
 }

		/* O F F I C I A L   A P I   R O U T I N E S */


/**************************************************************************************

Function:	PILGetBool
Description:	get (ask) for a value of boolean parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	This is an official API routine.

**************************************************************************************/

int	PILGetBool(const char *name, int *result)
 { return(PIL_get_bool(&PILDefaultPFile, name, result)); }


/**************************************************************************************

Function:	PILGetInt
Description:	get (ask) for a value of integer parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	This is an official API routine.

**************************************************************************************/

int	PILGetInt(const char *name, int *result)
 { return(PIL_get_int(&PILDefaultPFile, name, result)); }


/**************************************************************************************

Function:	PILGetReal
Description:	get (ask) for a value of double (real*8) parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	This is an official API routine.
	See also companion routine PILGetReal4

**************************************************************************************/

int	PILGetReal(const char *name, double *result)
 { return(PIL_get_real(&PILDefaultPFile, name, result)); }


/**************************************************************************************

Function:	PILGetReal4
Description:	get (ask) for a value of float (real*4) parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	This is an official API routine.
	Actually there is no REAL*4 parameter type, there is only REAL*8.
	This function calls PIL_get_real, then converts real*8 to real*4

**************************************************************************************/

int	PILGetReal4(const char *name, float *result)
 { double d;
   int r;

   if (NULL == result) return(PIL_NUL_PTR);
/* SPR 3302: prevent UMR if the parameter passed was not initialized, by
   removing the line which formerly read: d = *result; */
   r = PIL_get_real(&PILDefaultPFile, name, &d);
   if (PIL_OK == r) *result = d;
   return(r);
 }


/**************************************************************************************

Function:	PILGetString
Description:	get (ask) for a value of string parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	This is an official API routine.
	'result' variable has to be at least PIL_LINESIZE bytes long.

**************************************************************************************/

int	PILGetString(const char *name, char *result)
 { return(PIL_get_string(&PILDefaultPFile, name, result)); }


/**************************************************************************************

Function:	PILGetFname
Description:	get (ask) for a value of file name parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	This is an official API routine.
	'result' variable has to be at least PIL_LINESIZE bytes long.
	File access mode is checked (if required by parameter type)

**************************************************************************************/

int	PILGetFname(const char *name, char *result)
 { return(PIL_get_fname(&PILDefaultPFile, name, result)); }


/**************************************************************************************

Function:	PILGetDOL
Description:	get (ask) for a value of DOL parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	This is an official API routine.
	'result' variable has to be at least PIL_LINESIZE bytes long.
	Technically, DOL is a string parameter.

**************************************************************************************/

int	PILGetDOL(const char *name, char *result)
 { return(PIL_get_dol(&PILDefaultPFile, name, result)); }

/**************************************************************************************

Function:	PILGetAsString
Description:	get (ask) for a value of any parameter type, and return it as a string.
		This function was added to address SCREW 1497.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful value (or new value) is returned in 'result'
	This is an official API routine.
	'result' variable has to be at least PIL_LINESIZE bytes long.

**************************************************************************************/

int	PILGetAsString(const char *name, char *result)
 { return(PIL_get_as_string(&PILDefaultPFile, name, result)); }


/**************************************************************************************

Function:	PILPutBool
Description:	set new value of boolean parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
b (In)			returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine.

**************************************************************************************/

int     PILPutBool(const char *name, int b)
 { return(PIL_put_bool(&PILDefaultPFile, name, b)); }


/**************************************************************************************

Function:	PILPutInt
Description:	set new value of integer parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
i (In)			returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine.

**************************************************************************************/

int     PILPutInt(const char *name, int i)
 { return(PIL_put_int(&PILDefaultPFile, name, i)); }


/**************************************************************************************

Function:	PILPutReal
Description:	set new value of double (real*8) parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
d (In)			returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine.

**************************************************************************************/

int     PILPutReal(const char *name, double d)
 { return(PIL_put_real(&PILDefaultPFile, name, d)); }


/**************************************************************************************

Function:	PILPutString
Description:	set new value of string parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
s (In)			returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine.
	Function makes a copy of string s, so application is free to deallocate it
	after function terminates.

**************************************************************************************/

int     PILPutString(const char *name, const char *s)
 { return(PIL_put_string(&PILDefaultPFile, name, s)); }


/**************************************************************************************

Function:	PILPutFname
Description:	set new value of file name parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
s (In)			returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine.
	Function makes a copy of string s, so application is free to deallocate it
	after function terminates.

**************************************************************************************/

int     PILPutFname(const char *name, const char *s)
 { return(PIL_put_fname(&PILDefaultPFile, name, s)); }


/**************************************************************************************

Function:	PILSetModuleName
Description:	sets the name of module which uses PIL services.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name (ascii string)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is official API routine.

Warning:result is stored in global variable.

**************************************************************************************/

int	PILSetModuleName(const char *name)
 {
   if (NULL != name)
     { strncpy(PILModuleName, name, PIL_LINESIZE);
       PILModuleName[PIL_LINESIZE - 1] = 0;
     }

   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILSetModuleVersion
Description:	sets the version of module which uses PIL services
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
version (In)		version (ascii string)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	if NULL pointer is passed, then version is set to "version unspecified"
	This is an official API routine.

Warning:result is stored in global variable. Version is an ascii string.

**************************************************************************************/

int	PILSetModuleVersion(const char *version)
 {
   if (NULL != version)
     { strncpy(PILModuleVersion, version, PIL_LINESIZE);
       PILModuleVersion[PIL_LINESIZE - 1] = 0;
     }
   else strcpy(PILModuleVersion, "version unspecified");

   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILInit
Description:	initialize PIL library
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
argc (In)		number of items in parameter list (usually from main)
argv (In)		parameter list (usually from main)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This function should be called before any other PIL services are accessed.
	It does the following:
	Based on argv[0] and/or PILModuleName it calculates name of parameter file,
	then scans directories specified in PFILES environment variable for matching
	file. If not found error is reported.
	Then it calls PIL_pfile_open to open and read parameter file. Then it browses
	through list of read-in parameter, looking for parameter "ServerMode".
	If found and if its value is "yes" or "y" it sets global variable
	PILRunMode to ISDC_SERVER_MODE. Otherwise it sets PILrunMode to
	ISDC_SINGLE_MODE.
	This is an official API routine.

**************************************************************************************/

int	PILInit(int argc, char **argv)
 { static const int systries = 5;
   static const int usrtries = 10;
   /* SPR 3169: make sure tmpfile is constructed with reasonable values. */
   PIL_PFILE tmpfile = { 0, NULL, NULL, 0, NULL, 0, 0, NULL};
   time_t systime = 0;
   time_t usrtime = 0;
   char *p;
   int i, r, mode;

   if (argc <= 0) return(PIL_BAD_ARG);
   if (NULL == argv) return(PIL_NUL_PTR);

   PILLoggerFunction = NULL;

   if (0 == PILModuleName[0])
     { if (NULL != (p = strrchr(argv[0], sPathDelim))) { p++; }
       else { p = argv[0]; }	/* get base for param file name from argv[0] */
#ifdef WIN32
       { char *tmp_p;
         if ((NULL != (tmp_p = strrchr(p, '.'))) && !strcmp(tmp_p, ".exe")) *tmp_p = 0;
       }
#endif
       strncpy(PILModuleName, p, PIL_LINESIZE - 1);
       PILModuleName[PIL_LINESIZE - 1] = 0;
     }
   else
     { p = PILModuleName;	/* Use module as a base for param file name */
     }
   
   if (PIL_OK != (r = PIL_find_pfiles(p, PILUsrPFile, PILSysPFile)))
     { PILModuleName[0] = 0;
       return(r);
     }

   if (*PILSysPFile) /* If there is a system file, use it. */
     { if (PIL_OK != (r = PIL_get_mod_time(PILSysPFile, &systime)))
         { return(r); }

       r = PIL_get_mod_time(PILUsrPFile, &usrtime);

       /* If user parameter file does not exist or is out of date,
          copy the system file, then open the copy. */
       if (PIL_OK != r || systime >= usrtime)
         { if (PIL_OK != (r = PIL_file_copy(PILUsrPFile, PILSysPFile)))
             return(r);

           r = PIL_pfile_open(PILUsrPFile, argc, argv, PIL_RDWRITE, usrtries, &PILDefaultPFile);
         }
       else
         { r = PIL_pfile_open(PILSysPFile, argc, argv, PIL_RDONLY, systries, &PILDefaultPFile);
           if (PIL_OK != r) return(r);

           r = PIL_pfile_open(PILUsrPFile, argc, argv, PIL_RDWRITE, usrtries, &tmpfile);

           if (PIL_OK != r)
             { PIL_pfile_close(&PILDefaultPFile, -1);
               return(r);
             }

           r = PIL_merge_pfiles(&tmpfile, &PILDefaultPFile);

           /* Swap system and user par files; PILDefaultPFile will then
              contain the parameters from the system par file, possibly
              modified by the user parameter file, but the
              name, mode, and FILE * from the user par file. */
           if (PIL_OK == r)
             { FILE *fp = PILDefaultPFile.fp; /* Swap system and user par files. */
               int pil_access = PILDefaultPFile.pil_access;
               char *fname = PILDefaultPFile.fname;

               PILDefaultPFile.fp = tmpfile.fp;
               tmpfile.fp = fp;

               PILDefaultPFile.pil_access = tmpfile.pil_access;
               tmpfile.pil_access = pil_access;

               PILDefaultPFile.fname = tmpfile.fname;
               tmpfile.fname = fname;
             }
         }
     }
   else
     { r = PIL_pfile_open(PILUsrPFile, argc, argv, PIL_RDWRITE, usrtries, &PILDefaultPFile);
     }

   if (PIL_OK != r)
     { PIL_pfile_close(&PILDefaultPFile, r);
       return(r);
     }

   /* Close the tmpfile with a non-OK status, to make sure it is not written. */
   PIL_pfile_close(&tmpfile, -1);

   if (PIL_OK != r)
     { PILModuleName[0] = 0;
       PIL_pfile_close(&PILDefaultPFile, r);
       return(r);
     }

   mode = ISDC_SINGLE_MODE;

   if (PIL_OK == PILGetBool("ServerMode", &mode))
     { for (i=1; i< argc; i++)
	{ if (0 == strcmp("-s", argv[i]))  { mode = ISDC_SERVER_MODE; break; }
	  if (0 == strcmp("-c", argv[i]))  { mode = ISDC_SINGLE_MODE; break; }
	}
     }
   else
     { for (i=1; i< argc; i++)
	 { if (0 == strcmp("-c", argv[i]))  { mode = ISDC_SINGLE_MODE; break; }
	   if (0 == strcmp("-s", argv[i]))  { mode = ISDC_SERVER_MODE; break; }
	   if (0 == strcmp("ServerMode=n", argv[i]))  { mode = ISDC_SINGLE_MODE; break; }
	   if (0 == strcmp("ServerMode=y", argv[i]))  { mode = ISDC_SERVER_MODE; break; }
	   if (0 == strcmp("ServerMode=no", argv[i]))  { mode = ISDC_SINGLE_MODE; break; }
	   if (0 == strcmp("ServerMode=yes", argv[i]))  { mode = ISDC_SERVER_MODE; break; }
	 }
     }
   PILRunMode = mode;
   if (ISDC_SERVER_MODE == PILRunMode) PILServerState = ISDC_CONTINUE;
   else PILServerState = ISDC_TERMINATE;

   return(r);
 }


/**************************************************************************************

Function:	PILClose
Description:	shutdown PIL services
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
status (In)		status code to be reported to parent process

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If status is non-zero ( != PIL_OK) all changes made to parameters so far
	will be lost. Always pass value of zero (or PIL_OK) is you want parameter
	file to be updated (even if you subsequently return error code)
	Function also clears PILModuleName and PILModuleVersion global variables.
	This is official API routine.

Warning:This function does not terminate a process.

**************************************************************************************/

int	PILClose(int status)
 { int r;

   r = PIL_pfile_close(&PILDefaultPFile, status);
   PILModuleName[0] = 0;
   PILModuleVersion[0] = 0;
   PILLoggerFunction = NULL;
   return(r);
 }


/**************************************************************************************

Function:	PILF90AddArg
Description:	Add string to a list of F90 supplied parameters
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
arg (In)		argument in question (ascii string)

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine (used by F90 routines)
	This helper function provides F90 programs with C/C++ argc/argv support.
	There is a need for such a function when main module is written in F90, thus
	PIL would be unable to access argc/argv.
	This function is called from F90 code during F90 PIL initilization. F90 code
	subsequently calls PILF90Init.
	See also PILF90DeleteArgs

Warning:function modifies global variables PILF90argc and PILF90argv

**************************************************************************************/

int	PILF90AddArg(const char *arg)
 { char **p = 0;

   if (NULL == arg) return(PIL_BAD_ARG);
   if (PILF90argc < 0) return(PIL_BAD_ARG);
   
   if (0 == PILF90argc)
     { p = (char **)PIL_malloc(sizeof(char *));
     }
   if (PILF90argc > 0)
     { p = (char **)PIL_realloc((void *)PILF90argv, (PILF90argc + 1) * sizeof(char *));
     }
   if (NULL == p) return(PIL_NO_MEM);
   PILF90argv = p;
   PILF90argv[PILF90argc] = (char *)PIL_malloc(strlen(arg) + 1);
   if (NULL == PILF90argv[PILF90argc]) return(PIL_NO_MEM);
   strcpy(PILF90argv[PILF90argc], arg);
   PILF90argc++;
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILF90DeleteArgs
Description:	Delete F90 supplied parameter list
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine (used by F90 routines)
	This helper function provides F90 programs with C/C++ argc/argv support.
	There is a need for such a function when main module is written in F90, thus
	PIL would be unable to access argc/argv.
	This function is called from F90 code during F90 PIL shutdown.
	See also PILF90AddArg

Warning:function modifies global variables PILF90argc and PILF90argv

**************************************************************************************/

int	PILF90DeleteArgs(void)
 { int i;

   if (PILF90argc < 0) return(PIL_BAD_ARG);
   if (0 == PILF90argc) return(PIL_OK);
   if (NULL == PILF90argv) return(PIL_BAD_ARG);

   for (i=0; i<PILF90argc; i++)
    { if (NULL != PILF90argv[i]) PIL_free(PILF90argv[i]);
    }
   PIL_free(PILF90argv);
   PILF90argv = NULL;
   PILF90argc = 0;
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILF90Init
Description:	F90 equivalent for C/C++ PILInit
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine (used by F90 routines)
	This function is called from F90 code during F90 PIL initialization.
	It passes PILF90argc anf PILf90argv variables to PILIinit function.
	See also PILF90AddArg, PILIinit

Warning:function modifies global variables PILF90argc and PILF90argv

**************************************************************************************/

int	PILF90Init(void)
 {
   return(PILInit(PILF90argc, PILF90argv));
 }


/**************************************************************************************

Function:	PILF90Close
Description:	F90 equivalent to C/C++ PILClose
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
status (In)		status code to be reported to parent process

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine (used by F90 routines)
	This function is called from F90 code during F90 PIL shutdown.
	It call PILClose, passing it status parameter, then calls PILF90DeleteArgs
	to clean-up F90 parameter list.
	See also PILF90DeleteArgs, PILClose

Warning:this function does not terminate the process.

**************************************************************************************/

int	PILF90Close(int status)
 { int r, r2;

   r = PILClose(status);
   r2 = PILF90DeleteArgs();
   return((PIL_OK == r) ? r : r2);
 }


/**************************************************************************************

Function:	PILGetParFilename
Description:	get full path of used parameter file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
*fname (Out)		name/path of the parameter file

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	Absolute path is returned only when PFILES environment variable contains 
	absolute paths. If parameter file is taken from current dir then only
	filename is returned.
	This is an official API routine

Warning:Pointer returned points to a statically allocated buffer, applications should
	copy data from it using strcpy.

**************************************************************************************/

int	PILGetParFilename(char **fname)
 { if (NULL == fname) return(PIL_NUL_PTR);
   if (0 == PILUsrPFile[0]) return(PIL_NO_FILE);
   *fname = PILUsrPFile;
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILLockPFile
Description:	Lock parameter file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	lock parameter file for exclusive access.

Warning:Applications should not call this function directly. It is used internally
	by PIL_flush_parameters and PIL_reload_parameters.
	
**************************************************************************************/

int	PILLockPFile(void)
 { return(PIL_lock_pfile(&PILDefaultPFile, PIL_RDWRITE));
 }


/**************************************************************************************

Function:	PILUnlockPFile
Description:	
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	unlock parameter file.

Warning:Applications should not call this function directly. It is used internally
	by PIL_flush_parameters and PIL_reload_parameters.

**************************************************************************************/

int	PILUnlockPFile(void)
 { return(PIL_unlock_pfile(&PILDefaultPFile));
 }


/**************************************************************************************

Function:	PILReloadParameters
Description:	reload parameters from parameter file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	this function should be used by applications running in ISDC_SERVER_MODE
	to rescan parameter file and reload all parameters from it.
        This is an official API routine
	This function simply calls PIL_reload_parameters passing to it pointer
	to PILDefaultPFile structure.

Warning:Current parameter list in memory (including any modifications) is
	deleted.

**************************************************************************************/

int	PILReloadParameters(void)
 { return(PIL_reload_parameters(&PILDefaultPFile));
 }


/**************************************************************************************

Function:	PILFlushParameters
Description:	flush changes made to parameter list (in memory) to disk
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine
	This function simply calls PIL_flush_parameters passing to it pointer
	to PILDefaultPFile structure.

Warning:Current conents of parameter file is overwritten.
        
**************************************************************************************/

int	PILFlushParameters(void)
 { return(PIL_flush_parameters(&PILDefaultPFile));
 }


/**************************************************************************************

Function:	PILOverrideQueryMode
Description:	override quering user
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
newmode (In)		new value for query override mode

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	This is an official API routine
	When newmode is PIL_QUERY_OVERRIDE, prompting for new value of parameter
	is completely disabled. If value is bad or out of range PILGetXXX return
	error without asking user. No i/o in stdin/stdout is done by PIL in this
	mode. When newmode is PIL_QUERY_DEFAULT PIL reverts to default query
	mode.

Warning:When PIL_QUERY_OVERRIDE is in effect "query" mode value is ignored.
        
**************************************************************************************/

int	PILOverrideQueryMode(int newmode)
 { switch (newmode)
    { case PIL_QUERY_OVERRIDE : ;
      case PIL_QUERY_DEFAULT :	PILQueryMode = newmode;
				break;
      default : return(PIL_BAD_ARG);
    }
   return(PIL_OK);
 }

int	PILGetQueryMode(void)
 { return PILQueryMode; }

/**************************************************************************************

Function:	PILGetIntVector
Description:	get (ask) for a value of integer vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
nelem (In)		number of values to be read
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful values (or new values) is returned in 'result'
	This is an official API routine.

**************************************************************************************/

int	PILGetIntVector(const char *name, int nelem, int *result)
 { char buf[PIL_LINESIZE], buf2[2 * PIL_LINESIZE + 1000];
   int	r;

   if (NULL == result) return(PIL_NUL_PTR);
   PILVectorType = PIL_TYPE_INT4;
   PILVectorLen = nelem;
   r = PIL_get_string(&PILDefaultPFile, name, buf);
   PILVectorType = PIL_TYPE_STRING;
   PILVectorLen = 0;
   if (PIL_OK == r)  r = PILCheckIntVector(buf, nelem, result);
   if (PIL_OK != r)
     { sprintf(buf2, "PILGetIntVector(%s, %d)", name, nelem);
       PIL_log_error(buf2, r);
     }
   return(r);
 }


/**************************************************************************************

Function:	PILCheckIntVector
Description:	convert/check value of integer vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
buf (In)		input ascii string to check
nelem (In)		number of values to be read
result (Out)		returned value of parameter, may be NULL

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful, then values are returned in 'result'

**************************************************************************************/


int	PILCheckIntVector(const char *buf, int nelem, int *result)
 { int		i, ival, nbytes;
   char		dmy[PIL_LINESIZE];

   for (i = 0; i<nelem; i++)
    { 
      if (1 != sscanf(buf, "%d%n", &ival, &nbytes)) return(PIL_PFILES_FORMAT);
      if (NULL != result) result[i] = ival;
      buf += nbytes;
    }

   if (0 < sscanf(buf, "%s", dmy)) return(PIL_PFILES_FORMAT);
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILGetIntVarVector
Description:	get (ask) for a value of integer vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
nelem (In/Out)		max number of values to be read, out: values read
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful values (or new values) is returned in 'result'
	This is an official API routine.

**************************************************************************************/

int	PILGetIntVarVector(const char *name, int *nelem, int *result)
 { char buf[PIL_LINESIZE], buf2[2 * PIL_LINESIZE + 1000];
   int	r;

   if (NULL == result) return(PIL_NUL_PTR);
   if (NULL == nelem) return(PIL_NUL_PTR);
   PILVectorType = PIL_TYPE_INT4;
   PILVectorLen = *nelem;
   PILVectorVarMode = 1;
   r = PIL_get_string(&PILDefaultPFile, name, buf);
   PILVectorVarMode = 0;
   PILVectorType = PIL_TYPE_STRING;
   PILVectorLen = 0;
   if (PIL_OK == r)  r = PILCheckIntVarVector(buf, nelem, result);
   if (PIL_OK != r)
     { sprintf(buf2, "PILGetIntVarVector(%s)", name);
       PIL_log_error(buf2, r);
     }
   return(r);
 }


/**************************************************************************************

Function:	PILCheckIntVarVector
Description:	convert/check value of integer vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
buf (In)		input ascii string to check
nelem (In)		number of values to be read
result (Out)		returned value of parameter, may be NULL

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful, then values are returned in 'result'

**************************************************************************************/


int	PILCheckIntVarVector(const char *buf, int *nelem, int *result)
 { int		i, ival, nbytes, maxelem;
   char		dmy[PIL_LINESIZE];

   if (NULL == nelem) return(PIL_BAD_ARG);
   maxelem = *nelem;
   for (i = 0; i<maxelem; i++)
    { 
      if (1 != sscanf(buf, "%d %n", &ival, &nbytes))
         { if (sscanf(buf, "%s", dmy) > 0)  return(PIL_PFILES_FORMAT);
           break;
         }
      if (NULL != result) result[i] = ival;
      buf += nbytes;
    }

   *nelem = i;
   return(PIL_OK);
 }

/**************************************************************************************

Function:	PILGetRealVector
Description:	get (ask) for a value of real vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
nelem (In)		number of values to be read
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful values (or new values) is returned in 'result'
	This is an official API routine.

**************************************************************************************/

int	PILGetRealVector(const char *name, int nelem, double *result)
 { char buf[PIL_LINESIZE], buf2[2 * PIL_LINESIZE + 1000];
   int	r;

   if (NULL == result) return(PIL_NUL_PTR);
   PILVectorType = PIL_TYPE_REAL8;
   PILVectorLen = nelem;
   r = PIL_get_string(&PILDefaultPFile, name, buf);
   PILVectorType = PIL_TYPE_STRING;
   PILVectorLen = 0;
   if (PIL_OK == r)  r = PILCheckRealVector(buf, nelem, result);
   if (PIL_OK != r)
     { sprintf(buf2, "PILGetRealVector(%s, %d)", name, nelem);
       PIL_log_error(buf2, r);
     }
   return(r);
 }


/**************************************************************************************

Function:	PILCheckRealVector
Description:	convert/check value of real vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
buf (In)		input ascii string to check
nelem (In)		number of values to be read
result (Out)		returned value of parameter, may be NULL

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful, then values are returned in 'result'

**************************************************************************************/


int	PILCheckRealVector(const char *buf, int nelem, double *result)
 { double	ival;
   int		i, nbytes;
   char		dmy[PIL_LINESIZE];

   for (i = 0; i<nelem; i++)
    { 
      if (1 != sscanf(buf, "%lg%n", &ival, &nbytes)) return(PIL_PFILES_FORMAT);
      if (NULL != result) result[i] = ival;
      buf += nbytes;
    }

   if (0 < sscanf(buf, "%s", dmy)) return(PIL_PFILES_FORMAT);
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILGetRealVarVector
Description:	get (ask) for a value of integer vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
nelem (In/Out)		max number of values to be read, out: values read
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful values (or new values) is returned in 'result'
	This is an official API routine.

**************************************************************************************/

int	PILGetRealVarVector(const char *name, int *nelem, double *result)
 { char buf[PIL_LINESIZE], buf2[2 * PIL_LINESIZE + 1000];
   int	r;

   if (NULL == result) return(PIL_NUL_PTR);
   if (NULL == nelem) return(PIL_NUL_PTR);
   PILVectorType = PIL_TYPE_REAL4;
   PILVectorLen = *nelem;
   PILVectorVarMode = 1;
   r = PIL_get_string(&PILDefaultPFile, name, buf);
   PILVectorVarMode = 0;
   PILVectorType = PIL_TYPE_STRING;
   PILVectorLen = 0;
   if (PIL_OK == r)  r = PILCheckRealVarVector(buf, nelem, result);
   if (PIL_OK != r)
     { sprintf(buf2, "PILGetRealVarVector(%s)", name);
       PIL_log_error(buf2, r);
     }
   return(r);
 }


/**************************************************************************************

Function:	PILCheckRealVarVector
Description:	convert/check value of integer vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
buf (In)		input ascii string to check
nelem (In/Out)		max number of values to be read, out: actual number read
result (Out)		returned value of parameter, may be NULL

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful, then values are returned in 'result'

**************************************************************************************/


int	PILCheckRealVarVector(const char *buf, int *nelem, double *result)
 { int		i, nbytes, maxelem;
   char		dmy[PIL_LINESIZE];
   double	dval;

   if (NULL == nelem) return(PIL_BAD_ARG);
   maxelem = *nelem;
   for (i = 0; i<maxelem; i++)
    { 
      if (1 != sscanf(buf, "%lg %n", &dval, &nbytes))
        { if (sscanf(buf, "%s", dmy) > 0)  return(PIL_PFILES_FORMAT);
          break;
        }
      if (NULL != result) result[i] = dval;
      buf += nbytes;
    }

   *nelem = i;
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILGetReal4Vector
Description:	get (ask) for a value of real vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
nelem (In)		number of values to be read
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful values (or new values) is returned in 'result'
	This is an official API routine.

**************************************************************************************/

int	PILGetReal4Vector(const char *name, int nelem, float *result)
 { char buf[PIL_LINESIZE], buf2[2 * PIL_LINESIZE + 1000];
   int	r;

   if (NULL == result) return(PIL_NUL_PTR);
   PILVectorType = PIL_TYPE_REAL4;
   PILVectorLen = nelem;
   r = PIL_get_string(&PILDefaultPFile, name, buf);
   PILVectorType = PIL_TYPE_STRING;
   PILVectorLen = 0;
   if (PIL_OK == r)  r = PILCheckReal4Vector(buf, nelem, result);
   if (PIL_OK != r)
     { sprintf(buf2, "PILGetReal4Vector(%s, %d)", name, nelem);
       PIL_log_error(buf2, r);
     }
   return(r);
 }


/**************************************************************************************

Function:	PILCheckReal4Vector
Description:	convert/check value of real vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
buf (In)		input ascii string to check
nelem (In)		number of values to be read
result (Out)		returned value of parameter, may be NULL

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful, then values are returned in 'result'

**************************************************************************************/


int	PILCheckReal4Vector(const char *buf, int nelem, float *result)
 { float	ival;
   int		i, nbytes;
   char		dmy[PIL_LINESIZE];

   for (i = 0; i<nelem; i++)
    { 
      if (1 != sscanf(buf, "%g%n", &ival, &nbytes)) return(PIL_PFILES_FORMAT);
      if (NULL != result) result[i] = ival;
      buf += nbytes;
    }

   if (0 < sscanf(buf, "%s", dmy)) return(PIL_PFILES_FORMAT);
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILGetReal4VarVector
Description:	get (ask) for a value of integer vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)		name of parameter in question
nelem (In/Out)		max number of values to be read, out: values read
result (Out)		returned value of parameter

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful values (or new values) is returned in 'result'
	This is an official API routine.

**************************************************************************************/

int	PILGetReal4VarVector(const char *name, int *nelem, float *result)
 { char buf[PIL_LINESIZE], buf2[2 * PIL_LINESIZE + 1000];
   int	r;

   if (NULL == result) return(PIL_NUL_PTR);
   if (NULL == nelem) return(PIL_NUL_PTR);
   PILVectorType = PIL_TYPE_REAL4;
   PILVectorLen = *nelem;
   PILVectorVarMode = 1;
   r = PIL_get_string(&PILDefaultPFile, name, buf);
   PILVectorVarMode = 0;
   PILVectorType = PIL_TYPE_STRING;
   PILVectorLen = 0;
   if (PIL_OK == r)  r = PILCheckReal4VarVector(buf, nelem, result);
   if (PIL_OK != r)
     { sprintf(buf2, "PILGetReal4VarVector(%s)", name);
       PIL_log_error(buf2, r);
     }
   return(r);
 }


/**************************************************************************************

Function:	PILCheckReal4VarVector
Description:	convert/check value of integer vector parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
buf (In)		input ascii string to check
nelem (In/Out)		max number of values to be read, out: actual number read
result (Out)		returned value of parameter, may be NULL

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful, then values are returned in 'result'

**************************************************************************************/


int	PILCheckReal4VarVector(const char *buf, int *nelem, float *result)
 { int		i, nbytes, maxelem;
   char		dmy[PIL_LINESIZE];
   float	fval;

   if (NULL == nelem) return(PIL_BAD_ARG);
   maxelem = *nelem;
   for (i = 0; i<maxelem; i++)
    { 
      if (1 != sscanf(buf, "%g %n", &fval, &nbytes))
        { if (sscanf(buf, "%s", dmy) > 0)  return(PIL_PFILES_FORMAT);
          break;
        }
      if (NULL != result) result[i] = fval;
      buf += nbytes;
    }

   *nelem = i;
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILGetNumParameters
Description:	Return number of parameters in the parameter file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
parnum (Out)		returned number of parameters

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	If successful, then values are returned in 'parnum'

**************************************************************************************/


int	PILGetNumParameters(int *parnum)
{
  if (NULL == parnum) return(PIL_NUL_PTR);
  *parnum = PILDefaultPFile.pcnt;
  return(PIL_OK);
}


/**************************************************************************************

Function:	PILGetParameter
Description:	Return number of parameters in the parameter file
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
idx (In)	index of parameter to return
pp (Out)	returned parameter's data
minmaxok (Out)	flag whether min/max(1) or enum(2) values are defined (and returned) for that parameter
vmin (Out)	min value for returned parameter (converted to proper type), only when minmaxok==1
vmax (Out)	max value for returned parameter (converted to proper type), only when minmaxok==1

Return code:	standard PIL error code or PIL_OK (zero)

Notes : one can use the following program to list parameters :

PILGetNumParameters(&parcnt);
for (i=0; i<parcnt; i++)
 { if (PIL_OK != PILGetParameter(i, &pardata, &minmaxflag, &minval, &maxval)) break;
   if (PIL_FORMAT_OK != pp->format) continue;

   minmaxstr[0] = 0;
   if (1 == minmaxflag) sprintf(minmaxstr, "min=%s, max=%s ", pp->strmin, pp->strmax);
   if (2 == minmaxflag) sprintf(minmaxstr, "enum=%s ", pp->strmin);

   printf("%-16.16s 0x%02x 0x%02x %20s %s// %s\n", pp->strname, pp->type, pp->mode,
			pp->strvalue, minmaxstr, pp->strprompt);
 }

Symbolic values defined in pil.h for pp->type are :

#define	PIL_TYPE_BOOL		(1)
#define	PIL_TYPE_INT4		(2)
#define	PIL_TYPE_REAL8		(3)
#define	PIL_TYPE_STRING		(4)
#define	PIL_TYPE_FNAME		(5)

and for pp->mode are (xor-ed when more than one active at a time) :

#define	PIL_MODE_AUTO		(0x1)
#define	PIL_MODE_HIDDEN		(0x2)
#define	PIL_MODE_LEARN		(0x4)
#define	PIL_MODE_QUERY		(0x8)

and for pp->format :

#define	PIL_FORMAT_OK		(0x0)
#define	PIL_FORMAT_BLANK	(0x1)
#define	PIL_FORMAT_COMMENT	(0x2)
#define	PIL_FORMAT_ERR		(0x4)

**************************************************************************************/


int	PILGetParameter(int idx, PIL_PARAM *pp, int *minmaxok, PIL_VALUE *vmin, PIL_VALUE *vmax)
{ PIL_VALUE	*venum;
  int		venumcount;

  if (NULL == pp) return(PIL_NUL_PTR);
  if (NULL == minmaxok) return(PIL_NUL_PTR);
  if (NULL == vmin) return(PIL_NUL_PTR);
  if (NULL == vmax) return(PIL_NUL_PTR);
  if ((idx < 0) || (idx >= PILDefaultPFile.pcnt)) return(PIL_BAD_ARG);
  *pp = PILDefaultPFile.pp[idx];
  *minmaxok = 0;				/* NO min/max/enum value constraints */
  if (PIL_FORMAT_OK == pp->format)
    { if (PIL_OK == PIL_get_enum_range(&PILDefaultPFile, idx, vmin, vmax, &venum, &venumcount))
	{ if (-1 == venumcount)
	    { *minmaxok = 1;			/* min/max defined */
	    }
	  else
	    { *minmaxok = 2;			/* enum list defined (in min field) */
	      if (NULL != venum)  { PIL_free((void *)venum); venum = NULL; }
	    }
	}
    }
  return(PIL_OK);
}


/**************************************************************************************

Function:	PILVerifyCmdLine
Description:	return value of given parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
Return code:	standard PIL error code or PIL_OK (zero)

Notes:	PILVerifyCmdLine scans argument list (given by argc and argv parameters) and
	for all parameters in format Name=Value checks whether there is parameter with
	such name in parameter table (in memory). If this is not the case it returns
	with an error (meaning: bogus parameters specified in command line).

**************************************************************************************/

int	PILVerifyCmdLine(void)
 { char	buf[PIL_LINESIZE], buf2[PIL_LINESIZE + 100 + 65000];
   int	i, j, r, k, ok, delta;

   if (PILDefaultPFile.argc <= 0) return(PIL_NOT_FOUND);
   if (NULL == PILDefaultPFile.argv) return(PIL_NUL_PTR);

   switch (PILCmndArgMode)
    { case PIL_CMNDARG_ALLOW_SPACES:

	   for (i = 0; i < PILDefaultPFile.argc; i++)
	    { if (NULL == PILDefaultPFile.argv[i]) continue;
	      strncpy(buf, PILDefaultPFile.argv[i], PIL_LINESIZE - 1);
	      buf[PIL_LINESIZE - 1] = 0;
	      j = PIL_find_equal(buf);		/* find '=' in the string, neg. offset means not found */
	      if (j < 0) continue;
	      if ((0 == j) && (0 == i)) continue;
	      ok = 0;
	      if (j == 0)			/* solitary '=' */
	        { strncpy(buf, PILDefaultPFile.argv[i-1], PIL_LINESIZE-1);
		  buf[PIL_LINESIZE - 1] = 0;
	        }
	      else
	        { buf[j] = 0; }
              /* following 3 lines inserted by WDP */
              if (strchr(buf, '[')) continue;	/* this is a positional parameter */
						/* consisting of a virtual filename appended */
						/* with a qualifying filter containing an '=' */
	      for (k = 0; k < PILDefaultPFile.pcnt; k++)
	       { if (PIL_FORMAT_OK != PILDefaultPFile.pp[k].format) continue;
		 if (0 == strcmp(buf,PILDefaultPFile.pp[k].strname))
		   { ok = 1;
		     break;
		   }
	       }
	      if (ok == 0)
	        { r = PIL_BOGUS_CMDLINE;
		  strcpy(buf2, "PILVerifyCmdLine :");
						/* a long command line may overflow buf2[] !!! */
		  for (k = 0; k < PILDefaultPFile.argc; k++)
		    { strcat(buf2, " ");
		      strcat(buf2, PILDefaultPFile.argv[k]);
		    }
		  PIL_log_error(buf2, r);
		  strcat(buf2, "PILVerifyCmdLine : bogus parameter name = ");
		  strcat(buf2, buf);
		  PIL_log_error(buf2, r);
		  return(r);
	      }
	    }
	   break;
      default:
	   for (i = 1; i < PILDefaultPFile.argc; i++)
	    { if (NULL == PILDefaultPFile.argv[i]) continue;
	      strncpy(buf, PILDefaultPFile.argv[i], PIL_LINESIZE - 1);
	      buf[PIL_LINESIZE - 1] = 0;
	
	      if (0 == strcmp("ServerMode=n", buf)) continue;
	      if (0 == strcmp("ServerMode=y", buf)) continue;
	      if (0 == strcmp("ServerMode=no", buf)) continue;
	      if (0 == strcmp("ServerMode=yes", buf)) continue;
	
	      delta = 0;
	      k = PIL_find_equal(buf);
	      if (k < 0)
	        {
		  if ((i + 1) >= PILDefaultPFile.argc)  continue;
		  if (NULL == PILDefaultPFile.argv[i + 1])  continue;
		  if ('=' != PILDefaultPFile.argv[i + 1][0])  continue;
		  if (0 != PILDefaultPFile.argv[i + 1][1])  { delta = 1; }
		  else  { delta = 2; }
	        }
	      else
	        {
	          buf[k] = 0;		/* name is extracted now */
	          if (0 == buf[k + 1])  delta = 1;  /* argv[i+1] = value, skip it for next iteration */
                  if (strchr(buf, '[')) continue;  /* this is a positional parameter consisting of a */
                     /* virtual filename appended with a qualifying filter containing an '=' */
	        }
	
	      ok = 0;		/* go, find matching entry in the table */
	      for (j=0; j<PILDefaultPFile.pcnt; j++)
	       { if (PIL_FORMAT_OK != PILDefaultPFile.pp[j].format) continue;
		 if (0 == strcmp(buf, PILDefaultPFile.pp[j].strname))
		   { ok = 1;
		     break;
		   }
	       }
	      if (0 == ok)
	        { r = PIL_BOGUS_CMDLINE;
		  sprintf(buf2, "PILVerifyCmdLine : bogus parameter name = %s", buf);
		  PIL_log_error(buf2, r);
	          return(r);
	        }
	      i += delta;
	    }
	   break;
    }
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILSetRootNameFunction
Description:	return value of given parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
func (In)		Pointer to the function

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	PILSetRootNameFunction instructs PIL to use specified function during
	conversion from URL syntax to the file name. Passing NULL disables URL syntax
	checking, meaning any URL string is treated verbatim as a file name. This
	functionality is used whenever parameter of type 'fr', 'fw', 'fn' or 'fe'
	is validated by the PIL.

	Unless this function is called PIL_cfitsio_root_name() function is used
	by PIL.

**************************************************************************************/

int	PILSetRootNameFunction(int (*func)(char *s))
 {
   PILRootFnameFunction = func;
   return(PIL_OK);
 }


/**************************************************************************************

Function:	PILSetReadlinePromptMode
Description:	set READLINE prompt mode
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
mode (In)		New mode

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	PILSetReadlinePromptMode sets PIL's prompt mode when using READLINE
	library. Currently 2 modes are defined :
		PIL_RL_PROMPT_PIL - old mode, compatible with previous PIL versions.
			Default value is printed before ':'. Edit buffer is
			initially empty. So to change, say last character in a long 
			string/filename one must retype (or copy/paste with a mouse)
			whole string. To enter empty string, some tricks are needed :
			i.e. PIL_EMPTY_STRING env.var. must be defined.
		PIL_RL_PROMPT_IEB - new mode. Default value is _NOT_ printed 
			before ':'. Edit buffer is initially set to default value.
			So, to change last character in a long string/filename one
			simply needs to press Del/BACKSPACE then the new character.
			Furthermore to enter empty string one simply needs to
			press Ctrl-W (or ctrl-A then ctrl-R).

	Upon startup PIL mode is PIL_RL_PROMPT_PIL, unless STARTUP_PIL_RL_PROMPT_IEB 
	was defined during compile time.

**************************************************************************************/

int	PILSetReadlinePromptMode(int mode)
 { int	r;

   r = PIL_OK;

   switch (mode)
    { case PIL_RL_PROMPT_PIL: ;
      case PIL_RL_PROMPT_IEB :
		PILReadlinePromptMode = mode;
		break;
      default:
		r = PIL_BAD_ARG;
		break;
    }

   PIL_log_error("PILSetReadlinePromptMode", r);
   return(r);
 }


/**************************************************************************************

Function:	PILSetCmndArgMode
Description:	set command line arguments mode
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
mode (In)		New mode

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	PILSetCmndArgMode sets PIL's command line arguments mode (syntax).
	Currently there is 1 mode defined :
		PIL_CMNDARG_PIL - old mode, compatible with previous PIL versions.
			Command line argument should have one of the following formats:
				Name=Value
				Value (positional)

**************************************************************************************/

int	PILSetCmndArgMode(int mode)
 { int	r;

   r = PIL_OK;

   switch (mode)
    { case PIL_CMNDARG_PIL: ;
      case PIL_CMNDARG_ALLOW_SPACES: ;
		PILCmndArgMode = mode;
		break;
      default:
		r = PIL_BAD_ARG;
		break;
    }

   PIL_log_error("PILSetCmndArgMode", r);
   return(r);
 }



/**************************************************************************************

Function:	PILSetLoggerFunction
Description:	return value of given parameter
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
func (In)		Pointer to the function

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	PILSetLoggerFunction instructs PIL to use specified function for
	logging error messages. Passing NULL disables logging. After PILInit
	the setup is as PILSetLoggerFunction(NULL) has been called.

**************************************************************************************/

int	PILSetLoggerFunction(int (*func)(char *s))
 {
   PILLoggerFunction = func;
   return(PIL_OK);
 }



/**************************************************************************************

Function:	PILSetReprompt
Description:	SCREW 1496. Puts parameter in a state such that if one then calls any
		of the PILGet* functions, a prompt will occur even if the parameter was
		supplied on the command line, or was already prompted for, or is
		hidden. This does not trump override query mode, however.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
par_name (In)		Name of the parameter
reprompt (In)		Whether or not to reprompt (0 - do not reprompt, 1 - reprompt)

Return code:	standard PIL error code or PIL_OK (zero)

**************************************************************************************/

int	PILSetReprompt(const char *par_name, int reprompt)
 { int r = PIL_OK;
   int idx = 0;

   r = PIL_find_name(&PILDefaultPFile, par_name, &idx); 
   if (PIL_OK == r) PILDefaultPFile.pp[idx].reprompt = reprompt;
   return r;
 }
