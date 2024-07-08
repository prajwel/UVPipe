/**************************************************************************************

	I N T E G R A L   S C I E N C E   D A T A   C E N T E R

	P A R A M E T E R   I N T E R F A C E   L I B R A R Y

Copyright:	(c) 1998-2002 ISDC http://isdc.unige.ch
File:		pil.h
Description:	Parameter Interface Library - core & API header file
Authors:	Jerzy.Borkowski@obs.unige.ch (JB)
History:	10-May-98 (JB) : PILGetReal4 function added 
		13-May-98 (JB) : version 0.7 released. Code cleanup, unnecesary
			dependencies removed. Uses/links only PIL library. Changes
			to Makefile.
		13-Aug-98 (JB) : version 1.0 released. More compatible with Common
			library some new functions: PILReloadParameters,
			PILFlushParameters, locking of parameter file during
			read/write implemented.
		17-Aug-98 (JB) : reformatted, more comments added
		25-Jun-99 (JB) : added "" -> empty string substitution.
			added vector parameters support (3 new API functions)
		27-Jan-00 (JB) : char * were const-ized to support ANSI C++
		22-Aug-00 (JB) : added PILVerifyCmdLine().
		27-Sep-00 (JB) : version 1.5.1: xxxVarVector() added
		09-Feb-01 (JB) : version 1.6.3: enumerated values implemented
		20-Apr-01 (JB) : version 1.7.0 : readline support added
		03-Jul-01 (JB) : version 1.7.1 : access renamed to pil_access
		14-Aug-01 (JB) : version 1.7.2 : more access variables renamed to pil_access
		09-Mar-02 (JB) : version 1.8.0 : added root_name support

**************************************************************************************/


#ifndef	PIL_H_INCLUDED
#define	PIL_H_INCLUDED

#ifndef EXPSYM
#ifdef WIN32
#define EXPSYM __declspec(dllexport)
#else
#define EXPSYM
#endif
#endif

#ifndef IMPSYM
#ifdef WIN32
#define IMPSYM __declspec(dllimport)
#else
#define IMPSYM
#endif
#endif

/* SCREW 1145: Only include unix system headers on Unix. */
#ifndef WIN32
#include <sys/types.h>
#include <sys/stat.h>
#endif

/* SCREW 1145: Needed for time_t; it must have been included "by accident"
   on Unix by some other header file. */
#include <time.h>
#include <ctype.h>
#include <stdio.h>

#ifdef HAVE_READLINE_LIBS
#include <readline/readline.h>
#include <readline/history.h>
#endif

/* SCREW 1326: Add version identifier. */
#ifndef NO_PIL_VERSION_STRING
static const char * pilversion[] = {"ISDC component pil 2.0.0"};
#endif

#if defined(__cplusplus) && ! defined(__CINT__)
extern "C" {
#endif

		/* this is sometimes needed, avoids necessity to include isdc.h */
#ifndef	ISDC_OK
#define	ISDC_OK			(0)
#endif

			/* manifest constants */

#define	PIL_LINESIZE		(2000)
#define	PIL_CURLYSIZE		(10000)
#define	PIL_PATH_MAX		(1024)
#define	PIL_COPYCHUNK		(8192)

#define	PIL_RDONLY		(0)
#define	PIL_RDWRITE		(1)

#define	PIL_WHAT_VALUE		(0)
#define	PIL_WHAT_MIN		(1)
#define	PIL_WHAT_MAX		(2)

/* bugfix 25-Jun-99 (JB) - conditional definitions, now defined in common.h (?) */

#ifndef ISDC_SINGLE_MODE
#define	ISDC_SINGLE_MODE	(0)
#endif 

#ifndef ISDC_SERVER_MODE
#define	ISDC_SERVER_MODE	(1)
#endif 

#ifndef ISDC_CONTINUE
#define	ISDC_CONTINUE		(0)
#endif 

#ifndef ISDC_TERMINATE
#define	ISDC_TERMINATE		(1)
#endif 

	/* constants used for type */

#define	PIL_TYPE_UNKNOWN	(0)
#define	PIL_TYPE_BOOL		(1)
#define	PIL_TYPE_INT4		(2)
#define	PIL_TYPE_REAL8		(3)
#define	PIL_TYPE_STRING		(4)
#define	PIL_TYPE_FNAME		(5)
#define	PIL_TYPE_STRUCT		(6)
#define	PIL_TYPE_GCUR		(7)
#define	PIL_TYPE_IMCUR		(8)
#define	PIL_TYPE_D		(9)		/* extension, found in ftools, pobably fortran real, i.e -> 1.234D-12 */
#define	PIL_TYPE_G		(10)		/* this is probably list of ranges i.e: 1.2-3.4,5.6-7.8 */
#define	PIL_TYPE_REAL4		(11)		/* used only for vector of REAL*4 */
#define	PIL_TYPE_DOL		(12)		/* extra checks/expansion by readline - technically its STRING */

	/* constants used for attrib */

#define	PIL_FILECHECK_EXIST	(0x1)
#define	PIL_FILECHECK_NOEXIST	(0x2)
#define	PIL_FILECHECK_READ	(0x4)		/* open for read, fail if no file */
#define	PIL_FILECHECK_WRITE	(0x8)		/* open for update, fail if no file present */
#define	PIL_ATTRIB_REDIRECTION	(0x10)		/* parameter redirection */
#define	PIL_VALUE_INVALID	(0x20)		/* value does not follow type standard */
#define	PIL_VALUE_FROM_CMDLINE	(0x40)		/* value overriden by cmd line arg - do not ask if query */
#define	PIL_FILECHECK_CREATE	(0x80)		/* create new file, erase old, if any */

	/* constants used for format */

#define	PIL_FORMAT_OK		(0x0)		/* format seems to be ok, it is 7 tokens are present in line */
#define	PIL_FORMAT_BLANK	(0x1)		/* line is empty */
#define	PIL_FORMAT_COMMENT	(0x2)		/* line is a comment (1st char is a #) */
#define	PIL_FORMAT_ERR		(0x4)		/* generic format error, bad type, mode */
#define	PIL_FORMAT_TOO_MANY	(0x8)		/* too many tokens found ( >7 ) */
#define	PIL_FORMAT_TOO_FEW	(0x10)		/* not enough tokens found ( <7 ) */
#define	PIL_FORMAT_QUOTE	(0x20)		/* unbalanced (double)quote found */
#define	PIL_FORMAT_NO_LF	(0x40)		/* line without terminating LF character */
#define	PIL_FORMAT_EXTRA_SPACES	(0x80)		/* line with extra space(s)/TAB(s) after trailing (double)quote */

	/* constants used for mode */

#define	PIL_MODE_AUTO		(0x1)
#define	PIL_MODE_HIDDEN		(0x2)
#define	PIL_MODE_LEARN		(0x4)
#define	PIL_MODE_QUERY		(0x8)

	/* constants used for modified */

#define	PIL_MODIFIED_NAME	(0x1)
#define	PIL_MODIFIED_TYPE	(0x2)
#define	PIL_MODIFIED_MODE	(0x4)
#define	PIL_MODIFIED_VALUE	(0x8)
#define	PIL_MODIFIED_MIN	(0x10)
#define	PIL_MODIFIED_MAX	(0x20)
#define	PIL_MODIFIED_PROMPT	(0x40)
#define	PIL_MODIFIED_ANY	(0x7F)

#define	PIL_PARFILE_MODIFIED	(0x100)		/* for fp->mode if not set PIL_flush_parameters does nothing */
#define	PIL_SUBST_ENV		(0x200)		/* for fp->mode, honor $ENVVAR expansions */
#define	PIL_SUBST_INDIR		(0x400)		/* for fp->mode, honor )name expansions */
#define	PIL_PSET_MODE		(0x800)		/* pset/pil_set mode - _ALWAYS_ update, also pos. args shifted by 1 */
#define	PIL_NO_POSITIONAL	(0x1000)	/* disable positional args - for i.e. pget */
#define	PIL_BYPASS_NAMING	(0x2000)	/* bypass naming resolution, use directly module name */
#define	PIL_SUBST_EMPTYSTR	(0x4000)	/* for fp->mode, don't honor "" substitution for empty string */
#define	PIL_OPEN_RDONLY		(0x8000)	/* open always in read only mode, dont update on exit */
#define	PIL_SUBST_CURLY		(0x10000)	/* for fp->mode, don't honor ${envname} expansions */
#define	PIL_FLUSH_PERFORMED	(0x20000)	/* set after flush, should reset PARFILE_MODIFIED/CMDLINE */
#define	PIL_LINT_MODE		(0x40000)	/* use very strict syntax checking (for pil_lint) */

#define	PIL_MODIFIED_YES	(1)
#define	PIL_CHECK_RANGE_YES	(2)

#define	PIL_MAX_INDIR		(10)

#define	PIL_malloc(x)		malloc(x)
#define	PIL_realloc(x,y)	realloc(x,y)
#define	PIL_free(x)		free(x)

#define	PIL_PFILE_SUFFIX	".par"

#define	PIL_QUERY_OVERRIDE	(0)
#define	PIL_QUERY_DEFAULT	(1)

#define	PIL_ROOTNAME_FILE	(0)
#define	PIL_ROOTNAME_NOTFILE	(1)
#define	PIL_ROOTNAME_STDIN	(2)
#define	PIL_ROOTNAME_STDOUT	(3)
#define	PIL_ROOTNAME_STDINOUT	(4)

#define	PIL_RL_PROMPT_PIL	(0)
#define	PIL_RL_PROMPT_IEB	(1)		/* IEB = Init Edit Buffer */

#define	PIL_CMNDARG_PIL		(0)
#define	PIL_CMNDARG_ALLOW_SPACES (1)		/* Name [ Space ] Equal [ Space ] Argument */

			/* type definitions */

typedef	union	PIL_VALUE_UNION
      { int	b;		/* boolean, 0 - false, true otherwise */
	int	i;		/* integer*4 or larger */
	double	d;		/* real*8 or larger */
	char	s[PIL_LINESIZE];/* string, filename, struct (?), URL (?) */
      } PIL_VALUE;


typedef	struct	PIL_PARAM_STRUCT
      {	char	*strline;	/* line as read from file */
	char	*strname;	/* name token from line */
	char	*strtype;	/* type token from line */
	char	*strmode;	/* mode token from line */
	char	*strvalue	/* value token from line */;
	char	*strmin;	/* min token from line */
	char	*strmax;	/* max token from line */
	char	*strprompt;	/* prompt token from line */
	int	type;		/* parameter type */
	int	mode;		/* auto/learn/hidden/query */
	int	attrib;		/* for FNAME check flags,  also indirection */
	int	modified;	/* various flags to signal modification of fields took place */
	int	format;		/* is the format ok ? */
	int	reprompt;	/* are reprompts enabled? */
      } PIL_PARAM;


typedef struct PIL_PFILE_STRUCT
      {	int	argc;
	char	**argv;
	FILE	*fp;
	int 	pil_access;
	PIL_PARAM *pp;
	int	pcnt;
	int	mode;
	char	*fname; /* SCREW 1145: Filename is needed for portable file truncation etc. */
      } PIL_PFILE;

typedef struct	PIL_TYPE_TABLE_STRUCT
      {	int	type;
	const char *name;
      } PIL_TYPE_TABLE;

			/* variables */

extern	PIL_TYPE_TABLE	PIL_ttab[];			/* table of data types */

extern
#ifndef __CINT__
IMPSYM
#endif
PIL_PFILE		PILDefaultPFile;		/* handle of default PFILE for applic */
extern	int		PILRunMode;			/* mode of operation */
extern	int		PILF90argc;			/* argc from F90 main() */
extern	char		**PILF90argv;			/* argv from F90 main() */
extern	char		PILModuleName[PIL_LINESIZE];	/* name of module */
extern	char		PILModuleVersion[PIL_LINESIZE];	/* version of module */
extern	int		PILServerState;			/* state of server mode */
extern	int		PILSpecialMode;			/* extras, like PSET mode */
extern	char		*PILEmptyStringTemplate;	/* value for empty string substitution, if NULL -> \"\" */
extern	int		PILVectorType;			/* used to verify vector values */
extern	int		PILVectorLen;			/* used to verify vector values */
extern	int		PILVectorVarMode;		/* nonzero - variable length vector */
extern	char		PILLastNotFoundName[PIL_LINESIZE]; /* name of last not found parameter */
extern	int		(*PILRootFnameFunction)(char *s);  /* pointer to a URL->filename conversion function */
extern	int		PILReadlinePromptMode;		/* Prompting mode when using READLINE */
extern	int		PILCmndArgMode;			/* Command line arguments mode */
extern	int		(*PILLoggerFunction)(char *msg);

			/* function prototypes */

int	PIL_log_error(char *s, int errcode);
int	PIL_root_name(char *s);
int	PIL_cfitsio_root_name(char *s);
int	PIL_file_copy(const char *dest, const char *src);
int	PIL_allow_spaces_cmndarg2value(PIL_PFILE * fp);
int	PIL_cmndarg2value(const char *name, int numpar, char **result, int argc, char **argv);
int	PIL_verify_access(const char *name, const char *path, int mode, time_t *modtime);
int	PIL_traverse_token(const char *name, char *path, char **result, int mode);
int	PIL_build_names(char *name, int argc, char **argv);
int	PIL_find_pfiles(char *name, char *userpfile, char *syspfile);
int	PIL_stricmp(const char *s1, const char *s2);
int	PIL_trim_spaces(char *s);
int	PIL_dup(char **dst, const char *src);
#ifndef __CINT__
EXPSYM
#endif
const char *PIL_type2string(int type);
#ifndef __CINT__
EXPSYM
#endif
int	PIL_mode2string(int mode, char *buf);
int	PIL_indir2string(PIL_PFILE *fp, int what, char **dst, char *src);
int	PIL_string2value(const char *buf, int type, PIL_VALUE *v);
int	PIL_value2string(PIL_VALUE *v, int type, char *buf, int maxlen);
int	PIL_check_mode(char *ptab, int *mode);
int	PIL_check_pfile(PIL_PFILE *fp);
int	PIL_check_par(PIL_PFILE *fp, int idx);
int	PIL_merge_pfiles(PIL_PFILE *spfile, PIL_PFILE *upfile);
int	PIL_init_param(PIL_PARAM *pp);
int	PIL_free_entry(PIL_PARAM *pp);
int	PIL_delete_param(PIL_PFILE *fp, int idx);
int	PIL_append_param(PIL_PFILE *fp, PIL_PARAM *pp);
int	PIL_reload_parameters(PIL_PFILE *fp);
int	PIL_sort_parameters(PIL_PFILE *fp, PIL_PARAM ***parray);
int	PIL_pfile_open(char *fname, int argc, char **argv, int pfile_acc_mode, int numtries, PIL_PFILE *fp);
int	PIL_flush_parameters(PIL_PFILE *fp);
int	PIL_pfile_close(PIL_PFILE *fp, int status);
int	PIL_read_line(PIL_PFILE *fp, char *buf, int lineno);
int	PIL_mark_eos(char **p, char term);	/* scan up to 1st [term\n\r\0], mark EOS with \0 */
int	PIL_parse_line(PIL_PFILE *fp, char *buf, int numpar, PIL_PARAM *pp, int argc, char **argv, char *buf2);
int	PIL_find_name(PIL_PFILE *fp, const char *name, int *idx);
#ifndef __CINT__
EXPSYM
#endif
int	PIL_get_range(PIL_PFILE *fp, int idx, PIL_VALUE *vmin, PIL_VALUE *vmax);
int	PIL_get_enum_range(PIL_PFILE *fp, int idx, PIL_VALUE *vmin, PIL_VALUE *vmax, PIL_VALUE **enumlist, int *enumcount);
int	PIL_put_value(PIL_PFILE *fp, int idx, int type, PIL_VALUE *vs, const char *strv, int mode);
int     PIL_put_by_value(PIL_PFILE *fp, const char *name, int type, void *p);
int     PIL_put_by_string(PIL_PFILE *fp, const char *name, const char *src);
int     PIL_put_bool(PIL_PFILE *fp, const char *name, int b);
int     PIL_put_int(PIL_PFILE *fp, const char *name, int i);
int     PIL_put_real(PIL_PFILE *fp, const char *name, double d);
int     PIL_put_string(PIL_PFILE *fp, const char *name, const char *s);
int     PIL_put_fname(PIL_PFILE *fp, const char *name, const char *s);
int	PIL_get_by_value(PIL_PFILE *fp, const char *name, PIL_VALUE *v, int type);
#ifndef __CINT__
EXPSYM
#endif
int	PIL_get_bool(PIL_PFILE *fp, const char *name, int *result);
#ifndef __CINT__
EXPSYM
#endif
int	PIL_get_int(PIL_PFILE *fp, const char *name, int *result);
#ifndef __CINT__
EXPSYM
#endif
int	PIL_get_real(PIL_PFILE *fp, const char *name, double *result);
#ifndef __CINT__
EXPSYM
#endif
int	PIL_get_string(PIL_PFILE *fp, const char *name, char *result);
#ifndef __CINT__
EXPSYM
#endif
int	PIL_get_fname(PIL_PFILE *fp, const char *name, char *result);
#ifndef __CINT__
EXPSYM
#endif
int	PIL_get_dol(PIL_PFILE *fp, const char *name, char *result);

/* SCREW 1145: Added encapsulations for OS-dependent functions. */
		/* OS-dependent functions. */

int	PIL_get_file_size(const char *fp, size_t *fs);
int	PIL_get_mod_time(const char *name, time_t *modtime);
int	PIL_file_exists(const char *name);
int	PIL_file_readable(const char *name);
int	PIL_file_writable(const char *name);
int	PIL_lock_pfile(PIL_PFILE *fp, int mode);
int	PIL_unlock_pfile(PIL_PFILE *fp);
int	PIL_truncate_file(PIL_PFILE *fp);
void	PIL_sleep(int t);
char	*PIL_uniq_fname(const char *orig);

		/* official API stuff */

#ifndef __CINT__
EXPSYM
#endif
int	PILGetBool(const char *name, int *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetInt(const char *name, int *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetReal(const char *name, double *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetReal4(const char *name, float *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetString(const char *name, char *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetFname(const char *name, char *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetDOL(const char *name, char *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetAsString(const char *name, char *result);
#ifndef __CINT__
EXPSYM
#endif
int     PILPutBool(const char *name, int b);
#ifndef __CINT__
EXPSYM
#endif
int     PILPutInt(const char *name, int i);
#ifndef __CINT__
EXPSYM
#endif
int     PILPutReal(const char *name, double d);
#ifndef __CINT__
EXPSYM
#endif
int     PILPutString(const char *name, const char *s);
#ifndef __CINT__
EXPSYM
#endif
int     PILPutFname(const char *name, const char *s);
#ifndef __CINT__
EXPSYM
#endif
int	PILSetModuleName(const char *name);
#ifndef __CINT__
EXPSYM
#endif
int	PILSetModuleVersion(const char *version);
#ifndef __CINT__
EXPSYM
#endif
int	PILInit(int argc, char **argv);
#ifndef __CINT__
EXPSYM
#endif
int	PILClose(int status);
int	PILF90AddArg(const char *arg);
int	PILF90DeleteArgs(void);
int	PILF90Init(void);
int	PILF90Close(int status);
int	PILGetParFilename(char **fname);
int	PILLockPFile(void);
int	PILUnlockPFile(void);
#ifndef __CINT__
EXPSYM
#endif
int	PILReloadParameters(void);
#ifndef __CINT__
EXPSYM
#endif
int	PILFlushParameters(void);
#ifndef __CINT__
EXPSYM
#endif
int	PILOverrideQueryMode(int newmode);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetQueryMode(void);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetIntVector(const char *name, int nelem, int *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILCheckIntVector(const char *buf, int nelem, int *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetIntVarVector(const char *name, int *nelem, int *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILCheckIntVarVector(const char *buf, int *nelem, int *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetRealVector(const char *name, int nelem, double *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILCheckRealVector(const char *buf, int nelem, double *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetRealVarVector(const char *name, int *nelem, double *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILCheckRealVarVector(const char *buf, int *nelem, double *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetReal4Vector(const char *name, int nelem, float *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILCheckReal4Vector(const char *buf, int nelem, float *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetReal4VarVector(const char *name, int *nelem, float *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILCheckReal4VarVector(const char *buf, int *nelem, float *result);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetNumParameters(int *parnum);
#ifndef __CINT__
EXPSYM
#endif
int	PILGetParameter(int idx, PIL_PARAM *pp, int *minmaxok, PIL_VALUE *vmin, PIL_VALUE *vmax);
#ifndef __CINT__
EXPSYM
#endif
int	PILVerifyCmdLine(void);
#ifndef __CINT__
EXPSYM
#endif
int	PILSetRootNameFunction(int (*func)(char *s));
#ifndef __CINT__
EXPSYM
#endif
int	PILSetReadlinePromptMode(int mode);
#ifndef __CINT__
EXPSYM
#endif
int	PILSetCmndArgMode(int mode);
#ifndef __CINT__
EXPSYM
#endif
int	PILSetLoggerFunction(int (*func)(char *s));
#ifndef __CINT__
EXPSYM
#endif
/* SCREW 1498: Add capability for client code to provide a custom function with which to check
   whether a file exists. */
int	PILSetFileAccessFunction(int (*func)(const char *file_name, const char *open_mode));
#ifndef __CINT__
EXPSYM
#endif
/* SCREW 1496: Allow for reprompting. */
int	PILSetReprompt(const char *par_name, int reprompt);

#if defined(__cplusplus) && ! defined(__CINT__)
}
#endif
		/* local include files */

#include "pil_error.h"

#endif
