/**************************************************************************************

	I N T E G R A L   S C I E N C E   D A T A   C E N T E R

	P A R A M E T E R   I N T E R F A C E   L I B R A R Y

Copyright:	(c) 1998-2002 ISDC http://isdc.unige.ch
File:		pil_sys.c
Description:	Parameter Interface Library - OS-dependent functions
Authors:	James Peachey, peachey@lheamail.gsfc.nasa.gov (JP)
History:	14-Apr-03 (JP) : separating OS-specific functions

**************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef WIN32
#include <windows.h>
static const DWORD sFILE_NOT_FOUND = 0xFFFFFFFF;
#else

#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>

#endif

#include "pil.h"				/* include our stuff */

/* SCREW 1498: Pointer to a function which determines whether or not the file name may
   be opened with the given mode. */
static int (*sPILFileAccessFunction)(const char *file_name, const char *open_mode) = 0;

/**************************************************************************************
Function:	PIL_get_file_size
Description:	uses operating system call(s) to get file size.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)	name of file whose size to get
fs (Out)	size of the file

Return code:	standard PIL error code or PIL_OK (zero)
**************************************************************************************/

int	PIL_get_file_size(const char *name, size_t *fs)
 { int r;

   if (NULL == name) return(PIL_NUL_PTR);
   if (NULL == fs) return(PIL_NUL_PTR);
   if (0 == name[0]) return(PIL_NO_FILE);	/* name must be specified */

   r = PIL_OK;

#ifdef WIN32
   { HANDLE fh;
     WIN32_FIND_DATA fd;

     fh = FindFirstFile(name, &fd);
     if (INVALID_HANDLE_VALUE == fh) { r = PIL_NO_FILE; }
     else
      /* fd contains a 64-bit file size. Assume all par files will
         be less than 4 Gb long, and take the low part. */
      { *fs = fd.nFileSizeLow;
        FindClose(fh);
      }
   }
#else
   { struct stat fst;

     if (0 == stat(name, &fst)) { *fs = (size_t)fst.st_size; }
     else { r = PIL_NO_FILE; }
   }
#endif

   return(r);
 }


/**************************************************************************************
Function:	PIL_get_mod_time
Description:	uses operating system call(s) to get file modification time
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)	name of file
modtime (Out)	time of last modification

Return code:	standard PIL error code or PIL_OK (zero)
**************************************************************************************/

int	PIL_get_mod_time(const char *name, time_t *modtime)
 { int r;

   if (NULL == name) return(PIL_NUL_PTR);
   if (NULL == modtime) return(PIL_NUL_PTR);
   if (0 == name[0]) return(PIL_NO_FILE);	/* name must be specified */

   r = PIL_OK;

#ifdef WIN32
   { HANDLE fh;
     WIN32_FIND_DATA fd;

     fh = FindFirstFile(name, &fd);
     if (INVALID_HANDLE_VALUE == fh) { r = PIL_NO_FILE; }
     else
      /* fd contains a 64-bit time, in a struct containing two
         unsigned 32-bit numbers, ordered low then high. The time
         is nanoseconds, so ignore the low order number and take
         just the high order part. */
      { *modtime = fd.ftLastWriteTime.dwHighDateTime;
        FindClose(fh);
      }
   }
#else
   { struct stat statbuf;
     if (0 == stat(name, &statbuf)) { *modtime = statbuf.st_mtime; } 
     else { r = PIL_NO_FILE; }
   }
#endif

   return(r);
 }

/**************************************************************************************
Function:	PIL_file_testable
Description:	determines whether a file can be tested for existence etc.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)	name of file

Return code:	1 if file testable, 0 if file is not
**************************************************************************************/
static int	PIL_file_testable(const char *name)
   /* SCREW 1498: Check whether there is a custom file access checking function. If
      so, any file type should be considered "testable". */
 { if (0 != sPILFileAccessFunction) return 1;
   while (('\0' != *name) && (0 != isspace(*name))) ++name;
   if (!strncmp(name, "ftp:",   4) ||		/* various URLs */
       !strncmp(name, "http:",  5) ||
       !strncmp(name, "root:",  5) ||
       !strncmp(name, "mem:",   4) ||
       !strncmp(name, "shmem:", 6)) return(0);
   return(1);
 }

/**************************************************************************************
Function:	PIL_file_exists
Description:	uses operating system call(s) to find out if file exists
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)	name of file

Return code:	1 if file exists, 0 if file does not exist.
**************************************************************************************/

int	PIL_file_exists(const char *name)
   /* SCREW 1498: Check whether there is a custom file access checking function. */
 { if (0 != sPILFileAccessFunction) return (*sPILFileAccessFunction)(name, 0);
   if (NULL == name) return(PIL_NUL_PTR);
   if (0 == name[0]) return(PIL_NO_FILE);	/* name must be specified */
   if (!PIL_file_testable(name)) return 1;	/* assume non-testable file exists */
#ifdef WIN32
   if (0xFFFFFFFF == GetFileAttributes(name)) { return(0); }
   return(1);
#else
   return (!access(name, F_OK));
#endif
 }

/**************************************************************************************
Function:	PIL_file_readable
Description:	uses operating system call(s) to find out if file is readable
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)	name of file

Return code:	1 if file readable, 0 if not.
**************************************************************************************/

int	PIL_file_readable(const char *name)
   /* SCREW 1498: Check whether there is a custom file access checking function. */
 { if (0 != sPILFileAccessFunction) return (*sPILFileAccessFunction)(name, "r");
   if (NULL == name) return(PIL_NUL_PTR);
   if (0 == name[0]) return(PIL_NO_FILE);	/* name must be specified */
   if (!PIL_file_testable(name)) return 1;	/* assume non-testable file readable */
#ifdef WIN32
   { DWORD access;

     access = GetFileAttributes(name);
     if (sFILE_NOT_FOUND == access) { return(0); }
     /* More code here? */
   }
   return(1);
#else
   return (!access(name, R_OK));
#endif
 }

/**************************************************************************************
Function:	PIL_file_writable
Description:	uses operating system call(s) to find out if file is writable
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
name (In)	name of file

Return code:	1 if file writable, 0 if not.
**************************************************************************************/

int	PIL_file_writable(const char *name)
   /* SCREW 1498: Check whether there is a custom file access checking function. */
 { if (0 != sPILFileAccessFunction) return (*sPILFileAccessFunction)(name, "w");
   if (NULL == name) return(PIL_NUL_PTR);
   if (0 == name[0]) return(PIL_NO_FILE);	/* name must be specified */
   if (!PIL_file_testable(name)) return 0;	/* assume non-testable file not writable */
#ifdef WIN32
   { DWORD access;

     access = GetFileAttributes(name);
     if (sFILE_NOT_FOUND == access) { return(0); }
     else if (FILE_ATTRIBUTE_READONLY & access) { return(0); }
   }
   return(1);
#else
   return (!access(name, W_OK));
#endif
 }

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

int	PIL_lock_pfile(PIL_PFILE *fp, int mode)
#if(defined(WIN32) || !defined(PRE_SPR_1369))
 {
   /* File lock not yet implemented. The version in Windows requires
      a handle, not a FILE * */
   return(PIL_OK);
 }
#else
 { struct flock flk;

   if (NULL == fp) return(PIL_NUL_PTR);

   if (PIL_RDWRITE & mode)  { flk.l_type = F_WRLCK; }
   else  { flk.l_type = F_RDLCK; } 
   flk.l_whence = 0;
   flk.l_start = 0;
   flk.l_len = 0;
   if (-1 == fcntl(fileno(fp->fp), F_SETLKW, &flk))  return(PIL_LOCK_FAILED);

   return(PIL_OK);
 }
#endif


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

int	PIL_unlock_pfile(PIL_PFILE *fp)
#if(defined(WIN32) || !defined(PRE_SPR_1369))
 {
   /* File unlock not yet implemented. The version in Windows requires
      a handle, not a FILE * */
   return(PIL_OK);
 }
#else
 { struct flock flk;

   if (NULL == fp) return(PIL_NUL_PTR);

   flk.l_type = F_UNLCK;
   flk.l_whence = 0;
   flk.l_start = 0;
   flk.l_len = 0;
   if (-1 == fcntl(fileno(fp->fp), F_SETLKW, &flk))  return(PIL_LOCK_FAILED);

   return(PIL_OK);
 }
#endif


/**************************************************************************************
Function:	PIL_truncate_file
Description:	uses operating system call(s) to truncate file to size 1
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
file (In)	file pointer to truncate

Return code:	standard PIL error code or PIL_OK (zero)
**************************************************************************************/

int	PIL_truncate_file(PIL_PFILE *fp)
 { int r;

   if (NULL == fp) return(PIL_NUL_PTR);
   if (NULL == fp->fp) return(PIL_NUL_PTR);
#ifdef WIN32
   r = PIL_check_pfile(fp);
   if (PIL_OK == r)
     { fp->fp = freopen(fp->fname, "w", fp->fp);
       if (NULL == fp->fp) r = PIL_ERR_FWRITE;
     }
#else
   r = PIL_OK;
   if (ftruncate(fileno(fp->fp), 1)) { r = PIL_ERR_FWRITE; }
#endif
   return(r);
 }


/**************************************************************************************
Function:	PIL_sleep
Description:	uses operating system call(s) to sleep for the given number of seconds
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
t (In)		time in seconds to sleep

Return code:	standard PIL error code or PIL_OK (zero)
**************************************************************************************/

void	PIL_sleep(int t)
 {
#ifdef WIN32
   /* The function on Windows uses milliseconds */
   Sleep(1000 * t);
#else
   sleep(t);
#endif
 }


/**************************************************************************************
Function:	PIL_uniq_fname
Description:	generate a file name which is guaranteed unique,
		based on the input file name.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
orig (In)		original file name

Return:			tmpfile name based on input
**************************************************************************************/

#ifndef HOST_NAME_MAX
#define HOST_NAME_MAX 1023
#endif

char	*PIL_uniq_fname(const char *orig)
 { char *r = NULL;
   char hostname[HOST_NAME_MAX + 1] = "";
   char pid[32] = "";
   if (orig)
     {
#ifdef WIN32
       sprintf(pid, "%lu", GetCurrentProcessId());
#else
/* On Cygwin, gethostname tries to go out to the network and hangs if firewall installed. */
#ifndef __CYGWIN__
       gethostname(hostname, HOST_NAME_MAX);
#endif
       sprintf(pid, "%lu", (unsigned long) getpid());
#endif
       r = (char *) calloc(strlen(orig) + strlen(hostname) + strlen(pid) + 1, sizeof(char));
       if (NULL != r)
         { strcat(r, orig);
           strcat(r, hostname);
           strcat(r, pid);
           if (!strcmp(r, orig))
             { free(r);
               r = NULL;
             }
         }
     }

   return(r);
 }



/**************************************************************************************

Function:	PILSetFileAccessFunction
Description:	Supply PIL with a function to call to check for file accessibility.
		SCREW 1498.
ParameterName (I/O)	ParameterDescription:
---------------------------------------------------------------------------------------
func (In)		Pointer to the function

Return code:	standard PIL error code or PIL_OK (zero)

Notes:	PILSetFileAccessFunction instructs PIL to use specified function for
	logging error messages. Passing NULL disables logging. After PILInit
	the setup is as PILSetFileAccessFunction(NULL) has been called.

**************************************************************************************/

int	PILSetFileAccessFunction(int (*func)(const char *file_name, const char *open_mode))
 {
   sPILFileAccessFunction = func;
   return(PIL_OK);
 }
