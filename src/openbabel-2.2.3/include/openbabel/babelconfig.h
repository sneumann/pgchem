/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.in by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* Used by other #define statements for symbol exporting
 Use symbol hiding in GCC 4.0 and later where available */
#ifdef HAVE_GCC_VISIBILITY
  #define OB_EXPORT __attribute__ ((visibility("default")))
  #define OB_IMPORT __attribute__ ((visibility("default")))
  #define OB_HIDDEN __attribute__ ((visibility("hidden")))
#elif defined(WIN32)
  #define OB_EXPORT __declspec(dllexport)
  #define OB_IMPORT __declspec(dllimport)
  #define OB_HIDDEN
#else
  #define OB_EXPORT
  #define OB_IMPORT
  #define OB_HIDDEN
#endif


/* Where the data files are located */
#define BABEL_DATADIR "c:/tigress/obdata"

/* The version of Open Babel */
#define BABEL_VERSION "2.2.3"

/* Used to export symbols for DLL / shared library builds */
#if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define EXTERN OB_IMPORT extern
#else
  #define EXTERN OB_EXPORT extern
#endif


/* define if the Boost library is available */
#define HAVE_BOOST /**/

/* Define to 1 if the system has the type `clock_t'. */
/* #undef HAVE_CLOCK_T */

/* Define to 1 if you have the <conio.h> header file. */
#define HAVE_CONIO_H 1

/* Define to 1 if you have the <ctype.h> header file. */
#define HAVE_CTYPE_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
/* #undef HAVE_DLFCN_H */

/* Define to 1 if you have the <fstream> header file. */
#define HAVE_FSTREAM 1

/* Define to 1 if you have the <fstream.h> header file. */
/* #undef HAVE_FSTREAM_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <iostream> header file. */
#define HAVE_IOSTREAM 1

/* Define to 1 if you have the <iostream.h> header file. */
/* #undef HAVE_IOSTREAM_H */

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if you have the `z' library (-lz). */
#define HAVE_LIBZ 1

/* Define to 1 if you have the <locale.h> header file. */
#define HAVE_LOCALE_H 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the <regex.h> header file. */
/* #undef HAVE_REGEX_H */

/* Define to 1 if you have the `rint' function. */
#define HAVE_RINT 1

/* Define to 1 if you have the <rpc/types.h> header file. */
/* #undef HAVE_RPC_TYPES_H */

/* Define to 1 if you have the <rpc/xdr.h> header file. */
/* #undef HAVE_RPC_XDR_H */

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define to 1 if you have the `sranddev' function. */
/* #undef HAVE_SRANDDEV */

/* Define to 1 if you have the <sstream> header file. */
#define HAVE_SSTREAM 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strncasecmp' function. */
#define HAVE_STRNCASECMP 1

/* Define to 1 if you have the <strstream> header file. */
#define HAVE_STRSTREAM 1

/* Define to 1 if you have the <strstream.h> header file. */
/* #undef HAVE_STRSTREAM_H */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the <tr1/memory> header file. */
#define HAVE_TR1_MEMORY 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the `uselocale' function. */
/* #undef HAVE_USELOCALE */

/* Define to 1 if you have the <xlocale.h> header file. */
/* #undef HAVE_XLOCALE_H */

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#define LT_OBJDIR ".libs/"

/* The file extension used for shared modules */
#define MODULE_EXTENSION ".dll"

/* Used to export symbols for DLL / shared library builds */
#if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBAPI OB_IMPORT
#else
  #define OBAPI OB_EXPORT
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBCOMMON OB_IMPORT
#else
  #define OBCOMMON OB_EXPORT
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBCONV OB_IMPORT
#else
  #define OBCONV OB_EXPORT
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBERROR OB_IMPORT
#else
  #define OBERROR OB_EXPORT
#endif


/* Used to export symbols for DLL / shared library builds */
#if defined(USING_OBDLL) // e.g. in src/main.cpp
  #define OBFPRT OB_IMPORT
#else
  #define OBFPRT OB_EXPORT
#endif


/* Define to the address where bug reports for this package should be sent. */

/* Define to the full name of this package. */

/* Define to the full name and version of this package. */

/* Define to the one symbol short name of this package. */

/* Define to the version of this package. */

/* set if scandir needs a const */
#define SCANDIR_CONST const

/* set if scandir needs a const */
#define SCANDIR_T (int (*)(const dirent *))


#if !HAVE_SNPRINTF
extern "C" int snprintf( char *, size_t, const char *, /* args */ ...);
#endif


/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you can safely include both <sys/time.h> and <time.h>. */
#define TIME_WITH_SYS_TIME 1

/* Using BOOST for shared_pointer */
/* #undef USE_BOOST */

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
