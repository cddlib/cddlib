/*************************************************************************

        Copyright 1986 through 2000 by Wolfram Research Inc.
        All rights reserved

*************************************************************************/

#ifndef _MATHLINK_H
#define _MATHLINK_H

#if __BORLANDC__ && ! __BCPLUSPLUS__
#pragma warn -stu
#endif

#ifndef _MLVERS_H
#define _MLVERS_H

#ifndef _MLPLATFM_H
#define _MLPLATFM_H

#if ! MACINTOSH_MATHLINK && ! WINDOWS_MATHLINK && ! UNIX_MATHLINK && ! OS2_MATHLINK
#	if 0
#	if __BEOS__
#		define BE_MATHLINK 1
#	elif macintosh || Macintosh || THINK_C || defined(_MAC) || defined(__MRC__)
#		define MACINTOSH_MATHLINK 1
#	elif defined(WIN16) || defined(_WIN16)
#		define WINDOWS_MATHLINK 1
#	elif defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#		define WINDOWS_MATHLINK 1
#	elif unix || __unix || __unix__ ||_AIX
#		define UNIX_MATHLINK 1
#	endif
#	endif
#	define UNIX_MATHLINK 1
#endif

#if MACINTOSH_MATHLINK
#	if (powerc || __powerc || __powerc__)
#		define POWERMACINTOSH_MATHLINK 1
#	else
#		define M68KMACINTOSH_MATHLINK 1
#		if defined(__CFM68K__)
#			define CFM68K_MACINTOSH_MATHLINK 1
#		else
#			define CLASSIC68K_MACINTOSH_MATHLINK 1
#		endif
#	endif
#elif WINDOWS_MATHLINK
#	if defined(WIN32) || defined(__WIN32__) || defined(__NT__) || defined(_WIN32)
#		define WIN32_MATHLINK 1
#		if( _M_IX86 || __i386 || __i386__ || i386)
#			define I86_WIN32_MATHLINK 1
#		elif _M_ALPHA || __alpha || __alpha__ || alpha
#			define ALPHA_WIN32_MATHLINK 1
#		else
#		endif
#	else
#		define WIN16_MATHLINK 1
#	endif
#elif UNIX_MATHLINK
#	if (__sun || __sun__ || sun) && !defined(SUN_MATHLINK)
#		define SUN_MATHLINK 1
#		if __SVR4 || __svr4__
#			define SOLARIS_MATHLINK 1
#		else
#			define SUNOS_MATHLINK 1
#		endif
#		if __sparc || __sparc__ || sparc
#			define SPARC_SUN_MATHLINK 1
#		elif __i386 || __i386__ || i386
#			define I86_SUN_MATHLINK 1
#		else
			unknown platform
#		endif
#       elif (__MACH || __MACH__ || MACH) && !defined(DARWIN_MATHLINK)
#               define DARWIN_MATHLINK 1
#               if __ppc || __ppc__ || ppc
#                       define PPC_DARWIN_MATHLINK 1
#               else
                        not yet implemented
#               endif
#	elif (__linux || __linux__ || linux) && !defined(LINUX_MATHLINK)
#		define LINUX_MATHLINK 1
#		if __i386 || __i386__ || i386
#			define I86_LINUX_MATHLINK 1
#		elif __PPC || __PPC__ || PPC
#			define PPC_LINUX_MATHLINK 1
#		elif __alpha || __alpha__ || alpha
#			define AXP_LINUX_MATHLINK 1
#		else
			not yet implemented
#		endif
#	elif (__osf || __osf__ || osf || OSF1) && !defined(DIGITAL_MATHLINK)
#		define DIGITAL_MATHLINK 1
#		if __alpha || __alpha__ || alpha
#			define ALPHA_DIGITAL_MATHLINK 1
#		else
			unknown platform
#		endif
#	elif (_AIX || _IBMR2 || __xlC__) && !defined(AIX_MATHLINK)
#		define AIX_MATHLINK 1
#	elif (__sgi || __sgi__ || sgi || mips) && !defined(IRIX_MATHLINK)
#		define IRIX_MATHLINK 1
#	elif (hpux || __hpux) && !defined(HPUX_MATHLINK)
#		define HPUX_MATHLINK 1
#	elif (M_I386 || _SCO_DS || SCO) && !defined(SCO_MATHLINK)
#		define SCO_MATHLINK 1
#	elif (__NetBSD__) && !defined(NETBSD_MATHLINK)
#		define NETBSD_MATHLINK 1
#	elif (__FreeBSD__) && !defined(FREEBSD_MATHLINK)
#		define FREEBSD_MATHLINK 1
#	elif (bsdi || __bsdi__) && !defined(BSDI_MATHLINK)
#		define BSDI_MATHLINK 1
#	else
#	endif
#else
#	if defined(__amigaos__) || defined(AMIGA)
#		define AMIGA_MATHLINK 1
#	endif
#	if __BEOS__
#		define BE_MATHLINK 1
#	endif
#endif



#ifndef NO_GLOBAL_DATA
#	if M68KMACINTOSH_MATHLINK
#		define NO_GLOBAL_DATA 1
#	else
#		define NO_GLOBAL_DATA 0
#	endif
#endif

#if WINDOWS_MATHLINK || __i386 || __i386__ || i386 || _M_IX86
#	define LITTLEENDIAN_NUMERIC_TYPES 1
#else
#	define BIGENDIAN_NUMERIC_TYPES 1
#endif

#endif /* _MLPLATFM_H */

#ifndef MLVERSION
#	define MLVERSION 3
#endif

#if !OLD_VERSIONING


/*
 * MathLink adopts a simple versioning strategy that can be adapted to many
 * compile-time and run-time environments.  In particular, it is amenable to
 * the various shared library facilities in use.  (Although certain of these
 * facilities provide more sophisticated mechanisms than are required by the
 * following simple strategy.)
 * 
 * MathLink evolves by improving its implementation and by improving its
 * interface.  The values of MLREVISION or MLINTERFACE defined here are 
 * incremented whenever an improvement is made and released.
 * 
 * MLREVISION is the current revision number. It is incremented every time 
 * a change is made to the source and MathLink is rebuilt and distributed
 * on any platform.  (Bug fixes, optimizations, or other improvements
 * transparent to the interface increment only this number.)
 * 
 * MLINTERFACE is a name for a documented interface to MathLink.  This
 * number is incremented whenever a named constant or function is added,
 * removed, or its behavior is changed in a way that could break existing
 * correct* client programs.  It is expected that the interface to MathLink
 * is improved over time so that implemenations of higher numbered 
 * interfaces are more complete or more convenient to use for writing
 * effective client programs.  In particular, a specific interface provides
 * all the useful functionality of an earlier interface.
 * 
 *     *(It is possible that an incorrect MathLink program still works
 *     because it relies on some undocumented detail of a particular
 *     revision.  It may not always be possible to change the interface
 *     number when such a detail changes.  For example, one program may
 *     be relying on a bug in MathLink that a great many other programs
 *     need fixed.  In this case, we would likely choose to potentially
 *     break the incorrect program in order to fix the correct programs
 *     by incrementing the revision number leaving the interface number
 *     unchanged.  It is possible to bind to a particular revision of a
 *     MathLink interface if that is important for some programs.  One
 *     could use a statically linked version of the library, make use of
 *     the search algorithm used by the runtime loader, or dynamically
 *     load the MathLink library manually.)
 * 
 * 
 * If a distributed MathLink implmentation were labeled with its revision 
 * and interface numbers in dotted notation so that, say, ML.1.6 means the
 * sixth revision of interface one, then the following may represent the
 * distribution history of MathLink.
 * 
 *     first distribution
 *         ML.1.5   (Perhaps earlier revisions were never
 *                   distributed for this platform.)
 * 
 *     second distribution
 *         ML.1.6   (Bug fixes or other improvements were
 *                   made that don't affect the interface.)
 * 
 *     third distribution
 *         ML.2.7   (Perhaps some new functions were added.)
 *         
 *         ML.1.7   (And improvements were made that don't
 *                   affect the old interface.)
 * 
 *     fourth distribution
 *         ML.3.8   (Perhaps the return values of an existing
 *                   function changed.)
 *         ML.2.8   (Revision 8 also adds improvements transparent
 *                   to interface 2.)
 *         ML.1.7   (Clients of interface 1 see no improvements
 *                   in this eighth revision.)
 * 
 * Note that the distribution history may not be the same on different
 * platforms.  But revision numbers represent a named body of source code
 * across all platforms.
 * 
 * The mechanism for deploying this strategy differs between platforms
 * because of differing platform-specific facilities and conventions.
 * The interface and revision numbers may form part of the filename of
 * the MathLink library, or they may not.  This information is always
 * available in some conventional form so that it is easy and natural for
 * client programs to bind with and use the best available implementation
 * of a particular MathLink interface.  The details are described in the
 * MathLink Developer's Guide for each platform.
 */

#define MLREVISION 9

#define MLAPI1REVISION 1 /* the first revision to support interface 1 */
#define MLAPI2REVISION 6 /* the first revision to support interface 2 */


#ifndef MLINTERFACE
#	define MLINTERFACE 2
#	define MLAPIREVISION MLAPI2REVISION
    	/*
		 * Interface 2 adds the following exported functions:
		 *      MLGetBinaryNumberArray0
		 *      MLTransfer0
		 *      MLNextCharacter0
		 * And, for WINDOWS_MATHLINK, some constants in "mlntypes.h"
		 * were changed in a way that causes MLGetRawType to return
		 * different values.
		 *
		 *      MLPutNullSequence and MLEGETENDEXPR pushed to interface 3
		 */
#else
#	if MLINTERFACE == 1
#		define MLAPIREVISION MLAPI1REVISION
#	elif MLINTERFACE == 2
#		define MLAPIREVISION MLAPI2REVISION
#	else
/* syntax error */ )
#	endif
#endif


/* It may be possible for an implementation of one MathLink interface to
 * fully support an earlier interface.  MLNewParameters() may succeed when
 * passed an interface number less than the value of MLAPIREVISION when the
 * library was built.  This would happen, if the newer interface is a proper
 * superset of the older interface, or if the implementation can adjust its
 * behavior at runtime to conform to the older requested interface.
 */

#ifndef MLOLDDEFINITION
#	if WINDOWS_MATHLINK
#		if MLINTERFACE == 1
#			define MLOLDDEFINITION MLAPI1REVISION
#		elif MLINTERFACE == 2
#			define MLOLDDEFINITION MLAPI2REVISION
#		else
/* syntax error */ )
#		endif
#	else
#		define MLOLDDEFINITION MLAPI1REVISION
#	endif
#endif


#if 0
MLParameters s;
MLNewParameters( s, MLREVISION, MLAPIREVISION); or MLNewParameters( s, 0, MLAPIREVISION);
MLSetAllocParameter( s, allocator, deallocator);
MLIntialize(s);
#endif








#else
/* syntax error */ )
#endif

#endif /* _MLVERS_H */


#ifndef ML_EXTERN_C

#if defined(__cplusplus)
#	define ML_C "C"
#	define ML_EXTERN_C extern "C" {
#	define ML_END_EXTERN_C }
#else
#	define ML_C
#	define ML_EXTERN_C
#	define ML_END_EXTERN_C
#endif

#endif



#if WINDOWS_MATHLINK && (MPREP_REVISION || !defined(APIENTRY) || !defined(FAR))

#if defined(WIN32_LEAN_AND_MEAN) && defined(WIN32_EXTRA_LEAN)
#	include <windows.h>
#elif defined( WIN32_LEAN_AND_MEAN)
#	define WIN32_EXTRA_LEAN
#	include <windows.h>
#	undef WIN32_EXTRA_LEAN
#elif defined( WIN32_EXTRA_LEAN)
#	define WIN32_LEAN_AND_MEAN
#	include <windows.h>
#	undef WIN32_LEAN_AND_MEAN
#else
#	define WIN32_EXTRA_LEAN
#	define WIN32_LEAN_AND_MEAN
#	include <windows.h>
#	undef WIN32_EXTRA_LEAN
#	undef WIN32_LEAN_AND_MEAN
#endif

#endif

#ifndef _MLCFM_H
#define _MLCFM_H


#ifndef GENERATINGCFM
#	ifdef USESROUTINEDESCRIPTORS
#		define GENERATINGCFM USESROUTINEDESCRIPTORS
#	elif MACINTOSH_MATHLINK
#		include <ConditionalMacros.h>
#		ifndef GENERATINGCFM
#			define GENERATINGCFM USESROUTINEDESCRIPTORS
#		endif
#	else
#		define GENERATINGCFM 0
#	endif
#endif


#if MACINTOSH_MATHLINK
#	include <MixedMode.h>
#elif DARWIN_MATHLINK && defined(__CARBON__)
#	if defined(GENERATINGCFM)
#		undef GENERATINGCFM
#		define GENERATINGCFM 	0
#	endif /* defined(GENERATINGCFM) */
#	if defined(GENERATING68K)
#		undef GENERATING68K
#		define GENERATING68K	0
#	endif /* defined(GENERATING68K) */
#else
	enum {
		kPascalStackBased = 0,
		kCStackBased = 0,
		kThinkCStackBased = 0
	};
#	define SIZE_CODE(size) (0)
#	define RESULT_SIZE(sizeCode) (0)
#	define STACK_ROUTINE_PARAMETER(whichParam, sizeCode) (0)
#endif


#endif /* _MLCFM_H */


#ifdef __CFM68K__
#pragma import on
#endif


#ifndef _MLDEVICE_H
#define _MLDEVICE_H


#ifndef P

#  ifndef MLPROTOTYPES
#    define MLPROTOTYPES 1
#  endif

#  if MLPROTOTYPES || __STDC__ || defined(__cplusplus) || ! UNIX_MATHLINK
#    define P(s) s
#	 undef MLPROTOTYPES
#	 define MLPROTOTYPES 1
#  else
#    define P(s) ()
#	 undef MLPROTOTYPES
#	 define MLPROTOTYPES 0
#  endif
#endif
#ifndef _MLFAR_H
#define _MLFAR_H

#ifndef FAR

#if WINDOWS_MATHLINK
#	ifndef FAR
/* syntax error */ )
#	endif
#else
#	define FAR
#endif


#endif

/* //rename this file mlfarhuge.h */

#ifndef MLHUGE
#  if WINDOWS_MATHLINK && ! WIN32_MATHLINK
#    define MLHUGE huge
#  else
#    define MLHUGE
#  endif
#endif

#endif /* _MLFAR_H */

#ifndef _MLTYPES_H
#define _MLTYPES_H


#if WINDOWS_MATHLINK
#	ifndef	APIENTRY
#		define APIENTRY far pascal
#	endif
#	ifndef CALLBACK
#		define CALLBACK APIENTRY
#	endif
#	if WIN32_MATHLINK
 /* try this #define MLEXPORT __declspec(dllexport) */
#		define MLEXPORT
#	else
#		define MLEXPORT __export
#	endif
#	define MLCB APIENTRY MLEXPORT
#	define MLAPI APIENTRY

#elif OS2_MATHLINK
#	include <os2def.h>
#	define MLEXPORT
#	define MLCB APIENTRY
#	define MLAPI APIENTRY
#elif CLASSIC68K_MACINTOSH_MATHLINK
#	define MLAPI pascal
#	define MLEXPORT
#	if defined(__MWERKS__)
#		if !__fourbyteints__
#			define __uint_ct__ unsigned long
#			define __int_ct__ long
#		endif
#	elif defined(THINK_C) || defined(SYMANTEC_C) || defined(SYMANTEC_CPLUS)
#		if !__option(int_4)
#			define __uint_ct__ unsigned long
#			define __int_ct__ long
#		endif
#	endif
#else
#	define MLCB
#	define MLAPI
#	define MLEXPORT
#endif

#define MLAPI_ MLAPI


#ifndef MLDEFN
#	define MLDEFN( rtype, name, params) extern rtype MLAPI MLEXPORT name params
#endif
#ifndef MLDECL
#	define MLDECL( rtype, name, params) extern rtype MLAPI name P(params)
#endif

#ifndef ML_DEFN
#	define ML_DEFN( rtype, name, params) extern rtype MLAPI_ MLEXPORT name params
#endif
#ifndef ML_DECL
#	define ML_DECL( rtype, name, params) extern ML_C rtype MLAPI_ name P(params)
#endif



#if MACINTOSH_MATHLINK

#ifndef MLCBPROC
#	define MLCBPROC( rtype, name, params) typedef pascal rtype (* name) P(params)
#endif
#ifndef MLCBDECL
#	define MLCBDECL( rtype, name, params) extern pascal rtype name P(params)
#endif
#ifndef MLCBDEFN
#	define MLCBDEFN( rtype, name, params) extern pascal rtype name params
#endif

#elif OS2_MATHLINK

#ifndef MLCBPROC
#	define MLCBPROC( rtype, name, params) typedef rtype (* MLCB name) P(params)
#endif
#ifndef MLCBDECL
#	define MLCBDECL( rtype, name, params) extern rtype MLCB name P(params)
#endif
#ifndef MLCBDEFN
#	define MLCBDEFN( rtype, name, params) extern rtype MLCB name params
#endif

#else

#ifndef MLCBPROC
#	define MLCBPROC( rtype, name, params) typedef rtype (MLCB * name) P(params)
#endif
#ifndef MLCBDECL
#	define MLCBDECL( rtype, name, params) extern rtype MLCB name P(params)
#endif
#ifndef MLCBDEFN
#	define MLCBDEFN( rtype, name, params) extern rtype MLCB name params
#endif

#endif




/* move into mlalert.h */
#ifndef MLDPROC
#	define MLDPROC MLCBPROC
#endif
#ifndef MLDDECL
#	define MLDDECL MLCBDECL
#endif
#ifndef MLDDEFN
#	define MLDDEFN MLCBDEFN
#endif




/* move into ml3state.h or mlstrenv.h */
#ifndef MLTPROC
#	define MLTPROC MLCBPROC
#endif
#ifndef MLTDECL
#	define MLTDECL MLCBDECL
#endif
#ifndef MLTDEFN
#	define MLTDEFN MLCBDEFN
#endif


/* move into mlnumenv.h */
#ifndef MLNPROC
#	define MLNPROC MLCBPROC
#endif
#ifndef MLNDECL
#	define MLNDECL MLCBDECL
#endif
#ifndef MLNDEFN
#	define MLNDEFN MLCBDEFN
#endif


/* move into mlalloc.h */
#ifndef MLAPROC
#	define MLAPROC MLCBPROC
#endif
#ifndef MLADECL
#	define MLADECL MLCBDECL
#endif
#ifndef MLADEFN
#	define MLADEFN MLCBDEFN
#endif
#ifndef MLFPROC
#	define MLFPROC MLCBPROC
#endif
#ifndef MLFDECL
#	define MLFDECL MLCBDECL
#endif
#ifndef MLFDEFN
#	define MLFDEFN MLCBDEFN
#endif




/* move into mlstddev.h */
#ifndef MLYPROC
#	define MLYPROC MLCBPROC
#endif
#ifndef MLYDECL
#	define MLYDECL MLCBDECL
#endif
#ifndef MLYDEFN
#	define MLYDEFN MLCBDEFN
#endif
#ifndef MLMPROC
#	define MLMPROC MLCBPROC
#endif
#ifndef MLMDECL
#	define MLMDECL MLCBDECL
#endif
#ifndef MLMDEFN
#	define MLMDEFN MLCBDEFN
#endif


/* move into mlmake.h */
#ifndef MLUPROC
#	define MLUPROC MLCBPROC
#endif
#ifndef MLUDECL
#	define MLUDECL MLCBDECL
#endif
#ifndef MLUDEFN
#	define MLUDEFN MLCBDEFN
#endif


/* move into mlmake.h */
#ifndef MLBPROC
#	define MLBPROC MLCBPROC
#endif
#ifndef MLBDECL
#	define MLBDECL MLCBDECL
#endif
#ifndef MLBDEFN
#	define MLBDEFN MLCBDEFN
#endif

#ifndef MLDMPROC
#	define MLDMPROC MLCBPROC
#endif
#ifndef MLDMDECL
#	define MLDMDECL MLCBDECL
#endif
#ifndef MLDMDEFN
#	define MLDMDEFN MLCBDEFN
#endif


#ifndef __uint_ct__
#define __uint_ct__ unsigned int
#endif
#ifndef __int_ct__
#define __int_ct__ int
#endif


typedef unsigned char        uchar_ct;
typedef uchar_ct       FAR * ucharp_ct;
typedef ucharp_ct      FAR * ucharpp_ct;
typedef ucharpp_ct     FAR * ucharppp_ct;
typedef unsigned short       ushort_ct;
typedef ushort_ct      FAR * ushortp_ct;
typedef ushortp_ct     FAR * ushortpp_ct;
typedef ushortpp_ct    FAR * ushortppp_ct;
typedef __uint_ct__          uint_ct;
typedef __int_ct__           int_ct;
typedef void           FAR * voidp_ct;
typedef voidp_ct       FAR * voidpp_ct;
typedef char           FAR * charp_ct;
typedef charp_ct       FAR * charpp_ct;
typedef charpp_ct      FAR * charppp_ct;
typedef long           FAR * longp_ct;
typedef longp_ct       FAR * longpp_ct;
typedef unsigned long        ulong_ct;
typedef ulong_ct       FAR * ulongp_ct;




#ifndef MLCONST
#	if MLPROTOTYPES
#		define MLCONST const
#	else
#		define MLCONST
#	endif
#endif

typedef MLCONST unsigned short FAR * kushortp_ct;
typedef MLCONST unsigned short FAR * FAR * kushortpp_ct;
typedef MLCONST unsigned char FAR * kucharp_ct;
typedef MLCONST unsigned char FAR * FAR * kucharpp_ct;
typedef MLCONST char FAR * kcharp_ct;
typedef MLCONST char FAR * FAR * kcharpp_ct;
typedef MLCONST void FAR * kvoidp_ct;


typedef void FAR * MLPointer;

#ifndef __MLENV__
	typedef struct ml_environment FAR *MLENV;
	typedef MLENV MLEnvironment;
#	define __MLENV__
#endif

#ifndef __MLINK__
	typedef struct MLink FAR *MLINK;
#	define __MLINK__
#endif

#ifndef __MLMARK__
	typedef struct MLinkMark FAR *MLMARK;
	typedef MLMARK MLINKMark;
#	define __MLMARK__
#endif

#ifndef __mlapi_token__
#define __mlapi_token__ int_ct
#endif
typedef __mlapi_token__   mlapi_token;


typedef unsigned long      mlapi__token;
typedef mlapi__token FAR * mlapi__tokenp;

#ifndef __mlapi_packet__
#define __mlapi_packet__ int_ct
#endif
typedef __mlapi_packet__  mlapi_packet;


typedef long mlapi_error;
typedef long mlapi__error;

typedef long      long_st;
typedef longp_ct  longp_st;
typedef longpp_ct longpp_st;

typedef long long_et;


#ifndef __mlapi_result__
#define __mlapi_result__ int_ct
#endif
typedef __mlapi_result__ mlapi_result;


#define MLSUCCESS (1) /*bugcheck:  this stuff doesnt belong where it can be seen at MLAPI_ layer */
#define MLFAILURE (0)

ML_EXTERN_C

#if WINDOWS_MATHLINK
typedef int (CALLBACK *__MLProcPtr__)();
#else
typedef long (*__MLProcPtr__)();
#endif

ML_END_EXTERN_C

#endif /* _MLTYPES_H */


#if WINDOWS_MATHLINK
#	ifndef	APIENTRY
#		define	APIENTRY far pascal
#	endif
#	define MLBN APIENTRY /* bottleneck function: upper layer calls lower layer */
#else
#	define MLBN
#endif

#define BN MLBN



ML_EXTERN_C



typedef void FAR * dev_voidp;
typedef dev_voidp dev_type;
typedef dev_type FAR * dev_typep;
typedef long devproc_error;
typedef unsigned long devproc_selector;


#define MLDEV_WRITE_WINDOW  0
#define MLDEV_WRITE         1
#define MLDEV_HAS_DATA      2
#define MLDEV_READ          3
#define MLDEV_READ_COMPLETE 4
#define MLDEV_ACKNOWLEDGE   5

#define T_DEV_WRITE_WINDOW  MLDEV_WRITE_WINDOW
#define T_DEV_WRITE         MLDEV_WRITE
#define T_DEV_HAS_DATA      MLDEV_HAS_DATA
#define T_DEV_READ          MLDEV_READ
#define T_DEV_READ_COMPLETE MLDEV_READ_COMPLETE


#ifndef SCATTERED
#define SCATTERED 0
#undef NOT_SCATTERED
#define NOT_SCATTERED 1
#endif


#if powerc
#pragma options align=mac68k
#endif

typedef struct read_buf {
	unsigned short length;
	unsigned char* ptr;
} read_buf;

typedef read_buf FAR * read_bufp;
typedef read_bufp FAR * read_bufpp;

#if powerc
#pragma options align=reset
#endif



MLDMPROC( devproc_error, MLDeviceProcPtr, ( dev_type dev, devproc_selector selector, dev_voidp p1, dev_voidp p2));
MLDMDECL( devproc_error, MLDeviceMain, ( dev_type dev, devproc_selector selector, dev_voidp p1, dev_voidp p2));

enum {
	uppMLDeviceProcInfo = kPascalStackBased
		 | RESULT_SIZE(SIZE_CODE(sizeof(devproc_error)))
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(dev_type)))
		 | STACK_ROUTINE_PARAMETER(2, SIZE_CODE(sizeof(devproc_selector)))
		 | STACK_ROUTINE_PARAMETER(3, SIZE_CODE(sizeof(dev_voidp)))
		 | STACK_ROUTINE_PARAMETER(4, SIZE_CODE(sizeof(dev_voidp)))
};

#if GENERATINGCFM

	typedef UniversalProcPtr MLDeviceUPP;
#	define CallMLDeviceProc(userRoutine, thing, selector, p1, p2) \
		CallUniversalProc((userRoutine), uppMLDeviceProcInfo, (thing), (selector), (p1), (p2))
#	define NewMLDeviceProc(userRoutine) \
		NewRoutineDescriptor((ProcPtr)(userRoutine), uppMLDeviceProcInfo, GetCurrentArchitecture())

#else

	typedef MLDeviceProcPtr MLDeviceUPP;
#	define CallMLDeviceProc(userRoutine, thing, selector, p1, p2) (*(userRoutine))((thing), (selector), (dev_voidp)(p1), (dev_voidp)(p2))
#	define NewMLDeviceProc(userRoutine) (userRoutine)

#endif

typedef MLDeviceUPP dev_main_type;
typedef dev_main_type FAR * dev_main_typep;

ML_END_EXTERN_C


#endif /* _MLDEVICE_H */


#ifndef _MLAPI_H
#define _MLAPI_H


ML_EXTERN_C

#ifndef _MLALLOC_H
#define _MLALLOC_H




MLAPROC( MLPointer, MLAllocatorProcPtr, (unsigned long));

enum {
	uppMLAllocatorProcInfo = kPascalStackBased
		 | RESULT_SIZE(SIZE_CODE(sizeof(MLPointer)))
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(unsigned long)))
};

#if GENERATINGCFM
	typedef UniversalProcPtr MLAllocatorUPP;
#	define CallMLAllocatorProc(userRoutine, size) \
		(MLPointer)CallUniversalProc((userRoutine), uppMLAllocatorProcInfo, (size))
#	define NewMLAllocatorProc(userRoutine) \
		NewRoutineDescriptor(MLAllocatorCast((userRoutine)), uppMLAllocatorProcInfo, GetCurrentArchitecture())
#else
	typedef MLAllocatorProcPtr MLAllocatorUPP;
#	define CallMLAllocatorProc(userRoutine, size) (*(userRoutine))((size))
#	define NewMLAllocatorProc(userRoutine) (userRoutine)
#endif




MLFPROC( void, MLDeallocatorProcPtr, (MLPointer));

enum {
	uppMLDeallocatorProcInfo = kPascalStackBased
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(MLPointer)))
};

#if GENERATINGCFM
	typedef UniversalProcPtr MLDeallocatorUPP;
#	define CallMLDeallocatorProc(userRoutine, p) \
		CallUniversalProc((userRoutine), uppMLDeallocatorProcInfo, (p))
#	define NewMLDeallocatorProc(userRoutine) \
		NewRoutineDescriptor(MLDeallocatorCast((userRoutine)), uppMLDeallocatorProcInfo, GetCurrentArchitecture())
#else
	typedef MLDeallocatorProcPtr MLDeallocatorUPP;
#	define CallMLDeallocatorProc(userRoutine, p) (*(userRoutine))((p))
#	define NewMLDeallocatorProc(userRoutine) (userRoutine)
#endif



#endif /* _MLALLOC_H */


/* explicitly not protected by _MLALLOC_H in case MLDECL is redefined for multiple inclusion */


/* just some type-safe casts */
MLDECL( __MLProcPtr__, MLAllocatorCast,   ( MLAllocatorProcPtr f));
MLDECL( __MLProcPtr__, MLDeallocatorCast, ( MLDeallocatorProcPtr f));

ML_END_EXTERN_C

typedef MLAllocatorUPP MLAllocator;
typedef MLAllocator FAR * MLAllocatorp;
#define MLCallAllocator CallMLAllocatorProc
#define MLNewAllocator NewMLAllocatorProc

typedef MLDeallocatorUPP MLDeallocator;
typedef MLDeallocator FAR * MLDeallocatorp;
#define MLCallDeallocator CallMLDeallocatorProc
#define MLNewDeallocator NewMLDeallocatorProc

#define MLallocator MLAllocator
#define MLdeallocator MLDeallocator

#endif /* _MLAPI_H */


#ifndef _MLNTYPES_H
#define _MLNTYPES_H


#ifndef _MLNUMENV_H
#define _MLNUMENV_H


/* mlne__s2 must convert empty strings to zero */



ML_EXTERN_C


#define REALBIT 4
#define REAL_MASK (1 << REALBIT)
#define XDRBIT 5
#define XDR_MASK (1 << XDRBIT)
#define BINARYBIT 7
#define BINARY_MASK (1 << BINARYBIT)
#define SIZEVARIANTBIT 6
#define SIZEVARIANT_MASK (1 << SIZEVARIANTBIT)



#define MLNE__IMPLIED_SIZE( tok, num_dispatch) ((tok) & XDR_MASK || !((tok) & SIZEVARIANT_MASK) \
		? (tok) & 0x08 ? (tok) & (0x0E + 2) : (1 << ((tok)>>1 & 0x03)) \
		: call_num_dispatch( (num_dispatch), MLNE__SIZESELECTOR((tok)), 0,0,0))

/* Range[-128, 127] */
#define	MLTK_8BIT_SIGNED_2sCOMPLEMENT_INTEGER                 160 /* ((unsigned char)'\240') */
/* Range[0, 255] */
#define	MLTK_8BIT_UNSIGNED_2sCOMPLEMENT_INTEGER               161 /* ((unsigned char)'\241') */
#define MLTK_8BIT_UNSIGNED_INTEGER MLTK_8BIT_UNSIGNED_2sCOMPLEMENT_INTEGER

/* Range[-32768, 32767] */
#define	MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER      162 /* ((unsigned char)'\242') */
/* Range[0, 65535] */
#define	MLTK_16BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER    163 /* ((unsigned char)'\243') */
#define	MLTK_16BIT_UNSIGNED_BIGENDIAN_INTEGER MLTK_16BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
/* Range[-2147483648, 2147483647] */
#define	MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER      164 /* ((unsigned char)'\244') */
/* Range[0, 4294967295] */
#define	MLTK_32BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER    165 /* ((unsigned char)'\245') */
#define	MLTK_32BIT_UNSIGNED_BIGENDIAN_INTEGER MLTK_32BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
/* Range[-9223372036854775808, 9223372036854775807] */
#define	MLTK_64BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER      166 /* ((unsigned char)'\246') */
/* Range[0, 18446744073709551615] */
#define	MLTK_64BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER    167 /* ((unsigned char)'\247') */
#define	MLTK_64BIT_UNSIGNED_BIGENDIAN_INTEGER MLTK_64BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER


/* Range[-32768, 32767] */
#define	MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER   226 /* ((unsigned char)'\342') */
/* Range[0, 65535] */
#define	MLTK_16BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER 227 /* ((unsigned char)'\343') */
#define	MLTK_16BIT_UNSIGNED_LITTLEENDIAN_INTEGER MLTK_16BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
/* Range[-2147483648, 2147483647] */
#define	MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER   228 /* ((unsigned char)'\344') */
/* Range[0, 4294967295] */
#define	MLTK_32BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER 229 /* ((unsigned char)'\345') */
#define	MLTK_32BIT_UNSIGNED_LITTLEENDIAN_INTEGER MLTK_32BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
/* Range[-9223372036854775808, 9223372036854775807] */
#define	MLTK_64BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER   230 /* ((unsigned char)'\346') */
/* Range[0, 18446744073709551615] */
#define	MLTK_64BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER 231 /* ((unsigned char)'\347') */
#define	MLTK_64BIT_UNSIGNED_LITTLEENDIAN_INTEGER MLTK_64BIT_UNSIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER

/* Interval[{-3.402823e+38, 3.402823e+38}] */
#define	MLTK_BIGENDIAN_IEEE754_SINGLE	                      180 /* ((unsigned char)'\264') */
/* Interval[{-1.79769313486232e+308, 1.79769313486232e+308}] */
#define	MLTK_BIGENDIAN_IEEE754_DOUBLE	                      182 /* ((unsigned char)'\266') */

/* Interval[{-3.402823e+38, 3.402823e+38}] */
#define	MLTK_LITTLEENDIAN_IEEE754_SINGLE	                  244 /* ((unsigned char)'\364') */
/* Interval[{-1.79769313486232e+308, 1.79769313486232e+308}] */
#define	MLTK_LITTLEENDIAN_IEEE754_DOUBLE	                  246 /* ((unsigned char)'\366') */

/* Note, if the future brings...
 * #define MLTK_128BIT_UNSIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER   ((unsigned char)'\257')
 * with  Range[0, 340282366920938463463374607431768211456 (*approximately 3.40282e+38*)]
 * the dynamic range is still a monotonically increasing function of the token value.
 * An implementation might choose to set the high varient bit to mainain this property
 * and dispatch more efficiently by avoiding overflow checks
 */

#define MLNE__SELECTOR( dtok, stok) \
	(((dtok) << 8) | (stok)) /* maybe should mask of high word and cast stok */

#define MLNE__SIZESELECTOR( tok) MLNE__SELECTOR( 0, tok) 
#define MLNE__INITSELECTOR (0)
#define MLNE__TOSTRINGSELECTOR( tok) MLNE__SELECTOR( MLNE__IS_REAL(tok) ? MLTKREAL : MLTKINT, tok)
#define MLNE__FROMSTRINGSELECTOR( dtok, stok) MLNE__SELECTOR( dtok, stok) 

#define MLNE__STOK( selector) ( (selector) & 0x000000FF)
#define MLNE__DTOK( selector) ( ((selector) & 0x0000FF00)>>8)

#define MLNE__IS_BINARY( tok) ((tok) & BINARY_MASK)
#define MLNE__IS_REAL( tok) ((tok) & REAL_MASK)
#define MLNE__TEXT_TOKEN( tok) (MLNE__IS_REAL( tok) ? MLTKREAL : MLTKINT)



MLNDECL( long_et, mlne__dispatch, ( unsigned long selector, void* dptr, void* sptr, long* countp));
MLNPROC( long_et, dispatch_procptr_mlnet, ( unsigned long selector, void* dptr, void* sptr, long* countp));

/* will null terminate  strings only if countp is null */

enum {
	uppNumDispatchProcInfo = kPascalStackBased
		 | RESULT_SIZE(SIZE_CODE(sizeof(long_et)))
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(unsigned long)))
		 | STACK_ROUTINE_PARAMETER(2, SIZE_CODE(sizeof(void*)))
		 | STACK_ROUTINE_PARAMETER(3, SIZE_CODE(sizeof(void*)))
		 | STACK_ROUTINE_PARAMETER(4, SIZE_CODE(sizeof(long*)))
};


#if GENERATINGCFM && GENERATING68K /* on the powerpc use standard native shared library calling conventions */
	typedef UniversalProcPtr dispatch_function_mlnet;
#	define call_num_dispatch( num_dispatch, selector, dptr, sptr, countp) \
		(long_et)CallUniversalProc( (num_dispatch), uppNumDispatchProcInfo, (selector), (dptr), (sptr), (countp))
#	define new_num_dispatch(num_dispatch) \
		NewRoutineDescriptor((ProcPtr)(num_dispatch), uppNumDispatchProcInfo, GetCurrentArchitecture())
#else
	typedef dispatch_procptr_mlnet dispatch_function_mlnet;
#	define call_num_dispatch(num_dispatch, selector, dptr, sptr, countp) \
		((*(num_dispatch))( (selector), (dptr), (sptr), (countp)))
#	define new_num_dispatch(num_dispatch) (num_dispatch)
#endif



ML_END_EXTERN_C


#endif /* _MLNUMENV_H */

#ifndef MLINTERFACE
/* syntax error */ )
#endif

#define MLTK_CSHORT_P       (( BINARY_MASK | SIZEVARIANT_MASK | 1))              /* 193 */
#define MLTK_CINT_P         (( BINARY_MASK | SIZEVARIANT_MASK | 2))              /* 194 */
#define MLTK_CLONG_P        (( BINARY_MASK | SIZEVARIANT_MASK | 3))              /* 195 */
#define MLTK_CFLOAT_P       (( BINARY_MASK | SIZEVARIANT_MASK | REAL_MASK | 1))  /* 209 */
#define MLTK_CDOUBLE_P      (( BINARY_MASK | SIZEVARIANT_MASK | REAL_MASK | 2))  /* 210 */
#define MLTK_CLONGDOUBLE_P  (( BINARY_MASK | SIZEVARIANT_MASK | REAL_MASK | 3))  /* 211 */


#define MLTK_CUCHAR  MLTK_8BIT_UNSIGNED_2sCOMPLEMENT_INTEGER
#define MLTK_MLUCHAR MLTK_8BIT_UNSIGNED_2sCOMPLEMENT_INTEGER

#if MACINTOSH_MATHLINK
	/* two private tokens */
	/* Interval[{-1.189731495357231765e+4932, 1.189731495357231765e+4932}] */
#	define MLTK_80BIT_SANE_EXTENDED  152 /* ((unsigned char)'\230') */
#	define MLTK_96BIT_68881_EXTENDED 154 /* ((unsigned char)'\232') */
#endif

#if POWERMACINTOSH_MATHLINK
	/* one private token */
#	define MLTK_128BIT_LONGDOUBLE  158 /* ((unsigned char)'\236') */
#endif

#if WINDOWS_MATHLINK
	/* one private token */
	/* Interval[{-1.189731495357231765e+4932, 1.189731495357231765e+4932}] */
#	define MLTK_INTEL_80BIT_EXTENDED  216 /* ((unsigned char)'\330') */
#	if MLINTERFACE > 1
#		define NEW_WIN32_NUMENV 1
#	endif
#endif





#if M68KMACINTOSH_MATHLINK
#	define MATHLINK_NUMERICS_ENVIRONMENT_ID "Sep 16 1996, 23:14:20"

#	define MLTK_CSHORT MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#	define MLTK_CLONG  MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER

#	define MLTK_CFLOAT MLTK_BIGENDIAN_IEEE754_SINGLE

#	if defined(__MWERKS__)
#		if __fourbyteints__
#			define MLTK_CINT MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		else
#			define MLTK_CINT MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		endif
#		if __MC68881__
#			define MLTK_CLONGDOUBLE  MLTK_96BIT_68881_EXTENDED
#		else
#			define MLTK_CLONGDOUBLE  MLTK_80BIT_SANE_EXTENDED
#		endif
#		if __IEEEdoubles__ || __ieeedoubles__
#			define MLTK_CDOUBLE  MLTK_BIGENDIAN_IEEE754_DOUBLE
#		else
#			define MLTK_CDOUBLE  MLTK_CLONGDOUBLE
#		endif
#	elif defined(THINK_C) || defined(SYMANTEC_C) || defined(SYMANTEC_CPLUS)
#		if __option(int_4)
#			define MLTK_CINT MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		else
#			define MLTK_CINT MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		endif
#		if __option(native_fp) && !__option(mc68881)
#			define MLTK_CLONGDOUBLE  MLTK_80BIT_SANE_EXTENDED
#		else
#			define MLTK_CLONGDOUBLE  MLTK_96BIT_68881_EXTENDED
#		endif
#		if __option(double_8)
#			define MLTK_CDOUBLE  MLTK_BIGENDIAN_IEEE754_DOUBLE
#		else
#			define MLTK_CDOUBLE  MLTK_CLONGDOUBLE
#		endif
#	else /* applec */
#		define MLTK_CINT MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_CDOUBLE  MLTK_BIGENDIAN_IEEE754_DOUBLE
#		if mc68881
#			define MLTK_CLONGDOUBLE  MLTK_96BIT_68881_EXTENDED
#		else
#			define MLTK_CLONGDOUBLE  MLTK_80BIT_SANE_EXTENDED
#		endif
#	endif

#	if 0 /* no more statically linked Macintosh libraries STATICALLY_LINKED_MATHLINK */
#		define MLTK_MLSHORT       MLTK_CSHORT
#		define MLTK_MLINT         MLTK_CINT
#		define MLTK_MLLONG        MLTK_CLONG
#		define MLTK_MLFLOAT       MLTK_CFLOAT
#		define MLTK_MLDOUBLE      MLTK_CDOUBLE
#		define MLTK_MLLONGDOUBLE  MLTK_CLONGDOUBLE
#	else
#		define MLTK_MLSHORT       MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_MLINT         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_MLLONG        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_MLFLOAT       MLTK_BIGENDIAN_IEEE754_SINGLE
#		define MLTK_MLDOUBLE      MLTK_BIGENDIAN_IEEE754_DOUBLE
#		define MLTK_MLLONGDOUBLE  MLTK_80BIT_SANE_EXTENDED
#	endif

#elif POWERMACINTOSH_MATHLINK
#	define MATHLINK_NUMERICS_ENVIRONMENT_ID "newdog"

#	define MLTK_CSHORT MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#	define MLTK_CINT   MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#	define MLTK_CLONG  MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER

#	define MLTK_CFLOAT MLTK_BIGENDIAN_IEEE754_SINGLE
#	define MLTK_CDOUBLE  MLTK_BIGENDIAN_IEEE754_DOUBLE


#	ifndef MLTK_CLONGDOUBLE
#		if defined(__MWERKS__) || defined(SYMANTEC_C) || defined(SYMANTEC_CPLUS)
#			define MLTK_CLONGDOUBLE  MLTK_BIGENDIAN_IEEE754_DOUBLE
#		elif defined(__MRC__) && __MRC__ >= 0x0300 && __MRC__ != 0x0800 
			/* MrC version 1.0 defined __MRC__ to be 0x0800 presumably because of its Symantec heritage */
			/* One cannot querry value of -ldsize with old MrC or PPCC -- assume -ldsize 128 */
#			if __option(ldsize128)
#				define MLTK_CLONGDOUBLE  MLTK_128BIT_LONGDOUBLE
#			else
#				define MLTK_CLONGDOUBLE  MLTK_BIGENDIAN_IEEE754_DOUBLE
#			endif
#		else
#			define MLTK_CLONGDOUBLE  MLTK_128BIT_LONGDOUBLE
#		endif
#	endif

#	define MLTK_MLSHORT       MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#	define MLTK_MLINT         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#	define MLTK_MLLONG        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#	define MLTK_MLFLOAT       MLTK_BIGENDIAN_IEEE754_SINGLE
#	define MLTK_MLDOUBLE      MLTK_BIGENDIAN_IEEE754_DOUBLE
#	define MLTK_MLLONGDOUBLE  MLTK_128BIT_LONGDOUBLE

#elif SUN_MATHLINK


#	if __sparc || __sparc__ || sparc

#		define MLTK_CSHORT       MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_CINT         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_CLONG        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_CFLOAT       MLTK_BIGENDIAN_IEEE754_SINGLE
#		define MLTK_CDOUBLE      MLTK_BIGENDIAN_IEEE754_DOUBLE
#		ifndef MLTK_CLONGDOUBLE
#			if __SUNPRO_C >= 0x301
				/* one private token */
#				define MLTK_128BIT_LONGDOUBLE  158 /* ((unsigned char)'\236') */
#				define MLTK_CLONGDOUBLE  MLTK_128BIT_LONGDOUBLE
#			elif defined(__GNUC__) || defined(__GNUG__)
#				define MLTK_CLONGDOUBLE  MLTK_CDOUBLE
#			else
				/* no error directive here as the user may be
				 * using a different compiler.  Some macros
				 * simply won't be available.
				 */
#			endif
#		endif

#		define MLTK_MLSHORT       MLTK_16BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_MLINT         MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_MLLONG        MLTK_32BIT_SIGNED_2sCOMPLEMENT_BIGENDIAN_INTEGER
#		define MLTK_MLFLOAT       MLTK_BIGENDIAN_IEEE754_SINGLE
#		define MLTK_MLDOUBLE      MLTK_BIGENDIAN_IEEE754_DOUBLE
#		define MLTK_MLLONGDOUBLE  MLTK_128BIT_LONGDOUBLE

#	elif __i386 || __i386__ || i386

/* syntax error */ )

#		if __SUNPRO_C >= 0x301
			/* one private token */
			/* Interval[{-1.189731495357231765e+4932, 1.189731495357231765e+4932}] */
#			define MLTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED 218 /* ((unsigned char)'\332') */
#			define MLTK_CSHORT       MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#			define MLTK_CINT         MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#			define MLTK_CLONG        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#			define MLTK_CFLOAT       MLTK_LITTLEENDIAN_IEEE754_SINGLE
#			define MLTK_CDOUBLE      MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#			define MLTK_CLONGDOUBLE  MLTK_96BIT_HIGHPADDED_INTEL_80BIT_EXTENDED
#		elif defined(__GNUC__) || defined(__GNUG__)
			/* no error directive here as the user may be
			 * using a different compiler.  Some macros
			 * simply won't be available.
			 */
#		else
			/* no error directive here as the user may be
			 * using a different compiler.  Some macros
			 * simply won't be available.
			 */
#		endif


#	else
/* syntax error */ )
#	endif


#elif WIN16_MATHLINK || (WIN32_MATHLINK && NEW_WIN32_NUMENV)
#	if WIN16_MATHLINK
#		define MATHLINK_NUMERICS_ENVIRONMENT_ID "poodle"
#	elif WIN32_MATHLINK
#		define MATHLINK_NUMERICS_ENVIRONMENT_ID "setter"
#	endif

#	define MLTK_CSHORT       MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#	define MLTK_CLONG        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#	define MLTK_CFLOAT       MLTK_LITTLEENDIAN_IEEE754_SINGLE
#	define MLTK_CDOUBLE      MLTK_LITTLEENDIAN_IEEE754_DOUBLE

#	if WIN16_MATHLINK
#		define MLTK_CINT         MLTK_CSHORT
#	elif WIN32_MATHLINK
#		define MLTK_CINT         MLTK_CLONG
#	endif

#	if __WATCOMC__ || __SC__
#		define MLTK_CLONGDOUBLE  MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#	elif __BORLANDC__ || __BCPLUSPLUS__ || __TURBOC__ || __TCPLUSPLUS__
#		define MLTK_CLONGDOUBLE  MLTK_INTEL_80BIT_EXTENDED
#	elif _MSC_VER
#		if WIN16_MATHLINK
#			define MLTK_CLONGDOUBLE  MLTK_INTEL_80BIT_EXTENDED
#		elif WIN32_MATHLINK
#			define MLTK_CLONGDOUBLE  MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#		endif
#	endif


#	define MLTK_MLSHORT       MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#	if WIN16_MATHLINK
#		define MLTK_MLINT     MLTK_16BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#	elif WIN32_MATHLINK
#		define MLTK_MLINT     MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#	endif
#	define MLTK_MLLONG        MLTK_32BIT_SIGNED_2sCOMPLEMENT_LITTLEENDIAN_INTEGER
#	define MLTK_MLFLOAT       MLTK_LITTLEENDIAN_IEEE754_SINGLE
#	define MLTK_MLDOUBLE      MLTK_LITTLEENDIAN_IEEE754_DOUBLE
#	define MLTK_MLLONGDOUBLE  MLTK_INTEL_80BIT_EXTENDED

#else

#	if WIN32_MATHLINK
#		define MATHLINK_NUMERICS_ENVIRONMENT_ID "Sep 13 1996, 13:46:34"
#	endif

#	define MLTK_CSHORT       MLTK_CSHORT_P
#	define MLTK_CINT         MLTK_CINT_P
#	define MLTK_CLONG        MLTK_CLONG_P
#	define MLTK_CFLOAT       MLTK_CFLOAT_P
#	define MLTK_CDOUBLE      MLTK_CDOUBLE_P
#	define MLTK_CLONGDOUBLE  MLTK_CLONGDOUBLE_P

#	define MLTK_MLSHORT       MLTK_CSHORT_P
#	define MLTK_MLINT         MLTK_CINT_P
#	define MLTK_MLLONG        MLTK_CLONG_P
#	define MLTK_MLFLOAT       MLTK_CFLOAT_P
#	define MLTK_MLDOUBLE      MLTK_CDOUBLE_P
#	define MLTK_MLLONGDOUBLE  MLTK_CLONGDOUBLE_P

#endif

/* Objects of these numeric types exist in MathLink only in the numerics
 * environment and, unfortunately, in the "stack frames" of the functions that
 * put atomic numbers like MLPutInteger.  These C types are used by client
 * programs solely for type-checked access to the BinaryNumber functions.
 */
typedef unsigned char uchar_nt;
typedef uchar_nt     FAR * ucharp_nt;
typedef ucharp_nt    FAR * ucharpp_nt;

typedef short              short_nt;
typedef short_nt     FAR * shortp_nt;
typedef shortp_nt    FAR * shortpp_nt;

typedef int                int_nt;
typedef int_nt       FAR * intp_nt;
typedef intp_nt      FAR * intpp_nt;

typedef long               long_nt;
typedef long_nt      FAR * longp_nt;
typedef longp_nt     FAR * longpp_nt;

typedef float              float_nt;
typedef float_nt     FAR * floatp_nt;
typedef floatp_nt    FAR * floatpp_nt;

typedef double             double_nt;
typedef double_nt    FAR * doublep_nt;
typedef doublep_nt   FAR * doublepp_nt;

#ifndef CC_SUPPORTS_LONG_DOUBLE
#	if defined( __STDC__) || defined(__cplusplus) || ! UNIX_MATHLINK
#		define CC_SUPPORTS_LONG_DOUBLE 1
#	else
#		define CC_SUPPORTS_LONG_DOUBLE MLPROTOTYPES
#	endif
#endif

struct _i87extended_nt { unsigned short w[5];};

#if CC_SUPPORTS_LONG_DOUBLE
#	ifndef __extended_nt__
#		if WINDOWS_MATHLINK && (MLTK_CLONGDOUBLE != MLTK_MLLONGDOUBLE) /* subtle predicate that works for old and new windows numenvs */
#			define __extended_nt__ struct _i87extended_nt
#		else
#			define __extended_nt__ long double
#		endif
#	endif
	typedef __extended_nt__    extended_nt;
	typedef extended_nt  FAR * extendedp_nt;
	typedef extendedp_nt FAR * extendedpp_nt;
#endif

#endif /* _MLNTYPES_H */

#ifndef _ML0TYPES_H
#define _ML0TYPES_H


#if USING_OLD_TYPE_NAMES
typedef charp_ct ml_charp;
typedef charpp_ct ml_charpp;
typedef charppp_ct ml_charppp;
typedef ucharp_ct ml_ucharp;
typedef longp_ct ml_longp;
typedef longpp_ct ml_longpp;
typedef ulongp_ct ml_ulongp;
typedef shortp_nt ml_shortp;
typedef shortpp_nt ml_shortpp;
typedef intp_nt ml_intp;
typedef intpp_nt ml_intpp;
typedef floatp_nt ml_floatp;
typedef floatpp_nt ml_floatpp;
typedef doublep_nt ml_doublep;
typedef doublepp_nt ml_doublepp;
#if CC_SUPPORTS_LONG_DOUBLE
typedef extended_nt ml_extended;
typedef extendedp_nt ml_extendedp;
typedef extendedpp_nt ml_extendedpp;
#endif
typedef charp_ct MLBuffer;
typedef kcharp_ct MLKBuffer;
typedef charpp_ct MLBufferArray;

#endif

#endif /* _ML0TYPES_H */

ML_EXTERN_C

#ifndef _MLSTDDEV_H
#define _MLSTDDEV_H


#if WINDOWS_MATHLINK
#endif

#if OS2_MATHLINK
#	include <os2def.h>
#endif






typedef void FAR * dev_world;
typedef MLINK dev_cookie;

typedef dev_world FAR * dev_worldp;
typedef dev_cookie FAR * dev_cookiep;

typedef  MLAllocatorUPP dev_allocator;
#define call_dev_allocator CallMLAllocatorProc
#define new_dev_allocator NewMLAllocatorProc

typedef  MLDeallocatorUPP dev_deallocator;
#define call_dev_deallocator CallMLDeallocatorProc
#define new_dev_deallocator NewMLDeallocatorProc


typedef dev_main_type world_main_type;

#define MLSTDWORLD_INIT        16
#define MLSTDWORLD_DEINIT      17
#define MLSTDWORLD_MAKE        18
#define MLSTDDEV_CONNECT_READY 19
#define MLSTDDEV_CONNECT       20
#define MLSTDDEV_DESTROY       21

#define MLSTDDEV_SET_YIELDER   22
#define MLSTDDEV_GET_YIELDER   23

#define MLSTDDEV_WRITE_MSG     24
#define MLSTDDEV_HAS_MSG       25
#define MLSTDDEV_READ_MSG      26
#define MLSTDDEV_SET_HANDLER   27
#define MLSTDDEV_GET_HANDLER   28



#define T_WORLD_INIT        MLSTDWORLD_INIT
#define T_WORLD_DEINIT      MLSTDWORLD_DEINIT
#define T_WORLD_MAKE        MLSTDWORLD_MAKE
#define T_DEV_CONNECT_READY MLSTDDEV_CONNECT_READY
#define T_DEV_CONNECT       MLSTDDEV_CONNECT
#define T_DEV_DESTROY       MLSTDDEV_DESTROY

#define T_DEV_SET_YIELDER   MLSTDDEV_SET_YIELDER
#define T_DEV_GET_YIELDER   MLSTDDEV_GET_YIELDER

#define T_DEV_WRITE_MSG     MLSTDDEV_WRITE_MSG
#define T_DEV_HAS_MSG       MLSTDDEV_HAS_MSG
#define T_DEV_READ_MSG      MLSTDDEV_READ_MSG
#define T_DEV_SET_HANDLER   MLSTDDEV_SET_HANDLER
#define T_DEV_GET_HANDLER   MLSTDDEV_GET_HANDLER




typedef unsigned long dev_mode;
/* edit here and in mathlink.r */
#define NOMODE           ((dev_mode)0x0000)
#define LOOPBACKBIT      ((dev_mode)0x0001)
#define LISTENBIT        ((dev_mode)0x0002)
#define CONNECTBIT       ((dev_mode)0x0004)
#define LAUNCHBIT        ((dev_mode)0x0008)
#define PARENTCONNECTBIT ((dev_mode)0x0010)
#define READBIT          ((dev_mode)0x0020)
#define WRITEBIT         ((dev_mode)0x0040)
#define SERVERBIT        ((dev_mode)0x0080)
#define ANYMODE          (~(dev_mode)0)

typedef dev_mode FAR * dev_modep;





typedef unsigned long dev_options;

#define _DefaultOptions      ((dev_options)0x00000000)

#define _NetworkVisibleMask  ((dev_options)0x00000003)
#define _BrowseMask          ((dev_options)0x00000010)
#define _NonBlockingMask     ((dev_options)0x00000020)
#define _InteractMask        ((dev_options)0x00000100)
#define _VersionMask         ((dev_options)0x0F000000)

#define _NetworkVisible      ((dev_options)0x00000000)
#define _LocallyVisible      ((dev_options)0x00000001)
#define _InternetVisible     ((dev_options)0x00000002)

#define _Browse              ((dev_options)0x00000000)
#define _DontBrowse          ((dev_options)0x00000010)

#define _NonBlocking         ((dev_options)0x00000000)
#define _Blocking            ((dev_options)0x00000020)

#define _Interact            ((dev_options)0x00000000)
#define _DontInteract        ((dev_options)0x00000100)




/* values returned by selector DEVICE_TYPE */
#define UNREGISTERED_TYPE	0
#define UNIXPIPE_TYPE	1
#define UNIXSOCKET_TYPE 2
#define PPC_TYPE	3
#define MACTCP_TYPE	4
#define LOOPBACK_TYPE	5
#define COMMTB_TYPE	6
#define ADSP_TYPE	7
#define LOCAL_TYPE	8
#define WINLOCAL_TYPE	9
#define WINFMAP_TYPE	10

/* selectors */
#define DEVICE_TYPE 0                                       /* long */
#define DEVICE_NAME 1                                       /* char */
#define DEVICE_WORLD_ID 5                                   /* char */
#define PIPE_FD                (UNIXPIPE_TYPE * 256 + 0)    /* int */
#define PIPE_CHILD_PID         (UNIXPIPE_TYPE * 256 + 1)    /* int */
#define SOCKET_FD              (UNIXSOCKET_TYPE * 256 + 0)  /* int */
#define SOCKET_PARTNER_ADDR    (UNIXSOCKET_TYPE * 256 + 1)  /* unsigned long */
#define SOCKET_PARTNER_PORT    (UNIXSOCKET_TYPE * 256 + 2)  /* unsigned short */
#define PPC_SESS_REF_NUM       (PPC_TYPE * 256 + 0)         /* PPCSessRefNum */
#define PPC_PARTNER_PSN        (PPC_TYPE * 256 + 1)         /* ProcessSerialNumber */
#define PPC_PARTNER_LOCATION   (PPC_TYPE * 256 + 2)         /* LocationNameRec */
#define PPC_PARTNER_PORT       (PPC_TYPE * 256 + 3)         /* PPCPortRec */
#define MACTCP_STREAM          (MACTCP_TYPE * 256 + 0)      /* StreamPtr */
#define MACTCP_PARTNER_ADDR    (MACTCP_TYPE * 256 + 1)      /* ip_addr */
#define MACTCP_PARTNER_PORT    (MACTCP_TYPE * 256 + 2)      /* tcp_port */
#define MACTCP_IPDRIVER        (MACTCP_TYPE * 256 + 3)      /* short */
#define MACTCP_SETSIMPLESOCKET (MACTCP_TYPE * 256 + 9)      /* buf, buflen ignored */
#define COMMTB_CONNHANDLE      (COMMTB_TYPE * 256 + 0)      /* ConnHandle */
#define ADSP_CCBREFNUM         (ADSP_TYPE * 256 + 0)        /* short */
#define ADSP_IOCREFNUM         (ADSP_TYPE * 256 + 3)        /* short */

#define	WINDOWS_SET_NOTIFY_WINDOW     2330 /* HWND */
#define	WINDOWS_REMOVE_NOTIFY_WINDOW  2331 /* HWND */

/* info selectors */
#define WORLD_THISLOCATION 1        /* char */
#define WORLD_MODES 2               /* dev_mode */
#define WORLD_PROTONAME 3           /* char */
#define WORLD_STREAMCAPACITY 4      /* long */ /*this belongs in mlolddev.h*/
#define WORLD_ID DEVICE_WORLD_ID    /* char */



#ifndef MATHLINK_DEVICE_WORLD_ID
#define MATHLINK_DEVICE_WORLD_ID (__DATE__ ", " __TIME__)
#endif




#define YIELDVERSION 1

typedef long devyield_result;
typedef long devyield_place;
typedef long devyield_count;
typedef unsigned long devyield_sleep;

#define INTERNAL_YIELDING 0
#define MAKE_YIELDING 1
#define CONNECT_YIELDING 2
#define READ_YIELDING 3
#define WRITE_YIELDING 4
#define DESTROY_YIELDING 5
#define READY_YIELDING 6


typedef struct MLYieldParams FAR * MLYieldParameters;


#define MAX_SLEEP (600)
typedef struct MLYieldData{
	union {long l; double d; void FAR * p;} private_data[8];
} FAR * MLYieldDataPointer;

void MLNewYieldData P(( MLYieldDataPointer ydp   /* , dev_allocator, dev_deallocator */));
void MLFreeYieldData P(( MLYieldDataPointer ydp));
MLYieldParameters MLResetYieldData P(( MLYieldDataPointer ydp, devyield_place func_id));
mlapi_result   MLSetYieldParameter P(( MLYieldParameters yp, unsigned long selector, void* data, unsigned long* len));
mlapi_result   MLYieldParameter P(( MLYieldParameters yp, unsigned long selector, void* data, unsigned long* len));
devyield_sleep MLSetSleepYP P(( MLYieldParameters yp, devyield_sleep sleep));
devyield_count MLSetCountYP P(( MLYieldParameters yp, devyield_count count));


enum { MLSleepParameter = 1, MLCountParameter, MLPlaceParameter};





MLYPROC( devyield_result, MLYielderProcPtr, (MLINK mlp, MLYieldParameters yp));
typedef	MLYielderProcPtr MLDeviceYielderProcPtr;

enum {
	uppMLYielderProcInfo = kPascalStackBased
		 | RESULT_SIZE(SIZE_CODE(sizeof(devyield_result)))
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(MLINK)))
		 | STACK_ROUTINE_PARAMETER(2, SIZE_CODE(sizeof(MLYieldParameters))),
	uppMLDeviceYielderProcInfo = uppMLYielderProcInfo
};

#if GENERATINGCFM
	typedef UniversalProcPtr MLYielderUPP, MLDeviceYielderUPP;
#	define NewMLYielderProc(userRoutine) \
		NewRoutineDescriptor(MLYielderCast((userRoutine)), uppMLYielderProcInfo, GetCurrentArchitecture())
#elif WIN16_MATHLINK
	typedef __MLProcPtr__ MLYielderUPP, MLDeviceYielderUPP;
#	define NewMLYielderProc( userRoutine) \
		(MakeProcInstance( MLYielderCast(userRoutine), MLInstance))
#else
	typedef MLYielderProcPtr MLYielderUPP, MLDeviceYielderUPP;
#	define NewMLYielderProc(userRoutine) (userRoutine)
#endif

#define NewMLDeviceYielderProc NewMLYielderProc

typedef  MLYielderUPP MLYieldFunctionType;
typedef unsigned long MLYieldFunctionObject; /* bugcheck should I change this back to void* for 64 bit machines */

typedef  MLYieldFunctionObject dev_yielder;
typedef dev_yielder FAR* dev_yielderp;








typedef unsigned long dev_message;
typedef dev_message FAR * dev_messagep;


MLMPROC( void, MLHandlerProcPtr, (MLINK mlp, dev_message m, dev_message n));
typedef MLHandlerProcPtr MLDeviceHandlerProcPtr;


enum {
	uppMLHandlerProcInfo = kPascalStackBased
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(MLINK)))
		 | STACK_ROUTINE_PARAMETER(2, SIZE_CODE(sizeof(dev_message)))
		 | STACK_ROUTINE_PARAMETER(3, SIZE_CODE(sizeof(dev_message))),
	uppMLDeviceHandlerProcInfo = uppMLHandlerProcInfo
};

#if GENERATINGCFM
	typedef UniversalProcPtr MLHandlerUPP, MLDeviceHandlerUPP;
#	define NewMLHandlerProc(userRoutine) \
		NewRoutineDescriptor(MLHandlerCast((userRoutine)), uppMLHandlerProcInfo, GetCurrentArchitecture())
#elif WIN16_MATHLINK
	typedef __MLProcPtr__ MLHandlerUPP, MLDeviceHandlerUPP;
#	define NewMLHandlerProc( userRoutine) \
		(MakeProcInstance( MLHandlerCast(userRoutine), MLInstance))
#else
	typedef MLHandlerProcPtr MLHandlerUPP, MLDeviceHandlerUPP;
#	define NewMLHandlerProc(userRoutine) (userRoutine)
#endif

#define NewMLDeviceHandlerProc NewMLHandlerProc

typedef  MLHandlerUPP MLMessageHandlerType;
typedef unsigned long MLMessageHandlerObject;

typedef  MLMessageHandlerObject dev_msghandler;
typedef dev_msghandler FAR* dev_msghandlerp;



#endif /* _MLSTDDEV_H */



/* explicitly not protected by _MLSTDDEV_H in case MLDECL is redefined for multiple inclusion */

/*bugcheck //should the rest of YP stuff be exported? */
MLDECL( devyield_sleep,         MLSleepYP,               ( MLYieldParameters yp));
MLDECL( devyield_count,         MLCountYP,               ( MLYieldParameters yp));
MLDECL( MLYieldFunctionObject,  MLCreateYieldFunction,   ( MLEnvironment ep, MLYieldFunctionType yf, MLPointer reserved)); /* reserved must be 0 */
#ifndef MLINTERFACE
/* syntax error */ )
#endif
#if MLINTERFACE > 1
MLDECL( MLYieldFunctionObject,  MLCreateYieldFunction0,   ( MLEnvironment ep, MLYieldFunctionType yf, MLPointer reserved)); /* reserved must be 0 */
#endif
MLDECL( MLYieldFunctionType,    MLDestroyYieldFunction,  ( MLYieldFunctionObject yfo));
MLDECL( devyield_result,        MLCallYieldFunction,     ( MLYieldFunctionObject yfo, MLINK mlp, MLYieldParameters p));
MLDECL( MLMessageHandlerObject, MLCreateMessageHandler,  ( MLEnvironment ep, MLMessageHandlerType mh, MLPointer reserved)); /* reserved must be 0 */
#if MLINTERFACE > 1
MLDECL( MLMessageHandlerObject, MLCreateMessageHandler0,  ( MLEnvironment ep, MLMessageHandlerType mh, MLPointer reserved)); /* reserved must be 0 */
#endif
MLDECL( MLMessageHandlerType,   MLDestroyMessageHandler, ( MLMessageHandlerObject mho));
MLDECL( void,                   MLCallMessageHandler,    ( MLMessageHandlerObject mho, MLINK mlp, dev_message m, dev_message n));


/* just some type-safe casts */
MLDECL( __MLProcPtr__, MLYielderCast, ( MLYielderProcPtr yp));
MLDECL( __MLProcPtr__, MLHandlerCast, ( MLHandlerProcPtr mh));

ML_END_EXTERN_C



#ifndef _MLMAKE_H
#define _MLMAKE_H

/* --binding layer-- */
/*************** Starting MathLink ***************/

#define MLPARAMETERSIZE_R1 256
#define MLPARAMETERSIZE 256
typedef char FAR * MLParametersPointer;
typedef char MLParameters[MLPARAMETERSIZE];

#define MLLoopBackOpen MLLoopbackOpen



ML_EXTERN_C
MLUPROC( void, MLUserProcPtr, (MLINK));
ML_END_EXTERN_C

enum {
	uppMLUserFunctionProcInfo = kPascalStackBased
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(MLINK)))
};

#if GENERATINGCFM
	typedef UniversalProcPtr MLUserUPP;
#	define NewMLUserProc(userRoutine) \
		NewRoutineDescriptor(MLUserCast((userRoutine)), uppMLUserFunctionProcInfo, GetCurrentArchitecture())
#else
	typedef MLUserProcPtr MLUserUPP;
#	define NewMLUserProc(userRoutine) (userRoutine)
#endif

typedef MLUserUPP MLUserFunctionType;
typedef MLUserFunctionType FAR * MLUserFunctionTypePointer;





/* The following defines are
 * currently for internal use only.
 */


/* edit here and in mldevice.h and mathlink.r */
#define MLNetworkVisibleMask ((unsigned long)0x00000003)
#define MLBrowseMask         ((unsigned long)0x00000010)
#define MLNonBlockingMask    ((unsigned long)0x00000020)
#define MLInteractMask       ((unsigned long)0x00000100)
#define MLVersionMask        ((unsigned long)0x0000F000)

#define MLDefaultOptions     ((unsigned long)0x00000000)
#define MLNetworkVisible     ((unsigned long)0x00000000)
#define MLLocallyVisible     ((unsigned long)0x00000001)
#define MLInternetVisible    ((unsigned long)0x00000002)

#define MLBrowse             ((unsigned long)0x00000000)
#define MLDontBrowse         ((unsigned long)0x00000010)

#define MLNonBlocking        ((unsigned long)0x00000000)
#define MLBlocking           ((unsigned long)0x00000020)

#define MLInteract           ((unsigned long)0x00000000)
#define MLDontInteract       ((unsigned long)0x00000100)

#endif /* _MLMAKE_H */


/* explicitly not protected by _MLMAKE_H in case MLDECL is redefined for multiple inclusion */


ML_EXTERN_C
MLDECL( ulong_ct, MLNewParameters,     ( MLParametersPointer p, ulong_ct rev, ulong_ct apirev));
MLDECL( void,     MLSetAllocParameter, ( MLParametersPointer p, MLAllocator allocator, MLDeallocator deallocator));
#ifndef MLINTERFACE
/* syntax error */ )
#endif
#if MLINTERFACE > 1
MLDECL( long,     MLAllocParameter,       (MLParametersPointer p, MLAllocatorp allocatorp, MLDeallocatorp deallocatorp));
MLDECL( long,     MLSetResourceParameter, (MLParametersPointer p, kcharp_ct path));
MLDECL( long,     MLSetDeviceParameter,   (MLParametersPointer p, kcharp_ct devspec));
#endif
MLDECL( long,     MLErrorParameter,    ( MLParametersPointer p));

MLDECL( MLEnvironment, MLInitialize,   ( MLParametersPointer p)); /* pass in NULL */
MLDECL( void,          MLDeinitialize, ( MLEnvironment env));

/* or, if you use MLOpenArgv, ...*/

MLDECL( MLEnvironment, MLBegin, ( MLParametersPointer p)); /* pass in NULL */
MLDECL( void,          MLEnd,   ( MLEnvironment env));

#if MLNTESTPOINTS < 1
#undef MLNTESTPOINTS
#define MLNTESTPOINTS 1
#endif
MLDECL( long_et, MLTestPoint1, ( MLEnvironment ep, ulong_ct selector, voidp_ct p1, voidp_ct p2, longp_ct np));

#ifndef MLINTERFACE
/* syntax error */ )
#endif
#if MLINTERFACE > 1

#if MLNTESTPOINTS < 2
#undef MLNTESTPOINTS
#define MLNTESTPOINTS 2
#endif
MLDECL( void,    MLTestPoint2,     ( MLINK mlp));

#if MLNTESTPOINTS < 3
#undef MLNTESTPOINTS
#define MLNTESTPOINTS 3
#endif
MLDECL( ulong_ct,    MLTestPoint3,     ( MLINK mlp));

#if MLNTESTPOINTS < 4
#undef MLNTESTPOINTS
#define MLNTESTPOINTS 4
#endif
MLDECL( ulong_ct,    MLTestPoint4,     ( MLINK mlp));

MLDECL( long_et, MLNumberControl0, ( MLEnvironment ep, ulong_ct selector, voidp_ct p1, voidp_ct p2, longp_ct np));
#else
extern long_et MLNumberControl0( MLEnvironment ep, ulong_ct selector, voidp_ct p1, voidp_ct p2, longp_ct np);
#endif


/*************** Connection interface ***************/
MLDECL( MLINK,         MLCreate0,       ( MLEnvironment ep, dev_type dev, dev_main_type dev_main, longp_ct errp));
MLDECL( MLINK,         MLMake,          ( MLPointer ep, dev_type dev, dev_main_type dev_main, longp_ct errp));
MLDECL( void,          MLDestroy,       ( MLINK mlp, dev_typep devp, dev_main_typep dev_mainp));

MLDECL( long,          MLFeatureString, ( MLINK mlp, charp_ct buf, long buffsize));
MLDECL( charpp_ct,     MLFilterArgv0,   ( MLEnvironment ep, charpp_ct argv, charpp_ct argv_end));
MLDECL( MLINK,         MLOpenArgv,      ( MLEnvironment ep, charpp_ct argv, charpp_ct argv_end, longp_ct errp));
MLDECL( MLINK,         MLOpenString,    ( MLEnvironment ep, kcharp_ct command_line, longp_ct errp));
MLDECL( MLINK,         MLLoopbackOpen,  ( MLEnvironment ep, longp_ct errp));
MLDECL( MLINK,         MLLoopbackOpen0, ( MLEnvironment ep, kcharp_ct features, longp_ct errp));
MLDECL( int_ct,        MLStringToArgv,  ( kcharp_ct commandline, charp_ct buf, charpp_ct argv, int_ct len));
MLDECL( long,          MLScanString,    ( charpp_ct argv, charppp_ct argv_end, charpp_ct commandline, charpp_ct buf));
MLDECL( long,          MLPrintArgv,     ( charp_ct buf, charpp_ct buf_endp, charppp_ct argvp, charpp_ct argv_end));

MLDECL( kcharp_ct,     MLErrorMessage,  ( MLINK mlp));
MLDECL( kcharp_ct,     MLErrorString,   ( MLEnvironment env, long err));

MLDECL( MLINK,         MLOpen,          ( int_ct argc, charpp_ct argv));
MLDECL( MLINK,         MLOpenInEnv,     ( MLEnvironment env, int_ct argc, charpp_ct argv, longp_ct errp));
MLDECL( MLINK,         MLOpenS,         ( kcharp_ct command_line));

MLDECL( MLINK,         MLDuplicateLink,   ( MLINK parentmlp, kcharp_ct name, longp_ct errp ));
MLDECL( mlapi_result,  MLConnect,         ( MLINK mlp));
#define MLActivate MLConnect

#ifndef __feature_setp__
#define __feature_setp__
typedef struct feature_set* feature_setp;
#endif
MLDECL( mlapi_result,  MLEstablish,       ( MLINK mlp, feature_setp features));

MLDECL( mlapi_result,  MLEstablishString, ( MLINK mlp, kcharp_ct features));
MLDECL( void,          MLClose,           ( MLINK mlp));

MLDECL( void,          MLSetUserData,   ( MLINK mlp, MLPointer data, MLUserFunctionType f));
MLDECL( MLPointer,     MLUserData,      ( MLINK mlp, MLUserFunctionTypePointer fp));
MLDECL( void,          MLSetUserBlock,  ( MLINK mlp, MLPointer userblock));
MLDECL( MLPointer,     MLUserBlock,     ( MLINK mlp));

/* just a type-safe cast */
MLDECL( __MLProcPtr__, MLUserCast, ( MLUserProcPtr f));





/* MLName returns a pointer to the link's name.
 * Links are generally named when they are created
 * and are based on information that is potentially
 * useful and is available at that time.
 * Do not attempt to deallocate the name's storage
 * through this pointer.  The storage should be
 * considered in read-only memory.
 */
MLDECL( kcharp_ct, MLName,    ( MLINK mlp));
MLDECL( long,      MLNumber,  ( MLINK mlp));
#if MLINTERFACE > 1
MLDECL( long,  MLToLinkID,  ( MLINK mlp));
MLDECL( MLINK, MLFromLinkID, ( MLEnvironment ep, long n));
#else
extern MLINK MLFromLinkID( MLEnvironment ep, long n);
#endif
MLDECL( charp_ct,  MLSetName, ( MLINK mlp, kcharp_ct name));



/* The following functions are
 * currently for internal use only.
 */

MLDECL( MLPointer, MLInit,   ( MLallocator alloc, MLdeallocator dealloc, MLPointer enclosing_environment));
MLDECL( void,      MLDeinit, ( MLPointer env));
MLDECL( MLPointer, MLEnclosingEnvironment, ( MLPointer ep));
MLDECL( MLPointer, MLinkEnvironment, ( MLINK mlp));

/* the following two functions are for internal use only */
MLDECL( MLYieldFunctionObject, MLDefaultYieldFunction,    ( MLEnvironment env));
MLDECL( mlapi_result,          MLSetDefaultYieldFunction, ( MLEnvironment env, MLYieldFunctionObject yf));

ML_END_EXTERN_C


#ifndef _MLERRORS_H
#define _MLERRORS_H


/*************** MathLink errors ***************/
/*
 * When some problem is detected within MathLink, routines
 * will return a simple indication of failure and store
 * an error code internally. (For routines that have nothing
 * else useful to return, success is indicated by returning
 * non-zero and failure by returning 0.)  MLerror() returns
 * the current error code;  MLErrorMessage returns an English
 * language description of the error.
 * The error MLEDEAD is irrecoverable.  For the others, MLClearError()
 * will reset the error code to MLEOK.
 */



#ifndef _MLERRNO_H
#define _MLERRNO_H

/* edit here and in mlerrstr.h */

#define MLEUNKNOWN       -1
#define MLEOK             0
#define MLEDEAD           1
#define MLEGBAD           2
#define MLEGSEQ           3
#define MLEPBTK           4
#define MLEPSEQ           5
#define MLEPBIG           6
#define MLEOVFL           7
#define MLEMEM            8
#define MLEACCEPT         9
#define MLECONNECT       10
#define MLECLOSED        11
#define MLEDEPTH         12  /* internal error */
#define MLENODUPFCN      13  /* stream cannot be duplicated */

#define MLENOACK         15  /* */
#define MLENODATA        16  /* */
#define MLENOTDELIVERED  17  /* */
#define MLENOMSG         18  /* */
#define MLEFAILED        19  /* */

#define MLEGETENDEXPR    20
#define MLEPUTENDPACKET  21 /* unexpected call of MLEndPacket */
                            /* currently atoms aren't
                             * counted on the way out so this error is raised only when
                             * MLEndPacket is called in the midst of an atom
                             */
#define MLENEXTPACKET    22
#define MLEUNKNOWNPACKET 23
#define MLEGETENDPACKET  24
#define MLEABORT         25
#define MLEMORE          26 /* internal error */
#define MLENEWLIB        27
#define MLEOLDLIB        28
#define MLEBADPARAM      29
#define MLENOTIMPLEMENTED 30


#define MLEINIT          32
#define MLEARGV          33
#define MLEPROTOCOL      34
#define MLEMODE          35
#define MLELAUNCH        36
#define MLELAUNCHAGAIN   37
#define MLELAUNCHSPACE   38
#define MLENOPARENT      39
#define MLENAMETAKEN     40
#define MLENOLISTEN      41
#define MLEBADNAME       42
#define MLEBADHOST       43
#define MLERESOURCE      44  /* a required resource was missing */
#define MLELAUNCHFAILED  45
#define MLELAUNCHNAME    46
#define MLELAST MLELAUNCHNAME /* for internal use only */

#define MLETRACEON      996  /* */
#define MLETRACEOFF     997  /* */
#define MLEDEBUG        998  /* */
#define MLEASSERT       999  /* an internal assertion failed */
#define MLEUSER        1000  /* start of user defined errors */


#endif /* _MLERRNO_H */


#endif /* _MLERRORS_H */

/* explicitly not protected by _MLERRORS_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( mlapi_error,   MLError,        ( MLINK mlp));
MLDECL( mlapi_result,  MLClearError,   ( MLINK mlp));
MLDECL( mlapi_result,  MLSetError,     ( MLINK mlp, mlapi_error err));
ML_END_EXTERN_C


#ifndef _MLYLDMSG_H
#define _MLYLDMSG_H



enum {	MLTerminateMessage = 1, MLInterruptMessage, MLAbortMessage,
	MLEndPacketMessage, MLSynchronizeMessage, MLImDyingMessage,
	MLWaitingAcknowledgment, MLMarkTopLevelMessage,
	MLFirstUserMessage = 128, MLLastUserMessage = 255 };

typedef unsigned long devinfo_selector;


#endif /* _MLYLDMSG_H */

/* explicitly not protected by _MLYLDMSG_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( mlapi_result,   MLPutMessage,   ( MLINK mlp, dev_message  msg));
MLDECL( mlapi_result,   MLMessageReady, ( MLINK mlp));
MLDECL( mlapi_result,   MLGetMessage,   ( MLINK mlp, dev_messagep mp, dev_messagep np));

MLDECL( MLMessageHandlerObject, MLMessageHandler,    ( MLINK mlp));
MLDECL( mlapi_result,           MLSetMessageHandler, ( MLINK mlp, MLMessageHandlerObject h));
MLDECL( MLYieldFunctionObject,  MLYieldFunction,     ( MLINK mlp));
MLDECL( mlapi_result,           MLSetYieldFunction,  ( MLINK mlp, MLYieldFunctionObject yf));
#ifndef MLINTERFACE
/* syntax error */ )
#endif
#if MLINTERFACE > 1
MLDECL( mlapi_result,  MLSetYieldFunction0,  ( MLINK mlp, MLYieldFunctionObject yf, MLINK cookie));
MLDECL( mlapi_result,  MLSetMessageHandler0, ( MLINK mlp, MLMessageHandlerObject func, MLINK cookie));
#endif


MLDECL( mlapi_result, MLDeviceInformation, ( MLINK mlp, devinfo_selector selector, MLPointer buf, longp_st buflen));
ML_END_EXTERN_C

/*************** Textual interface ***************/


#ifndef _MLGET_H
#define _MLGET_H


#endif /* _MLGET_H */

/* explicitly not protected by _MLGET_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( mlapi_token,    MLGetNext,          ( MLINK mlp));
MLDECL( mlapi_token,    MLGetNextRaw,       ( MLINK mlp));
MLDECL( mlapi_token,    MLGetType,          ( MLINK mlp));
MLDECL( mlapi_token,    MLGetRawType,       ( MLINK mlp));
MLDECL( mlapi_result,   MLGetRawData,       ( MLINK mlp, ucharp_ct data, long_st size, longp_st gotp));
MLDECL( mlapi_result,   MLGetData,          ( MLINK mlp, charp_ct data, long_st size, longp_st gotp));
MLDECL( mlapi_result,   MLGetArgCount,      ( MLINK mlp, longp_st countp));
MLDECL( mlapi_result,   MLGetRawArgCount,   ( MLINK mlp, longp_st countp));
MLDECL( mlapi_result,   MLBytesToGet,       ( MLINK mlp, longp_st leftp));
MLDECL( mlapi_result,   MLRawBytesToGet,    ( MLINK mlp, longp_st leftp));
MLDECL( mlapi_result,   MLExpressionsToGet, ( MLINK mlp, longp_st countp));
MLDECL( mlapi_result,   MLNewPacket,        ( MLINK mlp));
MLDECL( mlapi_result,   MLTakeLast,         ( MLINK mlp, long_st eleft));
MLDECL( mlapi_result,   MLReady,            ( MLINK mlp));
MLDECL( mlapi_result,   MLFill,             ( MLINK mlp));
ML_END_EXTERN_C


#ifndef _MLPUT_H
#define _MLPUT_H


#define MLPutExpression is obsolete, use MLPutComposite

#endif /* _MLPUT_H */

/* explicitly not protected by _MLPUT_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( mlapi_result,   MLPutNext,      ( MLINK mlp, mlapi_token tok));
MLDECL( mlapi_result,   MLPutType,      ( MLINK mlp, mlapi__token tok));
MLDECL( mlapi_result,   MLPutRawSize,   ( MLINK mlp, long_st size));
MLDECL( mlapi_result,   MLPutRawData,   ( MLINK mlp, kucharp_ct data, long_st len));
MLDECL( mlapi_result,   MLPutArgCount,  ( MLINK mlp, long_st argc));
MLDECL( mlapi_result,   MLPutComposite, ( MLINK mlp, long_st argc));
MLDECL( mlapi_result,   MLBytesToPut,   ( MLINK mlp, longp_st leftp));
MLDECL( mlapi_result,   MLEndPacket,    ( MLINK mlp));
MLDECL( mlapi_result,   MLFlush,        ( MLINK mlp));
ML_END_EXTERN_C


#ifndef _MLTK_H
#define _MLTK_H


#define	MLTKOLDINT     'I'		/* 73 Ox49 01001001 */ /* integer leaf node */
#define	MLTKOLDREAL    'R'		/* 82 Ox52 01010010 */ /* real leaf node */


#define	MLTKFUNC    'F'		/* 70 Ox46 01000110 */ /* non-leaf node */

#define	MLTKERROR   (0)		/* bad token */
#define	MLTKERR     (0)		/* bad token */

/* text token bit patterns: 0010x01x --exactly 2 bits worth chosen to make things somewhat readable */
#define MLTK__IS_TEXT( tok) ( (tok & 0x00F6) == 0x0022)

#define	MLTKSTR     '"'  /* 00100010 */
#define	MLTKSYM     '\043'  /* 00100011 */ /* octal here as hash requires a trigraph */

#define	MLTKREAL    '*'  /* 00101010 */
#define	MLTKINT     '+'  /* 00101011 */



/* The following defines are for internal use only */
#define	MLTKPCTEND  ']'     /* at end of top level expression */
#define	MLTKAPCTEND '\n'    /* at end of top level expression */
#define	MLTKEND     '\n'
#define	MLTKAEND    '\r'
#define	MLTKSEND    ','

#define	MLTKCONT    '\\'
#define	MLTKELEN    ' '

#define	MLTKNULL    '.'
#define	MLTKOLDSYM  'Y'     /* 89 0x59 01011001 */
#define	MLTKOLDSTR  'S'     /* 83 0x53 01010011 */


typedef unsigned long decoder_mask;
#define	MLTKPACKED	'P'     /* 80 0x50 01010000 */
#define	MLTKARRAY	'A'     /* 65 0x41 01000001 */
#define	MLTKDIM		'D'     /* 68 0x44 01000100 */

#define MLLENGTH_DECODER        ((decoder_mask) 1<<16)
#define MLTKPACKED_DECODER      ((decoder_mask) 1<<17)
#define MLTKARRAY_DECODER	    ((decoder_mask) 1<<18)
#define MLTKMODERNCHARS_DECODER ((decoder_mask) 1<<19)
#if 0
#define MLTKNULLSEQUENCE_DECODER ((decoder_mask) 1<<20)
#else
#define MLTKNULLSEQUENCE_DECODER ((decoder_mask) 0)
#endif
#define MLTKALL_DECODERS (MLLENGTH_DECODER | MLTKPACKED_DECODER | MLTKARRAY_DECODER | MLTKMODERNCHARS_DECODER | MLTKNULLSEQUENCE_DECODER)

#define MLTK_FIRSTUSER '\x30' /* user token */
#define MLTK_LASTUSER  '\x3F'



#endif /* _MLTK_H */

/*************** Native C types interface ***************/


#ifndef _MLCGET_H
#define _MLCGET_H



#define MLGetReal MLGetDouble

#endif /* _MLCGET_H */


/* explicitly not protected by _MLCGET_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( mlapi_result,   MLGetBinaryNumber, ( MLINK mlp, voidp_ct np, long type));
MLDECL( mlapi_result,   MLGetShortInteger, ( MLINK mlp, shortp_nt hp));
MLDECL( mlapi_result,   MLGetInteger,      ( MLINK mlp, intp_nt ip));
MLDECL( mlapi_result,   MLGetLongInteger,  ( MLINK mlp, longp_nt lp));
MLDECL( mlapi_result,   MLGetFloat,        ( MLINK mlp, floatp_nt fp));
MLDECL( mlapi_result,   MLGetDouble,       ( MLINK mlp, doublep_nt dp));
#if CC_SUPPORTS_LONG_DOUBLE
MLDECL( mlapi_result,   MLGetLongDouble,   ( MLINK mlp, extendedp_nt xp));
#endif

MLDECL( mlapi_result,   MLGet16BitCharacters,  ( MLINK mlp, longp_st chars_left, ushortp_ct buf, long_st cardof_buf, longp_st got));
MLDECL( mlapi_result,   MLGet8BitCharacters,   ( MLINK mlp, longp_st chars_left, ucharp_ct  buf, long_st cardof_buf, longp_st got, long missing));
MLDECL( mlapi_result,   MLGet7BitCharacters,   ( MLINK mlp, longp_st chars_left, charp_ct   buf, long_st cardof_buf, longp_st got));

MLDECL( mlapi_result,   MLGetUnicodeString,    ( MLINK mlp, kushortpp_ct sp, longp_st lenp));
MLDECL( mlapi_result,   MLGetByteString,       ( MLINK mlp, kucharpp_ct  sp, longp_st lenp, long missing));
MLDECL( mlapi_result,   MLGetString,           ( MLINK mlp, kcharpp_ct   sp));
#ifndef MLINTERFACE
/* syntax error */ )
#endif
#if MLINTERFACE > 1
MLDECL( mlapi_result,   MLGetUnicodeString0,   ( MLINK mlp, kushortpp_ct sp, longp_st lenp));
MLDECL( mlapi_result,   MLGetByteString0,      ( MLINK mlp, kucharpp_ct  sp, longp_st lenp, long missing));
MLDECL( mlapi_result,   MLGetString0,          ( MLINK mlp, kcharpp_ct   sp));
#endif

MLDECL( mlapi_result,   MLGetUnicodeSymbol,    ( MLINK mlp, kushortpp_ct sp, longp_st lenp));
MLDECL( mlapi_result,   MLGetByteSymbol,       ( MLINK mlp, kucharpp_ct  sp, longp_st lenp, long missing));
MLDECL( mlapi_result,   MLGetSymbol,           ( MLINK mlp, kcharpp_ct   sp));

MLDECL( void,           MLDisownUnicodeString, ( MLINK mlp, kushortp_ct s,   long_st len));
MLDECL( void,           MLDisownByteString,    ( MLINK mlp, kucharp_ct  s,   long_st len));
MLDECL( void,           MLDisownString,        ( MLINK mlp, kcharp_ct   s));

MLDECL( void,           MLDisownUnicodeSymbol, ( MLINK mlp, kushortp_ct s,   long_st len));
MLDECL( void,           MLDisownByteSymbol,    ( MLINK mlp, kucharp_ct  s,   long_st len));
MLDECL( void,           MLDisownSymbol,        ( MLINK mlp, kcharp_ct   s));



MLDECL( mlapi_result,   MLCheckString,   ( MLINK mlp, kcharp_ct name));
MLDECL( mlapi_result,   MLCheckSymbol,   ( MLINK mlp, kcharp_ct name));
MLDECL( mlapi_result,   MLGetFunction,   ( MLINK mlp, kcharpp_ct sp, longp_st countp));
MLDECL( mlapi_result,   MLCheckFunction, ( MLINK mlp, kcharp_ct s, longp_st countp));
MLDECL( mlapi_result,   MLCheckFunctionWithArgCount, ( MLINK mlp, kcharp_ct s, longp_st countp));
ML_END_EXTERN_C


#ifndef _MLCPUT_H
#define _MLCPUT_H


#define MLPutReal MLPutDouble

#endif /* _MLCPUT_H */

/* explicitly not protected by _MLCPUT_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( mlapi_result,   MLPutBinaryNumber, ( MLINK mlp, voidp_ct np, long type));
MLDECL( mlapi_result,   MLPutShortInteger, ( MLINK mlp, int_nt h));
MLDECL( mlapi_result,   MLPutInteger,      ( MLINK mlp, int_nt i));
MLDECL( mlapi_result,   MLPutLongInteger,  ( MLINK mlp, long_nt l));
MLDECL( mlapi_result,   MLPutFloat,        ( MLINK mlp, double_nt f));
MLDECL( mlapi_result,   MLPutDouble,       ( MLINK mlp, double_nt d));
#if CC_SUPPORTS_LONG_DOUBLE
MLDECL( mlapi_result,   MLPutLongDouble,   ( MLINK mlp, extended_nt x));
#endif

MLDECL( mlapi_result,   MLPut16BitCharacters, ( MLINK mlp, long_st chars_left, kushortp_ct codes, long_st ncodes));
MLDECL( mlapi_result,   MLPut8BitCharacters,  ( MLINK mlp, long_st chars_left, kucharp_ct bytes, long_st nbytes));
MLDECL( mlapi_result,   MLPut7BitCount,       ( MLINK mlp, long_st count, long_st size));
MLDECL( mlapi_result,   MLPut7BitCharacters,  ( MLINK mlp, long_st chars_left, kcharp_ct bytes, long_st nbytes, long_st nchars_now));

MLDECL( mlapi_result,   MLPutUnicodeString, ( MLINK mlp, kushortp_ct s, long_st len));
MLDECL( mlapi_result,   MLPutByteString,    ( MLINK mlp, kucharp_ct  s, long_st len));
MLDECL( mlapi_result,   MLPutString,        ( MLINK mlp, kcharp_ct   s));
#ifndef MLINTERFACE
/* syntax error */ )
#endif
#if MLINTERFACE > 1
MLDECL( mlapi_result,   MLPutRealUnicodeString0, ( MLINK mlp, ushortp_ct s));
MLDECL( mlapi_result,   MLPutRealByteString0,    ( MLINK mlp, ucharp_ct  s));
#endif

MLDECL( mlapi_result,   MLPutUnicodeSymbol, ( MLINK mlp, kushortp_ct s, long_st len));
MLDECL( mlapi_result,   MLPutByteSymbol,    ( MLINK mlp, kucharp_ct  s, long_st len));
MLDECL( mlapi_result,   MLPutSymbol,        ( MLINK mlp, kcharp_ct   s));

MLDECL( mlapi_result,   MLPutFunction,      ( MLINK mlp, kcharp_ct s, long_st argc));


MLDECL( mlapi_result,   MLPutSize, ( MLINK mlp, long_st size));
MLDECL( mlapi_result,   MLPutData, ( MLINK mlp, kcharp_ct buff, long_st len));
ML_END_EXTERN_C



#ifndef _MLSTRING_H
#define _MLSTRING_H



#define MAX_BYTES_PER_OLD_CHARACTER 3
#define MAX_BYTES_PER_NEW_CHARACTER 6

#define ML_MAX_BYTES_PER_CHARACTER MAX_BYTES_PER_NEW_CHARACTER

/* for source code compatibility with earlier versions of MathLink */

#if POWERMACINTOSH_MATHLINK
#pragma options align=mac68k
#endif

typedef struct {
	kcharp_ct str;
	kcharp_ct end;
} MLStringPosition;

#if POWERMACINTOSH_MATHLINK
#pragma options align=reset
#endif

typedef MLStringPosition FAR * MLStringPositionPointer;

#define MLStringFirstPos(s,pos) MLStringFirstPosFun( s, &(pos))

#define MLforString( s, pos) \
	for( MLStringFirstPos(s,pos); MLStringCharacter( (pos).str, (pos).end) >= 0; MLNextCharacter(&(pos).str, (pos).end))

#define MLStringChar( pos) MLStringCharacter( (pos).str, (pos).end)

#define MLPutCharToString MLConvertCharacter


/* for internal use only */

#if POWERMACINTOSH_MATHLINK
#pragma options align=mac68k
#endif

typedef struct {
	ucharp_ct cc;
	int_ct  mode;
	int_ct  more;
	ucharp_ct head;
} MLOldStringPosition;

#if POWERMACINTOSH_MATHLINK
#pragma options align=reset
#endif

typedef MLOldStringPosition FAR * MLOldStringPositionPointer;


#define MLOldforString( s, pos) \
  for ( MLOldStringFirstPos( s, pos); (pos).more; MLOldStringNextPos( pos))

#define MLOldStringChar(pos) \
  ( ((pos).mode <= 1) ? (uint_ct)(*(ucharp_ct)((pos).cc)) : MLOldStringCharFun( &pos) )


#define MLOldStringFirstPos(s,pos) MLOldStringFirstPosFun( s, &(pos))

#define MLOldStringNextPos(pos)  ( \
	((pos).mode == 0) \
		? ((*(*(pos).cc ? ++(pos).cc : (pos).cc) ? 0 : ((pos).more = 0)), (pos).cc) \
		: MLOldStringNextPosFun( &pos) )





#endif /* _MLSTRING_H */




/* explicitly not protected by _MLXDATA_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
/* assumes *startp aligned on char boundary, if n == -1, returns ~(char_count) */
MLDECL( long, MLCharacterOffset,           ( kcharpp_ct startp, kcharp_ct end, long n));
MLDECL( long, MLStringCharacter,           ( kcharp_ct  start,  kcharp_ct end));
MLDECL( long, MLNextCharacter,             ( kcharpp_ct startp, kcharp_ct end));
#ifndef MLINTERFACE
/* syntax error */ )
#endif
#if MLINTERFACE > 1
MLDECL( long, MLNextCharacter0,            ( kcharp_ct str, longp_ct indexp, long len));
#endif

MLDECL( long, MLConvertNewLine,            ( charpp_ct sp));
MLDECL( long, MLConvertCharacter,          ( ulong_ct ch, charpp_ct sp));
MLDECL( long, MLConvertByteString,         ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end));
MLDECL( long, MLConvertByteStringNL,       ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end, ulong_ct nl));
MLDECL( long, MLConvertUnicodeString,      ( ushortp_ct codes, long len, charpp_ct strp, charp_ct str_end));
MLDECL( long, MLConvertUnicodeStringNL,    ( ushortp_ct codes, long len, charpp_ct strp, charp_ct str_end, ulong_ct nl));
MLDECL( long, MLConvertDoubleByteString,   ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end));
MLDECL( long, MLConvertDoubleByteStringNL, ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end, ulong_ct nl));







/* for source code compatibility with earlier versions of MathLink */
MLDECL( kcharp_ct,     MLStringFirstPosFun,  ( kcharp_ct s, MLStringPositionPointer p));

/* for internal use only */
MLDECL( mlapi_result, MLOldPutCharToString,      ( uint_ct ch, charpp_ct sp));
MLDECL( ucharp_ct,    MLOldStringNextPosFun,     ( MLOldStringPositionPointer p));
MLDECL( ucharp_ct,    MLOldStringFirstPosFun,    ( charp_ct s, MLOldStringPositionPointer p));
MLDECL( uint_ct,      MLOldStringCharFun,        ( MLOldStringPositionPointer p));
MLDECL( long,         MLOldConvertByteString,    ( ucharp_ct  codes, long len, charpp_ct strp, charp_ct str_end));
MLDECL( long,         MLOldConvertUnicodeString, ( ushortp_ct codes, long len, charpp_ct strp, charp_ct str_end));

ML_END_EXTERN_C


#ifndef _MLCAPUT_H
#define _MLCAPUT_H

#ifndef MLINTERFACE
/* syntax error */ )
#endif

#ifndef __array_meterp__
#define __array_meterp__
typedef struct array_meter FAR * array_meterp;
typedef array_meterp FAR * array_meterpp;
#endif


#define MLPutRealArray MLPutDoubleArray

#endif /* _MLCAPUT_H */


/* explicitly not protected by _MLCAPUT_H in case MLDECL is redefined for multiple inclusion */

/*bugcheck: bugcheck need FAR here */
ML_EXTERN_C
MLDECL( mlapi_result,   MLPutArray,                 ( MLINK mlp, array_meterp meterp));
MLDECL( mlapi_result,   MLPutBinaryNumberArrayData, ( MLINK mlp, array_meterp meterp, voidp_ct     datap, long_st count, long type));
MLDECL( mlapi_result,   MLPutByteArrayData,         ( MLINK mlp, array_meterp meterp, ucharp_nt    datap, long_st count));
MLDECL( mlapi_result,   MLPutShortIntegerArrayData, ( MLINK mlp, array_meterp meterp, shortp_nt    datap, long_st count));
MLDECL( mlapi_result,   MLPutIntegerArrayData,      ( MLINK mlp, array_meterp meterp, intp_nt      datap, long_st count));
MLDECL( mlapi_result,   MLPutLongIntegerArrayData,  ( MLINK mlp, array_meterp meterp, longp_nt     datap, long_st count));
MLDECL( mlapi_result,   MLPutFloatArrayData,        ( MLINK mlp, array_meterp meterp, floatp_nt    datap, long_st count));
MLDECL( mlapi_result,   MLPutDoubleArrayData,       ( MLINK mlp, array_meterp meterp, doublep_nt   datap, long_st count));
#if CC_SUPPORTS_LONG_DOUBLE
MLDECL( mlapi_result,   MLPutLongDoubleArrayData,   ( MLINK mlp, array_meterp meterp, extendedp_nt datap, long_st count));
#endif

MLDECL( mlapi_result,   MLPutBinaryNumberArray, ( MLINK mlp, voidp_ct     data, longp_st dimp, charpp_ct heads, long_st depth, long type));
MLDECL( mlapi_result,   MLPutByteArray,         ( MLINK mlp, ucharp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( mlapi_result,   MLPutShortIntegerArray, ( MLINK mlp, shortp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( mlapi_result,   MLPutIntegerArray,      ( MLINK mlp, intp_nt      data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( mlapi_result,   MLPutLongIntegerArray,  ( MLINK mlp, longp_nt     data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( mlapi_result,   MLPutDoubleArray,       ( MLINK mlp, doublep_nt   data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( mlapi_result,   MLPutFloatArray,        ( MLINK mlp, floatp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
#if CC_SUPPORTS_LONG_DOUBLE
MLDECL( mlapi_result,   MLPutLongDoubleArray,   ( MLINK mlp, extendedp_nt data, longp_st dims, charpp_ct heads, long_st depth));
#endif

MLDECL( mlapi_result,   MLPutBinaryNumberList, ( MLINK mlp, voidp_ct   data, long_st count, long type));
MLDECL( mlapi_result,   MLPutIntegerList,      ( MLINK mlp, intp_nt    data, long_st count));
MLDECL( mlapi_result,   MLPutRealList,         ( MLINK mlp, doublep_nt data, long_st count));

#if MLINTERFACE > 1
MLDECL( mlapi_result, MLPutArrayType0,             ( MLINK mlp, MLINK heads, long depth, array_meterpp meterpp));
MLDECL( mlapi_result, MLPutBinaryNumberArrayData0, ( MLINK mlp, MLINK heads, array_meterp meterp, voidp_ct datap, long_st count, long type));
MLDECL( mlapi_result, MLReleasePutArrayState0,     ( MLINK mlp, MLINK heads, array_meterp meterp));
MLDECL( mlapi_result, MLPutArrayLeaves0,           ( MLINK mlp, MLINK heads, array_meterp meterp, MLINK leaves, long_st count));
#endif

ML_END_EXTERN_C


#ifndef _MLCAGET_H
#define _MLCAGET_H

#ifndef MLINTERFACE
/* syntax error */ )
#endif

#ifndef __array_meterp__
#define __array_meterp__
typedef struct array_meter FAR * array_meterp;
typedef array_meterp FAR * array_meterpp;
#endif


#define MLGetRealArray    MLGetDoubleArray
#define MLDisownRealArray MLDisownDoubleArray

#endif /* _MLCAGET_H */



/* explicitly not protected by _MLCAGET_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( mlapi_result,  MLGetBinaryNumberList, ( MLINK mlp, voidpp_ct   datap, longp_st countp, long type));
MLDECL( mlapi_result,  MLGetIntegerList,      ( MLINK mlp, intpp_nt    datap, longp_st countp));
MLDECL( mlapi_result,  MLGetRealList,         ( MLINK mlp, doublepp_nt datap, longp_st countp));

MLDECL( void, MLDisownBinaryNumberList, ( MLINK mlp, voidp_ct   data, long_st count, long type));
MLDECL( void, MLDisownIntegerList,      ( MLINK mlp, intp_nt    data, long_st count));
MLDECL( void, MLDisownRealList,         ( MLINK mlp, doublep_nt data, long_st count));

MLDECL( mlapi_token,    MLGetArrayType,             ( MLINK mlp, array_meterp meterp));
MLDECL( mlapi_result,   MLGetArrayDimensions,       ( MLINK mlp, array_meterp meterp));

MLDECL( mlapi_result,   MLGetBinaryNumberArrayData, ( MLINK mlp, array_meterp meterp, voidp_ct     datap, long_st count, long type));
MLDECL( mlapi_result,   MLGetByteArrayData,         ( MLINK mlp, array_meterp meterp, ucharp_nt    datap, long_st count));
MLDECL( mlapi_result,   MLGetShortIntegerArrayData, ( MLINK mlp, array_meterp meterp, shortp_nt    datap, long_st count));
MLDECL( mlapi_result,   MLGetIntegerArrayData,      ( MLINK mlp, array_meterp meterp, intp_nt      datap, long_st count));
MLDECL( mlapi_result,   MLGetLongIntegerArrayData,  ( MLINK mlp, array_meterp meterp, longp_nt     datap, long_st count));
MLDECL( mlapi_result,   MLGetFloatArrayData,        ( MLINK mlp, array_meterp meterp, floatp_nt    datap, long_st count));
MLDECL( mlapi_result,   MLGetDoubleArrayData,       ( MLINK mlp, array_meterp meterp, doublep_nt   datap, long_st count));
#if CC_SUPPORTS_LONG_DOUBLE
MLDECL( mlapi_result,   MLGetLongDoubleArrayData,   ( MLINK mlp, array_meterp meterp, extendedp_nt datap, long_st count));
#endif

#if MLINTERFACE > 1
MLDECL( mlapi_result,   MLGetArrayType0,             ( MLINK mlp, MLINK heads, array_meterpp meterpp, longp_st depthp, mlapi__tokenp leaf_tokp));
MLDECL( mlapi_result,   MLGetBinaryNumberArrayData0, ( MLINK mlp, MLINK heads, array_meterp  meterp, voidp_ct datap, longp_st countp, long type));
MLDECL( void,           MLReleaseGetArrayState0,     ( MLINK mlp, MLINK heads, array_meterp  meterp));

MLDECL( mlapi_result,   MLGetBinaryNumberArray0,   ( MLINK mlp, voidpp_ct     datap, longpp_st dimpp, charppp_ct headsp, longp_st depthp, long type, mlapi__tokenp leaf_tokp));
#endif
MLDECL( mlapi_result,   MLGetBinaryNumberArray,    ( MLINK mlp, voidpp_ct     datap, longpp_st dimpp, charppp_ct headsp, longp_st depthp, long type));
MLDECL( mlapi_result,   MLGetByteArray,            ( MLINK mlp, ucharpp_nt    datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
MLDECL( mlapi_result,   MLGetShortIntegerArray,    ( MLINK mlp, shortpp_nt    datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
MLDECL( mlapi_result,   MLGetIntegerArray,         ( MLINK mlp, intpp_nt      datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
MLDECL( mlapi_result,   MLGetLongIntegerArray,     ( MLINK mlp, longpp_nt     datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
MLDECL( mlapi_result,   MLGetDoubleArray,          ( MLINK mlp, doublepp_nt   datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
MLDECL( mlapi_result,   MLGetFloatArray,           ( MLINK mlp, floatpp_nt    datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
#if CC_SUPPORTS_LONG_DOUBLE
MLDECL( mlapi_result,   MLGetLongDoubleArray,      ( MLINK mlp, extendedpp_nt datap, longpp_st dimsp, charppp_ct headsp, longp_st depthp));
#endif

MLDECL( void,           MLDisownBinaryNumberArray, ( MLINK mlp, voidp_ct     data, longp_st dimp, charpp_ct heads, long_st len, long type));
MLDECL( void,           MLDisownByteArray,         ( MLINK mlp, ucharp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( void,           MLDisownShortIntegerArray, ( MLINK mlp, shortp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( void,           MLDisownIntegerArray,      ( MLINK mlp, intp_nt      data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( void,           MLDisownLongIntegerArray,  ( MLINK mlp, longp_nt     data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( void,           MLDisownFloatArray,        ( MLINK mlp, floatp_nt    data, longp_st dims, charpp_ct heads, long_st depth));
MLDECL( void,           MLDisownDoubleArray,       ( MLINK mlp, doublep_nt   data, longp_st dims, charpp_ct heads, long_st depth));
#if CC_SUPPORTS_LONG_DOUBLE
MLDECL( void,           MLDisownLongDoubleArray,   ( MLINK mlp, extendedp_nt data, longp_st dims, charpp_ct heads, long_st depth));
#endif
ML_END_EXTERN_C


/*************** seeking, transfering  and synchronization ***************/

#ifndef _MLMARK_H
#define _MLMARK_H


#endif /* _MLMARK_H */

/* explicitly not protected by _MLMARK_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( MLINKMark,  MLCreateMark,  ( MLINK mlp));
MLDECL( MLINKMark,  MLSeekToMark,  ( MLINK mlp, MLINKMark mark, long index));
MLDECL( MLINKMark,  MLSeekMark,    ( MLINK mlp, MLINKMark mark, long index));
MLDECL( void,       MLDestroyMark, ( MLINK mlp, MLINKMark mark));
ML_END_EXTERN_C


#ifndef _MLXFER_H
#define _MLXFER_H


#endif /* _MLXFER_H */

/* explicitly not protected by _MLXFER_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( mlapi_result, MLTransferExpression, ( MLINK dmlp, MLINK smlp));
MLDECL( mlapi_result, MLTransferToEndOfLoopbackLink, ( MLINK dmlp, MLINK smlp));
#ifndef MLINTERFACE
/* syntax error */ )
#endif
#if MLINTERFACE > 1
MLDECL( mlapi_result, MLTransfer0, ( MLINK dmlp, MLINK smlp, ulong_ct sequence_no));
#endif
ML_END_EXTERN_C


#ifndef _MLSYNC_H
#define _MLSYNC_H


/* export mls__wait and mls__align(mlsp) */

#endif /* _MLSYNC_H */

/* explicitly not protected by _MLSYNC_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
/* in response to a reset message */
MLDECL( mlapi_result, MLForwardReset, ( MLINK mlp, ulong_ct marker));
MLDECL( mlapi_result, MLAlign,        ( MLINK lmlp, MLINK rmlp));
ML_END_EXTERN_C

/*************************************************************/


#ifndef _MLPKT_H
#define _MLPKT_H

/*************** Mathematica packet interface ***************/

			/* MLNextPacket returns one of... */


/* edit here and in mlpktstr.h */

#ifndef _MLPKTNO_H
#define _MLPKTNO_H

#define ILLEGALPKT      0

#define CALLPKT         7
#define EVALUATEPKT    13
#define RETURNPKT       3

#define INPUTNAMEPKT    8
#define ENTERTEXTPKT   14
#define ENTEREXPRPKT   15
#define OUTPUTNAMEPKT   9
#define RETURNTEXTPKT   4
#define RETURNEXPRPKT  16

#define DISPLAYPKT     11
#define DISPLAYENDPKT  12

#define MESSAGEPKT      5
#define TEXTPKT         2

#define INPUTPKT        1
#define INPUTSTRPKT    21
#define MENUPKT         6
#define SYNTAXPKT      10

#define SUSPENDPKT     17
#define RESUMEPKT      18

#define BEGINDLGPKT    19
#define ENDDLGPKT      20

#define FIRSTUSERPKT  128
#define LASTUSERPKT   255


#endif /* _MLPKTNO_H */

#endif /* _MLPKT_H */

/* explicitly not protected by _MLPKT_H in case MLDECL is redefined for multiple inclusion */

ML_EXTERN_C
MLDECL( mlapi_packet,  MLNextPacket, ( MLINK mlp));
ML_END_EXTERN_C


#ifndef _MLALERT_H
#define _MLALERT_H



ML_EXTERN_C
/*************** User interaction--for internal use only ***************/
typedef long mldlg_result;

MLDPROC( mldlg_result, MLAlertProcPtr,             ( MLEnvironment env, kcharp_ct message));
MLDPROC( mldlg_result, MLRequestProcPtr,           ( MLEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response));
MLDPROC( mldlg_result, MLConfirmProcPtr,           ( MLEnvironment env, kcharp_ct question, mldlg_result default_answer));
MLDPROC( mldlg_result, MLRequestArgvProcPtr,       ( MLEnvironment env, charpp_ct argv, long cardof_argv, charp_ct buf, long sizeof_buf));
MLDPROC( mldlg_result, MLRequestToInteractProcPtr, ( MLEnvironment env, mldlg_result wait_for_permission));
MLDPROC( mldlg_result, MLDialogProcPtr,            ( MLEnvironment env));

enum {
	uppMLAlertFunctionProcInfo = kPascalStackBased
		 | RESULT_SIZE(SIZE_CODE(sizeof(mldlg_result)))
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(MLEnvironment)))
		 | STACK_ROUTINE_PARAMETER(2, SIZE_CODE(sizeof(kcharp_ct))),
	uppMLRequestFunctionProcInfo = kPascalStackBased
		 | RESULT_SIZE(SIZE_CODE(sizeof(mldlg_result)))
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(MLEnvironment)))
		 | STACK_ROUTINE_PARAMETER(2, SIZE_CODE(sizeof(kcharp_ct)))
		 | STACK_ROUTINE_PARAMETER(3, SIZE_CODE(sizeof(charp_ct)))
		 | STACK_ROUTINE_PARAMETER(4, SIZE_CODE(sizeof(long))),
	uppMLConfirmFunctionProcInfo = kPascalStackBased
		 | RESULT_SIZE(SIZE_CODE(sizeof(mldlg_result)))
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(MLEnvironment)))
		 | STACK_ROUTINE_PARAMETER(2, SIZE_CODE(sizeof(kcharp_ct)))
		 | STACK_ROUTINE_PARAMETER(3, SIZE_CODE(sizeof(mldlg_result))),
	uppMLRequestArgvFunctionProcInfo = kPascalStackBased
		 | RESULT_SIZE(SIZE_CODE(sizeof(mldlg_result)))
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(MLEnvironment)))
		 | STACK_ROUTINE_PARAMETER(2, SIZE_CODE(sizeof(charpp_ct)))
		 | STACK_ROUTINE_PARAMETER(3, SIZE_CODE(sizeof(long)))
		 | STACK_ROUTINE_PARAMETER(4, SIZE_CODE(sizeof(charp_ct)))
		 | STACK_ROUTINE_PARAMETER(5, SIZE_CODE(sizeof(long))),
	uppMLRequestToInteractFunctionProcInfo = kPascalStackBased
		 | RESULT_SIZE(SIZE_CODE(sizeof(mldlg_result)))
		 | STACK_ROUTINE_PARAMETER(1, SIZE_CODE(sizeof(MLEnvironment)))
		 | STACK_ROUTINE_PARAMETER(2, SIZE_CODE(sizeof(mldlg_result)))
};



#if GENERATINGCFM

	typedef UniversalProcPtr MLDialogUPP;
	typedef UniversalProcPtr MLAlertUPP;
	typedef UniversalProcPtr MLRequestUPP;
	typedef UniversalProcPtr MLConfirmUPP;
	typedef UniversalProcPtr MLRequestArgvUPP;
	typedef UniversalProcPtr MLRequestToInteractUPP;

#	define NewMLAlertProc(userRoutine) \
		NewRoutineDescriptor((ProcPtr)MLAlertCast((userRoutine)), \
			uppMLAlertFunctionProcInfo, GetCurrentArchitecture())
#	define NewMLRequestProc(userRoutine) \
		NewRoutineDescriptor((ProcPtr)MLRequestCast((userRoutine)), \
			uppMLRequestFunctionProcInfo, GetCurrentArchitecture())
#	define NewMLConfirmProc(userRoutine) \
		NewRoutineDescriptor((ProcPtr)MLConfirmCast((userRoutine)), \
			uppMLConfirmFunctionProcInfo, GetCurrentArchitecture())
#	define NewMLRequestArgvProc(userRoutine) \
		NewRoutineDescriptor((ProcPtr)MLRequestArgvCast((userRoutine)), \
			uppMLRequestArgvFunctionProcInfo, GetCurrentArchitecture())
#	define NewMLRequestToInteractProc(userRoutine) \
		NewRoutineDescriptor((ProcPtr)MLRequestToInteractCast((userRoutine)), \
			uppMLRequestToInteractFunctionProcInfo, GetCurrentArchitecture())

#else

	typedef MLDialogProcPtr MLDialogUPP;
	typedef MLAlertProcPtr MLAlertUPP;
	typedef MLRequestProcPtr MLRequestUPP;
	typedef MLConfirmProcPtr MLConfirmUPP;
	typedef MLRequestArgvProcPtr MLRequestArgvUPP;
	typedef MLRequestToInteractProcPtr MLRequestToInteractUPP;

#	define NewMLAlertProc(userRoutine) MLAlertCast((userRoutine))
#	define NewMLRequestProc(userRoutine) MLRequestCast((userRoutine))
#	define NewMLConfirmProc(userRoutine) MLConfirmCast((userRoutine))
#	define NewMLRequestArgvProc(userRoutine) MLRequestArgvCast((userRoutine))
#	define NewMLRequestToInteractProc(userRoutine) MLRequestToInteractCast((userRoutine))

#endif


typedef MLAlertUPP MLAlertFunctionType;
typedef MLRequestUPP MLRequestFunctionType;
typedef MLConfirmUPP MLConfirmFunctionType;
typedef MLRequestArgvUPP MLRequestArgvFunctionType;
typedef MLRequestToInteractUPP MLRequestToInteractFunctionType;
typedef MLDialogUPP MLDialogFunctionType;



/* 
	MLDDECL( mldlg_result, alert_user, ( MLEnvironment env, kcharp_ct message));
	MLDDEFN( mldlg_result, alert_user, ( MLEnvironment env, kcharp_ct message))
	{
		fprintf( stderr, "%s\n", message);
	}


	...
	MLDialogFunctionType f = NewMLAlertProc(alert_user);
	MLSetDialogFunction( ep, MLAlertFunction, f);
	...
	or
	...
	MLSetDialogFunction( ep, MLAlertFunction, NewMLAlertProc(alert_user));
	...
*/



enum {	MLAlertFunction = 1, MLRequestFunction, MLConfirmFunction,
	MLRequestArgvFunction, MLRequestToInteractFunction };


#define ML_DEFAULT_DIALOG ( (MLDialogFunctionType) 1)
#define ML_IGNORE_DIALOG ( (MLDialogFunctionType) 0)
#define ML_SUPPRESS_DIALOG ML_IGNORE_DIALOG




#if MACINTOSH_MATHLINK

#ifndef _MLMAC_H
#define _MLMAC_H


ML_EXTERN_C

MLDDECL( mldlg_result, MLPermit_application,  ( MLEnvironment env, mldlg_result wait_for_permission));

MLDDECL( mldlg_result, MLAlert_application,   ( MLEnvironment env, kcharp_ct message));
MLDDECL( mldlg_result, MLAlert_tool,          ( MLEnvironment env, kcharp_ct message));
MLDDECL( mldlg_result, MLAlert_siow,          ( MLEnvironment env, kcharp_ct message));
MLDDECL( mldlg_result, MLAlert_console,       ( MLEnvironment env, kcharp_ct message));

MLDDECL( mldlg_result, MLRequest_application, ( MLEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response));
MLDDECL( mldlg_result, MLRequest_console,     ( MLEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response));
MLDDECL( mldlg_result, MLRequest_tool,        ( MLEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response));
MLDDECL( mldlg_result, MLRequest_siow,        ( MLEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response));

MLDDECL( mldlg_result, MLConfirm_application, ( MLEnvironment env, kcharp_ct question, mldlg_result default_answer));
MLDDECL( mldlg_result, MLConfirm_tool,        ( MLEnvironment env, kcharp_ct question, mldlg_result default_answer));
MLDDECL( mldlg_result, MLConfirm_siow,        ( MLEnvironment env, kcharp_ct question, mldlg_result default_answer));
MLDDECL( mldlg_result, MLConfirm_console,     ( MLEnvironment env, kcharp_ct question, mldlg_result default_answer));

ML_END_EXTERN_C

#endif /* _MLMAC_H */
#define MLALERT  	MLAlert_application
#define MLREQUEST	MLRequest_application
#define MLCONFIRM	MLConfirm_application
#define MLPERMIT 	MLPermit_application
#define MLREQUESTARGV	default_request_argv
#endif

#if WINDOWS_MATHLINK

#ifndef _MLWIN_H
#define _MLWIN_H



ML_EXTERN_C
MLDDECL( mldlg_result, MLAlert_win,   ( MLEnvironment ep, kcharp_ct alertstr));
MLDDECL( mldlg_result, MLRequest_win, ( MLEnvironment ep, kcharp_ct prompt, charp_ct response, long n));
MLDDECL( mldlg_result, MLConfirm_win, ( MLEnvironment ep, kcharp_ct okcancelquest, mldlg_result default_answer));
MLDDECL( mldlg_result, MLPermit_win,  ( MLEnvironment ep, mldlg_result wait));
ML_END_EXTERN_C

/* edit here and in mlwin.rc -- in both places because of command-line length limitations in dos */
#define DLG_LINKNAME                101
#define DLG_TEXT                    102
#define RIDOK                       1
#define RIDCANCEL                   104

#endif /* _MLWIN_H */
#define MLALERT         MLAlert_win
#define MLREQUEST       MLRequest_win
#define MLCONFIRM       MLConfirm_win
#define MLPERMIT        MLPermit_win
#define MLREQUESTARGV	default_request_argv
#endif

#if UNIX_MATHLINK

#ifndef _MLUNIX_H
#define _MLUNIX_H


ML_EXTERN_C

MLDDECL( mldlg_result, MLAlert_unix,   ( MLEnvironment env, kcharp_ct message));
MLDDECL( mldlg_result, MLRequest_unix, ( MLEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response));
MLDDECL( mldlg_result, MLConfirm_unix, ( MLEnvironment env, kcharp_ct question, mldlg_result default_answer));
MLDDECL( mldlg_result, MLPermit_unix,  ( MLEnvironment env, mldlg_result wait_for_permission));

ML_END_EXTERN_C

#endif /* _MLUNIX_H */
#define MLALERT  	MLAlert_unix
#define MLREQUEST	MLRequest_unix
#define MLCONFIRM	MLConfirm_unix
#define MLPERMIT 	MLPermit_unix
#define MLREQUESTARGV	default_request_argv
#endif

#if OS2_MATHLINK

#ifndef _MLOS2_H
#define _MLOS2_H



ML_EXTERN_C
MLDDECL( mldlg_result, MLAlert_os2,   ( MLEnvironment ep, kcharp_ct alertstr));
MLDDECL( mldlg_result, MLRequest_os2, ( MLEnvironment ep, kcharp_ct prompt, charp_ct response, long n));
MLDDECL( mldlg_result, MLConfirm_os2, ( MLEnvironment ep, kcharp_ct okcancelquest, mldlg_result default_answer));
MLDDECL( mldlg_result, MLPermit_os2,  ( MLEnvironment ep, mldlg_result wait));
ML_END_EXTERN_C

#endif /* _MLOS2_H */
#define MLALERT         MLAlert_os2
#define MLREQUEST       MLRequest_os2
#define MLCONFIRM       MLConfirm_os2
#define MLPERMIT        MLPermit_os2
#define MLREQUESTARGV	default_request_argv
#endif

MLDDECL( mldlg_result, default_request_argv, ( MLEnvironment ep, charpp_ct argv, long len, charp_ct buff, long size));
ML_END_EXTERN_C

#endif /* _MLALERT_H */


/* explicitly not protected by _MLXDATA_H in case MLDECL is redefined for multiple inclusion */
ML_EXTERN_C
MLDECL( mldlg_result,  MLAlert,             ( MLEnvironment env, kcharp_ct message));
MLDECL( mldlg_result,  MLRequest,           ( MLEnvironment env, kcharp_ct prompt, charp_ct response, long sizeof_response)); /* initialize response with default*/
MLDECL( mldlg_result,  MLConfirm,           ( MLEnvironment env, kcharp_ct question, mldlg_result default_answer));
MLDECL( mldlg_result,  MLRequestArgv,       ( MLEnvironment env, charpp_ct argv, long cardof_argv, charp_ct buff, long size));
MLDECL( mldlg_result,  MLRequestToInteract, ( MLEnvironment env, mldlg_result wait_for_permission));
MLDECL( mlapi_result,  MLSetDialogFunction, ( MLEnvironment env, long funcnum, MLDialogFunctionType func));

/* just some type-safe casts */
MLDECL( MLDialogProcPtr, MLAlertCast, ( MLAlertProcPtr f));
MLDECL( MLDialogProcPtr, MLRequestCast, ( MLRequestProcPtr f));
MLDECL( MLDialogProcPtr, MLConfirmCast, ( MLConfirmProcPtr f));
MLDECL( MLDialogProcPtr, MLRequestArgvCast, ( MLRequestArgvProcPtr f));
MLDECL( MLDialogProcPtr, MLRequestToInteractCast, ( MLRequestToInteractProcPtr f));
ML_END_EXTERN_C

/*************************************************************/

#ifdef __CFM68K__
#pragma import off
#endif


#ifndef _MLTM_H
#define _MLTM_H


/*************** Template interface ***************/

/* The following are useful only when using template files as
 * their definitions are produced by mprep.
 */

extern MLINK stdlink;
extern MLEnvironment stdenv;
extern MLYieldFunctionObject stdyielder;
extern MLMessageHandlerObject stdhandler;
extern int	MLMain P((int, charpp_ct)); /* pass in argc and argv */
extern int  MLMainString P(( charp_ct commandline));
extern int  MLMainArgv P(( char** argv, char** argv_end)); /* note not FAR pointers */
            
extern int	MLInstall P((MLINK));
extern mlapi_packet	MLAnswer P((MLINK));
extern int	MLDoCallPacket P((MLINK));
extern int	MLEvaluate P(( MLINK, charp_ct));
extern int	MLEvaluateString P(( MLINK, charp_ct));
MLMDECL( void, MLDefaultHandler, ( MLINK, unsigned long, unsigned long));
MLYDECL( devyield_result, MLDefaultYielder, ( MLINK, MLYieldParameters));

#if WINDOWS_MATHLINK
extern HWND MLInitializeIcon P(( HINSTANCE hinstCurrent, int nCmdShow));
extern HANDLE MLInstance;
extern HWND MLIconWindow;
#endif
extern int	MLAbort, MLDone;
extern long MLSpecialCharacter;

#endif /* _MLTM_H */



#endif /* _MATHLINK_H */
