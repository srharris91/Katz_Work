#ifndef __FORTTYPE


#ifdef CRAY

#define __REAL	REAL
#define __float	float
#define __NCHPWD	8

#define __INTEGER	INTEGER
#define __int	int
#define __NCHINT	8

#define __POINTER	INTEGER
#define __NCHPTR	8

#define __FORTTYPE
#endif


#if ( ( sgi && mips ) || DEC_ALPHA || PTR64 )

#ifdef D_PRECISION
#define __REAL	REAL*8
#define __float	double
#define __NCHPWD	8

#else
#define __REAL	REAL
#define __float	float
#define __NCHPWD	4
#endif

#define __INTEGER	INTEGER
#define __int	int
#define __NCHINT	4

#if ( MIPS4 || _MIPS_SZPTR==64 || DEC_ALPHA || PTR64 )
#define __POINTER	INTEGER*8
#define __NCHPTR	8

#else
#define __POINTER	INTEGER
#define __NCHPTR	4
#endif

#define __FORTTYPE
#endif


#if ( CONVEX_PD8 )

#define __REAL	REAL
#define __float	double
#define __NCHPWD	8

#define __INTEGER	INTEGER
#define __int	long long int
#define __NCHINT	8

#define __POINTER	INTEGER*4
#define __NCHPTR	8

#define __FORTTYPE
#endif


#ifndef __FORTTYPE

#ifdef D_PRECISION
#define __REAL	REAL*8
#define __float	double
#define __NCHPWD	8

#else
#define __REAL	REAL
#define __float	float
#define __NCHPWD	4
#endif

#define __INTEGER	INTEGER
#define __int	int
#define __NCHINT	4

#define __POINTER	INTEGER
#define __NCHPTR	4

#define __FORTTYPE
#endif

#endif
