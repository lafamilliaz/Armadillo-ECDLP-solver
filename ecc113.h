/*
	Armadillo ECDLP solver coded by ME!
*/

#pragma once

#include <Windows.h>
#include <emmintrin.h> 

#include "src_x64\miracl.h"
#include "matrix.h"

#define WORMS 256  //Play with this to see different speeds
#define MAX_SQUARING 128

typedef __m128i poly128;

typedef struct
{
	DWORD hamming;
	BYTE A;           //A to define the curve
	char *szOrder;    //Order of the public point P
	char *szCofactor; //Cofactor

	char *szPx;       //Public point x value
	char *szPy;       //Public point y value
	char *szQx;       //Fixed point Qx where Q = x*P
	char *szQy;       //Fixed point Qy where Q = x*P
} ECC_CURVE_DATA;


typedef struct
{
	poly128 x;
	poly128 y;
	BOOL point_at_infinity;
} EPOINT;

typedef struct
{
	DWORD64 iter;
	DWORD dist_cnt;
	DWORD64 square_cnt[MAX_SQUARING];
	DWORD64 neg_cnt; //For negation after point add like X = X + o^l*X. Resulting X is using lowest y coord

	big c, d;
} WORM;

typedef struct
{
	miracl *mip;
	big t1, t2, t3; //big temps
	DWORD A;
} ECC2M;

typedef struct
{
	BOOL status;
	big *mul_array;
} FROBENIUS;

typedef struct
{
	BOOL status;
	big mul_vgal;
} NEGATION;

//
// wdcnt is worm distinguish count. This is how many points this worm have found
//
typedef BOOL (__cdecl *FOUNDCALLBACK)(miracl *mip, DWORD totcnt, DWORD wdcnt, DWORD64 witer, big c, big d, epoint *X, void *ud);

typedef struct
{
	MATRIX *ma;  //Matrix for base convertion
	
	BYTE A;
	BYTE B;
	DWORD m, a, b, c;
	DWORD hamming;

	BOOL collision;

	big order; // Order of P
	epoint *P; // Base point
	epoint *Q; // where Q = x*P
	
	FROBENIUS frobenius;
	NEGATION negation;

	volatile LONG64 iter;
	DWORD restart_cnt;
	DWORD dist_cnt;

	CRITICAL_SECTION crit;
	miracl *mip;
	
	FOUNDCALLBACK f;
	void *ud; //Userdata for callback
} ECCTX;

typedef struct
{
	ECCTX *ecctx;
	DWORD array_id; //NOT Thread ID!
} THREAD_DATA;

void print_pol(poly128 pol);
void print_point(EPOINT *ep);

DWORD WINAPI rho_113(void *args);
ECCTX *ecc_create(const ECC_CURVE_DATA *curve, FOUNDCALLBACK f, void *ud);
void ecc_destroy(ECCTX **_ecctx);

