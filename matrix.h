/*
	Armadillo ECDLP solver coded by ME!
*/

#pragma once

#include <Windows.h>
#include <emmintrin.h> 

#include "src_x64\miracl.h"

typedef struct
{
	DWORD degree; //Degree of poly
	DWORD size8;  //Size of matrix 8bit
	BYTE Beta;    //Usually t or t+1
	big *m8;      //The matrix itself for 8bit
	big *m8nb;    //NOT inversed matrix to go from Normal base to polynomial base

	__m128i ma128_8[((113/8)+1)*256];
	__m128i ma128_nb_8[((113/8)+1)*256];
} MATRIX;

MATRIX * matrix_create(miracl *mip, DWORD size, BYTE Beta, BOOL do_inverse);
void matrix_destroy(MATRIX **_m);
__m128i change_base_113_poly_8(__m128i *ma, __m128i a);
