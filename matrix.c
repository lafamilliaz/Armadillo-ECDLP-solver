/*
	Armadillo ECDLP solver coded by ME!
*/

#include "matrix.h"

#include <math.h>

#define XOR(v, n) _mm_xor_si128(v, n) 

typedef struct
{
	big a;
	big b;
} AB;


static BOOL find_line(miracl *mip, AB *ab, DWORD size, DWORD pos, DWORD *line)
{
	DWORD i;
	for(i = 0; i < size; i++)
	{
		if(1 == getdig(mip, ab[i].a, pos))
		{
			(*line) = i;
			return TRUE;
		}
	}
	return FALSE;
}

static void swap_line(AB *ab, DWORD l1, DWORD l2)
{
	AB tmp;
	memcpy(&tmp, &ab[l1], sizeof(AB));
	memcpy(&ab[l1], &ab[l2], sizeof(AB));
	memcpy(&ab[l2], &tmp, sizeof(AB));
}

static void reverse_array(AB *ab, DWORD size)
{
	DWORD i;
	AB *abtmp;
	abtmp = (AB*)calloc(size, sizeof(AB));
	for(i = 0; i < size; i++)
	{
		memcpy(&abtmp[i], &ab[size-1-i], sizeof(AB));
	}
	memcpy(ab, abtmp, size*sizeof(AB));
	free(abtmp);
}

static void get_pol(miracl *mip, DWORD size, AB *ab, big res, DWORD pol, BOOL use_a)
{
	DWORD i;
	big t;
	t = mirvar(mip, 0);
	for(i = 0; i < size; i++)
	{
		if(pol&1)
		{
			if(use_a)
			{
				add2(ab[i].a, t, t);
			}
			else
			{
				add2(ab[i].b, t, t);
			}
		}
		pol = pol >> 1;
	}
	copy(t, res);
	mirkill(t);
}

static void prepare_matrix_x_bit(big *mm, miracl *mip, DWORD size, AB *ab, DWORD bit, BOOL use_a)
{
	DWORD i, j, kk, totbits;
	totbits = (DWORD)pow((double)2, (double)bit);
	for(i = 0, kk = 0; i < size; i+=bit)
	{
		for(j = 0; j < totbits; j++)
		{
			if(size > i + (DWORD)(log((double)j) / log(2.0)))
			{
				get_pol(mip, bit, &ab[i], mm[kk+j], j, use_a);
			}
			if(i+1 >= size && j >= 1)
			{
				return;
			}
		}
		kk+=totbits;
	}
}

static void prepare_nb_matrix_8bit(MATRIX *m, miracl *mip, DWORD size, AB *ab)
{
	prepare_matrix_x_bit(m->m8nb, mip, size, ab, 8, TRUE);
}

static void prepare_matrix_8bit(MATRIX *m, miracl *mip, DWORD size, AB *ab)
{
	prepare_matrix_x_bit(m->m8, mip, size, ab, 8, FALSE);
}


//
// Using Gaussian elimination to find the invers of a matrix
// Miracl library must be set up with a curve for this to work
// Beta migh be t or t+1
//
// Will looks something like this before solving. Result will be in B matrix
//     A        B
// |10110101|00000001|
// |11100100|00000010|
// |11110111|00000100|
// |01010100|00001000|
// |01011101|00010000|
// |10000110|00100000|
// |01101010|01000000|
// |11100011|10000000|
//
// A matrix will contain Beta^2^i
//

static BOOL matrix_set_and_inverse(MATRIX *m, miracl *mip, DWORD degree)
{
	BOOL retVal = FALSE;
	DWORD i, line, current;
	AB *ab;
	big t;
	
	if(NULL == (ab = (AB*)calloc(degree, sizeof(AB))))
	{
		return FALSE;
	}

	t = mirvar(mip, m->Beta); //Beta is usually t+1

	for(i = 0; i < degree; i++)
	{
		ab[i].a = mirvar(mip, 0);
		ab[i].b = mirvar(mip, 0);
		copy(t, ab[i].a);             //Store B, B^2, B^2^2 ...  in A matrix
		putdig(mip, 1, ab[i].b, i+1); //Store 1, t, t^2, t^3 ... in B matrix
		modsquare2(mip, t, t);
	}

	//
	// Make a copy of A matrix before we sole. This will be the matrix to go from NB to poly base
	//
	
	prepare_nb_matrix_8bit(m, mip, degree, ab);

	//
	// First take the low side if the matrix
	//
	for(current = 0; current < degree; current++)
	{
		if(!find_line(mip, &ab[current], degree-current, current+1, &line))
		{
			goto Cleanup;
		}

		line += current;

		if(line != current)
		{
			swap_line(ab, current, line);
		}

		for(i = current+1; i < degree; i++)
		{
			if(1 == getdig(mip, ab[i].a, current+1))
			{
				add2(ab[current].a,  ab[i].a, ab[i].a);
				add2(ab[current].b,  ab[i].b, ab[i].b);	
			}
		}
	}

	reverse_array(ab, degree);

	//
	// Now the upper side
	//

	for(current = 0; current < degree; current++)
	{
		if(!find_line(mip, &ab[current], degree-current, degree-current, &line))
		{
			goto Cleanup;
		}

		line += current;

		if(line != current)
		{
			goto Cleanup;
		}

		for(i = current+1; i < degree; i++)
		{
			if(1 == getdig(mip, ab[i].a, degree-current))
			{
				add2(ab[current].a,  ab[i].a, ab[i].a);
				add2(ab[current].b,  ab[i].b, ab[i].b);	
			}
		}
	}
	
	reverse_array(ab, degree);

	prepare_matrix_8bit(m, mip, degree, ab);

	retVal = TRUE;

Cleanup:
	for(i = 0; i < degree; i++)
	{
		mirkill(ab[i].a);
		mirkill(ab[i].b);
	}
	free(ab);
	mirkill(t);
	return retVal;
}


//
// Change from t^k to B^2^k where B is t or t+1
//

//
// 8 bit at a time. Loop removed for speedup
//

__m128i change_base_113_poly_8(__m128i *ma, __m128i a)
{
	DWORD off;
	DWORD64 bits;
	__m128i res128 = {0};

	//
	// First 64 bits
	//

	bits = a.m128i_u64[0];

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;
		
	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;
	
	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;
		
	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;


	//
	// Next 64 bit
	//
	
	bits = a.m128i_u64[1];
	
	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	off = bits & 255;
	res128 = XOR(ma[off], res128);
	bits = bits >> 8;
	ma += 256;

	//
	// last bit
	//

	off = bits & 1;
	res128 = XOR(ma[off], res128);

	return res128;
}



//
// CTOR
//
MATRIX * matrix_create(miracl *mip, DWORD degree, BYTE Beta, BOOL do_inverse)
{
	DWORD i;
	MATRIX *m;
	DWORD max_size8 = ((degree/8)+1)*256;
	
	big b, t;
	BOOL error = FALSE;

	if(113 != degree)
	{
		return NULL;
	}

	b = mirvar(mip, Beta);
	t = mirvar(mip, 0);

	for(i = 0; i < degree; i++)
	{
		modsquare2(mip, b, b);
		add2(t, b, t);
	}

	if(1 != size(t))
	{
		printf("!WARNING! Matrix can't be solved\n");
		error = TRUE;
	}

	mirkill(b);
	mirkill(t);

	if(error)
	{
		return NULL;
	}

	if(NULL == (m = (MATRIX*)calloc(1, sizeof(MATRIX))))
	{
		return NULL;
	}

	m->degree = degree;
	m->size8 = max_size8;
	m->Beta = Beta;
	
	if(NULL == (m->m8 = (big*)calloc(max_size8, sizeof(big))))
	{
		matrix_destroy(&m);
		return NULL;
	}
	
	if(NULL == (m->m8nb = (big*)calloc(max_size8, sizeof(big))))
	{
		matrix_destroy(&m);
		return NULL;
	}
	
	for(i = 0; i < max_size8; i++)
	{
		m->m8[i] = mirvar(mip, 0);
		m->m8nb[i] = mirvar(mip, 0);
	}

	if(do_inverse)
	{
		if(!matrix_set_and_inverse(m, mip, degree))
		{
			matrix_destroy(&m);
			return NULL;
		}

		for(i = 0; i < max_size8; i++)
		{
			m->ma128_8[i] = _mm_lddqu_si128((__m128i*)m->m8[i]->w);
			m->ma128_nb_8[i] = _mm_lddqu_si128((__m128i*)m->m8nb[i]->w);
		}
	}

	return m;
}

//
// DTOR
//
void matrix_destroy(MATRIX **_m)
{
	DWORD i;
	MATRIX *m;
	if(NULL == _m || NULL == *_m)
	{
		return;
	}
	m = (*_m);

	for(i = 0; i < m->size8; i++)
	{
		mirkill(m->m8[i]);
		mirkill(m->m8nb[i]);
	}
	free(m->m8);
	free(m->m8nb);

	free(m);
	(*_m) = NULL;
}

