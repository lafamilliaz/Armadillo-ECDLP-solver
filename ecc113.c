/*
	Coded by ME!
	Based on the code from:
	*
	* > ecdl2K-108.c <
	* Purpose: Fast arithmetic for computing discrete logs on elliptic curves.
	* Copyright: Robert J. Harley, 1997-1999.
	* Contact: Robert.Harley@inria.fr
	*
*/

#include <time.h>
#include <emmintrin.h> 

#include "ecc113.h"

#define XOR(v, n) _mm_xor_si128(v, n) 
#define OR(v, n) _mm_or_si128(v, n)

//
// Private functions
//

//
// The following function are called all the time and are time critical. Can you make them faster than me?? :)
//

//
// Multiply two polynomial. Source: https://www.intel.com/content/dam/www/public/us/en/documents/white-papers/polynomial-multiplication-instructions-paper.pdf
//
static void __inline GF2m_mul_2x2_xmm(__m128i *c1, __m128i *c0, __m128i b, __m128i a)
{
	__m128i t1, t2;
	*c0 = _mm_clmulepi64_si128(a, b, 0x00);
	*c1 = _mm_clmulepi64_si128(a, b, 0x11);
	t1 = _mm_shuffle_epi32(a, 0xEE);
	t1 = _mm_xor_si128(a, t1);
	t2 = _mm_shuffle_epi32(b, 0xEE);
	t2 = _mm_xor_si128(b, t2);
	t1 = _mm_clmulepi64_si128(t1, t2, 0x00);
	t1 = _mm_xor_si128(*c0, t1);
	t1 = _mm_xor_si128(*c1, t1);
	t2 = t1;
	t1 = _mm_slli_si128(t1, 8);
	t2 = _mm_srli_si128(t2, 8);
	*c0 = _mm_xor_si128(*c0, t1);
	*c1 = _mm_xor_si128(*c1, t2);
}

//
// Shift left
//
static __m128i __inline shl128(__m128i v, DWORD n)
{
    __m128i v1, v2;

	if (n >= 64)
	{
		v1 = _mm_slli_si128(v, 8);
		v1 = _mm_slli_epi64(v1, (n) - 64);
	}
	else
	{
		v1 = _mm_slli_epi64(v, n);
		v2 = _mm_slli_si128(v, 8);
		v2 = _mm_srli_epi64(v2, 64 - (n));
		v1 = _mm_or_si128(v1, v2);
	}
	return v1;
}

//
// Shift right
//
static __m128i __inline shr128(__m128i v, DWORD n)
{
	__m128i v1, v2;
	if (n >= 64)
	{
		v1 = _mm_srli_si128(v, 8);
		v1 = _mm_srli_epi64(v1, (n) - 64);
	}
	else
	{
		v1 = _mm_srli_epi64(v, n);
		v2 = _mm_srli_si128(v, 8);
		v2 = _mm_slli_epi64(v2, 64 - (n));
		v1 = _mm_or_si128(v1, v2);
	}
	return v1;
}

//
// Reduce modulo t^113+t^9+1
//

static __m128i __inline poly_reduce_113_9(poly128 hi00, poly128 lo00)
{
	__m128i tmp;
	lo00 = XOR(lo00, OR(shl128(hi00, 15), shr128(hi00, 113)));
	lo00 = XOR(lo00, OR(shl128(hi00, 24), shr128(hi00, 104)));
	tmp = shr128(lo00, 113);
	lo00 = XOR(lo00, tmp);
	lo00 = XOR(lo00, shl128(tmp, 113));
	lo00 = XOR(lo00, shl128(tmp, 9));
	return lo00;
}

//
// Square a polynomial
//
static void __inline GF2m_square(__m128i *c1, __m128i *c0, __m128i a0)
{
	*c0 = _mm_clmulepi64_si128(a0, a0, 0x00);
	*c1 = _mm_clmulepi64_si128(a0, a0, 0x11);
}

//
// Polynomial mul in the base t^113 + t^9 + 1
//
static poly128 __inline sse_poly_product(poly128 p1, poly128 p2)
{
	poly128 hi, lo;
	GF2m_mul_2x2_xmm(&hi, &lo, p1, p2);
	return poly_reduce_113_9(hi, lo);
}

static poly128 __inline sse_poly_square(poly128 p1)
{
	poly128 hi, lo;
	GF2m_square(&hi, &lo, p1);
	return poly_reduce_113_9(hi, lo);	
}

static poly128 __inline sse_poly_square_2(poly128 p1)
{
	poly128 hi, lo, res;
	GF2m_square(&hi, &lo, p1);
	res = poly_reduce_113_9(hi, lo);
	GF2m_square(&hi, &lo, res);
	return poly_reduce_113_9(hi, lo);
}

static poly128 __inline poly_add(poly128 p1, poly128 p2)
{
	return _mm_xor_si128(p1, p2);
}

static BOOL __inline poly_equal(poly128 p1, poly128 p2)
{
	return !(p1.m128i_u64[0] ^ p2.m128i_u64[0] | p1.m128i_u64[1] ^ p2.m128i_u64[1]);
}

//
// Find the hamming weight of a poly
//
static DWORD __inline hamming_popcnt_poly(poly128 n)
{
	DWORD res;
	res = (DWORD)_mm_popcnt_u64(n.m128i_u64[0]);
	res += (DWORD)_mm_popcnt_u64(n.m128i_u64[1]);
	return res;
}

//
// End time critical
//

//
// Following functions are not so time critical
//
static poly128 poly_inverse(poly128 a)
{
	poly128 shift, temp, res;
	DWORD m, r, rsft;
	int s;

	s = 5;
	res = a;
	m = 113 - 1;	

	while (s >= 0) 
	{
		r = m >> s;
		shift = res;
		for (rsft = 0; rsft < (r>>1); rsft++)
		{
			shift = sse_poly_square(shift);
		}

		temp = sse_poly_product(res, shift);
		
		if (r&1) 
		{
			temp = sse_poly_square(temp);
			res = sse_poly_product(temp, a);
		} 
		else
		{
			res = temp;
		}
		
		--s;
	};
	return (sse_poly_square(res));
}


static __inline poly128 poly_quotient(poly128 p1, poly128 p2)
{
	return sse_poly_product(p1, poly_inverse(p2));
}

static void poly_multi_inverse(poly128 *p, poly128 *res)
{
	int i;
	poly128 c[WORMS], u, t;

	if(1 == WORMS)
	{
		(*res) = poly_inverse(*p);
		return;
	}

	c[0] = p[0];
	
	for(i = 1; i < (int)WORMS; i++)
	{
		c[i] = sse_poly_product(c[i-1], p[i]);
	}
	
	u = poly_inverse(c[WORMS-1]);
	
	for(i = (int)(WORMS-1); i >= 1; i--)
	{
		res[i] = sse_poly_product(u, c[i-1]);
		t = u;
		u = sse_poly_product(p[i], t);
	}

	res[0] = u;
}

//
// This works for 64 bit's miracl only but should be a quick fix to get it to work with x86
//
static void EP_set(EPOINT *E, epoint *X)
{
	E->x.m128i_u64[0] = X->X->w[0];
	E->x.m128i_u64[1] = X->X->w[1];
	
	E->y.m128i_u64[0] = X->Y->w[0];
	E->y.m128i_u64[1] = X->Y->w[1];
}

static void EP_get(EPOINT *E, epoint *X)
{
	X->X->w[0] = E->x.m128i_u64[0];
	X->X->w[1] = E->x.m128i_u64[1];
	X->Y->w[0] = E->y.m128i_u64[0];
	X->Y->w[1] = E->y.m128i_u64[1];
	X->X->len = 2;
	X->Y->len = 2;
}

static poly128 __inline poly_squareNTimes_i(poly128 p1, DWORD n)
{
	int i = (int)n;
	for(; i >= 2; i -= 2)
	{
		p1 = sse_poly_square_2(p1);
	}
	if(i)
	{
		p1 = sse_poly_square(p1);
	}
	return p1;
}

static poly128 __inline poly_squareNTimes_nb_i(poly128 p1, DWORD n)
{
	__m128i v1, v2, mask;
	mask = _mm_xor_si128(mask, mask);
	mask.m128i_u64[1] = 0xFFFE000000000000;
	v1 = _mm_slli_epi64(p1, n);
	v2 = _mm_slli_si128(p1, 8);
	v2 = _mm_srli_epi64(v2, 64 - (n));
	v1 = _mm_or_si128(v1, v2);
	v2 = _mm_and_si128(v1, mask);
	v2 = _mm_srli_si128(v2, 8);
	v2 = _mm_srli_epi64(v2, (113) - 64);
	v1 = _mm_or_si128(v1, v2);
	v1 = _mm_andnot_si128(mask, v1);
	return v1;
}


static BOOL ellipticMultiSum(DWORD A, EPOINT *p1, EPOINT *p2, poly128 *denT, EPOINT *res)
{
	DWORD i;
	poly128 lam, t;
	poly128 inv_x1x2[WORMS];
	poly128 a;

	a.m128i_u64[0] = A;
	a.m128i_u64[1] = 0;
		
	poly_multi_inverse(denT, inv_x1x2);
	
	for(i = 0; i < WORMS; i++)
	{
		lam = sse_poly_product(poly_add(p1[i].y, p2[i].y), inv_x1x2[i]);
		t = poly_add(poly_add(poly_add(sse_poly_square(lam), lam), p2[i].x), a);
		res[i].x = poly_add(t, p1[i].x);
		res[i].y = poly_add(poly_add(sse_poly_product(lam, t), res[i].x), p1[i].y);
	}

	return TRUE;
}

static void __inline poly_point_copy(EPOINT *P, EPOINT *R)
{
	R->point_at_infinity = P->point_at_infinity;
	R->x = P->x;
	R->y = P->y;
}

static BOOL __inline poly_point_neg_i(EPOINT *E)
{
	poly128 neg;
	neg = _mm_xor_si128(E->x, E->y);

	if(neg.m128i_u64[1] < E->y.m128i_u64[1])
	{
		E->y = neg;
		return TRUE;
	}
	if(neg.m128i_u64[1] > E->y.m128i_u64[1])
	{
		return FALSE;
	}
	if(neg.m128i_u64[0] < E->y.m128i_u64[0])
	{
		E->y = neg;
		return TRUE;
	}
	return FALSE;
}

static ECC2M *ecc113_create(miracl *mip, DWORD A)
{
	ECC2M *ecc;
	if(NULL == (ecc = (ECC2M*)calloc(1, sizeof(ECC2M))))
	{
		return NULL;
	}

	ecc->mip = mip;
	ecc->A = A;
	ecc->t1 = mirvar(mip, 0);
	ecc->t2 = mirvar(mip, 0);
	ecc->t3 = mirvar(mip, 0);

	return ecc;
}

static BOOL ecc113_destroy(ECC2M **_ecc)
{
	ECC2M *ecc;
	if(NULL == _ecc || NULL == *_ecc)
	{
		return FALSE;
	}

	ecc = (*_ecc);
	mirkill(ecc->t1);
	mirkill(ecc->t2);
	mirkill(ecc->t3);
	free(ecc);
	(*_ecc) = NULL;
	return TRUE;
}

static void reset_square_cnt(DWORD64 *a, DWORD size)
{
	DWORD i;
	for(i = 0; i < size; i++)
	{
		a[i] = 0;
	}
}

static BOOL get_neg_point(epoint *X, big neg)
{
	add2(X->X, X->Y, neg);
	if(mr_compare(neg, X->Y) < 0)
	{
		copy(neg, X->Y);
		return TRUE;
	}
	return FALSE;
}

static void get_worm_startpos(ECCTX *ecctx, miracl *mip, WORM *wp, epoint *X)
{
	big neg;
	EnterCriticalSection(&ecctx->crit);

	neg = mirvar(mip, 0);
	reset_square_cnt(wp->square_cnt, MAX_SQUARING);
	bigrand(mip, ecctx->order, wp->c);
	bigrand(mip, ecctx->order, wp->d);
	ecurve2_mult2(mip, wp->c, ecctx->Q, wp->d, ecctx->P, X); //X = c*Q + d*P
	if(ecctx->negation.status)
	{
		if(get_neg_point(X, neg))
		{
			mad(mip, ecctx->negation.mul_vgal, wp->c, ecctx->negation.mul_vgal, ecctx->order, ecctx->order, wp->c);
			mad(mip, ecctx->negation.mul_vgal, wp->d, ecctx->negation.mul_vgal, ecctx->order, ecctx->order, wp->d);
		}
	}
	wp->dist_cnt = 0;
	wp->iter = 0;
	wp->neg_cnt = 0;
	mirkill(neg);

	LeaveCriticalSection(&ecctx->crit);
}

static void update_cd(ECCTX *ecctx, WORM *wp)
{
	DWORD i;
	big f, exp;
	miracl *mip = mirsys(100, 0);

	f = mirvar(mip, 0);
	exp = mirvar(mip, 1);
	
	if(ecctx->negation.status && wp->neg_cnt)
	{
		copy(ecctx->negation.mul_vgal, f);
		exp->w[0] = wp->neg_cnt;
		powmod(mip, f, exp, ecctx->order, f);
		mad(mip, f, wp->c, f, ecctx->order, ecctx->order, wp->c);
		mad(mip, f, wp->d, f, ecctx->order, ecctx->order, wp->d);
		wp->neg_cnt = 0;
	}

	if(ecctx->frobenius.status)
	{
		for(i = 0; i < MAX_SQUARING; i++)
		{
			if(wp->square_cnt[i])
			{
				copy(ecctx->frobenius.mul_array[i-1], f);
				incr(mip, f, 1, f);
				exp->w[0] = wp->square_cnt[i];
				powmod(mip, f, exp, ecctx->order, f);
				mad(mip, f, wp->c, f, ecctx->order, ecctx->order, wp->c);
				mad(mip, f, wp->d, f, ecctx->order, ecctx->order, wp->d);
			}

			wp->square_cnt[i] = 0;
		}
	}

	mirkill(f);
	mirkill(exp);
	mirexit(mip);
}

static BOOL verify_point(ECCTX *ecctx, miracl *mip, big c, big d, epoint *X)
{
	BOOL retVal = FALSE;
	epoint *T;
	T = epoint_init(mip);
	ecurve2_mult2(mip, c, ecctx->Q, d, ecctx->P, T);
	if(epoint2_comp(mip, T, X))
	{
		retVal = TRUE;
	}
	epoint_free(T);
	return retVal;
}

static BOOL found_point(ECCTX *ecctx, miracl *mip, WORM *wp, epoint *X, BOOL *restart, DWORD worm_id)
{
	BOOL retVal;
	EnterCriticalSection(&ecctx->crit);
	
	retVal = FALSE;
	(*restart) = FALSE;

	if(point_at_infinity(X))
	{
		printf("\n* Hit point at infinity. Restarting Worm.. (TID:%d WORM:%d) *\n", GetCurrentThreadId(), worm_id+1);
		get_worm_startpos(ecctx, mip, wp, X);
		ecctx->restart_cnt++;
		(*restart) = TRUE;
		goto Cleanup;
	}

	update_cd(ecctx, wp);
				
	if(!verify_point(ecctx, mip, wp->c, wp->d, X))
	{
		ecctx->restart_cnt++;
		printf("\n*** Point calc ERROR. Restarting Worm.. (TID:%d WORM:%d ErrorCnt:%d) ***\n", GetCurrentThreadId(), worm_id+1, ecctx->restart_cnt);
		get_worm_startpos(ecctx, mip, wp, X);
		(*restart) = TRUE;
		goto Cleanup;
	}

	wp->dist_cnt++;
	ecctx->dist_cnt++;

	//
	// We have a distinguished point and we can collect it. Call user callback
	//

	if(ecctx->f)
	{
		if(retVal = ecctx->f(mip, ecctx->dist_cnt, wp->dist_cnt, wp->iter, wp->c, wp->d, X, ecctx->ud))
		{
			ecctx->collision = TRUE;
		}
	}

	wp->iter = 0;

Cleanup:
	LeaveCriticalSection(&ecctx->crit);

	return retVal;
}


//
// This is the rho ecdlp loop itself
//
DWORD WINAPI rho_113(void *args)
{
	BOOL restart_worm = FALSE;
	DWORD i, ri, ham;
	big A, B;
	miracl *mip;

	THREAD_DATA *tdata = (THREAD_DATA*)args;
	ECCTX *ecctx = tdata->ecctx;

	epoint *X;

	WORM *wp;
	ECC2M *ecc;
	WORM w[WORMS] = {0};
	poly128 xnb, denT[WORMS];
	EPOINT EX[WORMS], ER[WORMS];

	if(NULL == (mip = mirsys(256, 2)))
	{
		return 0;
	}

	mip->IOBASE = 16;
	
	if(NULL == (ecc = ecc113_create(mip, ecctx->A)))
	{
		return 0;
	}

	A = mirvar(mip, ecctx->A);
	B = mirvar(mip, ecctx->B);

	irand(mip, GetTickCount()*GetCurrentThreadId());
	ecurve2_init(mip, ecctx->m, ecctx->a, ecctx->b, ecctx->c, A, B, TRUE, MR_AFFINE);
	
	X = epoint_init(mip);

	for(i = 0; i < WORMS; i++)
	{
		w[i].c = mirvar(mip, 0);
		w[i].d = mirvar(mip, 0);		

		//
		// Get random starting points or read them from save state
		//

		//
		// TODO : Read points from save state file so we can continue
		//

		get_worm_startpos(ecctx, mip, &w[i], X);
		EP_set(&EX[i], X);
		EX[i].point_at_infinity = FALSE;
	}


	for(;!ecctx->collision;)
	{
		for(i = 0; (i < WORMS) && !ecctx->collision; i++)
		{
			wp = &w[i];
			
			xnb = change_base_113_poly_8(ecctx->ma->ma128_8, EX[i].x);
			ham = hamming_popcnt_poly(xnb);

			if(ecctx->hamming == ham)
			{	
				EP_get(&EX[i], X);
				if(found_point(ecctx, mip, wp, X, &restart_worm, i))
				{
					goto Cleanup;
				}
				EP_set(&EX[i], X);
				if(restart_worm)
				{
					i -= 1;
					continue;
				}
			}

			//
			// Get a random walk based on the hamming weight. RND walk need to be the same on all clients
			//

			ri = (ham % 7) + 1;
						
			wp->square_cnt[ri]++;

			if(ri > 2)
			{
				//
				// Do squaring in normal base
				//

				ER[i].x = change_base_113_poly_8(ecctx->ma->ma128_nb_8, poly_squareNTimes_nb_i(xnb, ri));
				denT[i] = _mm_xor_si128(EX[i].x, ER[i].x);
				ER[i].y = change_base_113_poly_8(ecctx->ma->ma128_nb_8, poly_squareNTimes_nb_i(change_base_113_poly_8(ecctx->ma->ma128_8, EX[i].y), ri));

			}
			else
			{	
				//
				// Do squaring in polynomial base with reduction
				//

				ER[i].x = poly_squareNTimes_i(EX[i].x, ri);
				denT[i] = _mm_xor_si128(EX[i].x, ER[i].x);
				ER[i].y = poly_squareNTimes_i(EX[i].y, ri);
			}

			wp->iter++;
		}


		//
		// Using Montgomery tric to reduce number of field inversion
		//
		
		ellipticMultiSum(ecc->A, EX, ER, denT, EX);

		for(i = 0; i < WORMS; i++)
		{
			wp = &w[i];

			//
			// Check negation. We only use the point with the lowest Y value inspected as an integer
			//

			if(poly_point_neg_i(&EX[i]))
			{
				wp->neg_cnt++;
			}				
		}

		InterlockedAdd64(&ecctx->iter, WORMS);
	}

Cleanup:
	ecc113_destroy(&ecc);

	mirkill(A);
	mirkill(B);
	epoint_free(X);
	for(i = 0; i < WORMS; i++)
	{
		mirkill(w[i].c);
		mirkill(w[i].d);
	}
	mirexit(mip);
	return 0;
}

//
// Create frobenius multipliers
//
static BOOL make_frobenius(ECCTX *ecctx)
{
	BOOL retVal = FALSE;
	DWORD i;
	big t, m, x, y;
	epoint *P1, *P2;
	miracl *mip = ecctx->mip;
	miracl *mip2 = mirsys(100, 0);

	if(!isprime(mip, ecctx->order))
	{
		return FALSE;
	}

	t = mirvar(mip, 0);
	m = mirvar(mip, 0);
	x = mirvar(mip, 0);
	y = mirvar(mip, 0);
	P1 = epoint_init(mip);
	P2 = epoint_init(mip);
	

	epoint2_copy(ecctx->P, P1);
	modsquare2(mip, P1->X, P1->X);
	modsquare2(mip, P1->Y, P1->Y);

	decr(mip, ecctx->order, 1, t);
	convert(mip, ecctx->m, m);
	divide(mip, t, m, t);

	do
	{
		bigrand(mip, ecctx->order, m);
		powmod(mip2, m, t, ecctx->order, m);
	} while(1 == size(m));

	for(i = 1; i < ecctx->m; i++)
	{
		power(mip, m, i, ecctx->order, x);
		ecurve2_mult(mip, x, ecctx->P, P2);
		if(epoint2_comp(mip, P1, P2))
		{
			retVal = TRUE;
			break;
		}
	}

	if(retVal)
	{
		if(NULL == (ecctx->frobenius.mul_array = (big*)calloc(ecctx->m, sizeof(big))))
		{
			retVal = FALSE;
			goto Cleanup;
		}
		for(i = 0; i < ecctx->m; i++)
		{
			power(mip, x, i+1, ecctx->order, m);
			ecctx->frobenius.mul_array[i] = mirvar(mip, 0);
			copy(m, ecctx->frobenius.mul_array[i]);
		}
	}

	ecctx->frobenius.status = TRUE;

Cleanup:
	mirkill(t);
	mirkill(m);
	mirkill(x);
	mirkill(y);
	epoint_free(P1);
	epoint_free(P2);
	mirexit(mip2);
	return retVal;
}

static BOOL make_negation(ECCTX *ecctx)
{
	BOOL retVal = FALSE;
	big t, m;
	miracl *mip = ecctx->mip;
	miracl *mip2 = mirsys(100, 0);

	if(!isprime(mip, ecctx->order))
	{
		return FALSE;
	}

	t = mirvar(mip, 0);
	m = mirvar(mip, 0);

	decr(mip, ecctx->order, 1, t);
	convert(mip, 2, m);
	divide(mip, t, m, t);

	do
	{
		bigrand(mip, ecctx->order, m);
		powmod(mip2, m, t, ecctx->order, m);
	} while(1 == size(m));

	ecctx->negation.mul_vgal = mirvar(mip, 0);
	copy(m, ecctx->negation.mul_vgal);
	ecctx->negation.status = TRUE;
	retVal = TRUE;

	mirkill(t);
	mirkill(m);
	mirexit(mip2);
	return retVal;
}


//
// Public functions
//


void print_pol(poly128 pol)
{
	printf("%016I64X%016I64X\n", pol.m128i_u64[1], pol.m128i_u64[0]);
}

void print_point(EPOINT *ep)
{
	print_pol(ep->x);
	print_pol(ep->y);
}

void ecc_destroy(ECCTX **_ecctx)
{
	ECCTX *ecctx = (*_ecctx);

	matrix_destroy(&ecctx->ma);
	DeleteCriticalSection(&ecctx->crit);
	if(ecctx->mip)
	{

		epoint_free(ecctx->P);
		epoint_free(ecctx->Q);
		mirkill(ecctx->order);
		mirexit(ecctx->mip);
	}
	free(ecctx);
	(*_ecctx) = NULL;
}

ECCTX *ecc_create(const ECC_CURVE_DATA *curve, FOUNDCALLBACK f, void *ud)
{
	ECCTX *ecctx;
	ECCTX *retVal = NULL;

	miracl *mip;
	big x, y;
	epoint *T;
	big A, B;

	if(NULL == (ecctx = (ECCTX*)calloc(1, sizeof(ECCTX))))
	{
		return NULL;
	}

	InitializeCriticalSection(&ecctx->crit);
	ecctx->f = f;
	ecctx->ud = ud;

	if(NULL == (mip = mirsys(256, 2)))
	{
		ecc_destroy(&ecctx);
		return NULL;
	}

	mip->IOBASE = 16;	
	irand(mip, GetTickCount());
	
	
	ecctx->hamming = curve->hamming;
	ecctx->A = curve->A;
	ecctx->B = 1;
	ecctx->m = 113;
	ecctx->a = 9;
	ecctx->b = 0;
	ecctx->c = 0;

	ecctx->mip = mip;

	ecctx->order = mirvar(mip, 0);

	A = mirvar(mip, curve->A);
	B = mirvar(mip, 1);

	if(!ecurve2_init(mip, 113, 9, 0, 0, A, B, TRUE, MR_AFFINE))
	{
		mirkill(A);
		mirkill(B);
		ecc_destroy(&ecctx);
		return NULL;
	}
	mirkill(A);
	mirkill(B);

	ecctx->P = epoint_init(mip);
	ecctx->Q = epoint_init(mip);
	
	if(NULL == (ecctx->ma = matrix_create(mip, 113, 3, TRUE)))
	{
		ecc_destroy(&ecctx);
		return NULL;
	}

	x = mirvar(mip, 0);
	y = mirvar(mip, 0);
	T = epoint_init(mip);

	cinstr(mip, ecctx->order, curve->szOrder);	
	cinstr(mip, x, curve->szPx);
	cinstr(mip, y, curve->szPy);

	if(1 == (DWORD)strlen(curve->szPy) && '0' == curve->szPy[0])
	{
		if(!epoint2_set(mip, x, x, 0, ecctx->P))
		{
			if(!epoint2_set(mip, x, x, 1, ecctx->P))
			{
				ecc_destroy(&ecctx);
				goto Cleanup;
			}
		}
	}
	else
	{
		if(!epoint2_set(mip, x, y, 0, ecctx->P))
		{
			ecc_destroy(&ecctx);
			goto Cleanup;
		}
	}


	ecurve2_mult(mip, ecctx->order, ecctx->P, T);

	if(!point_at_infinity(T))
	{
		ecc_destroy(&ecctx);
		goto Cleanup;
	}


	make_frobenius(ecctx);
	make_negation(ecctx);

	retVal = ecctx;

Cleanup:
	epoint_free(T);
	mirkill(x);
	mirkill(y);
	return retVal;
}