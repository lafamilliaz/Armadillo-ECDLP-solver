/*
	Armadillo ECDLP solver coded by ME!
*/

#include <stdio.h>
#include <windows.h>
#include <conio.h>

#include "ecc113.h"

#ifdef _WIN64
	#include "src_x64/miracl.h"
#else
	This does NOT work :)
#endif


#define THREADCNT 1 //How many threads to use
#define MAXTHREADS 128

typedef enum
{
	ECC_TESTCURVE = 0,
	ECC_ARMACURVE
};

ECC_CURVE_DATA ecc_curve_array[] =
{
	//Test curev with order 2^74
	{28, 0, "31F1F2A998BD31AF391", "A40522F35C", "1F60D9876DB4712343C8CA7AE8256", "0FA15CFCE7C0BE8776B607388D618", "1B4CF7B27FBED64D518BCF7262AD2", "6BCC7C019026BD2A2636BDD558F2"},

	//Real curve with 112 bits order with order 2^112. This was the armadillo parameters we just solved :)
	{25, 1, "FFFFFFFFFFFFFFFDBF91AF6DEA73", "2", "112576594A7C5775C3EACF722DB3D", "1F95740EF56C833AAE2A77B92CD69", "1B7FAB8B73F679CB3431ED56F0225", "1B7F89783E73BEC1B86172A07C03E"},
};


//
// return TRUE to stop it
//
BOOL __cdecl found_point(miracl *mip, DWORD totcnt, DWORD wcnt, DWORD64 witer, big c, big d, epoint *X, void *ud)
{
	//
	// Collect points here
	//
	printf("Worm iter: %08I64X - Total cnt:%d Worm cnt:%d          \n", witer, totcnt, wcnt);
	printf("c: ");
	cotnum(mip, c, stdout);
	printf("d: ");
	cotnum(mip, d, stdout);
	printf("x: ");
	cotnum(mip, X->X, stdout);
	printf("y: ");
	cotnum(mip, X->Y, stdout);	
	printf("\n");
	//if(10 < totcnt)
	//{
	//	return TRUE;
	//}
	return FALSE;
}

static BOOL start_threads(THREAD_DATA *tdata, ECCTX *ecctx, HANDLE *hThread, DWORD thread_cnt)
{
	DWORD i, tid;
	for(i = 0; i < thread_cnt; i++)
	{
		tdata[i].ecctx = ecctx;
		tdata[i].array_id = i;
		if(INVALID_HANDLE_VALUE == (hThread[i] = CreateThread(NULL, 0, rho_113, &tdata[i], 0, &tid)))
		{
			printf("\nCreateThread error\n");
			return FALSE;
		}
		printf("Thread ID %d Started\n", tid);
		//
		// We will be nice to you OS :)
		//
		SetThreadPriority(hThread[i], THREAD_PRIORITY_IDLE);
		Sleep(200);
	}
	printf("\n");
	return TRUE;
}

static void wait_for_thread_exit(HANDLE *hThread, DWORD thread_cnt)
{
	DWORD i;
	printf("\nWaiting for all threads to stop\n");
	WaitForMultipleObjects(THREADCNT, hThread, TRUE, INFINITE);
	for(i = 0; i < thread_cnt; i++)
	{
		CloseHandle(hThread[i]);
		hThread[i] = NULL;
	}
	printf("All threads stopped\n");
}

int main (int argc, char *argv[])
{
	ECCTX *ecctx = NULL;
	DWORD t_start;
	float speed, t_used;
	LARGE_INTEGER li;
	HANDLE hThread[MAXTHREADS];
	THREAD_DATA tdata[MAXTHREADS] = {0};

	printf("ECDLP solver by ME!\n");

	if(NULL == (ecctx = ecc_create(&ecc_curve_array[ECC_TESTCURVE], found_point, NULL)))
	{
		printf("ecc create failed\n");
		return 1;
	}

	printf("Base point order 2^%d\n", logb2(ecctx->mip, ecctx->order));

	if(!start_threads(tdata, ecctx, hThread, THREADCNT))
	{
		return 1;
	}

	//
	// Reset counter so we get a more correct speed reading
	//

	InterlockedExchange64(&ecctx->iter, 0);
	t_start = GetTickCount();

	for(;;)
	{
		if(ecctx->collision)
		{
			wait_for_thread_exit(hThread, THREADCNT);
			break;
		}

		t_used = (float)((float)(GetTickCount()-t_start) / 1000.0);

		if(t_used)
		{
			speed = (float)((double)ecctx->iter / t_used);
			li.QuadPart = ecctx->iter;
			printf("ECDLP i:0x%X:%08X Dist:%d speed:%.02f             \r", li.HighPart, li.LowPart, ecctx->dist_cnt, speed);
		}

		Sleep(1000);
	}

	ecc_destroy(&ecctx);
	printf("DONE\n");
	return 0;
}

