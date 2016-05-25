/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	main program for Monte Carlo simulation of photon
 *	distribution in multi-layered turbid media.
 *
 ****/

/****
 *	THINKCPROFILER is defined to generate profiler calls in
 *	Think C. If 1, remember to turn on "Generate profiler
 *	calls" in the options menu.
 ****/
#define THINKCPROFILER 0

/* GNU cc does not support difftime() and CLOCKS_PER_SEC.*/
#define GNUCC 0
#define PROFILE 1

#if THINKCPROFILER
#include <profile.h>
#include <console.h>
#endif

#include <omp.h>
#include "mcml.h"
#include <sys/time.h>
#include <time.h>
#include <mkl.h>

/*	Declare before they are used in main(). */
FILE *GetFile(char *);
short ReadNumRuns(FILE* );
void ReadParm(FILE* , InputStruct * );
void CheckParm(FILE* , InputStruct * );
void InitOutputData(InputStruct, OutStruct *);
void FreeData(InputStruct, OutStruct *);
double Rspecular(LayerStruct * );
void Launch8Photon(double, LayerStruct *, Photon8Struct *);
void HopDropSpin(InputStruct  *,Photon8Struct *,OutStruct *, VSLStreamStatePtr *);
void SumScaleResult(InputStruct, OutStruct *);
void WriteResult(InputStruct, OutStruct, char *);
void CollectResult(InputStruct , OutStruct *, OutStruct *);
void Spin8AndRoul(InputStruct  * In_Ptr, Photon8Struct * Photon_Ptr, VSLStreamStatePtr  stream, double * result, short * count, double *);
static double dtime()
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday(&mytime,(struct timezone *) 0);
    tseconds = (double) (mytime.tv_sec + (double)mytime.tv_usec * 1.0e-6);
    return (tseconds) ;
}

/***********************************************************
 *	If F = 0, reset the clock and return 0.
 *
 *	If F = 1, pass the user time to Msg and print Msg on
 *	screen, return the real time since F=0.
 *
 *	If F = 2, same as F=1 except no printing.
 *
 *	Note that clock() and time() return user time and real
 *	time respectively.
 *	User time is whatever the system allocates to the
 *	running of the program;
 *	real time is wall-clock time.  In a time-shared system,
 *	they need not be the same.
 *
 *	clock() only hold 16 bit integer, which is about 32768
 *	clock ticks.
 ****/
time_t PunchTime(char F, char *Msg)
{
#if GNUCC
    return(0);
#else
    static clock_t ut0;	/* user time reference. */
    static time_t  rt0;	/* real time reference. */
    double secs;
    char s[STRLEN];

    if(F==0) {
        ut0 = clock();
        rt0 = time(NULL);
        return(0);
    } else if(F==1)  {
        secs = (clock() - ut0)/(double)CLOCKS_PER_SEC;
        if (secs<0) secs=0;	/* clock() can overflow. */
        sprintf(s, "User time: %8.0lf sec = %8.2lf hr.  %s\n",
                secs, secs/3600.0, Msg);
        puts(s);
        strcpy(Msg, s);
        return(difftime(time(NULL), rt0));
    } else if(F==2) return(difftime(time(NULL), rt0));
    else return(0);
#endif
}

/***********************************************************
 *	Print the current time and the estimated finishing time.
 *
 *	P1 is the number of computed photon packets.
 *	Pt is the total number of photon packets.
 ****/
void PredictDoneTime(long P1, long Pt)
{
    time_t now, done_time;
    struct tm *date;
    char s[80];

    now = time(NULL);
    date = localtime(&now);
    strftime(s, 80, "%H:%M %x", date);
    printf("Now %s, ", s);

    done_time = now +
                (time_t) (PunchTime(2,"")/(double)P1*(Pt-P1));
    date = localtime(&done_time);
    strftime(s, 80, "%H:%M %x", date);
    printf("End %s\n", s);
}

/***********************************************************
 *	Report time and write results.
 ****/
void ReportResult(InputStruct In_Parm, OutStruct Out_Parm)
{
    char time_report[STRLEN];

    strcpy(time_report, " Simulation time of this run.");
    PunchTime(1, time_report);

    SumScaleResult(In_Parm, &Out_Parm);
    WriteResult(In_Parm, Out_Parm, time_report);
}

/***********************************************************
 *	Get the file name of the input data file from the
 *	argument to the command line.
 ****/
void GetFnameFromArgv(int argc,
                      char * argv[],
                      char * input_filename)
{
    if(argc>=2) {			/* filename in command line */
        strcpy(input_filename, argv[1]);
    } else
        input_filename[0] = '\0';
}


/***********************************************************
 *	Execute Monte Carlo simulation for one independent run.
 ****/
void DoOneRun(short NumRuns, InputStruct *In_Ptr, int num_threads)
{
    double t_0, t_1, t_2, t_3;
    // time 0
    t_0 = dtime();

    int i;
    register long i_photon;
    /* index to photon. register for speed.*/
    // Declare out_parm[num_threads]. (out_parm[0] for suming up all results)
    OutStruct * out_parm;		/* distribution of photons.*/
    out_parm = (OutStruct *)malloc(sizeof(OutStruct) * num_threads);

    Photon8Struct photon;
    long num_photons = In_Ptr->num_photons, photon_rep=10;

#if THINKCPROFILER
    InitProfile(200,200);
    cecho2file("prof.rpt",0, stdout);
#endif
    //Initial Output Data
    for( i = 0 ; i < num_threads ; i++ ) {
        InitOutputData(*In_Ptr, &out_parm[i]);
        out_parm[i].Rsp = Rspecular(In_Ptr->layerspecs);
    }

    //Initial Random Seed
    unsigned int *rand_seed;
    rand_seed = (unsigned int *)malloc(sizeof(unsigned int) * num_threads);
    for (i = 0 ; i < num_threads ; i++ ) {
        rand_seed[i] = (unsigned int) (time(NULL) ^ i);
    }
    // time 1
    t_1 = dtime();

    PunchTime(0, "");
    #pragma omp parallel private(photon)
    {
        int tid = omp_get_thread_num(), j;
	    VSLStreamStatePtr stream;
        Launch8Photon(out_parm[tid].Rsp, In_Ptr->layerspecs, &photon);
    	double result[1024];
    	short count = 0;
        long num_photons_thread = num_photons/omp_get_num_threads();
        double* result_t;
        result_t = (double *)malloc(sizeof(double)*16);
//    	vslNewStream( &stream, VSL_BRNG_MT19937, rand_seed[tid] );
        vslNewStream( &stream, VSL_BRNG_MT19937, 777 );
        while(num_photons_thread > 0)
        {
            HopDropSpin(In_Ptr, &photon, &out_parm[tid], stream); 
            // May launch new photons
            for(j = 0; j < 8; j++)
            {
                int vec_num = j;
                if(photon.dead[vec_num] == 1)
                {
                    LaunchPhoton(out_parm[tid].Rsp, In_Ptr->layerspecs, &photon, vec_num);
                    num_photons_thread--;
                    //printf("%ddead !!  \n",num_photons/omp_get_num_threads() - num_photons_thread);
                }    
            }
        }
	    vslDeleteStream( &stream );
    }

    // time 2
    t_2 = dtime();

#if THINKCPROFILER
    exit(0);
#endif

    for( i = 1 ; i < num_threads ; i++ )
        CollectResult(*In_Ptr, &out_parm[0], &out_parm[i]);
    // time 3
    t_3 = dtime();
    printf("Initial time    : %3.3f (s)\n", t_1 - t_0 );
    printf("Kernel  time    : %3.3f (s)\n", t_2 - t_1 );
    printf("Collection time : %3.3f (s)\n", t_3 - t_2 );

    ReportResult(*In_Ptr, out_parm[0]);
//    for( i = 0 ; i < num_threads ; i++ )
//        FreeData(*In_Ptr, &out_parm[i]);
}

/***********************************************************
 *	The argument to the command line is filename, if any.
 *	Macintosh does not support command line.
 ****/
char main(int argc, char *argv[])
{
    char input_filename[STRLEN];
    FILE *input_file_ptr;
    short num_runs;	/* number of independent runs. */
    InputStruct in_parm;
    int num_threads;

    ShowVersion("Version 1.2, 1993");
    GetFnameFromArgv(argc, argv, input_filename);
    input_file_ptr = GetFile(input_filename);
    CheckParm(input_file_ptr, &in_parm);
    num_runs = ReadNumRuns(input_file_ptr);

    printf("Input the number of threads : ");
    scanf("%d", &num_threads);
    omp_set_num_threads(num_threads);

    #pragma omp parallel
    #pragma omp master
    {
        printf("%d threads start... \n", omp_get_num_threads());
    }

    while(num_runs--)  {
        ReadParm(input_file_ptr, &in_parm);
        DoOneRun(num_runs, &in_parm, num_threads);
    }

    fclose(input_file_ptr);
    return(0);
}
