/* --------------------------------------------------------------------------
 * File: master.c
 * Version 1.0.6
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
 * Copyright IBM Corporation, 2012, 2017. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 * --------------------------------------------------------------------------
 */
 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* ********************************************************************** *
 *                                                                        *
 *    S t u f f   t h a t   i s   c o m m o n   t o   b o t h             *
 *                                                                        *
 * ********************************************************************** */

enum {
	USERACTION_ADDCALLBACK,    /**< Install the info callback in the worker. */
	USERACTION_REMOVECALLBACK, /**< Remove the info callback from the worker. */
	USERACTION_CHANGEOBJDIFF   /**< Change the minimal difference for objective
								*   changes. Consecutive objective values must
								*   differ by at least this amount to be
								*   considered different. */
};

enum {
	INFO_NEWDUAL,   /**< Reports a new dual bound. */
	INFO_NEWPRIMAL, /**< Reports a new primal bound. */
	INFO_DETTIME    /**< Reports the current deterministic time stamp. */
};

/* ********************************************************************** *
 *                                                                        *
 *    S e r v e r   s i d e   i m p l e m e n t a t i o n                 *
 *                                                                        *
 * ********************************************************************** */

#if defined(COMPILE_MASTER)

#include <ilcplex/cplexremotemasterx.h>



#include <limits.h>
#include <unistd.h>
#define MAX_PATH_LEN PATH_MAX
#define millisleep(x) usleep((x) * 1000)

#define LOGFILE                 "./log/log.txt"
static double b_time = 0;
static double g_time = 0;

/* Time get func */\
static double
printtime (int fg, FILE *ptlog){
	struct timeval s_time;
	double time;
	gettimeofday(&s_time, 0);
	time = (double)(s_time.tv_sec*1000000 + s_time.tv_usec);
	time = time/1000000 - b_time;
    if ( fg == 1 ){
	    printf("The clock time is %f sec.\n", time);
        fseek(ptlog, 0, SEEK_END);
	    fprintf(ptlog, "The clock time is %f sec.\n", time);
    }
	return time;
}

struct paramvalue {
	int num;  /**< Parameter number. */
	int type; /**< Parameter type. */
	struct {
		CPXINT  i;
		CPXLONG l;
		double  d;
	} value;  /**< The value for the parameter. This should better be a
				*   union but some compilers do not support a simple way of
				*   statically initializing a union so we use a struct.
				*/
};
#define INTPARAM(num,val)    { num, CPX_PARAMTYPE_INT, { val, 0, 0 } }
#define LONGPARAM(num,val)   { num, CPX_PARAMTYPE_LONG, { 0, val, 0 } }
#define DOUBLEPARAM(num,val) { num, CPX_PARAMTYPE_DOUBLE, { 0, 0, val } }
#define ENDPARAM             { -1,  CPX_PARAMTYPE_NONE, { 0, 0, 0 } }
#define GOLBALGAP				(primal.bound - dual.bound)*100/primal.bound
#define LOCALGAP				(rs->primal - rs->dual)*100/rs->primal

/** Parameter settings to have CPLEX work mainly on the primal bound. */
static struct paramvalue primalbound[] = {
	LONGPARAM (CPXPARAM_MIP_Strategy_HeuristicFreq, 1), /* Run heuristics at each node. */
	LONGPARAM (CPXPARAM_MIP_Strategy_RINSHeur, 1), /* Run RINS every two nodes. */
	//LONGPARAM (CPXPARAM_MIP_Limits_CutPasses, -1), /* Disable cuts. */
	ENDPARAM
};

static struct paramvalue defaultbound[] = {
	ENDPARAM
};

/** Predefined parameter settings. */
static struct {
	char const        *name;
	struct paramvalue const *values;
} const settings[] = {
	{ "default", defaultbound },
	{ "primal only", primalbound }
};
#define NUMSETTINGS ((int)(sizeof (settings) / sizeof (settings[0])) - 1)

/** Apply a predefined parameter settings.
 * \param env     The environment to which the predefined setting is applied.
 * \param setting Index of the setting to apply.
 */
static void
applySettings (CPXENVptr env, int setting)
{
	int status = 0;
    struct paramvalue const *vals;
    if ( setting == 0 )
	    vals = settings[0].values;
    else
	    vals = settings[1 + setting % NUMSETTINGS].values;

	status = CPXXsetintparam (env, CPXPARAM_RandomSeed, setting);
	if ( status ) {
		abort ();
	}

    status = CPXXsetintparam (env, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
	if ( status ) {
		abort ();
	}

	while (vals->type != CPX_PARAMTYPE_NONE) {
		switch (vals->type) {
		case CPX_PARAMTYPE_INT:
			status = CPXXsetintparam (env, vals->num, vals->value.i);
			break;
		case CPX_PARAMTYPE_LONG:
			status = CPXXsetlongparam (env, vals->num, vals->value.l);
			break;
		case CPX_PARAMTYPE_DOUBLE:
			status = CPXXsetdblparam (env, vals->num, vals->value.d);
			break;
		}
		if ( status ) {
			abort ();
		}
		++vals;
	}
}


/** Install/reset the info callback at the remote end.
 */
static int installCallback (CPXENVptr env)
{
	return CPXXuserfunction (env, USERACTION_ADDCALLBACK, 0, NULL, 0, 0, NULL);
}

/** Remove the info callback at the remote end.
 */
static int removeCallback (CPXENVptr env)
{
	return CPXXuserfunction (env, USERACTION_REMOVECALLBACK, 0, NULL, 0, 0, NULL);
}

/** Change the tolerance that controls whether consecutive objective
 * function values are considered different.
 */
static int changeObjdiff (CPXENVptr env, double newdiff)
{
	int status;
	CPXSERIALIZERptr s = NULL;

	CPXXserializercreate (&s);
	s->adddouble (s, newdiff);
	status = CPXXuserfunction (env, USERACTION_CHANGEOBJDIFF,
								CPXXserializerlength (s),
								CPXXserializerpayload (s),
								0, 0, NULL);
	CPXXserializerdestroy (s);

	return status;
}

/** Description of best known primal bound. */
static struct {
	int volatile     valid; /**< True if the bound/env fields are valid. */
	double volatile  bound; /**< Best known primal bound. */
	int volatile     idx;   /**< The environment that reported bound. */
} primal = { 0, 0.0, -1 };

/** Description of best known dual bound. */
static struct {
	int volatile     valid; /**< True if the bound/env fields are valid. */
	double volatile  bound; /**< Best known dual bound. */
	int volatile     idx;   /**< The environment that reported bound. */
} dual = { 0, 0.0, -1 };

static int objsen;

/** Current status of a remote worker. */
typedef struct remotestat {
	CPXENVptr env;
	double    primal;
	double    dual;
	double    gap;
	int       idx;
} RmtStat, *RmtStatPtr;
static RmtStatPtr remotestats = NULL;

static int jobs = 0;

static int GapCompare (const void *x, const void *y){
    double gap_x = ((RmtStatPtr)x)->gap;
    double gap_y = ((RmtStatPtr)y)->gap;
    return (gap_x > gap_y) ? 1 : -1;
}

/* This function is invoked whenever a remote worker reports new
 * information.
 * The information reported by remote workers are:
 * - new dual bounds,
 * - new primal bounds,
 * The handler functions updates the remotestats[] entry for the respective
 * remote worker and also updates that global best known primal/dual bound
 * if necessary.
 */
static void CPXPUBLIC
infohandler (CPXENVptr xenv, CPXINFOTYPE type, int tag, CPXLONG length,
             void const *data, void *handle)
{
	struct remotestat *const rs = handle;
	double d;
    FILE *iflog = fopen(LOGFILE, "r+");
	(void)xenv;
	(void)type;
	(void)length;

	switch (tag) {
	case INFO_NEWDUAL:
		assert (type == CPXINFO_DOUBLE);
		assert (length == 1);
		d = *(double const *)data;
		printf ("[%d] New dual bound: %e\n", rs->idx, d);
        fseek(iflog, 0, SEEK_END);
		fprintf (iflog, "[%d] New dual bound: %e\n", rs->idx, d);
		fflush (stderr);
		rs->dual = d;
		rs->gap = LOCALGAP;
		if ( !dual.valid ||
			(objsen == CPX_MIN && d > dual.bound) ||
			(objsen == CPX_MAX && d < dual.bound) )
		{
			dual.valid = 1;
			dual.bound = *(double const *)data;
			dual.idx = rs->idx;
			printf ("\033[1m\033[33m[%d] Bound updated: (\033[32m%e\033[33m, %e),\033[0m", 
                    rs->idx, dual.bound, primal.bound);
            printf ("\033[1m\033[33m gap: %3.2f%%. \033[0m\n", GOLBALGAP);
            fseek(iflog, 0, SEEK_END);
	        fprintf (iflog, "[%d] Bound updated: (%e, %e), gap: %3.2f%%.\n", 
                     rs->idx, dual.bound, primal.bound, GOLBALGAP);
		}
		break;
	case INFO_NEWPRIMAL:
		assert (type == CPXINFO_DOUBLE);
		assert (length == 1);
		d = *(double const *)data;
		printf ("[%d] New primal bound: %e\n", rs->idx, d);
        fseek(iflog, 0, SEEK_END);
		fprintf (iflog, "[%d] New primal bound: %e\n", rs->idx, d);
		fflush (stderr);
		rs->primal = d;
		rs->gap = LOCALGAP;
		if ( !primal.valid ||
			(objsen == CPX_MIN && d < primal.bound) ||
			(objsen == CPX_MAX && d > primal.bound) )
		{
			primal.valid = 1;
			primal.bound = *(double const *)data;
			primal.idx = rs->idx;
			printf ("\033[1m\033[33m[%d] Bound updated: (%e, \033[32m%e\033[33m),\033[0m",
					rs->idx, dual.bound, primal.bound);
            printf ("\033[1m\033[33m gap: %3.2f%%. \033[0m\n", GOLBALGAP);
            fseek(iflog, 0, SEEK_END);
	        fprintf (iflog, "[%d] Bound updated: (%e, %e), gap: %3.2f%%.\n", 
                     rs->idx, dual.bound, primal.bound, GOLBALGAP);
		}
		break;
	default:
		fflush (stderr);
	}
	printtime(1, iflog);
    fclose(iflog);
}

/** Function used for prefixed output.
 *  The function just prints \a message prefixed by the value of \a handle.
 */
static void CPXPUBLIC
printer (void *handle, char const *message)
{
   printf ("[%p] %s", handle, message);
}


int
main(int argc, char **argv){
    int status;
    CPXENVptr *env, l_env, t_env;
    CPXLPptr *lp, l_lp;
    char const *modelfile = NULL;
    CPXASYNCptr *handle;
    CPXENVGROUPptr group;
    int i;
    int max_jobs = 0;
    int active;
    int *finished;
    char const **machine;
    char cwd[MAX_PATH_LEN];
    char usrfunc[MAX_PATH_LEN];
    double time_rc;
    double gap = 0.01, diff_obj = 0;
    FILE *log = fopen(LOGFILE, "w+");
    int bestidx;
    CPXDIM c, cols;
    double *x;
    enum {
      OUTPUT_SILENT, OUTPUT_PREFIXED, OUTPUT_LOG
    } output = OUTPUT_SILENT;

#if defined(USE_MPI)
    int numprocs, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    if ( numprocs < 3 ) {
        fprintf (stderr, "Invalid number of processors (%d)\n", numprocs);
        abort ();
    }
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if ( rank != 0 ) {
        fprintf (stderr, "Master must have rank 0!\n");
        MPI_Finalize ();
        abort ();
    }
    machine = malloc (sizeof (*machine) * numprocs);
    if ( machine == NULL ) {
        fprintf (stderr, "Out of memory!\n");
        abort ();
    }
    for (i = 0; i < numprocs; ++i)
        machine[i] = "mpimachine";
    jobs = numprocs - 1;
#elif defined(USE_PROCESS)
    char const *bin = "./cplex";

    machine = malloc (sizeof (*machine) * argc);
    if ( machine == NULL ) {
        fprintf (stderr, "Out of memory!\n");
        abort ();
    }
#elif defined(USE_TCPIP)
   machine = malloc (sizeof (*machine) * argc);
   if ( machine == NULL ) {
      fprintf (stderr, "Out of memory!\n");
      abort ();
   }
#else
#   error "No transport type selected"
#endif

    /* Parse the command line. */
    for (i = 1; i < argc; ++i) {
        if ( strncmp (argv[i], "-model=", 7) == 0 )
            modelfile = argv[i] + 7;
#if defined(USE_TCPIP)
        else if ( strncmp (argv[i], "-address=", 9) == 0 )
            machine[max_jobs++] = argv[i];
#endif
        else if ( strncmp (argv[i], "-gap=", 5) == 0 )
            gap = strtod (argv[i] + 5, NULL);
    }

    /* Get the base time. */
    b_time = printtime(0, log);

    /* Validate arguments. */
    if ( modelfile == NULL ) {
        fprintf (stderr, "No model file specified with -model=<modelfile>\n");
        abort ();
    }

    if ( max_jobs < 1 ) {
        fprintf (stderr, "Invalid job count %d\n", max_jobs);
        abort ();
    }

	/* Dynamic allocate worker env. */
    char const *transport = "tcpiptransport";
    char const *args[16];
	for (i = 0; i < max_jobs; ++i){
		*args = machine[i];
		t_env = CPXXopenCPLEXremote(transport, 1, args, &status);
        if ( (status == 0) && (t_env != NULL) ) {
			printf ("Creating env on %s\n", machine[i]);
			fseek(log, 0, SEEK_END);	
			fprintf (log, "Creating env on %s\n", machine[i]);
            jobs++;
			if ( (env = realloc (env, sizeof(*env) * jobs)) == NULL ){
				fprintf (stderr, "Out of memory!\n");
				abort ();
			}
			env[jobs-1] = t_env;
        }
	}

    /* Allocate working arrays. */
    if ( (handle = malloc (sizeof (*handle) * jobs)) == NULL ||
        (lp = malloc (sizeof (*lp) * jobs)) == NULL ||
        (finished = calloc (jobs, sizeof (*finished))) == NULL ||
        (remotestats = calloc (jobs, sizeof (*remotestats))) == NULL )
    {
        fprintf (stderr, "Out of memory!\n");
        abort ();
    }

    /* Find the place at which to find the shared object that implements
     * the user function. On Windows the path to the current directory is
     * likely to contain blanks, so better quote it.
     */
    getcwd (cwd, sizeof (cwd));
    usrfunc[0] = 0;

    strcat (usrfunc, "-libpath=");
    strcat (usrfunc, cwd);

    /* Create a remote object instances. */
    for (i = 0; i < jobs; ++i) {
        /* Create empty problem object for this remote solver. */
        printf ("Creating LP %d\n", i);
	    fseek(log, 0, SEEK_END);	
        fprintf (log, "Creating LP %d\n", i);
        lp[i] = CPXXcreateprob(env[i], &status, "problem");
        if ( status || lp[i] == NULL ) {
            fprintf (stderr, "CPXXcreateprob: %d\n", status);
            abort ();
        }

        /* Install and configure callbacks. */
        remotestats[i].env = env[i];
        remotestats[i].idx = i;
        if ( (status = CPXXsetinfohandler (env[i], infohandler, &remotestats[i])) != 0 ) {
            fprintf (stderr, "CPXXsetinfohandler: %d\n", status);
            abort ();
        }

        if ( (status = changeObjdiff (env[i], 1e-5)) != 0 ) {
            fprintf (stderr, "changeObjdiff: %d\n", status);
            abort ();
        }

        if ( (status = installCallback (env[i])) != 0 ) {
            fprintf (stderr, "installCallback: %d\n", status);
            abort ();
        }

        /* Apply predefined perameter settings for this solver. */
        applySettings (env[i], i);
    }

    /* Create the local env and lp to presolve the problem. */
    l_env = CPXXopenCPLEX(&status);
    printf("Creating local env...\n");
    if ( status || l_env == NULL ) {
        fprintf (stderr, "CPXXopenCPLEX: %d\n", status);
        abort ();
    }

    printf("Creating local problem...\n");
    l_lp = CPXXcreateprob(l_env, &status, "presolve");
    if ( status || l_lp == NULL ){
        fprintf (stderr, "CPXXcreateprob: %d\n", status);
        abort ();       
    }

    printf("Opening the problem %s...\n", modelfile);
    if ( (status = CPXXreadcopyprob(l_env, l_lp, modelfile, NULL)) != 0 ){
        fprintf (stderr, "CPXXreadcopyprob: %d\n", status);
        abort ();       
    }

    printf("Presolving the proble...\n");
    if ( (status = CPXXpresolve(l_env, l_lp, CPX_ALG_NONE)) != 0 ){
        fprintf (stderr, "CPXXpresolve: %d\n", status);
        abort ();
    }

    printf("Rewriting the presolved problem...\n");
    if ( (status = CPXXpreslvwrite(l_env, l_lp, "./data/preslvprob.gz", &diff_obj)) != 0 ){
        fprintf (stderr, "CPXXpreslvwrite: %d\n", status);
        abort ();
    }

    CPXXfreeprob (l_env, l_lp);
    CPXXcloseCPLEX (&l_env);
    free (l_env);

    /* Put all environments into one group. */
    if ( (status = CPXXcreateenvgroup (&group, jobs, env)) != 0 ) {
        fprintf (stderr, "CPXXcreateenvgroup: %d\n", status);
        abort ();
    }

    /* Read the model into all remote solver. */
    printf ("Reading and broadcasting the model...\n");
	fseek(log, 0, SEEK_END);	
    fprintf (log, "Reading and broadcasting the model.\n");
    if ( (status = CPXXreadcopyprob_multicast (group, "./data/preslvprob.gz", NULL)) != 0 ) {
        fprintf (stderr, "CPXXreadcopyprob_multicast: %d\n", status);
        abort ();
    }
    objsen = CPXXgetobjsen (env[0], lp[0]);
   

    /* Set the threads of each machine */
    printf ("Setting the params.\n");
    fseek(log, 0, SEEK_END);
    fprintf (log, "Setting the params.\n");
    if ( (status = CPXXsetintparam_multicast (group, CPXPARAM_Threads, 2)) != 0 ) {
        fprintf (stderr, "CPXXsetintparam_multicast: %d\n", status);
        abort ();
    }

    /* Start an asynchronous solve on each remote solver. */
    for (i = 0; i < jobs; ++i) {
        printf ("Solving %d\n", i);
	    fseek(log, 0, SEEK_END);	
        fprintf (log, "Solving %d\n", i); 
        if ( (status = CPXXmipopt_async (env[i], lp[i], &handle[i])) != 0 ) {
            fprintf (stderr, "CPXXmipopt_async: %d\n", status);
            abort ();
        }
    }

    /* All solves are started. Loop until the stopping criterion is met. */
    active = jobs;
    int cnt = 0;
    while (active > 0) {
        int running = 0;
        /* We stop them if the absolute mipgap is reached. */
        if (  primal.valid && dual.valid &&
            ((objsen == CPX_MIN && (primal.bound - dual.bound)/dual.bound <= gap) ||
             (objsen == CPX_MAX && (dual.bound - primal.bound)/primal.bound <= gap)) ){
            printf ("Stopping criterion reached. Stopping all pending solves.\n");
	        fseek(log, 0, SEEK_END);	
            fprintf (log, "Stopping criterion reached. Stopping all pending solves.\n");
		    time_rc = printtime(1, log);
            for (i = 0; i < jobs; ++i) {
                if ( !finished[i] )
                CPXXasynckill (handle[i]);
            }
            break;
        }

        /* Loop over all solvers and test if they are still running. */
        for (i = 0; i < jobs; ++i) {
            if ( finished[i] )
            continue;
            CPXXasynctest (handle[i], &running);
            if ( !running ) {
            /* The job is finished. Kill all others. */
                int j;
                --active;
                finished[i] = 1;
                printf ("First job (%d) is finished, killing the rest\n", i);
	            fseek(log, 0, SEEK_END);	
                fprintf (log, "First job (%d) is finished, killing the rest\n", i);
                time_rc = printtime(1, log);
                for (j = 0; j < jobs; ++j) {
                    if ( j != i )
                        CPXXasynckill (handle[j]);
                }
                break;
            }
        }

		/* Resort the remotestat according to the local gap. */
        if ( (cnt++) == 1000 ){
            qsort(remotestats, jobs, sizeof(struct remotestat), GapCompare);
            for (i = 0; i < jobs; ++i) 
                printf("[%d] worker %d, gap: %f\n", 
                       i, remotestats[i].idx, remotestats[i].gap);
            cnt = 0;
        }

        millisleep (10);
    }

    /* All solves have finished. Join them. */
    for (i = 0; i < jobs; ++i) {
        double obj = -CPX_INFBOUND;
        int stat;
      
        status = CPXXmipopt_join (&handle[i]);
        if ( status ) {
            fprintf (stderr, "CPXXmipopt_join: %d\n", status);
            abort ();
        }

        status = CPXXgetobjval (env[i], lp[i], &obj);
        if ( status == CPXERR_NO_SOLN ) {
            /* No feasible solution found (yet) on this machine. */
            obj = (objsen == CPX_MIN) ? CPX_INFBOUND : -CPX_INFBOUND;
        }
        else if ( status ) {
            fprintf (stderr, "CPXXgetobjval: %d\n", status);
            abort ();
        }
        stat = CPXXgetstat (env[i], lp[i]);
        printf ("Job %d: %f, stat %d\n", i, obj, stat);
        printf ("\t%e, %e\n", remotestats[i].dual, remotestats[i].primal);
	    fseek(log, 0, SEEK_END);	
        fprintf (log, "Job %d: %f, stat %d\n", i, obj, stat);
        fprintf (log, "\t%e, %e\n", remotestats[i].dual, remotestats[i].primal);
        if ( (status = removeCallback (env[i])) != 0 ) {
            fprintf (stderr, "removeCallback: %d\n", status);
            abort ();
        }
    }

	/* Fetch the x vector from the solver that produced the best
	* primal bound. */
	/*
	bestidx = primal.idx;
	cols = CPXXgetnumcols (env[bestidx], lp[bestidx]);
	if ( (x = malloc (cols * sizeof (*x))) == NULL ) {
		fprintf (stderr, "Out of memory!\n");
		abort ();
	}
	status = CPXXgetx (env[bestidx], lp[bestidx], x, 0, cols - 1);
	if ( status ) {
		fprintf (stderr, "CPXXgetx: %d\n", status);
		abort ();
	}
	printf ("Optimal solution:\n");
	for (c = 0; c < cols; ++c)
		printf ("x[%5d]: %f\n", c, x[c]);
	free (x);
	*/
    CPXXfreeenvgroup (&group);

    /* Close the CPLEX objects in _reverse_ order. */
    for (i = jobs - 1; i >= 0; --i)
        CPXXcloseCPLEX (&env[i]);
    free (remotestats);
    free (env);
    free (lp);
    free (handle);
    free (finished);
    free ((char **)machine);

    fclose(log);
    return 0;
}

#endif /* COMPILE_MASTER */
