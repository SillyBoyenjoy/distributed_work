/* --------------------------------------------------------------------------
 * File: master.c
 * Version 1.1.2
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
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <keyboard.h>

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
#include <pthread.h>
#define MAX_PATH_LEN PATH_MAX
#define MAX_CMD_LEN 64
#define millisleep(x) usleep((x) * 1000)

#define LOGFILE                 "./log/log.txt"
static double b_time = 0;

/* Time get func */
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

static int lflag = 1;   /* Loop Flag */

/* Command input func */
static void
commandinput (void){
    char command[MAX_CMD_LEN] = {0};
    init_keyboard();
    while ( lflag == 1 ){
        /* Detect the keyboard command. */
        if ( (kbhit() == 1) && (readch() == 10) ){
            printf(">");
            memset(command, 0, MAX_CMD_LEN);
            fgets(command, MAX_CMD_LEN, stdin);
            if ( strncmp(command, "quit", 4) == 0 ){
                lflag = 0;
                break;
            }
        }
    }
    close_keyboard();
    pthread_exit(NULL);
}

static int jobs = 0;    /* Current jobs */

/* Current status of a remote worker. */
typedef struct remotestat {
    int       idx;
	CPXENVptr env;
	double    primal;
	double    dual;
	double    gap;
    char const  *paramname;
} RmtStat, *RmtStatPtr;
static RmtStatPtr remotestats = NULL;

/* The list of best Integers and the list of best bounds. */
typedef struct boundlist {
    int idx;
    double bound;
    double gap;
} BoundList, *BoundListPtr;
static BoundListPtr BestPrimals = NULL, BestDuals = NULL;


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
#define GOLBALGAP			(primal.bound - dual.bound)*100/primal.bound
#define LOCALGAP			(rs->primal - rs->dual)*100/rs->primal
#define MAXDIFF			    (BestPrimals[jobs-1].bound - BestPrimals[0].bound)*100/BestPrimals[jobs-1].bound

static struct paramvalue defaultbound[] = {
	ENDPARAM
};

static struct paramvalue primalbound[] = {
	LONGPARAM (CPXPARAM_MIP_Strategy_HeuristicFreq, 1), /* Run heuristics at each node. */
	LONGPARAM (CPXPARAM_MIP_Strategy_RINSHeur, 1), /* Run RINS every two nodes. */
	//LONGPARAM (CPXPARAM_Emphasis_MIP, CPX_MIPEMPHASIS_BESTBOUND), /* Emphasis_bestbound. */
	//LONGPARAM (CPXPARAM_MIP_Limits_CutPasses, -1), /* Disable cuts. */
	ENDPARAM
};

static struct paramvalue dualbound[] = {
   LONGPARAM (CPXPARAM_MIP_Strategy_HeuristicFreq, -1), /* Heuristics off. */
   INTPARAM (CPXPARAM_MIP_Cuts_Cliques, 3),
   INTPARAM (CPXPARAM_MIP_Cuts_Covers, 3),
   INTPARAM (CPXPARAM_MIP_Cuts_Disjunctive, 3),
   INTPARAM (CPXPARAM_MIP_Cuts_FlowCovers, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_Gomory, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_GUBCovers, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_Implied, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_MIRCut, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_ZeroHalfCut, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_MCFCut, 2),
   ENDPARAM
};


/** Predefined parameter settings. */
static struct {
	char const        *name;
	struct paramvalue const *values;
} const settings[] = {
	{ "default", defaultbound },
	{ "primal only", primalbound },
    { "dual only", dualbound }
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
    if ( setting == 0 ){
	    vals = settings[2].values;
        remotestats[setting].paramname = settings[2].name;
        printf("Set worker %d param: %s.\n", 
               setting, remotestats[setting].paramname);
    }
    else{
	    vals = settings[1].values;
        remotestats[setting].paramname = settings[1].name;
        printf("Set worker %d param: %s.\n", 
               setting, remotestats[setting].paramname);
    }

    if( (status = CPXXsetintparam (env, CPXPARAM_RandomSeed, 10*setting)) != 0){
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
} primal = { 0, DBL_MAX, -1 };

/** Description of best known dual bound. */
static struct {
	int volatile     valid; /**< True if the bound/env fields are valid. */
	double volatile  bound; /**< Best known dual bound. */
	int volatile     idx;   /**< The environment that reported bound. */
} dual = { 0, 0.0, -1 };

static int objsen;

static int PriCompare (const void *x, const void *y){
    double primal_x = ((BoundListPtr)x)->bound;
    double primal_y = ((BoundListPtr)y)->bound;
    return (primal_x > primal_y) ? 1 : -1;
}

static int DualCompare (const void *x, const void *y){
    double dual_x = ((BoundListPtr)x)->bound;
    double dual_y = ((BoundListPtr)y)->bound;
    return (dual_x < dual_y) ? 1 : -1;
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
    int i;
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
        for ( i = 0; i < jobs; i++ ){
            if ( BestDuals[i].idx == rs->idx ){
                BestDuals[i].bound = rs->dual;
                BestDuals[i].gap = rs->gap;
                break;
            }
        }
		if ( (!dual.valid) || (d > dual.bound) )
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
        qsort(BestDuals, jobs, sizeof(struct boundlist), DualCompare);
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
        for ( i = 0; i < jobs; i++ ){
            if ( BestPrimals[i].idx == rs->idx ){
                BestPrimals[i].bound = rs->primal;
                BestPrimals[i].gap = rs->gap;
                break;
            }
        }
		if ( (!primal.valid) || (d < primal.bound) )
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
        qsort(BestPrimals, jobs, sizeof(struct boundlist), PriCompare);
		break;
	default:
		fflush (stderr);
	}
	printtime(1, iflog);

    /*
    for ( i = 0; i < jobs; i++)
        printf("worker %d: primal %e, gap %3.2f%%.\tworker %d: dual %e, gap %3.2f%%.\n",
              BestPrimals[i].idx, BestPrimals[i].bound, BestPrimals[i].gap,
              BestDuals[i].idx, BestDuals[i].bound, BestDuals[i].gap);
    */
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
    char const *modelfile = NULL, *startfile = NULL;
    char incumbent[MAX_PATH_LEN] = "./sol/icb_";
    CPXASYNCptr *handle;
    CPXENVGROUPptr group;
    int i, j, max_jobs = 0;
    int warm_start = 0, iteration = 0, cut = 0, polish = 0;
    int active, *conn;
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

    machine = malloc (sizeof (*machine) * argc);
    if ( machine == NULL ) {
        fprintf (stderr, "Out of memory!\n");
        abort ();
    }

    /* Parse the command line. */
    for (i = 1; i < argc; ++i) {
        if ( strncmp (argv[i], "-model=./data/", 14) == 0 ){
            modelfile = argv[i] + 7;
            printf("Modelfile is %s\n", modelfile);
            fseek(log, 0, SEEK_END);
            fprintf(log, "Modelfile is %s\n", modelfile);
            strcat(incumbent, argv[i] + 14);
            strcat(incumbent, ".sol.gz");
        }
        else if ( strncmp (argv[i], "-address=", 9) == 0 )
            machine[max_jobs++] = argv[i];
        else if ( strncmp (argv[i], "-gap=", 5) == 0 ){
            gap = strtod (argv[i] + 5, NULL);
            printf("Solve gap is %f\n", gap);
            fseek(log, 0, SEEK_END);
            fprintf(log, "Solve gap is %f\n", gap);
        }
        else if ( strncmp (argv[i], "-warmstart=./sol/", 17) == 0 ){
            warm_start = 1;
            startfile = argv[i] + 11;
            printf("Startfile is %s\n", startfile);
            fseek(log, 0, SEEK_END);
            fprintf(log, "Startfile is %s\n", startfile);
        }
        else if ( strncmp (argv[i], "-iteration=on", 13) == 0 ){
            iteration = 1;
            printf("Iteration: on.\n");
            fseek(log, 0, SEEK_END);
            fprintf(log, "Iteration: on.\n");
        }
        else if ( strncmp (argv[i], "-cut=on", 7) == 0 ){
            cut = 1;
            printf("Branch&Cut: on.\n");
            fseek(log, 0, SEEK_END);
            fprintf(log, "Branch&Cut: on.\n");
        }
        else if ( strncmp (argv[i], "-polish=on", 10) == 0 ){
            polish = 1;
            printf("Polish: on.\n");
            fseek(log, 0, SEEK_END);
            fprintf(log, "Polish: on.\n");
        }
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
            if ( jobs == 1 ){
                if ( (env = malloc (sizeof(*env))) == NULL ){
                    fprintf (stderr, "Out of memory!\n");
                    abort ();
                }
            }
            else{
                if ( (env = realloc (env, sizeof(*env) * jobs)) == NULL ){
                    fprintf (stderr, "Out of memory!\n");
                    abort ();
                }
            }
			env[jobs-1] = t_env;
        }
	}

    /* Allocate working arrays. */
    if ( (handle = malloc (sizeof (*handle) * jobs)) == NULL ||
        (lp = malloc (sizeof (*lp) * jobs)) == NULL ||
        (conn = calloc (jobs, sizeof (*conn))) == NULL ||
        (remotestats = calloc (jobs, sizeof (*remotestats))) == NULL ||
        (BestPrimals = calloc (jobs, sizeof (*BestPrimals))) == NULL ||
        (BestDuals = calloc (jobs, sizeof (*BestDuals))) == NULL )
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
        lp[i] = CPXXcreateprob(env[i], &status, "problem");
        if ( status || lp[i] == NULL ) {
            fprintf (stderr, "CPXXcreateprob: %d\n", status);
            abort ();
        }
        else
            conn[i] = 1;

        /* Install and configure callbacks. */
        remotestats[i].env = env[i];
        remotestats[i].idx = i;
        BestPrimals[i].idx = i;
        BestDuals[i].idx = i;
        remotestats[i].primal = DBL_MAX;
        BestPrimals[i].bound = DBL_MAX;
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

    CPXfreeprob (l_env, l_lp);
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
    printf ("Setting the params...\n");
    fseek(log, 0, SEEK_END);
    fprintf (log, "Setting the params.\n");
    if ( (status = CPXXsetintparam_multicast (group, CPXPARAM_Threads, 4)) != 0 ) {
        fprintf (stderr, "CPXXsetintparam_multicast: %d\n", status);
        abort ();
    }

    if ( (status = CPXXsetintparam_multicast (group, CPXPARAM_Advance, 1)) != 0 ) {
        fprintf (stderr, "CPXXsetintparam_multicast: %d\n", status);
        abort ();
    }

    if ( warm_start == 1 ){
        printf("Warm starting...\n");
        if ( (status = CPXXreadcopymipstarts_multicast (group, startfile)) != 0)
		    printf ("Warm start failed, error: %d. Normally starting...\n", status);
        if ( polish == 1 ){
            printf("Setting polish...\n");
            if ( (status = CPXXsetdblparam (env[0], CPXPARAM_MIP_PolishAfter_Time, 0)) != 0 )
                printf("Polish failed, error: %d.", status);
        }
    } 

    /* Start an asynchronous solve on each remote solver. */
    printf ("Solving...\n");
    fseek(log, 0, SEEK_END);	
    fprintf (log, "Solving...\n");
    for (i = 0; i < jobs; ++i) {
        if ( (status = CPXXmipopt_async (env[i], lp[i], &handle[i])) != 0 ) {
            fprintf (stderr, "CPXXmipopt_async: %d\n", status);
            abort ();
        }
    }
	
    /* Crete a new thread for command input. */
    pthread_t tfd;
    if ( (status = pthread_create(&tfd, NULL, (void*)commandinput, NULL)) != 0 ){
        fprintf(stderr, "Thread creation failed: %d\n", status);
        abort ();
    }

    /* All solves are started. Loop until the stopping criterion is met. */
	int cnt = 0;
    int running = 0, r1 = 0, r2 = 0; /* Running test bit */
    /* The main loop. */
    while ( 1 ) {
        /* QUIT */
        if ( lflag == 0 ){
            break;
        }
        
        /* We stop them if the absolute mipgap is reached. */
        if ( (primal.valid && dual.valid && GOLBALGAP <= 100*gap) ){
            printf ("Stopping criterion reached. Stopping all pending solves.\n");
		    time_rc = printtioe(1, log);
            lflag = 0;
            break;
        }

        /* Restart some died worker. */
        for (i = 0; i < jobs; ++i) {
            CPXXasynctest (handle[i], &running);
            if ( (!running) && (conn[i] == 1) && GOLBALGAP > 100*gap ){
                if ( (status = CPXXmipopt_async (env[i], lp[i], &handle[i])) != 0 ){
                    printf("Restart failed, worker %d disconnected.\n", i);
                    conn[i] = 0;
                }
            }
        }

        /* Branch and cut if golbalgap < 3.00%. */
        if ( (GOLBALGAP < 3.00) && (cut == 1) ){
            printf ("Branch and cut begin...\n");
            
            for ( i = 0; i < jobs; i++ ){
                CPXXasynckill (handle[i]);
                if ( (status = CPXXmipopt_join (&handle[i])) != 0 ){
                    fprintf (stderr, "CPXXmipopt_join job:%d: %d\n", i, status);
                    abort ();
                }
            }

            i = BestDuals[0].idx;
            printf("Worker %d is writing the best dual bound...\n", i);
            if ( (status = CPXXwriteprob(env[i], lp[i], "./data/BestDual.lp.gz", NULL)) != 0 ){
                fprintf (stderr, "CPXXwriteprob: %d\n", status);
                abort ();
            }
            printf ("Worker writing completed.\n");

            i = BestPrimals[0].idx;
            printf("Worker %d is writing the best primal bound...\n", i);
            if ( (status = CPXXsolwrite(env[i], lp[i], "./data/BestPrimal.sol.gz")) != 0 ){
                fprintf (stderr, "CPXXsolwrite: %d\n", status);
                abort ();
            }
            printf ("Worker writing completed.\n");

            printf ("Workers restarting...\n");
            for ( i = 0; i < jobs; i++ ){
                printf("Worker %d is restart...\n", i);
                if ( i == 0 ){                /* A dual worker every 5 */
                    status = CPXXreadcopyprob(env[i], lp[i], "./data/BestDual.lp.gz", NULL);
                    if ( status != 0 ){
                        fprintf (stderr, "CPXXreadcopyprob: %d\n", status);
                        abort ();
                    }
                }
                else{
                    status = CPXXreadcopymipstarts(env[i], lp[i], "./data/BestPrimal.sol.gz");
                    if ( status != 0 ){
                        fprintf (stderr, "CPXXreadcopymipstarts: %d\n", status);
                        abort ();
                    }
                }
                if ( (status = CPXXsetdblparam(env[i], CPXPARAM_MIP_PolishAfter_Time, 0)) != 0 ){
                    fprintf (stderr, "CPXX set polish time: %d\n", status);
                    abort ();
                }
                do {
                    CPXXmipopt_async (env[i], lp[i], &handle[i]);
                    millisleep (500);
                    CPXXasynctest (handle[i], &running);                   
                }
                while ( !running );
            }
            printf("Workers restart completed.\n");
            cut = 0;            /* Turn off branch and cut */
            iteration = 0;      /* Turn off iteration */
        }

        /* Iterate the worst integers solution every 5 min. */
		if ( ((cnt++) > 30000) && (MAXDIFF > 200.00) && 
             (iteration == 1) && (GOLBALGAP < 10.00) ){
			printf ("Iteration begin...\n");

			i = BestPrimals[0].idx;
            printf ("Worker %d is writing sol and restarting...\n", i);
	        fseek(log, 0, SEEK_END);	
            fprintf (log, "Worker %d is writing sol and restarting...\n", i);
            CPXXasynckill (handle[i]);
            if ( (status = CPXXmipopt_join (&handle[i])) != 0 ){
				fprintf (stderr, "CPXXmipopt_join: %d\n", status);
				abort ();
			} 
            if ( (status = CPXXsolwrite( env[i], lp[i], incumbent )) != 0 ){
				fprintf (stderr, "CPXXsolwrite: %d\n", status);
				abort ();
			}
            if ( (status = CPXXmipopt_async (env[i], lp[i], &handle[i])) != 0 ) {
				fprintf (stderr, "CPXXmipopt_async: %d\n", status);
				abort ();
			}
            printf ("Worker writing completed.\n");

			j = BestPrimals[jobs-1].idx;
            printf ("Worker %d is reading sol and restarting...\n", j);
	        fseek(log, 0, SEEK_END);	
            fprintf (log, "Worker %d is reading sol and restarting...\n", j);
			CPXXasynckill (handle[j]);
			if ( (status = CPXXmipopt_join (&handle[j])) != 0 ){
				fprintf (stderr, "CPXXmipopt_join: %d\n", status);
				abort ();
			}
            if ( (status = CPXXreadcopymipstarts(env[j], lp[j], incumbent)) != 0){
				fprintf (stderr, "CPXXreadcopymipstarts: %d\n", status);
				abort ();
			}
            if ( (status = CPXXsetdblparam(env[j], CPXPARAM_MIP_PolishAfter_Time, 0)) != 0){
				fprintf (stderr, "CPXX set polish time: %d\n", status);
				abort ();
			}
			if ( (status = CPXXmipopt_async (env[j], lp[j], &handle[j])) != 0 ) {
				fprintf (stderr, "CPXXmipopt_async: %d\n", status);
				abort ();
			}
            do{
                CPXXasynctest (handle[i], &r1);
                CPXXasynctest (handle[j], &r2);
                if ( !r1 )
                    CPXXmipopt_async (env[i], lp[i], &handle[i]);
                if ( !r2 )
                    CPXXmipopt_async (env[j], lp[j], &handle[j]);
                millisleep (500);
            }
            while ( (!r1) || (!r2) );
            printf ("Worker reading completed.\n");

			printf ("Iteration completed.\n");
			cnt = 0;
		}
        millisleep (10);
    }

    /* Join the thread */
    if ( (status = pthread_join(tfd, NULL)) != 0 ){
        fprintf(stderr, "Pthread join error: %d\n", status);
        abort ();
    }

    /* Kill all workers. */
    for ( i = 0; i < jobs; ++i ){
        printf("Stopping worker %d...\n", i);
        do{
            CPXXasynckill (handle[i]);
            millisleep (1000);
            CPXXasynctest (handle[i], &running);
        }
        while ( running );
    }

    /* All solves have finished. Join them. */
    for (i = 0; i < jobs; ++i) {
        if ( (status = CPXXmipopt_join (&handle[i])) != 0 ) {
            fprintf (stderr, "CPXXmipopt_join: %d\n", status);
            abort ();
        }
    }

    i = BestPrimals[0].idx;
    if ( (status = CPXXsolwrite( env[i], lp[i], "./sol/bestsol.sol" )) != 0 ){
		fprintf (stderr, "CPXXsolwrite: %d\n", status);
		abort ();
	}

    CPXXfreeenvgroup (&group);

    /* Close the CPLEX objects in _reverse_ order. */
    for (i = jobs - 1; i >= 0; --i)
        CPXXcloseCPLEX (&env[i]);
    free (remotestats);
    free (BestPrimals);
    free (BestDuals);
    free (env);
    free (lp);
    free (handle);
    free (conn);
    free ((char **)machine);

    fclose(log);
    return 0;
}

#endif /* COMPILE_MASTER */
