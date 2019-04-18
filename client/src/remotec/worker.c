/* --------------------------------------------------------------------------
 * File: worker.c
 * Version 1.0.2
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
 * Copyright IBM Corporation, 2012, 2019. All Rights Reserved.
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
 *    C l i e n t   s i d e   u s e r f u n c t i o n                     *
 *                                                                        *
 * ********************************************************************** */

#if defined(COMPILE_USERFUNCTION)

#include <ilcplex/cplexremoteworkerx.h>

#include <pthread.h>
#define MUTEX pthread_mutex_t
#define MUTEX_INIT(mtx)    pthread_mutex_init ((mtx), NULL)
#define MUTEX_DESTROY(mtx) pthread_mutex_destroy (mtx)
#define MUTEX_LOCK(mtx)    pthread_mutex_lock (mtx)
#define MUTEX_UNLOCK(mtx)  pthread_mutex_unlock (mtx)

#define BESTGAP            (best.primal - best.dual)*100/best.primal
/** Best known primal and dual bounds in this remote worker.
 *  We need to keep track of these values so that we only report bound
 *  improvements to the master process.
 */
static struct {
    MUTEX           objmutex;   /**< Mutex for synchronizing access to
                                 *   this structure. */
    int             haveDual;   /**< Is the 'dual' field valid? */
    int             havePrimal; /**< Is the 'primal' field valid? */
    double          dual;       /**< Best known dual bound in this solver. */
    double          primal;     /**< Best known primal bound in this solver. */
    double          gap;        /**< Best known gap in this solver. */
    double          objdiff;    /**< Minimal delta between consecutive bounds. */
} best;

/** MIP info callback that is registered with CPLEX.
 *  This callback picks up primal and dual bounds as well as the current
 *  deterministic time. If bounds changed then updated bounds are send
 *  to the master. Deterministic time stamps are always send to the master.
 */
static int CPXPUBLIC
infocallback (CPXCENVptr cbenv, void *cbdata, int wherefrom, void *cbhandle)
{
    CPXCENVptr env = cbhandle;
    double dual, primal, ts;
    double BeginGap = 100;

    /* Test if we have improved the primal bound and report if so. */
    if ( CPXXgetcallbackinfo (cbenv, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &primal) == 0 ) {
        MUTEX_LOCK (&best.objmutex);
        if ( !best.havePrimal || fabs (best.primal - primal) > best.objdiff ) {
            best.havePrimal = 1;
            best.primal = primal;
            if ( best.havePrimal && best.haveDual ){
                best.gap = BESTGAP;
                if ( best.gap <= BeginGap )
                    (void)CPXXsendinfodouble (env, INFO_NEWPRIMAL, 1, &primal);
            }
        }
        MUTEX_UNLOCK (&best.objmutex);
    }

   /* Test if we have improved the dual bound and report if so. */
    if ( CPXXgetcallbackinfo (cbenv, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &dual) == 0 ) {
        MUTEX_LOCK (&best.objmutex);
        if ( !best.haveDual || fabs (best.dual - dual) > best.objdiff ) {
            best.haveDual = 1;
            best.dual = dual;
            if ( best.havePrimal && best.haveDual ){
                best.gap = BESTGAP;
                if ( best.gap <= BeginGap )
                    (void)CPXXsendinfodouble (env, INFO_NEWDUAL, 1, &dual);
            }
        }
        MUTEX_UNLOCK (&best.objmutex);
    }

    return 0;
}

/** User function implementation.
 *  This function is executed when the master invokes a user function
 *  on this remote solver.
 */
static int CPXPUBLIC
userfunction (CPXENVptr env, int id, CPXLONG inlen,
              void const *indata, CPXLONG maxout,
              CPXLONG *outlen_p, void *outdata, void *handle)
{
    CPXDESERIALIZERptr d = NULL;
    int status = 0;

    (void)maxout;
    (void)outdata;
    (void)handle;

    *outlen_p = 0;
    CPXXdeserializercreate (&d, inlen, indata);

    switch (id) {
        case USERACTION_ADDCALLBACK:
            best.havePrimal = 0;
            best.haveDual = 0;
            if ( MUTEX_INIT (&best.objmutex) )
                status = -1;
            else
                (void)CPXXsetinfocallbackfunc (env, infocallback, env);
            break;
        case USERACTION_REMOVECALLBACK:
            (void)CPXXsetinfocallbackfunc (env, NULL, NULL);
            status = MUTEX_DESTROY (&best.objmutex);
            break;
        case USERACTION_CHANGEOBJDIFF:{
            double dbl;
            dbl = best.objdiff;
            d->getdouble (d, &dbl);
            best.objdiff = dbl;
            break;
	  }
    }

    CPXXdeserializerdestroy (d);

    return status;
}

/** Register user function handler.
 *  This function is invoked by the remote worker at startup to allow us
 * to register a user function handler.
 */
CPXEXPORT void CPXPUBLIC REGISTER_USERFUNCTION (struct messagehandler *handler);
CPXEXPORT void CPXPUBLIC
REGISTER_USERFUNCTION (struct messagehandler *handler)
{
    CPXXsetuserfunction (handler, userfunction, NULL);
}

#endif /* COMPILE_USERFUNCTION */
