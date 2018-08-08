#define OLD_CODE
#ifdef NEW_CODE

/**************************************************************************
PBEAMPP.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions 
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.           
Copyright (c) 2000-2002 ZIB & Loebel.  
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:04 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pbeampp.c,v 1.10 2005/02/17 19:42:32 bzfloebe Exp $  */





#include "pbeampp.h"
#include <processthreadsapi.h>
#include <omp.h>





#ifdef _PROTO_
int bea_is_dual_infeasible( arc_t *arc, cost_t red_cost )
#else
int bea_is_dual_infeasible( arc, red_cost )
    arc_t *arc;
    cost_t red_cost;
#endif
{
    return(    (red_cost < 0 && arc->ident == AT_LOWER)
            || (red_cost > 0 && arc->ident == AT_UPPER) );
}







typedef struct basket
{
    arc_t *a;
    cost_t cost;
    cost_t abs_cost;
} BASKET;

static long basket_size;
static BASKET basket[B+K+1];
static BASKET *perm[B+K+1];



#ifdef _PROTO_
void sort_basket( long min, long max )
#else
void sort_basket( min, max )
    long min, max;
#endif
{
    long l, r;
    cost_t cut;
    BASKET *xchange;

    l = min; r = max;

    cut = perm[ (long)( (l+r) / 2 ) ]->abs_cost;

    do
    {
        while( perm[l]->abs_cost > cut )
            l++;
        while( cut > perm[r]->abs_cost )
            r--;
            
        if( l < r )
        {
            xchange = perm[l];
            perm[l] = perm[r];
            perm[r] = xchange;
        }
        if( l <= r )
        {
            l++; r--;
        }

    }
    while( l <= r );

    if( min < r )
        sort_basket( min, r );
    if( l < max && l <= B )
        sort_basket( l, max ); 
}

static long nr_group;
static long group_pos;

int max_basket_size = 0;
static long initialize = 1;





#ifdef _PROTO_
arc_t *primal_bea_mpp( long m,  arc_t *arcs, arc_t *stop_arcs, 
                              cost_t *red_cost_of_bea )
#else
arc_t *primal_bea_mpp( m, arcs, stop_arcs, red_cost_of_bea )
    long m;
    arc_t *arcs;
    arc_t *stop_arcs;
    cost_t *red_cost_of_bea;
#endif
{
    long i, next, old_group_pos;
    arc_t *arc;
    cost_t red_cost;
    resHolder resArr[MAX_BASKET_SIZE];
    for (int j = 0; j < MAX_BASKET_SIZE; ++j) {
        resArr[j].arrIndex = 0;
        resArr[j].red_cost = 0;
    }

    if( initialize )
    {
        for( i=1; i < K+B+1; i++ )
            perm[i] = &(basket[i]);
        nr_group = ( (m-1) / K ) + 1;
        group_pos = 0;
        basket_size = 0;
        initialize = 0;
    }

    else
    {
        //  Used to perserve the red cost condition result
        int ConditionResultsArray[MAX_BASKET_SIZE];
        for (int l = 0; l < MAX_BASKET_SIZE; ++l) {
            ConditionResultsArray[l] = 0;
        }
        //  Used to keep track of where the results should be written into
        //  For example if their were 3 results for the second thread, then
        //  we know that in the results array we will write 3 serial results
        //  in the position right after the first X results of the first thread.
        int ResultsPerThread[NUM_OF_THREADS];
        for (int k = 0; k < NUM_OF_THREADS; ++k) {
            ResultsPerThread[k] = 0;
        }

        //  Set the minimum size that we will iterate, max size of it is B
        //  Due to the iteration size in the original for loop  
        int minimum_size_of_iteration;
        if(B < basket_size)
            minimum_size_of_iteration = B;
        else
            minimum_size_of_iteration = basket_size;

        //  Calculate the samples per thread, according to the original for loop
        //  its the minimum between B and basket size + 1 minus the 2 since
        //  the iteration index started at 2
        int samples_per_thread = (minimum_size_of_iteration - 2 + 1) / NUM_OF_THREADS;

//        printf("Starting parralized first loop\n");
        #pragma omp parallel for num_threads(NUM_OF_THREADS)
        for (int j = 0; j < NUM_OF_THREADS; j++)
        {

            //  Now each thread needs to run on his own samples.
            //  The original loop started its iteration at 2, so we
            //  need to calculate each thread start and end interval.
            int thread_sample_start_index = j * samples_per_thread + 2;

            // TODO (Eli): think about the following issue
            // Here at the end index, their might be left overs that the last
            //  thread might not take care of.
            int thread_sample_end_index = thread_sample_start_index + samples_per_thread;

//            printf("thread_sample_start_index: %d, thread_sample_end_index: %d\n", thread_sample_start_index, thread_sample_end_index);

            for (int k = thread_sample_start_index; k < thread_sample_end_index; k++) {
                //  Calculate arc and red_cost as it was in the original loop
                //  with the relative k
                arc = perm[k]->a;
                red_cost = arc->cost - arc->tail->potential + arc->head->potential;

                if((red_cost < 0 && arc->ident == AT_LOWER)
                   || (red_cost > 0 && arc->ident == AT_UPPER)) {
                    resArr[k].red_cost = red_cost;
                    resArr[k].arc      = arc;
                    //  Keep the if condition result
                    ConditionResultsArray[k] = 1;
                    //  Sum up the number of times the condition was met
                    ResultsPerThread[j] += 1;
                }
            }
        }

        #pragma omp parallel for num_threads(NUM_OF_THREADS)
        for (int j = 0; j < NUM_OF_THREADS; j++)
        {
            int CurrentThreadResultPositionTracker = 0;

            //  This variable is used to point where the thread starts to write
            //  results. For example if their was 3 good results for the first thread,
            //  the second thread should start placing results in position 3 (starting from 0)
            int ThreadResultStartingPosition = 0;
            for (int n = 0; n < j; ++n) {
                ThreadResultStartingPosition += ResultsPerThread[n];
            }

            //  Lets divide the samples to process between the
            // threads as we did in the previos loop
            int thread_sample_start_index = j * samples_per_thread + 2;
            // TODO (Eli): This still might suffer from the same issue as in the previos loop
            int thread_sample_end_index = thread_sample_start_index + samples_per_thread;

            for (int k = thread_sample_start_index; k < thread_sample_end_index; k++) {

                //  Each thread will check if in his relative position the red_cost
                //  condition was met
                if(1 == ConditionResultsArray[k]) {
                    //  We increase by 1 first since in the original code they started placing results
                    //  after the next++ (line 150 in the original code)
                    CurrentThreadResultPositionTracker++;
                    //  Result position is the starting position +
                    int ResultPosition = ThreadResultStartingPosition + CurrentThreadResultPositionTracker;
//                    if(j==1)
//                        printf("ResultPosition: %d, CurrentThreadResultPositionTracker: %d \n",ResultPosition,CurrentThreadResultPositionTracker);
                    //  Write the results
                    perm[ResultPosition]->a = resArr[k].arc;
                    perm[ResultPosition]->cost = resArr[k].red_cost;
                    perm[ResultPosition]->abs_cost = ABS(resArr[k].red_cost);
                }
            }

            //  Now set the basket size to be the sum of all results
            for (int l = 0; l < NUM_OF_THREADS; ++l) {
                basket_size = ResultsPerThread[l];
            }
        }
    }


NEXT:
    /* price next group */
    arc = arcs + group_pos;
    for( ; arc < stop_arcs; arc += nr_group )
    {
        if( arc->ident > BASIC )
        {
            /* red_cost = bea_compute_red_cost( arc ); */
            red_cost = arc->cost - arc->tail->potential + arc->head->potential;
            if( bea_is_dual_infeasible( arc, red_cost ) )
            {
                basket_size++;
                perm[basket_size]->a = arc;
                perm[basket_size]->cost = red_cost;
                perm[basket_size]->abs_cost = ABS(red_cost);
            }
        }

    }

    if( ++group_pos == nr_group )
        group_pos = 0;

    if( basket_size < B && group_pos != old_group_pos )
        goto NEXT;

    if( basket_size == 0 )
    {
        initialize = 1;
        *red_cost_of_bea = 0;
        return NULL;
    }

    sort_basket( 1, basket_size );

    *red_cost_of_bea = perm[1]->cost;
    return( perm[1]->a );
}
#endif

#ifdef OLD_CODE
/**************************************************************************
PBEAMPP.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.
Copyright (c) 2000-2002 ZIB & Loebel.
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:04 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pbeampp.c,v 1.10 2005/02/17 19:42:32 bzfloebe Exp $  */



#define K 300
#define B  50



#include "pbeampp.h"




#ifdef _PROTO_
int bea_is_dual_infeasible( arc_t *arc, cost_t red_cost )
#else
int bea_is_dual_infeasible( arc, red_cost )
    arc_t *arc;
    cost_t red_cost;
#endif
{
    return(    (red_cost < 0 && arc->ident == AT_LOWER)
            || (red_cost > 0 && arc->ident == AT_UPPER) );
}







typedef struct basket
{
    arc_t *a;
    cost_t cost;
    cost_t abs_cost;
} BASKET;

static long basket_size;
static BASKET basket[B+K+1];
static BASKET *perm[B+K+1];



#ifdef _PROTO_
void sort_basket( long min, long max )
#else
void sort_basket( min, max )
    long min, max;
#endif
{
    long l, r;
    cost_t cut;
    BASKET *xchange;

    l = min; r = max;

    cut = perm[ (long)( (l+r) / 2 ) ]->abs_cost;

    do
    {
        while( perm[l]->abs_cost > cut )
            l++;
        while( cut > perm[r]->abs_cost )
            r--;

        if( l < r )
        {
            xchange = perm[l];
            perm[l] = perm[r];
            perm[r] = xchange;
        }
        if( l <= r )
        {
            l++; r--;
        }

    }
    while( l <= r );

    if( min < r )
        sort_basket( min, r );
    if( l < max && l <= B )
        sort_basket( l, max );
}






static long nr_group;
static long group_pos;


static long initialize = 1;


#ifdef _PROTO_
arc_t *primal_bea_mpp( long m,  arc_t *arcs, arc_t *stop_arcs,
                              cost_t *red_cost_of_bea )
#else
arc_t *primal_bea_mpp( m, arcs, stop_arcs, red_cost_of_bea )
    long m;
    arc_t *arcs;
    arc_t *stop_arcs;
    cost_t *red_cost_of_bea;
#endif
{
    long i, next, old_group_pos;
    arc_t *arc;
    cost_t red_cost;

    if( initialize )
    {
        for( i=1; i < K+B+1; i++ )
            perm[i] = &(basket[i]);
        nr_group = ( (m-1) / K ) + 1;
        group_pos = 0;
        basket_size = 0;
        initialize = 0;
    }
    else
    {
        for( i = 2, next = 0; i <= B && i <= basket_size; i++ )
        {
            arc = perm[i]->a;
            red_cost = arc->cost - arc->tail->potential + arc->head->potential;
            if( (red_cost < 0 && arc->ident == AT_LOWER)
                || (red_cost > 0 && arc->ident == AT_UPPER) )
            {
                next++;
                perm[next]->a = arc;
                perm[next]->cost = red_cost;
                perm[next]->abs_cost = ABS(red_cost);
            }
                }
        basket_size = next;
        }

    old_group_pos = group_pos;

NEXT:
    /* price next group */
    arc = arcs + group_pos;
    for( ; arc < stop_arcs; arc += nr_group )
    {
        if( arc->ident > BASIC )
        {
            /* red_cost = bea_compute_red_cost( arc ); */
            red_cost = arc->cost - arc->tail->potential + arc->head->potential;
            if( bea_is_dual_infeasible( arc, red_cost ) )
            {
                basket_size++;
                perm[basket_size]->a = arc;
                perm[basket_size]->cost = red_cost;
                perm[basket_size]->abs_cost = ABS(red_cost);
            }
        }

    }

    if( ++group_pos == nr_group )
        group_pos = 0;

    if( basket_size < B && group_pos != old_group_pos )
        goto NEXT;

    if( basket_size == 0 )
    {
        initialize = 1;
        *red_cost_of_bea = 0;
        return NULL;
    }

    sort_basket( 1, basket_size );

    *red_cost_of_bea = perm[1]->cost;
    return( perm[1]->a );
}

#endif







