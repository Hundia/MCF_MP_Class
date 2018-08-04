NEXT:
    /* price next group */
    #define P 8
    int i,s,myar[P];    //shared arry to hold partial sums
    
    for(i=0;i<P;i++)
        myar[i]=0;
    sarc=arc = arcs + group_pos;
    s = (stop_arcs-arc)/P;
 #pragma omp parallel for
 for(j=arc ; j < stop_arcs; j += s)
 { int arct=arc;
    for(arct=j ; arct < j+s; arct += nr_group )
    {
        if( arct->ident > BASIC )
        {
            /* red_cost = bea_compute_red_cost( arc ); */
            red_cost = arct->cost - arct->tail->potential + arct->head->potential;
            if( bea_is_dual_infeasible( arct, red_cost ) )
            { int count;
                isone[count]=1;
                isarc[count]=arct;
                isred_cost[count]=red_cost;
                count++;
            }
        }
    }
    myar[getthreadid()]=count;
  }
    for(i=1;i<P;i++)myar[i]+=myar[i-1];
 #pragma omp parallel for
 for(j=arc ; j < stop_arcs; j += s)
 {
    int arct=arc;
    int basket_size=myar[getthreadid()];
    for(arct=j ; arct < j+s; arct += nr_group )
    {
        int count=0;
        perm[basket_size]->a = isarc[count];
        perm[basket_size]->cost = red_cost;
        perm[basket_size]->abs_cost = ABS(red_cost);
        basket_size++;
     }
        
 }
