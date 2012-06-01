#ifndef __VEGAS__
#define __VEGAS__


typedef  struct vegasGrid 
{
   
int ndim,     /* number of dimensions */ 
    ndmx;   
    long nCubs;
    double * x_grid;
    double * c_grid;
    float  * fMax;
} vegasGrid;
        
extern vegasGrid *  vegas_init
(int dim,  /* number of dimensions */
 int ndmx   /* size of grid */
);

extern void vegas_finish( vegasGrid * vegPtr);

extern int vegas_int(vegasGrid * vegPtr, 
 long ncall0,                       /* number of integrand calls */
 double alph,                       /* rate of grid improvement  */
 double(*fxn)(double *,double),     /* integrand */
 double *ti,                        /* integral estimation */ 
 double *tsi                        /* standard deviation */
);


extern int vegas_max(
vegasGrid * vegPtr, 
long  nCubs, 
long nRandom,
long nSimplex,
double (*fxn)( double *,double), 
double * eff
);


extern int vegas_events(
vegasGrid * vegPtr, 
long  nEvents,
double gmax,
double milk,
int nSimplex, 
double (*fxn)( double *,double), 
void (*out)(long ,int,char*),
double * eff,  /* efficiency */
double * rmax, /* max reached */
double * mult, /* partion of multiple events */
double * neg   /* partion of events with negative weght */
);

#endif
