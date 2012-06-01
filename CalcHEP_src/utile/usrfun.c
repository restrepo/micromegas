/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include<stdlib.h>
#include<stdio.h>

/*                FOR USER 

The given file contains a dummy version of 
       double usrfun(char * name)
function which should be  replaced by other ones written by the user.
We expect that user version will be attached to calchep numerical code
as  via 'Libraries' model file, while the current dummy version will 
be kept unchanged here. 
  If one uses U<txt> for  n_calchep  cuts and histograms then usrfun(<txt>) 
will be called. Say 'Uabs'  corresponds to usrfun("abs").
 
Below we present tools which can be used for usrfun.       
*/

extern double calcPhysVal(char key,char * lv);
/*   This function allow to reproduce built-in CalcHEP functions.
Here 'key' is one-character  function identifier ( see $CALCHEP/help/n_cut.txt)   
lv presents particle numbers terminates by zero. For example 
     calcPhysVal('M',"\3\4");  
return joint mass of particle 3 and 4 in reaction list. 
   One can extend list of observables using
*/

extern int nin_int, nout_int; /* numbers of incoming and outgoing particles */
 
extern double pvect[400];     /* momenta of particles 

             q[k]=pvect[4*(I-1)+k]  k=0,1,2,3 - momenta of I^{th} particle; 
                  1<=I<=nin_int          - incoming particles;
             nin_int<I<=nin_int+nout_int - outgoing particles;
             Energy of all particle are positive  pvect[4*(I-1)]>0;
             Axis of collision k=3.    
                               */ 

extern char * (*pinf_int)(int nsub, int nprtcl, double * pmass, long*pnum);
                               /* Input parameters are 
                      nsub - current subprocess number,*/
extern int Nsub;               /* returns number of current subprocess
                      nprtcl - number of particle in reaction.

                                 Return value 
                      name of particle, reference of static object.           

                                 Outgoing parameters are
                      pmass - particle mass;
                      pnum - PDG code. One can substitute NULL to ignore
                             this information.
                               */
                                                                
double usrfun(char * name)
{   
   fprintf(stdout," usrfun(char* name)  called with parameter %s\n"
                  " But is not  defined!\n",name);
   sortie(54);
   return 0.;
}
