/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "f_c.h"
#include "math.h"


double pow_dl(double ap, long bp)
{
double pow, x;
long n;
unsigned long u;

pow = 1;
x = ap;
n = bp;

if(n != 0)
        {
        if(n < 0)
                {
                n = -n;
                x = 1/x;
                }
        for(u = n; ; )
                {
                if(u & 01)
                        pow *= x;
                if(u >>= 1)
                        x *= x;
                else
                        break;
                }
        }
return(pow);
}
