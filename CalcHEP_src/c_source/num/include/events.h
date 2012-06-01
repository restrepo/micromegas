#ifndef __EVENTS__
#define __EVENTS__
#include"vegas.h"

extern void generateEvents(vegasGrid * vegPtr, double (*func)(double*,double), char* fname, FILE* iprt);

int saveEventSettings(FILE * f);
int readEventSettings(FILE * f);

#endif
