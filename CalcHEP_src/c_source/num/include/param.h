#ifndef __PARAM__
#define __PARAM__

extern int checkParam(void);

extern int selectParam(int x,int y,  char*mess, void ** pscrPrt,
   int for22, int vars_on,  int func_on, double ** varPos, char * varName,int*nPos);

extern int change_parameter(int x,int y, int for22);

extern void show_depend(int x, int y);

extern double Pcm22;

#endif
