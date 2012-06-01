#ifndef __HIST__
#define __HIST__

#include"file_scr.h"

extern table histTab;
extern void editHist(void);
extern void showHist(int X, int Y);
extern int wrt_hist(FILE *nchan);
extern int rdr_hist(FILE *nchan);
extern int wrt_hist2(FILE *nchan, char * comment);
extern int rdr_hist2(FILE *nchan,char * comment);
extern int clearHists(void);
extern void fillHists(double w);
extern int correctHistList(void);
extern int add_hist(FILE *f, char *procname);

#endif
