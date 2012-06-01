#ifndef __READ_MDL_
#define __READ_MDL_

#include"file_scr.h"
extern table modelTab[5];

#define vars_tab modelTab[0]
#define func_tab modelTab[1]
#define prtcls_tab modelTab[2]
#define lgrng_tab modelTab[3]


extern void  readModelFiles(char * path, int l);
extern int  loadModel(int check, int ugForce);
extern void  readEXTLIB(void);
extern int read2VarsParticles(void);
#endif
