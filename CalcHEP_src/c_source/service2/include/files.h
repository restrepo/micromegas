#ifndef __FILES_
#define __FILES_
#include"syst.h"
    extern char * outputDir;
    extern char pathtocomphep[STRSIZ];
    extern char  version[36];
    extern void  copyfile(char *   namefrom, char *  nameto);
    extern void nextFileName(char* f_name,char * firstname,char * ext);
#define f_slash '/'
#endif
