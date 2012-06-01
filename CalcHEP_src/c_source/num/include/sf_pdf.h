#ifndef __SF_PDF__
#define __SF_PDF__

extern int     p_pdf(long pNum);
extern void    n_pdf(int i, char *name);
extern int     r_pdf(int i, char *name);
extern int     m_pdf(int i);
extern int     mc_pdf(int i);
extern int     be_pdf(int i,double * be, double * mass);
extern double  c_pdf(int i, double x,double q);

#endif
