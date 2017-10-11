/*
 Copyright (C) 1997, Victor Edneral
*/
/* Modified by S.A. May 11 1991     */
#include <math.h>
#include <unistd.h>
#include "syst.h"
#include "crt.h"
#include "crt_util.h"
#include "tex_util.h"
#include "files.h"
#include "plot.h"

 static int X1, Y1, X2, Y2;

 static  struct{double  xmin, xmax, ymin, ymax;} grafminmax;
 static  double xscale, yscale;
    

static   int pictureX=300;
static   int pictureY=200;
static   char  letterSize[14]="normalsize";

static int  logScale=1;
static int  fgcolor = FGmain,  bkcolor = BGmain; 


                           
static int nlog10(double x)
{ double lg=log10(x);
  if (lg<0) return lg-1; else return lg;
}


static void axisDesign(double xmin,double xmax, int islog,
              double * xfirst, double * step,int * nsub)
{
  double dx,dx0,dx100;
  int n,n0,n1;
  char xmintxt[50];

  if (islog) 
  { 
    n=1+log10(xmax/xmin)/10;
    *step=pow((double)10,(double)n);     
    *xfirst=pow((double)10,(double)nlog10(xmin));
    *nsub=9; 
  }else
  {      
    dx=xmax-xmin;
    n=nlog10(dx); 
    dx0=pow((double)10,(double)n);  

    dx100= 10*dx/dx0;   
    if (dx100<15.0) { *step= 0.3*dx0;  *nsub=3; } else
    if (dx100<24.0) { *step= 0.5*dx0;  *nsub=5; } else 
    if (dx100<30.0) { *step= 0.8*dx0;  *nsub=8; } else 
    if (dx100<45.0) { *step=   1*dx0;  *nsub=10;} else 
    if (dx100<90.0) { *step=   2*dx0;  *nsub=2; } else 
                    { *step=   3*dx0;  *nsub=3; }
    if( fabs(xmin)<=(*step)*10 ) *xfirst=0;else
    {                  
       n0=nlog10(*step);     
       n1=nlog10(fabs(xmin)); 
    
       sprintf(xmintxt,"%.*E",MAX(0,n1-n0-1),xmin);
       trim(xmintxt);
       sscanf(xmintxt,"%lf",xfirst);         
    }       
    while(*xfirst>xmin+(*step)/(*nsub)  ) *xfirst -= *step;
    while(*xfirst+(*step) < xmin) *xfirst += *step;
  }
}

static long  chpround(double  x)
{  if(x>0) x+=0.5; else x-=0.5;
   return (long)x;
}

static void  gminmax(double *f ,double * df,int dim, double *ymin,double*ymax)
{
   int i;
   if (dim==0) return;

   *ymin=f[0];
   *ymax=f[0];
   if(df) { *ymax += df[0]; *ymin -= df[0];}
 
   if(df) for(i=1;i<dim;i++)
   { *ymin=MIN(*ymin, f[i]-df[i]); 
     *ymax=MAX(*ymax, f[i]+df[i]);
   } else for(i=1;i<dim;i++)
   { *ymin=MIN(*ymin, f[i]); 
     *ymax=MAX(*ymax, f[i]);
   } 
}


static int scx(double  x)
{ 
   x -= grafminmax.xmin; 
   return X1 + chpround(xscale * x);
} 


static int  scy(double  x)
{  
   if (logScale)   x = log10(x / grafminmax.ymin); 
   else            x -= grafminmax.ymin; 
   return Y1 - chpround(yscale * x);
} 


static double  dscx(double x)
{ return X1 + (xscale * (x-grafminmax.xmin));} 


static double  dscy(double  x)
{ 
   if (logScale) 
   {  
     if(x>0)  x = log10(x / grafminmax.ymin);
     else     x=   -10000*log10(grafminmax.ymax/grafminmax.ymin) ;
   }
   else            x -= grafminmax.ymin; 
   return Y1 - (yscale * x);
} 

static double xPhys(void)
{
  return       (mouse_info.x -X1)/xscale + grafminmax.xmin;
}

static double yPhys(void)
{  if(logScale)
   return pow((Y1-mouse_info.y)/yscale,10.)*grafminmax.ymin;  
   else   return -(mouse_info.y -Y1)/yscale + grafminmax.ymin;
}


static void doubleToStr(double x,double step,char * s)
{ 
   int    n1,n0;
   char * emark;
   char s1[200],s2[200];
   
   if(fabs(x) < 1.E-2 * step) strcpy(s,"0.0"); else
   {  int d;
      n1=nlog10(fabs(x));
      n0=nlog10(step);

      sprintf(s1,"%.*E",MAX(0,n1-n0),x);
      emark=strstr(s1,"E");
      emark[0]=0;
      sscanf(emark+1,"%d",&d);
      sprintf(s2,"%sE%d",s1,d);
      sprintf(s1,"%.*lf",MAX(0,-n0),x);
      if (strlen(s1)<=strlen(s2)) strcpy(s,s1);else strcpy(s,s2);
   }   
}    

static int colList[]={Black,Blue, Green,Cyan,Red, Magenta,Brown,LightGray,
                 DarkGray, LightBlue, LightGreen,LightCyan,
                 LightRed, LightMagenta, Yellow, White}; 
static char*colListTxt[]={"Black","Blue", "Green","Cyan","Red", "Magenta","Brown","LightGray",
                 "DarkGray", "LightBlue", "LightGreen","LightCyan",
                                  "LightRed", "LightMagenta", "Yellow", "White"};


static void  gaxes(char* upstr, char* xstr, int N, char**Y)
{double  xmax,xmin,ymax,ymin,aa,step; 
 int     m,naxis,islog,n_sub,i;

 int th = tg_textheight("0");
 int tw = tg_textwidth("0");
 int  hash = (th+tw)/2; 
 int  texhash=(th*texyscale+tw*texxscale+1)/2;
   xmin = grafminmax.xmin; 
   xmax = grafminmax.xmax; 
   ymin = grafminmax.ymin; 
   ymax = grafminmax.ymax;
   if ( 1.E-4 * (fabs(ymax) + fabs(ymin)) >= fabs(ymax - ymin)) 
   {  if (ymin == 0.0) {  ymin = -0.5; ymax = 0.5; } 
        else  if (ymin < 0.0) 
              {  
                 ymin *= 2.0;
                 if (ymax < 0.0) ymax *= 0.5;  else ymax *= 2.0; 
              }  else { ymax *= 2.0; ymin *= 0.5;} 
   }           
   grafminmax.ymin = ymin; 
   grafminmax.ymax = ymax; 
  
   if (logScale && ( ymin <= 0.0 || log10(ymax/ymin)<1 )) logScale=0;

   tg_setlinestyle(SolidLn,NormWidth);

   X1 = 10*tw+hash; X2 = tg_getmaxx() - 2*tw;  
   Y1 = tg_getmaxy() - 5*th; Y2 = 5*th/2; 

   xscale = (X2 - X1) / (xmax - xmin);
   yscale = (Y1 - Y2) / (logScale ? log10(ymax / ymin) : ymax - ymin); 
   tg_settextjustify(CenterText,TopText);
   tg_outtextxy((X1+X2)/2,CenterText,upstr);  
    
   for(naxis=0;naxis<2;naxis++)
   { double zmin,zmax,zmin1,zmax1;
     char xy;
     int xend,yend;
     int len;

     if(naxis==0)
     {  
        xy='X';
        if(texflag){ hash=texhash/texyscale;texhash= -texhash;}      
        islog = 0; 
        tg_settextjustify(CenterText,TopText);
        zmin=xmin;
        zmax=xmax;
        xend=X2;
        yend=Y1;              
     } else
     { 
        xy='Y';
        if(texflag){texhash= -texhash; hash=texhash/texxscale;}      
        islog = logScale; 
        tg_settextjustify(RightText,CenterText);
        zmin=ymin;
        zmax=ymax;
        xend=X1;
        yend=Y2;              
     }
     
     if(texflag) f_printf(out_tex,"%% ====================   %c-axis =============\n",xy );
     
     zmax1 = zmax + fabs(zmax)*1.E-7;
     zmin1 = zmin - fabs(zmin)*1.E-7;
     axisDesign(zmin,zmax, islog,&aa, &step,&n_sub);
     if(texflag)
     { double Nd,offset;
       char axis[10];
       char d[10];
       
       if(islog) {  strcpy(axis,"LogAxis"); 
                    Nd=log10(zmax/zmin)/log10(step);
                    offset=10*zmin/(step*aa); 
                    strcpy(d,"");
                 }                              
        else     {  strcpy(axis,"LinAxis");
                    Nd=(zmax-zmin)/step;       
                    offset=-n_sub*(aa-zmin)/step; 
                    sprintf(d,"%d,",n_sub); 
                    if (fabs(offset) >fabs(offset-n_sub) ) offset -=n_sub;                             
                 }
                 
      f_printf(out_tex,"\\%s(%.2f,%.2f)(%.2f,%.2f)(%.3f,%s%d,%.3f,1.5)\n",
        axis,texX(X1),texY(Y1),texX(xend),texY(yend),Nd,d,(2*texhash)/3,offset);                 
     }
     else tg_line(X1,Y1,xend,yend);
     len=0;
     while (aa <= zmax1) 
     {  double da;
        char snum[30];      
        int i;     
        if(aa >= zmin1)
        {  
           if(naxis==0){m=scx(aa);if(!texflag)tg_line(m,Y1,m,Y1+hash);}
           else        {m=scy(aa);if(!texflag)tg_line(X1-hash,m,X1,m);}
           if (islog)  doubleToStr(aa,aa,snum); else doubleToStr(aa,step,snum);
           len=MAX(len,tg_textwidth(snum));
           if(naxis==0) tg_outtextxy(m,Y1+hash,snum);
           else         tg_outtextxy(X1-hash,m,snum);     
        }
        if (islog) da = aa*step -aa ;else da=step;
        da=da/n_sub;
        for(i=1;i<n_sub;i++)
        {  aa +=da;    
           if(!texflag && aa >=zmin1 && aa<=zmax1)
           {
             if (naxis==0) { m = scx(aa);tg_line(m,Y1,m,Y1 + 2*hash/3); }
             else          { m = scy(aa);tg_line(X1 - 2*hash/3, m,X1,m);}
           }    
        }
        aa += da;                                                                                   
     }
     
     if(naxis==0) 
     {  tg_settextjustify(RightText,TopText); 
        tg_outtextxy(X2,Y1 + hash + th ,xstr);
     }else
/*     if(texflag)
     {
         
        f_printf(out_tex,"\\rText(%.1lf,%.1lf)[tr][l]{%s}\n",
        texX(X1-len-3*hash/2) ,texY(Y2),Y[0]); 
     
     }else
*/          
     {  tg_settextjustify(LeftText,CenterText);
        goto_xy(10,2);
        for(i=0;i<N;i++) 
        {  if(texflag) 
           { fprintf(out_tex,"\\SetColor{%s}\n",colListTxt[i]);
             f_printf(out_tex,"\\Text(%.1lf,%.1lf)[rb]{%s}\n",
             texX(X1) ,texY(Y2-i*hash),Y[i]);
             fprintf(out_tex,"\\Line(%.1lf,%.1lf)(%.1lf,%.1lf)\n",
             texX(X1+hash ),
             texY(Y2-i*hash),   texX(X1 +6*hash), texY(Y2-i*hash));
              
           }  
            else  
          {  scrcolor(colList[i],bkcolor);
             print("%s; ",Y[i]);
          } 
        }      
//        tg_outtextxy(X1,2*th ,ystr);  
     }
   }
   if(texflag) f_printf(out_tex,"%% ============== end of axis ============\n");  
     
}  /* GAxes */ 


static void plot_curve(double xMin,double xMax,int dim, double *f)
{  double   x,y,xx,yy;
   int i;
   double step=(xMax-xMin)/dim;
   double ymax = dscy(grafminmax.ymin);
   double ymin = dscy(grafminmax.ymax);
          
   for(i=1;i<dim;i++)
   {
      x =dscx(xMin+(i-0.5)*step);
      xx=dscx(xMin+(i+0.5)*step); 
      y =dscy(f[i-1]);
      yy=dscy(f[i]);
      
      if ( yy < y) 
      { double z;
        z=yy;yy=y;y=z;
        z=xx;xx=x;x=z;
      }
   
      if (yy>ymin &&  y<ymax)         
      {
        if(yy>ymax){ xx= xx-((xx-x)*(yy-ymax))/(yy-y); yy=ymax;}
        if(y<ymin) {  x=  x-((x-xx)*(y -ymin))/(y-yy); y=ymin;}
        tg_line((int)x,(int)y,(int)xx,(int)yy);
      }
   }
      
}


static void plot_hist(double xMin,double xMax,int dim, double *f,double *df)
{  double   x,y,yy;
   int i;

   double ymax = dscy(grafminmax.ymin);
   double ymin = dscy(grafminmax.ymax);
   double step=(xMax-xMin)/dim;

   for(i=0;i<dim;i++)
   {  
     y =dscy(f[i]);
     if(y<ymax && y>ymin) tg_line((int)dscx(xMin+i*step),  (int)y,
                                    (int)dscx(xMin+(i+1)*step),(int)y);
     y =MIN(dscy(f[i]-df[i]),ymax);
     yy=MAX(dscy(f[i]+df[i]),ymin);
     x=(dscx(xMin+i*step)+dscx(xMin+(i+1)*step))/2;
     if(y>yy) tg_line((int)x,(int)y, (int)x,(int)yy );        
   }
}



static void  writetable1(double xMin,double xMax,int dim,char*x_str,char*upstr, 
 int N, double**f,double**df,char**Y)  
{  char       filename[100], buff[STRSIZ],command[100];
   FILE *     outfile;
   int        i,k;
   int ID;
   
   double dx=xMax-xMin;
   double x1=xMin+0.35*dx, x2=xMin+0.43*dx, x3=xMin+0.45*dx;
   double ymin=grafminmax.ymin, ymax=grafminmax.ymax, dy=ymax-ymin, yl;
   int ilab=0;
   double step=(xMax-xMin)/dim;
   nextFileName(filename,"plot_",".txt");
   for(i=strlen(filename); filename[i]!='_';i--) continue;
   
   sscanf(filename+i+1,"%d",&ID);  
   outfile=fopen(filename,"w");
   fprintf(outfile,"#title %s\n",upstr);
   fprintf(outfile,"#yName ");
   for(k=0;k<N;k++)
   {  fprintf(outfile,"%s",Y[k]);
      if(df[k])  fprintf(outfile,"{h}"); else fprintf(outfile,"{c} ");
   } 
   fprintf(outfile,"\n");
   fprintf(outfile,"#xName %s\n",x_str);
   fprintf(outfile,"#xMin %E\n",xMin);
   fprintf(outfile,"#xMax %E\n",xMax);
   fprintf(outfile,"#xDim %d\n",dim);

   fprintf(outfile,"#--- GNUPLOT section ---\n");
//   fprintf(outfile,"#GNUPLOT set key off\n");
   fprintf(outfile,"#GNUPLOT set title  '%s'\n",upstr);
   fprintf(outfile,"#GNUPLOT set xlabel '%s'\n",x_str);
   fprintf(outfile,"#GNUPLOT plot[%G:%G] ",xMin,xMax);
   
   for(i=1,k=0;k<N;k++,i++)
   { 
   if(df[k]==NULL)fprintf(outfile,
       "'%s' using (%G +$0*%G):%d w l ti '%s'",
                 filename,xMin+0.5*step,step,i,Y[k]);
   else {        fprintf(outfile,
       "'%s' using (%G +$0*%G):%d:%d  w error ti '%s'",
                 filename,xMin,step,i,i+1,Y[k]);
         i++;
        }         
   if(k<N-1) fprintf(outfile,", ");else fprintf(outfile,"\n");               
   }      
         
   fprintf(outfile,"#--- PAW section ---\n");
   fprintf(outfile,"#PAW  TITLE '%s'\n",upstr);
   fprintf(outfile,"#PAW  vector/Create XX(%d)  \n",dim);
   fprintf(outfile,"#PAW  sigma XX=ARRAY(%d,%G#%G)\n",dim,xMin+0.5*step,xMax-0.5*step);
   sprintf(buff,"");

   for(k=0;k<N;k++)
   {  
     fprintf(outfile,"#PAW  vector/Create Y%d(%d)\n",k+1,dim);
     if(df[k])fprintf(outfile,"#PAW  vector/Create dY%d(%d)\n",k+1,dim);  
   }

   fprintf(outfile,"#PAW vector/Read ");
   for(k=0;k<N;k++)
   { if(k) fprintf(outfile,",");
     fprintf(outfile,"Y%d",k+1);   
     if(df[k]) fprintf(outfile,",dY%d",k+1);  
   }   
   fprintf(outfile," '%s' ' ' 'OC' '-/#/' \n",filename); 
   fprintf(outfile,"#PAW 1d 10 ! 1 %G %G;  min 10 %G; max 10 %G\n",
              xMin,xMax,grafminmax.ymin*0.9,grafminmax.ymax*1.1);
   fprintf(outfile,"#PAW his/plo 10\n"); 
   
   fprintf(outfile,"#PAW set chhe 0.3; set lwidt 5; set hwidt 5 10\n"); 
  
   for(k=0;k<N;k++)
     {  
     fprintf(outfile,"#PAW * %s\n",Y[k]);
     fprintf(outfile,"#PAW  set dmod %d\n",k+1);   

     ilab=ilab+1;
     fprintf(outfile,"#PAW  set plci %i ; set txci %i  \n",ilab,ilab);
     
     yl=ymax*1.1-0.05*ilab*dy;
     
     fprintf(outfile,"#PAW  dline %G %G %G %G ; itx %G %G '%s' \n",
     x1,x2, yl, yl, x3, yl, Y[k]);
 
     if(df[k]==NULL)
     { 
        fprintf(outfile,"#PAW  GRAPH %d XX Y%d l \n",dim,k+1);
     }else 
     {
       fprintf(outfile,"#PAW   histogram/create/1dhisto %d '%s' %d %G %G\n",
                         k+1,upstr,dim,xMin,xMax);                
       fprintf(outfile,"#PAW  histogram/put/contents %d Y%d\n",k+1,k+1);
       fprintf(outfile,"#PAW  histogram/put/errors   %d dY%d\n",k+1,k+1);   
       fprintf(outfile,"#PAW  histogram/plot %d  s\n",k+1);


//       fprintf(outfile,"#PAW  atitle  '%s' '%s' \n",x_str,y_str);
     }
   }       
       fprintf(outfile,"#PAW  set txci 1 \n" );
       fprintf(outfile,"#PAW  atitle  '%s'\n",x_str);
        
       
   fprintf(outfile,"#---   starting of data ---\n#");
   for(k=0;k<N;k++) 
   {  fprintf(outfile,"  %-12.12s",Y[k]);
      if(df[k])
      { char txt[100];
        sprintf(txt,"d(%s)",Y[k]);
        fprintf(outfile,"  %-12.12s",txt);
      }
   }
   
   fprintf(outfile,"\n");        
   for(i=0;i<dim;i++)
   { for(k=0;k<N;k++)  
     { fprintf(outfile,"%-12E  ",f[k][i]);
       if(df[k]) fprintf(outfile,"  %-12E",df[k][i]);
     } 
     fprintf(outfile,"\n");
   }
   fclose(outfile);
   sprintf(buff," You can find results in the file\n%s",filename);
   messanykey(10,12,buff);

//printf("    ID=%d\n", ID);
   
   sprintf(command,"grep \\#PAW plot_%d.txt |sed s/#PAW// > plot_%d.kumac",ID,ID);
   system(command);
   sprintf(command,"grep \\#GNUPLOT plot_%d.txt |sed s/#GNUPLOT// > plot_%d.gnu",ID,ID);
   system(command);
}

static void  writetable2(double xMin, double xMax, int dimX, 
                         double yMin,double yMax, int dimY, 
double * f, double *df, char*  upstr,  char*  x_str, char*  y_str)
{  char       filename[STRSIZ], buff[STRSIZ];
   FILE *     outfile;
   int        i;
   int ID;
   
   nextFileName(filename,"plot_",".txt");
   sscanf(filename,"plot_%d",&ID);
   
   outfile=fopen(filename,"w");
   fprintf(outfile,"#type 2  %%2d-plot\n");
   fprintf(outfile,"#title %s\n",upstr);
   fprintf(outfile,"#xName %s\n",x_str);
   fprintf(outfile,"#xMin %E\n",xMin);
   fprintf(outfile,"#xMax %E\n",xMax);
   fprintf(outfile,"#xDim %d\n",dimX);
   fprintf(outfile,"#yName %s\n",y_str);
   fprintf(outfile,"#yMin %E\n",yMin);
   fprintf(outfile,"#yMax %E\n",yMax);
   fprintf(outfile,"#yDim %d\n",dimY);

   fprintf(outfile,"#--- PAW section ---\n");
   fprintf(outfile,"#PAW   histogram/create/2dhisto %d '%s' %d %G %G %d %G %G\n",
                         ID,upstr,dimY,yMin,yMax, dimX,xMin,xMax);
    fprintf(outfile,"#PAW  vector/Create  Y%d(%d,%d)\n",ID,dimY,dimX);
    fprintf(outfile,"#PAW  vector/Create dY%d(%d,%d)\n",ID,dimY,dimX);
    fprintf(outfile,"#PAW  vector/Read Y%d,dY%d '%s' ' ' 'OC' '-/#/' \n",ID,ID,filename);                      
    fprintf(outfile,"#PAW  histogram/put/contents %d Y%d\n",ID,ID);
    fprintf(outfile,"#PAW  histogram/put/errors   %d dY%d\n",ID,ID);   
    fprintf(outfile,"#PAW  histogram/2d/lego %d\n",ID);
    fprintf(outfile,"#PAW  atitle  '%s' '%s' \n",y_str,x_str);
                                                                                
     
   fprintf(outfile,"#---   starting of data ---\n");
   fprintf(outfile,"#--- F ---    --- dF ---\n");
   for(i=0;i<dimX*dimY;i++) fprintf(outfile,"%-12E  %-12E\n",f[i],df[i]); 
 
   fclose(outfile);
   sprintf(buff," You can find results in the file\n%s",filename);
   messanykey(10,12,buff);
}

static void  writeMath1(double xMin, double xMax, int dim, double * f, 
double *df, char*  upstr,  char*  x_str, char*  y_str)
{  char       filename[STRSIZ], buff[STRSIZ];
   FILE *     outfile;
   int        i;
   int ID;
   int expon, expon2;
   double mant, mant2;
   nextFileName(filename,"plot_",".math");
   sscanf(filename,"plot_%d",&ID);
   
   outfile=fopen(filename,"w");

   getcwd(buff,1000);

   fprintf(outfile,"(*******************************************************\n");
   fprintf(outfile,"Use the following command to load this data:\n");
   fprintf(outfile,"Get[\"%s%c%s\"]\n",buff,f_slash,filename);
   fprintf(outfile,"(If you move this file, then change the path and/or filename accordingly.)\n");
   fprintf(outfile,"********************************************************)\n\n");
   
   mant=frexp(f[0], &expon);
   fprintf(outfile,"dataCH={\n\t{%-12f,%-12f*2^%-2d}",xMin+(1.0/2.0)*(xMax-xMin)/dim,mant,expon);
   for(i=1;i<dim;i++) {
     mant=frexp(f[i], &expon);
     fprintf(outfile,",\n\t{%-12f,%-12f*2^%-2d}",xMin+(1.0/2.0+i)*(xMax-xMin)/dim,mant,expon);
   }
   fprintf(outfile,"\n};\n\n");

   if(df){
     mant=frexp(df[0], &expon);
     fprintf(outfile,"uncCH={\n\t{%-12f,%-12f*2^%-2d}",xMin+(1.0/2.0)*(xMax-xMin)/dim,mant,expon);
     for(i=1;i<dim;i++) {
       mant=frexp(df[i], &expon);
       fprintf(outfile,",\n\t{%-12f,%-12f*2^%-2d}",xMin+(1.0/2.0+i)*(xMax-xMin)/dim,mant,expon);
     }
     fprintf(outfile,"\n};\n\n");
   }

   mant=frexp(f[0], &expon);
   fprintf(outfile,"histDataCH={\n\t{%-12f,%-12f*2^%-2d},{%-12f,%-12f*2^%-2d}",xMin,mant,expon,xMin+(xMax-xMin)/dim,mant,expon);
   for(i=1;i<dim;i++) {
     mant=frexp(f[i], &expon);
     fprintf(outfile,",\n\t{%-12f,%-12f*2^%-2d},{%-12f,%-12f*2^%-2d}",xMin+(xMax-xMin)/dim*i,mant,expon,xMin+(xMax-xMin)/dim*(i+1.0),mant,expon);
   }
   fprintf(outfile,"\n};\n\n");

   if(df){
     mant=frexp(f[0]-df[0], &expon);
     mant2=frexp(f[0]+df[0], &expon2);
     fprintf(outfile,"histUncCH={\n\t{%-12f,%-12f,%-12f*2^%-2d,%-12f*2^%-2d}",xMin,xMin+(xMax-xMin)/dim,mant,expon,mant2,expon2);
     for(i=1;i<dim;i++) {
       mant=frexp(f[i]-df[i], &expon);
       mant2=frexp(f[i]+df[i],&expon2);
       fprintf(outfile,",\n\t{%-12f,%-12f,%-12f*2^%-2d,%-12f*2^%-2d}",xMin+(xMax-xMin)/dim*i,xMin+(xMax-xMin)/dim*(i+1.0),mant,expon,mant2,expon2);
     }
     fprintf(outfile,"\n};\n\n");
   }

   fprintf(outfile,"(*******************************************************\n");
   fprintf(outfile,"Below is an example of how to plot this data:\n");
   fprintf(outfile,"********************************************************)\n\n");
   
   fprintf(outfile,"histCH=ListPlot[histDataCH,Joined->True, Frame -> True, PlotRange -> Full, Axes->False, FrameLabel -> {\"%s\", \"%s\", \"%s\", \"\"}, PlotStyle->Black];\n\n",x_str,y_str,upstr);

   fprintf(outfile,"ErrBarCH[dt_] := Block[\n");
   fprintf(outfile,"  {w = dt[[2]] - dt[[1]],\n");
   fprintf(outfile,"   c = (dt[[1]] + dt[[2]])/2,\n");
   fprintf(outfile,"   x0 = c - w/4,\n");
   fprintf(outfile,"   x1 = c + w/4,\n");
   fprintf(outfile,"   y0 = dt[[3]],\n");
   fprintf(outfile,"   y1 = dt[[4]]},\n");
   fprintf(outfile,"   Graphics[{Blue, Line[{\n");
   fprintf(outfile,"       {{x0, y0}, {x1, y0}},\n");
   fprintf(outfile,"       {{x0, y1}, {x1, y1}},\n");
   fprintf(outfile,"       {{c, y0}, {c, y1}}\n");
   fprintf(outfile,"   }]}]\n");
   fprintf(outfile,"];\n\n");
   
   fprintf(outfile,"errHistCH = Show[ErrBarCH /@ histUncCH];\n\n");

   fprintf(outfile,"histCombCH=Show[{errHistCH, histCH}, Frame -> True, Axes -> False, FrameLabel -> {\"%s\", \"%s\", \"%s\", \"\"}, AspectRatio -> 2/3]\n\n",x_str,y_str,upstr);
   
     
   fclose(outfile);
   sprintf(buff," You can find the results in the file\n%s",filename);
   messanykey(10,12,buff);
}


static void  writeMath2(double xMin, double xMax, int dimX, 
                         double yMin,double yMax, int dimY, 
double * f, double *df, char*  upstr,  char*  x_str, char*  y_str)
{  char       filename[STRSIZ], buff[STRSIZ];
   FILE *     outfile;
   int        i,j;
   int ID;
   int expon;double mant;
   nextFileName(filename,"plot_",".math");
   sscanf(filename,"plot_%d",&ID);
   
   outfile=fopen(filename,"w");

   getcwd(buff,1000);


   fprintf(outfile,"(*******************************************************\n");
   fprintf(outfile,"Use the following command to load this data:\n");
   fprintf(outfile,"Get[\"%s%c%s\"]\n",buff,f_slash,filename);
   fprintf(outfile,"(If you move this file, then change the path and/or filename accordingly.)\n");
   fprintf(outfile,"********************************************************)\n\n");



   fprintf(outfile,"dataCH={");
   for(i=0;i<dimX;i++) 
     for(j=0;j<dimY;j++){
       mant=frexp(f[i*dimX+j], &expon);
       if(i!=0||j!=0)fprintf(outfile,",");
       fprintf(outfile,"\n\t{%-12f,%-12f,%-12f*2^%-2d}",
	       xMin+(1.0/2.0+i)*(xMax-xMin)/dimX,yMin+(1.0/2.0+j)*(yMax-yMin)/dimY,mant,expon);
   }
   fprintf(outfile,"\n};\n\n");
   
   if(df){
     fprintf(outfile,"uncCH={");
     for(i=0;i<dimX;i++) 
       for(j=0;j<dimY;j++){
	 mant=frexp(df[i*dimX+j], &expon);
	 if(i!=0||j!=0)fprintf(outfile,",");
	 fprintf(outfile,"\n\t{%-12f,%-12f,%-12f*2^%-2d}",
		 xMin+(1.0/2.0+i)*(xMax-xMin)/dimX,yMin+(1.0/2.0+j)*(yMax-yMin)/dimY,mant,expon);
       }
     fprintf(outfile,"\n};\n\n");
   }

   fprintf(outfile,"(*******************************************************\n");
   fprintf(outfile,"Below is an example of how to plot this data:\n");
   fprintf(outfile,"********************************************************)\n\n");
   
   fprintf(outfile,"histCH=ListContourPlot[dataCH, FrameLabel -> {\"%s\", \"%s\", \"%s\", \"\"}]\n\n",x_str,y_str,upstr);

        
   fclose(outfile);
   sprintf(buff," You can find the results in the file\n%s",filename);
   messanykey(10,12,buff);
}


   
void  plot_Nar( char*  title, double xMin, double xMax, char*xName,  int dim, 
               int N, double **f,double**ff,char**Y)
{   
  void *     prtscr;  
  int i,k,nCol0,nRow0,key;
  double      ymin, ymax;
  char        f_name[STRSIZ], menustr[STRSIZ], buff[STRSIZ];

#define MESS "Esc - exit; Mouse - function value; Other keys - menu"
   get_text(1,1,maxCol(),maxRow(),&prtscr);

   grafminmax.xmax=xMax;
   grafminmax.xmin=xMin;

   
   gminmax(f[0],ff[0],dim,&grafminmax.ymin,&grafminmax.ymax);
   for(i=1;i<N;i++)
   { gminmax(f[i],ff[i],dim, &ymin,&ymax);
     if(grafminmax.ymin>ymin) grafminmax.ymin=ymin;
     if(grafminmax.ymax<ymax) grafminmax.ymax=ymax;   
   }
   ymin=grafminmax.ymin;
   ymax=grafminmax.ymax;
   
   logScale=(ymin >0 && grafminmax.ymax/grafminmax.ymin >10);   
   k = 0;      
REDRAW:
   nCol0=maxCol();
   nRow0=maxRow();
   clr_scr(fgcolor,bkcolor);
   
   gaxes(title,xName,N,Y);

   for(i=0;i<N;i++)
   {   fColor=colList[i];
       if(texflag) fprintf(out_tex,"\\SetColor{%s}\n", colListTxt[i]);
       if(ff[i]) plot_hist(xMin,xMax,dim,f[i],ff[i]); 
          else   plot_curve(xMin,xMax,dim,f[i]);
   }                    
   if (texflag)
   {  f_printf(out_tex,"\\end{picture}\n");
      texFinish();
      sprintf(buff,"%s\n%s","LaTeX output is saved in file",f_name);
      messanykey(35,15,buff);
      goto contin;              
   }
   tg_settextjustify(BottomText,LeftText);
   scrcolor(Red,bkcolor);
   tg_outtextxy(0,tg_getmaxy(), MESS);
   for(;;) 
   { double x;
     int i;
     key=inkey();
     if(key!=KB_MOUSE) break;
     if(mouse_info.but1 != 2) continue;
     x=xPhys();

     if(x<=xMin||x>=xMax) continue;
     if(ff[0])
     { 
        i= (x-xMin)/((xMax-xMin)/dim);
        goto_xy(1,maxRow()-1); scrcolor(Red,bkcolor); print("Mouse Info:"); 
        scrcolor(Black,bkcolor);
        print(" X=%.3E  F=%.2E+/-%.1E  ",x,f[0][i],ff[0][i]);
     } else
     { double al;
       i= (x-xMin)/((xMax-xMin)/(dim-1));
       al=(x-xMin)/((xMax-xMin)/(dim-1))-i;
       goto_xy(1,maxRow()-1); scrcolor(Red,bkcolor); print("Mouse Info:"); 
       scrcolor(Black,bkcolor);
       print(" X=%.3E ",x);
       for(k=0;k<N;k++) 
       {  scrcolor(colList[k],bkcolor);
          print("%s",Y[k]);
          scrcolor(fgcolor,bkcolor);
          print("=%.2E ",f[k][i]*(1-al)+f[k][i+1]*al);
       }
     }  
   }
   if( nCol0 != maxCol() && nRow0 != maxRow() ) goto REDRAW; 
   scrcolor(bkcolor,bkcolor);
   tg_outtextxy(0,tg_getmaxy(),MESS);
   if(key==KB_ESC) goto exi;     
contin:   
   do
   {  char sScale[20];
      void * pscr = NULL;

      if(logScale) strcpy(sScale,"Log.   "); else strcpy(sScale,"Lin.   ");  
   
      sprintf(menustr,"%c Y-max = %-9.3G Y-min = %-9.3G Y-scale = %s"
      " Save plot in file"
      " Math file        "
      " LaTeX file       ",
      18,grafminmax.ymax,grafminmax.ymin,sScale);

      menu1(nCol0-20, 2 ,"",menustr,"n_plot_*",&pscr,&k);

      switch (k)
      {   
           case 1:   
             correctDouble(33,11,"Y-max = ",&grafminmax.ymax,1 );
             break;
           case 2:
             correctDouble(33,11,"Y-min = ",&grafminmax.ymin,1 );
             break;
           case 3:  logScale=!logScale; break;
           case 4:  writetable1(xMin,xMax,dim,xName,title,N, f,ff,Y); 
             break;
           case 0:
           case 6:  
             if(grafminmax.ymin >=ymax|| grafminmax.ymax <=ymin ||
               grafminmax.ymin >= grafminmax.ymax)
             { messanykey(10,10," Wrong Y-range");
               break;
             } 
             if(logScale && (grafminmax.ymin <=0 ||
                   grafminmax.ymax/grafminmax.ymin <=10 ) )
             { messanykey(10,10," To use the logarithmic scale,\n"
                                " please, set Ymin and Ymax limits\n"
                                " such that Ymax > 10*Ymin > 0"); 
               logScale = 0; 
             } 
               
             if(k==6)
             { 
               if( !texmenu(&pictureX,&pictureY,letterSize)) break;
               nextFileName(f_name,"plot_",".tex"); 
               texStart(f_name,title,letterSize);
               texPicture(0,0,tg_getmaxx(),tg_getmaxy()-tg_textheight("0"),
                                                  pictureX,pictureY);
               f_printf(out_tex,"\\begin{picture}(%d,%d)(0,0)\n",pictureX,pictureY);  
               del_text(&pscr);
             }      
             goto REDRAW;
           case 5: writeMath1(xMin,xMax,dim,f[0],ff[0],title,xName,Y[0]); 
      }   
      if( nCol0 != maxCol() && nRow0 != maxRow() ) goto REDRAW;
   }  while (1); 
      
exi:
   clr_scr(FGmain,BGmain);
   put_text(&prtscr);
}

void  plot_N(char*title, double xMin, double xMax,  char*  xName,   int dim, int N,...)
{ int         k,i;
                                                               
  double **f; double**ff; char**Y;
  va_list ap;

   if (  1.e-4 * (fabs(xMax) + fabs(xMin)) >= fabs(xMax - xMin))
   {  messanykey(10,10, "Too short interval in X axis !");
      return;
   }   
     
   f= malloc(N*sizeof(double*));
   ff=malloc(N*sizeof(double*));
   Y= malloc(N*sizeof(char*));
   va_start(ap,N);
   
   for(i=0;i<N;i++) 
   { f[i]=va_arg(ap,double*);
     ff[i]=va_arg(ap,double*);
     Y[i]=va_arg(ap,char*);
   }  
   va_end(ap);

   plot_Nar(title,xMin,xMax,xName,dim,N,f,ff,Y);

   free(f);free(ff);free(Y);
}

void   plot_1(double xMin, double xMax, int dim,
                    double *f, double *ff,char* title, char* xstr, char* ystr)
{ plot_Nar(title,xMin, xMax,xstr,dim, 1,&f,&ff, &ystr); } 
                    

void plot_2D(double hMin1,double hMax1,int nBin1,double hMin2,double hMax2,int nBin2,
            double * f,double *df,char *upstr,char* xstr,char * ystr) 
{
  int        k,nCol0,nRow0,key,i,j;
  char       f_name[STRSIZ],  buff[STRSIZ];
  double     fmax,DX,DY;
  void *     prtscr;
  double     P=1;

   logScale = 0;
   get_text(1,1,maxCol(),maxRow(),&prtscr);

   grafminmax.xmax=hMax1;
   grafminmax.xmin=hMin1;
   
   grafminmax.ymax=hMax2;
   grafminmax.ymin=hMin2;
   
   fmax=0;
   for(i=0;i<nBin1;i++) for(j=0;j<nBin2;j++)
   {  double ff;
      ff=fabs(f[i*nBin2+j]+df[i*nBin2+j]);
      if(fmax<ff) fmax=ff;
   } 
   fmax*=1.2; 
   DX=(hMax1-hMin1)/nBin1;
   DY=(hMax2-hMin2)/nBin2;
   
   k = 0;      
REDRAW:
   nCol0=maxCol();
   nRow0=maxRow();
   clr_scr(fgcolor,bkcolor);

   gaxes(upstr,xstr,1,&ystr);
 
   { int i,j;
     for(i=0;i<nBin1;i++) for(j=0;j<nBin2;j++)
     { double x= hMin1+((double)i+0.5)*DX;
       double y= hMin2+((double)j+0.5)*DY;
       double g=pow(f[i*nBin2+j]/fmax,0.5*P);
       double g1=pow((f[i*nBin2+j]+df[i*nBin2+j])/fmax,0.5*P);
       double g2=pow(fabs(f[i*nBin2+j]-df[i*nBin2+j])/fmax,0.5*P);
       if(f[i*nBin2+j]<df[i*nBin2+j]) g2=0;
       if(g && fabs((dscx(x-g*DX/2)-dscx(x+g*DX/2))* 
       (dscy(y+g*DY/2)-dscy(y-g*DY/2)))>1 )
       { double dx=0.9*DX,dy=0.9*DY;
         bColor=Black;
         tg_bar( scx(x-g*dx/2), scy(y+g*dy/2),scx(x+g*dx/2),scy(y-g*dy/2));

         fColor=Black; 
         tg_line(scx(x-g*dx/2)-1,scy(y),       scx(x-g1*dx/2),scy(y)); 
         tg_line(scx(x+g*dx/2)+1,scy(y),       scx(x+g1*dx/2),scy(y));
         tg_line(scx(x),         scy(y+g*dy/2)+1,scx(x),        scy(y+g1*dy/2));
         tg_line(scx(x),         scy(y-g*dy/2)-1,scx(x),        scy(y-g1*dy/2)); 

         fColor=White;
         tg_line(scx(x-g2*dx/2),scy(y),       scx(x-g*dx/2),scy(y));
         tg_line(scx(x+g2*dx/2),scy(y),       scx(x+g*dx/2),scy(y));
         tg_line(scx(x),       scy(y+g2*dy/2),scx(x),        scy(y+g*dy/2));
         tg_line(scx(x),       scy(y-g2*dy/2),scx(x),        scy(y-g*dy/2)); 
       }                
     } 
   }
   if (texflag)
   {  f_printf(out_tex,"\\end{picture}\n");
      texFinish();
      sprintf(buff,"%s\n%s","LaTeX output is saved in file",f_name);
      messanykey(35,15,buff);
      goto contin;              
   }
   tg_settextjustify(BottomText,LeftText);
   scrcolor(Red,bkcolor);
   tg_outtextxy(0,tg_getmaxy(), MESS);
   for(;;) 
   { double x,y;
     int i,j;
     key=inkey();
     if(key!=KB_MOUSE) break;
     if(mouse_info.but1 != 2) continue;
     x=xPhys();
     y=yPhys();
     if(x<=hMin1||x>=hMax1) continue;
     if(y<=hMin2||y>=hMax2) continue;
     i= (x-hMin1)/DX;
     j= (y-hMin2)/DY;
     
     goto_xy(1,maxRow()-1); scrcolor(Red,bkcolor); print("Mouse Info:"); 
     scrcolor(Black,bkcolor);
     print(" X=%.3E Y=%.3E F=%.2E+/-%.1E  ",x,y,f[i*nBin2+j],df[i*nBin2+j]);
   }
   if( nCol0 != maxCol() && nRow0 != maxRow() ) goto REDRAW; 
   scrcolor(bkcolor,bkcolor);
   tg_outtextxy(0,tg_getmaxy(),MESS);
   if(key==KB_ESC) goto exi; 
contin:   
   do
   {  char menustr[]="\022"
      " S = F^(power     "
      " Save plot        "
      " Math  file       "
      " LaTeX file       ";

      void * pscr = NULL;
         
      improveStr(menustr,"power", "%.2f)",P);      
      menu1(nCol0-20, 2 ,"",menustr,"n_plot2_*",&pscr,&k);

      switch (k)
      {   
           case 0:
             goto REDRAW;
           case 1: correctDouble(33,11,"S=F^(P), P= ",&P,1 ); break;
           case 2: writetable2(hMin1,hMax1,nBin1,hMin2,hMax2,nBin2,
                               f,df,upstr,xstr,ystr);
                   break;
            case 3: writeMath2(hMin1,hMax1,nBin1,hMin2,hMax2,nBin2,
                                           f,df,upstr,xstr,ystr); 
                   break;
           case 4:  
             if( !texmenu(&pictureX,&pictureY,letterSize)) break;
             nextFileName(f_name,"plot_",".tex"); 
             texStart(f_name,upstr,letterSize);
             texPicture(0,0,tg_getmaxx(),tg_getmaxy()-tg_textheight("0"),
                                                  pictureX,pictureY);
             f_printf(out_tex,"\\begin{picture}(%d,%d)(0,0)\n",pictureX,pictureY);  
             del_text(&pscr);            
             goto REDRAW; 
      }   
      if( nCol0 != maxCol() && nRow0 != maxRow() ) goto REDRAW;
   }  while (1); 
exi:
   clr_scr(FGmain,BGmain);
   put_text(&prtscr);
}
