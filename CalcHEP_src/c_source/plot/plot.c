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
#include "../ntools/include/polint.h"
#include "plot.h"


static    int logX=0;



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
  if (lg<0) return lg-0.999; else return lg;
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
   if(logX) x=log10(x/ grafminmax.xmin);
   else x -= grafminmax.xmin; 
   return X1 + chpround(xscale * x);
} 


static int  scy(double  x)
{  
   if (logScale)   x = log10(x / grafminmax.ymin); 
   else            x -= grafminmax.ymin; 
   return Y1 - chpround(yscale * x);
} 


static double  dscx(double x)
{ if(logX)  x = log10(x / grafminmax.xmin); 
   else     x-=grafminmax.xmin;
   return X1 + (xscale *x);
} 


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
{ if(logX) return   pow(10.,(mouse_info.x-X1)/xscale)*grafminmax.xmin;  
  else     return  (mouse_info.x-X1)/xscale + grafminmax.xmin;
}

static double yPhys(void)
{  if(logScale)
   return pow(10,(Y1-mouse_info.y)/yscale)*grafminmax.ymin;  
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


static int colList[]={Black,Blue,Red,Yellow,Magenta, Green,Cyan, Brown,LightGray,
               DarkGray, LightBlue, LightGreen,LightCyan, LightRed, LightMagenta,  White}; 

static char*colListTxt[]={"Black","Blue","Red","Yellow","Magenta", "Green","Cyan", "Brown","LightGray",
              "DarkGray", "LightBlue","LightGreen","LightCyan","LightRed","LightMagenta","White"};


static void  gaxes(char* upstr, char* xstr, int N, char**Y)
{

   double  xmax,xmin,ymax,ymin,aa,step; 
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

   xscale = (X2 - X1) / ( logX    ? log10(xmax/xmin) :  xmax - xmin);
   yscale = (Y1 - Y2) / (logScale ? log10(ymax/ymin) :  ymax - ymin); 

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
        islog = logX; 
        tg_settextjustify(CenterText,TopText);
        zmin=xmin;
        zmax=xmax;
        xend=X2;
        yend=Y1;              
     } else
     { 
        xy='Y';
        islog = logScale; 
        tg_settextjustify(RightText,CenterText);
        zmin=ymin;
        zmax=ymax;
        xend=X1;
        yend=Y2;              
     }
     zmax1 = zmax + fabs(zmax-zmin)*1.E-4;
     zmin1 = zmin - fabs(zmax-zmin)*1.E-4;
     axisDesign(zmin,zmax, islog,&aa, &step,&n_sub);

     tg_line(X1,Y1,xend,yend);
     len=0;
     while (aa <= zmax1) 
     {  double da;
        char snum[30];      
        int i;     
        if(aa >= zmin1)
        {  
           if(naxis==0){m=scx(aa);tg_line(m,Y1,m,Y1+hash);}
           else        {m=scy(aa);tg_line(X1-hash,m,X1,m);}
           if (islog)  doubleToStr(aa,aa,snum); else doubleToStr(aa,step,snum);
           len=MAX(len,tg_textwidth(snum));
           if(naxis==0) tg_outtextxy(m,Y1+hash,snum);
           else         tg_outtextxy(X1-hash,m,snum);     
        }
        if (islog) da = aa*step -aa ;else da=step;
        da=da/n_sub;
        for(i=1;i<n_sub;i++)
        {  aa +=da;    
           if( aa >=zmin1 && aa<=zmax1)
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
     {  tg_settextjustify(LeftText,CenterText);
        goto_xy(10,2);
        for(i=0;i<N;i++)   
        {  scrcolor(colList[i],bkcolor);
           print("%s; ",Y[i]);
        }       
     }
   }
     
}  /* GAxes */ 


static void plot_curve(double xMin,double xMax,int dim, double *f)
{  double   x,y,xx,yy;
   int i;

   double step= logX? pow(xMax/xMin,1./dim) : (xMax-xMin)/dim;
   double ymax = dscy(grafminmax.ymin);
   double ymin = dscy(grafminmax.ymax);
   
   for(i=1;i<dim;i++)
   {
      if(logX)
      {
        x =dscx(xMin*pow(step,i-0.5));
        xx=dscx(xMin*pow(step,i+0.5));
      } else
      {  x =dscx(xMin+(i-0.5)*step);
         xx=dscx(xMin+(i+0.5)*step);
      }    
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

static void plot_spline(double xMin,double xMax,int dim, int lScale, double *f)
{  double   x,y,xx,yy;
   int i;
   double step=(xMax-xMin)/dim;
   double ymax = dscy(grafminmax.ymin);
   double ymin = dscy(grafminmax.ymax);
   double xgrid[300];
   double sdim=300;
   
   if(lScale)for(i=0;i<dim;i++) xgrid[i]=xMin*pow(xMax/xMin, (i+0.5)/(double)dim);
   else      for(i=0;i<dim;i++) xgrid[i]=xMin+(i+0.5)*step; 

   step=(xMax-xMin)/sdim;
          
   for(i=1;i<sdim;i++)
   {  double x1,x2;
      if(lScale)
      { x1=xMin*pow(xMax/xMin,(i+0.5)/sdim);
        x2=xMin*pow(xMax/xMin,(i+1.5)/sdim);
      } else  
      { x1=xMin+(i+0.5)*step;
        x2=x1+step;
      }  
      x =dscx(x1);
      xx=dscx(x2); 
      y =dscy(polint3(x1,dim,xgrid,f));
      yy=dscy(polint3(x2,dim,xgrid,f));
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


static int writetable1(char * filename, double xMin,double xMax, double yMin, double yMax, int logScale, 
     int dim,char*x_str,char*upstr, 
 int N, double**f,double**df,char**Y)  
{ 
   FILE *     outfile;
   int        i,k;
   
   double dx=xMax-xMin;
   double x1=xMin+0.35*dx, x2=xMin+0.43*dx, x3=xMin+0.45*dx;
   double  yl;
   int ilab=0;
   double step=(xMax-xMin)/dim;

   outfile=fopen(filename,"w");
   if(outfile==NULL) return 1;
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
   fprintf(outfile,"#xScale %d\n",logX);
       
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
     { fprintf(outfile," %-12E  ",f[k][i]);
       if(df[k]) fprintf(outfile,"  %-12E",df[k][i]);
     } 
     fprintf(outfile,"\n");
   }
   fclose(outfile);
   return 0;
}

static void  writeGNU(char*fname, double xMin,double xMax, double yMin, double yMax, int logScale, 
     int dim,char*x_str,char*upstr, 
 int N, char* type,char**Y)  
{  char       filename[100], buff[STRSIZ],command[100];
   FILE *     outfile;
   int        i,k;
   
   double dx=xMax-xMin;
   double x1=xMin+0.35*dx, x2=xMin+0.43*dx, x3=xMin+0.45*dx;
   double  yl;
   int ilab=0;
   double step=(xMax-xMin)/dim;
   strcpy(filename,fname);
   for(i=strlen(filename)-1; i>=0;i--) if(filename[i]=='.') { filename[i]=0; break; }
   strcat(filename,".gnu");
   outfile=fopen(filename,"w");


   fprintf(outfile," set title  '%s'\n",upstr);
   fprintf(outfile," set xlabel '%s'\n",x_str);
   
   fprintf(outfile," xMin=%E\n",xMin);
   fprintf(outfile," xMax=%E\n",xMax);
   fprintf(outfile," yMin=%E\n",yMin);
   fprintf(outfile," yMax=%E\n",yMax);
   fprintf(outfile," dim=%d\n",dim);
   
   if(logScale)  fprintf(outfile," set logscale y\n");
   if(logX)
   {  fprintf(outfile," set logscale x\n"); 
      fprintf(outfile," X(n)=xMin*(xMax/xMin)**((n+0.5)/dim)\n");
   } 
   else fprintf(outfile," X(n)=xMin+(n+0.5)/dim*(xMax-xMin)\n");

   fprintf(outfile," plot[xMin:xMax][yMin:yMax]");
   
   for(i=1,k=0;k<N;k++,i++)
   { 
//      printf("type[%d]=%c\n",k,type[k]);
      if(type[k]=='h')              
      {  fprintf(outfile,"'%s' using (X($0)):%d:%d  w error ti '%s'",
                 fname,i,i+1,Y[k]);
         i++;
      } else 
      {
         fprintf(outfile,
           "'%s' using (X($0)):%d w l ti '%s'",
                  fname,i,Y[k]);
         if(type[k]=='s')  fprintf(outfile," smooth csplines");
   
      }      
       if(k<N-1) fprintf(outfile,", ");else fprintf(outfile,"\n");               
   } 
        
   fclose(outfile);
   
   messanykey(10,12,filename);
}



static void  writePAW(char *fname,  double xMin,double xMax, double yMin, double yMax, int logScale, 
     int dim,char*x_str,char*upstr, 
 int N, char*type,char**Y)  
{  char       filename[100], buff[STRSIZ],command[100];
   FILE *     outfile;
   int        i,k;
   
   double dx=xMax-xMin;
   double x1=xMin+0.35*dx, x2=xMin+0.43*dx, x3=xMin+0.45*dx;
   double  yl;
   int ilab=0;
   double step=(xMax-xMin)/dim;

   strcpy(filename,fname);

   for(i=strlen(filename)-1; i>=0;i--) if(filename[i]=='.') { filename[i]=0; break; }
   strcat(filename,".kumac");
   
   outfile=fopen(filename,"w");
          
   fprintf(outfile,"*  http://paw.web.cern.ch/paw/allfaqs.html\n");
   fprintf(outfile," opt zfl\n");
   fprintf(outfile,"set *FON -60\n");
   if(logScale)  fprintf(outfile,"opt  logy\n");
   if(logX) fprintf(outfile,"opt  logx\n");   
   fprintf(outfile," TITLE '%s'\n",upstr);
   fprintf(outfile," vector/Create XX(%d)  \n",dim);
   if(logX)
   {
     fprintf(outfile,   " sigma XX=ARRAY(%d,0#%d)\n",dim,dim-1);
     fprintf(outfile,   " sigma XX=%E*(%E)**((0.5+XX)/%d)\n",xMin,xMax/xMin,dim);
   }
   else fprintf(outfile," sigma XX=ARRAY(%d,%G#%G)\n",dim,xMin+0.5*step,xMax-0.5*step);
   sprintf(buff,"");

   for(k=0;k<N;k++)
   {  
     fprintf(outfile," vector/Create Y%d(%d)\n",k+1,dim);
     if(type[k]=='h')fprintf(outfile," vector/Create dY%d(%d)\n",k+1,dim);  
   }

   fprintf(outfile," vector/Read ");
   for(k=0;k<N;k++)
   { if(k) fprintf(outfile,",");
     fprintf(outfile,"Y%d",k+1);   
     if(type[k]=='h') fprintf(outfile,",dY%d",k+1);  
   }   
   fprintf(outfile," '%s' ' ' 'OC' '-/#/' \n",fname); 
   fprintf(outfile," 1d 10 ! 1 %G %G;  min 10 %G; max 10 %G\n",
              xMin,xMax, yMin, yMax);
   fprintf(outfile,"  his/plo 10\n"); 
   
   fprintf(outfile,"  set chhe 0.3; set lwidt 5; set hwidt 5 10\n"); 
  
   for(k=0;k<N;k++)
   {  
     fprintf(outfile,"* %s\n",Y[k]);
     fprintf(outfile," set dmod %d\n",k+1);   
     if(type[k]=='h' ) fprintf(outfile," set HCOL %d\n",k+1); 
     fprintf(outfile," set plci %i ; set txci %i  \n",k+1,k+1);
     
     if(logScale) yl=exp(log(yMax)-0.05*(k+1)*log(yMax/yMin));
     else  yl=yMax -0.05*(k+1)*(yMax-yMin);
             
     
     fprintf(outfile," dline %G %G %G %G ; itx %G %G '%s' \n",
     x1,x2, yl, yl, x3, yl, Y[k]);
 
     if(type[k]=='h')
     {
       fprintf(outfile," set dmod 1\n");
       fprintf(outfile," histogram/create/1dhisto %d '%s' %d %G %G\n",
                         k+1,upstr,dim,xMin,xMax);                
       fprintf(outfile," histogram/put/contents %d Y%d\n",k+1,k+1);
       fprintf(outfile," histogram/put/errors   %d dY%d\n",k+1,k+1);   
       fprintf(outfile," histogram/plot %d  s\n",k+1);
//       fprintf(outfile,"#PAW  atitle  '%s' '%s' \n",x_str,y_str);
     }
     else 
     { fprintf(outfile," GRAPH %d XX Y%d ",dim,k+1);
       if(type[k]=='l')  fprintf(outfile," l\n");
       else fprintf(outfile," c\n");
     }  
   }
          
   fprintf(outfile," set txci 1 \n" );
   fprintf(outfile," atitle  '%s'\n",x_str);
   fprintf(outfile," pic/print  %s.eps\n",fname); 

   fclose(outfile);
   
   messanykey(10,12,filename);
}


static void  writeROOT(char *fname,  double xMin,double xMax, double yMin, double yMax, int logScale, 
     int dim,char*x_str,char*upstr, 
 int N, char*type,char**Y)  
{  char       filename[100],basename[100], buff[STRSIZ];
   FILE *     outfile;
   int        i,k;
   char  format[150]={}, read[150]={};   
   
   strcpy(basename,fname);

   for(i=strlen(basename)-1; i>=0;i--) if(basename[i]=='.') { basename[i]=0; break; }
   sprintf(filename,"%s.C",basename);
   
   outfile=fopen(filename,"w");
fprintf(outfile,          
"#include <string>\n"
"#include <fstream>\n"
"#include <iostream>\n"
"#include <TMath.h>\n"
"#include <TGraph.h>\n"
"#include <TGraphErrors.h>\n"
"#include <TLegend.h>\n"
"#include <TStyle.h>\n"
"#include <TCanvas.h>\n"
"#include <TH1.h>\n"
"#include <TH2.h>\n"
"#include <TLatex.h>\n"
"#include <TLine.h>\n");

   for(i=strlen(basename)-1; i>=0 && basename[i]!='/';i--);
   if(i<0) i=0;
   fprintf(outfile,"void %s(void)\n{\n",basename+i);
   
   fprintf(outfile," TCanvas *Can = new TCanvas(\"c\",\"Plot\",200,10,600,400);\n");
   fprintf(outfile," double xMin=%E, xMax=%E, yMin=%E, yMax=%E;\n",xMin,xMax,yMin,yMax);
   fprintf(outfile," char buff[2000];\n");          
   fprintf(outfile," double  X[%d];\n",dim);
   if(logX)
   { 
        fprintf(outfile," Can->SetLogx();\n"); 
        fprintf(outfile," for(int i=0;i<%d;i++) X[i]=xMin*pow(Xmax/Xmin,(i+0.5)/%d);\n", dim,dim);
   } 
   else fprintf(outfile," for(int i=0;i<%d;i++) X[i]=xMin+(i+0.5)*(xMax-xMin)/%d;\n", dim,dim); 

   fprintf(outfile," double Xtext = xMin+0.3*(xMax-xMin);\n");  
   fprintf(outfile," double Ytext= yMax;\n");

   if(logScale)   
   {
     fprintf(outfile," Can->SetLogy();\n");
     fprintf(outfile," double dYtext=pow(yMax/yMin,0.05);\n");
   } else fprintf(outfile," double dYtext= (yMax-yMin)/20.;\n");
   
       
   fprintf(outfile," TH1F *hr = Can->DrawFrame(xMin,yMin,xMax,yMax);\n");
   fprintf(outfile," hr->SetXTitle(\"%s\");\n",x_str);   
   

   for(k=0;k<N;k++)
   {  
     fprintf(outfile," double  Y%d[%d];\n",k+1,dim);
     sprintf(format+strlen(format)," %%lf");
     sprintf(read+strlen(read)," ,Y%d+i",k+1);
     if(type[k]=='h')
     {  fprintf(outfile," double  dY%d[%d];\n",k+1,dim);  
        sprintf(format+strlen(format)," %%lf");
        sprintf(read+strlen(read)," ,dY%d+i",k+1);
     }
   }
   
   fprintf(outfile, "  FILE*f=fopen(\"%s\",\"r\");\n",fname);
   fprintf(outfile, "  for(int i=0;i<%d;)\n {\n",dim);
   fprintf(outfile, "  fscanf(f,\"%%[^\\n]%%*c\",buff);\n");
   fprintf(outfile, "  if(buff[0]!='#')\n  {\n");
   fprintf(outfile, "  sscanf(buff,\"%s\" %s); i++;\n", format, read);
   fprintf(outfile,"    }\n }\n");
   fprintf(outfile,"   fclose(f);\n");

   fprintf(outfile,"   double dX[%d];\n",dim);
   fprintf(outfile,"   for(int i=0;i<%d;i++) dX[i]=(xMax-xMin)/%d/3.;\n",dim,dim);
      
   for(k=1;k<=N;k++)
   { if(type[k-1]=='h') 
     {  fprintf(outfile,"   TGraphErrors *gr%d =new TGraphErrors(%d,X,Y%d,dX,dY%d);\n",k,dim,k,k);
        fprintf(outfile,"   gr%d->SetLineColor(%d);\n",k,k); 
        fprintf(outfile,"   gr%d->Draw(\"P\");\n",k);
     } else 
     {   
       fprintf(outfile, "   TGraph *gr%d = new TGraph (%d,X,Y%d);\n",k,dim,k); 
       fprintf(outfile, "   gr%d->SetLineColor(%d);\n",k,k);
       if(type[k-1]=='l') fprintf(outfile, "   gr%d->Draw(\"L\");\n",k);
       else fprintf(outfile, "   gr%d->Draw(\"C\");\n",k);
     }
     if(logScale) fprintf(outfile,"  Ytext/=dYtext;\n"); else  fprintf(outfile,"  Ytext-=dYtext;\n");
     fprintf(outfile,"  TText *t%d = new TText(Xtext,Ytext,\"%s\");\n",k,Y[k-1]);
     fprintf(outfile,"  t%d->SetTextColor(%d);\n",k,k);
     fprintf(outfile,"  t%d->Draw();\n",k); 
   }  
   fprintf(outfile,"  Can->Print(\"%s.pdf\");\n",basename);
   
   fprintf(outfile,"}\n");
   fclose(outfile);
   
   messanykey(10,12,filename);
}


static void  writeROOT2D(char *fname,  double xMin,double xMax, double yMin, double yMax,
     int dimX, int dimY, char*x_str,char*y_str,char*upstr)   
{  char       filename[100],basename[100], buff[STRSIZ];
   FILE *     outfile;
   int        i,k;
   char  format[150]={}, read[150]={};   
   
   strcpy(basename,fname);

   for(i=strlen(basename)-1; i>=0;i--) if(basename[i]=='.') { basename[i]=0; break; }
   sprintf(filename,"%s.C",basename);
   
   outfile=fopen(filename,"w");
fprintf(outfile,          
"#include <string>\n"
"#include <fstream>\n"
"#include <iostream>\n"
"#include <TMath.h>\n"
"#include <TGraph.h>\n"
"#include <TGraphErrors.h>\n"
"#include <TLegend.h>\n"
"#include <TStyle.h>\n"
"#include <TCanvas.h>\n"
"#include <TH1.h>\n"
"#include <TH2.h>\n"
"#include <TLatex.h>\n"
"#include <TLine.h>\n");

   for(i=strlen(basename)-1; i>=0 && basename[i]!='/';i--);
   if(i<0) i=0;
   fprintf(outfile,"void %s(void)\n{\n",basename+i);
   
   fprintf(outfile," TCanvas *Can = new TCanvas(\"c\",\"Plot\",200,10,600,400);\n");
   fprintf(outfile," double xMin=%E, xMax=%E, yMin=%E, yMax=%E;\n",xMin,xMax,yMin,yMax);
   fprintf(outfile," int i,j,xDim=%d, yDim=%d;\n",dimX,dimY);
   fprintf(outfile," double data[%d][%d], dData[%d][%d];\n",dimX,dimY,dimX,dimY);
   fprintf(outfile," char buff[100];\n");          
   
   fprintf(outfile," FILE*f=fopen(\"%s\",\"r\");\n",fname);
   fprintf(outfile," for(i=0;i<xDim;i++) for(j=0;j<yDim;)\n{\n");
   fprintf(outfile," fscanf(f,\"%%[^\\n]%%*c\",buff);\n");
   fprintf(outfile," if(buff[0]!='#')\n  {\n");
   fprintf(outfile," sscanf(buff,\" %%lf %%lf \",&(data[i][j]),&(dData[i][j])); j++;\n");
   fprintf(outfile,"    }\n }\n");
   fprintf(outfile," fclose(f);\n");

   fprintf(outfile," TH2D * h2 = new TH2D(\"h2\",\"%s\",xDim,xMin,xMax,yDim,yMin,yMax);\n",upstr);

   fprintf(outfile," for(int i=0;i<xDim;i++) for(j=0;j<yDim;j++) h2->SetBinContent(i+1,j+1,data[i][j]);\n");
   fprintf(outfile," h2->Draw(\"LEGO2\");\n");      
   fprintf(outfile," Can->Print(\"%s.pdf\");\n",basename);
      
         fprintf(outfile,"}\n");
            fclose(outfile);
               
                  messanykey(10,12,filename);
}                  

#ifdef OLD
   fprintf(outfile," double Xtext = xMin+0.3*(xMax-xMin);\n");  
   fprintf(outfile," double Ytext= yMax;\n");

   if(logScale)   
   {
     fprintf(outfile," Can->SetLogy();\n");
     fprintf(outfile," double dYtext=pow(yMax/yMin,0.05);\n");
   } else fprintf(outfile," double dYtext= (yMax-yMin)/20.;\n");
   
       
   fprintf(outfile," TH1F *hr = Can->DrawFrame(xMin,yMin,xMax,yMax);\n");
   fprintf(outfile," hr->SetXTitle(\"%s\");\n",x_str);   
   

   

   fprintf(outfile,"   double dX[%d];\n",dim);
   fprintf(outfile,"   for(int i=0;i<%d;i++) dX[i]=(xMax-xMin)/%d/3.;\n",dim,dim);
      
   for(k=1;k<=N;k++)
   { if(type[k-1]=='h') 
     {  fprintf(outfile,"   TGraphErrors *gr%d =new TGraphErrors(%d,X,Y%d,dX,dY%d);\n",k,dim,k,k);
        fprintf(outfile,"   gr%d->SetLineColor(%d);\n",k,k); 
        fprintf(outfile,"   gr%d->Draw(\"P\");\n",k);
     } else 
     {   
       fprintf(outfile, "   TGraph *gr%d = new TGraph (%d,X,Y%d);\n",k,dim,k); 
       fprintf(outfile, "   gr%d->SetLineColor(%d);\n",k,k);
       if(type[k-1]=='l') fprintf(outfile, "   gr%d->Draw(\"L\");\n",k);
       else fprintf(outfile, "   gr%d->Draw(\"C\");\n",k);
     }
     if(logScale) fprintf(outfile,"  Ytext/=dYtext;\n"); else  fprintf(outfile,"  Ytext-=dYtext;\n");
     fprintf(outfile,"  TText *t%d = new TText(Xtext,Ytext,\"%s\");\n",k,Y[k-1]);
     fprintf(outfile,"  t%d->SetTextColor(%d);\n",k,k);
     fprintf(outfile,"  t%d->Draw();\n",k); 
   }  
   fprintf(outfile,"  Can->Print(\"%s.pdf\");\n",basename);
   
   fprintf(outfile,"}\n");
   fclose(outfile);
   
   messanykey(10,12,filename);
}

#endif

static int  writetable2(char*baseF, double xMin, double xMax, int dimX, 
                         double yMin,double yMax, int dimY, 
double * f, double *df, char*  upstr,  char*  x_str, char*  y_str)
{
   FILE *     outfile;
   int        i;
   int ID;
   
   outfile=fopen(baseF,"w");
   if(outfile==NULL) return 1;
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

//   fprintf(outfile,"#--- PAW section ---\n");
//   fprintf(outfile,"#PAW   histogram/create/2dhisto %d '%s' %d %G %G %d %G %G\n",
//                         ID,upstr,dimY,yMin,yMax, dimX,xMin,xMax);
//    fprintf(outfile,"#PAW  vector/Create  Y%d(%d,%d)\n",ID,dimY,dimX);
//    fprintf(outfile,"#PAW  vector/Create dY%d(%d,%d)\n",ID,dimY,dimX);
//    fprintf(outfile,"#PAW  vector/Read Y%d,dY%d '%s' ' ' 'OC' '-/#/' \n",ID,ID,filename);                      
//    fprintf(outfile,"#PAW  histogram/put/contents %d Y%d\n",ID,ID);
//    fprintf(outfile,"#PAW  histogram/put/errors   %d dY%d\n",ID,ID);   
//    fprintf(outfile,"#PAW  histogram/2d/lego %d\n",ID);
//    fprintf(outfile,"#PAW  atitle  '%s' '%s' \n",y_str,x_str);
                                                                                
     
   fprintf(outfile,"#---   starting of data ---\n");
   fprintf(outfile,"#--- F ---    --- dF ---\n");
   for(i=0;i<dimX*dimY;i++) fprintf(outfile,"%-12E  %-12E\n",f[i],df[i]); 
   fclose(outfile);
   return 0;
}

static void  writePAW2D(char *fname,  double xMin,double xMax, double yMin, double yMax,
     int dimX, int dimY, char*x_str,char*y_str,char*upstr)  

{  char       filename[100], buff[STRSIZ],command[100];
   FILE *     outfile;
   int        i,k;
   
   double dx=xMax-xMin;
   double x1=xMin+0.35*dx, x2=xMin+0.43*dx, x3=xMin+0.45*dx;
   double  yl;
   int ilab=0;

   strcpy(filename,fname);
   for(i=strlen(filename)-1; i>=0;i--) if(filename[i]=='.') { filename[i]=0; break; }
   strcat(filename,".kumac");
   
   outfile=fopen(filename,"w");
          
   fprintf(outfile,"*  http://paw.web.cern.ch/paw/allfaqs.html\n");
   fprintf(outfile," opt zfl\n");
   fprintf(outfile,"set *FON -60\n");
   fprintf(outfile,"TITLE='%s'\n",upstr);
   
   fprintf(outfile,"NBINX=%d\n",dimX);
   fprintf(outfile,"NBINY=%d\n",dimY);
   fprintf(outfile,"XMIN=%E\n",xMin);
   fprintf(outfile,"XMAX=%E\n",xMax);
   fprintf(outfile,"YMIN=%E\n",yMin);
   fprintf(outfile,"YMAX=%E\n",yMax);
   fprintf(outfile,"ID=1\n");
   fprintf(outfile,"FILE='%s'\n",fname);
   fprintf(outfile,"2D [ID]  [TITLE]  [NBINX] [XMIN] [XMAX] [NBINY] [YMIN] [YMAX]\n");
   fprintf(outfile,"vector/Create H(%d,%d),dH(%d,%d)\n",dimX,dimY,dimX,dimY);
   fprintf(outfile,"vector/Read H,DH [FILE] ' ' 'OC' '-/#/'\n");
   fprintf(outfile,"PUT_VECT/CONTENTS [ID] H\n");
   fprintf(outfile,"* 2D styles:  BOX COL Z SURF SURF1 SURF2 SURF3 SURF4 LEGO LEGO1 LEGO2\n");
   
   fprintf(outfile,"HIS/PLOT [ID]  LEGO2\n");
          
//   fprintf(outfile," set txci 1 \n" );
   fprintf(outfile," atitle  '%s' '%s'\n",x_str,y_str);
   fprintf(outfile," pic/print  %s.eps\n",fname); 

   fclose(outfile);
   
   messanykey(10,12,filename);
}
                                                                                  


static void  writeMath1(char*fname, double xMin, double xMax, int dim, double * f, 
double *df, char*  upstr,  char*  x_str, char*  y_str)
{  char       filename[STRSIZ], buff[STRSIZ],basename[STRSIZ];
   FILE *     outfile;
   int        i;
   int expon, expon2;
   double mant, mant2;

   
   strcpy(basename,fname); 
   for(i=strlen(basename)-1; i>=0;i--) if(basename[i]=='.') { basename[i]=0; break; }
   sprintf(filename,"%s.math",basename);
        

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
   messanykey(10,12,filename);
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


   
int  plot_Nar( char*file, char*  title, double xMin, double xMax, char*xName,  int dim, int xScale,
               int N, double **f,double**ff,char**Y)
{   
  void *     prtscr;  
  int i,k,nCol0,nRow0,key;
  double      ymin, ymax;
  char        f_name[STRSIZ], menustr[STRSIZ], buff[STRSIZ];
  char * type=malloc(sizeof(char)*N);
  static int mouseXY=1;
  
  for(i=0;i<N;i++) 
  {            for(k=0;k<dim;k++) if(!isfinite( f[i][k])) return 1;
     if(ff[i]) for(k=0;k<dim;k++) if(!isfinite(ff[i][k])) return 1; 
  }
  
#define MESS "Esc - exit; Mouse - position; Other keys - menu"
   get_text(1,1,maxCol(),maxRow(),&prtscr);

   grafminmax.xmax=xMax;
   grafminmax.xmin=xMin;

      
   gminmax(f[0],ff[0],dim,&grafminmax.ymin,&grafminmax.ymax);
   for(i=1;i<N;i++)
   { gminmax(f[i],ff[i],dim, &ymin,&ymax);
     if(grafminmax.ymin>ymin) grafminmax.ymin=ymin;
     if(grafminmax.ymax<ymax) grafminmax.ymax=ymax;   
   }
   for(i=0;i<N;i++) if(ff[i]==NULL) type[i]='l'; else type[i]='h';
   
   ymin=grafminmax.ymin;
   ymax=grafminmax.ymax;
   
   logX=xScale;
   logScale=0;//(ymin >0 && grafminmax.ymax/grafminmax.ymin >10);   
   k = 0;      
REDRAW:
   nCol0=maxCol();
   nRow0=maxRow();
   clr_scr(fgcolor,bkcolor);
   
   gaxes(title,xName,N,Y);

   for(i=0;i<N;i++)
   {   fColor=colList[i];
       if(ff[i]) plot_hist(xMin,xMax,dim,f[i],ff[i]); 
          else  if( type[i]=='l') plot_curve(xMin,xMax,dim,f[i]);
          else  plot_spline(xMin,xMax,dim,xScale,f[i]);
   }                    
   tg_settextjustify(BottomText,LeftText);
   scrcolor(Red,bkcolor);
   tg_outtextxy(0,tg_getmaxy(), MESS);
   for(;;) 
   { double x,y;
     int i;
     key=inkey();
     if(key!=KB_MOUSE) break;
     if(mouse_info.but1 != 2) continue;
     x=xPhys();
     if(x<=xMin||x>=xMax) continue;
     y=yPhys();
//     if(y<=yMin||y>=yMax) continue;
     int  w=2+log10(fabs((xMax+xMin)/(xMax-xMin) ));    
     goto_xy(1,maxRow()-1); scrcolor(Blue,bkcolor); print("Mouse: X=%.*E ",w,x);
     if(mouseXY) print(" Y=%.2E ",yPhys()); else
     { double x0,h,al;
       int i0;
       if(logX)
       {
         h=log(xMax/xMin)/dim;
         i= log(x/xMin)/h;
         i0=log(x/xMin)/h -0.5;
         if(i0<0) i0=0; if(i0>=dim-1) i0=dim-2;
         x0=xMin*exp(h*(0.5+i0*h));
         al=log(x/x0)/h;
       }else
       {
         h=(xMax-xMin)/dim;
         i= (x-xMin)/h;
         i0=(x-xMin)/h -0.5;
         if(i0<0) i0=0; if(i0>=dim-1) i0=dim-2;
         x0=xMin+h/2+i0*h;
         al=(x-x0)/h;
       }
       print(" Func =");    
       for(k=0;k<N;k++) 
       {  scrcolor(colList[k],bkcolor);
          if(ff[k]) print(" %.2E+/-%.1E ",f[k][i],ff[k][i]);
          else print(" %.2E ",f[k][i0]*(1-al)+f[k][i0+1]*al);
       }
     }
   }
   if( nCol0 != maxCol() && nRow0 != maxRow() ) goto REDRAW; 
   scrcolor(bkcolor,bkcolor);
   tg_outtextxy(0,tg_getmaxy(),MESS);
   if(key==KB_ESC) goto exi;     
   do
   {  char sScale[20];
      void * pscr = NULL;

      if(logScale) strcpy(sScale,"Log.   "); else strcpy(sScale,"Lin.   ");  
   
      sprintf(menustr,"%c Y-max = %-9.3G Y-min = %-9.3G Y-scale = %s"
      " Spline           "
      " Mouse: show XY   "
      " Save plot in file",
      18,grafminmax.ymax,grafminmax.ymin,sScale);
 
      if(!mouseXY) improveStr(menustr,"XY","Func"); 
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
           case 4: 
           {  char * splineMenu=malloc((2+N*24)*sizeof(char));
              int m=1;
              splineMenu[0]=24;
              for(;;)
              { int l=0; 
                splineMenu[1]=0;
                for(i=0;i<N;i++)
                { if(type[i]=='l')sprintf(splineMenu+strlen(splineMenu)," %-19.19s OFF",Y[i]);
                  if(type[i]=='s')sprintf(splineMenu+strlen(splineMenu)," %-19.19s ON ",Y[i]);           
                }               
                if(strlen(splineMenu)>1)
                {             
                  menu1(nCol0-26, 6, "Spline for",splineMenu,"",NULL,&m);
                  if(m==0) break;

                  for(i=0,l=0; ;l++,i++) { if(type[i]=='h')i++;  if(l==m-1) break;  }
                  if(type[i]=='l') type[i]='s'; else type[i]='l';
                } else 
                {  messanykey(10,12,"Only histograms are presented in this plot.");
                   break;
                }       
              }  
              free( splineMenu);
           }
           break;
           case 5: mouseXY=!mouseXY; break; 
           case 6:
           { 
             char buff[150];
             char  baseF[100];
             if(file)
             {
                sprintf(buff,"This plot is already stored in the file\n%s",file);
                messanykey(10,12,buff);
                strcpy(baseF,file);   
             } else   
             {  int err;
                nextFileName(baseF,"plot_",".tab");
                correctStr(10,17,"File name: ",baseF,29,1);
                trim(baseF);
                err=writetable1(baseF,xMin,xMax,grafminmax.ymin,grafminmax.ymax,logScale,
                    dim,xName,title,N, f,ff,Y);
                if(err) 
                {         
                   sprintf(buff,"Can not create file '%s'\n",baseF);
                   messanykey(10,12,buff);
                   break;
                }   
             }
                   
             char  fileMenu[]="\011"
             " PAW     "
             " Root    "
             " Gnuplot "
//             " Math    "
              ;
             int m=1;
contin:                                                   
             
             while(m)
             {  menu1(nCol0-15, 6, "For  Package",fileMenu,"n_graph_*",NULL,&m);             
                switch(m)        
                {  case 1: writePAW(baseF,xMin,xMax,grafminmax.ymin,grafminmax.ymax,logScale,
                                    dim,xName,title,N, type,Y); 
                   case 2: writeROOT(baseF,xMin,xMax,grafminmax.ymin,grafminmax.ymax,logScale,
                                                        dim,xName,title,N, type,Y);
                   break;          
                   break;
                   case 3: writeGNU(baseF,xMin,xMax,grafminmax.ymin,grafminmax.ymax,logScale,           
                                     dim,xName,title,N, type,Y);
                   break;
//                   case 4:
//                     writeMath1(baseF,xMin,xMax,dim,f[0],ff[0],title,xName,Y[0]);  
                }
             }   
           }  
             break;
           case 0:  
             if(grafminmax.ymin >=ymax|| grafminmax.ymax <=ymin ||
               grafminmax.ymin >= grafminmax.ymax)
             { messanykey(10,10," Wrong Y-range");
               break;
             } 
             if(logScale && (grafminmax.ymax <=0)) 
             {  messanykey(10,10," Can not implement log_Y scale because Ymax<0");
                logScale = 0;
             }
                 
             if(logScale && (grafminmax.ymin <=0)) grafminmax.ymin=grafminmax.ymax*1E-4;
             
              
             if(logScale &&  grafminmax.ymax/grafminmax.ymin <=10 ) grafminmax.ymin=grafminmax.ymax/10.;
                            
             goto REDRAW;
      }   
      if( nCol0 != maxCol() && nRow0 != maxRow() ) goto REDRAW;
   }  while (1); 
      
exi:
   free(type);
   clr_scr(FGmain,BGmain);
   put_text(&prtscr);
   return 0;
}

int  plot_N(char*title, double xMin, double xMax,  char*  xName,   int dim, int xScale, int N,...)
{ int err, k,i;
                                                               
  double **f; double**ff; char**Y;
  
  va_list ap;

   if (  1.e-4 * (fabs(xMax) + fabs(xMin)) >= fabs(xMax - xMin))
   {  messanykey(10,10, "Too short interval in X axis !");
      return 2;
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

   err=plot_Nar(NULL,title,xMin,xMax,xName,dim,xScale, N,f,ff,Y);

   free(f);free(ff);free(Y);
   return err;
}

void   plot_1(double xMin, double xMax, int dim,
                    double *f, double *ff,char* title, char* xstr, char* ystr)
{ plot_Nar(NULL,title,xMin, xMax,xstr,dim,0, 1,&f,&ff, &ystr); } 
                    

void plot_2D(double hMin1,double hMax1,int nBin1,double hMin2,double hMax2,int nBin2,
            double * f,double *df,char *upstr,char* xstr,char * ystr) 
{
  int        k,nCol0,nRow0,key,i,j;
  char       f_name[STRSIZ],  buff[STRSIZ];
  double     fmax,DX,DY;
  void *     prtscr;
  double     P=1;

   logScale = 0;
   logX=0;
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
//      " LaTeX file       "
      ;

      void * pscr = NULL;
         
      improveStr(menustr,"power", "%.2f)",P);      
      menu1(nCol0-20, 2 ,"",menustr,"n_plot2_*",&pscr,&k);

      switch (k)
      {   
           case 0:
             goto REDRAW;
           case 1: correctDouble(33,11,"S=F^(P), P= ",&P,1 ); break;
           case 2: 
             {  char baseF[100];
                nextFileName(baseF,"plot2d_",".tab");
                correctStr(10,17,"File name: ",baseF,29,1);
                trim(baseF);
                writetable2(baseF,hMin1,hMax1,nBin1,hMin2,hMax2,nBin2,
                               f,df,upstr,xstr,ystr);
                char  fileMenu[]="\011"
                                 " PAW     "
                                 " ROOT    ";
                int m=1;  
                while(m)
                {  menu1(nCol0-15, 6, "For  Package",fileMenu,"n_graph_*",NULL,&m);
                   switch(m)        
                   {   
                      case 1: writePAW2D(baseF,hMin1,hMax1,hMin2,hMax2, nBin1,nBin2,
                                          xstr,ystr,upstr); break;
                      case 2: writeROOT2D(baseF,hMin1,hMax1,hMin2,hMax2, nBin1,nBin2,
                                           xstr,ystr,upstr); break;                
                      break;
                   }                                                                                                                                              
                               
                }                                     
                break;
                case 3: writeMath2(hMin1,hMax1,nBin1,hMin2,hMax2,nBin2,
                                           f,df,upstr,xstr,ystr); 
                break;
             }   
      }   
      if( nCol0 != maxCol() && nRow0 != maxRow() ) goto REDRAW;
   }  while (1); 
exi:
   clr_scr(FGmain,BGmain);
   put_text(&prtscr);
}
