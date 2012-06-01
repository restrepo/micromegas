/*
Restored version of standard drand48 and seed48 procedures 
*/

#include<stdlib.h>
#include<stdio.h>
#include"drandXX.h" 

static double float48=1./(((double) 0x10000)*((double) 0x10000)*((double) 0x10000));

static unsigned long Xlong=0x1234ABCD;
static unsigned long Xshort=0x330E;
#define Along   0x5DEEC
#define Ashort  0xE66D
#define Along16 0xDEEC0000
#define  Clong  0
#define  Cshort 0xB   
#define FIRST16 0xFFFF


double drandXX(void)
{
  unsigned  long  bot=Xshort*Ashort; 
  unsigned  long  top= Clong+ (bot>>16);

  bot=(bot&FIRST16)+Cshort;

  Xlong=top + (bot>>16)+Along*Xshort+Ashort*Xlong+((Along16*Xlong));
  Xshort=bot&FIRST16;

  return ( (double)Xshort + (Xlong&0xFFFFFFFF)*(double)0x10000 )*float48; 
}

char * seedXX(char * init)
{
  unsigned long Xlong_,Xshort_;
  static char cbuff[13];
  
  sprintf(cbuff,"%08X%04X",Xlong,Xshort);
  
  if(init)
  { if(2==sscanf(init,"%8lX%4lX",&Xlong_,&Xshort_))
     { Xlong=Xlong_;
       Xshort=Xshort_;
     } else return NULL;
  }     
  return cbuff;
}
/*
#include<limits.h>
int main(void)
{

  unsigned long i,j;
  unsigned short W[3]={10,11,12};
  char w[13];
  double q;
    
      strcpy(w,seedXX(NULL)); printf("%s\n",w); 
q=drandXX();
   
  for(i=0;i<ULONG_MAX;i++)
  { printf("%u\n",i);
    for(j=0;j< USHRT_MAX;j++)
    {if(q==drandXX()) goto exi;}
  }    
exi: printf("%u %u\n",i,j);

}
*/
