#ifndef SORTARR
#define SORTARR(arr,len) {int I=1;long L;  while (I < (len)) \
if (arr[I-1] <= arr[I]) I++; else {L=arr[I-1];arr[I-1]=arr[I];arr[I]=L;if(I>1)I--;}} 
#endif
