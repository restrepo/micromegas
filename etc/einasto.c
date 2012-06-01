#include"../sources/micromegas.h"

int main(int argc,char** argv)
{
  displayFunc(hProfileStd,0.001, 30,"STD");
  displayFunc(hProfileEinasto,0.001, 30,"Einasto");
}