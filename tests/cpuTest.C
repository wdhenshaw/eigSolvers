#include <stdio.h>
#include <math.h>
#include <sys/time.h>
// #include <sys/resource.h>

#define Real double


/* This measures wall clock time */
/* Here is the version that works in serial */
void second( Real *time )
{
  static struct timeval firstTime;
  static int called=0;
  
  struct timeval theTime;
  gettimeofday( &theTime, 0 );
  if (!called)
  {
    firstTime = theTime;
    called=1;
  }
  
  *time=(theTime.tv_sec - firstTime.tv_sec) + (theTime.tv_usec - firstTime.tv_usec)*1.e-6;
  /*  printf("time inside second_: %e\n", *time); */
}


//=============================================================
//  Return the current value of the amount of CPU time used
//============================================================
Real getCPU()
{
  Real time;
  second( &time );

  return time;
}    



// ==================================================================================
// ================================ MAIN ============================================
// ==================================================================================
int main(int argc,char **argv)
{

  Real cpu0 = getCPU();

  int maxIt=10000; 
  for( int i=0; i<maxIt ; i++ )
  {
    Real s = sin( i*.1 );
  }

  Real cpu = getCPU() - cpu0;

  printf("maxIt=%6d, total time = %9.2e (s)\n",maxIt,cpu);

  cpu0 = getCPU();
  maxIt=100000; 
  for( int i=0; i<maxIt ; i++ )
  {
    Real s = sin( i*.1 );
  }

  cpu = getCPU() -cpu0;
  printf("maxIt=%6d, total time = %9.2e (s)\n",maxIt,cpu);



  return 0;
}