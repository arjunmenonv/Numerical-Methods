#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>

double r;
double x;
float delta = 0.01;
double r_beg = 2.4;
double r_end = 4.0;
double x_beg;
double *states;
int max_states = 100;			//limits number of oscillatory states stored to 100

void main()
{
 FILE *fptr;
 int count = 0;
 srand(time(0));
 x_beg = (double)(((rand()%99)+1)*0.01);//The beginning value of x is generated randomly within the range of (0, 1) {open interval}
 fptr = fopen("logistic.txt","w");
 states = (double *)malloc(max_states*(sizeof(double)));
 r= r_beg;
 while(r < r_end) 
 {
   x = x_beg;
    for(int i = 0; i< 10000; i++)
    {
      x = r*x*(1-x);
//The iterative relation for Logistic Map is run 750 times so that it approaches steady state. 750 was chosen as it was observed that there
//wasn't any significant difference in the maps for 750 iterations, 1000 iterations and 10000 iterations.
    }   
//
  int isosc=0;   
  int i = 0;
  int match;			     //Block of code to check for oscillatory states of x
  while(!isosc)
  {
      x = r*x*(1-x);
    if(i == 0)
    {
      *(states)  = x;		    //The first value of x is stored to the pointer states
    }
    else
     {
       match = 0;
       for(int j = 0; j<i; j++)
       {
        if(*(states + j)==x)
        {
         isosc=1;
	 match=j;		  //The number of oscillatory states is equal to j
	 break;
	}
       }
	if(!isosc)
		*(states+i)=x;
    }	
      i++;
      if(i == max_states)
      {
        break;
      }
   }
  for(int j = match; j<i-1; j++)
  {
    fprintf(fptr, "\n%lf  %lf", r,*(states+j));//The oscillating values of x are written to the file
    count++;
  }
  r += delta;
 }
 FILE *gnuplot_ptr = popen("gnuplot - persist", "w");
//Opens gnuplot. 'Persistent' is used so that the gnuplot window remains open until the program is stopped
 for(int i = 0; i<count; i++)
 {
   fprintf(gnuplot_ptr, "%s \n", "set title \"Logistic map\"");//Sets title of the gnuplot plot as Logistic Map
   fprintf(gnuplot_ptr, "%s \n", "plot 'logistic.txt' with p");//Plots values held in the file with points
 }
 pclose(gnuplot_ptr);
 fclose(fptr);
 free(states);
}

