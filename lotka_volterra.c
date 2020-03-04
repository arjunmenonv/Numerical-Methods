#include<stdlib.h>
#include<stdio.h>
#include<math.h>

double a, b, c, d;			//Constants appearing in the Differential Equation
double x = 100, y = 5, t = 0;		//x  = number of rabbits; y =number of foxes
double h = 0.005;			//h = stepsize
int size_store = 4000;			//number of iterations
//
//
 FILE *fptr;
//
double Euler(double x, double y);	//Function to run Euler's method of solvinf ODEs
double phase_plot();			//Function that generates values of y versus x and stores in a file
//
double xt(double x, double y)
{
  return a*x - b*x*y;			//  dx/dt
}
double yt(double x, double y)
{
  return d*x*y - c*y;			//  dy/dt
}
//
int main(int argc, char** argv)
{
 a = atof(argv[1]);			//The values of a, b, c and d are entered by the user during run-time
 b = atof(argv[2]);			//For stable x and y, the following are the approximate values of the constants
 c = atof(argv[3]);			// a= 10; b = 0.75; c = 1; d = 0.1
 d = atof(argv[4]);
 int opt;				//switch variable for plotting the graphs one at a time
 fptr = fopen("volterra.txt", "w");
 Euler(x, y);
 fclose(fptr);
 //
 phase_plot();
 FILE *gnuplot_ptr = popen("gnuplot", "w");
 FILE *gnu2_ptr = popen("gnuplot", "w");
 //
 printf("\n1. Press 1 to print Lotka-Volterra Map\n2. Press 2 to print Phase-Trajectory\n");
 scanf("%d", &opt);
 switch(opt)
 {
  case 1: 
    printf("\nPress Ctrl+Z to stop\n");	//Plotting Lotka-Volterra Map
    for(int i = 0; i<10000000; i++)
    {
     fprintf(gnuplot_ptr, "%s \n", "set title \"Lotka-Volterra map\"");
     fprintf(gnuplot_ptr, "%s \n", "plot 'volterra.txt' using 1:2 title \"x(t)\" with linesp, 'volterra.txt' using 1:3 title \"y(t)\" with linesp");
    }
    break;
  case 2:
    printf("\nPress Ctrl+Z to pause\n");
    for(int i = 0; i<10000000; i++)	//Plotting phase trajectory
    {
     fprintf(gnu2_ptr, "%s \n", "set title\"Phase Space\"");
     fprintf(gnu2_ptr, "%s \n", "plot 'Phase_space.txt' using 2:3 with lines");
    }
   break;
  default: printf("\nInvalid Entry");
    break;
 }
 pclose(gnu2_ptr);
 pclose(gnuplot_ptr);
 return 0;
}
//
double Euler(double x, double y)
{
 for(int i = 0; i<size_store; i++)
 {
  x += h*xt(x,y);
  y += h*yt(x,y);
  t += h;
  fprintf(fptr, "\n%lf %lf %lf", t, x, y);
 }
 return 0;
}
double phase_plot()
{
  fptr = fopen("Phase_space.txt", "w");
 x *= 2.56;			//Taking different initial values for generating different phase plots
 y *= 2.56;
 for(int i = 0; i<2; i++)
 {
  t = 0;
  x *= 0.625;
  y *= 0.625;
  Euler(x, y);
 }
 for(int i = 0; i<2; i++)
 {
  t = 0;
  x *= 0.3;
  y *= 0.3;
  Euler(x, y);
 }
 fclose(fptr);
 return 0;
}
