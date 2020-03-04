/*********************************************************************************************************************************************
																	     
				QUIZ 2- EE1103 (OCTOBER 2018)
Group: 9
Members: Anirudh R (EE18B103), Arjun Menon Vadakkeveedu (EE18B104), Aswin Rajesh (EE18B106)
Date: 21st October 2018
Description: The following C code models the Landau-Lifshitz equation for a magnetic dipole placed in a magnetic field. In this program, the  		     magnetic field is assumed to be uniform wrt space and time. The C code solves the ODE for all three components of M in cartesian  		     coordinates using Euler's method and Fourth-Order Runge-Kutta method.
Version: 2.0
Commands used on TERMINAL to run the code: "gcc -o a.out Quiz2_v2.c -lm", followed by "./a.out"

**********************************************************************************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>  			//arguments to functions Euler and RK are of type bool. stdbool.h defines bool data type
#include<time.h>

//Declaration of constants and variables
double alpha = -0.05;			
double gm = 1.0;			//factor that causes precession of the magnetic dipole
int h = 1;				//step-size in milliseconds
int t = 0;				//Note that both h and t here are taken in milliseconds so that the computations 
					//involving them may be done using integers
double x = -0.99, y = 0, z = 0;		//initial values of M_x, M_y, M_z
long int no_points = (long int)(200.0*1000);
double sigma = 3*0.00001;		//standard deviation for generating random numbers belonging to a Gaussian distribution
long double error;			
int opt;				//switch variable in main()
//Declaration of differential variables for RK 
double k_0, k_1, k_2, k_3;
double l_0, l_1, l_2, l_3;
double m_0, m_1, m_2, m_3;
//Noise for x, y and z
double n_x, n_y, n_z;
//
struct file_var			
{
  int t_var;
  double x_var;
  double y_var;
  double z_var;
};
struct file_var *var_ptrE;	//These structure pointers(var_ptr and store_ptr) are used to compare values obtained in the functions  
struct file_var *var_ptrR;	//error_find(), correlate_E() and correlate_R
struct file_var *startingE;	//stores starting position of the pointers. 
struct file_var *startingR;	//This allows us to navigate through the array held by the pointer	
struct file_var *store_ptr;
struct file_var *starting_store;
//Declaration of functions
double xt(double x, double y, double z)		//Differential equations obtained for M_x, M_y and M_z after simplifying the vector-DE
{
  double xt_val;
  xt_val = alpha*x*x + 0.01*x*y*alpha + 0.01*gm*z - alpha*(x*x + y*y + z*z);
  return xt_val;
}
double yt(double x, double y, double z)
{
  double yt_val;
  yt_val = alpha*x*y + 0.01*alpha*y*y - gm*z - alpha*(x*x + y*y + z*z)*0.01;
  return yt_val;
}
double zt(double x, double y, double z)
{
  double zt_val;
  zt_val = -0.01*gm*x + gm*y +alpha*x*z + 0.01*alpha*y*z;
  return zt_val;
}
//
void Euler(bool n);		
void RK(bool n);
void ALVAR();			//Varies alpha, plots time taken for x to switch versus alpha
void normalise();	//Normalises value of x, y and z so that the points lie on a sphere always(as mandated by the Differential Eqn)
void Plot_Eu();			//Functions to plot the results using gnuplot 
void Plot_RK();
void vary_h();		//varies step-size, runs both Euler and RK4 and calls function that finds error(error_find()) wrt data with h = 0.001
//
void error_find(char c);		
void define_varptr();	//initialises pointer that stores values of M_x, M_y and M_z obtained for h = 0.001
void Noise_gen(double s);	//Generates 3 random noises belonging to a Gaussian distribution
void correlate_E();		//finds correlation factor for Euler's method
void correlate_R();		//finds correlation factor for RK4 method
//Declaration of file pointers
FILE *fptr_E, *fptr_RK, *err_fptr;;
//
int main()
{
 int choice;
//
 var_ptrE = (struct file_var *)malloc(no_points*(sizeof(struct file_var)));  
 startingE = var_ptrE;
 var_ptrR = (struct file_var *)malloc(no_points*(sizeof(struct file_var))); 
 startingR = var_ptrR;
 err_fptr = fopen("Error.txt", "w"); //Resetting Error.txt to empty file. Later, it will be opened in append mode
 fprintf(err_fptr, "\n%s", " ");
 fclose(err_fptr);

 printf("\nChoose one of the following options: \n1. Euler's method (Forward Difference)\n2. 4th Order Runge-Kutta method\n3. Vary alpha\n");
 printf("4. Vary h\n5. Add noise to Euler\n6. Add noise to RK\n7. Correlation(Euler)\n8. Correlation(RK4)\n");
 scanf("%d", &choice);
 switch(choice)
 {
  case 1: fptr_E = fopen("Euler_data.txt", "w");
          Euler(0);	//argument = 0: no noise is added
          Plot_Eu();
          fclose(fptr_E);
	  break;
  case 2: fptr_RK = fopen("RK_data.txt", "w");
	  RK(0);
          Plot_RK();
          fclose(fptr_RK);
	  break;
  case 3: ALVAR();
          break;
  case 4: define_varptr();
          vary_h();
          break;
  case 5: fptr_E = fopen("Euler_data.txt", "w");
          Euler(1);	//argument = 1: adds noise values returned by Noise_gen() to x, y and z in each iteration of computing the DE
          Plot_Eu();
          fclose(fptr_E);
	  break;
  case 6: fptr_RK = fopen("RK_data.txt", "w");
	  RK(1);
          Plot_RK();
          fclose(fptr_RK);
	  break;
  case 7: //define_varptr();
	  printf("COMPUTING... PLEASE WAIT>>>\n");
          correlate_E();	
          break;
  case 8: //define_varptr();
	  printf("COMPUTING... PLEASE WAIT>>>\n");
          correlate_R();
          break;

  default: printf("\nInvalid Entry\n");
	  break;
 }
 var_ptrE = startingE;
 free(var_ptrE);
 var_ptrR = startingR;
 free(var_ptrR);
 return 0;
}
void Euler(bool n)	//if n = 0, noise is not added; else <n = 1>, noise is added
{
 t = 0;
 x = -0.99;
 y = 0;
 z = 0;
 if(n == 0)
 {
  while(t<=200000)	//t goes till 200 s = 200*1000ms
  {
   fprintf(fptr_E, "\n%d ", t);
   fprintf(fptr_E, "%lf ", x);
   fprintf(fptr_E, "%lf ", y);
   fprintf(fptr_E, "%lf ", z);
 //
   double temp_x = h*xt(x, y, z)*0.001;
   double temp_y = h*yt(x, y, z)*0.001;
   double temp_z = h*zt(x, y, z)*0.001;
   x += temp_x;
   y += temp_y;
   z += temp_z;
   t += h;
   normalise();
  }
 }
//
 else
 {
  while(t<=200000)
  {
   fprintf(fptr_E, "\n%d ", t);
   fprintf(fptr_E, "%lf ", x);
   fprintf(fptr_E, "%lf ", y);
   fprintf(fptr_E, "%lf ", z);
  //
   Noise_gen(sigma);
   double temp_x = h*xt(x + n_x*0.1, y + n_y*0.1, z + n_z*0.1)*0.001;	//noise is added to the arguments of functions computing the DE
   double temp_y = h*yt(x + n_x*0.1, y + n_y*0.1, z + n_z*0.1)*0.001;
   double temp_z = h*zt(x + n_x*0.1, y + n_y*0.1, z + n_z*0.1)*0.001;
   Noise_gen(sigma);
   x += temp_x + n_x;	//noise is again added to the solved values. this value is taken as the x, y, z corresponding to the next iteration
   y += temp_y + n_y;
   z += temp_z + n_z;
   t += h;
   normalise();
  }
 }
}
//
void RK(bool n)
{
 t = 0;
 x = -0.99;
 y = 0;
 z = 0;
 //normalise();
 if(n == 0)
 {
  while(t<=200000)
  {
   fprintf(fptr_RK, "\n%d ", t);
   fprintf(fptr_RK, "%lf ", x);
   fprintf(fptr_RK, "%lf ", y);
   fprintf(fptr_RK, "%lf ", z);
 //
   k_0 = h*xt(x, y, z)*0.001;		
   l_0 = h*yt(x, y, z)*0.001;
   m_0 = h*zt(x, y, z)*0.001;
 //
   k_1 = h*xt(x + 0.5*k_0, y + 0.5*l_0, z + 0.5*m_0)*0.001;	//for better accuracy, the DEs are computed at intermediate points
   l_1 = h*yt(x + 0.5*k_0, y + 0.5*l_0, z + 0.5*m_0)*0.001; 	//in RK4 method.k0, k1, k2, k3 correspond to those values
   m_1 = h*zt(x + 0.5*k_0, y + 0.5*l_0, z + 0.5*m_0)*0.001;
 //
   k_2 = h*xt(x + 0.5*k_1, y + 0.5*l_1, z + 0.5*m_1)*0.001;
   l_2 = h*yt(x + 0.5*k_1, y + 0.5*l_1, z + 0.5*m_1)*0.001; 
   m_2 = h*zt(x + 0.5*k_1, y + 0.5*l_1, z + 0.5*m_1)*0.001;
 //
   k_3 = h*xt(x + k_2, y + l_2, z + m_2)*0.001;
   l_3 = h*yt(x + k_2, y + l_2, z + m_2)*0.001; 
   m_3 = h*zt(x + k_2, y + l_2, z + m_2)*0.001;
 //
   x += (k_0 + 2*k_1 + 2*k_2 + k_3)/6;
   y += (l_0 + 2*l_1 + 2*l_2 + l_3)/6; 
   z += (m_0 + 2*m_1 + 2*m_2 + m_3)/6;
   t += h;
   normalise();
 //
   
  }
 }
//
 else
 {
  while(t<=200000)
  {
   fprintf(fptr_RK, "\n%d ", t);	
   fprintf(fptr_RK, "%lf ", x);
   fprintf(fptr_RK, "%lf ", y);
   fprintf(fptr_RK, "%lf ", z);
 //
   Noise_gen(sigma);
   k_0 = h*xt(x + n_x*0.1, y + n_y*0.1, z + n_z*0.1)*0.001;	//adding noise
   l_0 = h*yt(x + n_x*0.1, y + n_y*0.1, z + n_z*0.1)*0.001;
   m_0 = h*zt(x + n_x*0.1, y + n_y*0.1, z + n_z*0.1)*0.001;
 //
   k_1 = h*xt(x + 0.5*k_0 + n_x*0.1, y + 0.5*l_0 + n_y*0.1, z + 0.5*m_0 + n_z*0.1)*0.001;
   l_1 = h*yt(x + 0.5*k_0 + n_x*0.1, y + 0.5*l_0 + n_y*0.1, z + 0.5*m_0 + n_z*0.1)*0.001; 
   m_1 = h*zt(x + 0.5*k_0 + n_x*0.1, y + 0.5*l_0 + n_y*0.1, z + 0.5*m_0 + n_z*0.1)*0.001;
 //
   k_2 = h*xt(x + 0.5*k_1 + n_x*0.1, y + 0.5*l_1 + n_y*0.1, z + 0.5*m_1 + n_z*0.1)*0.001;
   l_2 = h*yt(x + 0.5*k_1 + n_x*0.1, y + 0.5*l_1 + n_y*0.1, z + 0.5*m_1 + n_z*0.1)*0.001; 
   m_2 = h*zt(x + 0.5*k_1 + n_x*0.1, y + 0.5*l_1 + n_y*0.1, z + 0.5*m_1 + n_z*0.1)*0.001;
 //
   k_3 = h*xt(x + k_2 + n_x*0.1, y + l_2 + n_y*0.1, z + m_2 + n_z*0.1)*0.001;
   l_3 = h*yt(x + k_2 + n_x*0.1, y + l_2 + n_y*0.1, z + m_2 + n_z*0.1)*0.001; 
   m_3 = h*zt(x + k_2 + n_x*0.1, y + l_2 + n_y*0.1, z + m_2 + n_z*0.1)*0.001;
 //
   Noise_gen(sigma);
   x += (k_0 + 2*k_1 + 2*k_2 + k_3)/6 + n_x;
   y += (l_0 + 2*l_1 + 2*l_2 + l_3)/6 + n_y; 
   z += (m_0 + 2*m_1 + 2*m_2 + m_3)/6 + n_z;
   t += h;
   normalise();
 //
  }
 } 
}
//
void normalise()
{
  double coeff;			//when noise is added, the points may move out of the sphere. In order to restrict them to fall on the sphere
  coeff = (sqrt((x*x) + (y*y) + (z*z)));//always, they are normalised by dividing each value of x, y and z obtained by the current modulus
  x = x/coeff;
  y = y/coeff;
  z = z/coeff;
}
//
void Plot_Eu()			
{
 FILE *gnu_Eptr = popen("gnuplot -persistent", "w");
 printf("\nPlot the curves: \n1.Plot trajectory \n2.Plot M_x v/s t \n3.Plot M_y v/s t \n4.Plot M_z v/s t\n5.Exit\n");
 scanf("%d", &opt);
 switch(opt)
 {
  case 1:		//gnuplot commands are written using file pointers opened using "popen"
   	   fprintf(gnu_Eptr, "%s \n", "set title \"Landau-Lifshitz trajectory\"");
           fprintf(gnu_Eptr, "%s \n", "set view equal xyz");	//The points are plotted with equal ratio of x, y and z so that the resulting 
	   fprintf(gnu_Eptr, "%s \n", "set xlabel \"M_x\"");	//plot is not stretched or squeezed
           fprintf(gnu_Eptr, "%s \n", "set ylabel \"M_y\"");
     	   fprintf(gnu_Eptr, "%s \n", "set zlabel \"M_z\"");
     	   fprintf(gnu_Eptr, "%s \n", "splot 'Euler_data.txt' using 2:3:4 title \"Euler\" with dots");
           break;
  case 2:
 	   fprintf(gnu_Eptr, "%s \n", "set title \"M_x v/s t\"");	
	   fprintf(gnu_Eptr, "%s \n", "set xlabel \"time in milliseconds\"");
           fprintf(gnu_Eptr, "%s \n", "set ylabel \"value of M_x\"");
     	   fprintf(gnu_Eptr, "%s \n", "plot 'Euler_data.txt' using 1:2 title \"Euler\" with dots");
           break;
  case 3:
   	   fprintf(gnu_Eptr, "%s \n", "set title \"M_y v/s t\"");
	   fprintf(gnu_Eptr, "%s \n", "set xlabel \"time in milliseconds\"");
           fprintf(gnu_Eptr, "%s \n", "set ylabel \"value of M_y\"");
     	   fprintf(gnu_Eptr, "%s \n", "plot 'Euler_data.txt' using 1:3 title \"Euler\" with dots");
           break;
  case 4:
 	   fprintf(gnu_Eptr, "%s \n", "set title \"M_z v/s t\"");
	   fprintf(gnu_Eptr, "%s \n", "set xlabel \"time in milliseconds\"");
           fprintf(gnu_Eptr, "%s \n", "set ylabel \"value of M_z\"");
     	   fprintf(gnu_Eptr, "%s \n", "plot 'Euler_data.txt' using 1:4 title \"Euler\" with dots");
           break;
  case 5: exit(0);
 	  break;
  default: printf("\nInvalid Option\n");
	   break;
 }
 pclose(gnu_Eptr);
}
//
void Plot_RK()
{
 FILE *gnu_RKptr = popen("gnuplot -persistent", "w");
 printf("\nPlot the curves: \n1.Plot trajectory \n2.Plot M_x v/s t \n3.Plot M_y v/s t \n4.Plot M_z v/s t\n5.Exit\n");
 scanf("%d", &opt);
 switch(opt)
 {
  case 1: 
 	   fprintf(gnu_RKptr, "%s \n", "set title \"Landau-Lifshitz trajectory\"");
           fprintf(gnu_RKptr, "%s \n", "set view equal xyz");
      	   fprintf(gnu_RKptr, "%s \n", "splot 'RK_data.txt' using 2:3:4 title \"4th Order Runge-Kutta\" with dots");
           break;
  case 2: 
 	   fprintf(gnu_RKptr, "%s \n", "set title \"M_x v/s t\"");
           fprintf(gnu_RKptr, "%s \n", "set xlabel \"time in milliseconds\"");
           fprintf(gnu_RKptr, "%s \n", "set ylabel \"value of M_x\"");
     	   fprintf(gnu_RKptr, "%s \n", "plot 'RK_data.txt' using 1:2 title \"4th Order Runge-Kutta\" with dots");
           break;
  case 3:
   	   fprintf(gnu_RKptr, "%s \n", "set title \"M_y v/s t\"");
	   fprintf(gnu_RKptr, "%s \n", "set xlabel \"time in milliseconds\"");
           fprintf(gnu_RKptr, "%s \n", "set ylabel \"value of M_y\"");
     	   fprintf(gnu_RKptr, "%s \n", "plot 'RK_data.txt' using 1:3 title \"4th Order Runge-Kutta\" with dots");
           break;
  case 4:
 	   fprintf(gnu_RKptr, "%s \n", "set title \"M_z v/s t\"");
     	   fprintf(gnu_RKptr, "%s \n", "set xlabel \"time in milliseconds\"");
           fprintf(gnu_RKptr, "%s \n", "set ylabel \"value of M_z\"");
	   fprintf(gnu_RKptr, "%s \n", "plot 'RK_data.txt' using 1:4 title \"4th Order Runge-Kutta\" with dots");
           break;
  case 5: exit(0);
 	  break;
  default: printf("\nInvalid Option\n");
	   break;
 }
 pclose(gnu_RKptr);
}
//
void ALVAR()
{
 int a;			//stores time
 double b, c, d;	//stores x, y and z values respectively
 FILE *fptr1, *fptr2;
 fptr2 = fopen("Euler_dat_comp.txt", "w");	//This file stores values of time taken for switching of the dipole
 fprintf(fptr2, "\n%s", " ");
 fclose(fptr2);
 for(int i=10;i>=0;i--)
 {
  alpha = (-1)*pow(10,-1*0.1*i);		//varying alpha by powers of ten, changing the index by -0.1 in every iteration
  fptr_E = fopen("Euler_data.txt", "w");
  Euler(0);
  fclose(fptr_E);
  fptr1 = fopen("Euler_data.txt", "r");
  while(fscanf(fptr1,"%d %lf %lf %lf",&a,&b,&c,&d)!=EOF)
  { 
   if(b>=0)	//the time at minimum positive value of x is noted(closest to switching time) and written to a file
    {
     fptr2 = fopen("Euler_dat_comp.txt", "a");
     fprintf(fptr2, "%lf\t%d\n ", alpha ,a);
     fclose(fptr2);
     break;
    }
  }
 fclose(fptr1);
 }
 
 FILE *E_alpha = popen("gnuplot -persistent", "w");
 fprintf(E_alpha, "%s \n", "set title \"t_0 v/s alpha\"");		//plots alpha v/s t_switch
 fprintf(E_alpha, "%s \n", "plot 'Euler_dat_comp.txt' using 1:2 with linesp");
//
 fptr2 = fopen("RK_dat_comp.txt", "w");		//The same algorithm is repeated for RK4 method
 fprintf(fptr2, "\n%s ", " ");
 fclose(fptr2);
 for(int j=11;j>=0;j--)
 {
  alpha = (-1)*pow(10,-1*0.1*j);
  fptr_RK = fopen("RK_data.txt", "w");
  RK(0);
  fclose(fptr_RK);
  fptr1 = fopen("RK_data.txt","r");
  while(fscanf(fptr1,"%d %lf %lf %lf",&a,&b,&c,&d)!=EOF)
  {
   if(b>=0)
    {
     fptr2 = fopen("RK_dat_comp.txt", "a");
     fprintf(fptr2, "%lf\t%d\n ",alpha ,a);
     fclose(fptr2);
     break;
    }
  }
  fclose(fptr1);
 }
 FILE *RK_alpha = popen("gnuplot -persistent", "w");
 fprintf(RK_alpha, "%s \n", "set title \"t_0 v/s alpha\"");
 fprintf(RK_alpha, "%s \n", "plot 'RK_dat_comp.txt' using 1:2 with linesp");
}
//
void define_varptr()
{
  FILE *nE_fptr, *nR_fptr;
  h = 1;			//setting h back to initial value. If vary_h function has been run before, the value of h would have modified
  no_points = (long int)(200000/h);
  int tmp_t;
  double tmp_x, tmp_y, tmp_z;
//
  fptr_E = fopen("Euler_data.txt", "w");
  Euler(0);			//writes values of Euler(without noise) to Euler_data.txt
  fclose(fptr_E);
//
  fptr_RK = fopen("RK_data.txt", "w");
  RK(0);			//RK4 without noise is written to RK_data.txt
  fclose(fptr_RK);
//
  nE_fptr = fopen("Euler_data.txt", "r+");
  for(int i = 0; i<no_points; i++)
 {
  fscanf(nE_fptr, "%d %lf %lf %lf", &(var_ptrE->t_var), &(var_ptrE->x_var), &(var_ptrE->y_var), &(var_ptrE->z_var));
  var_ptrE++;				//reads data from Euler_data.txt and stores it in memory pointed to by the structure pointer
 }
  var_ptrE = startingE;
//
  nR_fptr = fopen("RK_data.txt", "r+");
   for(int i = 0; i<no_points; i++)
 {
  fscanf(nR_fptr, "%d %lf %lf %lf", &(var_ptrR->t_var), &(var_ptrR->x_var), &(var_ptrR->y_var), &(var_ptrR->z_var));
  var_ptrR++;				//corresponding process for RK_data.txt
 }
  var_ptrR = startingR;
}
//
//Function to find error between different step-sizes considering values corresponding to h = 0.001 as "accurate". 
void vary_h()
{ 
  FILE *gnu_h;
  int opt_h;
  h = 1;
  no_points = (200000/h);
  int h_new[10] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000};
  h = h_new[0];
  
  printf("Vary h for:\n1. Euler's method\n2. Runge-Kutta method\n");
  scanf("%d", &opt_h);
  switch(opt_h)
  {
   case 1: for(int i = 0; i<10; i++)
              {
  	        h = h_new[i];
    		fptr_E = fopen("Euler_data.txt", "w");
    		Euler(0);
   		fclose(fptr_E);
    		error_find('E');
              }
                gnu_h = popen("gnuplot -persistent", "w");
 		fprintf(gnu_h, "%s \n", "set title \"Error v/s h\"");
 		fprintf(gnu_h, "%s \n", "plot 'Error.txt' using 1:2 with linesp");
		break;
   case 2: for(int i = 0; i<10; i++)
              {
  	        h = h_new[i];
    		fptr_RK = fopen("RK_data.txt", "w");
    		RK(0);
   		fclose(fptr_RK);
    		error_find('R');
              }
                gnu_h = popen("gnuplot -persistent", "w");
 		fprintf(gnu_h, "%s \n", "set title \"Error v/s h\"");
 		fprintf(gnu_h, "%s \n", "plot 'Error.txt' using 1:2 with linesp");
		break;
   default: printf("\nInvalid!\n");
	    break;
  }
  h =1;
}
//
void error_find(char c)			
//if c = 'E', calculates error for Euler. else if c = 'R', calculates error for RK4. This is matched with the method of solving used in vary_h 
{
 FILE *new_ptr;
//
 no_points = (200000/h);
 store_ptr = (struct file_var *)malloc(no_points*(sizeof(struct file_var)));  	//stores values corresponding to varied values of h
 starting_store = store_ptr;
 if(c == 'E')
 {
  new_ptr = fopen("Euler_data.txt", "r+");
  for(int i = 0; i<no_points; i++)
  {
   fscanf(new_ptr, "%d %lf %lf %lf", &(store_ptr->t_var), &(store_ptr->x_var), &(store_ptr->y_var), &(store_ptr->z_var));
   store_ptr++;
  }
  fclose(new_ptr);
  store_ptr = starting_store;
  double e_x = 0, e_y = 0, e_z = 0;
  var_ptrE = startingE;
  error = 0;
  t = 0;
  while(t < 200000)
  {
   e_x = (store_ptr->x_var)-(var_ptrE->x_var);
   e_y = (store_ptr->y_var)-(var_ptrE->y_var);
   e_z = (store_ptr->z_var)-(var_ptrE->z_var);
   error += sqrt(e_x*e_x + e_y*e_y + e_z*e_z);
   store_ptr++;
   for(int i = 0; i<h; i++)
   {
    var_ptrE++;
   }
   t += h;
   }
  }
//
 if(c == 'R')
 {
  new_ptr = fopen("RK_data.txt", "r+");
  for(int i = 0; i<no_points; i++)
  {
   fscanf(new_ptr, "%d %lf %lf %lf", &(store_ptr->t_var), &(store_ptr->x_var), &(store_ptr->y_var), &(store_ptr->z_var));
   store_ptr++;
  }
  fclose(new_ptr);
  store_ptr = starting_store;
  double e_x = 0, e_y = 0, e_z = 0;
  var_ptrR = startingR;
  error = 0;
  t = 0;
  while(t < 200000)
  {
   e_x = (store_ptr->x_var)-(var_ptrR->x_var);
   e_y = (store_ptr->y_var)-(var_ptrR->y_var);
   e_z = (store_ptr->z_var)-(var_ptrR->z_var);
   error += sqrt(e_x*e_x + e_y*e_y + e_z*e_z);
   store_ptr++;
   for(int i = 0; i<h; i++)
   {
    var_ptrR++;
   }
   t += h;
   }
  }
  error /= no_points;
  store_ptr = starting_store;
  free(store_ptr);
  err_fptr = fopen("Error.txt", "a");
  fprintf(err_fptr, "\n%d %Lf", h, error);
  fclose(err_fptr);
//
}
//
void Noise_gen(double s)
{
 #define pi 3.1415
 double u_x, u_y, u_z;		//a set of uniformly distributed random numbers are generated using rand()
 double v_x, v_y, v_z;		//these numbers are then used to generate a Gaussian distribution
/******************************************************************************************************************************************

	For a normal distribution with mean p and standard deviation s, the equation of each point in the distribution
	may be represented as:
	y = (N*).exp((-(x-p)^2)/(2s^2))
	where y is a number belonging to the normal distribution
	      x is a random number from a set that is uniformly distributed
	      N* = normalising factor = (1/sqrt(2.pi.s^2)) 
	This will generate positive values of y. For a random walk, there may be a kick in the positive direction or in 
	the negative direction, with almost equal probability of the kick going either way.

	In order to take that into account, we are multiplying the obtained y by a "toss" variable, a number that may be 
	1 or -1 at random, with equal probability. 	

*******************************************************************************************************************************************/
 double toss;
 u_x = (rand()%(1000))/10000;		//u_x, u_y, u_z belongs to (0.0000, 0.0999) 
 u_y = (rand()%(1000))/10000;		//larger values of the random variables will cause the path to deviate significantly
 u_z = (rand()%(1000))/10000;
 //
 v_x = ((exp((-1*u_x*u_x)/(2*s*s)))/(sqrt(2*pi*s*s)))*pow(10, -7);	//The generated noise values are further scaled down.
 v_y = ((exp((-1*u_y*u_y)/(2*s*s)))/(sqrt(2*pi*s*s)))*pow(10, -7);	//The factor for scaling down was determined through trial and error
 v_z = ((exp((-1*u_z*u_z)/(2*s*s)))/(sqrt(2*pi*s*s)))*pow(10, -7);	//by checking the plots for different values of the factor
 //
 toss = rand()%2;
 n_x = v_x*(pow(-1, toss));
 toss = rand()%2;
 n_y = v_y*(pow(-1, toss));
 toss = rand()%2;
 n_z = v_z*(pow(-1, toss));
}
//
void correlate_E()
{
 FILE *new_ptr, *org_ptr, *corr_ptr;		
 FILE *gnu_corrptr;
 int choice_corr;
 h = 1;
 double corr_x;		//cross-correlated values of x, y and z
 double corr_y;
 double corr_z;
 //
 double acor_x;		//auto-correlated values of x, y and z
 double acor_y;
 double acor_z; 
 //
 double coeff_x;	//calculates correlation coefficient of x, y and z. <coeff = cross/auto>
 double coeff_y;
 double coeff_z;
 no_points = (200000/h);
 //
 fptr_E = fopen("Euler_data.txt", "w");
 Euler(0);		//writes noise-free data to Euler_data.txt
 fclose(fptr_E);
 sigma = pow(10, -5);		//sigma is varied, starting from 10^-5
 corr_ptr = fopen("Correlate.txt", "w");
 for(int i = 0; i < 50 ; i++)
 {
   fptr_E = fopen("Euler_noise.txt", "w");
   Euler(1);		//writes noisy data to Euler_noise.txt
   fclose(fptr_E);
   org_ptr = fopen("Euler_data.txt", "r+");
   new_ptr = fopen("Euler_noise.txt", "r+");
//
   store_ptr = starting_store;
   corr_x = 0;
   corr_y = 0;
   corr_z = 0;
   acor_x = 0;
   acor_y = 0;
   acor_z = 0;
   for(int p = 0; p<200000; p++)
    {
     double t_o, x_o, y_o, z_o;
     double t_n, x_n, y_n, z_n;
     fscanf(org_ptr, "%lf %lf %lf %lf", &(t_o), &(x_o), &(y_o), &(z_o));	//reading noise-free data
     fscanf(new_ptr, "%lf %lf %lf %lf", &(t_n), &(x_n), &(y_n), &(z_n));	//reading noisy data
/************************************************************************************************************************************

	Cross-correlation y as a function of lag L is defined as:
	y(L) = integral{x_1(t).x_2(t+L)dt};
	the limits of integration are from -inf to +inf. x_1(t) and x_2(t) are two different 
	functions of time. Auto-correlation is done for the same function. We have computed the correlations for lag = 0
	Numerically, it can replaced by a summation, as we have done here.  

************************************************************************************************************************************/
     corr_x += (x_o)*(x_n)*h;
     corr_y += (y_o)*(y_n)*h;
     corr_z += (z_o)*(z_n)*h;
//
     acor_x += (x_o)*(x_o)*h;
     acor_y += (y_o)*(y_o)*h;
     acor_z += (z_o)*(z_o)*h;
    }
   fclose(org_ptr);
   fclose(new_ptr);
     if((acor_x == 0 || acor_y == 0) || (acor_z == 0))
     {
      printf("\nDIVISION BY ZERO! ABORTING");
      exit(0);
     }
     coeff_x = corr_x/acor_x;		//correlation ratio for x, y and z
     coeff_y = corr_y/acor_y;
     coeff_z = corr_z/acor_z;
     fprintf(corr_ptr, "\n%lf ", sigma*pow(10, 5));	//scaling up sigma while storing in file for ease of plotting
     fprintf(corr_ptr, "%lf ", coeff_x);
     fprintf(corr_ptr, "%lf ", coeff_y);
     fprintf(corr_ptr, "%lf ", coeff_z);
   sigma += 0.000005;
 }
 fclose(corr_ptr);
 
 gnu_corrptr = popen("gnuplot -persistent", "w");
 printf("Choose an option to plot the correlation:\n1.of M_x v/s sigma\n2.of M_y v/s sigma\n3.of M_z v/s sigma\n");
 scanf("%d", &choice_corr);
 switch(choice_corr)
 {
  case 1: 
 	   fprintf(gnu_corrptr, "%s \n", "set title \"Corr(M_x) v/s sigma - Euler\"");
           fprintf(gnu_corrptr, "%s \n", "set xlabel \"sigma scaled up by 10^5 for convenience of plotting\"");
           fprintf(gnu_corrptr, "%s \n", "plot 'Correlate.txt' using 1:2 with linesp");
           break;
  case 2:
   	   fprintf(gnu_corrptr, "%s \n", "set title \"Corr(M_y) v/s sigma - Euler\"");
           fprintf(gnu_corrptr, "%s \n", "set xlabel \"sigma scaled up by 10^5 for convenience of plotting\"");
     	   fprintf(gnu_corrptr, "%s \n", "plot 'Correlate.txt' using 1:3 with linesp");
           break;
  case 3:
 	   fprintf(gnu_corrptr, "%s \n", "set title \"Corr(M_z) v/s sigma - Euler\"");
           fprintf(gnu_corrptr, "%s \n", "set xlabel \"sigma scaled up by 10^5 for convenience of plotting\"");
     	   fprintf(gnu_corrptr, "%s \n", "plot 'Correlate.txt' using 1:4 with linesp");
           break;
  default: printf("\nInvalid!\n");
	   break;
 }
 pclose(gnu_corrptr);
}
//
void correlate_R()		//same function as correlate_E(), with only the pointers and files that are being dealt with changing
{
 FILE *new_ptr, *org_ptr, *corr_ptr;
 FILE *gnu_corrptr;
 int choice_corr;
 h = 1;
 double corr_x;
 double corr_y;
 double corr_z;
 //
 double acor_x;
 double acor_y;
 double acor_z; 
 //
 double coeff_x;
 double coeff_y;
 double coeff_z;
 no_points = (200000/h);
 fptr_RK = fopen("RK_data.txt", "w");
 RK(0);
 fclose(fptr_RK);
   
 sigma = pow(10, -5);
 corr_ptr = fopen("Correlate.txt", "w");
 for(int i = 0; i < 50 ; i++)
 {
   fptr_RK = fopen("RK_noise.txt", "w");
   RK(1);
   fclose(fptr_RK);
   org_ptr = fopen("RK_data.txt", "r+");
   new_ptr = fopen("RK_noise.txt", "r+");
//
   store_ptr = starting_store;
   corr_x = 0;
   corr_y = 0;
   corr_z = 0;
   acor_x = 0;
   acor_y = 0;
   acor_z = 0;
   for(int p = 0; p<200000; p++)
    {
     double t_o, x_o, y_o, z_o;
     double t_n, x_n, y_n, z_n;
  //
     fscanf(org_ptr, "%lf %lf %lf %lf", &(t_o), &(x_o), &(y_o), &(z_o));
     fscanf(new_ptr, "%lf %lf %lf %lf", &(t_n), &(x_n), &(y_n), &(z_n));
     corr_x += (x_o)*(x_n)*h;
     corr_y += (y_o)*(y_n)*h;
     corr_z += (z_o)*(z_n)*h;
//
     acor_x += (x_o)*(x_o)*h;
     acor_y += (y_o)*(y_o)*h;
     acor_z += (z_o)*(z_o)*h;
    }
     if((acor_x == 0 || acor_y == 0) || (acor_z == 0))
     {
      printf("\nDIVISION BY ZERO! ABORTING");
      exit(0);
     }
     coeff_x = corr_x/acor_x;
     coeff_y = corr_y/acor_y;
     coeff_z = corr_z/acor_z;
     fprintf(corr_ptr, "\n%lf ", sigma*pow(10, 5));
     fprintf(corr_ptr, "%lf ", coeff_x);
     fprintf(corr_ptr, "%lf ", coeff_y);
     fprintf(corr_ptr, "%lf ", coeff_z);
   sigma += 0.000005;
 }
 fclose(new_ptr);
 fclose(org_ptr);
 fclose(corr_ptr);
 
 gnu_corrptr = popen("gnuplot -persistent", "w");
 printf("Choose an option to plot the correlation:\n1.of M_x v/s sigma\n2.of M_y v/s sigma\n3.of M_z v/s sigma\n");
 scanf("%d", &choice_corr);
 switch(choice_corr)
 {
  case 1: 
 	   fprintf(gnu_corrptr, "%s \n", "set title \"Corr(M_x) v/s sigma - RK4\"");
           fprintf(gnu_corrptr, "%s \n", "set xlabel \"sigma scaled up by 10^5 for convenience of plotting\"");
           fprintf(gnu_corrptr, "%s \n", "plot 'Correlate.txt' using 1:2 with linesp");
           break;
  case 2:
   	   fprintf(gnu_corrptr, "%s \n", "set title \"Corr(M_y) v/s sigma - RK4\"");
           fprintf(gnu_corrptr, "%s \n", "set xlabel \"sigma scaled up by 10^5 for convenience of plotting\"");
     	   fprintf(gnu_corrptr, "%s \n", "plot 'Correlate.txt' using 1:3 with linesp");
           break;
  case 3:
 	   fprintf(gnu_corrptr, "%s \n", "set title \"Corr(M_z) v/s sigma - RK4\"");
           fprintf(gnu_corrptr, "%s \n", "set xlabel \"sigma scaled up by 10^5 for convenience of plotting\"");
     	   fprintf(gnu_corrptr, "%s \n", "plot 'Correlate.txt' using 1:4 with linesp");
           break;
  default: printf("\nInvalid!\n");
	   break;
 }
 pclose(gnu_corrptr);
 store_ptr = starting_store;
 free(store_ptr);
}

