/*
 ASSIGNMENT 4: INTEGRATION USING TRAPEZOID AND ROMBERG QUADRATURE.

 GROUP 9: ANIRUDH R(EE18B103), ARJUN MENON VADAKKEVEEDU (EE18B104), ASWIN RAJESH (EE18B106)

 C PROGRAM FOR THE PROBLEM STATEMENT IS ATTACHED BELOW. 

 **PLEASE NOTE THAT WE HAVE NOT ADDED MIDPOINT AND SIMPSON'S INTEGRATION METHODS EXPLICITLY, ALTHOUGH R(n, 1) OF ROMBERG'S METHOD IS
 ESSENTIALLY SIMPSON'S 1/3rd RULE (THE ORIGINAL PROBLEM STATEMENT DID NOT CONTAIN THE PART OF THE QUESTION ON MIDPOINT OR SIMPSON'S RULE). 

 THE ALGORITHM FOR ROMBERG QUADRATURE SLOWS DOWN TREMENDOUSLY AS THE NUMBER OF ITERATIONS OF HALVING THE INTERVAL INCREASES (BEYOND 9 THERE IS
 A SIGNIFICANT DELAY). THIS IS BECAUSE FOR THE nth ITERATION OF HALVING, THERE ARE 2^n COMPUTATIONS OF THE CUBIC SPLINE INTERPOLATION THAT IS
 REQUIRED. THE SPEED MAY BE IMPROVED BY SAVING THE INTERPOLATED VALUES INTO AN ARRAY CREATED USING HEAP MEMORY.   
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

int size_file = 0;
//
struct data_xy
{
  double index;
  long double y_val;
};
//
struct data_xy *org_ptr;
struct data_xy *down_ptr;
struct data_xy *starting;
struct data_xy *last;
//
int down_count = 5;
//
double romberg();
long double C_spline(double x);
//
int main()
{
  FILE *file_ptr;
  double temp_var;
  double trapezoid = 0;
  file_ptr = fopen("out1_test0.txt", "r");
  while(fscanf(file_ptr, "%lf", &temp_var)!=EOF)
  {
    size_file++;
  }
  fclose(file_ptr);
//
  org_ptr = (struct data_xy *)malloc(size_file*(sizeof(struct data_xy)));
  starting = org_ptr;
  file_ptr = fopen("out1_test0.txt", "r");
  org_ptr->index = 0;
  for(int i =0; i<size_file; i++)
  {
    org_ptr->y_val = temp_var;
    org_ptr->index = i;
    if(i == size_file - 1)
    {
      break;
    }
    org_ptr++;
  }
  fclose(file_ptr);
  last = org_ptr;
  org_ptr = starting;
  down_ptr = org_ptr + 1;
//Integrating using Trapezoid Rule
  for(int i = 0; i<size_file; i++)
  {
    trapezoid += ((org_ptr->y_val) + (down_ptr->y_val))*0.5;
    org_ptr++;
    down_ptr++;
  }
  org_ptr = starting;
  printf("\nThe value of the integral estimating by trapezoid rule is: %lf", trapezoid);
 //
  romberg();
  org_ptr = starting;
  free(org_ptr);
  return 0;
}

double romberg()
{
 double h = (double)size_file;
 int n;
 double x_val;
 printf("\nChoose an integer value of 'n',  the number of iterations for which the integration is to be done: ");
 scanf("%d", &n);
 double num_samples = pow(2,n);
 struct p_array
 {
  long double integral[n][n];
 };
 struct p_array *arr_ptr;
 arr_ptr = (struct p_array *)malloc(sizeof(struct p_array));
//Initialising the array
 for(int i = 0; i< n; i++)
 {
  for(int j = 0; j< n; j++)
  {
   arr_ptr->integral[i][j] = 0;
  }
 }
//
 for(int i = 0; i<n; i++)
 {
   x_val = starting->index;
   do
   {
    arr_ptr->integral[i][0] += (C_spline(x_val) + C_spline(x_val + h - 1))*h*0.5; 
    x_val += h;
   }
   while(x_val + h - 2 <= last->index);
   h /= 2;
 }
 printf("\nR(n, 0): %llf", arr_ptr->integral[n-1][0]);
//
 for(int j = 1; j<n; j++)
 { 
   for(int i = j; i<n; i++)
   {
     arr_ptr->integral[i][j] = arr_ptr->integral[i][j-1] + (arr_ptr->integral[i][j-1] - arr_ptr->integral[i-1][j-1])/(pow(4, j) - 1);
   }
  printf("\nR(n, %d", j);
  printf("): %llf", arr_ptr->integral[n-1][j]);
 }
 free(arr_ptr);
 return 0;
}
//
long double C_spline(double x)
{
 int size_array;
 double pos = 0;
 int quotient;
 double t1, t2, t3, t4, t5, t6;
 double CSval;
 size_array = ceil(size_file/down_count);
 struct inverse
 {
  long double inv_mat[size_array-2][size_array-2];
 };
 static struct inverse *inv_ptr;
 static struct inverse *start_inv;
 inv_ptr = (struct inverse *)malloc(sizeof(struct inverse));
 start_inv = inv_ptr;
//
 struct pseudo_arr
 {
  long double d_der[size_array];
  long double q[size_array];
  long double rhs[size_array-2];
 };
 struct pseudo_arr *arr_ptr1;
 static struct pseudo_arr *start_arr;
 arr_ptr1 = (struct pseudo_arr *)malloc(sizeof(struct pseudo_arr));
 start_arr = arr_ptr1;
 //
 for(int i = 0; i<size_array-2; i++)
  {
    for(int j = 0; j<size_array-2; j++)
    {
     if(i == j)
      {
       inv_ptr->inv_mat[i][j] = 1;
      }
     else
      {
       inv_ptr->inv_mat[i][j] = 0;
      }
    }
  }
//
 org_ptr = starting;
 down_ptr = starting + down_count;
 for(int i = 0; i< size_array; i++)
 {
  arr_ptr1->q[i] = 6*(down_ptr->y_val - org_ptr->y_val)/(down_count);
  down_ptr += down_count;
  org_ptr += down_count;
 }
 for(int i = 0; i<size_array-2; i++)
 {
  arr_ptr1->rhs[i] = arr_ptr1->q[i+1] - arr_ptr1->q[i];
 }
 //
 arr_ptr1->d_der[0] = 0;
 arr_ptr1->d_der[size_array-1] = 0;
 //Matrix Inversion Algorithm
 for(int i = 0; i<size_array-2; i++)
 {
  for(int j = 0; j<size_array-2; j++)
  {
   if(i!=j)
   {
    if((j == i+1) || (j == i-1))
    {
      pos = 0.25;
    }
    for(int k = 0; k< size_array-2; k++)
    {
     inv_ptr->inv_mat[j][k] -= inv_ptr->inv_mat[i][k]*pos;
    }
   }
  }
 }
//
 for(int i = 0; i<size_array-2; i++)
 {
  double sum = 0;
  for(int j = 0; j<size_array-2; j++)
  {
    sum += (inv_ptr->inv_mat[i][j])*(arr_ptr1->rhs[j]);  
  }
  arr_ptr1->d_der[i+1] = sum;
 }
//
 quotient = (double)(x/down_count);
 double temp0, temp1, temp2, temp3;
 org_ptr = starting + quotient;
 down_ptr = starting + quotient + down_count;
 temp0 = (x - (quotient+1)*(down_count));
 temp1 = pow(temp0, 3);
 temp2 = (x - quotient*down_count);
 temp3 = pow(temp2, 3);
 t1 = (-1)*((arr_ptr1->d_der[quotient])/(6*down_count))*temp1;
 t2 = ((arr_ptr1->d_der[quotient+1])/(6*down_count))*temp3;
 t3 = ((down_ptr->y_val)/(down_count))*(temp2);
 t4 = (-1)*((arr_ptr1->d_der[quotient+1])/6)*(down_count)*temp2;
 t5 = (-1)*((org_ptr->y_val)/(down_count))*temp0;
 t6 = ((arr_ptr1->d_der[quotient])/6)*(down_count)*(temp0);
 CSval = t1+t2+t3+t4+t5+t6;
//
inv_ptr = start_inv;
arr_ptr1 = start_arr;
free(start_inv);
free(arr_ptr1);
return CSval;
}

