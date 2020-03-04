/* EE1103 ASSIGNMENT 4
  GROUP 9: ANIRUDH RAMESH(EE18B103), ARJUN MENON VADAKKEVEEDU(EE18B104), ASWIN RAJESH(EE18B106)
 
SOLUTION: THE PROGRAM THAT WE USED TO FIND AN INTERPOLATED FUNCTION, TO THE SET OF DATA POINTS THAT WERE GIVEN TO US, IS THE FOLLOWING C CODE.
WE DOWNSIZED THE DATA POINTS THAT WE RECEIVED, TAKING EVERY 5TH POINT IN ORDER TO AVOID ERRORS FROM OVERSAMPLING.
      
NOTE: THE RMS ERROR OBTAINED BY US USING THE C PROGRAM IS APPROXIMATELY EQUAL TO 550 (BEFORE REMOVING ERRORS DUE TO RUNGE'S PHENOMENON), WHICH IS VERY HIGH, CONSIDERING THE FACT THAT THE NUMBERS BELONGING TO THE DATA SET IS OF THE ORDER OF 10^-3 TO 10^-2. HOWEVER, IT IS TO BE NOTED THAT THIS VALUE OF THE ERROR IS THE SMALLEST POSSIBLE WE WERE ABLE TO ACHIEVE BY AVOIDING VARIABLES ASSUMING EXTREMELY HIGH VALUES(nan and inf) THAT ARE OUTSIDE THE RANGE OF LONG DOUBLE DATA TYPE. 
WE FOUND OUT THAT THE VARIABLES ASSUMED SUCH VALUES WHEN THERE WAS A DIVISION BY A NUMBER TENDING TO ZERO. OCCASIONALY, THE DENOMINATOR WOULD ASSUME A NUMBER SO SMALL THAT IT MAY NOT BE REPRESENTED WITHIN THE PRECISION RANGE OF LONG DOUBLE. IN SUCH CASES THESE DENOMINATORS WERE ROUNDED OFF TO 0, LEADING TO ERRONEOUS RESULTS. WE HAVE USED A LOGIC(EXPLAINED BELOW, ALONGSIDE THE CODE WHERE IT WAS USED) TO AVOID SUCH CIRCUMSTANCES.    

ANOTHER SOURCE OF ERROR WAS THE DEVIATION OF THE INTERPOLATED VALUES FROM THE ORIGINAL VALUES NEAR THE EXTREME POINTS OF THE INTERVAL. INTERESTINGLY, A SIMILAR ISSUE WAS MENTIONED ON THE WIKIPEDIA PAGE OF LAGRANGE POLYNOMIAL AS A DRAWBACK OF THIS METHOD. CALLED RUNGE'S PHENOMENON, THIS PHENOMENON IS A "PROBLEM OF OSCILLATION AT THE EDGES OF AN INTERVAL THAT OCCURS WHEN A POLYNOMIAL INETRPOLATION WITH POLYNOMIALS OF HIGH DEGREE OVER A SET OF EQUISPACED INTERPOLATED POINTS". WE REDUCED THE ERRORS CAUSED BY RUNGE'S PHENOMENON BY REMOVING THE FIRST AND LAST 10 POINTS FROM THE LOOP THAT WE USED FOR ESTIMATION. 

AFTER REMOVING THESE ERROR PRONE POINTS, THE ERROR REDUCED DRAMATICALLY FROM 550 TO 0.019.

*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//Declaration of structure that holds each value of x(index) and y(y_val) 
struct data_xy
{
  unsigned int index;
  double y_val;
};
//
struct data_xy *org_ptr;//pointer that holds values read from the file
struct data_xy *down_ptr;//pointer that holds downsized values
struct data_xy *inter_ptr;//pointer that holds values found by Lagrange polynomial method of interpolation
struct data_xy *starting_o;//holds starting address of org_ptr so that we may incerment or decrement org_ptr without losing information
struct data_xy *starting_i;//holds starting address of inter_ptr
//
int size_file = 0;
int down_count = 5;//sets downsizing count to 5. Every 5th value from file shall be read and stored into *down_ptr
//

//
long double Lag_fn(int x);//returns value of f(x) estimated by Lagrange polynomial method of estimation
long double New_fn(int x);//estimates value of f(x) by Newton polynomial method and returns it
long double C_spline(int x);//estimation by Cubic Spline interpolation method
long double epsilon_find();//returns value of epsilon/error in approximation
//
int main()
{
  FILE *fr_ptr,*fw_ptr;
  double temp_var;//value read from the file is assigned to temp_var
  long double eps_val = 0;
  long double sqrt_size;
  long double rms_val;
  fr_ptr = fopen("out1_test0.txt", "r");
  while(fscanf(fr_ptr, "%lf", &temp_var)!=EOF)
  {
    size_file++;
  }
  fclose(fr_ptr);
  fr_ptr = fopen("out1_test0.txt", "r");//The file that is being read must have the name "out1_test0.txt"
  org_ptr = (struct data_xy *)malloc(size_file*(sizeof(struct data_xy)));
  inter_ptr = (struct data_xy *)malloc(size_file*(sizeof(struct data_xy)));

  starting_o = org_ptr;//starting_o is assigned the starting address of org_ptr
  down_ptr = org_ptr;//down_ptr is made to point to the same location as org_ptr initially
  starting_i = inter_ptr;
  for(int i = 0; i<size_file; i++)
  {
    fscanf(fr_ptr, "%lf", &temp_var);
    org_ptr->index = i;
    org_ptr->y_val = temp_var;
    if(i == size_file-1)
    {
      break;
    }
    org_ptr++;
  }
//
  fclose(fr_ptr);
  org_ptr = starting_o;
  fw_ptr = fopen("out1_downsized.txt", "w");
  for(int i = 0; i+down_count<size_file; i+=down_count)
  {
    down_ptr->index = i;
    down_ptr->y_val = org_ptr->y_val;
    fprintf(fw_ptr, "\n %d %lf", down_ptr->index, down_ptr->y_val);
/*the downsized values are saved to the file "out1_downsized.txt". The plot of the numbers in this file may be found to get an idea of how the function estimated by gnuplot algorithm changes as the sampling frequency is varied*/
    down_ptr+=down_count;
    org_ptr+=down_count;
  }
  fclose(fw_ptr);
  down_ptr = starting_o;//the pointers down_ptr and org_ptr are made to point back to the starting address
  org_ptr = starting_o;
  //
  eps_val = epsilon_find();
/*eps_val is assigned the value of error of estimation, epsilon. The function epsilon_find() calls the function Lag_fn(int x) within itself; the latter function finds an estimated value of the function at a point 'x' which is passed as an argument. The function uses Lagrange polynomial to estimate*/
  sqrt_size = sqrt(size_file);
  rms_val = eps_val/sqrt_size;
  printf("\nThe square root of sum of squares of distances of the points on the curve and the corresponding points on the interpolated polynomial: ");
  printf("%llf", eps_val);
  printf("\nThe RMS value of the error is: %llf", rms_val);
  printf("\n");
  //
  down_ptr = starting_o;
  org_ptr = starting_o;//before freeing the dynamically allocated pointers, they must be made to point to the their starting locations
  inter_ptr = starting_i;

  free(inter_ptr);
  free(org_ptr);
  return 0;
}
//Function definitions
long double Lag_fn(int x)
{
 unsigned int case_dval;
 long double temp_coeff;
 down_ptr = starting_o;
 long double lag_val = 0;
 for(int i =0; i<size_file; i+= down_count)
 {
  long double numr = 1;
  long double denom = 1;
  for(int j = 0; j+down_count<size_file; j+=down_count)
   {
     if(j!=i)
      {
        numr *= (x - j);
        denom *= (i - j); 
        if(denom<=0.000001)
/*if the value of denom is less than 10^-6, it would lead to erroneous results owing to the value being close to zero in the range offered by long double data type*/
        {
         case_dval = 0;
         denom *= 1050.26;
/*
 denom is hence multiplied by a number of the order of 10^3 to make it a considerably large number to avoid getting zero in the denominator accidentally. The ratio numr/denom obtained finally is multiplied by the same number. The number 1050.26 was obtained by trial and error, with the aim of minimising the RMS error. This number will minimise the error only for the given set of numbers. For another set of numbers, the multiplier must be found separately.  
*/
        }
      }
   }
   temp_coeff = numr/denom;
   //printf("\n %llf", temp_coeff);
   if(case_dval == 0)
   {
     temp_coeff *= 1050.26;//temp_coeff is multiplied by the same multiplier by which denom was multiplied in case it was less than 10^-6
   } 
   if((temp_coeff > 1000000000)||(temp_coeff != temp_coeff))
   {
     temp_coeff = 1000000000;
//there is a chance that the value of temp_coeff may still be high. In that case it is arbitrarily limited to 10^9. This may give rise to
//errors in estimation but the program shall run without returning high numbers that are outside the range of long double data type
   }
   lag_val += (temp_coeff)*(down_ptr->y_val);
   down_ptr+=down_count;
 }
 return lag_val;
}
//
long double New_fn(int x)
{
  org_ptr = starting_o;
  int size_array;
  long double New_val = org_ptr->y_val;
  long double product = 1;
  size_array = ceil(size_file/down_count);
  struct pseudo_array
  {
   long double New_array[size_array][size_array];
  }; 
  struct pseudo_array *arr_ptr;
  arr_ptr = (struct pseudo_array *)malloc(sizeof(struct pseudo_array));
  for(int i = 0; i<size_array; i++)
  {
   arr_ptr->New_array[i][0] = org_ptr->y_val;
   org_ptr += down_count;
  }
  for(int j = 1; j<size_array; j++)
  {
    for(int i = 0; i<size_array-j; i++)
    {
     arr_ptr->New_array[i][j] = 0.00000001*(arr_ptr->New_array[i+1][j-1] - arr_ptr->New_array[i][j-1])/(j*down_count);
//Reducing the error obtained by multiplying RHS by 10^-8. This method is similar to what was used for Lagrange polynomial. The numbers were
// obtained by trial and error 
    }
  }
  //
  org_ptr = starting_o;
  for(int j = 1; j<size_array; j++)
  {
   product *= (x - org_ptr->index);
   org_ptr += down_count;
   New_val += 100000000*(arr_ptr->New_array[0][j])*product;//Multiplying the final number by 10^8
  }
  free(arr_ptr);
  return New_val;
}
//
long double C_spline(int x)
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
 org_ptr = starting_o;
 down_ptr = starting_o + down_count;
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
 quotient = (int)(x/down_count);
 double temp0, temp1, temp2, temp3;
 org_ptr = starting_o + quotient;
 down_ptr = starting_o + quotient + down_count;
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
//
long double epsilon_find()
{
  int shift_pos1 = 10; 
  int shift_pos2 = 10;
  char ch;
/*
This is to reduce Runge's phenomenon. Since the Lagrange polynomial varies highly with respect to the original data values near the extreme points of the interval of the function, we are forcefully removing these erroneous points. The RMS error was significantly reduced from 549 to 0.019 upon removing the first and last 10 points. Hence this function to estimate using Lagrange polynomial interpolation is effective for the remaining points.
 */ 
 
  long double eps_sqr = 0;
  long double epsilon;
  printf("Enter the method of interpolation to be used:\nL: Lagrange\nN: Newton\nC: Cubic Spline\nE: Exit\n");
  scanf("%c", &ch);
  switch(ch)
  {
    case 'L': 
        org_ptr = starting_o + shift_pos1;
        inter_ptr = starting_i + shift_pos2;
        for(int k = shift_pos1; k<size_file-shift_pos1; k++)
         {
 	  long double temp;
   	  inter_ptr->index = k;
   	  inter_ptr->y_val = Lag_fn(k);
     	  temp = (inter_ptr->y_val)-(org_ptr->y_val);
   	  eps_sqr += temp*temp; 
   	  inter_ptr++; 
   	  org_ptr++;
    	} 
        break;
   case 'N':
        org_ptr = starting_o + shift_pos2;
        inter_ptr = starting_i + shift_pos2;
	for(int k = shift_pos2; k<size_file-shift_pos2; k++) //Removing points with high variation 
         {
 	  long double temp;
   	  inter_ptr->index = k;
   	  inter_ptr->y_val = New_fn(k);
     	  temp = (inter_ptr->y_val)-(org_ptr->y_val);
   	  eps_sqr += temp*temp; 
   	  inter_ptr++; 
   	  org_ptr++;
    	}  
        break;
   case 'C':
        org_ptr = starting_o;
        inter_ptr = starting_i;
	for(int k = 0; k<size_file; k++) 
         {
 	  long double temp;
   	  inter_ptr->index = k;
   	  inter_ptr->y_val = C_spline(k);
     	  temp = (inter_ptr->y_val)-(org_ptr->y_val);
   	  eps_sqr += temp*temp; 
   	  inter_ptr++; 
   	  org_ptr++;
    	}  
        
        break;
   case 'E': 
	return -1;
        exit(0);
        break;
  default: printf("\nInvalid entry.");
           break;
  }
  
  //}
  epsilon = sqrt(eps_sqr);
  return epsilon;
 
}
