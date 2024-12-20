#include "pwsim.h"


/* Function to generate additive white Gaussian Noise samples with zero mean and a standard deviation of 1. 
 * This code was developed by "Cagri Tanriover" and published in https://www.embeddedrelated.com/
*/
double AWGN_generator()
{
	double temp1;
	double temp2;
	double result;
	
	while(1)
	{
		temp2 = ( rand() / ( (double)RAND_MAX ) );
		if ( temp2 != 0 )
			break;
	}
	
	temp1 = cos( ( 2.0 * M_PI ) * rand() / ( (double)RAND_MAX ) );
	result = sqrt( -2.0 * log( temp2 ) ) * temp1;
	
	return result;
  
}


/* Function to remove linear trends.
 * It receives an array with its size and return a pointer to a new 
 * array with the linear trend removed.
 * */
double *detrend(double a[],int a_size){
	double x_m = (a_size-1.)/2.;
	double y_m = 0.0;
	int i;
	for(i=0;i<a_size;i++)
		y_m = y_m + a[i]/a_size;
	
	double num = 0.0; 
	double den = 0.0;
	for(i=0;i<a_size;i++)
	{
		num = num + (i-x_m)*(a[i]-y_m);
		den = den + (i-x_m)*(i-x_m);
	}
	
	double m = num/den;
	double b = y_m - m*x_m;
	
	double *det=(double*)malloc(a_size*sizeof(double));
	
	for(i=0;i<a_size;i++)
		det[i]=a[i]-(m*i+b);
	
	return det;
}



/* Receives an integer and return the higher next pow of 2.
 * */
int nextpow2(int x)
{
	return (int)pow(2, ceil(log(x)/log(2)));
}




/*////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
 * A standard implementation of mergesort algorithm and iterative
 * binary search. The three function was mainly taken from 
 * geeksforgeeks.org
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////*/

void merge(double arr[], int l, int m, int r) 
{
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 =  r - m; 
    double L[n1], R[n2]; 

    for (i = 0; i < n1; i++) 
        L[i] = arr[l + i]; 
    for (j = 0; j < n2; j++) 
        R[j] = arr[m + 1+ j]; 
  
    i = 0; 
    j = 0; 
    k = l; 
    while (i < n1 && j < n2)
    { 
		
        if (L[i] <= R[j])
        { 
            arr[k] = L[i]; 
            i++; 
        } 
        else
        { 
            arr[k] = R[j]; 
            j++; 
        } 
        
        k++; 
    }
  
    while (i < n1)
    { 
        arr[k] = L[i]; 
        i++; 
        k++; 
    } 

    while (j < n2)
    { 
        arr[k] = R[j]; 
        j++; 
        k++; 
    }
    
    return;
}
  
void mergeSort(double arr[], int l, int r)
{ 
    if (l < r)
    {
        int m = l+(r-l)/2; 
        mergeSort(arr, l, m); 
        mergeSort(arr, m+1, r); 
        merge(arr, l, m, r); 
    } 
    
    return;
} 

int binarySearch(double arr[], int l, int r, double x) 
{
  while (l <= r) 
  {
    int m = l + (r-l)/2; 

    if (arr[m] == x)  
        return m;   

    if (arr[m] < x)  
        l = m + 1;  

    else 
         r = m - 1;  
  } 

  return -1;  
} 
/*/////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////*/

/* Function that moves acceleration waveform in a time t0 adding padding zeroes.
 * Receives a time t0, an acceleration array, the sample frequency, the size
 * of the acceleration array, and the final size of the array with zeroes.
 * */
double *offset(double t0,double a[],double Fs,int a_size,int final_size)
{
	
	double *S = (double*)calloc(final_size,sizeof(double)); 
	
	int ind_offset = (int)(t0/(1./Fs));
	
	for(int i=0;i<final_size;i++)
	{
		if(i<ind_offset || i>=ind_offset+a_size)
			S[i]=0.0;
		else
			S[i]=a[i-ind_offset];
	}
	
	return S;
	
}


/* Smooth a single vector value depending on its index, using a moving average.
 * The index of the value and the order determinates the amount of neighbors
 * to use.
 * */
double complex smooth(double complex *A,int A_size, int ind, double complex *cum_sum, int a_order)
{
	int order = (a_order%2 ==0)? a_order-1:a_order;
	
	int aux;
	int j;
	double complex aux_sum=0.0;
	if(ind<(int)((order-1)/2) || (A_size-1)-ind<(int)((order-1)/2) )
	{
		
		aux = ind<(A_size-1)-ind? ind:(A_size-1)-ind;
		for(j=ind-aux;j<=ind+aux;j++)
		{
			aux_sum = aux_sum + A[j]/(2.*aux+1.);
		}
		
		return aux_sum;
	}
	
	else
	{
		if(*cum_sum==0.0)
		{
			for(j=ind-(int)((order-1)/2);j<=ind+(int)((order-1)/2);j++)
			{
				aux_sum = aux_sum + A[j]/order;
			}
			*cum_sum=aux_sum;
			return aux_sum;
		}
		else
		{
			aux_sum = *cum_sum - A[ind-(int)((order-1)/2)-1]/order;
			aux_sum =  aux_sum + A[ind+(int)((order-1)/2)]/order;
			*cum_sum = aux_sum;
			return aux_sum;
		}
	}
	
}


/* interpolate transfers functions amplitudes provided manually by user
 * receives the filename with frequencies, and horizontal and vertical
 * amplifications coefficients.
 * */
double complex ** interpolateFA(char *filename)
{
	
	int F=countLines(filename);
    double* ff =(double*)malloc((F+2)*sizeof(double));
	double* Ah =(double*)malloc((F+2)*sizeof(double));
	double* Av =(double*)malloc((F+2)*sizeof(double));

    
    
    FILE *archivo = fopen(filename,"rt"); 
    if(archivo == NULL)
    {
		printf("Frequency-amplifications file %s cannot be opened\n", filename);
		exit(EXIT_FAILURE);
	}
	int i;
    int C=1; 
    int lee;
    double A1,A2,A3;
    do
    {
		fscanf(archivo,"%lf  %lf  %lf", &A1,&A2,&A3);
        lee=feof(archivo);
        if (lee==1) break;  //stops if cannot read
        if(A1>f[nnyq-1])
        {
			printf("Defined frequencies exceeds nyquist frequency\n");
			exit(EXIT_FAILURE);
		}
        ff[C]=A1;
		Ah[C]=A2;
		Av[C]=A3;

		
        C++;

    }while(1);
    fclose(archivo);
    
    
    ff[0]=-0.001;
	Ah[0]=Ah[1];
	Av[0]=Av[1];
	ff[F+1]=f[nnyq-1]+0.001;
	Ah[F+1]=Ah[F];
	Av[F+1]=Av[F];
    
    double complex **Famps = (double complex**)malloc(2*sizeof(double complex*));
	
	for(i=0;i<2;i++)
		Famps[i]  = (double complex *)malloc(nnyq*sizeof(double complex));
		
	C=0;
	
	for(i=0;i<nnyq;i++)
	{
		if( !(f[i]>=ff[C] && f[i]<=ff[C+1])	)
		{
			C++; //  :O !!!!!
		}


		Famps[0][i] = Ah[C] + ((Ah[C+1]-Ah[C])/(ff[C+1]-ff[C])) * (f[i] - ff[C]);
		Famps[1][i] = Av[C] + ((Av[C+1]-Av[C])/(ff[C+1]-ff[C])) * (f[i] - ff[C]);
	
	}
	
	return Famps;

}

