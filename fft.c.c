// Please go through the READ- ME file attached with the code

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <tgmath.h>
#include <stdlib.h>

/* ---------------------- MACROS DESCRIPTION--------------------------

   The macros defined below are used to convert the Floating point data to the Fixed point data.
   Whenever required the macros are used to convert floating to fixed point and vice versa.
   Separate macros are used to perform multiplication and complex multiplication on Fixed Point data.
*/

#define pi 3.14
#define DoubleToFixed(x) ((x) * (double)(1<<scale))
#define FixedToDouble(x) ((double)(x) /(double)(1<<scale))
#define CmplxToFixed(x) ((DoubleToFixed(creal(x))+I*(DoubleToFixed(cimag(x)));
#define FixedToCmplx(x) ((FixedToDouble(creal(x))+I*(FixedToDouble(cimag(x)));
//#define IntToFixed(x) ((x)<<scale)
//#define FixedToInt(x) ((x)>>scale)
#define MUL (x,y) ((((x)>>8)*((y)>>8))>>0)
#define CMUL(x,y) ((MUL(creal(x),creal(y))-MUL(cimag(x),cimag(y)))+I*(MUL(creal(x),cimag(y))-MUL(cimag(x),creal(y))))

const int scale = 16;

typedef int complex cmplx;

int x[1][100];

/* ---------------- FUNCTION DESCRIPTION ------------------------

   This fuction is a C- program adaptation of Cooley - Tukey algortihm for the compuation of FFT of input data.
   The complex array of input values and the value of N which decides the accuracy of the output ( no.  of splittings) 
   is given as input to the function

*/


cmplx** fft(int x[1][100], int L, int N)
{
	if(N<L)
	{
		printf("N should be  greater than L");
		//break;
		
	}
	
	//to ensure that the subsequent splittings occur without any problem

	int xn[1][N];
	
	for(int i=0;i<L;i++)
	{
		xn[1][i] = x[i];
	}
	
    for (int i=L;i<N;i++)
	{
		xn[1][i]=0;
	}
    
	cmplx xo[1][N/2],transpose_xo[1][N/2];
	cmplx xe[1][N/2],transpose_xe[1][N/2];
	cmplx Xe[N][N/2] ;
	cmplx Xo[N][N/2];
	
	int j = 0;
    for(int i=0;i<N;i+=2)
    {
        xo[0][j]=xn[0][i+1]+I*0;
        xe[0][j]=xn[0][i]+I*0;
        j+=1;
    }
	
	//splitting the input matrix to sub parts even and odd(Cooley Tukey algorithm) 
	
	cmplx Wn[N/2];
	
	for(int k=0;k<(N-1);k++)
	{
		for(int r=0;r<((N/2)-1);r++)
		{
			int j=0;
			Wn[j]=CmplxtoFixed(cexp(-(4*I*pi*k*r/N)));
			Xe[k+1][r+1]=Wn[j];
			Xo[k+1][r+1]=Wn[j];
			j=j+1;
		}
	}
	
	cmplx Wn1[N][1];
	
	for(int k=0;k<(N-1);k++)
	{
		Wn1[k][0] = CmplxtoFixed(cexp(-(I*2*pi*(k/N))));
	}
	
	static cmplx p[N][1];
    for(int i=0; i<1; ++i)
	{
		for(int j=0; j<N/2; ++j)
		{
            transpose_xe[j][i] = xe[i][j];
			transpose_xo[j][i] = xo[i][j];
        }
	}
    
	cmplx m1[N][1],m2[N][1],m3[N][1];
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < 1; j++)
		{
            double tmp1 = 0.0 ,tmp2 =0.0;

            for (int k = 0; k < N/2; k++)
            {
                tmp1 += CMUL(Xe[i][k],transpose_xe[k][j]);
                tmp2 += CMUL([i][k],transpose_xo[k][j]);
            }

             m1[i][j]=tmp1;
             m2[i][j]=tmp2;
        }
		
		for (int i = 0; i < N; i++)
		{
			m3[i][0]= CMUL(Wn1[i][0],m2[i][0]);
		}
		
		for(int i=0;i<N;i++)
		{
			p[i][0]= CMUL(m1[i][0],m3[i][0]);
		}
		
	}

return p;

}


/* --------------------------- MAIN PROGRAM DESCRIPTION -------------------------
   
   In the main program similar to the MATLAB program conisders double tone sinusoidal siganal . We are computing the Fourier transform 
   of the signal using the function defined above. In the output we compute the maximum value in the FFT array returned. This helps us 
   to know the maximum and second maximum frequency component values in the output for further analysis

*/

int main()
{
	
	int f1 =  100;
	int f2 =  400;
	int N;
	printf("Enter N value");
	scanf("%d\n",N);
	
	double T = 1/3000;
	
	int len = 0;
	
	for(int i=0;i<=(1-T);i+=T)
	{
		x[0][i] = DoubleToFixed(sin(2*pi*f1*i)+sin(2*pi*f2*i));
		len++;
	}
    
	cmplx z[N][1];
	
	cmplx** y;
    y = fft(x,len,N);
	y = &z;
	
	for(int i=0;i<N;i++)
	{
		z[i][0] = cabs(FixedToCmplx(z[i][0]));
	}
	
	double max1, max2 ;
	
	max1 = z[0][0];
    for (int i = 0; i < N; i++)
	{
        if (z[i][0]>max1) 
		{
            max1 = z[i][0];
        }
    }
	
    max2 = z[0][0];
    for (int i = 1; i < N; i++) 
	{
        if (z[i][0] > max2 && z[i][0] < max1)
            max2 = z[i][0];
    }
	
	printf("The frequencies in the input signal are %lf and %lf ",max1,max2);
	
}