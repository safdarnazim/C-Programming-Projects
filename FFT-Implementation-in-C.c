#include<stdio.h>
#include<math.h>
#include<complex.h>
#define pi 3.14159265359
#define N 4


int n_point_fft_dit(int,double complex *);
int n_point_ifft_dit(int,double complex *);
int bit_reverse(int,int);
int no_of_bits(int);
int n_point_fft_dif(int,double complex *);
int n_point_ifft_dif(int,double complex *);

int bit_reverse(int n,int bits)
 {
    int reversed = 0;
    int i;

    for (i = 0; i < bits;i++)
    {
        reversed=reversed|(((n >> i) & 1) << (bits - 1 - i));
    }
    return reversed;
}

int no_of_bits(int n)
{
   int count=0;
   while(n>0)
   {
       count=count+1;
       n=n/2;
   }
   return count-1;
}

int main()
{
    printf("%d-POINT FFT\n",N);
    double complex x[N]={0,1,2,3,4,5,6,7};
    printf("Input sequence for DIT\n\n");
    for(int i=0;i<N;i++)
        printf("%f+j%f\t",creal(x[i]),cimag(x[i]));
    printf("\n\n\n");
    double complex X[N];
    double complex Y[N];
    int n=no_of_bits(N);
    for(int i=0;i<N;i++)
    {
        int j=bit_reverse(i,n);
        X[i]=x[j];
    }
    printf("DFT using DIT:\n\n");
    n_point_fft_dit(N,X);
    for(int i=0;i<N;i++)
    {
        printf("%f+j%f\t",creal(X[i]),cimag(X[i]));
    }
    printf("\n\n\n");
    for(int i=0;i<N;i++)
    {
        int j=bit_reverse(i,n);
        Y[i]=X[j];
    }
    printf("Bit reversed DFT\n\n");
    for(int i=0;i<N;i++)
        {
            printf("%f+j%f\t",creal(Y[i]),cimag(Y[i]));
        }
    printf("\n\n\n");
    printf("IFFT using DIT\n\n");
    n_point_ifft_dit(N,Y);
    for(int i=0;i<N;i++)
    {
        Y[i]=Y[i]/N;
        printf("%f+j%f\t",creal(Y[i]),cimag(Y[i]));
    }
    printf("\n\n\n");

    printf("Input sequence for DIF\n\n");
    for(int i=0;i<N;i++)
    {
        X[i]=x[i];
    }
    for(int i=0;i<N;i++)
        printf("%f+j%f\t",creal(X[i]),cimag(X[i]));
    printf("\n\n\n");


    printf("DFT using DIF:\n\n");
    n_point_fft_dif(N,X);
    for(int i=0;i<N;i++)
    {
        int j=bit_reverse(i,n);
        Y[i]=X[j];
    }
    for(int i=0;i<N;i++)
    {
        printf("%f+j%f\t",creal(Y[i]),cimag(Y[i]));
    }
    printf("\n\n\n");

    printf("IFFT using DIF\n\n");
    n_point_ifft_dif(N,Y);
    for(int i=0;i<N;i++)
    {
        int j=bit_reverse(i,n);
        X[i]=Y[j];
    }
    for(int i=0;i<N;i++)
    {
        X[i]=X[i]/N;
        printf("%f+j%f\t",creal(X[i]),cimag(X[i]));
    }
    printf("\n\n\n");




}

int n_point_fft_dit(int n,double complex * X)
{
    double complex Y[n];
    int count=no_of_bits(n);
    while(count>0)
    {
        int n1=n/pow(2,count-1);
        double angle=(2*pi*(-1))/n1;
        double complex w=cexp(I*angle);
        for(int i=0;i<n;i=i+n1)
        {
            int k=0;
            for(int j=0;j<n1;j++)
            {
                if(j==n1/2)
                {
                    k=0;
                }
                if(j>=n1/2)
                {
                    *(Y+i+j)=*(X+i+k)-((*(X+i+k+(n1/2)))*(cpow(w,k)));
                    k=k+1;
                }
                else
                {
                    *(Y+i+j)=*(X+i+k)+((*(X+i+k+(n1/2)))*cpow(w,k));
                    k=k+1;
                }
            }
        }
        for(int i=0;i<n;i++)
            {
                *(X+i)=*(Y+i);
            }
        count=count-1;

    }
    return 0;
}

int n_point_ifft_dit(int n,double complex * X)
{
    double complex Y[n];
    int count=no_of_bits(n);
    while(count>0)
    {
        int n1=n/pow(2,count-1);
        double angle=(2*pi*(-1))/n1;
        double complex w=cexp(I*angle);
        for(int i=0;i<n;i=i+n1)
        {
            int k=0;
            for(int j=0;j<n1;j++)
            {
                if(j==n1/2)
                {
                    k=0;
                }
                if(j>=n1/2)
                {
                    *(Y+i+j)=*(X+i+k)-((*(X+i+k+(n1/2)))*(cpow(w,-k)));
                    k=k+1;
                }
                else
                {
                    *(Y+i+j)=*(X+i+k)+((*(X+i+k+(n1/2)))*cpow(w,-k));
                    k=k+1;
                }
            }
        }
        for(int i=0;i<n;i++)
            {
                *(X+i)=*(Y+i);
            }
        count=count-1;
    }
    return 0;
}

int n_point_fft_dif(int n,double complex * X)
{
    double complex Y[n];
    int count=no_of_bits(n);
    int looper=0;
    while(count>0)
    {
        int n1=n/pow(2,looper);
        double angle=(2*pi*(-1))/n1;
        double complex w=cexp(I*angle);
        for(int i=0;i<n;i=i+n1)
        {
            int k=0;
            for(int j=0;j<n1;j++)
            {
                if(j==n1/2)
                {
                    k=0;
                }
                if(j>=n1/2)
                {
                    *(Y+i+j)=(*(X+i+k)-(*(X+i+k+(n1/2))))*(cpow(w,k));
                    k=k+1;
                }
                else
                {
                    *(Y+i+j)=*(X+i+k)+((*(X+i+k+(n1/2))));
                    k=k+1;
                }
            }
        }
        for(int i=0;i<n;i++)
            {
                *(X+i)=*(Y+i);
            }
        count=count-1;
        looper++;
    }
    return 0;
}

int n_point_ifft_dif(int n,double complex * X)
{
    double complex Y[n];
    int count=no_of_bits(n);
    int looper=0;
    while(count>0)
    {
        int n1=n/pow(2,looper);
        double angle=(2*pi*(-1))/n1;
        double complex w=cexp(I*angle);
        for(int i=0;i<n;i=i+n1)
        {
            int k=0;
            for(int j=0;j<n1;j++)
            {
                if(j==n1/2)
                {
                    k=0;
                }
                if(j>=n1/2)
                {
                    *(Y+i+j)=(*(X+i+k)-(*(X+i+k+(n1/2))))*(cpow(w,-k));
                    k=k+1;
                }
                else
                {
                    *(Y+i+j)=*(X+i+k)+((*(X+i+k+(n1/2))));
                    k=k+1;
                }
            }
        }
        for(int i=0;i<n;i++)
            {
                *(X+i)=*(Y+i);
            }
        count=count-1;
        looper++;
    }
    return 0;
}
