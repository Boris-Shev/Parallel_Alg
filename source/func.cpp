#include "func.h"
#include "matrix.h"

int solve(double *a,double *b,double *x,int n,double *d, double *l,  ARGS * arg)
{
    int thread_num = arg->thread_num;
    int thread_count = arg->thread_count;
    int i,j,k;
    double tmp;
    int ret;
    double crutch;
    long double time1=get_time();
    long double time2=get_full_time();
    double *curr;
    double norm = 0, max = 0;
    double eps = 1.737274e-14;


    if(thread_num==0)
    {
        for(i = 0;i < n; i++)
        {
            norm = 0;
            for(j = 0;j < n;j++)
            {
                norm += fabs(a[i*n + j]);
            }
            if(max < norm)
                max = norm;
        }
    }
    ret = synchronize(arg);
    if (ret < 0)
        return -1;
    for(i = 0;i < n; ++i)
    {
        if(thread_num==0)
        {
            crutch = a[i * n + i];
            for(j = 0; j < i; ++j)
            {

                crutch -= l[i*n+j] * (x[j] = l[i * n + j]*d[j]);

            }
            //printf("crutch = %lf\n",crutch);
            if(fabs(crutch) <= eps*norm)
                arg->ret = -1;

            l[i*n + i] = tmp = sqrt(fabs(crutch));// r(i,i)

            //d(i,i)
            if(crutch > 0)
            {
                d[i] = 1.0;
            }
            else
            {
                d[i] = -1.0;
                tmp = -tmp;
            }
        }

        ret = synchronize(arg);
        if (arg[thread_num].ret < 0)
            return -1;
        //printf("Ret = %d\n",ret);
        tmp = l[i*n+i]*d[i];
        if (fabs(tmp) <= eps*norm)
            arg->ret = -1;
        ret = synchronize(arg);
        if (ret < 0)
            return -1;

        int start=(n-i-1)*thread_num;
    		start=start/thread_count+i+1;
    		int end=(n-i-1)*(thread_num+1);
    		end=end/thread_count+i+1;

        for(j = start; j < end; ++j)
        {
            crutch = a[j*n + i];
            for(k = 0; k < i; ++k)
            {
                //printf("X[%d] = %lf, %d\n",k,x[k],thread_num); --Первые nan это l
                //printf("L(%d) = %lf I = %d J =%d Num = %d\n",j*n+k,l[j * n + k],i,j,thread_num);
                //printf("L = %lf x = %lf num = %d\n",l[j * n + k],x[k],thread_num);
                crutch -= l[j * n + k] * x[k];
            }
            //printf("TMP = %lf\n",tmp);
            l[j*n + i] = crutch / tmp; // r(j,i) tmp - r(i,i)*d(i)
          //  printf("TMP = %lf\n",tmp);
            //printf("L[] = %lf j = %d i = %d num = %d\n",l[j*n + i],j,i,thread_num); //--- это падает раньше
        }
        ret = synchronize(arg);
        if (ret < 0)
            return -1;

    }
////////////////////////////////////////////////////////////////////
    if(thread_num==0){
    for (i = 0; i < n; ++i)
    {
        curr = l + i * n;
        tmp = b[i];
        for (j = 0; j < i; ++j)
            tmp -= x[j] * curr[j];
        x[i] = tmp / curr[i];
    }
    //print_vector(d,n,n);

    for (i = n - 1; i >= 0; --i)
    {
        curr = l + i;
        tmp = x[i];
        if (d[i] > 0)
        {
            for (j = n - 1; j > i; --j)
                tmp -= d[j] * curr[j * n];
            d[i] = tmp / curr[i * n];
        }
        else
        {
            for (j = n - 1; j > i; --j)
                tmp += d[j] * curr[j * n];
            d[i] = -tmp / curr[i * n];
        }
    }
    for (i = 0; i < n; ++i)
        x[i] = d[i];
    }
    ret = synchronize(arg);
    if (ret < 0)
        return -1;
    time1-=get_time();
    time2-=get_full_time();
    printf("Thread %d: Processor time: %Lf sec. Astronomic time: %Lf sec.\n", thread_num+1, -time1, -time2);


    return 0;

}
