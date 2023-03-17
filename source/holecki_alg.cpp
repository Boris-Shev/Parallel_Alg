#include "header.hpp"

using namespace std;

void* HolecAlgParallel(void* arg) {
  ARGS* args = (ARGS*)arg;

  int n = args->n;
  args->time_ast = currentTime();
  args->time_cpu = cpu_time();
  args->err = HolecAlg(args->n, args->id, args->total_threads,
                      args->A, args->b, args->x, args->R, args->R + n*n);
  args->time_ast = currentTime() - args->time_ast;
  args->time_cpu = cpu_time() - args->time_cpu;

  pthread_exit(nullptr);
}

int HolecAlg (int n, int id, int total_threads, double* A,
              double* b, double* x, double* R, double* d) {
  double sum;
  double eps = 1e-14;
  double norm = 1;
  double tmp;
  int left = (n*n / total_threads) * (id - 1);
  int right = (n*n / total_threads) * id;
  if (id == total_threads)
    right = n*n;

  for(int i = left; i < right; i++)
    R[i] = 0;

  synchronize(total_threads);
  for (int i = 0; i < n; i++) {
    x[i] = 0;
    d[i] = 0;
  }

  if(id == 1) {
    norm = 0;
    for(int i = 0; i < n; i++)
      for(int j = 0; j < n; j++)
        norm += fabs(A[i*n + j]);
  }
  synchronize(total_threads);
  ///////////////////////////#Построение Холецкого#////////////////////////////
  for (int i = 0; i < n; i++) {
    synchronize(total_threads);
    if(id == 1) {
      sum = A[i*n + i];
      for (int k = 0; k < i; k++) {
        sum -= R[k*n + i] * (x[k] =  R[k*n + i] * d[k]);
      }

      if (fabs(sum) < eps * norm)
        d[0] = 0;

      d[i] = (double(0) < sum) - (sum < double(0));
      R[i*n + i] = sqrt(fabs(sum));
    }
    synchronize(total_threads);
    if (fabs(d[0]) < eps)
      return -1;

    tmp = R[i*n + i] * d[i];
    left = ((n-i-1) / total_threads) * (id - 1) + (i + 1);
    right = ((n-i-1) / total_threads) * id + (i + 1);
    if (id == total_threads)
      right = n;
    for (int j = left; j < right; j++) {
      sum = A[i*n + j];
      for (int k = 0; k < i; k++) {
        sum -= x[k] * R[k*n + j];
      }
      R[i*n + j] = sum / tmp;
    }
  }
  // if (id == 1) {
  //   PrintMat(A,n,n);
  //   printf("\n" );
  //   PrintMat(R,n,n);
  //   printf("\n" );
  //   PrintMat(d,1,n);
  //   printf("\n" );
  // }
  synchronize(total_threads);

  //////////////////////////#Обращение матрицы R#//////////////////////////////
  left = (n*n / total_threads) * (id - 1);
  right = (n*n / total_threads) * id;
  if (id == total_threads)
    right = n*n;

  for(int i = left; i < right; i++)
    A[i] = 0;
  synchronize(total_threads);
  for (int i = 0; i < n; i++)
    A[i*n + i] = 1;




  for (int i = n - 1; i >= 0; i--) {
    tmp = R[i*n + i];
    if(id == 1) {
      for (int j = i; j < n; j++) {
        A[i*n + j] /= R[i*n + i];
      }
    }
    synchronize(total_threads);

    //NEW
    left = ((n-i) / total_threads) * (id - 1) + i;
    right = ((n-i) / total_threads) * id + i;
    if (id == total_threads)
      right = n;
    //printf("%d: %d, %d\n", id, left, right);

    for (int k = 0; k < i; k++) {
      tmp = R[k*n + i];
      for (int j = left; j < right; j++) {
        A[k*n + j] -= tmp * A[i*n + j];
      }
    }
    synchronize(total_threads);
    // //OLD
    // left = (i / total_threads) * (id - 1);
    // right = (i / total_threads) * id;
    // if (id == total_threads)
    //   right = i;
    //
    // for (int k = left; k < right; k++) {
    //   for (int j = i; j < n; j++) {
    //     RR[k*n + j] = RR[k*n + j] - R[k*n + i] * RR[i*n + j];
    //   }
    // }
    // synchronize(total_threads);
  }
  // if (id == 1) {
  //   printf("\n" );
  //   PrintMat(RR,n,n);
  //   printf("\n" );
  // }
  ////////////////////////#Перемножение матриц#///////////////////////////////
  synchronize(total_threads);
  left = (n / total_threads) * (id - 1);
  right = (n / total_threads) * id;
  if (id == total_threads)
    right = n;

  //// NEW
  // for (int i = left; i < right; i++) {
  //   double tmp1 = 0;
  //   for (int k = 0; k < n; k++) {
  //     tmp = RR[i*n + k] * d[k];
  //     for (int q = 0; q < n; q++) {
  //       tmp1 += tmp * RR[q*n + k] * b[q];
  //     }
  //   }
  //   x[i] = tmp1;
  // }

  // OLD
  for (int i = left; i < right; i++) {
    x[i] = 0;
    for (int k = 0; k < n; k++) {
      for (int q = 0; q < n; q++) {
        x[i] += A[i*n + k] * d[k] * A[q*n + k] * b[q];
      }
    }
  }

   return 0;
}


template <typename T>
int Sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
