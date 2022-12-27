#include "header.hpp"

using namespace std;

void* HolecAlgParallel(void* arg) {
  ARGS* args = (ARGS*)arg;

  HolecAlg(args->n, args->id, args->total_threads,
                      args->A, args->b, args->x, args->ExtraMem);

  return nullptr;
}

int HolecAlg (int n, int id, int total_threads, double* A,
              double* b, double* x, double* R) {
id=id;total_threads=total_threads;

  double sum = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      sum = 0;
      for (int k = 0; k < j; k++){                  // Распараллелить хз как
        sum += R[i*n + k] * R[j*n + k];
      }
      if (i == j)
        R[i*n + j] = sqrt(fabs(A[i*n + i] - sum));

      else {
        if (fabs(R[j*n + j]) < std::numeric_limits<double>::epsilon())
          return -1; // Матрица вырождена

        R[i*n + j] = (1.0 / R[j*n + j] * (A[i*n + j] - sum));
      }
    }
    // PrintMat(R,n,n);
    // printf("\n");
  }
   // PrintMat(R,n,n);
   // printf("\n");

  for (int k = 0; k < n; k++) {
    for (int i = 0; i < k; i++)
      R[k*n + i] /= R[k*n + k];

    R[k*n + k] = 1 / R[k*n + k];

    for(int i = k + 1; i < n; i++)
      for (int j = 0; j < k; j++)                 // Распараллелить можно
        R[i*n + j] += - R[k*n + j] * R[i*n + k];

    for (int i = k + 1; i < n; i++)
      R[i*n + k] = - R[k*n + k] * R[i*n + k];
  }
  // PrintMat(R, n, n, n);
  // printf("\n");

  for (int m = 0; m < n; m++) {
    double sum2 = 0;
    for (int z = 0; z < m; z++) {
      sum = 0;
      for (int k = m; k < n; k++) {               // Распараллелить
        sum += R[k*n + z] * R[k*n + m];
      }
      sum2 += sum * b[z];
    }
    for (int z = m; z < n; z++) {
      sum = 0;
      for (int k = z; k < n; k++) {                 // Распараллелить
        sum += R[k*n + z] * R[k*n + m];
      }
      sum2 += sum * b[z];
    }
    x[m] = sum2;
  }
  // PrintMat(x, n, 1, n);
  // printf("\n");
  return 0;
}


template <typename T>
int Sgn(T val) {
  return (T(0) < val) - (val < T(0));
}
