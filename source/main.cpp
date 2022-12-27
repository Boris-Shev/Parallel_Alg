#include "header.hpp"

int main(int argc, char* argv[]) {
  int n, m, total_threads, k, err;
  err = TestInitArg(argc, argv, &n, &m, &total_threads, &k);
  switch (err) {
    case -1:
      printf("Недостаточное количество аргументов\n");
      return -1;

    case -2:
      printf("Некорректные аргументы\n");
      return -2;

    default:
      break;
  }

  if (n <= 0) {
    printf("Некорректный размер матрицы\n");
    return -3;
  }
  double* matrx = new double[n*n];
  if (k == 0)
    err = InMat(n, k, matrx, argv[4]);
  else
    err = InMat(n, k, matrx, nullptr);
  switch (err) {
    case -1:
      std::cout << "Ошибка открытия файла" << std::endl;
      delete[] matrx;
      return -4;
    case -2:
      std::cout << "Неверный формат данных" << std::endl;
      delete[] matrx;
      return -5;
    case -3:
      std::cout << "Ошибка ввода-вывода при чтении" << std::endl;
      delete[] matrx;
      return -6;
    case -4:
      std::cout << "Недостаточное количество элементов" << std::endl;
      delete[] matrx;
      return -7;
    default:
      break;
  }

  double* b = new double[n];
  double* x = new double[n];
  for (int i = 0; i < n; i++) {
    double sum = 0;
    for (int j = 0; j < (n + 1)/ 2; j++)
      sum += matrx[i*n + 2*j];
    b[i] = sum;
  }

  double* extra_mem = new double[n*n];
  ARGS* args = new ARGS[total_threads];
  pthread_t* threads = new pthread_t[total_threads];
  for (int i = 0; i < total_threads; i++)
	{
		args[i].n = n;
		args[i].A = matrx;
    args[i].b = b;
		args[i].x = x;
		args[i].id = i;
    args[i].ExtraMem = extra_mem;
		args[i].total_threads = total_threads;
	}

  err = 0;
  double time = (double)clock();
  ////////////////////////////////////////////////////////
  for (int i = 0; i < total_threads; i++)
		if (pthread_create(threads + i, nullptr, HolecAlgParallel, args + i)) {
			printf("Не получилось создать поток %d!\n", i);
      delete[] matrx;
      delete[] b;
      delete[] x;
      delete[] extra_mem;
      delete[] args;
      delete[] threads;

			return -10;
		}

	for (int i = 0; i < total_threads; i++)
		if (pthread_join(threads[i], nullptr))
		{
			printf("Не получилось подождать поток %d!\n", i);
      delete[] matrx;
      delete[] b;
      delete[] x;
      delete[] extra_mem;
      delete[] args;
      delete[] threads;
			return -11;
		}
    ////////////////////////////////////////////////////
  time = (double)(clock() - time) / CLOCKS_PER_SEC;
  if (err == -1) {
    printf("Матрица вырождена или некорректна\n");
    delete[] matrx;
    delete[] b;
    delete[] x;
    delete[] extra_mem;
    delete[] args;
    delete[] threads;
    return -8;
  }

  PrintMat(x, n, 1, m);
  printf("Время алгоритма: %.3lf\n", time);
  printf("Норма невязки: %lf\n", Residual(matrx, n, b, x));
  printf("Норма погрешности: %lf\n", Inaccuracy(x, n));

  delete[] matrx;
  delete[] b;
  delete[] x;
  delete[] extra_mem;
  delete[] args;
  delete[] threads;
  return 0;
}
