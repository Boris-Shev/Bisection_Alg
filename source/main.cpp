#include "header.hpp"

int main(int argc, char* argv[]) {
  int n, m, k, err;
  double eps;
  err = TestInitArg(argc, argv, &n, &m, &k, &eps);
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
    err = InMat(n, k, matrx, argv[5]);
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

  //PrintMat(matrx, n, n, n);
///////////// Конец проверок

  double* x = new double[n];
  double* extra_mem = new double[n*n];
  err = 0;
  double time = (double)clock();
  err = BisecAlg (n, matrx, x, eps, extra_mem);
  time = (double)(clock() - time) / CLOCKS_PER_SEC;
  if (err == -1) {
    printf("Матрица вырождена или некорректна\n");
    delete[] matrx;
    delete[] x;
    delete[] extra_mem;
    return -8;
  }
  PrintMat(x, n, 1, m);
  printf("Время алгоритма: %.3lf\n", time);
  // printf("Норма невязки: %lf\n", Residual(matrx, n, b, x));
  // printf("Норма погрешности: %lf\n", Inaccuracy(x, n));

  delete[] matrx;
  delete[] x;
  delete[] extra_mem;
  return 0;
}
